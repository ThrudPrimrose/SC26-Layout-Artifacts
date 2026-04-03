/*
 * Layout-conflict benchmark (CPU) — NUMA-aware, per-thread PMC
 * =============================================================
 * C[i,j] += A[i,j] + B[i,j]
 *
 * Per-thread perf_event_open: each OMP thread opens its own counters
 * in the initial parallel region.  Main thread enable/disable/reads
 * all fds around the kernel call.  Sum across threads = true total.
 *
 * Compile: g++ -O3 -march=native -mtune=native -fopenmp -ffast-math
 *          -fvect-cost-model=cheap -std=c++17 -o bench_cpu bench_cpu.cpp -lnuma
 * Run:     ./bench_cpu <csv_file> [nthreads]
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cerrno>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>

#include <omp.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <linux/perf_event.h>

#if __has_include(<numaif.h>)
  #include <numaif.h>
  #define HAS_NUMA 1
#else
  #define HAS_NUMA 0
  #define MPOL_BIND 2
  static int mbind(void*,size_t,int,const unsigned long*,unsigned long,unsigned){return 0;}
#endif

#ifndef M_DIM
#define M_DIM 16384
#endif
#ifndef N_DIM
#define N_DIM 16384
#endif
static const int M = M_DIM;
static const int N = N_DIM;
#define NREP    100
#define NWARMUP 5
#define MAX_THREADS 512


/* ════════════════════════════════════════════════════════════════════
 *  PER-THREAD PMC
 *
 *  Each OMP thread calls perf_event_open(pid=0) to create counters
 *  that track THAT thread only.  Fds stored in g_pmc_fds[tid][evt].
 *  Main thread can ioctl/read any fd (fds are process-wide).
 *
 *  Flow:
 *    pmc_open_all()    — called inside #pragma omp parallel at startup
 *    pmc_enable_all()  — main thread, before kernel
 *    fn(A, B, C)       — kernel (own #pragma omp parallel for)
 *    pmc_read_all()    — main thread, after kernel; returns sum
 * ════════════════════════════════════════════════════════════════════ */

static long sys_perf_event_open(struct perf_event_attr *a, pid_t p,
                                int cpu, int grp, unsigned long fl) {
    return syscall(__NR_perf_event_open, a, p, cpu, grp, fl);
}

struct PmcDef { const char *name; uint32_t type; uint64_t config; };

static const PmcDef PMC_EVENTS[] = {
    {"cycles",        PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
    {"instructions",  PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS},
    {"cache-refs",    PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_REFERENCES},
    {"cache-misses",  PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES},
    {"L1d-load-miss", PERF_TYPE_HW_CACHE,
        (PERF_COUNT_HW_CACHE_L1D) |
        ((uint64_t)PERF_COUNT_HW_CACHE_OP_READ << 8) |
        ((uint64_t)PERF_COUNT_HW_CACHE_RESULT_MISS << 16)},
};
static constexpr int N_PMC = sizeof(PMC_EVENTS) / sizeof(PMC_EVENTS[0]);

static int  g_pmc_fds[MAX_THREADS][N_PMC];
static bool g_pmc_avail[N_PMC];    /* per-event availability */
static bool g_pmc_active = false;
static int  g_nthreads = 0;

/* Called once inside #pragma omp parallel at startup */
static void pmc_open_all() {
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nth = omp_get_num_threads();

        #pragma omp master
        {
            g_nthreads = nth;
            for (int i = 0; i < N_PMC; i++) g_pmc_avail[i] = true;
        }
        #pragma omp barrier

        for (int i = 0; i < N_PMC; i++) {
            struct perf_event_attr pe = {};
            pe.size = sizeof(pe);
            pe.type = PMC_EVENTS[i].type;
            pe.config = PMC_EVENTS[i].config;
            pe.disabled = 1;
            pe.exclude_kernel = 1;
            pe.exclude_hv = 1;
            /* pid=0 = THIS thread, cpu=-1 = any cpu */
            int fd = (int)sys_perf_event_open(&pe, 0, -1, -1, 0);
            g_pmc_fds[tid][i] = fd;
            if (fd < 0) {
                #pragma omp critical
                {
                    if (g_pmc_avail[i]) {
                        fprintf(stderr, "  [PMC] '%s' unavailable (errno=%d) -- skipped\n",
                                PMC_EVENTS[i].name, errno);
                        g_pmc_avail[i] = false;
                    }
                }
            }
        }
    }

    /* Check if at least one event works */
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) { g_pmc_active = true; break; }

    if (g_pmc_active) {
        int ncounters = 0;
        for (int i = 0; i < N_PMC; i++) if (g_pmc_avail[i]) ncounters++;
        printf("  PMC: ACTIVE — %d events x %d threads = %d counters\n",
               ncounters, g_nthreads, ncounters * g_nthreads);
    } else {
        printf("  PMC: DISABLED (no events available)\n");
    }
}

static void pmc_close_all() {
    for (int t = 0; t < g_nthreads; t++)
        for (int i = 0; i < N_PMC; i++)
            if (g_pmc_fds[t][i] >= 0) close(g_pmc_fds[t][i]);
}

/* Main thread calls: reset + enable all fds across all threads */
static void pmc_enable_all() {
    if (!g_pmc_active) return;
    for (int t = 0; t < g_nthreads; t++)
        for (int i = 0; i < N_PMC; i++)
            if (g_pmc_fds[t][i] >= 0) {
                ioctl(g_pmc_fds[t][i], PERF_EVENT_IOC_RESET, 0);
                ioctl(g_pmc_fds[t][i], PERF_EVENT_IOC_ENABLE, 0);
            }
}

/* Main thread calls: disable all, read all, sum per event */
static void pmc_read_all(uint64_t out[N_PMC]) {
    memset(out, 0, N_PMC * sizeof(uint64_t));
    if (!g_pmc_active) return;
    for (int t = 0; t < g_nthreads; t++)
        for (int i = 0; i < N_PMC; i++) {
            int fd = g_pmc_fds[t][i];
            if (fd < 0) continue;
            ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
            uint64_t val = 0;
            if (read(fd, &val, sizeof(val)) == sizeof(val))
                out[i] += val;
        }
}

/* ── helpers ──────────────────────────────────────────────────────── */
static inline int idx_rm(int i, int j) { return i * N + j; }
static inline int idx_cm(int i, int j) { return j * M + i; }
template<int SB> static inline int idx_blk_rm(int i, int j) {
    return (i/SB*(N/SB)+j/SB)*SB*SB+(i%SB)*SB+(j%SB); }
template<int SB> static inline int idx_blk_cm(int i, int j) {
    return (i/SB*(N/SB)+j/SB)*SB*SB+(j%SB)*SB+(i%SB); }

using Clock = std::chrono::high_resolution_clock;
static inline double esec(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b-a).count(); }

/* ════════════════════════════════════════════════════════════════════
 *  NUMA alloc + first-touch
 * ════════════════════════════════════════════════════════════════════ */
static long PAGE_SZ;
static int get_numa_node(){unsigned c,n;syscall(__NR_getcpu,&c,&n,nullptr);return(int)n;}

static void *numa_alloc(size_t bytes) {
    void *p = mmap(nullptr,bytes,PROT_READ|PROT_WRITE,
                   MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE,-1,0);
    if (p==MAP_FAILED){perror("mmap");std::abort();}
#if NOHUGEPAGE
    madvise(p,bytes,MADV_NOHUGEPAGE);
#else
    madvise(p,bytes,MADV_HUGEPAGE);
#endif
    return p;
}
static void numa_free(void *p, size_t b) { munmap(p,b); }

static void bind_and_touch(void *base, size_t total_bytes) {
    int64_t np=(total_bytes+PAGE_SZ-1)/PAGE_SZ;
    #pragma omp parallel
    {
        int tid=omp_get_thread_num(),nt=omp_get_num_threads();
        int64_t ch=(np+nt-1)/nt,lo=tid*ch,hi=std::min(lo+ch,np);
        if(lo<np&&lo<hi){
            uintptr_t pb=(uintptr_t)base+lo*PAGE_SZ;
            uintptr_t pe=(uintptr_t)base+hi*PAGE_SZ;
            if(pe>(uintptr_t)base+total_bytes) pe=(uintptr_t)base+total_bytes;
            int nd=get_numa_node();
            unsigned long mask[16]={};
            mask[nd/(8*sizeof(unsigned long))]=1UL<<(nd%(8*sizeof(unsigned long)));
            mbind((void*)pb,pe-pb,MPOL_BIND,mask,sizeof(mask)*8+1,0);
            for(uintptr_t a=pb;a<pe;a+=PAGE_SZ) *(volatile char*)a=0;
        }
    }
}

static void ft_init_rm(double *b,size_t tot,int M_,int N_,bool isB){
    bind_and_touch(b,tot*sizeof(double));
    #pragma omp parallel for schedule(static)
    for(size_t k=0;k<tot;k++){int i=(int)(k/N_),j=(int)(k%N_);
        b[k]=isB?(double)((i*13+j*37)%1000)/100.0:(double)((i*17+j*31)%1000)/100.0;}}

static void ft_init_cm(double *b,size_t tot,int M_,int N_){
    bind_and_touch(b,tot*sizeof(double));
    #pragma omp parallel for schedule(static)
    for(int j=0;j<N_;j++) for(int i=0;i<M_;i++)
        b[j*M_+i]=(double)((i*13+j*37)%1000)/100.0;}

template<int SB> static void ft_init_blk_rm(double *buf,size_t tot,int M_,int N_,bool isB){
    bind_and_touch(buf,tot*sizeof(double)); int NB=N_/SB;
    #pragma omp parallel for schedule(static) collapse(2)
    for(int br=0;br<M_/SB;br++) for(int bc=0;bc<NB;bc++){
        double *bl=buf+(br*NB+bc)*SB*SB;
        for(int lr=0;lr<SB;lr++) for(int lc=0;lc<SB;lc++){
            int i=br*SB+lr,j=bc*SB+lc;
            bl[lr*SB+lc]=isB?(double)((i*13+j*37)%1000)/100.0
                             :(double)((i*17+j*31)%1000)/100.0;}}}

template<int SB> static void ft_init_blk_cm(double *buf,size_t tot,int M_,int N_){
    bind_and_touch(buf,tot*sizeof(double)); int NB=N_/SB;
    #pragma omp parallel for schedule(static) collapse(2)
    for(int br=0;br<M_/SB;br++) for(int bc=0;bc<NB;bc++){
        double *bl=buf+(br*NB+bc)*SB*SB;
        for(int lr=0;lr<SB;lr++) for(int lc=0;lc<SB;lc++){
            int i=br*SB+lr,j=bc*SB+lc;
            bl[lc*SB+lr]=(double)((i*13+j*37)%1000)/100.0;}}}

static constexpr int64_t FLUSH_N=1<<27;
static double *g_flush_buf=nullptr;
static void cache_flush(){
    if(!g_flush_buf){g_flush_buf=(double*)numa_alloc(FLUSH_N*sizeof(double));
        bind_and_touch(g_flush_buf,FLUSH_N*sizeof(double));}
    #pragma omp parallel for schedule(static)
    for(int64_t i=0;i<FLUSH_N;i++) g_flush_buf[i]=(double)i;
    volatile double s=g_flush_buf[FLUSH_N-1];(void)s;}

/* ════════════════════════════════════════════════════════════════════
 *  Kernels
 * ════════════════════════════════════════════════════════════════════ */
static void kernel_row_major(const double*__restrict__ A,const double*__restrict__ B,double*__restrict__ C){
    #pragma omp parallel for schedule(static)
    for(int i=0;i<M;i++){
        #pragma omp simd 
        for(int j=0;j<N;j++) C[idx_rm(i,j)]+=A[idx_rm(i,j)]+B[idx_cm(i,j)];}}

static void kernel_col_major(const double*__restrict__ A,const double*__restrict__ B,double*__restrict__ C){
    #pragma omp parallel for schedule(static)
    for(int j=0;j<N;j++){
        #pragma omp simd
        for(int i=0;i<M;i++) C[idx_rm(i,j)]+=A[idx_rm(i,j)]+B[idx_cm(i,j)];}}

template<int T>
static void kernel_tiled(const double*__restrict__ A,const double*__restrict__ B,double*__restrict__ C){
    int nti=(M+T-1)/T,ntj=(N+T-1)/T;
    #pragma omp parallel for schedule(static) collapse(2)
    for(int ti=0;ti<nti;ti++) for(int tj=0;tj<ntj;tj++){
        int ii=ti*T,jj=tj*T,ie=std::min(ii+T,M),je=std::min(jj+T,N);
        double bl[T*T];
        for(int j=jj;j<je;j++){
            #pragma omp simd
            for(int i=ii;i<ie;i++) bl[(i-ii)*T+(j-jj)]=B[idx_cm(i,j)];}
        for(int i=ii;i<ie;i++){
            #pragma omp simd
            for(int j=jj;j<je;j++) C[idx_rm(i,j)]+=A[idx_rm(i,j)]+bl[(i-ii)*T+(j-jj)];}}}

static void kernel_all_rowmajor(const double*__restrict__ A,const double*__restrict__ Br,double*__restrict__ C){
    #pragma omp parallel for schedule(static)
    for(int i=0;i<M;i++){
        #pragma omp simd
        for(int j=0;j<N;j++) C[idx_rm(i,j)]+=A[idx_rm(i,j)]+Br[idx_rm(i,j)];}}

template<int SB>
static void kernel_blk_all_rm(const double*__restrict__ A,const double*__restrict__ B,double*__restrict__ C){
    int NBi=M/SB,NBj=N/SB,NB=NBj;
    #pragma omp parallel for schedule(static) collapse(2)
    for(int br=0;br<NBi;br++) for(int bc=0;bc<NBj;bc++){
        const double*a=A+(br*NB+bc)*SB*SB,*b=B+(br*NB+bc)*SB*SB;
        double*c=C+(br*NB+bc)*SB*SB;
        for(int lr=0;lr<SB;lr++){
            #pragma omp simd
            for(int lc=0;lc<SB;lc++) c[lr*SB+lc]+=a[lr*SB+lc]+b[lr*SB+lc];}}}

template<int SB>
static void kernel_blk_conflict(const double*__restrict__ A,const double*__restrict__ B,double*__restrict__ C){
    int NBi=M/SB,NBj=N/SB,NB=NBj;
    #pragma omp parallel for schedule(static) collapse(2)
    for(int br=0;br<NBi;br++) for(int bc=0;bc<NBj;bc++){
        const double*a=A+(br*NB+bc)*SB*SB,*b=B+(br*NB+bc)*SB*SB;
        double*c=C+(br*NB+bc)*SB*SB;
        for(int lr=0;lr<SB;lr++){
            #pragma omp simd
            for(int lc=0;lc<SB;lc++) c[lr*SB+lc]+=a[lr*SB+lc]+b[lc*SB+lr];}}}

/* ════════════════════════════════════════════════════════════════════
 *  Benchmark with per-thread PMC
 * ════════════════════════════════════════════════════════════════════ */
using KernelFn = void(*)(const double*,const double*,double*);
struct IterRecord { std::string v; int t,nth,r; double ts,bw,cs; const char *st; };
struct PmcRow { std::string v; int t; double med_gbs; uint64_t tot[N_PMC]; };

static std::vector<IterRecord> g_rec;
static std::vector<PmcRow> g_pmcr;

static void bench(const char *vn, int ts, KernelFn fn,
                  const double *A, const double *B, double *C,
                  double rcs, int nth) {
    size_t total=(size_t)M*N;
    double dbytes=(double)total*sizeof(double)*4.0;  /* C += A + B: 2R + 1RW */
    auto zC=[&](){
        #pragma omp parallel for schedule(static)
        for(size_t k=0;k<total;k++) C[k]=0.0;};

    /* warmup (PMC off) */
    for(int w=0;w<NWARMUP;w++){zC();cache_flush();fn(A,B,C);}

    /* correctness */
    zC(); fn(A,B,C);
    double cs=0;
    #pragma omp parallel for schedule(static) reduction(+:cs)
    for(size_t k=0;k<total;k++) cs+=C[k];
    bool ok=(std::fabs(cs-rcs)<=1e-3*std::fabs(rcs));
    const char *st=ok?"PASS":"FAIL";
    if(!ok) fprintf(stderr,"  [%s T=%d] MISMATCH: %.6e vs %.6e\n",vn,ts,cs,rcs);

    PmcRow pr; pr.v=vn; pr.t=ts; memset(pr.tot,0,sizeof(pr.tot));
    uint64_t cnt[N_PMC];

    for(int r=0;r<NREP;r++){
        //zC();
        cache_flush();

        pmc_enable_all();            /* reset+enable all per-thread fds */
        auto t0=Clock::now();
        fn(A,B,C);                   /* kernel — own #pragma omp parallel for */
        auto t1=Clock::now();
        pmc_read_all(cnt);           /* disable+read all, sum across threads */

        double dt=esec(t0,t1), bw=dbytes/dt/1e9;
        g_rec.push_back({vn,ts,nth,r,dt,bw,cs,st});
        for(int i=0;i<N_PMC;i++) pr.tot[i]+=cnt[i];
    }

    std::vector<double> bws;
    for(int i=(int)g_rec.size()-NREP;i<(int)g_rec.size();i++) bws.push_back(g_rec[i].bw);
    std::sort(bws.begin(),bws.end());
    pr.med_gbs=bws[NREP/2];
    g_pmcr.push_back(pr);

    printf("  %-24s T=%-4d  med %7.1f GB/s  [%7.1f..%7.1f]  %s\n",
           vn,ts,bws[NREP/2],bws[0],bws[NREP-1],st);
}

template<int SB> static void bench_blocked(int nth){
    static_assert(M_DIM%SB==0&&N_DIM%SB==0,"");
    size_t total=(size_t)M*N,bytes=total*sizeof(double);
    printf("\n  --- Blocked SB=%d ---\n",SB);
    double *Ab=(double*)numa_alloc(bytes),*Brm=(double*)numa_alloc(bytes);
    double *Bcm=(double*)numa_alloc(bytes),*Cb=(double*)numa_alloc(bytes);
    ft_init_blk_rm<SB>(Ab,total,M,N,false);ft_init_blk_rm<SB>(Brm,total,M,N,true);
    ft_init_blk_cm<SB>(Bcm,total,M,N);
    bind_and_touch(Cb,bytes);
    #pragma omp parallel for schedule(static)
    for(size_t k=0;k<total;k++) Cb[k]=0;
    kernel_blk_all_rm<SB>(Ab,Brm,Cb);
    double rc=0;
    #pragma omp parallel for schedule(static) reduction(+:rc)
    for(size_t k=0;k<total;k++) rc+=Cb[k];
    bench("blk_all_rm",SB,kernel_blk_all_rm<SB>,Ab,Brm,Cb,rc,nth);
    bench("blk_conflict",SB,kernel_blk_conflict<SB>,Ab,Bcm,Cb,rc,nth);
    numa_free(Ab,bytes);numa_free(Brm,bytes);numa_free(Bcm,bytes);numa_free(Cb,bytes);
}

/* ════════════════════════════════════════════════════════════════════
 *  PMC summary table
 * ════════════════════════════════════════════════════════════════════ */
static void print_pmc_table() {
    if(g_pmcr.empty()||!g_pmc_active) return;
    double elems=(double)M*(double)N*(double)NREP;

    printf("\n");
    printf("=============================================================================================================\n");
    printf("  PMC Summary (per-element = sum_all_threads / (M*N*NREP),  M=%d  N=%d  NREP=%d  threads=%d)\n",
           M,N,NREP,g_nthreads);
    printf("=============================================================================================================\n");
    printf("  %-20s %4s %7s","variant","T","GB/s");
    for(int i=0;i<N_PMC;i++) if(g_pmc_avail[i]) printf(" | %12s/el",PMC_EVENTS[i].name);
    printf("\n  %-20s %4s %7s","--------------------","----","-------");
    for(int i=0;i<N_PMC;i++) if(g_pmc_avail[i]) printf(" | %15s","---------------");
    printf("\n");

    for(auto &p : g_pmcr){
        char ts[16]; if(p.t>0) snprintf(ts,sizeof(ts),"%d",p.t); else snprintf(ts,sizeof(ts),"-");
        printf("  %-20s %4s %7.1f",p.v.c_str(),ts,p.med_gbs);
        for(int i=0;i<N_PMC;i++) if(g_pmc_avail[i])
            printf(" | %15.4f",(double)p.tot[i]/elems);
        printf("\n");
    }
    printf("=============================================================================================================\n");

    /* key comparisons */
    auto find=[&](const char *n,int t)->PmcRow*{
        for(auto &p:g_pmcr) if(p.v==n&&p.t==t) return &p; return nullptr;};

    PmcRow *best_ba=nullptr;
    for(auto &p:g_pmcr) if(p.v=="blk_all_rm"&&(!best_ba||p.med_gbs>best_ba->med_gbs)) best_ba=&p;
    auto *rm=find("all_rowmajor",0);
    auto *rw=find("row_major",0);

    auto print_ratios=[&](const char *title, PmcRow *a, PmcRow *b){
        if(!a||!b) return;
        printf("\n  %s vs %s(T=%d):  %.0f vs %.0f GB/s\n",
               a->v.c_str(),b->v.c_str(),b->t,a->med_gbs,b->med_gbs);
        for(int i=0;i<N_PMC;i++){
            if(!g_pmc_avail[i]||b->tot[i]==0) continue;
            double r=(double)a->tot[i]/(double)b->tot[i];
            const char *tag = r>1.1 ? "  <-- MORE" : (r<0.9 ? "  <-- LESS" : "");
            printf("    %-16s  %.2fx%s\n",PMC_EVENTS[i].name,r,tag);
        }
    };

    if(rm && best_ba)
        print_ratios("Cache set conflicts", rm, best_ba);
    if(rw && rm)
        print_ratios("Layout conflict cost", rw, rm);

    /* Sanity check: expected instructions per element */
    if(rm) {
        double ipe = 0;
        for(int i=0;i<N_PMC;i++)
            if(g_pmc_avail[i] && strcmp(PMC_EVENTS[i].name,"instructions")==0)
                ipe = (double)rm->tot[i] / elems;
        if(ipe > 0)
            printf("\n  Sanity: all_rowmajor instructions/elem = %.1f "
                   "(expect ~4-8 for vectorized add+store)\n", ipe);
    }
}

/* ════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv) {
    if(argc<2){fprintf(stderr,"Usage: %s <csv> [nthreads]\n",argv[0]);return 1;}
    const char *csv=argv[1];
    int rth=(argc>2)?atoi(argv[2]):0;
    if(rth>0) omp_set_num_threads(rth);
    PAGE_SZ=sysconf(_SC_PAGESIZE);

    /* Force thread pool creation + open per-thread PMC fds */
    int nth=0;
    #pragma omp parallel
    { nth=omp_get_num_threads(); }

    printf("Layout-conflict benchmark (CPU) -- per-thread PMC\n");
    printf("  M=%d  N=%d  reps=%d  warmup=%d  threads=%d  page=%ld\n",
           M,N,NREP,NWARMUP,nth,PAGE_SZ);
    printf("  NUMA=%s\n", HAS_NUMA?"yes":"fallback");

    /* Open PMC — must happen INSIDE a parallel region so each thread
       calls perf_event_open(pid=0) for itself */
    pmc_open_all();
    printf("\n");

    size_t total=(size_t)M*N, bytes=total*sizeof(double);
    double *A=(double*)numa_alloc(bytes),*Bcm=(double*)numa_alloc(bytes);
    double *Brm=(double*)numa_alloc(bytes),*C=(double*)numa_alloc(bytes);
    printf("Init arrays ...\n");
    ft_init_rm(A,total,M,N,false); ft_init_cm(Bcm,total,M,N);
    ft_init_rm(Brm,total,M,N,true);
    bind_and_touch(C,bytes);
    #pragma omp parallel for schedule(static)
    for(size_t k=0;k<total;k++) C[k]=0;

    /* ref checksums */
    #pragma omp parallel for schedule(static)
    for(size_t k=0;k<total;k++) C[k]=0;
    kernel_row_major(A,Bcm,C);
    double rcs=0;
    #pragma omp parallel for schedule(static) reduction(+:rcs)
    for(size_t k=0;k<total;k++) rcs+=C[k];

    #pragma omp parallel for schedule(static)
    for(size_t k=0;k<total;k++) C[k]=0;
    kernel_all_rowmajor(A,Brm,C);
    double rcs_rm=0;
    #pragma omp parallel for schedule(static) reduction(+:rcs_rm)
    for(size_t k=0;k<total;k++) rcs_rm+=C[k];

    printf("\n=== Naive schedules ===\n");
    bench("row_major",0,kernel_row_major,A,Bcm,C,rcs,nth);
    bench("col_major",0,kernel_col_major,A,Bcm,C,rcs,nth);
    bench("all_rowmajor",0,kernel_all_rowmajor,A,Brm,C,rcs_rm,nth);

    printf("\n=== Tiled ===\n");
    bench("tiled",8,kernel_tiled<8>,A,Bcm,C,rcs,nth);
    bench("tiled",16,kernel_tiled<16>,A,Bcm,C,rcs,nth);
    bench("tiled",32,kernel_tiled<32>,A,Bcm,C,rcs,nth);
    bench("tiled",64,kernel_tiled<64>,A,Bcm,C,rcs,nth);
    bench("tiled",128,kernel_tiled<128>,A,Bcm,C,rcs,nth);

    numa_free(A,bytes);numa_free(Bcm,bytes);numa_free(Brm,bytes);numa_free(C,bytes);

    printf("\n=== Blocked storage ===\n");
    bench_blocked<8>(nth);bench_blocked<16>(nth);bench_blocked<32>(nth);
    bench_blocked<64>(nth);bench_blocked<128>(nth);

    /* ── PMC table ── */
    print_pmc_table();

    /* ── CSV ── */
    FILE *fp=fopen(csv,"w");
    if(fp){
        fprintf(fp,"variant,M,N,tile,nthreads,rep,time_s,bw_gbs,checksum,status\n");
        for(auto &r:g_rec)
            fprintf(fp,"%s,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e,%s\n",
                    r.v.c_str(),M,N,r.t,r.nth,r.r,r.ts,r.bw,r.cs,r.st);
        fclose(fp);
        printf("\nWrote %zu records to %s\n",g_rec.size(),csv);
    }

    /* ── BW summary ── */
    auto med=[&](const char *n,int t)->double{
        std::vector<double> v;
        for(auto &r:g_rec) if(r.v==n&&r.t==t) v.push_back(r.bw);
        if(v.empty()) return 0; std::sort(v.begin(),v.end()); return v[v.size()/2];};

    double bt=0; for(int t:{8,16,32,64,128}) bt=std::max(bt,med("tiled",t));
    double br=med("row_major",0),bc=med("col_major",0),bctrl=med("all_rowmajor",0);
    double bba=0,bbc=0; int bas=0,bcs=0;
    for(int s:{8,16,32,64,128}){double v;
        v=med("blk_all_rm",s);if(v>bba){bba=v;bas=s;}
        v=med("blk_conflict",s);if(v>bbc){bbc=v;bcs=s;}}

    printf("\n=== Summary ===\n");
    printf("  Row-major:  tiled=%.0f  row=%.0f(%.1fx)  col=%.0f(%.1fx)  ctrl=%.0f GB/s\n",
           bt,br,bt/std::max(br,.1),bc,bt/std::max(bc,.1),bctrl);
    printf("  Blocked:    all_rm=%.0f(SB=%d)  conflict=%.0f(SB=%d)  ratio=%.2f\n",
           bba,bas,bbc,bcs,bbc/std::max(bba,.1));

    printf("\n=== Model (B_eff=8, C+=A+B -> 4 streams) ===\n");
    printf("  row_major:  cost=3+8=11   predicted %.1fx\n",(3+8.0)/4.0);
    printf("  col_major:  cost=1+3*8=25  predicted %.1fx\n",(1+24.0)/4.0);

    pmc_close_all();
    if(g_flush_buf) numa_free(g_flush_buf,FLUSH_N*sizeof(double));
    return 0;
}