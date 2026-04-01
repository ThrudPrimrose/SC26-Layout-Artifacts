/*
 * cost_metrics.cpp -- Cost metrics for z_v_grad_w stencil
 *
 * Four metrics only:
 *   μ              = avg new-block count per step
 *   μ·Δ            = avg (new-blocks × min-distance), unweighted
 *   μ·Δ_numa       = avg (new-blocks × NUMA-weighted distance)
 *   Δ_max          = avg block distance using max (worst-case reuse)
 *
 * Supports 5 cell distributions: uniform, normal_var1, normal_var4,
 * sequential, exact (from ICON serialised p_patch files).
 *
 * Environment:
 *   ICON_DATA_PATH  - directory with p_patch.*.data (default: primrose path)
 *   Timestep fixed to 9 for the connectivity load.
 *
 * Compile: g++ -O3 -std=c++17 -fopenmp -o cost_metrics cost_metrics.cpp
 * Run:     ./cost_metrics [N] [nlev] [beta] [alpha] [gamma] [P_NUMA]
 */
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <climits>
#include <cstring>
#include <unistd.h>

#include "icon_data_loader.h"

static int    BETA = 1;
static double ALPHA = 0.012;
static double GAMMA = 1.8;
static int    P_NUMA = 4;
static int    BLOCK_BYTES = 64;
#pragma omp threadprivate(BLOCK_BYTES)

/* ---- Layout index functions ---- */
inline int IC(int V, int je, int jk, int N, int nlev) {
    return (V <= 2) ? je + jk * N : jk + je * nlev;
}
inline int IN(int V, int je, int n, int N) {
    return (V == 1 || V == 3) ? je + n * N : n + je * 2;
}

/* ---- Arrays ---- */
enum ArrID {
    A_OUT=0, A_VN_IE, A_W, A_Z_VT_IE, A_Z_W_V,
    A_INV_DUAL, A_INV_PRIMAL, A_TANGENT,
    A_CELL_IDX, A_VERT_IDX, NUM_ARR
};
static const int ebytes[] = {8,8,8,8,8, 8,8,8, 4,4};

inline int64_t blk(int elem_idx, int arr) {
    return (int64_t)elem_idx * ebytes[arr] / BLOCK_BYTES;
}

/* ---- Schedule / loop order ---- */
enum Schedule  { SCHED_OMP_FOR=0, SCHED_OMP_COLLAPSE2=1 };
static const char* sched_name[] = {"omp_for","omp_collapse2"};
enum LoopOrder { KLON_FIRST=0, KLEV_FIRST=1 };

inline LoopOrder iteration_order(Schedule sched, int V) {
    if (sched == SCHED_OMP_COLLAPSE2) return KLON_FIRST;
    return (V <= 2) ? KLON_FIRST : KLEV_FIRST;
}

/* ---- NUMA placement ---- */
static Schedule G_SCHED;
#pragma omp threadprivate(G_SCHED)

inline int numa_2d(Schedule sched, int V, int elem, int N, int nlev) {
    int64_t total = (int64_t)N * nlev;
    if (sched == SCHED_OMP_COLLAPSE2) {
        int64_t lin;
        if (V <= 2) lin = elem;
        else { int jk = elem % nlev, je = elem / nlev; lin = (int64_t)jk * N + je; }
        int64_t chunk = (total + P_NUMA - 1) / P_NUMA;
        int d = (int)(lin / chunk);
        return (d < P_NUMA) ? d : P_NUMA - 1;
    }
    if (V <= 2) { int jk=elem/N; int chunk=(nlev+P_NUMA-1)/P_NUMA; int d=jk/chunk; return (d<P_NUMA)?d:P_NUMA-1; }
    else        { int je=elem/nlev; int chunk=(N+P_NUMA-1)/P_NUMA; int d=je/chunk; return (d<P_NUMA)?d:P_NUMA-1; }
}
inline int numa_1d(int elem, int N) {
    int chunk=(N+P_NUMA-1)/P_NUMA; if(chunk==0)return 0;
    int d=elem/chunk; return (d<P_NUMA)?d:P_NUMA-1;
}
inline int numa_idx(int V, int elem, int N) {
    int je=(V==1||V==3)?elem%N:elem/2; return numa_1d(je,N);
}
inline int nu_elem(int arr, int V, int elem, int N, int nlev) {
    switch(arr){
    case A_OUT:case A_VN_IE:case A_W:case A_Z_VT_IE:case A_Z_W_V:
        return numa_2d(G_SCHED,V,elem,N,nlev);
    case A_INV_DUAL:case A_INV_PRIMAL:case A_TANGENT: return numa_1d(elem,N);
    case A_CELL_IDX:case A_VERT_IDX: return numa_idx(V,elem,N);
    default: return 0;
    }
}
inline int nu_block(int arr, int V, int64_t baddr, int N, int nlev) {
    int elem=(int)(baddr*BLOCK_BYTES/ebytes[arr]); return nu_elem(arr,V,elem,N,nlev);
}

static int NU_HOME = 0;
#pragma omp threadprivate(NU_HOME)

inline double w_numa(int arr, int V, int64_t baddr, int64_t raw_dist, int N, int nlev) {
    if(nu_block(arr,V,baddr,N,nlev)!=NU_HOME) return GAMMA;
    return (raw_dist<BETA)?ALPHA:1.0;
}

/* ---- References per stencil iteration ---- */
static constexpr int NREFS = 14;
struct Ref { int arr; int64_t block; };

inline void collect_refs(int V, int je, int jk, int N, int nlev,
                         const int* cidx, const int* vidx, Ref refs[NREFS]) {
    int c2d=IC(V,je,jk,N,nlev);
    int ci0=cidx[IN(V,je,0,N)],ci1=cidx[IN(V,je,1,N)];
    int vi0=vidx[IN(V,je,0,N)],vi1=vidx[IN(V,je,1,N)];
    refs[ 0]={A_OUT,     blk(c2d,A_OUT)};
    refs[ 1]={A_VN_IE,   blk(c2d,A_VN_IE)};
    refs[ 2]={A_INV_DUAL,blk(je,A_INV_DUAL)};
    refs[ 3]={A_W,       blk(IC(V,ci0,jk,N,nlev),A_W)};
    refs[ 4]={A_W,       blk(IC(V,ci1,jk,N,nlev),A_W)};
    refs[ 5]={A_Z_VT_IE, blk(c2d,A_Z_VT_IE)};
    refs[ 6]={A_INV_PRIMAL,blk(je,A_INV_PRIMAL)};
    refs[ 7]={A_TANGENT, blk(je,A_TANGENT)};
    refs[ 8]={A_Z_W_V,   blk(IC(V,vi0,jk,N,nlev),A_Z_W_V)};
    refs[ 9]={A_Z_W_V,   blk(IC(V,vi1,jk,N,nlev),A_Z_W_V)};
    refs[10]={A_CELL_IDX, blk(IN(V,je,0,N),A_CELL_IDX)};
    refs[11]={A_CELL_IDX, blk(IN(V,je,1,N),A_CELL_IDX)};
    refs[12]={A_VERT_IDX, blk(IN(V,je,0,N),A_VERT_IDX)};
    refs[13]={A_VERT_IDX, blk(IN(V,je,1,N),A_VERT_IDX)};
}

/* ---- Block set ---- */
struct BlockSet {
    std::vector<int64_t> a[NUM_ARR];
    void clear(){for(int i=0;i<NUM_ARR;i++)a[i].clear();}
    void add(int arr,int64_t b){a[arr].push_back(b);}
    void finalize(){for(int i=0;i<NUM_ARR;i++){auto&v=a[i];std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());}}
    int total()const{int n=0;for(int i=0;i<NUM_ARR;i++)n+=(int)a[i].size();return n;}
    bool has(int arr,int64_t b)const{return std::binary_search(a[arr].begin(),a[arr].end(),b);}
};

/* ---- Metrics: only 4 ---- */
struct Metrics {
    double mu;             /* avg new-block count per step             */
    double mu_delta;       /* avg sum of (min-distance) over new blocks */
    double mu_delta_numa;  /* avg sum of (NUMA-weighted distance)       */
    double delta_max;      /* avg (max-distance) over new blocks        */
    int64_t T;
};

/* ---- Process one vectorised step ---- */
static inline void process_step(
    int V, int W_vec, int N, int nlev,
    LoopOrder loop, int je0, int jk0,
    const int* cidx, const int* vidx,
    BlockSet& prev, BlockSet& curr,
    int64_t& T,
    double& s_mu, double& s_md, double& s_md_numa, double& s_dmax)
{
    curr.clear();
    for(int w=0;w<W_vec;w++){
        int je=(loop==KLON_FIRST)?je0+w:je0;
        int jk=(loop==KLON_FIRST)?jk0:jk0+w;
        Ref refs[NREFS];
        collect_refs(V,je,jk,N,nlev,cidx,vidx,refs);
        for(int r=0;r<NREFS;r++) curr.add(refs[r].arr,refs[r].block);
    }
    curr.finalize();

    if(T==0){
        /* First step: everything is new, distance = 1 (sentinel) */
        int nb=curr.total();
        s_mu+=nb;
        s_md+=nb;              /* 1.0 per new block */
        s_md_numa+=nb;
        s_dmax+=1.0;
    } else {
        int nc=0;
        double sum_min=0, sum_dmax=0, sum_numa=0;

        for(int a=0;a<NUM_ARR;a++){
          for(int64_t b:curr.a[a]){
            if(prev.has(a,b)) continue;
            nc++;

            int64_t best_min=INT64_MAX, best_max=0;
            for(int64_t bp:prev.a[a]){
                int64_t d=std::abs(b-bp);
                if(d<best_min) best_min=d;
                if(d>best_max) best_max=d;
            }

            if(best_min==INT64_MAX){
                /* No previous block for this array */
                sum_min+=1.0; sum_dmax+=1.0;
                sum_numa+=w_numa(a,V,b,INT64_MAX,N,nlev);
            } else {
                double rmin=(double)best_min;
                double rmax=(double)best_max;
                double wn=w_numa(a,V,b,best_min,N,nlev);
                sum_min+=rmin;
                sum_dmax+=rmax;
                sum_numa+=wn*rmin;
            }
          }
        }
        s_mu+=nc;
        s_md+=sum_min;
        s_md_numa+=sum_numa;
        if(nc>0) s_dmax+=sum_dmax/nc;
    }
    T++;
    std::swap(prev,curr);
}

/* ---- Compute metrics over one NUMA slice ---- */
Metrics compute_slice(int V, int W_vec, int N, int nlev,
                      Schedule sched, const int* cidx, const int* vidx,
                      int64_t range_lo, int64_t range_hi) {
    G_SCHED=sched;
    LoopOrder loop=iteration_order(sched,V);
    BlockSet prev,curr;

    int64_t T=0;
    double s_mu=0, s_md=0, s_md_numa=0, s_dmax=0;

    if(sched==SCHED_OMP_COLLAPSE2){
        int64_t lin=range_lo;
        while(lin<range_hi){
            int jk=(int)(lin/N); int je_start=(int)(lin%N);
            int64_t row_end=std::min((int64_t)(jk+1)*N,range_hi);
            int je_end=(int)(row_end-(int64_t)jk*N);
            for(int je0=je_start;je0+W_vec<=je_end;je0+=W_vec)
                process_step(V,W_vec,N,nlev,KLON_FIRST,je0,jk,
                             cidx,vidx,prev,curr,T,
                             s_mu,s_md,s_md_numa,s_dmax);
            lin=row_end;
        }
    } else {
        int inner_n=(loop==KLON_FIRST)?N:nlev;
        for(int64_t outer=range_lo;outer<range_hi;outer++)
            for(int inner0=0;inner0+W_vec<=inner_n;inner0+=W_vec){
                int je=(loop==KLON_FIRST)?inner0:(int)outer;
                int jk=(loop==KLON_FIRST)?(int)outer:inner0;
                process_step(V,W_vec,N,nlev,loop,je,jk,
                             cidx,vidx,prev,curr,T,
                             s_mu,s_md,s_md_numa,s_dmax);
            }
    }
    if(T==0) return {};
    double dT=(double)T;
    return { s_mu/dT, s_md/dT, s_md_numa/dT, s_dmax/dT, T };
}

/* ---- Average over all P_NUMA domains ---- */
Metrics compute_metrics(int V, int W_vec, int N, int nlev,
                        Schedule sched, const int* cidx, const int* vidx) {
    LoopOrder loop=iteration_order(sched,V);
    int outer_n=(loop==KLON_FIRST)?nlev:N;
    Metrics acc={}; int n_domains=0;

    auto accumulate=[&](Metrics& m){
        if(m.T==0) return;
        acc.mu+=m.mu; acc.mu_delta+=m.mu_delta;
        acc.mu_delta_numa+=m.mu_delta_numa; acc.delta_max+=m.delta_max;
        acc.T+=m.T; n_domains++;
    };

    if(sched==SCHED_OMP_COLLAPSE2){
        int64_t total=(int64_t)outer_n*((loop==KLON_FIRST)?N:nlev);
        int64_t chunk=(total+P_NUMA-1)/P_NUMA;
        for(int d=0;d<P_NUMA;d++){
            int64_t lo=d*chunk,hi=std::min((d+1)*chunk,total);
            if(lo>=hi)continue; NU_HOME=d;
            auto m=compute_slice(V,W_vec,N,nlev,sched,cidx,vidx,lo,hi);
            accumulate(m);
        }
    } else {
        int chunk=(outer_n+P_NUMA-1)/P_NUMA;
        for(int d=0;d<P_NUMA;d++){
            int lo=d*chunk,hi=std::min((d+1)*chunk,outer_n);
            if(lo>=hi)continue; NU_HOME=d;
            auto m=compute_slice(V,W_vec,N,nlev,sched,cidx,vidx,lo,hi);
            accumulate(m);
        }
    }
    if(n_domains==0) return {};
    double nd=(double)n_domains;
    return { acc.mu/nd, acc.mu_delta/nd, acc.mu_delta_numa/nd,
             acc.delta_max/nd, acc.T };
}

/* ================================================================ */
/*  Cell-index distributions (synthetic)                            */
/* ================================================================ */
enum CellDist {UNIFORM=0,NORMAL1=1,NORMAL4=2,SEQUENTIAL=3};
static const char* dist_name[]={"uniform","normal_var1","normal_var4","sequential","exact"};

static void gen_cell_idx(int* dst,int V,int N,CellDist dist,std::mt19937& rng){
    std::vector<int> L(N*2);
    switch(dist){
    case UNIFORM:{std::uniform_int_distribution<int>u(0,N-1);for(int i=0;i<N;i++){L[2*i]=u(rng);L[2*i+1]=u(rng);}break;}
    case NORMAL1:{std::normal_distribution<double>nd(0,1);for(int i=0;i<N;i++){int v0=i+1+(int)std::round(nd(rng)),v1=i-1+(int)std::round(nd(rng));L[2*i]=((v0%N)+N)%N;L[2*i+1]=((v1%N)+N)%N;}break;}
    case NORMAL4:{std::normal_distribution<double>nd(0,2);for(int i=0;i<N;i++){int v0=i+1+(int)std::round(nd(rng)),v1=i-1+(int)std::round(nd(rng));L[2*i]=((v0%N)+N)%N;L[2*i+1]=((v1%N)+N)%N;}break;}
    case SEQUENTIAL:for(int i=0;i<N;i++){L[2*i]=(i+1)%N;L[2*i+1]=(i+1)%N;}break;
    }
    for(int je=0;je<N;je++){dst[IN(V,je,0,N)]=L[2*je];dst[IN(V,je,1,N)]=L[2*je+1];}
}
static void gen_vert_idx(int* dst,int V,int N,std::mt19937& rng){
    std::vector<int>p(N);
    std::iota(p.begin(),p.end(),0);std::shuffle(p.begin(),p.end(),rng);
    for(int je=0;je<N;je++)dst[IN(V,je,0,N)]=p[je];
    std::iota(p.begin(),p.end(),0);std::shuffle(p.begin(),p.end(),rng);
    for(int je=0;je<N;je++)dst[IN(V,je,1,N)]=p[je];
}

/* ================================================================ */
/*  Target configurations                                           */
/* ================================================================ */
struct Target{const char*name;const char*csv_target;int block_bytes;int vec_width;};
static const Target targets[]={
    {"CPU_scalar","cpu_scalar",64,1},
    {"CPU_AVX512","cpu_avx512",64,8},
    {"GPU_scalar","gpu_scalar",128,1},
    {"GPU_Warp32","gpu_warp32",128,32},
    {"GPU_Wave64","gpu_wave64",128,64},
};
static constexpr int N_TGT=sizeof(targets)/sizeof(targets[0]);

/* ================================================================ */
/*  Result struct                                                   */
/* ================================================================ */
struct Result{
    int V,ti,si,di;
    const char*target_name,*csv_target,*sched,*dist;
    int block_bytes,vec_width;
    Metrics m;
};

/* ================================================================ */
/*  main                                                            */
/* ================================================================ */
int main(int argc,char**argv){
    int N=(argc>1)?atoi(argv[1]):81920;
    int nlev=(argc>2)?atoi(argv[2]):90;
    BETA=(argc>3)?atoi(argv[3]):1;
    ALPHA=(argc>4)?atof(argv[4]):0.012;
    GAMMA=(argc>5)?atof(argv[5]):1.8;
    P_NUMA=(argc>6)?atoi(argv[6]):4;

    fprintf(stderr,"Cost metrics: N=%d nlev=%d beta=%d alpha=%.4f gamma=%.3f P=%d\n\n",
            N,nlev,BETA,ALPHA,GAMMA,P_NUMA);

    /* ---- Load ICON exact data ---- */
    std::string patch_path = icon_patch_path(9);
    IconEdgeData icon_ed;
    bool have_exact = icon_load_patch(patch_path.c_str(), icon_ed);
    if (have_exact) {
        fprintf(stderr, "ICON exact: n_edges_valid=%d  n_cells=%d  n_verts=%d\n\n",
                icon_ed.n_edges_valid, icon_ed.n_cells, icon_ed.n_verts);
    } else {
        fprintf(stderr, "WARNING: no ICON data at '%s', skipping exact dist\n\n",
                patch_path.c_str());
    }

    /* ---- CSV ---- */
    FILE*csv=fopen("metrics.csv","w");
    fprintf(csv,"nlev,variant,cell_dist,schedule,target,block_bytes,vector_width,"
                "mu,mu_delta,mu_delta_numa,delta_max,"
                "beta,alpha,gamma,P_NUMA,T\n");

    /* ---- Synthetic distributions: 4 variants × 2 schedules × 4 dists × N_TGT ---- */
    int n_synth = 4*2*4*N_TGT;
    std::vector<Result> results(n_synth);

    #pragma omp parallel for schedule(dynamic,1)
    for(int job=0;job<n_synth;job++){
        int rem=job;
        int ti=rem%N_TGT;rem/=N_TGT;
        int d=rem%4;rem/=4;
        int si=rem%2;rem/=2;
        int V=rem+1;

        Schedule sched=(si==0)?SCHED_OMP_FOR:SCHED_OMP_COLLAPSE2;
        LoopOrder loop=iteration_order(sched,V);
        int Wv=targets[ti].vec_width;
        int inner=(loop==KLON_FIRST)?N:nlev;

        Result&r=results[job];
        r.V=V;r.ti=ti;r.si=si;r.di=d;
        r.target_name=targets[ti].name;r.csv_target=targets[ti].csv_target;
        r.sched=sched_name[si];r.dist=dist_name[d];
        r.block_bytes=targets[ti].block_bytes;r.vec_width=Wv;
        r.m={};
        if(Wv>inner)continue;

        std::mt19937 rng(42+d);
        std::vector<int>cidx(N*2),vidx(N*2);
        gen_cell_idx(cidx.data(),V,N,(CellDist)d,rng);
        gen_vert_idx(vidx.data(),V,N,rng);
        BLOCK_BYTES=targets[ti].block_bytes;
        r.m=compute_metrics(V,Wv,N,nlev,sched,cidx.data(),vidx.data());
    }

    /* ---- Exact distribution: 4 variants × 2 schedules × N_TGT ---- */
    int n_exact = have_exact ? 4*2*N_TGT : 0;
    std::vector<Result> exact_results(n_exact);

    if(have_exact){
        int Ne = icon_ed.n_edges_valid;

        #pragma omp parallel for schedule(dynamic,1)
        for(int job=0;job<n_exact;job++){
            int rem=job;
            int ti=rem%N_TGT;rem/=N_TGT;
            int si=rem%2;rem/=2;
            int V=rem+1;

            Schedule sched=(si==0)?SCHED_OMP_FOR:SCHED_OMP_COLLAPSE2;
            LoopOrder loop=iteration_order(sched,V);
            int Wv=targets[ti].vec_width;
            int inner=(loop==KLON_FIRST)?Ne:nlev;

            Result&r=exact_results[job];
            r.V=V;r.ti=ti;r.si=si;r.di=4; /* 4 = exact */
            r.target_name=targets[ti].name;r.csv_target=targets[ti].csv_target;
            r.sched=sched_name[si];r.dist="exact";
            r.block_bytes=targets[ti].block_bytes;r.vec_width=Wv;
            r.m={};
            if(Wv>inner)continue;

            /* Build cidx/vidx in variant-V layout from flat ICON connectivity */
            std::vector<int>cidx(Ne*2),vidx(Ne*2);
            for(int je=0;je<Ne;je++){
                cidx[IN(V,je,0,Ne)] = icon_ed.cell_idx[je*2+0];
                cidx[IN(V,je,1,Ne)] = icon_ed.cell_idx[je*2+1];
                vidx[IN(V,je,0,Ne)] = icon_ed.vert_idx[je*2+0];
                vidx[IN(V,je,1,Ne)] = icon_ed.vert_idx[je*2+1];
            }

            BLOCK_BYTES=targets[ti].block_bytes;
            r.m=compute_metrics(V,Wv,Ne,nlev,sched,cidx.data(),vidx.data());
        }
    }

    /* ---- Write CSV rows ---- */
    auto write_csv=[&](const Result&r){
        if(r.m.T==0)return;
        fprintf(csv,"%d,%d,%s,%s,%s,%d,%d,"
                "%.6f,%.6f,%.6f,%.6f,"
                "%d,%.4f,%.4f,%d,%ld\n",
                nlev,r.V,r.dist,r.sched,r.csv_target,r.block_bytes,r.vec_width,
                r.m.mu,r.m.mu_delta,r.m.mu_delta_numa,r.m.delta_max,
                BETA,ALPHA,GAMMA,P_NUMA,(long)r.m.T);
    };
    for(auto&r:results) write_csv(r);
    for(auto&r:exact_results) write_csv(r);
    fclose(csv);

    /* ---- Stdout table ---- */
    auto fmt=[](double v,char*buf){
        if(v<100)snprintf(buf,16,"%8.3f",v);
        else if(v<10000)snprintf(buf,16,"%8.1f",v);
        else snprintf(buf,16,"%8.0f",v);
    };

    printf("\n  Cost Metrics: N=%d nlev=%d beta=%d alpha=%.4f gamma=%.3f P=%d\n\n",
           N,nlev,BETA,ALPHA,GAMMA,P_NUMA);
    printf("  %-4s %-11s %-8s %-11s %8s %8s %8s %8s\n",
           "V","Target","Schedule","Dist",
           "mu","mu*D","mu*Dn","D_max");
    printf("  %-4s %-11s %-8s %-11s %8s %8s %8s %8s\n",
           "----","-----------","--------","-----------",
           "--------","--------","--------","--------");

    auto print_row=[&](const Result&r){
        if(r.m.T==0)return;
        char b[4][16];
        fmt(r.m.mu,b[0]); fmt(r.m.mu_delta,b[1]);
        fmt(r.m.mu_delta_numa,b[2]); fmt(r.m.delta_max,b[3]);
        printf("  V%-3d %-11s %-8s %-11s %s %s %s %s\n",
               r.V,r.target_name,r.sched,r.dist,
               b[0],b[1],b[2],b[3]);
    };

    int prev_V=-1;
    for(auto&r:results){
        if(r.m.T==0)continue;
        if(prev_V!=-1&&r.V!=prev_V)printf("\n");
        prev_V=r.V;
        print_row(r);
    }
    if(!exact_results.empty()){
        printf("\n  --- Exact (ICON) distribution ---\n\n");
        prev_V=-1;
        for(auto&r:exact_results){
            if(r.m.T==0)continue;
            if(prev_V!=-1&&r.V!=prev_V)printf("\n");
            prev_V=r.V;
            print_row(r);
        }
    }

    printf("\n");
    fprintf(stderr,"\nWritten: metrics.csv\n");
    return 0;
}