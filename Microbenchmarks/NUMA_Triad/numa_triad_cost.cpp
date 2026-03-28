/*
 * numa_triad_cost.cpp -- Three-tier cost metrics for NUMA triad benchmark
 *
 * Same metric definitions as cost_metrics.cpp (μ, Δ, σ).
 * Triad kernel: C[i] = A[i] + s * B[i], 3 arrays, stride-1.
 *
 * Compile: g++ -O3 -std=c++17 -o numa_triad_cost numa_triad_cost.cpp
 * Run:     ./numa_triad_cost [N] [P_threads] [β] [α] [γ] [P_NUMA]
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <climits>
#include <unistd.h>

static int    BETA=64;
static double ALPHA=0.16, GAMMA=2.0;
static int    P_NUMA=4, P_THR=96;
#define CL_BYTES 64
#define CL_DOUBLES (CL_BYTES/(int)sizeof(double))

/* ---- Schedules ---- */
enum SchedType { SCHED_BLOCK, SCHED_SHIFT1, SCHED_INTERLEAVE };
static const char* sched_label[] = {"block","shift1","interleave"};

inline int owner_thread(SchedType s, int64_t i, int64_t N, int P, int chunk) {
    switch(s){
    case SCHED_BLOCK: { int64_t b=(N+P-1)/P; int t=(int)(i/b);
        return (t<P)?t:P-1; }
    case SCHED_SHIFT1: { int64_t b=(N+P-1)/P; int t=(int)(i/b);
        t=(t<P)?t:P-1; return (t+1)%P; }
    case SCHED_INTERLEAVE: { int c=(chunk>0)?chunk:1; return (int)(i/c)%P; }
    } return 0;
}

/* ---- Placement function ν ---- */
inline int thread_to_numa(int tid) {
    int per=(P_THR+P_NUMA-1)/P_NUMA;
    int d=tid/per; return (d<P_NUMA)?d:P_NUMA-1;
}
inline int elem_numa(SchedType touch, int64_t i, int64_t N, int tc) {
    return thread_to_numa(owner_thread(touch,i,N,P_THR,tc));
}

/* ---- Weight functions ---- */
static constexpr int NU_HOME = 0;

inline double w_uma(int64_t raw_dist) {
    return (raw_dist < BETA) ? ALPHA : 1.0;
}
inline double w_numa(int64_t baddr, int64_t raw_dist,
                     SchedType touch, int64_t N, int tc) {
    int64_t elem = baddr * CL_DOUBLES;
    if (elem_numa(touch, elem, N, tc) != NU_HOME) return GAMMA;
    return (raw_dist < BETA) ? ALPHA : 1.0;
}

/* ---- Metrics ---- */
static constexpr int N_ARR = 3;
struct Metrics { double mu, delta, sigma; int64_t T; double frac_remote; };

Metrics compute_metrics(int64_t N, int W,
                        SchedType touch, int tc,
                        SchedType compute, int cc) {
    std::vector<int64_t> elems;
    switch(compute){
    case SCHED_BLOCK: {
        int64_t hi=N/P_THR;
        for(int64_t i=0;i<hi;i++) elems.push_back(i); break;}
    case SCHED_SHIFT1: {
        int64_t lo=N/P_THR, hi=2*(N/P_THR);
        if(hi>N)hi=N;
        for(int64_t i=lo;i<hi;i++) elems.push_back(i); break;}
    case SCHED_INTERLEAVE: {
        int c=(cc>0)?cc:1;
        for(int64_t base=0;base<N;base+=(int64_t)P_THR*c){
            int64_t end=base+c; if(end>N)end=N;
            for(int64_t i=base;i<end;i++) elems.push_back(i);
        } break;}
    }
    int64_t ni=(int64_t)elems.size();
    if(ni==0) return {0,0,0,0,0};

    int64_t T=0, remote=0, total_new=0;
    double s_mu=0, s_delta=0, s_sigma=0;
    std::vector<int64_t> prev[N_ARR], cur[N_ARR];

    for(int64_t it=0; it+W<=ni; it+=W){
        for(int a=0;a<N_ARR;a++) cur[a].clear();
        for(int w=0;w<W;w++){
            int64_t cl=elems[it+w]/CL_DOUBLES;
            for(int a=0;a<N_ARR;a++) cur[a].push_back(cl);
        }
        for(int a=0;a<N_ARR;a++){
            std::sort(cur[a].begin(),cur[a].end());
            cur[a].erase(std::unique(cur[a].begin(),cur[a].end()),cur[a].end());
        }
        if(T==0){
            int nb=0; for(int a=0;a<N_ARR;a++) nb+=(int)cur[a].size();
            s_mu+=nb; s_delta+=1.0; s_sigma+=nb; total_new+=nb;
            for(int a=0;a<N_ARR;a++)
                for(int64_t b:cur[a])
                    if(elem_numa(touch,b*CL_DOUBLES,N,tc)!=NU_HOME) remote++;
        } else {
            int nc=0; double d_uma=0, d_numa=0;
            for(int a=0;a<N_ARR;a++){
              for(int64_t b:cur[a]){
                if(std::binary_search(prev[a].begin(),prev[a].end(),b)) continue;
                nc++; total_new++;
                int64_t best=INT64_MAX;
                for(int64_t bp:prev[a]){ int64_t d=std::abs(b-bp); if(d<best)best=d; }
                if(best==INT64_MAX){ d_uma+=1.0; d_numa+=1.0; }
                else {
                    d_uma  += w_uma(best)*(double)best;
                    d_numa += w_numa(b,best,touch,N,tc)*(double)best;
                }
                if(elem_numa(touch,b*CL_DOUBLES,N,tc)!=NU_HOME) remote++;
              }
            }
            s_mu+=nc; s_delta+=(nc>0)?d_uma/nc:0.0; s_sigma+=d_numa;
        }
        T++;
        for(int a=0;a<N_ARR;a++) std::swap(prev[a],cur[a]);
    }
    return { s_mu/T, s_delta/T, s_sigma/T, T,
             total_new>0?(double)remote/total_new:0.0 };
}

/* ---- Configurations ---- */
struct Config { const char* name; SchedType touch,compute; int tp,cp; };
static const Config configs[] = {
    {"local",              SCHED_BLOCK,      SCHED_BLOCK,      0,  0},
    {"all_remote",         SCHED_BLOCK,      SCHED_SHIFT1,     0,  0},
    {"touch_interl_1pg",   SCHED_INTERLEAVE, SCHED_BLOCK,      1,  0},
    {"touch_interl_16pg",  SCHED_INTERLEAVE, SCHED_BLOCK,      16, 0},
    {"compute_interl_1pg", SCHED_BLOCK,      SCHED_INTERLEAVE, 0,  1},
    {"compute_interl_16pg",SCHED_BLOCK,      SCHED_INTERLEAVE, 0,  16},
};
static constexpr int N_CFG = sizeof(configs)/sizeof(configs[0]);

int main(int argc, char** argv) {
    int64_t N=(argc>1)?atol(argv[1]):268435456;
    P_THR=(argc>2)?atoi(argv[2]):96;
    BETA =(argc>3)?atoi(argv[3]):(int)(sysconf(_SC_PAGESIZE)/CL_BYTES);
    ALPHA=(argc>4)?atof(argv[4]):0.16;
    GAMMA=(argc>5)?atof(argv[5]):2.0;
    P_NUMA=(argc>6)?atoi(argv[6]):4;

    long ps=sysconf(_SC_PAGESIZE);
    int pcl=(int)(ps/CL_BYTES), pe=pcl*CL_DOUBLES;
    int Ws[]={1,8,32};

    fprintf(stderr,"NUMA Triad Cost: N=%ld P=%d β=%d α=%.3f γ=%.3f P_NUMA=%d\n\n",
            (long)N,P_THR,BETA,ALPHA,GAMMA,P_NUMA);
    printf("config,touch,compute,touch_chunk_cls,compute_chunk_cls,"
           "W,T,mu,delta,sigma,frac_remote,N,P,beta,alpha,gamma,P_NUMA\n");

    for(int ci=0;ci<N_CFG;ci++){
        auto& c=configs[ci];
        int tc_cls=c.tp*pcl, cc_cls=c.cp*pcl;
        int tc_el=c.tp*pe,   cc_el=c.cp*pe;
        for(int wi=0;wi<3;wi++){
            int W=Ws[wi];
            fprintf(stderr,"  %-24s W=%-3d ...",c.name,W);
            auto m=compute_metrics(N,W,c.touch,tc_el,c.compute,cc_el);
            fprintf(stderr,"  μ=%.2f Δ=%.1f σ=%.1f remote=%.0f%%\n",
                    m.mu,m.delta,m.sigma,m.frac_remote*100);
            /*
            printf("%s,%s,%s,%d,%d,%d,%ld,%.6f,%.6f,%.6f,%.6f,"
                   "%ld,%d,%d,%.4f,%.4f,%d\n",
                   c.name,sched_label[c.touch],sched_label[c.compute],
                   tc_cls,cc_cls,W,(long)m.T,m.mu,m.delta,m.sigma,
                   m.frac_remote,(long)N,P_THR,BETA,ALPHA,GAMMA,P_NUMA);
            */
        }
    }
    return 0;
}