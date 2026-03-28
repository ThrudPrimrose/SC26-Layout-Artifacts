/*
 * cost_metrics.cpp -- Cost metrics for z_v_grad_w stencil
 *
 * Metrics:
 *   Block-set: mu, W, delta_raw, delta_max, delta_uma, delta_numa,
 *              mu_delta_raw, mu_delta_uma, mu_delta_numa,
 *              sigma, cost_uma, cost_numa, span, W_delta
 *   Per-ref:   ref_SL, ref_RL, ref_SR, ref_RR, ref_cost
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

static int    BETA = 4;
static double ALPHA = 1.0/8.0;
static double GAMMA = 8.0;
static int    P_NUMA = 4;
static int    BLOCK_BYTES = 64;
#pragma omp threadprivate(BLOCK_BYTES)

inline int IC(int V, int je, int jk, int N, int nlev) {
    return (V <= 2) ? je + jk * N : jk + je * nlev;
}
inline int IN(int V, int je, int n, int N) {
    return (V == 1 || V == 3) ? je + n * N : n + je * 2;
}

enum ArrID {
    A_OUT=0, A_VN_IE, A_W, A_Z_VT_IE, A_Z_W_V,
    A_INV_DUAL, A_INV_PRIMAL, A_TANGENT,
    A_CELL_IDX, A_VERT_IDX, NUM_ARR
};
static const int ebytes[] = {8,8,8,8,8, 8,8,8, 4,4};

inline int64_t blk(int elem_idx, int arr) {
    return (int64_t)elem_idx * ebytes[arr] / BLOCK_BYTES;
}

enum Schedule  { SCHED_OMP_FOR=0, SCHED_OMP_COLLAPSE2=1 };
static const char* sched_name[] = {"omp_for","omp_collapse2"};
enum LoopOrder { KLON_FIRST=0, KLEV_FIRST=1 };
static const char* loop_name[] = {"klon_first","klev_first"};

inline LoopOrder iteration_order(Schedule sched, int V) {
    if (sched == SCHED_OMP_COLLAPSE2) return KLON_FIRST;
    return (V <= 2) ? KLON_FIRST : KLEV_FIRST;
}

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
    if (V <= 2) { int jk = elem/N; int chunk=(nlev+P_NUMA-1)/P_NUMA; int d=jk/chunk; return (d<P_NUMA)?d:P_NUMA-1; }
    else { int je = elem/nlev; int chunk=(N+P_NUMA-1)/P_NUMA; int d=je/chunk; return (d<P_NUMA)?d:P_NUMA-1; }
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

inline double w_uma(int64_t raw_dist) { return (raw_dist<BETA)?ALPHA:1.0; }
inline double w_numa(int arr, int V, int64_t baddr, int64_t raw_dist, int N, int nlev) {
    if(nu_block(arr,V,baddr,N,nlev)!=NU_HOME) return GAMMA;
    return (raw_dist<BETA)?ALPHA:1.0;
}

static constexpr int NREFS = 14;
struct Ref { int arr; int64_t block; };

inline void collect_refs(int V, int je, int jk, int N, int nlev,
                         const int* cidx, const int* vidx, Ref refs[NREFS]) {
    int c2d=IC(V,je,jk,N,nlev);
    int ci0=cidx[IN(V,je,0,N)],ci1=cidx[IN(V,je,1,N)];
    int vi0=vidx[IN(V,je,0,N)],vi1=vidx[IN(V,je,1,N)];
    refs[0]={A_OUT,blk(c2d,A_OUT)}; refs[1]={A_VN_IE,blk(c2d,A_VN_IE)};
    refs[2]={A_INV_DUAL,blk(je,A_INV_DUAL)};
    refs[3]={A_W,blk(IC(V,ci0,jk,N,nlev),A_W)}; refs[4]={A_W,blk(IC(V,ci1,jk,N,nlev),A_W)};
    refs[5]={A_Z_VT_IE,blk(c2d,A_Z_VT_IE)};
    refs[6]={A_INV_PRIMAL,blk(je,A_INV_PRIMAL)}; refs[7]={A_TANGENT,blk(je,A_TANGENT)};
    refs[8]={A_Z_W_V,blk(IC(V,vi0,jk,N,nlev),A_Z_W_V)}; refs[9]={A_Z_W_V,blk(IC(V,vi1,jk,N,nlev),A_Z_W_V)};
    refs[10]={A_CELL_IDX,blk(IN(V,je,0,N),A_CELL_IDX)}; refs[11]={A_CELL_IDX,blk(IN(V,je,1,N),A_CELL_IDX)};
    refs[12]={A_VERT_IDX,blk(IN(V,je,0,N),A_VERT_IDX)}; refs[13]={A_VERT_IDX,blk(IN(V,je,1,N),A_VERT_IDX)};
}

struct BlockSet {
    std::vector<int64_t> a[NUM_ARR];
    void clear(){for(int i=0;i<NUM_ARR;i++)a[i].clear();}
    void add(int arr,int64_t b){a[arr].push_back(b);}
    void finalize(){for(int i=0;i<NUM_ARR;i++){auto&v=a[i];std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());}}
    int total()const{int n=0;for(int i=0;i<NUM_ARR;i++)n+=(int)a[i].size();return n;}
    bool has(int arr,int64_t b)const{return std::binary_search(a[arr].begin(),a[arr].end(),b);}
};

struct Metrics {
    double mu, W;
    double delta_raw, delta_max, delta_uma, delta_numa;
    double mu_delta_raw, mu_delta_uma, mu_delta_numa;
    double sigma, cost_uma, cost_numa;
    double ref_SL, ref_RL, ref_SR, ref_RR, ref_cost;
    double span, W_delta;
    int64_t T;
};

static inline void process_step(
    int V, int W_vec, int N, int nlev,
    LoopOrder loop, int je0, int jk0,
    const int* cidx, const int* vidx,
    BlockSet& prev, BlockSet& curr,
    int64_t prev_ref[NREFS], bool& have_prev_ref,
    int64_t& T,
    double& s_mu, double& s_W,
    double& s_d_raw, double& s_d_max, double& s_d_uma, double& s_d_numa,
    double& s_md_raw, double& s_md_uma, double& s_md_numa,
    double& s_sigma, double& s_cost_uma, double& s_cost_numa,
    double& s_ref_SL, double& s_ref_RL, double& s_ref_SR, double& s_ref_RR,
    double& s_span)
{
    Ref all_refs[NREFS];
    curr.clear();
    for(int w=0;w<W_vec;w++){
        int je=(loop==KLON_FIRST)?je0+w:je0;
        int jk=(loop==KLON_FIRST)?jk0:jk0+w;
        Ref refs[NREFS];
        collect_refs(V,je,jk,N,nlev,cidx,vidx,refs);
        for(int r=0;r<NREFS;r++){curr.add(refs[r].arr,refs[r].block);if(w==0)all_refs[r]=refs[r];}
    }
    curr.finalize();

    /* W and span */
    {
        int64_t bmin=INT64_MAX, bmax=0; int W_t=0;
        for(int a=0;a<NUM_ARR;a++){
            W_t+=(int)curr.a[a].size();
            for(int64_t b:curr.a[a]){if(b<bmin)bmin=b;if(b>bmax)bmax=b;}
        }
        s_W += W_t;
        s_span += (bmin<=bmax)?(double)(bmax-bmin):0.0;
    }

    /* Per-reference stride */
    {
        int nSL=0,nRL=0,nSR=0,nRR=0;
        for(int r=0;r<NREFS;r++){
            int64_t b=all_refs[r].block; int a=all_refs[r].arr;
            bool remote=(nu_block(a,V,b,N,nlev)!=NU_HOME);
            if(!have_prev_ref){if(remote)nSR++;else nSL++;}
            else{
                int64_t stride=std::abs(b-prev_ref[r]); bool streaming=(stride<BETA);
                if(!remote&&streaming)nSL++;else if(!remote&&!streaming)nRL++;
                else if(remote&&streaming)nSR++;else nRR++;
            }
            prev_ref[r]=b;
        }
        s_ref_SL+=nSL;s_ref_RL+=nRL;s_ref_SR+=nSR;s_ref_RR+=nRR;
        have_prev_ref=true;
    }

    /* Block-set distance metrics */
    if(T==0){
        int nb=curr.total();
        s_mu+=nb; s_d_raw+=1.0; s_d_max+=1.0; s_d_uma+=1.0; s_d_numa+=1.0;
        s_md_raw+=nb; s_md_uma+=nb; s_md_numa+=nb;
        s_sigma+=nb; s_cost_uma+=nb*ALPHA; s_cost_numa+=nb*ALPHA;
    } else {
        int nc=0;
        double sum_raw=0,sum_dmax=0,sum_uma=0,sum_numa=0;
        double sum_addr=0,sum_w_uma=0,sum_w_numa=0;

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
                sum_raw+=1.0; sum_dmax+=1.0; sum_uma+=1.0;
                double wn=w_numa(a,V,b,INT64_MAX,N,nlev);
                sum_numa+=wn; sum_addr+=1.0;
                sum_w_uma+=1.0; sum_w_numa+=wn;
            } else {
                double rmin=(double)best_min;
                double rmax=(double)best_max;
                double wu=w_uma(best_min);
                double wn=w_numa(a,V,b,best_min,N,nlev);
                sum_raw+=rmin; sum_dmax+=rmax;
                sum_uma+=wu*rmin; sum_numa+=wn*rmin;
                sum_addr+=rmin*((double)BLOCK_BYTES/ebytes[a]);
                sum_w_uma+=wu; sum_w_numa+=wn;
            }
          }
        }
        s_mu+=nc;
        if(nc>0){s_d_raw+=sum_raw/nc; s_d_max+=sum_dmax/nc; s_d_uma+=sum_uma/nc; s_d_numa+=sum_numa/nc;}
        s_md_raw+=sum_raw; s_md_uma+=sum_uma; s_md_numa+=sum_numa;
        s_sigma+=sum_addr; s_cost_uma+=sum_w_uma; s_cost_numa+=sum_w_numa;
    }
    T++;
    std::swap(prev,curr);
}

Metrics compute_slice(int V, int W_vec, int N, int nlev,
                      Schedule sched, const int* cidx, const int* vidx,
                      int64_t range_lo, int64_t range_hi) {
    G_SCHED=sched;
    LoopOrder loop=iteration_order(sched,V);
    BlockSet prev,curr;
    int64_t prev_ref[NREFS]; bool have_prev_ref=false;
    memset(prev_ref,0,sizeof(prev_ref));

    int64_t T=0; double s_mu=0,s_W=0;
    double s_d_raw=0,s_d_max=0,s_d_uma=0,s_d_numa=0;
    double s_md_raw=0,s_md_uma=0,s_md_numa=0;
    double s_sigma=0,s_cost_uma=0,s_cost_numa=0;
    double s_ref_SL=0,s_ref_RL=0,s_ref_SR=0,s_ref_RR=0;
    double s_span=0;

    if(sched==SCHED_OMP_COLLAPSE2){
        int64_t lin=range_lo;
        while(lin<range_hi){
            int jk=(int)(lin/N); int je_start=(int)(lin%N);
            int64_t row_end=std::min((int64_t)(jk+1)*N,range_hi);
            int je_end=(int)(row_end-(int64_t)jk*N);
            for(int je0=je_start;je0+W_vec<=je_end;je0+=W_vec)
                process_step(V,W_vec,N,nlev,KLON_FIRST,je0,jk,cidx,vidx,
                    prev,curr,prev_ref,have_prev_ref,T,
                    s_mu,s_W,s_d_raw,s_d_max,s_d_uma,s_d_numa,
                    s_md_raw,s_md_uma,s_md_numa,
                    s_sigma,s_cost_uma,s_cost_numa,
                    s_ref_SL,s_ref_RL,s_ref_SR,s_ref_RR,s_span);
            lin=row_end;
        }
    } else {
        int inner_n=(loop==KLON_FIRST)?N:nlev;
        for(int64_t outer=range_lo;outer<range_hi;outer++)
            for(int inner0=0;inner0+W_vec<=inner_n;inner0+=W_vec){
                int je=(loop==KLON_FIRST)?inner0:(int)outer;
                int jk=(loop==KLON_FIRST)?(int)outer:inner0;
                process_step(V,W_vec,N,nlev,loop,je,jk,cidx,vidx,
                    prev,curr,prev_ref,have_prev_ref,T,
                    s_mu,s_W,s_d_raw,s_d_max,s_d_uma,s_d_numa,
                    s_md_raw,s_md_uma,s_md_numa,
                    s_sigma,s_cost_uma,s_cost_numa,
                    s_ref_SL,s_ref_RL,s_ref_SR,s_ref_RR,s_span);
            }
    }
    if(T==0) return {};
    double dT=(double)T;
    double rSL=s_ref_SL/dT,rRL=s_ref_RL/dT,rSR=s_ref_SR/dT,rRR=s_ref_RR/dT;
    double Wavg=s_W/dT, dmin=s_d_raw/dT;
    return {
        s_mu/dT, Wavg,
        dmin, s_d_max/dT, s_d_uma/dT, s_d_numa/dT,
        s_md_raw/dT, s_md_uma/dT, s_md_numa/dT,
        s_sigma/dT, s_cost_uma/dT, s_cost_numa/dT,
        rSL, rRL, rSR, rRR,
        ALPHA*rSL + 1.0*rRL + ALPHA*GAMMA*rSR + GAMMA*rRR,
        s_span/dT, Wavg*dmin,
        T
    };
}

Metrics compute_metrics(int V, int W_vec, int N, int nlev,
                        Schedule sched, const int* cidx, const int* vidx) {
    LoopOrder loop=iteration_order(sched,V);
    int outer_n=(loop==KLON_FIRST)?nlev:N;
    Metrics acc={}; int n_domains=0;

    auto accumulate=[&](Metrics& m){
        if(m.T==0) return;
        acc.mu+=m.mu; acc.W+=m.W;
        acc.delta_raw+=m.delta_raw; acc.delta_max+=m.delta_max;
        acc.delta_uma+=m.delta_uma; acc.delta_numa+=m.delta_numa;
        acc.mu_delta_raw+=m.mu_delta_raw; acc.mu_delta_uma+=m.mu_delta_uma;
        acc.mu_delta_numa+=m.mu_delta_numa; acc.sigma+=m.sigma;
        acc.cost_uma+=m.cost_uma; acc.cost_numa+=m.cost_numa;
        acc.ref_SL+=m.ref_SL; acc.ref_RL+=m.ref_RL;
        acc.ref_SR+=m.ref_SR; acc.ref_RR+=m.ref_RR;
        acc.ref_cost+=m.ref_cost;
        acc.span+=m.span; acc.W_delta+=m.W_delta;
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
    return {
        acc.mu/nd, acc.W/nd,
        acc.delta_raw/nd, acc.delta_max/nd, acc.delta_uma/nd, acc.delta_numa/nd,
        acc.mu_delta_raw/nd, acc.mu_delta_uma/nd, acc.mu_delta_numa/nd,
        acc.sigma/nd, acc.cost_uma/nd, acc.cost_numa/nd,
        acc.ref_SL/nd, acc.ref_RL/nd, acc.ref_SR/nd, acc.ref_RR/nd,
        acc.ref_cost/nd, acc.span/nd, acc.W_delta/nd, acc.T
    };
}

enum CellDist {UNIFORM=0,NORMAL1=1,NORMAL4=2,SEQUENTIAL=3};
static const char* dist_name[]={"uniform","normal_var1","normal_var4","sequential"};

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

struct Target{const char*name;const char*csv_target;int block_bytes;int vec_width;};
static const Target targets[]={
    {"CPU scalar","cpu_scalar",64,1},{"CPU AVX-512","cpu_avx512",64,8},
    {"GPU scalar","gpu_scalar",128,1},{"GPU Warp32","gpu_warp32",128,32},
    {"GPU Wave64","gpu_wave64",128,64},
};
static constexpr int N_TGT=sizeof(targets)/sizeof(targets[0]);

struct Result{
    int V,ti,si,di;
    const char*target_name,*csv_target,*sched,*loop,*dist;
    int block_bytes,vec_width; Metrics m;
};

int main(int argc,char**argv){
    int N=(argc>1)?atoi(argv[1]):81920;
    int nlev=(argc>2)?atoi(argv[2]):96;
    BETA=(argc>3)?atoi(argv[3]):4;
    ALPHA=(argc>4)?atof(argv[4]):1.0/8.0;
    GAMMA=(argc>5)?atof(argv[5]):8.0;
    P_NUMA=(argc>6)?atoi(argv[6]):4;

    fprintf(stderr,"Cost metrics: N=%d nlev=%d beta=%d alpha=%.4f gamma=%.3f P=%d\n\n",
            N,nlev,BETA,ALPHA,GAMMA,P_NUMA);

    FILE*csv=fopen("metrics.csv","w");
    fprintf(csv,"nlev,variant,cell_dist,loop_order,target,block_bytes,vector_width,"
                "mu,W,delta_raw,delta_max,delta_uma,delta_numa,"
                "mu_delta_raw,mu_delta_uma,mu_delta_numa,"
                "sigma,cost_uma,cost_numa,"
                "ref_SL,ref_RL,ref_SR,ref_RR,ref_cost,"
                "span,W_delta,"
                "beta,alpha,gamma,P_NUMA,T\n");

    int n_jobs=4*2*4*N_TGT;
    std::vector<Result> results(n_jobs);

    #pragma omp parallel for schedule(dynamic,1)
    for(int job=0;job<n_jobs;job++){
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
        r.sched=sched_name[si];r.loop=loop_name[loop];r.dist=dist_name[d];
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

    for(int job=0;job<n_jobs;job++){
        const Result&r=results[job];
        if(r.m.T==0)continue;
        fprintf(csv,"%d,%d,%s,%s,%s,%d,%d,"
                "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,%.6f,%.6f,"
                "%.6f,%.6f,"
                "%d,%.4f,%.4f,%d,%ld\n",
                nlev,r.V,r.dist,r.loop,r.csv_target,r.block_bytes,r.vec_width,
                r.m.mu,r.m.W,r.m.delta_raw,r.m.delta_max,r.m.delta_uma,r.m.delta_numa,
                r.m.mu_delta_raw,r.m.mu_delta_uma,r.m.mu_delta_numa,
                r.m.sigma,r.m.cost_uma,r.m.cost_numa,
                r.m.ref_SL,r.m.ref_RL,r.m.ref_SR,r.m.ref_RR,r.m.ref_cost,
                r.m.span,r.m.W_delta,
                BETA,ALPHA,GAMMA,P_NUMA,(long)r.m.T);
    }
    fclose(csv);

    /* ── Tables ── */
    auto fmt=[](double v,char*buf){
        if(v<100)snprintf(buf,16,"%8.3f",v);
        else if(v<10000)snprintf(buf,16,"%8.1f",v);
        else snprintf(buf,16,"%8.0f",v);
    };

    printf("\n  Cost Metrics: N=%d nlev=%d beta=%d alpha=%.4f gamma=%.3f P=%d\n\n",
           N,nlev,BETA,ALPHA,GAMMA,P_NUMA);
    printf("  %-4s %-11s %-8s %-11s %-10s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
           "V","Target","Schedule","Dist","Loop",
           "mu","W","D_min","D_max","D_uma","D_numa",
           "muD_raw","muD_uma","muD_numa",
           "sigma","w_uma","w_numa","span","W*D_min");
    printf("  %-4s %-11s %-8s %-11s %-10s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
           "----","-----------","--------","-----------","----------",
           "--------","--------","--------","--------","--------","--------",
           "--------","--------","--------","--------","--------","--------","--------","--------");

    int prev_V=-1;
    for(int job=0;job<n_jobs;job++){
        const Result&r=results[job];
        if(r.m.T==0)continue;
        if(prev_V!=-1&&r.V!=prev_V)printf("\n");
        prev_V=r.V;
        char b[19][16];
        fmt(r.m.mu,b[0]);fmt(r.m.W,b[1]);
        fmt(r.m.delta_raw,b[2]);fmt(r.m.delta_max,b[3]);
        fmt(r.m.delta_uma,b[4]);fmt(r.m.delta_numa,b[5]);
        fmt(r.m.mu_delta_raw,b[6]);fmt(r.m.mu_delta_uma,b[7]);fmt(r.m.mu_delta_numa,b[8]);
        fmt(r.m.sigma,b[9]);fmt(r.m.cost_uma,b[10]);fmt(r.m.cost_numa,b[11]);
        fmt(r.m.span,b[12]);fmt(r.m.W_delta,b[13]);
        printf("  V%-3d %-11s %-8s %-11s %-10s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
               r.V,r.target_name,r.sched,r.dist,r.loop,
               b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9],b[10],b[11],b[12],b[13]);
    }

    printf("\n\n  Per-Reference Stride Metrics:\n\n");
    printf("  %-4s %-11s %-8s %-11s %-10s %8s %8s %8s %8s %10s\n",
           "V","Target","Schedule","Dist","Loop","ref_SL","ref_RL","ref_SR","ref_RR","ref_cost");
    printf("  %-4s %-11s %-8s %-11s %-10s %8s %8s %8s %8s %10s\n",
           "----","-----------","--------","-----------","----------",
           "--------","--------","--------","--------","----------");
    prev_V=-1;
    for(int job=0;job<n_jobs;job++){
        const Result&r=results[job];
        if(r.m.T==0)continue;
        if(prev_V!=-1&&r.V!=prev_V)printf("\n");
        prev_V=r.V;
        printf("  V%-3d %-11s %-8s %-11s %-10s %8.3f %8.3f %8.3f %8.3f %10.3f\n",
               r.V,r.target_name,r.sched,r.dist,r.loop,
               r.m.ref_SL,r.m.ref_RL,r.m.ref_SR,r.m.ref_RR,r.m.ref_cost);
    }
    printf("\n");
    fprintf(stderr,"\nWritten: metrics.csv\n");
    return 0;
}