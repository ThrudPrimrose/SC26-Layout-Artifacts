// benchmark_local.cu — ICON z_ekinh with LOCAL neighbor lists + CSV output
// Neighbors: self, jc-1, jc+1 (all same block) — best-case locality

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

static const int NBLKS = 2;
static const int NITER = 50;

#define CK(e) do { cudaError_t _r=(e); if(_r!=cudaSuccess){ fprintf(stderr,"CUDA %s @L%d: %s\n",#e,__LINE__,cudaGetErrorString(_r)); exit(1); } } while(0)

#define IA(jc,jk,jb,NP,NL)      ((jc)  + (long)(jk)*(NP)    + (long)(jb)*(NP)*(NL))
#define IB(jc,jk,jb,NP,NL)      ((jk)  + (long)(jc)*(NL)    + (long)(jb)*(NL)*(NP))
#define IBLN_A(jc,n,jb,NP)      ((jc)  + (long)(n)*(NP)     + (long)(jb)*(NP)*3)
#define IMAP_A(jc,jb,n,NP,NB)   ((jc)  + (long)(jb)*(NP)    + (long)(n)*(NP)*(NB))
#define IBLN_B(n,jc,jb,NP)      ((n)   + 3L*(jc)            + 3L*(NP)*(jb))
#define IMAP_B(n,jc,jb,NP)      ((n)   + 3L*(jc)            + 3L*(NP)*(jb))

template<int ELEMS, int BSZJC, int BSZJK = 1, int ELEMC = 1>
__global__ void ker_A(float* __restrict__ out, const float* __restrict__ bln,
                      const float* __restrict__ kin, const int* __restrict__ ix,
                      const int* __restrict__ ib, int NP, int NL, int NB) {
    int jc_base=(blockIdx.x*BSZJC+threadIdx.x)*ELEMC, jk_base=(blockIdx.y*BSZJK+threadIdx.y)*ELEMS, jb=blockIdx.z;
    if(jc_base>=NP||jb>=NB) return;
    float b0[ELEMC],b1[ELEMC],b2[ELEMC]; int i0[ELEMC],i1[ELEMC],i2[ELEMC],bl0[ELEMC],bl1[ELEMC],bl2[ELEMC];
    #pragma unroll
    for(int c=0;c<ELEMC;c++){int jc=jc_base+c;if(jc>=NP)break;
        b0[c]=bln[IBLN_A(jc,0,jb,NP)];b1[c]=bln[IBLN_A(jc,1,jb,NP)];b2[c]=bln[IBLN_A(jc,2,jb,NP)];
        i0[c]=ix[IMAP_A(jc,jb,0,NP,NB)];bl0[c]=ib[IMAP_A(jc,jb,0,NP,NB)];
        i1[c]=ix[IMAP_A(jc,jb,1,NP,NB)];bl1[c]=ib[IMAP_A(jc,jb,1,NP,NB)];
        i2[c]=ix[IMAP_A(jc,jb,2,NP,NB)];bl2[c]=ib[IMAP_A(jc,jb,2,NP,NB)];}
    #pragma unroll
    for(int e=0;e<ELEMS;e++){if(jk_base+e>=NL)break;
        #pragma unroll
        for(int c=0;c<ELEMC;c++){int jc=jc_base+c;if(jc>=NP)break;
            out[IA(jc,jk_base+e,jb,NP,NL)]=
                b0[c]*kin[IA(i0[c],jk_base+e,bl0[c],NP,NL)]+
                b1[c]*kin[IA(i1[c],jk_base+e,bl1[c],NP,NL)]+
                b2[c]*kin[IA(i2[c],jk_base+e,bl2[c],NP,NL)];}}
}

template<int ELEMS, int BSZJK, int BSZJC>
__global__ void ker_B(float* __restrict__ out, const float* __restrict__ bln,
                      const float* __restrict__ kin, const int* __restrict__ ix,
                      const int* __restrict__ ib, int NP, int NL, int NB) {
    int jk=(blockIdx.x*BSZJK+threadIdx.x)*ELEMS, jc=blockIdx.y*BSZJC+threadIdx.y, jb=blockIdx.z;
    if(jc>=NP||jb>=NB) return;
    float b0=bln[IBLN_B(0,jc,jb,NP)],b1=bln[IBLN_B(1,jc,jb,NP)],b2=bln[IBLN_B(2,jc,jb,NP)];
    int i0=ix[IMAP_B(0,jc,jb,NP)],bl0=ib[IMAP_B(0,jc,jb,NP)];
    int i1=ix[IMAP_B(1,jc,jb,NP)],bl1=ib[IMAP_B(1,jc,jb,NP)];
    int i2=ix[IMAP_B(2,jc,jb,NP)],bl2=ib[IMAP_B(2,jc,jb,NP)];
    #pragma unroll
    for(int e=0;e<ELEMS;e++){if(jk+e>=NL)break;
        out[IB(jc,jk+e,jb,NP,NL)]=b0*kin[IB(i0,jk+e,bl0,NP,NL)]+b1*kin[IB(i1,jk+e,bl1,NP,NL)]+b2*kin[IB(i2,jk+e,bl2,NP,NL)];}
}
template<int ELEMS, int BSZJC, int BSZJK = 1, int ELEMC = 1>
__global__ void ker_C(float* __restrict__ out, const float* __restrict__ bln,
                      const float* __restrict__ kin, const int* __restrict__ ix,
                      const int* __restrict__ ib, int NP, int NL, int NB) {
    int jc_base = (blockIdx.x * BSZJC + threadIdx.x) * ELEMC;
    int jk_base = (blockIdx.y * BSZJK + threadIdx.y) * ELEMS;
    int jb      =  blockIdx.z;
    if (jc_base >= NP || jb >= NB) return;
    float b0[ELEMC], b1[ELEMC], b2[ELEMC];
    int i0[ELEMC], i1[ELEMC], i2[ELEMC], bl0[ELEMC], bl1[ELEMC], bl2[ELEMC];
    #pragma unroll
    for (int c = 0; c < ELEMC; c++) {
        int jc = jc_base + c; if (jc >= NP) break;
        b0[c]=bln[IBLN_B(0,jc,jb,NP)]; b1[c]=bln[IBLN_B(1,jc,jb,NP)]; b2[c]=bln[IBLN_B(2,jc,jb,NP)];
        i0[c]=ix[IMAP_B(0,jc,jb,NP)]; bl0[c]=ib[IMAP_B(0,jc,jb,NP)];
        i1[c]=ix[IMAP_B(1,jc,jb,NP)]; bl1[c]=ib[IMAP_B(1,jc,jb,NP)];
        i2[c]=ix[IMAP_B(2,jc,jb,NP)]; bl2[c]=ib[IMAP_B(2,jc,jb,NP)];
    }
    #pragma unroll
    for (int e = 0; e < ELEMS; e++) { if (jk_base+e >= NL) break;
        #pragma unroll
        for (int c = 0; c < ELEMC; c++) { int jc=jc_base+c; if (jc>=NP) break;
            out[IA(jc,jk_base+e,jb,NP,NL)] =
                b0[c]*kin[IA(i0[c],jk_base+e,bl0[c],NP,NL)] +
                b1[c]*kin[IA(i1[c],jk_base+e,bl1[c],NP,NL)] +
                b2[c]*kin[IA(i2[c],jk_base+e,bl2[c],NP,NL)];
    }}
}
// ─── Local neighbor generation ─────────────────────────────────────────────
static void gen_local_A(int* ix, int* ib, int NP, int NB) {
    for(int n=0;n<3;n++) for(int jb=0;jb<NB;jb++) for(int jc=0;jc<NP;jc++){
        int nb = (n==0)?jc : (n==1)?((jc>0)?jc-1:jc) : ((jc<NP-1)?jc+1:jc);
        ix[IMAP_A(jc,jb,n,NP,NB)]=nb; ib[IMAP_A(jc,jb,n,NP,NB)]=jb;}
}
static void gen_local_B(int* ix, int* ib, int NP, int NB) {
    for(int jb=0;jb<NB;jb++) for(int jc=0;jc<NP;jc++) for(int n=0;n<3;n++){
        int nb = (n==0)?jc : (n==1)?((jc>0)?jc-1:jc) : ((jc<NP-1)?jc+1:jc);
        ix[IMAP_B(n,jc,jb,NP)]=nb; ib[IMAP_B(n,jc,jb,NP)]=jb;}
}

// ─── CPU ref + correctness ─────────────────────────────────────────────────
static void ref_cpu_A(float* out, const float* bln, const float* kin, const int* ix, const int* ib, int NP, int NL, int NB) {
    for(int jb=0;jb<NB;jb++) for(int jk=0;jk<NL;jk++) for(int jc=0;jc<NP;jc++){
        float b0=bln[IBLN_A(jc,0,jb,NP)],b1=bln[IBLN_A(jc,1,jb,NP)],b2=bln[IBLN_A(jc,2,jb,NP)];
        int i0=ix[IMAP_A(jc,jb,0,NP,NB)],bl0=ib[IMAP_A(jc,jb,0,NP,NB)];
        int i1=ix[IMAP_A(jc,jb,1,NP,NB)],bl1=ib[IMAP_A(jc,jb,1,NP,NB)];
        int i2=ix[IMAP_A(jc,jb,2,NP,NB)],bl2=ib[IMAP_A(jc,jb,2,NP,NB)];
        out[IA(jc,jk,jb,NP,NL)]=b0*kin[IA(i0,jk,bl0,NP,NL)]+b1*kin[IA(i1,jk,bl1,NP,NL)]+b2*kin[IA(i2,jk,bl2,NP,NL)];}
}
struct CheckResult{float max_abs,max_rel;int mismatches;};
static CheckResult check_A(const float*g,const float*c,int NP,int NL,int NB,float tol=1e-5f){
    size_t n=(size_t)NP*NL*NB;CheckResult r={0,0,0};
    for(size_t i=0;i<n;i++){float ae=fabsf(g[i]-c[i]),d=fabsf(c[i]),re=(d>1e-10f)?ae/d:0.f;
        if(ae>r.max_abs)r.max_abs=ae;if(re>r.max_rel)r.max_rel=re;if(re>tol)r.mismatches++;}return r;}
static CheckResult check_B_vs_A(const float*gB,const float*cA,int NP,int NL,int NB,float tol=1e-5f){
    CheckResult r={0,0,0};
    for(int jb=0;jb<NB;jb++) for(int jk=0;jk<NL;jk++) for(int jc=0;jc<NP;jc++){
        float gv=gB[IB(jc,jk,jb,NP,NL)],cv=cA[IA(jc,jk,jb,NP,NL)],ae=fabsf(gv-cv),d=fabsf(cv),re=(d>1e-10f)?ae/d:0.f;
        if(ae>r.max_abs)r.max_abs=ae;if(re>r.max_rel)r.max_rel=re;if(re>tol)r.mismatches++;}return r;}

static void verify_correctness(){
    const int NP=1024,NL=90,NB=NBLKS;const float TOL=1e-5f;
    size_t sm=(size_t)NP*NL*NB,sa=(size_t)NP*3*NB;
    float*hbA=new float[sa],*hkA=new float[sm],*hoG=new float[sm],*hoC=new float[sm];
    int*hiA=new int[sa],*hbiA=new int[sa];
    float*hbB=new float[sa],*hkB=new float[sm];int*hiB=new int[sa],*hbiB=new int[sa];float*hoB=new float[sm];
    srand(42);
    for(size_t i=0;i<sa;i++) hbA[i]=(float)rand()/RAND_MAX;
    for(size_t i=0;i<sm;i++) hkA[i]=(float)rand()/RAND_MAX;
    gen_local_A(hiA,hbiA,NP,NB); gen_local_B(hiB,hbiB,NP,NB);
    for(int jb=0;jb<NB;jb++) for(int jc=0;jc<NP;jc++) for(int n=0;n<3;n++)
        hbB[IBLN_B(n,jc,jb,NP)]=hbA[IBLN_A(jc,n,jb,NP)];
    for(int jb=0;jb<NB;jb++) for(int jk=0;jk<NL;jk++) for(int jc=0;jc<NP;jc++)
        hkB[IB(jc,jk,jb,NP,NL)]=hkA[IA(jc,jk,jb,NP,NL)];
    ref_cpu_A(hoC,hbA,hkA,hiA,hbiA,NP,NL,NB);

    float*d_o,*d_b,*d_k;int*d_i,*d_ib;
    CK(cudaMalloc(&d_o,sm*4));CK(cudaMalloc(&d_b,sa*4));CK(cudaMalloc(&d_k,sm*4));CK(cudaMalloc(&d_i,sa*4));CK(cudaMalloc(&d_ib,sa*4));
    CK(cudaMemcpy(d_b,hbA,sa*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(d_k,hkA,sm*4,cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_i,hiA,sa*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(d_ib,hbiA,sa*4,cudaMemcpyHostToDevice));
    ker_A<1,128><<<dim3((NP+127)/128,NL,NB),128>>>(d_o,d_b,d_k,d_i,d_ib,NP,NL,NB);
    CK(cudaDeviceSynchronize());CK(cudaMemcpy(hoG,d_o,sm*4,cudaMemcpyDeviceToHost));
    auto rA=check_A(hoG,hoC,NP,NL,NB,TOL);
    printf("  ker_A vs CPU: max_abs=%.3e max_rel=%.3e mismatches=%d %s\n",rA.max_abs,rA.max_rel,rA.mismatches,rA.mismatches==0?"PASS":"FAIL");
    CK(cudaFree(d_o));CK(cudaFree(d_b));CK(cudaFree(d_k));CK(cudaFree(d_i));CK(cudaFree(d_ib));

    CK(cudaMalloc(&d_o,sm*4));CK(cudaMalloc(&d_b,sa*4));CK(cudaMalloc(&d_k,sm*4));CK(cudaMalloc(&d_i,sa*4));CK(cudaMalloc(&d_ib,sa*4));
    CK(cudaMemcpy(d_b,hbB,sa*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(d_k,hkB,sm*4,cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_i,hiB,sa*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(d_ib,hbiB,sa*4,cudaMemcpyHostToDevice));
    ker_B<1,32,8><<<dim3((NL+31)/32,(NP+7)/8,NB),dim3(32,8,1)>>>(d_o,d_b,d_k,d_i,d_ib,NP,NL,NB);
    CK(cudaDeviceSynchronize());CK(cudaMemcpy(hoB,d_o,sm*4,cudaMemcpyDeviceToHost));
    auto rB=check_B_vs_A(hoB,hoC,NP,NL,NB,TOL);
    printf("  ker_B vs CPU: max_abs=%.3e max_rel=%.3e mismatches=%d %s\n",rB.max_abs,rB.max_rel,rB.mismatches,rB.mismatches==0?"PASS":"FAIL");
    CK(cudaFree(d_o));CK(cudaFree(d_b));CK(cudaFree(d_k));CK(cudaFree(d_i));CK(cudaFree(d_ib));
    delete[]hbA;delete[]hkA;delete[]hiA;delete[]hbiA;delete[]hoG;delete[]hoC;
    delete[]hbB;delete[]hkB;delete[]hiB;delete[]hbiB;delete[]hoB;
}

// ─── CSV ───────────────────────────────────────────────────────────────────
struct TimingRecord{std::string variant;int nlev,iter;float time_ms;};
static std::vector<TimingRecord> g_records;
static void record(const char* label,int nlev,const float* times){
    float s=0; for(int i=0;i<NITER;i++){g_records.push_back({label,nlev,i,times[i]});s+=times[i];}
    printf("  %-52s  avg=%8.3f ms\n",label,s/NITER);
}
static void write_csv(const char* path){
    FILE*f=fopen(path,"w"); fprintf(f,"variant,nlev,iteration,time_ms\n");
    for(auto&r:g_records) fprintf(f,"%s,%d,%d,%.6f\n",r.variant.c_str(),r.nlev,r.iter,r.time_ms);
    fclose(f); printf("Wrote %zu records to %s\n",g_records.size(),path);
}

// ─── Bench ─────────────────────────────────────────────────────────────────
struct Arrays{float*out_A,*kin_A,*bln_A;int*ix_A,*ib_A;float*out_B,*kin_B,*bln_B;int*ix_B,*ib_B;};

template<int ELEMS,int BSZJC,int BSZJK=1,int ELEMC=1>
static void bench_A(const char*label,Arrays&a,int NP,int NL){
    dim3 block(BSZJC,BSZJK,1),grid((NP+BSZJC*ELEMC-1)/(BSZJC*ELEMC),(NL+BSZJK*ELEMS-1)/(BSZJK*ELEMS),NBLKS);
    ker_A<ELEMS,BSZJC,BSZJK,ELEMC><<<grid,block>>>(a.out_A,a.bln_A,a.kin_A,a.ix_A,a.ib_A,NP,NL,NBLKS);
    CK(cudaDeviceSynchronize());
    float times[NITER];cudaEvent_t e0,e1;CK(cudaEventCreate(&e0));CK(cudaEventCreate(&e1));
    for(int i=0;i<NITER;i++){CK(cudaEventRecord(e0));
        ker_A<ELEMS,BSZJC,BSZJK,ELEMC><<<grid,block>>>(a.out_A,a.bln_A,a.kin_A,a.ix_A,a.ib_A,NP,NL,NBLKS);
        CK(cudaEventRecord(e1));CK(cudaEventSynchronize(e1));CK(cudaEventElapsedTime(&times[i],e0,e1));}
    CK(cudaEventDestroy(e0));CK(cudaEventDestroy(e1));record(label,NL,times);
}
template<int ELEMS,int BSZJK,int BSZJC>
static void bench_B(const char*label,Arrays&a,int NP,int NL){
    dim3 block(BSZJK,BSZJC,1),grid((NL+BSZJK*ELEMS-1)/(BSZJK*ELEMS),(NP+BSZJC-1)/BSZJC,NBLKS);
    ker_B<ELEMS,BSZJK,BSZJC><<<grid,block>>>(a.out_B,a.bln_B,a.kin_B,a.ix_B,a.ib_B,NP,NL,NBLKS);
    CK(cudaDeviceSynchronize());
    float times[NITER];cudaEvent_t e0,e1;CK(cudaEventCreate(&e0));CK(cudaEventCreate(&e1));
    for(int i=0;i<NITER;i++){CK(cudaEventRecord(e0));
        ker_B<ELEMS,BSZJK,BSZJC><<<grid,block>>>(a.out_B,a.bln_B,a.kin_B,a.ix_B,a.ib_B,NP,NL,NBLKS);
        CK(cudaEventRecord(e1));CK(cudaEventSynchronize(e1));CK(cudaEventElapsedTime(&times[i],e0,e1));}
    CK(cudaEventDestroy(e0));CK(cudaEventDestroy(e1));record(label,NL,times);
}
template<int ELEMS, int BSZJC, int BSZJK = 1, int ELEMC = 1>
static void bench_C(const char* label, Arrays& a, int NP, int NL) {
    dim3 block(BSZJC, BSZJK, 1);
    dim3 grid((NP+BSZJC*ELEMC-1)/(BSZJC*ELEMC), (NL+BSZJK*ELEMS-1)/(BSZJK*ELEMS), NBLKS);
    // Uses A's out/kin/bln but B's ix/ib
    ker_C<ELEMS,BSZJC,BSZJK,ELEMC><<<grid,block>>>(a.out_A,a.bln_A,a.kin_A,a.ix_B,a.ib_B,NP,NL,NBLKS);
    CK(cudaDeviceSynchronize());
    float times[NITER]; cudaEvent_t e0,e1;
    CK(cudaEventCreate(&e0)); CK(cudaEventCreate(&e1));
    for (int i=0;i<NITER;i++) {
        CK(cudaEventRecord(e0));
        ker_C<ELEMS,BSZJC,BSZJK,ELEMC><<<grid,block>>>(a.out_A,a.bln_A,a.kin_A,a.ix_B,a.ib_B,NP,NL,NBLKS);
        CK(cudaEventRecord(e1)); CK(cudaEventSynchronize(e1));
        CK(cudaEventElapsedTime(&times[i],e0,e1));
    }
    CK(cudaEventDestroy(e0)); CK(cudaEventDestroy(e1));
    record(label, NL, times);
}
int main(int argc,char**argv){
    const char*csv_path="timings_local.csv"; if(argc>1) csv_path=argv[1];
    int dev;CK(cudaGetDevice(&dev));cudaDeviceProp prop;CK(cudaGetDeviceProperties(&prop,dev));
    printf("GPU: %s (SM %d.%d), NITER=%d\n\n",prop.name,prop.major,prop.minor,NITER);

    printf("══ CORRECTNESS (local neighbors) ═══════════════════════════════════\n");
    verify_correctness(); printf("\n");

    const int NP=81920, NL_CONFIGS[]={90,96};
    for(int ci=0;ci<2;ci++){
        int NL=NL_CONFIGS[ci];
        size_t sm=(size_t)NP*NL*NBLKS, sa=(size_t)NP*3*NBLKS;

        float*hbA=new float[sa],*hkA=new float[sm]; int*hiA=new int[sa],*hbiA=new int[sa];
        float*hbB=new float[sa],*hkB=new float[sm]; int*hiB=new int[sa],*hbiB=new int[sa];
        srand(0xDEADBEEF);
        for(size_t i=0;i<sa;i++) hbA[i]=(float)rand()/RAND_MAX;
        for(size_t i=0;i<sm;i++) hkA[i]=(float)rand()/RAND_MAX;
        gen_local_A(hiA,hbiA,NP,NBLKS); gen_local_B(hiB,hbiB,NP,NBLKS);
        for(int jb=0;jb<NBLKS;jb++) for(int jc=0;jc<NP;jc++) for(int n=0;n<3;n++)
            hbB[IBLN_B(n,jc,jb,NP)]=hbA[IBLN_A(jc,n,jb,NP)];
        for(int jb=0;jb<NBLKS;jb++) for(int jk=0;jk<NL;jk++) for(int jc=0;jc<NP;jc++)
            hkB[IB(jc,jk,jb,NP,NL)]=hkA[IA(jc,jk,jb,NP,NL)];

        Arrays a;
        CK(cudaMalloc(&a.out_A,sm*4));CK(cudaMalloc(&a.kin_A,sm*4));CK(cudaMalloc(&a.bln_A,sa*4));
        CK(cudaMalloc(&a.ix_A,sa*4));CK(cudaMalloc(&a.ib_A,sa*4));
        CK(cudaMalloc(&a.out_B,sm*4));CK(cudaMalloc(&a.kin_B,sm*4));CK(cudaMalloc(&a.bln_B,sa*4));
        CK(cudaMalloc(&a.ix_B,sa*4));CK(cudaMalloc(&a.ib_B,sa*4));
        CK(cudaMemcpy(a.bln_A,hbA,sa*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(a.kin_A,hkA,sm*4,cudaMemcpyHostToDevice));
        CK(cudaMemcpy(a.ix_A,hiA,sa*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(a.ib_A,hbiA,sa*4,cudaMemcpyHostToDevice));
        CK(cudaMemcpy(a.bln_B,hbB,sa*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(a.kin_B,hkB,sm*4,cudaMemcpyHostToDevice));
        CK(cudaMemcpy(a.ix_B,hiB,sa*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(a.ib_B,hbiB,sa*4,cudaMemcpyHostToDevice));
        delete[]hbA;delete[]hkA;delete[]hiA;delete[]hbiA;delete[]hbB;delete[]hkB;delete[]hiB;delete[]hbiB;

        printf("══ nlev=%-3d nproma=%-6d nblks=%d (LOCAL jc±1) ══════════════════\n",NL,NP,NBLKS);

        bench_A<1,128>    ("A1-01  ELEMC=1  ELEMS=1  1D blk=128",       a,NP,NL);
        bench_A<1,256>    ("A1-02  ELEMC=1  ELEMS=1  1D blk=256",       a,NP,NL);
        bench_A<1,128,4>  ("A1-03  ELEMC=1  ELEMS=1  2D blk=128x4",     a,NP,NL);
        bench_A<1,32,8>   ("A1-04  ELEMC=1  ELEMS=1  2D blk=32x8",      a,NP,NL);
        bench_A<2,128>    ("A1-05  ELEMC=1  ELEMS=2  1D blk=128",       a,NP,NL);
        bench_A<2,256>    ("A1-06  ELEMC=1  ELEMS=2  1D blk=256",       a,NP,NL);
        bench_A<2,64,4>   ("A1-07  ELEMC=1  ELEMS=2  2D blk=64x4",      a,NP,NL);
        bench_A<2,32,8>   ("A1-08  ELEMC=1  ELEMS=2  2D blk=32x8",      a,NP,NL);
        bench_A<4,128>    ("A1-09  ELEMC=1  ELEMS=4  1D blk=128",       a,NP,NL);
        bench_A<4,256>    ("A1-10  ELEMC=1  ELEMS=4  1D blk=256",       a,NP,NL);
        bench_A<4,64,4>   ("A1-11  ELEMC=1  ELEMS=4  2D blk=64x4",      a,NP,NL);
        bench_A<4,32,8>   ("A1-12  ELEMC=1  ELEMS=4  2D blk=32x8",      a,NP,NL);
        bench_A<1,128,1,2>("A2-03  ELEMS=1  ELEMC=2  1D blk=128",       a,NP,NL);
        bench_A<1,256,1,2>("A2-04  ELEMS=1  ELEMC=2  1D blk=256",       a,NP,NL);
        bench_A<1,64,4,2> ("A2-05  ELEMS=1  ELEMC=2  2D blk=64x4",      a,NP,NL);
        bench_A<1,32,8,2> ("A2-06  ELEMS=1  ELEMC=2  2D blk=32x8",      a,NP,NL);
        bench_A<1,128,1,4>("A2-07  ELEMS=1  ELEMC=4  1D blk=128",       a,NP,NL);
        bench_A<1,256,1,4>("A2-08  ELEMS=1  ELEMC=4  1D blk=256",       a,NP,NL);
        bench_A<1,64,4,4> ("A2-09  ELEMS=1  ELEMC=4  2D blk=64x4",      a,NP,NL);
        bench_A<1,32,8,4> ("A2-10  ELEMS=1  ELEMC=4  2D blk=32x8",      a,NP,NL);
        bench_A<4,128,1,1>("A2-11  ELEMS=4  ELEMC=1  1D blk=128",       a,NP,NL);
        bench_A<4,256,1,1>("A2-12  ELEMS=4  ELEMC=1  1D blk=256",       a,NP,NL);
        bench_A<4,64,4,1> ("A2-13  ELEMS=4  ELEMC=1  2D blk=64x4",      a,NP,NL);
        bench_A<4,32,8,1> ("A2-14  ELEMS=4  ELEMC=1  2D blk=32x8",      a,NP,NL);
        bench_B<1,96,2>   ("B-01  LayoutB  ELEMS=1  blk=96x2",          a,NP,NL);
        bench_B<1,32,8>   ("B-02  LayoutB  ELEMS=1  blk=32x8",          a,NP,NL);
        bench_B<1,64,8>   ("B-03  LayoutB  ELEMS=1  blk=64x8",          a,NP,NL);
        bench_B<1,32,16>  ("B-04  LayoutB  ELEMS=1  blk=32x16",         a,NP,NL);
        bench_B<3,32,8>   ("B-05  LayoutB  ELEMS=3  blk=32x8",          a,NP,NL);
        bench_B<3,64,8>   ("B-06  LayoutB  ELEMS=3  blk=64x8",          a,NP,NL);
        bench_B<3,32,16>  ("B-07  LayoutB  ELEMS=3  blk=32x16",         a,NP,NL);
        printf("\n");
        bench_C<1,256>    ("C1-02  ELEMC=1  ELEMS=1  1D blk=256",       a,NP,NL);

        CK(cudaFree(a.out_A));CK(cudaFree(a.kin_A));CK(cudaFree(a.bln_A));CK(cudaFree(a.ix_A));CK(cudaFree(a.ib_A));
        CK(cudaFree(a.out_B));CK(cudaFree(a.kin_B));CK(cudaFree(a.bln_B));CK(cudaFree(a.ix_B));CK(cudaFree(a.ib_B));
    }
    write_csv(csv_path); return 0;
}
