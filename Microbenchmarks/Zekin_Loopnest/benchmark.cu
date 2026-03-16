// benchmark.cu — ICON z_ekinh kernel variants
// z_ekinh(jc,jk,jb) = sum_{n=0..2} e_bln_c_s(...,n,...) * z_kin_hor_e(ieidx(...,n,...), jk, ieblk(...,n,...))
// All arrays column-major (Fortran style), 0-based indexing in CUDA.

#include <cuda_runtime.h>
#include <curand.h>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

static const int NBLKS = 2;
static const int NITER = 20;

#define CK(e)  do { cudaError_t _r=(e);      if(_r!=cudaSuccess)          { fprintf(stderr,"CUDA %s @L%d: %s\n",       #e,__LINE__,cudaGetErrorString(_r)); exit(1); } } while(0)
#define CRK(e) do { curandStatus_t _r=(e);   if(_r!=CURAND_STATUS_SUCCESS){ fprintf(stderr,"cuRAND %s @L%d: err %d\n", #e,__LINE__,(int)_r);              exit(1); } } while(0)

// ─── index macros ──────────────────────────────────────────────────────────
// Main arrays
//   Layout-A : [nproma, nlev, nblks]   nproma-major (original Fortran order)
#define IA(jc,jk,jb,NP,NL)      ((jc)  + (long)(jk)*(NP)    + (long)(jb)*(NP)*(NL))
//   Layout-B : [nlev, nproma, nblks]   nlev-major (loop-exchanged)
#define IB(jc,jk,jb,NP,NL)      ((jk)  + (long)(jc)*(NL)    + (long)(jb)*(NL)*(NP))

// ker_A auxiliary arrays — match original Fortran layout
//   e_bln_c_s : [nproma, 3, nblks]
#define IBLN_A(jc,n,jb,NP)      ((jc)  + (long)(n)*(NP)     + (long)(jb)*(NP)*3)
//   ieidx/ieblk : [nproma, nblks, 3]
#define IMAP_A(jc,jb,n,NP,NB)   ((jc)  + (long)(jb)*(NP)    + (long)(n)*(NP)*(NB))

// ker_B auxiliary arrays — permuted: n is innermost (stride-1 access for 0,1,2 per column)
//   e_bln_c_s : [3, nproma, nblks]
#define IBLN_B(n,jc,jb,NP)      ((n)   + 3L*(jc)            + 3L*(NP)*(jb))
//   ieidx/ieblk : [3, nproma, nblks]
#define IMAP_B(n,jc,jb,NP)      ((n)   + 3L*(jc)            + 3L*(NP)*(jb))

// ─── helper: uniform float [0,1) → bounded int ─────────────────────────────
__global__ void float_to_int(int* __restrict__ dst, const float* __restrict__ src,
                              int n, int maxval) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) dst[i] = min((int)(src[i] * maxval), maxval - 1);
}

// ═══════════════════════════════════════════════════════════════════════════
// ker_A — layout [nproma, nlev, nblks], aux in [nproma,3,nblks]/[nproma,nblks,3]
//
//   Block: (BSZJC, BSZJK, 1)   threadIdx.x → jc  (fast/coalesced dim)
//                               threadIdx.y → jk  (BSZJK=1 → 1-D block)
//   Grid : (⌈NP/(BSZJC·ELEMC)⌉ , ⌈NL/(BSZJK·ELEMS)⌉ , NBLKS)
//
//   Each thread: covers ELEMC consecutive jc columns × ELEMS consecutive jk levels.
//   ELEMC=1, BSZJK=1 → original 1-D block, 1-element-per-thread behaviour.
//
//   jc layout: thread t covers jc = t*ELEMC, t*ELEMC+1, ..., t*ELEMC+ELEMC-1
//   (consecutive → coalesced reads/writes since nproma is stride-1 in layout A)
// ═══════════════════════════════════════════════════════════════════════════
template<int ELEMS, int BSZJC, int BSZJK = 1, int ELEMC = 1>
__global__ void ker_A(float* __restrict__ out,
                      const float* __restrict__ bln,
                      const float* __restrict__ kin,
                      const int*   __restrict__ ix,
                      const int*   __restrict__ ib,
                      int NP, int NL, int NB)
{
    // Each thread's base jc: threads tile in groups of ELEMC consecutive columns
    int jc_base = (blockIdx.x * BSZJC + threadIdx.x) * ELEMC;
    int jk_base = (blockIdx.y * BSZJK + threadIdx.y) * ELEMS;
    int jb      =  blockIdx.z;
    if (jc_base >= NP || jb >= NB) return;

    // Per-column prefetch for each of the ELEMC jc values.
    // Scalars are reused across all ELEMS jk levels → good register reuse.
    float b0[ELEMC], b1[ELEMC], b2[ELEMC];
    int   i0[ELEMC], i1[ELEMC], i2[ELEMC];
    int  bl0[ELEMC],bl1[ELEMC],bl2[ELEMC];

    #pragma unroll
    for (int c = 0; c < ELEMC; c++) {
        int jc = jc_base + c;
        if (jc >= NP) break;                  // guard last partial tile
        b0[c]  = bln[IBLN_A(jc,0,jb,NP)];
        b1[c]  = bln[IBLN_A(jc,1,jb,NP)];
        b2[c]  = bln[IBLN_A(jc,2,jb,NP)];
        i0[c]  = ix [IMAP_A(jc,jb,0,NP,NB)]; bl0[c] = ib[IMAP_A(jc,jb,0,NP,NB)];
        i1[c]  = ix [IMAP_A(jc,jb,1,NP,NB)]; bl1[c] = ib[IMAP_A(jc,jb,1,NP,NB)];
        i2[c]  = ix [IMAP_A(jc,jb,2,NP,NB)]; bl2[c] = ib[IMAP_A(jc,jb,2,NP,NB)];
    }

    // Inner loops fully unrolled: ELEMC jc values × ELEMS jk levels
    #pragma unroll
    for (int e = 0; e < ELEMS; e++) {
        if (jk_base + e >= NL) break;
        #pragma unroll
        for (int c = 0; c < ELEMC; c++) {
            int jc = jc_base + c;
            if (jc >= NP) break;
            out[IA(jc, jk_base+e, jb, NP, NL)] =
                b0[c] * kin[IA(i0[c], jk_base+e, bl0[c], NP, NL)] +
                b1[c] * kin[IA(i1[c], jk_base+e, bl1[c], NP, NL)] +
                b2[c] * kin[IA(i2[c], jk_base+e, bl2[c], NP, NL)];
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// ker_B — layout [nlev, nproma, nblks], aux in [3, nproma, nblks]
//
//   Block: (BSZJK, BSZJC, 1)   threadIdx.x → jk,  threadIdx.y → jc
//   Grid : (⌈NL/(BSZJK·ELEMS)⌉ , ⌈NP/BSZJC⌉ , NBLKS)
//
//   Each thread: covers 1 jc (from threadIdx.y) × ELEMS consecutive jk levels.
//   Loading b0/b1/b2 and indices hits [3,jc,jb] — n=0,1,2 are stride-1.
// ═══════════════════════════════════════════════════════════════════════════
template<int ELEMS, int BSZJK, int BSZJC>
__global__ void ker_B(float* __restrict__ out,
                      const float* __restrict__ bln,
                      const float* __restrict__ kin,
                      const int*   __restrict__ ix,
                      const int*   __restrict__ ib,
                      int NP, int NL, int NB)
{
    int jk = (blockIdx.x * BSZJK + threadIdx.x) * ELEMS;   // base level for this thread
    int jc =  blockIdx.y * BSZJC + threadIdx.y;             // column index
    int jb =  blockIdx.z;
    if (jc >= NP || jb >= NB) return;

    // n=0,1,2 are stride-1 in IBLN_B / IMAP_B → 3 consecutive loads per warp-column
    float b0 = bln[IBLN_B(0, jc, jb, NP)];
    float b1 = bln[IBLN_B(1, jc, jb, NP)];
    float b2 = bln[IBLN_B(2, jc, jb, NP)];
    int i0  = ix[IMAP_B(0, jc, jb, NP)],  bl0 = ib[IMAP_B(0, jc, jb, NP)];
    int i1  = ix[IMAP_B(1, jc, jb, NP)],  bl1 = ib[IMAP_B(1, jc, jb, NP)];
    int i2  = ix[IMAP_B(2, jc, jb, NP)],  bl2 = ib[IMAP_B(2, jc, jb, NP)];

    #pragma unroll
    for (int e = 0; e < ELEMS; e++) {
        if (jk + e >= NL) break;
        out[IB(jc, jk+e, jb, NP, NL)] =
            b0 * kin[IB(i0, jk+e, bl0, NP, NL)] +
            b1 * kin[IB(i1, jk+e, bl1, NP, NL)] +
            b2 * kin[IB(i2, jk+e, bl2, NP, NL)];
    }
}

// ─── timing ────────────────────────────────────────────────────────────────
static cudaEvent_t ev0, ev1;
static void tstart() { CK(cudaEventRecord(ev0)); }
static float tstop()  {
    CK(cudaEventRecord(ev1)); CK(cudaEventSynchronize(ev1));
    float ms; CK(cudaEventElapsedTime(&ms, ev0, ev1));
    return ms / NITER;
}

// ─── arrays ────────────────────────────────────────────────────────────────
struct Arrays {
    float *out_A, *kin_A, *bln_A;   int *ix_A, *ib_A;  // Layout-A with A-aux
    float *out_B, *kin_B, *bln_B;   int *ix_B, *ib_B;  // Layout-B with B-aux
};

// ─── bench helpers ─────────────────────────────────────────────────────────
// BSZJK=1, ELEMC=1 → original 1-D, 1-elem/thr behaviour
// BSZJK>1           → 2-D block, threadIdx.y covers jk dimension
// ELEMC>1           → each thread computes ELEMC consecutive jc columns
template<int ELEMS, int BSZJC, int BSZJK = 1, int ELEMC = 1>
static void bench_A(const char* label, Arrays& a, int NP, int NL) {
    dim3 block(BSZJC, BSZJK, 1);
    dim3 grid((NP + BSZJC*ELEMC - 1) / (BSZJC*ELEMC),   // each block covers BSZJC*ELEMC columns
              (NL + BSZJK*ELEMS - 1) / (BSZJK*ELEMS),   // each block covers BSZJK*ELEMS levels
              NBLKS);
    ker_A<ELEMS,BSZJC,BSZJK,ELEMC><<<grid,block>>>(a.out_A, a.bln_A, a.kin_A, a.ix_A, a.ib_A, NP, NL, NBLKS);
    CK(cudaDeviceSynchronize());
    tstart();
    for (int i = 0; i < NITER; i++)
        ker_A<ELEMS,BSZJC,BSZJK,ELEMC><<<grid,block>>>(a.out_A, a.bln_A, a.kin_A, a.ix_A, a.ib_A, NP, NL, NBLKS);
    printf("  %-52s  %8.3f ms\n", label, tstop());
}

template<int ELEMS, int BSZJK, int BSZJC>
static void bench_B(const char* label, Arrays& a, int NP, int NL) {
    // x → ⌈NL / (BSZJK * ELEMS)⌉  blocks covering levels
    // y → ⌈NP / BSZJC⌉             blocks covering columns
    // z → NBLKS
    dim3 block(BSZJK, BSZJC, 1);
    dim3 grid((NL + BSZJK*ELEMS - 1) / (BSZJK*ELEMS),
              (NP + BSZJC       - 1) /  BSZJC,
              NBLKS);
    ker_B<ELEMS,BSZJK,BSZJC><<<grid,block>>>(a.out_B, a.bln_B, a.kin_B, a.ix_B, a.ib_B, NP, NL, NBLKS);
    CK(cudaDeviceSynchronize());
    tstart();
    for (int i = 0; i < NITER; i++)
        ker_B<ELEMS,BSZJK,BSZJC><<<grid,block>>>(a.out_B, a.bln_B, a.kin_B, a.ix_B, a.ib_B, NP, NL, NBLKS);
    printf("  %-52s  %8.3f ms\n", label, tstop());
}

// ─── main ──────────────────────────────────────────────────────────────────
int main()
{
    int dev; CK(cudaGetDevice(&dev));
    cudaDeviceProp prop; CK(cudaGetDeviceProperties(&prop, dev));
    printf("GPU: %s  (SM %d.%d)\n\n", prop.name, prop.major, prop.minor);

    CK(cudaEventCreate(&ev0)); CK(cudaEventCreate(&ev1));

    const int NP = 81920;
    const int NL_CONFIGS[] = {90, 96};

    for (int ci = 0; ci < 2; ci++) {
        int NL = NL_CONFIGS[ci];

        size_t sz_main = (size_t)NP * NL * NBLKS;   // main/kin arrays
        size_t sz_aux  = (size_t)NP * 3  * NBLKS;   // bln/ix/ib (same total, different layout)

        Arrays a;
        CK(cudaMalloc(&a.out_A, sz_main * sizeof(float)));
        CK(cudaMalloc(&a.kin_A, sz_main * sizeof(float)));
        CK(cudaMalloc(&a.bln_A, sz_aux  * sizeof(float)));
        CK(cudaMalloc(&a.ix_A,  sz_aux  * sizeof(int)));
        CK(cudaMalloc(&a.ib_A,  sz_aux  * sizeof(int)));
        CK(cudaMalloc(&a.out_B, sz_main * sizeof(float)));
        CK(cudaMalloc(&a.kin_B, sz_main * sizeof(float)));
        CK(cudaMalloc(&a.bln_B, sz_aux  * sizeof(float)));
        CK(cudaMalloc(&a.ix_B,  sz_aux  * sizeof(int)));
        CK(cudaMalloc(&a.ib_B,  sz_aux  * sizeof(int)));

        // ── cuRAND XORWOW — fully parallel on-GPU RNG ─────────────
        curandGenerator_t gen;
        CRK(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_XORWOW));
        CRK(curandSetPseudoRandomGeneratorSeed(gen, 0xDEADBEEFULL));

        CRK(curandGenerateUniform(gen, a.bln_A, sz_aux));
        CRK(curandGenerateUniform(gen, a.kin_A, sz_main));
        CRK(curandGenerateUniform(gen, a.bln_B, sz_aux));
        CRK(curandGenerateUniform(gen, a.kin_B, sz_main));

        float* tmp; CK(cudaMalloc(&tmp, sz_aux * sizeof(float)));

        CRK(curandGenerateUniform(gen, tmp, sz_aux));
        float_to_int<<<(sz_aux+255)/256, 256>>>(a.ix_A, tmp, sz_aux, NP);  // [0, NP-1]
        CRK(curandGenerateUniform(gen, tmp, sz_aux));
        float_to_int<<<(sz_aux+255)/256, 256>>>(a.ix_B, tmp, sz_aux, NP);

        CRK(curandGenerateUniform(gen, tmp, sz_aux));
        float_to_int<<<(sz_aux+255)/256, 256>>>(a.ib_A, tmp, sz_aux, 2);   // {0, 1}
        CRK(curandGenerateUniform(gen, tmp, sz_aux));
        float_to_int<<<(sz_aux+255)/256, 256>>>(a.ib_B, tmp, sz_aux, 2);

        CK(cudaDeviceSynchronize());
        CK(cudaFree(tmp));
        CRK(curandDestroyGenerator(gen));

        printf("══════════════════════════════════════════════════════════════════════\n");
        printf("  nlev=%-3d  nproma=%-6d  nblks=%d\n", NL, NP, NBLKS);
        printf("══════════════════════════════════════════════════════════════════════\n");
        printf("  %-52s  %8s\n", "Variant", "Avg ms");
        printf("  %s\n", "──────────────────────────────────────────────────────────────");

        // ════════════════════════════════════════════════════════════
        // Layout-A  —  <ELEMS, BSZJC, BSZJK, ELEMC>
        // ════════════════════════════════════════════════════════════

        // ── A1: vary ELEMS (nlev/jk per thread), ELEMC=1 fixed ───────────
        //        isolates the effect of processing more levels per thread
        printf("  --- A1: ELEMC=1 fixed, ELEMS varies ---\n");
        bench_A<1,128>    ("A1-01  ELEMC=1  ELEMS=1  1D blk=128",       a, NP, NL);
        bench_A<1,256>    ("A1-02  ELEMC=1  ELEMS=1  1D blk=256",       a, NP, NL);
        bench_A<1,128, 4> ("A1-03  ELEMC=1  ELEMS=1  2D blk=128x4",     a, NP, NL);
        bench_A<1, 32, 8> ("A1-04  ELEMC=1  ELEMS=1  2D blk=32x8",      a, NP, NL);

        bench_A<2,128>    ("A1-05  ELEMC=1  ELEMS=2  1D blk=128",       a, NP, NL);
        bench_A<2,256>    ("A1-06  ELEMC=1  ELEMS=2  1D blk=256",       a, NP, NL);
        bench_A<2, 64, 4> ("A1-07  ELEMC=1  ELEMS=2  2D blk=64x4",      a, NP, NL);
        bench_A<2, 32, 8> ("A1-08  ELEMC=1  ELEMS=2  2D blk=32x8",      a, NP, NL);

        bench_A<4,128>    ("A1-09  ELEMC=1  ELEMS=4  1D blk=128",       a, NP, NL);
        bench_A<4,256>    ("A1-10  ELEMC=1  ELEMS=4  1D blk=256",       a, NP, NL);
        bench_A<4, 64, 4> ("A1-11  ELEMC=1  ELEMS=4  2D blk=64x4",      a, NP, NL);
        bench_A<4, 32, 8> ("A1-12  ELEMC=1  ELEMS=4  2D blk=32x8",      a, NP, NL);

        // ── A2: vary ELEMC (nproma/jc per thread), ELEMS=1 fixed ─────────
        //        isolates the effect of processing more columns per thread
        printf("  --- A2: ELEMS=1 fixed, ELEMC varies ---\n");
        bench_A<1,128,1,1>("A2-01  ELEMS=1  ELEMC=1  1D blk=128",       a, NP, NL);
        bench_A<1,256,1,1>("A2-02  ELEMS=1  ELEMC=1  1D blk=256",       a, NP, NL);

        bench_A<1,128,1,2>("A2-03  ELEMS=1  ELEMC=2  1D blk=128",       a, NP, NL);
        bench_A<1,256,1,2>("A2-04  ELEMS=1  ELEMC=2  1D blk=256",       a, NP, NL);
        bench_A<1, 64,4,2>("A2-05  ELEMS=1  ELEMC=2  2D blk=64x4",      a, NP, NL);
        bench_A<1, 32,8,2>("A2-06  ELEMS=1  ELEMC=2  2D blk=32x8",      a, NP, NL);

        bench_A<1,128,1,4>("A2-07  ELEMS=1  ELEMC=4  1D blk=128",       a, NP, NL);
        bench_A<1,256,1,4>("A2-08  ELEMS=1  ELEMC=4  1D blk=256",       a, NP, NL);
        bench_A<1, 64,4,4>("A2-09  ELEMS=1  ELEMC=4  2D blk=64x4",      a, NP, NL);
        bench_A<1, 32,8,4>("A2-10  ELEMS=1  ELEMC=4  2D blk=32x8",      a, NP, NL);

        bench_A<4,128,1,1>("A2-07  ELEMS=4  ELEMC=1  1D blk=128",       a, NP, NL);
        bench_A<4,256,1,1>("A2-08  ELEMS=4  ELEMC=1  1D blk=256",       a, NP, NL);
        bench_A<4, 64,4,1>("A2-09  ELEMS=4  ELEMC=1  2D blk=64x4",      a, NP, NL);
        bench_A<4, 32,8,1>("A2-10  ELEMS=4  ELEMC=1  2D blk=32x8",      a, NP, NL);

        // ── ker_B  (loop-exchanged reference, unchanged) ──────────────────
        printf("  --- B: LayoutB [nlev,nproma,nblks] reference ---\n");
        bench_B<1,96, 2>("B-01  LayoutB  ELEMS=1  blk=96x1",            a, NP, NL);
        bench_B<1,32, 8>("B-01  LayoutB  ELEMS=1  blk=32x8",            a, NP, NL);
        bench_B<1,64, 8>("B-02  LayoutB  ELEMS=1  blk=64x8",            a, NP, NL);
        bench_B<1,32,16>("B-03  LayoutB  ELEMS=1  blk=32x16",           a, NP, NL);
        bench_B<3,32, 8>("B-04  LayoutB  ELEMS=3  blk=32x8",            a, NP, NL);
        bench_B<3,64, 8>("B-05  LayoutB  ELEMS=3  blk=64x8",            a, NP, NL);
        bench_B<3,32,16>("B-06  LayoutB  ELEMS=3  blk=32x16",           a, NP, NL);
        printf("\n");

        CK(cudaFree(a.out_A)); CK(cudaFree(a.kin_A)); CK(cudaFree(a.bln_A));
        CK(cudaFree(a.ix_A));  CK(cudaFree(a.ib_A));
        CK(cudaFree(a.out_B)); CK(cudaFree(a.kin_B)); CK(cudaFree(a.bln_B));
        CK(cudaFree(a.ix_B));  CK(cudaFree(a.ib_B));
    }

    CK(cudaEventDestroy(ev0)); CK(cudaEventDestroy(ev1));
    return 0;
}
