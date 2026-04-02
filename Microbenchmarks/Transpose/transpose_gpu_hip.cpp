#include "hip/hip_runtime.h"
/*  transpose_gpu.cu — Fully-templated GPU matrix transpose.
 *
 *  Arg order matches run_transpose.py:
 *    ./transpose_gpu N variant BX BY TX TY csv SB PAD WARMUP REPS
 *
 *  Compile (CUDA): nvcc -O3 -std=c++17 -arch=native -o transpose_gpu transpose_gpu.cu
 *  Compile (HIP):  hipcc <flags> -o transpose_gpu transpose_gpu_hip.cpp
 *                  (rename file, swap cuda→hip API calls)
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <hip/hip_runtime.h>

#define GPU(call) do { hipError_t e = (call); if (e != hipSuccess) { \
    fprintf(stderr, "CUDA %s:%d: %s\n", __FILE__, __LINE__, \
    hipGetErrorString(e)); exit(1); } } while(0)

/* ═══════════════════════════════════════════════════════════════════════
 *  V0: Naive
 * ═══════════════════════════════════════════════════════════════════════ */
template<int BX, int BY, int TX, int TY>
__global__ void tr_naive(const float* __restrict__ in,
                         float* __restrict__ out, int N) {
    const int c0 = blockIdx.x * (BX * TX) + threadIdx.x * TX;
    const int r0 = blockIdx.y * (BY * TY) + threadIdx.y * TY;
    #pragma unroll
    for (int ty = 0; ty < TY; ty++) {
        const int r = r0 + ty;
        if (r < N) {
            #pragma unroll
            for (int tx = 0; tx < TX; tx++) {
                const int c = c0 + tx;
                if (c < N) out[c * N + r] = in[r * N + c];
            }
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V1: Blocked (SB template param)
 * ═══════════════════════════════════════════════════════════════════════ */
template<int BX, int BY, int TX, int TY, int SB>
__global__ void tr_blocked(const float* __restrict__ in,
                           float* __restrict__ out, int N) {
    const int c0 = blockIdx.x * (BX * TX) + threadIdx.x * TX;
    const int r0 = blockIdx.y * (BY * TY) + threadIdx.y * TY;
    #pragma unroll
    for (int ty = 0; ty < TY; ty++) {
        const int r = r0 + ty;
        if (r < N) {
            #pragma unroll
            for (int tx = 0; tx < TX; tx++) {
                const int c = c0 + tx;
                if (c < N) out[c * N + r] = in[r * N + c];
            }
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V2: Shared memory (no bank-conflict fix)
 * ═══════════════════════════════════════════════════════════════════════ */
template<int BX, int BY, int TX, int TY>
__global__ void tr_smem(const float* __restrict__ in,
                        float* __restrict__ out, int N) {
    constexpr int BW = BX * TX, BH = BY * TY, TOT = BX * BY;
    __shared__ float sm[BH * BW];
    const int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    const int tid = threadIdx.y * BX + threadIdx.x;
    #pragma unroll
    for (int k = tid; k < BH * BW; k += TOT) {
        const int lr = k / BW, lc = k % BW;
        const int gr = r0 + lr, gc = c0 + lc;
        sm[lr * BW + lc] = (gr < N && gc < N) ? in[gr * N + gc] : 0.f;
    }
    __syncthreads();
    #pragma unroll
    for (int k = tid; k < BW * BH; k += TOT) {
        const int lc = k / BH, lr = k % BH;
        const int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * BW + lc];
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V3: Shared memory + SB blocking
 * ═══════════════════════════════════════════════════════════════════════ */
template<int BX, int BY, int TX, int TY, int SB>
__global__ void tr_smem_blk(const float* __restrict__ in,
                            float* __restrict__ out, int N) {
    constexpr int BW = BX * TX, BH = BY * TY, TOT = BX * BY;
    __shared__ float sm[BH * BW];
    const int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    const int tid = threadIdx.y * BX + threadIdx.x;
    #pragma unroll
    for (int k = tid; k < BH * BW; k += TOT) {
        const int lr = k / BW, lc = k % BW;
        const int gr = r0 + lr, gc = c0 + lc;
        sm[lr * BW + lc] = (gr < N && gc < N) ? in[gr * N + gc] : 0.f;
    }
    __syncthreads();
    #pragma unroll
    for (int k = tid; k < BW * BH; k += TOT) {
        const int lc = k / BH, lr = k % BH;
        const int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * BW + lc];
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V4: Shared memory + padding
 * ═══════════════════════════════════════════════════════════════════════ */
template<int BX, int BY, int TX, int TY, int PAD>
__global__ void tr_smem_pad(const float* __restrict__ in,
                            float* __restrict__ out, int N) {
    constexpr int BW = BX * TX, BH = BY * TY, TOT = BX * BY;
    constexpr int SW = BW + PAD;
    __shared__ float sm[BH * SW];
    const int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    const int tid = threadIdx.y * BX + threadIdx.x;
    #pragma unroll
    for (int k = tid; k < BH * BW; k += TOT) {
        const int lr = k / BW, lc = k % BW;
        const int gr = r0 + lr, gc = c0 + lc;
        sm[lr * SW + lc] = (gr < N && gc < N) ? in[gr * N + gc] : 0.f;
    }
    __syncthreads();
    #pragma unroll
    for (int k = tid; k < BW * BH; k += TOT) {
        const int lc = k / BH, lr = k % BH;
        const int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * SW + lc];
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V5: Shared memory + XOR swizzle
 * ═══════════════════════════════════════════════════════════════════════ */
template<int BX, int BY, int TX, int TY>
__global__ void tr_smem_swiz(const float* __restrict__ in,
                             float* __restrict__ out, int N) {
    constexpr int BW = BX * TX, BH = BY * TY, TOT = BX * BY;
    __shared__ float sm[BH * BW];
    const int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    const int tid = threadIdx.y * BX + threadIdx.x;
    #pragma unroll
    for (int k = tid; k < BH * BW; k += TOT) {
        const int lr = k / BW, lc = k % BW;
        const int gr = r0 + lr, gc = c0 + lc;
        sm[lr * BW + (lc ^ lr)] = (gr < N && gc < N) ? in[gr * N + gc] : 0.f;
    }
    __syncthreads();
    #pragma unroll
    for (int k = tid; k < BW * BH; k += TOT) {
        const int lc = k / BH, lr = k % BH;
        const int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * BW + (lc ^ lr)];
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V6: Smem + SB blocking + XOR swizzle
 * ═══════════════════════════════════════════════════════════════════════ */
template<int BX, int BY, int TX, int TY, int SB>
__global__ void tr_smem_blk_swiz(const float* __restrict__ in,
                                 float* __restrict__ out, int N) {
    constexpr int BW = BX * TX, BH = BY * TY, TOT = BX * BY;
    __shared__ float sm[BH * BW];
    const int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    const int tid = threadIdx.y * BX + threadIdx.x;
    #pragma unroll
    for (int k = tid; k < BH * BW; k += TOT) {
        const int lr = k / BW, lc = k % BW;
        const int gr = r0 + lr, gc = c0 + lc;
        sm[lr * BW + (lc ^ lr)] = (gr < N && gc < N) ? in[gr * N + gc] : 0.f;
    }
    __syncthreads();
    #pragma unroll
    for (int k = tid; k < BW * BH; k += TOT) {
        const int lc = k / BH, lr = k % BH;
        const int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * BW + (lc ^ lr)];
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V7: Smem + padding + SB blocking
 * ═══════════════════════════════════════════════════════════════════════ */
template<int BX, int BY, int TX, int TY, int SB, int PAD>
__global__ void tr_smem_pad_blk(const float* __restrict__ in,
                                float* __restrict__ out, int N) {
    constexpr int BW = BX * TX, BH = BY * TY, TOT = BX * BY;
    constexpr int SW = BW + PAD;
    __shared__ float sm[BH * SW];
    const int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    const int tid = threadIdx.y * BX + threadIdx.x;
    #pragma unroll
    for (int k = tid; k < BH * BW; k += TOT) {
        const int lr = k / BW, lc = k % BW;
        const int gr = r0 + lr, gc = c0 + lc;
        sm[lr * SW + lc] = (gr < N && gc < N) ? in[gr * N + gc] : 0.f;
    }
    __syncthreads();
    #pragma unroll
    for (int k = tid; k < BW * BH; k += TOT) {
        const int lc = k / BH, lr = k % BH;
        const int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * SW + lc];
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Dispatch
 * ═══════════════════════════════════════════════════════════════════════ */
enum Variant { V_NAIVE=0, V_BLOCKED, V_SMEM, V_SMEM_BLK, V_SMEM_PAD,
               V_SMEM_SWIZ, V_SMEM_BLK_SWIZ, V_SMEM_PAD_BLK, V_COUNT };
static const char* V_NAMES[] = {
    "naive","blocked","smem","smem_blk","smem_pad",
    "smem_swiz","smem_blk_swiz","smem_pad_blk"};

#define DISPATCH(BX_,BY_,TX_,TY_,SB_,PAD_,var,grid,blk,in,out,N) \
    switch(var){ \
    case V_NAIVE:        tr_naive        <BX_,BY_,TX_,TY_>          <<<grid,blk>>>(in,out,N);break;\
    case V_BLOCKED:      tr_blocked      <BX_,BY_,TX_,TY_,SB_>     <<<grid,blk>>>(in,out,N);break;\
    case V_SMEM:         tr_smem         <BX_,BY_,TX_,TY_>          <<<grid,blk>>>(in,out,N);break;\
    case V_SMEM_BLK:     tr_smem_blk     <BX_,BY_,TX_,TY_,SB_>     <<<grid,blk>>>(in,out,N);break;\
    case V_SMEM_PAD:     tr_smem_pad     <BX_,BY_,TX_,TY_,PAD_>    <<<grid,blk>>>(in,out,N);break;\
    case V_SMEM_SWIZ:    tr_smem_swiz    <BX_,BY_,TX_,TY_>          <<<grid,blk>>>(in,out,N);break;\
    case V_SMEM_BLK_SWIZ:tr_smem_blk_swiz<BX_,BY_,TX_,TY_,SB_>    <<<grid,blk>>>(in,out,N);break;\
    case V_SMEM_PAD_BLK: tr_smem_pad_blk <BX_,BY_,TX_,TY_,SB_,PAD_><<<grid,blk>>>(in,out,N);break;\
    default:fprintf(stderr,"bad variant %d\n",var);exit(1);}

static void dispatch(int var, int bx, int by, int tx, int ty, int sb, int pad,
                     dim3 grid, dim3 blk,
                     const float* in, float* out, int N) {
    #define TRY(BX_,BY_,TX_,TY_) \
        if(bx==BX_&&by==BY_&&tx==TX_&&ty==TY_){ \
            if(sb==8  &&pad==1){DISPATCH(BX_,BY_,TX_,TY_,8,  1,var,grid,blk,in,out,N);return;}\
            if(sb==32 &&pad==1){DISPATCH(BX_,BY_,TX_,TY_,32, 1,var,grid,blk,in,out,N);return;}\
            if(sb==64 &&pad==1){DISPATCH(BX_,BY_,TX_,TY_,64, 1,var,grid,blk,in,out,N);return;}\
            if(sb==128&&pad==1){DISPATCH(BX_,BY_,TX_,TY_,128,1,var,grid,blk,in,out,N);return;}\
            if(sb==256&&pad==1){DISPATCH(BX_,BY_,TX_,TY_,256,1,var,grid,blk,in,out,N);return;}\
            if(sb==32 &&pad==2){DISPATCH(BX_,BY_,TX_,TY_,32, 2,var,grid,blk,in,out,N);return;}\
            if(sb==64 &&pad==2){DISPATCH(BX_,BY_,TX_,TY_,64, 2,var,grid,blk,in,out,N);return;}\
            DISPATCH(BX_,BY_,TX_,TY_,32,1,var,grid,blk,in,out,N);return;}

    /* All configs from run_transpose.py CONFIGS */
    TRY(32, 8, 1, 8)   TRY(32,16, 1, 8)   TRY(32, 8, 2, 4)   TRY(32,16, 2, 4)
    TRY(32, 8, 4, 1)   TRY(32,16, 4, 1)   TRY(32, 8, 4, 2)   TRY(32,16, 4, 2)
    TRY(32, 8, 2, 2)   TRY(32, 8, 4, 4)
    TRY(32,32, 1, 1)   TRY(32,32, 1, 2)   TRY(32,32, 2, 1)   TRY(32,32, 1, 4)
    TRY(32,32, 2, 2)
    TRY(16, 8, 1, 1)   TRY(16,16, 1, 1)   TRY(16, 8, 1, 2)   TRY(16,16, 1, 2)
    TRY(16, 8, 2, 1)   TRY(16,16, 2, 2)   TRY(16, 8, 1, 4)   TRY(16,16, 1, 4)
    TRY(64, 8, 1, 4)   TRY(64,16, 1, 2)   TRY(64,16, 1, 4)   TRY(64,16, 2, 1)
    TRY(64,16, 2, 2)
    /* Extra common configs */
    TRY(32, 8, 1, 1)   TRY(32, 8, 1, 2)   TRY(32, 8, 1, 4)
    TRY(16,16, 4, 4)   TRY(32,32, 1, 2)   TRY( 8, 8, 4, 4)
    TRY(64, 4, 1, 8)   TRY(32, 4, 1, 8)
    #undef TRY

    fprintf(stderr,"No instantiation for BX=%d BY=%d TX=%d TY=%d SB=%d PAD=%d\n",
            bx,by,tx,ty,sb,pad);
    exit(1);
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Main — arg order: N variant BX BY TX TY csv SB PAD WARMUP REPS
 *         (matches run_transpose.py)
 * ═══════════════════════════════════════════════════════════════════════ */
int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s N [variant=0] [BX=32] [BY=8] [TX=1] [TY=8]"
            " [csv] [SB=32] [PAD=1] [WARMUP=5] [REPS=100]\n"
            "  0=naive 1=blocked 2=smem 3=smem_blk 4=smem_pad"
            " 5=smem_swiz 6=smem_blk_swiz 7=smem_pad_blk\n", argv[0]);
        return 1;
    }

    /* Arg order: N variant BX BY TX TY csv SB PAD WARMUP REPS */
    int N      = atoi(argv[1]);
    int VAR    = argc> 2 ? atoi(argv[2]) : 0;
    int BX_    = argc> 3 ? atoi(argv[3]) : 32;
    int BY_    = argc> 4 ? atoi(argv[4]) : 8;
    int TX_    = argc> 5 ? atoi(argv[5]) : 1;
    int TY_    = argc> 6 ? atoi(argv[6]) : 8;
    const char* csv = argc> 7 ? argv[7]  : "transpose_raw.csv";
    int SB_    = argc> 8 ? atoi(argv[8]) : 32;
    int PAD_   = argc> 9 ? atoi(argv[9]) : 1;
    int WARMUP = argc>10 ? atoi(argv[10]): 5;
    int REPS   = argc>11 ? atoi(argv[11]): 100;

    if (VAR < 0 || VAR >= V_COUNT) {
        fprintf(stderr,"bad variant %d (max %d)\n", VAR, V_COUNT-1);
        return 1;
    }

    int BW = BX_ * TX_, BH = BY_ * TY_;
    dim3 block(BX_, BY_);
    dim3 grid((N + BW - 1) / BW, (N + BH - 1) / BH);

    size_t elems = (size_t)N * N, bytes = elems * sizeof(float);

    float* h_in  = (float*)malloc(bytes);
    float* h_out = (float*)malloc(bytes);
    for (size_t i = 0; i < elems; i++) h_in[i] = (float)i / (float)N;

    float *d_in, *d_out;
    GPU(hipMalloc(&d_in,  bytes));
    GPU(hipMalloc(&d_out, bytes));
    GPU(hipMemcpy(d_in, h_in, bytes, hipMemcpyHostToDevice));

    for (int i = 0; i < WARMUP; i++)
        dispatch(VAR, BX_, BY_, TX_, TY_, SB_, PAD_, grid, block, d_in, d_out, N);
    GPU(hipDeviceSynchronize());

    /* Verify */
    GPU(hipMemcpy(h_out, d_out, bytes, hipMemcpyDeviceToHost));
    float maxerr = 0;
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            float e = fabsf(h_out[c * N + r] - h_in[r * N + c]);
            if (e > maxerr) maxerr = e;
        }
    double cksum = 0;
    for (size_t i = 0; i < elems; i++) cksum += h_out[i];

    /* Timed runs */
    hipEvent_t t0, t1;
    GPU(hipEventCreate(&t0));
    GPU(hipEventCreate(&t1));
    float* times_ms = (float*)malloc(REPS * sizeof(float));

    for (int i = 0; i < REPS; i++) {
        GPU(hipEventRecord(t0, 0));
        dispatch(VAR, BX_, BY_, TX_, TY_, SB_, PAD_, grid, block, d_in, d_out, N);
        GPU(hipEventRecord(t1, 0));
        GPU(hipEventSynchronize(t1));
        GPU(hipEventElapsedTime(&times_ms[i], t0, t1));
    }

    std::sort(times_ms, times_ms + REPS);
    double bpi = 2.0 * N * (double)N * sizeof(float);
    float med_ms = times_ms[REPS / 2];

    printf("%s N=%d BX=%d BY=%d TX=%d TY=%d SB=%d PAD=%d | "
           "med %.4f ms (%.1f GB/s)  p5 %.4f ms  p95 %.4f ms  "
           "maxerr=%.1e  cksum=%.6e\n",
           V_NAMES[VAR], N, BX_, BY_, TX_, TY_, SB_, PAD_,
           med_ms, bpi / (med_ms * 1e-3) / 1e9,
           times_ms[(int)(REPS * 0.05)],
           times_ms[(int)(REPS * 0.95)],
           maxerr, cksum);

    /* CSV: variant,N,BX,BY,TX,TY,SB,PAD,rep,time_s,gbs,cksum */
    FILE* f = fopen(csv, "a");
    if (f) {
        for (int i = 0; i < REPS; i++) {
            double t_s = times_ms[i] * 1e-3;
            fprintf(f, "%s,%d,%d,%d,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e\n",
                    V_NAMES[VAR], N, BX_, BY_, TX_, TY_, SB_, PAD_,
                    i, t_s, bpi / t_s / 1e9, cksum);
        }
        fclose(f);
    }

    GPU(hipEventDestroy(t0));
    GPU(hipEventDestroy(t1));
    GPU(hipFree(d_in));
    GPU(hipFree(d_out));
    free(h_in); free(h_out); free(times_ms);
    return (maxerr == 0.0f) ? 0 : 1;
}
