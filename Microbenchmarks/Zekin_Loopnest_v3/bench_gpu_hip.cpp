#include "hip/hip_runtime.h"
/*
 * bench_gpu.cu -- GPU-only z_v_grad_w stencil benchmark
 *
 * V1/V2 (je,jk) layout: threadIdx.x -> je (stride-1), threadIdx.y -> jk
 * V3/V4 (jk,je) layout: threadIdx.x -> jk (stride-1), threadIdx.y -> je
 *
 * Now supports a 5th cell distribution "exact" loaded from ICON serialised
 * p_patch files.
 *
 * nlev  = array stride (padded, e.g. 96)
 * nlev_end = actual compute bound (e.g. 90)
 * For synthetic distributions nlev_end == nlev.
 *
 * Environment:
 *   ICON_DATA_PATH  - directory containing p_patch.*.data files
 *                     default: /home/primrose/Work/icon-artifacts/velocity/data_r02b05
 *   Timestep defaults to 9.  Pass a different step as argv[1] (integer).
 *
 * Compile:
 *   nvcc -O3 -arch=sm_80 -std=c++17 -Xcompiler -fopenmp bench_gpu.cu -o bench_gpu
 *
 * Run:
 *   ./bench_gpu                  # synthetic + exact (step 9)
 *   ./bench_gpu 43               # synthetic + exact (step 43)
 *   ICON_DATA_PATH=/my/data ./bench_gpu 1
 */

#include "bench_common.h"
#include "icon_data_loader.h"
#include <ctime>

#if __HIP_PLATFORM_AMD__
#include "hip/hip_runtime.h"
#endif

/* ================================================================ */
/*  CUDA helpers                                                     */
/* ================================================================ */
#define CUDA_CHECK(call) do {                                              \
    hipError_t e = (call);                                                \
    if (e != hipSuccess) {                                                \
        fprintf(stderr, "CUDA error %s:%d: %s\n",                         \
                __FILE__, __LINE__, hipGetErrorString(e));                 \
        exit(1);                                                           \
    }                                                                      \
} while(0)

#define CUDA_LAUNCH_CHECK() do {                                           \
    hipError_t e = hipGetLastError();                                     \
    if (e != hipSuccess) {                                                \
        fprintf(stderr, "CUDA launch error %s:%d: %s\n",                   \
                __FILE__, __LINE__, hipGetErrorString(e));                 \
        exit(1);                                                           \
    }                                                                      \
} while(0)

/* ================================================================ */
/*  CPU reference (serial, for verification)                         */
/*  nlev = stride, nlev_end = loop bound                             */
/* ================================================================ */
template<int V>
static void cpu_reference(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev, int nlev_end)
{
    for (int jk = 0; jk < nlev_end; jk++)
        for (int je = 0; je < N; je++) { STENCIL_BODY(V) }
}

static void cpu_reference_v(int V,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev, int nlev_end)
{
    switch (V) {
        case 1: cpu_reference<1>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
        case 2: cpu_reference<2>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
        case 3: cpu_reference<3>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
        case 4: cpu_reference<4>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
    }
}

/* ================================================================ */
/*  numerical verification                                           */
/* ================================================================ */
static bool verify(const double* got, const double* ref, size_t n,
                   double rtol, double atol,
                   int* n_fail, double* max_rel,
                   size_t* first_fail_idx)
{
    *n_fail  = 0;
    *max_rel = 0.0;
    *first_fail_idx = 0;
    for (size_t i = 0; i < n; i++) {
        double diff = std::abs(got[i] - ref[i]);
        double denom = std::max(std::abs(ref[i]), 1e-300);
        double rel = diff / denom;
        if (rel > *max_rel) *max_rel = rel;
        if (diff > atol + rtol * std::abs(ref[i])) {
            if (*n_fail == 0) *first_fail_idx = i;
            (*n_fail)++;
        }
    }
    return *n_fail == 0;
}

/* ================================================================ */
/*  GPU kernel -- V1/V2: threadIdx.x = je, threadIdx.y = jk          */
/*  nlev = stride, nlev_end = loop bound                             */
/* ================================================================ */

template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_je_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev, int nlev_end)
{
    const int je_base = ((int)blockIdx.x * BX + (int)threadIdx.x) * TX;
    const int jk_base = ((int)blockIdx.y * BY + (int)threadIdx.y) * TY;

    int    ci0_a[TX], ci1_a[TX], vi0_a[TX], vi1_a[TX];
    double id_a[TX],  ip_a[TX],  tg_a[TX];

    #pragma unroll
    for (int tx = 0; tx < TX; tx++) {
        int je = je_base + tx;
        if (je < N) {
            ci0_a[tx] = cell_idx[IN<V>(je, 0, N)];
            ci1_a[tx] = cell_idx[IN<V>(je, 1, N)];
            vi0_a[tx] = vert_idx[IN<V>(je, 0, N)];
            vi1_a[tx] = vert_idx[IN<V>(je, 1, N)];
            id_a[tx]  = inv_dual[je];
            ip_a[tx]  = inv_primal[je];
            tg_a[tx]  = tangent[je];
        }
    }

    #pragma unroll
    for (int ty = 0; ty < TY; ty++) {
        int jk = jk_base + ty;
        if (jk >= nlev_end) continue;
        #pragma unroll
        for (int tx = 0; tx < TX; tx++) {
            int je = je_base + tx;
            if (je >= N) continue;
            int c2d = IC<V>(je, jk, N, nlev);
            out[c2d] =
                vn_ie[c2d] * id_a[tx] *
                    (w[IC<V>(ci0_a[tx], jk, N, nlev)] -
                     w[IC<V>(ci1_a[tx], jk, N, nlev)])
              + z_vt_ie[c2d] * ip_a[tx] * tg_a[tx] *
                    (z_w_v[IC<V>(vi0_a[tx], jk, N, nlev)] -
                     z_w_v[IC<V>(vi1_a[tx], jk, N, nlev)]);
        }
    }
}

/* ================================================================ */
/*  GPU kernel -- V3/V4: threadIdx.x = jk, threadIdx.y = je          */
/*  nlev = stride, nlev_end = loop bound                             */
/* ================================================================ */

template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_jk_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev, int nlev_end)
{
    const int jk_base = ((int)blockIdx.x * BX + (int)threadIdx.x) * TX;
    const int je_base = ((int)blockIdx.y * BY + (int)threadIdx.y) * TY;

    int    ci0_a[TY], ci1_a[TY], vi0_a[TY], vi1_a[TY];
    double id_a[TY],  ip_a[TY],  tg_a[TY];

    #pragma unroll
    for (int ty = 0; ty < TY; ty++) {
        int je = je_base + ty;
        if (je < N) {
            ci0_a[ty] = cell_idx[IN<V>(je, 0, N)];
            ci1_a[ty] = cell_idx[IN<V>(je, 1, N)];
            vi0_a[ty] = vert_idx[IN<V>(je, 0, N)];
            vi1_a[ty] = vert_idx[IN<V>(je, 1, N)];
            id_a[ty]  = inv_dual[je];
            ip_a[ty]  = inv_primal[je];
            tg_a[ty]  = tangent[je];
        }
    }

    #pragma unroll
    for (int ty = 0; ty < TY; ty++) {
        int je = je_base + ty;
        if (je >= N) continue;
        #pragma unroll
        for (int tx = 0; tx < TX; tx++) {
            int jk = jk_base + tx;
            if (jk >= nlev_end) continue;
            int c2d = IC<V>(je, jk, N, nlev);
            out[c2d] =
                vn_ie[c2d] * id_a[ty] *
                    (w[IC<V>(ci0_a[ty], jk, N, nlev)] -
                     w[IC<V>(ci1_a[ty], jk, N, nlev)])
              + z_vt_ie[c2d] * ip_a[ty] * tg_a[ty] *
                    (z_w_v[IC<V>(vi0_a[ty], jk, N, nlev)] -
                     z_w_v[IC<V>(vi1_a[ty], jk, N, nlev)]);
        }
    }
}

/* ================================================================ */
/*  GPU config table  (69 configs)                                   */
/* ================================================================ */
struct GpuCfg { int tx, ty, bx, by; const char* label; };
static constexpr GpuCfg GCFG[] = {
    /* ---- 1-D blocks (BY=1) ---- */
    {1, 1, 256,  1, "1x1_256x1"},
    {1, 1, 128,  1, "1x1_128x1"},
    {1, 1,  64,  1, "1x1_64x1"},
    {2, 1, 256,  1, "2x1_256x1"},
    {2, 1, 128,  1, "2x1_128x1"},
    {4, 1, 256,  1, "4x1_256x1"},
    {4, 1, 128,  1, "4x1_128x1"},
    {4, 1,  64,  1, "4x1_64x1"},
    {8, 1, 128,  1, "8x1_128x1"},
    {8, 1,  64,  1, "8x1_64x1"},
    {1, 2, 256,  1, "1x2_256x1"},
    {1, 2, 128,  1, "1x2_128x1"},
    {1, 4, 256,  1, "1x4_256x1"},
    {1, 4, 128,  1, "1x4_128x1"},
    {1, 8, 128,  1, "1x8_128x1"},
    {2, 2, 256,  1, "2x2_256x1"},
    {2, 2, 128,  1, "2x2_128x1"},
    {2, 4, 128,  1, "2x4_128x1"},
    {4, 2, 128,  1, "4x2_128x1"},
    {4, 2,  64,  1, "4x2_64x1"},
    {4, 4,  64,  1, "4x4_64x1"},
    /* ---- 2-D blocks: 32xBY ---- */
    {1, 1,  32, 32, "1x1_32x32"},
    {1, 1,  32, 16, "1x1_32x16"},
    {1, 1,  32,  8, "1x1_32x8"},
    {1, 1,  32,  4, "1x1_32x4"},
    {2, 1,  32, 16, "2x1_32x16"},
    {2, 1,  32,  8, "2x1_32x8"},
    {2, 1,  32,  4, "2x1_32x4"},
    {4, 1,  32, 16, "4x1_32x16"},
    {4, 1,  32,  8, "4x1_32x8"},
    {4, 1,  32,  4, "4x1_32x4"},
    {1, 2,  32, 16, "1x2_32x16"},
    {1, 2,  32,  8, "1x2_32x8"},
    {1, 2,  32,  4, "1x2_32x4"},
    {1, 4,  32, 16, "1x4_32x16"},
    {1, 4,  32,  8, "1x4_32x8"},
    {1, 4,  32,  4, "1x4_32x4"},
    {2, 2,  32, 16, "2x2_32x16"},
    {2, 2,  32,  8, "2x2_32x8"},
    {2, 2,  32,  4, "2x2_32x4"},
    {2, 4,  32,  8, "2x4_32x8"},
    {2, 4,  32,  4, "2x4_32x4"},
    {4, 2,  32,  8, "4x2_32x8"},
    {4, 2,  32,  4, "4x2_32x4"},
    /* ---- 2-D blocks: 64xBY ---- */
    {1, 1,  64,  8, "1x1_64x8"},
    {1, 1,  64,  4, "1x1_64x4"},
    {1, 1,  64,  2, "1x1_64x2"},
    {2, 1,  64,  8, "2x1_64x8"},
    {2, 1,  64,  4, "2x1_64x4"},
    {2, 1,  64,  2, "2x1_64x2"},
    {4, 1,  64,  4, "4x1_64x4"},
    {4, 1,  64,  2, "4x1_64x2"},
    {1, 2,  64,  8, "1x2_64x8"},
    {1, 2,  64,  4, "1x2_64x4"},
    {1, 4,  64,  4, "1x4_64x4"},
    {2, 2,  64,  4, "2x2_64x4"},
    {2, 2,  64,  2, "2x2_64x2"},
    /* ---- 2-D blocks: 128xBY ---- */
    {1, 1, 128,  4, "1x1_128x4"},
    {1, 1, 128,  2, "1x1_128x2"},
    {2, 1, 128,  4, "2x1_128x4"},
    {2, 1, 128,  2, "2x1_128x2"},
    {1, 2, 128,  4, "1x2_128x4"},
    {1, 2, 128,  2, "1x2_128x2"},
    /* ---- 2-D blocks: 16xBY ---- */
    {1, 1,  16, 16, "1x1_16x16"},
    {2, 1,  16, 16, "2x1_16x16"},
    {1, 2,  16, 16, "1x2_16x16"},
    {2, 2,  16, 16, "2x2_16x16"},
    {4, 1,  16, 16, "4x1_16x16"},
    {4, 2,  16, 16, "4x2_16x16"},
};
static constexpr int N_GCFG = sizeof(GCFG) / sizeof(GCFG[0]);

/* ================================================================ */
/*  GPU launch dispatch                                              */
/*  nlev = stride, nlev_end = compute bound                          */
/* ================================================================ */

template<int V>
static void launch_gpu(int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev, int nlev_end)
{
    #define LG_JE(TX_,TY_,BX_,BY_) do {                                   \
        dim3 blk(BX_, BY_);                                                \
        dim3 grd(((unsigned)N        + (BX_)*(TX_) - 1) / ((BX_)*(TX_)),  \
                 ((unsigned)nlev_end + (BY_)*(TY_) - 1) / ((BY_)*(TY_))); \
        gpu_kernel_je_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>(              \
            out, vn_ie, inv_dual, w, cell_idx,                             \
            z_vt_ie, inv_primal, tangent, z_w_v, vert_idx,                 \
            N, nlev, nlev_end);                                            \
        CUDA_LAUNCH_CHECK();                                               \
    } while(0)

    #define LG_JK(TX_,TY_,BX_,BY_) do {                                       \
        if ((BX_)*(TX_) > nlev_end) break;                                     \
        dim3 blk(BX_, BY_);                                                    \
        dim3 grd(((unsigned)nlev_end + (BX_)*(TX_) - 1) / ((BX_)*(TX_)),      \
                 ((unsigned)N        + (BY_)*(TY_) - 1) / ((BY_)*(TY_)));     \
        if (grd.x > 65535u || grd.y > 65535u) break;                          \
        gpu_kernel_jk_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>(                  \
            out, vn_ie, inv_dual, w, cell_idx,                                 \
            z_vt_ie, inv_primal, tangent, z_w_v, vert_idx,                    \
            N, nlev, nlev_end);                                                \
        CUDA_LAUNCH_CHECK();                                                   \
    } while(0)

    #define LG(TX_,TY_,BX_,BY_) do {                                       \
        if constexpr (V <= 2) { LG_JE(TX_,TY_,BX_,BY_); }                 \
        else                  { LG_JK(TX_,TY_,BX_,BY_); }                  \
    } while(0)

    switch (cfg) {
        /* ---- 1-D blocks (BY=1) ---- */
        case  0: LG(1,1,256,1);  break;
        case  1: LG(1,1,128,1);  break;
        case  2: LG(1,1,64,1);   break;
        case  3: LG(2,1,256,1);  break;
        case  4: LG(2,1,128,1);  break;
        case  5: LG(4,1,256,1);  break;
        case  6: LG(4,1,128,1);  break;
        case  7: LG(4,1,64,1);   break;
        case  8: LG(8,1,128,1);  break;
        case  9: LG(8,1,64,1);   break;
        case 10: LG(1,2,256,1);  break;
        case 11: LG(1,2,128,1);  break;
        case 12: LG(1,4,256,1);  break;
        case 13: LG(1,4,128,1);  break;
        case 14: LG(1,8,128,1);  break;
        case 15: LG(2,2,256,1);  break;
        case 16: LG(2,2,128,1);  break;
        case 17: LG(2,4,128,1);  break;
        case 18: LG(4,2,128,1);  break;
        case 19: LG(4,2,64,1);   break;
        case 20: LG(4,4,64,1);   break;
        /* ---- 2-D blocks: 32xBY ---- */
        case 21: LG(1,1,32,32);  break;
        case 22: LG(1,1,32,16);  break;
        case 23: LG(1,1,32,8);   break;
        case 24: LG(1,1,32,4);   break;
        case 25: LG(2,1,32,16);  break;
        case 26: LG(2,1,32,8);   break;
        case 27: LG(2,1,32,4);   break;
        case 28: LG(4,1,32,16);  break;
        case 29: LG(4,1,32,8);   break;
        case 30: LG(4,1,32,4);   break;
        case 31: LG(1,2,32,16);  break;
        case 32: LG(1,2,32,8);   break;
        case 33: LG(1,2,32,4);   break;
        case 34: LG(1,4,32,16);  break;
        case 35: LG(1,4,32,8);   break;
        case 36: LG(1,4,32,4);   break;
        case 37: LG(2,2,32,16);  break;
        case 38: LG(2,2,32,8);   break;
        case 39: LG(2,2,32,4);   break;
        case 40: LG(2,4,32,8);   break;
        case 41: LG(2,4,32,4);   break;
        case 42: LG(4,2,32,8);   break;
        case 43: LG(4,2,32,4);   break;
        /* ---- 2-D blocks: 64xBY ---- */
        case 44: LG(1,1,64,8);   break;
        case 45: LG(1,1,64,4);   break;
        case 46: LG(1,1,64,2);   break;
        case 47: LG(2,1,64,8);   break;
        case 48: LG(2,1,64,4);   break;
        case 49: LG(2,1,64,2);   break;
        case 50: LG(4,1,64,4);   break;
        case 51: LG(4,1,64,2);   break;
        case 52: LG(1,2,64,8);   break;
        case 53: LG(1,2,64,4);   break;
        case 54: LG(1,4,64,4);   break;
        case 55: LG(2,2,64,4);   break;
        case 56: LG(2,2,64,2);   break;
        /* ---- 2-D blocks: 128xBY ---- */
        case 57: LG(1,1,128,4);  break;
        case 58: LG(1,1,128,2);  break;
        case 59: LG(2,1,128,4);  break;
        case 60: LG(2,1,128,2);  break;
        case 61: LG(1,2,128,4);  break;
        case 62: LG(1,2,128,2);  break;
        /* ---- 2-D blocks: 16xBY ---- */
        case 63: LG(1,1,16,16);  break;
        case 64: LG(2,1,16,16);  break;
        case 65: LG(1,2,16,16);  break;
        case 66: LG(2,2,16,16);  break;
        case 67: LG(4,1,16,16);  break;
        case 68: LG(4,2,16,16);  break;
    }
    #undef LG
    #undef LG_JE
    #undef LG_JK
}

static void launch_gpu_v(int V, int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev, int nlev_end)
{
    switch (V) {
        case 1: launch_gpu<1>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
        case 2: launch_gpu<2>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
        case 3: launch_gpu<3>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
        case 4: launch_gpu<4>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
    }
}

/* ================================================================ */
/*  GPU cache flush -- persistent device buffers                     */
/* ================================================================ */

__global__ void flush_stencil_step(const double* __restrict__ A,
                                   double* __restrict__ B, int N)
{
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= 1 && i < N-1 && j >= 1 && j < N-1)
        B[i * N + j] = 0.25 * (
            A[(i-1)*N + j] + A[(i+1)*N + j] +
            A[i*N + (j-1)] + A[i*N + (j+1)]);
}

static constexpr int FLUSH_N     = 8192*4;
static constexpr int FLUSH_STEPS = 3;

struct GpuFlush {
    double *d_A, *d_B;
    bool inited;

    GpuFlush() : d_A(nullptr), d_B(nullptr), inited(false) {}

    void init() {
        if (inited) return;
        size_t n = (size_t)FLUSH_N * FLUSH_N;
        size_t bytes = n * sizeof(double);

        double* h = new double[n];
        for (size_t i = 0; i < n; i++) {
            uint64_t v = splitmix64(12345ULL + i);
            h[i] = (double)(v >> 11) / (double)(1ULL << 53);
        }

        CUDA_CHECK(hipMalloc(&d_A, bytes));
        CUDA_CHECK(hipMalloc(&d_B, bytes));
        CUDA_CHECK(hipMemcpy(d_A, h, bytes, hipMemcpyHostToDevice));
        CUDA_CHECK(hipMemcpy(d_B, h, bytes, hipMemcpyHostToDevice));
        delete[] h;
        inited = true;
    }

    void flush() {
        init();
        dim3 block(16, 16);
        dim3 grid((FLUSH_N + block.x - 1) / block.x,
                  (FLUSH_N + block.y - 1) / block.y);
        for (int s = 0; s < FLUSH_STEPS; s++) {
            flush_stencil_step<<<grid, block>>>(d_A, d_B, FLUSH_N);
            std::swap(d_A, d_B);
        }
        CUDA_CHECK(hipDeviceSynchronize());

        int ri = rand() % (FLUSH_N * FLUSH_N);
        double val;
        CUDA_CHECK(hipMemcpy(&val, d_A + ri, sizeof(double),
                              hipMemcpyDeviceToHost));
    }

    void destroy() {
        if (d_A) hipFree(d_A);
        if (d_B) hipFree(d_B);
        d_A = d_B = nullptr;
        inited = false;
    }
};

static GpuFlush g_flush;

/* ================================================================ */
/*  Helper: run all GPU configs for one (N, nlev, dist, V) combo     */
/*  nlev = stride (padded), nlev_end = compute bound                 */
/* ================================================================ */
static void run_variant_configs(
    FILE* fcsv,
    int V, int N, int nlev, int nlev_end, const char* dist_label,
    double* h_ref,
    double* h_vn_ie, double* inv_dual,
    double* h_w,     int*    h_cidx,
    double* h_z_vt_ie, double* inv_primal,
    double* tangent_o, double* h_z_w_v,
    int*    h_vidx,
    size_t  sz2d,
    double* d_vn_ie, double* d_inv_dual,
    double* d_w,     int*    d_cidx,
    double* d_z_vt_ie, double* d_inv_primal,
    double* d_tangent, double* d_z_w_v,
    int*    d_vidx, double* d_out,
    hipEvent_t ev0, hipEvent_t ev1,
    double* h_gpu_out)
{
    /* CPU reference -- zero-init so padding matches */
    memset(h_ref, 0, sz2d * sizeof(double));
    cpu_reference_v(V,
        h_ref, h_vn_ie, inv_dual,
        h_w, h_cidx, h_z_vt_ie, inv_primal,
        tangent_o, h_z_w_v, h_vidx, N, nlev, nlev_end);

    /* upload variant data */
    CUDA_CHECK(hipMemcpy(d_vn_ie,      h_vn_ie,    sz2d*sizeof(double), hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_w,          h_w,        sz2d*sizeof(double), hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_z_vt_ie,    h_z_vt_ie,  sz2d*sizeof(double), hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_z_w_v,      h_z_w_v,    sz2d*sizeof(double), hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_cidx,       h_cidx,     N*2*sizeof(int),     hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_vidx,       h_vidx,     N*2*sizeof(int),     hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_inv_dual,   inv_dual,   N*sizeof(double),    hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_inv_primal, inv_primal, N*sizeof(double),    hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_tangent,    tangent_o,  N*sizeof(double),    hipMemcpyHostToDevice));

    for (int ci = 0; ci < N_GCFG; ci++) {
        /* zero d_out so padding matches h_ref */
        CUDA_CHECK(hipMemset(d_out, 0, sz2d * sizeof(double)));

        /* warmup */
        for (int r = 0; r < WARMUP; r++) {
            g_flush.flush();
            launch_gpu_v(V, ci, d_out, d_vn_ie, d_inv_dual,
                d_w, d_cidx, d_z_vt_ie, d_inv_primal,
                d_tangent, d_z_w_v, d_vidx, N, nlev, nlev_end);
            CUDA_CHECK(hipDeviceSynchronize());
        }

        /* verify */
        CUDA_CHECK(hipMemcpy(h_gpu_out, d_out,
            sz2d * sizeof(double), hipMemcpyDeviceToHost));
        int n_fail = 0;
        double max_rel = 0.0;
        size_t first_fail = 0;
        bool ok = verify(h_gpu_out, h_ref, sz2d,
                         1e-8, 1e-12, &n_fail, &max_rel, &first_fail);
        if (!ok) {
            int ff_je, ff_jk;
            if (V <= 2) { ff_je = first_fail % N; ff_jk = first_fail / N; }
            else         { ff_jk = first_fail % nlev; ff_je = first_fail / nlev; }
            printf("VERIFY FAIL: nlev=%d(%d) dist=%-12s V=%d cfg=%-14s  "
                   "fails=%d/%zu max_rel=%.3e\n"
                   "  first_fail: idx=%zu (je=%d,jk=%d) got=%.6e ref=%.6e\n",
                   nlev, nlev_end, dist_label, V, GCFG[ci].label,
                   n_fail, sz2d, max_rel,
                   first_fail, ff_je, ff_jk,
                   h_gpu_out[first_fail], h_ref[first_fail]);
            continue;
        } else if (ci == 0) {
            printf("VERIFY OK:   nlev=%d(%d) dist=%-12s V=%d max_rel=%.3e\n",
                   nlev, nlev_end, dist_label, V, max_rel);
        }

        /* timed runs */
        for (int r = 0; r < NRUNS; r++) {
            g_flush.flush();
            CUDA_CHECK(hipEventRecord(ev0));
            launch_gpu_v(V, ci, d_out, d_vn_ie, d_inv_dual,
                d_w, d_cidx, d_z_vt_ie, d_inv_primal,
                d_tangent, d_z_w_v, d_vidx, N, nlev, nlev_end);
            CUDA_CHECK(hipEventRecord(ev1));
            CUDA_CHECK(hipEventSynchronize(ev1));
            float ms = 0.0f;
            CUDA_CHECK(hipEventElapsedTime(&ms, ev0, ev1));
            fprintf(fcsv,
                "gpu,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%.6f\n",
                V, nlev_end, N, dist_label,
                GCFG[ci].label,
                GCFG[ci].tx, GCFG[ci].ty,
                GCFG[ci].bx, GCFG[ci].by,
                r, (double)ms);
            g_flush.flush();
        }
    }

    printf("Done: nlev=%d(%d)  dist=%-12s  V=%d\n", nlev, nlev_end, dist_label, V);
}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main(int argc, char* argv[]) {
    FILE* fcsv = fopen("z_v_grad_w_gpu.csv", "w");
    if (!fcsv) { perror("fopen"); return 1; }
    fprintf(fcsv,
        "backend,variant,nlev,nproma,cell_dist,"
        "config_label,TX,TY,BX,BY,run_id,time_ms\n");

    const int N = NPROMA;
    std::mt19937 rng(42);

    VertData vd;
    vd.init(N, rng);

    int* cell_logical = new int[N * 2];

    /* ---- load ICON exact data ---- */
    int icon_step = 9;
    if (argc >= 2) icon_step = atoi(argv[1]);

    std::string global_path = icon_global_path(icon_step);
    std::string patch_path  = icon_patch_path(icon_step);

    int icon_nproma = icon_read_nproma(global_path.c_str());
    if (icon_nproma <= 0) {
        fprintf(stderr, "WARNING: could not read nproma from '%s'\n",
                global_path.c_str());
    }
    printf("Loading ICON data from: %s  (nproma=%d)\n",
           patch_path.c_str(), icon_nproma);

    IconEdgeData icon_ed;
    bool have_exact = (icon_nproma > 0) &&
                      icon_load_patch(patch_path.c_str(), icon_nproma, icon_ed);
    if (have_exact) {
        printf("ICON exact data loaded: nproma=%d  "
               "n_edges=%d (valid=%d)  n_cells=%d  n_verts=%d\n",
               icon_ed.nproma, icon_ed.n_edges, icon_ed.n_edges_valid,
               icon_ed.n_cells, icon_ed.n_verts);
    } else {
        fprintf(stderr, "WARNING: failed to load ICON data, "
                "will run synthetic distributions only\n");
    }

    /* print GPU info */
    hipDeviceProp_t prop;
    CUDA_CHECK(hipGetDeviceProperties(&prop, 0));
    printf("GPU: %s  SM count: %d\n", prop.name, prop.multiProcessorCount);
    printf("Max grid dims: (%d, %d, %d)\n",
           prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
    printf("Configs: %d\n", N_GCFG);

    srand((unsigned)time(NULL));

    g_flush.init();
    printf("GPU flush buffers initialized (2 x %.0f MB)\n",
           (double)FLUSH_N * FLUSH_N * 8 / 1e6);

    for (int nlev_i = 0; nlev_i < N_NLEVS; nlev_i++) {
        int nlev = NLEVS[nlev_i];

        /* ============================================================ */
        /*  Synthetic distributions -- nlev_end == nlev (no padding)    */
        /* ============================================================ */
        {
            BenchData bd;
            bd.alloc(N, nlev);
            bd.fill(nlev);

            double* h_ref     = new double[bd.sz2d];
            double* h_gpu_out = new double[bd.sz2d];

            double *d_vn_ie, *d_w, *d_z_vt_ie, *d_z_w_v, *d_out;
            double *d_inv_dual, *d_inv_primal, *d_tangent;
            int    *d_cidx, *d_vidx;

            CUDA_CHECK(hipMalloc(&d_vn_ie,      bd.sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_w,          bd.sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_z_vt_ie,    bd.sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_z_w_v,      bd.sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_out,        bd.sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_inv_dual,   N * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_inv_primal, N * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_tangent,    N * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_cidx,       N * 2 * sizeof(int)));
            CUDA_CHECK(hipMalloc(&d_vidx,       N * 2 * sizeof(int)));

            hipEvent_t ev0, ev1;
            CUDA_CHECK(hipEventCreate(&ev0));
            CUDA_CHECK(hipEventCreate(&ev1));

            for (int di = 0; di < 4; di++) {
                CellDist dist = (CellDist)di;
                gen_cell_idx_logical(cell_logical, N, dist, rng);

                for (int V = 1; V <= 4; V++) {
                    bd.set_variant(V, cell_logical, vd.logical);

                    run_variant_configs(fcsv, V, N, nlev, nlev, dist_name[di],
                        h_ref,
                        bd.h_vn_ie, bd.inv_dual,
                        bd.h_w, bd.h_cidx,
                        bd.h_z_vt_ie, bd.inv_primal,
                        bd.tangent_o, bd.h_z_w_v,
                        bd.h_vidx, bd.sz2d,
                        d_vn_ie, d_inv_dual,
                        d_w, d_cidx,
                        d_z_vt_ie, d_inv_primal,
                        d_tangent, d_z_w_v,
                        d_vidx, d_out,
                        ev0, ev1, h_gpu_out);

                    fflush(fcsv);
                }
            }

            hipFree(d_vn_ie);      hipFree(d_w);
            hipFree(d_z_vt_ie);    hipFree(d_z_w_v);    hipFree(d_out);
            hipFree(d_inv_dual);   hipFree(d_inv_primal);
            hipFree(d_tangent);    hipFree(d_cidx);      hipFree(d_vidx);
            CUDA_CHECK(hipEventDestroy(ev0));
            CUDA_CHECK(hipEventDestroy(ev1));
            delete[] h_ref;
            delete[] h_gpu_out;
            bd.free_all();
        }

        /* ============================================================ */
        /*  EXACT distribution from ICON grid data                      */
        /*  nlev_end = original nlev, nlev = padded to 32               */
        /* ============================================================ */
        if (have_exact) {
            const int Ne = icon_ed.n_edges;
            const int nlev_end = nlev;
            const int nlev_padded = icon_pad_nlev(nlev);

            printf("\n=== EXACT distribution: nlev_end=%d nlev_padded=%d  "
                   "Ne=%d (nproma=%d * nblks_e=%d) ===\n",
                   nlev_end, nlev_padded, Ne, icon_ed.nproma, icon_ed.nblks_e);

            if (icon_ed.n_cells > Ne || icon_ed.n_verts > Ne) {
                fprintf(stderr, "WARNING: n_cells=%d or n_verts=%d > Ne=%d; "
                        "skipping exact for nlev=%d\n",
                        icon_ed.n_cells, icon_ed.n_verts, Ne, nlev_end);
                continue;
            }

            BenchData bd;
            bd.alloc(Ne, nlev_padded);
            bd.fill(nlev_padded);

            size_t sz2d = bd.sz2d;
            double* h_ref     = new double[sz2d];
            double* h_gpu_out = new double[sz2d];

            double *d_vn_ie, *d_w, *d_z_vt_ie, *d_z_w_v, *d_out;
            double *d_inv_dual, *d_inv_primal, *d_tangent;
            int    *d_cidx, *d_vidx;

            CUDA_CHECK(hipMalloc(&d_vn_ie,      sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_w,          sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_z_vt_ie,    sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_z_w_v,      sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_out,        sz2d * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_inv_dual,   Ne * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_inv_primal, Ne * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_tangent,    Ne * sizeof(double)));
            CUDA_CHECK(hipMalloc(&d_cidx,       Ne * 2 * sizeof(int)));
            CUDA_CHECK(hipMalloc(&d_vidx,       Ne * 2 * sizeof(int)));

            hipEvent_t ev0, ev1;
            CUDA_CHECK(hipEventCreate(&ev0));
            CUDA_CHECK(hipEventCreate(&ev1));

            /* Build logical connectivity arrays from icon_ed */
            int* exact_cell_logical = new int[Ne * 2];
            int* exact_vert_logical = new int[Ne * 2];
            for (int i = 0; i < Ne * 2; i++) {
                exact_cell_logical[i] = icon_ed.cell_idx[i];
                exact_vert_logical[i] = icon_ed.vert_idx[i];
            }

            /* Override 1-D geometry with real ICON values */
            for (int je = 0; je < Ne; je++) {
                bd.inv_dual[je]   = icon_ed.inv_dual[je];
                bd.inv_primal[je] = icon_ed.inv_primal[je];
                bd.tangent_o[je]  = icon_ed.tangent_o[je];
            }

            for (int V = 1; V <= 4; V++) {
                bd.set_variant(V, exact_cell_logical, exact_vert_logical);

                run_variant_configs(fcsv, V, Ne, nlev_padded, nlev_end, "exact",
                    h_ref,
                    bd.h_vn_ie, bd.inv_dual,
                    bd.h_w, bd.h_cidx,
                    bd.h_z_vt_ie, bd.inv_primal,
                    bd.tangent_o, bd.h_z_w_v,
                    bd.h_vidx, sz2d,
                    d_vn_ie, d_inv_dual,
                    d_w, d_cidx,
                    d_z_vt_ie, d_inv_primal,
                    d_tangent, d_z_w_v,
                    d_vidx, d_out,
                    ev0, ev1, h_gpu_out);

                fflush(fcsv);
            }

            delete[] exact_cell_logical;
            delete[] exact_vert_logical;

            hipFree(d_vn_ie);      hipFree(d_w);
            hipFree(d_z_vt_ie);    hipFree(d_z_w_v);    hipFree(d_out);
            hipFree(d_inv_dual);   hipFree(d_inv_primal);
            hipFree(d_tangent);    hipFree(d_cidx);      hipFree(d_vidx);
            CUDA_CHECK(hipEventDestroy(ev0));
            CUDA_CHECK(hipEventDestroy(ev1));
            delete[] h_ref;
            delete[] h_gpu_out;
            bd.free_all();
        }
    }

    g_flush.destroy();
    vd.free_all();
    delete[] cell_logical;
    if (have_exact) icon_ed.free_all();
    fclose(fcsv);

    printf("\nResults written to z_v_grad_w_gpu.csv\n");
    return 0;
}
