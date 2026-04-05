/*
 * bench_gpu.cu -- GPU-only z_v_grad_w stencil benchmark
 *
 * V1/V2 (je,jk) layout: threadIdx.x -> je (stride-1), threadIdx.y -> jk
 * V3/V4 (jk,je) layout: threadIdx.x -> jk (stride-1), threadIdx.y -> je
 * V5    same as V4 but nlev padded to next multiple of 32
 *
 * nlev     = array stride (== nlev_end for V1-V4, padded for V5)
 * nlev_end = actual compute bound
 *
 * Environment:
 *   ICON_DATA_PATH  - directory containing p_patch.*.data files
 *   Timestep defaults to 9.  Pass a different step as argv[1].
 *
 * Compile:
 *   nvcc -O3 -arch=sm_80 -std=c++17 -Xcompiler -fopenmp bench_gpu.cu -o bench_gpu
 */

#include "bench_common.h"
#include "icon_data_loader.h"
#include <ctime>
#include <cassert>

#if __HIP_PLATFORM_AMD__
#include "hip/hip_runtime.h"
#endif

/* ================================================================ */
/*  CUDA helpers                                                     */
/* ================================================================ */
#define CUDA_CHECK(call) do {                                              \
    cudaError_t e = (call);                                                \
    if (e != cudaSuccess) {                                                \
        fprintf(stderr, "CUDA error %s:%d: %s\n",                         \
                __FILE__, __LINE__, cudaGetErrorString(e));                 \
        exit(1);                                                           \
    }                                                                      \
} while(0)

#define CUDA_LAUNCH_CHECK() do {                                           \
    cudaError_t e = cudaGetLastError();                                     \
    if (e != cudaSuccess) {                                                \
        fprintf(stderr, "CUDA launch error %s:%d: %s\n",                   \
                __FILE__, __LINE__, cudaGetErrorString(e));                 \
        exit(1);                                                           \
    }                                                                      \
} while(0)

/* ================================================================ */
/*  CPU reference  (nlev = stride, nlev_end = loop bound)            */
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
    int kV = kern_v(V);
    switch (kV) {
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
                   int* n_fail, double* max_rel, size_t* first_fail_idx)
{
    *n_fail = 0; *max_rel = 0.0; *first_fail_idx = 0;
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
/*  GPU kernels  (nlev = stride, nlev_end = loop bound)              */
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
    int ci0_a[TX],ci1_a[TX],vi0_a[TX],vi1_a[TX];
    double id_a[TX],ip_a[TX],tg_a[TX];
    #pragma unroll
    for (int tx=0;tx<TX;tx++){int je=je_base+tx;if(je<N){
        ci0_a[tx]=cell_idx[IN<V>(je,0,N)]; ci1_a[tx]=cell_idx[IN<V>(je,1,N)];
        vi0_a[tx]=vert_idx[IN<V>(je,0,N)]; vi1_a[tx]=vert_idx[IN<V>(je,1,N)];
        id_a[tx]=inv_dual[je]; ip_a[tx]=inv_primal[je]; tg_a[tx]=tangent[je];
    }}
    #pragma unroll
    for (int ty=0;ty<TY;ty++){int jk=jk_base+ty;if(jk>=nlev_end)continue;
    #pragma unroll
    for (int tx=0;tx<TX;tx++){int je=je_base+tx;if(je>=N)continue;
        int c2d=IC<V>(je,jk,N,nlev);
        out[c2d]=vn_ie[c2d]*id_a[tx]*
            (w[IC<V>(ci0_a[tx],jk,N,nlev)]-w[IC<V>(ci1_a[tx],jk,N,nlev)])
          +z_vt_ie[c2d]*ip_a[tx]*tg_a[tx]*
            (z_w_v[IC<V>(vi0_a[tx],jk,N,nlev)]-z_w_v[IC<V>(vi1_a[tx],jk,N,nlev)]);
    }}
}

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
    int ci0_a[TY],ci1_a[TY],vi0_a[TY],vi1_a[TY];
    double id_a[TY],ip_a[TY],tg_a[TY];
    #pragma unroll
    for (int ty=0;ty<TY;ty++){int je=je_base+ty;if(je<N){
        ci0_a[ty]=cell_idx[IN<V>(je,0,N)]; ci1_a[ty]=cell_idx[IN<V>(je,1,N)];
        vi0_a[ty]=vert_idx[IN<V>(je,0,N)]; vi1_a[ty]=vert_idx[IN<V>(je,1,N)];
        id_a[ty]=inv_dual[je]; ip_a[ty]=inv_primal[je]; tg_a[ty]=tangent[je];
    }}
    #pragma unroll
    for (int ty=0;ty<TY;ty++){int je=je_base+ty;if(je>=N)continue;
    #pragma unroll
    for (int tx=0;tx<TX;tx++){int jk=jk_base+tx;if(jk>=nlev_end)continue;
        int c2d=IC<V>(je,jk,N,nlev);
        out[c2d]=vn_ie[c2d]*id_a[ty]*
            (w[IC<V>(ci0_a[ty],jk,N,nlev)]-w[IC<V>(ci1_a[ty],jk,N,nlev)])
          +z_vt_ie[c2d]*ip_a[ty]*tg_a[ty]*
            (z_w_v[IC<V>(vi0_a[ty],jk,N,nlev)]-z_w_v[IC<V>(vi1_a[ty],jk,N,nlev)]);
    }}
}

/* ================================================================ */
/*  GPU config table  (69 configs)                                   */
/* ================================================================ */
struct GpuCfg { int tx, ty, bx, by; const char* label; };
/* ================================================================ */
/*  GPU config table                                                 */
/*  Convention: TX tiles the stride-1 dim, TY tiles the other dim.   */
/*  je_first (V<=2): TX->je, TY->jk                                 */
/*  jk_first (V>=3): TX->jk, TY->je                                 */
/*                                                                   */
/*  Groups:                                                          */
/*    A) TX >= TY  (original bias)                                   */
/*    B) TY >  TX  (new: favours NPROMA-heavy tiling for V3-V5)     */
/* ================================================================ */
static constexpr GpuCfg GCFG[] = {
    /* ---- A: TX >= TY, 1D blocks (BY=1) ---- */
    /* 0  */ {1,1,256,1,"1x1_256x1"},
    /* 1  */ {1,1,128,1,"1x1_128x1"},
    /* 2  */ {1,1,64,1,"1x1_64x1"},
    /* 3  */ {2,1,256,1,"2x1_256x1"},
    /* 4  */ {2,1,128,1,"2x1_128x1"},
    /* 5  */ {4,1,256,1,"4x1_256x1"},
    /* 6  */ {4,1,128,1,"4x1_128x1"},
    /* 7  */ {4,1,64,1,"4x1_64x1"},
    /* 8  */ {8,1,128,1,"8x1_128x1"},
    /* 9  */ {8,1,64,1,"8x1_64x1"},
    /* ---- A: TX >= TY, 1D blocks (BY=1) with small TY ---- */
    /* 10 */ {1,2,256,1,"1x2_256x1"},
    /* 11 */ {1,2,128,1,"1x2_128x1"},
    /* 12 */ {1,4,256,1,"1x4_256x1"},
    /* 13 */ {1,4,128,1,"1x4_128x1"},
    /* 14 */ {1,8,128,1,"1x8_128x1"},
    /* 15 */ {2,2,256,1,"2x2_256x1"},
    /* 16 */ {2,2,128,1,"2x2_128x1"},
    /* 17 */ {2,4,128,1,"2x4_128x1"},
    /* 18 */ {4,2,128,1,"4x2_128x1"},
    /* 19 */ {4,2,64,1,"4x2_64x1"},
    /* 20 */ {4,4,64,1,"4x4_64x1"},
    /* ---- A: 2D blocks, square-ish ---- */
    /* 21 */ {1,1,32,32,"1x1_32x32"},
    /* 22 */ {1,1,32,16,"1x1_32x16"},
    /* 23 */ {1,1,32,8,"1x1_32x8"},
    /* 24 */ {1,1,32,4,"1x1_32x4"},
    /* 25 */ {2,1,32,16,"2x1_32x16"},
    /* 26 */ {2,1,32,8,"2x1_32x8"},
    /* 27 */ {2,1,32,4,"2x1_32x4"},
    /* 28 */ {4,1,32,16,"4x1_32x16"},
    /* 29 */ {4,1,32,8,"4x1_32x8"},
    /* 30 */ {4,1,32,4,"4x1_32x4"},
    /* 31 */ {1,2,32,16,"1x2_32x16"},
    /* 32 */ {1,2,32,8,"1x2_32x8"},
    /* 33 */ {1,2,32,4,"1x2_32x4"},
    /* 34 */ {1,4,32,16,"1x4_32x16"},
    /* 35 */ {1,4,32,8,"1x4_32x8"},
    /* 36 */ {1,4,32,4,"1x4_32x4"},
    /* 37 */ {2,2,32,16,"2x2_32x16"},
    /* 38 */ {2,2,32,8,"2x2_32x8"},
    /* 39 */ {2,2,32,4,"2x2_32x4"},
    /* 40 */ {2,4,32,8,"2x4_32x8"},
    /* 41 */ {2,4,32,4,"2x4_32x4"},
    /* 42 */ {4,2,32,8,"4x2_32x8"},
    /* 43 */ {4,2,32,4,"4x2_32x4"},
    /* ---- A: 2D blocks, moderate ---- */
    /* 44 */ {1,1,64,8,"1x1_64x8"},
    /* 45 */ {1,1,64,4,"1x1_64x4"},
    /* 46 */ {1,1,64,2,"1x1_64x2"},
    /* 47 */ {2,1,64,8,"2x1_64x8"},
    /* 48 */ {2,1,64,4,"2x1_64x4"},
    /* 49 */ {2,1,64,2,"2x1_64x2"},
    /* 50 */ {4,1,64,4,"4x1_64x4"},
    /* 51 */ {4,1,64,2,"4x1_64x2"},
    /* 52 */ {1,2,64,8,"1x2_64x8"},
    /* 53 */ {1,2,64,4,"1x2_64x4"},
    /* 54 */ {1,4,64,4,"1x4_64x4"},
    /* 55 */ {2,2,64,4,"2x2_64x4"},
    /* 56 */ {2,2,64,2,"2x2_64x2"},
    /* 57 */ {1,1,128,4,"1x1_128x4"},
    /* 58 */ {1,1,128,2,"1x1_128x2"},
    /* 59 */ {2,1,128,4,"2x1_128x4"},
    /* 60 */ {2,1,128,2,"2x1_128x2"},
    /* 61 */ {1,2,128,4,"1x2_128x4"},
    /* 62 */ {1,2,128,2,"1x2_128x2"},
    /* ---- A: 2D blocks, 16x16 ---- */
    /* 63 */ {1,1,16,16,"1x1_16x16"},
    /* 64 */ {2,1,16,16,"2x1_16x16"},
    /* 65 */ {1,2,16,16,"1x2_16x16"},
    /* 66 */ {2,2,16,16,"2x2_16x16"},
    /* 67 */ {4,1,16,16,"4x1_16x16"},
    /* 68 */ {4,2,16,16,"4x2_16x16"},
    /* ---- A: 1D blocks, 96-wide ---- */
    /* 69 */ {1,1,96,1,"1x1_96x1"},
    /* 70 */ {1,2,96,1,"1x2_96x1"},
    /* 71 */ {1,4,96,1,"1x4_96x1"},
    /* 72 */ {1,8,96,1,"1x8_96x1"},
    /* 73 */ {1,1,96,2,"1x1_96x2"},
    /* 74 */ {1,2,96,2,"1x2_96x2"},
    /* 75 */ {1,4,96,2,"1x4_96x2"},
    /* 76 */ {1,8,96,2,"1x8_96x2"},
    /* 77 */ {1,1,96,4,"1x1_96x4"},
    /* 78 */ {1,2,96,4,"1x2_96x4"},
    /* 79 */ {1,4,96,4,"1x4_96x4"},
    /* 80 */ {1,8,96,4,"1x8_96x4"},
 
    /* ============================================================ */
    /*  B: TY > TX  (NPROMA-heavy tiling for jk_first kernels)      */
    /* ============================================================ */
 
    /* ---- B: 1D blocks (BY=1), big TY ---- */
    /* 81  */ {1,2,64,1,"1x2_64x1"},
    /* 82  */ {1,4,64,1,"1x4_64x1"},
    /* 83  */ {1,8,64,1,"1x8_64x1"},
    /* 84  */ {1,2,32,1,"1x2_32x1"},
    /* 85  */ {1,4,32,1,"1x4_32x1"},
    /* 86  */ {1,8,32,1,"1x8_32x1"},
    /* 87  */ {1,16,32,1,"1x16_32x1"},
    /* 88  */ {1,2,16,1,"1x2_16x1"},
    /* 89  */ {1,4,16,1,"1x4_16x1"},
    /* 90  */ {1,8,16,1,"1x8_16x1"},
    /* 91  */ {1,16,16,1,"1x16_16x1"},
 
    /* ---- B: 2D blocks, TY > TX, BX=32 ---- */
    /* 92  */ {1,2,32,32,"1x2_32x32"},
    /* 93  */ {1,4,32,16,"1x4_32x16b"},
    /* 94  */ {1,4,32,8,"1x4_32x8b"},
    /* 95  */ {1,8,32,8,"1x8_32x8"},
    /* 96  */ {1,8,32,4,"1x8_32x4"},
    /* 97  */ {1,8,32,16,"1x8_32x16"},
    /* 98  */ {1,16,32,4,"1x16_32x4"},
    /* 99  */ {1,16,32,8,"1x16_32x8"},
 
    /* ---- B: 2D blocks, TY > TX, BX=16 ---- */
    /* 100 */ {1,2,16,16,"1x2_16x16b"},
    /* 101 */ {1,4,16,16,"1x4_16x16"},
    /* 102 */ {1,8,16,16,"1x8_16x16"},
    /* 103 */ {1,4,16,8,"1x4_16x8"},
    /* 104 */ {1,8,16,8,"1x8_16x8"},
    /* 105 */ {1,16,16,4,"1x16_16x4"},
    /* 106 */ {1,16,16,8,"1x16_16x8"},
 
    /* ---- B: 2D blocks, TY > TX, BX=64 ---- */
    /* 107 */ {1,2,64,4,"1x2_64x4b"},
    /* 108 */ {1,4,64,2,"1x4_64x2"},
    /* 109 */ {1,8,64,2,"1x8_64x2"},
    /* 110 */ {1,8,64,4,"1x8_64x4"},
    /* 111 */ {1,4,64,8,"1x4_64x8b"},
 
    /* ---- B: mixed TX=2 with TY=4,8 ---- */
    /* 112 */ {2,4,32,16,"2x4_32x16"},
    /* 113 */ {2,8,32,8,"2x8_32x8"},
    /* 114 */ {2,8,32,4,"2x8_32x4"},
    /* 115 */ {2,4,64,4,"2x4_64x4"},
    /* 116 */ {2,8,64,2,"2x8_64x2"},
    /* 117 */ {2,4,16,16,"2x4_16x16"},
    /* 118 */ {2,8,16,8,"2x8_16x8"},
    /* 119 */ {2,8,16,16,"2x8_16x16"},
 
    /* ---- B: 1D blocks, 128/256-wide, bigger TY ---- */
    /* 120 */ {1,4,128,2,"1x4_128x2"},
    /* 121 */ {1,8,128,2,"1x8_128x2"},
    /* 122 */ {1,4,256,1,"1x4_256x1b"},
    /* 123 */ {1,8,256,1,"1x8_256x1"},
    /* 124 */ {1,16,128,1,"1x16_128x1"},
    /* 125 */ {1,16,64,1,"1x16_64x1"},
 
    /* ---- B: 96-wide, bigger TY ---- */
    /* 126 */ {1,16,96,1,"1x16_96x1"},
    /* 127 */ {1,16,96,2,"1x16_96x2"},
    /* 128 */ {1,8,96,8,"1x8_96x8"},
};
static constexpr int N_GCFG = sizeof(GCFG)/sizeof(GCFG[0]);

/* ================================================================ */
/*  GPU launch dispatch                                              */
/* ================================================================ */
template<int V>
static bool launch_gpu(int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev, int nlev_end)
{
    bool launched = false;

    #define LG_JE(TX_,TY_,BX_,BY_) do {                                   \
        dim3 blk(BX_,BY_);                                                 \
        dim3 grd(((unsigned)N       +(BX_)*(TX_)-1)/((BX_)*(TX_)),        \
                 ((unsigned)nlev_end+(BY_)*(TY_)-1)/((BY_)*(TY_)));       \
        if (grd.x>65535u||grd.y>65535u) break;                             \
        gpu_kernel_je_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>(              \
            out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,      \
            z_w_v,vert_idx,N,nlev,nlev_end);                               \
        CUDA_LAUNCH_CHECK(); launched=true;                                \
    } while(0)

    #define LG_JK(TX_,TY_,BX_,BY_) do {                                   \
        if ((BX_)*(TX_) > (unsigned)nlev_end) break;                       \
        dim3 blk(BX_,BY_);                                                 \
        dim3 grd(((unsigned)nlev_end+(BX_)*(TX_)-1)/((BX_)*(TX_)),        \
                 ((unsigned)N       +(BY_)*(TY_)-1)/((BY_)*(TY_)));       \
        if (grd.x>65535u||grd.y>65535u) break;                             \
        gpu_kernel_jk_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>(              \
            out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,      \
            z_w_v,vert_idx,N,nlev,nlev_end);                               \
        CUDA_LAUNCH_CHECK(); launched=true;                                \
    } while(0)

    #define LG(TX_,TY_,BX_,BY_) do {                                      \
        if constexpr (V<=2) { LG_JE(TX_,TY_,BX_,BY_); }                   \
        else                { LG_JK(TX_,TY_,BX_,BY_); }                   \
    } while(0)


    switch(cfg){
    /* ---- A: original configs ---- */
    case  0:LG(1,1,256,1);break; case  1:LG(1,1,128,1);break; case  2:LG(1,1,64,1);break;
    case  3:LG(2,1,256,1);break; case  4:LG(2,1,128,1);break;
    case  5:LG(4,1,256,1);break; case  6:LG(4,1,128,1);break; case  7:LG(4,1,64,1);break;
    case  8:LG(8,1,128,1);break; case  9:LG(8,1,64,1);break;
    case 10:LG(1,2,256,1);break; case 11:LG(1,2,128,1);break;
    case 12:LG(1,4,256,1);break; case 13:LG(1,4,128,1);break; case 14:LG(1,8,128,1);break;
    case 15:LG(2,2,256,1);break; case 16:LG(2,2,128,1);break; case 17:LG(2,4,128,1);break;
    case 18:LG(4,2,128,1);break; case 19:LG(4,2,64,1);break;  case 20:LG(4,4,64,1);break;
    case 21:LG(1,1,32,32);break; case 22:LG(1,1,32,16);break;
    case 23:LG(1,1,32,8);break;  case 24:LG(1,1,32,4);break;
    case 25:LG(2,1,32,16);break; case 26:LG(2,1,32,8);break;  case 27:LG(2,1,32,4);break;
    case 28:LG(4,1,32,16);break; case 29:LG(4,1,32,8);break;  case 30:LG(4,1,32,4);break;
    case 31:LG(1,2,32,16);break; case 32:LG(1,2,32,8);break;  case 33:LG(1,2,32,4);break;
    case 34:LG(1,4,32,16);break; case 35:LG(1,4,32,8);break;  case 36:LG(1,4,32,4);break;
    case 37:LG(2,2,32,16);break; case 38:LG(2,2,32,8);break;  case 39:LG(2,2,32,4);break;
    case 40:LG(2,4,32,8);break;  case 41:LG(2,4,32,4);break;
    case 42:LG(4,2,32,8);break;  case 43:LG(4,2,32,4);break;
    case 44:LG(1,1,64,8);break;  case 45:LG(1,1,64,4);break;  case 46:LG(1,1,64,2);break;
    case 47:LG(2,1,64,8);break;  case 48:LG(2,1,64,4);break;  case 49:LG(2,1,64,2);break;
    case 50:LG(4,1,64,4);break;  case 51:LG(4,1,64,2);break;
    case 52:LG(1,2,64,8);break;  case 53:LG(1,2,64,4);break;  case 54:LG(1,4,64,4);break;
    case 55:LG(2,2,64,4);break;  case 56:LG(2,2,64,2);break;
    case 57:LG(1,1,128,4);break; case 58:LG(1,1,128,2);break;
    case 59:LG(2,1,128,4);break; case 60:LG(2,1,128,2);break;
    case 61:LG(1,2,128,4);break; case 62:LG(1,2,128,2);break;
    case 63:LG(1,1,16,16);break; case 64:LG(2,1,16,16);break; case 65:LG(1,2,16,16);break;
    case 66:LG(2,2,16,16);break; case 67:LG(4,1,16,16);break; case 68:LG(4,2,16,16);break;
    case 69:LG(1,1,96,1);break;  case 70:LG(1,2,96,1);break;
    case 71:LG(1,4,96,1);break;  case 72:LG(1,8,96,1);break;
    case 73:LG(1,1,96,2);break;  case 74:LG(1,2,96,2);break;
    case 75:LG(1,4,96,2);break;  case 76:LG(1,8,96,2);break;
    case 77:LG(1,1,96,4);break;  case 78:LG(1,2,96,4);break;
    case 79:LG(1,4,96,4);break;  case 80:LG(1,8,96,4);break;
 
    /* ---- B: TY > TX configs ---- */
    case  81:LG(1,2,64,1);break;   case  82:LG(1,4,64,1);break;
    case  83:LG(1,8,64,1);break;
    case  84:LG(1,2,32,2);break;   case  85:LG(1,4,32,2);break;
    case  86:LG(1,8,32,2);break;   case  87:LG(1,16,32,2);break;
    case  88:LG(1,2,16,4);break;   case  89:LG(1,4,16,4);break;
    case  90:LG(1,8,16,4);break;   case  91:LG(1,16,16,4);break;
    case  92:LG(1,2,32,32);break;  case  93:LG(1,4,32,16);break;
    case  94:LG(1,4,32,8);break;   case  95:LG(1,8,32,8);break;
    case  96:LG(1,8,32,4);break;   case  97:LG(1,8,32,16);break;
    case  98:LG(1,16,32,4);break;  case  99:LG(1,16,32,8);break;
    case 100:LG(1,2,16,16);break;  case 101:LG(1,4,16,16);break;
    case 102:LG(1,8,16,16);break;  case 103:LG(1,4,16,8);break;
    case 104:LG(1,8,16,8);break;   case 105:LG(1,16,16,4);break;
    case 106:LG(1,16,16,8);break;
    case 107:LG(1,2,64,4);break;   case 108:LG(1,4,64,2);break;
    case 109:LG(1,8,64,2);break;   case 110:LG(1,8,64,4);break;
    case 111:LG(1,4,64,8);break;
    case 112:LG(2,4,32,16);break;  case 113:LG(2,8,32,8);break;
    case 114:LG(2,8,32,4);break;   case 115:LG(2,4,64,4);break;
    case 116:LG(2,8,64,2);break;   case 117:LG(2,4,16,16);break;
    case 118:LG(2,8,16,8);break;   case 119:LG(2,8,16,16);break;
    case 120:LG(1,4,128,2);break;  case 121:LG(1,8,128,2);break;
    case 122:LG(1,4,256,1);break;  case 123:LG(1,8,256,1);break;
    case 124:LG(1,16,128,1);break; case 125:LG(1,16,64,1);break;
    case 126:LG(1,16,96,1);break;  case 127:LG(1,16,96,2);break;
    case 128:LG(1,8,96,8);break;
    }
    #undef LG
    #undef LG_JE
    #undef LG_JK
    return launched;
}

static bool launch_gpu_v(int V, int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev, int nlev_end)
{
    int kV = kern_v(V);
    switch(kV){
    case 1: return launch_gpu<1>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
    case 2: return launch_gpu<2>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
    case 3: return launch_gpu<3>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
    case 4: return launch_gpu<4>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev,nlev_end); break;
    }
    return false;
}

/* ================================================================ */
/*  GPU cache flush                                                  */
/* ================================================================ */
__global__ void flush_stencil_step(const double* __restrict__ A,
                                   double* __restrict__ B, int N) {
    int i=blockIdx.y*blockDim.y+threadIdx.y;
    int j=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=1&&i<N-1&&j>=1&&j<N-1)
        B[i*N+j]=0.25*(A[(i-1)*N+j]+A[(i+1)*N+j]+A[i*N+(j-1)]+A[i*N+(j+1)]);
}
static constexpr int FLUSH_N=8192*4, FLUSH_STEPS=3;
struct GpuFlush {
    double *d_A=nullptr, *d_B=nullptr; bool inited=false;
    void init(){if(inited)return; size_t n=(size_t)FLUSH_N*FLUSH_N;
        double*h=new double[n]; for(size_t i=0;i<n;i++){uint64_t v=splitmix64(12345ULL+i);h[i]=(double)(v>>11)/(double)(1ULL<<53);}
        CUDA_CHECK(cudaMalloc(&d_A,n*8)); CUDA_CHECK(cudaMalloc(&d_B,n*8));
        CUDA_CHECK(cudaMemcpy(d_A,h,n*8,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_B,h,n*8,cudaMemcpyHostToDevice));
        delete[]h; inited=true;}
    void flush(){init(); dim3 bl(16,16),gr((FLUSH_N+15)/16,(FLUSH_N+15)/16);
        for(int s=0;s<FLUSH_STEPS;s++){flush_stencil_step<<<gr,bl>>>(d_A,d_B,FLUSH_N);std::swap(d_A,d_B);}
        CUDA_CHECK(cudaDeviceSynchronize());
        int ri=rand()%(FLUSH_N*FLUSH_N); double val;
        CUDA_CHECK(cudaMemcpy(&val,d_A+ri,8,cudaMemcpyDeviceToHost));}
    void destroy(){if(d_A)cudaFree(d_A);if(d_B)cudaFree(d_B);d_A=d_B=nullptr;inited=false;}
};
static GpuFlush g_flush;

/* ================================================================ */
/*  run_variant_configs                                              */
/*  V = label (1-5), nlev = stride, nlev_end = compute bound         */
/* ================================================================ */
static void run_variant_configs(
    FILE* fcsv,
    int V, int N, int nlev, int nlev_end, const char* dist_label,
    double* h_ref,
    double* h_vn_ie, double* inv_dual,
    double* h_w, int* h_cidx,
    double* h_z_vt_ie, double* inv_primal,
    double* tangent_o, double* h_z_w_v,
    int* h_vidx, size_t sz2d,
    double* d_vn_ie, double* d_inv_dual,
    double* d_w, int* d_cidx,
    double* d_z_vt_ie, double* d_inv_primal,
    double* d_tangent, double* d_z_w_v,
    int* d_vidx, double* d_out,
    cudaEvent_t ev0, cudaEvent_t ev1,
    double* h_gpu_out)
{
    int kV = kern_v(V);

    /* CPU reference -- zero-init so padding matches */
    memset(h_ref, 0, sz2d * sizeof(double));
    cpu_reference_v(V,
        h_ref, h_vn_ie, inv_dual,
        h_w, h_cidx, h_z_vt_ie, inv_primal,
        tangent_o, h_z_w_v, h_vidx, N, nlev, nlev_end);

    /* upload */
    CUDA_CHECK(cudaMemcpy(d_vn_ie,      h_vn_ie,    sz2d*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_w,          h_w,        sz2d*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_vt_ie,    h_z_vt_ie,  sz2d*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_w_v,      h_z_w_v,    sz2d*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cidx,       h_cidx,     N*2*4,  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vidx,       h_vidx,     N*2*4,  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_dual,   inv_dual,   N*8,    cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_primal, inv_primal, N*8,    cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_tangent,    tangent_o,  N*8,    cudaMemcpyHostToDevice));

    for (int ci = 0; ci < N_GCFG; ci++) {
        CUDA_CHECK(cudaMemset(d_out, 0, sz2d*8));
        
        bool launched = true;
        for (int r = 0; r < WARMUP; r++) {
            g_flush.flush();
            launched = launch_gpu_v(V, ci, d_out, d_vn_ie, d_inv_dual,
                d_w, d_cidx, d_z_vt_ie, d_inv_primal,
                d_tangent, d_z_w_v, d_vidx, N, nlev, nlev_end);
            if (!launched) break;
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        if (!launched) {
            printf("SKIP:        nlev=%d(%d) dist=%-12s V=%d cfg=%-14s  (does not fit)\n",
                   nlev, nlev_end, dist_label, V, GCFG[ci].label);
            continue;
        }

        CUDA_CHECK(cudaMemcpy(h_gpu_out, d_out, sz2d*8, cudaMemcpyDeviceToHost));
        int n_fail=0; double max_rel=0; size_t first_fail=0;
        bool ok = verify(h_gpu_out, h_ref, sz2d, 1e-8, 1e-12, &n_fail, &max_rel, &first_fail);
        if (!ok) {
            int ff_je, ff_jk;
            if (kV<=2){ff_je=first_fail%N;ff_jk=first_fail/N;}
            else      {ff_jk=first_fail%nlev;ff_je=first_fail/nlev;}
            printf("VERIFY FAIL: nlev=%d(%d) dist=%-12s V=%d cfg=%-14s  "
                   "fails=%d/%zu max_rel=%.3e\n"
                   "  first_fail: idx=%zu (je=%d,jk=%d) got=%.6e ref=%.6e\n",
                   nlev, nlev_end, dist_label, V, GCFG[ci].label,
                   n_fail, sz2d, max_rel, first_fail, ff_je, ff_jk,
                   h_gpu_out[first_fail], h_ref[first_fail]);
            continue;
        } else if (ci == 0) {
            printf("VERIFY OK:   nlev=%d(%d) dist=%-12s V=%d max_rel=%.3e\n",
                   nlev, nlev_end, dist_label, V, max_rel);
        }

        for (int r = 0; r < NRUNS; r++) {
            g_flush.flush();
            CUDA_CHECK(cudaEventRecord(ev0));
            launch_gpu_v(V, ci, d_out, d_vn_ie, d_inv_dual,
                d_w, d_cidx, d_z_vt_ie, d_inv_primal,
                d_tangent, d_z_w_v, d_vidx, N, nlev, nlev_end);
            CUDA_CHECK(cudaEventRecord(ev1));
            CUDA_CHECK(cudaEventSynchronize(ev1));
            float ms=0; CUDA_CHECK(cudaEventElapsedTime(&ms, ev0, ev1));
            fprintf(fcsv, "gpu,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%.6f\n",
                V, nlev_end, N, dist_label, GCFG[ci].label,
                GCFG[ci].tx, GCFG[ci].ty, GCFG[ci].bx, GCFG[ci].by,
                r, (double)ms);
            g_flush.flush();
        }
    }
    printf("Done: nlev=%d(%d)  dist=%-12s  V=%d\n", nlev, nlev_end, dist_label, V);
}

/* ================================================================ */
/*  run_dist_block -- allocate, run V_start..V_end, free             */
/* ================================================================ */
static void run_dist_block(
    FILE* fcsv, int N, int nlev, int nlev_end,
    const char* dist_label,
    int V_start, int V_end,
    int* cell_logical, int* vert_logical,
    double* icon_inv_dual, double* icon_inv_primal, double* icon_tangent)
{
    BenchData bd;
    bd.alloc(N, nlev);
    bd.fill(nlev);

    if (icon_inv_dual)
        for (int je = 0; je < N; je++) {
            bd.inv_dual[je]   = icon_inv_dual[je];
            bd.inv_primal[je] = icon_inv_primal[je];
            bd.tangent_o[je]  = icon_tangent[je];
        }

    size_t sz2d = bd.sz2d;
    double* h_ref     = new double[sz2d];
    double* h_gpu_out = new double[sz2d];

    double *d_vn_ie, *d_w, *d_z_vt_ie, *d_z_w_v, *d_out;
    double *d_inv_dual, *d_inv_primal, *d_tangent;
    int *d_cidx, *d_vidx;
    CUDA_CHECK(cudaMalloc(&d_vn_ie,      sz2d*8));
    CUDA_CHECK(cudaMalloc(&d_w,          sz2d*8));
    CUDA_CHECK(cudaMalloc(&d_z_vt_ie,    sz2d*8));
    CUDA_CHECK(cudaMalloc(&d_z_w_v,      sz2d*8));
    CUDA_CHECK(cudaMalloc(&d_out,        sz2d*8));
    CUDA_CHECK(cudaMalloc(&d_inv_dual,   N*8));
    CUDA_CHECK(cudaMalloc(&d_inv_primal, N*8));
    CUDA_CHECK(cudaMalloc(&d_tangent,    N*8));
    CUDA_CHECK(cudaMalloc(&d_cidx,       N*2*4));
    CUDA_CHECK(cudaMalloc(&d_vidx,       N*2*4));

    cudaEvent_t ev0, ev1;
    CUDA_CHECK(cudaEventCreate(&ev0));
    CUDA_CHECK(cudaEventCreate(&ev1));

    for (int V = V_start; V <= V_end; V++) {
        int kV = kern_v(V);
        bd.set_variant(kV, cell_logical, vert_logical);

        run_variant_configs(fcsv, V, N, nlev, nlev_end, dist_label,
            h_ref, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
            bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
            bd.h_vidx, sz2d,
            d_vn_ie, d_inv_dual, d_w, d_cidx,
            d_z_vt_ie, d_inv_primal, d_tangent, d_z_w_v,
            d_vidx, d_out, ev0, ev1, h_gpu_out);
        fflush(fcsv);
    }

    cudaFree(d_vn_ie); cudaFree(d_w); cudaFree(d_z_vt_ie);
    cudaFree(d_z_w_v); cudaFree(d_out);
    cudaFree(d_inv_dual); cudaFree(d_inv_primal); cudaFree(d_tangent);
    cudaFree(d_cidx); cudaFree(d_vidx);
    CUDA_CHECK(cudaEventDestroy(ev0));
    CUDA_CHECK(cudaEventDestroy(ev1));
    delete[] h_ref;
    delete[] h_gpu_out;
    bd.free_all();
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
    if (icon_nproma <= 0)
        fprintf(stderr, "WARNING: could not read nproma from '%s'\n", global_path.c_str());
    printf("Loading ICON data from: %s  (nproma=%d)\n", patch_path.c_str(), icon_nproma);

    IconEdgeData icon_ed;
    bool have_exact = (icon_nproma > 0) &&
                      icon_load_patch(patch_path.c_str(), icon_nproma, icon_ed);
    assert(have_exact);
    if (have_exact)
        printf("ICON exact data loaded: nproma=%d  n_edges=%d (valid=%d)  n_cells=%d  n_verts=%d\n",
               icon_ed.nproma, icon_ed.n_edges, icon_ed.n_edges_valid, icon_ed.n_cells, icon_ed.n_verts);
    else
        fprintf(stderr, "WARNING: failed to load ICON data, synthetic only\n");

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    printf("GPU: %s  SM count: %d\n", prop.name, prop.multiProcessorCount);
    printf("Configs: %d\n", N_GCFG);

    srand((unsigned)time(NULL));
    g_flush.init();
    printf("GPU flush buffers initialized (2 x %.0f MB)\n",
           (double)FLUSH_N*FLUSH_N*8/1e6);

    for (int nlev_i = 0; nlev_i < N_NLEVS; nlev_i++) {
        int nlev_base   = NLEVS[nlev_i];
        int nlev_padded = icon_pad_nlev(nlev_base);

        /* ============================================================ */
        /*  EXACT distribution                                          */
        /* ============================================================ */
        if (have_exact) {
            const int Ne = icon_ed.n_edges;
            if (icon_ed.n_cells > Ne || icon_ed.n_verts > Ne) {
                fprintf(stderr, "WARNING: n_cells/n_verts > Ne; skipping exact nlev=%d\n", nlev_base);
                continue;
            }

            int* ecl = new int[Ne*2];
            int* evl = new int[Ne*2];
            for (int i = 0; i < Ne*2; i++) {
                ecl[i] = icon_ed.cell_idx[i];
                evl[i] = icon_ed.vert_idx[i];
            }

            /* V1-V4: nlev = nlev_end = nlev_base */
            run_dist_block(fcsv, Ne, nlev_base, nlev_base,
                "exact", 1, 4,
                ecl, evl,
                icon_ed.inv_dual.data(), icon_ed.inv_primal.data(),
                icon_ed.tangent_o.data());

            /* V5: nlev = padded, nlev_end = nlev_base */
            if (nlev_padded != nlev_base) {
                run_dist_block(fcsv, Ne, nlev_padded, nlev_base,
                    "exact", 5, 5,
                    ecl, evl,
                    icon_ed.inv_dual.data(), icon_ed.inv_primal.data(),
                    icon_ed.tangent_o.data());
            }

            delete[] ecl;
            delete[] evl;
        }

        /* ============================================================ */
        /*  Synthetic distributions                                     */
        /* ============================================================ */
        for (int di = 0; di < 2; di++) {
            gen_cell_idx_logical(cell_logical, N, (CellDist)di, rng);

            /* V1-V4: nlev = nlev_end = nlev_base (no padding) */
            run_dist_block(fcsv, N, nlev_base, nlev_base,
                dist_name[di], 1, 4,
                cell_logical, vd.logical,
                nullptr, nullptr, nullptr);

            /* V5: nlev = padded, nlev_end = nlev_base */
            if (nlev_padded != nlev_base) {
                run_dist_block(fcsv, N, nlev_padded, nlev_base,
                    dist_name[di], 5, 5,
                    cell_logical, vd.logical,
                    nullptr, nullptr, nullptr);
            }

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
