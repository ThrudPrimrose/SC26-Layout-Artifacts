/*
 * bench_gpu.cu -- GPU z_v_grad_w stencil benchmark
 *
 * Separate array dimensions: N_e (edges), N_c (cells), N_v (verts).
 * Grid-stride loops; TX*TY <= 8.
 *
 * Compile (NVIDIA): nvcc -O3 -arch=sm_80 -std=c++17 -Xcompiler -fopenmp bench_gpu.cu -o bench_gpu
 * Compile (AMD):    hipcc -O3 -std=c++17 -fopenmp bench_gpu.cu -o bench_gpu
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
/*  CPU reference                                                    */
/* ================================================================ */
template<int V>
static void cpu_reference(
    double* __restrict__ out,
    const double* __restrict__ vn_ie, const double* __restrict__ inv_dual,
    const double* __restrict__ w, const int* __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int* __restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    for (int jk = 0; jk < nlev_end; jk++)
        for (int je = 0; je < N_e; je++) { STENCIL_BODY(V) }
}
static void cpu_reference_v(int V,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    int kV = kern_v(V);
    switch (kV) {
    case 1: cpu_reference<1>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end); break;
    case 2: cpu_reference<2>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end); break;
    case 3: cpu_reference<3>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end); break;
    case 4: cpu_reference<4>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end); break;
    }
}

/* ================================================================ */
/*  Verify                                                           */
/* ================================================================ */
static bool verify(const double* got, const double* ref, size_t n,
                   double rtol, double atol,
                   int* n_fail, double* max_rel, size_t* first_fail_idx)
{
    *n_fail=0; *max_rel=0.0; *first_fail_idx=0;
    for (size_t i=0;i<n;i++) {
        double diff=std::abs(got[i]-ref[i]);
        double denom=std::max(std::abs(ref[i]),1e-300);
        double rel=diff/denom;
        if (rel>*max_rel) *max_rel=rel;
        if (diff>atol+rtol*std::abs(ref[i])) {
            if (*n_fail==0) *first_fail_idx=i;
            (*n_fail)++;
        }
    }
    return *n_fail==0;
}

/* ================================================================ */
/*  GPU kernels — grid-stride, N_e/N_c/N_v                          */
/* ================================================================ */
template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_je_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie, const double* __restrict__ inv_dual,
    const double* __restrict__ w, const int* __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int* __restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    const int je_stride = (int)gridDim.x * BX * TX;
    const int jk_stride = (int)gridDim.y * BY * TY;
    for (int je_base0 = ((int)blockIdx.x*BX+(int)threadIdx.x)*TX;
         je_base0 < N_e; je_base0 += je_stride)
    {
        int ci0_a[TX],ci1_a[TX],vi0_a[TX],vi1_a[TX];
        double id_a[TX],ip_a[TX],tg_a[TX];
        #pragma unroll
        for (int tx=0;tx<TX;tx++){int je=je_base0+tx;if(je<N_e){
            ci0_a[tx]=cell_idx[IN<V>(je,0,N_e)]; ci1_a[tx]=cell_idx[IN<V>(je,1,N_e)];
            vi0_a[tx]=vert_idx[IN<V>(je,0,N_e)]; vi1_a[tx]=vert_idx[IN<V>(je,1,N_e)];
            id_a[tx]=inv_dual[je]; ip_a[tx]=inv_primal[je]; tg_a[tx]=tangent[je];
        }}
        for (int jk_base0 = ((int)blockIdx.y*BY+(int)threadIdx.y)*TY;
             jk_base0 < nlev_end; jk_base0 += jk_stride)
        {
            #pragma unroll
            for (int ty=0;ty<TY;ty++){int jk=jk_base0+ty;if(jk>=nlev_end)continue;
            #pragma unroll
            for (int tx=0;tx<TX;tx++){int je=je_base0+tx;if(je>=N_e)continue;
                int c2d=IC<V>(je,jk,N_e,nlev);
                out[c2d]=vn_ie[c2d]*id_a[tx]*
                    (w[IC<V>(ci0_a[tx],jk,N_c,nlev)]-w[IC<V>(ci1_a[tx],jk,N_c,nlev)])
                  +z_vt_ie[c2d]*ip_a[tx]*tg_a[tx]*
                    (z_w_v[IC<V>(vi0_a[tx],jk,N_v,nlev)]-z_w_v[IC<V>(vi1_a[tx],jk,N_v,nlev)]);
            }}
        }
    }
}

template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_jk_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie, const double* __restrict__ inv_dual,
    const double* __restrict__ w, const int* __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int* __restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    const int jk_stride = (int)gridDim.x * BX * TX;
    const int je_stride = (int)gridDim.y * BY * TY;
    for (int je_base0 = ((int)blockIdx.y*BY+(int)threadIdx.y)*TY;
         je_base0 < N_e; je_base0 += je_stride)
    {
        int ci0_a[TY],ci1_a[TY],vi0_a[TY],vi1_a[TY];
        double id_a[TY],ip_a[TY],tg_a[TY];
        #pragma unroll
        for (int ty=0;ty<TY;ty++){int je=je_base0+ty;if(je<N_e){
            ci0_a[ty]=cell_idx[IN<V>(je,0,N_e)]; ci1_a[ty]=cell_idx[IN<V>(je,1,N_e)];
            vi0_a[ty]=vert_idx[IN<V>(je,0,N_e)]; vi1_a[ty]=vert_idx[IN<V>(je,1,N_e)];
            id_a[ty]=inv_dual[je]; ip_a[ty]=inv_primal[je]; tg_a[ty]=tangent[je];
        }}
        #pragma unroll
        for (int ty=0;ty<TY;ty++){int je=je_base0+ty;if(je>=N_e)continue;
        #pragma unroll
        for (int tx=0;tx<TX;tx++){int jk=jk_base0+tx;if(jk>=nlev_end)continue;
            int c2d=IC<V>(je,jk,N_e,nlev);
            out[c2d]=vn_ie[c2d]*id_a[ty]*
                (w[IC<V>(ci0_a[ty],jk,N_c,nlev)]-w[IC<V>(ci1_a[ty],jk,N_c,nlev)])
                +z_vt_ie[c2d]*ip_a[ty]*tg_a[ty]*
                (z_w_v[IC<V>(vi0_a[ty],jk,N_v,nlev)]-z_w_v[IC<V>(vi1_a[ty],jk,N_v,nlev)]);
        }}
    }
}

/* ================================================================ */
/*  V6: Flat 1D schedule, jk-first layout (V4)                      */
/* ================================================================ */
template<unsigned BX, int V>
__global__ void gpu_kernel_v6_flat(
    double* __restrict__ out,
    const double* __restrict__ vn_ie, const double* __restrict__ inv_dual,
    const double* __restrict__ w, const int* __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int* __restrict__ vert_idx,
    unsigned N_e, unsigned N_c, unsigned N_v,
    unsigned nlev, unsigned nlev_end, unsigned nlev_padded)
{
    unsigned bid = (unsigned)blockIdx.y * gridDim.x + (unsigned)blockIdx.x;
    unsigned gid = bid * BX + threadIdx.x;
    unsigned jk = gid % nlev_padded;
    unsigned je = gid / nlev_padded;
    if (je >= N_e || jk >= nlev_end) return;
    int c2d = IC<V>((int)je,(int)jk,(int)N_e,(int)nlev);
    int ci0 = cell_idx[IN<V>((int)je,0,(int)N_e)];
    int ci1 = cell_idx[IN<V>((int)je,1,(int)N_e)];
    int vi0 = vert_idx[IN<V>((int)je,0,(int)N_e)];
    int vi1 = vert_idx[IN<V>((int)je,1,(int)N_e)];
    out[c2d] = vn_ie[c2d]*inv_dual[je]*
        (w[IC<V>(ci0,(int)jk,(int)N_c,(int)nlev)]-w[IC<V>(ci1,(int)jk,(int)N_c,(int)nlev)])
      + z_vt_ie[c2d]*inv_primal[je]*tangent[je]*
        (z_w_v[IC<V>(vi0,(int)jk,(int)N_v,(int)nlev)]-z_w_v[IC<V>(vi1,(int)jk,(int)N_v,(int)nlev)]);
}

/* ================================================================ */
/*  V7: Flat 1D with TX elements per thread, jk-first, UNPADDED     */
/*  gid = (bid*BX + tid)*TX + t  →  jk = gid % nlev, je = gid / nlev */
/*  Consecutive gids sharing same je reuse index lookups.            */
/*  No nlev padding — gid maps 1:1 to memory address in V4 layout.  */
/* ================================================================ */
template<unsigned BX, int TX, int V>
__global__ void gpu_kernel_v7_flat(
    double* __restrict__ out,
    const double* __restrict__ vn_ie, const double* __restrict__ inv_dual,
    const double* __restrict__ w, const int* __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int* __restrict__ vert_idx,
    unsigned N_e, unsigned N_c, unsigned N_v, unsigned nlev, unsigned total)
{
    unsigned bid=(unsigned)blockIdx.y*gridDim.x+(unsigned)blockIdx.x;
    unsigned gid_base=(bid*BX+threadIdx.x)*(unsigned)TX;
    #pragma unroll
    for(int t=0;t<TX;t++){
        unsigned gid=gid_base+(unsigned)t;
        if(gid>=total) return;
        unsigned jk=gid%nlev, je=gid/nlev;
        int c2d=(int)gid;
        int ci0=cell_idx[IN<V>((int)je,0,(int)N_e)],ci1=cell_idx[IN<V>((int)je,1,(int)N_e)];
        int vi0=vert_idx[IN<V>((int)je,0,(int)N_e)],vi1=vert_idx[IN<V>((int)je,1,(int)N_e)];
        out[c2d]=vn_ie[c2d]*inv_dual[je]*
            (w[IC<V>(ci0,(int)jk,(int)N_c,(int)nlev)]-w[IC<V>(ci1,(int)jk,(int)N_c,(int)nlev)])
          +z_vt_ie[c2d]*inv_primal[je]*tangent[je]*
            (z_w_v[IC<V>(vi0,(int)jk,(int)N_v,(int)nlev)]-z_w_v[IC<V>(vi1,(int)jk,(int)N_v,(int)nlev)]);
    }
}

/* ================================================================ */
/*  Config table — TX*TY <= 8, pruned                                */
/* ================================================================ */
struct GpuCfg { int tx,ty,bx,by; const char* label; };
static constexpr GpuCfg GCFG[] = {
    /* 1D, TX>=TY */
    {1,1,256,1,"1x1_256x1"}, {1,1,128,1,"1x1_128x1"}, {1,1,64,1,"1x1_64x1"},
    {2,1,128,1,"2x1_128x1"}, {4,1,128,1,"4x1_128x1"}, {4,1,64,1,"4x1_64x1"},
    {8,1,128,1,"8x1_128x1"}, {8,1,64,1,"8x1_64x1"},
    /* 1D, TY>TX */
    {1,2,128,1,"1x2_128x1"}, {1,4,128,1,"1x4_128x1"}, {1,8,128,1,"1x8_128x1"},
    {1,2,64,1,"1x2_64x1"}, {1,4,64,1,"1x4_64x1"}, {1,8,64,1,"1x8_64x1"},
    {1,2,32,1,"1x2_32x1"}, {1,4,32,1,"1x4_32x1"}, {1,8,32,1,"1x8_32x1"},
    {1,4,16,1,"1x4_16x1"}, {1,8,16,1,"1x8_16x1"},
    {2,2,128,1,"2x2_128x1"}, {2,4,128,1,"2x4_128x1"},
    /* 2D, BX=32 */
    {1,1,32,32,"1x1_32x32"}, {1,1,32,16,"1x1_32x16"}, {1,1,32,8,"1x1_32x8"},
    {2,1,32,8,"2x1_32x8"}, {4,1,32,8,"4x1_32x8"}, {4,1,32,4,"4x1_32x4"},
    {1,2,32,16,"1x2_32x16"}, {1,4,32,8,"1x4_32x8"}, {1,8,32,4,"1x8_32x4"},
    {2,2,32,8,"2x2_32x8"}, {2,4,32,4,"2x4_32x4"},
    /* 2D, BX=64 */
    {1,1,64,8,"1x1_64x8"}, {1,1,64,4,"1x1_64x4"}, {1,1,64,2,"1x1_64x2"},
    {2,1,64,4,"2x1_64x4"}, {4,1,64,2,"4x1_64x2"},
    {1,2,64,4,"1x2_64x4"}, {1,4,64,4,"1x4_64x4"}, {1,8,64,2,"1x8_64x2"},
    {2,2,64,4,"2x2_64x4"},
    /* 2D, BX=128 */
    {1,1,128,4,"1x1_128x4"}, {1,1,128,2,"1x1_128x2"},
    {1,2,128,2,"1x2_128x2"}, {1,4,128,2,"1x4_128x2"}, {1,8,128,2,"1x8_128x2"},
    /* 2D, BX=16 */
    {1,1,16,16,"1x1_16x16"}, {2,1,16,16,"2x1_16x16"}, {1,2,16,16,"1x2_16x16"},
    {1,4,16,16,"1x4_16x16"}, {1,8,16,16,"1x8_16x16"}, {4,2,16,16,"4x2_16x16"},
    /* 96-wide */
    {1,1,96,1,"1x1_96x1"}, {1,2,96,1,"1x2_96x1"}, {1,4,96,1,"1x4_96x1"}, {1,8,96,1,"1x8_96x1"},
    {1,1,96,2,"1x1_96x2"}, {1,2,96,2,"1x2_96x2"}, {1,4,96,2,"1x4_96x2"},
    {1,1,96,4,"1x1_96x4"}, {1,8,96,4,"1x8_96x4"}, {1,8,96,8,"1x8_96x8"},
    /* Remaining useful */
    {1,8,256,1,"1x8_256x1"}, {2,4,16,16,"2x4_16x16"},

    {4,1,16,16,"4x1_16x16"}, {2,2,16,16,"2x2_16x16"},
};
static constexpr int N_GCFG = sizeof(GCFG)/sizeof(GCFG[0]);

/* V6 configs */
struct GpuCfgV6 { unsigned bx; const char* label; };
static constexpr GpuCfgV6 GCFG_V6[] = {
    {32,"v6_32"}, {64,"v6_64"}, {96,"v6_96"}, {128,"v6_128"},
    {192,"v6_192"}, {256,"v6_256"}, {384,"v6_384"}, {512,"v6_512"}, {1024,"v6_1024"},
};
static constexpr int N_GCFG_V6 = sizeof(GCFG_V6)/sizeof(GCFG_V6[0]);

/* V7 configs: flat + TX elements per thread, unpadded nlev */
struct GpuCfgV7 { unsigned bx; int tx; const char* label; };
static constexpr GpuCfgV7 GCFG_V7[] = {
    { 64,2,"v7_64_tx2"},  { 64,4,"v7_64_tx4"},
    {128,2,"v7_128_tx2"}, {128,4,"v7_128_tx4"}, {128,8,"v7_128_tx8"},
    {256,2,"v7_256_tx2"}, {256,4,"v7_256_tx4"}, {256,8,"v7_256_tx8"},
    {512,2,"v7_512_tx2"}, {512,4,"v7_512_tx4"},
};
static constexpr int N_GCFG_V7 = sizeof(GCFG_V7)/sizeof(GCFG_V7[0]);

/* ================================================================ */
/*  Launch dispatch V1-V5                                            */
/* ================================================================ */
template<int V>
static bool launch_gpu(int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v, const int* vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    bool launched = false;
    #define CLAMP65535(x) ((x)>65535u ? 65535u : (x))
    #define LG_JE(TX_,TY_,BX_,BY_) do { \
        dim3 blk(BX_,BY_); \
        dim3 grd(CLAMP65535(((unsigned)N_e+(BX_)*(TX_)-1)/((BX_)*(TX_))), \
                 CLAMP65535(((unsigned)nlev_end+(BY_)*(TY_)-1)/((BY_)*(TY_)))); \
        gpu_kernel_je_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>( \
            out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent, \
            z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end); \
        CUDA_LAUNCH_CHECK(); launched=true; \
    } while(0)
    #define LG_JK(TX_,TY_,BX_,BY_) do { \
        dim3 blk(BX_,BY_); \
        dim3 grd(CLAMP65535(((unsigned)nlev_end+(BX_)*(TX_)-1)/((BX_)*(TX_))), \
                 CLAMP65535(((unsigned)N_e+(BY_)*(TY_)-1)/((BY_)*(TY_)))); \
        gpu_kernel_jk_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>( \
            out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent, \
            z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end); \
        CUDA_LAUNCH_CHECK(); launched=true; \
    } while(0)
    #define LG(TX_,TY_,BX_,BY_) do { \
        if constexpr (V<=2) { LG_JE(TX_,TY_,BX_,BY_); } \
        else                { LG_JK(TX_,TY_,BX_,BY_); } \
    } while(0)

    switch(cfg){
    case  0: LG(1,1,256,1); break;
    case  1: LG(1,1,128,1); break;
    case  2: LG(1,1,64,1); break;
    case  3: LG(2,1,128,1); break;
    case  4: LG(4,1,128,1); break;
    case  5: LG(4,1,64,1); break;
    case  6: LG(8,1,128,1); break;
    case  7: LG(8,1,64,1); break;
    case  8: LG(1,2,128,1); break;
    case  9: LG(1,4,128,1); break;
    case 10: LG(1,8,128,1); break;
    case 11: LG(1,2,64,1); break;
    case 12: LG(1,4,64,1); break;
    case 13: LG(1,8,64,1); break;
    case 14: LG(1,2,32,1); break;
    case 15: LG(1,4,32,1); break;
    case 16: LG(1,8,32,1); break;
    case 17: LG(1,4,16,1); break;
    case 18: LG(1,8,16,1); break;
    case 19: LG(2,2,128,1); break;
    case 20: LG(2,4,128,1); break;
    case 21: LG(1,1,32,32); break;
    case 22: LG(1,1,32,16); break;
    case 23: LG(1,1,32,8); break;
    case 24: LG(2,1,32,8); break;
    case 25: LG(4,1,32,8); break;
    case 26: LG(4,1,32,4); break;
    case 27: LG(1,2,32,16); break;
    case 28: LG(1,4,32,8); break;
    case 29: LG(1,8,32,4); break;
    case 30: LG(2,2,32,8); break;
    case 31: LG(2,4,32,4); break;
    case 32: LG(1,1,64,8); break;
    case 33: LG(1,1,64,4); break;
    case 34: LG(1,1,64,2); break;
    case 35: LG(2,1,64,4); break;
    case 36: LG(4,1,64,2); break;
    case 37: LG(1,2,64,4); break;
    case 38: LG(1,4,64,4); break;
    case 39: LG(1,8,64,2); break;
    case 40: LG(2,2,64,4); break;
    case 41: LG(1,1,128,4); break;
    case 42: LG(1,1,128,2); break;
    case 43: LG(1,2,128,2); break;
    case 44: LG(1,4,128,2); break;
    case 45: LG(1,8,128,2); break;
    case 46: LG(1,1,16,16); break;
    case 47: LG(2,1,16,16); break;
    case 48: LG(1,2,16,16); break;
    case 49: LG(1,4,16,16); break;
    case 50: LG(1,8,16,16); break;
    case 51: LG(4,2,16,16); break;
    case 52: LG(1,1,96,1); break;
    case 53: LG(1,2,96,1); break;
    case 54: LG(1,4,96,1); break;
    case 55: LG(1,8,96,1); break;
    case 56: LG(1,1,96,2); break;
    case 57: LG(1,2,96,2); break;
    case 58: LG(1,4,96,2); break;
    case 59: LG(1,1,96,4); break;
    case 60: LG(1,8,96,4); break;
    case 61: LG(1,8,96,8); break;
    case 62: LG(1,8,256,1); break;
    case 63: LG(2,4,16,16); break;
    case 64: LG(4,1,16,16); break;
    case 65: LG(4,2,16,16); break;
    }
    #undef LG
    #undef LG_JE
    #undef LG_JK
    #undef CLAMP65535
    return launched;
}

static bool launch_gpu_v(int V, int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v, const int* vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    int kV = kern_v(V);
    switch(kV){
    case 1: return launch_gpu<1>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end);
    case 2: return launch_gpu<2>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end);
    case 3: return launch_gpu<3>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end);
    case 4: return launch_gpu<4>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev,nlev_end);
    }
    return false;
}

/* ================================================================ */
/*  V6 launch                                                        */
/* ================================================================ */
template<int V>
static bool launch_gpu_v6(int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v, const int* vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end, unsigned nlev_padded)
{
    if (cfg<0||cfg>=N_GCFG_V6) return false;
    unsigned BX=GCFG_V6[cfg].bx;
    unsigned total=nlev_padded*(unsigned)N_e;
    unsigned grd_linear=(total+BX-1u)/BX;
    dim3 blk(BX); dim3 grd;
    if (grd_linear<=65535u) grd=dim3(grd_linear,1);
    else { grd.y=(grd_linear+65534u)/65535u; grd.x=(grd_linear+grd.y-1u)/grd.y; }
    #define LAUNCH_V6(BX_) gpu_kernel_v6_flat<BX_,V><<<grd,blk>>>( \
        out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx, \
        (unsigned)N_e,(unsigned)N_c,(unsigned)N_v,(unsigned)nlev,(unsigned)nlev_end,nlev_padded); \
        CUDA_LAUNCH_CHECK()
    switch(BX){
    case   32: LAUNCH_V6(  32); break; case   64: LAUNCH_V6(  64); break;
    case   96: LAUNCH_V6(  96); break; case  128: LAUNCH_V6( 128); break;
    case  192: LAUNCH_V6( 192); break; case  256: LAUNCH_V6( 256); break;
    case  384: LAUNCH_V6( 384); break; case  512: LAUNCH_V6( 512); break;
    case 1024: LAUNCH_V6(1024); break; default: return false;
    }
    #undef LAUNCH_V6
    return true;
}

/* ================================================================ */
/*  V7 launch — flat with TX elements per thread, unpadded           */
/* ================================================================ */
template<int V>
static bool launch_gpu_v7(int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v, const int* vert_idx,
    int N_e, int N_c, int N_v, int nlev)
{
    if (cfg<0||cfg>=N_GCFG_V7) return false;
    unsigned BX=GCFG_V7[cfg].bx;
    int TX=GCFG_V7[cfg].tx;
    unsigned total=(unsigned)nlev*(unsigned)N_e;
    unsigned nthreads=(total+(unsigned)TX-1u)/(unsigned)TX;
    unsigned grd_linear=(nthreads+BX-1u)/BX;
    dim3 blk(BX); dim3 grd;
    if (grd_linear<=65535u) grd=dim3(grd_linear,1);
    else { grd.y=(grd_linear+65534u)/65535u; grd.x=(grd_linear+grd.y-1u)/grd.y; }

    #define LAUNCH_V7(BX_,TX_) gpu_kernel_v7_flat<BX_,TX_,V><<<grd,blk>>>( \
        out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx, \
        (unsigned)N_e,(unsigned)N_c,(unsigned)N_v,(unsigned)nlev,total); \
        CUDA_LAUNCH_CHECK()

    /* dispatch on (BX, TX) pairs from config table */
    if      (BX== 64&&TX==2) { LAUNCH_V7( 64,2); }
    else if (BX== 64&&TX==4) { LAUNCH_V7( 64,4); }
    else if (BX==128&&TX==2) { LAUNCH_V7(128,2); }
    else if (BX==128&&TX==4) { LAUNCH_V7(128,4); }
    else if (BX==128&&TX==8) { LAUNCH_V7(128,8); }
    else if (BX==256&&TX==2) { LAUNCH_V7(256,2); }
    else if (BX==256&&TX==4) { LAUNCH_V7(256,4); }
    else if (BX==256&&TX==8) { LAUNCH_V7(256,8); }
    else if (BX==512&&TX==2) { LAUNCH_V7(512,2); }
    else if (BX==512&&TX==4) { LAUNCH_V7(512,4); }
    else { return false; }
    #undef LAUNCH_V7
    return true;
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
/*  run_variant_configs  (V1-V5)                                     */
/* ================================================================ */
static void run_variant_configs(
    FILE* fcsv, int V, int N_e, int N_c, int N_v, int nlev, int nlev_end,
    const char* dist_label,
    double* h_ref, BenchData& bd, size_t sz_e,
    double* d_vn_ie, double* d_inv_dual, double* d_w, int* d_cidx,
    double* d_z_vt_ie, double* d_inv_primal, double* d_tangent,
    double* d_z_w_v, int* d_vidx, double* d_out,
    cudaEvent_t ev0, cudaEvent_t ev1, double* h_gpu_out)
{
    memset(h_ref, 0, sz_e*sizeof(double));
    cpu_reference_v(V, h_ref, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
        bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v, bd.h_vidx,
        N_e, N_c, N_v, nlev, nlev_end);

    CUDA_CHECK(cudaMemcpy(d_vn_ie,   bd.h_vn_ie,   bd.sz_e*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_w,       bd.h_w,       bd.sz_c*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_vt_ie, bd.h_z_vt_ie, bd.sz_e*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_w_v,   bd.h_z_w_v,   bd.sz_v*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cidx,    bd.h_cidx,    N_e*2*4,   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vidx,    bd.h_vidx,    N_e*2*4,   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_dual,   bd.inv_dual,   N_e*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_primal, bd.inv_primal, N_e*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_tangent,    bd.tangent_o,  N_e*8, cudaMemcpyHostToDevice));

    for (int ci=0; ci<N_GCFG; ci++) {
        CUDA_CHECK(cudaMemset(d_out, 0, sz_e*8));
        bool launched=true;
        for (int r=0;r<WARMUP;r++) {
            g_flush.flush();
            launched=launch_gpu_v(V,ci,d_out,d_vn_ie,d_inv_dual,d_w,d_cidx,
                d_z_vt_ie,d_inv_primal,d_tangent,d_z_w_v,d_vidx,N_e,N_c,N_v,nlev,nlev_end);
            if(!launched)break;
            CUDA_CHECK(cudaDeviceSynchronize());
        }
        if(!launched){
            printf("SKIP: nlev=%d(%d) dist=%-12s V=%d cfg=%-14s\n",
                   nlev,nlev_end,dist_label,V,GCFG[ci].label);
            continue;
        }
        CUDA_CHECK(cudaMemcpy(h_gpu_out,d_out,sz_e*8,cudaMemcpyDeviceToHost));
        int nf=0;double mr=0;size_t ff=0;
        bool ok=verify(h_gpu_out,h_ref,sz_e,1e-8,1e-12,&nf,&mr,&ff);
        if(!ok){
            printf("FAIL: V=%d cfg=%-14s fails=%d max_rel=%.3e\n",V,GCFG[ci].label,nf,mr);
            continue;
        } else if(ci==0) printf("OK:   nlev=%d(%d) dist=%-12s V=%d mr=%.3e\n",nlev,nlev_end,dist_label,V,mr);
        for (int r=0;r<NRUNS;r++){
            g_flush.flush();
            CUDA_CHECK(cudaEventRecord(ev0));
            launch_gpu_v(V,ci,d_out,d_vn_ie,d_inv_dual,d_w,d_cidx,
                d_z_vt_ie,d_inv_primal,d_tangent,d_z_w_v,d_vidx,N_e,N_c,N_v,nlev,nlev_end);
            CUDA_CHECK(cudaEventRecord(ev1));
            CUDA_CHECK(cudaEventSynchronize(ev1));
            float ms=0; CUDA_CHECK(cudaEventElapsedTime(&ms,ev0,ev1));
            fprintf(fcsv,"gpu,%d,%d,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%.6f\n",
                V,nlev_end,N_e,N_c,N_v,dist_label,GCFG[ci].label,
                GCFG[ci].tx,GCFG[ci].ty,GCFG[ci].bx,GCFG[ci].by,r,(double)ms);
            g_flush.flush();
        }
    }
    printf("Done: nlev=%d(%d) dist=%-12s V=%d\n",nlev,nlev_end,dist_label,V);
}

/* ================================================================ */
/*  run_v6_configs                                                   */
/* ================================================================ */
static void run_v6_configs(
    FILE* fcsv, int N_e, int N_c, int N_v, int nlev, int nlev_end, unsigned nlev_padded,
    const char* dist_label,
    double* h_ref, BenchData& bd, size_t sz_e,
    double* d_vn_ie, double* d_inv_dual, double* d_w, int* d_cidx,
    double* d_z_vt_ie, double* d_inv_primal, double* d_tangent,
    double* d_z_w_v, int* d_vidx, double* d_out,
    cudaEvent_t ev0, cudaEvent_t ev1, double* h_gpu_out)
{
    memset(h_ref,0,sz_e*sizeof(double));
    cpu_reference<4>(h_ref,bd.h_vn_ie,bd.inv_dual,bd.h_w,bd.h_cidx,
        bd.h_z_vt_ie,bd.inv_primal,bd.tangent_o,bd.h_z_w_v,bd.h_vidx,
        N_e,N_c,N_v,nlev,nlev_end);
    CUDA_CHECK(cudaMemcpy(d_vn_ie,bd.h_vn_ie,bd.sz_e*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_w,bd.h_w,bd.sz_c*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_vt_ie,bd.h_z_vt_ie,bd.sz_e*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_w_v,bd.h_z_w_v,bd.sz_v*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cidx,bd.h_cidx,N_e*2*4,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vidx,bd.h_vidx,N_e*2*4,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_dual,bd.inv_dual,N_e*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_primal,bd.inv_primal,N_e*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_tangent,bd.tangent_o,N_e*8,cudaMemcpyHostToDevice));
    for (int ci=0;ci<N_GCFG_V6;ci++){
        CUDA_CHECK(cudaMemset(d_out,0,sz_e*8));
        bool launched=true;
        for(int r=0;r<WARMUP;r++){g_flush.flush();
            launched=launch_gpu_v6<4>(ci,d_out,d_vn_ie,d_inv_dual,d_w,d_cidx,
                d_z_vt_ie,d_inv_primal,d_tangent,d_z_w_v,d_vidx,N_e,N_c,N_v,nlev,nlev_end,nlev_padded);
            if(!launched)break; CUDA_CHECK(cudaDeviceSynchronize());}
        if(!launched){printf("SKIP: V=6 cfg=%s\n",GCFG_V6[ci].label);continue;}
        CUDA_CHECK(cudaMemcpy(h_gpu_out,d_out,sz_e*8,cudaMemcpyDeviceToHost));
        int nf=0;double mr=0;size_t ff=0;
        bool ok=verify(h_gpu_out,h_ref,sz_e,1e-8,1e-12,&nf,&mr,&ff);
        if(!ok){printf("FAIL: V=6 cfg=%s fails=%d mr=%.3e\n",GCFG_V6[ci].label,nf,mr);continue;}
        else if(ci==0) printf("OK:   V=6 dist=%-12s mr=%.3e\n",dist_label,mr);
        for(int r=0;r<NRUNS;r++){
            g_flush.flush(); CUDA_CHECK(cudaEventRecord(ev0));
            launch_gpu_v6<4>(ci,d_out,d_vn_ie,d_inv_dual,d_w,d_cidx,
                d_z_vt_ie,d_inv_primal,d_tangent,d_z_w_v,d_vidx,N_e,N_c,N_v,nlev,nlev_end,nlev_padded);
            CUDA_CHECK(cudaEventRecord(ev1)); CUDA_CHECK(cudaEventSynchronize(ev1));
            float ms=0; CUDA_CHECK(cudaEventElapsedTime(&ms,ev0,ev1));
            fprintf(fcsv,"gpu,6,%d,%d,%d,%d,%s,%s,1,1,%d,1,%d,%.6f\n",
                nlev_end,N_e,N_c,N_v,dist_label,GCFG_V6[ci].label,(int)GCFG_V6[ci].bx,r,(double)ms);
            g_flush.flush();
        }
    }
    printf("Done: nlev=%d(%d) dist=%-12s V=6\n",nlev,nlev_end,dist_label);
}

/* ================================================================ */
/*  run_dist_block                                                   */
/* ================================================================ */
static void run_dist_block(
    FILE* fcsv, int N_e, int N_c, int N_v, int nlev, int nlev_end,
    const char* dist_label, int V_start, int V_end,
    int* cell_logical, int* vert_logical,
    double* icon_inv_dual, double* icon_inv_primal, double* icon_tangent)
{
    BenchData bd; bd.alloc(N_e,N_c,N_v,nlev); bd.fill(nlev);
    if(icon_inv_dual) for(int je=0;je<N_e;je++){
        bd.inv_dual[je]=icon_inv_dual[je]; bd.inv_primal[je]=icon_inv_primal[je]; bd.tangent_o[je]=icon_tangent[je];}

    double* h_ref=new double[bd.sz_e]; double* h_gpu_out=new double[bd.sz_e];
    double *d_vn,*d_w,*d_vt,*d_zw,*d_out,*d_id,*d_ip,*d_tg; int *d_ci,*d_vi;
    CUDA_CHECK(cudaMalloc(&d_vn,bd.sz_e*8)); CUDA_CHECK(cudaMalloc(&d_w,bd.sz_c*8));
    CUDA_CHECK(cudaMalloc(&d_vt,bd.sz_e*8)); CUDA_CHECK(cudaMalloc(&d_zw,bd.sz_v*8));
    CUDA_CHECK(cudaMalloc(&d_out,bd.sz_e*8));
    CUDA_CHECK(cudaMalloc(&d_id,N_e*8)); CUDA_CHECK(cudaMalloc(&d_ip,N_e*8)); CUDA_CHECK(cudaMalloc(&d_tg,N_e*8));
    CUDA_CHECK(cudaMalloc(&d_ci,N_e*2*4)); CUDA_CHECK(cudaMalloc(&d_vi,N_e*2*4));
    cudaEvent_t ev0,ev1; CUDA_CHECK(cudaEventCreate(&ev0)); CUDA_CHECK(cudaEventCreate(&ev1));

    for(int V=V_start;V<=V_end;V++){
        if(V==2||V==3) continue;
        bd.set_variant(kern_v(V),cell_logical,vert_logical);
        run_variant_configs(fcsv,V,N_e,N_c,N_v,nlev,nlev_end,dist_label,
            h_ref,bd,bd.sz_e,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,d_out,ev0,ev1,h_gpu_out);
        fflush(fcsv);
    }
    cudaFree(d_vn);cudaFree(d_w);cudaFree(d_vt);cudaFree(d_zw);cudaFree(d_out);
    cudaFree(d_id);cudaFree(d_ip);cudaFree(d_tg);cudaFree(d_ci);cudaFree(d_vi);
    CUDA_CHECK(cudaEventDestroy(ev0)); CUDA_CHECK(cudaEventDestroy(ev1));
    delete[]h_ref; delete[]h_gpu_out; bd.free_all();
}

static void run_dist_block_v6(
    FILE* fcsv, int N_e, int N_c, int N_v, int nlev, int nlev_end, unsigned nlev_padded,
    const char* dist_label, int* cell_logical, int* vert_logical,
    double* icon_inv_dual, double* icon_inv_primal, double* icon_tangent)
{
    BenchData bd; bd.alloc(N_e,N_c,N_v,nlev); bd.fill(nlev);
    if(icon_inv_dual) for(int je=0;je<N_e;je++){
        bd.inv_dual[je]=icon_inv_dual[je]; bd.inv_primal[je]=icon_inv_primal[je]; bd.tangent_o[je]=icon_tangent[je];}
    bd.set_variant(4,cell_logical,vert_logical);

    double* h_ref=new double[bd.sz_e]; double* h_gpu_out=new double[bd.sz_e];
    double *d_vn,*d_w,*d_vt,*d_zw,*d_out,*d_id,*d_ip,*d_tg; int *d_ci,*d_vi;
    CUDA_CHECK(cudaMalloc(&d_vn,bd.sz_e*8)); CUDA_CHECK(cudaMalloc(&d_w,bd.sz_c*8));
    CUDA_CHECK(cudaMalloc(&d_vt,bd.sz_e*8)); CUDA_CHECK(cudaMalloc(&d_zw,bd.sz_v*8));
    CUDA_CHECK(cudaMalloc(&d_out,bd.sz_e*8));
    CUDA_CHECK(cudaMalloc(&d_id,N_e*8)); CUDA_CHECK(cudaMalloc(&d_ip,N_e*8)); CUDA_CHECK(cudaMalloc(&d_tg,N_e*8));
    CUDA_CHECK(cudaMalloc(&d_ci,N_e*2*4)); CUDA_CHECK(cudaMalloc(&d_vi,N_e*2*4));
    cudaEvent_t ev0,ev1; CUDA_CHECK(cudaEventCreate(&ev0)); CUDA_CHECK(cudaEventCreate(&ev1));

    run_v6_configs(fcsv,N_e,N_c,N_v,nlev,nlev_end,nlev_padded,dist_label,
        h_ref,bd,bd.sz_e,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,d_out,ev0,ev1,h_gpu_out);
    fflush(fcsv);
    cudaFree(d_vn);cudaFree(d_w);cudaFree(d_vt);cudaFree(d_zw);cudaFree(d_out);
    cudaFree(d_id);cudaFree(d_ip);cudaFree(d_tg);cudaFree(d_ci);cudaFree(d_vi);
    CUDA_CHECK(cudaEventDestroy(ev0)); CUDA_CHECK(cudaEventDestroy(ev1));
    delete[]h_ref; delete[]h_gpu_out; bd.free_all();
}

/* ================================================================ */
/*  run_v7_configs — flat + TX, unpadded                             */
/* ================================================================ */
static void run_v7_configs(
    FILE* fcsv, int N_e, int N_c, int N_v, int nlev,
    const char* dist_label,
    double* h_ref, BenchData& bd, size_t sz_e,
    double* d_vn_ie, double* d_inv_dual, double* d_w, int* d_cidx,
    double* d_z_vt_ie, double* d_inv_primal, double* d_tangent,
    double* d_z_w_v, int* d_vidx, double* d_out,
    cudaEvent_t ev0, cudaEvent_t ev1, double* h_gpu_out)
{
    memset(h_ref,0,sz_e*sizeof(double));
    cpu_reference<4>(h_ref,bd.h_vn_ie,bd.inv_dual,bd.h_w,bd.h_cidx,
        bd.h_z_vt_ie,bd.inv_primal,bd.tangent_o,bd.h_z_w_v,bd.h_vidx,
        N_e,N_c,N_v,nlev,nlev);
    CUDA_CHECK(cudaMemcpy(d_vn_ie,bd.h_vn_ie,bd.sz_e*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_w,bd.h_w,bd.sz_c*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_vt_ie,bd.h_z_vt_ie,bd.sz_e*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_w_v,bd.h_z_w_v,bd.sz_v*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cidx,bd.h_cidx,N_e*2*4,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vidx,bd.h_vidx,N_e*2*4,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_dual,bd.inv_dual,N_e*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_primal,bd.inv_primal,N_e*8,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_tangent,bd.tangent_o,N_e*8,cudaMemcpyHostToDevice));
    for(int ci=0;ci<N_GCFG_V7;ci++){
        CUDA_CHECK(cudaMemset(d_out,0,sz_e*8));
        bool launched=true;
        for(int r=0;r<WARMUP;r++){g_flush.flush();
            launched=launch_gpu_v7<4>(ci,d_out,d_vn_ie,d_inv_dual,d_w,d_cidx,
                d_z_vt_ie,d_inv_primal,d_tangent,d_z_w_v,d_vidx,N_e,N_c,N_v,nlev);
            if(!launched)break; CUDA_CHECK(cudaDeviceSynchronize());}
        if(!launched){printf("SKIP: V=7 cfg=%s\n",GCFG_V7[ci].label);continue;}
        CUDA_CHECK(cudaMemcpy(h_gpu_out,d_out,sz_e*8,cudaMemcpyDeviceToHost));
        int nf=0;double mr=0;size_t ff=0;
        bool ok=verify(h_gpu_out,h_ref,sz_e,1e-8,1e-12,&nf,&mr,&ff);
        if(!ok){printf("FAIL: V=7 cfg=%s fails=%d mr=%.3e\n",GCFG_V7[ci].label,nf,mr);continue;}
        else if(ci==0) printf("OK:   V=7 dist=%-12s mr=%.3e\n",dist_label,mr);
        for(int r=0;r<NRUNS;r++){
            g_flush.flush(); CUDA_CHECK(cudaEventRecord(ev0));
            launch_gpu_v7<4>(ci,d_out,d_vn_ie,d_inv_dual,d_w,d_cidx,
                d_z_vt_ie,d_inv_primal,d_tangent,d_z_w_v,d_vidx,N_e,N_c,N_v,nlev);
            CUDA_CHECK(cudaEventRecord(ev1)); CUDA_CHECK(cudaEventSynchronize(ev1));
            float ms=0; CUDA_CHECK(cudaEventElapsedTime(&ms,ev0,ev1));
            fprintf(fcsv,"gpu,7,%d,%d,%d,%d,%s,%s,%d,1,%d,1,%d,%.6f\n",
                nlev,N_e,N_c,N_v,dist_label,GCFG_V7[ci].label,
                GCFG_V7[ci].tx,(int)GCFG_V7[ci].bx,r,(double)ms);
            g_flush.flush();
        }
    }
    printf("Done: nlev=%d dist=%-12s V=7\n",nlev,dist_label);
}

static void run_dist_block_v7(
    FILE* fcsv, int N_e, int N_c, int N_v, int nlev,
    const char* dist_label, int* cell_logical, int* vert_logical,
    double* icon_inv_dual, double* icon_inv_primal, double* icon_tangent)
{
    BenchData bd; bd.alloc(N_e,N_c,N_v,nlev); bd.fill(nlev);
    if(icon_inv_dual) for(int je=0;je<N_e;je++){
        bd.inv_dual[je]=icon_inv_dual[je]; bd.inv_primal[je]=icon_inv_primal[je]; bd.tangent_o[je]=icon_tangent[je];}
    bd.set_variant(4,cell_logical,vert_logical);
    double* h_ref=new double[bd.sz_e]; double* h_gpu_out=new double[bd.sz_e];
    double *d_vn,*d_w,*d_vt,*d_zw,*d_out,*d_id,*d_ip,*d_tg; int *d_ci,*d_vi;
    CUDA_CHECK(cudaMalloc(&d_vn,bd.sz_e*8)); CUDA_CHECK(cudaMalloc(&d_w,bd.sz_c*8));
    CUDA_CHECK(cudaMalloc(&d_vt,bd.sz_e*8)); CUDA_CHECK(cudaMalloc(&d_zw,bd.sz_v*8));
    CUDA_CHECK(cudaMalloc(&d_out,bd.sz_e*8));
    CUDA_CHECK(cudaMalloc(&d_id,N_e*8)); CUDA_CHECK(cudaMalloc(&d_ip,N_e*8)); CUDA_CHECK(cudaMalloc(&d_tg,N_e*8));
    CUDA_CHECK(cudaMalloc(&d_ci,N_e*2*4)); CUDA_CHECK(cudaMalloc(&d_vi,N_e*2*4));
    cudaEvent_t ev0,ev1; CUDA_CHECK(cudaEventCreate(&ev0)); CUDA_CHECK(cudaEventCreate(&ev1));
    run_v7_configs(fcsv,N_e,N_c,N_v,nlev,dist_label,
        h_ref,bd,bd.sz_e,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,d_out,ev0,ev1,h_gpu_out);
    fflush(fcsv);
    cudaFree(d_vn);cudaFree(d_w);cudaFree(d_vt);cudaFree(d_zw);cudaFree(d_out);
    cudaFree(d_id);cudaFree(d_ip);cudaFree(d_tg);cudaFree(d_ci);cudaFree(d_vi);
    CUDA_CHECK(cudaEventDestroy(ev0)); CUDA_CHECK(cudaEventDestroy(ev1));
    delete[]h_ref; delete[]h_gpu_out; bd.free_all();
}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main(int argc, char* argv[]) {
    FILE* fcsv=fopen("z_v_grad_w_gpu.csv","w");
    if(!fcsv){perror("fopen");return 1;}
    fprintf(fcsv,"backend,variant,nlev,N_e,N_c,N_v,cell_dist,"
                 "config_label,TX,TY,BX,BY,run_id,time_ms\n");

    int icon_step=(argc>=2)?atoi(argv[1]):9;
    std::string gp=icon_global_path(icon_step), pp=icon_patch_path(icon_step);
    int icon_nproma=icon_read_nproma(gp.c_str());
    printf("Loading ICON: %s (nproma=%d)\n",pp.c_str(),icon_nproma);

    IconEdgeData ied;
    bool have_exact=(icon_nproma>0)&&icon_load_patch(pp.c_str(),icon_nproma,ied);
    assert(have_exact);
    if(have_exact) icon_print_locality_metrics(ied);

    const int N_e=ied.n_edges_valid, N_c=ied.n_cells, N_v=ied.n_verts;
    printf("Dimensions: N_e=%d N_c=%d N_v=%d\n",N_e,N_c,N_v);

    std::mt19937 rng(42);
    int* cell_logical=new int[N_e*2];
    int* vert_logical=new int[N_e*2];

    cudaDeviceProp prop; CUDA_CHECK(cudaGetDeviceProperties(&prop,0));
    printf("GPU: %s  SM=%d  Configs: V1-V5=%d V6=%d V7=%d\n",prop.name,prop.multiProcessorCount,N_GCFG,N_GCFG_V6,N_GCFG_V7);
    srand((unsigned)time(NULL)); g_flush.init();

    for(int ni=0;ni<N_NLEVS;ni++){
        int nlev_base=NLEVS[ni], nlev_padded=icon_pad_nlev(nlev_base);

        /* EXACT */
        if(have_exact){
            int *ecl=new int[N_e*2], *evl=new int[N_e*2];
            for(int i=0;i<N_e*2;i++){ecl[i]=ied.cell_idx[i]; evl[i]=ied.vert_idx[i];}
            run_dist_block(fcsv,N_e,N_c,N_v,nlev_base,nlev_base,"exact",1,4,ecl,evl,
                ied.inv_dual.data(),ied.inv_primal.data(),ied.tangent_o.data());
            run_dist_block(fcsv,N_e,N_c,N_v,nlev_padded,nlev_base,"exact",5,5,ecl,evl,
                ied.inv_dual.data(),ied.inv_primal.data(),ied.tangent_o.data());
            run_dist_block_v6(fcsv,N_e,N_c,N_v,nlev_padded,nlev_base,(unsigned)nlev_padded,
                "exact",ecl,evl,ied.inv_dual.data(),ied.inv_primal.data(),ied.tangent_o.data());
            /* V7: flat+TX, unpadded — always runs */
            run_dist_block_v7(fcsv,N_e,N_c,N_v,nlev_base,"exact",ecl,evl,
                ied.inv_dual.data(),ied.inv_primal.data(),ied.tangent_o.data());
            delete[]ecl; delete[]evl;
        }

        /* Synthetic */
        for(int di=0;di<2;di++){
            gen_idx_logical(cell_logical,N_e,N_c,(CellDist)di,rng);
            gen_idx_logical(vert_logical,N_e,N_v,(CellDist)di,rng);
            run_dist_block(fcsv,N_e,N_c,N_v,nlev_base,nlev_base,dist_name[di],1,4,
                cell_logical,vert_logical,nullptr,nullptr,nullptr);
            run_dist_block(fcsv,N_e,N_c,N_v,nlev_padded,nlev_base,dist_name[di],5,5,
                cell_logical,vert_logical,nullptr,nullptr,nullptr);
            run_dist_block_v6(fcsv,N_e,N_c,N_v,nlev_padded,nlev_base,(unsigned)nlev_padded,
                dist_name[di],cell_logical,vert_logical,nullptr,nullptr,nullptr);
            /* V7: flat+TX, unpadded — always runs */
            run_dist_block_v7(fcsv,N_e,N_c,N_v,nlev_base,dist_name[di],
                cell_logical,vert_logical,nullptr,nullptr,nullptr);
        }
    }

    g_flush.destroy();
    delete[]cell_logical; delete[]vert_logical;
    if(have_exact) ied.free_all();
    fclose(fcsv);
    printf("\nWritten: z_v_grad_w_gpu.csv\n");
    return 0;
}
