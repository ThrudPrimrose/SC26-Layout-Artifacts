/*
 * bench_gpu_oldstyle.cu -- OLD-STYLE GPU benchmark, adapted for new headers
 *
 * Flat kernels (no grid-stride), cross-mapped jk_first.
 * Minimal changes from the original fast version: only API compatibility
 * fixes for new bench_common.h / icon_data_loader.h (N_e/N_c/N_v).
 *
 * Compile (NVIDIA): nvcc -O3 -arch=sm_90 -std=c++17 -Xcompiler -fopenmp bench_gpu_oldstyle.cu -o bench_gpu_old
 * Compile (AMD):    hipcc -O3 -std=c++17 -fopenmp bench_gpu_oldstyle.cu -o bench_gpu_old
 */
#include "bench_common.h"
#include "icon_data_loader.h"
#include <ctime>
#include <cassert>

#if __HIP_PLATFORM_AMD__
#include "hip/hip_runtime.h"
#endif

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
        fprintf(stderr, "CUDA launch warn %s:%d: %s\n",                    \
                __FILE__, __LINE__, cudaGetErrorString(e));                 \
        return false;                                                      \
    }                                                                      \
} while(0)

/* ================================================================ */
/*  CPU reference                                                    */
/* ================================================================ */
template<int V>
static void cpu_reference(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev)
{
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N_e; je++) { STENCIL_BODY(V) }
}

static void cpu_reference_v(int V,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N_e, int N_c, int N_v, int nlev)
{
    switch (V) {
        case 1: cpu_reference<1>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev); break;
        case 2: cpu_reference<2>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev); break;
        case 3: cpu_reference<3>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev); break;
        case 4: cpu_reference<4>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev); break;
    }
}

/* ================================================================ */
/*  Verify                                                           */
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
/*  GPU kernel -- V1/V2: je-first (UNCHANGED from original)          */
/*  threadIdx.x = je (stride-1), threadIdx.y = jk                   */
/* ================================================================ */
template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_je_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev)
{
    const int je_base = ((int)blockIdx.x * BX + (int)threadIdx.x) * TX;
    const int jk_base = ((int)blockIdx.y * BY + (int)threadIdx.y) * TY;

    int    ci0_a[TX], ci1_a[TX], vi0_a[TX], vi1_a[TX];
    double id_a[TX],  ip_a[TX],  tg_a[TX];

    #pragma unroll
    for (int tx = 0; tx < TX; tx++) {
        int je = je_base + tx;
        if (je < N_e) {
            ci0_a[tx] = cell_idx[IN<V>(je, 0, N_e)];
            ci1_a[tx] = cell_idx[IN<V>(je, 1, N_e)];
            vi0_a[tx] = vert_idx[IN<V>(je, 0, N_e)];
            vi1_a[tx] = vert_idx[IN<V>(je, 1, N_e)];
            id_a[tx]  = inv_dual[je];
            ip_a[tx]  = inv_primal[je];
            tg_a[tx]  = tangent[je];
        }
    }

    #pragma unroll
    for (int ty = 0; ty < TY; ty++) {
        int jk = jk_base + ty;
        if (jk >= nlev) continue;
        #pragma unroll
        for (int tx = 0; tx < TX; tx++) {
            int je = je_base + tx;
            if (je >= N_e) continue;
            int c2d = IC<V>(je, jk, N_e, nlev);
            out[c2d] =
                vn_ie[c2d] * id_a[tx] *
                    (w[IC<V>(ci0_a[tx], jk, N_c, nlev)] -
                     w[IC<V>(ci1_a[tx], jk, N_c, nlev)])
              + z_vt_ie[c2d] * ip_a[tx] * tg_a[tx] *
                    (z_w_v[IC<V>(vi0_a[tx], jk, N_v, nlev)] -
                     z_w_v[IC<V>(vi1_a[tx], jk, N_v, nlev)]);
        }
    }
}

/* ================================================================ */
/*  GPU kernel -- V3/V4: jk-first, cross-mapped (UNCHANGED)          */
/*  threadIdx.x = jk (stride-1), threadIdx.y = je                   */
/*  blockIdx.x = je (large), blockIdx.y = jk (small)                */
/* ================================================================ */
template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_jk_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev)
{
    const int jk_base = ((int)blockIdx.y * BX + (int)threadIdx.x) * TX;
    const int je_base = ((int)blockIdx.x * BY + (int)threadIdx.y) * TY;

    int    ci0_a[TY], ci1_a[TY], vi0_a[TY], vi1_a[TY];
    double id_a[TY],  ip_a[TY],  tg_a[TY];

    #pragma unroll
    for (int ty = 0; ty < TY; ty++) {
        int je = je_base + ty;
        if (je < N_e) {
            ci0_a[ty] = cell_idx[IN<V>(je, 0, N_e)];
            ci1_a[ty] = cell_idx[IN<V>(je, 1, N_e)];
            vi0_a[ty] = vert_idx[IN<V>(je, 0, N_e)];
            vi1_a[ty] = vert_idx[IN<V>(je, 1, N_e)];
            id_a[ty]  = inv_dual[je];
            ip_a[ty]  = inv_primal[je];
            tg_a[ty]  = tangent[je];
        }
    }

    #pragma unroll
    for (int ty = 0; ty < TY; ty++) {
        int je = je_base + ty;
        if (je >= N_e) continue;
        #pragma unroll
        for (int tx = 0; tx < TX; tx++) {
            int jk = jk_base + tx;
            if (jk >= nlev) continue;
            int c2d = IC<V>(je, jk, N_e, nlev);
            out[c2d] =
                vn_ie[c2d] * id_a[ty] *
                    (w[IC<V>(ci0_a[ty], jk, N_c, nlev)] -
                     w[IC<V>(ci1_a[ty], jk, N_c, nlev)])
              + z_vt_ie[c2d] * ip_a[ty] * tg_a[ty] *
                    (z_w_v[IC<V>(vi0_a[ty], jk, N_v, nlev)] -
                     z_w_v[IC<V>(vi1_a[ty], jk, N_v, nlev)]);
        }
    }
}

/* ================================================================ */
/*  Config table (UNCHANGED from original)                           */
/* ================================================================ */
struct GpuCfg { int tx, ty, bx, by; const char* label; };
static constexpr GpuCfg GCFG[] = {
    {1, 1, 256, 1, "1x1_256x1"},
    {1, 1, 128, 1, "1x1_128x1"},
    {2, 1, 128, 1, "2x1_128x1"},
    {4, 1, 128, 1, "4x1_128x1"},
    {4, 1, 64,  1, "4x1_64x1"},
    {1, 2, 256, 1, "1x2_256x1"},
    {1, 4, 256, 1, "1x4_256x1"},
    {2, 2, 128, 1, "2x2_128x1"},
    {2, 4, 128, 1, "2x4_128x1"},
    {4, 2, 64,  1, "4x2_64x1"},
    {1, 1, 32, 16, "1x1_32x16"},
    {1, 1, 32,  8, "1x1_32x8"},
    {1, 1, 32,  4, "1x1_32x4"},
    {2, 1, 32, 16, "2x1_32x16"},
    {2, 1, 32,  8, "2x1_32x8"},
    {2, 1, 32,  4, "2x1_32x4"},
    {1, 2, 32, 16, "1x2_32x16"},
    {1, 2, 32,  8, "1x2_32x8"},
    {1, 2, 32,  4, "1x2_32x4"},
    {2, 2, 32,  8, "2x2_32x8"},
    {2, 2, 32,  4, "2x2_32x4"},
    {4, 2, 32,  4, "4x2_32x4"},
};
static constexpr int N_GCFG = sizeof(GCFG) / sizeof(GCFG[0]);

/* ================================================================ */
/*  Launch dispatch (UNCHANGED structure, added N_c/N_v + skip>65k)  */
/* ================================================================ */
template<int V>
static bool launch_gpu(int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N_e, int N_c, int N_v, int nlev)
{
    #define LG_JE(TX_,TY_,BX_,BY_) do {                                   \
        unsigned gx = ((unsigned)N_e  + (BX_)*(TX_)-1) / ((BX_)*(TX_));   \
        unsigned gy = ((unsigned)nlev + (BY_)*(TY_)-1) / ((BY_)*(TY_));   \
        if (gx>65535u||gy>65535u) return false;                            \
        dim3 blk(BX_, BY_); dim3 grd(gx, gy);                             \
        gpu_kernel_je_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>(              \
            out,vn_ie,inv_dual,w,cell_idx,                                 \
            z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev);  \
        CUDA_LAUNCH_CHECK(); return true;                                  \
    } while(0)

    #define LG_JK(TX_,TY_,BX_,BY_) do {                                   \
        unsigned gx = ((unsigned)N_e  + (BY_)*(TY_)-1) / ((BY_)*(TY_));   \
        unsigned gy = ((unsigned)nlev + (BX_)*(TX_)-1) / ((BX_)*(TX_));   \
        if (gx>65535u||gy>65535u) return false;                            \
        dim3 blk(BX_, BY_); dim3 grd(gx, gy);                             \
        gpu_kernel_jk_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>(              \
            out,vn_ie,inv_dual,w,cell_idx,                                 \
            z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev);  \
        CUDA_LAUNCH_CHECK(); return true;                                  \
    } while(0)

    #define LG(TX_,TY_,BX_,BY_) do {                                       \
        if constexpr (V <= 2) { LG_JE(TX_,TY_,BX_,BY_); }                 \
        else                  { LG_JK(TX_,TY_,BX_,BY_); }                  \
    } while(0)

    switch (cfg) {
        case  0: LG(1,1,256,1);  break;
        case  1: LG(1,1,128,1);  break;
        case  2: LG(2,1,128,1);  break;
        case  3: LG(4,1,128,1);  break;
        case  4: LG(4,1,64,1);   break;
        case  5: LG(1,2,256,1);  break;
        case  6: LG(1,4,256,1);  break;
        case  7: LG(2,2,128,1);  break;
        case  8: LG(2,4,128,1);  break;
        case  9: LG(4,2,64,1);   break;
        case 10: LG(1,1,32,16);  break;
        case 11: LG(1,1,32,8);   break;
        case 12: LG(1,1,32,4);   break;
        case 13: LG(2,1,32,16);  break;
        case 14: LG(2,1,32,8);   break;
        case 15: LG(2,1,32,4);   break;
        case 16: LG(1,2,32,16);  break;
        case 17: LG(1,2,32,8);   break;
        case 18: LG(1,2,32,4);   break;
        case 19: LG(2,2,32,8);   break;
        case 20: LG(2,2,32,4);   break;
        case 21: LG(4,2,32,4);   break;
    }
    #undef LG
    #undef LG_JE
    #undef LG_JK
    return false;
}

static bool launch_gpu_v(int V, int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N_e, int N_c, int N_v, int nlev)
{
    switch (V) {
        case 1: return launch_gpu<1>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev);
        case 2: return launch_gpu<2>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev);
        case 3: return launch_gpu<3>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev);
        case 4: return launch_gpu<4>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N_e,N_c,N_v,nlev);
    }
    return false;
}

/* ================================================================ */
/*  GPU cache flush (UNCHANGED)                                      */
/* ================================================================ */
__global__ void flush_stencil_step(const double* __restrict__ A,
                                   double* __restrict__ B, int N) {
    int i = blockIdx.y*blockDim.y+threadIdx.y;
    int j = blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=1&&i<N-1&&j>=1&&j<N-1)
        B[i*N+j]=0.25*(A[(i-1)*N+j]+A[(i+1)*N+j]+A[i*N+(j-1)]+A[i*N+(j+1)]);
}
static constexpr int FLUSH_N=8192*4, FLUSH_STEPS=3;
struct GpuFlush {
    double *d_A=nullptr,*d_B=nullptr; bool inited=false;
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
/*  run_variant_configs (adapted for N_e/N_c/N_v)                    */
/* ================================================================ */
static void run_variant_configs(
    FILE* fcsv, int V, int N_e, int N_c, int N_v, int nlev,
    const char* dist_label,
    double* h_ref, BenchData& bd, size_t sz_e,
    double* d_vn_ie, double* d_inv_dual, double* d_w, int* d_cidx,
    double* d_z_vt_ie, double* d_inv_primal, double* d_tangent,
    double* d_z_w_v, int* d_vidx, double* d_out,
    cudaEvent_t ev0, cudaEvent_t ev1, double* h_gpu_out)
{
    cpu_reference_v(V, h_ref, bd.h_vn_ie, bd.inv_dual,
        bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
        bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N_e, N_c, N_v, nlev);

    CUDA_CHECK(cudaMemcpy(d_vn_ie,   bd.h_vn_ie,   bd.sz_e*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_w,       bd.h_w,       bd.sz_c*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_vt_ie, bd.h_z_vt_ie, bd.sz_e*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z_w_v,   bd.h_z_w_v,   bd.sz_v*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cidx,    bd.h_cidx,    N_e*2*4,   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vidx,    bd.h_vidx,    N_e*2*4,   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_dual,   bd.inv_dual,   N_e*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_inv_primal, bd.inv_primal, N_e*8, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_tangent,    bd.tangent_o,  N_e*8, cudaMemcpyHostToDevice));

    for (int ci = 0; ci < N_GCFG; ci++) {
        CUDA_CHECK(cudaMemset(d_out, 0, sz_e*8));
        bool launched = true;
        for (int r = 0; r < WARMUP; r++) {
            g_flush.flush();
            launched = launch_gpu_v(V, ci, d_out, d_vn_ie, d_inv_dual,
                d_w, d_cidx, d_z_vt_ie, d_inv_primal,
                d_tangent, d_z_w_v, d_vidx, N_e, N_c, N_v, nlev);
            if (!launched) break;
            CUDA_CHECK(cudaDeviceSynchronize());
        }
        if (!launched) {
            printf("SKIP: nlev=%d dist=%-12s V=%d cfg=%-14s\n",
                   nlev, dist_label, V, GCFG[ci].label);
            continue;
        }

        CUDA_CHECK(cudaMemcpy(h_gpu_out, d_out, sz_e*8, cudaMemcpyDeviceToHost));
        int nf=0; double mr=0; size_t ff=0;
        bool ok = verify(h_gpu_out, h_ref, sz_e, 1e-8, 1e-12, &nf, &mr, &ff);
        if (!ok) {
            printf("FAIL: V=%d cfg=%-14s fails=%d max_rel=%.3e\n",
                   V, GCFG[ci].label, nf, mr);
            continue;
        } else if (ci == 0) {
            printf("OK:   nlev=%d dist=%-12s V=%d mr=%.3e\n",
                   nlev, dist_label, V, mr);
        }

        for (int r = 0; r < NRUNS; r++) {
            g_flush.flush();
            CUDA_CHECK(cudaEventRecord(ev0));
            launch_gpu_v(V, ci, d_out, d_vn_ie, d_inv_dual,
                d_w, d_cidx, d_z_vt_ie, d_inv_primal,
                d_tangent, d_z_w_v, d_vidx, N_e, N_c, N_v, nlev);
            CUDA_CHECK(cudaEventRecord(ev1));
            CUDA_CHECK(cudaEventSynchronize(ev1));
            float ms = 0.0f;
            CUDA_CHECK(cudaEventElapsedTime(&ms, ev0, ev1));
            fprintf(fcsv, "gpu,%d,%d,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%.6f\n",
                V, nlev, N_e, N_c, N_v, dist_label,
                GCFG[ci].label,
                GCFG[ci].tx, GCFG[ci].ty,
                GCFG[ci].bx, GCFG[ci].by,
                r, (double)ms);
            g_flush.flush();
        }
    }
    printf("Done: nlev=%d dist=%-12s V=%d\n", nlev, dist_label, V);
}

/* ================================================================ */
/*  main — uses new ICON loader + gen_idx_logical                    */
/* ================================================================ */
int main(int argc, char* argv[]) {
    FILE* fcsv = fopen("z_v_grad_w_gpu_old.csv", "w");
    if (!fcsv) { perror("fopen"); return 1; }
    fprintf(fcsv, "backend,variant,nlev,N_e,N_c,N_v,cell_dist,"
                  "config_label,TX,TY,BX,BY,run_id,time_ms\n");

    int icon_step = (argc >= 2) ? atoi(argv[1]) : 9;
    std::string gp = icon_global_path(icon_step), pp = icon_patch_path(icon_step);
    int icon_nproma = icon_read_nproma(gp.c_str());
    printf("Loading ICON: %s (nproma=%d)\n", pp.c_str(), icon_nproma);

    IconEdgeData ied;
    bool have_exact = (icon_nproma > 0) && icon_load_patch(pp.c_str(), icon_nproma, ied);
    assert(have_exact);
    if (have_exact) icon_print_locality_metrics(ied);

    const int N_e = ied.n_edges_valid, N_c = ied.n_cells, N_v = ied.n_verts;
    printf("Dimensions: N_e=%d N_c=%d N_v=%d\n", N_e, N_c, N_v);

    std::mt19937 rng(42);
    int* cell_logical = new int[N_e * 2];
    int* vert_logical = new int[N_e * 2];

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    printf("GPU: %s  SM=%d  Configs: %d (old-style)\n",
           prop.name, prop.multiProcessorCount, N_GCFG);

    srand((unsigned)time(NULL));
    g_flush.init();

    for (int ni = 0; ni < N_NLEVS; ni++) {
        int nlev = NLEVS[ni];

        /* ---- EXACT ---- */
        if (have_exact) {
            int* ecl = new int[N_e * 2];
            int* evl = new int[N_e * 2];
            for (int i = 0; i < N_e * 2; i++) {
                ecl[i] = ied.cell_idx[i];
                evl[i] = ied.vert_idx[i];
            }

            BenchData bd;
            bd.alloc(N_e, N_c, N_v, nlev);
            bd.fill(nlev);
            for (int je = 0; je < N_e; je++) {
                bd.inv_dual[je]   = ied.inv_dual[je];
                bd.inv_primal[je] = ied.inv_primal[je];
                bd.tangent_o[je]  = ied.tangent_o[je];
            }

            double* h_ref     = new double[bd.sz_e];
            double* h_gpu_out = new double[bd.sz_e];
            double *d_vn,*d_w,*d_vt,*d_zw,*d_out,*d_id,*d_ip,*d_tg;
            int *d_ci,*d_vi;
            CUDA_CHECK(cudaMalloc(&d_vn, bd.sz_e*8));
            CUDA_CHECK(cudaMalloc(&d_w,  bd.sz_c*8));
            CUDA_CHECK(cudaMalloc(&d_vt, bd.sz_e*8));
            CUDA_CHECK(cudaMalloc(&d_zw, bd.sz_v*8));
            CUDA_CHECK(cudaMalloc(&d_out,bd.sz_e*8));
            CUDA_CHECK(cudaMalloc(&d_id, N_e*8));
            CUDA_CHECK(cudaMalloc(&d_ip, N_e*8));
            CUDA_CHECK(cudaMalloc(&d_tg, N_e*8));
            CUDA_CHECK(cudaMalloc(&d_ci, N_e*2*4));
            CUDA_CHECK(cudaMalloc(&d_vi, N_e*2*4));
            cudaEvent_t ev0, ev1;
            CUDA_CHECK(cudaEventCreate(&ev0));
            CUDA_CHECK(cudaEventCreate(&ev1));

            for (int V = 1; V <= 4; V++) {
                bd.set_variant(V, ecl, evl);
                run_variant_configs(fcsv, V, N_e, N_c, N_v, nlev, "exact",
                    h_ref, bd, bd.sz_e,
                    d_vn, d_id, d_w, d_ci, d_vt, d_ip, d_tg, d_zw, d_vi, d_out,
                    ev0, ev1, h_gpu_out);
                fflush(fcsv);
            }

            cudaFree(d_vn); cudaFree(d_w); cudaFree(d_vt);
            cudaFree(d_zw); cudaFree(d_out);
            cudaFree(d_id); cudaFree(d_ip); cudaFree(d_tg);
            cudaFree(d_ci); cudaFree(d_vi);
            CUDA_CHECK(cudaEventDestroy(ev0));
            CUDA_CHECK(cudaEventDestroy(ev1));
            delete[] h_ref; delete[] h_gpu_out;
            delete[] ecl; delete[] evl;
            bd.free_all();
        }

        /* ---- Synthetic ---- */
        for (int di = 0; di < 2; di++) {
            gen_idx_logical(cell_logical, N_e, N_c, (CellDist)di, rng);
            gen_idx_logical(vert_logical, N_e, N_v, (CellDist)di, rng);

            BenchData bd;
            bd.alloc(N_e, N_c, N_v, nlev);
            bd.fill(nlev);

            double* h_ref     = new double[bd.sz_e];
            double* h_gpu_out = new double[bd.sz_e];
            double *d_vn,*d_w,*d_vt,*d_zw,*d_out,*d_id,*d_ip,*d_tg;
            int *d_ci,*d_vi;
            CUDA_CHECK(cudaMalloc(&d_vn, bd.sz_e*8));
            CUDA_CHECK(cudaMalloc(&d_w,  bd.sz_c*8));
            CUDA_CHECK(cudaMalloc(&d_vt, bd.sz_e*8));
            CUDA_CHECK(cudaMalloc(&d_zw, bd.sz_v*8));
            CUDA_CHECK(cudaMalloc(&d_out,bd.sz_e*8));
            CUDA_CHECK(cudaMalloc(&d_id, N_e*8));
            CUDA_CHECK(cudaMalloc(&d_ip, N_e*8));
            CUDA_CHECK(cudaMalloc(&d_tg, N_e*8));
            CUDA_CHECK(cudaMalloc(&d_ci, N_e*2*4));
            CUDA_CHECK(cudaMalloc(&d_vi, N_e*2*4));
            cudaEvent_t ev0, ev1;
            CUDA_CHECK(cudaEventCreate(&ev0));
            CUDA_CHECK(cudaEventCreate(&ev1));

            for (int V = 1; V <= 4; V++) {
                bd.set_variant(V, cell_logical, vert_logical);
                run_variant_configs(fcsv, V, N_e, N_c, N_v, nlev, dist_name[di],
                    h_ref, bd, bd.sz_e,
                    d_vn, d_id, d_w, d_ci, d_vt, d_ip, d_tg, d_zw, d_vi, d_out,
                    ev0, ev1, h_gpu_out);
                fflush(fcsv);
            }

            cudaFree(d_vn); cudaFree(d_w); cudaFree(d_vt);
            cudaFree(d_zw); cudaFree(d_out);
            cudaFree(d_id); cudaFree(d_ip); cudaFree(d_tg);
            cudaFree(d_ci); cudaFree(d_vi);
            CUDA_CHECK(cudaEventDestroy(ev0));
            CUDA_CHECK(cudaEventDestroy(ev1));
            delete[] h_ref; delete[] h_gpu_out;
            bd.free_all();
        }
    }

    g_flush.destroy();
    delete[] cell_logical;
    delete[] vert_logical;
    if (have_exact) ied.free_all();
    fclose(fcsv);
    printf("\nWritten: z_v_grad_w_gpu_old.csv\n");
    return 0;
}
