/*
 * bench_gpu.cu -- GPU-only z_v_grad_w stencil benchmark
 *
 * V1/V2 (je,jk) layout: threadIdx.x -> je (stride-1), threadIdx.y -> jk
 * V3/V4 (jk,je) layout: threadIdx.x -> jk (stride-1), threadIdx.y -> je
 *
 * Compile:
 *   nvcc -O3 -arch=sm_80 -std=c++17 -Xcompiler -fopenmp bench_gpu.cu -o bench_gpu
 */

#include "bench_common.h"

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

/* ================================================================ */
/*  CPU reference (serial, for verification)                         */
/* ================================================================ */
template<int V>
static void cpu_reference(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev)
{
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) { STENCIL_BODY(V) }
}

static void cpu_reference_v(int V,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev)
{
    switch (V) {
        case 1: cpu_reference<1>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 2: cpu_reference<2>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 3: cpu_reference<3>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 4: cpu_reference<4>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
    }
}

/* ================================================================ */
/*  numerical verification                                           */
/* ================================================================ */
static bool verify(const double* got, const double* ref, size_t n,
                   double rtol, double atol,
                   int* n_fail, double* max_rel)
{
    *n_fail  = 0;
    *max_rel = 0.0;
    for (size_t i = 0; i < n; i++) {
        double diff = std::abs(got[i] - ref[i]);
        double denom = std::max(std::abs(ref[i]), 1e-300);
        double rel = diff / denom;
        if (rel > *max_rel) *max_rel = rel;
        if (diff > atol + rtol * std::abs(ref[i]))
            (*n_fail)++;
    }
    return *n_fail == 0;
}

/* ================================================================ */
/*  GPU kernel -- V1/V2: threadIdx.x = je, threadIdx.y = jk          */
/*  TX tiles je, TY tiles jk                                         */
/* ================================================================ */

template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_je_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev)
{
    const int je_base = ((int)blockIdx.x * BX + (int)threadIdx.x) * TX;
    const int jk_base = ((int)blockIdx.y * BY + (int)threadIdx.y) * TY;

    int    ci0_a[TX], ci1_a[TX], vi0_a[TX], vi1_a[TX];
    double id_a[TX],  ip_a[TX],  tg_a[TX];

    if (je_base >= N || jk_base >= nlev) return;

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
        if (jk >= nlev) continue;
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
/*  TX tiles jk (stride-1), TY tiles je                              */
/* ================================================================ */

template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_jk_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev)
{
    const int jk_base = ((int)blockIdx.x * BX + (int)threadIdx.x) * TX;
    const int je_base = ((int)blockIdx.y * BY + (int)threadIdx.y) * TY;

    int    ci0_a[TY], ci1_a[TY], vi0_a[TY], vi1_a[TY];
    double id_a[TY],  ip_a[TY],  tg_a[TY];

    if (je_base >= N || jk_base >= nlev) return;

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
            if (jk >= nlev) continue;
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
/*  GPU config table                                                 */
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
/*  GPU launch dispatch                                              */
/* ================================================================ */

template<int V>
static void launch_gpu(int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev)
{
    #define LG_JE(TX_,TY_,BX_,BY_) do {                                   \
        dim3 blk(BX_, BY_);                                                \
        dim3 grd(((unsigned)N    + (BX_)*(TX_) - 1) / ((BX_)*(TX_)),      \
                 ((unsigned)nlev + (BY_)*(TY_) - 1) / ((BY_)*(TY_)));     \
        gpu_kernel_je_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>(              \
            out, vn_ie, inv_dual, w, cell_idx,                             \
            z_vt_ie, inv_primal, tangent, z_w_v, vert_idx, N, nlev);       \
    } while(0)

    #define LG_JK(TX_,TY_,BX_,BY_) do {                                   \
        dim3 blk(BX_, BY_);                                                \
        dim3 grd(((unsigned)nlev + (BX_)*(TX_) - 1) / ((BX_)*(TX_)),      \
                 ((unsigned)N    + (BY_)*(TY_) - 1) / ((BY_)*(TY_)));     \
        gpu_kernel_jk_first<TX_,TY_,BX_,BY_,V><<<grd,blk>>>(              \
            out, vn_ie, inv_dual, w, cell_idx,                             \
            z_vt_ie, inv_primal, tangent, z_w_v, vert_idx, N, nlev);       \
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
}

static void launch_gpu_v(int V, int cfg,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev)
{
    switch (V) {
        case 1: launch_gpu<1>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 2: launch_gpu<2>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 3: launch_gpu<3>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 4: launch_gpu<4>(cfg,out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
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

static constexpr int FLUSH_N     = 8192;
static constexpr int FLUSH_STEPS = 3;

struct GpuFlush {
    double *d_A, *d_B;
    bool inited;

    GpuFlush() : d_A(nullptr), d_B(nullptr), inited(false) {}

    void init() {
        if (inited) return;
        size_t n = (size_t)FLUSH_N * FLUSH_N;
        size_t bytes = n * sizeof(double);

        /* host init */
        double* h = new double[n];
        for (size_t i = 0; i < n; i++) {
            uint64_t v = splitmix64(12345ULL + i);
            h[i] = (double)(v >> 11) / (double)(1ULL << 53);
        }

        CUDA_CHECK(cudaMalloc(&d_A, bytes));
        CUDA_CHECK(cudaMalloc(&d_B, bytes));
        CUDA_CHECK(cudaMemcpy(d_A, h, bytes, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_B, h, bytes, cudaMemcpyHostToDevice));
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
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    void destroy() {
        if (d_A) cudaFree(d_A);
        if (d_B) cudaFree(d_B);
        d_A = d_B = nullptr;
        inited = false;
    }
};

static GpuFlush g_flush;

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main() {
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

    /* print GPU info */
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    printf("GPU: %s  SM count: %d\n", prop.name, prop.multiProcessorCount);
    printf("Configs: %d\n", N_GCFG);

    /* init flush buffers once */
    g_flush.init();
    printf("GPU flush buffers initialized (2 x %.0f MB)\n",
           (double)FLUSH_N * FLUSH_N * 8 / 1e6);

    for (int nlev_i = 0; nlev_i < N_NLEVS; nlev_i++) {
        int nlev = NLEVS[nlev_i];

        BenchData bd;
        bd.alloc(N, nlev);
        bd.fill(nlev);

        /* host reference buffer (computed once per nlev,dist,V) */
        double* h_ref = new double[bd.sz2d];
        /* host buffer to readback GPU output for verification */
        double* h_gpu_out = new double[bd.sz2d];

        /* device arrays */
        double *d_vn_ie, *d_w, *d_z_vt_ie, *d_z_w_v, *d_out;
        double *d_inv_dual, *d_inv_primal, *d_tangent;
        int    *d_cidx, *d_vidx;

        CUDA_CHECK(cudaMalloc(&d_vn_ie,      bd.sz2d * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_w,          bd.sz2d * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_z_vt_ie,    bd.sz2d * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_z_w_v,      bd.sz2d * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_out,        bd.sz2d * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_inv_dual,   N * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_inv_primal, N * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_tangent,    N * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_cidx,       N * 2 * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_vidx,       N * 2 * sizeof(int)));

        CUDA_CHECK(cudaMemcpy(d_inv_dual,   bd.inv_dual,
                              N*sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_inv_primal, bd.inv_primal,
                              N*sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_tangent,    bd.tangent_o,
                              N*sizeof(double), cudaMemcpyHostToDevice));

        cudaEvent_t ev0, ev1;
        CUDA_CHECK(cudaEventCreate(&ev0));
        CUDA_CHECK(cudaEventCreate(&ev1));

        for (int di = 0; di < 4; di++) {
            CellDist dist = (CellDist)di;
            gen_cell_idx_logical(cell_logical, N, dist, rng);

            for (int V = 1; V <= 4; V++) {
                bd.set_variant(V, cell_logical, vd.logical);

                /* ---- compute CPU reference ONCE per (nlev, dist, V) ---- */
                cpu_reference_v(V,
                    h_ref, bd.h_vn_ie, bd.inv_dual,
                    bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
                    bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);

                /* upload variant data to device */
                CUDA_CHECK(cudaMemcpy(d_vn_ie,   bd.h_vn_ie,
                    bd.sz2d*sizeof(double), cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_w,       bd.h_w,
                    bd.sz2d*sizeof(double), cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_z_vt_ie, bd.h_z_vt_ie,
                    bd.sz2d*sizeof(double), cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_z_w_v,   bd.h_z_w_v,
                    bd.sz2d*sizeof(double), cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_cidx,    bd.h_cidx,
                    N*2*sizeof(int), cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_vidx,    bd.h_vidx,
                    N*2*sizeof(int), cudaMemcpyHostToDevice));

                for (int ci = 0; ci < N_GCFG; ci++) {
                    /* warmup */
                    for (int r = 0; r < WARMUP; r++) {
                        g_flush.flush();
                        launch_gpu_v(V, ci, d_out, d_vn_ie, d_inv_dual,
                            d_w, d_cidx, d_z_vt_ie, d_inv_primal,
                            d_tangent, d_z_w_v, d_vidx, N, nlev);
                        CUDA_CHECK(cudaDeviceSynchronize());
                    }

                    /* ---- verify first config run against CPU ref ---- */
                    CUDA_CHECK(cudaMemcpy(h_gpu_out, d_out,
                        bd.sz2d * sizeof(double), cudaMemcpyDeviceToHost));
                    int n_fail = 0;
                    double max_rel = 0.0;
                    bool ok = verify(h_gpu_out, h_ref, bd.sz2d,
                                     1e-12, 1e-15, &n_fail, &max_rel);
                    if (!ok) {
                        printf("VERIFY FAIL: nlev=%d dist=%-12s V=%d "
                               "cfg=%-14s  fails=%d max_rel=%.3e\n",
                               nlev, dist_name[di], V,
                               GCFG[ci].label, n_fail, max_rel);
                    } else if (ci == 0) {
                        /* print once per (nlev, dist, V) to confirm */
                        printf("VERIFY OK:   nlev=%d dist=%-12s V=%d "
                               "max_rel=%.3e\n",
                               nlev, dist_name[di], V, max_rel);
                    }

                    /* timed runs */
                    for (int r = 0; r < NRUNS; r++) {
                        g_flush.flush();
                        CUDA_CHECK(cudaEventRecord(ev0));
                        launch_gpu_v(V, ci, d_out, d_vn_ie, d_inv_dual,
                            d_w, d_cidx, d_z_vt_ie, d_inv_primal,
                            d_tangent, d_z_w_v, d_vidx, N, nlev);
                        CUDA_CHECK(cudaEventRecord(ev1));
                        CUDA_CHECK(cudaEventSynchronize(ev1));
                        float ms = 0.0f;
                        CUDA_CHECK(cudaEventElapsedTime(&ms, ev0, ev1));
                        fprintf(fcsv,
                            "gpu,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%.6f\n",
                            V, nlev, N, dist_name[di],
                            GCFG[ci].label,
                            GCFG[ci].tx, GCFG[ci].ty,
                            GCFG[ci].bx, GCFG[ci].by,
                            r, (double)ms);
                    }
                }

                printf("Done: nlev=%d  dist=%-12s  V=%d\n",
                       nlev, dist_name[di], V);
                fflush(fcsv);
            }
        }

        cudaFree(d_vn_ie);      cudaFree(d_w);
        cudaFree(d_z_vt_ie);    cudaFree(d_z_w_v);    cudaFree(d_out);
        cudaFree(d_inv_dual);   cudaFree(d_inv_primal);
        cudaFree(d_tangent);    cudaFree(d_cidx);      cudaFree(d_vidx);
        CUDA_CHECK(cudaEventDestroy(ev0));
        CUDA_CHECK(cudaEventDestroy(ev1));
        delete[] h_ref;
        delete[] h_gpu_out;
        bd.free_all();
    }

    g_flush.destroy();
    vd.free_all();
    delete[] cell_logical;
    fclose(fcsv);

    printf("\nResults written to z_v_grad_w_gpu.csv\n");
    printf("Total rows: 2 nlevs x 4 dists x 4 variants x %d cfgs x %d = %d\n",
           N_GCFG, NRUNS, 2 * 4 * 4 * N_GCFG * NRUNS);
    return 0;
}
