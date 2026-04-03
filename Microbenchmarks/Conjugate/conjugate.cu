#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

/*  GPU out-of-place: negate all K fields.
 *  TOTAL = K * N_base = const.   BW = 2·K·N_base·8.
 *  VL sweep: 2 … 64.                                                 */
constexpr int64_t TOTAL = 1LL << 27;          // ~1 GB per side on GPU
constexpr int64_t BLK   = 256;
constexpr int64_t RUNS  = 100;
constexpr int64_t GPU_MAX_VL = 64;

/* ═══ Layout structs ═══ */

template<int K>
struct AoS { double f[K]; };

template<int K, int VL>
struct AoSoA { double f[K][VL]; };

/* ═══ GPU kernels ═══ */

template<int K>
__global__ void g_aos(const AoS<K> *__restrict__ in,
                      AoS<K> *__restrict__ out, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n)
        for (int k = 0; k < K; k++)
            out[i].f[k] = -in[i].f[k];
}

/* SoA: K separate arrays packed into contiguous in/out buffers.
 * in_ptrs[k]  = base_in  + k * n_base
 * out_ptrs[k] = base_out + k * n_base                                */
template<int K>
__global__ void g_soa(const double *__restrict__ base_in,
                      double *__restrict__ base_out,
                      int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n)
        for (int k = 0; k < K; k++)
            base_out[(int64_t)k * n + i] = -base_in[(int64_t)k * n + i];
}

template<int K, int VL>
__global__ void g_aosoa(const AoSoA<K,VL> *__restrict__ in,
                        AoSoA<K,VL> *__restrict__ out, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n) {
        int blk  = i / VL;
        int lane = i % VL;
        for (int k = 0; k < K; k++)
            out[blk].f[k][lane] = -in[blk].f[k][lane];
    }
}

/* ═══ Reporting ═══ */

static FILE *csv;

static double bw_oop(int K, int64_t n, double ms) {
    return 2.0 * K * n * sizeof(double) / (ms * 1e6);
}

#define GPU_BENCH(K_val, n_base, label, call) do { \
    for (int w = 0; w < 5; w++) { call; } \
    cudaDeviceSynchronize(); \
    for (int r = 0; r < RUNS; r++) { \
        cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b); \
        cudaEventRecord(a); \
        call; \
        cudaEventRecord(b); cudaEventSynchronize(b); \
        float ms; cudaEventElapsedTime(&ms, a, b); \
        fprintf(csv, "%d,%s,%d,%.6f,%.2f\n", \
                K_val, label, r, (double)ms, bw_oop(K_val, n_base, ms)); \
        cudaEventDestroy(a); cudaEventDestroy(b); \
    } \
    { /* print avg */ \
        float tot = 0; \
        for (int r = 0; r < 1; r++) { \
            cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b); \
            cudaEventRecord(a); call; \
            cudaEventRecord(b); cudaEventSynchronize(b); \
            float ms; cudaEventElapsedTime(&ms, a, b); tot += ms; \
            cudaEventDestroy(a); cudaEventDestroy(b); \
        } \
        printf("  %-14s  (see csv)\n", label); \
    } \
} while (0)

/* ═══ Per-K bench functions ═══ */

template<int K>
static void bench_aos(int64_t n, double *di, double *dout) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(K, n, "AoS",
              (g_aos<K><<<grid, BLK>>>((AoS<K>*)di, (AoS<K>*)dout, n)));
}

template<int K>
static void bench_soa(int64_t n, double *di, double *dout) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(K, n, "SoA",
              (g_soa<K><<<grid, BLK>>>(di, dout, n)));
}

template<int K, int VL>
static void bench_aosoa(int64_t n, double *di, double *dout, const char *label) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(K, n, label,
              (g_aosoa<K,VL><<<grid, BLK>>>(
                  (AoSoA<K,VL>*)di, (AoSoA<K,VL>*)dout, n)));
}

template<int K>
static void run_all(double *di, double *dout) {
    int64_t n = (TOTAL / K / GPU_MAX_VL) * GPU_MAX_VL;
    printf("\n── K=%d  N_base=%lld  total=%.1f GB (per side) ──\n",
           K, (long long)n, (double)K * n * 8 / 1e9);

    bench_aos<K>(n, di, dout);
    bench_soa<K>(n, di, dout);
    bench_aosoa<K,  2>(n, di, dout, "AoSoA-2");
    bench_aosoa<K,  4>(n, di, dout, "AoSoA-4");
    bench_aosoa<K,  8>(n, di, dout, "AoSoA-8");
    bench_aosoa<K, 16>(n, di, dout, "AoSoA-16");
    bench_aosoa<K, 32>(n, di, dout, "AoSoA-32");
    bench_aosoa<K, 64>(n, di, dout, "AoSoA-64");
}

int main() {
    csv = fopen("results_gpu_oop.csv", "w");
    fprintf(csv, "K,layout,run,ms,gbps\n");

    printf("conj OOP (GPU): TOTAL=%lldM  runs=%d  dtype=double\n",
           (long long)(TOTAL >> 20), (int)RUNS);

    /* Allocate max needed: K=3 has largest N_base, needs K*N_base*8 per side.
     * But K*N_base = TOTAL (const), so each side = TOTAL*8 bytes.           */
    size_t bytes = TOTAL * sizeof(double);
    double *di, *dout;
    if (cudaMalloc(&di, bytes) != cudaSuccess ||
        cudaMalloc(&dout, bytes) != cudaSuccess) {
        printf("cudaMalloc failed\n"); return 1;
    }

    /* Init with some pattern */
    {
        double *h = (double *)malloc(bytes);
        for (int64_t i = 0; i < TOTAL; i++) h[i] = (double)(i % 997) * 0.001;
        cudaMemcpy(di, h, bytes, cudaMemcpyHostToDevice);
        free(h);
    }

    run_all< 3>(di, dout);
    run_all< 6>(di, dout);
    run_all< 9>(di, dout);
    run_all<12>(di, dout);
    run_all<15>(di, dout);

    cudaFree(di); cudaFree(dout);
    fclose(csv);
    printf("\nwrote results_gpu_oop.csv\n");
}