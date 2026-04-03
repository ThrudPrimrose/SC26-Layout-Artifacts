#include "hip/hip_runtime.h"
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

/*  GPU in-place: negate all K fields.
 *  TOTAL = K * N_base = const.   BW = 2·K·N_base·8 (read+writeback). */
constexpr int64_t TOTAL = 1LL << 27;
constexpr int64_t BLK   = 256;
constexpr int64_t RUNS  = 100;
constexpr int64_t GPU_MAX_VL = 64;

template<int K>
struct AoS { double f[K]; };

template<int K, int VL>
struct AoSoA { double f[K][VL]; };

/* ═══ GPU kernels (in-place) ═══ */

template<int K>
__global__ void g_aos(AoS<K> *__restrict__ buf, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n)
        for (int k = 0; k < K; k++)
            buf[i].f[k] = -buf[i].f[k];
}

template<int K>
__global__ void g_soa(double *__restrict__ base, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n)
        for (int k = 0; k < K; k++)
            base[(int64_t)k * n + i] = -base[(int64_t)k * n + i];
}

template<int K, int VL>
__global__ void g_aosoa(AoSoA<K,VL> *__restrict__ buf, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n) {
        int blk  = i / VL;
        int lane = i % VL;
        for (int k = 0; k < K; k++)
            buf[blk].f[k][lane] = -buf[blk].f[k][lane];
    }
}

/* ═══ Reporting ═══ */

static FILE *csv;

static double bw_ip(int K, int64_t n, double ms) {
    return 2.0 * K * n * sizeof(double) / (ms * 1e6);
}

#define GPU_BENCH(K_val, n_base, label, call) do { \
    for (int w = 0; w < 5; w++) { call; } \
    hipDeviceSynchronize(); \
    for (int r = 0; r < RUNS; r++) { \
        hipEvent_t a, b; hipEventCreate(&a); hipEventCreate(&b); \
        hipEventRecord(a); \
        call; \
        hipEventRecord(b); hipEventSynchronize(b); \
        float ms; hipEventElapsedTime(&ms, a, b); \
        fprintf(csv, "%d,%s,%d,%.6f,%.2f\n", \
                K_val, label, r, (double)ms, bw_ip(K_val, n_base, ms)); \
        hipEventDestroy(a); hipEventDestroy(b); \
    } \
    printf("  %-14s  (see csv)\n", label); \
} while (0)

template<int K>
static void bench_aos(int64_t n, double *dbuf) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(K, n, "AoS",
              (g_aos<K><<<grid, BLK>>>((AoS<K>*)dbuf, n)));
}

template<int K>
static void bench_soa(int64_t n, double *dbuf) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(K, n, "SoA",
              (g_soa<K><<<grid, BLK>>>(dbuf, n)));
}

template<int K, int VL>
static void bench_aosoa(int64_t n, double *dbuf, const char *label) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(K, n, label,
              (g_aosoa<K,VL><<<grid, BLK>>>(
                  (AoSoA<K,VL>*)dbuf, n)));
}

template<int K>
static void run_all(double *dbuf) {
    int64_t n = (TOTAL / K / GPU_MAX_VL) * GPU_MAX_VL;
    printf("\n── K=%d  N_base=%lld  total=%.1f GB ──\n",
           K, (long long)n, (double)K * n * 8 / 1e9);

    bench_aos<K>(n, dbuf);
    bench_soa<K>(n, dbuf);
    bench_aosoa<K,  2>(n, dbuf, "AoSoA-2");
    bench_aosoa<K,  4>(n, dbuf, "AoSoA-4");
    bench_aosoa<K,  8>(n, dbuf, "AoSoA-8");
    bench_aosoa<K, 16>(n, dbuf, "AoSoA-16");
    bench_aosoa<K, 32>(n, dbuf, "AoSoA-32");
    bench_aosoa<K, 64>(n, dbuf, "AoSoA-64");
}

int main() {
    csv = fopen("results_gpu_inplace.csv", "w");
    fprintf(csv, "K,layout,run,ms,gbps\n");

    printf("conj IN-PLACE (GPU): TOTAL=%lldM  runs=%d  dtype=double\n",
           (long long)(TOTAL >> 20), (int)RUNS);

    size_t bytes = TOTAL * sizeof(double);
    double *dbuf;
    if (hipMalloc(&dbuf, bytes) != hipSuccess) {
        printf("hipMalloc failed\n"); return 1;
    }

    {
        double *h = (double *)malloc(bytes);
        for (int64_t i = 0; i < TOTAL; i++) h[i] = (double)(i % 997) * 0.001;
        hipMemcpy(dbuf, h, bytes, hipMemcpyHostToDevice);
        free(h);
    }

    run_all< 3>(dbuf);
    run_all< 6>(dbuf);
    run_all< 9>(dbuf);
    run_all<12>(dbuf);
    run_all<15>(dbuf);

    hipFree(dbuf);
    fclose(csv);
    printf("\nwrote results_gpu_inplace.csv\n");
}