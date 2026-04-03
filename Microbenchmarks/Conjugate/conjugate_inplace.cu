#include <cstdio>
#include <cstdlib>
#include <cstdint>

#ifdef __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>
#define cudaMalloc       hipMalloc
#define cudaFree         hipFree
#define cudaMemcpy       hipMemcpy
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaDeviceSynchronize  hipDeviceSynchronize
#define cudaGetDeviceCount     hipGetDeviceCount
#define cudaSetDevice          hipSetDevice
#define cudaEvent_t      hipEvent_t
#define cudaEventCreate  hipEventCreate
#define cudaEventDestroy hipEventDestroy
#define cudaEventRecord  hipEventRecord
#define cudaEventSynchronize hipEventSynchronize
#define cudaEventElapsedTime hipEventElapsedTime
#define cudaSuccess      hipSuccess
#else
#include <cuda_runtime.h>
#endif

/*  Conjugate P complex arrays in-place on GPU.
 *      buf.im_p = -buf.im_p   (negate im only)
 *  Total doubles = 2·P·N_base = const.
 *  Useful BW = 2·P·N_base·8  (read im + write im).                  */

constexpr int64_t TOTAL_DOUBLES = 1LL << 27;
constexpr int64_t BLK   = 256;
constexpr int64_t RUNS  = 100;
constexpr int64_t GPU_MAX_VL = 64;

template<int P>
struct AoS {
    struct { double re, im; } c[P];
};

template<int P, int VL>
struct AoSoA {
    struct { double re[VL], im[VL]; } c[P];
};

/* ═══ GPU kernels (in-place): negate im only ═══ */

template<int P>
__global__ void g_aos(AoS<P> *__restrict__ buf, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n)
        for (int p = 0; p < P; p++)
            buf[i].c[p].im = -buf[i].c[p].im;
}

/* SoA: re at offset p*n, im at (P+p)*n — only touch im */
template<int P>
__global__ void g_soa(double *__restrict__ base, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n)
        for (int p = 0; p < P; p++)
            base[(int64_t)(P + p) * n + i] = -base[(int64_t)(P + p) * n + i];
}

template<int P, int VL>
__global__ void g_aosoa(AoSoA<P,VL> *__restrict__ buf, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n) {
        int blk  = i / VL;
        int lane = i % VL;
        for (int p = 0; p < P; p++)
            buf[blk].c[p].im[lane] = -buf[blk].c[p].im[lane];
    }
}

/* ═══ Reporting ═══ */

static FILE *csv;

static double bw_ip(int P, int64_t n, double ms) {
    return 2.0 * P * n * sizeof(double) / (ms * 1e6);
}

#define GPU_BENCH(P_val, n_base, label, call) do { \
    for (int w = 0; w < 5; w++) { call; } \
    cudaDeviceSynchronize(); \
    for (int r = 0; r < RUNS; r++) { \
        cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b); \
        cudaEventRecord(a); \
        call; \
        cudaEventRecord(b); cudaEventSynchronize(b); \
        float ms; cudaEventElapsedTime(&ms, a, b); \
        fprintf(csv, "%d,%s,%d,%.6f,%.2f\n", \
                P_val, label, r, (double)ms, bw_ip(P_val, n_base, ms)); \
        cudaEventDestroy(a); cudaEventDestroy(b); \
    } \
    printf("  %-14s  (see csv)\n", label); \
} while (0)

template<int P>
static void bench_aos(int64_t n, double *dbuf) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(P, n, "AoS",
              (g_aos<P><<<grid, BLK>>>((AoS<P>*)dbuf, n)));
}

template<int P>
static void bench_soa(int64_t n, double *dbuf) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(P, n, "SoA",
              (g_soa<P><<<grid, BLK>>>(dbuf, n)));
}

template<int P, int VL>
static void bench_aosoa(int64_t n, double *dbuf, const char *label) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(P, n, label,
              (g_aosoa<P,VL><<<grid, BLK>>>(
                  (AoSoA<P,VL>*)dbuf, n)));
}

template<int P>
static void run_all(double *dbuf) {
    int64_t n = (TOTAL_DOUBLES / (2 * P) / GPU_MAX_VL) * GPU_MAX_VL;
    printf("\n── P=%d complex pairs  (%d im streams)  N_base=%lld  "
           "total=%.1f GB ──\n",
           P, P, (long long)n, 2.0 * P * n * 8 / 1e9);

    bench_aos<P>(n, dbuf);
    bench_soa<P>(n, dbuf);
    bench_aosoa<P,  2>(n, dbuf, "AoSoA-2");
    bench_aosoa<P,  4>(n, dbuf, "AoSoA-4");
    bench_aosoa<P,  8>(n, dbuf, "AoSoA-8");
    bench_aosoa<P, 16>(n, dbuf, "AoSoA-16");
    bench_aosoa<P, 32>(n, dbuf, "AoSoA-32");
    bench_aosoa<P, 64>(n, dbuf, "AoSoA-64");
}

int main() {
    csv = fopen("results_gpu_inplace.csv", "w");
    fprintf(csv, "P,layout,run,ms,gbps\n");

    printf("conjugate IN-PLACE (GPU): TOTAL_DOUBLES=%lldM  runs=%d\n",
           (long long)(TOTAL_DOUBLES >> 20), (int)RUNS);

    size_t bytes = TOTAL_DOUBLES * sizeof(double);
    double *dbuf;
    if (cudaMalloc(&dbuf, bytes) != cudaSuccess) {
        printf("cudaMalloc failed\n"); return 1;
    }

    {
        double *h = (double *)malloc(bytes);
        for (int64_t i = 0; i < TOTAL_DOUBLES; i++)
            h[i] = (double)(i % 997) * 0.001;
        cudaMemcpy(dbuf, h, bytes, cudaMemcpyHostToDevice);
        free(h);
    }

    run_all< 3>(dbuf);
    run_all< 6>(dbuf);
    run_all< 9>(dbuf);
    run_all<12>(dbuf);
    run_all<15>(dbuf);
    run_all<18>(dbuf);
    run_all<21>(dbuf);

    cudaFree(dbuf);
    fclose(csv);
    printf("\nwrote results_gpu_inplace.csv\n");
}