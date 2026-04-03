#include "hip/hip_runtime.h"
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#ifdef __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>
#define hipMalloc       hipMalloc
#define hipFree         hipFree
#define hipMemcpy       hipMemcpy
#define hipMemcpyHostToDevice hipMemcpyHostToDevice
#define hipDeviceSynchronize  hipDeviceSynchronize
#define hipGetDeviceCount     hipGetDeviceCount
#define hipSetDevice          hipSetDevice
#define hipEvent_t      hipEvent_t
#define hipEventCreate  hipEventCreate
#define hipEventDestroy hipEventDestroy
#define hipEventRecord  hipEventRecord
#define hipEventSynchronize hipEventSynchronize
#define hipEventElapsedTime hipEventElapsedTime
#define hipSuccess      hipSuccess
#else
#include <hip/hip_runtime.h>
#endif

/*  Conjugate P complex arrays out-of-place on GPU.
 *  out.re_p = in.re_p  (copy),  out.im_p = -in.im_p  (negate).
 *  Total doubles = 2·P·N_base = const.
 *  OOP BW = 4·P·N_base·8.
 *  VL sweep: 2 … 64.                                                 */

constexpr int64_t TOTAL_DOUBLES = 1LL << 29;   // ~1 GB per side
constexpr int64_t BLK   = 256;
constexpr int64_t RUNS  = 100;
constexpr int64_t GPU_MAX_VL = 64;

/* ═══ Layout structs ═══ */

template<int P>
struct AoS {
    struct { double re, im; } c[P];
};

template<int P, int VL>
struct AoSoA {
    struct { double re[VL], im[VL]; } c[P];
};

/* ═══ GPU kernels ═══ */

template<int P>
__global__ void g_aos(const AoS<P> *__restrict__ in,
                      AoS<P> *__restrict__ out, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n)
        for (int p = 0; p < P; p++) {
            out[i].c[p].re =  in[i].c[p].re;
            out[i].c[p].im = -in[i].c[p].im;
        }
}

/* SoA: 2P arrays packed contiguously.  re_p at offset p*n, im_p at (P+p)*n */
template<int P>
__global__ void g_soa(const double *__restrict__ base_in,
                      double *__restrict__ base_out, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n)
        for (int p = 0; p < P; p++) {
            base_out[(int64_t)p * n + i]       =  base_in[(int64_t)p * n + i];         /* re */
            base_out[(int64_t)(P + p) * n + i] = -base_in[(int64_t)(P + p) * n + i];   /* im */
        }
}

template<int P, int VL>
__global__ void g_aosoa(const AoSoA<P,VL> *__restrict__ in,
                        AoSoA<P,VL> *__restrict__ out, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n) {
        int blk  = i / VL;
        int lane = i % VL;
        for (int p = 0; p < P; p++) {
            out[blk].c[p].re[lane] =  in[blk].c[p].re[lane];
            out[blk].c[p].im[lane] = -in[blk].c[p].im[lane];
        }
    }
}

/* ═══ Reporting ═══ */

static FILE *csv;

static double bw_oop(int P, int64_t n, double ms) {
    return 4.0 * P * n * sizeof(double) / (ms * 1e6);
}

#define GPU_BENCH(P_val, n_base, label, call) do { \
    for (int w = 0; w < 5; w++) { call; } \
    hipDeviceSynchronize(); \
    for (int r = 0; r < RUNS; r++) { \
        hipEvent_t a, b; hipEventCreate(&a); hipEventCreate(&b); \
        hipEventRecord(a); \
        call; \
        hipEventRecord(b); hipEventSynchronize(b); \
        float ms; hipEventElapsedTime(&ms, a, b); \
        fprintf(csv, "%d,%s,%d,%.6f,%.2f\n", \
                P_val, label, r, (double)ms, bw_oop(P_val, n_base, ms)); \
        hipEventDestroy(a); hipEventDestroy(b); \
    } \
    printf("  %-14s  (see csv)\n", label); \
} while (0)

template<int P>
static void bench_aos(int64_t n, double *di, double *dout) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(P, n, "AoS",
              (g_aos<P><<<grid, BLK>>>((AoS<P>*)di, (AoS<P>*)dout, n)));
}

template<int P>
static void bench_soa(int64_t n, double *di, double *dout) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(P, n, "SoA",
              (g_soa<P><<<grid, BLK>>>(di, dout, n)));
}

template<int P, int VL>
static void bench_aosoa(int64_t n, double *di, double *dout, const char *label) {
    int64_t grid = (n + BLK - 1) / BLK;
    GPU_BENCH(P, n, label,
              (g_aosoa<P,VL><<<grid, BLK>>>(
                  (AoSoA<P,VL>*)di, (AoSoA<P,VL>*)dout, n)));
}

template<int P>
static void run_all(double *di, double *dout) {
    int64_t n = (TOTAL_DOUBLES / (2 * P) / GPU_MAX_VL) * GPU_MAX_VL;
    printf("\n── P=%d complex pairs  (%d SoA streams)  N_base=%lld  "
           "total=%.1f GB/side ──\n",
           P, 4*P, (long long)n, 2.0 * P * n * 8 / 1e9);

    bench_aos<P>(n, di, dout);
    bench_soa<P>(n, di, dout);
    bench_aosoa<P,  2>(n, di, dout, "AoSoA-2");
    bench_aosoa<P,  4>(n, di, dout, "AoSoA-4");
    bench_aosoa<P,  8>(n, di, dout, "AoSoA-8");
    bench_aosoa<P, 16>(n, di, dout, "AoSoA-16");
    bench_aosoa<P, 32>(n, di, dout, "AoSoA-32");
    bench_aosoa<P, 64>(n, di, dout, "AoSoA-64");
}

int main() {
    csv = fopen("results_gpu_oop.csv", "w");
    fprintf(csv, "P,layout,run,ms,gbps\n");

    printf("conjugate OOP (GPU): TOTAL_DOUBLES=%lldM  runs=%d\n",
           (long long)(TOTAL_DOUBLES >> 20), (int)RUNS);

    size_t bytes = TOTAL_DOUBLES * sizeof(double);
    double *di, *dout;
    if (hipMalloc(&di, bytes) != hipSuccess ||
        hipMalloc(&dout, bytes) != hipSuccess) {
        printf("hipMalloc failed\n"); return 1;
    }

    {
        double *h = (double *)malloc(bytes);
        for (int64_t i = 0; i < TOTAL_DOUBLES; i++)
            h[i] = (double)(i % 997) * 0.001;
        hipMemcpy(di, h, bytes, hipMemcpyHostToDevice);
        free(h);
    }

    run_all< 3>(di, dout);
    run_all< 6>(di, dout);
    run_all< 9>(di, dout);
    run_all<12>(di, dout);
    run_all<15>(di, dout);
    run_all<18>(di, dout);
    run_all<21>(di, dout);

    hipFree(di); hipFree(dout);
    fclose(csv);
    printf("\nwrote results_gpu_oop.csv\n");
}