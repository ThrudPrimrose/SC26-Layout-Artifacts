#include "hip/hip_runtime.h"
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <type_traits>

constexpr int64_t N    = 1 << 24;
constexpr int64_t BLK  = 256;
constexpr int64_t GRID = (N + BLK - 1) / BLK;
constexpr int64_t RUNS = 100;

struct C2 { double re, im; };

template<int VL>
struct C2V { double re[VL], im[VL]; };

static_assert(std::is_trivially_copyable<C2>::value, "C2 must be trivially copyable");
static_assert(sizeof(C2) == 2 * sizeof(double), "C2 has unexpected padding");

#define CHECK_C2V(VL) \
    static_assert(std::is_trivially_copyable<C2V<VL>>::value, "C2V<" #VL "> must be trivially copyable"); \
    static_assert(sizeof(C2V<VL>) == 2 * (VL) * sizeof(double), "C2V<" #VL "> has unexpected padding");
CHECK_C2V(2) CHECK_C2V(4) CHECK_C2V(8)
CHECK_C2V(16) CHECK_C2V(32) CHECK_C2V(64)
#undef CHECK_C2V

/* ═══ GPU Kernels ═══ */

__global__ void g_aos(const C2* __restrict__ in, C2* __restrict__ out, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n) out[i] = {in[i].re, -in[i].im};
}

__global__ void g_soa(const double* __restrict__ ri, const double* __restrict__ ii,
                      double* __restrict__ ro, double* __restrict__ io, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n) { ro[i] = ri[i]; io[i] = -ii[i]; }
}

template<int VL>
__global__ void g_aosoa(const C2V<VL>* __restrict__ in, C2V<VL>* __restrict__ out, int n) {
    int i = blockIdx.x * BLK + threadIdx.x;
    if (i < n) {
        int blk = i / VL, lane = i % VL;
        out[blk].re[lane] =  in[blk].re[lane];
        out[blk].im[lane] = -in[blk].im[lane];
    }
}


/* ═══ Reporting ═══ */

static FILE *csv;
static double bw(double ms) { return 4.0 * N * sizeof(double) / (ms * 1e6); }

static void emit(const char *dev, const char *layout, double ms) {
    double g = bw(ms);
    printf("  %-4s %-14s %8.4f ms  %7.1f GB/s\n", dev, layout, ms, g);
    fprintf(csv, "%s,%s,%.6f,%.2f\n", dev, layout, ms, g);
}

/* ═══ GPU bench (CUDA events) ═══ */

#define GPU_BENCH(label, call) do { \
    call; hipDeviceSynchronize(); \
    hipEvent_t a, b; hipEventCreate(&a); hipEventCreate(&b); \
    hipEventRecord(a); \
    for (int r = 0; r < RUNS; r++) { call; } \
    hipEventRecord(b); hipEventSynchronize(b); \
    float ms; hipEventElapsedTime(&ms, a, b); \
    emit("GPU", label, ms / RUNS); \
    hipEventDestroy(a); hipEventDestroy(b); \
} while (0)


int main() {
    size_t bytes = 2ULL * N * sizeof(double);

    double *hi = (double*)malloc(bytes), *ho = (double*)malloc(bytes);
    for (int i = 0; i < 2 * N; i++) hi[i] = (double)(i % 997) * 0.001;

    csv = fopen("results_gpu_oop.csv", "w");
    fprintf(csv, "device,layout,ms,gbps\n");

    printf("conj: N=%lldM  block=%ld  runs=%d dtype=double\n\n",
           (long long)(N >> 20), (long)BLK, (int)RUNS);

    /* --- GPU --- */
    int devcount = 0;
    hipGetDeviceCount(&devcount);
    if (devcount > 0 && hipSetDevice(0) == hipSuccess) {
        double *di, *dout;
        if (hipMalloc(&di, bytes) != hipSuccess ||
            hipMalloc(&dout, bytes) != hipSuccess) {
            printf("[GPU] hipMalloc failed, skipping\n\n");
        } else {
            hipMemcpy(di, hi, bytes, hipMemcpyHostToDevice);
            printf("[GPU]\n");
            GPU_BENCH("AoS",      (g_aos<<<GRID,BLK>>>((C2*)di, (C2*)dout, N)));
            GPU_BENCH("SoA",      (g_soa<<<GRID,BLK>>>(di, di+N, dout, dout+N, N)));
            GPU_BENCH("AoSoA-2",  (g_aosoa< 2><<<GRID,BLK>>>((C2V< 2>*)di, (C2V< 2>*)dout, N)));
            GPU_BENCH("AoSoA-4",  (g_aosoa< 4><<<GRID,BLK>>>((C2V< 4>*)di, (C2V< 4>*)dout, N)));
            GPU_BENCH("AoSoA-8",  (g_aosoa< 8><<<GRID,BLK>>>((C2V< 8>*)di, (C2V< 8>*)dout, N)));
            GPU_BENCH("AoSoA-16", (g_aosoa<16><<<GRID,BLK>>>((C2V<16>*)di, (C2V<16>*)dout, N)));
            GPU_BENCH("AoSoA-32", (g_aosoa<32><<<GRID,BLK>>>((C2V<32>*)di, (C2V<32>*)dout, N)));
            GPU_BENCH("AoSoA-64", (g_aosoa<64><<<GRID,BLK>>>((C2V<64>*)di, (C2V<64>*)dout, N)));
            hipFree(di); hipFree(dout);
        }
    } else {
        printf("[GPU] no device available, skipping\n\n");
    }


    fclose(csv);
    free(hi); free(ho);
    printf("\nwrote results_gpu_oop.csv\n");
}
