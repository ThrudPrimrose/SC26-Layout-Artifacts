/*
 * gpu_triad.cu / gpu_triad.hip -- GPU STREAM Triad benchmark
 *
 * Triad: a[i] = b[i] + scalar * c[i]
 *
 * Build (NVIDIA):
 *   nvcc -O3 -o gpu_triad gpu_triad.cu
 *
 * Build (AMD / HIP):
 *   hipcc -O3 -o gpu_triad gpu_triad.cu
 *
 * Run:
 *   ./gpu_triad [N]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ---- HIP / CUDA portability ---- */
#ifdef __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>
#define GPU(call) do {                                          \
    hipError_t err = (call);                                    \
    if (err != hipSuccess) {                                    \
        fprintf(stderr, "HIP error %s:%d: %s\n",               \
                __FILE__, __LINE__, hipGetErrorString(err));     \
        exit(1);                                                \
    }                                                           \
} while(0)
#define gpuMalloc          hipMalloc
#define gpuFree            hipFree
#define gpuMemcpy          hipMemcpy
#define gpuMemset          hipMemset
#define gpuMemcpyH2D       hipMemcpyHostToDevice
#define gpuMemcpyD2H       hipMemcpyDeviceToHost
#define gpuDeviceSynchronize hipDeviceSynchronize
#define gpuEventCreate     hipEventCreate
#define gpuEventRecord     hipEventRecord
#define gpuEventSynchronize hipEventSynchronize
#define gpuEventElapsedTime hipEventElapsedTime
#define gpuEventDestroy    hipEventDestroy
#define gpuEvent_t         hipEvent_t
#define gpuGetDeviceProperties hipGetDeviceProperties
#define gpuDeviceProp      hipDeviceProp_t
#else
#include <cuda_runtime.h>
#define GPU(call) do {                                          \
    cudaError_t err = (call);                                   \
    if (err != cudaSuccess) {                                   \
        fprintf(stderr, "CUDA error %s:%d: %s\n",              \
                __FILE__, __LINE__, cudaGetErrorString(err));    \
        exit(1);                                                \
    }                                                           \
} while(0)
#define gpuMalloc          cudaMalloc
#define gpuFree            cudaFree
#define gpuMemcpy          cudaMemcpy
#define gpuMemset          cudaMemset
#define gpuMemcpyH2D       cudaMemcpyHostToDevice
#define gpuMemcpyD2H       cudaMemcpyDeviceToHost
#define gpuDeviceSynchronize cudaDeviceSynchronize
#define gpuEventCreate     cudaEventCreate
#define gpuEventRecord     cudaEventRecord
#define gpuEventSynchronize cudaEventSynchronize
#define gpuEventElapsedTime cudaEventElapsedTime
#define gpuEventDestroy    cudaEventDestroy
#define gpuEvent_t         cudaEvent_t
#define gpuGetDeviceProperties cudaGetDeviceProperties
#define gpuDeviceProp      cudaDeviceProp
#endif

/* ================================================================ */
/*  Triad kernels                                                    */
/* ================================================================ */

__global__ void triad_kernel(double* __restrict__ a,
                             const double* __restrict__ b,
                             const double* __restrict__ c,
                             double scalar, size_t N)
{
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = (size_t)blockDim.x * gridDim.x;
    for (; i < N; i += stride)
        a[i] = b[i] + scalar * c[i];
}

/* ================================================================ */
/*  Constants                                                        */
/* ================================================================ */

#ifndef N_DEFAULT
#define N_DEFAULT  (1L << 30)    /* ~1B doubles = 8 GiB per array */
#endif
#define WARMUP     5
#define NTIMES     100
#define SCALAR     3.0

/* ================================================================ */
/*  Benchmark one block-size configuration                           */
/* ================================================================ */

static void bench_config(double* d_a, double* d_b, double* d_c,
                         size_t N, int block_size, double bytes)
{
    int grid_size = (int)((N + block_size - 1) / block_size);

    /* cap grid to avoid huge launch overhead on small N */
    int max_grid = 65535 * 16;
    if (grid_size > max_grid) grid_size = max_grid;

    gpuEvent_t start, stop;
    GPU(gpuEventCreate(&start));
    GPU(gpuEventCreate(&stop));

    /* warm-up */
    for (int r = 0; r < WARMUP; r++) {
        triad_kernel<<<grid_size, block_size>>>(d_a, d_b, d_c, SCALAR, N);
    }
    GPU(gpuDeviceSynchronize());

    double tmin = 1e30, tmax = 0.0, tavg = 0.0;

    for (int r = 0; r < NTIMES; r++) {
        /* reset a */
        GPU(gpuMemset(d_a, 0, N * sizeof(double)));
        GPU(gpuDeviceSynchronize());

        GPU(gpuEventRecord(start, 0));
        triad_kernel<<<grid_size, block_size>>>(d_a, d_b, d_c, SCALAR, N);
        GPU(gpuEventRecord(stop, 0));
        GPU(gpuEventSynchronize(stop));

        float ms;
        GPU(gpuEventElapsedTime(&ms, start, stop));
        double dt = (double)ms / 1e3;  /* seconds */

        if (dt < tmin) tmin = dt;
        if (dt > tmax) tmax = dt;
        tavg += dt;
    }
    tavg /= NTIMES;

    printf("  Block %-4d | Grid %-8d | "
           "Best: %7.3f ms -> %8.2f GB/s | "
           "Avg: %7.3f ms -> %8.2f GB/s | "
           "Worst: %7.3f ms -> %8.2f GB/s\n",
           block_size, grid_size,
           tmin * 1e3, bytes / tmin / 1e9,
           tavg * 1e3, bytes / tavg / 1e9,
           tmax * 1e3, bytes / tmax / 1e9);

    GPU(gpuEventDestroy(start));
    GPU(gpuEventDestroy(stop));
}

/* ================================================================ */
/*  Validation                                                       */
/* ================================================================ */

static int validate(double* d_a, size_t N)
{
    /* only check a subset to avoid huge D2H transfer */
    size_t check_n = (N < 1024 * 1024) ? N : 1024 * 1024;
    double* h_check = (double*)malloc(check_n * sizeof(double));
    GPU(gpuMemcpy(h_check, d_a, check_n * sizeof(double), gpuMemcpyD2H));

    double expected = 1.0 + SCALAR * 2.0;
    int ok = 1;
    for (size_t i = 0; i < check_n; i++) {
        double err = h_check[i] - expected;
        if (err > 1e-12 || err < -1e-12) {
            fprintf(stderr, "VALIDATION FAILED i=%zu: got %.15f expected %.15f\n",
                    i, h_check[i], expected);
            ok = 0;
            break;
        }
    }
    free(h_check);
    return ok;
}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */

int main(int argc, char** argv)
{
    size_t N = N_DEFAULT;
    if (argc > 1) N = atol(argv[1]);

    /* query device */
    gpuDeviceProp prop;
    GPU(gpuGetDeviceProperties(&prop, 0));

    printf("GPU Triad Benchmark\n");
    printf("  Device       : %s\n", prop.name);
    printf("  CUs/SMs      : %d\n", prop.multiProcessorCount);
    printf("  Clock        : %d MHz\n", prop.clockRate / 1000);
    printf("  Global mem   : %.2f GiB\n",
           (double)prop.totalGlobalMem / (1L << 30));
    printf("  Elements     : %zu (%.2f GiB per array)\n",
           N, (double)(N * sizeof(double)) / (1L << 30));
    printf("  Warmup       : %d\n", WARMUP);
    printf("  Repetitions  : %d\n\n", NTIMES);

    /* allocate device arrays */
    double *d_a, *d_b, *d_c;
    GPU(gpuMalloc(&d_a, N * sizeof(double)));
    GPU(gpuMalloc(&d_b, N * sizeof(double)));
    GPU(gpuMalloc(&d_c, N * sizeof(double)));

    /* init on host, copy to device */
    double* h_tmp = (double*)malloc(N * sizeof(double));
    if (!h_tmp) { perror("malloc"); return 1; }

    for (size_t i = 0; i < N; i++) h_tmp[i] = 1.0;
    GPU(gpuMemcpy(d_b, h_tmp, N * sizeof(double), gpuMemcpyH2D));

    for (size_t i = 0; i < N; i++) h_tmp[i] = 2.0;
    GPU(gpuMemcpy(d_c, h_tmp, N * sizeof(double), gpuMemcpyH2D));

    GPU(gpuMemset(d_a, 0, N * sizeof(double)));
    free(h_tmp);

    /* 2 reads + 1 write = 3 streams */
    double bytes = 3.0 * (double)N * sizeof(double);

    printf("Results (Triad: a = b + s*c)\n");

    /* bench both block sizes */
    int block_sizes[] = {128, 256};
    int n_configs = sizeof(block_sizes) / sizeof(block_sizes[0]);

    for (int i = 0; i < n_configs; i++) {
        bench_config(d_a, d_b, d_c, N, block_sizes[i], bytes);
    }

    /* validate last run */
    printf("\n");
    if (validate(d_a, N))
        printf("VALIDATION PASSED\n");

    GPU(gpuFree(d_a));
    GPU(gpuFree(d_b));
    GPU(gpuFree(d_c));

    return 0;
}
