/*
 * gpu_triad.cu -- GPU STREAM Triad benchmark
 *
 * Triad: a[i] = b[i] + scalar * c[i]
 *
 * Build (AMD / HIP):
 *   hipcc -O3 -o gpu_triad gpu_triad.cu
 *
 * Build (NVIDIA):
 *   nvcc -O3 -o gpu_triad gpu_triad.cu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <hip/hip_runtime.h>


#define HIP_CHECK(call) do {                                        \
    hipError_t err = (call);                                        \
    if (err != hipSuccess) {                                        \
        fprintf(stderr, "GPU error at %s:%d: %s\n",                 \
                __FILE__, __LINE__, hipGetErrorString(err));         \
        exit(1);                                                    \
    }                                                               \
} while(0)

/* Also check for async errors after kernel launches */
#define KERNEL_CHECK() do {                                         \
    HIP_CHECK(hipGetLastError());                                   \
    HIP_CHECK(hipDeviceSynchronize());                              \
} while(0)

/* ================================================================ */
/*  Triad kernel                                                     */
/* ================================================================ */

__global__ void triad_kernel(double* __restrict__ a,
                             const double* __restrict__ b,
                             const double* __restrict__ c,
                             const double scalar, const long N)
{
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long stride = (long)blockDim.x * gridDim.x;
    for (; i < N; i += stride)
        a[i] = b[i] + scalar * c[i];
}

/* ================================================================ */
/*  Init kernel (avoid huge host->device copy)                       */
/* ================================================================ */

__global__ void init_kernel(double* __restrict__ a,
                            double* __restrict__ b, double* __restrict__ c,
                            double va, double vb, double vc,
                            const long N)
{
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long stride = (long)blockDim.x * gridDim.x;
    for (; i < N; i += stride) {
        a[i] = va;
        b[i] = vb;
        c[i] = vc;
    }
}

/* ================================================================ */
/*  Constants                                                        */
/* ================================================================ */

#ifndef N_DEFAULT
#define N_DEFAULT  (1L << 30)
#endif
#define WARMUP     5
#define NTIMES     100
#define SCALAR     3.0

/* ================================================================ */
/*  Benchmark one block-size configuration                           */
/* ================================================================ */

static void bench_config(double* d_a, double* d_b, double* d_c,
                         long N, int block_size, double bytes)
{
    unsigned int grid_size = ((unsigned long)N + block_size - 1) / block_size;

    /* reasonable cap -- grid-stride loop handles the rest */
    if (grid_size > 65536u) grid_size = 65536u;

    dim3 grid(grid_size);
    dim3 block(block_size);

    /* verify launch config works */
    printf("  Testing block=%d grid=%u ... ", block_size, grid_size);
    fflush(stdout);

    triad_kernel<<<grid, block>>>(d_a, d_b, d_c, SCALAR, N);
    HIP_CHECK(hipGetLastError());
    HIP_CHECK(hipDeviceSynchronize());
    printf("launch OK\n");

    hipEvent_t start, stop;
    HIP_CHECK(hipEventCreate(&start));
    HIP_CHECK(hipEventCreate(&stop));

    /* warm-up */
    for (int r = 0; r < WARMUP; r++) {
        triad_kernel<<<grid, block>>>(d_a, d_b, d_c, SCALAR, N);
        HIP_CHECK(hipGetLastError());
    }
    HIP_CHECK(hipDeviceSynchronize());

    double tmin = 1e30, tmax = 0.0, tavg = 0.0;

    for (int r = 0; r < NTIMES; r++) {
        /* reset a[] on device */
        HIP_CHECK(hipMemset(d_a, 0, (size_t)N * sizeof(double)));
        HIP_CHECK(hipDeviceSynchronize());

        HIP_CHECK(hipEventRecord(start, 0));
        triad_kernel<<<grid, block>>>(d_a, d_b, d_c, SCALAR, N);
        HIP_CHECK(hipGetLastError());
        HIP_CHECK(hipEventRecord(stop, 0));
        HIP_CHECK(hipEventSynchronize(stop));

        float ms = 0.0f;
        HIP_CHECK(hipEventElapsedTime(&ms, start, stop));
        double dt = (double)ms / 1e3;

        if (dt < tmin) tmin = dt;
        if (dt > tmax) tmax = dt;
        tavg += dt;
    }
    tavg /= NTIMES;

    printf("  Block %-4d | Grid %-8u | "
           "Best: %7.3f ms -> %8.2f GB/s | "
           "Avg: %7.3f ms -> %8.2f GB/s | "
           "Worst: %7.3f ms -> %8.2f GB/s\n",
           block_size, grid_size,
           tmin * 1e3, bytes / tmin / 1e9,
           tavg * 1e3, bytes / tavg / 1e9,
           tmax * 1e3, bytes / tmax / 1e9);

    HIP_CHECK(hipEventDestroy(start));
    HIP_CHECK(hipEventDestroy(stop));
}

/* ================================================================ */
/*  Validation                                                       */
/* ================================================================ */

static int validate(double* d_a, long N)
{
    size_t check_n = (N < 1024 * 1024) ? (size_t)N : 1024 * 1024;

    /* check head */
    double* h_check = (double*)malloc(check_n * sizeof(double));
    HIP_CHECK(hipMemcpy(h_check, d_a, check_n * sizeof(double),
                        hipMemcpyDeviceToHost));

    double expected = 1.0 + SCALAR * 2.0;
    int ok = 1;
    for (size_t i = 0; i < check_n; i++) {
        double err = h_check[i] - expected;
        if (err > 1e-12 || err < -1e-12) {
            fprintf(stderr, "VALIDATION FAILED at head i=%zu: "
                    "got %.15f expected %.15f\n",
                    i, h_check[i], expected);
            ok = 0;
            break;
        }
    }

    /* also check tail */
    if (ok && N > (long)check_n) {
        HIP_CHECK(hipMemcpy(h_check, d_a + (N - (long)check_n),
                            check_n * sizeof(double),
                            hipMemcpyDeviceToHost));
        for (size_t i = 0; i < check_n; i++) {
            double err = h_check[i] - expected;
            if (err > 1e-12 || err < -1e-12) {
                fprintf(stderr, "VALIDATION FAILED at tail offset=%ld i=%zu: "
                        "got %.15f expected %.15f\n",
                        N - (long)check_n, i, h_check[i], expected);
                ok = 0;
                break;
            }
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
    long N = N_DEFAULT;
    if (argc > 1) N = atol(argv[1]);

    hipDeviceProp_t prop;
    HIP_CHECK(hipGetDeviceProperties(&prop, 0));

    printf("GPU Triad Benchmark\n");
    printf("  Device       : %s\n", prop.name);
    printf("  CUs/SMs      : %d\n", prop.multiProcessorCount);
    printf("  Clock        : %d MHz\n", prop.clockRate / 1000);
    printf("  Global mem   : %.2f GiB\n",
           (double)prop.totalGlobalMem / (1L << 30));
    printf("  Elements     : %ld (%.2f GiB per array)\n",
           N, (double)(N * sizeof(double)) / (1L << 30));
    printf("  Warmup       : %d\n", WARMUP);
    printf("  Repetitions  : %d\n\n", NTIMES);

    /* allocate device arrays */
    double *d_a, *d_b, *d_c;
    printf("Allocating 3 x %.2f GiB ...\n",
           (double)(N * sizeof(double)) / (1L << 30));
    HIP_CHECK(hipMalloc(&d_a, (size_t)N * sizeof(double)));
    HIP_CHECK(hipMalloc(&d_b, (size_t)N * sizeof(double)));
    HIP_CHECK(hipMalloc(&d_c, (size_t)N * sizeof(double)));
    printf("Allocation OK\n\n");

    /* init on device (avoid 8 GiB host->device copies) */
    {
        unsigned int ig = 65536;
        int ib = 256;
        printf("Initializing arrays on device ...\n");
        init_kernel<<<dim3(ig), dim3(ib)>>>(d_a, d_b, d_c,
                                             0.0, 1.0, 2.0, N);
        HIP_CHECK(hipGetLastError());
        HIP_CHECK(hipDeviceSynchronize());
        printf("Init OK\n\n");
    }

    /* quick sanity: read back one element */
    {
        double tmp[3];
        HIP_CHECK(hipMemcpy(&tmp[0], d_b, sizeof(double),
                            hipMemcpyDeviceToHost));
        HIP_CHECK(hipMemcpy(&tmp[1], d_c, sizeof(double),
                            hipMemcpyDeviceToHost));
        HIP_CHECK(hipMemcpy(&tmp[2], d_a, sizeof(double),
                            hipMemcpyDeviceToHost));
        printf("Sanity check: b[0]=%.1f c[0]=%.1f a[0]=%.1f\n\n",
               tmp[0], tmp[1], tmp[2]);
    }

    /* 2 reads + 1 write = 3 streams */
    double bytes = 3.0 * (double)N * sizeof(double);

    printf("Results (Triad: a = b + s*c)\n");

    int block_sizes[] = {128, 256};
    int n_configs = sizeof(block_sizes) / sizeof(block_sizes[0]);

    for (int i = 0; i < n_configs; i++) {
        bench_config(d_a, d_b, d_c, N, block_sizes[i], bytes);
    }

    printf("\n");
    if (validate(d_a, N))
        printf("VALIDATION PASSED\n");

    HIP_CHECK(hipFree(d_a));
    HIP_CHECK(hipFree(d_b));
    HIP_CHECK(hipFree(d_c));

    return 0;
}
