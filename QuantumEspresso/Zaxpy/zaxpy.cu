// zaxpy_indirect_sweep.cu
// Sweep threadblock sizes × coarsening factors for indirect zaxpy
//
// Four kernel variants (2 layouts × 2 index orders):
//   aos_scatter  — cuDoubleComplex (interleaved re/im), original random indices
//   aos_sorted   — cuDoubleComplex, sorted y_index_map for write locality
//   soa_scatter  — separate real[] / imag[] arrays, original random indices
//   soa_sorted   — separate real[] / imag[] arrays, sorted y_index_map
//
// Each individual repetition is logged to CSV.
//
// Compile:
//   nvcc -O3 -std=c++17 zaxpy_indirect_sweep.cu -o zaxpy_indirect_sweep

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <numeric>
#include <cuComplex.h>
#include <array>
#include <cmath>

#define CUDA_CHECK(call)                                     \
    do                                                       \
    {                                                        \
        cudaError_t err = (call);                            \
        if (err != cudaSuccess)                              \
        {                                                    \
            std::cerr << "CUDA error at " << __FILE__ << ":" \
                      << __LINE__ << " -> "                  \
                      << cudaGetErrorString(err)             \
                      << std::endl;                          \
            std::exit(EXIT_FAILURE);                         \
        }                                                    \
    } while (0)

/* ================================================================ */
/*  Cache-flush helper (write a big buffer to evict L2)             */
/* ================================================================ */
static double *d_flush = nullptr;
static constexpr size_t FLUSH_SIZE = 256ULL * 1024 * 1024; // 256 MiB

__global__ void flush_kernel(double *buf, size_t n)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) buf[i] = 1.0;
}

static void init_flush()
{
    if (!d_flush)
        CUDA_CHECK(cudaMalloc(&d_flush, FLUSH_SIZE));
}

static void flush_l2()
{
    size_t n = FLUSH_SIZE / sizeof(double);
    flush_kernel<<<(n + 255) / 256, 256>>>(d_flush, n);
    CUDA_CHECK(cudaDeviceSynchronize());
}

/* ================================================================ */
/*  AoS kernels (cuDoubleComplex — interleaved re,im)               */
/* ================================================================ */

template <int COARSEN>
__global__ void kern_aos_scatter(
    int nx,
    const int *__restrict__ y_index_map,
    const cuDoubleComplex *__restrict__ x,
    cuDoubleComplex *__restrict__ y)
{
    int base = (blockIdx.x * blockDim.x + threadIdx.x) * COARSEN;
    #pragma unroll
    for (int c = 0; c < COARSEN; c++)
    {
        int gid = base + c;
        if (gid < nx)
        {
            int yi = y_index_map[gid];
            y[yi] = cuCadd(x[gid], y[yi]);
        }
    }
}

template <int COARSEN>
__global__ void kern_aos_sorted(
    int nx,
    const int *__restrict__ y_index_map_sorted,
    const int *__restrict__ x_index_map,
    const cuDoubleComplex *__restrict__ x,
    cuDoubleComplex *__restrict__ y)
{
    int base = (blockIdx.x * blockDim.x + threadIdx.x) * COARSEN;
    #pragma unroll
    for (int c = 0; c < COARSEN; c++)
    {
        int gid = base + c;
        if (gid < nx)
        {
            int yi = y_index_map_sorted[gid];
            y[yi] = cuCadd(x[x_index_map[gid]], y[yi]);
        }
    }
}

/* ================================================================ */
/*  SoA kernels (separate double *re, double *im arrays)            */
/* ================================================================ */

template <int COARSEN>
__global__ void kern_soa_scatter(
    int nx,
    const int *__restrict__ y_index_map,
    const double *__restrict__ x_re,
    const double *__restrict__ x_im,
    double *__restrict__ y_re,
    double *__restrict__ y_im)
{
    int base = (blockIdx.x * blockDim.x + threadIdx.x) * COARSEN;
    #pragma unroll
    for (int c = 0; c < COARSEN; c++)
    {
        int gid = base + c;
        if (gid < nx)
        {
            int yi = y_index_map[gid];
            y_re[yi] += x_re[gid];
            y_im[yi] += x_im[gid];
        }
    }
}

template <int COARSEN>
__global__ void kern_soa_sorted(
    int nx,
    const int *__restrict__ y_index_map_sorted,
    const int *__restrict__ x_index_map,
    const double *__restrict__ x_re,
    const double *__restrict__ x_im,
    double *__restrict__ y_re,
    double *__restrict__ y_im)
{
    int base = (blockIdx.x * blockDim.x + threadIdx.x) * COARSEN;
    #pragma unroll
    for (int c = 0; c < COARSEN; c++)
    {
        int gid = base + c;
        if (gid < nx)
        {
            int yi = y_index_map_sorted[gid];
            int xi = x_index_map[gid];
            y_re[yi] += x_re[xi];
            y_im[yi] += x_im[xi];
        }
    }
}

/* ================================================================ */
/*  Dispatch tables                                                  */
/* ================================================================ */

using aos_scatter_fn = void (*)(int, const int*, const cuDoubleComplex*, cuDoubleComplex*);
using aos_sorted_fn  = void (*)(int, const int*, const int*, const cuDoubleComplex*, cuDoubleComplex*);
using soa_scatter_fn = void (*)(int, const int*, const double*, const double*, double*, double*);
using soa_sorted_fn  = void (*)(int, const int*, const int*, const double*, const double*, double*, double*);

#define ENTRY(K, C) K<C>

aos_scatter_fn aos_scatter_tbl[] = {nullptr, ENTRY(kern_aos_scatter,1), ENTRY(kern_aos_scatter,2), nullptr, ENTRY(kern_aos_scatter,4), nullptr,nullptr,nullptr, ENTRY(kern_aos_scatter,8)};
aos_sorted_fn  aos_sorted_tbl[]  = {nullptr, ENTRY(kern_aos_sorted,1),  ENTRY(kern_aos_sorted,2),  nullptr, ENTRY(kern_aos_sorted,4),  nullptr,nullptr,nullptr, ENTRY(kern_aos_sorted,8)};
soa_scatter_fn soa_scatter_tbl[] = {nullptr, ENTRY(kern_soa_scatter,1), ENTRY(kern_soa_scatter,2), nullptr, ENTRY(kern_soa_scatter,4), nullptr,nullptr,nullptr, ENTRY(kern_soa_scatter,8)};
soa_sorted_fn  soa_sorted_tbl[]  = {nullptr, ENTRY(kern_soa_sorted,1),  ENTRY(kern_soa_sorted,2),  nullptr, ENTRY(kern_soa_sorted,4),  nullptr,nullptr,nullptr, ENTRY(kern_soa_sorted,8)};

/* ================================================================ */
/*  Verification helpers                                             */
/* ================================================================ */

bool verify_aos(int n, cuDoubleComplex *gpu, cuDoubleComplex *ref, double tol = 1e-6)
{
    for (int i = 0; i < n; i++)
    {
        double diff = cuCabs(cuCsub(gpu[i], ref[i]));
        double rel  = diff / (cuCabs(ref[i]) + 1e-10);
        if (rel > tol) { printf("AoS mismatch at %d: rel=%e\n", i, rel); return false; }
    }
    return true;
}

bool verify_soa(int n, double *gpu_re, double *gpu_im, double *ref_re, double *ref_im, double tol = 1e-6)
{
    for (int i = 0; i < n; i++)
    {
        double dr = gpu_re[i] - ref_re[i];
        double di = gpu_im[i] - ref_im[i];
        double diff = std::sqrt(dr * dr + di * di);
        double mag  = std::sqrt(ref_re[i] * ref_re[i] + ref_im[i] * ref_im[i]) + 1e-10;
        if (diff / mag > tol) { printf("SoA mismatch at %d: rel=%e\n", i, diff / mag); return false; }
    }
    return true;
}

/* ================================================================ */
/*  Core profiling: single (tpb, coarsen) config, all 4 variants    */
/* ================================================================ */

void profile_config(
    FILE *csv,
    const char *experiment_name,
    int ny, int nx,
    int *h_ymap,
    int *h_ymap_sorted,
    int *h_xmap,
    int tpb, int coarsen,
    int iters, int warmup,
    unsigned int seed)
{
    int elems_per_block = tpb * coarsen;
    int num_blocks = (nx + elems_per_block - 1) / elems_per_block;

    // ---- Host data ----
    cuDoubleComplex *h_y_aos   = (cuDoubleComplex *)malloc(ny * sizeof(cuDoubleComplex));
    cuDoubleComplex *h_x_aos   = (cuDoubleComplex *)malloc(nx * sizeof(cuDoubleComplex));
    cuDoubleComplex *h_out_aos = (cuDoubleComplex *)malloc(ny * sizeof(cuDoubleComplex));
    cuDoubleComplex *h_ref_aos = (cuDoubleComplex *)malloc(ny * sizeof(cuDoubleComplex));

    double *h_y_re  = (double *)malloc(ny * sizeof(double));
    double *h_y_im  = (double *)malloc(ny * sizeof(double));
    double *h_x_re  = (double *)malloc(nx * sizeof(double));
    double *h_x_im  = (double *)malloc(nx * sizeof(double));
    double *h_out_re = (double *)malloc(ny * sizeof(double));
    double *h_out_im = (double *)malloc(ny * sizeof(double));
    double *h_ref_re = (double *)malloc(ny * sizeof(double));
    double *h_ref_im = (double *)malloc(ny * sizeof(double));

    // ---- Device data ----
    cuDoubleComplex *d_y_aos, *d_x_aos;
    double *d_y_re, *d_y_im, *d_x_re, *d_x_im;
    int *d_ymap, *d_ymap_s, *d_xmap;

    CUDA_CHECK(cudaMalloc(&d_y_aos, ny * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_x_aos, nx * sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(&d_y_re,  ny * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_y_im,  ny * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_x_re,  nx * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_x_im,  nx * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_ymap,   nx * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_ymap_s, nx * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_xmap,   nx * sizeof(int)));

    // ---- Fill random data ----
    srand(seed);
    for (int i = 0; i < ny; i++)
    {
        double re = (double)rand() / RAND_MAX;
        double im = (double)rand() / RAND_MAX;
        h_y_aos[i]   = make_cuDoubleComplex(re, im);
        h_ref_aos[i] = h_y_aos[i];
        h_y_re[i] = re;  h_y_im[i] = im;
        h_ref_re[i] = re; h_ref_im[i] = im;
    }
    for (int i = 0; i < nx; i++)
    {
        double re = (double)rand() / RAND_MAX;
        double im = (double)rand() / RAND_MAX;
        h_x_aos[i] = make_cuDoubleComplex(re, im);
        h_x_re[i] = re;  h_x_im[i] = im;
    }

    // ---- CPU reference ----
    for (int i = 0; i < nx; i++)
    {
        int yi = h_ymap[i];
        h_ref_aos[yi] = cuCadd(h_x_aos[i], h_ref_aos[yi]);
        h_ref_re[yi] += h_x_re[i];
        h_ref_im[yi] += h_x_im[i];
    }

    // ---- Upload index maps ----
    CUDA_CHECK(cudaMemcpy(d_ymap,   h_ymap,        nx * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ymap_s, h_ymap_sorted, nx * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_xmap,   h_xmap,        nx * sizeof(int), cudaMemcpyHostToDevice));

    cudaEvent_t e0, e1;
    CUDA_CHECK(cudaEventCreate(&e0));
    CUDA_CHECK(cudaEventCreate(&e1));
    float ms;

    auto run_variant = [&](const char *vname, auto launch_fn, auto verify_fn)
    {
        if (!verify_fn())
        {
            fprintf(stderr, "FAIL %s tpb=%d coarsen=%d\n", vname, tpb, coarsen);
            return;
        }
        for (int w = 0; w < warmup; w++) { launch_fn(); CUDA_CHECK(cudaDeviceSynchronize()); }
        for (int r = 0; r < iters; r++)
        {
            flush_l2();
            CUDA_CHECK(cudaEventRecord(e0));
            launch_fn();
            CUDA_CHECK(cudaEventRecord(e1));
            CUDA_CHECK(cudaEventSynchronize(e1));
            CUDA_CHECK(cudaEventElapsedTime(&ms, e0, e1));
            fprintf(csv, "%s,%s,%d,%d,%d,%d,%d,%.6f\n",
                    vname, experiment_name, ny, nx, tpb, coarsen, r, (double)ms);
        }
    };

    // =========== AoS scatter ===========
    {
        auto kfn = aos_scatter_tbl[coarsen];
        CUDA_CHECK(cudaMemcpy(d_x_aos, h_x_aos, nx * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
        auto launch = [&]() { kfn<<<num_blocks, tpb>>>(nx, d_ymap, d_x_aos, d_y_aos); };
        auto verify = [&]() -> bool {
            CUDA_CHECK(cudaMemcpy(d_y_aos, h_y_aos, ny * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
            launch(); CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_out_aos, d_y_aos, ny * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
            return verify_aos(ny, h_out_aos, h_ref_aos);
        };
        CUDA_CHECK(cudaMemcpy(d_y_aos, h_y_aos, ny * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
        run_variant("aos_scatter", launch, verify);
    }

    // =========== AoS sorted ===========
    {
        auto kfn = aos_sorted_tbl[coarsen];
        CUDA_CHECK(cudaMemcpy(d_x_aos, h_x_aos, nx * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
        auto launch = [&]() { kfn<<<num_blocks, tpb>>>(nx, d_ymap_s, d_xmap, d_x_aos, d_y_aos); };
        auto verify = [&]() -> bool {
            CUDA_CHECK(cudaMemcpy(d_y_aos, h_y_aos, ny * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
            launch(); CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_out_aos, d_y_aos, ny * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
            return verify_aos(ny, h_out_aos, h_ref_aos);
        };
        CUDA_CHECK(cudaMemcpy(d_y_aos, h_y_aos, ny * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
        run_variant("aos_sorted", launch, verify);
    }

    // =========== SoA scatter ===========
    {
        auto kfn = soa_scatter_tbl[coarsen];
        CUDA_CHECK(cudaMemcpy(d_x_re, h_x_re, nx * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_x_im, h_x_im, nx * sizeof(double), cudaMemcpyHostToDevice));
        auto launch = [&]() { kfn<<<num_blocks, tpb>>>(nx, d_ymap, d_x_re, d_x_im, d_y_re, d_y_im); };
        auto verify = [&]() -> bool {
            CUDA_CHECK(cudaMemcpy(d_y_re, h_y_re, ny * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_y_im, h_y_im, ny * sizeof(double), cudaMemcpyHostToDevice));
            launch(); CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_out_re, d_y_re, ny * sizeof(double), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_out_im, d_y_im, ny * sizeof(double), cudaMemcpyDeviceToHost));
            return verify_soa(ny, h_out_re, h_out_im, h_ref_re, h_ref_im);
        };
        CUDA_CHECK(cudaMemcpy(d_y_re, h_y_re, ny * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_y_im, h_y_im, ny * sizeof(double), cudaMemcpyHostToDevice));
        run_variant("soa_scatter", launch, verify);
    }

    // =========== SoA sorted ===========
    {
        auto kfn = soa_sorted_tbl[coarsen];
        CUDA_CHECK(cudaMemcpy(d_x_re, h_x_re, nx * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_x_im, h_x_im, nx * sizeof(double), cudaMemcpyHostToDevice));
        auto launch = [&]() { kfn<<<num_blocks, tpb>>>(nx, d_ymap_s, d_xmap, d_x_re, d_x_im, d_y_re, d_y_im); };
        auto verify = [&]() -> bool {
            CUDA_CHECK(cudaMemcpy(d_y_re, h_y_re, ny * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_y_im, h_y_im, ny * sizeof(double), cudaMemcpyHostToDevice));
            launch(); CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_out_re, d_y_re, ny * sizeof(double), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_out_im, d_y_im, ny * sizeof(double), cudaMemcpyDeviceToHost));
            return verify_soa(ny, h_out_re, h_out_im, h_ref_re, h_ref_im);
        };
        CUDA_CHECK(cudaMemcpy(d_y_re, h_y_re, ny * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_y_im, h_y_im, ny * sizeof(double), cudaMemcpyHostToDevice));
        run_variant("soa_sorted", launch, verify);
    }

    // ---- Cleanup ----
    CUDA_CHECK(cudaEventDestroy(e0));
    CUDA_CHECK(cudaEventDestroy(e1));
    CUDA_CHECK(cudaFree(d_y_aos)); CUDA_CHECK(cudaFree(d_x_aos));
    CUDA_CHECK(cudaFree(d_y_re));  CUDA_CHECK(cudaFree(d_y_im));
    CUDA_CHECK(cudaFree(d_x_re));  CUDA_CHECK(cudaFree(d_x_im));
    CUDA_CHECK(cudaFree(d_ymap));  CUDA_CHECK(cudaFree(d_ymap_s));
    CUDA_CHECK(cudaFree(d_xmap));

    free(h_y_aos); free(h_x_aos); free(h_out_aos); free(h_ref_aos);
    free(h_y_re);  free(h_y_im);  free(h_x_re);    free(h_x_im);
    free(h_out_re); free(h_out_im); free(h_ref_re); free(h_ref_im);
}

/* ================================================================ */
/*  Helpers: sort, sample, binary I/O                               */
/* ================================================================ */

void sortWithIndices(int N, const int *input, int *sortedValues, int *sortedIndices)
{
    std::iota(sortedIndices, sortedIndices + N, 0);
    std::sort(sortedIndices, sortedIndices + N,
              [&input](int a, int b) { return input[a] < input[b]; });
    for (int i = 0; i < N; ++i)
        sortedValues[i] = input[sortedIndices[i]];
}

std::vector<int> uniformSample(int nx, int ny, std::mt19937 &rng)
{
    std::vector<int> pool(ny);
    std::iota(pool.begin(), pool.end(), 0);
    for (int i = 0; i < nx; ++i)
    {
        std::uniform_int_distribution<int> dist(i, ny - 1);
        std::swap(pool[i], pool[dist(rng)]);
    }
    return std::vector<int>(pool.begin(), pool.begin() + nx);
}

std::vector<int> normalSample(int nx, int ny, std::mt19937 &rng)
{
    const double mu = (ny - 1) / 2.0;
    const double sigma = ny / 6.0;
    std::normal_distribution<double> dist(mu, sigma);
    std::unordered_set<int> seen;
    std::vector<int> result;
    result.reserve(nx);
    while ((int)result.size() < nx)
    {
        int s = static_cast<int>(std::round(dist(rng)));
        if (s < 0 || s >= ny) continue;
        if (!seen.insert(s).second) continue;
        result.push_back(s);
    }
    return result;
}

int *readBinaryFile(const std::string &filepath, int &N)
{
    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + filepath);
    std::streamsize sz = file.tellg();
    if (sz % sizeof(int) != 0)
        throw std::runtime_error("File size not multiple of sizeof(int).");
    N = sz / sizeof(int);
    file.seekg(0, std::ios::beg);
    int *data = new int[N];
    if (!file.read(reinterpret_cast<char *>(data), sz))
    {
        delete[] data;
        throw std::runtime_error("Failed to read file: " + filepath);
    }
    return data;
}

/* ================================================================ */
/*  CSV header                                                       */
/* ================================================================ */
static const char *CSV_HEADER = "variant,experiment,ny,nx,tpb,coarsen,rep,time_ms\n";

/* ================================================================ */
/*  Main                                                             */
/* ================================================================ */

int main(void)
{
    constexpr int ITERS  = 100;
    constexpr int WARMUP = 5;
    constexpr int NUM_SAMPLES = 5;

    std::array<int, 7> tpb_values     = {32, 64, 128, 256, 384, 512, 1024};
    std::array<int, 4> coarsen_values = {1, 2, 4, 8};
    std::array<int, 5> ny_multiples   = {1, 2, 4, 8, 16};

    std::mt19937 rng(42);
    init_flush();

    // ============================================================
    //  QE real-world data
    // ============================================================
    struct QECase { const char *name; int ny; int nx; const char *binfile; };
    QECase qe_cases[] = {
        {"qe_nat05", 64001,  27609,  "./bin/vexx_k_gpu__dfftt__nl_nat05.bin"},
        {"qe_nat10", 120001, 55191,  "./bin/vexx_k_gpu__dfftt__nl_nat10.bin"},
        {"qe_nat20", 225001, 110273, "./bin/vexx_k_gpu__dfftt__nl_nat20.bin"},
    };

    FILE *csv = fopen("zaxpy_sweep_qe.csv", "w");
    fprintf(csv, "%s", CSV_HEADER);

    for (auto &qe : qe_cases)
    {
        printf("\n>>> %s  ny=%d nx=%d\n", qe.name, qe.ny, qe.nx);
        int N = 0;
        int *ymap = readBinaryFile(qe.binfile, N);
        assert(N == qe.nx);

        int *ymap_s = new int[qe.nx];
        int *xmap   = new int[qe.nx];
        sortWithIndices(qe.nx, ymap, ymap_s, xmap);

        for (int tpb : tpb_values)
            for (int cf : coarsen_values)
            {
                printf("  tpb=%d coarsen=%d\n", tpb, cf);
                profile_config(csv, qe.name, qe.ny, qe.nx,
                               ymap, ymap_s, xmap, tpb, cf, ITERS, WARMUP, 0);
                fflush(csv);
            }

        delete[] ymap; delete[] ymap_s; delete[] xmap;
    }
    fclose(csv);

    // ============================================================
    //  Synthetic: uniform
    // ============================================================
    int nx = 1024 * 1024;
    int *ymap_s = new int[nx];
    int *xmap   = new int[nx];

    csv = fopen("zaxpy_sweep_uniform.csv", "w");
    fprintf(csv, "%s", CSV_HEADER);

    for (int mult : ny_multiples)
    {
        int ny = mult * nx;
        printf("\n>>> uniform_mult%d  ny=%d nx=%d\n", mult, ny, nx);

        for (int s = 0; s < NUM_SAMPLES; s++)
        {
            auto ymap_vec = uniformSample(nx, ny, rng);
            sortWithIndices(nx, ymap_vec.data(), ymap_s, xmap);

            char ename[80];
            snprintf(ename, sizeof(ename), "uniform_mult%d_s%d", mult, s);

            for (int tpb : tpb_values)
                for (int cf : coarsen_values)
                {
                    printf("  %s tpb=%d coarsen=%d\n", ename, tpb, cf);
                    profile_config(csv, ename, ny, nx,
                                   ymap_vec.data(), ymap_s, xmap,
                                   tpb, cf, ITERS, WARMUP, 0);
                    fflush(csv);
                }
        }
    }
    fclose(csv);

    // ============================================================
    //  Synthetic: normal
    // ============================================================
    csv = fopen("zaxpy_sweep_normal.csv", "w");
    fprintf(csv, "%s", CSV_HEADER);

    for (int mult : ny_multiples)
    {
        int ny = mult * nx;
        printf("\n>>> normal_mult%d  ny=%d nx=%d\n", mult, ny, nx);

        for (int s = 0; s < NUM_SAMPLES; s++)
        {
            auto ymap_vec = normalSample(nx, ny, rng);
            sortWithIndices(nx, ymap_vec.data(), ymap_s, xmap);

            char ename[80];
            snprintf(ename, sizeof(ename), "normal_mult%d_s%d", mult, s);

            for (int tpb : tpb_values)
                for (int cf : coarsen_values)
                {
                    printf("  %s tpb=%d coarsen=%d\n", ename, tpb, cf);
                    profile_config(csv, ename, ny, nx,
                                   ymap_vec.data(), ymap_s, xmap,
                                   tpb, cf, ITERS, WARMUP, 0);
                    fflush(csv);
                }
        }
    }
    fclose(csv);

    delete[] ymap_s; delete[] xmap;
    if (d_flush) cudaFree(d_flush);

    printf("\nDone. CSVs: zaxpy_sweep_qe.csv, zaxpy_sweep_uniform.csv, zaxpy_sweep_normal.csv\n");
    return 0;
}
