#include <cub/cub.cuh>
#include <cuda_runtime.h>
#include <cstdio>
#include <cfloat>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

#define CHECK_CUDA(call)                                                       \
    do {                                                                        \
        cudaError_t err = (call);                                               \
        if (err != cudaSuccess) {                                               \
            fprintf(stderr, "CUDA error at %s:%d — %s\n", __FILE__, __LINE__,  \
                    cudaGetErrorString(err));                                    \
            exit(1);                                                            \
        }                                                                       \
    } while (0)

// ---------------------------------------------------------------------------
//  Constants
// ---------------------------------------------------------------------------
constexpr int NROWS   = 81920;
constexpr int NCOLS   = 90;
constexpr int N       = NROWS * NCOLS;

// 3 warps = 96 threads.  Threads 0..89 each own one column; 90..95 idle.
constexpr int BLOCK_THREADS  = 96;
// How many rows each block chews through.  Tune this.
constexpr int ROWS_PER_BLOCK = 1024;
constexpr int NUM_BLOCKS     = (NROWS + ROWS_PER_BLOCK - 1) / ROWS_PER_BLOCK; // 80

// ---------------------------------------------------------------------------
//  Pass 1 — partial column-max
//  Grid : NUM_BLOCKS blocks  (each 96 threads)
//  Each block handles rows [block_start, block_end) for ALL 90 columns.
//  Thread tid (tid < NCOLS) walks its column with stride NCOLS — coalesced
//  because adjacent threads read adjacent addresses within the same row.
//
//  Output: partial[blockIdx.x * NCOLS + tid] = max over the block's row range
// ---------------------------------------------------------------------------
__global__ void column_max_partial(const float* __restrict__ data,
                                   float*       __restrict__ partial,
                                   int nrows, int ncols)
{
    const int tid   = threadIdx.x;
    const int bid   = blockIdx.x;
    const int r0    = bid * ROWS_PER_BLOCK;
    const int r1    = min(r0 + ROWS_PER_BLOCK, nrows);

    if (tid >= ncols) return;           // threads 90..95 bail out

    float mx = -FLT_MAX;
    // Walk down the column; consecutive threads hit consecutive floats
    for (int r = r0; r < r1; ++r) {
        mx = fmaxf(mx, data[r * ncols + tid]);
    }
    partial[bid * ncols + tid] = mx;
}

// ---------------------------------------------------------------------------
//  Pass 2 — final reduction across NUM_BLOCKS partial maxes per column
//  One block of 96 threads, each thread reduces its column across all blocks.
// ---------------------------------------------------------------------------
__global__ void column_max_final(const float* __restrict__ partial,
                                 float*       __restrict__ out,
                                 int num_blocks, int ncols)
{
    const int tid = threadIdx.x;
    if (tid >= ncols) return;

    float mx = -FLT_MAX;
    for (int b = 0; b < num_blocks; ++b) {
        mx = fmaxf(mx, partial[b * ncols + tid]);
    }
    out[tid] = mx;
}

// ---------------------------------------------------------------------------
//  CUB baselines (same as before, inlined for self-contained file)
// ---------------------------------------------------------------------------
struct ColumnAccessOp {
    const float* data;
    int num_cols, seg_len;
    __host__ __device__ __forceinline__
    float operator()(int vi) const {
        return data[(vi % seg_len) * num_cols + vi / seg_len];
    }
};

static void make_offsets(int num_segs, int seg_len, int** d_off) {
    std::vector<int> h(num_segs + 1);
    for (int i = 0; i <= num_segs; ++i) h[i] = i * seg_len;
    CHECK_CUDA(cudaMalloc(d_off, h.size() * sizeof(int)));
    CHECK_CUDA(cudaMemcpy(*d_off, h.data(), h.size() * sizeof(int),
                           cudaMemcpyHostToDevice));
}

// ---------------------------------------------------------------------------
//  Benchmark harness
// ---------------------------------------------------------------------------
struct Stats { float mean, median, stddev, mn, mx; };

static Stats compute_stats(std::vector<float>& v) {
    std::sort(v.begin(), v.end());
    Stats s;
    s.mn = v.front(); s.mx = v.back(); s.median = v[v.size()/2];
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    s.mean = (float)(sum / v.size());
    double sq = 0;
    for (auto x : v) sq += (x - s.mean) * (x - s.mean);
    s.stddev = (float)std::sqrt(sq / v.size());
    return s;
}

static float bw_gbs(float ms) {
    return ((size_t)N * sizeof(float) / 1e9) / (ms / 1e3);
}

static void print_stats(const char* label, Stats& s) {
    printf("  %-50s\n", label);
    printf("  %-8s: %8.4f ms   (%7.2f GB/s)\n", "Mean",   s.mean,   bw_gbs(s.mean));
    printf("  %-8s: %8.4f ms   (%7.2f GB/s)\n", "Median", s.median, bw_gbs(s.median));
    printf("  %-8s: %8.4f ms\n",                 "Stddev", s.stddev);
    printf("  %-8s: %8.4f ms   (%7.2f GB/s)\n", "Min",    s.mn,     bw_gbs(s.mn));
    printf("  %-8s: %8.4f ms   (%7.2f GB/s)\n", "Max",    s.mx,     bw_gbs(s.mx));
}

// ---------------------------------------------------------------------------
int main() {
    constexpr int TOTAL = 100, WARMUP = 5, MEASURED = TOTAL - WARMUP;

    cudaDeviceProp prop;
    CHECK_CUDA(cudaGetDeviceProperties(&prop, 0));
    printf("GPU: %s  (SM %d.%d, %d SMs)\n", prop.name, prop.major,
           prop.minor, prop.multiProcessorCount);
    printf("Array: %d × %d = %d elements (%.2f MB)\n",
           NROWS, NCOLS, N, N * sizeof(float) / (1024.0 * 1024.0));
    printf("Custom kernel: %d threads/block, %d rows/block, %d blocks\n\n",
           BLOCK_THREADS, ROWS_PER_BLOCK, NUM_BLOCKS);

    // ---- Allocate & init ----
    std::vector<float> h_data(N);
    srand(42);
    for (auto& v : h_data) v = (float)rand() / RAND_MAX;

    float *d_data;
    CHECK_CUDA(cudaMalloc(&d_data, N * sizeof(float)));
    CHECK_CUDA(cudaMemcpy(d_data, h_data.data(), N * sizeof(float),
                           cudaMemcpyHostToDevice));

    float *d_out;
    CHECK_CUDA(cudaMalloc(&d_out, NCOLS * sizeof(float)));

    float *d_partial;
    CHECK_CUDA(cudaMalloc(&d_partial, NUM_BLOCKS * NCOLS * sizeof(float)));

    // ---- Verify custom kernel ----
    column_max_partial<<<NUM_BLOCKS, BLOCK_THREADS>>>(d_data, d_partial, NROWS, NCOLS);
    column_max_final<<<1, BLOCK_THREADS>>>(d_partial, d_out, NUM_BLOCKS, NCOLS);
    CHECK_CUDA(cudaDeviceSynchronize());

    std::vector<float> gpu_custom(NCOLS);
    CHECK_CUDA(cudaMemcpy(gpu_custom.data(), d_out, NCOLS * sizeof(float),
                           cudaMemcpyDeviceToHost));

    std::vector<float> ref(NCOLS, -FLT_MAX);
    for (int r = 0; r < NROWS; ++r)
        for (int c = 0; c < NCOLS; ++c)
            ref[c] = fmaxf(ref[c], h_data[r * NCOLS + c]);

    int errs = 0;
    for (int i = 0; i < NCOLS; ++i)
        if (gpu_custom[i] != ref[i]) errs++;
    printf("Custom kernel verification: %s (%d mismatches)\n\n",
           errs ? "FAIL" : "PASS", errs);

    // ==================================================================
    //  Benchmark 1: CUB strided (TransformIterator)
    // ==================================================================
    {
        int* d_off; make_offsets(NCOLS, NROWS, &d_off);
        cub::CountingInputIterator<int> cnt(0);
        ColumnAccessOp op{d_data, NCOLS, NROWS};
        cub::TransformInputIterator<float, ColumnAccessOp,
                                    cub::CountingInputIterator<int>> it(cnt, op);
        void* tmp = nullptr; size_t nb = 0;
        cub::DeviceSegmentedReduce::Max(tmp, nb, it, d_out, NCOLS,
                                        d_off, d_off + 1);
        cudaMalloc(&tmp, nb);

        cudaEvent_t t0, t1;
        cudaEventCreate(&t0); cudaEventCreate(&t1);
        std::vector<float> times;

        for (int i = 0; i < TOTAL; ++i) {
            cudaEventRecord(t0);
            cub::DeviceSegmentedReduce::Max(tmp, nb, it, d_out, NCOLS,
                                            d_off, d_off + 1);
            cudaEventRecord(t1);
            cudaEventSynchronize(t1);
            if (i >= WARMUP) { float ms; cudaEventElapsedTime(&ms, t0, t1); times.push_back(ms); }
        }
        auto s = compute_stats(times);
        printf("========================================================\n");
        print_stats("CUB DeviceSegmentedReduce (strided TransformIter)", s);
        printf("========================================================\n");

        cudaFree(tmp); cudaFree(d_off);
        cudaEventDestroy(t0); cudaEventDestroy(t1);
    }

    // ==================================================================
    //  Benchmark 2: CUB contiguous (90×81920 row-reduce, as reference)
    // ==================================================================
    {
        // Lay out same data as 90×81920 row-major for a fair comparison
        float* d_data_row;
        CHECK_CUDA(cudaMalloc(&d_data_row, N * sizeof(float)));

        // Transpose on host: row-major 81920×90 → 90×81920
        std::vector<float> h_row(N);
        for (int r = 0; r < NROWS; ++r)
            for (int c = 0; c < NCOLS; ++c)
                h_row[c * NROWS + r] = h_data[r * NCOLS + c];
        CHECK_CUDA(cudaMemcpy(d_data_row, h_row.data(), N * sizeof(float),
                               cudaMemcpyHostToDevice));

        float* d_out_row;
        CHECK_CUDA(cudaMalloc(&d_out_row, NCOLS * sizeof(float)));

        int* d_off; make_offsets(NCOLS, NROWS, &d_off);
        void* tmp = nullptr; size_t nb = 0;
        cub::DeviceSegmentedReduce::Max(tmp, nb, d_data_row, d_out_row, NCOLS,
                                        d_off, d_off + 1);
        cudaMalloc(&tmp, nb);

        cudaEvent_t t0, t1;
        cudaEventCreate(&t0); cudaEventCreate(&t1);
        std::vector<float> times;

        for (int i = 0; i < TOTAL; ++i) {
            cudaEventRecord(t0);
            cub::DeviceSegmentedReduce::Max(tmp, nb, d_data_row, d_out_row, NCOLS,
                                            d_off, d_off + 1);
            cudaEventRecord(t1);
            cudaEventSynchronize(t1);
            if (i >= WARMUP) { float ms; cudaEventElapsedTime(&ms, t0, t1); times.push_back(ms); }
        }
        auto s = compute_stats(times);
        printf("========================================================\n");
        print_stats("CUB DeviceSegmentedReduce (contiguous, 90×81920)", s);
        printf("========================================================\n");

        cudaFree(tmp); cudaFree(d_off); cudaFree(d_data_row); cudaFree(d_out_row);
        cudaEventDestroy(t0); cudaEventDestroy(t1);
    }

    // ==================================================================
    //  Benchmark 3: Custom coalesced column-reduce (two-pass)
    // ==================================================================
    {
        cudaEvent_t t0, t1;
        cudaEventCreate(&t0); cudaEventCreate(&t1);
        std::vector<float> times;

        for (int i = 0; i < TOTAL; ++i) {
            cudaEventRecord(t0);
            column_max_partial<<<NUM_BLOCKS, BLOCK_THREADS>>>(
                d_data, d_partial, NROWS, NCOLS);
            column_max_final<<<1, BLOCK_THREADS>>>(
                d_partial, d_out, NUM_BLOCKS, NCOLS);
            cudaEventRecord(t1);
            cudaEventSynchronize(t1);
            if (i >= WARMUP) { float ms; cudaEventElapsedTime(&ms, t0, t1); times.push_back(ms); }
        }
        auto s = compute_stats(times);
        printf("========================================================\n");
        print_stats("Custom 3-warp coalesced column-reduce (two-pass)", s);
        printf("========================================================\n");

        cudaEventDestroy(t0); cudaEventDestroy(t1);
    }

    // ==================================================================
    //  Summary
    // ==================================================================
    printf("\n");

    CHECK_CUDA(cudaFree(d_data));
    CHECK_CUDA(cudaFree(d_out));
    CHECK_CUDA(cudaFree(d_partial));

    return 0;
}
