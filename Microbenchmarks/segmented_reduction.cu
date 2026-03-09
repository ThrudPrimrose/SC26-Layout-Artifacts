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
// For the 81920×90 case (reduce along dim-0, i.e. each column → 1 value),
// elements of column j are strided: data[0*90+j], data[1*90+j], ...
// CUB needs contiguous virtual segments, so we build a transform iterator
// that remaps virtual linear index → physical strided index.
// ---------------------------------------------------------------------------
struct ColumnAccessOp {
    const float* data;
    int num_cols;   // 90
    int seg_len;    // 81920

    __host__ __device__ __forceinline__
    float operator()(int virtual_idx) const {
        int seg = virtual_idx / seg_len;   // which column
        int pos = virtual_idx % seg_len;   // which row
        return data[pos * num_cols + seg];
    }
};

// ---------------------------------------------------------------------------
// Build host-side offset array: [0, seg_len, 2*seg_len, ..., num_segs*seg_len]
// ---------------------------------------------------------------------------
static void make_offsets(int num_segs, int seg_len, int** d_offsets) {
    std::vector<int> h(num_segs + 1);
    for (int i = 0; i <= num_segs; ++i) h[i] = i * seg_len;
    CHECK_CUDA(cudaMalloc(d_offsets, h.size() * sizeof(int)));
    CHECK_CUDA(cudaMemcpy(*d_offsets, h.data(), h.size() * sizeof(int),
                           cudaMemcpyHostToDevice));
}

// ---------------------------------------------------------------------------
// Case A: 81920×90 row-major, reduce each column (strided) → 90 outputs
// ---------------------------------------------------------------------------
static void bench_column_reduce(const float* d_data, float* d_out,
                                int num_rows, int num_cols, int iters,
                                int warmup, std::vector<float>& times_ms) {
    const int num_segs = num_cols;  // 90
    const int seg_len  = num_rows;  // 81920

    int* d_offsets;
    make_offsets(num_segs, seg_len, &d_offsets);

    // Build the virtual iterator: virtual_idx → data[row * num_cols + col]
    cub::CountingInputIterator<int> counting(0);
    ColumnAccessOp op{d_data, num_cols, seg_len};
    cub::TransformInputIterator<float, ColumnAccessOp,
                                cub::CountingInputIterator<int>>
        col_iter(counting, op);

    // Query temp storage
    void*  d_temp = nullptr;
    size_t temp_bytes = 0;
    CHECK_CUDA(cub::DeviceSegmentedReduce::Max(
        d_temp, temp_bytes, col_iter, d_out,
        num_segs, d_offsets, d_offsets + 1));
    CHECK_CUDA(cudaMalloc(&d_temp, temp_bytes));

    cudaEvent_t t0, t1;
    CHECK_CUDA(cudaEventCreate(&t0));
    CHECK_CUDA(cudaEventCreate(&t1));

    for (int i = 0; i < warmup + iters; ++i) {
        CHECK_CUDA(cudaEventRecord(t0));
        CHECK_CUDA(cub::DeviceSegmentedReduce::Max(
            d_temp, temp_bytes, col_iter, d_out,
            num_segs, d_offsets, d_offsets + 1));
        CHECK_CUDA(cudaEventRecord(t1));
        CHECK_CUDA(cudaEventSynchronize(t1));
        if (i >= warmup) {
            float ms;
            CHECK_CUDA(cudaEventElapsedTime(&ms, t0, t1));
            times_ms.push_back(ms);
        }
    }

    CHECK_CUDA(cudaEventDestroy(t0));
    CHECK_CUDA(cudaEventDestroy(t1));
    CHECK_CUDA(cudaFree(d_temp));
    CHECK_CUDA(cudaFree(d_offsets));
}

// ---------------------------------------------------------------------------
// Case B: 90×81920 row-major, reduce each row (contiguous) → 90 outputs
// ---------------------------------------------------------------------------
static void bench_row_reduce(const float* d_data, float* d_out,
                             int num_rows, int num_cols, int iters,
                             int warmup, std::vector<float>& times_ms) {
    const int num_segs = num_rows;  // 90
    const int seg_len  = num_cols;  // 81920

    int* d_offsets;
    make_offsets(num_segs, seg_len, &d_offsets);

    // Direct pointer — rows are already contiguous
    void*  d_temp = nullptr;
    size_t temp_bytes = 0;
    CHECK_CUDA(cub::DeviceSegmentedReduce::Max(
        d_temp, temp_bytes, d_data, d_out,
        num_segs, d_offsets, d_offsets + 1));
    CHECK_CUDA(cudaMalloc(&d_temp, temp_bytes));

    cudaEvent_t t0, t1;
    CHECK_CUDA(cudaEventCreate(&t0));
    CHECK_CUDA(cudaEventCreate(&t1));

    for (int i = 0; i < warmup + iters; ++i) {
        CHECK_CUDA(cudaEventRecord(t0));
        CHECK_CUDA(cub::DeviceSegmentedReduce::Max(
            d_temp, temp_bytes, d_data, d_out,
            num_segs, d_offsets, d_offsets + 1));
        CHECK_CUDA(cudaEventRecord(t1));
        CHECK_CUDA(cudaEventSynchronize(t1));
        if (i >= warmup) {
            float ms;
            CHECK_CUDA(cudaEventElapsedTime(&ms, t0, t1));
            times_ms.push_back(ms);
        }
    }

    CHECK_CUDA(cudaEventDestroy(t0));
    CHECK_CUDA(cudaEventDestroy(t1));
    CHECK_CUDA(cudaFree(d_temp));
    CHECK_CUDA(cudaFree(d_offsets));
}

// ---------------------------------------------------------------------------
// Statistics helper
// ---------------------------------------------------------------------------
struct Stats {
    float mean, median, stddev, min, max;
};

static Stats compute_stats(std::vector<float>& v) {
    Stats s{};
    std::sort(v.begin(), v.end());
    s.min = v.front();
    s.max = v.back();
    s.median = v[v.size() / 2];
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    s.mean = static_cast<float>(sum / v.size());
    double sq = 0;
    for (auto x : v) sq += (x - s.mean) * (x - s.mean);
    s.stddev = static_cast<float>(std::sqrt(sq / v.size()));
    return s;
}

// ---------------------------------------------------------------------------
// Correctness check: compare both reductions against each other
// ---------------------------------------------------------------------------
static void verify(const float* d_data_col, const float* d_data_row,
                   int N, int S) {
    // Both arrays contain the same logical data (N total elements, S segments
    // of length N/S).  Reduce on host and compare GPU outputs.
    std::vector<float> h(N);

    // Verify column-reduce (81920×90)
    CHECK_CUDA(cudaMemcpy(h.data(), d_data_col, N * sizeof(float),
                           cudaMemcpyDeviceToHost));
    std::vector<float> ref_col(S, -FLT_MAX);
    for (int r = 0; r < 81920; ++r)
        for (int c = 0; c < 90; ++c)
            ref_col[c] = fmaxf(ref_col[c], h[r * 90 + c]);

    // Verify row-reduce (90×81920)
    CHECK_CUDA(cudaMemcpy(h.data(), d_data_row, N * sizeof(float),
                           cudaMemcpyDeviceToHost));
    std::vector<float> ref_row(S, -FLT_MAX);
    for (int r = 0; r < 90; ++r)
        for (int c = 0; c < 81920; ++c)
            ref_row[r] = fmaxf(ref_row[r], h[r * 81920 + c]);

    // Check GPU results
    float* d_out_col;
    float* d_out_row;
    CHECK_CUDA(cudaMalloc(&d_out_col, S * sizeof(float)));
    CHECK_CUDA(cudaMalloc(&d_out_row, S * sizeof(float)));

    // Run once each
    int* d_off_col; make_offsets(90, 81920, &d_off_col);
    int* d_off_row; make_offsets(90, 81920, &d_off_row);

    { // column reduce
        cub::CountingInputIterator<int> cnt(0);
        ColumnAccessOp op{d_data_col, 90, 81920};
        cub::TransformInputIterator<float, ColumnAccessOp,
                                    cub::CountingInputIterator<int>> it(cnt, op);
        void* tmp = nullptr; size_t nb = 0;
        cub::DeviceSegmentedReduce::Max(tmp, nb, it, d_out_col, 90,
                                        d_off_col, d_off_col + 1);
        cudaMalloc(&tmp, nb);
        cub::DeviceSegmentedReduce::Max(tmp, nb, it, d_out_col, 90,
                                        d_off_col, d_off_col + 1);
        cudaDeviceSynchronize();
        cudaFree(tmp);
    }
    { // row reduce
        void* tmp = nullptr; size_t nb = 0;
        cub::DeviceSegmentedReduce::Max(tmp, nb, d_data_row, d_out_row, 90,
                                        d_off_row, d_off_row + 1);
        cudaMalloc(&tmp, nb);
        cub::DeviceSegmentedReduce::Max(tmp, nb, d_data_row, d_out_row, 90,
                                        d_off_row, d_off_row + 1);
        cudaDeviceSynchronize();
        cudaFree(tmp);
    }

    std::vector<float> gpu_col(S), gpu_row(S);
    cudaMemcpy(gpu_col.data(), d_out_col, S * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(gpu_row.data(), d_out_row, S * sizeof(float), cudaMemcpyDeviceToHost);

    int errs = 0;
    for (int i = 0; i < S; ++i) {
        if (gpu_col[i] != ref_col[i]) { errs++; if (errs <= 3) printf("  col[%d]: gpu=%.6f ref=%.6f\n", i, gpu_col[i], ref_col[i]); }
        if (gpu_row[i] != ref_row[i]) { errs++; if (errs <= 3) printf("  row[%d]: gpu=%.6f ref=%.6f\n", i, gpu_row[i], ref_row[i]); }
    }
    printf("Verification: %s (%d mismatches)\n\n", errs == 0 ? "PASS" : "FAIL", errs);

    cudaFree(d_out_col); cudaFree(d_out_row);
    cudaFree(d_off_col); cudaFree(d_off_row);
}

// ---------------------------------------------------------------------------
int main() {
    constexpr int ROWS_A = 81920, COLS_A = 90;   // Case A: column reduce
    constexpr int ROWS_B = 90,    COLS_B = 81920; // Case B: row reduce
    constexpr int N = ROWS_A * COLS_A;            // same total elements
    constexpr int NUM_SEGMENTS = 90;
    constexpr int TOTAL_ITERS = 100;
    constexpr int WARMUP      = 5;
    constexpr int MEASURED     = TOTAL_ITERS - WARMUP;

    // Print GPU info
    cudaDeviceProp prop;
    CHECK_CUDA(cudaGetDeviceProperties(&prop, 0));
    printf("GPU: %s  (SM %d.%d, %d SMs)\n\n", prop.name, prop.major,
           prop.minor, prop.multiProcessorCount);

    // Allocate and init with random data
    std::vector<float> h_data(N);
    srand(42);
    for (auto& v : h_data) v = static_cast<float>(rand()) / RAND_MAX;

    float *d_data_A, *d_data_B;
    CHECK_CUDA(cudaMalloc(&d_data_A, N * sizeof(float)));
    CHECK_CUDA(cudaMalloc(&d_data_B, N * sizeof(float)));
    CHECK_CUDA(cudaMemcpy(d_data_A, h_data.data(), N * sizeof(float),
                           cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_data_B, h_data.data(), N * sizeof(float),
                           cudaMemcpyHostToDevice));

    float *d_out_A, *d_out_B;
    CHECK_CUDA(cudaMalloc(&d_out_A, NUM_SEGMENTS * sizeof(float)));
    CHECK_CUDA(cudaMalloc(&d_out_B, NUM_SEGMENTS * sizeof(float)));

    // Verify correctness first
    printf("--- Correctness ---\n");
    verify(d_data_A, d_data_B, N, NUM_SEGMENTS);

    // Benchmark
    size_t total_bytes = (size_t)N * sizeof(float);
    printf("Array size      : %d × %d = %d elements (%.2f MB)\n",
           ROWS_A, COLS_A, N, total_bytes / (1024.0 * 1024.0));
    printf("Segments        : %d  (each %d elements)\n", NUM_SEGMENTS,
           N / NUM_SEGMENTS);
    printf("Iterations      : %d total, %d warmup, %d measured\n\n",
           TOTAL_ITERS, WARMUP, MEASURED);

    std::vector<float> times_A, times_B;

    printf("Running Case A: 81920×90 column-reduce (strided via TransformIterator)...\n");
    bench_column_reduce(d_data_A, d_out_A, ROWS_A, COLS_A,
                        MEASURED, WARMUP, times_A);

    printf("Running Case B: 90×81920 row-reduce (contiguous segments)...\n");
    bench_row_reduce(d_data_B, d_out_B, ROWS_B, COLS_B,
                     MEASURED, WARMUP, times_B);

    // Report
    auto sA = compute_stats(times_A);
    auto sB = compute_stats(times_B);

    auto bw = [&](float ms) {
        return (total_bytes / 1e9) / (ms / 1e3);  // GB/s
    };

    printf("\n========================================================\n");
    printf("  Case A — 81920×90  column-reduce  (strided access)\n");
    printf("--------------------------------------------------------\n");
    printf("  Mean   : %8.4f ms   (%6.2f GB/s)\n", sA.mean, bw(sA.mean));
    printf("  Median : %8.4f ms   (%6.2f GB/s)\n", sA.median, bw(sA.median));
    printf("  Stddev : %8.4f ms\n", sA.stddev);
    printf("  Min    : %8.4f ms   (%6.2f GB/s)\n", sA.min, bw(sA.min));
    printf("  Max    : %8.4f ms   (%6.2f GB/s)\n", sA.max, bw(sA.max));
    printf("========================================================\n");
    printf("  Case B — 90×81920  row-reduce     (contiguous access)\n");
    printf("--------------------------------------------------------\n");
    printf("  Mean   : %8.4f ms   (%6.2f GB/s)\n", sB.mean, bw(sB.mean));
    printf("  Median : %8.4f ms   (%6.2f GB/s)\n", sB.median, bw(sB.median));
    printf("  Stddev : %8.4f ms\n", sB.stddev);
    printf("  Min    : %8.4f ms   (%6.2f GB/s)\n", sB.min, bw(sB.min));
    printf("  Max    : %8.4f ms   (%6.2f GB/s)\n", sB.max, bw(sB.max));
    printf("========================================================\n");
    printf("\n  Speedup (contiguous / strided) : %.2fx\n\n",
           sA.median / sB.median);

    CHECK_CUDA(cudaFree(d_data_A));
    CHECK_CUDA(cudaFree(d_data_B));
    CHECK_CUDA(cudaFree(d_out_A));
    CHECK_CUDA(cudaFree(d_out_B));

    return 0;
}
