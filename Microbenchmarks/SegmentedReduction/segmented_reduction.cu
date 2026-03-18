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

constexpr int NROWS = 81920;
constexpr int NCOLS = 90;
constexpr int N     = NROWS * NCOLS;
constexpr int BT    = 96;  // 3 warps

// ===================================================================
//  NEW — Single-pass column-max.
//
//  Template TILE_C: columns handled per block.
//  96 / TILE_C workers share each column, striding through ALL rows.
//  Grid = ceil(NCOLS / TILE_C) blocks → output is final (no pass 2).
//
//  Coalescing:  within a warp, consecutive threads (tid, tid+1, …)
//  map to consecutive columns  (tid % TILE_C increments), so loads
//  within each row-group hit adjacent 4-byte addresses.
// ===================================================================
template <int TILE_C>
__global__ void col_max(const float* __restrict__ data,
                        float*       __restrict__ out,
                        int nrows, int ncols)
{
    static_assert(BT % TILE_C == 0, "TILE_C must divide block size");
    constexpr int WORKERS = BT / TILE_C;

    const int col  = blockIdx.x * TILE_C + threadIdx.x % TILE_C;
    const int rank = threadIdx.x / TILE_C;

    float mx = -FLT_MAX;
    if (col < ncols)
        for (int r = rank; r < nrows; r += WORKERS)
            mx = fmaxf(mx, data[r * ncols + col]);

    __shared__ float smem[BT];
    smem[threadIdx.x] = mx;
    __syncthreads();

    #pragma unroll
    for (int s = WORKERS >> 1; s > 0; s >>= 1) {
        if (rank < s)
            smem[threadIdx.x] = fmaxf(smem[threadIdx.x],
                                       smem[threadIdx.x + s * TILE_C]);
        __syncthreads();
    }

    if (rank == 0 && col < ncols)
        out[col] = smem[threadIdx.x % TILE_C];
}

// ===================================================================
//  V1 two-pass reference (from previous benchmark)
// ===================================================================
template <int RPB>
__global__ void v1_partial(const float* __restrict__ data,
                           float* __restrict__ part,
                           int nrows, int ncols)
{
    const int tid = threadIdx.x, bid = blockIdx.x;
    if (tid >= ncols) return;
    const int r0 = bid * RPB, r1 = min(r0 + RPB, nrows);
    float mx = -FLT_MAX;
    for (int r = r0; r < r1; ++r) mx = fmaxf(mx, data[r * ncols + tid]);
    part[bid * ncols + tid] = mx;
}

__global__ void v1_final(const float* __restrict__ part,
                         float* __restrict__ out,
                         int nblocks, int ncols)
{
    const int tid = threadIdx.x;
    if (tid >= ncols) return;
    float mx = -FLT_MAX;
    for (int b = 0; b < nblocks; ++b) mx = fmaxf(mx, part[b * ncols + tid]);
    out[tid] = mx;
}

// ===================================================================
//  CUB strided accessor
// ===================================================================
struct ColumnAccessOp {
    const float* data; int nc, nr;
    __host__ __device__ __forceinline__
    float operator()(int vi) const { return data[(vi % nr) * nc + vi / nr]; }
};

static void make_offsets(int nsegs, int seglen, int** d_off) {
    std::vector<int> h(nsegs + 1);
    for (int i = 0; i <= nsegs; ++i) h[i] = i * seglen;
    CHECK_CUDA(cudaMalloc(d_off, h.size() * sizeof(int)));
    CHECK_CUDA(cudaMemcpy(*d_off, h.data(), h.size() * sizeof(int),
                           cudaMemcpyHostToDevice));
}

// ===================================================================
//  Bench helpers
// ===================================================================
struct Stats { float mean, median, sd, mn, mx; };

static Stats compute(std::vector<float>& v) {
    std::sort(v.begin(), v.end());
    Stats s{};
    s.mn = v.front(); s.mx = v.back(); s.median = v[v.size()/2];
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    s.mean = (float)(sum / v.size());
    double sq = 0; for (auto x : v) sq += (x - s.mean)*(x - s.mean);
    s.sd = (float)std::sqrt(sq / v.size());
    return s;
}

static float gbps(float ms) { return (N * 4.0 / 1e9) / (ms / 1e3); }

static void report(const char* tag, Stats& s) {
    printf("  %-55s\n", tag);
    printf("    Mean %8.4f ms  (%7.2f GB/s)  Median %8.4f ms  (%7.2f GB/s)\n",
           s.mean, gbps(s.mean), s.median, gbps(s.median));
    printf("    Std  %8.4f ms   Min %8.4f ms   Max %8.4f ms\n",
           s.sd, s.mn, s.mx);
}

template <typename Fn>
static Stats bench(Fn fn, int total = 100, int warm = 5) {
    cudaEvent_t t0, t1;
    cudaEventCreate(&t0); cudaEventCreate(&t1);
    std::vector<float> times;
    for (int i = 0; i < total; ++i) {
        cudaEventRecord(t0); fn(); cudaEventRecord(t1);
        cudaEventSynchronize(t1);
        if (i >= warm) { float ms; cudaEventElapsedTime(&ms, t0, t1); times.push_back(ms); }
    }
    cudaEventDestroy(t0); cudaEventDestroy(t1);
    return compute(times);
}

static bool verify(float* d_out, const std::vector<float>& ref) {
    std::vector<float> g(ref.size());
    cudaMemcpy(g.data(), d_out, ref.size()*4, cudaMemcpyDeviceToHost);
    for (size_t i = 0; i < ref.size(); ++i) if (g[i] != ref[i]) return false;
    return true;
}

// ===================================================================
int main() {
    cudaDeviceProp prop;
    CHECK_CUDA(cudaGetDeviceProperties(&prop, 0));
    printf("GPU: %s  (%d SMs)\n", prop.name, prop.multiProcessorCount);
    printf("Array: %d×%d = %d elems (%.2f MB), block = %d threads (3 warps)\n\n",
           NROWS, NCOLS, N, N*4.0/(1<<20), BT);

    std::vector<float> h(N);
    srand(42);
    for (auto& v : h) v = (float)rand() / RAND_MAX;

    float *d_data, *d_out;
    CHECK_CUDA(cudaMalloc(&d_data, N * 4));
    CHECK_CUDA(cudaMemcpy(d_data, h.data(), N * 4, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMalloc(&d_out, NCOLS * 4));

    // Host reference
    std::vector<float> ref(NCOLS, -FLT_MAX);
    for (int r = 0; r < NROWS; ++r)
        for (int c = 0; c < NCOLS; ++c)
            ref[c] = fmaxf(ref[c], h[r * NCOLS + c]);

    printf("==========================================================\n");

    // ---- 1. CUB strided ----
    {
        int* d_off; make_offsets(NCOLS, NROWS, &d_off);
        cub::CountingInputIterator<int> cnt(0);
        ColumnAccessOp op{d_data, NCOLS, NROWS};
        cub::TransformInputIterator<float, ColumnAccessOp,
                                    cub::CountingInputIterator<int>> it(cnt, op);
        void* tmp = nullptr; size_t nb = 0;
        cub::DeviceSegmentedReduce::Max(tmp, nb, it, d_out, NCOLS, d_off, d_off+1);
        cudaMalloc(&tmp, nb);

        auto s = bench([&]{ cub::DeviceSegmentedReduce::Max(tmp, nb, it, d_out, NCOLS, d_off, d_off+1); });
        report("CUB strided (TransformIterator)", s);
        printf("    Verify: %s\n", verify(d_out, ref) ? "PASS" : "FAIL");
        cudaFree(tmp); cudaFree(d_off);
    }
    printf("----------------------------------------------------------\n");

    // ---- 2. CUB contiguous (transposed layout) ----
    {
        float* d_row;
        CHECK_CUDA(cudaMalloc(&d_row, N * 4));
        std::vector<float> ht(N);
        for (int r = 0; r < NROWS; ++r)
            for (int c = 0; c < NCOLS; ++c)
                ht[c * NROWS + r] = h[r * NCOLS + c];
        cudaMemcpy(d_row, ht.data(), N * 4, cudaMemcpyHostToDevice);

        int* d_off; make_offsets(NCOLS, NROWS, &d_off);
        void* tmp = nullptr; size_t nb = 0;
        cub::DeviceSegmentedReduce::Max(tmp, nb, d_row, d_out, NCOLS, d_off, d_off+1);
        cudaMalloc(&tmp, nb);

        auto s = bench([&]{ cub::DeviceSegmentedReduce::Max(tmp, nb, d_row, d_out, NCOLS, d_off, d_off+1); });
        report("CUB contiguous (transposed 90×81920)", s);

        std::vector<float> ref_row(NCOLS, -FLT_MAX);
        for (int c = 0; c < NCOLS; ++c)
            for (int r = 0; r < NROWS; ++r)
                ref_row[c] = fmaxf(ref_row[c], ht[c * NROWS + r]);
        printf("    Verify: %s\n", verify(d_out, ref_row) ? "PASS" : "FAIL");
        cudaFree(tmp); cudaFree(d_off); cudaFree(d_row);
    }
    printf("----------------------------------------------------------\n");

    // ---- 3. Custom single-pass at various TILE_C ----
    auto run = [&](auto tag, const char* label) {
        constexpr int TC = decltype(tag)::value;
        const int grid = (NCOLS + TC - 1) / TC;

        col_max<TC><<<grid, BT>>>(d_data, d_out, NROWS, NCOLS);
        cudaDeviceSynchronize();
        bool ok = verify(d_out, ref);

        auto s = bench([&]{ col_max<TC><<<grid, BT>>>(d_data, d_out, NROWS, NCOLS); });
        report(label, s);
        printf("    Grid: %d blocks, %d workers/col — Verify: %s\n",
               grid, BT / TC, ok ? "PASS" : "FAIL");
    };

    run(std::integral_constant<int,2>{},  "Custom TILE_C=2  (48 workers/col, 45 blks)");
    printf("----------------------------------------------------------\n");
    run(std::integral_constant<int,3>{},  "Custom TILE_C=3  (32 workers/col, 30 blks)");
    printf("----------------------------------------------------------\n");
    run(std::integral_constant<int,6>{},  "Custom TILE_C=6  (16 workers/col, 15 blks)");
    printf("----------------------------------------------------------\n");
    run(std::integral_constant<int,12>{}, "Custom TILE_C=12 ( 8 workers/col,  8 blks)");
    printf("----------------------------------------------------------\n");

    // ---- 4. Two-pass v1 (thread-per-column, 80+1 blocks) ----
    {
        constexpr int RPB = 1024;
        constexpr int NB  = (NROWS + RPB - 1) / RPB;  // 80
        float* d_part;
        CHECK_CUDA(cudaMalloc(&d_part, NB * NCOLS * 4));

        v1_partial<RPB><<<NB, BT>>>(d_data, d_part, NROWS, NCOLS);
        v1_final<<<1, BT>>>(d_part, d_out, NB, NCOLS);
        cudaDeviceSynchronize();
        bool ok = verify(d_out, ref);

        auto s = bench([&]{
            v1_partial<RPB><<<NB, BT>>>(d_data, d_part, NROWS, NCOLS);
            v1_final<<<1, BT>>>(d_part, d_out, NB, NCOLS);
        });
        report("Two-pass v1 (80+1 blks, 1024 rows/blk, thread-per-col)", s);
        printf("    Partial: %d×%d — Verify: %s\n", NB, NCOLS, ok ? "PASS" : "FAIL");
        cudaFree(d_part);
    }

    printf("==========================================================\n");
    cudaFree(d_data); cudaFree(d_out);
    return 0;
}

