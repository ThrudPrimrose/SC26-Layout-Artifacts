/*
 * elementwise_16x16.cu
 *
 * B[:] = alpha * (A[:] + B[:])
 *
 * Fixed 16x16 thread block (flat variants use 256x1).
 * 12 kernel variants x 4 layout permutations = 48 configs.
 *
 * Direct (no shared memory):
 *   0  direct       2D 16x16, 1x1 elem/thread
 *   1  flat_row     1D x256, row-major traversal (idx/N, idx%N)
 *   2  flat_col     1D x256, col-major traversal (idx%M, idx/M)
 *   3  direct_2x2   2D 16x16, 2x2 elem/thread
 *   4  direct_4x1   2D 16x16, 4 cols per thread
 *   5  direct_1x4   2D 16x16, 4 rows per thread
 *
 * Shared memory (cooperative load -> smem compute -> cooperative store):
 *   6  smem         16x16 tile, 1x1 elem/thread
 *   7  smem_flat_row  16x16 tile via 256x1 block, row-major tile order
 *   8  smem_flat_col  16x16 tile via 256x1 block, col-major tile order
 *   9  smem_2x2     32x32 tile, 2x2 elem/thread
 *  10  smem_4x1     64x16 tile, 4 cols per thread
 *  11  smem_1x4     16x64 tile, 4 rows per thread
 *
 * Usage: ./elementwise_16x16 <csv_file> [M] [N]
 *        defaults: M=4096, N=4096
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cuda_runtime.h>

static constexpr int BX = 16;
static constexpr int BY = 16;
static constexpr double ALPHA = 1.5;

#define NUM_WARMUP  5
#define NUM_MEASURE 100

#define CUDA_CHECK(call) do {                                       \
    cudaError_t _e = (call);                                        \
    if (_e != cudaSuccess) {                                        \
        fprintf(stderr, "CUDA error %s:%d: %s\n",                  \
                __FILE__, __LINE__, cudaGetErrorString(_e));        \
        exit(1);                                                    \
    }                                                               \
} while (0)

/* ================================================================
 *  Layout functors
 * ================================================================ */
struct RowMajor {
    int N;
    static constexpr bool is_col = false;
    __device__ __host__ __forceinline__
    int operator()(int i, int j) const { return i * N + j; }
};

struct ColMajor {
    int M;
    static constexpr bool is_col = true;
    __device__ __host__ __forceinline__
    int operator()(int i, int j) const { return j * M + i; }
};

/* ================================================================
 *  Cooperative tile load/store helpers
 *
 *  TILE_R x TILE_C tile starting at (bi, bj) in global coords.
 *  nthd threads cooperatively iterate over TILE_R*TILE_C elements.
 *
 *  Global memory: unit stride guaranteed.
 *    row-major source: consecutive threads -> consecutive j -> stride-1
 *    col-major source: consecutive threads -> consecutive i -> stride-1
 *
 *  Shared memory: row-major with padded stride (STRIDE = TILE_C + 1)
 *    to eliminate bank conflicts when the global iteration order
 *    walks down columns (col-major load/store).
 *    smem[li * STRIDE + lj]
 *
 *  Callers must declare smem as:  double smem[TILE_R * STRIDE]
 *  where STRIDE = TILE_C + 1.
 * ================================================================ */

template<int TILE_R, int TILE_C, int STRIDE, typename Layout>
__device__ __forceinline__
void coop_load(double* __restrict__ smem,
               const double* __restrict__ gmem,
               Layout lay, int bi, int bj,
               int M, int N, int tid, int nthd)
{
    constexpr int TELEMS = TILE_R * TILE_C;
    for (int k = tid; k < TELEMS; k += nthd) {
        int li, lj;
        if constexpr (Layout::is_col) {
            /* col-major global: consecutive k -> consecutive i -> unit stride */
            li = k % TILE_R;
            lj = k / TILE_R;
        } else {
            /* row-major global: consecutive k -> consecutive j -> unit stride */
            li = k / TILE_C;
            lj = k % TILE_C;
        }
        int gi = bi + li, gj = bj + lj;
        smem[li * STRIDE + lj] = (gi < M && gj < N) ? gmem[lay(gi, gj)] : 0.0;
    }
}

template<int TILE_R, int TILE_C, int STRIDE, typename Layout>
__device__ __forceinline__
void coop_store(const double* __restrict__ smem,
                double* __restrict__ gmem,
                Layout lay, int bi, int bj,
                int M, int N, int tid, int nthd)
{
    constexpr int TELEMS = TILE_R * TILE_C;
    for (int k = tid; k < TELEMS; k += nthd) {
        int li, lj;
        if constexpr (Layout::is_col) {
            li = k % TILE_R;
            lj = k / TILE_R;
        } else {
            li = k / TILE_C;
            lj = k % TILE_C;
        }
        int gi = bi + li, gj = bj + lj;
        if (gi < M && gj < N)
            gmem[lay(gi, gj)] = smem[li * STRIDE + lj];
    }
}

/* ================================================================
 *  PRNG
 * ================================================================ */
__device__ __host__ __forceinline__
uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

__device__ __host__ __forceinline__
uint64_t xorshift64(uint64_t s) {
    s ^= s << 13; s ^= s >> 7; s ^= s << 17;
    return s;
}

/* ================================================================
 *  Init kernels
 * ================================================================ */
__global__ void init_row(double* __restrict__ buf, uint64_t base_seed,
                         int M, int N)
{
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
         idx < M * N; idx += blockDim.x * gridDim.x) {
        int i = idx / N, j = idx % N;
        uint64_t seed = (uint64_t)i * N + j;
        buf[i * N + j] = (double)(xorshift64(splitmix64(base_seed + seed)) % 10000) / 100.0;
    }
}

__global__ void init_col(double* __restrict__ buf, uint64_t base_seed,
                         int M, int N)
{
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
         idx < M * N; idx += blockDim.x * gridDim.x) {
        int i = idx / N, j = idx % N;
        uint64_t seed = (uint64_t)i * N + j;
        buf[j * M + i] = (double)(xorshift64(splitmix64(base_seed + seed)) % 10000) / 100.0;
    }
}

/* ================================================================
 *  DIRECT KERNELS (no shared memory)
 * ================================================================ */

/* 0: direct -- 2D 16x16, 1 elem/thread */
template<typename LA, typename LB>
__global__ void __launch_bounds__(BX * BY)
kernel_direct(const double* __restrict__ A, double* __restrict__ B,
              LA la, LB lb, int M, int N)
{
    int j = blockIdx.x * BX + threadIdx.x;
    int i = blockIdx.y * BY + threadIdx.y;
    if (i < M && j < N) {
        int ai = la(i, j), bi = lb(i, j);
        B[bi] = ALPHA * (A[ai] + B[bi]);
    }
}

/* 1: flat_row -- 1D x256, row-major traversal */
template<typename LA, typename LB>
__global__ void __launch_bounds__(256)
kernel_flat_row(const double* __restrict__ A, double* __restrict__ B,
                LA la, LB lb, int M, int N)
{
    int total = M * N;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
         idx < total; idx += blockDim.x * gridDim.x) {
        int i = idx / N, j = idx % N;
        int ai = la(i, j), bi = lb(i, j);
        B[bi] = ALPHA * (A[ai] + B[bi]);
    }
}

/* 2: flat_col -- 1D x256, col-major traversal */
template<typename LA, typename LB>
__global__ void __launch_bounds__(256)
kernel_flat_col(const double* __restrict__ A, double* __restrict__ B,
                LA la, LB lb, int M, int N)
{
    int total = M * N;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
         idx < total; idx += blockDim.x * gridDim.x) {
        int i = idx % M, j = idx / M;
        int ai = la(i, j), bi = lb(i, j);
        B[bi] = ALPHA * (A[ai] + B[bi]);
    }
}

/* 3: direct_2x2 -- 2D 16x16, 2x2 elem/thread */
template<typename LA, typename LB>
__global__ void __launch_bounds__(BX * BY)
kernel_direct_2x2(const double* __restrict__ A, double* __restrict__ B,
                  LA la, LB lb, int M, int N)
{
    int base_j = blockIdx.x * (BX * 2) + threadIdx.x;
    int base_i = blockIdx.y * (BY * 2) + threadIdx.y;
    #pragma unroll
    for (int dy = 0; dy < 2; dy++) {
        int i = base_i + dy * BY;
        #pragma unroll
        for (int dx = 0; dx < 2; dx++) {
            int j = base_j + dx * BX;
            if (i < M && j < N) {
                int ai = la(i, j), bi = lb(i, j);
                B[bi] = ALPHA * (A[ai] + B[bi]);
            }
        }
    }
}

/* 4: direct_4x1 -- 2D 16x16, 4 cols per thread */
template<typename LA, typename LB>
__global__ void __launch_bounds__(BX * BY)
kernel_direct_4x1(const double* __restrict__ A, double* __restrict__ B,
                  LA la, LB lb, int M, int N)
{
    int base_j = blockIdx.x * (BX * 4) + threadIdx.x;
    int i      = blockIdx.y * BY + threadIdx.y;
    if (i < M) {
        #pragma unroll
        for (int dx = 0; dx < 4; dx++) {
            int j = base_j + dx * BX;
            if (j < N) {
                int ai = la(i, j), bi = lb(i, j);
                B[bi] = ALPHA * (A[ai] + B[bi]);
            }
        }
    }
}

/* 5: direct_1x4 -- 2D 16x16, 4 rows per thread */
template<typename LA, typename LB>
__global__ void __launch_bounds__(BX * BY)
kernel_direct_1x4(const double* __restrict__ A, double* __restrict__ B,
                  LA la, LB lb, int M, int N)
{
    int j      = blockIdx.x * BX + threadIdx.x;
    int base_i = blockIdx.y * (BY * 4) + threadIdx.y;
    if (j < N) {
        #pragma unroll
        for (int dy = 0; dy < 4; dy++) {
            int i = base_i + dy * BY;
            if (i < M) {
                int ai = la(i, j), bi = lb(i, j);
                B[bi] = ALPHA * (A[ai] + B[bi]);
            }
        }
    }
}

/* ================================================================
 *  SHARED MEMORY KERNELS
 *  Each mirrors a direct variant above.
 *  Pattern: cooperative load -> sync -> compute in smem -> sync -> store
 * ================================================================ */

/* 6: smem -- 16x16 tile, 1x1 elem/thread, 16x16 block */
template<typename LA, typename LB>
__global__ void __launch_bounds__(BX * BY)
kernel_smem(const double* __restrict__ A, double* __restrict__ B,
            LA la, LB lb, int M, int N)
{
    constexpr int S = BX + 1;  /* padded stride */
    __shared__ double A_s[BY * S];
    __shared__ double B_s[BY * S];

    const int bj = blockIdx.x * BX, bi = blockIdx.y * BY;
    const int tid  = threadIdx.y * BX + threadIdx.x;
    const int nthd = BX * BY;

    coop_load<BY, BX, S>(A_s, A, la, bi, bj, M, N, tid, nthd);
    coop_load<BY, BX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
    __syncthreads();

    const int si = threadIdx.y * S + threadIdx.x;
    B_s[si] = ALPHA * (A_s[si] + B_s[si]);
    __syncthreads();

    coop_store<BY, BX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
}

/* 7: smem_flat_row -- 16x16 tile via 256x1 block, row-major tile order */
template<typename LA, typename LB>
__global__ void __launch_bounds__(256)
kernel_smem_flat_row(const double* __restrict__ A, double* __restrict__ B,
                     LA la, LB lb, int M, int N)
{
    constexpr int S = BX + 1;
    __shared__ double A_s[BY * S];
    __shared__ double B_s[BY * S];

    const int tiles_j = (N + BX - 1) / BX;
    const int tiles_i = (M + BY - 1) / BY;
    const int total_tiles = tiles_i * tiles_j;
    const int tid = threadIdx.x;
    const int nthd = 256;

    for (int t = blockIdx.x; t < total_tiles; t += gridDim.x) {
        int ti = t / tiles_j;
        int tj = t % tiles_j;
        int bi = ti * BY, bj = tj * BX;

        coop_load<BY, BX, S>(A_s, A, la, bi, bj, M, N, tid, nthd);
        coop_load<BY, BX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
        __syncthreads();

        /* 256 threads, 256 elements: remap tid -> (li,lj) in padded smem */
        int cli = tid / BX;
        int clj = tid % BX;
        B_s[cli * S + clj] = ALPHA * (A_s[cli * S + clj] + B_s[cli * S + clj]);
        __syncthreads();

        coop_store<BY, BX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
        __syncthreads();
    }
}

/* 8: smem_flat_col -- 16x16 tile via 256x1 block, col-major tile order */
template<typename LA, typename LB>
__global__ void __launch_bounds__(256)
kernel_smem_flat_col(const double* __restrict__ A, double* __restrict__ B,
                     LA la, LB lb, int M, int N)
{
    constexpr int S = BX + 1;
    __shared__ double A_s[BY * S];
    __shared__ double B_s[BY * S];

    const int tiles_j = (N + BX - 1) / BX;
    const int tiles_i = (M + BY - 1) / BY;
    const int total_tiles = tiles_i * tiles_j;
    const int tid = threadIdx.x;
    const int nthd = 256;

    for (int t = blockIdx.x; t < total_tiles; t += gridDim.x) {
        int tj = t / tiles_i;
        int ti = t % tiles_i;
        int bi = ti * BY, bj = tj * BX;

        coop_load<BY, BX, S>(A_s, A, la, bi, bj, M, N, tid, nthd);
        coop_load<BY, BX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
        __syncthreads();

        int cli = tid / BX;
        int clj = tid % BX;
        B_s[cli * S + clj] = ALPHA * (A_s[cli * S + clj] + B_s[cli * S + clj]);
        __syncthreads();

        coop_store<BY, BX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
        __syncthreads();
    }
}

/* 9: smem_2x2 -- 32x32 tile, 2x2 elem/thread, 16x16 block */
template<typename LA, typename LB>
__global__ void __launch_bounds__(BX * BY)
kernel_smem_2x2(const double* __restrict__ A, double* __restrict__ B,
                LA la, LB lb, int M, int N)
{
    constexpr int TX = BX * 2;
    constexpr int TY = BY * 2;
    constexpr int S  = TX + 1;  /* padded stride */
    __shared__ double A_s[TY * S];
    __shared__ double B_s[TY * S];

    const int bj = blockIdx.x * TX, bi = blockIdx.y * TY;
    const int tid  = threadIdx.y * BX + threadIdx.x;
    const int nthd = BX * BY;

    coop_load<TY, TX, S>(A_s, A, la, bi, bj, M, N, tid, nthd);
    coop_load<TY, TX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
    __syncthreads();

    #pragma unroll
    for (int dy = 0; dy < 2; dy++) {
        int li = threadIdx.y + dy * BY;
        #pragma unroll
        for (int dx = 0; dx < 2; dx++) {
            int lj = threadIdx.x + dx * BX;
            int si = li * S + lj;
            B_s[si] = ALPHA * (A_s[si] + B_s[si]);
        }
    }
    __syncthreads();

    coop_store<TY, TX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
}

/* 10: smem_4x1 -- 64x16 tile, 4 cols per thread, 16x16 block */
template<typename LA, typename LB>
__global__ void __launch_bounds__(BX * BY)
kernel_smem_4x1(const double* __restrict__ A, double* __restrict__ B,
                LA la, LB lb, int M, int N)
{
    constexpr int TX = BX * 4;
    constexpr int TY = BY;
    constexpr int S  = TX + 1;
    __shared__ double A_s[TY * S];
    __shared__ double B_s[TY * S];

    const int bj = blockIdx.x * TX, bi = blockIdx.y * TY;
    const int tid  = threadIdx.y * BX + threadIdx.x;
    const int nthd = BX * BY;

    coop_load<TY, TX, S>(A_s, A, la, bi, bj, M, N, tid, nthd);
    coop_load<TY, TX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
    __syncthreads();

    int li = threadIdx.y;
    #pragma unroll
    for (int dx = 0; dx < 4; dx++) {
        int lj = threadIdx.x + dx * BX;
        int si = li * S + lj;
        B_s[si] = ALPHA * (A_s[si] + B_s[si]);
    }
    __syncthreads();

    coop_store<TY, TX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
}

/* 11: smem_1x4 -- 16x64 tile, 4 rows per thread, 16x16 block */
template<typename LA, typename LB>
__global__ void __launch_bounds__(BX * BY)
kernel_smem_1x4(const double* __restrict__ A, double* __restrict__ B,
                LA la, LB lb, int M, int N)
{
    constexpr int TX = BX;
    constexpr int TY = BY * 4;
    constexpr int S  = TX + 1;
    __shared__ double A_s[TY * S];
    __shared__ double B_s[TY * S];

    const int bj = blockIdx.x * TX, bi = blockIdx.y * TY;
    const int tid  = threadIdx.y * BX + threadIdx.x;
    const int nthd = BX * BY;

    coop_load<TY, TX, S>(A_s, A, la, bi, bj, M, N, tid, nthd);
    coop_load<TY, TX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
    __syncthreads();

    int lj = threadIdx.x;
    #pragma unroll
    for (int dy = 0; dy < 4; dy++) {
        int li = threadIdx.y + dy * BY;
        int si = li * S + lj;
        B_s[si] = ALPHA * (A_s[si] + B_s[si]);
    }
    __syncthreads();

    coop_store<TY, TX, S>(B_s, B, lb, bi, bj, M, N, tid, nthd);
}

/* ================================================================
 *  Host: layout + kernel tables
 * ================================================================ */
struct LayoutConfig {
    int a_layout;
    int b_layout;
    const char* name;
};

static const LayoutConfig LAYOUTS[4] = {
    {0, 0, "Ar_Br"},
    {0, 1, "Ar_Bc"},
    {1, 0, "Ac_Br"},
    {1, 1, "Ac_Bc"},
};

struct KernelInfo {
    int id;
    const char* tag;
    int uses_smem;
};

static const KernelInfo KERNELS[] = {
    { 0, "direct",        0},
    { 1, "flat_row",      0},
    { 2, "flat_col",      0},
    { 3, "direct_2x2",    0},
    { 4, "direct_4x1",    0},
    { 5, "direct_1x4",    0},
    { 6, "smem",          1},
    { 7, "smem_flat_row", 1},
    { 8, "smem_flat_col", 1},
    { 9, "smem_2x2",      1},
    {10, "smem_4x1",      1},
    {11, "smem_1x4",      1},
};
static constexpr int NUM_KERNELS = sizeof(KERNELS) / sizeof(KERNELS[0]);

/* ================================================================
 *  Host: init helper
 * ================================================================ */
static void init_array(double* d_buf, int layout, uint64_t seed,
                       int M, int N)
{
    int nblk = (M * N + 255) / 256;
    if (layout == 0)
        init_row<<<nblk, 256>>>(d_buf, seed, M, N);
    else
        init_col<<<nblk, 256>>>(d_buf, seed, M, N);
}

/* ================================================================
 *  Host: launch dispatch
 * ================================================================ */
#define DISPATCH(K, GRID, BLOCK) do {                                \
    switch (key) {                                                    \
    case 0: K<<<GRID,BLOCK>>>(d_A,d_B,RowMajor{N},RowMajor{N},M,N); break; \
    case 1: K<<<GRID,BLOCK>>>(d_A,d_B,RowMajor{N},ColMajor{M},M,N); break; \
    case 2: K<<<GRID,BLOCK>>>(d_A,d_B,ColMajor{M},RowMajor{N},M,N); break; \
    case 3: K<<<GRID,BLOCK>>>(d_A,d_B,ColMajor{M},ColMajor{M},M,N); break; \
    }                                                                \
} while (0)

static void launch(const double* d_A, double* d_B,
                   int a_layout, int b_layout, int kernel_id,
                   int M, int N)
{
    int key = a_layout * 2 + b_layout;

    switch (kernel_id) {
    case 0: { /* direct */
        dim3 g((N+BX-1)/BX, (M+BY-1)/BY);  dim3 b(BX, BY);
        DISPATCH(kernel_direct, g, b);
        break;
    }
    case 1: { /* flat_row */
        int nblk = (M*N+255)/256;
        dim3 g(nblk); dim3 b(256);
        DISPATCH(kernel_flat_row, g, b);
        break;
    }
    case 2: { /* flat_col */
        int nblk = (M*N+255)/256;
        dim3 g(nblk); dim3 b(256);
        DISPATCH(kernel_flat_col, g, b);
        break;
    }
    case 3: { /* direct_2x2 */
        dim3 g((N+BX*2-1)/(BX*2), (M+BY*2-1)/(BY*2));  dim3 b(BX, BY);
        DISPATCH(kernel_direct_2x2, g, b);
        break;
    }
    case 4: { /* direct_4x1 */
        dim3 g((N+BX*4-1)/(BX*4), (M+BY-1)/BY);  dim3 b(BX, BY);
        DISPATCH(kernel_direct_4x1, g, b);
        break;
    }
    case 5: { /* direct_1x4 */
        dim3 g((N+BX-1)/BX, (M+BY*4-1)/(BY*4));  dim3 b(BX, BY);
        DISPATCH(kernel_direct_1x4, g, b);
        break;
    }
    case 6: { /* smem */
        dim3 g((N+BX-1)/BX, (M+BY-1)/BY);  dim3 b(BX, BY);
        DISPATCH(kernel_smem, g, b);
        break;
    }
    case 7: { /* smem_flat_row */
        int tiles = ((M+BY-1)/BY) * ((N+BX-1)/BX);
        int nblk = (tiles < 1024) ? tiles : 1024;
        dim3 g(nblk); dim3 b(256);
        DISPATCH(kernel_smem_flat_row, g, b);
        break;
    }
    case 8: { /* smem_flat_col */
        int tiles = ((M+BY-1)/BY) * ((N+BX-1)/BX);
        int nblk = (tiles < 1024) ? tiles : 1024;
        dim3 g(nblk); dim3 b(256);
        DISPATCH(kernel_smem_flat_col, g, b);
        break;
    }
    case 9: { /* smem_2x2 */
        dim3 g((N+BX*2-1)/(BX*2), (M+BY*2-1)/(BY*2));  dim3 b(BX, BY);
        DISPATCH(kernel_smem_2x2, g, b);
        break;
    }
    case 10: { /* smem_4x1 */
        dim3 g((N+BX*4-1)/(BX*4), (M+BY-1)/BY);  dim3 b(BX, BY);
        DISPATCH(kernel_smem_4x1, g, b);
        break;
    }
    case 11: { /* smem_1x4 */
        dim3 g((N+BX-1)/BX, (M+BY*4-1)/(BY*4));  dim3 b(BX, BY);
        DISPATCH(kernel_smem_1x4, g, b);
        break;
    }
    }
}

#undef DISPATCH

/* ================================================================
 *  Host: checksum
 * ================================================================ */
static double checksum_host(const double* buf, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += buf[i];
    return s;
}

/* ================================================================
 *  main
 * ================================================================ */
int main(int argc, char** argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <csv_file> [M] [N]\n", argv[0]);
        return 1;
    }
    const char* csv_path = argv[1];
    int M = (argc >= 3) ? atoi(argv[2]) : 4096;
    int N = (argc >= 4) ? atoi(argv[3]) : 4096;

    printf("elementwise_16x16: M=%d N=%d  kernels=%d  layouts=4  "
           "warmup=%d  measure=%d  total_configs=%d\n",
           M, N, NUM_KERNELS, NUM_WARMUP, NUM_MEASURE, NUM_KERNELS * 4);

    const int total = M * N;
    const size_t bytes = (size_t)total * sizeof(double);
    const double data_bytes = (double)total * sizeof(double) * 3.0;

    double* h_B = (double*)malloc(bytes);
    if (!h_B) { fprintf(stderr, "malloc failed\n"); return 1; }

    double *d_A, *d_B, *d_B_orig;
    CUDA_CHECK(cudaMalloc(&d_A, bytes));
    CUDA_CHECK(cudaMalloc(&d_B, bytes));
    CUDA_CHECK(cudaMalloc(&d_B_orig, bytes));

    cudaEvent_t ev_start, ev_stop;
    CUDA_CHECK(cudaEventCreate(&ev_start));
    CUDA_CHECK(cudaEventCreate(&ev_stop));

    FILE* fp = fopen(csv_path, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", csv_path); return 1; }
    fprintf(fp, "kernel_name,layout,use_smem,M,N,run,time_ms,bw_gbps,checksum\n");

    /* print header to stdout too */
    printf("kernel_name,layout,use_smem,M,N,run,time_ms,bw_gbps,checksum\n");

    /* ---- 12 kernels x 4 layouts = 48 configs ---- */
    for (int li = 0; li < 4; li++) {
        const LayoutConfig& lc = LAYOUTS[li];

        for (int ki = 0; ki < NUM_KERNELS; ki++) {
            const KernelInfo& kinfo = KERNELS[ki];
            char kname[64];
            snprintf(kname, sizeof(kname), "%s_%s", lc.name, kinfo.tag);

            fprintf(stderr, "--- %s ---\n", kname);

            /* init */
            init_array(d_A, lc.a_layout, 0x123456789ABCDEF0ULL, M, N);
            init_array(d_B, lc.b_layout, 0xFEDCBA9876543210ULL, M, N);
            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(d_B_orig, d_B, bytes, cudaMemcpyDeviceToDevice));

            /* reference checksum */
            launch(d_A, d_B, lc.a_layout, lc.b_layout, kinfo.id, M, N);
            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_B, d_B, bytes, cudaMemcpyDeviceToHost));
            double ref_cs = checksum_host(h_B, total);

            /* warmup */
            for (int w = 0; w < NUM_WARMUP; w++) {
                CUDA_CHECK(cudaMemcpy(d_B, d_B_orig, bytes, cudaMemcpyDeviceToDevice));
                launch(d_A, d_B, lc.a_layout, lc.b_layout, kinfo.id, M, N);
            }
            CUDA_CHECK(cudaDeviceSynchronize());

            /* measured runs */
            for (int r = 0; r < NUM_MEASURE; r++) {
                CUDA_CHECK(cudaMemcpy(d_B, d_B_orig, bytes, cudaMemcpyDeviceToDevice));
                CUDA_CHECK(cudaDeviceSynchronize());

                CUDA_CHECK(cudaEventRecord(ev_start));
                launch(d_A, d_B, lc.a_layout, lc.b_layout, kinfo.id, M, N);
                CUDA_CHECK(cudaEventRecord(ev_stop));
                CUDA_CHECK(cudaEventSynchronize(ev_stop));

                float ms;
                CUDA_CHECK(cudaEventElapsedTime(&ms, ev_start, ev_stop));
                double bw = data_bytes / ((double)ms * 1e-3) / 1e9;

                fprintf(fp, "%s,%s,%d,%d,%d,%d,%.6f,%.3f,%.6f\n",
                        kname, lc.name, kinfo.uses_smem, M, N,
                        r, (double)ms, bw, ref_cs);
                printf("%s,%s,%d,%d,%d,%d,%.6f,%.3f,%.6f\n",
                       kname, lc.name, kinfo.uses_smem, M, N,
                       r, (double)ms, bw, ref_cs);
            }

            fflush(fp);
            fflush(stdout);
            fprintf(stderr, "  checksum=%.6f\n", ref_cs);
        }
    }

    fclose(fp);
    fprintf(stderr, "\nDone. Results: %s\n", csv_path);

    CUDA_CHECK(cudaEventDestroy(ev_start));
    CUDA_CHECK(cudaEventDestroy(ev_stop));
    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_B));
    CUDA_CHECK(cudaFree(d_B_orig));
    free(h_B);

    return 0;
}
