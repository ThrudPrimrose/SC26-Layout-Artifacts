/*
 * Elementwise kernel benchmark: B[:] = alpha * (A[:] + B[:])
 *
 * Compile-time parameters (set via -D):
 *   M_DIM, N_DIM  - matrix dimensions
 *   LAYOUT_A/B     - 0=row-major, 1=col-major, 2=blocked-row, 3=blocked-col
 *   BLK_X, BLK_Y   - thread block dimensions (always 2D)
 *   TILE_X, TILE_Y  - elements per thread in X (col) / Y (row) direction
 *   USE_SMEM        - 0=direct global, 1=shared memory
 *   BSIZE           - block size for blocked layouts (default 16)
 *
 * Usage: ./exe <csv_file>
 * Appends one summary line per run to the CSV.
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <cuda_runtime.h>

/* ---------- defaults ---------- */
#ifndef M_DIM
#define M_DIM 4096
#endif
#ifndef N_DIM
#define N_DIM 4096
#endif
#ifndef LAYOUT_A
#define LAYOUT_A 0
#endif
#ifndef LAYOUT_B
#define LAYOUT_B 0
#endif
#ifndef BLK_X
#define BLK_X 128
#endif
#ifndef BLK_Y
#define BLK_Y 1
#endif
#ifndef TILE_X
#define TILE_X 1
#endif
#ifndef TILE_Y
#define TILE_Y 1
#endif
#ifndef USE_SMEM
#define USE_SMEM 0
#endif
#ifndef BSIZE
#define BSIZE 16
#endif

static constexpr double ALPHA = 1.5;
#define NUM_WARMUP  5
#define NUM_MEASURE 100

/* ---------- tile geometry (compile-time) ---------- */
#define TILE_ROWS (BLK_Y * TILE_Y)
#define TILE_COLS (BLK_X * TILE_X)
#define TILE_ELEMS (TILE_ROWS * TILE_COLS)

/* ---------- error check ---------- */
#define CUDA_CHECK(call) do {                                       \
    cudaError_t _e = (call);                                        \
    if (_e != cudaSuccess) {                                        \
        fprintf(stderr, "CUDA error %s:%d: %s\n",                  \
                __FILE__, __LINE__, cudaGetErrorString(_e));        \
        exit(1);                                                    \
    }                                                               \
} while (0)

/* ================================================================
 *  Layout index macros
 *  (i,j) with i in [0,M_DIM), j in [0,N_DIM)
 * ================================================================ */

/* --- A layout --- */
#if LAYOUT_A == 0
  #define A_IDX(i, j)  ((i) * N_DIM + (j))
#elif LAYOUT_A == 1
  #define A_IDX(i, j)  ((j) * M_DIM + (i))
#elif LAYOUT_A == 2
  /* blocked, blocks in row-major order, elements within block row-major */
  #define A_IDX(i, j) \
      (((i) / BSIZE * (N_DIM / BSIZE) + (j) / BSIZE) * (BSIZE * BSIZE) \
       + ((i) % BSIZE) * BSIZE + (j) % BSIZE)
#elif LAYOUT_A == 3
  /* blocked, blocks in col-major order, elements within block row-major */
  #define A_IDX(i, j) \
      (((j) / BSIZE * (M_DIM / BSIZE) + (i) / BSIZE) * (BSIZE * BSIZE) \
       + ((i) % BSIZE) * BSIZE + (j) % BSIZE)
#endif

/* --- B layout --- */
#if LAYOUT_B == 0
  #define B_IDX(i, j)  ((i) * N_DIM + (j))
#elif LAYOUT_B == 1
  #define B_IDX(i, j)  ((j) * M_DIM + (i))
#elif LAYOUT_B == 2
  #define B_IDX(i, j) \
      (((i) / BSIZE * (N_DIM / BSIZE) + (j) / BSIZE) * (BSIZE * BSIZE) \
       + ((i) % BSIZE) * BSIZE + (j) % BSIZE)
#elif LAYOUT_B == 3
  #define B_IDX(i, j) \
      (((j) / BSIZE * (M_DIM / BSIZE) + (i) / BSIZE) * (BSIZE * BSIZE) \
       + ((i) % BSIZE) * BSIZE + (j) % BSIZE)
#endif

/* ================================================================
 *  PRNG (deterministic, parallelizable)
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
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
}


template<int TILE_R, int TILE_C, int STRIDE, int Layout>
__device__ __forceinline__
void coop_load(double* __restrict__ smem,
               const double* __restrict__ gmem,
               int bi, int bj,
               int M, int N, int tid, int nthd)
{
    constexpr int TELEMS = TILE_R * TILE_C;
    for (int k = tid; k < TELEMS; k += nthd) {
        int li, lj;
        if constexpr (Layout == 0) {
            /* col-major global: consecutive k -> consecutive i -> unit stride */
            li = k % TILE_R;
            lj = k / TILE_R;
        } else {
            /* row-major global: consecutive k -> consecutive j -> unit stride */
            li = k / TILE_C;
            lj = k % TILE_C;
        }
        int gi = bi + li, gj = bj + lj;
        if constexpr (Layout == 0) {
            /* row-major global: consecutive k -> consecutive j -> unit stride */
            smem[li * STRIDE + lj] = (gi < M && gj < N) ? gmem[gi * N + gj] : 0.0;
        } else {
            /* col-major global: consecutive k -> consecutive i -> unit stride */   
            smem[li * STRIDE + lj] = (gi < M && gj < N) ? gmem[gj * M + gi] : 0.0;
        }
    }
}

template<int TILE_R, int TILE_C, int STRIDE, int Layout>
__device__ __forceinline__
void coop_store(const double* __restrict__ smem,
                double* __restrict__ gmem,
                int bi, int bj,
                int M, int N, int tid, int nthd)
{
    constexpr int TELEMS = TILE_R * TILE_C;
    for (int k = tid; k < TELEMS; k += nthd) {
        int li, lj;
        if constexpr (Layout == 0) {
            li = k % TILE_R;
            lj = k / TILE_R;
        } else {
            li = k / TILE_C;
            lj = k % TILE_C;
        }
        int gi = bi + li, gj = bj + lj;
        if (gi < M && gj < N){
            if constexpr (Layout == 0) {
                /* row-major global: consecutive k -> consecutive j -> unit stride */
                gmem[gi * N + gj] = smem[li * STRIDE + lj];
            } else {
                /* col-major global: consecutive k -> consecutive i -> unit stride */
                gmem[gj * M + gi] = smem[li * STRIDE + lj];
            }
        }
    }
}

/* ================================================================
 *  Init kernel -- layout-aware, deterministic per (i,j)
 * ================================================================ */
__global__ void init_kernel(double* __restrict__ A,
                            double* __restrict__ B)
{
    const int total = M_DIM * N_DIM;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
         idx < total;
         idx += blockDim.x * gridDim.x)
    {
        int i = idx / N_DIM;
        int j = idx % N_DIM;
        uint64_t seed = (uint64_t)i * N_DIM + j;
        uint64_t sa = splitmix64(0x123456789ABCDEF0ULL + seed);
        uint64_t sb = splitmix64(0xFEDCBA9876543210ULL + seed);
        double va = (double)(xorshift64(sa) % 10000) / 100.0;
        double vb = (double)(xorshift64(sb) % 10000) / 100.0;
        A[A_IDX(i, j)] = va;
        B[B_IDX(i, j)] = vb;
    }
}

/* ================================================================
 *  Kernel 0: direct global memory (no shared memory)
 *  blockIdx.x -> columns (j),  blockIdx.y -> rows (i)
 *  threadIdx.x -> j,  threadIdx.y -> i
 * ================================================================ */
__global__ void __launch_bounds__(BLK_X * BLK_Y)
kernel_direct(const double* __restrict__ A,
              double*       __restrict__ B)
{
    int base_j = blockIdx.x * TILE_COLS + threadIdx.x;
    int base_i = blockIdx.y * TILE_ROWS + threadIdx.y;

    #pragma unroll
    for (int dy = 0; dy < TILE_Y; dy++) {
        int i = base_i + dy * BLK_Y;
        #pragma unroll
        for (int dx = 0; dx < TILE_X; dx++) {
            int j = base_j + dx * BLK_X;
            if (i < M_DIM && j < N_DIM) {
                int ai = A_IDX(i, j);
                int bi = B_IDX(i, j);
                B[bi] = ALPHA * (A[ai] + B[bi]);
            }
        }
    }
}


/* ================================================================
 *  Kernel 1: shared-memory tile
 *  Load tile cooperatively (coalesced for source layout),
 *  compute in shared mem, store cooperatively.
 * ================================================================ */
__global__ void __launch_bounds__(BLK_X * BLK_Y)
kernel_smem(const double* __restrict__ A,
            double*       __restrict__ B)
{
    extern __shared__ double smem[];
    double* __restrict__ A_s = smem;
    double* __restrict__ B_s = smem + TILE_ELEMS;

    const int block_j = blockIdx.x * TILE_COLS;
    const int block_i = blockIdx.y * TILE_ROWS;
    const int tid     = threadIdx.y * BLK_X + threadIdx.x;
    const int nthds   = BLK_X * BLK_Y;
    constexpr int S  = (BLK_X * TILE_COLS) + 1;  /* padded stride */
    constexpr int nthd = BLK_X * BLK_Y;


    /* ---- cooperative load A ---- */
    coop_load<TILE_ROWS, TILE_COLS, S, LAYOUT_A>(A_s, A, block_i, block_j, M_DIM, N_DIM, tid, nthd);
    coop_load<TILE_ROWS, TILE_COLS, S, LAYOUT_B>(B_s, B, block_i, block_j, M_DIM, N_DIM, tid, nthd);

    __syncthreads();

    /* ---- per-thread compute ---- */
    #pragma unroll
    for (int dy = 0; dy < TILE_Y; dy++) {
        int li = threadIdx.y + dy * BLK_Y;
        #pragma unroll
        for (int dx = 0; dx < TILE_X; dx++) {
            int lj = threadIdx.x + dx * BLK_X;
            int si = li * TILE_COLS + lj;
            B_s[si] = ALPHA * (A_s[si] + B_s[si]);
        }
    }

    __syncthreads();

    /* ---- cooperative store B ---- */
    coop_store<TILE_ROWS, TILE_COLS, S, LAYOUT_B>(B_s, B, block_i, block_j, M_DIM, N_DIM, tid, nthd);

}

/* ================================================================
 *  Host helpers
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
        fprintf(stderr, "Usage: %s <csv_file>\n", argv[0]);
        return 1;
    }
    const char* csv_path = argv[1];

    const int total = M_DIM * N_DIM;
    const size_t bytes = (size_t)total * sizeof(double);

    /* allocate device */
    double *d_A, *d_B, *d_B_orig;
    CUDA_CHECK(cudaMalloc(&d_A, bytes));
    CUDA_CHECK(cudaMalloc(&d_B, bytes));
    CUDA_CHECK(cudaMalloc(&d_B_orig, bytes));

    /* allocate host (for checksum) */
    double* h_B = (double*)malloc(bytes);
    if (!h_B) { fprintf(stderr, "malloc failed\n"); return 1; }

    /* init arrays on GPU */
    {
        int nblk = (total + 255) / 256;
        init_kernel<<<nblk, 256>>>(d_A, d_B);
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    /* backup B_orig */
    CUDA_CHECK(cudaMemcpy(d_B_orig, d_B, bytes, cudaMemcpyDeviceToDevice));

    /* kernel launch geometry */
    const int grid_x = (N_DIM + TILE_COLS - 1) / TILE_COLS;
    const int grid_y = (M_DIM + TILE_ROWS - 1) / TILE_ROWS;
    const dim3 grid(grid_x, grid_y);
    const dim3 block(BLK_X, BLK_Y);
    const size_t smem_bytes = 2 * TILE_ELEMS * sizeof(double);

    /* macro to invoke the correct kernel */
    #define RUN_KERNEL() do {                                         \
        if (USE_SMEM)                                                 \
            kernel_smem<<<grid, block, smem_bytes>>>(d_A, d_B);       \
        else                                                          \
            kernel_direct<<<grid, block>>>(d_A, d_B);                 \
    } while (0)

    /* reference: run once, grab checksum */
    RUN_KERNEL();
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(h_B, d_B, bytes, cudaMemcpyDeviceToHost));
    double ref_cs = checksum_host(h_B, total);

    /* warmup */
    for (int w = 0; w < NUM_WARMUP; w++) {
        CUDA_CHECK(cudaMemcpy(d_B, d_B_orig, bytes, cudaMemcpyDeviceToDevice));
        RUN_KERNEL();
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    /* timed runs -- write each run individually to CSV */
    cudaEvent_t ev_start, ev_stop;
    CUDA_CHECK(cudaEventCreate(&ev_start));
    CUDA_CHECK(cudaEventCreate(&ev_stop));

    FILE* fp = fopen(csv_path, "a");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", csv_path); return 1; }

    double data_bytes = (double)total * sizeof(double) * 3.0;

    for (int r = 0; r < NUM_MEASURE; r++) {
        CUDA_CHECK(cudaMemcpy(d_B, d_B_orig, bytes, cudaMemcpyDeviceToDevice));
        CUDA_CHECK(cudaDeviceSynchronize());

        CUDA_CHECK(cudaEventRecord(ev_start));
        RUN_KERNEL();
        CUDA_CHECK(cudaEventRecord(ev_stop));
        CUDA_CHECK(cudaEventSynchronize(ev_stop));

        float elapsed_ms;
        CUDA_CHECK(cudaEventElapsedTime(&elapsed_ms, ev_start, ev_stop));

        double bw_gbps = data_bytes / ((double)elapsed_ms * 1e-3) / 1e9;

        /* use_smem,layout_a,layout_b,blk_x,blk_y,tile_x,tile_y,bsize,M,N,
           run,time_ms,bw_gbps,checksum */
        fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,"
                    "%d,%.6f,%.3f,%.6f\n",
                USE_SMEM, LAYOUT_A, LAYOUT_B, BLK_X, BLK_Y,
                TILE_X, TILE_Y, BSIZE, M_DIM, N_DIM,
                r, (double)elapsed_ms, bw_gbps, ref_cs);
    }

    fclose(fp);

    /* stdout: just a one-line summary */
    printf("smem=%d la=%d lb=%d blk=%dx%d tile=%dx%d bsize=%d M=%d N=%d  done (%d runs)\n",
           USE_SMEM, LAYOUT_A, LAYOUT_B, BLK_X, BLK_Y,
           TILE_X, TILE_Y, BSIZE, M_DIM, N_DIM, NUM_MEASURE);

    /* cleanup */
    CUDA_CHECK(cudaEventDestroy(ev_start));
    CUDA_CHECK(cudaEventDestroy(ev_stop));
    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_B));
    CUDA_CHECK(cudaFree(d_B_orig));
    free(h_B);

    return 0;
}
