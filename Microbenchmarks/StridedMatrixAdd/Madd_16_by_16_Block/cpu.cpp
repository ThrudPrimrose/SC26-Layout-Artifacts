/*
 * elementwise_cpu_tiled.cpp
 *
 * CPU analogue of CUDA experiment:
 * - direct (no tiling)
 * - tiled (stack "shared memory")
 * - inner tiling sizes: 1,2,4,8,16,32,64
 * - layout-aware load patterns (row/col)
 *
 * Usage:
 *   ./cpu_elem <csv_file> [M] [N]
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

static constexpr double ALPHA = 1.5;
#define NUM_WARMUP 5
#define NUM_MEASURE 100

/* ================================================================
 * Layouts
 * ================================================================ */
struct RowMajor {
    int N;
    static constexpr bool is_col = false;
    inline int operator()(int i, int j) const { return i * N + j; }
};

struct ColMajor {
    int M;
    static constexpr bool is_col = true;
    inline int operator()(int i, int j) const { return j * M + i; }
};

/* ================================================================
 * PRNG (same as CUDA)
 * ================================================================ */
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

static inline uint64_t xorshift64(uint64_t s) {
    s ^= s << 13; s ^= s >> 7; s ^= s << 17;
    return s;
}

/* ================================================================
 * Init
 * ================================================================ */
template<typename Layout>
void init(double* buf, Layout lay, uint64_t seed, int M, int N)
{
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            uint64_t s = (uint64_t)i * N + j;
            buf[lay(i,j)] =
                (double)(xorshift64(splitmix64(seed + s)) % 10000) / 100.0;
        }
    }
}

/* ================================================================
 * DIRECT (no tiling)
 * ================================================================ */
template<typename LA, typename LB>
void kernel_direct(const double* A, double* B,
                   LA la, LB lb, int M, int N)
{
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            int ai = la(i,j), bi = lb(i,j);
            B[bi] = ALPHA * (A[ai] + B[bi]);
        }
    }
}

/* ================================================================
 * TILED (stack "shared memory")
 * Layout-aware load pattern
 * ================================================================ */
template<int INNER, typename LA, typename LB>
void kernel_tiled(const double* A, double* B,
                  LA la, LB lb, int M, int N)
{
    constexpr int TILE_R = 16;
    constexpr int TILE_C = INNER;

    const int tiles_i = (M + TILE_R - 1) / TILE_R;
    const int tiles_j = (N + TILE_C - 1) / TILE_C;

    #pragma omp parallel for collapse(2) schedule(static)
    for (int ti = 0; ti < tiles_i; ti++) {
        for (int tj = 0; tj < tiles_j; tj++) {

            int bi = ti * TILE_R;
            int bj = tj * TILE_C;

            int rows = std::min(TILE_R, M - bi);
            int cols = std::min(TILE_C, N - bj);

            alignas(64) double A_s[TILE_R * TILE_C];
            alignas(64) double B_s[TILE_R * TILE_C];

            /* ---------------------------
             * LOAD (layout-optimized)
             * --------------------------- */
            if constexpr (LA::is_col) {
                // unit stride in i
                for (int j = 0; j < cols; j++)
                    for (int i = 0; i < rows; i++)
                        A_s[i*TILE_C + j] =
                            A[la(bi+i, bj+j)];
            } else {
                // unit stride in j
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        A_s[i*TILE_C + j] =
                            A[la(bi+i, bj+j)];
            }

            if constexpr (LB::is_col) {
                for (int j = 0; j < cols; j++)
                    for (int i = 0; i < rows; i++)
                        B_s[i*TILE_C + j] =
                            B[lb(bi+i, bj+j)];
            } else {
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        B_s[i*TILE_C + j] =
                            B[lb(bi+i, bj+j)];
            }

            /* ---------------------------
             * COMPUTE
             * --------------------------- */
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    B_s[i*TILE_C + j] =
                        ALPHA * (A_s[i*TILE_C + j] +
                                 B_s[i*TILE_C + j]);

            /* ---------------------------
             * STORE (layout-optimized)
             * --------------------------- */
            if constexpr (LB::is_col) {
                for (int j = 0; j < cols; j++)
                    for (int i = 0; i < rows; i++)
                        B[lb(bi+i, bj+j)] =
                            B_s[i*TILE_C + j];
            } else {
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        B[lb(bi+i, bj+j)] =
                            B_s[i*TILE_C + j];
            }
        }
    }
}

/* ================================================================
 * Dispatch
 * ================================================================ */
template<typename LA, typename LB>
void launch(int variant,
            const double* A, double* B,
            LA la, LB lb,
            int M, int N)
{
    switch (variant) {
        case 0: kernel_direct(A,B,la,lb,M,N); break;
        case 1: kernel_tiled<1 >(A,B,la,lb,M,N); break;
        case 2: kernel_tiled<2 >(A,B,la,lb,M,N); break;
        case 3: kernel_tiled<4 >(A,B,la,lb,M,N); break;
        case 4: kernel_tiled<8 >(A,B,la,lb,M,N); break;
        case 5: kernel_tiled<16>(A,B,la,lb,M,N); break;
        case 6: kernel_tiled<32>(A,B,la,lb,M,N); break;
        case 7: kernel_tiled<64>(A,B,la,lb,M,N); break;
    }
}

/* ================================================================
 * Checksum
 * ================================================================ */
double checksum(const double* buf, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += buf[i];
    return s;
}

/* ================================================================
 * MAIN
 * ================================================================ */
int main(int argc, char** argv)
{
    if (argc < 2) {
        printf("Usage: %s <csv> [M] [N]\n", argv[0]);
        return 1;
    }

    const char* csv = argv[1];
    int M = (argc>2)?atoi(argv[2]):4096;
    int N = (argc>3)?atoi(argv[3]):4096;

    size_t bytes = (size_t)M*N*sizeof(double);
    double* A = (double*)malloc(bytes);
    double* B = (double*)malloc(bytes);
    double* B0= (double*)malloc(bytes);

    FILE* fp = fopen(csv,"w");
    fprintf(fp,"kernel,layout,M,N,run,time_ms,bw,checksum\n");

    for (int layout = 0; layout < 4; layout++) {

        bool a_col = (layout==2||layout==3);
        bool b_col = (layout==1||layout==3);

        for (int k = 0; k < 8; k++) {

            printf("layout=%d kernel=%d\n",layout,k);

            if (!a_col && !b_col) {
                init(A, RowMajor{N}, 123, M,N);
                init(B, RowMajor{N}, 456, M,N);
            } else if (!a_col && b_col) {
                init(A, RowMajor{N}, 123, M,N);
                init(B, ColMajor{M}, 456, M,N);
            } else if (a_col && !b_col) {
                init(A, ColMajor{M}, 123, M,N);
                init(B, RowMajor{N}, 456, M,N);
            } else {
                init(A, ColMajor{M}, 123, M,N);
                init(B, ColMajor{M}, 456, M,N);
            }

            memcpy(B0,B,bytes);

            auto run_once = [&](double& ms){
                auto t0 = std::chrono::high_resolution_clock::now();

                if (!a_col && !b_col)
                    launch(k,A,B,RowMajor{N},RowMajor{N},M,N);
                else if (!a_col && b_col)
                    launch(k,A,B,RowMajor{N},ColMajor{M},M,N);
                else if (a_col && !b_col)
                    launch(k,A,B,ColMajor{M},RowMajor{N},M,N);
                else
                    launch(k,A,B,ColMajor{M},ColMajor{M},M,N);

                auto t1 = std::chrono::high_resolution_clock::now();
                ms = std::chrono::duration<double, std::milli>(t1-t0).count();
            };

            double ref_ms;
            run_once(ref_ms);
            double ref = checksum(B,M*N);

            for(int w=0;w<NUM_WARMUP;w++){
                memcpy(B0,B,bytes);
                run_once(ref_ms);
            }

            for(int r=0;r<NUM_MEASURE;r++){
                memcpy(B,B0,bytes);
                double ms;
                run_once(ms);

                double bw = (double)bytes*3 / (ms*1e-3) / 1e9;

                fprintf(fp,"%d,%d,%d,%d,%d,%.4f,%.3f,%.6f\n",
                        k,layout,M,N,r,ms,bw,ref);
            }
        }
    }

    fclose(fp);
    printf("Done -> %s\n",csv);
}