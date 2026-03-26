/*
 * bench_cpu_simple.cpp -- CPU z_v_grad_w: omp parallel for, V1-V4, no collapse
 *
 * Compile:  g++ -O3 -fopenmp -march=native -std=c++17 bench_cpu_simple.cpp -o bench_cpu_simple
 */

#include "bench_common.h"
#include <omp.h>

/* ================================================================ */
/*  CPU kernels -- omp parallel for only                             */
/* ================================================================ */

template<int V>
static void cpu_kernel(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev)
{
    if constexpr (V <= 2) {
        #pragma omp parallel for schedule(static)
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) { STENCIL_BODY(V) }
    } else {
        #pragma omp parallel for schedule(static)
        for (int je = 0; je < N; je++)
            for (int jk = 0; jk < nlev; jk++) { STENCIL_BODY(V) }
    }
}

typedef void (*kern_t)(double*, const double*, const double*,
    const double*, const int*, const double*, const double*,
    const double*, const double*, const int*, int, int);

static kern_t kern_tbl[] = {
    cpu_kernel<1>, cpu_kernel<2>, cpu_kernel<3>, cpu_kernel<4>
};

/* ================================================================ */
/*  Cache flush                                                      */
/* ================================================================ */

static constexpr int FLUSH_N = 8192;
static constexpr int FLUSH_STEPS = 3;

static void flush_caches()
{
    double *A = new double[FLUSH_N * FLUSH_N];
    double *B = new double[FLUSH_N * FLUSH_N];

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < FLUSH_N * FLUSH_N; i++) {
        A[i] = (double)(i % 1000) * 0.001;
        B[i] = A[i];
    }

    for (int s = 0; s < FLUSH_STEPS; s++) {
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < FLUSH_N - 1; i++)
            for (int j = 1; j < FLUSH_N - 1; j++)
                B[i * FLUSH_N + j] = 0.25 * (
                    A[(i-1)*FLUSH_N + j] + A[(i+1)*FLUSH_N + j] +
                    A[i*FLUSH_N + (j-1)] + A[i*FLUSH_N + (j+1)]);
        std::swap(A, B);
    }

    volatile double sink = A[FLUSH_N + 1];
    (void)sink;

    delete[] A;
    delete[] B;
}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main()
{
    const int N = NPROMA;
    std::mt19937 rng(42);

    VertData vd;
    vd.init(N, rng);

    int* cell_logical = new int[N * 2];

    printf("OMP threads: %d\n", omp_get_max_threads());
    printf("NPROMA:      %d\n", N);

    int nlev = 96;
    size_t bytes_2d = (size_t)N * nlev * sizeof(double);
    double total_bytes = (double)(
        4 * bytes_2d          /* vn_ie, w, z_vt_ie, z_w_v read */
        + 1 * bytes_2d          /* out write */
        + 3 * N * sizeof(double)/* inv_dual, inv_primal, tangent */
        + 2 * N * sizeof(int)   /* cell_idx */
        + 2 * N * sizeof(int)   /* vert_idx */
    );

    BenchData bd;
    bd.alloc(N, nlev);
    bd.fill(nlev);

    for (int di = 0; di < 4; di++) {
        CellDist dist = (CellDist)di;
        gen_cell_idx_logical(cell_logical, N, dist, rng);

        for (int V = 1; V <= 4; V++) {
            bd.set_variant(V, cell_logical, vd.logical);

            /* warmup */
            for (int r = 0; r < WARMUP; r++) {
                flush_caches();
                kern_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
                    bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
                    bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
            }

            /* timed runs */
            double sum_ms = 0;
            for (int r = 0; r < NRUNS; r++) {
                flush_caches();
                auto t0 = std::chrono::high_resolution_clock::now();
                kern_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
                    bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
                    bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
                auto t1 = std::chrono::high_resolution_clock::now();
                double dt = std::chrono::duration<double, std::milli>(t1-t0).count();
                sum_ms += dt;
            }

            double avg_ms = sum_ms / NRUNS;
            double bw_gbs = (total_bytes / (avg_ms * 1e-3)) / 1e9;
            printf("nlev=%d  dist=%-12s  V%d  avg=%.3f ms  BW=%.1f GB/s\n",
                    nlev, dist_name[di], V, avg_ms, bw_gbs);
        }
    }
    bd.free_all();

    vd.free_all();
    delete[] cell_logical;
    return 0;
}