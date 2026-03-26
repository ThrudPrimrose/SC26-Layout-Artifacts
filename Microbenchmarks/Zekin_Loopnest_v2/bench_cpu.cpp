/*
 * bench_cpu.cpp -- CPU-only z_v_grad_w stencil benchmark
 *
 * Compile:  g++ -O3 -fopenmp -march=native -std=c++17 bench_cpu.cpp -o bench_cpu
 */

#include "bench_common.h"
#include <omp.h>
#include <cstdlib>
#include <cstring>
#include <utility>

/* ================================================================ */
/*  CPU kernels                                                      */
/* ================================================================ */

template<int V>
static void cpu_par_for(
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

template<int V>
static void cpu_collapse2(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev)
{
    if constexpr (V <= 2) {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) { STENCIL_BODY(V) }
    } else {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int je = 0; je < N; je++)
            for (int jk = 0; jk < nlev; jk++) { STENCIL_BODY(V) }
    }
}

/* dispatch tables */
typedef void (*cpu_fn_t)(double*, const double*, const double*,
    const double*, const int*, const double*, const double*,
    const double*, const double*, const int*, int, int);

static cpu_fn_t cpu_par_tbl[] = {
    cpu_par_for<1>, cpu_par_for<2>, cpu_par_for<3>, cpu_par_for<4>
};
static cpu_fn_t cpu_col_tbl[] = {
    cpu_collapse2<1>, cpu_collapse2<2>, cpu_collapse2<3>, cpu_collapse2<4>
};

/* ================================================================ */
/*  Cache-flush: 2-D Jacobi stencil                                  */
/* ================================================================ */

static constexpr int FLUSH_N = 2048;          /* ~32 MB per grid */
static constexpr int FLUSH_STEPS = 3;

static double flush_buf0[FLUSH_N * FLUSH_N];
static double flush_buf1[FLUSH_N * FLUSH_N];

static void flush_caches()
{
    static bool inited = false;
    double *A = flush_buf0, *B = flush_buf1;

    if (!inited) {
        srand(12345);
        for (int i = 0; i < FLUSH_N * FLUSH_N; i++)
            A[i] = (double)rand() / RAND_MAX;
        std::memcpy(B, A, sizeof(flush_buf0));
        inited = true;
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

    /* print one random element so the compiler cannot elide the work */
    int ri = rand() % (FLUSH_N * FLUSH_N);
    printf("  [flush] A[%d] = %.12e\n", ri, A[ri]);
}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main() {
    FILE* fcsv = fopen("z_v_grad_w_cpu.csv", "w");
    if (!fcsv) { perror("fopen"); return 1; }
    fprintf(fcsv,
        "backend,variant,nlev,nproma,cell_dist,"
        "parallelization,run_id,time_ms\n");

    const int N = NPROMA;
    std::mt19937 rng(42);

    VertData vd;
    vd.init(N, rng);

    int* cell_logical = new int[N * 2];

    printf("OMP threads: %d\n", omp_get_max_threads());

    for (int nlev_i = 0; nlev_i < N_NLEVS; nlev_i++) {
        int nlev = NLEVS[nlev_i];

        BenchData bd;
        bd.alloc(N, nlev);
        bd.fill(nlev);

        for (int di = 0; di < 4; di++) {
            CellDist dist = (CellDist)di;
            gen_cell_idx_logical(cell_logical, N, dist, rng);

            for (int V = 1; V <= 4; V++) {
                bd.set_variant(V, cell_logical, vd.logical);

                /* ---- omp parallel for ---- */
                for (int r = 0; r < WARMUP; r++)
                    cpu_par_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
                        bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
                        bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);

                for (int r = 0; r < NRUNS; r++) {
                    flush_caches();
                    auto t0 = std::chrono::high_resolution_clock::now();
                    cpu_par_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
                        bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
                        bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    double dt = std::chrono::duration<double, std::milli>(t1-t0).count();
                    fprintf(fcsv, "cpu,%d,%d,%d,%s,omp_for,%d,%.9f\n",
                            V, nlev, N, dist_name[di], r, dt);
                }

                /* ---- omp collapse(2) ---- */
                for (int r = 0; r < WARMUP; r++)
                    cpu_col_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
                        bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
                        bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);

                for (int r = 0; r < NRUNS; r++) {
                    flush_caches();
                    auto t0 = std::chrono::high_resolution_clock::now();
                    cpu_col_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
                        bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
                        bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    double dt = std::chrono::duration<double, std::milli>(t1-t0).count();
                    fprintf(fcsv, "cpu,%d,%d,%d,%s,omp_collapse2,%d,%.9f\n",
                            V, nlev, N, dist_name[di], r, dt);
                }

                printf("Done: nlev=%d  dist=%-12s  V=%d\n",
                       nlev, dist_name[di], V);
                fflush(fcsv);
            }
        }
        bd.free_all();
    }

    vd.free_all();
    delete[] cell_logical;
    fclose(fcsv);

    printf("\nResults written to z_v_grad_w_cpu.csv\n");
    printf("Total rows: 2 nlevs x 4 dists x 4 variants x 2 omp x %d = %d\n",
           NRUNS, 2 * 4 * 4 * 2 * NRUNS);
    return 0;
}