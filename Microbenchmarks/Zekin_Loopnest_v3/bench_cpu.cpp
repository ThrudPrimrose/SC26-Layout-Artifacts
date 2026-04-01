/*
 * bench_cpu.cpp -- NUMA-aware CPU-only z_v_grad_w stencil benchmark
 *
 * Now supports a 5th cell distribution "exact" loaded from ICON serialised
 * p_patch files.
 *
 * Environment:
 *   ICON_DATA_PATH  - directory containing p_patch.*.data files
 *                     default: /home/primrose/Work/icon-artifacts/velocity/data_r02b05
 *   Timestep defaults to 9.  Pass a different step as argv[1] (integer).
 *
 * Compile:  g++ -O3 -fopenmp -march=native -std=c++17 bench_cpu.cpp -o bench_cpu
 *
 * Run:
 *   ./bench_cpu                  # synthetic + exact (step 9)
 *   ./bench_cpu 43               # synthetic + exact (step 43)
 *   ICON_DATA_PATH=/my/data ./bench_cpu 1
 */

#include "bench_common.h"
#include "icon_data_loader.h"
#include <omp.h>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <utility>

/* ================================================================ */
/*  Map variant + parallelization -> SchedKind                       */
/* ================================================================ */
static SchedKind sched_for_par_for(int V) {
    return (V <= 2) ? SCHED_JK_OUTER : SCHED_JE_OUTER;
}

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
/*  CPU reference (serial, for verification)                         */
/* ================================================================ */
template<int V>
static void cpu_reference(
    double* __restrict__ out,
    const double* __restrict__ vn_ie,   const double* __restrict__ inv_dual,
    const double* __restrict__ w,       const int*    __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int*    __restrict__ vert_idx, int N, int nlev)
{
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) { STENCIL_BODY(V) }
}

static void cpu_reference_v(int V,
    double* out, const double* vn_ie, const double* inv_dual,
    const double* w, const int* cell_idx,
    const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v,
    const int* vert_idx, int N, int nlev)
{
    switch (V) {
        case 1: cpu_reference<1>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 2: cpu_reference<2>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 3: cpu_reference<3>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
        case 4: cpu_reference<4>(out,vn_ie,inv_dual,w,cell_idx,z_vt_ie,inv_primal,tangent,z_w_v,vert_idx,N,nlev); break;
    }
}

/* ================================================================ */
/*  numerical verification                                           */
/* ================================================================ */
static bool verify(const double* got, const double* ref, size_t n,
                   double rtol, double atol,
                   int* n_fail, double* max_rel)
{
    *n_fail  = 0;
    *max_rel = 0.0;
    for (size_t i = 0; i < n; i++) {
        double diff = std::abs(got[i] - ref[i]);
        double denom = std::max(std::abs(ref[i]), 1e-300);
        double rel = diff / denom;
        if (rel > *max_rel) *max_rel = rel;
        if (diff > atol + rtol * std::abs(ref[i]))
            (*n_fail)++;
    }
    return *n_fail == 0;
}

/* ================================================================ */
/*  Cache-flush: 2-D Jacobi on NUMA-aware persistent buffers         */
/* ================================================================ */

static constexpr int FLUSH_N = 8192 * 4;
static constexpr int FLUSH_STEPS = 3;

static double* flush_buf0 = nullptr;
static double* flush_buf1 = nullptr;

static void flush_caches()
{
    static bool inited = false;

    if (!inited) {
        size_t n = (size_t)FLUSH_N * FLUSH_N;
        flush_buf0 = numa_alloc_unfaulted<double>(n);
        flush_buf1 = numa_alloc_unfaulted<double>(n);
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < n; i++) {
            uint64_t h = splitmix64(12345ULL + (uint64_t)i);
            flush_buf0[i] = (double)(h >> 11) / (double)(1ULL << 53);
            flush_buf1[i] = flush_buf0[i];
        }
        inited = true;
    }

    double *A = flush_buf0;
    double *B = flush_buf1;

    for (int s = 0; s < FLUSH_STEPS; s++) {
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < FLUSH_N - 1; i++)
            for (int j = 1; j < FLUSH_N - 1; j++)
                B[i * FLUSH_N + j] = 0.25 * (
                    A[(i-1)*FLUSH_N + j] + A[(i+1)*FLUSH_N + j] +
                    A[i*FLUSH_N + (j-1)] + A[i*FLUSH_N + (j+1)]);
        std::swap(A, B);
    }

    int ri = rand() % (FLUSH_N * FLUSH_N);
    printf("  [flush] A[%d] = %.12e\n", ri, A[ri]);
}

/* ================================================================ */
/*  Helper: run omp_for + collapse2 for one (N, nlev, dist, V)      */
/* ================================================================ */
static void run_variant_cpu(
    FILE* fcsv,
    int V, int N, int nlev, const char* dist_label,
    BenchData& bd, double* h_ref)
{
    /* ---- omp parallel for ---- */
    SchedKind par_sched = sched_for_par_for(V);
    bd.change_schedule(par_sched);

    /* serial reference */
    cpu_reference_v(V,
        h_ref, bd.h_vn_ie, bd.inv_dual,
        bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
        bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
    flush_caches();

    /* warmup omp_for */
    for (int r = 0; r < WARMUP; r++) {
        flush_caches();
        cpu_par_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
            bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
            bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
    }
    flush_caches();

    /* verify omp_for */
    {
        int n_fail = 0; double max_rel = 0.0;
        bool ok = verify(bd.h_out, h_ref, bd.sz2d,
                         1e-12, 1e-15, &n_fail, &max_rel);
        if (!ok)
            printf("VERIFY FAIL: nlev=%d dist=%-12s V=%d "
                   "omp_for  fails=%d max_rel=%.3e\n",
                   nlev, dist_label, V, n_fail, max_rel);
        else
            printf("VERIFY OK:   nlev=%d dist=%-12s V=%d "
                   "omp_for  max_rel=%.3e\n",
                   nlev, dist_label, V, max_rel);
    }

    /* timed omp_for */
    for (int r = 0; r < NRUNS; r++) {
        flush_caches();
        auto t0 = std::chrono::high_resolution_clock::now();
        cpu_par_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
            bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
            bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
        auto t1 = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration<double, std::milli>(t1-t0).count();
        fprintf(fcsv, "cpu,%d,%d,%d,%s,omp_for,%d,%.9f\n",
                V, nlev, N, dist_label, r, dt);
        flush_caches();
    }

    /* ---- omp collapse(2) ---- */
    bd.change_schedule(SCHED_COLLAPSE2);

    /* warmup collapse2 */
    for (int r = 0; r < WARMUP; r++) {
        flush_caches();
        cpu_col_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
            bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
            bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
        flush_caches();
    }

    /* verify collapse2 */
    {
        int n_fail = 0; double max_rel = 0.0;
        bool ok = verify(bd.h_out, h_ref, bd.sz2d,
                         1e-12, 1e-15, &n_fail, &max_rel);
        if (!ok)
            printf("VERIFY FAIL: nlev=%d dist=%-12s V=%d "
                   "collapse2  fails=%d max_rel=%.3e\n",
                   nlev, dist_label, V, n_fail, max_rel);
        else
            printf("VERIFY OK:   nlev=%d dist=%-12s V=%d "
                   "collapse2  max_rel=%.3e\n",
                   nlev, dist_label, V, max_rel);
    }
    flush_caches();

    /* timed collapse2 */
    for (int r = 0; r < NRUNS; r++) {
        flush_caches();
        auto t0 = std::chrono::high_resolution_clock::now();
        cpu_col_tbl[V-1](bd.h_out, bd.h_vn_ie, bd.inv_dual,
            bd.h_w, bd.h_cidx, bd.h_z_vt_ie, bd.inv_primal,
            bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
        auto t1 = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration<double, std::milli>(t1-t0).count();
        fprintf(fcsv, "cpu,%d,%d,%d,%s,omp_collapse2,%d,%.9f\n",
                V, nlev, N, dist_label, r, dt);
        flush_caches();
    }
    flush_caches();

    printf("Done: nlev=%d  dist=%-12s  V=%d\n",
           nlev, dist_label, V);
}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main(int argc, char* argv[]) {
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

    /* ---- load ICON exact data ----
     *
     * ICON_DATA_PATH env var -> directory (default: the primrose path)
     * argv[1] -> timestep (default: 9)
     */
    int icon_step = 9;
    if (argc >= 2) icon_step = atoi(argv[1]);

    std::string patch_path = icon_patch_path(icon_step);
    printf("Loading ICON data from: %s\n", patch_path.c_str());

    IconEdgeData icon_ed;
    bool have_exact = icon_load_patch(patch_path.c_str(), icon_ed);
    if (have_exact) {
        printf("ICON exact data loaded: n_edges_valid=%d  n_edges_alloc=%d  "
               "n_cells=%d  n_verts=%d  nproma=%d\n",
               icon_ed.n_edges_valid, icon_ed.n_edges_alloc,
               icon_ed.n_cells, icon_ed.n_verts, icon_ed.nproma);
    } else {
        fprintf(stderr, "WARNING: failed to load '%s', "
                "will run synthetic distributions only\n",
                patch_path.c_str());
    }

    printf("OMP threads: %d\n", omp_get_max_threads());

    srand((unsigned)time(NULL));

    /* fault in flush buffers before any timing */
    flush_caches();
    printf("Flush buffers initialized (2 x %.0f MB)\n",
           (double)FLUSH_N * FLUSH_N * 8 / 1e6);

    for (int nlev_i = 0; nlev_i < N_NLEVS; nlev_i++) {
        int nlev = NLEVS[nlev_i];

        /* ============================================================ */
        /*  Synthetic distributions (uniform, normal1, normal4, ...)    */
        /* ============================================================ */
        {
            BenchData bd;
            bd.alloc(N, nlev);
            bd.fill(nlev);

            double* h_ref = numa_alloc_unfaulted<double>(bd.sz2d);
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < bd.sz2d; i++)
                h_ref[i] = 0.0;

            for (int di = 0; di < 4; di++) {
                CellDist dist = (CellDist)di;
                gen_cell_idx_logical(cell_logical, N, dist, rng);

                for (int V = 1; V <= 4; V++) {
                    SchedKind par_sched = sched_for_par_for(V);
                    bd.set_variant(V, cell_logical, vd.logical, par_sched);

                    run_variant_cpu(fcsv, V, N, nlev, dist_name[di],
                                   bd, h_ref);
                    fflush(fcsv);
                }
            }
            numa_dealloc(h_ref, bd.sz2d);
            bd.free_all();
        }

        /* ============================================================ */
        /*  EXACT distribution from ICON grid data                      */
        /* ============================================================ */
        if (have_exact) {
            const int Ne = icon_ed.n_edges_valid;

            printf("\n=== EXACT distribution: nlev=%d  Ne=%d ===\n", nlev, Ne);

            if (icon_ed.n_cells > Ne || icon_ed.n_verts > Ne) {
                fprintf(stderr, "WARNING: n_cells=%d or n_verts=%d > Ne=%d; "
                        "skipping exact for nlev=%d\n",
                        icon_ed.n_cells, icon_ed.n_verts, Ne, nlev);
                continue;
            }

            BenchData bd;
            bd.alloc(Ne, nlev);
            bd.fill(nlev);

            double* h_ref = numa_alloc_unfaulted<double>(bd.sz2d);
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < bd.sz2d; i++)
                h_ref[i] = 0.0;

            /* Build logical connectivity arrays from icon_ed */
            int* exact_cell_logical = new int[Ne * 2];
            int* exact_vert_logical = new int[Ne * 2];
            for (int i = 0; i < Ne * 2; i++) {
                exact_cell_logical[i] = icon_ed.cell_idx[i];
                exact_vert_logical[i] = icon_ed.vert_idx[i];
            }

            /* Override 1-D geometry with real ICON values */
            for (int je = 0; je < Ne; je++) {
                bd.inv_dual[je]   = icon_ed.inv_dual[je];
                bd.inv_primal[je] = icon_ed.inv_primal[je];
                bd.tangent_o[je]  = icon_ed.tangent_o[je];
            }

            for (int V = 1; V <= 4; V++) {
                SchedKind par_sched = sched_for_par_for(V);
                bd.set_variant(V, exact_cell_logical, exact_vert_logical,
                               par_sched);

                run_variant_cpu(fcsv, V, Ne, nlev, "exact", bd, h_ref);
                fflush(fcsv);
            }

            delete[] exact_cell_logical;
            delete[] exact_vert_logical;
            numa_dealloc(h_ref, bd.sz2d);
            bd.free_all();
        }
    }

    /* clean up flush buffers */
    numa_dealloc(flush_buf0, (size_t)FLUSH_N * FLUSH_N);
    numa_dealloc(flush_buf1, (size_t)FLUSH_N * FLUSH_N);

    vd.free_all();
    delete[] cell_logical;
    if (have_exact) icon_ed.free_all();
    fclose(fcsv);

    printf("\nResults written to z_v_grad_w_cpu.csv\n");
    return 0;
}