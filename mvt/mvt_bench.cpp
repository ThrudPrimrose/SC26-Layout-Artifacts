// mvt_bench.cpp — MVT layout/schedule/blocking/tiling explorer (OpenMP integrated)
// PolyBench mvt: y1[i] += A[i][j]*x1[j],  y2[j] += A[i][j]*x2[i]  (A is N×N)
//
// Build: g++ -O3 -fopenmp -std=c++17 -march=native -ffast-math \
//        -DSZ_N=1024 -DSZ_B=32 -DSZ_T=32 -DNRUNS=50 -DNWARM=5 \
//        -o mvt_bench mvt_bench.cpp -lopenblas
//
// Run:   ./mvt_bench <mode:1-6> [output.csv]
//  1 = A layout variants (row-major, col-major)            — 2 variants
//  2 = loop orderings (4 combined y1+y2 orderings)         — 4 variants
//  3 = blocked A (4 outer×inner layout combos) + blocked   — 4 variants
//  4 = tiled loops, row-major A                            — 1 variant
//  5 = OpenBLAS (2× cblas_dgemv)                           — 1 variant
//  6 = FUSED: single pass over A, two variants                — 5 variants
//      fused_tiled             : row-major A, SZ_T×SZ_T tiling, heap-private y2
//      fused_blk_{RR,RC,CR,CC} : blocked A (4 outer×inner combos),
//                                heap-private y2
//
// Fused strategy: each thread accumulates y1[i] into a scalar and y2[j]
// into a thread-private heap buffer y2_priv[N].  After the parallel for a
// second parallel for reduces all per-thread copies into y2[].  No atomics.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include <cblas.h>

// ── compile-time parameters ──────────────────────────────────────────────
#ifndef SZ_N
#define SZ_N 1024
#endif
#ifndef SZ_B           // block size (square) for mode 3
#define SZ_B 32
#endif
#ifndef SZ_T           // tile size (square) for mode 4
#define SZ_T 32
#endif
#ifndef NRUNS
#define NRUNS 50
#endif
#ifndef NWARM
#define NWARM 5
#endif

static constexpr int N     = SZ_N;
static constexpr double alpha = 1.0;
using Clock = std::chrono::high_resolution_clock;

// ── fast parallel deterministic init ────────────────────────────────────
// xorshift32, seeded per-thread: reproducible across runs, parallelisable
static void init_rand(double *a, long n) {
    #pragma omp parallel
    {
        uint32_t s = 12345u + (uint32_t)omp_get_thread_num() * 2654435761u;
        #pragma omp for schedule(static)
        for (long i = 0; i < n; i++) {
            s ^= s << 13; s ^= s >> 17; s ^= s << 5;
            a[i] = (s & 0xFFFFu) * (1.0 / 65536.0);
        }
    }
}

// ── verification ─────────────────────────────────────────────────────────
static double max_relerr(const double *got, const double *ref, int n) {
    double mx = 0.0;
    for (int i = 0; i < n; i++)
        mx = std::max(mx, std::fabs(got[i] - ref[i]) / (1.0 + std::fabs(ref[i])));
    return mx;
}

// ── flat row/col-major index helpers ────────────────────────────────────
enum Lay { R, C };
template<Lay L>
inline int fidx(int r, int c, int nr, int nc) {
    return (L == R) ? r * nc + c : c * nr + r;
}

// convert row-major source to target layout
static void to_lay(const double *rm, double *dst, int nr, int nc, Lay lay) {
    if (lay == R) { memcpy(dst, rm, (size_t)nr * nc * sizeof(double)); return; }
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++)
            dst[j * nr + i] = rm[i * nc + j];
}
static void to_rm(const double *src, double *dst, int nr, int nc, Lay lay) {
    if (lay == R) { memcpy(dst, src, (size_t)nr * nc * sizeof(double)); return; }
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++)
            dst[i * nc + j] = src[j * nr + i];
}

// ── blocked layout index ─────────────────────────────────────────────────
// outer_col: block grid stored col-major  (block (bi,bj) offset = bj*NB+bi)
// inner_col: within-block elements stored col-major (row r, col c → c*B+r)
static constexpr int NB = (N + SZ_B - 1) / SZ_B;        // number of blocks per dim
static constexpr int BSIZ = NB * NB * SZ_B * SZ_B;       // padded allocation size

template<bool outer_col, bool inner_col>
inline int blkidx(int i, int j) {
    const int bi = i / SZ_B, bj = j / SZ_B;
    const int r  = i % SZ_B, c  = j % SZ_B;
    const int block_off = outer_col ? (bj * NB + bi) : (bi * NB + bj);
    const int inner_off = inner_col ? (c * SZ_B + r)    : (r * SZ_B + c);
    return block_off * SZ_B * SZ_B + inner_off;
}

static void to_blk_rr(const double *rm, double *blk) {
    memset(blk, 0, BSIZ * sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            blk[blkidx<false,false>(i,j)] = rm[i*N+j];
}
static void to_blk_rc(const double *rm, double *blk) {
    memset(blk, 0, BSIZ * sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            blk[blkidx<false,true>(i,j)] = rm[i*N+j];
}
static void to_blk_cr(const double *rm, double *blk) {
    memset(blk, 0, BSIZ * sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            blk[blkidx<true,false>(i,j)] = rm[i*N+j];
}
static void to_blk_cc(const double *rm, double *blk) {
    memset(blk, 0, BSIZ * sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            blk[blkidx<true,true>(i,j)] = rm[i*N+j];
}

// ── MODE 1: A layout variants ────────────────────────────────────────────
// y1[i] += A[i][j]*x1[j]   y2[j] += A[i][j]*x2[i]
// Parallel: y1 parallel on i (row owner), y2 parallel on j (col owner)

template<Lay LA>
static void kern_lay(const double *Ad,
                     const double *x1, const double *x2,
                     double *y1, double *y2)
{
    // y1 = A * x1  — parallel on i, inner j is read-only on y1
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double acc = 0.0;
        for (int j = 0; j < N; j++)
            acc += alpha * Ad[fidx<LA>(i, j, N, N)] * x1[j];
        y1[i] += acc;
    }
    // y2 = A^T * x2 — parallel on j (= column of A), inner i read-only on y2
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++) {
        double acc = 0.0;
        for (int i = 0; i < N; i++)
            acc += alpha * Ad[fidx<LA>(i, j, N, N)] * x2[i];
        y2[j] += acc;
    }
}

static void run_layout(Lay lay,
                       const double *A_rm,
                       const double *x1, const double *x2,
                       const double *y10, const double *y20,
                       const double *y1ref, const double *y2ref,
                       FILE *fp, const char *tag)
{
    double *Ad  = new double[N * N];
    double *y1w = new double[N];
    double *y2w = new double[N];
    to_lay(A_rm, Ad, N, N, lay);

    for (int r = 0; r < NRUNS; r++) {
        memcpy(y1w, y10, N * sizeof(double));
        memcpy(y2w, y20, N * sizeof(double));
        auto t0 = Clock::now();
        if (lay == R) kern_lay<R>(Ad, x1, x2, y1w, y2w);
        else          kern_lay<C>(Ad, x1, x2, y1w, y2w);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,1,%s,0,0,%d,%.2f\n",
                    N, tag, r - NWARM,
                    std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    double e1 = max_relerr(y1w, y1ref, N);
    double e2 = max_relerr(y2w, y2ref, N);
    fprintf(stderr, "%s %s: y1_err=%.2e y2_err=%.2e\n",
            (e1 > 1e-9 || e2 > 1e-9) ? "FAIL" : "OK", tag, e1, e2);
    delete[] Ad; delete[] y1w; delete[] y2w;
}

// ── MODE 2: loop orderings (row-major A) ─────────────────────────────────
// Four combined orderings for (y1 loop, y2 loop)
// y1_ij: for i { for j { y1[i] += A[i][j]*x1[j] } }  → parallel i (safe)
// y1_ji: for j { for i { y1[i] += A[i][j]*x1[j] } }  → serial outer
// y2_ij: for i { for j { y2[j] += A[i][j]*x2[i] } }  → serial outer
// y2_ji: for j { for i { y2[j] += A[i][j]*x2[i] } }  → parallel j (safe)

#define RM(a,r,c) (a)[(r)*N+(c)]

static void kern_ij_ji(const double *A,
                       const double *x1, const double *x2,
                       double *y1, double *y2)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double acc = 0.0;
        for (int j = 0; j < N; j++) acc += alpha * RM(A,i,j) * x1[j];
        y1[i] += acc;
    }
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++) {
        double acc = 0.0;
        for (int i = 0; i < N; i++) acc += alpha * RM(A,i,j) * x2[i];
        y2[j] += acc;
    }
}

static void kern_ji_ij(const double *A,
                       const double *x1, const double *x2,
                       double *y1, double *y2)
{
    // y1: outer j serial (write to multiple y1[i] inside)
    for (int j = 0; j < N; j++)
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++)
            y1[i] += alpha * RM(A,i,j) * x1[j];
    // y2: outer i serial (write to multiple y2[j] inside)
    for (int i = 0; i < N; i++)
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < N; j++)
            y2[j] += alpha * RM(A,i,j) * x2[i];
}

static void kern_ij_ij(const double *A,
                       const double *x1, const double *x2,
                       double *y1, double *y2)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        double acc = 0.0;
        for (int j = 0; j < N; j++) acc += alpha * RM(A,i,j) * x1[j];
        y1[i] += acc;
    }
    for (int i = 0; i < N; i++)
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < N; j++)
            y2[j] += alpha * RM(A,i,j) * x2[i];
}

static void kern_ji_ji(const double *A,
                       const double *x1, const double *x2,
                       double *y1, double *y2)
{
    for (int j = 0; j < N; j++)
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++)
            y1[i] += alpha * RM(A,i,j) * x1[j];
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++) {
        double acc = 0.0;
        for (int i = 0; i < N; i++) acc += alpha * RM(A,i,j) * x2[i];
        y2[j] += acc;
    }
}

using LoopKern = void(*)(const double*, const double*, const double*, double*, double*);

static void run_loop(LoopKern fn,
                     const double *A,
                     const double *x1, const double *x2,
                     const double *y10, const double *y20,
                     const double *y1ref, const double *y2ref,
                     FILE *fp, const char *tag)
{
    double *y1w = new double[N];
    double *y2w = new double[N];
    for (int r = 0; r < NRUNS; r++) {
        memcpy(y1w, y10, N * sizeof(double));
        memcpy(y2w, y20, N * sizeof(double));
        auto t0 = Clock::now();
        fn(A, x1, x2, y1w, y2w);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,2,%s,0,0,%d,%.2f\n",
                    N, tag, r - NWARM,
                    std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    double e1 = max_relerr(y1w, y1ref, N);
    double e2 = max_relerr(y2w, y2ref, N);
    fprintf(stderr, "%s %s: y1_err=%.2e y2_err=%.2e\n",
            (e1 > 1e-9 || e2 > 1e-9) ? "FAIL" : "OK", tag, e1, e2);
    delete[] y1w; delete[] y2w;
}

// ── MODE 3: blocked A + blocked loops ───────────────────────────────────
// Four storage variants: outer {R,C} × inner {R,C}
// y1: parallel on block-row bi (each bi owns y1[bi*SZ_B .. bi*SZ_B+SZ_B-1])
// y2: parallel on block-col bj (each bj owns y2[bj*SZ_B .. bj*SZ_B+SZ_B-1])

template<bool outer_col, bool inner_col>
static void kern_blocked(const double *Ab,
                         const double *x1, const double *x2,
                         double *y1, double *y2)
{
    // y1 = A * x1
    #pragma omp parallel for schedule(static)
    for (int bi = 0; bi < NB; bi++) {
        const int ie = std::min(bi * SZ_B + SZ_B, N);
        for (int bj = 0; bj < NB; bj++) {
            const int je = std::min(bj * SZ_B + SZ_B, N);
            for (int i = bi * SZ_B; i < ie; i++) {
                double acc = 0.0;
                for (int j = bj * SZ_B; j < je; j++)
                    acc += alpha * Ab[blkidx<outer_col, inner_col>(i, j)] * x1[j];
                y1[i] += acc;
            }
        }
    }
    // y2 = A^T * x2
    #pragma omp parallel for schedule(static)
    for (int bj = 0; bj < NB; bj++) {
        const int je = std::min(bj * SZ_B + SZ_B, N);
        for (int bi = 0; bi < NB; bi++) {
            const int ie = std::min(bi * SZ_B + SZ_B, N);
            for (int j = bj * SZ_B; j < je; j++) {
                double acc = 0.0;
                for (int i = bi * SZ_B; i < ie; i++)
                    acc += alpha * Ab[blkidx<outer_col, inner_col>(i, j)] * x2[i];
                y2[j] += acc;
            }
        }
    }
}

struct BlkVariant {
    const char *tag;
    void (*to_blk)(const double*, double*);
    void (*kern)(const double*, const double*, const double*, double*, double*);
};

static const BlkVariant BLK_VARIANTS[] = {
    { "blk_RR", to_blk_rr, kern_blocked<false,false> },
    { "blk_RC", to_blk_rc, kern_blocked<false,true>  },
    { "blk_CR", to_blk_cr, kern_blocked<true, false> },
    { "blk_CC", to_blk_cc, kern_blocked<true, true>  },
};

static void run_blocked(const BlkVariant &v,
                        const double *A_rm,
                        const double *x1, const double *x2,
                        const double *y10, const double *y20,
                        const double *y1ref, const double *y2ref,
                        FILE *fp)
{
    double *Ab  = new double[BSIZ];
    double *y1w = new double[N];
    double *y2w = new double[N];
    v.to_blk(A_rm, Ab);

    for (int r = 0; r < NRUNS; r++) {
        memcpy(y1w, y10, N * sizeof(double));
        memcpy(y2w, y20, N * sizeof(double));
        auto t0 = Clock::now();
        v.kern(Ab, x1, x2, y1w, y2w);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,3,%s,%d,0,%d,%.2f\n",
                    N, v.tag, SZ_B, r - NWARM,
                    std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    double e1 = max_relerr(y1w, y1ref, N);
    double e2 = max_relerr(y2w, y2ref, N);
    fprintf(stderr, "%s %s SZ_B=%d: y1_err=%.2e y2_err=%.2e\n",
            (e1 > 1e-9 || e2 > 1e-9) ? "FAIL" : "OK", v.tag, SZ_B, e1, e2);
    delete[] Ab; delete[] y1w; delete[] y2w;
}

// ── MODE 4: tiled loops, row-major A ────────────────────────────────────
static void kern_tiled(const double *A,
                       const double *x1, const double *x2,
                       double *y1, double *y2)
{
    // y1 = A * x1  — parallel on tile-rows
    #pragma omp parallel for schedule(static)
    for (int ii = 0; ii < N; ii += SZ_T) {
        const int ie = std::min(ii + SZ_T, N);
        for (int jj = 0; jj < N; jj += SZ_T) {
            const int je = std::min(jj + SZ_T, N);
            for (int i = ii; i < ie; i++) {
                double acc = 0.0;
                for (int j = jj; j < je; j++)
                    acc += alpha * RM(A,i,j) * x1[j];
                y1[i] += acc;
            }
        }
    }
    // y2 = A^T * x2 — parallel on tile-cols
    #pragma omp parallel for schedule(static)
    for (int jj = 0; jj < N; jj += SZ_T) {
        const int je = std::min(jj + SZ_T, N);
        for (int ii = 0; ii < N; ii += SZ_T) {
            const int ie = std::min(ii + SZ_T, N);
            for (int j = jj; j < je; j++) {
                double acc = 0.0;
                for (int i = ii; i < ie; i++)
                    acc += alpha * RM(A,i,j) * x2[i];
                y2[j] += acc;
            }
        }
    }
}

static void run_tiled(const double *A,
                      const double *x1, const double *x2,
                      const double *y10, const double *y20,
                      const double *y1ref, const double *y2ref,
                      FILE *fp)
{
    double *y1w = new double[N];
    double *y2w = new double[N];
    for (int r = 0; r < NRUNS; r++) {
        memcpy(y1w, y10, N * sizeof(double));
        memcpy(y2w, y20, N * sizeof(double));
        auto t0 = Clock::now();
        kern_tiled(A, x1, x2, y1w, y2w);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,4,tiled,0,%d,%d,%.2f\n",
                    N, SZ_T, r - NWARM,
                    std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    double e1 = max_relerr(y1w, y1ref, N);
    double e2 = max_relerr(y2w, y2ref, N);
    fprintf(stderr, "%s tiled SZ_T=%d: y1_err=%.2e y2_err=%.2e\n",
            (e1 > 1e-9 || e2 > 1e-9) ? "FAIL" : "OK", SZ_T, e1, e2);
    delete[] y1w; delete[] y2w;
}

// ── MODE 5: OpenBLAS (2× dgemv) ─────────────────────────────────────────
// OpenBLAS has no dedicated MVT; use dgemv twice (row-major A)
static void run_openblas(const double *A,
                         const double *x1, const double *x2,
                         const double *y10, const double *y20,
                         const double *y1ref, const double *y2ref,
                         FILE *fp)
{
    double *y1w = new double[N];
    double *y2w = new double[N];
    for (int r = 0; r < NRUNS; r++) {
        memcpy(y1w, y10, N * sizeof(double));
        memcpy(y2w, y20, N * sizeof(double));
        auto t0 = Clock::now();
        // y1 += alpha * A * x1
        cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N,
                    alpha, A, N, x1, 1, 1.0, y1w, 1);
        // y2 += alpha * A^T * x2
        cblas_dgemv(CblasRowMajor, CblasTrans, N, N,
                    alpha, A, N, x2, 1, 1.0, y2w, 1);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,5,openblas,0,0,%d,%.2f\n",
                    N, r - NWARM,
                    std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    double e1 = max_relerr(y1w, y1ref, N);
    double e2 = max_relerr(y2w, y2ref, N);
    fprintf(stderr, "%s openblas: y1_err=%.2e y2_err=%.2e\n",
            (e1 > 1e-9 || e2 > 1e-9) ? "FAIL" : "OK", e1, e2);
    delete[] y1w; delete[] y2w;
}

// ── MODE 6: fused kernels — single pass over A ──────────────────────────
//
// Both variants read each A[i][j] exactly once and accumulate into y1[i]
// (via a scalar) and y2[j] (via a thread-private heap buffer) simultaneously.
//
// Heap-private buffer discipline
// ──────────────────────────────
// y2_priv[nt][N] is allocated once in run_fused before the timing loop and
// re-zeroed between runs.  Kernels receive it as a parameter — no allocation
// or probing inside the hot path.  After each kernel call a parallel reduce
// accumulates the per-thread copies into y2[].
//
// Parallelism: collapse(2) schedule(static) over the two outer tile/block
// loops.  ie/je (or ie/je for blocks) are computed inside the loop body so
// that the collapsed loop nest remains canonical.  y2p is looked up by
// thread id inside each iteration.
//
// Two variants
// ─────────────
//  fused_tiled           : row-major A, SZ_T×SZ_T tiling
//  fused_blk_{RR,RC,CR,CC}: blocked A (4 outer×inner layout combos)

// ── helpers: manage heap-private y2 buffers ─────────────────────────────
static int y2_priv_alloc(double **&priv) {
    int nt = 0;
    #pragma omp parallel reduction(max:nt)
    { nt = std::max(nt, omp_get_num_threads()); }
    priv = new double*[nt];
    for (int t = 0; t < nt; t++)
        priv[t] = new double[N];
    return nt;
}
static void y2_priv_zero(double **priv, int nt) {
    for (int t = 0; t < nt; t++)
        memset(priv[t], 0, N * sizeof(double));
}
static void y2_priv_reduce(double **priv, int nt, double *y2) {
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++) {
        double s = 0.0;
        for (int t = 0; t < nt; t++) s += priv[t][j];
        y2[j] += s;
    }
}
static void y2_priv_free(double **priv, int nt) {
    for (int t = 0; t < nt; t++) delete[] priv[t];
    delete[] priv;
}

// ── 6a: fused tiled, row-major A ─────────────────────────────────────────
// collapse(2) distributes (N/SZ_T)² tile pairs across threads.
// ie and je are derived inside the body from the collapsed loop variables.
static void kern_fused_tiled(const double *A,
                              const double *x1, const double *x2,
                              double *y1, double **y2_priv)
{
    #pragma omp parallel for collapse(2) schedule(static)
    for (int ii = 0; ii < N; ii += SZ_T) {
        for (int jj = 0; jj < N; jj += SZ_T) {
            const int ie = std::min(ii + SZ_T, N);
            const int je = std::min(jj + SZ_T, N);
            double *y2p = y2_priv[omp_get_thread_num()];
            for (int i = ii; i < ie; i++) {
                const double *Ai = A + (long)i * N + jj;
                double acc = 0.0;
                const double xi2 = alpha * x2[i];
                for (int j = 0; j < je - jj; j++) {
                    const double aij = Ai[j];
                    acc         += alpha * aij * x1[jj + j];
                    y2p[jj + j] += aij * xi2;
                }
                y1[i] += acc;
            }
        }
    }
}

// ── 6b: fused blocked A (outer×inner = {R,C}²) ───────────────────────────
// collapse(2) distributes NB² block pairs across threads.
// ie and je are derived inside the body.
template<bool outer_col, bool inner_col>
static void kern_fused_blk(const double *Ab,
                            const double *x1, const double *x2,
                            double *y1, double **y2_priv)
{
    #pragma omp parallel for collapse(2) schedule(static)
    for (int bi = 0; bi < NB; bi++) {
        for (int bj = 0; bj < NB; bj++) {
            const int ie = std::min(bi * SZ_B + SZ_B, N);
            const int je = std::min(bj * SZ_B + SZ_B, N);
            double *y2p = y2_priv[omp_get_thread_num()];
            for (int i = bi * SZ_B; i < ie; i++) {
                double acc = 0.0;
                const double xi2 = alpha * x2[i];
                for (int j = bj * SZ_B; j < je; j++) {
                    const double aij = Ab[blkidx<outer_col, inner_col>(i, j)];
                    acc     += alpha * aij * x1[j];
                    y2p[j] += aij * xi2;
                }
                y1[i] += acc;
            }
        }
    }
}

// ── runner for mode 6 ────────────────────────────────────────────────────
struct FusedVariant {
    const char *tag;
    bool  uses_blk;
    void (*to_blk_fn)(const double*, double*);
    void (*kern_plain)(const double*, const double*, const double*, double*, double**);
    void (*kern_blk  )(const double*, const double*, const double*, double*, double**);
};

static const FusedVariant FUSED_VARIANTS[] = {
    { "fused_tiled",  false, nullptr,    kern_fused_tiled,            nullptr                    },
    { "fused_blk_RR", true,  to_blk_rr, nullptr, kern_fused_blk<false,false> },
    { "fused_blk_RC", true,  to_blk_rc, nullptr, kern_fused_blk<false,true>  },
    { "fused_blk_CR", true,  to_blk_cr, nullptr, kern_fused_blk<true, false> },
    { "fused_blk_CC", true,  to_blk_cc, nullptr, kern_fused_blk<true, true>  },
};

static void run_fused(const FusedVariant &v,
                      const double *A_rm,
                      const double *x1, const double *x2,
                      const double *y10, const double *y20,
                      const double *y1ref, const double *y2ref,
                      FILE *fp)
{
    // Pack A once into required layout (outside timing loop).
    double *Ab = nullptr;
    const double *Aptr = A_rm;
    if (v.uses_blk) {
        Ab = new double[BSIZ];
        v.to_blk_fn(A_rm, Ab);
        Aptr = Ab;
    }

    // Allocate per-thread y2 buffers once for all NRUNS; re-zero each run.
    double **y2_priv; const int nt = y2_priv_alloc(y2_priv);

    double *y1w = new double[N];
    double *y2w = new double[N];

    for (int r = 0; r < NRUNS; r++) {
        memcpy(y1w, y10, N * sizeof(double));
        memcpy(y2w, y20, N * sizeof(double));
        y2_priv_zero(y2_priv, nt);
        auto t0 = Clock::now();
        if (v.uses_blk)
            v.kern_blk  (Aptr, x1, x2, y1w, y2_priv);
        else
            v.kern_plain(Aptr, x1, x2, y1w, y2_priv);
        y2_priv_reduce(y2_priv, nt, y2w);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,6,%s,%d,%d,%d,%.2f\n",
                    N, v.tag,
                    v.uses_blk ? SZ_B : 0,
                    v.uses_blk ? 0 : SZ_T,
                    r - NWARM,
                    std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    double e1 = max_relerr(y1w, y1ref, N);
    double e2 = max_relerr(y2w, y2ref, N);
    fprintf(stderr, "%s %s: y1_err=%.2e y2_err=%.2e\n",
            (e1 > 1e-9 || e2 > 1e-9) ? "FAIL" : "OK", v.tag, e1, e2);

    y2_priv_free(y2_priv, nt);
    delete[] y1w; delete[] y2w;
    delete[] Ab;
}

// ── main ────────────────────────────────────────────────────────────────
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <mode:1-6> [output.csv]\n", argv[0]);
        fprintf(stderr, "  1 = A layout (R,C)\n");
        fprintf(stderr, "  2 = loop orderings\n");
        fprintf(stderr, "  3 = blocked A, SZ_B=%d, variants RR/RC/CR/CC\n", SZ_B);
        fprintf(stderr, "  4 = tiled loops, SZ_T=%d\n", SZ_T);
        fprintf(stderr, "  5 = OpenBLAS 2x dgemv\n");
        fprintf(stderr, "  6 = fused single-pass: tiled (SZ_T=%d) + blk_{RR,RC,CR,CC} (SZ_B=%d)\n", SZ_T, SZ_B);
        return 1;
    }
    int mode = atoi(argv[1]);
    const char *csv = (argc > 2) ? argv[2] : "results_mvt.csv";

    // ── allocate and initialise ──
    double *A   = new double[(long)N * N];
    double *x1  = new double[N];
    double *x2  = new double[N];
    double *y10 = new double[N];   // initial y vectors (kept read-only)
    double *y20 = new double[N];
    double *y1ref = new double[N];
    double *y2ref = new double[N];

    init_rand(A,   (long)N * N);
    init_rand(x1,  N);
    init_rand(x2,  N);
    init_rand(y10, N);
    init_rand(y20, N);

    // ── reference: serial ij_ji, row-major A ──
    fprintf(stderr, "computing reference (N=%d)...\n", N);
    memcpy(y1ref, y10, N * sizeof(double));
    memcpy(y2ref, y20, N * sizeof(double));
    {
        int saved = omp_get_max_threads();
        omp_set_num_threads(1);
        kern_ij_ji(A, x1, x2, y1ref, y2ref);
        omp_set_num_threads(saved);
    }

    // ── open CSV (append; write header once) ──
    FILE *fp = fopen(csv, "a");
    if (!fp) { perror("fopen"); return 1; }
    fseek(fp, 0, SEEK_END);
    if (ftell(fp) == 0)
        fprintf(fp, "n,mode,variant,block_size,tile_size,run,time_us\n");

    switch (mode) {
    case 1:
        fprintf(stderr, "mode 1: A layout variants\n");
        run_layout(R, A, x1, x2, y10, y20, y1ref, y2ref, fp, "A_row");
        run_layout(C, A, x1, x2, y10, y20, y1ref, y2ref, fp, "A_col");
        break;
    case 2:
        fprintf(stderr, "mode 2: loop orderings\n");
        run_loop(kern_ij_ji, A, x1, x2, y10, y20, y1ref, y2ref, fp, "ij_ji");
        run_loop(kern_ji_ji, A, x1, x2, y10, y20, y1ref, y2ref, fp, "ji_ji");
        run_loop(kern_ij_ij, A, x1, x2, y10, y20, y1ref, y2ref, fp, "ij_ij");
        run_loop(kern_ji_ij, A, x1, x2, y10, y20, y1ref, y2ref, fp, "ji_ij");
        break;
    case 3:
        fprintf(stderr, "mode 3: blocked A SZ_B=%d\n", SZ_B);
        for (auto &v : BLK_VARIANTS)
            run_blocked(v, A, x1, x2, y10, y20, y1ref, y2ref, fp);
        break;
    case 4:
        fprintf(stderr, "mode 4: tiled loops SZ_T=%d\n", SZ_T);
        run_tiled(A, x1, x2, y10, y20, y1ref, y2ref, fp);
        break;
    case 5:
        fprintf(stderr, "mode 5: OpenBLAS 2x dgemv\n");
        run_openblas(A, x1, x2, y10, y20, y1ref, y2ref, fp);
        break;
    case 6:
        fprintf(stderr, "mode 6: fused single-pass variants\n");
        for (auto &v : FUSED_VARIANTS)
            run_fused(v, A, x1, x2, y10, y20, y1ref, y2ref, fp);
        break;
    default:
        fprintf(stderr, "unknown mode %d\n", mode);
        fclose(fp); return 1;
    }

    fclose(fp);
    delete[] A; delete[] x1; delete[] x2;
    delete[] y10; delete[] y20; delete[] y1ref; delete[] y2ref;
    fprintf(stderr, "done. results appended to %s\n", csv);
    return 0;
}