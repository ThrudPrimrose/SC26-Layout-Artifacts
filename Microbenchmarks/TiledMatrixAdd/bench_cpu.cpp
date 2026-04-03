/*
 * Layout-conflict benchmark (CPU) — NUMA-aware, per-iteration CSV
 * ================================================================
 * C[i,j] = A[i,j] + B[i,j]
 *
 * Row-major schedules:
 *   row_major     — for i, for j  (A,C stream; B stride-M)
 *   col_major     — for j, for i  (B streams; A,C stride-N)
 *   tiled_T       — T×T tiles, stack-local B transpose
 *   all_rowmajor  — control (B also row-major)
 *
 * Blocked-storage schedules (SB×SB contiguous blocks):
 *   blk_all_rm    — A,B,C all row-major-in-block  (control)
 *   blk_conflict  — A,C row-major-in-block; B col-major-in-block
 *
 * NUMA: mmap + MADV_NOHUGEPAGE + per-thread mbind(MPOL_BIND) + first-touch
 * CSV:  one row per iteration for violin plots
 *
 * Compile: g++ -O3 -march=native -std=c++17 -fopenmp -o bench_cpu bench_cpu.cpp
 * Run:     ./bench_cpu <csv_file> [nthreads]
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>

#include <omp.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>

#if __has_include(<numaif.h>)
  #include <numaif.h>
  #define HAS_NUMA 1
#else
  #define HAS_NUMA 0
  #define MPOL_BIND 2
  static int mbind(void*, size_t, int, const unsigned long*, unsigned long, unsigned) { return 0; }
#endif

#ifndef M_DIM
#define M_DIM 16384
#endif
#ifndef N_DIM
#define N_DIM 16384
#endif

static const int M = M_DIM;
static const int N = N_DIM;

#define NREP    100
#define NWARMUP 5

/* ── index helpers ────────────────────────────────────────────────── */
static inline int idx_rm(int i, int j) { return i * N + j; }
static inline int idx_cm(int i, int j) { return j * M + i; }

/* Blocked storage: block (br,bc) is contiguous SB*SB.
 * Within block: row-major = lr*SB+lc, col-major = lc*SB+lr. */
template <int SB>
static inline int idx_blk_rm(int i, int j) {
    int NB = N / SB;
    return (i / SB * NB + j / SB) * SB * SB + (i % SB) * SB + (j % SB);
}
template <int SB>
static inline int idx_blk_cm(int i, int j) {
    int NB = N / SB;
    return (i / SB * NB + j / SB) * SB * SB + (j % SB) * SB + (i % SB);
}

/* ── timing ───────────────────────────────────────────────────────── */
using Clock = std::chrono::high_resolution_clock;
static inline double elapsed_sec(Clock::time_point t0, Clock::time_point t1) {
    return std::chrono::duration<double>(t1 - t0).count();
}

/* ════════════════════════════════════════════════════════════════════
 *  NUMA-aware allocation
 * ════════════════════════════════════════════════════════════════════ */
static long PAGE_SZ;

static int get_numa_node() {
    unsigned cpu, node;
    syscall(__NR_getcpu, &cpu, &node, nullptr);
    return (int)node;
}

static void *numa_alloc(size_t bytes) {
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return p;
}

static void numa_free(void *p, size_t bytes) { munmap(p, bytes); }

static void bind_and_touch(void *base, size_t total_bytes) {
    int64_t n_pages = (total_bytes + PAGE_SZ - 1) / PAGE_SZ;
    #pragma omp parallel
    {
        int tid  = omp_get_thread_num();
        int nthr = omp_get_num_threads();
        int64_t chunk = (n_pages + nthr - 1) / nthr;
        int64_t lo = tid * chunk;
        int64_t hi = std::min(lo + chunk, n_pages);
        if (lo < n_pages && lo < hi) {
            uintptr_t pbeg = (uintptr_t)base + lo * PAGE_SZ;
            uintptr_t pend = (uintptr_t)base + hi * PAGE_SZ;
            if (pend > (uintptr_t)base + total_bytes)
                pend = (uintptr_t)base + total_bytes;

            int node = get_numa_node();
            unsigned long mask[16] = {};
            mask[node / (8 * sizeof(unsigned long))] =
                1UL << (node % (8 * sizeof(unsigned long)));
            mbind((void *)pbeg, pend - pbeg, MPOL_BIND, mask,
                  sizeof(mask) * 8 + 1, 0);

            for (uintptr_t addr = pbeg; addr < pend; addr += PAGE_SZ)
                *(volatile char *)addr = 0;
        }
    }
}

/* First-touch init for row-major */
static void ft_init_rm(double *buf, size_t total, int M_, int N_, bool is_B) {
    bind_and_touch(buf, total * sizeof(double));
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++) {
        int i = (int)(k / N_), j = (int)(k % N_);
        buf[k] = is_B ? (double)((i * 13 + j * 37) % 1000) / 100.0
                       : (double)((i * 17 + j * 31) % 1000) / 100.0;
    }
}

/* First-touch init for col-major B[j*M+i] */
static void ft_init_cm(double *buf, size_t total, int M_, int N_) {
    bind_and_touch(buf, total * sizeof(double));
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N_; j++)
        for (int i = 0; i < M_; i++)
            buf[j * M_ + i] = (double)((i * 13 + j * 37) % 1000) / 100.0;
}

/* First-touch init for blocked-rowmajor */
template <int SB>
static void ft_init_blk_rm(double *buf, size_t total, int M_, int N_, bool is_B) {
    bind_and_touch(buf, total * sizeof(double));
    int NB = N_ / SB;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int br = 0; br < M_ / SB; br++)
        for (int bc = 0; bc < NB; bc++) {
            double *blk = buf + (br * NB + bc) * SB * SB;
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++) {
                    int i = br * SB + lr, j = bc * SB + lc;
                    blk[lr * SB + lc] = is_B
                        ? (double)((i * 13 + j * 37) % 1000) / 100.0
                        : (double)((i * 17 + j * 31) % 1000) / 100.0;
                }
        }
}

/* First-touch init for blocked-colmajor (B with col-major inner layout) */
template <int SB>
static void ft_init_blk_cm(double *buf, size_t total, int M_, int N_) {
    bind_and_touch(buf, total * sizeof(double));
    int NB = N_ / SB;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int br = 0; br < M_ / SB; br++)
        for (int bc = 0; bc < NB; bc++) {
            double *blk = buf + (br * NB + bc) * SB * SB;
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++) {
                    int i = br * SB + lr, j = bc * SB + lc;
                    /* col-major within block: index = lc*SB + lr */
                    blk[lc * SB + lr] = (double)((i * 13 + j * 37) % 1000) / 100.0;
                }
        }
}

/* ════════════════════════════════════════════════════════════════════
 *  Cache flush
 * ════════════════════════════════════════════════════════════════════ */
static constexpr int64_t FLUSH_N = 1 << 27;
static double *g_flush_buf = nullptr;

static void cache_flush() {
    if (!g_flush_buf) {
        g_flush_buf = (double *)numa_alloc(FLUSH_N * sizeof(double));
        bind_and_touch(g_flush_buf, FLUSH_N * sizeof(double));
    }
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++)
        g_flush_buf[i] = (double)i;
    volatile double sink = g_flush_buf[FLUSH_N - 1];
    (void)sink;
}

/* ════════════════════════════════════════════════════════════════════
 *  Row-major kernels (unchanged)
 * ════════════════════════════════════════════════════════════════════ */

static void kernel_row_major(const double *__restrict__ A,
                             const double *__restrict__ B,
                             double *__restrict__ C) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < M; i++)
        #pragma omp simd 
        for (int j = 0; j < N; j++)
            C[idx_rm(i, j)] += A[idx_rm(i, j)] + B[idx_cm(i, j)];
}

static void kernel_col_major(const double *__restrict__ A,
                             const double *__restrict__ B,
                             double *__restrict__ C) {
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++)
        #pragma omp simd 
        for (int i = 0; i < M; i++)
            C[idx_rm(i, j)] += A[idx_rm(i, j)] + B[idx_cm(i, j)];
}

template <int T>
static void kernel_tiled(const double *__restrict__ A,
                         const double *__restrict__ B,
                         double *__restrict__ C) {
    const int ntiles_i = (M + T - 1) / T;
    const int ntiles_j = (N + T - 1) / T;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int ti = 0; ti < ntiles_i; ti++) {
        for (int tj = 0; tj < ntiles_j; tj++) {
            const int ii = ti * T, jj = tj * T;
            const int iend = (ii + T < M) ? ii + T : M;
            const int jend = (jj + T < N) ? jj + T : N;

            double b_local[T * T];
            for (int j = jj; j < jend; j++)
                #pragma omp simd
                for (int i = ii; i < iend; i++)
                    b_local[(i - ii) * T + (j - jj)] = B[idx_cm(i, j)];

            for (int i = ii; i < iend; i++)
                #pragma omp simd 
                for (int j = jj; j < jend; j++)
                    C[idx_rm(i, j)] += A[idx_rm(i, j)]
                                    + b_local[(i - ii) * T + (j - jj)];
        }
    }
}

static void kernel_all_rowmajor(const double *__restrict__ A,
                                const double *__restrict__ B_rm,
                                double *__restrict__ C) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < M; i++)
        #pragma omp simd 
        for (int j = 0; j < N; j++)
            C[idx_rm(i, j)] += A[idx_rm(i, j)] + B_rm[idx_rm(i, j)];
}

/* ════════════════════════════════════════════════════════════════════
 *  Blocked-storage kernels
 *
 *  Block (br,bc) for each array is contiguous SB*SB doubles.
 *  Within each block, arrays use either row-major or col-major layout.
 *  The block fits in L1 cache for SB ≤ ~64, so even conflicting
 *  inner layouts (A row-major, B col-major) incur no DRAM-level
 *  cache line waste — the "conflict cost" is absorbed by L1.
 * ════════════════════════════════════════════════════════════════════ */

/* blk_all_rm: A,B,C all row-major-in-block.  Control — no conflict. */
template <int SB>
static void kernel_blk_all_rm(const double *__restrict__ A,
                              const double *__restrict__ B,
                              double *__restrict__ C) {
    const int NB_i = M / SB;
    const int NB_j = N / SB;
    const int NB = NB_j;  /* blocks per row */
    #pragma omp parallel for schedule(static) collapse(2)
    for (int br = 0; br < NB_i; br++) {
        for (int bc = 0; bc < NB_j; bc++) {
            const double *a_blk = A + (br * NB + bc) * SB * SB;
            const double *b_blk = B + (br * NB + bc) * SB * SB;
            double       *c_blk = C + (br * NB + bc) * SB * SB;
            #pragma unroll
            for (int lr = 0; lr < SB; lr++)
                #pragma omp simd
                for (int lc = 0; lc < SB; lc++)
                    c_blk[lr * SB + lc] += a_blk[lr * SB + lc]
                                        + b_blk[lr * SB + lc];
        }
    }
}

/* blk_conflict: A,C row-major-in-block; B col-major-in-block.
 * B is stored as blk[lc*SB + lr] instead of blk[lr*SB + lc]. */
template <int SB>
static void kernel_blk_conflict(const double *__restrict__ A,
                                const double *__restrict__ B,
                                double *__restrict__ C) {
    const int NB_i = M / SB;
    const int NB_j = N / SB;
    const int NB = NB_j;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int br = 0; br < NB_i; br++) {
        for (int bc = 0; bc < NB_j; bc++) {
            const double *a_blk = A + (br * NB + bc) * SB * SB;
            const double *b_blk = B + (br * NB + bc) * SB * SB;
            double       *c_blk = C + (br * NB + bc) * SB * SB;
            #pragma unroll
            for (int lr = 0; lr < SB; lr++)
                #pragma omp simd
                for (int lc = 0; lc < SB; lc++)
                    c_blk[lr * SB + lc] += a_blk[lr * SB + lc]
                                        + b_blk[lc * SB + lr];  /* col-major B */
        }
    }
}

/* ════════════════════════════════════════════════════════════════════
 *  Benchmark infrastructure
 * ════════════════════════════════════════════════════════════════════ */
using KernelFn = void (*)(const double *, const double *, double *);

struct IterRecord {
    std::string variant;
    int tile;
    int nthreads;
    int rep;
    double time_s;
    double bw_gbs;
    double checksum;
    const char *status;
};

static std::vector<IterRecord> g_records;

static void bench(const char *variant_name, int tile_sz, KernelFn fn,
                  const double *A, const double *B, double *C,
                  double ref_checksum, int nthreads) {
    const size_t total = (size_t)M * N;
    /* C += A + B:  read A + read B + read C + write C = 4 array transfers */
    const double data_bytes = (double)total * sizeof(double) * 4.0;

    /* Zero C before warmup */
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++)
        C[k] = 0.0;

    for (int w = 0; w < NWARMUP; w++) {
        cache_flush();
        fn(A, B, C);
    }

    /* Zero C and run once for correctness check */
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++)
        C[k] = 0.0;
    fn(A, B, C);

    double cs = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:cs)
    for (size_t k = 0; k < total; k++)
        cs += C[k];
    bool ok = (std::fabs(cs - ref_checksum) <= 1e-3 * std::fabs(ref_checksum));
    const char *status = ok ? "PASS" : "FAIL";

    if (!ok)
        fprintf(stderr, "  [%s T=%d] CHECKSUM MISMATCH: got %.6e expected %.6e\n",
                variant_name, tile_sz, cs, ref_checksum);

    for (int r = 0; r < NREP; r++) {
        cache_flush();
        auto t0 = Clock::now();
        fn(A, B, C);
        auto t1 = Clock::now();
        double dt = elapsed_sec(t0, t1);
        double bw = data_bytes / dt / 1e9;
        g_records.push_back({variant_name, tile_sz, nthreads, r, dt, bw, cs, status});
    }

    std::vector<double> bws;
    for (int i = (int)g_records.size() - NREP; i < (int)g_records.size(); i++)
        bws.push_back(g_records[i].bw_gbs);
    std::sort(bws.begin(), bws.end());
    printf("  %-28s T=%-4d  median %7.1f GB/s  [%7.1f .. %7.1f]  %s\n",
           variant_name, tile_sz, bws[NREP / 2], bws[0], bws[NREP - 1], status);
}

/* ════════════════════════════════════════════════════════════════════
 *  Blocked-storage benchmark wrapper
 *
 *  Allocates separate blocked arrays, inits with NUMA first-touch,
 *  computes reference checksum, runs the benchmark, frees arrays.
 * ════════════════════════════════════════════════════════════════════ */

template <int SB>
static void bench_blocked(int nthreads) {
    static_assert(M_DIM % SB == 0 && N_DIM % SB == 0,
                  "M and N must be divisible by SB");

    const size_t total = (size_t)M * N;
    const size_t bytes = total * sizeof(double);

    printf("\n  --- Blocked SB=%d ---\n", SB);

    double *A_blk = (double *)numa_alloc(bytes);
    double *B_blk_rm = (double *)numa_alloc(bytes);  /* B row-major-in-block */
    double *B_blk_cm = (double *)numa_alloc(bytes);  /* B col-major-in-block */
    double *C_blk = (double *)numa_alloc(bytes);

    ft_init_blk_rm<SB>(A_blk, total, M, N, false);
    ft_init_blk_rm<SB>(B_blk_rm, total, M, N, true);
    ft_init_blk_cm<SB>(B_blk_cm, total, M, N);
    bind_and_touch(C_blk, bytes);
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++)
        C_blk[k] = 0.0;

    /* Reference checksum (blk_all_rm) */
    kernel_blk_all_rm<SB>(A_blk, B_blk_rm, C_blk);
    double ref_cs = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:ref_cs)
    for (size_t k = 0; k < total; k++)
        ref_cs += C_blk[k];

    bench("blk_all_rm", SB, kernel_blk_all_rm<SB>, A_blk, B_blk_rm, C_blk,
          ref_cs, nthreads);
    bench("blk_conflict", SB, kernel_blk_conflict<SB>, A_blk, B_blk_cm, C_blk,
          ref_cs, nthreads);

    numa_free(A_blk, bytes);
    numa_free(B_blk_rm, bytes);
    numa_free(B_blk_cm, bytes);
    numa_free(C_blk, bytes);
}

/* ════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr,
                "Usage: %s <csv_file> [nthreads]\n"
                "  csv_file:  output CSV (one row per iteration)\n"
                "  nthreads:  override OMP_NUM_THREADS (0 = use env)\n",
                argv[0]);
        return 1;
    }

    const char *csv_path = argv[1];
    int req_threads = (argc > 2) ? atoi(argv[2]) : 0;
    if (req_threads > 0)
        omp_set_num_threads(req_threads);

    PAGE_SZ = sysconf(_SC_PAGESIZE);
    int nthreads;
    #pragma omp parallel
    { nthreads = omp_get_num_threads(); }

    const size_t total = (size_t)M * N;
    const size_t bytes = total * sizeof(double);

    printf("Layout-conflict benchmark (CPU)\n");
    printf("  M=%d  N=%d  reps=%d  warmup=%d\n", M, N, NREP, NWARMUP);
    printf("  threads=%d  page=%ld B  NUMA=%s\n",
           nthreads, PAGE_SZ, HAS_NUMA ? "yes" : "fallback");
    printf("  A: row-major   B: col-major   C: row-major\n");
    printf("  cache line = 64 B  ->  %.0f doubles/line\n\n",
           64.0 / sizeof(double));

    /* ── Row-major arrays ──────────────────────────────────────────── */
    double *A = (double *)numa_alloc(bytes);
    double *B_cm = (double *)numa_alloc(bytes);
    double *B_rm = (double *)numa_alloc(bytes);
    double *C = (double *)numa_alloc(bytes);

    printf("Initializing row-major arrays (NUMA bind+touch) ...\n");
    ft_init_rm(A, total, M, N, false);
    ft_init_cm(B_cm, total, M, N);
    ft_init_rm(B_rm, total, M, N, true);
    bind_and_touch(C, bytes);
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++)
        C[k] = 0.0;

    /* Reference checksums — zero C before each to get single-pass result */
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++)
        C[k] = 0.0;
    kernel_row_major(A, B_cm, C);
    double ref_cs = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:ref_cs)
    for (size_t k = 0; k < total; k++)
        ref_cs += C[k];

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++)
        C[k] = 0.0;
    kernel_all_rowmajor(A, B_rm, C);
    double ref_cs_rm = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:ref_cs_rm)
    for (size_t k = 0; k < total; k++)
        ref_cs_rm += C[k];

    /* ── Row-major benchmarks ──────────────────────────────────────── */
    printf("\n=== Naive schedules ===\n");
    bench("row_major", 0, kernel_row_major, A, B_cm, C, ref_cs, nthreads);
    bench("col_major", 0, kernel_col_major, A, B_cm, C, ref_cs, nthreads);
    bench("all_rowmajor", 0, kernel_all_rowmajor, A, B_rm, C, ref_cs_rm, nthreads);

    printf("\n=== Tiled schedules (stack-local B transpose) ===\n");
    bench("tiled", 8, kernel_tiled<8>, A, B_cm, C, ref_cs, nthreads);
    bench("tiled", 16, kernel_tiled<16>, A, B_cm, C, ref_cs, nthreads);
    bench("tiled", 32, kernel_tiled<32>, A, B_cm, C, ref_cs, nthreads);
    bench("tiled", 64, kernel_tiled<64>, A, B_cm, C, ref_cs, nthreads);
    bench("tiled", 128, kernel_tiled<128>, A, B_cm, C, ref_cs, nthreads);

    numa_free(A, bytes);
    numa_free(B_cm, bytes);
    numa_free(B_rm, bytes);
    numa_free(C, bytes);

    /* ── Blocked-storage benchmarks ────────────────────────────────── */
    printf("\n=== Blocked storage (SB×SB contiguous blocks) ===\n");
    bench_blocked<8>(nthreads);
    bench_blocked<16>(nthreads);
    bench_blocked<32>(nthreads);
    bench_blocked<64>(nthreads);
    bench_blocked<128>(nthreads);

    /* ── Write CSV ─────────────────────────────────────────────────── */
    FILE *fp = fopen(csv_path, "w");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", csv_path);
        return 1;
    }
    fprintf(fp, "variant,M,N,tile,nthreads,rep,time_s,bw_gbs,checksum,status\n");
    for (auto &r : g_records)
        fprintf(fp, "%s,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e,%s\n",
                r.variant.c_str(), M, N, r.tile, r.nthreads,
                r.rep, r.time_s, r.bw_gbs, r.checksum, r.status);
    fclose(fp);
    printf("\nWrote %zu records to %s\n", g_records.size(), csv_path);

    /* ── Summary ───────────────────────────────────────────────────── */
    auto median_bw = [&](const char *name, int tile) -> double {
        std::vector<double> bws;
        for (auto &r : g_records)
            if (r.variant == name && r.tile == tile)
                bws.push_back(r.bw_gbs);
        if (bws.empty())
            return 0;
        std::sort(bws.begin(), bws.end());
        return bws[bws.size() / 2];
    };

    double best_tiled = 0;
    for (int t : {8, 16, 32, 64, 128})
        best_tiled = std::max(best_tiled, median_bw("tiled", t));
    double bw_row = median_bw("row_major", 0);
    double bw_col = median_bw("col_major", 0);
    double bw_ctrl = median_bw("all_rowmajor", 0);

    double best_blk_all = 0, best_blk_conflict = 0;
    int best_blk_all_sb = 0, best_blk_conflict_sb = 0;
    for (int sb : {8, 16, 32, 64, 128}) {
        double bw;
        bw = median_bw("blk_all_rm", sb);
        if (bw > best_blk_all) { best_blk_all = bw; best_blk_all_sb = sb; }
        bw = median_bw("blk_conflict", sb);
        if (bw > best_blk_conflict) { best_blk_conflict = bw; best_blk_conflict_sb = sb; }
    }

    double B_eff = 64.0 / sizeof(double);
    printf("\n=== Summary ===\n");
    printf("  Row-major schedules:\n");
    printf("    best tiled:   %7.1f GB/s\n", best_tiled);
    printf("    row_major:    %7.1f GB/s  (%.1fx vs tiled)\n",
           bw_row, best_tiled / std::max(bw_row, 0.1));
    printf("    col_major:    %7.1f GB/s  (%.1fx vs tiled)\n",
           bw_col, best_tiled / std::max(bw_col, 0.1));
    printf("    all_rowmajor: %7.1f GB/s  (peak control)\n", bw_ctrl);
    printf("  Blocked storage:\n");
    printf("    blk_all_rm:   %7.1f GB/s  SB=%d  (%.1f%% of all_rowmajor)\n",
           best_blk_all, best_blk_all_sb,
           100.0 * best_blk_all / std::max(bw_ctrl, 0.1));
    printf("    blk_conflict: %7.1f GB/s  SB=%d  (%.1f%% of blk_all_rm)\n",
           best_blk_conflict, best_blk_conflict_sb,
           100.0 * best_blk_conflict / std::max(best_blk_all, 0.1));
    printf("  Conflict cost absorbed by blocking:\n");
    printf("    row_major conflict:  %.1fx slowdown vs all_rowmajor\n",
           bw_ctrl / std::max(bw_row, 0.1));
    printf("    blocked conflict:    %.1fx slowdown vs blk_all_rm\n",
           best_blk_all / std::max(best_blk_conflict, 0.1));

    printf("\n=== Model (B_eff=%.0f) ===\n", B_eff);
    printf("  row_major predicted:  %.1fx  (cost %g vs 3)\n",
           (2 + B_eff) / 3.0, 2 + B_eff);
    printf("  col_major predicted:  %.1fx  (cost %g vs 3)\n",
           (1 + 2 * B_eff) / 3.0, 1 + 2 * B_eff);

    if (g_flush_buf)
        numa_free(g_flush_buf, FLUSH_N * sizeof(double));

    return 0;
}