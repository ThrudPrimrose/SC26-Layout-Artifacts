/*
 * Layout-conflict benchmark (CPU) — NUMA-aware, per-iteration CSV
 * ================================================================
 * C[i,j] = A[i,j] + B[i,j]
 *
 * Layouts:  A row-major, B col-major, C row-major
 *
 * Schedules:
 *   row_major   — for i, for j  (A,C stream; B stride-M)
 *   col_major   — for j, for i  (B streams; A,C stride-N)
 *   tiled_T     — T×T tiles, stack-local B transpose (swept T)
 *   all_rowmajor — control (B also row-major)
 *
 * NUMA: mmap + MADV_NOHUGEPAGE + per-thread mbind(MPOL_BIND) + first-touch
 * CSV:  one row per iteration for violin plots / distribution analysis
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

/* NUMA headers — graceful fallback if unavailable */
#if __has_include(<numaif.h>)
  #include <numaif.h>
  #define HAS_NUMA 1
#else
  #define HAS_NUMA 0
  #define MPOL_BIND 2
  static int mbind(void*, size_t, int, const unsigned long*, unsigned long, unsigned) { return 0; }
#endif

/* ---- dimensions ---- */
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

/* ---- index helpers ---- */
static inline int idx_rm(int i, int j) { return i * N + j; }
static inline int idx_cm(int i, int j) { return j * M + i; }

/* ---- timing ---- */
using Clock = std::chrono::high_resolution_clock;
static inline double elapsed_sec(Clock::time_point t0, Clock::time_point t1) {
    return std::chrono::duration<double>(t1 - t0).count();
}


/* ================================================================
 *  NUMA-aware allocation
 *
 *  1. mmap with MAP_ANONYMOUS | MAP_NORESERVE (no physical pages yet)
 *  2. MADV_NOHUGEPAGE to prevent THP from defeating first-touch
 *  3. bind_and_touch: each thread mbinds its schedule(static) chunk
 *     to its local NUMA node, then faults each page
 *  4. Data init with same schedule(static) → pages on local node
 * ================================================================ */

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
    madvise(p, bytes, MADV_NOHUGEPAGE);
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

/* First-touch init for row-major array */
static void ft_init_rm(double *buf, size_t total, int M_, int N_,
                       bool is_B_values) {
    bind_and_touch(buf, total * sizeof(double));
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++) {
        int i = (int)(k / N_), j = (int)(k % N_);
        if (is_B_values)
            buf[k] = (double)((i * 13 + j * 37) % 1000) / 100.0;
        else
            buf[k] = (double)((i * 17 + j * 31) % 1000) / 100.0;
    }
}

/* First-touch init for col-major array B[j*M+i] */
static void ft_init_cm(double *buf, size_t total, int M_, int N_) {
    bind_and_touch(buf, total * sizeof(double));
    /* Touch in column-major order so pages are distributed
       proportionally to the access pattern */
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N_; j++)
        for (int i = 0; i < M_; i++)
            buf[j * M_ + i] = (double)((i * 13 + j * 37) % 1000) / 100.0;
}


/* ================================================================
 *  Cache flush (touch a large buffer to evict working set)
 * ================================================================ */
static constexpr int64_t FLUSH_N = 1 << 24;  /* 128 MiB */
static double *g_flush_buf = nullptr;

static void cache_flush() {
    if (!g_flush_buf) {
        g_flush_buf = (double *)numa_alloc(FLUSH_N * sizeof(double));
        bind_and_touch(g_flush_buf, FLUSH_N * sizeof(double));
    }
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++)
        g_flush_buf[i] = (double)i;
    /* prevent dead-code elimination */
    volatile double sink = g_flush_buf[FLUSH_N - 1];
    (void)sink;
}


/* ================================================================
 *  Kernels
 * ================================================================ */

static void kernel_row_major(const double* __restrict__ A,
                             const double* __restrict__ B,
                             double*       __restrict__ C) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < M; i++)
#pragma omp simd nontemporal(C)
        for (int j = 0; j < N; j++)
            C[idx_rm(i, j)] = A[idx_rm(i, j)] + B[idx_cm(i, j)];
}

static void kernel_col_major(const double* __restrict__ A,
                             const double* __restrict__ B,
                             double*       __restrict__ C) {
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++)
#pragma omp simd nontemporal(C)
        for (int i = 0; i < M; i++)
            C[idx_rm(i, j)] = A[idx_rm(i, j)] + B[idx_cm(i, j)];
}

template<int T>
static void kernel_tiled(const double* __restrict__ A,
                         const double* __restrict__ B,
                         double*       __restrict__ C) {
    const int ntiles_i = (M + T - 1) / T;
    const int ntiles_j = (N + T - 1) / T;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int ti = 0; ti < ntiles_i; ti++) {
        for (int tj = 0; tj < ntiles_j; tj++) {
            const int ii = ti * T;
            const int jj = tj * T;
            const int iend = (ii + T < M) ? ii + T : M;
            const int jend = (jj + T < N) ? jj + T : N;

            /* Phase 1: load B tile col-major → row-major local buffer */
            double b_local[T * T];
            for (int j = jj; j < jend; j++)
#pragma omp simd
                for (int i = ii; i < iend; i++)
                    b_local[(i - ii) * T + (j - jj)] = B[idx_cm(i, j)];

            /* Phase 2: compute with unit stride on all arrays */
            for (int i = ii; i < iend; i++)
#pragma omp simd nontemporal(C)
                for (int j = jj; j < jend; j++)
                    C[idx_rm(i, j)] = A[idx_rm(i, j)]
                                    + b_local[(i - ii) * T + (j - jj)];
        }
    }
}

static void kernel_all_rowmajor(const double* __restrict__ A,
                                const double* __restrict__ B_rm,
                                double*       __restrict__ C) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < M; i++)
#pragma omp simd nontemporal(C)
        for (int j = 0; j < N; j++)
            C[idx_rm(i, j)] = A[idx_rm(i, j)] + B_rm[idx_rm(i, j)];
}


/* ================================================================
 *  Benchmark infrastructure — per-iteration CSV dump
 * ================================================================ */
using KernelFn = void(*)(const double*, const double*, double*);

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

static void bench(const char* variant_name, int tile_sz, KernelFn fn,
                  const double* A, const double* B, double* C,
                  double ref_checksum, int nthreads)
{
    const size_t total = (size_t)M * N;
    const double data_bytes = (double)total * sizeof(double) * 3.0;

    /* warmup */
    for (int w = 0; w < NWARMUP; w++) {
        cache_flush();
        fn(A, B, C);
    }

    /* correctness check */
    double cs = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:cs)
    for (size_t k = 0; k < total; k++) cs += C[k];
    bool ok = (std::fabs(cs - ref_checksum) <= 1e-3 * std::fabs(ref_checksum));
    const char *status = ok ? "PASS" : "FAIL";

    if (!ok)
        fprintf(stderr, "  [%s T=%d] CHECKSUM MISMATCH: got %.6e expected %.6e\n",
                variant_name, tile_sz, cs, ref_checksum);

    /* timed runs — record every iteration */
    for (int r = 0; r < NREP; r++) {
        cache_flush();
        auto t0 = Clock::now();
        fn(A, B, C);
        auto t1 = Clock::now();
        double dt = elapsed_sec(t0, t1);
        double bw = data_bytes / dt / 1e9;
        g_records.push_back({variant_name, tile_sz, nthreads,
                             r, dt, bw, cs, status});
    }

    /* console summary */
    std::vector<double> bws;
    for (int i = (int)g_records.size() - NREP; i < (int)g_records.size(); i++)
        bws.push_back(g_records[i].bw_gbs);
    std::sort(bws.begin(), bws.end());
    printf("  %-28s T=%-4d  median %7.1f GB/s  [%7.1f .. %7.1f]  %s\n",
           variant_name, tile_sz, bws[NREP/2], bws[0], bws[NREP-1], status);
}


/* ================================================================ */
int main(int argc, char** argv)
{
    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s <csv_file> [nthreads]\n"
            "  csv_file:  output CSV (one row per iteration)\n"
            "  nthreads:  override OMP_NUM_THREADS (0 = use env)\n",
            argv[0]);
        return 1;
    }

    const char* csv_path = argv[1];
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

    /* ---- NUMA-aware allocation + first-touch ---- */
    double* A    = (double*)numa_alloc(bytes);
    double* B_cm = (double*)numa_alloc(bytes);
    double* B_rm = (double*)numa_alloc(bytes);
    double* C    = (double*)numa_alloc(bytes);

    printf("Initializing arrays (NUMA bind+touch) ...\n");
    ft_init_rm(A,    total, M, N, false);
    ft_init_cm(B_cm, total, M, N);
    ft_init_rm(B_rm, total, M, N, true);   /* row-major copy of B values */
    bind_and_touch(C, bytes);
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++) C[k] = 0.0;

    /* ---- reference checksum ---- */
    kernel_row_major(A, B_cm, C);
    double ref_cs = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:ref_cs)
    for (size_t k = 0; k < total; k++) ref_cs += C[k];

    kernel_all_rowmajor(A, B_rm, C);
    double ref_cs_rm = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:ref_cs_rm)
    for (size_t k = 0; k < total; k++) ref_cs_rm += C[k];

    /* ---- benchmarks ---- */
    printf("\n=== Naive schedules ===\n");
    bench("row_major",     0, kernel_row_major,    A, B_cm, C, ref_cs,    nthreads);
    bench("col_major",     0, kernel_col_major,    A, B_cm, C, ref_cs,    nthreads);
    bench("all_rowmajor",  0, kernel_all_rowmajor, A, B_rm, C, ref_cs_rm, nthreads);

    printf("\n=== Tiled schedules (stack-local B transpose) ===\n");
    bench("tiled",  8,  kernel_tiled<8>,   A, B_cm, C, ref_cs, nthreads);
    bench("tiled", 16,  kernel_tiled<16>,  A, B_cm, C, ref_cs, nthreads);
    bench("tiled", 32,  kernel_tiled<32>,  A, B_cm, C, ref_cs, nthreads);
    bench("tiled", 64,  kernel_tiled<64>,  A, B_cm, C, ref_cs, nthreads);
    bench("tiled", 128, kernel_tiled<128>, A, B_cm, C, ref_cs, nthreads);

    /* ---- write CSV ---- */
    FILE* fp = fopen(csv_path, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", csv_path); return 1; }
    fprintf(fp, "variant,M,N,tile,nthreads,rep,time_s,bw_gbs,checksum,status\n");
    for (auto& r : g_records)
        fprintf(fp, "%s,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e,%s\n",
                r.variant.c_str(), M, N, r.tile, r.nthreads,
                r.rep, r.time_s, r.bw_gbs, r.checksum, r.status);
    fclose(fp);
    printf("\nWrote %zu records to %s\n", g_records.size(), csv_path);

    /* ---- summary ---- */
    auto median_bw = [&](const char* name, int tile) -> double {
        std::vector<double> bws;
        for (auto& r : g_records)
            if (r.variant == name && r.tile == tile)
                bws.push_back(r.bw_gbs);
        if (bws.empty()) return 0;
        std::sort(bws.begin(), bws.end());
        return bws[bws.size()/2];
    };

    double best_tiled = 0;
    for (int t : {8, 16, 32, 64, 128})
        best_tiled = std::max(best_tiled, median_bw("tiled", t));
    double bw_row  = median_bw("row_major", 0);
    double bw_col  = median_bw("col_major", 0);
    double bw_ctrl = median_bw("all_rowmajor", 0);

    double B_eff = 64.0 / sizeof(double);
    printf("\n=== Summary ===\n");
    printf("  Best tiled:   %7.1f GB/s\n", best_tiled);
    printf("  row_major:    %7.1f GB/s  (%.1fx vs tiled)\n",
           bw_row, best_tiled / bw_row);
    printf("  col_major:    %7.1f GB/s  (%.1fx vs tiled)\n",
           bw_col, best_tiled / bw_col);
    printf("  all_rowmajor: %7.1f GB/s  (peak control)\n", bw_ctrl);
    printf("\n=== Model (B_eff=%.0f) ===\n", B_eff);
    printf("  row_major predicted:  %.1fx  (cost %g vs 3)\n",
           (2 + B_eff) / 3.0, 2 + B_eff);
    printf("  col_major predicted:  %.1fx  (cost %g vs 3)\n",
           (1 + 2*B_eff) / 3.0, 1 + 2*B_eff);

    /* cleanup */
    numa_free(A,    bytes);
    numa_free(B_cm, bytes);
    numa_free(B_rm, bytes);
    numa_free(C,    bytes);
    if (g_flush_buf) numa_free(g_flush_buf, FLUSH_N * sizeof(double));

    return 0;
}