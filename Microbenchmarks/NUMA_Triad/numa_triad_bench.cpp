/*
 * numa_triad_bench.cpp -- NUMA-aware triad benchmark
 *
 * Runs C[i] = A[i] + s * B[i] under multiple (touch, compute) schedule
 * combinations.  The first-touch schedule determines NUMA page placement;
 * the compute schedule determines which thread accesses which elements.
 * When they mismatch, cross-NUMA traffic degrades bandwidth.
 *
 * Configurations:
 *   touch_block  + compute_block   → local (baseline)
 *   touch_block  + compute_shift1  → all remote (shifted by 1 NUMA domain)
 *   touch_block  + compute_interl  → mixed (interleaved chunks)
 *   touch_interl + compute_block   → mixed (reverse mismatch)
 *
 * Each configuration is parameterized by chunk size (in cache lines).
 *
 * Build: g++ -O3 -std=c++17 -march=native -fopenmp -o numa_triad_bench numa_triad_bench.cpp
 * Run:   OMP_PLACES=cores OMP_PROC_BIND=close ./numa_triad_bench [N_GiB]
 * Output: CSV to stdout, progress to stderr
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <sys/mman.h>
#include <unistd.h>
#include <omp.h>

#define CL_BYTES 64
#define CL_DOUBLES (CL_BYTES / (int)sizeof(double))

static void *alloc_unfaulted(size_t bytes) {
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); exit(1); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return p;
}

static void dealloc(void *p, size_t bytes) { if (p) munmap(p, bytes); }

/* ================================================================
 *  Schedule descriptors
 *
 *  Given N elements and P threads, thread_range(tid, P, N) returns
 *  the set of element indices that thread tid touches/computes.
 *  We encode this as (lo, hi, stride) triples per thread.
 *
 *  For simplicity, we implement the schedule inline in the loops
 *  using a "who owns element i" function.
 * ================================================================ */

enum SchedType {
    SCHED_BLOCK,      /* schedule(static): contiguous N/P blocks         */
    SCHED_SHIFT1,     /* like block, but shifted by 1 NUMA domain        */
    SCHED_INTERLEAVE, /* round-robin of chunk_size elements              */
};

static const char* sched_label[] = {"block", "shift1", "interleave"};

/* Return the thread that owns element i under the given schedule.
 * P = number of threads, chunk = interleave chunk in elements. */
inline int owner(SchedType s, int i, int N, int P, int chunk) {
    switch (s) {
    case SCHED_BLOCK: {
        int blk = (N + P - 1) / P;
        int t = i / blk;
        return (t < P) ? t : P - 1;
    }
    case SCHED_SHIFT1: {
        int blk = (N + P - 1) / P;
        int t = i / blk;
        t = (t < P) ? t : P - 1;
        return (t + 1) % P;  /* shift by 1 domain */
    }
    case SCHED_INTERLEAVE: {
        int c = (chunk > 0) ? chunk : 1;
        int group = i / c;
        return group % P;
    }
    }
    return 0;
}

/* ================================================================
 *  First-touch: thread tid writes 0.0 to all elements it owns.
 * ================================================================ */

static void first_touch(double *arr, size_t N, SchedType sched, int chunk) {
    int P = omp_get_max_threads();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (size_t i = 0; i < N; i++) {
            if (owner(sched, (int)i, (int)N, P, chunk) == tid)
                arr[i] = 0.0;
        }
    }
}

/* Faster: block first-touch (the common case) */
static void first_touch_block(double *arr, size_t N) {
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nth = omp_get_num_threads();
        size_t lo = (size_t)N * tid / nth;
        size_t hi = (size_t)N * (tid + 1) / nth;
        for (size_t i = lo; i < hi; i++) arr[i] = 0.0;
    }
}

static void first_touch_shift1(double *arr, size_t N) {
    int P = omp_get_max_threads();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        /* Thread tid touches the block that will be owned by (tid-1+P)%P */
        int src = (tid + 1) % P;  /* shift: thread tid touches src's block */
        size_t lo = (size_t)N * src / P;
        size_t hi = (size_t)N * (src + 1) / P;
        for (size_t i = lo; i < hi; i++) arr[i] = 0.0;
    }
}

static void first_touch_interleave(double *arr, size_t N, int chunk_elems) {
    int P = omp_get_max_threads();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (size_t i = 0; i < N; i++) {
            int group = (int)(i / chunk_elems);
            if (group % P == tid)
                arr[i] = 0.0;
        }
    }
}

static void do_first_touch(double *arr, size_t N, SchedType s, int chunk_elems) {
    switch (s) {
    case SCHED_BLOCK:      first_touch_block(arr, N); break;
    case SCHED_SHIFT1:     first_touch_shift1(arr, N); break;
    case SCHED_INTERLEAVE: first_touch_interleave(arr, N, chunk_elems); break;
    }
}

/* ================================================================
 *  Compute: triad C = A + s*B under a given schedule
 * ================================================================ */

static double run_triad(double *C, const double *A, const double *B,
                         size_t N, SchedType compute_sched,
                         int chunk_elems, int reps) {
    const double s = 3.0;

    /* Warmup */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nth = omp_get_num_threads();
        size_t lo = (size_t)N * tid / nth;
        size_t hi = (size_t)N * (tid + 1) / nth;
        for (size_t i = lo; i < hi; i++)
            C[i] = A[i] + s * B[i];
    }

    double best = 1e30;
    for (int r = 0; r < reps; r++) {
        double t0, t1;

        switch (compute_sched) {
        case SCHED_BLOCK:
            #pragma omp parallel
            {
                #pragma omp barrier
                #pragma omp master
                t0 = omp_get_wtime();
                #pragma omp barrier

                int tid = omp_get_thread_num();
                int nth = omp_get_num_threads();
                size_t lo = (size_t)N * tid / nth;
                size_t hi = (size_t)N * (tid + 1) / nth;
                for (size_t i = lo; i < hi; i++)
                    C[i] = A[i] + s * B[i];

                #pragma omp barrier
                #pragma omp master
                t1 = omp_get_wtime();
            }
            break;

        case SCHED_SHIFT1:
            #pragma omp parallel
            {
                #pragma omp barrier
                #pragma omp master
                t0 = omp_get_wtime();
                #pragma omp barrier

                int tid = omp_get_thread_num();
                int nth = omp_get_num_threads();
                /* Thread tid computes on the block of thread (tid+1)%P */
                int src = (tid + 1) % nth;
                size_t lo = (size_t)N * src / nth;
                size_t hi = (size_t)N * (src + 1) / nth;
                for (size_t i = lo; i < hi; i++)
                    C[i] = A[i] + s * B[i];

                #pragma omp barrier
                #pragma omp master
                t1 = omp_get_wtime();
            }
            break;

        case SCHED_INTERLEAVE:
            #pragma omp parallel
            {
                #pragma omp barrier
                #pragma omp master
                t0 = omp_get_wtime();
                #pragma omp barrier

                int tid = omp_get_thread_num();
                int nth = omp_get_num_threads();
                for (size_t base = (size_t)tid * chunk_elems;
                     base < N;
                     base += (size_t)nth * chunk_elems) {
                    size_t end = base + chunk_elems;
                    if (end > N) end = N;
                    for (size_t i = base; i < end; i++)
                        C[i] = A[i] + s * B[i];
                }

                #pragma omp barrier
                #pragma omp master
                t1 = omp_get_wtime();
            }
            break;
        }

        double dt = t1 - t0;
        if (dt < best) best = dt;
    }

    /* 2 reads + 1 write = 3N doubles */
    return 3.0 * (double)N * sizeof(double) / best / 1e9;
}

/* ================================================================
 *  Configuration table
 * ================================================================ */

struct Config {
    const char* name;
    SchedType   touch;
    SchedType   compute;
    int         touch_chunk_cls;   /* chunk in cache lines (0 = N/P) */
    int         compute_chunk_cls;
};

static const Config configs[] = {
    /* Baseline: touch=block, compute=block → all local */
    {"local",                SCHED_BLOCK,      SCHED_BLOCK,      0, 0},

    /* All remote: touch=block, compute=shifted → every access crosses NUMA */
    {"all_remote",           SCHED_BLOCK,      SCHED_SHIFT1,     0, 0},

    /* Touch interleaved (1 page), compute block → mixed */
    {"touch_interl_1pg",     SCHED_INTERLEAVE, SCHED_BLOCK,      0, 0},

    /* Touch interleaved (16 pages), compute block → less mixing */
    {"touch_interl_16pg",    SCHED_INTERLEAVE, SCHED_BLOCK,      0, 0},

    /* Touch block, compute interleaved (1 page) → mixed */
    {"compute_interl_1pg",   SCHED_BLOCK,      SCHED_INTERLEAVE, 0, 0},

    /* Touch block, compute interleaved (16 pages) → less mixing */
    {"compute_interl_16pg",  SCHED_BLOCK,      SCHED_INTERLEAVE, 0, 0},
};
static constexpr int N_CONFIGS = sizeof(configs) / sizeof(configs[0]);

/* ================================================================
 *  Main
 * ================================================================ */

int main(int argc, char **argv) {
    double gib = (argc > 1) ? atof(argv[1]) : 2.0;
    int reps   = (argc > 2) ? atoi(argv[2]) : 50;
    size_t N = (size_t)(gib * (1UL << 30) / sizeof(double));
    size_t bytes = N * sizeof(double);

    long page_bytes = sysconf(_SC_PAGESIZE);
    int page_cls = (int)(page_bytes / CL_BYTES);
    int P = omp_get_max_threads();

    fprintf(stderr, "NUMA Triad Benchmark\n");
    fprintf(stderr, "  N = %zu (%.2f GiB per array)\n", N, gib);
    fprintf(stderr, "  Threads = %d, Page = %ld B (%d CLs)\n",
            P, page_bytes, page_cls);
    fprintf(stderr, "  Reps = %d\n\n", reps);

    /* CSV header */
    printf("config,touch_sched,compute_sched,touch_chunk_cls,compute_chunk_cls,"
           "bw_gbs,N,P,page_cls\n");

    for (int ci = 0; ci < N_CONFIGS; ci++) {
        Config cfg = configs[ci];

        /* Resolve chunk sizes: 0 means N/P, otherwise scale by page */
        int t_chunk_cls = cfg.touch_chunk_cls;
        int c_chunk_cls = cfg.compute_chunk_cls;

        /* Set page-based chunk sizes for interleave configs */
        if (cfg.touch == SCHED_INTERLEAVE && t_chunk_cls == 0) {
            if (strstr(cfg.name, "1pg"))  t_chunk_cls = page_cls;
            if (strstr(cfg.name, "16pg")) t_chunk_cls = 16 * page_cls;
        }
        if (cfg.compute == SCHED_INTERLEAVE && c_chunk_cls == 0) {
            if (strstr(cfg.name, "1pg"))  c_chunk_cls = page_cls;
            if (strstr(cfg.name, "16pg")) c_chunk_cls = 16 * page_cls;
        }

        int t_chunk_elems = (t_chunk_cls > 0) ? t_chunk_cls * CL_DOUBLES : 0;
        int c_chunk_elems = (c_chunk_cls > 0) ? c_chunk_cls * CL_DOUBLES : 0;

        fprintf(stderr, "  %-24s touch=%-12s compute=%-12s tchunk=%d cchunk=%d ...",
                cfg.name, sched_label[cfg.touch], sched_label[cfg.compute],
                t_chunk_cls, c_chunk_cls);

        /* Allocate fresh arrays for each config (clean NUMA state) */
        double *A = (double *)alloc_unfaulted(bytes);
        double *B = (double *)alloc_unfaulted(bytes);
        double *C = (double *)alloc_unfaulted(bytes);

        /* First-touch */
        do_first_touch(A, N, cfg.touch, t_chunk_elems);
        do_first_touch(B, N, cfg.touch, t_chunk_elems);
        do_first_touch(C, N, cfg.touch, t_chunk_elems);

        /* Initialize A, B with values (same touch pattern) */
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nth = omp_get_num_threads();
            size_t lo = (size_t)N * tid / nth;
            size_t hi = (size_t)N * (tid + 1) / nth;
            for (size_t i = lo; i < hi; i++) {
                A[i] = 1.0;
                B[i] = 2.0;
            }
        }

        /* Run triad under the compute schedule */
        double bw = run_triad(C, A, B, N, cfg.compute, c_chunk_elems, reps);

        fprintf(stderr, " BW = %.2f GB/s\n", bw);

        printf("%s,%s,%s,%d,%d,%.4f,%zu,%d,%d\n",
               cfg.name, sched_label[cfg.touch], sched_label[cfg.compute],
               t_chunk_cls, c_chunk_cls, bw, N, P, page_cls);

        dealloc(A, bytes);
        dealloc(B, bytes);
        dealloc(C, bytes);
    }

    return 0;
}