/*
 * numa_calibrate.c -- Microbenchmark suite for cost-model parameter calibration
 *
 * Measures three parameters of the NUMA-aware average block distance:
 *   β  (page-locality radius)  -- detected from stride-sweep latency knee
 *   α  (page-locality discount) -- ratio of intra-page to inter-page latency
 *   γ  (NUMA penalty factor)    -- ratio of remote-NUMA to local-NUMA latency
 *
 * DESIGN:
 *   - Every access and every pointer-chase node sits on a cache-line boundary.
 *     The chase node is a 64-byte struct; the pointer occupies the first 8 bytes,
 *     the remaining 56 bytes are padding.  Stride is always a multiple of 64 B.
 *   - Arrays are sized to 4--8 GiB so they blow past LLC and TLB on all targets.
 *
 * Five benchmarks:
 *   1. stride_bw    -- Strided load bandwidth (1 CL per group, sweep stride).
 *   2. stride_lat   -- Pointer-chase latency at varying stride (β, α).
 *   3. numa_bw      -- Local vs remote NUMA triad bandwidth (γ_bw).
 *   4. numa_lat     -- Local vs remote NUMA pointer-chase latency (γ_lat, α).
 *   5. numa_matrix  -- Thread-to-owner bandwidth matrix.
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -o numa_calibrate numa_calibrate.c -lm
 *
 * Run:
 *   OMP_PLACES=cores OMP_PROC_BIND=close ./numa_calibrate all
 *
 * Related tools:
 *   - lmbench lat_mem_rd    (https://lmbench.sourceforge.net/man/lat_mem_rd.8.html)
 *   - torvalds/test-tlb     (https://github.com/torvalds/test-tlb)
 *   - pChase                (https://github.com/maleadt/pChase)
 *   - clamchowder/Microbenchmarks (https://github.com/clamchowder/Microbenchmarks)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sys/mman.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

/* ================================================================ */
/*  Constants                                                        */
/* ================================================================ */

#define CL_BYTES   64
#define CL_DOUBLES (CL_BYTES / (int)sizeof(double))   /* 8 */
#define PAGE_4K    4096
#define PAGE_64K   65536
#define PAGE_2M    (2*1024*1024)
#define GiB        (1UL << 32)

/* ================================================================ */
/*  Allocation                                                       */
/* ================================================================ */

static void *alloc_unfaulted(size_t bytes) {
    void *p = mmap(NULL, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); exit(1); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return p;
}

static void dealloc(void *p, size_t bytes) {
    if (p) munmap(p, bytes);
}

/* ================================================================ */
/*  Cache-line-aligned pointer-chase node                            */
/*                                                                   */
/*  Every node is exactly 64 bytes.  The pointer sits at offset 0.   */
/*  When we chase, each load fetches exactly one fresh cache line.   */
/* ================================================================ */

typedef struct __attribute__((aligned(CL_BYTES))) {
    void *next;
    char  pad[CL_BYTES - sizeof(void *)];
} cl_node_t;

_Static_assert(sizeof(cl_node_t) == CL_BYTES,
               "chase node must be exactly one cache line");

/* ================================================================ */
/*  Build a pointer chain with cache-line-granularity stride         */
/*                                                                   */
/*  arena         : base pointer (CL-aligned from mmap)              */
/*  arena_bytes   : total arena size                                 */
/*  stride_cls    : stride in cache lines (≥ 1)                      */
/*  Returns number of nodes in the chain.                            */
/* ================================================================ */

static size_t build_chase_chain(char *arena, size_t arena_bytes,
                                size_t stride_cls) {
    size_t stride_bytes = stride_cls * CL_BYTES;
    size_t n_nodes = arena_bytes / stride_bytes;
    if (n_nodes < 2) return 0;

    /* Sequential strided chain: node[i] → node[i+1] */
    for (size_t i = 0; i < n_nodes - 1; i++) {
        cl_node_t *cur  = (cl_node_t *)(arena + i * stride_bytes);
        cl_node_t *nxt  = (cl_node_t *)(arena + (i + 1) * stride_bytes);
        cur->next = nxt;
    }
    /* Wrap last → first */
    cl_node_t *last = (cl_node_t *)(arena + (n_nodes - 1) * stride_bytes);
    cl_node_t *first = (cl_node_t *)arena;
    last->next = first;

    return n_nodes;
}

/* Chase `count` pointers starting from `start`, return endpoint.
 * Every load is serialized by the asm barrier. */
static inline void *do_chase(void *start, int count) {
    void *p = start;
    for (int i = 0; i < count; i++) {
        p = *(void **)p;
        __asm__ volatile("" : "+r"(p));
    }
    return p;
}

/* Measure pointer-chase latency in ns/hop.  5 reps, return best. */
static double chase_latency_ns(char *arena, size_t arena_bytes,
                                size_t stride_cls, int chases) {
    size_t n = build_chase_chain(arena, arena_bytes, stride_cls);
    if (n == 0 || n < (size_t)chases) return -1.0;

    /* Warm up */
    void *p = do_chase(arena, chases);

    double best = 1e30;
    for (int rep = 0; rep < 5; rep++) {
        double t0 = omp_get_wtime();
        p = do_chase(p, chases);
        double t1 = omp_get_wtime();
        double ns_per = (t1 - t0) * 1e9 / chases;
        if (ns_per < best) best = ns_per;
    }
    /* Prevent p from being optimized away across calls */
    if (best < 0) printf("%p", p);
    return best;
}

/* ================================================================ */
/*  Benchmark 1: Strided Load Bandwidth                              */
/*                                                                   */
/*  For each stride S (in CLs), load one cache line (8 doubles),     */
/*  skip S-1 lines, repeat over a 4 GiB array.                      */
/*  TLB working set explodes at the page boundary → BW drops → β.   */
/* ================================================================ */

static void bench_stride_bw(void) {
    const size_t BYTES = 8UL * GiB;
    const size_t N = BYTES / sizeof(double);
    const int REPS = 100;
    double *arr = (double *)alloc_unfaulted(BYTES);

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < N; i++) arr[i] = 1.0;

    printf("=== Benchmark 1: Strided Load Bandwidth ===\n");
    printf("  Array: %.1f GiB, Reps: %d\n", (double)BYTES / GiB, REPS);
    printf("%12s %12s %12s %12s\n",
           "stride_CL", "stride_B", "eff_BW_GB/s", "ns_per_CL");

    for (int stride_cl = 1; stride_cl <= 32768; stride_cl *= 2) {
        size_t stride_elem = (size_t)stride_cl * CL_DOUBLES;
        size_t n_groups = N / stride_elem;
        if (n_groups < 256) break;
        size_t bytes_loaded = n_groups * CL_BYTES;

        double best = 1e30;
        for (int rep = 0; rep < REPS; rep++) {
            double t0 = omp_get_wtime();

            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                int nth = omp_get_num_threads();
                /* Each thread gets a CL-aligned range of groups */
                size_t g_lo = (n_groups * tid / nth);
                size_t g_hi = (n_groups * (tid + 1) / nth);
                double sink = 0.0;

                for (size_t g = g_lo; g < g_hi; g++) {
                    size_t base = g * stride_elem;
                    /* Load exactly 1 cache line = 8 doubles */
                    sink += arr[base+0] + arr[base+1]
                          + arr[base+2] + arr[base+3]
                          + arr[base+4] + arr[base+5]
                          + arr[base+6] + arr[base+7];
                }
                if (sink == -999.999) printf("%f", sink);
            }

            double dt = omp_get_wtime() - t0;
            if (dt < best) best = dt;
        }

        double bw = (double)bytes_loaded / best / 1e9;
        double ns = best / (double)n_groups * 1e9;
        printf("%12d %12zu %12.2f %12.2f\n",
               stride_cl, (size_t)stride_cl * CL_BYTES, bw, ns);
    }
    printf("\n");
    dealloc(arr, BYTES);
}

/* ================================================================ */
/*  Benchmark 2: Pointer-Chase Latency (stride sweep)                */
/*                                                                   */
/*  Single-thread chase through a 4 GiB arena.  Stride sweeps from   */
/*  1 CL to 64K CLs (= 4 MiB).  Every node is CL-aligned.          */
/*  Reveals L1/L2/LLC/DRAM/TLB tiers → α = same-page / cross-page.  */
/* ================================================================ */

static void bench_stride_lat(void) {
    const size_t ARENA = 4UL * GiB;
    const int CHASES = 1 << 22;   /* 4M chases: long chain for stable ns */
    const int REPS = 5;

    char *arena = (char *)alloc_unfaulted(ARENA);
    /* Single-thread first-touch for controlled NUMA placement */
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < ARENA / sizeof(double); i++)
        ((double *)arena)[i] = 0.0;

    printf("=== Benchmark 2: Pointer-Chase Latency (stride sweep) ===\n");
    printf("  Arena: %.1f GiB, Chases: %d per measurement\n",
           (double)ARENA / GiB, CHASES);
    printf("%12s %12s %12s %12s\n",
           "stride_CL", "stride_B", "lat_ns", "tier");

    /* Sweep stride from 1 CL (64B) to 64K CLs (4 MiB) */
    for (size_t scl = 1; scl <= 65536; scl *= 2) {
        size_t n_nodes = ARENA / (scl * CL_BYTES);
        if (n_nodes < (size_t)CHASES) break;

        double lat = chase_latency_ns(arena, ARENA, scl, CHASES);

        const char *tier = "???";
        if (lat < 3)        tier = "L1";
        else if (lat < 15)  tier = "L2";
        else if (lat < 60)  tier = "L3/LLC";
        else if (lat < 120) tier = "DRAM-inpg";
        else                tier = "DRAM-TLBmiss";

        printf("%12zu %12zu %12.2f %12s\n",
               scl, scl * CL_BYTES, lat, tier);
    }
    printf("\n");
    dealloc(arena, ARENA);
}

/* ================================================================ */
/*  Benchmark 3: Local vs Remote NUMA Bandwidth                      */
/*                                                                   */
/*  Thread 0 runs STREAM triad on arrays first-touched by itself     */
/*  (local) vs arrays first-touched by the last thread (remote).     */
/*  Each array: 1 GiB.  Ratio BW_local / BW_remote → γ_bw.          */
/* ================================================================ */

static double triad_bw_single_thread(double *a, double *b, double *c,
                                      size_t N, int reps) {
    const double scalar = 3.0;
    for (size_t i = 0; i < N; i++) a[i] = b[i] + scalar * c[i];

    double best = 1e30;
    for (int r = 0; r < reps; r++) {
        double t0 = omp_get_wtime();
        for (size_t i = 0; i < N; i++)
            a[i] = b[i] + scalar * c[i];
        double dt = omp_get_wtime() - t0;
        if (dt < best) best = dt;
    }
    return 3.0 * (double)N * sizeof(double) / best / 1e9;
}

static void bench_numa_bw(void) {
    const size_t BYTES = 1UL * GiB;
    const size_t N = BYTES / sizeof(double);
    const int REPS = 30;
    int nth = omp_get_max_threads();

    printf("=== Benchmark 3: Local vs Remote NUMA Bandwidth ===\n");
    printf("  Per-array: %.1f GiB, Threads: %d\n\n",
           (double)BYTES / GiB, nth);

    double *a_loc = (double *)alloc_unfaulted(BYTES);
    double *b_loc = (double *)alloc_unfaulted(BYTES);
    double *c_loc = (double *)alloc_unfaulted(BYTES);
    double *a_rem = (double *)alloc_unfaulted(BYTES);
    double *b_rem = (double *)alloc_unfaulted(BYTES);
    double *c_rem = (double *)alloc_unfaulted(BYTES);

    /* First-touch: local by thread 0, remote by last thread */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        if (tid == 0)
            for (size_t i = 0; i < N; i++) {
                a_loc[i] = 0; b_loc[i] = 1; c_loc[i] = 2;
            }
        if (tid == nth - 1)
            for (size_t i = 0; i < N; i++) {
                a_rem[i] = 0; b_rem[i] = 1; c_rem[i] = 2;
            }
    }

    double bw_local = 0, bw_remote = 0;
    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0) {
            bw_local  = triad_bw_single_thread(a_loc, b_loc, c_loc, N, REPS);
            bw_remote = triad_bw_single_thread(a_rem, b_rem, c_rem, N, REPS);
        }
    }

    printf("  Local  BW (thread 0 → own NUMA):    %8.2f GB/s\n", bw_local);
    printf("  Remote BW (thread 0 → remote NUMA): %8.2f GB/s\n", bw_remote);
    printf("  γ_bw = local / remote:              %8.3f\n",
           bw_local / bw_remote);
    printf("\n  NOTE: If γ ≈ 1, threads 0 and %d share a NUMA domain.\n"
           "        Ensure OMP_PLACES separates them.\n\n", nth - 1);

    dealloc(a_loc, BYTES); dealloc(b_loc, BYTES); dealloc(c_loc, BYTES);
    dealloc(a_rem, BYTES); dealloc(b_rem, BYTES); dealloc(c_rem, BYTES);
}

/* ================================================================ */
/*  Benchmark 4: Local vs Remote NUMA Pointer-Chase Latency          */
/*                                                                   */
/*  Thread 0 chases through a 4 GiB arena placed on its own vs       */
/*  remote NUMA node.  Chase stride = 1 page (64 or 1024 CLs).      */
/*  Also sweeps stride on local arena for α calibration.             */
/* ================================================================ */

static void bench_numa_lat(void) {
    const size_t ARENA = 4UL * GiB;
    const int CHASES = 1 << 22;    /* 4M chases */
    int nth = omp_get_max_threads();

    /* Detect page size to choose a sensible default stride */
    long ps = sysconf(_SC_PAGESIZE);
    if (ps <= 0) ps = PAGE_4K;
    size_t page_cls = (size_t)ps / CL_BYTES;

    printf("=== Benchmark 4: Local vs Remote NUMA Latency ===\n");
    printf("  Arena: %.1f GiB, Page: %ld B (%zu CLs), Chases: %d\n\n",
           (double)ARENA / GiB, ps, page_cls, CHASES);

    char *arena_local  = (char *)alloc_unfaulted(ARENA);
    char *arena_remote = (char *)alloc_unfaulted(ARENA);

    /* First-touch: thread 0 → local, last thread → remote */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        if (tid == 0) {
            #pragma omp parallel for schedule(static) num_threads(1)
            for (size_t i = 0; i < ARENA / sizeof(double); i++)
                ((double *)arena_local)[i] = 0.0;
        }
        #pragma omp barrier
        if (tid == nth - 1) {
            for (size_t i = 0; i < ARENA / sizeof(double); i++)
                ((double *)arena_remote)[i] = 0.0;
        }
    }

    double lat_local = 0, lat_remote = 0;

    /* Thread 0 chases on both arenas at page-stride */
    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0) {
            lat_local  = chase_latency_ns(arena_local,  ARENA,
                                           page_cls, CHASES);
            lat_remote = chase_latency_ns(arena_remote, ARENA,
                                           page_cls, CHASES);
        }
    }

    printf("  Chase stride: %zu CLs = %zu B (1 page)\n", page_cls,
           page_cls * CL_BYTES);
    printf("  Local  latency:  %8.2f ns\n", lat_local);
    printf("  Remote latency:  %8.2f ns\n", lat_remote);
    printf("  γ_lat = remote / local: %8.3f\n\n", lat_remote / lat_local);

    /* Stride sweep on local arena for α calibration */
    printf("  --- Stride sweep on local NUMA (α calibration) ---\n");
    printf("  %12s %12s %12s %12s\n",
           "stride_CL", "stride_B", "lat_ns", "ratio_vs_1CL");

    double lat_1cl = 0;
    for (size_t scl = 1; scl <= 65536; scl *= 2) {
        if (ARENA / (scl * CL_BYTES) < (size_t)CHASES) break;

        double lat;
        #pragma omp parallel
        {
            if (omp_get_thread_num() == 0)
                lat = chase_latency_ns(arena_local, ARENA, scl, CHASES);
        }

        if (scl == 1) lat_1cl = lat;
        printf("  %12zu %12zu %12.2f %12.3f\n",
               scl, scl * CL_BYTES, lat,
               lat_1cl > 0 ? lat / lat_1cl : 0.0);
    }

    printf("\n  α ≈ lat(stride < β) / lat(stride = β)\n");
    printf("  γ ≈ lat_remote / lat_local  (at stride = page)\n\n");

    dealloc(arena_local,  ARENA);
    dealloc(arena_remote, ARENA);
}

/* ================================================================ */
/*  Benchmark 5: NUMA BW Matrix (thread → owner)                     */
/*                                                                   */
/*  Per-thread arrays of 512 MiB.  Thread 0 triads on arrays owned   */
/*  by selected threads across the machine.                          */
/* ================================================================ */

static void bench_numa_matrix(void) {
    int nth = omp_get_max_threads();
    const size_t BYTES = 512UL * (1UL << 20);  /* 512 MiB per array */
    const size_t N = BYTES / sizeof(double);
    const int REPS = 20;

    double **as = (double **)calloc(nth, sizeof(double *));
    double **bs = (double **)calloc(nth, sizeof(double *));
    double **cs = (double **)calloc(nth, sizeof(double *));

    for (int t = 0; t < nth; t++) {
        as[t] = (double *)alloc_unfaulted(BYTES);
        bs[t] = (double *)alloc_unfaulted(BYTES);
        cs[t] = (double *)alloc_unfaulted(BYTES);
    }

    /* Each thread first-touches its own arrays */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (size_t i = 0; i < N; i++) {
            as[tid][i] = 0; bs[tid][i] = 1; cs[tid][i] = 2;
        }
    }

    printf("=== Benchmark 5: NUMA BW Matrix (thread → owner) ===\n");
    printf("  Per-array: %zu MiB, Threads: %d\n",
           BYTES >> 20, nth);
    printf("  Thread 0 accesses arrays owned by selected threads.\n\n");
    printf("  %8s %12s\n", "owner", "BW_GB/s");

    int step = (nth > 16) ? nth / 8 : 1;
    if (step < 1) step = 1;
    for (int owner = 0; owner < nth; owner += step) {
        double bw = 0;
        #pragma omp parallel
        {
            if (omp_get_thread_num() == 0)
                bw = triad_bw_single_thread(as[owner], bs[owner],
                                             cs[owner], N, REPS);
        }
        printf("  %8d %12.2f%s\n", owner, bw,
               owner == 0 ? "  (local)" : "");
    }
    printf("\n");

    for (int t = 0; t < nth; t++) {
        dealloc(as[t], BYTES); dealloc(bs[t], BYTES); dealloc(cs[t], BYTES);
    }
    free(as); free(bs); free(cs);
}

/* ================================================================ */
/*  Main                                                             */
/* ================================================================ */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s <benchmark>\n\n"
        "Benchmarks:\n"
        "  stride_bw   Strided load bandwidth (β detection)\n"
        "  stride_lat  Pointer-chase latency sweep (α calibration)\n"
        "  numa_bw     Local vs remote NUMA bandwidth (γ_bw)\n"
        "  numa_lat    Local vs remote NUMA latency + stride sweep (α, γ)\n"
        "  numa_matrix Thread-to-owner BW matrix\n"
        "  all         Run all benchmarks\n",
        prog);
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    const char *mode = argv[1];

    long ps = sysconf(_SC_PAGESIZE);
    printf("NUMA Cost-Model Calibration Suite\n");
    printf("  OMP threads : %d\n", omp_get_max_threads());
    printf("  Cache line  : %d B\n", CL_BYTES);
    printf("  OS page     : %ld B (%ld CLs)\n", ps, ps / CL_BYTES);
    printf("  Chase node  : %zu B (= 1 CL)\n\n", sizeof(cl_node_t));

    int run_all = (strcmp(mode, "all") == 0);

    if (run_all || strcmp(mode, "stride_bw") == 0)   bench_stride_bw();
    if (run_all || strcmp(mode, "stride_lat") == 0)   bench_stride_lat();
    if (run_all || strcmp(mode, "numa_bw") == 0)      bench_numa_bw();
    if (run_all || strcmp(mode, "numa_lat") == 0)      bench_numa_lat();
    if (run_all || strcmp(mode, "numa_matrix") == 0)   bench_numa_matrix();

    if (!run_all
        && strcmp(mode, "stride_bw") != 0
        && strcmp(mode, "stride_lat") != 0
        && strcmp(mode, "numa_bw") != 0
        && strcmp(mode, "numa_lat") != 0
        && strcmp(mode, "numa_matrix") != 0) {
        fprintf(stderr, "Unknown benchmark: %s\n", mode);
        usage(argv[0]);
        return 1;
    }
    return 0;
}