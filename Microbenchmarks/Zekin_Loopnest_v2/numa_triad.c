#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <omp.h>

/*
 * NUMA-aware STREAM Triad via first-touch only.
 * No libnuma -- placement is entirely driven by:
 *   1. mmap to reserve pages without faulting them
 *   2. parallel first-touch init so each thread faults its chunk
 *      on the NUMA node it is pinned to via OMP_PLACES
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -o numa_triad numa_triad.c
 *
 * Run (4x Grace, 72 cores/socket):
 *   OMP_NUM_THREADS=288                                             \
 *   OMP_PLACES="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"         \
 *   OMP_PROC_BIND=close                                             \
 *   ./numa_triad [N]
 *
 * Run (single Grace, 72 cores, 1 NUMA node):
 *   OMP_NUM_THREADS=72 OMP_PLACES=cores OMP_PROC_BIND=close        \
 *   numactl --cpunodebind=0 --membind=0 ./numa_triad [N]
 */

#ifndef N_DEFAULT
#define N_DEFAULT  (1L << 30)   /* ~1B doubles = 8 GiB per array */
#endif
#define NTIMES     100
#define SCALAR     3.0

/* Allocate virtual pages via mmap but do NOT fault them.
 * The first parallel touch decides NUMA placement. */
static double *alloc_unfaulted(size_t n)
{
    size_t bytes = n * sizeof(double);
    void *p = mmap(NULL, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); exit(1); }
    madvise(p, bytes, MADV_HUGEPAGE);   /* <-- add this */
    return (double *)p;
}

static void free_array(double *p, size_t n)
{
    munmap(p, n * sizeof(double));
}

int main(int argc, char **argv)
{
    size_t N = N_DEFAULT;
    if (argc > 1) N = atol(argv[1]);

    printf("NUMA Triad (first-touch, no libnuma)\n");
    printf("  Elements   : %zu (%.2f GiB per array)\n",
           N, (double)(N * sizeof(double)) / (1L << 30));
    printf("  OMP threads: %d\n", omp_get_max_threads());
    printf("  Repetitions: %d\n\n", NTIMES);

    /* Reserve virtual address space -- no pages faulted yet. */
    double *a = alloc_unfaulted(N);
    double *b = alloc_unfaulted(N);
    double *c = alloc_unfaulted(N);

    /* ----------------------------------------------------------
     * First-touch initialization.
     * schedule(static) guarantees that thread T always gets the
     * same iteration range, so the pages it faults here are the
     * same pages it will access in the triad loop.  Combined
     * with OMP_PROC_BIND=close this places each page on the
     * NUMA node where it will be consumed.
     * ---------------------------------------------------------- */
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < N; i++) {
        a[i] = 0.0;
        b[i] = 1.0;
        c[i] = 2.0;
    }

    /* warm-up */
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < N; i++)
        a[i] = b[i] + SCALAR * c[i];

    double tmin = 1e30, tmax = 0.0, tavg = 0.0;

    for (int t = 0; t < NTIMES; t++) {
        /* Reset a[] with the SAME static schedule so pages stay local. */
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < N; i++)
            a[i] = 0.0;

        double t0, t1;

        #pragma omp parallel
        {
            #pragma omp barrier
            #pragma omp master
            t0 = omp_get_wtime();
            #pragma omp barrier

            #pragma omp for schedule(static)
            for (size_t i = 0; i < N; i++)
                a[i] = b[i] + SCALAR * c[i];

            #pragma omp barrier
            #pragma omp master
            t1 = omp_get_wtime();
        }

        double dt = t1 - t0;
        if (dt < tmin) tmin = dt;
        if (dt > tmax) tmax = dt;
        tavg += dt;
    }
    tavg /= NTIMES;

    /* 2 reads + 1 write = 3 streams */
    double bytes = 3.0 * (double)N * sizeof(double);

    printf("Results (Triad: a = b + s*c)\n");
    printf("  Best  : %8.4f ms  ->  %8.2f GB/s\n",
           tmin * 1e3, bytes / tmin / 1e9);
    printf("  Avg   : %8.4f ms  ->  %8.2f GB/s\n",
           tavg * 1e3, bytes / tavg / 1e9);
    printf("  Worst : %8.4f ms  ->  %8.2f GB/s\n",
           tmax * 1e3, bytes / tmax / 1e9);

    /* validate */
    double expected = 1.0 + SCALAR * 2.0;
    int ok = 1;
    for (size_t i = 0; i < N; i++) {
        double err = a[i] - expected;
        if (err > 1e-12 || err < -1e-12) {
            fprintf(stderr, "VALIDATION FAILED i=%zu: got %.15f expected %.15f\n",
                    i, a[i], expected);
            ok = 0;
            break;
        }
    }
    if (ok) printf("\nVALIDATION PASSED\n");

    free_array(a, N);
    free_array(b, N);
    free_array(c, N);
    return 0;
}