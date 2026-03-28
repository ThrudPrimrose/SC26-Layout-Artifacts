/*
 * numa_rotation_bench.c
 *
 * NUMA round-robin bandwidth benchmark.
 *
 * Setup:
 *   P_N NUMA domains, each with arrays A,B,C first-touched locally.
 *   N_threads = P_N * threads_per_domain.
 *
 * Phases (0 .. P_N-1):
 *   Phase k: domain d's threads operate on arrays owned by domain (d+k) % P_N.
 *   Phase 0 = fully local.
 *   Phase P_N/2 = maximally remote (on 2-socket).
 *
 * Kernel: STREAM triad  A[i] = B[i] + scalar * C[i]
 *
 * Compile:
 *   gcc -O3 -fopenmp -march=native -o numa_rotation_bench numa_rotation_bench.c -lm -lnuma
 *   (without libnuma: compile with -DNO_LIBNUMA, uses first-touch instead)
 *
 * Run:
 *   OMP_PLACES=cores OMP_PROC_BIND=close ./numa_rotation_bench [N_elems] [iters] [csv_out]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#ifndef NO_LIBNUMA
#include <numa.h>
#endif

/* ─── defaults ─── */
#define DEFAULT_N      (1 << 24)   /* 16M doubles = 128 MiB per array */
#define DEFAULT_ITERS  50
#define DEFAULT_CSV    "numa_rotation.csv"
#define SCALAR         3.141592653589793
#define WARMUP         5

/* ─── NUMA helpers ─── */

static int detect_numa_domains(void) {
#ifndef NO_LIBNUMA
    if (numa_available() >= 0)
        return numa_max_node() + 1;
#endif
    /* Fallback: parse /sys */
    FILE *f = fopen("/sys/devices/system/node/online", "r");
    if (!f) return 1;
    int lo, hi;
    if (fscanf(f, "%d-%d", &lo, &hi) == 2) { fclose(f); return hi - lo + 1; }
    fclose(f);
    return 1;
}

static int thread_to_domain(int tid, int nth, int P_N) {
    /* Assumes OMP_PLACES=cores OMP_PROC_BIND=close → threads are packed
       contiguously, so threads [0..nth/P_N) → domain 0, etc. */
    int per = (nth + P_N - 1) / P_N;
    int d = tid / per;
    return (d < P_N) ? d : P_N - 1;
}

/* ─── Allocation with first-touch placement ─── */

typedef struct {
    double *A, *B, *C;
    size_t N;           /* elements per array */
    int    domain;      /* owning NUMA domain */
} DomainArrays;

static DomainArrays* alloc_domain_arrays(int P_N, size_t N, int nth) {
    DomainArrays *da = (DomainArrays *)calloc(P_N, sizeof(DomainArrays));
    for (int d = 0; d < P_N; d++) {
        da[d].N = N;
        da[d].domain = d;
        da[d].A = (double *)malloc(N * sizeof(double));
        da[d].B = (double *)malloc(N * sizeof(double));
        da[d].C = (double *)malloc(N * sizeof(double));
    }

    /* First-touch: each domain's arrays touched by that domain's threads */
    #pragma omp parallel num_threads(nth)
    {
        int tid = omp_get_thread_num();
        int my_domain = thread_to_domain(tid, nth, P_N);
        int tpd = (nth + P_N - 1) / P_N;  /* threads per domain */
        int local_tid = tid - my_domain * tpd;
        if (local_tid < 0) local_tid = 0;

        size_t chunk = (N + tpd - 1) / tpd;
        size_t lo = (size_t)local_tid * chunk;
        size_t hi = lo + chunk;
        if (hi > N) hi = N;

        double *A = da[my_domain].A;
        double *B = da[my_domain].B;
        double *C = da[my_domain].C;
        for (size_t i = lo; i < hi; i++) {
            A[i] = 0.0;
            B[i] = (double)(i & 0xFF) * 0.01;
            C[i] = (double)((i + 7) & 0xFF) * 0.01;
        }
    }
    return da;
}

static void free_domain_arrays(DomainArrays *da, int P_N) {
    for (int d = 0; d < P_N; d++) {
        free(da[d].A); free(da[d].B); free(da[d].C);
    }
    free(da);
}

/* ─── Triad kernel ─── */

static void triad_phase(DomainArrays *da, int P_N, int phase,
                        size_t N, int nth) {
    /*
     * Domain d's threads operate on arrays owned by domain (d + phase) % P_N.
     * Each thread computes its static chunk of the target domain's arrays.
     */
    #pragma omp parallel num_threads(nth)
    {
        int tid = omp_get_thread_num();
        int my_domain = thread_to_domain(tid, nth, P_N);
        int target_domain = (my_domain + phase) % P_N;

        int tpd = (nth + P_N - 1) / P_N;
        int local_tid = tid - my_domain * tpd;
        if (local_tid < 0) local_tid = 0;

        size_t chunk = (N + tpd - 1) / tpd;
        size_t lo = (size_t)local_tid * chunk;
        size_t hi = lo + chunk;
        if (hi > N) hi = N;

        double *__restrict__ A = da[target_domain].A;
        const double *__restrict__ B = da[target_domain].B;
        const double *__restrict__ C = da[target_domain].C;

        for (size_t i = lo; i < hi; i++)
            A[i] = B[i] + SCALAR * C[i];
    }
}

/* ─── Main ─── */

int main(int argc, char **argv) {
    size_t N     = (argc > 1) ? (size_t)atol(argv[1]) : DEFAULT_N;
    int    iters = (argc > 2) ? atoi(argv[2]) : DEFAULT_ITERS;
    const char *csv = (argc > 3) ? argv[3] : DEFAULT_CSV;

    int nth = omp_get_max_threads();
    int P_N = detect_numa_domains();

    printf("numa_rotation_bench\n");
    printf("  N = %zu (%.1f MiB/array), iters = %d\n",
           N, (double)N * sizeof(double) / (1 << 20), iters);
    printf("  threads = %d, NUMA domains = %d, threads/domain = %d\n\n",
           nth, P_N, (nth + P_N - 1) / P_N);

    DomainArrays *da = alloc_domain_arrays(P_N, N, nth);
    size_t bytes_per_iter = 3 * N * sizeof(double);  /* 2 reads + 1 write */

    FILE *fout = fopen(csv, "w");
    fprintf(fout, "phase,numa_shift,iter,time_s,bw_gbs\n");

    printf("  %-8s %-12s %-12s %-12s %-10s\n",
           "Phase", "NUMA shift", "Med BW GB/s", "Min BW", "Max BW");
    printf("  %s\n", "------------------------------------------------------");

    for (int phase = 0; phase < P_N; phase++) {
        /* Warmup */
        for (int w = 0; w < WARMUP; w++)
            triad_phase(da, P_N, phase, N, nth);

        double *times = (double *)malloc(iters * sizeof(double));

        for (int it = 0; it < iters; it++) {
            double t0 = omp_get_wtime();
            triad_phase(da, P_N, phase, N, nth);
            double t1 = omp_get_wtime();
            times[it] = t1 - t0;
            double bw = (double)bytes_per_iter / times[it] / 1e9;
            fprintf(fout, "%d,%d,%d,%.9f,%.4f\n", phase, phase, it, times[it], bw);
        }

        /* Sort for median/min/max */
        for (int i = 0; i < iters - 1; i++)
            for (int j = i + 1; j < iters; j++)
                if (times[j] < times[i]) {
                    double tmp = times[i]; times[i] = times[j]; times[j] = tmp;
                }

        double med_bw = (double)bytes_per_iter / times[iters / 2] / 1e9;
        double max_bw = (double)bytes_per_iter / times[0] / 1e9;
        double min_bw = (double)bytes_per_iter / times[iters - 1] / 1e9;

        const char *label;
        if (phase == 0) label = "local";
        else if (phase == P_N / 2) label = "max-remote";
        else label = "cross-NUMA";

        printf("  %-8d %-12s %-12.2f %-12.2f %-10.2f\n",
               phase, label, med_bw, min_bw, max_bw);

        free(times);
    }

    fclose(fout);
    free_domain_arrays(da, P_N);

    printf("\nWrote: %s\n", csv);

    /* Print gamma estimates */
    printf("\n=== Derived γ (NUMA penalty) ===\n");
    printf("  γ = BW_local / BW_remote for each phase shift.\n");
    printf("  Use phase 0 as baseline (local).\n");
    printf("  Rerun with: OMP_PLACES=cores OMP_PROC_BIND=close\n");

    return 0;
}