#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <numaif.h>
#include <unistd.h>

constexpr int64_t N    = 1 << 30;   // 1B doubles = 8 GB per array
constexpr int64_t RUNS = 10;
constexpr double  S    = 3.0;

static long PAGE_SZ;

static int get_numa_node() {
    unsigned cpu, node;
    syscall(__NR_getcpu, &cpu, &node, nullptr);
    return (int)node;
}

static void print_mems_allowed() {
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return;
    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "Mems_allowed:", 13) == 0 ||
            strncmp(line, "Mems_allowed_list:", 18) == 0 ||
            strncmp(line, "Cpus_allowed_list:", 18) == 0)
            printf("  %s", line);
    }
    fclose(f);
}

static void *numa_alloc(size_t bytes) {
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_NOHUGEPAGE);
    return p;
}
static void numa_free(void *p, size_t bytes) { munmap(p, bytes); }

/* Per-thread mbind(MPOL_BIND) + first-touch.
   Each OpenMP thread binds its schedule(static) chunk to its local NUMA node,
   then writes to fault pages locally. */
static double *alloc_ft(int64_t n, double val) {
    size_t bytes = n * sizeof(double);
    double *p = (double *)numa_alloc(bytes);

    #pragma omp parallel
    {
        int tid  = omp_get_thread_num();
        int nthr = omp_get_num_threads();
        /* Replicate schedule(static) partitioning */
        int64_t chunk = (n + nthr - 1) / nthr;
        int64_t lo = (int64_t)tid * chunk;
        int64_t hi = std::min(lo + chunk, n);
        if (lo < n && lo < hi) {
            /* Page-align the byte range */
            uintptr_t beg = (uintptr_t)p + lo * sizeof(double);
            uintptr_t end = (uintptr_t)p + hi * sizeof(double);
            uintptr_t pbeg = beg & ~(uintptr_t)(PAGE_SZ - 1);
            uintptr_t pend = (end + PAGE_SZ - 1) & ~(uintptr_t)(PAGE_SZ - 1);

            int node = get_numa_node();
            unsigned long mask = 1UL << node;
            mbind((void *)pbeg, pend - pbeg, MPOL_BIND, &mask, 37, 0);

            for (int64_t i = lo; i < hi; i++) p[i] = val;
        }
    }
    return p;
}

/* ═══ Cache flush ═══ */
constexpr int64_t FLUSH_N = 1 << 26;  // 512 MB
static double *flush_buf = nullptr;

static void flush_init() {
    constexpr size_t bytes = FLUSH_N * sizeof(double);
    flush_buf = (double *)numa_alloc(bytes);
    #pragma omp parallel
    {
        int tid  = omp_get_thread_num();
        int nthr = omp_get_num_threads();
        int64_t chunk = (FLUSH_N + nthr - 1) / nthr;
        int64_t lo = (int64_t)tid * chunk;
        int64_t hi = std::min(lo + chunk, FLUSH_N);
        if (lo < FLUSH_N && lo < hi) {
            uintptr_t beg = (uintptr_t)flush_buf + lo * sizeof(double);
            uintptr_t end = (uintptr_t)flush_buf + hi * sizeof(double);
            uintptr_t pbeg = beg & ~(uintptr_t)(PAGE_SZ - 1);
            uintptr_t pend = (end + PAGE_SZ - 1) & ~(uintptr_t)(PAGE_SZ - 1);
            int node = get_numa_node();
            unsigned long mask = 1UL << node;
            mbind((void *)pbeg, pend - pbeg, MPOL_BIND, &mask, 37, 0);
            for (int64_t i = lo; i < hi; i++) flush_buf[i] = 0.0;
        }
    }
}

static void flush_caches() {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] += 1.0;
}

static void flush_free() {
    numa_free(flush_buf, FLUSH_N * sizeof(double));
}

/* ═══ Store: a = b + s*c  →  3 elements/i ═══ */

static void triad_s1(double *__restrict__ a1, const double *__restrict__ b1, const double *__restrict__ c1,
                     double s, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) a1[i] = b1[i] + s * c1[i];
}

static void triad_s2(double *__restrict__ a1, const double *__restrict__ b1, const double *__restrict__ c1,
                     double *__restrict__ a2, const double *__restrict__ b2, const double *__restrict__ c2,
                     double s, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        a1[i] = b1[i] + s * c1[i];
        a2[i] = b2[i] + s * c2[i];
    }
}

static void triad_s3(double *__restrict__ a1, const double *__restrict__ b1, const double *__restrict__ c1,
                     double *__restrict__ a2, const double *__restrict__ b2, const double *__restrict__ c2,
                     double *__restrict__ a3, const double *__restrict__ b3, const double *__restrict__ c3,
                     double s, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        a1[i] = b1[i] + s * c1[i];
        a2[i] = b2[i] + s * c2[i];
        a3[i] = b3[i] + s * c3[i];
    }
}

static void triad_s5(double **__restrict__ a, const double **__restrict__ b, const double **__restrict__ c,
                     double s, int64_t n) {
    double *__restrict__ a0=a[0],*__restrict__ a1=a[1],*__restrict__ a2=a[2],*__restrict__ a3=a[3],*__restrict__ a4=a[4];
    const double *__restrict__ b0=b[0],*__restrict__ b1=b[1],*__restrict__ b2=b[2],*__restrict__ b3=b[3],*__restrict__ b4=b[4];
    const double *__restrict__ c0=c[0],*__restrict__ c1=c[1],*__restrict__ c2=c[2],*__restrict__ c3=c[3],*__restrict__ c4=c[4];
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        a0[i] = b0[i] + s * c0[i];
        a1[i] = b1[i] + s * c1[i];
        a2[i] = b2[i] + s * c2[i];
        a3[i] = b3[i] + s * c3[i];
        a4[i] = b4[i] + s * c4[i];
    }
}

static void triad_s7(double **__restrict__ a, const double **__restrict__ b, const double **__restrict__ c,
                     double s, int64_t n) {
    double *__restrict__ a0=a[0],*__restrict__ a1=a[1],*__restrict__ a2=a[2],*__restrict__ a3=a[3],
           *__restrict__ a4=a[4],*__restrict__ a5=a[5],*__restrict__ a6=a[6];
    const double *__restrict__ b0=b[0],*__restrict__ b1=b[1],*__restrict__ b2=b[2],*__restrict__ b3=b[3],
                 *__restrict__ b4=b[4],*__restrict__ b5=b[5],*__restrict__ b6=b[6];
    const double *__restrict__ c0=c[0],*__restrict__ c1=c[1],*__restrict__ c2=c[2],*__restrict__ c3=c[3],
                 *__restrict__ c4=c[4],*__restrict__ c5=c[5],*__restrict__ c6=c[6];
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        a0[i] = b0[i] + s * c0[i];
        a1[i] = b1[i] + s * c1[i];
        a2[i] = b2[i] + s * c2[i];
        a3[i] = b3[i] + s * c3[i];
        a4[i] = b4[i] + s * c4[i];
        a5[i] = b5[i] + s * c5[i];
        a6[i] = b6[i] + s * c6[i];
    }
}

/* ═══ Accum: a += b + s*c  →  4 elements/i ═══ */

static void triad_a1(double *__restrict__ a1, const double *__restrict__ b1, const double *__restrict__ c1,
                     double s, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) a1[i] += b1[i] + s * c1[i];
}

static void triad_a2(double *__restrict__ a1, const double *__restrict__ b1, const double *__restrict__ c1,
                     double *__restrict__ a2, const double *__restrict__ b2, const double *__restrict__ c2,
                     double s, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        a1[i] += b1[i] + s * c1[i];
        a2[i] += b2[i] + s * c2[i];
    }
}

static void triad_a3(double *__restrict__ a1, const double *__restrict__ b1, const double *__restrict__ c1,
                     double *__restrict__ a2, const double *__restrict__ b2, const double *__restrict__ c2,
                     double *__restrict__ a3, const double *__restrict__ b3, const double *__restrict__ c3,
                     double s, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        a1[i] += b1[i] + s * c1[i];
        a2[i] += b2[i] + s * c2[i];
        a3[i] += b3[i] + s * c3[i];
    }
}

static void triad_a5(double **__restrict__ a, const double **__restrict__ b, const double **__restrict__ c,
                     double s, int64_t n) {
    double *__restrict__ a0=a[0],*__restrict__ a1=a[1],*__restrict__ a2=a[2],*__restrict__ a3=a[3],*__restrict__ a4=a[4];
    const double *__restrict__ b0=b[0],*__restrict__ b1=b[1],*__restrict__ b2=b[2],*__restrict__ b3=b[3],*__restrict__ b4=b[4];
    const double *__restrict__ c0=c[0],*__restrict__ c1=c[1],*__restrict__ c2=c[2],*__restrict__ c3=c[3],*__restrict__ c4=c[4];
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        a0[i] += b0[i] + s * c0[i];
        a1[i] += b1[i] + s * c1[i];
        a2[i] += b2[i] + s * c2[i];
        a3[i] += b3[i] + s * c3[i];
        a4[i] += b4[i] + s * c4[i];
    }
}

static void triad_a7(double **__restrict__ a, const double **__restrict__ b, const double **__restrict__ c,
                     double s, int64_t n) {
    double *__restrict__ a0=a[0],*__restrict__ a1=a[1],*__restrict__ a2=a[2],*__restrict__ a3=a[3],
           *__restrict__ a4=a[4],*__restrict__ a5=a[5],*__restrict__ a6=a[6];
    const double *__restrict__ b0=b[0],*__restrict__ b1=b[1],*__restrict__ b2=b[2],*__restrict__ b3=b[3],
                 *__restrict__ b4=b[4],*__restrict__ b5=b[5],*__restrict__ b6=b[6];
    const double *__restrict__ c0=c[0],*__restrict__ c1=c[1],*__restrict__ c2=c[2],*__restrict__ c3=c[3],
                 *__restrict__ c4=c[4],*__restrict__ c5=c[5],*__restrict__ c6=c[6];
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        a0[i] += b0[i] + s * c0[i];
        a1[i] += b1[i] + s * c1[i];
        a2[i] += b2[i] + s * c2[i];
        a3[i] += b3[i] + s * c3[i];
        a4[i] += b4[i] + s * c4[i];
        a5[i] += b5[i] + s * c5[i];
        a6[i] += b6[i] + s * c6[i];
    }
}

/* ═══ Reporting ═══ */

static FILE *csv;
static double bw_store(int nt, double ms) { return (double)nt * 3.0 * N * sizeof(double) / (ms * 1e6); }
static double bw_accum(int nt, double ms) { return (double)nt * 4.0 * N * sizeof(double) / (ms * 1e6); }

#define BENCH(label, nt, call, bw_fn) do { \
    for (int w = 0; w < 5; w++) { flush_caches(); call; } \
    flush_caches(); \
    double times[RUNS]; \
    for (int r = 0; r < RUNS; r++) { \
        double t0 = omp_get_wtime(); \
        call; \
        times[r] = (omp_get_wtime() - t0) * 1e3; \
    } \
    double sum = 0; \
    for (int r = 0; r < RUNS; r++) { \
        fprintf(csv, "%s,%d,%.6f,%.2f\n", label, r, times[r], bw_fn(nt, times[r])); \
        sum += times[r]; \
    } \
    double avg = sum / RUNS; \
    printf("  %-22s  %8.4f ms  %7.1f GB/s\n", label, avg, bw_fn(nt, avg)); \
} while (0)

struct TriadSet { double *a, *b, *c; };
static TriadSet make_set(double av, double bv, double cv) {
    return { alloc_ft(N, av), alloc_ft(N, bv), alloc_ft(N, cv) };
}
static void free_set(TriadSet &s) {
    constexpr size_t bytes = N * sizeof(double);
    numa_free(s.a, bytes); numa_free(s.b, bytes); numa_free(s.c, bytes);
}

int main() {
    PAGE_SZ = sysconf(_SC_PAGESIZE);
    flush_init();
    constexpr size_t bytes = N * sizeof(double);
    csv = fopen("results_triad.csv", "w");
    fprintf(csv, "label,run,ms,gbps\n");

    printf("triad: N=%lldM (%.0f MB/array)  runs=%d  threads=%d  page=%ld\n",
           (long long)(N>>20), bytes/1e6, (int)RUNS, omp_get_max_threads(), PAGE_SZ);
    printf("store: 3N×8/triad  accum: 4N×8/triad  |  cache flush before each run\n");
    printf("NUMA: per-thread mbind(MPOL_BIND) + first-touch\n");
    print_mems_allowed();

    /* Verify thread spread */
    {
        int nodes[4] = {};
        #pragma omp parallel
        {
            int n = get_numa_node();
            if (n >= 0 && n < 4) {
                #pragma omp atomic
                nodes[n]++;
            }
        }
        printf("  thread spread: N0=%d N1=%d N2=%d N3=%d\n\n", nodes[0], nodes[1], nodes[2], nodes[3]);
    }

    { TriadSet s1 = make_set(0,1,2);
      BENCH("1×store  (3 arr)",  1, triad_s1(s1.a,s1.b,s1.c,S,N), bw_store);
      BENCH("1×accum  (3 arr)",  1, triad_a1(s1.a,s1.b,s1.c,S,N), bw_accum);
      free_set(s1); }

    { TriadSet s1=make_set(0,1,2), s2=make_set(0,3,4);
      BENCH("2×store  (6 arr)",  2, triad_s2(s1.a,s1.b,s1.c,s2.a,s2.b,s2.c,S,N), bw_store);
      BENCH("2×accum  (6 arr)",  2, triad_a2(s1.a,s1.b,s1.c,s2.a,s2.b,s2.c,S,N), bw_accum);
      free_set(s1); free_set(s2); }

    { TriadSet s1=make_set(0,1,2), s2=make_set(0,3,4), s3=make_set(0,5,6);
      BENCH("3×store  (9 arr)",  3, triad_s3(s1.a,s1.b,s1.c,s2.a,s2.b,s2.c,s3.a,s3.b,s3.c,S,N), bw_store);
      BENCH("3×accum  (9 arr)",  3, triad_a3(s1.a,s1.b,s1.c,s2.a,s2.b,s2.c,s3.a,s3.b,s3.c,S,N), bw_accum);
      free_set(s1); free_set(s2); free_set(s3); }

    { TriadSet sets[5]; double *ap[5],*bp[5],*cp[5];
      for (int k=0;k<5;k++) { sets[k]=make_set(0,(double)(2*k+1),(double)(2*k+2)); ap[k]=sets[k].a; bp[k]=sets[k].b; cp[k]=sets[k].c; }
      BENCH("5×store  (15 arr)", 5, triad_s5(ap,(const double**)bp,(const double**)cp,S,N), bw_store);
      BENCH("5×accum  (15 arr)", 5, triad_a5(ap,(const double**)bp,(const double**)cp,S,N), bw_accum);
      for (int k=0;k<5;k++) free_set(sets[k]); }

    { TriadSet sets[7]; double *ap[7],*bp[7],*cp[7];
      for (int k=0;k<7;k++) { sets[k]=make_set(0,(double)(2*k+1),(double)(2*k+2)); ap[k]=sets[k].a; bp[k]=sets[k].b; cp[k]=sets[k].c; }
      BENCH("7×store  (21 arr)", 7, triad_s7(ap,(const double**)bp,(const double**)cp,S,N), bw_store);
      BENCH("7×accum  (21 arr)", 7, triad_a7(ap,(const double**)bp,(const double**)cp,S,N), bw_accum);
      for (int k=0;k<7;k++) free_set(sets[k]); }

    fclose(csv);
    flush_free();
    printf("\nwrote results_triad.csv\n");
}