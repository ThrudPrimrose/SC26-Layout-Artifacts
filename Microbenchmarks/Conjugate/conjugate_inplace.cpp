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

/*  Conjugate P complex arrays in-place:
 *      buf[i].im_p = -buf[i].im_p      (negate im only)
 *
 *  Total doubles = 2·P·N_base = const across P.
 *
 *  Effective BW (useful bytes):
 *    SoA:   2·P·N_base·8   (only P im arrays read+written)
 *    AoS:   read full struct (2P fields) + write back ⇒ hw moves
 *           2·(2P)·N_base·8 but only P fields are useful.
 *  We report useful BW = 2·P·N_base·8 uniformly so all layouts
 *  are compared on the same "work done" basis.                       */

constexpr int64_t TOTAL_DOUBLES = 1LL << 28;
constexpr int64_t RUNS  = 100;
constexpr int64_t MAX_VL = 512;

template<int P>
struct AoS {
    struct { double re, im; } c[P];
};

template<int P, int VL>
struct AoSoA {
    struct { double re[VL], im[VL]; } c[P];
};

/* ═══ NUMA helpers ═══ */

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
            uintptr_t pend = std::min((uintptr_t)base + hi * PAGE_SZ,
                                      (uintptr_t)base + total_bytes);
            int node = get_numa_node();
            unsigned long mask = 1UL << node;
            mbind((void *)pbeg, pend - pbeg, MPOL_BIND, &mask, 64, 0);
            for (uintptr_t a = pbeg; a < pend; a += PAGE_SZ)
                *(volatile char *)a = 0;
        }
    }
}

/* ═══ Cache flush ═══ */

constexpr int64_t FLUSH_N = 1 << 27;   /* 1 GB — 4× per-CCX L3 on MI300A */
static double *flush_buf;

static void flush_init() {
    size_t b = FLUSH_N * sizeof(double);
    flush_buf = (double *)numa_alloc(b);
    bind_and_touch(flush_buf, b);
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] = 0.0;
}
static void flush_caches() {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] += 1.0;
    /* Full fence: ensure all stores (including write-backs) complete
       before we start timing the next kernel invocation.             */
    #pragma omp parallel
    { asm volatile("mfence" ::: "memory"); }
}
static void flush_free() { numa_free(flush_buf, FLUSH_N * sizeof(double)); }

/* ═══ Init helpers ═══ */

static void init_and_bind(double *buf, int64_t n) {
    bind_and_touch(buf, n * sizeof(double));
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        buf[i] = (double)(i % 997) * 0.001;
}

template<int P>
static void init_aos(AoS<P> *buf, int64_t n) {
    bind_and_touch(buf, n * sizeof(AoS<P>));
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        for (int p = 0; p < P; p++) {
            buf[i].c[p].re = (double)((i * P + p) % 997) * 0.001;
            buf[i].c[p].im = (double)((i * P + p) % 991) * 0.001;
        }
}

template<int P, int VL>
static void init_aosoa(AoSoA<P,VL> *buf, int64_t n_base) {
    int64_t nblks = n_base / VL;
    bind_and_touch(buf, nblks * sizeof(AoSoA<P,VL>));
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < nblks; b++)
        for (int p = 0; p < P; p++)
            for (int l = 0; l < VL; l++) {
                int64_t idx = b * VL + l;
                buf[b].c[p].re[l] = (double)((idx * P + p) % 997) * 0.001;
                buf[b].c[p].im[l] = (double)((idx * P + p) % 991) * 0.001;
            }
}

/* ═══ Kernels (in-place): negate im only ═══ */

template<int P>
static void kern_aos(AoS<P> *__restrict__ buf, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        for (int p = 0; p < P; p++)
            buf[i].c[p].im = -buf[i].c[p].im;
}

/* SoA: only the P im arrays are touched */
template<int P>
static void kern_soa(double *const *__restrict__ im, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        for (int p = 0; p < P; p++)
            im[p][i] = -im[p][i];
}

template<int P, int VL>
static void kern_aosoa(AoSoA<P,VL> *__restrict__ buf, int64_t n_base) {
    int64_t nblks = n_base / VL;
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < nblks; b++)
        for (int p = 0; p < P; p++)
            #pragma omp simd
            for (int l = 0; l < VL; l++)
                buf[b].c[p].im[l] = -buf[b].c[p].im[l];
}

/* ═══ Bench driver ═══ */

static FILE *csv;

/* Useful BW: P im fields read + P im fields written = 2·P·N_base·8 */
static double bw_ip(int P, int64_t n_base, double ms) {
    return 2.0 * P * n_base * sizeof(double) / (ms * 1e6);
}

#define CPU_BENCH(P_val, n_base, label, call) do { \
    for (int w = 0; w < 5; w++) { flush_caches(); call; } \
    double times[RUNS]; \
    for (int r = 0; r < RUNS; r++) { \
        flush_caches(); \
        double t0 = omp_get_wtime(); \
        call; \
        times[r] = (omp_get_wtime() - t0) * 1e3; \
    } \
    double sum = 0; \
    for (int r = 0; r < RUNS; r++) { \
        fprintf(csv, "%d,%s,%d,%.6f,%.2f\n", \
                P_val, label, r, times[r], bw_ip(P_val, n_base, times[r])); \
        sum += times[r]; \
    } \
    double avg = sum / RUNS; \
    printf("  %-14s %8.4f ms  %7.1f GB/s\n", label, avg, bw_ip(P_val, n_base, avg)); \
} while (0)

template<int P>
static void bench_aos(int64_t n) {
    size_t bytes = n * sizeof(AoS<P>);
    auto *buf = (AoS<P> *)numa_alloc(bytes); init_aos<P>(buf, n);
    CPU_BENCH(P, n, "AoS", (kern_aos<P>(buf, n)));
    numa_free(buf, bytes);
}

template<int P>
static void bench_soa(int64_t n) {
    size_t bytes = n * sizeof(double);
    double *re[P], *im[P];
    for (int p = 0; p < P; p++) {
        re[p] = (double *)numa_alloc(bytes); init_and_bind(re[p], n);
        im[p] = (double *)numa_alloc(bytes); init_and_bind(im[p], n);
    }
    CPU_BENCH(P, n, "SoA", (kern_soa<P>(im, n)));
    for (int p = 0; p < P; p++) {
        numa_free(re[p], bytes); numa_free(im[p], bytes);
    }
}

template<int P, int VL>
static void bench_aosoa(int64_t n, const char *label) {
    int64_t nblks = n / VL;
    size_t bytes = nblks * sizeof(AoSoA<P,VL>);
    auto *buf = (AoSoA<P,VL> *)numa_alloc(bytes); init_aosoa<P,VL>(buf, n);
    CPU_BENCH(P, n, label, (kern_aosoa<P,VL>(buf, n)));
    numa_free(buf, bytes);
}

template<int P>
static void run_all() {
    int64_t n = (TOTAL_DOUBLES / (2 * P) / MAX_VL) * MAX_VL;
    printf("\n── P=%d complex pairs  (%d im streams IP)  N_base=%lld  "
           "total=%.1f GB ──\n",
           P, P, (long long)n, 2.0 * P * n * 8 / 1e9);

    bench_aos<P>(n);
    bench_soa<P>(n);
    bench_aosoa<P,   2>(n, "AoSoA-2");
    bench_aosoa<P,   4>(n, "AoSoA-4");
    bench_aosoa<P,   8>(n, "AoSoA-8");
    bench_aosoa<P,  16>(n, "AoSoA-16");
    bench_aosoa<P,  32>(n, "AoSoA-32");
    bench_aosoa<P,  64>(n, "AoSoA-64");
    bench_aosoa<P, 128>(n, "AoSoA-128");
    bench_aosoa<P, 256>(n, "AoSoA-256");
    bench_aosoa<P, 512>(n, "AoSoA-512");
}

int main() {
    PAGE_SZ = sysconf(_SC_PAGESIZE);
    flush_init();

    csv = fopen("results_cpu_inplace.csv", "w");
    fprintf(csv, "P,layout,run,ms,gbps\n");

    printf("conjugate IN-PLACE (CPU): TOTAL_DOUBLES=%lldM  runs=%d  threads=%d\n",
           (long long)(TOTAL_DOUBLES >> 20), (int)RUNS, omp_get_max_threads());
    printf("BW = 2·P·N_base·8  (read P im + write P im)\n");

    {
        int nodes[8] = {};
        #pragma omp parallel
        { int nd = get_numa_node(); if (nd >= 0 && nd < 8) {
            #pragma omp atomic
            nodes[nd]++;
        }}
        printf("  thread spread:");
        for (int i = 0; i < 8; i++) if (nodes[i]) printf(" N%d=%d", i, nodes[i]);
        printf("\n");
    }

    run_all< 3>();
    run_all< 6>();
    run_all< 9>();
    run_all<12>();
    run_all<15>();
    run_all<18>();
    run_all<21>();

    fclose(csv);
    flush_free();
    printf("\nwrote results_cpu_inplace.csv\n");
}