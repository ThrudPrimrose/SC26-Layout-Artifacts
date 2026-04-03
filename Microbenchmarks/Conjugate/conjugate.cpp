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

/*  Total doubles across all K fields is constant.
 *  TOTAL = K * N_base   =>   same memory footprint for every K.
 *  OOP bandwidth = 2 * TOTAL * 8  (read all + write all).            */
constexpr int64_t TOTAL = 1LL << 28;          // ~2 GB per side
constexpr int64_t RUNS  = 100;
constexpr int64_t CPU_MAX_VL = 512;

/* ═══ Generic layout structs ═══ */

template<int K>
struct AoS { double f[K]; };

template<int K, int VL>
struct AoSoA { double f[K][VL]; };

/* ═══ NUMA helpers (same as original) ═══ */

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

constexpr int64_t FLUSH_N = 1 << 26;
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
}
static void flush_free() { numa_free(flush_buf, FLUSH_N * sizeof(double)); }

/* ═══ Init helpers ═══ */

static void init_and_bind(double *buf, int64_t n) {
    bind_and_touch(buf, n * sizeof(double));
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        buf[i] = (double)(i % 997) * 0.001;
}

template<int K>
static void init_aos(AoS<K> *buf, int64_t n) {
    bind_and_touch(buf, n * sizeof(AoS<K>));
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        for (int k = 0; k < K; k++)
            buf[i].f[k] = (double)((i * K + k) % 997) * 0.001;
}

template<int K, int VL>
static void init_aosoa(AoSoA<K,VL> *buf, int64_t n_base) {
    int64_t nblks = n_base / VL;
    bind_and_touch(buf, nblks * sizeof(AoSoA<K,VL>));
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < nblks; b++)
        for (int k = 0; k < K; k++)
            for (int l = 0; l < VL; l++)
                buf[b].f[k][l] = (double)(((b * VL + l) * K + k) % 997) * 0.001;
}

/* ═══ Kernels (out-of-place): negate all K fields ═══ */

template<int K>
static void kern_aos(const AoS<K> *__restrict__ in,
                     AoS<K> *__restrict__ out, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        for (int k = 0; k < K; k++)
            out[i].f[k] = -in[i].f[k];
}

template<int K>
static void kern_soa(double *const *__restrict__ in,
                     double *const *__restrict__ out, int64_t n) {
    /* Every thread touches all 2K streams → stresses TLB */
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        for (int k = 0; k < K; k++)
            out[k][i] = -in[k][i];
}

template<int K, int VL>
static void kern_aosoa(const AoSoA<K,VL> *__restrict__ in,
                       AoSoA<K,VL> *__restrict__ out, int64_t n_base) {
    int64_t nblks = n_base / VL;
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < nblks; b++)
        for (int k = 0; k < K; k++)
            #pragma omp simd
            for (int l = 0; l < VL; l++)
                out[b].f[k][l] = -in[b].f[k][l];
}

/* ═══ Bench driver ═══ */

static FILE *csv;

/* BW = 2 * K * n_base * 8   (read K fields + write K fields) */
static double bw_oop(int K, int64_t n_base, double ms) {
    return 2.0 * K * n_base * sizeof(double) / (ms * 1e6);
}

#define CPU_BENCH(K_val, n_base, label, call) do { \
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
                K_val, label, r, times[r], bw_oop(K_val, n_base, times[r])); \
        sum += times[r]; \
    } \
    double avg = sum / RUNS; \
    printf("  %-14s %8.4f ms  %7.1f GB/s\n", label, avg, bw_oop(K_val, n_base, avg)); \
} while (0)

template<int K>
static void bench_aos(int64_t n) {
    size_t bytes = n * sizeof(AoS<K>);
    auto *in  = (AoS<K> *)numa_alloc(bytes); init_aos<K>(in, n);
    auto *out = (AoS<K> *)numa_alloc(bytes); init_aos<K>(out, n);
    CPU_BENCH(K, n, "AoS", (kern_aos<K>(in, out, n)));
    numa_free(in, bytes); numa_free(out, bytes);
}

template<int K>
static void bench_soa(int64_t n) {
    size_t bytes = n * sizeof(double);
    double *in_arr[K], *out_arr[K];
    for (int k = 0; k < K; k++) {
        in_arr[k]  = (double *)numa_alloc(bytes); init_and_bind(in_arr[k], n);
        out_arr[k] = (double *)numa_alloc(bytes); init_and_bind(out_arr[k], n);
    }
    double *const *ip = in_arr, *const *op = out_arr;
    CPU_BENCH(K, n, "SoA", (kern_soa<K>(ip, op, n)));
    for (int k = 0; k < K; k++) {
        numa_free(in_arr[k], bytes); numa_free(out_arr[k], bytes);
    }
}

template<int K, int VL>
static void bench_aosoa(int64_t n, const char *label) {
    int64_t nblks = n / VL;
    size_t bytes = nblks * sizeof(AoSoA<K,VL>);
    auto *in  = (AoSoA<K,VL> *)numa_alloc(bytes); init_aosoa<K,VL>(in, n);
    auto *out = (AoSoA<K,VL> *)numa_alloc(bytes); init_aosoa<K,VL>(out, n);
    CPU_BENCH(K, n, label, (kern_aosoa<K,VL>(in, out, n)));
    numa_free(in, bytes); numa_free(out, bytes);
}

template<int K>
static void run_all(int64_t n) {
    printf("\n── K=%d  N_base=%lld  total=%.1f GB (per side) ──\n",
           K, (long long)n, (double)K * n * 8 / 1e9);
    bench_aos<K>(n);
    bench_soa<K>(n);
    bench_aosoa<K,   2>(n, "AoSoA-2");
    bench_aosoa<K,   4>(n, "AoSoA-4");
    bench_aosoa<K,   8>(n, "AoSoA-8");
    bench_aosoa<K,  16>(n, "AoSoA-16");
    bench_aosoa<K,  32>(n, "AoSoA-32");
    bench_aosoa<K,  64>(n, "AoSoA-64");
    bench_aosoa<K, 128>(n, "AoSoA-128");
    bench_aosoa<K, 256>(n, "AoSoA-256");
    bench_aosoa<K, 512>(n, "AoSoA-512");
}

static int64_t n_base(int K) {
    return (TOTAL / K / CPU_MAX_VL) * CPU_MAX_VL;
}

int main() {
    PAGE_SZ = sysconf(_SC_PAGESIZE);
    flush_init();

    csv = fopen("results_cpu_oop.csv", "w");
    fprintf(csv, "K,layout,run,ms,gbps\n");

    printf("conj OOP: TOTAL=%lldM  runs=%d  threads=%d  page=%ld\n",
           (long long)(TOTAL >> 20), (int)RUNS, omp_get_max_threads(), PAGE_SZ);
    printf("BW = 2·K·N_base·8  (read+write all K fields)\n");

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

    run_all< 3>(n_base( 3));
    run_all< 6>(n_base( 6));
    run_all< 9>(n_base( 9));
    run_all<12>(n_base(12));
    run_all<15>(n_base(15));

    fclose(csv);
    flush_free();
    printf("\nwrote results_cpu_oop.csv\n");
}