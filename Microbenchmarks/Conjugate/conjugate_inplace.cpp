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
#include <type_traits>

constexpr int64_t N    = 1 << 28;
constexpr int64_t RUNS = 100;

struct C2 { double re, im; };
template<int VL> struct C2V { double re[VL], im[VL]; };

static_assert(sizeof(C2) == 2 * sizeof(double), "");
#define CHECK_C2V(VL) static_assert(sizeof(C2V<VL>) == 2*(VL)*sizeof(double), "");
CHECK_C2V(2) CHECK_C2V(4) CHECK_C2V(8) CHECK_C2V(16) CHECK_C2V(32) CHECK_C2V(64) CHECK_C2V(128) CHECK_C2V(256) CHECK_C2V(512)
#undef CHECK_C2V

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
        if (strncmp(line, "Mems_allowed_list:", 18) == 0 ||
            strncmp(line, "Cpus_allowed_list:", 18) == 0)
            printf("  %s", line);
    }
    fclose(f);
}

/* ═══ NUMA-aware allocation: mmap + per-thread mbind + first-touch ═══ */

static void *numa_alloc(size_t bytes) {
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_NOHUGEPAGE);
    return p;
}

static void numa_free(void *p, size_t bytes) { munmap(p, bytes); }

/* Bind and first-touch: each thread mbinds its schedule(static) chunk to its NUMA node */
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
            unsigned long mask = 1UL << node;
            mbind((void *)pbeg, pend - pbeg, MPOL_BIND, &mask, 64, 0);

            volatile char *cp = (volatile char *)pbeg;
            for (uintptr_t addr = pbeg; addr < pend; addr += PAGE_SZ)
                *(volatile char *)addr = 0;
        }
    }
}

/* ═══ First-touch with data init (matching compute schedule) ═══ */

static void ft_aos(C2 *buf, int64_t n) {
    bind_and_touch(buf, n * sizeof(C2));
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        buf[i].re = (double)(i % 997) * 0.001;
        buf[i].im = (double)(i % 991) * 0.001;
    }
}

static void ft_soa(double *r, double *im, int64_t n) {
    bind_and_touch(r,  n * sizeof(double));
    bind_and_touch(im, n * sizeof(double));
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        r[i]  = (double)(i % 997) * 0.001;
        im[i] = (double)(i % 991) * 0.001;
    }
}

template<int VL>
static void ft_aosoa(C2V<VL> *buf, int64_t n) {
    bind_and_touch(buf, (n / VL) * sizeof(C2V<VL>));
    int64_t nblks = n / VL;
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < nblks; b++)
        for (int l = 0; l < VL; l++) {
            buf[b].re[l] = (double)((b * VL + l) % 997) * 0.001;
            buf[b].im[l] = (double)((b * VL + l) % 991) * 0.001;
        }
}

/* ═══ Cache flush ═══ */

constexpr int64_t FLUSH_N = 1 << 26;
static double *flush_buf = nullptr;

static void flush_init() {
    constexpr size_t bytes = FLUSH_N * sizeof(double);
    flush_buf = (double *)numa_alloc(bytes);
    bind_and_touch(flush_buf, bytes);
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] = 0.0;
}

static void flush_caches() {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] += 1.0;
}

static void flush_free() { numa_free(flush_buf, FLUSH_N * sizeof(double)); }

/* ═══ Kernels (in-place): im = -im ═══ */

static void c_aos(C2 *__restrict__ in, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) in[i].im = -in[i].im;
}

static void c_soa(double *__restrict__ ii, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) ii[i] = -ii[i];
}

template<int VL>
static void c_aosoa(C2V<VL> *__restrict__ in, int64_t n) {
    int64_t nblks = n / VL;
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < nblks; b++)
        #pragma omp simd
        for (int l = 0; l < VL; l++)
            in[b].im[l] = -in[b].im[l];
}

/* ═══ BW = N elements read + N elements written = 2N×8 ═══ */

static FILE *csv;
static double bw(double ms) { return 2.0 * N * sizeof(double) / (ms * 1e6); }

#define CPU_BENCH(label, call) do { \
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
        fprintf(csv, "%s,%d,%.6f,%.2f\n", label, r, times[r], bw(times[r])); \
        sum += times[r]; \
    } \
    double avg = sum / RUNS; \
    printf("  %-14s %8.4f ms  %7.1f GB/s\n", label, avg, bw(avg)); \
} while (0)

static void bench_aos() {
    constexpr size_t bytes = N * sizeof(C2);
    auto *in = (C2 *)numa_alloc(bytes);
    ft_aos(in, N);
    CPU_BENCH("AoS", c_aos(in, N));
    numa_free(in, bytes);
}

static void bench_soa() {
    constexpr size_t bytes = N * sizeof(double);
    double *ri = (double *)numa_alloc(bytes);
    double *ii = (double *)numa_alloc(bytes);
    ft_soa(ri, ii, N);
    CPU_BENCH("SoA", c_soa(ii, N));
    numa_free(ri, bytes);
    numa_free(ii, bytes);
}

template<int VL>
static void bench_aosoa(const char *label) {
    constexpr size_t bytes = (N / VL) * sizeof(C2V<VL>);
    auto *in = (C2V<VL> *)numa_alloc(bytes);
    ft_aosoa<VL>(in, N);
    CPU_BENCH(label, c_aosoa<VL>(in, N));
    numa_free(in, bytes);
}

int main() {
    PAGE_SZ = sysconf(_SC_PAGESIZE);
    flush_init();

    csv = fopen("results_cpu_inplace.csv", "w");
    fprintf(csv, "layout,run,ms,gbps\n");

    printf("conj (inplace): N=%lldM  runs=%d  threads=%d  page=%ld\n",
           (long long)(N >> 20), (int)RUNS, omp_get_max_threads(), PAGE_SZ);
    printf("BW = 2N×8 (read im + write im)  |  cache flush before each run\n");
    printf("NUMA: per-thread mbind(MPOL_BIND) + first-touch\n");
    print_mems_allowed();

    {
        int nodes[8] = {};
        #pragma omp parallel
        {
            int n = get_numa_node();
            if (n >= 0 && n < 8) {
                #pragma omp atomic
                nodes[n]++;
            }
        }
        printf("  thread spread:");
        for (int i = 0; i < 8; i++)
            if (nodes[i]) printf(" N%d=%d", i, nodes[i]);
        printf("\n\n");
    }

    printf("[CPU]\n");
    bench_aos();
    bench_soa();
    bench_aosoa< 2>("AoSoA-2");
    bench_aosoa< 4>("AoSoA-4");
    bench_aosoa< 8>("AoSoA-8");
    bench_aosoa<16>("AoSoA-16");
    bench_aosoa<32>("AoSoA-32");
    bench_aosoa<64>("AoSoA-64");
    bench_aosoa<128>("AoSoA-128");
    bench_aosoa<256>("AoSoA-256");
    bench_aosoa<512>("AoSoA-512");
    fclose(csv);
    flush_free();
    printf("\nwrote results_cpu_inplace.csv\n");
}