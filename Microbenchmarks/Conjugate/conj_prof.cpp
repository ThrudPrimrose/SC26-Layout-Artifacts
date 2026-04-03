/*
 * Conjugate in-place profiling benchmark (CPU) — NUMA-aware, per-thread PMC
 * ==========================================================================
 * Negates the imaginary part of P complex arrays in-place.
 * Sweeps P = 3,6,9,12,15,18,21 and layouts AoS/SoA/AoSoA-{2..512}.
 *
 * Per-thread perf_event_open: each OMP thread opens its own counters.
 * Main thread enable/disable/reads all fds around the kernel call.
 * Sum across threads = true total.
 *
 * PMC events include generic Linux + AMD Zen4 raw events:
 *   - cycles, instructions, cache-refs, cache-misses
 *   - L1d load misses, DTLB load misses  (generic HW_CACHE)
 *   - ls_l1_d_tlb_miss (AMD raw 0x45, umask 0xFF) — L1 DTLB misses all page sizes
 *
 * Compile:
 *   g++ -O3 -march=native -mtune=native -fopenmp -ffast-math \
 *       -std=c++17 -o conj_prof conj_prof.cpp -lnuma
 *
 * Run:
 *   OMP_NUM_THREADS=96 OMP_PROC_BIND=close ./conj_prof
 *
 * To discover AMD raw events on your system:
 *   perf list | grep -i tlb
 *   ls /sys/bus/event_source/devices/cpu/events/
 *   perf stat -e dTLB-load-misses,r00ff45 -a sleep 0.01
 */

#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cerrno>
#include <algorithm>
#include <vector>
#include <string>

#include <sys/mman.h>
#include <sys/syscall.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <unistd.h>

#if __has_include(<numaif.h>)
  #include <numaif.h>
  #define HAS_NUMA 1
#else
  #define HAS_NUMA 0
  #define MPOL_BIND 2
  static int mbind(void*,size_t,int,const unsigned long*,unsigned long,unsigned){return 0;}
#endif

/* ══════════════════════════════════════════════════════════════════════
 *  Config
 * ══════════════════════════════════════════════════════════════════════ */

constexpr int64_t TOTAL_DOUBLES = 1LL << 28;  /* ~2 GB */
constexpr int     RUNS   = 50;
constexpr int     WARMUP = 5;
constexpr int     MAX_VL = 512;
constexpr int     MAX_THREADS = 512;

/* ══════════════════════════════════════════════════════════════════════
 *  Layout structs
 * ══════════════════════════════════════════════════════════════════════ */

template<int P>
struct AoS {
    struct { double re, im; } c[P];
};

template<int P, int VL>
struct AoSoA {
    struct { double re[VL], im[VL]; } c[P];
};

/* ══════════════════════════════════════════════════════════════════════
 *  PMC infrastructure
 * ══════════════════════════════════════════════════════════════════════ */

static long sys_perf_event_open(struct perf_event_attr *a, pid_t p,
                                int cpu, int grp, unsigned long fl) {
    return syscall(__NR_perf_event_open, a, p, cpu, grp, fl);
}

struct PmcDef { const char *name; uint32_t type; uint64_t config; };

/*
 * Event catalog:
 *   [0..5] Generic Linux events (work on any x86)
 *   [6]    AMD raw: ls_l1_d_tlb_miss.all (event=0x45, umask=0xFF)
 *          → counts ALL L1 DTLB misses (4K + 2M + 1G pages)
 *   [7]    AMD raw: ls_l1_d_tlb_miss.4k_only (event=0x45, umask=0x01)
 *   [8]    AMD raw: ls_l1_d_tlb_miss.2m_only (event=0x45, umask=0x04)
 *   [9]    AMD raw: ls_l1_d_tlb_miss.1g_only (event=0x45, umask=0x08)
 *   [10]   AMD raw: ls_l1_d_tlb_miss.coalesced (event=0x45, umask=0x02)
 *
 *  On non-AMD, the raw events simply fail to open and are skipped.
 *
 *  AMD Zen4 PPR (Processor Programming Reference):
 *    PMCx045 [LS L1 DTLB Miss]
 *      UnitMask bits:
 *        [0] = TlbReload4kPage
 *        [1] = TlbReloadCoalescedPage
 *        [2] = TlbReload2mPage
 *        [3] = TlbReload1gPage
 *
 *  To encode raw:  config = (umask << 8) | event_select
 *    e.g. all L1 DTLB miss = (0xFF << 8) | 0x45 = 0xFF45
 *
 *  Discover more events:
 *    perf list | grep -i tlb
 *    perf list | grep -i cache
 *    cat /sys/bus/event_source/devices/cpu/events/*
 */

static const PmcDef PMC_EVENTS[] = {
    /* Generic Linux */
    {"cycles",          PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
    {"instructions",    PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS},
    {"cache-refs",      PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_REFERENCES},
    {"cache-misses",    PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES},
    {"L1d-load-miss",   PERF_TYPE_HW_CACHE,
        (PERF_COUNT_HW_CACHE_L1D) |
        ((uint64_t)PERF_COUNT_HW_CACHE_OP_READ  << 8) |
        ((uint64_t)PERF_COUNT_HW_CACHE_RESULT_MISS << 16)},
    {"dTLB-load-miss",  PERF_TYPE_HW_CACHE,
        (PERF_COUNT_HW_CACHE_DTLB) |
        ((uint64_t)PERF_COUNT_HW_CACHE_OP_READ  << 8) |
        ((uint64_t)PERF_COUNT_HW_CACHE_RESULT_MISS << 16)},

    /* AMD Zen4 raw: PMCx045 ls_l1_d_tlb_miss */
    {"amd:L1dTLB-miss-all",       PERF_TYPE_RAW, 0xFF45},  /* umask=0xFF */
    {"amd:L1dTLB-miss-4K",        PERF_TYPE_RAW, 0x0145},  /* umask=0x01 */
    {"amd:L1dTLB-miss-coalesced", PERF_TYPE_RAW, 0x0245},  /* umask=0x02 */
    {"amd:L1dTLB-miss-2M",        PERF_TYPE_RAW, 0x0445},  /* umask=0x04 */
    {"amd:L1dTLB-miss-1G",        PERF_TYPE_RAW, 0x0845},  /* umask=0x08 */
};
static constexpr int N_PMC = sizeof(PMC_EVENTS) / sizeof(PMC_EVENTS[0]);

static int  g_pmc_fds[MAX_THREADS][N_PMC];
static bool g_pmc_avail[N_PMC];
static bool g_pmc_active = false;
static int  g_nthreads = 0;

static void pmc_open_all() {
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp master
        {
            g_nthreads = omp_get_num_threads();
            for (int i = 0; i < N_PMC; i++) g_pmc_avail[i] = true;
        }
        #pragma omp barrier

        for (int i = 0; i < N_PMC; i++) {
            struct perf_event_attr pe = {};
            pe.size = sizeof(pe);
            pe.type = PMC_EVENTS[i].type;
            pe.config = PMC_EVENTS[i].config;
            pe.disabled = 1;
            pe.exclude_kernel = 1;
            pe.exclude_hv = 1;
            int fd = (int)sys_perf_event_open(&pe, 0, -1, -1, 0);
            g_pmc_fds[tid][i] = fd;
            if (fd < 0) {
                #pragma omp critical
                {
                    if (g_pmc_avail[i]) {
                        fprintf(stderr, "  [PMC] %-28s unavailable (errno=%d)\n",
                                PMC_EVENTS[i].name, errno);
                        g_pmc_avail[i] = false;
                    }
                }
            }
        }
    }

    g_pmc_active = false;
    int ncounters = 0;
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) { g_pmc_active = true; ncounters++; }

    printf("  PMC: %s — %d events × %d threads\n",
           g_pmc_active ? "ACTIVE" : "DISABLED", ncounters, g_nthreads);
    printf("  Available events:");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) printf(" %s", PMC_EVENTS[i].name);
    printf("\n");
}

static void pmc_close_all() {
    for (int t = 0; t < g_nthreads; t++)
        for (int i = 0; i < N_PMC; i++)
            if (g_pmc_fds[t][i] >= 0) close(g_pmc_fds[t][i]);
}

static void pmc_enable_all() {
    if (!g_pmc_active) return;
    for (int t = 0; t < g_nthreads; t++)
        for (int i = 0; i < N_PMC; i++)
            if (g_pmc_fds[t][i] >= 0) {
                ioctl(g_pmc_fds[t][i], PERF_EVENT_IOC_RESET, 0);
                ioctl(g_pmc_fds[t][i], PERF_EVENT_IOC_ENABLE, 0);
            }
}

static void pmc_read_all(uint64_t out[N_PMC]) {
    memset(out, 0, N_PMC * sizeof(uint64_t));
    if (!g_pmc_active) return;
    for (int t = 0; t < g_nthreads; t++)
        for (int i = 0; i < N_PMC; i++) {
            int fd = g_pmc_fds[t][i];
            if (fd < 0) continue;
            ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
            uint64_t val = 0;
            if (read(fd, &val, sizeof(val)) == sizeof(val))
                out[i] += val;
        }
}

/* ══════════════════════════════════════════════════════════════════════
 *  NUMA helpers
 * ══════════════════════════════════════════════════════════════════════ */

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

/* ══════════════════════════════════════════════════════════════════════
 *  Cache flush
 * ══════════════════════════════════════════════════════════════════════ */

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

/* ══════════════════════════════════════════════════════════════════════
 *  Init helpers
 * ══════════════════════════════════════════════════════════════════════ */

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

/* ══════════════════════════════════════════════════════════════════════
 *  Kernels (in-place): negate im only
 * ══════════════════════════════════════════════════════════════════════ */

template<int P>
static void kern_aos(AoS<P> *__restrict__ buf, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++)
        for (int p = 0; p < P; p++)
            buf[i].c[p].im = -buf[i].c[p].im;
}

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

/* ══════════════════════════════════════════════════════════════════════
 *  Bench driver with PMC
 * ══════════════════════════════════════════════════════════════════════ */

static FILE *csv_bw;   /* per-run BW */
static FILE *csv_pmc;  /* per-variant PMC totals */

/* Useful BW: P im fields read + P im fields written = 2·P·N_base·8 */
static double bw_ip(int P, int64_t n_base, double sec) {
    return 2.0 * P * n_base * sizeof(double) / (sec * 1e9);
}

struct PmcResult {
    int P;
    std::string layout;
    double med_gbps;
    int64_t n_base;
    uint64_t tot[N_PMC];
};
static std::vector<PmcResult> g_results;

template<typename Fn>
static void bench_with_pmc(int P, int64_t n_base, const char *layout, Fn &&kernel) {
    /* warmup */
    for (int w = 0; w < WARMUP; w++) { flush_caches(); kernel(); }

    PmcResult pr;
    pr.P = P;
    pr.layout = layout;
    pr.n_base = n_base;
    memset(pr.tot, 0, sizeof(pr.tot));

    std::vector<double> gbps_v;
    uint64_t cnt[N_PMC];

    for (int r = 0; r < RUNS; r++) {
        flush_caches();

        pmc_enable_all();
        double t0 = omp_get_wtime();
        kernel();
        double t1 = omp_get_wtime();
        pmc_read_all(cnt);

        double dt = t1 - t0;
        double gbs = bw_ip(P, n_base, dt);
        gbps_v.push_back(gbs);

        fprintf(csv_bw, "%d,%s,%d,%.6f,%.2f\n", P, layout, r, dt * 1e3, gbs);

        for (int i = 0; i < N_PMC; i++) pr.tot[i] += cnt[i];
    }

    std::sort(gbps_v.begin(), gbps_v.end());
    pr.med_gbps = gbps_v[RUNS / 2];
    g_results.push_back(pr);

    printf("  P=%-2d %-14s %7.1f GB/s  [%7.1f..%7.1f]\n",
           P, layout, pr.med_gbps, gbps_v.front(), gbps_v.back());
}

/* ══════════════════════════════════════════════════════════════════════
 *  Per-P bench orchestration
 * ══════════════════════════════════════════════════════════════════════ */

template<int P>
static void bench_aos(int64_t n) {
    size_t bytes = n * sizeof(AoS<P>);
    auto *buf = (AoS<P> *)numa_alloc(bytes); init_aos<P>(buf, n);
    bench_with_pmc(P, n, "AoS", [&]() { kern_aos<P>(buf, n); });
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
    bench_with_pmc(P, n, "SoA", [&]() { kern_soa<P>(im, n); });
    for (int p = 0; p < P; p++) {
        numa_free(re[p], bytes); numa_free(im[p], bytes);
    }
}

template<int P, int VL>
static void bench_aosoa(int64_t n, const char *label) {
    int64_t nblks = n / VL;
    size_t bytes = nblks * sizeof(AoSoA<P,VL>);
    auto *buf = (AoSoA<P,VL> *)numa_alloc(bytes); init_aosoa<P,VL>(buf, n);
    bench_with_pmc(P, n, label, [&]() { kern_aosoa<P,VL>(buf, n); });
    numa_free(buf, bytes);
}

template<int P>
static void run_all() {
    int64_t n = (TOTAL_DOUBLES / (2 * P) / MAX_VL) * MAX_VL;
    printf("\n── P=%d complex pairs  N_base=%lld  total=%.1f GB ──\n",
           P, (long long)n, 2.0 * P * n * 8 / 1e9);

    bench_aos<P>(n);
    bench_soa<P>(n);
    bench_aosoa<P,   8>(n, "AoSoA-8");
    bench_aosoa<P,  16>(n, "AoSoA-16");
    bench_aosoa<P,  32>(n, "AoSoA-32");
    bench_aosoa<P,  64>(n, "AoSoA-64");
    bench_aosoa<P, 128>(n, "AoSoA-128");
    bench_aosoa<P, 256>(n, "AoSoA-256");
    bench_aosoa<P, 512>(n, "AoSoA-512");
}

/* ══════════════════════════════════════════════════════════════════════
 *  PMC summary table
 * ══════════════════════════════════════════════════════════════════════ */

static void print_pmc_table() {
    if (g_results.empty() || !g_pmc_active) return;

    /* Header */
    printf("\n");
    printf("═══════════════════════════════════════════════════════════"
           "═══════════════════════════════════════════════════════\n");
    printf("  PMC Summary  (per-element = sum_threads / (N_base × RUNS),  "
           "RUNS=%d  threads=%d)\n", RUNS, g_nthreads);
    printf("═══════════════════════════════════════════════════════════"
           "═══════════════════════════════════════════════════════\n");

    printf("  %2s %-14s %7s", "P", "layout", "GB/s");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) printf("  %16s", PMC_EVENTS[i].name);
    printf("\n");

    printf("  %2s %-14s %7s", "--", "--------------", "-------");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) printf("  %16s", "----------------");
    printf("\n");

    int prev_P = -1;
    for (auto &pr : g_results) {
        if (pr.P != prev_P) {
            if (prev_P >= 0) printf("\n");
            prev_P = pr.P;
        }
        double elems = (double)pr.n_base * RUNS;
        printf("  %2d %-14s %7.1f", pr.P, pr.layout.c_str(), pr.med_gbps);
        for (int i = 0; i < N_PMC; i++)
            if (g_pmc_avail[i])
                printf("  %16.4f", (double)pr.tot[i] / elems);
        printf("\n");
    }

    printf("═══════════════════════════════════════════════════════════"
           "═══════════════════════════════════════════════════════\n");

    /* ── Key comparisons: SoA vs AoS DTLB miss ratio per P ── */
    printf("\n  DTLB miss ratio  (SoA / AoS)  per P:\n");
    printf("  %3s", "P");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i] &&
            (strstr(PMC_EVENTS[i].name, "TLB") || strstr(PMC_EVENTS[i].name, "tlb")))
            printf("  %16s", PMC_EVENTS[i].name);
    printf("\n");

    for (int P : {3, 6, 9, 12, 15, 18, 21}) {
        PmcResult *aos = nullptr, *soa = nullptr;
        for (auto &r : g_results) {
            if (r.P == P && r.layout == "AoS") aos = &r;
            if (r.P == P && r.layout == "SoA") soa = &r;
        }
        if (!aos || !soa) continue;
        printf("  %3d", P);
        for (int i = 0; i < N_PMC; i++) {
            if (!g_pmc_avail[i]) continue;
            if (!strstr(PMC_EVENTS[i].name, "TLB") &&
                !strstr(PMC_EVENTS[i].name, "tlb")) continue;
            if (aos->tot[i] == 0)
                printf("  %16s", "n/a");
            else
                printf("  %15.2fx", (double)soa->tot[i] / (double)aos->tot[i]);
        }
        printf("\n");
    }
}

static void write_pmc_csv() {
    if (g_results.empty()) return;

    fprintf(csv_pmc, "P,layout,med_gbps,n_base");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) fprintf(csv_pmc, ",%s", PMC_EVENTS[i].name);
    fprintf(csv_pmc, "\n");

    for (auto &pr : g_results) {
        double elems = (double)pr.n_base * RUNS;
        fprintf(csv_pmc, "%d,%s,%.2f,%lld",
                pr.P, pr.layout.c_str(), pr.med_gbps, (long long)pr.n_base);
        for (int i = 0; i < N_PMC; i++)
            if (g_pmc_avail[i])
                fprintf(csv_pmc, ",%.6f", (double)pr.tot[i] / elems);
        fprintf(csv_pmc, "\n");
    }
}

/* ══════════════════════════════════════════════════════════════════════
 *  Main
 * ══════════════════════════════════════════════════════════════════════ */

int main() {
    PAGE_SZ = sysconf(_SC_PAGESIZE);
    flush_init();

    int nth = 0;
    #pragma omp parallel
    { nth = omp_get_num_threads(); }

    printf("Conjugate in-place profiling (CPU)\n");
    printf("  TOTAL_DOUBLES=%lldM  runs=%d  warmup=%d  threads=%d  page=%ld\n",
           (long long)(TOTAL_DOUBLES >> 20), RUNS, WARMUP, nth, PAGE_SZ);
    printf("  NUMA=%s\n\n", HAS_NUMA ? "yes" : "fallback");

    pmc_open_all();

    csv_bw  = fopen("results_cpu_inplace_prof.csv", "w");
    csv_pmc = fopen("results_cpu_inplace_pmc.csv", "w");
    fprintf(csv_bw, "P,layout,run,ms,gbps\n");

    run_all< 3>();
    run_all< 6>();
    run_all< 9>();
    run_all<12>();
    run_all<15>();
    run_all<18>();
    run_all<21>();

    print_pmc_table();
    write_pmc_csv();

    fclose(csv_bw);
    fclose(csv_pmc);
    pmc_close_all();
    flush_free();

    printf("\nWrote results_cpu_inplace_prof.csv  (per-run BW)\n");
    printf("Wrote results_cpu_inplace_pmc.csv   (per-variant PMC summary)\n");
}