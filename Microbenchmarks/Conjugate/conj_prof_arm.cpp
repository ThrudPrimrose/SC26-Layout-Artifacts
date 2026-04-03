/*
 * Conjugate in-place profiling — ARM Neoverse V2 (GH200 Grace CPU)
 * =================================================================
 * Diagnose the high variance seen on Grace in the conjugate benchmark.
 *
 * ARM PMU events via PERF_TYPE_RAW:
 *   - L1D_TLB_REFILL (0x05)    — L1 DTLB miss (refill from L2 TLB)
 *   - L2D_TLB_REFILL (0x2D)    — L2 DTLB miss = page table walk
 *   - L1D_CACHE_REFILL (0x03)  — L1 DCache miss
 *   - L2D_CACHE_REFILL (0x17)  — L2 DCache miss (goes to SLC/DRAM)
 *   - L1D_CACHE_WB (0x15)      — L1 writebacks
 *   - MEM_ACCESS (0x13)         — total memory accesses
 *   - STALL_BACKEND (0x24)      — backend stall cycles
 *   - BUS_ACCESS (0x19)         — off-chip bus transactions
 *
 * Discover events on Grace:
 *   perf list | grep -i tlb
 *   perf list | grep -i cache
 *   perf list | grep -i stall
 *   ls /sys/bus/event_source/devices/armv8_pmuv3_0/events/
 *
 * Compile:
 *   g++ -O3 -march=native -mtune=native -fopenmp -ffast-math \
 *       -std=c++17 -o conj_prof_arm conj_prof_arm.cpp -lnuma
 *
 * Run:
 *   OMP_NUM_THREADS=72 OMP_PROC_BIND=close OMP_PLACES=threads ./conj_prof_arm
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
 *  PMC infrastructure — ARM Neoverse V2
 *
 *  ARM PMUv3 events are encoded as raw event numbers.
 *  Use PERF_TYPE_RAW with the event ID as config.
 *
 *  Key Neoverse V2 events (from ARM TRM):
 *    0x01  L1I_CACHE_REFILL     — L1 ICache miss
 *    0x03  L1D_CACHE_REFILL     — L1 DCache miss (refill)
 *    0x04  L1D_CACHE            — L1 DCache access
 *    0x05  L1D_TLB_REFILL       — L1 DTLB miss → served by L2 TLB
 *    0x08  INST_RETIRED         — instructions retired
 *    0x11  CPU_CYCLES           — CPU cycles
 *    0x13  MEM_ACCESS            — data memory accesses
 *    0x15  L1D_CACHE_WB         — L1 DCache writebacks
 *    0x16  L2D_CACHE            — L2 DCache access
 *    0x17  L2D_CACHE_REFILL     — L2 DCache miss → SLC/DRAM
 *    0x18  L2D_CACHE_WB         — L2 DCache writebacks
 *    0x19  BUS_ACCESS           — off-chip bus transactions
 *    0x24  STALL_BACKEND        — backend stall cycles (mem latency)
 *    0x25  L1D_TLB              — L1 DTLB access count
 *    0x2D  L2D_TLB_REFILL       — L2 DTLB miss = PAGE TABLE WALK
 *    0x2F  L2D_TLB              — L2 DTLB access count
 *
 *  Discover what's available:
 *    perf list                    # full list
 *    perf list | grep -i arm     # ARM-specific
 *    cat /sys/bus/event_source/devices/armv8_pmuv3_0/events/*
 * ══════════════════════════════════════════════════════════════════════ */

static long sys_perf_event_open(struct perf_event_attr *a, pid_t p,
                                int cpu, int grp, unsigned long fl) {
    return syscall(__NR_perf_event_open, a, p, cpu, grp, fl);
}

struct PmcDef { const char *name; uint32_t type; uint64_t config; };

static const PmcDef PMC_EVENTS[] = {
    /* Generic Linux (portable) */
    {"cycles",            PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
    {"instructions",      PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS},

    /* ARM Neoverse V2 raw events */
    {"L1D-access",        PERF_TYPE_RAW, 0x04},   /* L1D_CACHE             */
    {"L1D-refill",        PERF_TYPE_RAW, 0x03},   /* L1D_CACHE_REFILL      */
    {"L1D-writeback",     PERF_TYPE_RAW, 0x15},   /* L1D_CACHE_WB          */
    {"L2D-access",        PERF_TYPE_RAW, 0x16},   /* L2D_CACHE             */
    {"L2D-refill",        PERF_TYPE_RAW, 0x17},   /* L2D_CACHE_REFILL      */
    {"L2D-writeback",     PERF_TYPE_RAW, 0x18},   /* L2D_CACHE_WB          */
    {"mem-access",        PERF_TYPE_RAW, 0x13},   /* MEM_ACCESS            */
    {"bus-access",        PERF_TYPE_RAW, 0x19},   /* BUS_ACCESS (off-chip) */

    /* TLB hierarchy */
    {"L1D-TLB-access",   PERF_TYPE_RAW, 0x25},   /* L1D_TLB               */
    {"L1D-TLB-refill",   PERF_TYPE_RAW, 0x05},   /* L1D_TLB_REFILL        */
    {"L2D-TLB-access",   PERF_TYPE_RAW, 0x2F},   /* L2D_TLB               */
    {"L2D-TLB-refill",   PERF_TYPE_RAW, 0x2D},   /* L2D_TLB_REFILL = walk */

    /* Stalls — key for understanding variance */
    {"stall-backend",     PERF_TYPE_RAW, 0x24},   /* STALL_BACKEND         */
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
                        fprintf(stderr, "  [PMC] %-20s unavailable (errno=%d)\n",
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
    printf("  Available:");
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

static void print_numa_topo() {
    int nodes[8] = {};
    #pragma omp parallel
    {
        int nd = get_numa_node();
        if (nd >= 0 && nd < 8) {
            #pragma omp atomic
            nodes[nd]++;
        }
    }
    printf("  thread→NUMA:");
    for (int i = 0; i < 8; i++)
        if (nodes[i]) printf("  N%d=%d", i, nodes[i]);
    printf("\n");

    /* Check if threads span multiple NUMA nodes */
    int n_active = 0;
    for (int i = 0; i < 8; i++) if (nodes[i]) n_active++;
    if (n_active > 1)
        printf("  WARNING: threads span %d NUMA nodes — "
               "cross-node traffic may cause variance!\n", n_active);
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
 *
 *  Grace has 64KB L1D + 1MB L2 per core + ~117MB SLC.
 *  1<<27 doubles = 1 GB >> 117MB SLC, should be enough.
 *  We also add DC CIVAC (Clean+Invalidate by VA to PoC) for
 *  the actual benchmark buffers as an extra measure.
 * ══════════════════════════════════════════════════════════════════════ */

constexpr int64_t FLUSH_N = 1 << 27;   /* 1 GB */
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
    #pragma omp parallel
    { __sync_synchronize(); }   /* dmb ish on ARM */
}

/* Explicit DC CIVAC on a buffer — flushes through all cache levels */
static void arm_dc_civac(void *base, size_t bytes) {
    #pragma omp parallel
    {
        int tid  = omp_get_thread_num();
        int nthr = omp_get_num_threads();
        size_t chunk = (bytes + nthr - 1) / nthr;
        size_t lo = tid * chunk;
        size_t hi = std::min(lo + chunk, bytes);
        /* DC CIVAC: Clean and Invalidate by VA to Point of Coherency */
        for (size_t off = lo; off < hi; off += 64)
            asm volatile("dc civac, %0" :: "r"((char*)base + off) : "memory");
    }
    __sync_synchronize();
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
 *  Bench driver with PMC + per-run diagnostics
 * ══════════════════════════════════════════════════════════════════════ */

static FILE *csv_bw;
static FILE *csv_pmc;

static double bw_ip(int P, int64_t n_base, double sec) {
    return 2.0 * P * n_base * sizeof(double) / (sec * 1e9);
}

struct PmcResult {
    int P;
    std::string layout;
    double med_gbps;
    int64_t n_base;
    uint64_t tot[N_PMC];
    /* Per-run variance diagnostics */
    double gbps_min, gbps_max, gbps_p5, gbps_p95;
    std::vector<double> all_gbps;
};
static std::vector<PmcResult> g_results;

/* Flush the actual benchmark buffer via DC CIVAC before each run */
static void *g_active_buf = nullptr;
static size_t g_active_bytes = 0;

template<typename Fn>
static void bench_with_pmc(int P, int64_t n_base, const char *layout, Fn &&kernel) {
    for (int w = 0; w < WARMUP; w++) { flush_caches(); kernel(); }

    PmcResult pr;
    pr.P = P;
    pr.layout = layout;
    pr.n_base = n_base;
    memset(pr.tot, 0, sizeof(pr.tot));

    uint64_t cnt[N_PMC];

    for (int r = 0; r < RUNS; r++) {
        flush_caches();
        /* Also DC CIVAC the benchmark buffer for deterministic cold start */
        if (g_active_buf && g_active_bytes > 0)
            arm_dc_civac(g_active_buf, g_active_bytes);

        pmc_enable_all();
        double t0 = omp_get_wtime();
        kernel();
        double t1 = omp_get_wtime();
        pmc_read_all(cnt);

        double dt = t1 - t0;
        double gbs = bw_ip(P, n_base, dt);
        pr.all_gbps.push_back(gbs);

        fprintf(csv_bw, "%d,%s,%d,%.6f,%.2f\n", P, layout, r, dt * 1e3, gbs);
        for (int i = 0; i < N_PMC; i++) pr.tot[i] += cnt[i];
    }

    std::sort(pr.all_gbps.begin(), pr.all_gbps.end());
    pr.med_gbps = pr.all_gbps[RUNS / 2];
    pr.gbps_min = pr.all_gbps.front();
    pr.gbps_max = pr.all_gbps.back();
    pr.gbps_p5  = pr.all_gbps[(int)(0.05 * RUNS)];
    pr.gbps_p95 = pr.all_gbps[(int)(0.95 * RUNS)];
    g_results.push_back(pr);

    double spread = (pr.gbps_p95 - pr.gbps_p5) / pr.med_gbps * 100.0;
    printf("  P=%-2d %-14s med=%7.1f  [%7.1f..%7.1f]  P5-P95=%6.1f%%",
           P, layout, pr.med_gbps, pr.gbps_min, pr.gbps_max, spread);
    if (spread > 10.0)
        printf("  *** HIGH VARIANCE");
    printf("\n");
}

/* ══════════════════════════════════════════════════════════════════════
 *  Per-P bench orchestration
 * ══════════════════════════════════════════════════════════════════════ */

template<int P>
static void bench_aos(int64_t n) {
    size_t bytes = n * sizeof(AoS<P>);
    auto *buf = (AoS<P> *)numa_alloc(bytes); init_aos<P>(buf, n);
    g_active_buf = buf; g_active_bytes = bytes;
    bench_with_pmc(P, n, "AoS", [&]() { kern_aos<P>(buf, n); });
    g_active_buf = nullptr; g_active_bytes = 0;
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
    /* Can't DC CIVAC all SoA buffers easily, skip it */
    g_active_buf = nullptr; g_active_bytes = 0;
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
    g_active_buf = buf; g_active_bytes = bytes;
    bench_with_pmc(P, n, label, [&]() { kern_aosoa<P,VL>(buf, n); });
    g_active_buf = nullptr; g_active_bytes = 0;
    numa_free(buf, bytes);
}

template<int P>
static void run_all() {
    int64_t n = (TOTAL_DOUBLES / (2 * P) / MAX_VL) * MAX_VL;
    printf("\n── P=%d complex pairs  N_base=%lld  total=%.1f GB ──\n",
           P, (long long)n, 2.0 * P * n * 8 / 1e9);

    bench_aos<P>(n);
    bench_soa<P>(n);
    bench_aosoa<P,  16>(n, "AoSoA-16");
    bench_aosoa<P,  64>(n, "AoSoA-64");
    bench_aosoa<P, 128>(n, "AoSoA-128");
    bench_aosoa<P, 512>(n, "AoSoA-512");
}

/* ══════════════════════════════════════════════════════════════════════
 *  PMC summary + variance analysis
 * ══════════════════════════════════════════════════════════════════════ */

static void print_pmc_table() {
    if (g_results.empty() || !g_pmc_active) return;

    printf("\n");
    printf("═══════════════════════════════════════════════════════════"
           "═══════════════════════════════════════════════════════\n");
    printf("  PMC Summary  (per-element = sum / (N_base × RUNS),  "
           "RUNS=%d  threads=%d)\n", RUNS, g_nthreads);
    printf("═══════════════════════════════════════════════════════════"
           "═══════════════════════════════════════════════════════\n");

    printf("  %2s %-14s %7s %6s", "P", "layout", "GB/s", "var%");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) printf("  %14s", PMC_EVENTS[i].name);
    printf("\n");

    printf("  %2s %-14s %7s %6s", "--", "--------------", "-------", "------");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) printf("  %14s", "--------------");
    printf("\n");

    int prev_P = -1;
    for (auto &pr : g_results) {
        if (pr.P != prev_P) {
            if (prev_P >= 0) printf("\n");
            prev_P = pr.P;
        }
        double elems = (double)pr.n_base * RUNS;
        double spread = (pr.gbps_p95 - pr.gbps_p5) / pr.med_gbps * 100.0;
        printf("  %2d %-14s %7.1f %5.1f%%", pr.P, pr.layout.c_str(),
               pr.med_gbps, spread);
        for (int i = 0; i < N_PMC; i++)
            if (g_pmc_avail[i])
                printf("  %14.4f", (double)pr.tot[i] / elems);
        printf("\n");
    }

    printf("═══════════════════════════════════════════════════════════"
           "═══════════════════════════════════════════════════════\n");

    /* ── Variance diagnosis ── */
    printf("\n  Variance analysis (P5–P95 spread as %% of median):\n");
    printf("  %3s %-14s %7s %7s %7s %7s %6s  %s\n",
           "P", "layout", "P5", "med", "P95", "max", "sprd%", "diagnosis");
    printf("  %3s %-14s %7s %7s %7s %7s %6s  %s\n",
           "---", "--------------", "-------", "-------", "-------",
           "-------", "------", "-------------------");

    for (auto &pr : g_results) {
        double spread = (pr.gbps_p95 - pr.gbps_p5) / pr.med_gbps * 100.0;
        const char *diag = "OK";
        if (spread > 30.0)      diag = "SEVERE — likely NUMA or scheduler";
        else if (spread > 15.0) diag = "HIGH — check cache flush / NUMA";
        else if (spread > 5.0)  diag = "moderate — some OS jitter";

        printf("  %3d %-14s %7.1f %7.1f %7.1f %7.1f %5.1f%%  %s\n",
               pr.P, pr.layout.c_str(),
               pr.gbps_p5, pr.med_gbps, pr.gbps_p95, pr.gbps_max,
               spread, diag);
    }

    /* ── TLB comparison: SoA vs AoS ── */
    printf("\n  TLB miss ratio (SoA / AoS) per P:\n");
    printf("  %3s", "P");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i] && strstr(PMC_EVENTS[i].name, "TLB"))
            printf("  %14s", PMC_EVENTS[i].name);
    printf("  %14s  %14s\n", "stall-backend", "L2D-refill");

    for (int P : {3, 6, 9, 12, 15, 18, 21}) {
        PmcResult *aos = nullptr, *soa = nullptr;
        for (auto &r : g_results) {
            if (r.P == P && r.layout == "AoS") aos = &r;
            if (r.P == P && r.layout == "SoA") soa = &r;
        }
        if (!aos || !soa) continue;
        printf("  %3d", P);
        for (int i = 0; i < N_PMC; i++) {
            if (!g_pmc_avail[i] || !strstr(PMC_EVENTS[i].name, "TLB")) continue;
            if (aos->tot[i] == 0)
                printf("  %14s", "n/a");
            else
                printf("  %13.2fx", (double)soa->tot[i] / (double)aos->tot[i]);
        }
        /* Also stall-backend and L2D-refill ratios */
        for (const char *nm : {"stall-backend", "L2D-refill"}) {
            bool found = false;
            for (int i = 0; i < N_PMC; i++) {
                if (!g_pmc_avail[i] || strcmp(PMC_EVENTS[i].name, nm) != 0) continue;
                if (aos->tot[i] == 0)
                    printf("  %14s", "n/a");
                else
                    printf("  %13.2fx", (double)soa->tot[i] / (double)aos->tot[i]);
                found = true;
            }
            if (!found) printf("  %14s", "n/a");
        }
        printf("\n");
    }
}

static void write_pmc_csv() {
    if (g_results.empty()) return;

    fprintf(csv_pmc, "P,layout,med_gbps,p5_gbps,p95_gbps,spread_pct,n_base");
    for (int i = 0; i < N_PMC; i++)
        if (g_pmc_avail[i]) fprintf(csv_pmc, ",%s", PMC_EVENTS[i].name);
    fprintf(csv_pmc, "\n");

    for (auto &pr : g_results) {
        double elems = (double)pr.n_base * RUNS;
        double spread = (pr.gbps_p95 - pr.gbps_p5) / pr.med_gbps * 100.0;
        fprintf(csv_pmc, "%d,%s,%.2f,%.2f,%.2f,%.1f,%lld",
                pr.P, pr.layout.c_str(), pr.med_gbps,
                pr.gbps_p5, pr.gbps_p95, spread, (long long)pr.n_base);
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

    printf("Conjugate in-place profiling (ARM Neoverse V2)\n");
    printf("  TOTAL_DOUBLES=%lldM  runs=%d  warmup=%d  threads=%d  page=%ld\n",
           (long long)(TOTAL_DOUBLES >> 20), RUNS, WARMUP, nth, PAGE_SZ);
    printf("  NUMA=%s\n", HAS_NUMA ? "yes" : "fallback");

    print_numa_topo();
    printf("\n");

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
    printf("Wrote results_cpu_inplace_pmc.csv   (per-variant PMC + variance)\n");
}