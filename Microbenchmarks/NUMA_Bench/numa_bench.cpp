#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <numaif.h>
#include <unistd.h>
#include <sched.h>

/*
 *  Benchmark to empirically determine cost model parameters:
 *
 *  ═══ PART A: Latency-based (single-threaded) ═══
 *
 *    ALPHA  = latency(local, streaming) / latency(local, random)
 *    BETA   = stride threshold (in cache lines) below which prefetcher helps
 *    GAMMA  = latency(remote) / latency(local)  for same access pattern
 *
 *  Tests:
 *    1. Pointer chasing (random latency):  local NUMA  vs  remote NUMA
 *    2. Strided sequential (stride = 1,2,4,8,16,32 cache lines):  local vs remote
 *
 *  All single-threaded to measure pure memory latency without contention.
 *  Array size >> LLC to ensure DRAM-bound.
 *
 *  ═══ PART B: Bandwidth-based (multi-threaded, saturating) ═══
 *
 *    ALPHA_BW  = BW(local, streaming) / BW(local, random)
 *    BETA_BW   = stride threshold where BW drops below 50% of streaming
 *    GAMMA_BW  = BW(local) / BW(remote)  for same access pattern
 *
 *  Tests:
 *    1. Streaming read (sequential):  local vs remote, all threads
 *    2. Random gather (index-based):  local vs remote, all threads
 *    3. Strided read (stride = 1,2,4,8,16,32 CL): local vs remote, all threads
 *
 *  Multi-threaded to saturate memory controllers.
 */

static long PAGE_SZ;
static int  NUM_NODES;

static int get_numa_node() {
    unsigned cpu, node;
    syscall(__NR_getcpu, &cpu, &node, nullptr);
    return (int)node;
}

/* ═══ Allocate on a specific NUMA node ═══ */

static void *alloc_on_node(size_t bytes, int node) {
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_HUGEPAGE);

    unsigned long mask = 1UL << node;
    if (mbind(p, bytes, MPOL_BIND, &mask, 64, MPOL_MF_MOVE | MPOL_MF_STRICT) != 0)
        perror("mbind");

    /* First-touch to fault pages */
    volatile char *cp = (volatile char *)p;
    for (size_t off = 0; off < bytes; off += PAGE_SZ)
        cp[off] = 0;

    return p;
}

static void nfree(void *p, size_t bytes) { munmap(p, bytes); }

/* ═══════════════════════════════════════════════════════════════════
 *  PART A: LATENCY BENCHMARKS (single-threaded)
 * ═══════════════════════════════════════════════════════════════════ */

/* ═══ Pointer chasing: random permutation walk ═══ */

static double bench_pointer_chase(int node, int64_t n_lines) {
    size_t bytes = n_lines * 64;
    int64_t *buf = (int64_t *)alloc_on_node(bytes, node);

    std::vector<int64_t> perm(n_lines);
    std::iota(perm.begin(), perm.end(), 0);
    std::mt19937_64 rng(42);
    for (int64_t i = n_lines - 1; i > 0; i--) {
        std::uniform_int_distribution<int64_t> dist(0, i - 1);
        int64_t j = dist(rng);
        std::swap(perm[i], perm[j]);
    }

    for (int64_t i = 0; i < n_lines; i++)
        buf[i * 8] = perm[i] * 64;

    { int64_t idx = 0;
      for (int64_t i = 0; i < n_lines; i++)
          idx = *(int64_t *)((char *)buf + idx); }

    int reps = 3;
    double best = 1e30;
    for (int r = 0; r < reps; r++) {
        int64_t idx = 0;
        double t0 = omp_get_wtime();
        for (int64_t i = 0; i < n_lines; i++)
            idx = *(volatile int64_t *)((char *)buf + idx);
        double t1 = omp_get_wtime();
        double ns = (t1 - t0) * 1e9 / n_lines;
        if (ns < best) best = ns;
        if (idx == -999) printf("%ld\n", idx);
    }

    nfree(buf, bytes);
    return best;
}

/* ═══ Strided sequential access (latency, single-threaded) ═══ */

static double bench_strided(int node, int64_t n_lines, int stride_cls) {
    size_t bytes = n_lines * 64;
    int64_t *buf = (int64_t *)alloc_on_node(bytes, node);
    memset(buf, 0xAB, bytes);

    int64_t n_accesses = n_lines / stride_cls;

    { int64_t sum = 0;
      for (int64_t i = 0; i < n_accesses; i++)
          sum += buf[i * stride_cls * 8];
      if (sum == -999) printf("%ld\n", sum); }

    int reps = 5;
    double best = 1e30;
    for (int r = 0; r < reps; r++) {
        int64_t sum = 0;
        double t0 = omp_get_wtime();
        for (int64_t i = 0; i < n_accesses; i++)
            sum += *(volatile int64_t *)&buf[i * stride_cls * 8];
        double t1 = omp_get_wtime();
        double ns = (t1 - t0) * 1e9 / n_accesses;
        if (ns < best) best = ns;
        if (sum == -999) printf("%ld\n", sum);
    }

    nfree(buf, bytes);
    return best;
}

/* ═══════════════════════════════════════════════════════════════════
 *  PART B: BANDWIDTH BENCHMARKS (multi-threaded, saturating)
 * ═══════════════════════════════════════════════════════════════════ */

/*
 *  All BW benchmarks return GB/s of useful data touched.
 *  "Useful data" = number of distinct cache lines accessed × 64 bytes.
 *  Multi-threaded: each thread works on its own contiguous chunk.
 */

/* ═══ Streaming read bandwidth ═══
 *  Sequential scan, stride=1 CL.  Maximally prefetcher-friendly.
 */
static double bench_bw_stream(int node, int64_t n_lines, int nthreads) {
    size_t bytes = n_lines * 64;
    int64_t *buf = (int64_t *)alloc_on_node(bytes, node);
    memset(buf, 0xAB, bytes);

    int64_t lines_per_t = n_lines / nthreads;
    double data_gb = (double)(lines_per_t * nthreads * 64) / 1e9;

    /* Warmup */
    #pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        int64_t start = tid * lines_per_t;
        int64_t end   = start + lines_per_t;
        int64_t sum = 0;
        for (int64_t i = start; i < end; i++)
            sum += buf[i * 8];
        if (sum == -999) printf("%ld\n", sum);
    }

    int reps = 5;
    double best = 0.0;
    for (int r = 0; r < reps; r++) {
        double t0 = omp_get_wtime();
        #pragma omp parallel num_threads(nthreads)
        {
            int tid = omp_get_thread_num();
            int64_t start = tid * lines_per_t;
            int64_t end   = start + lines_per_t;
            int64_t sum = 0;
            for (int64_t i = start; i < end; i++)
                sum += *(volatile int64_t *)&buf[i * 8];
            if (sum == -999) printf("%ld\n", sum);
        }
        double t1 = omp_get_wtime();
        double gbs = data_gb / (t1 - t0);
        if (gbs > best) best = gbs;
    }

    nfree(buf, bytes);
    return best;
}

/* ═══ Strided read bandwidth ═══
 *  Each thread reads every stride-th cache line from its chunk.
 *  Returns GB/s of cache lines actually touched.
 */
static double bench_bw_strided(int node, int64_t n_lines, int stride_cls, int nthreads) {
    size_t bytes = n_lines * 64;
    int64_t *buf = (int64_t *)alloc_on_node(bytes, node);
    memset(buf, 0xAB, bytes);

    int64_t lines_per_t = n_lines / nthreads;
    int64_t accesses_per_t = lines_per_t / stride_cls;
    double data_gb = (double)(accesses_per_t * nthreads * 64) / 1e9;  /* CL touched */

    /* Warmup */
    #pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        int64_t base = tid * lines_per_t;
        int64_t sum = 0;
        for (int64_t i = 0; i < accesses_per_t; i++)
            sum += buf[(base + i * stride_cls) * 8];
        if (sum == -999) printf("%ld\n", sum);
    }

    int reps = 5;
    double best = 0.0;
    for (int r = 0; r < reps; r++) {
        double t0 = omp_get_wtime();
        #pragma omp parallel num_threads(nthreads)
        {
            int tid = omp_get_thread_num();
            int64_t base = tid * lines_per_t;
            int64_t sum = 0;
            for (int64_t i = 0; i < accesses_per_t; i++)
                sum += *(volatile int64_t *)&buf[(base + i * stride_cls) * 8];
            if (sum == -999) printf("%ld\n", sum);
        }
        double t1 = omp_get_wtime();
        double gbs = data_gb / (t1 - t0);
        if (gbs > best) best = gbs;
    }

    nfree(buf, bytes);
    return best;
}

/* ═══ Random gather bandwidth ═══
 *  Each thread reads cache lines at random (pre-shuffled) indices.
 *  Measures sustained BW under fully random access.
 */
static double bench_bw_random(int node, int64_t n_lines, int nthreads) {
    size_t bytes = n_lines * 64;
    int64_t *buf = (int64_t *)alloc_on_node(bytes, node);
    memset(buf, 0xAB, bytes);

    int64_t lines_per_t = n_lines / nthreads;
    double data_gb = (double)(lines_per_t * nthreads * 64) / 1e9;

    /* Build per-thread random index arrays */
    std::vector<std::vector<int64_t>> indices(nthreads);
    for (int t = 0; t < nthreads; t++) {
        indices[t].resize(lines_per_t);
        std::iota(indices[t].begin(), indices[t].end(), t * lines_per_t);
        std::mt19937_64 rng(42 + t);
        std::shuffle(indices[t].begin(), indices[t].end(), rng);
    }

    /* Warmup */
    #pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        const int64_t *idx = indices[tid].data();
        int64_t sum = 0;
        for (int64_t i = 0; i < lines_per_t; i++)
            sum += buf[idx[i] * 8];
        if (sum == -999) printf("%ld\n", sum);
    }

    int reps = 5;
    double best = 0.0;
    for (int r = 0; r < reps; r++) {
        double t0 = omp_get_wtime();
        #pragma omp parallel num_threads(nthreads)
        {
            int tid = omp_get_thread_num();
            const int64_t *idx = indices[tid].data();
            int64_t sum = 0;
            for (int64_t i = 0; i < lines_per_t; i++)
                sum += *(volatile int64_t *)&buf[idx[i] * 8];
            if (sum == -999) printf("%ld\n", sum);
        }
        double t1 = omp_get_wtime();
        double gbs = data_gb / (t1 - t0);
        if (gbs > best) best = gbs;
    }

    nfree(buf, bytes);
    return best;
}

/* ═══ Main ═══ */

int main(int argc, char **argv) {
    PAGE_SZ = sysconf(_SC_PAGESIZE);

    int nthreads_bw = omp_get_max_threads();
    if (argc > 1) nthreads_bw = atoi(argv[1]);

    /* Pin to core 0 for latency tests */
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);

    int local_node = get_numa_node();

    /* Detect available NUMA nodes */
    std::vector<int> nodes;
    for (int n = 0; n < 16; n++) {
        char path[64];
        snprintf(path, sizeof(path), "/sys/devices/system/node/node%d", n);
        if (access(path, F_OK) == 0) nodes.push_back(n);
    }
    NUM_NODES = (int)nodes.size();

    int remote_node = -1;
    for (int n : nodes)
        if (n != local_node) remote_node = n;
    if (remote_node < 0) {
        fprintf(stderr, "Only 1 NUMA node found, cannot measure remote latency.\n");
        remote_node = local_node;
    }

    constexpr int64_t N_LINES = 1 << 22;  /* 4M cache lines = 256 MB */
    constexpr int strides[] = {1, 2, 4, 8, 16, 32};
    constexpr int N_STRIDES = sizeof(strides) / sizeof(strides[0]);

    printf("NUMA latency + bandwidth probe\n");
    printf("  pinned to core 0 (node %d) for latency tests\n", local_node);
    printf("  local node  = %d\n", local_node);
    printf("  remote node = %d\n", remote_node);
    printf("  array       = %lld CL = %lld MB\n",
           (long long)N_LINES, (long long)(N_LINES * 64 / (1 << 20)));
    printf("  BW threads  = %d\n", nthreads_bw);
    printf("  NUMA nodes  = %d\n\n", NUM_NODES);

    FILE *csv = fopen("results_numa_latency.csv", "w");
    fprintf(csv, "test,node,stride_cl,value,unit,location\n");

    /* ══════════════════════════════════════════════════════════════
     *  PART A: LATENCY
     * ══════════════════════════════════════════════════════════════ */

    printf("  ╔══════════════════════════════════════════════════════════╗\n");
    printf("  ║  PART A: LATENCY (single-threaded, ns/CL)              ║\n");
    printf("  ╚══════════════════════════════════════════════════════════╝\n\n");

    printf("  %-40s %10s %10s %10s\n", "Test", "Local(ns)", "Remote(ns)", "Ratio");
    printf("  %-40s %10s %10s %10s\n", "----", "---------", "----------", "-----");

    double chase_local  = bench_pointer_chase(local_node, N_LINES);
    double chase_remote = bench_pointer_chase(remote_node, N_LINES);
    fprintf(csv, "lat_pointer_chase,%d,0,%.2f,ns,local\n", local_node, chase_local);
    fprintf(csv, "lat_pointer_chase,%d,0,%.2f,ns,remote\n", remote_node, chase_remote);
    printf("  %-40s %10.1f %10.1f %10.2f\n",
           "pointer chase (random)", chase_local, chase_remote, chase_remote / chase_local);

    printf("\n");
    double stride_local[N_STRIDES], stride_remote[N_STRIDES];
    for (int si = 0; si < N_STRIDES; si++) {
        stride_local[si]  = bench_strided(local_node, N_LINES, strides[si]);
        stride_remote[si] = bench_strided(remote_node, N_LINES, strides[si]);
        fprintf(csv, "lat_strided,%d,%d,%.2f,ns,local\n", local_node, strides[si], stride_local[si]);
        fprintf(csv, "lat_strided,%d,%d,%.2f,ns,remote\n", remote_node, strides[si], stride_remote[si]);

        char label[64];
        snprintf(label, sizeof(label), "stride %2d CL (%4d B)", strides[si], strides[si] * 64);
        printf("  %-40s %10.1f %10.1f %10.2f\n",
               label, stride_local[si], stride_remote[si],
               stride_remote[si] / stride_local[si]);
    }

    /* Latency-derived parameters */
    printf("\n  ═══ Latency-derived cost model parameters ═══\n\n");

    double lat_random_local  = chase_local;
    double lat_stream_local  = stride_local[0];
    double lat_random_remote = chase_remote;
    double lat_stream_remote = stride_remote[0];

    double alpha_lat = lat_stream_local / lat_random_local;
    double gamma_lat_random = lat_random_remote / lat_random_local;
    double gamma_lat_stream = lat_stream_remote / lat_stream_local;

    printf("  ALPHA_LAT = lat_stream / lat_random  (local)\n");
    printf("            = %.1f / %.1f = %.4f\n\n", lat_stream_local, lat_random_local, alpha_lat);

    printf("  GAMMA_LAT (random)  = %.1f / %.1f = %.4f\n",
           lat_random_remote, lat_random_local, gamma_lat_random);
    printf("  GAMMA_LAT (stream)  = %.1f / %.1f = %.4f\n\n",
           lat_stream_remote, lat_stream_local, gamma_lat_stream);

    double threshold_lat = 2.0 * lat_stream_local;
    int beta_lat = 0;
    printf("  BETA_LAT (stride where latency > 2× streaming):\n");
    for (int si = 0; si < N_STRIDES; si++) {
        bool over = stride_local[si] > threshold_lat;
        printf("    stride %2d CL: %6.1f ns  %s  (threshold=%.1f)\n",
               strides[si], stride_local[si], over ? "ABOVE" : "below", threshold_lat);
        if (!over) beta_lat = strides[si];
    }
    printf("  BETA_LAT = %d\n\n", beta_lat);

    /* ══════════════════════════════════════════════════════════════
     *  PART B: BANDWIDTH
     * ══════════════════════════════════════════════════════════════ */

    /* Unpin — let OpenMP spread threads across cores for BW tests */
    CPU_ZERO(&cpuset);
    for (int i = 0; i < CPU_SETSIZE; i++) CPU_SET(i, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);

    printf("  ╔══════════════════════════════════════════════════════════╗\n");
    printf("  ║  PART B: BANDWIDTH (%d threads, GB/s)                  ║\n", nthreads_bw);
    printf("  ╚══════════════════════════════════════════════════════════╝\n\n");

    printf("  %-40s %10s %10s %10s\n", "Test", "Local", "Remote", "Ratio");
    printf("  %-40s %10s %10s %10s\n", "----", "------", "------", "-----");

    /* Streaming */
    double bw_stream_local  = bench_bw_stream(local_node,  N_LINES, nthreads_bw);
    double bw_stream_remote = bench_bw_stream(remote_node, N_LINES, nthreads_bw);
    fprintf(csv, "bw_stream,%d,1,%.2f,GB/s,local\n",  local_node, bw_stream_local);
    fprintf(csv, "bw_stream,%d,1,%.2f,GB/s,remote\n", remote_node, bw_stream_remote);
    printf("  %-40s %10.1f %10.1f %10.2f\n",
           "streaming read (stride 1 CL)",
           bw_stream_local, bw_stream_remote, bw_stream_local / bw_stream_remote);

    /* Random gather */
    double bw_random_local  = bench_bw_random(local_node,  N_LINES, nthreads_bw);
    double bw_random_remote = bench_bw_random(remote_node, N_LINES, nthreads_bw);
    fprintf(csv, "bw_random,%d,0,%.2f,GB/s,local\n",  local_node, bw_random_local);
    fprintf(csv, "bw_random,%d,0,%.2f,GB/s,remote\n", remote_node, bw_random_remote);
    printf("  %-40s %10.1f %10.1f %10.2f\n",
           "random gather",
           bw_random_local, bw_random_remote, bw_random_local / bw_random_remote);

    /* Strided */
    printf("\n");
    double bw_stride_local[N_STRIDES], bw_stride_remote[N_STRIDES];
    for (int si = 0; si < N_STRIDES; si++) {
        bw_stride_local[si]  = bench_bw_strided(local_node,  N_LINES, strides[si], nthreads_bw);
        bw_stride_remote[si] = bench_bw_strided(remote_node, N_LINES, strides[si], nthreads_bw);
        fprintf(csv, "bw_strided,%d,%d,%.2f,GB/s,local\n",
                local_node, strides[si], bw_stride_local[si]);
        fprintf(csv, "bw_strided,%d,%d,%.2f,GB/s,remote\n",
                remote_node, strides[si], bw_stride_remote[si]);

        char label[64];
        snprintf(label, sizeof(label), "stride %2d CL (%4d B)", strides[si], strides[si] * 64);
        printf("  %-40s %10.1f %10.1f %10.2f\n",
               label, bw_stride_local[si], bw_stride_remote[si],
               bw_stride_local[si] / bw_stride_remote[si]);
    }

    /* BW-derived parameters */
    printf("\n  ═══ Bandwidth-derived cost model parameters ═══\n\n");

    double alpha_bw = bw_stream_local / bw_random_local;
    double gamma_bw_random = bw_stream_local / bw_stream_remote;  /* note: inverted vs latency */
    double gamma_bw_stream = bw_stream_local / bw_stream_remote;

    printf("  ALPHA_BW = BW_stream / BW_random  (local, higher = bigger penalty)\n");
    printf("           = %.1f / %.1f = %.4f\n\n", bw_stream_local, bw_random_local, alpha_bw);

    double gamma_bw_r = bw_random_local / bw_random_remote;
    printf("  GAMMA_BW (stream)  = BW_local / BW_remote\n");
    printf("           = %.1f / %.1f = %.4f\n", bw_stream_local, bw_stream_remote, gamma_bw_stream);
    printf("  GAMMA_BW (random)  = BW_local / BW_remote\n");
    printf("           = %.1f / %.1f = %.4f\n\n", bw_random_local, bw_random_remote, gamma_bw_r);

    double threshold_bw = 0.5 * bw_stream_local;
    int beta_bw = 0;
    printf("  BETA_BW (stride where BW < 50%% of streaming):\n");
    for (int si = 0; si < N_STRIDES; si++) {
        bool under = bw_stride_local[si] < threshold_bw;
        printf("    stride %2d CL: %6.1f GB/s  %s  (threshold=%.1f)\n",
               strides[si], bw_stride_local[si], under ? "BELOW" : "above", threshold_bw);
        if (!under) beta_bw = strides[si];
    }
    printf("  BETA_BW = %d\n\n", beta_bw);

    /* ══════════════════════════════════════════════════════════════
     *  SUMMARY: side-by-side comparison
     * ══════════════════════════════════════════════════════════════ */

    printf("  ╔══════════════════════════════════════════════════════════╗\n");
    printf("  ║  SUMMARY: Latency vs Bandwidth derived parameters      ║\n");
    printf("  ╚══════════════════════════════════════════════════════════╝\n\n");

    printf("  %-20s %12s %12s\n", "Parameter", "Latency", "Bandwidth");
    printf("  %-20s %12s %12s\n", "---------", "-------", "---------");
    printf("  %-20s %12.4f %12.4f\n", "ALPHA", alpha_lat, alpha_bw);
    printf("  %-20s %12d %12d\n",     "BETA  (CL)", beta_lat, beta_bw);
    printf("  %-20s %12.4f %12.4f\n", "GAMMA (stream)", gamma_lat_stream, gamma_bw_stream);
    printf("  %-20s %12.4f %12.4f\n", "GAMMA (random)", gamma_lat_random, gamma_bw_r);
    printf("  %-20s %12d %12d\n",     "P_NUMA", NUM_NODES, NUM_NODES);
    printf("\n");

    /* ══════════════════════════════════════════════════════════════
     *  Cross-NUMA latency matrix
     * ══════════════════════════════════════════════════════════════ */

    if (NUM_NODES > 1) {
        printf("  ═══ Cross-NUMA pointer chase latency matrix (ns) ═══\n\n");
        printf("  %8s", "from\\to");
        for (int to : nodes) printf("  node%-4d", to);
        printf("\n");

        for (int from : nodes) {
            char path[128];
            snprintf(path, sizeof(path), "/sys/devices/system/node/node%d/cpulist", from);
            FILE *f = fopen(path, "r");
            int first_cpu = 0;
            if (f) { fscanf(f, "%d", &first_cpu); fclose(f); }

            CPU_ZERO(&cpuset);
            CPU_SET(first_cpu, &cpuset);
            sched_setaffinity(0, sizeof(cpuset), &cpuset);

            printf("  node%-4d", from);
            for (int to : nodes) {
                double lat = bench_pointer_chase(to, N_LINES);
                fprintf(csv, "lat_matrix,%d->%d,0,%.2f,ns,%s\n",
                        from, to, lat, (from == to) ? "local" : "remote");
                printf("  %8.1f", lat);
            }
            printf("\n");
        }

        /* ── Cross-NUMA BW matrix ── */
        printf("\n  ═══ Cross-NUMA streaming BW matrix (GB/s, %d threads) ═══\n\n",
               nthreads_bw);
        printf("  %8s", "from\\to");
        for (int to : nodes) printf("  node%-4d", to);
        printf("\n");

        /* Unpin for BW */
        CPU_ZERO(&cpuset);
        for (int i = 0; i < CPU_SETSIZE; i++) CPU_SET(i, &cpuset);
        sched_setaffinity(0, sizeof(cpuset), &cpuset);

        for (int from : nodes) {
            /* Pin threads to 'from' node for BW matrix */
            char path[128];
            snprintf(path, sizeof(path), "/sys/devices/system/node/node%d/cpulist", from);
            FILE *f = fopen(path, "r");
            char cpulist[256] = {};
            if (f) { fgets(cpulist, sizeof(cpulist), f); fclose(f); }

            /* Parse cpulist "0-71" or "0,2,4" etc → set affinity */
            cpu_set_t node_cpuset;
            CPU_ZERO(&node_cpuset);
            /* Simple parse: first_cpu-last_cpu */
            int lo = 0, hi = 0;
            if (sscanf(cpulist, "%d-%d", &lo, &hi) == 2) {
                for (int c = lo; c <= hi; c++) CPU_SET(c, &node_cpuset);
            } else {
                /* Comma-separated fallback */
                char *tok = strtok(cpulist, ",\n");
                while (tok) {
                    int a, b;
                    if (sscanf(tok, "%d-%d", &a, &b) == 2)
                        for (int c = a; c <= b; c++) CPU_SET(c, &node_cpuset);
                    else
                        CPU_SET(atoi(tok), &node_cpuset);
                    tok = strtok(nullptr, ",\n");
                }
            }
            sched_setaffinity(0, sizeof(node_cpuset), &node_cpuset);

            printf("  node%-4d", from);
            for (int to : nodes) {
                double bw = bench_bw_stream(to, N_LINES, nthreads_bw);
                fprintf(csv, "bw_matrix,%d->%d,1,%.2f,GB/s,%s\n",
                        from, to, bw, (from == to) ? "local" : "remote");
                printf("  %8.1f", bw);
            }
            printf("\n");
        }
    }

    fclose(csv);

    /* Restore */
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);

    printf("\n  wrote results_numa_latency.csv\n");
    return 0;
}