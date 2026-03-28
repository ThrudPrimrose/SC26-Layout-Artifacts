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
    madvise(p, bytes, MADV_NOHUGEPAGE);

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

/* ═══ Pointer chasing: random permutation walk ═══
 *
 *  Build a random single-cycle permutation over N cache lines.
 *  Chase pointers, measure total time / N.
 *  Returns nanoseconds per cache line access.
 */

static double bench_pointer_chase(int node, int64_t n_lines) {
    size_t bytes = n_lines * 64;
    int64_t *buf = (int64_t *)alloc_on_node(bytes, node);

    /* Build random single-cycle permutation */
    std::vector<int64_t> perm(n_lines);
    std::iota(perm.begin(), perm.end(), 0);
    std::mt19937_64 rng(42);
    /* Sattolo's algorithm for a single cycle */
    for (int64_t i = n_lines - 1; i > 0; i--) {
        std::uniform_int_distribution<int64_t> dist(0, i - 1);
        int64_t j = dist(rng);
        std::swap(perm[i], perm[j]);
    }

    /* Write next-pointer: buf[line_i * 8] = byte offset of next line */
    for (int64_t i = 0; i < n_lines; i++)
        buf[i * 8] = perm[i] * 64;  /* 8 int64's per cache line, use slot 0 */

    /* Warmup */
    {
        int64_t idx = 0;
        for (int64_t i = 0; i < n_lines; i++)
            idx = *(int64_t *)((char *)buf + idx);
    }

    /* Timed run */
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
        /* Prevent dead-code elimination */
        if (idx == -999) printf("%ld\n", idx);
    }

    nfree(buf, bytes);
    return best;
}

/* ═══ Strided sequential access ═══
 *
 *  Read every stride-th cache line sequentially.
 *  Prefetcher can detect and run ahead for small strides.
 *  Returns nanoseconds per cache line touched.
 */

static double bench_strided(int node, int64_t n_lines, int stride_cls) {
    size_t bytes = n_lines * 64;
    int64_t *buf = (int64_t *)alloc_on_node(bytes, node);

    /* Fill with non-zero so reads are not special-cased */
    memset(buf, 0xAB, bytes);

    int64_t n_accesses = n_lines / stride_cls;

    /* Warmup */
    {
        int64_t sum = 0;
        for (int64_t i = 0; i < n_accesses; i++)
            sum += buf[i * stride_cls * 8];  /* 8 int64's per CL */
        if (sum == -999) printf("%ld\n", sum);
    }

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

/* ═══ Main ═══ */

int main() {
    PAGE_SZ = sysconf(_SC_PAGESIZE);

    /* Pin to core 0 for consistent single-thread results */
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

    /* Find a remote node (pick the one with highest node id, likely farthest) */
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

    printf("NUMA latency probe: pinned to core 0 (node %d)\n", local_node);
    printf("  local node  = %d\n", local_node);
    printf("  remote node = %d\n", remote_node);
    printf("  array       = %lld cache lines = %lld MB\n",
           (long long)N_LINES, (long long)(N_LINES * 64 / (1 << 20)));
    printf("  NUMA nodes  = %d (%d available)\n\n", NUM_NODES, (int)nodes.size());

    FILE *csv = fopen("results_numa_latency.csv", "w");
    fprintf(csv, "test,node,stride_cl,ns_per_cl,location\n");

    /* ── 1. Pointer chasing ── */
    printf("  %-40s %10s %10s %10s\n", "Test", "Local(ns)", "Remote(ns)", "Ratio");
    printf("  %-40s %10s %10s %10s\n", "----", "---------", "----------", "-----");

    double chase_local  = bench_pointer_chase(local_node, N_LINES);
    double chase_remote = bench_pointer_chase(remote_node, N_LINES);
    fprintf(csv, "pointer_chase,%d,0,%.2f,local\n", local_node, chase_local);
    fprintf(csv, "pointer_chase,%d,0,%.2f,remote\n", remote_node, chase_remote);
    printf("  %-40s %10.1f %10.1f %10.2f\n",
           "pointer chase (random)", chase_local, chase_remote, chase_remote / chase_local);

    /* ── 2. Strided sequential ── */
    printf("\n");
    double stride_local[N_STRIDES], stride_remote[N_STRIDES];
    for (int si = 0; si < N_STRIDES; si++) {
        stride_local[si]  = bench_strided(local_node, N_LINES, strides[si]);
        stride_remote[si] = bench_strided(remote_node, N_LINES, strides[si]);
        fprintf(csv, "strided,%d,%d,%.2f,local\n", local_node, strides[si], stride_local[si]);
        fprintf(csv, "strided,%d,%d,%.2f,remote\n", remote_node, strides[si], stride_remote[si]);

        char label[64];
        snprintf(label, sizeof(label), "stride %2d CL (%4d B)", strides[si], strides[si] * 64);
        printf("  %-40s %10.1f %10.1f %10.2f\n",
               label, stride_local[si], stride_remote[si],
               stride_remote[si] / stride_local[si]);
    }

    fclose(csv);

    /* ── 3. Derive parameters ── */
    printf("\n  ═══ Derived cost model parameters ═══\n\n");

    double lat_random_local = chase_local;
    double lat_stream_local = stride_local[0];  /* stride=1 CL */
    double lat_random_remote = chase_remote;
    double lat_stream_remote = stride_remote[0];

    double alpha = lat_stream_local / lat_random_local;
    double gamma_random = lat_random_remote / lat_random_local;
    double gamma_stream = lat_stream_remote / lat_stream_local;

    printf("  ALPHA  = lat_stream_local / lat_random_local\n");
    printf("         = %.1f / %.1f = %.4f\n\n", lat_stream_local, lat_random_local, alpha);

    printf("  GAMMA  (random)  = lat_random_remote / lat_random_local\n");
    printf("         = %.1f / %.1f = %.4f\n", lat_random_remote, lat_random_local, gamma_random);
    printf("  GAMMA  (stream)  = lat_stream_remote / lat_stream_local\n");
    printf("         = %.1f / %.1f = %.4f\n\n", lat_stream_remote, lat_stream_local, gamma_stream);

    /* BETA: find stride threshold where latency starts approaching random */
    printf("  BETA determination (stride where latency > 2× streaming):\n");
    double threshold = 2.0 * lat_stream_local;
    int beta = 0;
    for (int si = 0; si < N_STRIDES; si++) {
        bool over = stride_local[si] > threshold;
        printf("    stride %2d CL: %6.1f ns  %s  (threshold=%.1f)\n",
               strides[si], stride_local[si], over ? "ABOVE" : "below", threshold);
        if (!over) beta = strides[si];
    }
    printf("\n  BETA   = %d  (largest stride still below threshold)\n", beta);

    printf("\n  ═══ Summary for cost_metrics.cpp ═══\n");
    printf("  ALPHA = %.4f\n", alpha);
    printf("  BETA  = %d\n", beta);
    printf("  GAMMA = %.4f  (use random ratio)\n", gamma_random);
    printf("  P_NUMA = %d\n\n", NUM_NODES);

    /* ── 4. Full latency table for all node pairs ── */
    if (NUM_NODES > 1) {
        printf("  ═══ Cross-NUMA pointer chase latency matrix (ns) ═══\n\n");
        printf("  %8s", "from\\to");
        for (int to : nodes) printf("  node%-4d", to);
        printf("\n");

        for (int from : nodes) {
            /* Pin to a core on 'from' node */
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
                fprintf(csv, "matrix,%d→%d,0,%.2f,%s\n",
                        from, to, lat, (from == to) ? "local" : "remote");
                printf("  %8.1f", lat);
            }
            printf("\n");
        }
    }

    /* Restore pinning */
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);

    printf("\n  wrote results_numa_latency.csv\n");
    return 0;
}