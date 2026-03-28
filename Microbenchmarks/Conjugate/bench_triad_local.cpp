#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <sys/mman.h>

constexpr int64_t N    = 1 << 30;   // 1B doubles = 8 GB per array
constexpr int64_t RUNS = 10;
constexpr double  S    = 3.0;

/* ═══ Cache flush ═══ */
constexpr int64_t FLUSH_N = 1 << 26;  // 512 MB — larger than aggregate L3
static double *flush_buf = nullptr;

static void flush_init() {
    size_t bytes = FLUSH_N * sizeof(double);
    flush_buf = (double *)mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                               MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] = 0.0;
}

static void flush_caches() {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] += 1.0;
}

static void flush_free() {
    munmap(flush_buf, FLUSH_N * sizeof(double));
}

int main() {
    const int nthreads = omp_get_max_threads();
    const size_t bytes = N * sizeof(double);
    printf("simple triad: N=%lldM (%.0f MB/array)  runs=%d  threads=%d\n",
           (long long)(N >> 20), bytes / 1e6, (int)RUNS, nthreads);

    flush_init();

    double *a = (double *)malloc(bytes);
    double *b = (double *)malloc(bytes);
    double *c = (double *)malloc(bytes);
    if (!a || !b || !c) { perror("malloc"); return 1; }

    // init (also first-touch with default policy)
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < N; i++) { a[i] = 0.0; b[i] = 1.0; c[i] = 2.0; }

    // warmup
    for (int w = 0; w < 5; w++) {
        flush_caches();
        #pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < N; i++) a[i] = b[i] + S * c[i];
    }

    // benchmark
    double times[RUNS];
    for (int r = 0; r < RUNS; r++) {
        flush_caches();
        double t0 = omp_get_wtime();
        #pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < N; i++) a[i] = b[i] + S * c[i];
        times[r] = omp_get_wtime() - t0;
    }

    std::sort(times, times + RUNS);
    double best = times[0];
    double med  = times[RUNS / 2];
    double avg  = 0;
    for (int r = 0; r < RUNS; r++) avg += times[r];
    avg /= RUNS;

    auto bw = [&](double t) { return 3.0 * bytes / t / 1e9; };
    printf("  best:   %8.4f ms  %7.1f GB/s\n", best * 1e3, bw(best));
    printf("  median: %8.4f ms  %7.1f GB/s\n", med  * 1e3, bw(med));
    printf("  avg:    %8.4f ms  %7.1f GB/s\n", avg  * 1e3, bw(avg));

    // sanity check
    printf("  a[0] = %.1f (expect %.1f)\n", a[0], b[0] + S * c[0]);

    free(a); free(b); free(c);
    flush_free();
}