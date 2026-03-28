#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <sys/mman.h>

constexpr int64_t N_PER = 1 << 24;   // 16M doubles per stream
constexpr int64_t RUNS  = 50;
constexpr int     COUNTS[] = {1, 5, 10, 15, 20, 25, 30};
constexpr int     N_COUNTS = sizeof(COUNTS) / sizeof(COUNTS[0]);

static void *numa_alloc(size_t bytes) {
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return p;
}
static void numa_free(void *p, size_t bytes) { munmap(p, bytes); }

static void ft_fill(double *buf, int64_t n, int seed) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) buf[i] = (double)((i + seed) % 997) * 0.001;
}

/* ═══ Cache flush ═══ */
constexpr int64_t FLUSH_N = 1 << 26;  // 512 MB
static double *flush_buf = nullptr;

static void flush_init() {
    constexpr size_t bytes = FLUSH_N * sizeof(double);
    flush_buf = (double *)numa_alloc(bytes);
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] = 0.0;
}

static void flush_caches() {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < FLUSH_N; i++) flush_buf[i] += 1.0;
}

static void flush_free() {
    numa_free(flush_buf, FLUSH_N * sizeof(double));
}

/* Kernel: sequential per-array sweeps, 4 elements per i per array */
static void conj_multi(double **__restrict__ ri, double **__restrict__ ii,
                       double **__restrict__ ro, double **__restrict__ io,
                       int n_arrays, int64_t n) {
    for (int k = 0; k < n_arrays; k++) {
        const double *__restrict__ rk = ri[k];
        const double *__restrict__ ik = ii[k];
        double       *__restrict__ ok = ro[k];
        double       *__restrict__ wk = io[k];
        #pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < n; i++) {
            ok[i] =  rk[i];
            wk[i] = -ik[i];
        }
    }
}

/* BW = n_arrays × 4 × N_PER × 8 */
static FILE *csv;
static double bw(int n_arrays, double ms) {
    return (double)n_arrays * 4.0 * N_PER * sizeof(double) / (ms * 1e6);
}

static void bench(int n_arrays) {
    constexpr size_t bytes = N_PER * sizeof(double);
    double *ri[n_arrays], *ii[n_arrays], *ro[n_arrays], *io[n_arrays];

    for (int k = 0; k < n_arrays; k++) {
        ri[k] = (double *)numa_alloc(bytes); ft_fill(ri[k], N_PER, 2*k);
        ii[k] = (double *)numa_alloc(bytes); ft_fill(ii[k], N_PER, 2*k+1);
        ro[k] = (double *)numa_alloc(bytes); ft_fill(ro[k], N_PER, 0);
        io[k] = (double *)numa_alloc(bytes); ft_fill(io[k], N_PER, 0);
    }

    for (int w = 0; w < 5; w++) {
        flush_caches();
        conj_multi(ri, ii, ro, io, n_arrays, N_PER);
    }

    double times[RUNS];
    for (int r = 0; r < RUNS; r++) {
        flush_caches();
        double t0 = omp_get_wtime();
        conj_multi(ri, ii, ro, io, n_arrays, N_PER);
        times[r] = (omp_get_wtime() - t0) * 1e3;
    }

    double sum = 0;
    for (int r = 0; r < RUNS; r++) {
        fprintf(csv, "%d,%d,%.6f,%.2f\n", n_arrays, r, times[r], bw(n_arrays, times[r]));
        sum += times[r];
    }
    double avg = sum / RUNS;
    printf("  arrays=%-2d  %8.4f ms  %7.1f GB/s\n", n_arrays, avg, bw(n_arrays, avg));

    for (int k = 0; k < n_arrays; k++) {
        numa_free(ri[k], bytes); numa_free(ii[k], bytes);
        numa_free(ro[k], bytes); numa_free(io[k], bytes);
    }
}

int main() {
    flush_init();
    csv = fopen("results_streams.csv", "w");
    fprintf(csv, "n_arrays,run,ms,gbps\n");

    printf("streams: N_PER=%lldM  runs=%d  threads=%d\n",
           (long long)(N_PER>>20), (int)RUNS, omp_get_max_threads());
    printf("BW = n_arrays × 4 × N_PER × 8  |  cache flush before each run\n\n");

    for (int c = 0; c < N_COUNTS; c++)
        bench(COUNTS[c]);

    fclose(csv);
    flush_free();
    printf("\nwrote results_streams.csv\n");
}