#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <sys/mman.h>
#include <type_traits>

constexpr int64_t N    = 1 << 27;   // 128M elements
constexpr int64_t BLK  = 256;
constexpr int64_t GRID = (N + BLK - 1) / BLK;
constexpr int64_t RUNS = 100;

struct C2 { double re, im; };

template<int VL>
struct C2V { double re[VL], im[VL]; };

static_assert(std::is_trivially_copyable<C2>::value, "");
static_assert(sizeof(C2) == 2 * sizeof(double), "");
#define CHECK_C2V(VL) \
    static_assert(std::is_trivially_copyable<C2V<VL>>::value, ""); \
    static_assert(sizeof(C2V<VL>) == 2 * (VL) * sizeof(double), "");
CHECK_C2V(2) CHECK_C2V(4) CHECK_C2V(8)
CHECK_C2V(16) CHECK_C2V(32) CHECK_C2V(64)
#undef CHECK_C2V

/* ═══ NUMA-aware allocation ═══ */

static void *numa_alloc(size_t bytes) {
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    return p;
}

static void numa_free(void *p, size_t bytes) { munmap(p, bytes); }

/* ═══ First-touch — matches compute schedule(static) ═══ */

static void ft_aos(C2 *buf, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        buf[i].re = (double)((2*i)   % 997) * 0.001;
        buf[i].im = (double)((2*i+1) % 997) * 0.001;
    }
}

static void ft_soa(double *r, double *im, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) {
        r[i]  = (double)((2*i)   % 997) * 0.001;
        im[i] = (double)((2*i+1) % 997) * 0.001;
    }
}

template<int VL>
static void ft_aosoa(C2V<VL> *buf, int64_t n) {
    int64_t nblks = n / VL;
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < nblks; b++)
        for (int l = 0; l < VL; l++) {
            buf[b].re[l] = (double)((2*(b*VL+l))   % 997) * 0.001;
            buf[b].im[l] = (double)((2*(b*VL+l)+1) % 997) * 0.001;
        }
}

/* ═══ CPU Kernels (out-of-place) ═══ */

static void c_aos(const C2 *__restrict__ in, C2 *__restrict__ out, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) { out[i].re = in[i].re; out[i].im = -in[i].im; }
}

static void c_soa(const double *__restrict__ ri, const double *__restrict__ ii,
                  double *__restrict__ ro, double *__restrict__ io, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) { ro[i] = ri[i]; io[i] = -ii[i]; }
}

template<int VL>
static void c_aosoa(const C2V<VL> *__restrict__ in, C2V<VL> *__restrict__ out, int64_t n) {
    int64_t nblks = n / VL;
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < nblks; b++) {
        #pragma omp simd
        for (int l = 0; l < VL; l++) {
            out[b].re[l] =  in[b].re[l];
            out[b].im[l] = -in[b].im[l];
        }
    }
}

/* ═══ Reporting ═══ */

static FILE *csv;
static double bw(double ms) { return 4.0 * N * sizeof(double) / (ms * 1e6); }

static void emit_run(const char *dev, const char *layout, int run, double ms) {
    fprintf(csv, "%s,%s,%d,%.6f,%.2f\n", dev, layout, run, ms, bw(ms));
}

/* ═══ Bench macro ═══ */

#define CPU_BENCH(label, call) do { \
    for (int w = 0; w < 5; w++) { call; } \
    double times[RUNS]; \
    for (int r = 0; r < RUNS; r++) { \
        double t0 = omp_get_wtime(); \
        call; \
        times[r] = (omp_get_wtime() - t0) * 1e3; \
    } \
    double sum = 0; \
    for (int r = 0; r < RUNS; r++) { \
        emit_run("CPU", label, r, times[r]); \
        sum += times[r]; \
    } \
    double avg = sum / RUNS; \
    printf("  %-14s %8.4f ms  %7.1f GB/s\n", label, avg, bw(avg)); \
} while (0)

/* ═══ Per-layout alloc + first-touch + bench + free ═══ */

static void bench_aos() {
    constexpr size_t bytes = N * sizeof(C2);
    auto *in  = (C2 *)numa_alloc(bytes);
    auto *out = (C2 *)numa_alloc(bytes);
    ft_aos(in,  N);
    ft_aos(out, N);
    CPU_BENCH("AoS", (c_aos(in, out, N)));
    numa_free(in, bytes); numa_free(out, bytes);
}

static void bench_soa() {
    constexpr size_t bytes = N * sizeof(double);
    double *ri = (double *)numa_alloc(bytes);
    double *ii = (double *)numa_alloc(bytes);
    double *ro = (double *)numa_alloc(bytes);
    double *io = (double *)numa_alloc(bytes);
    ft_soa(ri, ii, N);
    ft_soa(ro, io, N);
    CPU_BENCH("SoA", (c_soa(ri, ii, ro, io, N)));
    numa_free(ri, bytes); numa_free(ii, bytes);
    numa_free(ro, bytes); numa_free(io, bytes);
}

template<int VL>
static void bench_aosoa(const char *label) {
    constexpr size_t bytes = (N / VL) * sizeof(C2V<VL>);
    auto *in  = (C2V<VL> *)numa_alloc(bytes);
    auto *out = (C2V<VL> *)numa_alloc(bytes);
    ft_aosoa<VL>(in,  N);
    ft_aosoa<VL>(out, N);
    CPU_BENCH(label, (c_aosoa<VL>(in, out, N)));
    numa_free(in, bytes); numa_free(out, bytes);
}

/* ═══ Main ═══ */

int main() {
    csv = fopen("results_cpu_oop.csv", "w");
    fprintf(csv, "device,layout,run,ms,gbps\n");

    printf("conj (oop): N=%lldM  block=%ld  runs=%d  omp_threads=%d  dtype=double\n",
           (long long)(N >> 20), (long)BLK, (int)RUNS, omp_get_max_threads());
    printf("NUMA first-touch: per-layout mmap + matching schedule(static)\n\n");

    printf("[CPU]\n");
    bench_aos();
    bench_soa();
    bench_aosoa< 2>("AoSoA-2");
    bench_aosoa< 4>("AoSoA-4");
    bench_aosoa< 8>("AoSoA-8");
    bench_aosoa<16>("AoSoA-16");
    bench_aosoa<32>("AoSoA-32");
    bench_aosoa<64>("AoSoA-64");

    fclose(csv);
    printf("\nwrote results_cpu_oop.csv\n");
}