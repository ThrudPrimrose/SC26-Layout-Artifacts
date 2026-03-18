#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <type_traits>

constexpr int64_t N    = 1 << 25;
constexpr int64_t BLK  = 256;
constexpr int64_t GRID = (N + BLK - 1) / BLK;
constexpr int64_t RUNS = 100;

struct C2 { double re, im; };

template<int VL>
struct C2V { double re[VL], im[VL]; };

static_assert(std::is_trivially_copyable<C2>::value, "C2 must be trivially copyable");
static_assert(sizeof(C2) == 2 * sizeof(double), "C2 has unexpected padding");

#define CHECK_C2V(VL) \
    static_assert(std::is_trivially_copyable<C2V<VL>>::value, "C2V<" #VL "> must be trivially copyable"); \
    static_assert(sizeof(C2V<VL>) == 2 * (VL) * sizeof(double), "C2V<" #VL "> has unexpected padding");
CHECK_C2V(2) CHECK_C2V(4) CHECK_C2V(8)
CHECK_C2V(16) CHECK_C2V(32) CHECK_C2V(64)
#undef CHECK_C2V

/* ═══ CPU Kernels (out-of-place) ═══ */

static void c_aos(const C2* __restrict__ in, C2* __restrict__ out, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) { out[i].re = in[i].re; out[i].im = -in[i].im; }
}

static void c_soa(const double* __restrict__ ri, const double* __restrict__ ii,
                  double* __restrict__ ro, double* __restrict__ io, int64_t n) {
    #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; i++) { ro[i] = ri[i]; io[i] = -ii[i]; }
}

template<int VL>
static void c_aosoa(const C2V<VL>* __restrict__ in, C2V<VL>* __restrict__ out, int64_t n) {
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

static void emit_header() {
    fprintf(csv, "device,layout,run,ms,gbps\n");
}

static void emit_run(const char *dev, const char *layout, int run, double ms) {
    fprintf(csv, "%s,%s,%d,%.6f,%.2f\n", dev, layout, run, ms, bw(ms));
}

/* ═══ CPU bench ═══ */

#define CPU_BENCH(label, call) do { \
    call; /* warmup */ \
    call; /* warmup */ \
    call; /* warmup */ \
    call; /* warmup */ \
    call; /* warmup */ \
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

int main() {
    size_t bytes = 2ULL * N * sizeof(double);

    double *hi = (double*)malloc(bytes), *ho = (double*)malloc(bytes);
    for (int64_t i = 0; i < 2 * N; i++) hi[i] = (double)(i % 997) * 0.001;

    csv = fopen("results_cpu.csv", "w");
    emit_header();

    printf("conj(oop): N=%lldM  runs=%d  omp_threads=%d  dtype=double\n\n",
           (long long)(N >> 20), (int)RUNS, omp_get_max_threads());

    printf("[CPU]\n");
    CPU_BENCH("AoS",      (c_aos((C2*)hi, (C2*)ho, N)));
    CPU_BENCH("SoA",      (c_soa(hi, hi+N, ho, ho+N, N)));
    CPU_BENCH("AoSoA-2",  (c_aosoa< 2>((C2V< 2>*)hi, (C2V< 2>*)ho, N)));
    CPU_BENCH("AoSoA-4",  (c_aosoa< 4>((C2V< 4>*)hi, (C2V< 4>*)ho, N)));
    CPU_BENCH("AoSoA-8",  (c_aosoa< 8>((C2V< 8>*)hi, (C2V< 8>*)ho, N)));
    CPU_BENCH("AoSoA-16", (c_aosoa<16>((C2V<16>*)hi, (C2V<16>*)ho, N)));
    CPU_BENCH("AoSoA-32", (c_aosoa<32>((C2V<32>*)hi, (C2V<32>*)ho, N)));
    CPU_BENCH("AoSoA-64", (c_aosoa<64>((C2V<64>*)hi, (C2V<64>*)ho, N)));

    fclose(csv);
    free(hi); free(ho);
    printf("\nwrote results_cpu.csv\n");
}
