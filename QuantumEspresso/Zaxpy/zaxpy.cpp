// zaxpy_indirect_sweep_cpu.cpp
// CPU indirect zaxpy benchmark with NUMA-aware first-touch allocation.
//
// Same experiments as GPU version (nat20 sizes, small + 1GB tiled).
// Four kernel variants: aos_scatter, aos_sorted, soa_scatter, soa_sorted.
// Arrays allocated via mmap, first-touched in parallel for NUMA placement.
//
// Compile:
//   g++ -O3 -std=c++17 -fopenmp -lnuma zaxpy_indirect_sweep_cpu.cpp -o zaxpy_indirect_sweep_cpu
//   # or with clang: clang++ -O3 -std=c++17 -fopenmp -lnuma ...

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include <omp.h>
#include <numaif.h>
#include <sys/mman.h>
#include <unistd.h>

/* ================================================================ */
/*  NUMA helpers                                                     */
/* ================================================================ */

static constexpr int NUMA_NODES = 4;

template <typename T>
static T *numa_alloc_unfaulted(size_t count)
{
    size_t bytes = count * sizeof(T);
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return static_cast<T *>(p);
}

template <typename T>
static void numa_dealloc(T *p, size_t count)
{
    if (p) munmap(p, count * sizeof(T));
}

/**
 * First-touch a 1-D array in parallel with schedule(static) so that
 * pages land on the NUMA node of the touching thread.
 * Relies on OMP_PROC_BIND=spread + OMP_PLACES=cores.
 */
template <typename T>
static void first_touch_zero(T *arr, size_t count)
{
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < count; i++)
        arr[i] = T(0);
}

template <typename T>
static void first_touch_copy(T *dst, const T *src, size_t count)
{
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < count; i++)
        dst[i] = src[i];
}

/**
 * Flush caches: write a large buffer to evict resident data.
 * Uses same 256 MiB approach as the GPU version.
 */
static constexpr size_t FLUSH_ELEMS = 256ULL * 1024 * 1024 / sizeof(double);
static double *g_flush = nullptr;

static void init_flush()
{
    if (!g_flush) {
        g_flush = numa_alloc_unfaulted<double>(FLUSH_ELEMS);
        first_touch_zero(g_flush, FLUSH_ELEMS);
    }
}

static void flush_caches()
{
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < FLUSH_ELEMS; i++)
        g_flush[i] = 1.0;
}

/* ================================================================ */
/*  Verification                                                     */
/* ================================================================ */

static bool verify_aos(int ny, const double *gpu, const double *ref, double tol = 1e-6)
{
    for (int i = 0; i < 2*ny; i += 2) {
        double dr = gpu[i]   - ref[i];
        double di = gpu[i+1] - ref[i+1];
        double d  = std::sqrt(dr*dr + di*di);
        double m  = std::sqrt(ref[i]*ref[i] + ref[i+1]*ref[i+1]) + 1e-10;
        if (d/m > tol) { printf("AoS mismatch at %d\n", i/2); return false; }
    }
    return true;
}

static bool verify_soa(int n, const double *gr, const double *gi,
                       const double *rr, const double *ri, double tol = 1e-6)
{
    for (int i = 0; i < n; i++) {
        double dr = gr[i] - rr[i];
        double di = gi[i] - ri[i];
        double d  = std::sqrt(dr*dr + di*di);
        double m  = std::sqrt(rr[i]*rr[i] + ri[i]*ri[i]) + 1e-10;
        if (d/m > tol) { printf("SoA mismatch at %d\n", i); return false; }
    }
    return true;
}

/* ================================================================ */
/*  CPU Kernels                                                      */
/* ================================================================ */

// AoS: y is interleaved [re0,im0,re1,im1,...], 2*ny doubles
// x is interleaved [re0,im0,...], 2*nx doubles

static void kern_aos_scatter(int nx, const int *ymap,
                             const double *x, double *y)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nx; i++) {
        int yi = ymap[i];
        y[2*yi]   += x[2*i];
        y[2*yi+1] += x[2*i+1];
    }
}

static void kern_aos_sorted(int nx, const int *ymap_s, const int *xmap,
                            const double *x, double *y)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nx; i++) {
        int yi = ymap_s[i];
        int xi = xmap[i];
        y[2*yi]   += x[2*xi];
        y[2*yi+1] += x[2*xi+1];
    }
}

static void kern_soa_scatter(int nx, const int *ymap,
                             const double *xr, const double *xi,
                             double *yr, double *yi)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nx; i++) {
        int y_ = ymap[i];
        yr[y_] += xr[i];
        yi[y_] += xi[i];
    }
}

static void kern_soa_sorted(int nx, const int *ymap_s, const int *xmap,
                            const double *xr, const double *xi,
                            double *yr, double *yi)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nx; i++) {
        int y_ = ymap_s[i];
        int x_ = xmap[i];
        yr[y_] += xr[x_];
        yi[y_] += xi[x_];
    }
}

/* ================================================================ */
/*  Core profiling                                                   */
/* ================================================================ */

void profile_config(
    FILE *csv, const char *ename,
    int ny, int nx,
    const int *h_ymap_orig,
    const int *h_ymap_s_orig,
    const int *h_xmap_orig,
    int nthreads, int iters, int warmup)
{
    omp_set_num_threads(nthreads);

    // ── NUMA-aware allocation via mmap + first-touch ────────────────────

    // Index maps (read-only, split by iteration i)
    int *ymap   = numa_alloc_unfaulted<int>(nx);
    int *ymap_s = numa_alloc_unfaulted<int>(nx);
    int *xmap   = numa_alloc_unfaulted<int>(nx);
    first_touch_copy(ymap,   h_ymap_orig,   (size_t)nx);
    first_touch_copy(ymap_s, h_ymap_s_orig, (size_t)nx);
    first_touch_copy(xmap,   h_xmap_orig,   (size_t)nx);

    // AoS arrays: interleaved [re,im] → 2*n doubles
    double *y_aos_init = numa_alloc_unfaulted<double>(2 * ny);
    double *x_aos      = numa_alloc_unfaulted<double>(2 * nx);
    double *y_aos      = numa_alloc_unfaulted<double>(2 * ny);
    double *ref_aos    = (double *)malloc(2 * ny * sizeof(double));

    // SoA arrays: separate re, im
    double *yr_init = numa_alloc_unfaulted<double>(ny);
    double *yi_init = numa_alloc_unfaulted<double>(ny);
    double *xr      = numa_alloc_unfaulted<double>(nx);
    double *xi      = numa_alloc_unfaulted<double>(nx);
    double *yr      = numa_alloc_unfaulted<double>(ny);
    double *yi      = numa_alloc_unfaulted<double>(ny);
    double *ref_re  = (double *)malloc(ny * sizeof(double));
    double *ref_im  = (double *)malloc(ny * sizeof(double));

    // ── Fill random data + first-touch ──────────────────────────────────
    srand(0);
    for (int i = 0; i < ny; i++) {
        double re = (double)rand() / RAND_MAX;
        double im = (double)rand() / RAND_MAX;
        ref_aos[2*i] = re; ref_aos[2*i+1] = im;
        ref_re[i] = re; ref_im[i] = im;
    }
    // First-touch y arrays — distribute by block across NUMA nodes
    // y is written indirectly, but first-touch by block is the best we can do
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < ny; i++) {
        y_aos_init[2*i]   = ref_aos[2*i];
        y_aos_init[2*i+1] = ref_aos[2*i+1];
        yr_init[i] = ref_re[i];
        yi_init[i] = ref_im[i];
    }

    for (int i = 0; i < nx; i++) {
        double re = (double)rand() / RAND_MAX;
        double im = (double)rand() / RAND_MAX;
        // Stored in temp, will first-touch copy below
        ref_re[0] = re; // reuse buffer temporarily
        // Direct write for now, first_touch_copy for x arrays below
        x_aos[2*i] = re; x_aos[2*i+1] = im;  // will be overwritten
        xr[i] = re; xi[i] = im;               // will be overwritten
    }
    // Regenerate x with proper first-touch
    {
        // Need stable random sequence — re-seed and skip ny values
        srand(0);
        for (int i = 0; i < ny; i++) { rand(); rand(); }  // skip y values
        std::vector<double> tmp_re(nx), tmp_im(nx);
        for (int i = 0; i < nx; i++) {
            tmp_re[i] = (double)rand() / RAND_MAX;
            tmp_im[i] = (double)rand() / RAND_MAX;
        }
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < nx; i++) {
            x_aos[2*i]   = tmp_re[i];
            x_aos[2*i+1] = tmp_im[i];
            xr[i] = tmp_re[i];
            xi[i] = tmp_im[i];
        }

        // ── CPU reference (serial) ─────────────────────────────────
        for (int i = 0; i < ny; i++) {
            ref_re[i]  = y_aos_init[2*i];       // re
            ref_im[i]  = y_aos_init[2*i+1];     // im
        }
        // ref_aos already has init values
        memcpy(ref_aos, y_aos_init, 2 * ny * sizeof(double));

        for (int i = 0; i < nx; i++) {
            int yi = h_ymap_orig[i];
            ref_aos[2*yi]   += tmp_re[i];
            ref_aos[2*yi+1] += tmp_im[i];
            ref_re[yi]      += tmp_re[i];
            ref_im[yi]      += tmp_im[i];
        }
    }

    double t0, t1, ms;

    auto run = [&](const char *vname,
                   auto reset_fn, auto launch_fn, auto verify_fn)
    {
        // Verify
        reset_fn();
        launch_fn();
        if (!verify_fn()) {
            fprintf(stderr, "FAIL %s %s nthreads=%d\n", vname, ename, nthreads);
            return;
        }
        // Warmup
        for (int w = 0; w < warmup; w++) { reset_fn(); launch_fn(); }
        // Timed
        for (int r = 0; r < iters; r++) {
            reset_fn();
            flush_caches();
            t0 = omp_get_wtime();
            launch_fn();
            t1 = omp_get_wtime();
            ms = (t1 - t0) * 1000.0;
            fprintf(csv, "%s,%s,%d,%d,%d,1,%d,%.6f\n",
                    vname, ename, ny, nx, nthreads, r, ms);
        }
    };

    // ── AoS scatter ─────────────────────────────────────────────────────
    run("aos_scatter",
        [&]{ memcpy(y_aos, y_aos_init, 2*ny*sizeof(double)); },
        [&]{ kern_aos_scatter(nx, ymap, x_aos, y_aos); },
        [&]{ return verify_aos(ny, y_aos, ref_aos); }
    );

    // ── AoS sorted ──────────────────────────────────────────────────────
    run("aos_sorted",
        [&]{ memcpy(y_aos, y_aos_init, 2*ny*sizeof(double)); },
        [&]{ kern_aos_sorted(nx, ymap_s, xmap, x_aos, y_aos); },
        [&]{ return verify_aos(ny, y_aos, ref_aos); }
    );

    // ── SoA scatter ─────────────────────────────────────────────────────
    run("soa_scatter",
        [&]{
            memcpy(yr, yr_init, ny*sizeof(double));
            memcpy(yi, yi_init, ny*sizeof(double));
        },
        [&]{ kern_soa_scatter(nx, ymap, xr, xi, yr, yi); },
        [&]{ return verify_soa(ny, yr, yi, ref_re, ref_im); }
    );

    // ── SoA sorted ──────────────────────────────────────────────────────
    run("soa_sorted",
        [&]{
            memcpy(yr, yr_init, ny*sizeof(double));
            memcpy(yi, yi_init, ny*sizeof(double));
        },
        [&]{ kern_soa_sorted(nx, ymap_s, xmap, xr, xi, yr, yi); },
        [&]{ return verify_soa(ny, yr, yi, ref_re, ref_im); }
    );

    // Cleanup
    numa_dealloc(ymap, nx); numa_dealloc(ymap_s, nx); numa_dealloc(xmap, nx);
    numa_dealloc(y_aos_init, 2*ny); numa_dealloc(x_aos, 2*nx);
    numa_dealloc(y_aos, 2*ny);
    numa_dealloc(yr_init, ny); numa_dealloc(yi_init, ny);
    numa_dealloc(xr, nx); numa_dealloc(xi, nx);
    numa_dealloc(yr, ny); numa_dealloc(yi, ny);
    free(ref_aos); free(ref_re); free(ref_im);
}

/* ================================================================ */
/*  Helpers: sort, sample, I/O, tile                                 */
/* ================================================================ */

void sortWithIndices(int N, const int *in, int *sv, int *si)
{
    std::iota(si, si+N, 0);
    std::sort(si, si+N, [&in](int a, int b){ return in[a]<in[b]; });
    for (int i = 0; i < N; i++) sv[i] = in[si[i]];
}

std::vector<int> uniformSample(int nx, int ny, std::mt19937 &rng)
{
    std::vector<int> pool(ny);
    std::iota(pool.begin(), pool.end(), 0);
    for (int i = 0; i < nx; i++) {
        std::uniform_int_distribution<int> d(i, ny-1);
        std::swap(pool[i], pool[d(rng)]);
    }
    return {pool.begin(), pool.begin()+nx};
}

std::vector<int> normalSample(int nx, int ny, std::mt19937 &rng)
{
    double mu = (ny-1)/2.0, sigma = ny/6.0;
    std::normal_distribution<double> dist(mu, sigma);
    std::unordered_set<int> seen;
    std::vector<int> res; res.reserve(nx);
    while ((int)res.size() < nx) {
        int s = (int)std::round(dist(rng));
        if (s<0 || s>=ny) continue;
        if (!seen.insert(s).second) continue;
        res.push_back(s);
    }
    return res;
}

int *readBinaryFile(const std::string &fp, int &N)
{
    std::ifstream f(fp, std::ios::binary|std::ios::ate);
    if (!f.is_open()) throw std::runtime_error("Cannot open "+fp);
    auto sz = f.tellg();
    if (sz % sizeof(int) != 0) throw std::runtime_error("Bad size");
    N = sz/sizeof(int); f.seekg(0);
    int *d = new int[N];
    if (!f.read((char*)d, sz)) { delete[]d; throw std::runtime_error("Read fail"); }
    return d;
}

std::vector<int> tileIndexMap(const int *base, int nx_base, int ny_base,
                              size_t target_bytes)
{
    size_t bpt = (size_t)(ny_base + nx_base) * 16;
    int K = std::max(1, (int)(target_bytes / bpt));
    int nx_t = K * nx_base;
    std::vector<int> out(nx_t);
    for (int k = 0; k < K; k++) {
        int yo = k*ny_base, xo = k*nx_base;
        for (int i = 0; i < nx_base; i++) out[xo+i] = base[i] + yo;
    }
    printf("    Tiled K=%d  nx=%d ny=%d  footprint=%.0f MB\n",
           K, nx_t, K*ny_base,
           (double)(K*(size_t)(ny_base+nx_base)*16)/(1024.0*1024.0));
    return out;
}

/* ================================================================ */
/*  Sweep driver                                                     */
/* ================================================================ */

static const char *CSV_HDR = "variant,experiment,ny,nx,tpb,coarsen,rep,time_ms\n";
// Note: tpb column holds nthreads for CPU, coarsen is always 1.

void sweep_all_configs(
    FILE *csv, const char *ename,
    int ny, int nx,
    int *ymap, int *ymap_s, int *xmap,
    const std::vector<int> &thread_counts,
    int iters, int warmup)
{
    for (int nt : thread_counts) {
        printf("  %s nthreads=%d\n", ename, nt);
        profile_config(csv, ename, ny, nx, ymap, ymap_s, xmap,
                       nt, iters, warmup);
        fflush(csv);
    }
}

/* ================================================================ */
/*  Main                                                             */
/* ================================================================ */

int main(int argc, char **argv)
{
    constexpr int ITERS  = 100;
    constexpr int WARMUP = 10;
    constexpr int NUM_SAMPLES = 5;
    constexpr size_t TARGET_1GB = 1ULL << 30;

    constexpr int NX_BASE = 110273;
    constexpr int NY_BASE = 225001;
    const char *QE_BIN = "./bin/vexx_k_gpu__dfftt__nl_nat20.bin";

    // Thread counts to sweep — adapt to your machine
    // MI300A: 96 Zen 4 cores (4 NUMA × 24)
    // Grace:  72 Neoverse V2 cores (2 NUMA × 36) — or 144 with 2 sockets
    int max_threads = omp_get_max_threads();
    std::vector<int> thread_counts;

    // Build sweep: 1 NUMA domain, 2, 3, 4, then full
    if (max_threads >= 96) {
        // MI300A: 4 NUMA × 24 cores
        thread_counts = {24, 48, 72, 96};
    } else if (max_threads >= 72) {
        // Grace: 72 cores
        thread_counts = {36, 72};
    } else {
        // Generic: half, full
        thread_counts = {std::max(1, max_threads/2), max_threads};
    }

    printf("Thread sweep: ");
    for (int t : thread_counts) printf("%d ", t);
    printf("\n");

    std::mt19937 rng(42);
    init_flush();

    // Load QE base index map
    int N = 0;
    int *qe_base = readBinaryFile(QE_BIN, N);
    assert(N == NX_BASE);

    // ================================================================
    //  SMALL scale
    // ================================================================
    FILE *csv = fopen("zaxpy_sweep_small_cpu.csv", "w");
    fprintf(csv, "%s", CSV_HDR);

    // QE
    {
        printf("\n>>> SMALL qe  ny=%d nx=%d\n", NY_BASE, NX_BASE);
        int *ys = new int[NX_BASE], *xm = new int[NX_BASE];
        sortWithIndices(NX_BASE, qe_base, ys, xm);
        sweep_all_configs(csv, "qe", NY_BASE, NX_BASE, qe_base, ys, xm,
                          thread_counts, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }

    // Uniform
    for (int s = 0; s < NUM_SAMPLES; s++) {
        char en[64]; snprintf(en, sizeof(en), "uniform_s%d", s);
        printf("\n>>> SMALL %s\n", en);
        auto yv = uniformSample(NX_BASE, NY_BASE, rng);
        int *ys = new int[NX_BASE], *xm = new int[NX_BASE];
        sortWithIndices(NX_BASE, yv.data(), ys, xm);
        sweep_all_configs(csv, en, NY_BASE, NX_BASE, yv.data(), ys, xm,
                          thread_counts, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }

    // Normal
    for (int s = 0; s < NUM_SAMPLES; s++) {
        char en[64]; snprintf(en, sizeof(en), "normal_s%d", s);
        printf("\n>>> SMALL %s\n", en);
        auto yv = normalSample(NX_BASE, NY_BASE, rng);
        int *ys = new int[NX_BASE], *xm = new int[NX_BASE];
        sortWithIndices(NX_BASE, yv.data(), ys, xm);
        sweep_all_configs(csv, en, NY_BASE, NX_BASE, yv.data(), ys, xm,
                          thread_counts, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }
    fclose(csv);

    // ================================================================
    //  1 GB tiled
    // ================================================================
    csv = fopen("zaxpy_sweep_1gb_cpu.csv", "w");
    fprintf(csv, "%s", CSV_HDR);

    // QE tiled
    {
        printf("\n>>> 1GB qe_tiled\n");
        auto tv = tileIndexMap(qe_base, NX_BASE, NY_BASE, TARGET_1GB);
        int nx_t = (int)tv.size(), K = nx_t/NX_BASE, ny_t = K*NY_BASE;
        int *ys = new int[nx_t], *xm = new int[nx_t];
        sortWithIndices(nx_t, tv.data(), ys, xm);
        sweep_all_configs(csv, "qe_tiled", ny_t, nx_t, tv.data(), ys, xm,
                          thread_counts, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }

    // Uniform tiled
    for (int s = 0; s < NUM_SAMPLES; s++) {
        char en[64]; snprintf(en, sizeof(en), "uniform_tiled_s%d", s);
        printf("\n>>> 1GB %s\n", en);
        auto ubase = uniformSample(NX_BASE, NY_BASE, rng);
        auto tv = tileIndexMap(ubase.data(), NX_BASE, NY_BASE, TARGET_1GB);
        int nx_t = (int)tv.size(), K = nx_t/NX_BASE, ny_t = K*NY_BASE;
        int *ys = new int[nx_t], *xm = new int[nx_t];
        sortWithIndices(nx_t, tv.data(), ys, xm);
        sweep_all_configs(csv, en, ny_t, nx_t, tv.data(), ys, xm,
                          thread_counts, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }

    // Normal tiled
    for (int s = 0; s < NUM_SAMPLES; s++) {
        char en[64]; snprintf(en, sizeof(en), "normal_tiled_s%d", s);
        printf("\n>>> 1GB %s\n", en);
        auto nbase = normalSample(NX_BASE, NY_BASE, rng);
        auto tv = tileIndexMap(nbase.data(), NX_BASE, NY_BASE, TARGET_1GB);
        int nx_t = (int)tv.size(), K = nx_t/NX_BASE, ny_t = K*NY_BASE;
        int *ys = new int[nx_t], *xm = new int[nx_t];
        sortWithIndices(nx_t, tv.data(), ys, xm);
        sweep_all_configs(csv, en, ny_t, nx_t, tv.data(), ys, xm,
                          thread_counts, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }
    fclose(csv);

    delete[] qe_base;
    if (g_flush) numa_dealloc(g_flush, FLUSH_ELEMS);
    printf("\nDone.  CSVs: zaxpy_sweep_small_cpu.csv  zaxpy_sweep_1gb_cpu.csv\n");
    return 0;
}