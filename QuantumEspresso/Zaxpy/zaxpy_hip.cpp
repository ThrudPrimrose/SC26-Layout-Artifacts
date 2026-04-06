#include "hip/hip_runtime.h"
// zaxpy_indirect_sweep.cu
// Indirect zaxpy benchmark: y[map[i]] += x[i]
//
// All experiments use QE nat20 dimensions (nx=110273, ny=225001):
//   - "small" : original nat20 sizes
//   - "1gb"   : tiled K times with offsets so array footprint ≈ 1 GiB
//
// Three index-map distributions per scale:
//   qe       — real Quantum ESPRESSO access pattern
//   uniform  — uniformly random subset of [0, ny)
//   normal   — normally distributed subset centred on ny/2
//
// Four kernel variants (2 layouts × 2 index orders):
//   aos_scatter, aos_sorted, soa_scatter, soa_sorted
//
// Sweep: 7 thread-block sizes × 4 coarsening factors, 100 reps each.
//
// Compile:
//   nvcc -O3 -std=c++17 zaxpy_indirect_sweep.cu -o zaxpy_indirect_sweep

#include <hip/hip_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <numeric>
#include <hip/hip_complex.h>
#include <array>
#include <cmath>

#define CUDA_CHECK(call)                                     \
    do                                                       \
    {                                                        \
        hipError_t err = (call);                            \
        if (err != hipSuccess)                              \
        {                                                    \
            std::cerr << "CUDA error at " << __FILE__ << ":" \
                      << __LINE__ << " -> "                  \
                      << hipGetErrorString(err)             \
                      << std::endl;                          \
            std::exit(EXIT_FAILURE);                         \
        }                                                    \
    } while (0)

/* ================================================================ */
/*  Cache-flush helper                                               */
/* ================================================================ */
static double *d_flush = nullptr;
static constexpr size_t FLUSH_SIZE = 256ULL * 1024 * 1024;

__global__ void flush_kernel(double *buf, size_t n)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) buf[i] = 1.0;
}
static void init_flush()
{
    if (!d_flush) CUDA_CHECK(hipMalloc(&d_flush, FLUSH_SIZE));
}
static void flush_l2()
{
    size_t n = FLUSH_SIZE / sizeof(double);
    flush_kernel<<<(n + 255) / 256, 256>>>(d_flush, n);
    CUDA_CHECK(hipDeviceSynchronize());
}

/* ================================================================ */
/*  AoS kernels (hipDoubleComplex)                                    */
/* ================================================================ */

template <int C>
__global__ void kern_aos_scatter(
    int nx, const int *__restrict__ ymap,
    const hipDoubleComplex *__restrict__ x, hipDoubleComplex *__restrict__ y)
{
    int base = (blockIdx.x * blockDim.x + threadIdx.x) * C;
    #pragma unroll
    for (int c = 0; c < C; c++) {
        int g = base + c;
        if (g < nx) { int yi = ymap[g]; y[yi] = hipCadd(x[g], y[yi]); }
    }
}

template <int C>
__global__ void kern_aos_sorted(
    int nx, const int *__restrict__ ymap_s, const int *__restrict__ xmap,
    const hipDoubleComplex *__restrict__ x, hipDoubleComplex *__restrict__ y)
{
    int base = (blockIdx.x * blockDim.x + threadIdx.x) * C;
    #pragma unroll
    for (int c = 0; c < C; c++) {
        int g = base + c;
        if (g < nx) { int yi = ymap_s[g]; y[yi] = hipCadd(x[xmap[g]], y[yi]); }
    }
}

/* ================================================================ */
/*  SoA kernels (separate re / im)                                   */
/* ================================================================ */

template <int C>
__global__ void kern_soa_scatter(
    int nx, const int *__restrict__ ymap,
    const double *__restrict__ xr, const double *__restrict__ xi,
    double *__restrict__ yr, double *__restrict__ yi)
{
    int base = (blockIdx.x * blockDim.x + threadIdx.x) * C;
    #pragma unroll
    for (int c = 0; c < C; c++) {
        int g = base + c;
        if (g < nx) { int y_ = ymap[g]; yr[y_] += xr[g]; yi[y_] += xi[g]; }
    }
}

template <int C>
__global__ void kern_soa_sorted(
    int nx, const int *__restrict__ ymap_s, const int *__restrict__ xmap,
    const double *__restrict__ xr, const double *__restrict__ xi,
    double *__restrict__ yr, double *__restrict__ yi)
{
    int base = (blockIdx.x * blockDim.x + threadIdx.x) * C;
    #pragma unroll
    for (int c = 0; c < C; c++) {
        int g = base + c;
        if (g < nx) { int y_ = ymap_s[g]; int x_ = xmap[g];
                       yr[y_] += xr[x_]; yi[y_] += xi[x_]; }
    }
}

/* ================================================================ */
/*  Dispatch tables                                                  */
/* ================================================================ */
using aos_scat_t = void(*)(int, const int*, const hipDoubleComplex*, hipDoubleComplex*);
using aos_sort_t = void(*)(int, const int*, const int*, const hipDoubleComplex*, hipDoubleComplex*);
using soa_scat_t = void(*)(int, const int*, const double*, const double*, double*, double*);
using soa_sort_t = void(*)(int, const int*, const int*, const double*, const double*, double*, double*);

#define E(K,C) K<C>
aos_scat_t aos_scat_tbl[] = {nullptr, E(kern_aos_scatter,1), E(kern_aos_scatter,2), nullptr, E(kern_aos_scatter,4), nullptr,nullptr,nullptr, E(kern_aos_scatter,8)};
aos_sort_t aos_sort_tbl[] = {nullptr, E(kern_aos_sorted,1),  E(kern_aos_sorted,2),  nullptr, E(kern_aos_sorted,4),  nullptr,nullptr,nullptr, E(kern_aos_sorted,8)};
soa_scat_t soa_scat_tbl[] = {nullptr, E(kern_soa_scatter,1), E(kern_soa_scatter,2), nullptr, E(kern_soa_scatter,4), nullptr,nullptr,nullptr, E(kern_soa_scatter,8)};
soa_sort_t soa_sort_tbl[] = {nullptr, E(kern_soa_sorted,1),  E(kern_soa_sorted,2),  nullptr, E(kern_soa_sorted,4),  nullptr,nullptr,nullptr, E(kern_soa_sorted,8)};

/* ================================================================ */
/*  Verification                                                     */
/* ================================================================ */
bool verify_aos(int n, hipDoubleComplex *g, hipDoubleComplex *r, double tol=1e-6) {
    for (int i = 0; i < n; i++) {
        double d = hipCabs(hipCsub(g[i], r[i]));
        if (d / (hipCabs(r[i]) + 1e-10) > tol) { printf("AoS mismatch %d\n",i); return false; }
    }
    return true;
}
bool verify_soa(int n, double *gr, double *gi, double *rr, double *ri, double tol=1e-6) {
    for (int i = 0; i < n; i++) {
        double d = std::sqrt((gr[i]-rr[i])*(gr[i]-rr[i])+(gi[i]-ri[i])*(gi[i]-ri[i]));
        double m = std::sqrt(rr[i]*rr[i]+ri[i]*ri[i]) + 1e-10;
        if (d/m > tol) { printf("SoA mismatch %d\n",i); return false; }
    }
    return true;
}

/* ================================================================ */
/*  Core profiling — one (tpb, coarsen), all 4 variants              */
/* ================================================================ */
void profile_config(
    FILE *csv, const char *ename,
    int ny, int nx,
    int *h_ymap, int *h_ymap_s, int *h_xmap,
    int tpb, int coarsen, int iters, int warmup, unsigned seed)
{
    int epb = tpb * coarsen;
    int nblk = (nx + epb - 1) / epb;

    hipDoubleComplex *h_ya, *h_xa, *h_oa, *h_ra;
    double *h_yr, *h_yi, *h_xr, *h_xi, *h_or, *h_oi, *h_rr, *h_ri;
    hipDoubleComplex *d_ya, *d_xa;
    double *d_yr, *d_yi, *d_xr, *d_xi;
    int *d_ym, *d_yms, *d_xm;

    h_ya = (hipDoubleComplex*)malloc(ny*sizeof(hipDoubleComplex));
    h_xa = (hipDoubleComplex*)malloc(nx*sizeof(hipDoubleComplex));
    h_oa = (hipDoubleComplex*)malloc(ny*sizeof(hipDoubleComplex));
    h_ra = (hipDoubleComplex*)malloc(ny*sizeof(hipDoubleComplex));
    h_yr = (double*)malloc(ny*sizeof(double)); h_yi = (double*)malloc(ny*sizeof(double));
    h_xr = (double*)malloc(nx*sizeof(double)); h_xi = (double*)malloc(nx*sizeof(double));
    h_or = (double*)malloc(ny*sizeof(double)); h_oi = (double*)malloc(ny*sizeof(double));
    h_rr = (double*)malloc(ny*sizeof(double)); h_ri = (double*)malloc(ny*sizeof(double));

    CUDA_CHECK(hipMalloc(&d_ya,  ny*sizeof(hipDoubleComplex)));
    CUDA_CHECK(hipMalloc(&d_xa,  nx*sizeof(hipDoubleComplex)));
    CUDA_CHECK(hipMalloc(&d_yr,  ny*sizeof(double)));
    CUDA_CHECK(hipMalloc(&d_yi,  ny*sizeof(double)));
    CUDA_CHECK(hipMalloc(&d_xr,  nx*sizeof(double)));
    CUDA_CHECK(hipMalloc(&d_xi,  nx*sizeof(double)));
    CUDA_CHECK(hipMalloc(&d_ym,  nx*sizeof(int)));
    CUDA_CHECK(hipMalloc(&d_yms, nx*sizeof(int)));
    CUDA_CHECK(hipMalloc(&d_xm,  nx*sizeof(int)));

    srand(seed);
    for (int i = 0; i < ny; i++) {
        double re = (double)rand()/RAND_MAX, im = (double)rand()/RAND_MAX;
        h_ya[i] = make_hipDoubleComplex(re,im); h_ra[i] = h_ya[i];
        h_yr[i] = re; h_yi[i] = im; h_rr[i] = re; h_ri[i] = im;
    }
    for (int i = 0; i < nx; i++) {
        double re = (double)rand()/RAND_MAX, im = (double)rand()/RAND_MAX;
        h_xa[i] = make_hipDoubleComplex(re,im);
        h_xr[i] = re; h_xi[i] = im;
    }
    for (int i = 0; i < nx; i++) {
        int yi = h_ymap[i];
        h_ra[yi] = hipCadd(h_xa[i], h_ra[yi]);
        h_rr[yi] += h_xr[i]; h_ri[yi] += h_xi[i];
    }

    CUDA_CHECK(hipMemcpy(d_ym,  h_ymap,   nx*sizeof(int), hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_yms, h_ymap_s, nx*sizeof(int), hipMemcpyHostToDevice));
    CUDA_CHECK(hipMemcpy(d_xm,  h_xmap,   nx*sizeof(int), hipMemcpyHostToDevice));

    hipEvent_t e0, e1;
    CUDA_CHECK(hipEventCreate(&e0)); CUDA_CHECK(hipEventCreate(&e1));
    float ms;

    auto run = [&](const char *vname, auto launch, auto verify_fn) {
        if (!verify_fn()) { fprintf(stderr,"FAIL %s %s tpb=%d cf=%d\n",vname,ename,tpb,coarsen); return; }
        for (int w=0;w<warmup;w++){launch();CUDA_CHECK(hipDeviceSynchronize());}
        for (int r=0;r<iters;r++){
            flush_l2();
            CUDA_CHECK(hipEventRecord(e0));
            launch();
            CUDA_CHECK(hipEventRecord(e1));
            CUDA_CHECK(hipEventSynchronize(e1));
            CUDA_CHECK(hipEventElapsedTime(&ms,e0,e1));
            fprintf(csv,"%s,%s,%d,%d,%d,%d,%d,%.6f\n",vname,ename,ny,nx,tpb,coarsen,r,(double)ms);
        }
    };

    // AoS scatter
    {auto kfn=aos_scat_tbl[coarsen];
     CUDA_CHECK(hipMemcpy(d_xa,h_xa,nx*sizeof(hipDoubleComplex),hipMemcpyHostToDevice));
     auto L=[&]{kfn<<<nblk,tpb>>>(nx,d_ym,d_xa,d_ya);};
     auto V=[&]()->bool{
        CUDA_CHECK(hipMemcpy(d_ya,h_ya,ny*sizeof(hipDoubleComplex),hipMemcpyHostToDevice));
        L();CUDA_CHECK(hipDeviceSynchronize());
        CUDA_CHECK(hipMemcpy(h_oa,d_ya,ny*sizeof(hipDoubleComplex),hipMemcpyDeviceToHost));
        return verify_aos(ny,h_oa,h_ra);};
     CUDA_CHECK(hipMemcpy(d_ya,h_ya,ny*sizeof(hipDoubleComplex),hipMemcpyHostToDevice));
     run("aos_scatter",L,V);}

    // AoS sorted
    {auto kfn=aos_sort_tbl[coarsen];
     CUDA_CHECK(hipMemcpy(d_xa,h_xa,nx*sizeof(hipDoubleComplex),hipMemcpyHostToDevice));
     auto L=[&]{kfn<<<nblk,tpb>>>(nx,d_yms,d_xm,d_xa,d_ya);};
     auto V=[&]()->bool{
        CUDA_CHECK(hipMemcpy(d_ya,h_ya,ny*sizeof(hipDoubleComplex),hipMemcpyHostToDevice));
        L();CUDA_CHECK(hipDeviceSynchronize());
        CUDA_CHECK(hipMemcpy(h_oa,d_ya,ny*sizeof(hipDoubleComplex),hipMemcpyDeviceToHost));
        return verify_aos(ny,h_oa,h_ra);};
     CUDA_CHECK(hipMemcpy(d_ya,h_ya,ny*sizeof(hipDoubleComplex),hipMemcpyHostToDevice));
     run("aos_sorted",L,V);}

    // SoA scatter
    {auto kfn=soa_scat_tbl[coarsen];
     CUDA_CHECK(hipMemcpy(d_xr,h_xr,nx*sizeof(double),hipMemcpyHostToDevice));
     CUDA_CHECK(hipMemcpy(d_xi,h_xi,nx*sizeof(double),hipMemcpyHostToDevice));
     auto L=[&]{kfn<<<nblk,tpb>>>(nx,d_ym,d_xr,d_xi,d_yr,d_yi);};
     auto V=[&]()->bool{
        CUDA_CHECK(hipMemcpy(d_yr,h_yr,ny*sizeof(double),hipMemcpyHostToDevice));
        CUDA_CHECK(hipMemcpy(d_yi,h_yi,ny*sizeof(double),hipMemcpyHostToDevice));
        L();CUDA_CHECK(hipDeviceSynchronize());
        CUDA_CHECK(hipMemcpy(h_or,d_yr,ny*sizeof(double),hipMemcpyDeviceToHost));
        CUDA_CHECK(hipMemcpy(h_oi,d_yi,ny*sizeof(double),hipMemcpyDeviceToHost));
        return verify_soa(ny,h_or,h_oi,h_rr,h_ri);};
     CUDA_CHECK(hipMemcpy(d_yr,h_yr,ny*sizeof(double),hipMemcpyHostToDevice));
     CUDA_CHECK(hipMemcpy(d_yi,h_yi,ny*sizeof(double),hipMemcpyHostToDevice));
     run("soa_scatter",L,V);}

    // SoA sorted
    {auto kfn=soa_sort_tbl[coarsen];
     CUDA_CHECK(hipMemcpy(d_xr,h_xr,nx*sizeof(double),hipMemcpyHostToDevice));
     CUDA_CHECK(hipMemcpy(d_xi,h_xi,nx*sizeof(double),hipMemcpyHostToDevice));
     auto L=[&]{kfn<<<nblk,tpb>>>(nx,d_yms,d_xm,d_xr,d_xi,d_yr,d_yi);};
     auto V=[&]()->bool{
        CUDA_CHECK(hipMemcpy(d_yr,h_yr,ny*sizeof(double),hipMemcpyHostToDevice));
        CUDA_CHECK(hipMemcpy(d_yi,h_yi,ny*sizeof(double),hipMemcpyHostToDevice));
        L();CUDA_CHECK(hipDeviceSynchronize());
        CUDA_CHECK(hipMemcpy(h_or,d_yr,ny*sizeof(double),hipMemcpyDeviceToHost));
        CUDA_CHECK(hipMemcpy(h_oi,d_yi,ny*sizeof(double),hipMemcpyDeviceToHost));
        return verify_soa(ny,h_or,h_oi,h_rr,h_ri);};
     CUDA_CHECK(hipMemcpy(d_yr,h_yr,ny*sizeof(double),hipMemcpyHostToDevice));
     CUDA_CHECK(hipMemcpy(d_yi,h_yi,ny*sizeof(double),hipMemcpyHostToDevice));
     run("soa_sorted",L,V);}

    CUDA_CHECK(hipEventDestroy(e0)); CUDA_CHECK(hipEventDestroy(e1));
    CUDA_CHECK(hipFree(d_ya)); CUDA_CHECK(hipFree(d_xa));
    CUDA_CHECK(hipFree(d_yr)); CUDA_CHECK(hipFree(d_yi));
    CUDA_CHECK(hipFree(d_xr)); CUDA_CHECK(hipFree(d_xi));
    CUDA_CHECK(hipFree(d_ym)); CUDA_CHECK(hipFree(d_yms)); CUDA_CHECK(hipFree(d_xm));
    free(h_ya);free(h_xa);free(h_oa);free(h_ra);
    free(h_yr);free(h_yi);free(h_xr);free(h_xi);
    free(h_or);free(h_oi);free(h_rr);free(h_ri);
}

/* ================================================================ */
/*  Helpers                                                          */
/* ================================================================ */

void sortWithIndices(int N, const int *in, int *sv, int *si)
{
    std::iota(si, si+N, 0);
    std::sort(si, si+N, [&in](int a, int b){return in[a]<in[b];});
    for (int i=0;i<N;i++) sv[i]=in[si[i]];
}

std::vector<int> uniformSample(int nx, int ny, std::mt19937 &rng)
{
    std::vector<int> pool(ny);
    std::iota(pool.begin(), pool.end(), 0);
    for (int i=0;i<nx;i++){
        std::uniform_int_distribution<int> d(i,ny-1);
        std::swap(pool[i], pool[d(rng)]);
    }
    return std::vector<int>(pool.begin(), pool.begin()+nx);
}

std::vector<int> normalSample(int nx, int ny, std::mt19937 &rng)
{
    double mu=(ny-1)/2.0, sigma=ny/6.0;
    std::normal_distribution<double> dist(mu,sigma);
    std::unordered_set<int> seen;
    std::vector<int> res; res.reserve(nx);
    while ((int)res.size()<nx){
        int s=(int)std::round(dist(rng));
        if(s<0||s>=ny) continue;
        if(!seen.insert(s).second) continue;
        res.push_back(s);
    }
    return res;
}

int *readBinaryFile(const std::string &fp, int &N)
{
    std::ifstream f(fp, std::ios::binary|std::ios::ate);
    if(!f.is_open()) throw std::runtime_error("Cannot open "+fp);
    auto sz=f.tellg();
    if(sz%sizeof(int)!=0) throw std::runtime_error("Bad size");
    N=sz/sizeof(int); f.seekg(0);
    int *d=new int[N];
    if(!f.read((char*)d,sz)){delete[]d; throw std::runtime_error("Read fail");}
    return d;
}

/**
 * Tile an index map K times with offsets.
 * Tile k: out[k*nx_base + i] = base[i] + k*ny_base
 */
std::vector<int> tileIndexMap(const int *base, int nx_base, int ny_base,
                              size_t target_bytes)
{
    size_t bpt = (size_t)(ny_base + nx_base) * 16;
    int K = std::max(1, (int)(target_bytes / bpt));
    int nx_t = K * nx_base;
    std::vector<int> out(nx_t);
    for (int k=0; k<K; k++){
        int yo = k*ny_base, xo = k*nx_base;
        for (int i=0; i<nx_base; i++) out[xo+i] = base[i] + yo;
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

void sweep_all_configs(
    FILE *csv, const char *ename,
    int ny, int nx,
    int *ymap, int *ymap_s, int *xmap,
    const std::array<int,7> &tpbs,
    const std::array<int,4> &cfs,
    int iters, int warmup)
{
    for (int tpb : tpbs)
        for (int cf : cfs) {
            printf("  %s tpb=%d cf=%d\n", ename, tpb, cf);
            profile_config(csv, ename, ny, nx, ymap, ymap_s, xmap,
                           tpb, cf, iters, warmup, 0);
            fflush(csv);
        }
}

/* ================================================================ */
/*  Main                                                             */
/* ================================================================ */

int main(void)
{
    constexpr int ITERS  = 100;
    constexpr int WARMUP = 10;
    constexpr int NUM_SAMPLES = 5;
    constexpr size_t TARGET_1GB = 1ULL << 30;

    // nat20 base dimensions
    constexpr int NX_BASE = 110273;
    constexpr int NY_BASE = 225001;
    const char *QE_BIN   = "./bin/vexx_k_gpu__dfftt__nl_nat20.bin";

    std::array<int,7> tpbs = {32, 64, 128, 256, 384, 512, 1024};
    std::array<int,4> cfs  = {1, 2, 4, 8};

    std::mt19937 rng(42);
    init_flush();

    // ── Load QE base index map ──────────────────────────────────────────
    int N = 0;
    int *qe_base = readBinaryFile(QE_BIN, N);
    assert(N == NX_BASE);

    // ================================================================
    //  SMALL scale (nat20 original sizes)
    // ================================================================
    FILE *csv = fopen("zaxpy_sweep_small.csv", "w");
    fprintf(csv, "%s", CSV_HDR);

    // -- QE --
    {
        printf("\n>>> SMALL qe  ny=%d nx=%d\n", NY_BASE, NX_BASE);
        int *ys = new int[NX_BASE], *xm = new int[NX_BASE];
        sortWithIndices(NX_BASE, qe_base, ys, xm);
        sweep_all_configs(csv, "qe", NY_BASE, NX_BASE, qe_base, ys, xm,
                          tpbs, cfs, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }

    // -- Uniform (NUM_SAMPLES samples) --
    for (int s = 0; s < NUM_SAMPLES; s++)
    {
        char en[64]; snprintf(en, sizeof(en), "uniform_s%d", s);
        printf("\n>>> SMALL %s  ny=%d nx=%d\n", en, NY_BASE, NX_BASE);
        auto yv = uniformSample(NX_BASE, NY_BASE, rng);
        int *ys = new int[NX_BASE], *xm = new int[NX_BASE];
        sortWithIndices(NX_BASE, yv.data(), ys, xm);
        sweep_all_configs(csv, en, NY_BASE, NX_BASE, yv.data(), ys, xm,
                          tpbs, cfs, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }

    // -- Normal (NUM_SAMPLES samples) --
    for (int s = 0; s < NUM_SAMPLES; s++)
    {
        char en[64]; snprintf(en, sizeof(en), "normal_s%d", s);
        printf("\n>>> SMALL %s  ny=%d nx=%d\n", en, NY_BASE, NX_BASE);
        auto yv = normalSample(NX_BASE, NY_BASE, rng);
        int *ys = new int[NX_BASE], *xm = new int[NX_BASE];
        sortWithIndices(NX_BASE, yv.data(), ys, xm);
        sweep_all_configs(csv, en, NY_BASE, NX_BASE, yv.data(), ys, xm,
                          tpbs, cfs, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }
    fclose(csv);

    // ================================================================
    //  1 GB scale (tiled with offsets)
    // ================================================================
    csv = fopen("zaxpy_sweep_1gb.csv", "w");
    fprintf(csv, "%s", CSV_HDR);

    // -- QE tiled --
    {
        printf("\n>>> 1GB qe_tiled\n");
        auto tv = tileIndexMap(qe_base, NX_BASE, NY_BASE, TARGET_1GB);
        int nx_t = (int)tv.size();
        int K = nx_t / NX_BASE;
        int ny_t = K * NY_BASE;
        int *ys = new int[nx_t], *xm = new int[nx_t];
        sortWithIndices(nx_t, tv.data(), ys, xm);
        sweep_all_configs(csv, "qe_tiled", ny_t, nx_t, tv.data(), ys, xm,
                          tpbs, cfs, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }

    // -- Uniform tiled --
    for (int s = 0; s < NUM_SAMPLES; s++)
    {
        char en[64]; snprintf(en, sizeof(en), "uniform_tiled_s%d", s);
        printf("\n>>> 1GB %s\n", en);
        auto ubase = uniformSample(NX_BASE, NY_BASE, rng);
        auto tv = tileIndexMap(ubase.data(), NX_BASE, NY_BASE, TARGET_1GB);
        int nx_t = (int)tv.size();
        int K = nx_t / NX_BASE;
        int ny_t = K * NY_BASE;
        int *ys = new int[nx_t], *xm = new int[nx_t];
        sortWithIndices(nx_t, tv.data(), ys, xm);
        sweep_all_configs(csv, en, ny_t, nx_t, tv.data(), ys, xm,
                          tpbs, cfs, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }

    // -- Normal tiled --
    for (int s = 0; s < NUM_SAMPLES; s++)
    {
        char en[64]; snprintf(en, sizeof(en), "normal_tiled_s%d", s);
        printf("\n>>> 1GB %s\n", en);
        auto nbase = normalSample(NX_BASE, NY_BASE, rng);
        auto tv = tileIndexMap(nbase.data(), NX_BASE, NY_BASE, TARGET_1GB);
        int nx_t = (int)tv.size();
        int K = nx_t / NX_BASE;
        int ny_t = K * NY_BASE;
        int *ys = new int[nx_t], *xm = new int[nx_t];
        sortWithIndices(nx_t, tv.data(), ys, xm);
        sweep_all_configs(csv, en, ny_t, nx_t, tv.data(), ys, xm,
                          tpbs, cfs, ITERS, WARMUP);
        delete[] ys; delete[] xm;
    }
    fclose(csv);

    delete[] qe_base;
    if (d_flush) hipFree(d_flush);
    printf("\nDone.  CSVs: zaxpy_sweep_small.csv  zaxpy_sweep_1gb.csv\n");
    return 0;
}
