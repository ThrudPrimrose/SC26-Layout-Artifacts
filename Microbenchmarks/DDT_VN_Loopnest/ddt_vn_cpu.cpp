/*
 * bench_ddt_vn_cpu_sweep.cpp -- CPU benchmark for ddt_vn_apc_pc stencil
 *
 * Compares three layouts:
 *   SoA     — baseline, 15 distinct array streams
 *   AoS     — grouped (grpA/grpC/grpD), 10 streams
 *   AoSoA16 — grouped with V=16 vector tiling, 10 streams
 *
 * V=16 means each tile holds 16 consecutive elements per field:
 *   16 doubles × 8 B = 128 B = 2 cache lines per field per tile.
 *   grpA tile: 5 fields × 128 B = 640 B (10 cache lines, fits L1).
 *   Vectoriser sees contiguous runs of 16 for each field.
 *
 * NUMA-aware first-touch allocation via mmap.
 * Cache flush: 32768² Jacobi stencil (8 GB per buffer), 3 sweeps.
 * 100 repetitions with per-run CSV output.
 *
 * Compile:
 *   g++ -O3 -fopenmp -march=native -std=c++17 bench_ddt_vn_cpu_sweep.cpp -o bench_ddt_vn_cpu_sweep
 *
 * Run:
 *   ./bench_ddt_vn_cpu_sweep          (defaults: nlev=65,90)
 *   ./bench_ddt_vn_cpu_sweep 90       (single nlev)
 */

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>
#include <omp.h>
#include <sys/mman.h>

/* ================================================================ */
/*  Constants                                                        */
/* ================================================================ */
static constexpr int NPROMA = 81920;
static constexpr int NRUNS  = 100;
static constexpr int WARMUP = 10;

/* Default nlev sweep (overridable via argv) */
static constexpr int DEFAULT_NLEVS[] = { 65, 90 };
static constexpr int DEFAULT_N_NLEVS = sizeof(DEFAULT_NLEVS)/sizeof(DEFAULT_NLEVS[0]);

/* AoSoA vector width for CPU */
static constexpr int VW = 16;

/* Group-A field indices (edge wind state, 3D + nlevp1) */
static constexpr int A_VN    = 0;
static constexpr int A_VT    = 1;
static constexpr int A_VN_IE = 2;
static constexpr int A_VT_IE = 3;
static constexpr int A_EKIN  = 4;
static constexpr int NA      = 5;

/* Group-C field indices (edge geometry, 2D) */
static constexpr int C_IDUAL = 0;
static constexpr int C_IPRIM = 1;
static constexpr int C_TANG  = 2;
static constexpr int C_FE    = 3;
static constexpr int NC      = 4;

/* Group-D field indices (cell gathered pair, 3D indirect) */
static constexpr int D_EKINH = 0;
static constexpr int D_WCON  = 1;
static constexpr int ND      = 2;

/* ================================================================ */
/*  NUMA allocation via mmap + first-touch                           */
/* ================================================================ */
template<typename T>
static T* numa_alloc(size_t n) {
    size_t bytes = n * sizeof(T);
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); exit(1); }
    return static_cast<T*>(p);
}

template<typename T>
static void numa_free(T *p, size_t n) {
    munmap(p, n * sizeof(T));
}

/* First-touch copy for 2D arrays: distribute along jk (= outer parallel dim) */
static double* numa_redistribute_2d(const double *src, int N, int K) {
    size_t n = (size_t)N * K;
    double *dst = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < K; jk++)
        memcpy(dst + (size_t)jk * N, src + (size_t)jk * N, N * sizeof(double));
    return dst;
}

/* First-touch copy for 1D arrays: stripe across threads */
static double* numa_redistribute_1d(const double *src, size_t n) {
    double *dst = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) dst[i] = src[i];
    return dst;
}

static int* numa_redistribute_idx(const int *src, int N) {
    size_t n = (size_t)N * 2;
    int *dst = numa_alloc<int>(n);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) dst[i] = src[i];
    return dst;
}

/* First-touch for AoS grpA: distribute along jk dimension */
static double* numa_redistribute_grpA(const double *src, int N, int nlevp1) {
    size_t n = (size_t)NA * N * nlevp1;
    double *dst = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlevp1; jk++) {
        size_t off = (size_t)NA * N * jk;
        memcpy(dst + off, src + off, (size_t)NA * N * sizeof(double));
    }
    return dst;
}

static double* numa_redistribute_grpC(const double *src, int N) {
    size_t n = (size_t)NC * N;
    double *dst = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) dst[i] = src[i];
    return dst;
}

static double* numa_redistribute_grpD(const double *src, int N_c, int nlev) {
    size_t n = (size_t)ND * N_c * nlev;
    double *dst = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++) {
        size_t off = (size_t)ND * N_c * jk;
        memcpy(dst + off, src + off, (size_t)ND * N_c * sizeof(double));
    }
    return dst;
}

/* First-touch for AoSoA grpA: distribute along jk */
static double* numa_redistribute_aosoa_grpA(const double *src, int N, int nlevp1) {
    int tiles = (N + VW - 1) / VW;
    size_t n = (size_t)tiles * NA * VW * nlevp1;
    double *dst = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlevp1; jk++) {
        size_t off = (size_t)jk * tiles * NA * VW;
        memcpy(dst + off, src + off, (size_t)tiles * NA * VW * sizeof(double));
    }
    return dst;
}

static double* numa_redistribute_aosoa_grpC(const double *src, int N) {
    int tiles = (N + VW - 1) / VW;
    size_t n = (size_t)tiles * NC * VW;
    double *dst = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) dst[i] = src[i];
    return dst;
}

static double* numa_redistribute_aosoa_grpD(const double *src, int N_c, int nlev) {
    int tiles = (N_c + VW - 1) / VW;
    size_t n = (size_t)tiles * ND * VW * nlev;
    double *dst = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++) {
        size_t off = (size_t)jk * tiles * ND * VW;
        memcpy(dst + off, src + off, (size_t)tiles * ND * VW * sizeof(double));
    }
    return dst;
}

/* ================================================================ */
/*  RNG                                                              */
/* ================================================================ */
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
static void fill_rand(double *a, size_t n, unsigned seed) {
    for (size_t i = 0; i < n; i++) {
        uint64_t h = splitmix64((uint64_t)seed * 2654435761ULL + i);
        a[i] = (double)(int64_t)(h & 0xFFFFF) / 100000.0 - 5.0;
    }
}

/* ================================================================ */
/*  Index helpers                                                    */
/* ================================================================ */
static inline int I2(int je, int jk, int N) { return je + jk * N; }
static inline int IX(int je, int n, int N)  { return je + n * N; }
static inline int IG(int je, int n, int N)  { return je + n * N; }

/* AoS group indices */
static inline int IA(int f, int je, int jk, int N) { return f + NA * (je + N * jk); }
static inline int IC_grp(int f, int je)              { return f + NC * je; }
static inline int ID(int f, int cell, int jk, int Nc) { return f + ND * (cell + Nc * jk); }

/* AoSoA group indices (V = VW = 16) */
static inline int IA_aosoa(int f, int je, int jk, int N) {
    int tiles = (N + VW - 1) / VW;
    return jk * tiles * NA * VW + (je / VW) * NA * VW + f * VW + (je % VW);
}
static inline int IC_aosoa(int f, int je) {
    return (je / VW) * NC * VW + f * VW + (je % VW);
}
static inline int ID_aosoa(int f, int cell, int jk, int Nc) {
    int tiles = (Nc + VW - 1) / VW;
    return jk * tiles * ND * VW + (cell / VW) * ND * VW + f * VW + (cell % VW);
}

/* AoSoA allocation sizes */
static inline size_t aosoa_grpA_size(int N, int nlevp1) {
    return (size_t)((N + VW - 1) / VW) * NA * VW * nlevp1;
}
static inline size_t aosoa_grpC_size(int N) {
    return (size_t)((N + VW - 1) / VW) * NC * VW;
}
static inline size_t aosoa_grpD_size(int Nc, int nlev) {
    return (size_t)((Nc + VW - 1) / VW) * ND * VW * nlev;
}

/* ================================================================ */
/*  Stencil body macros                                              */
/* ================================================================ */

/* SoA baseline — 15 streams */
#define SOA_BODY()                                                              \
    int c2 = I2(je, jk, N);                                                    \
    int ci0 = cell_idx[IX(je,0,N)], ci1 = cell_idx[IX(je,1,N)];                \
    int vi0 = vert_idx[IX(je,0,N)], vi1 = vert_idx[IX(je,1,N)];                \
    double cg0 = coeff_gradekin[IG(je,0,N)];                                   \
    double cg1 = coeff_gradekin[IG(je,1,N)];                                   \
    double cl0 = c_lin_e[IG(je,0,N)], cl1 = c_lin_e[IG(je,1,N)];              \
    out[c2] = -(                                                               \
        z_kin_hor_e[c2] * (cg0 - cg1)                                         \
        + cg1 * z_ekinh[I2(ci1, jk, N_c)]                                     \
        - cg0 * z_ekinh[I2(ci0, jk, N_c)]                                     \
        + vt[c2] * (f_e[je] + 0.5 *                                           \
            (zeta[I2(vi0, jk, N_v)] + zeta[I2(vi1, jk, N_v)]))                \
        + (cl0 * z_w_con_c_full[I2(ci0, jk, N_c)]                             \
         + cl1 * z_w_con_c_full[I2(ci1, jk, N_c)])                            \
          * (vn_ie[I2(je, jk, N)] - vn_ie[I2(je, jk+1, N)])                   \
          / ddqz_z_full_e[c2]                                                  \
    );

/* AoS grouped — 10 streams */
#define GRP_BODY()                                                              \
    int ci0 = cell_idx[IX(je,0,N)], ci1 = cell_idx[IX(je,1,N)];                \
    int vi0 = vert_idx[IX(je,0,N)], vi1 = vert_idx[IX(je,1,N)];                \
    double cg0 = coeff_gradekin[IG(je,0,N)];                                   \
    double cg1 = coeff_gradekin[IG(je,1,N)];                                   \
    double cl0 = c_lin_e[IG(je,0,N)], cl1 = c_lin_e[IG(je,1,N)];              \
    double ekin_e = grpA[IA(A_EKIN, je, jk, N)];                               \
    double vt_e   = grpA[IA(A_VT,   je, jk, N)];                               \
    double vnie_k = grpA[IA(A_VN_IE,je, jk, N)];                               \
    double vnie_k1= grpA[IA(A_VN_IE,je, jk+1, N)];                             \
    double fe_e   = grpC[IC_grp(C_FE, je)];                                    \
    double ekinh0 = grpD[ID(D_EKINH, ci0, jk, N_c)];                           \
    double ekinh1 = grpD[ID(D_EKINH, ci1, jk, N_c)];                           \
    double wcon0  = grpD[ID(D_WCON,  ci0, jk, N_c)];                           \
    double wcon1  = grpD[ID(D_WCON,  ci1, jk, N_c)];                           \
    out[I2(je, jk, N)] = -(                                                    \
        ekin_e * (cg0 - cg1)                                                   \
        + cg1 * ekinh1 - cg0 * ekinh0                                         \
        + vt_e * (fe_e + 0.5 *                                                \
            (zeta[I2(vi0, jk, N_v)] + zeta[I2(vi1, jk, N_v)]))                \
        + (cl0 * wcon0 + cl1 * wcon1)                                          \
          * (vnie_k - vnie_k1) / ddqz_z_full_e[I2(je, jk, N)]                \
    );

/* AoSoA-16 grouped — 10 streams, vector-tiled */
#define AOSOA_BODY()                                                            \
    int ci0 = cell_idx[IX(je,0,N)], ci1 = cell_idx[IX(je,1,N)];                \
    int vi0 = vert_idx[IX(je,0,N)], vi1 = vert_idx[IX(je,1,N)];                \
    double cg0 = coeff_gradekin[IG(je,0,N)];                                   \
    double cg1 = coeff_gradekin[IG(je,1,N)];                                   \
    double cl0 = c_lin_e[IG(je,0,N)], cl1 = c_lin_e[IG(je,1,N)];              \
    double ekin_e = grpA[IA_aosoa(A_EKIN, je, jk, N)];                         \
    double vt_e   = grpA[IA_aosoa(A_VT,   je, jk, N)];                         \
    double vnie_k = grpA[IA_aosoa(A_VN_IE,je, jk, N)];                         \
    double vnie_k1= grpA[IA_aosoa(A_VN_IE,je, jk+1, N)];                       \
    double fe_e   = grpC[IC_aosoa(C_FE, je)];                                  \
    double ekinh0 = grpD[ID_aosoa(D_EKINH, ci0, jk, N_c)];                     \
    double ekinh1 = grpD[ID_aosoa(D_EKINH, ci1, jk, N_c)];                     \
    double wcon0  = grpD[ID_aosoa(D_WCON,  ci0, jk, N_c)];                     \
    double wcon1  = grpD[ID_aosoa(D_WCON,  ci1, jk, N_c)];                     \
    out[I2(je, jk, N)] = -(                                                    \
        ekin_e * (cg0 - cg1)                                                   \
        + cg1 * ekinh1 - cg0 * ekinh0                                         \
        + vt_e * (fe_e + 0.5 *                                                \
            (zeta[I2(vi0, jk, N_v)] + zeta[I2(vi1, jk, N_v)]))                \
        + (cl0 * wcon0 + cl1 * wcon1)                                          \
          * (vnie_k - vnie_k1) / ddqz_z_full_e[I2(je, jk, N)]                \
    );

/* ================================================================ */
/*  Connectivity distributions                                       */
/* ================================================================ */
enum CellDist { UNIFORM=0, NORMAL1=1, NORMAL4=2, SEQUENTIAL=3 };
static const char *dist_names[] = {"uniform","normal_var1","normal_var4","sequential"};

static void gen_connectivity(int *L, int N, int N_target,
                             CellDist dist, std::mt19937 &rng) {
    switch (dist) {
    case UNIFORM: {
        std::uniform_int_distribution<int> u(0, N_target - 1);
        for (int i = 0; i < N; i++) { L[i*2] = u(rng); L[i*2+1] = u(rng); }
        break;
    }
    case NORMAL1: {
        std::normal_distribution<double> nd(0, 1);
        for (int i = 0; i < N; i++) {
            L[i*2]   = ((i+1+(int)std::round(nd(rng)))%N_target+N_target)%N_target;
            L[i*2+1] = ((i-1+(int)std::round(nd(rng)))%N_target+N_target)%N_target;
        }
        break;
    }
    case NORMAL4: {
        std::normal_distribution<double> nd(0, 2);
        for (int i = 0; i < N; i++) {
            L[i*2]   = ((i+1+(int)std::round(nd(rng)))%N_target+N_target)%N_target;
            L[i*2+1] = ((i-1+(int)std::round(nd(rng)))%N_target+N_target)%N_target;
        }
        break;
    }
    case SEQUENTIAL:
        for (int i = 0; i < N; i++) {
            L[i*2]   = (i+1) % N_target;
            L[i*2+1] = (i+1) % N_target;
        }
        break;
    }
}

/* ================================================================ */
/*  Cache flush — 32768² Jacobi, 3 sweeps (~8 GB per buffer)        */
/* ================================================================ */
static constexpr int FN = 8192 * 4;   /* 32768 */
static constexpr int FS = 3;
static double *fb0 = nullptr, *fb1 = nullptr;

static void flush_init() {
    size_t n = (size_t)FN * FN;
    fb0 = numa_alloc<double>(n);
    fb1 = numa_alloc<double>(n);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) {
        uint64_t h = splitmix64(12345ULL + i);
        fb0[i] = (double)(h >> 11) / (double)(1ULL << 53);
        fb1[i] = fb0[i];
    }
}
static void flush() {
    double *A = fb0, *B = fb1;
    for (int s = 0; s < FS; s++) {
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < FN - 1; i++)
            for (int j = 1; j < FN - 1; j++)
                B[i*FN+j] = 0.25 * (A[(i-1)*FN+j] + A[(i+1)*FN+j]
                                   + A[i*FN+(j-1)] + A[i*FN+(j+1)]);
        std::swap(A, B);
    }
    printf("  [flush] A[0]=%.6e\n", A[0]);
}
static void flush_destroy() {
    size_t n = (size_t)FN * FN;
    numa_free(fb0, n); numa_free(fb1, n);
}

/* ================================================================ */
/*  Verification                                                     */
/* ================================================================ */
static bool verify(const double *got, const double *ref, size_t n,
                   int &nfail, double &maxrel) {
    nfail = 0; maxrel = 0;
    for (size_t i = 0; i < n; i++) {
        double d = std::abs(got[i] - ref[i]);
        double dn = std::max(std::abs(ref[i]), 1e-300);
        double r = d / dn;
        if (r > maxrel) maxrel = r;
        if (d > 1e-12 + 1e-8 * std::abs(ref[i])) nfail++;
    }
    return nfail == 0;
}

/* ================================================================ */
/*  SoA data with NUMA redistribution                                */
/* ================================================================ */
struct SoADataNuma {
    int N, N_c, N_v, nlev;
    size_t sz_e, sz_c, sz_v;

    /* NUMA-distributed pointers (used in kernels) */
    double *vt, *vn_ie, *z_kin_hor_e, *ddqz_z_full_e;
    double *coeff_gradekin, *c_lin_e, *f_e;
    double *z_ekinh, *z_w_con_c_full, *zeta;
    int    *cell_idx, *vert_idx;
    double *out;

    /* Staging buffers (heap, for initial fill + packing) */
    double *_vt, *_vn_ie, *_ekin, *_ddqz;
    double *_cg, *_cl, *_fe;
    double *_ekinh, *_wcon, *_zeta;
    int    *_cidx, *_vidx;

    void alloc(int N_, int Nc, int Nv, int nlev_) {
        N=N_; N_c=Nc; N_v=Nv; nlev=nlev_;
        sz_e=(size_t)N*nlev; sz_c=(size_t)Nc*nlev; sz_v=(size_t)Nv*nlev;
        int nlevp1 = nlev + 1;
        _vt   = new double[sz_e];    _vn_ie = new double[(size_t)N*nlevp1];
        _ekin = new double[sz_e];    _ddqz  = new double[sz_e];
        _cg   = new double[N*2];     _cl    = new double[N*2];
        _fe   = new double[N];
        _ekinh= new double[sz_c];    _wcon  = new double[sz_c];
        _zeta = new double[sz_v];
        _cidx = new int[N*2];        _vidx  = new int[N*2];
    }
    void fill() {
        fill_rand(_vt, sz_e, 101);   fill_rand(_vn_ie, (size_t)N*(nlev+1), 102);
        fill_rand(_ekin, sz_e, 103);  fill_rand(_ddqz, sz_e, 104);
        fill_rand(_cg, N*2, 105);     fill_rand(_cl, N*2, 106);
        fill_rand(_fe, N, 107);
        fill_rand(_ekinh, sz_c, 108); fill_rand(_wcon, sz_c, 109);
        fill_rand(_zeta, sz_v, 110);
        for (size_t i = 0; i < sz_e; i++)
            if (std::abs(_ddqz[i]) < 1e-10) _ddqz[i] = 1.0;
    }
    void set_connectivity(const int *cell_log, const int *vert_log) {
        memcpy(_cidx, cell_log, N*2*sizeof(int));
        memcpy(_vidx, vert_log, N*2*sizeof(int));
    }
    void numa_distribute() {
        vt             = numa_redistribute_2d(_vt, N, nlev);
        vn_ie          = numa_redistribute_2d(_vn_ie, N, nlev+1);
        z_kin_hor_e    = numa_redistribute_2d(_ekin, N, nlev);
        ddqz_z_full_e  = numa_redistribute_2d(_ddqz, N, nlev);
        z_ekinh        = numa_redistribute_2d(_ekinh, N_c, nlev);
        z_w_con_c_full = numa_redistribute_2d(_wcon, N_c, nlev);
        zeta           = numa_redistribute_2d(_zeta, N_v, nlev);
        coeff_gradekin = numa_redistribute_1d(_cg, N*2);
        c_lin_e        = numa_redistribute_1d(_cl, N*2);
        f_e            = numa_redistribute_1d(_fe, N);
        cell_idx       = numa_redistribute_idx(_cidx, N);
        vert_idx       = numa_redistribute_idx(_vidx, N);
        out            = numa_alloc<double>(sz_e);
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < sz_e; i++) out[i] = 0.0;
    }
    void numa_free_all() {
        numa_free(vt, sz_e); numa_free(vn_ie, (size_t)N*(nlev+1));
        numa_free(z_kin_hor_e, sz_e); numa_free(ddqz_z_full_e, sz_e);
        numa_free(z_ekinh, sz_c); numa_free(z_w_con_c_full, sz_c);
        numa_free(zeta, sz_v);
        numa_free(coeff_gradekin, (size_t)N*2);
        numa_free(c_lin_e, (size_t)N*2);
        numa_free(f_e, (size_t)N);
        numa_free(cell_idx, (size_t)N*2);
        numa_free(vert_idx, (size_t)N*2);
        numa_free(out, sz_e);
    }
    void free_staging() {
        delete[] _vt; delete[] _vn_ie; delete[] _ekin; delete[] _ddqz;
        delete[] _cg; delete[] _cl; delete[] _fe;
        delete[] _ekinh; delete[] _wcon; delete[] _zeta;
        delete[] _cidx; delete[] _vidx;
    }
};

/* ================================================================ */
/*  AoS grouped data with NUMA redistribution                        */
/* ================================================================ */
struct GrpDataNuma {
    int N, N_c, N_v, nlev;
    size_t sz_e;

    double *grpA, *grpC, *grpD;
    double *ddqz_z_full_e, *coeff_gradekin, *c_lin_e, *zeta;
    int    *cell_idx, *vert_idx;
    double *out;

    void from_soa_and_distribute(const SoADataNuma &s) {
        N=s.N; N_c=s.N_c; N_v=s.N_v; nlev=s.nlev;
        int nlevp1 = nlev + 1;
        sz_e = (size_t)N * nlev;

        /* Pack grpA on heap */
        size_t szA = (size_t)NA * N * nlevp1;
        double *tmpA = new double[szA];
        memset(tmpA, 0, szA * sizeof(double));
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) {
                tmpA[IA(A_VT,    je, jk, N)] = s._vt[I2(je, jk, N)];
                tmpA[IA(A_VN_IE, je, jk, N)] = s._vn_ie[I2(je, jk, N)];
                tmpA[IA(A_EKIN,  je, jk, N)] = s._ekin[I2(je, jk, N)];
            }
        for (int je = 0; je < N; je++)
            tmpA[IA(A_VN_IE, je, nlev, N)] = s._vn_ie[I2(je, nlev, N)];

        /* Pack grpC */
        double *tmpC = new double[(size_t)NC * N];
        memset(tmpC, 0, NC * N * sizeof(double));
        for (int je = 0; je < N; je++)
            tmpC[IC_grp(C_FE, je)] = s._fe[je];

        /* Pack grpD */
        size_t szD = (size_t)ND * N_c * nlev;
        double *tmpD = new double[szD];
        for (int jk = 0; jk < nlev; jk++)
            for (int ic = 0; ic < N_c; ic++) {
                tmpD[ID(D_EKINH, ic, jk, N_c)] = s._ekinh[I2(ic, jk, N_c)];
                tmpD[ID(D_WCON,  ic, jk, N_c)] = s._wcon[I2(ic, jk, N_c)];
            }

        /* NUMA redistribute */
        grpA           = numa_redistribute_grpA(tmpA, N, nlevp1);
        grpC           = numa_redistribute_grpC(tmpC, N);
        grpD           = numa_redistribute_grpD(tmpD, N_c, nlev);
        ddqz_z_full_e  = numa_redistribute_2d(s._ddqz, N, nlev);
        coeff_gradekin = numa_redistribute_1d(s._cg, N * 2);
        c_lin_e        = numa_redistribute_1d(s._cl, N * 2);
        zeta           = numa_redistribute_2d(s._zeta, N_v, nlev);
        cell_idx       = numa_redistribute_idx(s._cidx, N);
        vert_idx       = numa_redistribute_idx(s._vidx, N);
        out            = numa_alloc<double>(sz_e);
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < sz_e; i++) out[i] = 0.0;

        delete[] tmpA; delete[] tmpC; delete[] tmpD;
    }
    void numa_free_all() {
        numa_free(grpA, (size_t)NA * N * (nlev + 1));
        numa_free(grpC, (size_t)NC * N);
        numa_free(grpD, (size_t)ND * N_c * nlev);
        numa_free(ddqz_z_full_e, sz_e);
        numa_free(coeff_gradekin, (size_t)N * 2);
        numa_free(c_lin_e, (size_t)N * 2);
        numa_free(zeta, (size_t)N_v * nlev);
        numa_free(cell_idx, (size_t)N * 2);
        numa_free(vert_idx, (size_t)N * 2);
        numa_free(out, sz_e);
    }
};

/* ================================================================ */
/*  AoSoA-16 data with NUMA redistribution                          */
/* ================================================================ */
struct AoSoADataNuma {
    int N, N_c, N_v, nlev;
    size_t sz_e;
    size_t szA, szC, szD;

    double *grpA, *grpC, *grpD;
    double *ddqz_z_full_e, *coeff_gradekin, *c_lin_e, *zeta;
    int    *cell_idx, *vert_idx;
    double *out;

    void from_soa_and_distribute(const SoADataNuma &s) {
        N=s.N; N_c=s.N_c; N_v=s.N_v; nlev=s.nlev;
        int nlevp1 = nlev + 1;
        sz_e = (size_t)N * nlev;
        szA = aosoa_grpA_size(N, nlevp1);
        szC = aosoa_grpC_size(N);
        szD = aosoa_grpD_size(N_c, nlev);

        /* Pack grpA into AoSoA layout on heap */
        double *tmpA = new double[szA];
        memset(tmpA, 0, szA * sizeof(double));
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) {
                tmpA[IA_aosoa(A_VT,    je, jk, N)] = s._vt[I2(je, jk, N)];
                tmpA[IA_aosoa(A_VN_IE, je, jk, N)] = s._vn_ie[I2(je, jk, N)];
                tmpA[IA_aosoa(A_EKIN,  je, jk, N)] = s._ekin[I2(je, jk, N)];
            }
        for (int je = 0; je < N; je++)
            tmpA[IA_aosoa(A_VN_IE, je, nlev, N)] = s._vn_ie[I2(je, nlev, N)];

        /* Pack grpC */
        double *tmpC = new double[szC];
        memset(tmpC, 0, szC * sizeof(double));
        for (int je = 0; je < N; je++)
            tmpC[IC_aosoa(C_FE, je)] = s._fe[je];

        /* Pack grpD */
        double *tmpD = new double[szD];
        memset(tmpD, 0, szD * sizeof(double));
        for (int jk = 0; jk < nlev; jk++)
            for (int ic = 0; ic < N_c; ic++) {
                tmpD[ID_aosoa(D_EKINH, ic, jk, N_c)] = s._ekinh[I2(ic, jk, N_c)];
                tmpD[ID_aosoa(D_WCON,  ic, jk, N_c)] = s._wcon[I2(ic, jk, N_c)];
            }

        /* NUMA redistribute */
        grpA           = numa_redistribute_aosoa_grpA(tmpA, N, nlevp1);
        grpC           = numa_redistribute_aosoa_grpC(tmpC, N);
        grpD           = numa_redistribute_aosoa_grpD(tmpD, N_c, nlev);
        ddqz_z_full_e  = numa_redistribute_2d(s._ddqz, N, nlev);
        coeff_gradekin = numa_redistribute_1d(s._cg, N * 2);
        c_lin_e        = numa_redistribute_1d(s._cl, N * 2);
        zeta           = numa_redistribute_2d(s._zeta, N_v, nlev);
        cell_idx       = numa_redistribute_idx(s._cidx, N);
        vert_idx       = numa_redistribute_idx(s._vidx, N);
        out            = numa_alloc<double>(sz_e);
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < sz_e; i++) out[i] = 0.0;

        delete[] tmpA; delete[] tmpC; delete[] tmpD;
    }
    void numa_free_all() {
        numa_free(grpA, szA);
        numa_free(grpC, szC);
        numa_free(grpD, szD);
        numa_free(ddqz_z_full_e, sz_e);
        numa_free(coeff_gradekin, (size_t)N * 2);
        numa_free(c_lin_e, (size_t)N * 2);
        numa_free(zeta, (size_t)N_v * nlev);
        numa_free(cell_idx, (size_t)N * 2);
        numa_free(vert_idx, (size_t)N * 2);
        numa_free(out, sz_e);
    }
};

/* ================================================================ */
/*  CPU kernels                                                      */
/* ================================================================ */

static void kernel_soa(SoADataNuma &d) {
    const int N=d.N, N_c=d.N_c, N_v=d.N_v, nlev=d.nlev;
    double       *__restrict__ out             = d.out;
    const double *__restrict__ vt              = d.vt;
    const double *__restrict__ vn_ie           = d.vn_ie;
    const double *__restrict__ z_kin_hor_e     = d.z_kin_hor_e;
    const double *__restrict__ ddqz_z_full_e   = d.ddqz_z_full_e;
    const double *__restrict__ coeff_gradekin  = d.coeff_gradekin;
    const double *__restrict__ c_lin_e         = d.c_lin_e;
    const double *__restrict__ f_e             = d.f_e;
    const double *__restrict__ z_ekinh         = d.z_ekinh;
    const double *__restrict__ z_w_con_c_full  = d.z_w_con_c_full;
    const double *__restrict__ zeta            = d.zeta;
    const int    *__restrict__ cell_idx        = d.cell_idx;
    const int    *__restrict__ vert_idx        = d.vert_idx;
    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) { SOA_BODY() }
}

static void kernel_grp(GrpDataNuma &d) {
    const int N=d.N, N_c=d.N_c, N_v=d.N_v, nlev=d.nlev;
    double       *__restrict__ out             = d.out;
    const double *__restrict__ grpA            = d.grpA;
    const double *__restrict__ grpC            = d.grpC;
    const double *__restrict__ grpD            = d.grpD;
    const double *__restrict__ ddqz_z_full_e   = d.ddqz_z_full_e;
    const double *__restrict__ coeff_gradekin  = d.coeff_gradekin;
    const double *__restrict__ c_lin_e         = d.c_lin_e;
    const double *__restrict__ zeta            = d.zeta;
    const int    *__restrict__ cell_idx        = d.cell_idx;
    const int    *__restrict__ vert_idx        = d.vert_idx;
    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) { GRP_BODY() }
}

static void kernel_aosoa(AoSoADataNuma &d) {
    const int N=d.N, N_c=d.N_c, N_v=d.N_v, nlev=d.nlev;
    double       *__restrict__ out             = d.out;
    const double *__restrict__ grpA            = d.grpA;
    const double *__restrict__ grpC            = d.grpC;
    const double *__restrict__ grpD            = d.grpD;
    const double *__restrict__ ddqz_z_full_e   = d.ddqz_z_full_e;
    const double *__restrict__ coeff_gradekin  = d.coeff_gradekin;
    const double *__restrict__ c_lin_e         = d.c_lin_e;
    const double *__restrict__ zeta            = d.zeta;
    const int    *__restrict__ cell_idx        = d.cell_idx;
    const int    *__restrict__ vert_idx        = d.vert_idx;
    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) { AOSOA_BODY() }
}

/* Serial SoA reference (for verification) */
static void kernel_soa_ref(SoADataNuma &d, double *ref_out) {
    const int N=d.N, N_c=d.N_c, N_v=d.N_v, nlev=d.nlev;
    double       *__restrict__ out             = ref_out;
    const double *__restrict__ vt              = d.vt;
    const double *__restrict__ vn_ie           = d.vn_ie;
    const double *__restrict__ z_kin_hor_e     = d.z_kin_hor_e;
    const double *__restrict__ ddqz_z_full_e   = d.ddqz_z_full_e;
    const double *__restrict__ coeff_gradekin  = d.coeff_gradekin;
    const double *__restrict__ c_lin_e         = d.c_lin_e;
    const double *__restrict__ f_e             = d.f_e;
    const double *__restrict__ z_ekinh         = d.z_ekinh;
    const double *__restrict__ z_w_con_c_full  = d.z_w_con_c_full;
    const double *__restrict__ zeta            = d.zeta;
    const int    *__restrict__ cell_idx        = d.cell_idx;
    const int    *__restrict__ vert_idx        = d.vert_idx;
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) { SOA_BODY() }
}

/* ================================================================ */
/*  Timed runner for a single layout                                 */
/* ================================================================ */
typedef void (*kern_fn_t)(void*);

/* Wrappers */
static SoADataNuma   *_soa_ptr;
static GrpDataNuma   *_grp_ptr;
static AoSoADataNuma *_aosoa_ptr;
static void _ksoa(void*)   { kernel_soa(*_soa_ptr); }
static void _kgrp(void*)   { kernel_grp(*_grp_ptr); }
static void _kaosoa(void*) { kernel_aosoa(*_aosoa_ptr); }

static void run_layout(FILE *csv, const char *layout, const char *dist_name,
                       int nlev, int N,
                       kern_fn_t kern, void * /*data*/, double *out_ptr,
                       const double *ref, size_t sz_e) {
    /* Verify */
    kern(nullptr);
    int nf; double mr;
    verify(out_ptr, ref, sz_e, nf, mr);
    printf("  %-8s %-12s verify: %s  max_rel=%.2e\n",
           layout, dist_name, nf ? "FAIL" : "OK", mr);

    /* Warmup */
    for (int w = 0; w < WARMUP; w++) { flush(); kern(nullptr); }

    /* Timed runs */
    for (int r = 0; r < NRUNS; r++) {
        flush();
        auto t0 = std::chrono::high_resolution_clock::now();
        kern(nullptr);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        fprintf(csv, "%s,%d,%d,%s,%d,%.9f\n", layout, nlev, N, dist_name, r, ms);
    }
    fflush(csv);
}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main(int argc, char *argv[]) {
    /* Determine nlev sweep */
    std::vector<int> nlevs;
    if (argc >= 2) {
        nlevs.push_back(atoi(argv[1]));
    } else {
        for (int i = 0; i < DEFAULT_N_NLEVS; i++)
            nlevs.push_back(DEFAULT_NLEVS[i]);
    }

    const int N   = NPROMA;
    const int N_c = N, N_v = N;

    printf("ddt_vn_apc_pc CPU sweep benchmark\n");
    printf("  N=%d  threads=%d  NRUNS=%d  WARMUP=%d  VW=%d\n",
           N, omp_get_max_threads(), NRUNS, WARMUP, VW);
    printf("  Flush buffers: 2 × %.1f GB  (%d² Jacobi, %d sweeps)\n",
           (double)FN * FN * 8.0 / 1e9, FN, FS);
    printf("  nlev sweep:");
    for (int lv : nlevs) printf(" %d", lv);
    printf("\n");

    FILE *csv = fopen("ddt_vn_apc_pc_cpu_sweep.csv", "w");
    fprintf(csv, "layout,nlev,nproma,cell_dist,run_id,time_ms\n");

    std::mt19937 rng(42);

    flush_init();
    srand((unsigned)time(NULL));
    flush();
    printf("Ready\n\n");

    for (int nlev : nlevs) {
        size_t sz_e = (size_t)N * nlev;

        for (int di = 0; di < 4; di++) {
            CellDist dist = (CellDist)di;
            printf("=== nlev=%d  dist=%s ===\n", nlev, dist_names[di]);

            /* Generate connectivity */
            int *cell_log = new int[N * 2], *vert_log = new int[N * 2];
            gen_connectivity(cell_log, N, N_c, dist, rng);
            gen_connectivity(vert_log, N, N_v, dist, rng);

            /* Build SoA */
            SoADataNuma soa;
            soa.alloc(N, N_c, N_v, nlev);
            soa.fill();
            soa.set_connectivity(cell_log, vert_log);
            soa.numa_distribute();

            /* Build AoS grouped */
            GrpDataNuma grp;
            grp.from_soa_and_distribute(soa);

            /* Build AoSoA-16 */
            AoSoADataNuma aosoa;
            aosoa.from_soa_and_distribute(soa);

            /* Serial reference */
            double *ref = numa_alloc<double>(sz_e);
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < sz_e; i++) ref[i] = 0.0;
            kernel_soa_ref(soa, ref);

            /* Run SoA */
            _soa_ptr = &soa;
            run_layout(csv, "soa", dist_names[di], nlev, N,
                       _ksoa, nullptr, soa.out, ref, sz_e);

            /* Run AoS grouped */
            _grp_ptr = &grp;
            run_layout(csv, "aos", dist_names[di], nlev, N,
                       _kgrp, nullptr, grp.out, ref, sz_e);

            /* Run AoSoA-16 */
            _aosoa_ptr = &aosoa;
            run_layout(csv, "aosoa16", dist_names[di], nlev, N,
                       _kaosoa, nullptr, aosoa.out, ref, sz_e);

            /* Quick console summary (medians from the CSV runs already written) */
            {
                std::vector<double> ts(NRUNS), tg(NRUNS), ta(NRUNS);
                for (int r = 0; r < NRUNS; r++) {
                    flush();
                    auto t0 = std::chrono::high_resolution_clock::now();
                    kernel_soa(soa);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    ts[r] = std::chrono::duration<double,std::milli>(t1-t0).count();
                }
                for (int r = 0; r < NRUNS; r++) {
                    flush();
                    auto t0 = std::chrono::high_resolution_clock::now();
                    kernel_grp(grp);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    tg[r] = std::chrono::duration<double,std::milli>(t1-t0).count();
                }
                for (int r = 0; r < NRUNS; r++) {
                    flush();
                    auto t0 = std::chrono::high_resolution_clock::now();
                    kernel_aosoa(aosoa);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    ta[r] = std::chrono::duration<double,std::milli>(t1-t0).count();
                }
                std::sort(ts.begin(), ts.end());
                std::sort(tg.begin(), tg.end());
                std::sort(ta.begin(), ta.end());
                double ms_s = ts[NRUNS/2], ms_g = tg[NRUNS/2], ms_a = ta[NRUNS/2];
                printf("  ** MEDIAN  SoA=%.3f  AoS=%.3f (%.2fx)  AoSoA16=%.3f (%.2fx) ms **\n\n",
                       ms_s, ms_g, ms_s/ms_g, ms_a, ms_s/ms_a);
            }

            /* Cleanup */
            numa_free(ref, sz_e);
            soa.numa_free_all();
            soa.free_staging();
            grp.numa_free_all();
            aosoa.numa_free_all();
            delete[] cell_log;
            delete[] vert_log;
        }
    }

    flush_destroy();
    fclose(csv);
    printf("Written: ddt_vn_apc_pc_cpu_sweep.csv\n");
    printf("  Total rows: %zu nlevs × 4 dists × 3 layouts × %d reps = %zu\n",
           nlevs.size(), NRUNS, nlevs.size() * 4 * 3 * NRUNS);
    return 0;
}