/*
 * bench_ddt_vn_gpu_sweep.cpp -- ddt_vn_apc_pc GPU stencil benchmark
 *
 * Compares three layouts:
 *   SoA   — baseline, 15 distinct array streams
 *   AoS   — grouped (grpA/grpC/grpD), 10 streams
 *   AoSoA — grouped with vector-width tiling (V=32 or V=64), 10 streams
 *
 * Sweeps:
 *   - Threadblock sizes: {32,64,128,256} × {1,2,4,8}
 *   - Coarsening factors: TX × TY elements per thread
 *   - Connectivity distributions: uniform, normal σ=1, normal σ=2, sequential
 *   - AoSoA vector widths: 32, 64
 *
 * 100 repetitions per configuration; each individual timing written to CSV.
 * L2 cache flush between every timed kernel launch.
 *
 * Compile (HIP / MI300A):
 *   hipcc -O3 -std=c++17 bench_ddt_vn_gpu_sweep.cpp -o bench_ddt_vn_gpu_sweep
 *
 * Compile (CUDA / GH200):
 *   nvcc -O3 -std=c++17 bench_ddt_vn_gpu_sweep.cpp -o bench_ddt_vn_gpu_sweep
 *
 * Run:
 *   ./bench_ddt_vn_gpu_sweep [nlev]    (default nlev=90)
 */

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

/* ================================================================ */
/*  GPU backend abstraction                                          */
/* ================================================================ */
#ifdef __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>
#define GPU_CHECK(call) do {                                        \
    hipError_t e = (call);                                          \
    if (e != hipSuccess) {                                          \
        fprintf(stderr, "HIP error %s:%d: %s\n",                   \
                __FILE__, __LINE__, hipGetErrorString(e));           \
        exit(1);                                                    \
    }                                                               \
} while(0)
#define GPU_MALLOC       hipMalloc
#define GPU_FREE         hipFree
#define GPU_MEMCPY_H2D   hipMemcpyHostToDevice
#define GPU_MEMCPY_D2H   hipMemcpyDeviceToHost
#define GPU_MEMCPY       hipMemcpy
#define GPU_SYNC         hipDeviceSynchronize
#define GPU_EVENT_T      hipEvent_t
#define GPU_EVENT_CREATE  hipEventCreate
#define GPU_EVENT_RECORD  hipEventRecord
#define GPU_EVENT_SYNC    hipEventSynchronize
#define GPU_EVENT_ELAPSED hipEventElapsedTime
#define GPU_EVENT_DESTROY hipEventDestroy
#define GPU_MEMSET       hipMemset
#define HD __host__ __device__ __forceinline__
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#define GPU_CHECK(call) do {                                        \
    cudaError_t e = (call);                                         \
    if (e != cudaSuccess) {                                         \
        fprintf(stderr, "CUDA error %s:%d: %s\n",                  \
                __FILE__, __LINE__, cudaGetErrorString(e));          \
        exit(1);                                                    \
    }                                                               \
} while(0)
#define GPU_MALLOC       cudaMalloc
#define GPU_FREE         cudaFree
#define GPU_MEMCPY_H2D   cudaMemcpyHostToDevice
#define GPU_MEMCPY_D2H   cudaMemcpyDeviceToHost
#define GPU_MEMCPY       cudaMemcpy
#define GPU_SYNC         cudaDeviceSynchronize
#define GPU_EVENT_T      cudaEvent_t
#define GPU_EVENT_CREATE  cudaEventCreate
#define GPU_EVENT_RECORD  cudaEventRecord
#define GPU_EVENT_SYNC    cudaEventSynchronize
#define GPU_EVENT_ELAPSED cudaEventElapsedTime
#define GPU_EVENT_DESTROY cudaEventDestroy
#define GPU_MEMSET       cudaMemset
#define HD __host__ __device__ __forceinline__
#else
#error "This file requires HIP or CUDA.  Use bench_ddt_vn.cpp for CPU."
#endif

/* ================================================================ */
/*  Constants                                                        */
/* ================================================================ */
static constexpr int NPROMA = 81920;
static constexpr int NRUNS  = 100;
static constexpr int WARMUP = 10;

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
/*  Index helpers — SoA and AoS                                      */
/* ================================================================ */
HD int I2(int je, int jk, int N) { return je + jk * N; }
HD int IX(int je, int n, int N)  { return je + n * N; }
HD int IG(int je, int n, int N)  { return je + n * N; }
HD int IA(int f, int je, int jk, int N) { return f + NA * (je + N * jk); }
HD int IC_grp(int f, int je)             { return f + NC * je; }
HD int ID(int f, int cell, int jk, int Nc) { return f + ND * (cell + Nc * jk); }

/* ================================================================ */
/*  Index helpers — AoSoA                                            */
/* ================================================================ */
/*
 * AoSoA with vector width V:  within each tile of V consecutive
 * elements the F fields are interleaved in SoA order, giving
 * coalesced warp/wavefront access to each field.
 *
 * Layout for a 3D array [F, N, K]:
 *   index(f, je, jk) = jk * tiles_N * F * V  +  (je / V) * F * V  +  f * V  +  (je % V)
 *   where tiles_N = ceil(N / V)
 *
 * Layout for a 2D array [F, N]:
 *   index(f, je) = (je / V) * F * V  +  f * V  +  (je % V)
 */
template<int V, int F>
HD int IAoSoA_3D(int f, int je, int jk, int N) {
    int tiles = (N + V - 1) / V;
    return jk * tiles * F * V + (je / V) * F * V + f * V + (je % V);
}
template<int V, int F>
HD int IAoSoA_2D(int f, int je) {
    return (je / V) * F * V + f * V + (je % V);
}

/* Convenience wrappers */
template<int V> HD int IA_aosoa(int f, int je, int jk, int N)
    { return IAoSoA_3D<V, NA>(f, je, jk, N); }

template<int V> HD int IC_aosoa(int f, int je)
    { return IAoSoA_2D<V, NC>(f, je); }

template<int V> HD int ID_aosoa(int f, int cell, int jk, int Nc)
    { return IAoSoA_3D<V, ND>(f, cell, jk, Nc); }

/* Total allocation sizes for AoSoA */
template<int V>
size_t aosoa_grpA_size(int N, int nlevp1)
    { return (size_t)((N + V - 1) / V) * NA * V * nlevp1; }
template<int V>
size_t aosoa_grpC_size(int N)
    { return (size_t)((N + V - 1) / V) * NC * V; }
template<int V>
size_t aosoa_grpD_size(int Nc, int nlev)
    { return (size_t)((Nc + V - 1) / V) * ND * V * nlev; }

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

/* AoSoA grouped — 10 streams, tiled vector access */
#define AOSOA_BODY(V)                                                           \
    int ci0 = cell_idx[IX(je,0,N)], ci1 = cell_idx[IX(je,1,N)];                \
    int vi0 = vert_idx[IX(je,0,N)], vi1 = vert_idx[IX(je,1,N)];                \
    double cg0 = coeff_gradekin[IG(je,0,N)];                                   \
    double cg1 = coeff_gradekin[IG(je,1,N)];                                   \
    double cl0 = c_lin_e[IG(je,0,N)], cl1 = c_lin_e[IG(je,1,N)];              \
    double ekin_e = grpA[IA_aosoa<V>(A_EKIN, je, jk, N)];                      \
    double vt_e   = grpA[IA_aosoa<V>(A_VT,   je, jk, N)];                      \
    double vnie_k = grpA[IA_aosoa<V>(A_VN_IE,je, jk, N)];                      \
    double vnie_k1= grpA[IA_aosoa<V>(A_VN_IE,je, jk+1, N)];                    \
    double fe_e   = grpC[IC_aosoa<V>(C_FE, je)];                               \
    double ekinh0 = grpD[ID_aosoa<V>(D_EKINH, ci0, jk, N_c)];                  \
    double ekinh1 = grpD[ID_aosoa<V>(D_EKINH, ci1, jk, N_c)];                  \
    double wcon0  = grpD[ID_aosoa<V>(D_WCON,  ci0, jk, N_c)];                  \
    double wcon1  = grpD[ID_aosoa<V>(D_WCON,  ci1, jk, N_c)];                  \
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
/*  L2 cache flush                                                   */
/* ================================================================ */
/*
 * MI300A:  256 MB shared L2.   GH200:  ~100 MB L2.
 * We use 384 MB (48 M doubles) to be safe on both.
 */
static constexpr size_t FLUSH_ELEMS = 48ULL * 1024 * 1024;  /* 384 MB */
static double *d_flush_buf = nullptr;

__global__ void flush_kernel(double *buf, size_t n) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) buf[i] = buf[i] * 1.00001 + 1e-12;
}

static void flush_init() {
    GPU_CHECK(GPU_MALLOC(&d_flush_buf, FLUSH_ELEMS * sizeof(double)));
    GPU_CHECK(GPU_MEMSET(d_flush_buf, 0, FLUSH_ELEMS * sizeof(double)));
    /* Warm the buffer once */
    int blk = 256;
    int grd = (int)((FLUSH_ELEMS + blk - 1) / blk);
    flush_kernel<<<grd, blk>>>(d_flush_buf, FLUSH_ELEMS);
    GPU_CHECK(GPU_SYNC());
}
static void flush() {
    int blk = 256;
    int grd = (int)((FLUSH_ELEMS + blk - 1) / blk);
    flush_kernel<<<grd, blk>>>(d_flush_buf, FLUSH_ELEMS);
    GPU_CHECK(GPU_SYNC());
}
static void flush_destroy() {
    GPU_FREE(d_flush_buf);
}

/* ================================================================ */
/*  Host data structures                                             */
/* ================================================================ */

struct SoAData {
    int N, N_c, N_v, nlev;
    size_t sz_e, sz_c, sz_v;
    double *vt, *vn_ie, *z_kin_hor_e, *ddqz_z_full_e;
    double *coeff_gradekin, *c_lin_e, *f_e;
    double *z_ekinh, *z_w_con_c_full, *zeta;
    int    *cell_idx, *vert_idx;
    double *out;

    void alloc(int N_, int Nc, int Nv, int nlev_) {
        N=N_; N_c=Nc; N_v=Nv; nlev=nlev_;
        sz_e=(size_t)N*nlev; sz_c=(size_t)Nc*nlev; sz_v=(size_t)Nv*nlev;
        vt=new double[sz_e]; vn_ie=new double[(size_t)N*(nlev+1)];
        z_kin_hor_e=new double[sz_e]; ddqz_z_full_e=new double[sz_e];
        coeff_gradekin=new double[N*2]; c_lin_e=new double[N*2]; f_e=new double[N];
        z_ekinh=new double[sz_c]; z_w_con_c_full=new double[sz_c]; zeta=new double[sz_v];
        cell_idx=new int[N*2]; vert_idx=new int[N*2]; out=new double[sz_e];
    }
    void fill() {
        fill_rand(vt,sz_e,101); fill_rand(vn_ie,(size_t)N*(nlev+1),102);
        fill_rand(z_kin_hor_e,sz_e,103); fill_rand(ddqz_z_full_e,sz_e,104);
        fill_rand(coeff_gradekin,N*2,105); fill_rand(c_lin_e,N*2,106);
        fill_rand(f_e,N,107);
        fill_rand(z_ekinh,sz_c,108); fill_rand(z_w_con_c_full,sz_c,109);
        fill_rand(zeta,sz_v,110);
        memset(out,0,sz_e*sizeof(double));
        for (size_t i=0; i<sz_e; i++)
            if (std::abs(ddqz_z_full_e[i])<1e-10) ddqz_z_full_e[i]=1.0;
    }
    void free_all() {
        delete[] vt; delete[] vn_ie; delete[] z_kin_hor_e; delete[] ddqz_z_full_e;
        delete[] coeff_gradekin; delete[] c_lin_e; delete[] f_e;
        delete[] z_ekinh; delete[] z_w_con_c_full; delete[] zeta;
        delete[] cell_idx; delete[] vert_idx; delete[] out;
    }
};

struct GrpData {
    int N, N_c, N_v, nlev;
    size_t sz_e;
    double *grpA, *grpC, *grpD;
    double *ddqz_z_full_e, *coeff_gradekin, *c_lin_e, *zeta;
    int *cell_idx, *vert_idx;
    double *out;

    void from_soa(const SoAData &s) {
        N=s.N; N_c=s.N_c; N_v=s.N_v; nlev=s.nlev;
        sz_e=(size_t)N*nlev;
        int nlevp1=nlev+1;
        grpA=new double[(size_t)NA*N*nlevp1];
        grpC=new double[(size_t)NC*N];
        grpD=new double[(size_t)ND*N_c*nlev];
        memset(grpA,0,(size_t)NA*N*nlevp1*sizeof(double));
        for (int jk=0; jk<nlev; jk++)
            for (int je=0; je<N; je++) {
                grpA[IA(A_VT,je,jk,N)]    = s.vt[I2(je,jk,N)];
                grpA[IA(A_VN_IE,je,jk,N)] = s.vn_ie[I2(je,jk,N)];
                grpA[IA(A_EKIN,je,jk,N)]  = s.z_kin_hor_e[I2(je,jk,N)];
            }
        for (int je=0; je<N; je++)
            grpA[IA(A_VN_IE,je,nlev,N)] = s.vn_ie[I2(je,nlev,N)];
        memset(grpC,0,(size_t)NC*N*sizeof(double));
        for (int je=0; je<N; je++)
            grpC[IC_grp(C_FE,je)] = s.f_e[je];
        for (int jk=0; jk<nlev; jk++)
            for (int ic=0; ic<N_c; ic++) {
                grpD[ID(D_EKINH,ic,jk,N_c)] = s.z_ekinh[I2(ic,jk,N_c)];
                grpD[ID(D_WCON,ic,jk,N_c)]  = s.z_w_con_c_full[I2(ic,jk,N_c)];
            }
        ddqz_z_full_e=new double[sz_e];
        coeff_gradekin=new double[N*2]; c_lin_e=new double[N*2];
        zeta=new double[(size_t)N_v*nlev];
        cell_idx=new int[N*2]; vert_idx=new int[N*2]; out=new double[sz_e];
        memcpy(ddqz_z_full_e,s.ddqz_z_full_e,sz_e*sizeof(double));
        memcpy(coeff_gradekin,s.coeff_gradekin,N*2*sizeof(double));
        memcpy(c_lin_e,s.c_lin_e,N*2*sizeof(double));
        memcpy(zeta,s.zeta,(size_t)N_v*nlev*sizeof(double));
        memcpy(cell_idx,s.cell_idx,N*2*sizeof(int));
        memcpy(vert_idx,s.vert_idx,N*2*sizeof(int));
        memset(out,0,sz_e*sizeof(double));
    }
    void free_all() {
        delete[] grpA; delete[] grpC; delete[] grpD;
        delete[] ddqz_z_full_e; delete[] coeff_gradekin; delete[] c_lin_e;
        delete[] zeta; delete[] cell_idx; delete[] vert_idx; delete[] out;
    }
};

template<int V>
struct AoSoAData {
    int N, N_c, N_v, nlev;
    size_t sz_e;
    double *grpA, *grpC, *grpD;
    double *ddqz_z_full_e, *coeff_gradekin, *c_lin_e, *zeta;
    int *cell_idx, *vert_idx;
    double *out;
    size_t szA, szC, szD;

    void from_soa(const SoAData &s) {
        N=s.N; N_c=s.N_c; N_v=s.N_v; nlev=s.nlev;
        sz_e=(size_t)N*nlev;
        int nlevp1=nlev+1;
        szA = aosoa_grpA_size<V>(N, nlevp1);
        szC = aosoa_grpC_size<V>(N);
        szD = aosoa_grpD_size<V>(N_c, nlev);
        grpA=new double[szA]; memset(grpA,0,szA*sizeof(double));
        grpC=new double[szC]; memset(grpC,0,szC*sizeof(double));
        grpD=new double[szD]; memset(grpD,0,szD*sizeof(double));
        for (int jk=0; jk<nlev; jk++)
            for (int je=0; je<N; je++) {
                grpA[IA_aosoa<V>(A_VT,je,jk,N)]    = s.vt[I2(je,jk,N)];
                grpA[IA_aosoa<V>(A_VN_IE,je,jk,N)] = s.vn_ie[I2(je,jk,N)];
                grpA[IA_aosoa<V>(A_EKIN,je,jk,N)]  = s.z_kin_hor_e[I2(je,jk,N)];
            }
        for (int je=0; je<N; je++)
            grpA[IA_aosoa<V>(A_VN_IE,je,nlev,N)] = s.vn_ie[I2(je,nlev,N)];
        for (int je=0; je<N; je++)
            grpC[IC_aosoa<V>(C_FE,je)] = s.f_e[je];
        for (int jk=0; jk<nlev; jk++)
            for (int ic=0; ic<N_c; ic++) {
                grpD[ID_aosoa<V>(D_EKINH,ic,jk,N_c)] = s.z_ekinh[I2(ic,jk,N_c)];
                grpD[ID_aosoa<V>(D_WCON,ic,jk,N_c)]   = s.z_w_con_c_full[I2(ic,jk,N_c)];
            }
        ddqz_z_full_e=new double[sz_e];
        coeff_gradekin=new double[N*2]; c_lin_e=new double[N*2];
        zeta=new double[(size_t)N_v*nlev];
        cell_idx=new int[N*2]; vert_idx=new int[N*2]; out=new double[sz_e];
        memcpy(ddqz_z_full_e,s.ddqz_z_full_e,sz_e*sizeof(double));
        memcpy(coeff_gradekin,s.coeff_gradekin,N*2*sizeof(double));
        memcpy(c_lin_e,s.c_lin_e,N*2*sizeof(double));
        memcpy(zeta,s.zeta,(size_t)N_v*nlev*sizeof(double));
        memcpy(cell_idx,s.cell_idx,N*2*sizeof(int));
        memcpy(vert_idx,s.vert_idx,N*2*sizeof(int));
        memset(out,0,sz_e*sizeof(double));
    }
    void free_all() {
        delete[] grpA; delete[] grpC; delete[] grpD;
        delete[] ddqz_z_full_e; delete[] coeff_gradekin; delete[] c_lin_e;
        delete[] zeta; delete[] cell_idx; delete[] vert_idx; delete[] out;
    }
};

/* ================================================================ */
/*  GPU device buffers (all layouts share standalone arrays)          */
/* ================================================================ */
struct DevSoA {
    double *vt, *vn_ie, *ekin, *ddqz, *cg, *cl, *fe;
    double *ekinh, *wcon, *zeta, *out;
    int *cidx, *vidx;
};
struct DevGrp {
    double *grpA, *grpC, *grpD;
    double *ddqz, *cg, *cl, *zeta, *out;
    int *cidx, *vidx;
};
struct DevAoSoA {
    double *grpA, *grpC, *grpD;
    double *ddqz, *cg, *cl, *zeta, *out;
    int *cidx, *vidx;
    size_t szA, szC, szD;
};

/* ================================================================ */
/*  GPU kernels — SoA, templated on coarsening TX × TY               */
/* ================================================================ */
template<int TX, int TY>
__global__ void kern_soa(
    double *__restrict__ out,
    const double *__restrict__ vt, const double *__restrict__ vn_ie,
    const double *__restrict__ z_kin_hor_e, const double *__restrict__ ddqz_z_full_e,
    const double *__restrict__ coeff_gradekin, const double *__restrict__ c_lin_e,
    const double *__restrict__ f_e,
    const double *__restrict__ z_ekinh, const double *__restrict__ z_w_con_c_full,
    const double *__restrict__ zeta,
    const int *__restrict__ cell_idx, const int *__restrict__ vert_idx,
    int N, int N_c, int N_v, int nlev)
{
    int je_base = ((int)blockIdx.x * (int)blockDim.x + (int)threadIdx.x) * TX;
    int jk_base = ((int)blockIdx.y * (int)blockDim.y + (int)threadIdx.y) * TY;
    #pragma unroll
    for (int dy = 0; dy < TY; dy++) {
        int jk = jk_base + dy;
        if (jk >= nlev) break;
        #pragma unroll
        for (int dx = 0; dx < TX; dx++) {
            int je = je_base + dx;
            if (je >= N) break;
            SOA_BODY()
        }
    }
}

/* ================================================================ */
/*  GPU kernels — AoS grouped, templated on TX × TY                  */
/* ================================================================ */
template<int TX, int TY>
__global__ void kern_grp(
    double *__restrict__ out,
    const double *__restrict__ grpA, const double *__restrict__ grpC,
    const double *__restrict__ grpD,
    const double *__restrict__ ddqz_z_full_e,
    const double *__restrict__ coeff_gradekin, const double *__restrict__ c_lin_e,
    const double *__restrict__ zeta,
    const int *__restrict__ cell_idx, const int *__restrict__ vert_idx,
    int N, int N_c, int N_v, int nlev)
{
    int je_base = ((int)blockIdx.x * (int)blockDim.x + (int)threadIdx.x) * TX;
    int jk_base = ((int)blockIdx.y * (int)blockDim.y + (int)threadIdx.y) * TY;
    #pragma unroll
    for (int dy = 0; dy < TY; dy++) {
        int jk = jk_base + dy;
        if (jk >= nlev) break;
        #pragma unroll
        for (int dx = 0; dx < TX; dx++) {
            int je = je_base + dx;
            if (je >= N) break;
            GRP_BODY()
        }
    }
}

/* ================================================================ */
/*  GPU kernels — AoSoA, templated on V, TX, TY                     */
/* ================================================================ */
template<int V, int TX, int TY>
__global__ void kern_aosoa(
    double *__restrict__ out,
    const double *__restrict__ grpA, const double *__restrict__ grpC,
    const double *__restrict__ grpD,
    const double *__restrict__ ddqz_z_full_e,
    const double *__restrict__ coeff_gradekin, const double *__restrict__ c_lin_e,
    const double *__restrict__ zeta,
    const int *__restrict__ cell_idx, const int *__restrict__ vert_idx,
    int N, int N_c, int N_v, int nlev)
{
    int je_base = ((int)blockIdx.x * (int)blockDim.x + (int)threadIdx.x) * TX;
    int jk_base = ((int)blockIdx.y * (int)blockDim.y + (int)threadIdx.y) * TY;
    #pragma unroll
    for (int dy = 0; dy < TY; dy++) {
        int jk = jk_base + dy;
        if (jk >= nlev) break;
        #pragma unroll
        for (int dx = 0; dx < TX; dx++) {
            int je = je_base + dx;
            if (je >= N) break;
            AOSOA_BODY(V)
        }
    }
}

/* ================================================================ */
/*  Verification                                                     */
/* ================================================================ */
static bool verify(const double *got, const double *ref, size_t n,
                   int &nfail, double &maxrel) {
    nfail=0; maxrel=0;
    for (size_t i=0; i<n; i++) {
        double d = std::abs(got[i]-ref[i]);
        double dn = std::max(std::abs(ref[i]),1e-300);
        double r = d/dn;
        if (r>maxrel) maxrel=r;
        if (d>1e-12+1e-8*std::abs(ref[i])) nfail++;
    }
    return nfail==0;
}

/* ================================================================ */
/*  Kernel launch helpers                                            */
/* ================================================================ */
struct KernelConfig {
    int bx, by;   /* block size */
    int tx, ty;   /* coarsening: elements per thread */
};

/* Configuration arrays */
static const KernelConfig configs[] = {
    /* bx   by  tx  ty */
    {  32,  1,  1,  1 },
    {  64,  1,  1,  1 },
    { 128,  1,  1,  1 },
    { 256,  1,  1,  1 },
    {  64,  2,  1,  1 },
    { 128,  2,  1,  1 },
    {  64,  4,  1,  1 },
    {  32,  4,  1,  1 },
    {  32,  8,  1,  1 },
    /* coarsening in je */
    {  64,  1,  2,  1 },
    { 128,  1,  2,  1 },
    {  64,  2,  2,  1 },
    {  64,  1,  4,  1 },
    { 128,  1,  4,  1 },
    /* coarsening in jk */
    {  64,  1,  1,  2 },
    { 128,  1,  1,  2 },
    {  64,  2,  1,  2 },
    {  64,  1,  1,  4 },
    { 128,  1,  1,  4 },
    /* coarsening in both */
    {  64,  1,  2,  2 },
    { 128,  1,  2,  2 },
    {  64,  2,  2,  2 },
    {  64,  1,  4,  2 },
    {  64,  1,  2,  4 },
};
static constexpr int NCONFIGS = sizeof(configs)/sizeof(configs[0]);

/* Macro to dispatch (TX,TY) at compile time — instantiates the needed combos */
#define DISPATCH_TXTY(KERN, tx, ty, grid, block, ...) do {                     \
    if      ((tx)==1 && (ty)==1) KERN<1,1><<<grid,block>>>(__VA_ARGS__);        \
    else if ((tx)==2 && (ty)==1) KERN<2,1><<<grid,block>>>(__VA_ARGS__);        \
    else if ((tx)==4 && (ty)==1) KERN<4,1><<<grid,block>>>(__VA_ARGS__);        \
    else if ((tx)==1 && (ty)==2) KERN<1,2><<<grid,block>>>(__VA_ARGS__);        \
    else if ((tx)==1 && (ty)==4) KERN<1,4><<<grid,block>>>(__VA_ARGS__);        \
    else if ((tx)==2 && (ty)==2) KERN<2,2><<<grid,block>>>(__VA_ARGS__);        \
    else if ((tx)==4 && (ty)==2) KERN<4,2><<<grid,block>>>(__VA_ARGS__);        \
    else if ((tx)==2 && (ty)==4) KERN<2,4><<<grid,block>>>(__VA_ARGS__);        \
    else { fprintf(stderr, "Unsupported TX=%d TY=%d\n", tx, ty); exit(1); }    \
} while(0)

/* Same but with an extra template parameter V for AoSoA */
#define DISPATCH_V_TXTY(KERN, V, tx, ty, grid, block, ...) do {                \
    if      ((tx)==1 && (ty)==1) KERN<V,1,1><<<grid,block>>>(__VA_ARGS__);      \
    else if ((tx)==2 && (ty)==1) KERN<V,2,1><<<grid,block>>>(__VA_ARGS__);      \
    else if ((tx)==4 && (ty)==1) KERN<V,4,1><<<grid,block>>>(__VA_ARGS__);      \
    else if ((tx)==1 && (ty)==2) KERN<V,1,2><<<grid,block>>>(__VA_ARGS__);      \
    else if ((tx)==1 && (ty)==4) KERN<V,1,4><<<grid,block>>>(__VA_ARGS__);      \
    else if ((tx)==2 && (ty)==2) KERN<V,2,2><<<grid,block>>>(__VA_ARGS__);      \
    else if ((tx)==4 && (ty)==2) KERN<V,4,2><<<grid,block>>>(__VA_ARGS__);      \
    else if ((tx)==2 && (ty)==4) KERN<V,2,4><<<grid,block>>>(__VA_ARGS__);      \
    else { fprintf(stderr, "Unsupported TX=%d TY=%d\n", tx, ty); exit(1); }    \
} while(0)

/* ================================================================ */
/*  Device allocation helpers                                        */
/* ================================================================ */
static void alloc_dev_soa(DevSoA &d, const SoAData &h) {
    int N=h.N, nlev=h.nlev;
    GPU_CHECK(GPU_MALLOC(&d.vt,    h.sz_e*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.vn_ie, (size_t)N*(nlev+1)*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.ekin,  h.sz_e*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.ddqz,  h.sz_e*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cg,    N*2*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cl,    N*2*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.fe,    N*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.ekinh, h.sz_c*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.wcon,  h.sz_c*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.zeta,  h.sz_v*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cidx,  N*2*sizeof(int)));
    GPU_CHECK(GPU_MALLOC(&d.vidx,  N*2*sizeof(int)));
    GPU_CHECK(GPU_MALLOC(&d.out,   h.sz_e*sizeof(double)));

    GPU_CHECK(GPU_MEMCPY(d.vt,    h.vt,             h.sz_e*sizeof(double),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.vn_ie, h.vn_ie,          (size_t)N*(nlev+1)*sizeof(double), GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.ekin,  h.z_kin_hor_e,    h.sz_e*sizeof(double),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.ddqz,  h.ddqz_z_full_e,  h.sz_e*sizeof(double),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cg,    h.coeff_gradekin, N*2*sizeof(double),                GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cl,    h.c_lin_e,        N*2*sizeof(double),                GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.fe,    h.f_e,            N*sizeof(double),                  GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.ekinh, h.z_ekinh,        h.sz_c*sizeof(double),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.wcon,  h.z_w_con_c_full, h.sz_c*sizeof(double),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.zeta,  h.zeta,           h.sz_v*sizeof(double),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cidx,  h.cell_idx,       N*2*sizeof(int),                   GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.vidx,  h.vert_idx,       N*2*sizeof(int),                   GPU_MEMCPY_H2D));
}
static void free_dev_soa(DevSoA &d) {
    GPU_FREE(d.vt); GPU_FREE(d.vn_ie); GPU_FREE(d.ekin); GPU_FREE(d.ddqz);
    GPU_FREE(d.cg); GPU_FREE(d.cl); GPU_FREE(d.fe);
    GPU_FREE(d.ekinh); GPU_FREE(d.wcon); GPU_FREE(d.zeta);
    GPU_FREE(d.cidx); GPU_FREE(d.vidx); GPU_FREE(d.out);
}

static void alloc_dev_grp(DevGrp &d, const GrpData &h) {
    int N=h.N, Nc=h.N_c, Nv=h.N_v, nlev=h.nlev;
    GPU_CHECK(GPU_MALLOC(&d.grpA, (size_t)NA*N*(nlev+1)*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.grpC, (size_t)NC*N*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.grpD, (size_t)ND*Nc*nlev*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.ddqz, h.sz_e*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cg,   N*2*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cl,   N*2*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.zeta, (size_t)Nv*nlev*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cidx, N*2*sizeof(int)));
    GPU_CHECK(GPU_MALLOC(&d.vidx, N*2*sizeof(int)));
    GPU_CHECK(GPU_MALLOC(&d.out,  h.sz_e*sizeof(double)));

    GPU_CHECK(GPU_MEMCPY(d.grpA, h.grpA, (size_t)NA*N*(nlev+1)*sizeof(double), GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.grpC, h.grpC, (size_t)NC*N*sizeof(double),          GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.grpD, h.grpD, (size_t)ND*Nc*nlev*sizeof(double),    GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.ddqz, h.ddqz_z_full_e, h.sz_e*sizeof(double),       GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cg,   h.coeff_gradekin, N*2*sizeof(double),          GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cl,   h.c_lin_e,        N*2*sizeof(double),          GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.zeta, h.zeta,           (size_t)Nv*nlev*sizeof(double), GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cidx, h.cell_idx,       N*2*sizeof(int),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.vidx, h.vert_idx,       N*2*sizeof(int),             GPU_MEMCPY_H2D));
}
static void free_dev_grp(DevGrp &d) {
    GPU_FREE(d.grpA); GPU_FREE(d.grpC); GPU_FREE(d.grpD);
    GPU_FREE(d.ddqz); GPU_FREE(d.cg); GPU_FREE(d.cl);
    GPU_FREE(d.zeta); GPU_FREE(d.cidx); GPU_FREE(d.vidx); GPU_FREE(d.out);
}

template<int V>
static void alloc_dev_aosoa(DevAoSoA &d, const AoSoAData<V> &h) {
    int N=h.N, Nc=h.N_c, Nv=h.N_v, nlev=h.nlev;
    d.szA=h.szA; d.szC=h.szC; d.szD=h.szD;
    GPU_CHECK(GPU_MALLOC(&d.grpA, h.szA*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.grpC, h.szC*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.grpD, h.szD*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.ddqz, h.sz_e*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cg,   N*2*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cl,   N*2*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.zeta, (size_t)Nv*nlev*sizeof(double)));
    GPU_CHECK(GPU_MALLOC(&d.cidx, N*2*sizeof(int)));
    GPU_CHECK(GPU_MALLOC(&d.vidx, N*2*sizeof(int)));
    GPU_CHECK(GPU_MALLOC(&d.out,  h.sz_e*sizeof(double)));

    GPU_CHECK(GPU_MEMCPY(d.grpA, h.grpA, h.szA*sizeof(double),                    GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.grpC, h.grpC, h.szC*sizeof(double),                    GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.grpD, h.grpD, h.szD*sizeof(double),                    GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.ddqz, h.ddqz_z_full_e, h.sz_e*sizeof(double),          GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cg,   h.coeff_gradekin, N*2*sizeof(double),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cl,   h.c_lin_e,        N*2*sizeof(double),             GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.zeta, h.zeta,           (size_t)Nv*nlev*sizeof(double), GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.cidx, h.cell_idx,       N*2*sizeof(int),                GPU_MEMCPY_H2D));
    GPU_CHECK(GPU_MEMCPY(d.vidx, h.vert_idx,       N*2*sizeof(int),                GPU_MEMCPY_H2D));
}
static void free_dev_aosoa(DevAoSoA &d) {
    GPU_FREE(d.grpA); GPU_FREE(d.grpC); GPU_FREE(d.grpD);
    GPU_FREE(d.ddqz); GPU_FREE(d.cg); GPU_FREE(d.cl);
    GPU_FREE(d.zeta); GPU_FREE(d.cidx); GPU_FREE(d.vidx); GPU_FREE(d.out);
}

/* ================================================================ */
/*  Benchmark runner                                                 */
/* ================================================================ */
static void bench_soa(FILE *csv, const char *dist_name,
                      const KernelConfig &cfg,
                      DevSoA &d, int N, int Nc, int Nv, int nlev,
                      GPU_EVENT_T ev0, GPU_EVENT_T ev1) {
    dim3 block(cfg.bx, cfg.by);
    dim3 grid((N  + block.x*cfg.tx - 1) / (block.x*cfg.tx),
              (nlev + block.y*cfg.ty - 1) / (block.y*cfg.ty));

    /* Warmup */
    for (int w = 0; w < WARMUP; w++) {
        DISPATCH_TXTY(kern_soa, cfg.tx, cfg.ty, grid, block,
            d.out, d.vt, d.vn_ie, d.ekin, d.ddqz,
            d.cg, d.cl, d.fe, d.ekinh, d.wcon, d.zeta,
            d.cidx, d.vidx, N, Nc, Nv, nlev);
    }
    GPU_CHECK(GPU_SYNC());

    /* Timed */
    for (int r = 0; r < NRUNS; r++) {
        flush();
        GPU_CHECK(GPU_EVENT_RECORD(ev0));
        DISPATCH_TXTY(kern_soa, cfg.tx, cfg.ty, grid, block,
            d.out, d.vt, d.vn_ie, d.ekin, d.ddqz,
            d.cg, d.cl, d.fe, d.ekinh, d.wcon, d.zeta,
            d.cidx, d.vidx, N, Nc, Nv, nlev);
        GPU_CHECK(GPU_EVENT_RECORD(ev1));
        GPU_CHECK(GPU_EVENT_SYNC(ev1));
        float ms; GPU_CHECK(GPU_EVENT_ELAPSED(&ms, ev0, ev1));
        fprintf(csv, "soa,%s,%d,%d,%d,%d,0,%d,%.6f\n",
                dist_name, cfg.bx, cfg.by, cfg.tx, cfg.ty, r, (double)ms);
    }
}

static void bench_grp(FILE *csv, const char *dist_name,
                      const KernelConfig &cfg,
                      DevGrp &d, int N, int Nc, int Nv, int nlev,
                      GPU_EVENT_T ev0, GPU_EVENT_T ev1) {
    dim3 block(cfg.bx, cfg.by);
    dim3 grid((N  + block.x*cfg.tx - 1) / (block.x*cfg.tx),
              (nlev + block.y*cfg.ty - 1) / (block.y*cfg.ty));

    for (int w = 0; w < WARMUP; w++) {
        DISPATCH_TXTY(kern_grp, cfg.tx, cfg.ty, grid, block,
            d.out, d.grpA, d.grpC, d.grpD,
            d.ddqz, d.cg, d.cl, d.zeta,
            d.cidx, d.vidx, N, Nc, Nv, nlev);
    }
    GPU_CHECK(GPU_SYNC());

    for (int r = 0; r < NRUNS; r++) {
        flush();
        GPU_CHECK(GPU_EVENT_RECORD(ev0));
        DISPATCH_TXTY(kern_grp, cfg.tx, cfg.ty, grid, block,
            d.out, d.grpA, d.grpC, d.grpD,
            d.ddqz, d.cg, d.cl, d.zeta,
            d.cidx, d.vidx, N, Nc, Nv, nlev);
        GPU_CHECK(GPU_EVENT_RECORD(ev1));
        GPU_CHECK(GPU_EVENT_SYNC(ev1));
        float ms; GPU_CHECK(GPU_EVENT_ELAPSED(&ms, ev0, ev1));
        fprintf(csv, "aos,%s,%d,%d,%d,%d,0,%d,%.6f\n",
                dist_name, cfg.bx, cfg.by, cfg.tx, cfg.ty, r, (double)ms);
    }
}

template<int V>
static void bench_aosoa(FILE *csv, const char *dist_name,
                        const KernelConfig &cfg,
                        DevAoSoA &d, int N, int Nc, int Nv, int nlev,
                        GPU_EVENT_T ev0, GPU_EVENT_T ev1) {
    dim3 block(cfg.bx, cfg.by);
    dim3 grid((N  + block.x*cfg.tx - 1) / (block.x*cfg.tx),
              (nlev + block.y*cfg.ty - 1) / (block.y*cfg.ty));

    for (int w = 0; w < WARMUP; w++) {
        DISPATCH_V_TXTY(kern_aosoa, V, cfg.tx, cfg.ty, grid, block,
            d.out, d.grpA, d.grpC, d.grpD,
            d.ddqz, d.cg, d.cl, d.zeta,
            d.cidx, d.vidx, N, Nc, Nv, nlev);
    }
    GPU_CHECK(GPU_SYNC());

    for (int r = 0; r < NRUNS; r++) {
        flush();
        GPU_CHECK(GPU_EVENT_RECORD(ev0));
        DISPATCH_V_TXTY(kern_aosoa, V, cfg.tx, cfg.ty, grid, block,
            d.out, d.grpA, d.grpC, d.grpD,
            d.ddqz, d.cg, d.cl, d.zeta,
            d.cidx, d.vidx, N, Nc, Nv, nlev);
        GPU_CHECK(GPU_EVENT_RECORD(ev1));
        GPU_CHECK(GPU_EVENT_SYNC(ev1));
        float ms; GPU_CHECK(GPU_EVENT_ELAPSED(&ms, ev0, ev1));
        fprintf(csv, "aosoa,%s,%d,%d,%d,%d,%d,%d,%.6f\n",
                dist_name, cfg.bx, cfg.by, cfg.tx, cfg.ty, V, r, (double)ms);
    }
}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main(int argc, char *argv[]) {
    int nlev = (argc >= 2) ? atoi(argv[1]) : 90;
    int N    = NPROMA;
    int N_c  = N, N_v = N;

    printf("ddt_vn_apc_pc GPU sweep benchmark\n");
    printf("  N=%d  nlev=%d  NRUNS=%d  WARMUP=%d  configs=%d\n",
           N, nlev, NRUNS, WARMUP, NCONFIGS);
    printf("  Flush buffer: %.0f MB\n", FLUSH_ELEMS*8.0/1e6);

    FILE *csv = fopen("ddt_vn_apc_pc_gpu_sweep.csv", "w");
    fprintf(csv, "layout,cell_dist,block_x,block_y,coarsen_x,coarsen_y,vec_width,run_id,time_ms\n");

    std::mt19937 rng(42);
    flush_init();

    GPU_EVENT_T ev0, ev1;
    GPU_CHECK(GPU_EVENT_CREATE(&ev0));
    GPU_CHECK(GPU_EVENT_CREATE(&ev1));

    for (int di = 0; di < 4; di++) {
        CellDist dist = (CellDist)di;
        printf("\n========== dist=%s ==========\n", dist_names[di]);

        /* ---- Build host data ---- */
        SoAData soa;
        soa.alloc(N, N_c, N_v, nlev);
        soa.fill();
        gen_connectivity(soa.cell_idx, N, N_c, dist, rng);
        gen_connectivity(soa.vert_idx, N, N_v, dist, rng);

        GrpData grp;
        grp.from_soa(soa);

        AoSoAData<32> aosoa32;
        aosoa32.from_soa(soa);

        AoSoAData<64> aosoa64;
        aosoa64.from_soa(soa);

        /* ---- CPU reference ---- */
        double *h_ref = new double[soa.sz_e];
        memset(soa.out, 0, soa.sz_e * sizeof(double));
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) {
                const double *vt = soa.vt, *vn_ie = soa.vn_ie;
                const double *z_kin_hor_e = soa.z_kin_hor_e;
                const double *ddqz_z_full_e = soa.ddqz_z_full_e;
                const double *coeff_gradekin = soa.coeff_gradekin;
                const double *c_lin_e = soa.c_lin_e;
                const double *f_e = soa.f_e;
                const double *z_ekinh = soa.z_ekinh;
                const double *z_w_con_c_full = soa.z_w_con_c_full;
                const double *zeta = soa.zeta;
                const int *cell_idx = soa.cell_idx, *vert_idx = soa.vert_idx;
                double *out = soa.out;
                SOA_BODY()
            }
        memcpy(h_ref, soa.out, soa.sz_e * sizeof(double));

        /* ---- Allocate device ---- */
        DevSoA   dev_soa;   alloc_dev_soa(dev_soa, soa);
        DevGrp   dev_grp;   alloc_dev_grp(dev_grp, grp);
        DevAoSoA dev_a32;   alloc_dev_aosoa<32>(dev_a32, aosoa32);
        DevAoSoA dev_a64;   alloc_dev_aosoa<64>(dev_a64, aosoa64);

        /* ---- Verify all layouts (baseline config) ---- */
        {
            dim3 block(64, 1);
            dim3 grid((N+63)/64, nlev);
            double *h_out = new double[soa.sz_e];
            int nf; double mr;

            kern_soa<1,1><<<grid,block>>>(
                dev_soa.out, dev_soa.vt, dev_soa.vn_ie, dev_soa.ekin, dev_soa.ddqz,
                dev_soa.cg, dev_soa.cl, dev_soa.fe, dev_soa.ekinh, dev_soa.wcon,
                dev_soa.zeta, dev_soa.cidx, dev_soa.vidx, N, N_c, N_v, nlev);
            GPU_CHECK(GPU_SYNC());
            GPU_CHECK(GPU_MEMCPY(h_out, dev_soa.out, soa.sz_e*sizeof(double), GPU_MEMCPY_D2H));
            verify(h_out, h_ref, soa.sz_e, nf, mr);
            printf("  SoA    verify: %s  max_rel=%.2e\n", nf?"FAIL":"OK", mr);

            kern_grp<1,1><<<grid,block>>>(
                dev_grp.out, dev_grp.grpA, dev_grp.grpC, dev_grp.grpD,
                dev_grp.ddqz, dev_grp.cg, dev_grp.cl, dev_grp.zeta,
                dev_grp.cidx, dev_grp.vidx, N, N_c, N_v, nlev);
            GPU_CHECK(GPU_SYNC());
            GPU_CHECK(GPU_MEMCPY(h_out, dev_grp.out, soa.sz_e*sizeof(double), GPU_MEMCPY_D2H));
            verify(h_out, h_ref, soa.sz_e, nf, mr);
            printf("  AoS    verify: %s  max_rel=%.2e\n", nf?"FAIL":"OK", mr);

            kern_aosoa<32,1,1><<<grid,block>>>(
                dev_a32.out, dev_a32.grpA, dev_a32.grpC, dev_a32.grpD,
                dev_a32.ddqz, dev_a32.cg, dev_a32.cl, dev_a32.zeta,
                dev_a32.cidx, dev_a32.vidx, N, N_c, N_v, nlev);
            GPU_CHECK(GPU_SYNC());
            GPU_CHECK(GPU_MEMCPY(h_out, dev_a32.out, soa.sz_e*sizeof(double), GPU_MEMCPY_D2H));
            verify(h_out, h_ref, soa.sz_e, nf, mr);
            printf("  AoSoA32 verify: %s  max_rel=%.2e\n", nf?"FAIL":"OK", mr);

            kern_aosoa<64,1,1><<<grid,block>>>(
                dev_a64.out, dev_a64.grpA, dev_a64.grpC, dev_a64.grpD,
                dev_a64.ddqz, dev_a64.cg, dev_a64.cl, dev_a64.zeta,
                dev_a64.cidx, dev_a64.vidx, N, N_c, N_v, nlev);
            GPU_CHECK(GPU_SYNC());
            GPU_CHECK(GPU_MEMCPY(h_out, dev_a64.out, soa.sz_e*sizeof(double), GPU_MEMCPY_D2H));
            verify(h_out, h_ref, soa.sz_e, nf, mr);
            printf("  AoSoA64 verify: %s  max_rel=%.2e\n", nf?"FAIL":"OK", mr);

            delete[] h_out;
        }

        /* ---- Sweep all configurations ---- */
        for (int ci = 0; ci < NCONFIGS; ci++) {
            const KernelConfig &cfg = configs[ci];
            printf("  config %2d/%d: block=(%3d,%d) coarsen=(%d,%d) ...",
                   ci+1, NCONFIGS, cfg.bx, cfg.by, cfg.tx, cfg.ty);
            fflush(stdout);

            bench_soa(csv, dist_names[di], cfg, dev_soa, N, N_c, N_v, nlev, ev0, ev1);
            bench_grp(csv, dist_names[di], cfg, dev_grp, N, N_c, N_v, nlev, ev0, ev1);
            bench_aosoa<32>(csv, dist_names[di], cfg, dev_a32, N, N_c, N_v, nlev, ev0, ev1);
            bench_aosoa<64>(csv, dist_names[di], cfg, dev_a64, N, N_c, N_v, nlev, ev0, ev1);

            fflush(csv);
            printf(" done\n");
        }

        /* ---- Cleanup ---- */
        free_dev_soa(dev_soa);
        free_dev_grp(dev_grp);
        free_dev_aosoa(dev_a32);
        free_dev_aosoa(dev_a64);
        soa.free_all();
        grp.free_all();
        aosoa32.free_all();
        aosoa64.free_all();
        delete[] h_ref;
    }

    GPU_CHECK(GPU_EVENT_DESTROY(ev0));
    GPU_CHECK(GPU_EVENT_DESTROY(ev1));
    flush_destroy();
    fclose(csv);
    printf("\nResults written to: ddt_vn_apc_pc_gpu_sweep.csv\n");
    printf("  Total rows: %d dist × %d configs × 4 layouts × %d reps = %d\n",
           4, NCONFIGS, NRUNS, 4 * NCONFIGS * 4 * NRUNS);
    return 0;
}