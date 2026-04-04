/*
 * bench_ddt_vn.cpp -- ddt_vn_apc_pc stencil benchmark
 *
 * Compares SoA baseline (15 array streams) vs AoS-grouped layout (10 streams).
 *
 * Groups applied:
 *   grpA[5, N, nlevp1]  = {vn, vt, vn_ie, z_vt_ie, z_kin_hor_e}  (edge wind state)
 *   grpC[4, N]           = {inv_dual, inv_primal, tangent, f_e}     (edge geometry 2D)
 *   grpD[2, N_c, nlev]   = {z_ekinh, z_w_con_c_full}               (cell gathered pair)
 *
 * Compile (CPU):
 *   g++ -O3 -fopenmp -march=native -std=c++17 bench_ddt_vn.cpp -o bench_ddt_vn_cpu
 *
 * Compile (GPU via HIP):
 *   hipcc -O3 -std=c++17 bench_ddt_vn.cpp -o bench_ddt_vn_gpu
 *
 * Run:
 *   ./bench_ddt_vn_cpu [nruns] [nlev]
 *   ./bench_ddt_vn_gpu [nruns] [nlev]
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
#define GPU_MALLOC      hipMalloc
#define GPU_FREE        hipFree
#define GPU_MEMCPY_H2D  hipMemcpyHostToDevice
#define GPU_MEMCPY_D2H  hipMemcpyDeviceToHost
#define GPU_MEMCPY      hipMemcpy
#define GPU_SYNC        hipDeviceSynchronize
#define GPU_EVENT_T     hipEvent_t
#define GPU_EVENT_CREATE hipEventCreate
#define GPU_EVENT_RECORD hipEventRecord
#define GPU_EVENT_SYNC   hipEventSynchronize
#define GPU_EVENT_ELAPSED hipEventElapsedTime
#define GPU_EVENT_DESTROY hipEventDestroy
#define HD __host__ __device__ __forceinline__
#define HAS_GPU 1
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
#define GPU_MALLOC      cudaMalloc
#define GPU_FREE        cudaFree
#define GPU_MEMCPY_H2D  cudaMemcpyHostToDevice
#define GPU_MEMCPY_D2H  cudaMemcpyDeviceToHost
#define GPU_MEMCPY      cudaMemcpy
#define GPU_SYNC        cudaDeviceSynchronize
#define GPU_EVENT_T     cudaEvent_t
#define GPU_EVENT_CREATE cudaEventCreate
#define GPU_EVENT_RECORD cudaEventRecord
#define GPU_EVENT_SYNC   cudaEventSynchronize
#define GPU_EVENT_ELAPSED cudaEventElapsedTime
#define GPU_EVENT_DESTROY cudaEventDestroy
#define HD __host__ __device__ __forceinline__
#define HAS_GPU 1
#else
/* CPU-only build */
#define HD inline
#define HAS_GPU 0
#include <omp.h>
#endif

/* ================================================================ */
/*  Constants                                                        */
/* ================================================================ */
static constexpr int NPROMA   = 81920;
static constexpr int WARMUP   = 5;

/* Group field indices */
static constexpr int A_VN     = 0;
static constexpr int A_VT     = 1;
static constexpr int A_VN_IE  = 2;
static constexpr int A_VT_IE  = 3;
static constexpr int A_EKIN   = 4;
static constexpr int NA       = 5;  /* fields in grpA */

static constexpr int C_IDUAL  = 0;
static constexpr int C_IPRIM  = 1;
static constexpr int C_TANG   = 2;
static constexpr int C_FE     = 3;
static constexpr int NC       = 4;  /* fields in grpC */

static constexpr int D_EKINH  = 0;
static constexpr int D_WCON   = 1;
static constexpr int ND       = 2;  /* fields in grpD */

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

/* 2D column-major: je fast, jk slow */
HD int I2(int je, int jk, int N) { return je + jk * N; }

/* Connectivity: 2 neighbors per edge, interleaved */
HD int IX(int je, int n, int N) { return je + n * N; }

/* Group A: field-first AoS, then je, then jk */
HD int IA(int f, int je, int jk, int N) { return f + NA * (je + N * jk); }

/* Group C: field-first AoS, then je (2D) */
HD int IC_grp(int f, int je) { return f + NC * je; }

/* Group D: field-first AoS, then cell, then jk */
HD int ID(int f, int cell, int jk, int N_c) { return f + ND * (cell + N_c * jk); }

/* coeff_gradekin: 2 coefficients per edge */
HD int IG(int je, int n, int N) { return je + n * N; }

/* ================================================================ */
/*  SoA stencil body (baseline — 15 array streams)                   */
/* ================================================================ */

/*
 * Arrays (SoA):
 *   vt[N,nlev], vn_ie[N,nlevp1], z_kin_hor_e[N,nlev]     — edge 3D
 *   ddqz_z_full_e[N,nlev]                                  — edge metric 3D
 *   coeff_gradekin[N,2], c_lin_e[N,2], f_e[N]              — edge 2D
 *   z_ekinh[N_c,nlev], z_w_con_c_full[N_c,nlev]            — cell 3D (indirect)
 *   zeta[N_v,nlev]                                          — vertex 3D (indirect)
 *   cell_idx[N,2], vert_idx[N,2]                            — connectivity
 *   out[N,nlev]                                             — output
 */
#define SOA_STENCIL_BODY()                                                     \
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

/* ================================================================ */
/*  Grouped stencil body (10 array streams)                          */
/* ================================================================ */

/*
 * Arrays (grouped):
 *   grpA[NA, N, nlevp1]              — vt=1, vn_ie=2, z_kin_hor_e=4
 *   grpC[NC, N]                      — f_e=3
 *   grpD[ND, N_c, nlev]              — ekinh=0, wcon_full=1  (indirect)
 *   ddqz_z_full_e[N, nlev]           — standalone edge metric
 *   coeff_gradekin[N, 2]             — standalone (already paired)
 *   c_lin_e[N, 2]                    — standalone (already paired)
 *   zeta[N_v, nlev]                  — standalone vertex (indirect)
 *   cell_idx[N, 2], vert_idx[N, 2]  — connectivity
 *   out[N, nlev]                     — output
 */
#define GRP_STENCIL_BODY()                                                     \
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


/* ================================================================ */
/*  Connectivity distributions                                       */
/* ================================================================ */
enum CellDist { UNIFORM=0, NORMAL1=1, NORMAL4=2, SEQUENTIAL=3 };
static const char *dist_names[] = {"uniform","normal_var1","normal_var4","sequential"};

static void gen_connectivity(int *L, int N, int N_target, CellDist dist, std::mt19937 &rng) {
    switch (dist) {
    case UNIFORM: {
        std::uniform_int_distribution<int> u(0, N_target - 1);
        for (int i = 0; i < N; i++) { L[i*2] = u(rng); L[i*2+1] = u(rng); }
        break;
    }
    case NORMAL1: {
        std::normal_distribution<double> nd(0, 1);
        for (int i = 0; i < N; i++) {
            L[i*2]   = ((i+1 + (int)std::round(nd(rng))) % N_target + N_target) % N_target;
            L[i*2+1] = ((i-1 + (int)std::round(nd(rng))) % N_target + N_target) % N_target;
        }
        break;
    }
    case NORMAL4: {
        std::normal_distribution<double> nd(0, 2);
        for (int i = 0; i < N; i++) {
            L[i*2]   = ((i+1 + (int)std::round(nd(rng))) % N_target + N_target) % N_target;
            L[i*2+1] = ((i-1 + (int)std::round(nd(rng))) % N_target + N_target) % N_target;
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
/*  SoA data                                                         */
/* ================================================================ */
struct SoAData {
    int N, N_c, N_v, nlev;
    size_t sz_e, sz_c, sz_v;  /* N*nlev, N_c*nlev, N_v*nlev */

    double *vt, *vn_ie, *z_kin_hor_e, *ddqz_z_full_e;
    double *coeff_gradekin, *c_lin_e, *f_e;
    double *z_ekinh, *z_w_con_c_full, *zeta;
    int    *cell_idx, *vert_idx;
    double *out;

    void alloc(int N_, int N_c_, int N_v_, int nlev_) {
        N = N_; N_c = N_c_; N_v = N_v_; nlev = nlev_;
        int nlevp1 = nlev + 1;
        sz_e = (size_t)N * nlev; sz_c = (size_t)N_c * nlev; sz_v = (size_t)N_v * nlev;
        vt              = new double[sz_e];
        vn_ie           = new double[(size_t)N * nlevp1];
        z_kin_hor_e     = new double[sz_e];
        ddqz_z_full_e   = new double[sz_e];
        coeff_gradekin  = new double[N * 2];
        c_lin_e         = new double[N * 2];
        f_e             = new double[N];
        z_ekinh         = new double[sz_c];
        z_w_con_c_full  = new double[sz_c];
        zeta            = new double[sz_v];
        cell_idx        = new int[N * 2];
        vert_idx        = new int[N * 2];
        out             = new double[sz_e];
    }
    void fill(std::mt19937 &rng) {
        fill_rand(vt, sz_e, 101); fill_rand(vn_ie, (size_t)N*(nlev+1), 102);
        fill_rand(z_kin_hor_e, sz_e, 103); fill_rand(ddqz_z_full_e, sz_e, 104);
        fill_rand(coeff_gradekin, N*2, 105); fill_rand(c_lin_e, N*2, 106);
        fill_rand(f_e, N, 107);
        fill_rand(z_ekinh, sz_c, 108); fill_rand(z_w_con_c_full, sz_c, 109);
        fill_rand(zeta, sz_v, 110);
        memset(out, 0, sz_e * sizeof(double));
        /* Ensure ddqz nonzero to avoid div-by-zero */
        for (size_t i = 0; i < sz_e; i++)
            if (std::abs(ddqz_z_full_e[i]) < 1e-10) ddqz_z_full_e[i] = 1.0;
    }
    void free_all() {
        delete[] vt; delete[] vn_ie; delete[] z_kin_hor_e; delete[] ddqz_z_full_e;
        delete[] coeff_gradekin; delete[] c_lin_e; delete[] f_e;
        delete[] z_ekinh; delete[] z_w_con_c_full; delete[] zeta;
        delete[] cell_idx; delete[] vert_idx; delete[] out;
    }
};

/* ================================================================ */
/*  Grouped data                                                     */
/* ================================================================ */
struct GrpData {
    int N, N_c, N_v, nlev;
    size_t sz_e;

    double *grpA;             /* NA * N * nlevp1 */
    double *grpC;             /* NC * N */
    double *grpD;             /* ND * N_c * nlev */
    double *ddqz_z_full_e;   /* N * nlev */
    double *coeff_gradekin;   /* N * 2 */
    double *c_lin_e;          /* N * 2 */
    double *zeta;             /* N_v * nlev */
    int    *cell_idx;         /* N * 2 */
    int    *vert_idx;         /* N * 2 */
    double *out;              /* N * nlev */

    void from_soa(const SoAData &s) {
        N = s.N; N_c = s.N_c; N_v = s.N_v; nlev = s.nlev;
        int nlevp1 = nlev + 1;
        sz_e = (size_t)N * nlev;

        grpA = new double[(size_t)NA * N * nlevp1];
        grpC = new double[(size_t)NC * N];
        grpD = new double[(size_t)ND * N_c * nlev];
        ddqz_z_full_e = new double[sz_e];
        coeff_gradekin = new double[N * 2];
        c_lin_e = new double[N * 2];
        zeta = new double[(size_t)N_v * nlev];
        cell_idx = new int[N * 2];
        vert_idx = new int[N * 2];
        out = new double[sz_e];

        /* Pack grpA: for each (je, jk), store {vn=0, vt, vn_ie, z_vt_ie=0, z_kin_hor_e} */
        memset(grpA, 0, (size_t)NA * N * nlevp1 * sizeof(double));
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) {
                grpA[IA(A_VT,       je, jk, N)] = s.vt[I2(je, jk, N)];
                grpA[IA(A_VN_IE,    je, jk, N)] = s.vn_ie[I2(je, jk, N)];
                grpA[IA(A_EKIN,     je, jk, N)] = s.z_kin_hor_e[I2(je, jk, N)];
            }
        /* vn_ie has nlevp1 */
        for (int je = 0; je < N; je++)
            grpA[IA(A_VN_IE, je, nlev, N)] = s.vn_ie[I2(je, nlev, N)];

        /* Pack grpC */
        for (int je = 0; je < N; je++) {
            grpC[IC_grp(C_IDUAL, je)] = 0.0;  /* not used in this stencil */
            grpC[IC_grp(C_IPRIM, je)] = 0.0;
            grpC[IC_grp(C_TANG,  je)] = 0.0;
            grpC[IC_grp(C_FE,    je)] = s.f_e[je];
        }

        /* Pack grpD */
        for (int jk = 0; jk < nlev; jk++)
            for (int ic = 0; ic < N_c; ic++) {
                grpD[ID(D_EKINH, ic, jk, N_c)] = s.z_ekinh[I2(ic, jk, N_c)];
                grpD[ID(D_WCON,  ic, jk, N_c)] = s.z_w_con_c_full[I2(ic, jk, N_c)];
            }

        /* Copy standalone arrays */
        memcpy(ddqz_z_full_e, s.ddqz_z_full_e, sz_e * sizeof(double));
        memcpy(coeff_gradekin, s.coeff_gradekin, N * 2 * sizeof(double));
        memcpy(c_lin_e, s.c_lin_e, N * 2 * sizeof(double));
        memcpy(zeta, s.zeta, (size_t)N_v * nlev * sizeof(double));
        memcpy(cell_idx, s.cell_idx, N * 2 * sizeof(int));
        memcpy(vert_idx, s.vert_idx, N * 2 * sizeof(int));
        memset(out, 0, sz_e * sizeof(double));
    }
    void free_all() {
        delete[] grpA; delete[] grpC; delete[] grpD;
        delete[] ddqz_z_full_e; delete[] coeff_gradekin; delete[] c_lin_e;
        delete[] zeta; delete[] cell_idx; delete[] vert_idx; delete[] out;
    }
};

/* ================================================================ */
/*  CPU kernels                                                      */
/* ================================================================ */
#if !HAS_GPU

static void cpu_soa(SoAData &d) {
    const int N = d.N, N_c = d.N_c, N_v = d.N_v, nlev = d.nlev;
    double *__restrict__ out = d.out;
    const double *__restrict__ vt = d.vt, *__restrict__ vn_ie = d.vn_ie;
    const double *__restrict__ z_kin_hor_e = d.z_kin_hor_e;
    const double *__restrict__ ddqz_z_full_e = d.ddqz_z_full_e;
    const double *__restrict__ coeff_gradekin = d.coeff_gradekin;
    const double *__restrict__ c_lin_e = d.c_lin_e;
    const double *__restrict__ f_e = d.f_e;
    const double *__restrict__ z_ekinh = d.z_ekinh;
    const double *__restrict__ z_w_con_c_full = d.z_w_con_c_full;
    const double *__restrict__ zeta = d.zeta;
    const int    *__restrict__ cell_idx = d.cell_idx;
    const int    *__restrict__ vert_idx = d.vert_idx;

    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) {
            SOA_STENCIL_BODY()
        }
}

static void cpu_grp(GrpData &d) {
    const int N = d.N, N_c = d.N_c, N_v = d.N_v, nlev = d.nlev;
    double *__restrict__ out = d.out;
    const double *__restrict__ grpA = d.grpA;
    const double *__restrict__ grpC = d.grpC;
    const double *__restrict__ grpD = d.grpD;
    const double *__restrict__ ddqz_z_full_e = d.ddqz_z_full_e;
    const double *__restrict__ coeff_gradekin = d.coeff_gradekin;
    const double *__restrict__ c_lin_e = d.c_lin_e;
    const double *__restrict__ zeta = d.zeta;
    const int    *__restrict__ cell_idx = d.cell_idx;
    const int    *__restrict__ vert_idx = d.vert_idx;

    #pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) {
            GRP_STENCIL_BODY()
        }
}

/* Serial reference for verification */
static void cpu_soa_ref(SoAData &d) {
    const int N = d.N, N_c = d.N_c, N_v = d.N_v, nlev = d.nlev;
    double *__restrict__ out = d.out;
    const double *__restrict__ vt = d.vt, *__restrict__ vn_ie = d.vn_ie;
    const double *__restrict__ z_kin_hor_e = d.z_kin_hor_e;
    const double *__restrict__ ddqz_z_full_e = d.ddqz_z_full_e;
    const double *__restrict__ coeff_gradekin = d.coeff_gradekin;
    const double *__restrict__ c_lin_e = d.c_lin_e;
    const double *__restrict__ f_e = d.f_e;
    const double *__restrict__ z_ekinh = d.z_ekinh;
    const double *__restrict__ z_w_con_c_full = d.z_w_con_c_full;
    const double *__restrict__ zeta = d.zeta;
    const int    *__restrict__ cell_idx = d.cell_idx;
    const int    *__restrict__ vert_idx = d.vert_idx;
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++) {
            SOA_STENCIL_BODY()
        }
}

#endif /* !HAS_GPU */

/* ================================================================ */
/*  GPU kernels                                                      */
/* ================================================================ */
#if HAS_GPU

__global__ void gpu_soa_kernel(
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
    int je = blockIdx.x * blockDim.x + threadIdx.x;
    int jk = blockIdx.y * blockDim.y + threadIdx.y;
    if (je >= N || jk >= nlev) return;
    SOA_STENCIL_BODY()
}

__global__ void gpu_grp_kernel(
    double *__restrict__ out,
    const double *__restrict__ grpA, const double *__restrict__ grpC,
    const double *__restrict__ grpD,
    const double *__restrict__ ddqz_z_full_e,
    const double *__restrict__ coeff_gradekin, const double *__restrict__ c_lin_e,
    const double *__restrict__ zeta,
    const int *__restrict__ cell_idx, const int *__restrict__ vert_idx,
    int N, int N_c, int N_v, int nlev)
{
    int je = blockIdx.x * blockDim.x + threadIdx.x;
    int jk = blockIdx.y * blockDim.y + threadIdx.y;
    if (je >= N || jk >= nlev) return;
    GRP_STENCIL_BODY()
}

#endif /* HAS_GPU */

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
/*  main                                                             */
/* ================================================================ */
int main(int argc, char *argv[]) {
    int NRUNS = (argc >= 2) ? atoi(argv[1]) : 20;
    int nlev  = (argc >= 3) ? atoi(argv[2]) : 90;
    int N     = NPROMA;
    int N_c   = N;  /* cells ~ edges for flat synthetic */
    int N_v   = N;  /* verts ~ edges */

    printf("ddt_vn_apc_pc benchmark\n");
    printf("  N=%d  nlev=%d  NRUNS=%d\n", N, nlev, NRUNS);

    FILE *fcsv = fopen("ddt_vn_apc_pc.csv", "w");
    fprintf(fcsv, "layout,cell_dist,run_id,time_ms\n");

    std::mt19937 rng(42);

    for (int di = 0; di < 4; di++) {
        CellDist dist = (CellDist)di;
        printf("\n=== dist=%s ===\n", dist_names[di]);

        /* Allocate and fill SoA data */
        SoAData soa;
        soa.alloc(N, N_c, N_v, nlev);
        soa.fill(rng);
        gen_connectivity(soa.cell_idx, N, N_c, dist, rng);
        gen_connectivity(soa.vert_idx, N, N_v, dist, rng);

        /* Build grouped data from SoA */
        GrpData grp;
        grp.from_soa(soa);

#if !HAS_GPU
        /* ---- CPU path ---- */

        /* Reference */
        double *ref = new double[soa.sz_e];
        memset(ref, 0, soa.sz_e * sizeof(double));
        cpu_soa_ref(soa);
        memcpy(ref, soa.out, soa.sz_e * sizeof(double));

        /* Verify grouped matches SoA */
        cpu_grp(grp);
        int nf; double mr;
        verify(grp.out, ref, soa.sz_e, nf, mr);
        printf("  GRP verify: %s  max_rel=%.2e  fails=%d\n",
               nf ? "FAIL" : "OK", mr, nf);

        /* Warmup */
        for (int r = 0; r < WARMUP; r++) { cpu_soa(soa); cpu_grp(grp); }

        /* Timed SoA */
        for (int r = 0; r < NRUNS; r++) {
            memset(soa.out, 0, soa.sz_e * sizeof(double));
            auto t0 = std::chrono::high_resolution_clock::now();
            cpu_soa(soa);
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            fprintf(fcsv, "soa,%s,%d,%.6f\n", dist_names[di], r, ms);
        }

        /* Timed Grouped */
        for (int r = 0; r < NRUNS; r++) {
            memset(grp.out, 0, grp.sz_e * sizeof(double));
            auto t0 = std::chrono::high_resolution_clock::now();
            cpu_grp(grp);
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            fprintf(fcsv, "grp,%s,%d,%.6f\n", dist_names[di], r, ms);
        }

        /* Compute median times */
        {
            std::vector<double> ts(NRUNS), tg(NRUNS);
            for (int r = 0; r < NRUNS; r++) {
                memset(soa.out, 0, soa.sz_e * sizeof(double));
                auto t0 = std::chrono::high_resolution_clock::now();
                cpu_soa(soa);
                auto t1 = std::chrono::high_resolution_clock::now();
                ts[r] = std::chrono::duration<double, std::milli>(t1 - t0).count();
            }
            for (int r = 0; r < NRUNS; r++) {
                memset(grp.out, 0, grp.sz_e * sizeof(double));
                auto t0 = std::chrono::high_resolution_clock::now();
                cpu_grp(grp);
                auto t1 = std::chrono::high_resolution_clock::now();
                tg[r] = std::chrono::duration<double, std::milli>(t1 - t0).count();
            }
            std::sort(ts.begin(), ts.end());
            std::sort(tg.begin(), tg.end());
            double ms_soa = ts[NRUNS/2], ms_grp = tg[NRUNS/2];
            printf("  SoA median: %.3f ms   GRP median: %.3f ms   speedup: %.2fx\n",
                   ms_soa, ms_grp, ms_soa / ms_grp);
        }

        delete[] ref;

#else /* HAS_GPU */
        /* ---- GPU path ---- */

        /* CPU reference */
        double *h_ref = new double[soa.sz_e];
        memset(soa.out, 0, soa.sz_e * sizeof(double));
        /* Serial reference */
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) {
                const int N_c_ = N_c, N_v_ = N_v;
                int c2 = I2(je, jk, N);
                int ci0 = soa.cell_idx[IX(je,0,N)], ci1 = soa.cell_idx[IX(je,1,N)];
                int vi0 = soa.vert_idx[IX(je,0,N)], vi1 = soa.vert_idx[IX(je,1,N)];
                double cg0 = soa.coeff_gradekin[IG(je,0,N)];
                double cg1 = soa.coeff_gradekin[IG(je,1,N)];
                double cl0 = soa.c_lin_e[IG(je,0,N)], cl1 = soa.c_lin_e[IG(je,1,N)];
                soa.out[c2] = -(
                    soa.z_kin_hor_e[c2] * (cg0 - cg1)
                    + cg1 * soa.z_ekinh[I2(ci1, jk, N_c_)]
                    - cg0 * soa.z_ekinh[I2(ci0, jk, N_c_)]
                    + soa.vt[c2] * (soa.f_e[je] + 0.5 *
                        (soa.zeta[I2(vi0, jk, N_v_)] + soa.zeta[I2(vi1, jk, N_v_)]))
                    + (cl0 * soa.z_w_con_c_full[I2(ci0, jk, N_c_)]
                     + cl1 * soa.z_w_con_c_full[I2(ci1, jk, N_c_)])
                      * (soa.vn_ie[I2(je, jk, N)] - soa.vn_ie[I2(je, jk+1, N)])
                      / soa.ddqz_z_full_e[c2]
                );
            }
        memcpy(h_ref, soa.out, soa.sz_e * sizeof(double));

        /* Allocate device — SoA */
        double *d_vt, *d_vn_ie, *d_ekin_e, *d_ddqz, *d_cg, *d_cl, *d_fe;
        double *d_ekinh, *d_wcon, *d_zeta, *d_out_s;
        int *d_cidx, *d_vidx;
        GPU_CHECK(GPU_MALLOC(&d_vt,     soa.sz_e * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_vn_ie,  (size_t)N*(nlev+1) * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_ekin_e, soa.sz_e * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_ddqz,   soa.sz_e * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_cg,     N*2 * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_cl,     N*2 * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_fe,     N * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_ekinh,  soa.sz_c * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_wcon,   soa.sz_c * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_zeta,   soa.sz_v * sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_cidx,   N*2 * sizeof(int)));
        GPU_CHECK(GPU_MALLOC(&d_vidx,   N*2 * sizeof(int)));
        GPU_CHECK(GPU_MALLOC(&d_out_s,  soa.sz_e * sizeof(double)));

        GPU_CHECK(GPU_MEMCPY(d_vt, soa.vt, soa.sz_e*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_vn_ie, soa.vn_ie, (size_t)N*(nlev+1)*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_ekin_e, soa.z_kin_hor_e, soa.sz_e*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_ddqz, soa.ddqz_z_full_e, soa.sz_e*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_cg, soa.coeff_gradekin, N*2*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_cl, soa.c_lin_e, N*2*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_fe, soa.f_e, N*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_ekinh, soa.z_ekinh, soa.sz_c*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_wcon, soa.z_w_con_c_full, soa.sz_c*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_zeta, soa.zeta, soa.sz_v*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_cidx, soa.cell_idx, N*2*sizeof(int), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_vidx, soa.vert_idx, N*2*sizeof(int), GPU_MEMCPY_H2D));

        /* Allocate device — Grouped */
        double *d_grpA, *d_grpC, *d_grpD, *d_ddqz_g, *d_cg_g, *d_cl_g;
        double *d_zeta_g, *d_out_g;
        int *d_cidx_g, *d_vidx_g;
        GPU_CHECK(GPU_MALLOC(&d_grpA,   (size_t)NA*N*(nlev+1)*sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_grpC,   (size_t)NC*N*sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_grpD,   (size_t)ND*N_c*nlev*sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_ddqz_g, grp.sz_e*sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_cg_g,   N*2*sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_cl_g,   N*2*sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_zeta_g, (size_t)N_v*nlev*sizeof(double)));
        GPU_CHECK(GPU_MALLOC(&d_cidx_g, N*2*sizeof(int)));
        GPU_CHECK(GPU_MALLOC(&d_vidx_g, N*2*sizeof(int)));
        GPU_CHECK(GPU_MALLOC(&d_out_g,  grp.sz_e*sizeof(double)));

        GPU_CHECK(GPU_MEMCPY(d_grpA, grp.grpA, (size_t)NA*N*(nlev+1)*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_grpC, grp.grpC, (size_t)NC*N*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_grpD, grp.grpD, (size_t)ND*N_c*nlev*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_ddqz_g, grp.ddqz_z_full_e, grp.sz_e*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_cg_g, grp.coeff_gradekin, N*2*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_cl_g, grp.c_lin_e, N*2*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_zeta_g, grp.zeta, (size_t)N_v*nlev*sizeof(double), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_cidx_g, grp.cell_idx, N*2*sizeof(int), GPU_MEMCPY_H2D));
        GPU_CHECK(GPU_MEMCPY(d_vidx_g, grp.vert_idx, N*2*sizeof(int), GPU_MEMCPY_H2D));

        dim3 block(64, 4);
        dim3 grid((N + block.x - 1) / block.x, (nlev + block.y - 1) / block.y);

        GPU_EVENT_T ev0, ev1;
        GPU_CHECK(GPU_EVENT_CREATE(&ev0));
        GPU_CHECK(GPU_EVENT_CREATE(&ev1));

        /* Verify SoA GPU */
        gpu_soa_kernel<<<grid, block>>>(d_out_s, d_vt, d_vn_ie, d_ekin_e, d_ddqz,
            d_cg, d_cl, d_fe, d_ekinh, d_wcon, d_zeta, d_cidx, d_vidx,
            N, N_c, N_v, nlev);
        GPU_CHECK(GPU_SYNC());
        double *h_gpu = new double[soa.sz_e];
        GPU_CHECK(GPU_MEMCPY(h_gpu, d_out_s, soa.sz_e*sizeof(double), GPU_MEMCPY_D2H));
        int nf; double mr;
        verify(h_gpu, h_ref, soa.sz_e, nf, mr);
        printf("  GPU SoA verify: %s  max_rel=%.2e\n", nf?"FAIL":"OK", mr);

        /* Verify GRP GPU */
        gpu_grp_kernel<<<grid, block>>>(d_out_g, d_grpA, d_grpC, d_grpD,
            d_ddqz_g, d_cg_g, d_cl_g, d_zeta_g, d_cidx_g, d_vidx_g,
            N, N_c, N_v, nlev);
        GPU_CHECK(GPU_SYNC());
        GPU_CHECK(GPU_MEMCPY(h_gpu, d_out_g, grp.sz_e*sizeof(double), GPU_MEMCPY_D2H));
        verify(h_gpu, h_ref, soa.sz_e, nf, mr);
        printf("  GPU GRP verify: %s  max_rel=%.2e\n", nf?"FAIL":"OK", mr);

        /* Warmup */
        for (int r = 0; r < WARMUP; r++) {
            gpu_soa_kernel<<<grid, block>>>(d_out_s, d_vt, d_vn_ie, d_ekin_e, d_ddqz,
                d_cg, d_cl, d_fe, d_ekinh, d_wcon, d_zeta, d_cidx, d_vidx,
                N, N_c, N_v, nlev);
            gpu_grp_kernel<<<grid, block>>>(d_out_g, d_grpA, d_grpC, d_grpD,
                d_ddqz_g, d_cg_g, d_cl_g, d_zeta_g, d_cidx_g, d_vidx_g,
                N, N_c, N_v, nlev);
        }
        GPU_CHECK(GPU_SYNC());

        /* Timed runs */
        std::vector<double> ts(NRUNS), tg(NRUNS);
        for (int r = 0; r < NRUNS; r++) {
            GPU_CHECK(GPU_EVENT_RECORD(ev0));
            gpu_soa_kernel<<<grid, block>>>(d_out_s, d_vt, d_vn_ie, d_ekin_e, d_ddqz,
                d_cg, d_cl, d_fe, d_ekinh, d_wcon, d_zeta, d_cidx, d_vidx,
                N, N_c, N_v, nlev);
            GPU_CHECK(GPU_EVENT_RECORD(ev1));
            GPU_CHECK(GPU_EVENT_SYNC(ev1));
            float ms; GPU_CHECK(GPU_EVENT_ELAPSED(&ms, ev0, ev1));
            ts[r] = ms;
            fprintf(fcsv, "soa,%s,%d,%.6f\n", dist_names[di], r, (double)ms);
        }
        for (int r = 0; r < NRUNS; r++) {
            GPU_CHECK(GPU_EVENT_RECORD(ev0));
            gpu_grp_kernel<<<grid, block>>>(d_out_g, d_grpA, d_grpC, d_grpD,
                d_ddqz_g, d_cg_g, d_cl_g, d_zeta_g, d_cidx_g, d_vidx_g,
                N, N_c, N_v, nlev);
            GPU_CHECK(GPU_EVENT_RECORD(ev1));
            GPU_CHECK(GPU_EVENT_SYNC(ev1));
            float ms; GPU_CHECK(GPU_EVENT_ELAPSED(&ms, ev0, ev1));
            tg[r] = ms;
            fprintf(fcsv, "grp,%s,%d,%.6f\n", dist_names[di], r, (double)ms);
        }

        std::sort(ts.begin(), ts.end());
        std::sort(tg.begin(), tg.end());
        printf("  GPU SoA median: %.3f ms   GRP median: %.3f ms   speedup: %.2fx\n",
               ts[NRUNS/2], tg[NRUNS/2], ts[NRUNS/2] / tg[NRUNS/2]);

        /* Cleanup GPU */
        GPU_FREE(d_vt); GPU_FREE(d_vn_ie); GPU_FREE(d_ekin_e); GPU_FREE(d_ddqz);
        GPU_FREE(d_cg); GPU_FREE(d_cl); GPU_FREE(d_fe);
        GPU_FREE(d_ekinh); GPU_FREE(d_wcon); GPU_FREE(d_zeta);
        GPU_FREE(d_cidx); GPU_FREE(d_vidx); GPU_FREE(d_out_s);
        GPU_FREE(d_grpA); GPU_FREE(d_grpC); GPU_FREE(d_grpD);
        GPU_FREE(d_ddqz_g); GPU_FREE(d_cg_g); GPU_FREE(d_cl_g);
        GPU_FREE(d_zeta_g); GPU_FREE(d_cidx_g); GPU_FREE(d_vidx_g); GPU_FREE(d_out_g);
        GPU_CHECK(GPU_EVENT_DESTROY(ev0));
        GPU_CHECK(GPU_EVENT_DESTROY(ev1));
        delete[] h_gpu;
        delete[] h_ref;
#endif

        soa.free_all();
        grp.free_all();
        fflush(fcsv);
    }

    fclose(fcsv);
    printf("\nResults: ddt_vn_apc_pc.csv\n");
    return 0;
}