/*
 * bench_ddt_vn_gpu_sweep.cpp -- ddt_vn_apc_pc GPU stencil benchmark
 *
 * Layouts:
 *   soa          — SoA baseline (15 streams)
 *   aos          — Full AoS
 *   grp          — Grouped (grpA[NA=2], grpD[ND=2])
 *   aosoa{32,64} — AoSoA with V=32/64
 *   soa_nf       — SoA nlev-first: arr[je * nlev + jk]
 *   blk{8,64}    — Blocked SoA: arr[N/BX, nlev, BX]
 *   tiled_hf_{64,256} — Tile-NF: arr[tile][jl][jk]
 *
 * Sweeps threadblock, coarsening, connectivity. 100 reps, L2 flush.
 *
 * hipcc -O3 -std=c++17 --offload-arch=gfx942 bench_ddt_vn_gpu_sweep.cpp -o bench
 * nvcc  -O3 -std=c++17 -arch=sm_90  -x cu     bench_ddt_vn_gpu_sweep.cpp -o bench
 */

#include "icon_data_loader.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

/* ================================================================ */
/*  HIP / CUDA portability                                           */
/* ================================================================ */

#ifdef __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>
#define GPU_CHECK(c) do { hipError_t e=(c); if(e!=hipSuccess){ \
    fprintf(stderr,"HIP %s:%d: %s\n",__FILE__,__LINE__,hipGetErrorString(e)); exit(1);} } while(0)
#define GM hipMalloc
#define GF hipFree
#define GC(d,s,n)  GPU_CHECK(hipMemcpy(d,s,n,hipMemcpyHostToDevice))
#define GCD(d,s,n) GPU_CHECK(hipMemcpy(d,s,n,hipMemcpyDeviceToHost))
#define GS hipDeviceSynchronize
#define GE hipEvent_t
#define GEC hipEventCreate
#define GER hipEventRecord
#define GES hipEventSynchronize
#define GEE hipEventElapsedTime
#define GED hipEventDestroy
#define GSET hipMemset
#define HD __host__ __device__ __forceinline__
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#define GPU_CHECK(c) do { cudaError_t e=(c); if(e!=cudaSuccess){ \
    fprintf(stderr,"CUDA %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); exit(1);} } while(0)
#define GM cudaMalloc
#define GF cudaFree
#define GC(d,s,n)  GPU_CHECK(cudaMemcpy(d,s,n,cudaMemcpyHostToDevice))
#define GCD(d,s,n) GPU_CHECK(cudaMemcpy(d,s,n,cudaMemcpyDeviceToHost))
#define GS cudaDeviceSynchronize
#define GE cudaEvent_t
#define GEC cudaEventCreate
#define GER cudaEventRecord
#define GES cudaEventSynchronize
#define GEE cudaEventElapsedTime
#define GED cudaEventDestroy
#define GSET cudaMemset
#define HD __host__ __device__ __forceinline__
#else
#error "Requires HIP or CUDA."
#endif

#define AL(p,n) GPU_CHECK(GM(&(p),(n)))

/* ================================================================ */
/*  Constants                                                        */
/* ================================================================ */

static constexpr int NPROMA = 81920, NRUNS = 100, WARMUP = 10;
static constexpr int A_VT = 0, A_VN_IE = 1, NA = 2;
static constexpr int C_CG0 = 0, C_CG1 = 1, C_CL0 = 2, C_CL1 = 3, NC = 4;
static constexpr int D_EKINH = 0, D_WCON = 1, ND = 2;
static constexpr int E3_VT = 0, E3_DDQZ = 1, NE3 = 2;
static constexpr int E2_CG0 = 0, E2_CG1 = 1, E2_CL0 = 2, E2_CL1 = 3, NE2 = 4;
static constexpr int CN_CI0 = 0, CN_CI1 = 1, CN_VI0 = 2, CN_VI1 = 3, NCONN = 4;

/* ================================================================ */
/*  RNG                                                              */
/* ================================================================ */

static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
static void fill_rand(double *a, size_t n, unsigned s) {
    for (size_t i = 0; i < n; i++) {
        uint64_t h = splitmix64((uint64_t)s * 2654435761ULL + i);
        a[i] = (double)(int64_t)(h & 0xFFFFF) / 100000.0 - 5.0;
    }
}

/* ================================================================ */
/*  Index helpers                                                    */
/* ================================================================ */

/* SoA je-first */
HD int I2(int je, int jk, int N)  { return je + jk * N; }
HD int IX(int je, int n, int N)   { return je + n * N; }
HD int IP(int je, int n)          { return n + 2 * je; }

/* Grouped */
HD int IA(int f, int je, int jk, int N)   { return f + NA * (je + N * jk); }
HD int ID(int f, int c, int jk, int Nc)   { return f + ND * (c + Nc * jk); }

/* Full AoS */
HD int IE3(int f, int je, int jk, int N)  { return f + NE3 * (je + N * jk); }
HD int IE2(int f, int je)                 { return f + NE2 * je; }
HD int ICN(int f, int je)                 { return f + NCONN * je; }
HD int IC(int f, int je)                  { return f + NC * je; }

/* AoSoA */
template <int V> HD int IAao(int f, int je, int jk, int N) {
    int t = (N + V - 1) / V;
    return jk * t * NA * V + (je / V) * NA * V + f * V + (je % V);
}
template <int V> HD int IDao(int f, int c, int jk, int Nc) {
    int t = (Nc + V - 1) / V;
    return jk * t * ND * V + (c / V) * ND * V + f * V + (c % V);
}
template <int V> size_t szAao(int N, int K) {
    return (size_t)((N + V - 1) / V) * NA * V * K;
}
template <int V> size_t szDao(int Nc, int K) {
    return (size_t)((Nc + V - 1) / V) * ND * V * K;
}

/* SoA nlev-first */
HD int I2nf(int je, int jk, int K) { return jk + K * je; }

/* Blocked SoA: arr[N/BX, K, BX] */
template <int BX> HD int I2b(int je, int jk, int K) {
    return (je / BX) * K * BX + jk * BX + (je % BX);
}
template <int BX> size_t szBlk(int M, int K) {
    return (size_t)((M + BX - 1) / BX) * K * BX;
}

/* Tile-NF: arr[N/T, T, K] */
template <int T> HD int I2ht(int je, int jk, int K) {
    return (je / T) * T * K + (je % T) * K + jk;
}
template <int T> size_t szHt(int M, int K) {
    return (size_t)((M + T - 1) / T) * T * K;
}

/* ================================================================ */
/*  Stencil body macros                                              */
/* ================================================================ */

#define STENCIL_EXPR(oi) \
    out[oi] = -(ekin_e * (cg0 - cg1) + cg1 * eh1 - cg0 * eh0 + \
                vt_e * (fe_e + 0.5 * (zeta0 + zeta1)) + \
                (cl0 * w0 + cl1 * w1) * (vk - vk1) / ddqz_e);

#define SOA_BODY() { \
    int c2 = I2(je,jk,N); \
    int ci0=cell_idx[IX(je,0,N)], ci1=cell_idx[IX(je,1,N)], \
        vi0=vert_idx[IX(je,0,N)], vi1=vert_idx[IX(je,1,N)]; \
    double cg0=coeff_gradekin[IX(je,0,N)], cg1=coeff_gradekin[IX(je,1,N)], \
           cl0=c_lin_e[IX(je,0,N)], cl1=c_lin_e[IX(je,1,N)]; \
    double ekin_e=z_kin_hor_e[c2], vt_e=vt[c2], ddqz_e=ddqz_z_full_e[c2], fe_e=f_e[je]; \
    double eh0=z_ekinh[I2(ci0,jk,N_c)], eh1=z_ekinh[I2(ci1,jk,N_c)], \
           w0=z_w_con_c_full[I2(ci0,jk,N_c)], w1=z_w_con_c_full[I2(ci1,jk,N_c)]; \
    double vk=vn_ie[I2(je,jk,N)], vk1=vn_ie[I2(je,jk+1,N)], \
           zeta0=zeta[I2(vi0,jk,N_v)], zeta1=zeta[I2(vi1,jk,N_v)]; \
    STENCIL_EXPR(c2) }

#define AOS_BODY() { \
    int ci0=aos_conn[ICN(CN_CI0,je)], ci1=aos_conn[ICN(CN_CI1,je)], \
        vi0=aos_conn[ICN(CN_VI0,je)], vi1=aos_conn[ICN(CN_VI1,je)]; \
    double cg0=aos_e2d[IE2(E2_CG0,je)], cg1=aos_e2d[IE2(E2_CG1,je)], \
           cl0=aos_e2d[IE2(E2_CL0,je)], cl1=aos_e2d[IE2(E2_CL1,je)]; \
    int eb = NE3*(je+N*jk); \
    double vt_e=aos_e3d[eb+E3_VT], ddqz_e=aos_e3d[eb+E3_DDQZ]; \
    double ekin_e=z_kin_hor_e[I2(je,jk,N)], fe_e=f_e[je]; \
    double eh0=aos_cell[ID(D_EKINH,ci0,jk,N_c)], eh1=aos_cell[ID(D_EKINH,ci1,jk,N_c)], \
           w0=aos_cell[ID(D_WCON,ci0,jk,N_c)], w1=aos_cell[ID(D_WCON,ci1,jk,N_c)]; \
    double vk=vn_ie[I2(je,jk,N)], vk1=vn_ie[I2(je,jk+1,N)], \
           zeta0=zeta[I2(vi0,jk,N_v)], zeta1=zeta[I2(vi1,jk,N_v)]; \
    STENCIL_EXPR(I2(je,jk,N)) }

#define GRP_BODY() { \
    int ci0=conn[ICN(CN_CI0,je)], ci1=conn[ICN(CN_CI1,je)], \
        vi0=conn[ICN(CN_VI0,je)], vi1=conn[ICN(CN_VI1,je)]; \
    double cg0=coef[IC(C_CG0,je)], cg1=coef[IC(C_CG1,je)], \
           cl0=coef[IC(C_CL0,je)], cl1=coef[IC(C_CL1,je)]; \
    double vt_e=grpA[IA(A_VT,je,jk,N)], \
           vk=grpA[IA(A_VN_IE,je,jk,N)], vk1=grpA[IA(A_VN_IE,je,jk+1,N)]; \
    double ekin_e=z_kin_hor_e[I2(je,jk,N)], ddqz_e=ddqz_z_full_e[I2(je,jk,N)], fe_e=f_e[je]; \
    double eh0=grpD[ID(D_EKINH,ci0,jk,N_c)], eh1=grpD[ID(D_EKINH,ci1,jk,N_c)], \
           w0=grpD[ID(D_WCON,ci0,jk,N_c)], w1=grpD[ID(D_WCON,ci1,jk,N_c)]; \
    double zeta0=zeta[I2(vi0,jk,N_v)], zeta1=zeta[I2(vi1,jk,N_v)]; \
    STENCIL_EXPR(I2(je,jk,N)) }

#define AOSOA_BODY(V) { \
    int ci0=conn[ICN(CN_CI0,je)], ci1=conn[ICN(CN_CI1,je)], \
        vi0=conn[ICN(CN_VI0,je)], vi1=conn[ICN(CN_VI1,je)]; \
    double cg0=coef[IC(C_CG0,je)], cg1=coef[IC(C_CG1,je)], \
           cl0=coef[IC(C_CL0,je)], cl1=coef[IC(C_CL1,je)]; \
    double vt_e=grpA[IAao<V>(A_VT,je,jk,N)], \
           vk=grpA[IAao<V>(A_VN_IE,je,jk,N)], vk1=grpA[IAao<V>(A_VN_IE,je,jk+1,N)]; \
    double ekin_e=z_kin_hor_e[I2(je,jk,N)], ddqz_e=ddqz_z_full_e[I2(je,jk,N)], fe_e=f_e[je]; \
    double eh0=grpD[IDao<V>(D_EKINH,ci0,jk,N_c)], eh1=grpD[IDao<V>(D_EKINH,ci1,jk,N_c)], \
           w0=grpD[IDao<V>(D_WCON,ci0,jk,N_c)], w1=grpD[IDao<V>(D_WCON,ci1,jk,N_c)]; \
    double zeta0=zeta[I2(vi0,jk,N_v)], zeta1=zeta[I2(vi1,jk,N_v)]; \
    STENCIL_EXPR(I2(je,jk,N)) }

/* ── New layout bodies ── */

#define SOA_NF_BODY() { \
    int c2 = I2nf(je,jk,nlev); int np = nlev+1; \
    int ci0=cell_idx[IX(je,0,N)], ci1=cell_idx[IX(je,1,N)], \
        vi0=vert_idx[IX(je,0,N)], vi1=vert_idx[IX(je,1,N)]; \
    double cg0=coeff_gradekin[IX(je,0,N)], cg1=coeff_gradekin[IX(je,1,N)], \
           cl0=c_lin_e[IX(je,0,N)], cl1=c_lin_e[IX(je,1,N)]; \
    double ekin_e=z_kin_hor_e[c2], vt_e=vt[c2], ddqz_e=ddqz_z_full_e[c2], fe_e=f_e[je]; \
    double eh0=z_ekinh[I2nf(ci0,jk,nlev)], eh1=z_ekinh[I2nf(ci1,jk,nlev)], \
           w0=z_w_con_c_full[I2nf(ci0,jk,nlev)], w1=z_w_con_c_full[I2nf(ci1,jk,nlev)]; \
    double vk=vn_ie[I2nf(je,jk,np)], vk1=vn_ie[I2nf(je,jk+1,np)], \
           zeta0=zeta[I2nf(vi0,jk,nlev)], zeta1=zeta[I2nf(vi1,jk,nlev)]; \
    STENCIL_EXPR(c2) }

#define BLK_GPU_BODY(BX) { \
    int ci0=cell_idx[IX(je,0,N)], ci1=cell_idx[IX(je,1,N)], \
        vi0=vert_idx[IX(je,0,N)], vi1=vert_idx[IX(je,1,N)]; \
    double cg0=coeff_gradekin[IX(je,0,N)], cg1=coeff_gradekin[IX(je,1,N)], \
           cl0=c_lin_e[IX(je,0,N)], cl1=c_lin_e[IX(je,1,N)]; \
    double ekin_e=z_kin_hor_e[I2b<BX>(je,jk,nlev)], vt_e=vt[I2b<BX>(je,jk,nlev)], \
           ddqz_e=ddqz_z_full_e[I2b<BX>(je,jk,nlev)], fe_e=f_e[je]; \
    double eh0=z_ekinh[I2b<BX>(ci0,jk,nlev)], eh1=z_ekinh[I2b<BX>(ci1,jk,nlev)], \
           w0=z_w_con_c_full[I2b<BX>(ci0,jk,nlev)], w1=z_w_con_c_full[I2b<BX>(ci1,jk,nlev)]; \
    int np_=nlev+1; \
    double vk=vn_ie[I2b<BX>(je,jk,np_)], vk1=vn_ie[I2b<BX>(je,jk+1,np_)], \
           zeta0=zeta[I2b<BX>(vi0,jk,nlev)], zeta1=zeta[I2b<BX>(vi1,jk,nlev)]; \
    STENCIL_EXPR(I2b<BX>(je,jk,nlev)) }

#define TILED_HF_GPU_BODY(T) { \
    int ci0=cell_idx[IX(je,0,N)], ci1=cell_idx[IX(je,1,N)], \
        vi0=vert_idx[IX(je,0,N)], vi1=vert_idx[IX(je,1,N)]; \
    double cg0=coeff_gradekin[IX(je,0,N)], cg1=coeff_gradekin[IX(je,1,N)], \
           cl0=c_lin_e[IX(je,0,N)], cl1=c_lin_e[IX(je,1,N)]; \
    double ekin_e=z_kin_hor_e[I2ht<T>(je,jk,nlev)], vt_e=vt[I2ht<T>(je,jk,nlev)], \
           ddqz_e=ddqz_z_full_e[I2ht<T>(je,jk,nlev)], fe_e=f_e[je]; \
    double eh0=z_ekinh[I2ht<T>(ci0,jk,nlev)], eh1=z_ekinh[I2ht<T>(ci1,jk,nlev)], \
           w0=z_w_con_c_full[I2ht<T>(ci0,jk,nlev)], w1=z_w_con_c_full[I2ht<T>(ci1,jk,nlev)]; \
    int np_=nlev+1; \
    double vk=vn_ie[I2ht<T>(je,jk,np_)], vk1=vn_ie[I2ht<T>(je,jk+1,np_)], \
           zeta0=zeta[I2ht<T>(vi0,jk,nlev)], zeta1=zeta[I2ht<T>(vi1,jk,nlev)]; \
    STENCIL_EXPR(I2ht<T>(je,jk,nlev)) }

/* ================================================================ */
/*  Connectivity                                                     */
/* ================================================================ */

enum CellDist { UNIFORM=0, NORMAL1=1, NORMAL4=2, SEQUENTIAL=3, EXACT=4 };
static const char *dist_names[] = {"uniform","normal_var1","normal_var4","sequential","exact"};

static void gen_conn(int *L, int N, int Nt, CellDist d, std::mt19937 &rng) {
    switch (d) {
    case UNIFORM: { std::uniform_int_distribution<int> u(0,Nt-1);
        for(int i=0;i<N;i++){L[i*2]=u(rng);L[i*2+1]=u(rng);} break; }
    case NORMAL1: { std::normal_distribution<double> nd(0,1);
        for(int i=0;i<N;i++){L[i*2]=((i+1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;
        L[i*2+1]=((i-1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;} break; }
    case NORMAL4: { std::normal_distribution<double> nd(0,2);
        for(int i=0;i<N;i++){L[i*2]=((i+1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;
        L[i*2+1]=((i-1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;} break; }
    case SEQUENTIAL: for(int i=0;i<N;i++){L[i*2]=(i+1)%Nt;L[i*2+1]=(i+1)%Nt;} break;
    case EXACT: break;
    }
}

/* ================================================================ */
/*  L2 flush                                                         */
/* ================================================================ */

static constexpr size_t FLUSH_N = 48ULL * 1024 * 1024;
static double *d_fb = nullptr;
__global__ void fk(double *b, size_t n) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) b[i] = b[i] * 1.00001 + 1e-12;
}
static void fi() { GPU_CHECK(GM(&d_fb, FLUSH_N*8)); GPU_CHECK(GSET(d_fb,0,FLUSH_N*8));
    fk<<<(int)((FLUSH_N+255)/256),256>>>(d_fb,FLUSH_N); GPU_CHECK(GS()); }
static void fl() { fk<<<(int)((FLUSH_N+255)/256),256>>>(d_fb,FLUSH_N); GPU_CHECK(GS()); }
static void fd() { GF(d_fb); }

static bool verify(const double *g, const double *r, size_t n, int &nf_, double &mr) {
    nf_=0; mr=0;
    for(size_t i=0;i<n;i++){double d=std::abs(g[i]-r[i]),dn=std::max(std::abs(r[i]),1e-300),rv=d/dn;
    if(rv>mr)mr=rv; if(d>1e-12+1e-8*std::abs(r[i]))nf_++;} return nf_==0; }

/* ================================================================ */
/*  GPU kernel launch macro                                          */
/* ================================================================ */

#define KL(TX,TY,B) \
    int jb=((int)blockIdx.x*(int)blockDim.x+(int)threadIdx.x)*TX, \
        kb=((int)blockIdx.y*(int)blockDim.y+(int)threadIdx.y)*TY; \
    _Pragma("unroll") for(int dy=0;dy<TY;dy++){int jk=kb+dy;if(jk>=nlev)break; \
    _Pragma("unroll") for(int dx=0;dx<TX;dx++){int je=jb+dx;if(je>=N)break; B }}

/* ================================================================ */
/*  GPU kernel templates                                             */
/* ================================================================ */

/* SoA */
template<int TX,int TY> __global__ void
ks(double*__restrict__ out, const double*__restrict__ vt, const double*__restrict__ vn_ie,
   const double*__restrict__ z_kin_hor_e, const double*__restrict__ ddqz_z_full_e,
   const double*__restrict__ coeff_gradekin, const double*__restrict__ c_lin_e,
   const double*__restrict__ f_e, const double*__restrict__ z_ekinh,
   const double*__restrict__ z_w_con_c_full, const double*__restrict__ zeta,
   const int*__restrict__ cell_idx, const int*__restrict__ vert_idx,
   int N, int N_c, int N_v, int nlev) { KL(TX,TY, SOA_BODY()) }

/* AoS */
template<int TX,int TY> __global__ void
ka(double*__restrict__ out, const double*__restrict__ aos_e3d, const double*__restrict__ vn_ie,
   const double*__restrict__ z_kin_hor_e, const double*__restrict__ f_e,
   const double*__restrict__ aos_e2d, const int*__restrict__ aos_conn,
   const double*__restrict__ aos_cell, const double*__restrict__ zeta,
   int N, int N_c, int N_v, int nlev) { KL(TX,TY, AOS_BODY()) }

/* Grouped */
template<int TX,int TY> __global__ void
kg(double*__restrict__ out, const double*__restrict__ grpA, const double*__restrict__ grpD,
   const double*__restrict__ z_kin_hor_e, const double*__restrict__ ddqz_z_full_e,
   const double*__restrict__ f_e, const double*__restrict__ coef, const double*__restrict__ zeta,
   const int*__restrict__ conn, int N, int N_c, int N_v, int nlev) { KL(TX,TY, GRP_BODY()) }

/* AoSoA */
template<int V,int TX,int TY> __global__ void
ko(double*__restrict__ out, const double*__restrict__ grpA, const double*__restrict__ grpD,
   const double*__restrict__ z_kin_hor_e, const double*__restrict__ ddqz_z_full_e,
   const double*__restrict__ f_e, const double*__restrict__ coef, const double*__restrict__ zeta,
   const int*__restrict__ conn, int N, int N_c, int N_v, int nlev) { KL(TX,TY, AOSOA_BODY(V)) }

/* SoA NF — all 3D arrays use I2nf */
template<int TX,int TY> __global__ void
ks_nf(double*__restrict__ out, const double*__restrict__ vt, const double*__restrict__ vn_ie,
      const double*__restrict__ z_kin_hor_e, const double*__restrict__ ddqz_z_full_e,
      const double*__restrict__ coeff_gradekin, const double*__restrict__ c_lin_e,
      const double*__restrict__ f_e, const double*__restrict__ z_ekinh,
      const double*__restrict__ z_w_con_c_full, const double*__restrict__ zeta,
      const int*__restrict__ cell_idx, const int*__restrict__ vert_idx,
      int N, int N_c, int N_v, int nlev) { KL(TX,TY, SOA_NF_BODY()) }

/* Blocked SoA */
template<int BX,int TX,int TY> __global__ void
ks_blk(double*__restrict__ out, const double*__restrict__ vt, const double*__restrict__ vn_ie,
       const double*__restrict__ z_kin_hor_e, const double*__restrict__ ddqz_z_full_e,
       const double*__restrict__ coeff_gradekin, const double*__restrict__ c_lin_e,
       const double*__restrict__ f_e, const double*__restrict__ z_ekinh,
       const double*__restrict__ z_w_con_c_full, const double*__restrict__ zeta,
       const int*__restrict__ cell_idx, const int*__restrict__ vert_idx,
       int N, int N_c, int N_v, int nlev) { KL(TX,TY, BLK_GPU_BODY(BX)) }

/* Tiled HF */
template<int T,int TX,int TY> __global__ void
ks_ht(double*__restrict__ out, const double*__restrict__ vt, const double*__restrict__ vn_ie,
      const double*__restrict__ z_kin_hor_e, const double*__restrict__ ddqz_z_full_e,
      const double*__restrict__ coeff_gradekin, const double*__restrict__ c_lin_e,
      const double*__restrict__ f_e, const double*__restrict__ z_ekinh,
      const double*__restrict__ z_w_con_c_full, const double*__restrict__ zeta,
      const int*__restrict__ cell_idx, const int*__restrict__ vert_idx,
      int N, int N_c, int N_v, int nlev) { KL(TX,TY, TILED_HF_GPU_BODY(T)) }

/* ================================================================ */
/*  Dispatch macros                                                  */
/* ================================================================ */

#define DT(K, tx, ty, g, b, ...) do { \
    if(tx==1&&ty==1) K<1,1><<<g,b>>>(__VA_ARGS__); \
    else if(tx==2&&ty==1) K<2,1><<<g,b>>>(__VA_ARGS__); \
    else if(tx==4&&ty==1) K<4,1><<<g,b>>>(__VA_ARGS__); \
    else if(tx==1&&ty==2) K<1,2><<<g,b>>>(__VA_ARGS__); \
    else if(tx==1&&ty==4) K<1,4><<<g,b>>>(__VA_ARGS__); \
    else if(tx==2&&ty==2) K<2,2><<<g,b>>>(__VA_ARGS__); \
    else if(tx==4&&ty==2) K<4,2><<<g,b>>>(__VA_ARGS__); \
    else if(tx==2&&ty==4) K<2,4><<<g,b>>>(__VA_ARGS__); \
} while(0)

#define DTV(K, V, tx, ty, g, b, ...) do { \
    if(tx==1&&ty==1) K<V,1,1><<<g,b>>>(__VA_ARGS__); \
    else if(tx==2&&ty==1) K<V,2,1><<<g,b>>>(__VA_ARGS__); \
    else if(tx==4&&ty==1) K<V,4,1><<<g,b>>>(__VA_ARGS__); \
    else if(tx==1&&ty==2) K<V,1,2><<<g,b>>>(__VA_ARGS__); \
    else if(tx==1&&ty==4) K<V,1,4><<<g,b>>>(__VA_ARGS__); \
    else if(tx==2&&ty==2) K<V,2,2><<<g,b>>>(__VA_ARGS__); \
    else if(tx==4&&ty==2) K<V,4,2><<<g,b>>>(__VA_ARGS__); \
    else if(tx==2&&ty==4) K<V,2,4><<<g,b>>>(__VA_ARGS__); \
} while(0)

/* Same as DTV but for ks_blk and ks_ht which take BX/T as first template arg */
#define DT_BLK(BX, tx, ty, g, b, ...) DTV(ks_blk, BX, tx, ty, g, b, __VA_ARGS__)
#define DT_HT(T, tx, ty, g, b, ...)   DTV(ks_ht, T, tx, ty, g, b, __VA_ARGS__)

/* ================================================================ */
/*  Host data structures                                             */
/* ================================================================ */

struct SoAH {
    int N, Nc, Nv, nl; size_t se, sc, sv;
    double *vt, *vni, *ek, *dq, *cg, *cl, *fe, *eh, *wc, *zt; int *ci, *vi; double *out;
    void alloc(int N_,int c,int v,int l) {
        N=N_;Nc=c;Nv=v;nl=l; se=(size_t)N*nl; sc=(size_t)c*nl; sv=(size_t)v*nl;
        vt=new double[se]; vni=new double[(size_t)N*(nl+1)]; ek=new double[se]; dq=new double[se];
        cg=new double[N*2]; cl=new double[N*2]; fe=new double[N];
        eh=new double[sc]; wc=new double[sc]; zt=new double[sv];
        ci=new int[N*2]; vi=new int[N*2]; out=new double[se]; }
    void fill() {
        fill_rand(vt,se,101); fill_rand(vni,(size_t)N*(nl+1),102);
        fill_rand(ek,se,103); fill_rand(dq,se,104);
        fill_rand(cg,N*2,105); fill_rand(cl,N*2,106); fill_rand(fe,N,107);
        fill_rand(eh,sc,108); fill_rand(wc,sc,109); fill_rand(zt,sv,110);
        memset(out,0,se*8);
        for(size_t i=0;i<se;i++) if(std::abs(dq[i])<1e-10) dq[i]=1.0; }
    void free_all() {
        delete[]vt;delete[]vni;delete[]ek;delete[]dq;delete[]cg;delete[]cl;
        delete[]fe;delete[]eh;delete[]wc;delete[]zt;delete[]ci;delete[]vi;delete[]out; }
};

/* Generic repacked host data — same pointer names for all SoA-like layouts */
struct RepackH {
    int N, Nc, Nv, nl; size_t se;
    double *vt, *vni, *ek, *dq, *cg, *cl, *fe, *eh, *wc, *zt; int *ci, *vi; double *out;
    size_t sz_vt, sz_vni, sz_ek, sz_dq, sz_eh, sz_wc, sz_zt, sz_out;

    /* Copy 2D + connectivity from SoAH (unchanged across layouts) */
    void copy_2d(const SoAH &s) {
        N=s.N; Nc=s.Nc; Nv=s.Nv; nl=s.nl; se=s.se;
        cg=new double[N*2]; memcpy(cg,s.cg,N*2*8);
        cl=new double[N*2]; memcpy(cl,s.cl,N*2*8);
        fe=new double[N]; memcpy(fe,s.fe,N*8);
        ci=new int[N*2]; memcpy(ci,s.ci,N*2*4);
        vi=new int[N*2]; memcpy(vi,s.vi,N*2*4);
    }
    void free_all() {
        delete[]vt;delete[]vni;delete[]ek;delete[]dq;delete[]cg;delete[]cl;
        delete[]fe;delete[]eh;delete[]wc;delete[]zt;delete[]ci;delete[]vi;delete[]out; }
};

/* Build NF host data */
static void build_nf(RepackH &h, const SoAH &s) {
    h.copy_2d(s);
    int N=s.N, Nc=s.Nc, Nv=s.Nv, nl=s.nl;
    auto tr = [](const double *src, int M, int K) -> double* {
        double *d = new double[(size_t)M*K];
        for(int jk=0;jk<K;jk++) for(int m=0;m<M;m++) d[jk+K*m]=src[m+jk*M];
        return d; };
    h.vt=tr(s.vt,N,nl); h.vni=tr(s.vni,N,nl+1); h.ek=tr(s.ek,N,nl); h.dq=tr(s.dq,N,nl);
    h.eh=tr(s.eh,Nc,nl); h.wc=tr(s.wc,Nc,nl); h.zt=tr(s.zt,Nv,nl);
    h.sz_vt=h.se; h.sz_vni=(size_t)N*(nl+1); h.sz_ek=h.se; h.sz_dq=h.se;
    h.sz_eh=(size_t)Nc*nl; h.sz_wc=(size_t)Nc*nl; h.sz_zt=(size_t)Nv*nl; h.sz_out=h.se;
    h.out=new double[h.se]; memset(h.out,0,h.se*8);
}

/* Build Blocked SoA host data */
template<int BX>
static void build_blk(RepackH &h, const SoAH &s) {
    h.copy_2d(s);
    int N=s.N, Nc=s.Nc, Nv=s.Nv, nl=s.nl;
    auto pk = [](const double *src, int M, int K) -> double* {
        size_t n = szBlk<BX>(M,K); double *d = new double[n]; memset(d,0,n*8);
        for(int jk=0;jk<K;jk++) for(int m=0;m<M;m++) d[I2b<BX>(m,jk,K)]=src[m+jk*M];
        return d; };
    h.vt=pk(s.vt,N,nl); h.vni=pk(s.vni,N,nl+1); h.ek=pk(s.ek,N,nl); h.dq=pk(s.dq,N,nl);
    h.eh=pk(s.eh,Nc,nl); h.wc=pk(s.wc,Nc,nl); h.zt=pk(s.zt,Nv,nl);
    h.sz_vt=szBlk<BX>(N,nl); h.sz_vni=szBlk<BX>(N,nl+1);
    h.sz_ek=szBlk<BX>(N,nl); h.sz_dq=szBlk<BX>(N,nl);
    h.sz_eh=szBlk<BX>(Nc,nl); h.sz_wc=szBlk<BX>(Nc,nl); h.sz_zt=szBlk<BX>(Nv,nl);
    h.sz_out=szBlk<BX>(N,nl);
    h.out=new double[h.sz_out]; memset(h.out,0,h.sz_out*8);
}

/* Build Tiled HF host data */
template<int T>
static void build_ht(RepackH &h, const SoAH &s) {
    h.copy_2d(s);
    int N=s.N, Nc=s.Nc, Nv=s.Nv, nl=s.nl;
    auto pk = [](const double *src, int M, int K) -> double* {
        size_t n = szHt<T>(M,K); double *d = new double[n]; memset(d,0,n*8);
        for(int jk=0;jk<K;jk++) for(int m=0;m<M;m++) d[I2ht<T>(m,jk,K)]=src[m+jk*M];
        return d; };
    h.vt=pk(s.vt,N,nl); h.vni=pk(s.vni,N,nl+1); h.ek=pk(s.ek,N,nl); h.dq=pk(s.dq,N,nl);
    h.eh=pk(s.eh,Nc,nl); h.wc=pk(s.wc,Nc,nl); h.zt=pk(s.zt,Nv,nl);
    h.sz_vt=szHt<T>(N,nl); h.sz_vni=szHt<T>(N,nl+1);
    h.sz_ek=szHt<T>(N,nl); h.sz_dq=szHt<T>(N,nl);
    h.sz_eh=szHt<T>(Nc,nl); h.sz_wc=szHt<T>(Nc,nl); h.sz_zt=szHt<T>(Nv,nl);
    h.sz_out=szHt<T>(N,nl);
    h.out=new double[h.sz_out]; memset(h.out,0,h.sz_out*8);
}

/* AoS / Grp / AoSoA host structs — same as original */
struct AoSH {
    int N,Nc,Nv,nl; size_t se,s3,s2,sn,sl;
    double *e3,*vni,*ek,*fe,*e2; int *cn; double *cd,*zt,*out;
    void from(const SoAH &s) {
        N=s.N;Nc=s.Nc;Nv=s.Nv;nl=s.nl;se=s.se;
        s3=(size_t)NE3*N*nl; s2=(size_t)NE2*N; sn=(size_t)NCONN*N; sl=(size_t)ND*Nc*nl;
        e3=new double[s3];
        for(int jk=0;jk<nl;jk++) for(int je=0;je<N;je++){
            e3[IE3(E3_VT,je,jk,N)]=s.vt[I2(je,jk,N)]; e3[IE3(E3_DDQZ,je,jk,N)]=s.dq[I2(je,jk,N)];}
        vni=new double[(size_t)N*(nl+1)]; memcpy(vni,s.vni,(size_t)N*(nl+1)*8);
        ek=new double[se]; memcpy(ek,s.ek,se*8);
        fe=new double[N]; memcpy(fe,s.fe,N*8);
        e2=new double[s2];
        for(int je=0;je<N;je++){e2[IE2(E2_CG0,je)]=s.cg[IX(je,0,N)];e2[IE2(E2_CG1,je)]=s.cg[IX(je,1,N)];
            e2[IE2(E2_CL0,je)]=s.cl[IX(je,0,N)];e2[IE2(E2_CL1,je)]=s.cl[IX(je,1,N)];}
        cn=new int[sn];
        for(int je=0;je<N;je++){cn[ICN(CN_CI0,je)]=s.ci[IX(je,0,N)];cn[ICN(CN_CI1,je)]=s.ci[IX(je,1,N)];
            cn[ICN(CN_VI0,je)]=s.vi[IX(je,0,N)];cn[ICN(CN_VI1,je)]=s.vi[IX(je,1,N)];}
        cd=new double[sl];
        for(int jk=0;jk<nl;jk++) for(int ic=0;ic<Nc;ic++){
            cd[ID(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)]; cd[ID(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}
        zt=new double[(size_t)Nv*nl]; memcpy(zt,s.zt,(size_t)Nv*nl*8);
        out=new double[se]; memset(out,0,se*8); }
    void free_all(){delete[]e3;delete[]vni;delete[]ek;delete[]fe;delete[]e2;delete[]cn;delete[]cd;delete[]zt;delete[]out;}
};

struct GrpH {
    int N,Nc,Nv,nl; size_t se;
    double *gA,*gD,*ek,*dq,*fe,*coef,*zt; int *conn; double *out;
    void from(const SoAH &s) {
        N=s.N;Nc=s.Nc;Nv=s.Nv;nl=s.nl;se=s.se; int np=nl+1;
        gA=new double[(size_t)NA*N*np]; memset(gA,0,(size_t)NA*N*np*8);
        for(int jk=0;jk<nl;jk++) for(int je=0;je<N;je++){
            gA[IA(A_VT,je,jk,N)]=s.vt[I2(je,jk,N)]; gA[IA(A_VN_IE,je,jk,N)]=s.vni[I2(je,jk,N)];}
        for(int je=0;je<N;je++) gA[IA(A_VN_IE,je,nl,N)]=s.vni[I2(je,nl,N)];
        gD=new double[(size_t)ND*Nc*nl];
        for(int jk=0;jk<nl;jk++) for(int ic=0;ic<Nc;ic++){
            gD[ID(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)]; gD[ID(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}
        ek=new double[se]; memcpy(ek,s.ek,se*8); dq=new double[se]; memcpy(dq,s.dq,se*8);
        fe=new double[N]; memcpy(fe,s.fe,N*8);
        zt=new double[(size_t)Nv*nl]; memcpy(zt,s.zt,(size_t)Nv*nl*8);
        coef=new double[(size_t)NC*N];
        for(int je=0;je<N;je++){coef[IC(C_CG0,je)]=s.cg[IX(je,0,N)];coef[IC(C_CG1,je)]=s.cg[IX(je,1,N)];
            coef[IC(C_CL0,je)]=s.cl[IX(je,0,N)];coef[IC(C_CL1,je)]=s.cl[IX(je,1,N)];}
        conn=new int[(size_t)NCONN*N];
        for(int je=0;je<N;je++){conn[ICN(CN_CI0,je)]=s.ci[IX(je,0,N)];conn[ICN(CN_CI1,je)]=s.ci[IX(je,1,N)];
            conn[ICN(CN_VI0,je)]=s.vi[IX(je,0,N)];conn[ICN(CN_VI1,je)]=s.vi[IX(je,1,N)];}
        out=new double[se]; memset(out,0,se*8); }
    void free_all(){delete[]gA;delete[]gD;delete[]ek;delete[]dq;delete[]fe;delete[]coef;delete[]zt;delete[]conn;delete[]out;}
};

template<int V> struct AoSoAH {
    int N,Nc,Nv,nl; size_t se,sA,sD;
    double *gA,*gD,*ek,*dq,*fe,*coef,*zt; int *conn; double *out;
    void from(const SoAH &s) {
        N=s.N;Nc=s.Nc;Nv=s.Nv;nl=s.nl;se=s.se; int np=nl+1;
        sA=szAao<V>(N,np); sD=szDao<V>(Nc,nl);
        gA=new double[sA]; memset(gA,0,sA*8);
        for(int jk=0;jk<nl;jk++) for(int je=0;je<N;je++){
            gA[IAao<V>(A_VT,je,jk,N)]=s.vt[I2(je,jk,N)];
            gA[IAao<V>(A_VN_IE,je,jk,N)]=s.vni[I2(je,jk,N)];}
        for(int je=0;je<N;je++) gA[IAao<V>(A_VN_IE,je,nl,N)]=s.vni[I2(je,nl,N)];
        gD=new double[sD]; memset(gD,0,sD*8);
        for(int jk=0;jk<nl;jk++) for(int ic=0;ic<Nc;ic++){
            gD[IDao<V>(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)];
            gD[IDao<V>(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}
        ek=new double[se]; memcpy(ek,s.ek,se*8); dq=new double[se]; memcpy(dq,s.dq,se*8);
        fe=new double[N]; memcpy(fe,s.fe,N*8);
        zt=new double[(size_t)Nv*nl]; memcpy(zt,s.zt,(size_t)Nv*nl*8);
        coef=new double[(size_t)NC*N];
        for(int je=0;je<N;je++){coef[IC(C_CG0,je)]=s.cg[IX(je,0,N)];coef[IC(C_CG1,je)]=s.cg[IX(je,1,N)];
            coef[IC(C_CL0,je)]=s.cl[IX(je,0,N)];coef[IC(C_CL1,je)]=s.cl[IX(je,1,N)];}
        conn=new int[(size_t)NCONN*N];
        for(int je=0;je<N;je++){conn[ICN(CN_CI0,je)]=s.ci[IX(je,0,N)];conn[ICN(CN_CI1,je)]=s.ci[IX(je,1,N)];
            conn[ICN(CN_VI0,je)]=s.vi[IX(je,0,N)];conn[ICN(CN_VI1,je)]=s.vi[IX(je,1,N)];}
        out=new double[se]; memset(out,0,se*8); }
    void free_all(){delete[]gA;delete[]gD;delete[]ek;delete[]dq;delete[]fe;delete[]coef;delete[]zt;delete[]conn;delete[]out;}
};

/* ================================================================ */
/*  Device structs + allocation                                      */
/* ================================================================ */

struct DS { double *vt,*vni,*ek,*dq,*cg,*cl,*fe,*eh,*wc,*zt,*out; int *ci,*vi; };
struct DA { double *e3,*vni,*ek,*fe,*e2,*cd,*zt,*out; int *cn; };
struct DG { double *gA,*gD,*ek,*dq,*fe,*coef,*zt,*out; int *conn; };
struct DO_ { double *gA,*gD,*ek,*dq,*fe,*coef,*zt,*out; int *conn; size_t sA,sD; };

static void ads(DS &d, const SoAH &h) {
    int N=h.N,nl=h.nl;
    AL(d.vt,h.se*8);AL(d.vni,(size_t)N*(nl+1)*8);AL(d.ek,h.se*8);AL(d.dq,h.se*8);
    AL(d.cg,N*2*8);AL(d.cl,N*2*8);AL(d.fe,N*8);AL(d.eh,h.sc*8);AL(d.wc,h.sc*8);
    AL(d.zt,h.sv*8);AL(d.ci,N*2*4);AL(d.vi,N*2*4);AL(d.out,h.se*8);
    GC(d.vt,h.vt,h.se*8);GC(d.vni,h.vni,(size_t)N*(nl+1)*8);GC(d.ek,h.ek,h.se*8);
    GC(d.dq,h.dq,h.se*8);GC(d.cg,h.cg,N*2*8);GC(d.cl,h.cl,N*2*8);GC(d.fe,h.fe,N*8);
    GC(d.eh,h.eh,h.sc*8);GC(d.wc,h.wc,h.sc*8);GC(d.zt,h.zt,h.sv*8);
    GC(d.ci,h.ci,N*2*4);GC(d.vi,h.vi,N*2*4); }
static void fds(DS &d){GF(d.vt);GF(d.vni);GF(d.ek);GF(d.dq);GF(d.cg);GF(d.cl);GF(d.fe);GF(d.eh);GF(d.wc);GF(d.zt);GF(d.ci);GF(d.vi);GF(d.out);}

/* Allocate DS from RepackH (used for NF, blocked, tiled) */
static void ads_rp(DS &d, const RepackH &h) {
    int N=h.N;
    AL(d.vt,h.sz_vt*8);AL(d.vni,h.sz_vni*8);AL(d.ek,h.sz_ek*8);AL(d.dq,h.sz_dq*8);
    AL(d.cg,N*2*8);AL(d.cl,N*2*8);AL(d.fe,N*8);AL(d.eh,h.sz_eh*8);AL(d.wc,h.sz_wc*8);
    AL(d.zt,h.sz_zt*8);AL(d.ci,N*2*4);AL(d.vi,N*2*4);AL(d.out,h.sz_out*8);
    GC(d.vt,h.vt,h.sz_vt*8);GC(d.vni,h.vni,h.sz_vni*8);GC(d.ek,h.ek,h.sz_ek*8);
    GC(d.dq,h.dq,h.sz_dq*8);GC(d.cg,h.cg,N*2*8);GC(d.cl,h.cl,N*2*8);GC(d.fe,h.fe,N*8);
    GC(d.eh,h.eh,h.sz_eh*8);GC(d.wc,h.wc,h.sz_wc*8);GC(d.zt,h.zt,h.sz_zt*8);
    GC(d.ci,h.ci,N*2*4);GC(d.vi,h.vi,N*2*4);GPU_CHECK(GSET(d.out,0,h.sz_out*8)); }

static void ada(DA &d, const AoSH &h) {
    int N=h.N,nl=h.nl;
    AL(d.e3,h.s3*8);AL(d.vni,(size_t)N*(nl+1)*8);AL(d.ek,h.se*8);AL(d.fe,N*8);
    AL(d.e2,h.s2*8);AL(d.cn,h.sn*4);AL(d.cd,h.sl*8);AL(d.zt,(size_t)h.Nv*nl*8);AL(d.out,h.se*8);
    GC(d.e3,h.e3,h.s3*8);GC(d.vni,h.vni,(size_t)N*(nl+1)*8);GC(d.ek,h.ek,h.se*8);
    GC(d.fe,h.fe,N*8);GC(d.e2,h.e2,h.s2*8);GC(d.cn,h.cn,h.sn*4);
    GC(d.cd,h.cd,h.sl*8);GC(d.zt,h.zt,(size_t)h.Nv*nl*8); }
static void fda(DA &d){GF(d.e3);GF(d.vni);GF(d.ek);GF(d.fe);GF(d.e2);GF(d.cn);GF(d.cd);GF(d.zt);GF(d.out);}

static void adg(DG &d, const GrpH &h) {
    int N=h.N,Nc=h.Nc,Nv=h.Nv,nl=h.nl;
    AL(d.gA,(size_t)NA*N*(nl+1)*8);AL(d.gD,(size_t)ND*Nc*nl*8);AL(d.ek,h.se*8);AL(d.dq,h.se*8);
    AL(d.fe,N*8);AL(d.coef,(size_t)NC*N*8);AL(d.zt,(size_t)Nv*nl*8);AL(d.conn,(size_t)NCONN*N*4);AL(d.out,h.se*8);
    GC(d.gA,h.gA,(size_t)NA*N*(nl+1)*8);GC(d.gD,h.gD,(size_t)ND*Nc*nl*8);
    GC(d.ek,h.ek,h.se*8);GC(d.dq,h.dq,h.se*8);GC(d.fe,h.fe,N*8);
    GC(d.coef,h.coef,(size_t)NC*N*8);GC(d.zt,h.zt,(size_t)Nv*nl*8);GC(d.conn,h.conn,(size_t)NCONN*N*4); }
static void fdg(DG &d){GF(d.gA);GF(d.gD);GF(d.ek);GF(d.dq);GF(d.fe);GF(d.coef);GF(d.zt);GF(d.conn);GF(d.out);}

template<int V> static void ado(DO_ &d, const AoSoAH<V> &h) {
    int N=h.N,Nv=h.Nv,nl=h.nl; d.sA=h.sA; d.sD=h.sD;
    AL(d.gA,h.sA*8);AL(d.gD,h.sD*8);AL(d.ek,h.se*8);AL(d.dq,h.se*8);
    AL(d.fe,N*8);AL(d.coef,(size_t)NC*N*8);AL(d.zt,(size_t)Nv*nl*8);AL(d.conn,(size_t)NCONN*N*4);AL(d.out,h.se*8);
    GC(d.gA,h.gA,h.sA*8);GC(d.gD,h.gD,h.sD*8);GC(d.ek,h.ek,h.se*8);GC(d.dq,h.dq,h.se*8);
    GC(d.fe,h.fe,N*8);GC(d.coef,h.coef,(size_t)NC*N*8);GC(d.zt,h.zt,(size_t)Nv*nl*8);
    GC(d.conn,h.conn,(size_t)NCONN*N*4); }
static void fdo(DO_ &d){GF(d.gA);GF(d.gD);GF(d.ek);GF(d.dq);GF(d.fe);GF(d.coef);GF(d.zt);GF(d.conn);GF(d.out);}

/* ================================================================ */
/*  Verification helpers for new layouts                             */
/* ================================================================ */

static void nf_to_soa(const double *nf, double *soa, int N, int nl) {
    for(int je=0;je<N;je++) for(int jk=0;jk<nl;jk++) soa[je+jk*N]=nf[jk+nl*je]; }
template<int BX>
static void blk_to_soa(const double *blk, double *soa, int N, int nl) {
    for(int je=0;je<N;je++) for(int jk=0;jk<nl;jk++) soa[je+jk*N]=blk[I2b<BX>(je,jk,nl)]; }
template<int T>
static void ht_to_soa(const double *ht, double *soa, int N, int nl) {
    for(int je=0;je<N;je++) for(int jk=0;jk<nl;jk++) soa[je+jk*N]=ht[I2ht<T>(je,jk,nl)]; }

/* ================================================================ */
/*  Config sweep + benchmark macro                                   */
/* ================================================================ */

struct KC { int bx, by, tx, ty; };
static const KC cfgs[] = {
    {32,1,1,1},{64,1,1,1},{128,1,1,1},{256,1,1,1},
    {64,2,1,1},{128,2,1,1},{64,4,1,1},{32,4,1,1},{32,8,1,1},
    {64,1,2,1},{128,1,2,1},{64,2,2,1},
    {64,1,4,1},{128,1,4,1},
    {64,1,1,2},{128,1,1,2},{64,2,1,2},
    {64,1,1,4},{128,1,1,4},
    {64,1,2,2},{128,1,2,2},{64,2,2,2},{64,1,4,2},{64,1,2,4}};
static constexpr int NC_ = sizeof(cfgs)/sizeof(cfgs[0]);

#define BENCH(lbl, K_CALL, vw) \
    for(int w=0;w<WARMUP;w++){K_CALL;} GPU_CHECK(GS()); \
    for(int r=0;r<NRUNS;r++){fl();GPU_CHECK(GER(e0));K_CALL;GPU_CHECK(GER(e1)); \
    GPU_CHECK(GES(e1));float ms;GPU_CHECK(GEE(&ms,e0,e1)); \
    fprintf(csv,lbl ",%s,%d,%d,%d,%d," #vw ",%d,%.6f\n",dn,c.bx,c.by,c.tx,c.ty,r,(double)ms);}

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */

int main(int argc, char *argv[]) {
    int nl = (argc >= 2) ? atoi(argv[1]) : 90;
    int N = NPROMA, Nc = N, Nv = N;
    printf("ddt_vn GPU sweep  N=%d nl=%d NRUNS=%d cfgs=%d\n", N, nl, NRUNS, NC_);
    printf("  Layouts: soa, aos, grp, aosoa{32,64}, soa_nf, blk{8,64}, tiled_hf_{64,256}\n");

    FILE *csv = fopen("ddt_vn_apc_pc_gpu_sweep.csv", "w");
    fprintf(csv, "layout,cell_dist,block_x,block_y,coarsen_x,coarsen_y,vec_width,run_id,time_ms\n");
    std::mt19937 rng(42);
    fi();
    GE e0, e1; GPU_CHECK(GEC(&e0)); GPU_CHECK(GEC(&e1));

    /* Load ICON connectivity */
    int icon_step = 9;
    std::string gp = icon_global_path(icon_step), pp = icon_patch_path(icon_step);
    int icon_nproma = icon_read_nproma(gp.c_str());
    if (icon_nproma <= 0) fprintf(stderr, "WARNING: could not read nproma\n");
    printf("Loading ICON: %s (nproma=%d)\n", pp.c_str(), icon_nproma);
    IconEdgeData ied;
    bool have_exact = (icon_nproma > 0) && icon_load_patch(pp.c_str(), icon_nproma, ied);
    if (have_exact)
        printf("ICON: nproma=%d n_edges=%d n_cells=%d n_verts=%d\n",
               ied.nproma, ied.n_edges, ied.n_cells, ied.n_verts);
    int *icon_ci = nullptr, *icon_vi = nullptr;
    if (have_exact) {
        icon_ci = new int[N*2]; icon_vi = new int[N*2];
        for (int i = 0; i < N; i++) {
            int src = i % ied.n_edges;
            icon_ci[i*2]=ied.cell_idx[src*2]%Nc; icon_ci[i*2+1]=ied.cell_idx[src*2+1]%Nc;
            icon_vi[i*2]=ied.vert_idx[src*2]%Nv; icon_vi[i*2+1]=ied.vert_idx[src*2+1]%Nv;
        }
        ied.free_all();
    }
    int n_dists = have_exact ? 5 : 4;

    for (int di = 0; di < n_dists; di++) {
        printf("\n=== dist=%s ===\n", dist_names[di]);

        /* Build SoA baseline */
        SoAH soa; soa.alloc(N, Nc, Nv, nl); soa.fill();
        if (di == EXACT) { memcpy(soa.ci, icon_ci, N*2*4); memcpy(soa.vi, icon_vi, N*2*4); }
        else { gen_conn(soa.ci, N, Nc, (CellDist)di, rng); gen_conn(soa.vi, N, Nv, (CellDist)di, rng); }

        /* Build all layouts */
        AoSH fao; fao.from(soa);
        GrpH grp; grp.from(soa);
        AoSoAH<32> a32; a32.from(soa);
        AoSoAH<64> a64; a64.from(soa);
        RepackH hnf; build_nf(hnf, soa);
        RepackH hb8; build_blk<8>(hb8, soa);
        RepackH hb64; build_blk<64>(hb64, soa);
        RepackH hht64; build_ht<64>(hht64, soa);
        RepackH hht256; build_ht<256>(hht256, soa);

        /* Serial reference */
        double *hr = new double[soa.se];
        memset(soa.out, 0, soa.se * 8);
        for (int jk = 0; jk < nl; jk++)
            for (int je = 0; je < N; je++) {
                const double *vt=soa.vt, *vn_ie=soa.vni, *z_kin_hor_e=soa.ek,
                    *ddqz_z_full_e=soa.dq, *coeff_gradekin=soa.cg, *c_lin_e=soa.cl,
                    *f_e=soa.fe, *z_ekinh=soa.eh, *z_w_con_c_full=soa.wc, *zeta=soa.zt;
                const int *cell_idx=soa.ci, *vert_idx=soa.vi;
                double *out=soa.out; int N_c=Nc, N_v=Nv;
                SOA_BODY()
            }
        memcpy(hr, soa.out, soa.se * 8);

        /* Allocate device */
        DS ds; ads(ds, soa);
        DA da; ada(da, fao);
        DG dg; adg(dg, grp);
        DO_ d32; ado<32>(d32, a32);
        DO_ d64; ado<64>(d64, a64);
        DS dnf; ads_rp(dnf, hnf);
        DS db8; ads_rp(db8, hb8);
        DS db64; ads_rp(db64, hb64);
        DS dht64; ads_rp(dht64, hht64);
        DS dht256; ads_rp(dht256, hht256);

        /* Verify all layouts */
        {
            dim3 bk(64,1), gr((N+63)/64, nl);
            double *ho = new double[soa.se], *ho2 = new double[soa.se];
            int nf; double mr;

            auto vfy = [&](const char *lbl, double *dout, size_t osz) {
                GCD(ho, dout, osz * 8);
                verify(ho, hr, soa.se, nf, mr);
                printf("  %-14s %s %.2e\n", lbl, nf ? "FAIL" : "OK", mr);
            };
            auto vfy_nf = [&](const char *lbl, double *dout) {
                GCD(ho, dout, soa.se * 8);
                nf_to_soa(ho, ho2, N, nl);
                verify(ho2, hr, soa.se, nf, mr);
                printf("  %-14s %s %.2e\n", lbl, nf ? "FAIL" : "OK", mr);
            };

            ks<1,1><<<gr,bk>>>(ds.out,ds.vt,ds.vni,ds.ek,ds.dq,ds.cg,ds.cl,ds.fe,ds.eh,ds.wc,ds.zt,ds.ci,ds.vi,N,Nc,Nv,nl); GPU_CHECK(GS());
            vfy("soa", ds.out, soa.se);

            ka<1,1><<<gr,bk>>>(da.out,da.e3,da.vni,da.ek,da.fe,da.e2,da.cn,da.cd,da.zt,N,Nc,Nv,nl); GPU_CHECK(GS());
            vfy("aos", da.out, soa.se);

            kg<1,1><<<gr,bk>>>(dg.out,dg.gA,dg.gD,dg.ek,dg.dq,dg.fe,dg.coef,dg.zt,dg.conn,N,Nc,Nv,nl); GPU_CHECK(GS());
            vfy("grp", dg.out, soa.se);

            ko<32,1,1><<<gr,bk>>>(d32.out,d32.gA,d32.gD,d32.ek,d32.dq,d32.fe,d32.coef,d32.zt,d32.conn,N,Nc,Nv,nl); GPU_CHECK(GS());
            vfy("aosoa32", d32.out, soa.se);

            ko<64,1,1><<<gr,bk>>>(d64.out,d64.gA,d64.gD,d64.ek,d64.dq,d64.fe,d64.coef,d64.zt,d64.conn,N,Nc,Nv,nl); GPU_CHECK(GS());
            vfy("aosoa64", d64.out, soa.se);

            ks_nf<1,1><<<gr,bk>>>(dnf.out,dnf.vt,dnf.vni,dnf.ek,dnf.dq,dnf.cg,dnf.cl,dnf.fe,dnf.eh,dnf.wc,dnf.zt,dnf.ci,dnf.vi,N,Nc,Nv,nl); GPU_CHECK(GS());
            vfy_nf("soa_nf", dnf.out);

            ks_blk<8,1,1><<<gr,bk>>>(db8.out,db8.vt,db8.vni,db8.ek,db8.dq,db8.cg,db8.cl,db8.fe,db8.eh,db8.wc,db8.zt,db8.ci,db8.vi,N,Nc,Nv,nl); GPU_CHECK(GS());
            { size_t osz=szBlk<8>(N,nl); double *tmp=new double[osz]; GCD(tmp,db8.out,osz*8);
              blk_to_soa<8>(tmp,ho2,N,nl); verify(ho2,hr,soa.se,nf,mr);
              printf("  %-14s %s %.2e\n","blk8",nf?"FAIL":"OK",mr); delete[]tmp; }

            ks_blk<64,1,1><<<gr,bk>>>(db64.out,db64.vt,db64.vni,db64.ek,db64.dq,db64.cg,db64.cl,db64.fe,db64.eh,db64.wc,db64.zt,db64.ci,db64.vi,N,Nc,Nv,nl); GPU_CHECK(GS());
            { size_t osz=szBlk<64>(N,nl); double *tmp=new double[osz]; GCD(tmp,db64.out,osz*8);
              blk_to_soa<64>(tmp,ho2,N,nl); verify(ho2,hr,soa.se,nf,mr);
              printf("  %-14s %s %.2e\n","blk64",nf?"FAIL":"OK",mr); delete[]tmp; }

            ks_ht<64,1,1><<<gr,bk>>>(dht64.out,dht64.vt,dht64.vni,dht64.ek,dht64.dq,dht64.cg,dht64.cl,dht64.fe,dht64.eh,dht64.wc,dht64.zt,dht64.ci,dht64.vi,N,Nc,Nv,nl); GPU_CHECK(GS());
            { size_t osz=szHt<64>(N,nl); double *tmp=new double[osz]; GCD(tmp,dht64.out,osz*8);
              ht_to_soa<64>(tmp,ho2,N,nl); verify(ho2,hr,soa.se,nf,mr);
              printf("  %-14s %s %.2e\n","tiled_hf_64",nf?"FAIL":"OK",mr); delete[]tmp; }

            ks_ht<256,1,1><<<gr,bk>>>(dht256.out,dht256.vt,dht256.vni,dht256.ek,dht256.dq,dht256.cg,dht256.cl,dht256.fe,dht256.eh,dht256.wc,dht256.zt,dht256.ci,dht256.vi,N,Nc,Nv,nl); GPU_CHECK(GS());
            { size_t osz=szHt<256>(N,nl); double *tmp=new double[osz]; GCD(tmp,dht256.out,osz*8);
              ht_to_soa<256>(tmp,ho2,N,nl); verify(ho2,hr,soa.se,nf,mr);
              printf("  %-14s %s %.2e\n","tiled_hf_256",nf?"FAIL":"OK",mr); delete[]tmp; }

            delete[] ho; delete[] ho2;
        }

        /* Config sweep */
        const char *dn = dist_names[di];
        for (int ci = 0; ci < NC_; ci++) {
            auto &c = cfgs[ci];
            printf("  cfg %2d/%d (%3d,%d) c(%d,%d)...", ci+1, NC_, c.bx, c.by, c.tx, c.ty);
            fflush(stdout);
            dim3 bk(c.bx, c.by);
            dim3 gr((N + bk.x*c.tx - 1) / (bk.x*c.tx), (nl + bk.y*c.ty - 1) / (bk.y*c.ty));

            BENCH("soa", DT(ks,c.tx,c.ty,gr,bk,ds.out,ds.vt,ds.vni,ds.ek,ds.dq,ds.cg,ds.cl,ds.fe,ds.eh,ds.wc,ds.zt,ds.ci,ds.vi,N,Nc,Nv,nl), 0)
            BENCH("aos", DT(ka,c.tx,c.ty,gr,bk,da.out,da.e3,da.vni,da.ek,da.fe,da.e2,da.cn,da.cd,da.zt,N,Nc,Nv,nl), 0)
            BENCH("grp", DT(kg,c.tx,c.ty,gr,bk,dg.out,dg.gA,dg.gD,dg.ek,dg.dq,dg.fe,dg.coef,dg.zt,dg.conn,N,Nc,Nv,nl), 0)
            BENCH("aosoa", DTV(ko,32,c.tx,c.ty,gr,bk,d32.out,d32.gA,d32.gD,d32.ek,d32.dq,d32.fe,d32.coef,d32.zt,d32.conn,N,Nc,Nv,nl), 32)
            BENCH("aosoa", DTV(ko,64,c.tx,c.ty,gr,bk,d64.out,d64.gA,d64.gD,d64.ek,d64.dq,d64.fe,d64.coef,d64.zt,d64.conn,N,Nc,Nv,nl), 64)

            BENCH("soa_nf", DT(ks_nf,c.tx,c.ty,gr,bk,dnf.out,dnf.vt,dnf.vni,dnf.ek,dnf.dq,dnf.cg,dnf.cl,dnf.fe,dnf.eh,dnf.wc,dnf.zt,dnf.ci,dnf.vi,N,Nc,Nv,nl), 0)
            BENCH("blk8", DT_BLK(8,c.tx,c.ty,gr,bk,db8.out,db8.vt,db8.vni,db8.ek,db8.dq,db8.cg,db8.cl,db8.fe,db8.eh,db8.wc,db8.zt,db8.ci,db8.vi,N,Nc,Nv,nl), 0)
            BENCH("blk64", DT_BLK(64,c.tx,c.ty,gr,bk,db64.out,db64.vt,db64.vni,db64.ek,db64.dq,db64.cg,db64.cl,db64.fe,db64.eh,db64.wc,db64.zt,db64.ci,db64.vi,N,Nc,Nv,nl), 0)
            BENCH("tiled_hf_64", DT_HT(64,c.tx,c.ty,gr,bk,dht64.out,dht64.vt,dht64.vni,dht64.ek,dht64.dq,dht64.cg,dht64.cl,dht64.fe,dht64.eh,dht64.wc,dht64.zt,dht64.ci,dht64.vi,N,Nc,Nv,nl), 0)
            BENCH("tiled_hf_256", DT_HT(256,c.tx,c.ty,gr,bk,dht256.out,dht256.vt,dht256.vni,dht256.ek,dht256.dq,dht256.cg,dht256.cl,dht256.fe,dht256.eh,dht256.wc,dht256.zt,dht256.ci,dht256.vi,N,Nc,Nv,nl), 0)

            fflush(csv);
            printf(" done\n");
        }

        /* Cleanup */
        fds(ds); fda(da); fdg(dg); fdo(d32); fdo(d64);
        fds(dnf); fds(db8); fds(db64); fds(dht64); fds(dht256);
        soa.free_all(); fao.free_all(); grp.free_all(); a32.free_all(); a64.free_all();
        hnf.free_all(); hb8.free_all(); hb64.free_all(); hht64.free_all(); hht256.free_all();
        delete[] hr;
    }

    if (icon_ci) { delete[] icon_ci; delete[] icon_vi; }
    GPU_CHECK(GED(e0)); GPU_CHECK(GED(e1));
    fd(); fclose(csv);
    printf("\nWritten: ddt_vn_apc_pc_gpu_sweep.csv  (%d dists × %d cfgs × 10 layouts × %d reps)\n",
           n_dists, NC_, NRUNS);
    return 0;
}