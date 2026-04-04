/*
 * bench_ddt_vn_cpu_sweep.cpp -- CPU benchmark for ddt_vn_apc_pc
 *
 * Group design (all sizes are 2/4/8):
 *   grpA[NA=4, N, nlevp1]  = {vn, vt, vn_ie, z_vt_ie}
 *   grpD[ND=2, N_c, nlev]  = {z_ekinh, z_w_con_c_full}
 *   Standalone: z_kin_hor_e, ddqz_z_full_e, f_e, zeta, out
 *   Paired (n-contiguous, IP(je,n)=n+2*je for AoS/grp/AoSoA):
 *     coeff_gradekin[2*N], c_lin_e[2*N], cell_idx[2*N], vert_idx[2*N]
 *
 * Full AoS:
 *   aos_e3d[NE3=2, N, nlev] = {vt, ddqz}
 *   aos_e2d[NE2=4, N]       = {cg0, cg1, cl0, cl1}
 *   aos_conn[NCONN=4, N]    = {ci0, ci1, vi0, vi1}
 *   aos_cell[ND=2, N_c, nl] = {ekinh, wcon}
 *   Standalone: z_kin_hor_e, vn_ie, f_e, zeta, out
 *
 * Layouts: soa, aos, grp, aosoa16  (each also with _nf = nlev-first)
 * Each × collapse(0/1) = 16 kernels per (nlev, dist).
 *
 * Compile: g++ -O3 -fopenmp -march=native -std=c++17 bench_ddt_vn_cpu_sweep.cpp -o bench
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
#include "icon_data_loader.h"

static constexpr int NPROMA=81920, NRUNS=100, WARMUP=10, VW=16;
static constexpr int DEFAULT_NLEVS[]={65,90};
static constexpr int DEFAULT_N_NLEVS=2;

/* Group A: NA=4 (was 5; A_EKIN removed) */
static constexpr int A_VN=0, A_VT=1, A_VN_IE=2, A_VT_IE=3, NA=4;
/* Group D: ND=2 */
static constexpr int D_EKINH=0, D_WCON=1, ND=2;
/* Full-AoS e3d: NE3=2 (z_kin_hor_e standalone) */
static constexpr int E3_VT=0, E3_DDQZ=1, NE3=2;
/* Full-AoS e2d: NE2=4 (f_e standalone) */
static constexpr int E2_CG0=0, E2_CG1=1, E2_CL0=2, E2_CL1=3, NE2=4;
/* Full-AoS conn: NCONN=4 */
static constexpr int CN_CI0=0, CN_CI1=1, CN_VI0=2, CN_VI1=3, NCONN=4;

/* ─── NUMA ────────────────────────────────────────────────────── */
template<typename T>static T*na(size_t n){void*p=mmap(0,n*sizeof(T),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);if(p==MAP_FAILED){perror("mmap");exit(1);}return(T*)p;}
template<typename T>static void nf_(T*p,size_t n){munmap(p,n*sizeof(T));}
/* jk-outer first-touch */
static double*nr2(const double*s,int N,int K){size_t n=(size_t)N*K;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int jk=0;jk<K;jk++)memcpy(d+(size_t)jk*N,s+(size_t)jk*N,N*8);return d;}
static double*nr1(const double*s,size_t n){double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(size_t i=0;i<n;i++)d[i]=s[i];return d;}
template<typename T>static T*nri(const T*s,size_t n){T*d=na<T>(n);_Pragma("omp parallel for schedule(static)")for(size_t i=0;i<n;i++)d[i]=s[i];return d;}
static double*nrA(const double*s,int N,int np){size_t n=(size_t)NA*N*np;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int jk=0;jk<np;jk++){size_t o=(size_t)NA*N*jk;memcpy(d+o,s+o,(size_t)NA*N*8);}return d;}
static double*nrD(const double*s,int Nc,int nl){size_t n=(size_t)ND*Nc*nl;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int jk=0;jk<nl;jk++){size_t o=(size_t)ND*Nc*jk;memcpy(d+o,s+o,(size_t)ND*Nc*8);}return d;}
static double*nre3(const double*s,int N,int nl){size_t n=(size_t)NE3*N*nl;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int jk=0;jk<nl;jk++){size_t o=(size_t)NE3*N*jk;memcpy(d+o,s+o,(size_t)NE3*N*8);}return d;}
/* je-outer first-touch (nlev-first) */
static double*nr2nf(const double*s,int N,int K){size_t n=(size_t)N*K;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int je=0;je<N;je++)memcpy(d+(size_t)je*K,s+(size_t)je*K,K*8);return d;}
static double*nrAnf(const double*s,int N,int np){size_t n=(size_t)NA*np*N;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int je=0;je<N;je++){size_t o=(size_t)NA*np*je;memcpy(d+o,s+o,(size_t)NA*np*8);}return d;}
static double*nrDnf(const double*s,int Nc,int nl){size_t n=(size_t)ND*nl*Nc;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int c=0;c<Nc;c++){size_t o=(size_t)ND*nl*c;memcpy(d+o,s+o,(size_t)ND*nl*8);}return d;}
static double*nre3nf(const double*s,int N,int nl){size_t n=(size_t)NE3*nl*N;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int je=0;je<N;je++){size_t o=(size_t)NE3*nl*je;memcpy(d+o,s+o,(size_t)NE3*nl*8);}return d;}
/* AoSoA NUMA */
static double*nraA(const double*s,int N,int np){int t=(N+VW-1)/VW;size_t n=(size_t)t*NA*VW*np;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int jk=0;jk<np;jk++){size_t o=(size_t)jk*t*NA*VW;memcpy(d+o,s+o,(size_t)t*NA*VW*8);}return d;}
static double*nraD(const double*s,int Nc,int nl){int t=(Nc+VW-1)/VW;size_t n=(size_t)t*ND*VW*nl;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int jk=0;jk<nl;jk++){size_t o=(size_t)jk*t*ND*VW;memcpy(d+o,s+o,(size_t)t*ND*VW*8);}return d;}
static double*nraAnf(const double*s,int N,int np){int t=(N+VW-1)/VW;size_t n=(size_t)t*NA*VW*np;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int ti=0;ti<t;ti++){size_t o=(size_t)ti*NA*VW*np;memcpy(d+o,s+o,(size_t)NA*VW*np*8);}return d;}
static double*nraDnf(const double*s,int Nc,int nl){int t=(Nc+VW-1)/VW;size_t n=(size_t)t*ND*VW*nl;double*d=na<double>(n);_Pragma("omp parallel for schedule(static)")for(int ti=0;ti<t;ti++){size_t o=(size_t)ti*ND*VW*nl;memcpy(d+o,s+o,(size_t)ND*VW*nl*8);}return d;}

/* ─── RNG ─────────────────────────────────────────────────────── */
static inline uint64_t splitmix64(uint64_t x){x+=0x9E3779B97F4A7C15ULL;x=(x^(x>>30))*0xBF58476D1CE4E5B9ULL;x=(x^(x>>27))*0x94D049BB133111EBULL;return x^(x>>31);}
static void fill_rand(double*a,size_t n,unsigned s){for(size_t i=0;i<n;i++){uint64_t h=splitmix64((uint64_t)s*2654435761ULL+i);a[i]=(double)(int64_t)(h&0xFFFFF)/100000.0-5.0;}}

/* ─── Index helpers ───────────────────────────────────────────── */
/* SoA: je-fast, pairs at stride N */
static inline int I2(int je,int jk,int N){return je+jk*N;}
static inline int IX(int je,int n,int N){return je+n*N;}      /* stride-N pairs (SoA) */
/* Pair-contiguous: n-fast (for AoS/grp/AoSoA) */
static inline int IP(int je,int n){return n+2*je;}

/* Grouped je-first */
static inline int IA(int f,int je,int jk,int N){return f+NA*(je+N*jk);}
static inline int ID(int f,int c,int jk,int Nc){return f+ND*(c+Nc*jk);}
/* Full AoS je-first */
static inline int IE3(int f,int je,int jk,int N){return f+NE3*(je+N*jk);}
static inline int IE2(int f,int je){return f+NE2*je;}
static inline int ICN(int f,int je){return f+NCONN*je;}

/* nlev-first: jk-fast */
static inline int I2nf(int je,int jk,int K){return jk+K*je;}
static inline int IAnf(int f,int je,int jk,int K){return f+NA*(jk+K*je);}
static inline int IDnf(int f,int c,int jk,int K){return f+ND*(jk+K*c);}
static inline int IE3nf(int f,int je,int jk,int K){return f+NE3*(jk+K*je);}

/* AoSoA je-first (jk outermost) */
static inline int IAao(int f,int je,int jk,int N){int t=(N+VW-1)/VW;return jk*t*NA*VW+(je/VW)*NA*VW+f*VW+(je%VW);}
static inline int IDao(int f,int c,int jk,int Nc){int t=(Nc+VW-1)/VW;return jk*t*ND*VW+(c/VW)*ND*VW+f*VW+(c%VW);}
/* AoSoA nlev-first (tile outermost) */
static inline int IAaonf(int f,int je,int jk,int K){return(je/VW)*NA*VW*K+jk*NA*VW+f*VW+(je%VW);}
static inline int IDaonf(int f,int c,int jk,int K){return(c/VW)*ND*VW*K+jk*ND*VW+f*VW+(c%VW);}

static inline size_t szAao(int N,int K){return(size_t)((N+VW-1)/VW)*NA*VW*K;}
static inline size_t szDao(int Nc,int K){return(size_t)((Nc+VW-1)/VW)*ND*VW*K;}

/* ─── Stencil body macros ─────────────────────────────────────── */
/* Common tail (shared computation once locals are set up) */
#define STENCIL_EXPR(out_idx) \
    out[out_idx]=-(ekin_e*(cg0-cg1)+cg1*eh1-cg0*eh0\
    +vt_e*(fe_e+0.5*(zeta0+zeta1))\
    +(cl0*w0+cl1*w1)*(vk-vk1)/ddqz_e);

/* ── SoA je-first ── */
#define SOA_BODY() { int c2=I2(je,jk,N);\
    int ci0=cell_idx[IX(je,0,N)],ci1=cell_idx[IX(je,1,N)];\
    int vi0=vert_idx[IX(je,0,N)],vi1=vert_idx[IX(je,1,N)];\
    double cg0=coeff_gradekin[IX(je,0,N)],cg1=coeff_gradekin[IX(je,1,N)];\
    double cl0=c_lin_e[IX(je,0,N)],cl1=c_lin_e[IX(je,1,N)];\
    double ekin_e=z_kin_hor_e[c2],vt_e=vt[c2],ddqz_e=ddqz_z_full_e[c2],fe_e=f_e[je];\
    double eh0=z_ekinh[I2(ci0,jk,N_c)],eh1=z_ekinh[I2(ci1,jk,N_c)];\
    double w0=z_w_con_c_full[I2(ci0,jk,N_c)],w1=z_w_con_c_full[I2(ci1,jk,N_c)];\
    double vk=vn_ie[I2(je,jk,N)],vk1=vn_ie[I2(je,jk+1,N)];\
    double zeta0=zeta[I2(vi0,jk,N_v)],zeta1=zeta[I2(vi1,jk,N_v)];\
    STENCIL_EXPR(c2) }

/* ── SoA nlev-first ── */
#define SOA_NF_BODY() { int c2=I2nf(je,jk,nlev);int np=nlev+1;\
    int ci0=cell_idx[IX(je,0,N)],ci1=cell_idx[IX(je,1,N)];\
    int vi0=vert_idx[IX(je,0,N)],vi1=vert_idx[IX(je,1,N)];\
    double cg0=coeff_gradekin[IX(je,0,N)],cg1=coeff_gradekin[IX(je,1,N)];\
    double cl0=c_lin_e[IX(je,0,N)],cl1=c_lin_e[IX(je,1,N)];\
    double ekin_e=z_kin_hor_e[c2],vt_e=vt[c2],ddqz_e=ddqz_z_full_e[c2],fe_e=f_e[je];\
    double eh0=z_ekinh[I2nf(ci0,jk,nlev)],eh1=z_ekinh[I2nf(ci1,jk,nlev)];\
    double w0=z_w_con_c_full[I2nf(ci0,jk,nlev)],w1=z_w_con_c_full[I2nf(ci1,jk,nlev)];\
    double vk=vn_ie[I2nf(je,jk,np)],vk1=vn_ie[I2nf(je,jk+1,np)];\
    double zeta0=zeta[I2nf(vi0,jk,nlev)],zeta1=zeta[I2nf(vi1,jk,nlev)];\
    STENCIL_EXPR(c2) }

/* ── Full AoS je-first ── */
#define AOS_BODY() { \
    int ci0=aos_conn[ICN(CN_CI0,je)],ci1=aos_conn[ICN(CN_CI1,je)];\
    int vi0=aos_conn[ICN(CN_VI0,je)],vi1=aos_conn[ICN(CN_VI1,je)];\
    double cg0=aos_e2d[IE2(E2_CG0,je)],cg1=aos_e2d[IE2(E2_CG1,je)];\
    double cl0=aos_e2d[IE2(E2_CL0,je)],cl1=aos_e2d[IE2(E2_CL1,je)];\
    int eb=NE3*(je+N*jk); double vt_e=aos_e3d[eb+E3_VT],ddqz_e=aos_e3d[eb+E3_DDQZ];\
    double ekin_e=z_kin_hor_e[I2(je,jk,N)],fe_e=f_e[je];\
    double eh0=aos_cell[ID(D_EKINH,ci0,jk,N_c)],eh1=aos_cell[ID(D_EKINH,ci1,jk,N_c)];\
    double w0=aos_cell[ID(D_WCON,ci0,jk,N_c)],w1=aos_cell[ID(D_WCON,ci1,jk,N_c)];\
    double vk=vn_ie[I2(je,jk,N)],vk1=vn_ie[I2(je,jk+1,N)];\
    double zeta0=zeta[I2(vi0,jk,N_v)],zeta1=zeta[I2(vi1,jk,N_v)];\
    STENCIL_EXPR(I2(je,jk,N)) }

/* ── Full AoS nlev-first ── */
#define AOS_NF_BODY() { int np=nlev+1;\
    int ci0=aos_conn[ICN(CN_CI0,je)],ci1=aos_conn[ICN(CN_CI1,je)];\
    int vi0=aos_conn[ICN(CN_VI0,je)],vi1=aos_conn[ICN(CN_VI1,je)];\
    double cg0=aos_e2d[IE2(E2_CG0,je)],cg1=aos_e2d[IE2(E2_CG1,je)];\
    double cl0=aos_e2d[IE2(E2_CL0,je)],cl1=aos_e2d[IE2(E2_CL1,je)];\
    double vt_e=aos_e3d[IE3nf(E3_VT,je,jk,nlev)],ddqz_e=aos_e3d[IE3nf(E3_DDQZ,je,jk,nlev)];\
    double ekin_e=z_kin_hor_e[I2nf(je,jk,nlev)],fe_e=f_e[je];\
    double eh0=aos_cell[IDnf(D_EKINH,ci0,jk,nlev)],eh1=aos_cell[IDnf(D_EKINH,ci1,jk,nlev)];\
    double w0=aos_cell[IDnf(D_WCON,ci0,jk,nlev)],w1=aos_cell[IDnf(D_WCON,ci1,jk,nlev)];\
    double vk=vn_ie[I2nf(je,jk,np)],vk1=vn_ie[I2nf(je,jk+1,np)];\
    double zeta0=zeta[I2nf(vi0,jk,nlev)],zeta1=zeta[I2nf(vi1,jk,nlev)];\
    STENCIL_EXPR(I2nf(je,jk,nlev)) }

/* ── Grouped je-first (pair-contiguous connectivity via IP) ── */
#define GRP_BODY() { \
    int ci0=cell_idx[IP(je,0)],ci1=cell_idx[IP(je,1)];\
    int vi0=vert_idx[IP(je,0)],vi1=vert_idx[IP(je,1)];\
    double cg0=coeff_gradekin[IP(je,0)],cg1=coeff_gradekin[IP(je,1)];\
    double cl0=c_lin_e[IP(je,0)],cl1=c_lin_e[IP(je,1)];\
    double vt_e=grpA[IA(A_VT,je,jk,N)],vk=grpA[IA(A_VN_IE,je,jk,N)],vk1=grpA[IA(A_VN_IE,je,jk+1,N)];\
    double ekin_e=z_kin_hor_e[I2(je,jk,N)],ddqz_e=ddqz_z_full_e[I2(je,jk,N)],fe_e=f_e[je];\
    double eh0=grpD[ID(D_EKINH,ci0,jk,N_c)],eh1=grpD[ID(D_EKINH,ci1,jk,N_c)];\
    double w0=grpD[ID(D_WCON,ci0,jk,N_c)],w1=grpD[ID(D_WCON,ci1,jk,N_c)];\
    double zeta0=zeta[I2(vi0,jk,N_v)],zeta1=zeta[I2(vi1,jk,N_v)];\
    STENCIL_EXPR(I2(je,jk,N)) }

/* ── Grouped nlev-first ── */
#define GRP_NF_BODY() { int np=nlev+1;\
    int ci0=cell_idx[IP(je,0)],ci1=cell_idx[IP(je,1)];\
    int vi0=vert_idx[IP(je,0)],vi1=vert_idx[IP(je,1)];\
    double cg0=coeff_gradekin[IP(je,0)],cg1=coeff_gradekin[IP(je,1)];\
    double cl0=c_lin_e[IP(je,0)],cl1=c_lin_e[IP(je,1)];\
    double vt_e=grpA[IAnf(A_VT,je,jk,np)],vk=grpA[IAnf(A_VN_IE,je,jk,np)],vk1=grpA[IAnf(A_VN_IE,je,jk+1,np)];\
    double ekin_e=z_kin_hor_e[I2nf(je,jk,nlev)],ddqz_e=ddqz_z_full_e[I2nf(je,jk,nlev)],fe_e=f_e[je];\
    double eh0=grpD[IDnf(D_EKINH,ci0,jk,nlev)],eh1=grpD[IDnf(D_EKINH,ci1,jk,nlev)];\
    double w0=grpD[IDnf(D_WCON,ci0,jk,nlev)],w1=grpD[IDnf(D_WCON,ci1,jk,nlev)];\
    double zeta0=zeta[I2nf(vi0,jk,nlev)],zeta1=zeta[I2nf(vi1,jk,nlev)];\
    STENCIL_EXPR(I2nf(je,jk,nlev)) }

/* ── AoSoA-16 je-first (pair-contiguous via IP) ── */
#define AOSOA_BODY() { \
    int ci0=cell_idx[IP(je,0)],ci1=cell_idx[IP(je,1)];\
    int vi0=vert_idx[IP(je,0)],vi1=vert_idx[IP(je,1)];\
    double cg0=coeff_gradekin[IP(je,0)],cg1=coeff_gradekin[IP(je,1)];\
    double cl0=c_lin_e[IP(je,0)],cl1=c_lin_e[IP(je,1)];\
    double vt_e=grpA[IAao(A_VT,je,jk,N)],vk=grpA[IAao(A_VN_IE,je,jk,N)],vk1=grpA[IAao(A_VN_IE,je,jk+1,N)];\
    double ekin_e=z_kin_hor_e[I2(je,jk,N)],ddqz_e=ddqz_z_full_e[I2(je,jk,N)],fe_e=f_e[je];\
    double eh0=grpD[IDao(D_EKINH,ci0,jk,N_c)],eh1=grpD[IDao(D_EKINH,ci1,jk,N_c)];\
    double w0=grpD[IDao(D_WCON,ci0,jk,N_c)],w1=grpD[IDao(D_WCON,ci1,jk,N_c)];\
    double zeta0=zeta[I2(vi0,jk,N_v)],zeta1=zeta[I2(vi1,jk,N_v)];\
    STENCIL_EXPR(I2(je,jk,N)) }

/* ── AoSoA-16 nlev-first ── */
#define AOSOA_NF_BODY() { int np=nlev+1;\
    int ci0=cell_idx[IP(je,0)],ci1=cell_idx[IP(je,1)];\
    int vi0=vert_idx[IP(je,0)],vi1=vert_idx[IP(je,1)];\
    double cg0=coeff_gradekin[IP(je,0)],cg1=coeff_gradekin[IP(je,1)];\
    double cl0=c_lin_e[IP(je,0)],cl1=c_lin_e[IP(je,1)];\
    double vt_e=grpA[IAaonf(A_VT,je,jk,np)],vk=grpA[IAaonf(A_VN_IE,je,jk,np)],vk1=grpA[IAaonf(A_VN_IE,je,jk+1,np)];\
    double ekin_e=z_kin_hor_e[I2nf(je,jk,nlev)],ddqz_e=ddqz_z_full_e[I2nf(je,jk,nlev)],fe_e=f_e[je];\
    double eh0=grpD[IDaonf(D_EKINH,ci0,jk,nlev)],eh1=grpD[IDaonf(D_EKINH,ci1,jk,nlev)];\
    double w0=grpD[IDaonf(D_WCON,ci0,jk,nlev)],w1=grpD[IDaonf(D_WCON,ci1,jk,nlev)];\
    double zeta0=zeta[I2nf(vi0,jk,nlev)],zeta1=zeta[I2nf(vi1,jk,nlev)];\
    STENCIL_EXPR(I2nf(je,jk,nlev)) }

/* ─── Connectivity ────────────────────────────────────────────── */
enum CellDist{UNIFORM=0,NORMAL1=1,NORMAL4=2,SEQUENTIAL=3,EXACT=4};
static const char*dist_names[]={"uniform","normal_var1","normal_var4","sequential","exact"};
static void gen_conn(int*L,int N,int Nt,CellDist d,std::mt19937&rng){
    switch(d){case UNIFORM:{std::uniform_int_distribution<int>u(0,Nt-1);for(int i=0;i<N;i++){L[i*2]=u(rng);L[i*2+1]=u(rng);}break;}
    case NORMAL1:{std::normal_distribution<double>nd(0,1);for(int i=0;i<N;i++){L[i*2]=((i+1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;L[i*2+1]=((i-1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;}break;}
    case NORMAL4:{std::normal_distribution<double>nd(0,2);for(int i=0;i<N;i++){L[i*2]=((i+1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;L[i*2+1]=((i-1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;}break;}
    case SEQUENTIAL:for(int i=0;i<N;i++){L[i*2]=(i+1)%Nt;L[i*2+1]=(i+1)%Nt;}break;}}

/* ─── Cache flush ─────────────────────────────────────────────── */
static constexpr int FN=8192*4,FS=3;
static double*fb0=nullptr,*fb1=nullptr;
static void flush_init(){size_t n=(size_t)FN*FN;fb0=na<double>(n);fb1=na<double>(n);_Pragma("omp parallel for schedule(static)")for(size_t i=0;i<n;i++){uint64_t h=splitmix64(12345ULL+i);fb0[i]=(double)(h>>11)/(double)(1ULL<<53);fb1[i]=fb0[i];}}
static void flush(){double*A=fb0,*B=fb1;for(int s=0;s<FS;s++){_Pragma("omp parallel for schedule(static)")for(int i=1;i<FN-1;i++)for(int j=1;j<FN-1;j++)B[i*FN+j]=0.25*(A[(i-1)*FN+j]+A[(i+1)*FN+j]+A[i*FN+(j-1)]+A[i*FN+(j+1)]);std::swap(A,B);}printf("  [flush] %.6e\n",A[0]);}
static void flush_destroy(){size_t n=(size_t)FN*FN;nf_(fb0,n);nf_(fb1,n);}
static bool verify(const double*g,const double*r,size_t n,int&nfail,double&mr){nfail=0;mr=0;for(size_t i=0;i<n;i++){double d=std::abs(g[i]-r[i]),dn=std::max(std::abs(r[i]),1e-300),rv=d/dn;if(rv>mr)mr=rv;if(d>1e-12+1e-8*std::abs(r[i]))nfail++;}return nfail==0;}

/* ─── Staging data (heap, not NUMA) ───────────────────────────── */
struct Stg{int N,Nc,Nv,nlev;size_t se,sc,sv;
    double*vt,*vn,*ek,*dq,*cg,*cl,*fe,*eh,*wc,*zt;int*ci,*vi;
    void alloc(int N_,int Nc_,int Nv_,int nl){N=N_;Nc=Nc_;Nv=Nv_;nlev=nl;se=(size_t)N*nl;sc=(size_t)Nc*nl;sv=(size_t)Nv*nl;
        vt=new double[se];vn=new double[(size_t)N*(nl+1)];ek=new double[se];dq=new double[se];
        cg=new double[N*2];cl=new double[N*2];fe=new double[N];eh=new double[sc];wc=new double[sc];zt=new double[sv];ci=new int[N*2];vi=new int[N*2];}
    void fill(){fill_rand(vt,se,101);fill_rand(vn,(size_t)N*(nlev+1),102);fill_rand(ek,se,103);fill_rand(dq,se,104);
        fill_rand(cg,N*2,105);fill_rand(cl,N*2,106);fill_rand(fe,N,107);fill_rand(eh,sc,108);fill_rand(wc,sc,109);fill_rand(zt,sv,110);
        for(size_t i=0;i<se;i++)if(std::abs(dq[i])<1e-10)dq[i]=1.0;}
    void set_conn(const int*c_,const int*v_){memcpy(ci,c_,N*2*4);memcpy(vi,v_,N*2*4);}
    void free_all(){delete[]vt;delete[]vn;delete[]ek;delete[]dq;delete[]cg;delete[]cl;delete[]fe;delete[]eh;delete[]wc;delete[]zt;delete[]ci;delete[]vi;}};

/* ─── Data structs (NUMA-distributed, per layout × dim_order) ── */

struct SoA_D{int N,Nc,Nv,nlev;size_t se;
    double*vt,*vn_ie,*z_kin_hor_e,*ddqz_z_full_e,*coeff_gradekin,*c_lin_e,*f_e,*z_ekinh,*z_w_con_c_full,*zeta;int*cell_idx,*vert_idx;double*out;
    void build(const Stg&s,bool nf_){N=s.N;Nc=s.Nc;Nv=s.Nv;nlev=s.nlev;se=s.se;
        auto r2=[&](const double*src,int M,int K)->double*{
            if(!nf_)return nr2(src,M,K);
            double*t=new double[(size_t)M*K];for(int jk=0;jk<K;jk++)for(int m=0;m<M;m++)t[I2nf(m,jk,K)]=src[I2(m,jk,M)];
            double*d=nr2nf(t,M,K);delete[]t;return d;};
        vt=r2(s.vt,N,nlev);vn_ie=r2(s.vn,N,nlev+1);z_kin_hor_e=r2(s.ek,N,nlev);ddqz_z_full_e=r2(s.dq,N,nlev);
        z_ekinh=r2(s.eh,Nc,nlev);z_w_con_c_full=r2(s.wc,Nc,nlev);zeta=r2(s.zt,Nv,nlev);
        coeff_gradekin=nr1(s.cg,N*2);c_lin_e=nr1(s.cl,N*2);f_e=nr1(s.fe,N);
        cell_idx=nri(s.ci,(size_t)N*2);vert_idx=nri(s.vi,(size_t)N*2);
        out=na<double>(se);_Pragma("omp parallel for schedule(static)")for(size_t i=0;i<se;i++)out[i]=0;}
    void nfree(){nf_(vt,se);nf_(vn_ie,(size_t)N*(nlev+1));nf_(z_kin_hor_e,se);nf_(ddqz_z_full_e,se);
        nf_(z_ekinh,(size_t)Nc*nlev);nf_(z_w_con_c_full,(size_t)Nc*nlev);nf_(zeta,(size_t)Nv*nlev);
        nf_(coeff_gradekin,(size_t)N*2);nf_(c_lin_e,(size_t)N*2);nf_(f_e,(size_t)N);
        nf_(cell_idx,(size_t)N*2);nf_(vert_idx,(size_t)N*2);nf_(out,se);}};

struct AoS_D{int N,Nc,Nv,nlev;size_t se,se3,se2,scn,scl;
    double*aos_e3d,*vn_ie,*z_kin_hor_e,*f_e,*aos_e2d;int*aos_conn;double*aos_cell,*zeta,*out;
    void build(const Stg&s,bool nf_){N=s.N;Nc=s.Nc;Nv=s.Nv;nlev=s.nlev;se=s.se;
        se3=(size_t)NE3*N*nlev;se2=(size_t)NE2*N;scn=(size_t)NCONN*N;scl=(size_t)ND*Nc*nlev;
        /* e3d: {vt, ddqz} */
        double*t3=new double[se3];
        if(!nf_){for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++){t3[IE3(E3_VT,je,jk,N)]=s.vt[I2(je,jk,N)];t3[IE3(E3_DDQZ,je,jk,N)]=s.dq[I2(je,jk,N)];}
            aos_e3d=nre3(t3,N,nlev);}
        else{for(int je=0;je<N;je++)for(int jk=0;jk<nlev;jk++){t3[IE3nf(E3_VT,je,jk,nlev)]=s.vt[I2(je,jk,N)];t3[IE3nf(E3_DDQZ,je,jk,nlev)]=s.dq[I2(je,jk,N)];}
            aos_e3d=nre3nf(t3,N,nlev);}
        delete[]t3;
        /* vn_ie, z_kin_hor_e, f_e standalone (same transpose as SoA) */
        if(!nf_){vn_ie=nr2(s.vn,N,nlev+1);z_kin_hor_e=nr2(s.ek,N,nlev);zeta=nr2(s.zt,Nv,nlev);}
        else{double*t;
            t=new double[(size_t)N*(nlev+1)];for(int jk=0;jk<=nlev;jk++)for(int je=0;je<N;je++)t[I2nf(je,jk,nlev+1)]=s.vn[I2(je,jk,N)];vn_ie=nr2nf(t,N,nlev+1);delete[]t;
            t=new double[se];for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++)t[I2nf(je,jk,nlev)]=s.ek[I2(je,jk,N)];z_kin_hor_e=nr2nf(t,N,nlev);delete[]t;
            t=new double[(size_t)Nv*nlev];for(int jk=0;jk<nlev;jk++)for(int iv=0;iv<Nv;iv++)t[I2nf(iv,jk,nlev)]=s.zt[I2(iv,jk,Nv)];zeta=nr2nf(t,Nv,nlev);delete[]t;}
        f_e=nr1(s.fe,N);
        /* e2d: {cg0,cg1,cl0,cl1} */
        double*t2=new double[se2];for(int je=0;je<N;je++){t2[IE2(E2_CG0,je)]=s.cg[IX(je,0,N)];t2[IE2(E2_CG1,je)]=s.cg[IX(je,1,N)];t2[IE2(E2_CL0,je)]=s.cl[IX(je,0,N)];t2[IE2(E2_CL1,je)]=s.cl[IX(je,1,N)];}
        aos_e2d=nr1(t2,se2);delete[]t2;
        /* conn: {ci0,ci1,vi0,vi1} */
        int*tc=new int[scn];for(int je=0;je<N;je++){tc[ICN(CN_CI0,je)]=s.ci[IX(je,0,N)];tc[ICN(CN_CI1,je)]=s.ci[IX(je,1,N)];tc[ICN(CN_VI0,je)]=s.vi[IX(je,0,N)];tc[ICN(CN_VI1,je)]=s.vi[IX(je,1,N)];}
        aos_conn=nri(tc,scn);delete[]tc;
        /* cell: {ekinh,wcon} */
        double*td=new double[scl];
        if(!nf_){for(int jk=0;jk<nlev;jk++)for(int ic=0;ic<Nc;ic++){td[ID(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)];td[ID(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}aos_cell=nrD(td,Nc,nlev);}
        else{for(int ic=0;ic<Nc;ic++)for(int jk=0;jk<nlev;jk++){td[IDnf(D_EKINH,ic,jk,nlev)]=s.eh[I2(ic,jk,Nc)];td[IDnf(D_WCON,ic,jk,nlev)]=s.wc[I2(ic,jk,Nc)];}aos_cell=nrDnf(td,Nc,nlev);}
        delete[]td;
        out=na<double>(se);_Pragma("omp parallel for schedule(static)")for(size_t i=0;i<se;i++)out[i]=0;}
    void nfree(){nf_(aos_e3d,se3);nf_(vn_ie,(size_t)N*(nlev+1));nf_(z_kin_hor_e,se);nf_(f_e,(size_t)N);
        nf_(aos_e2d,se2);nf_(aos_conn,scn);nf_(aos_cell,scl);nf_(zeta,(size_t)Nv*nlev);nf_(out,se);}};

/* Helper: pack pair-contiguous from stride-N staging */
static void pack_pairs(double*dst,const double*src,int N){
    for(int je=0;je<N;je++){dst[IP(je,0)]=src[IX(je,0,N)];dst[IP(je,1)]=src[IX(je,1,N)];}}
static void pack_pairs_i(int*dst,const int*src,int N){
    for(int je=0;je<N;je++){dst[IP(je,0)]=src[IX(je,0,N)];dst[IP(je,1)]=src[IX(je,1,N)];}}

struct Grp_D{int N,Nc,Nv,nlev;size_t se,szA,szD;
    double*grpA,*grpD,*z_kin_hor_e,*ddqz_z_full_e,*f_e,*coeff_gradekin,*c_lin_e,*zeta;int*cell_idx,*vert_idx;double*out;
    void build(const Stg&s,bool nf_){N=s.N;Nc=s.Nc;Nv=s.Nv;nlev=s.nlev;se=s.se;int np=nlev+1;
        szA=(size_t)NA*N*np;szD=(size_t)ND*Nc*nlev;
        double*tA=new double[szA];memset(tA,0,szA*8);double*tD=new double[szD];
        if(!nf_){
            for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++){tA[IA(A_VT,je,jk,N)]=s.vt[I2(je,jk,N)];tA[IA(A_VN_IE,je,jk,N)]=s.vn[I2(je,jk,N)];}
            for(int je=0;je<N;je++)tA[IA(A_VN_IE,je,nlev,N)]=s.vn[I2(je,nlev,N)];
            for(int jk=0;jk<nlev;jk++)for(int ic=0;ic<Nc;ic++){tD[ID(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)];tD[ID(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}
            grpA=nrA(tA,N,np);grpD=nrD(tD,Nc,nlev);
            z_kin_hor_e=nr2(s.ek,N,nlev);ddqz_z_full_e=nr2(s.dq,N,nlev);zeta=nr2(s.zt,Nv,nlev);
        }else{
            for(int je=0;je<N;je++)for(int jk=0;jk<nlev;jk++){tA[IAnf(A_VT,je,jk,np)]=s.vt[I2(je,jk,N)];tA[IAnf(A_VN_IE,je,jk,np)]=s.vn[I2(je,jk,N)];}
            for(int je=0;je<N;je++)tA[IAnf(A_VN_IE,je,nlev,np)]=s.vn[I2(je,nlev,N)];
            for(int ic=0;ic<Nc;ic++)for(int jk=0;jk<nlev;jk++){tD[IDnf(D_EKINH,ic,jk,nlev)]=s.eh[I2(ic,jk,Nc)];tD[IDnf(D_WCON,ic,jk,nlev)]=s.wc[I2(ic,jk,Nc)];}
            grpA=nrAnf(tA,N,np);grpD=nrDnf(tD,Nc,nlev);
            double*t;
            t=new double[se];for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++)t[I2nf(je,jk,nlev)]=s.ek[I2(je,jk,N)];z_kin_hor_e=nr2nf(t,N,nlev);delete[]t;
            t=new double[se];for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++)t[I2nf(je,jk,nlev)]=s.dq[I2(je,jk,N)];ddqz_z_full_e=nr2nf(t,N,nlev);delete[]t;
            t=new double[(size_t)Nv*nlev];for(int jk=0;jk<nlev;jk++)for(int iv=0;iv<Nv;iv++)t[I2nf(iv,jk,nlev)]=s.zt[I2(iv,jk,Nv)];zeta=nr2nf(t,Nv,nlev);delete[]t;}
        delete[]tA;delete[]tD;
        f_e=nr1(s.fe,N);
        /* Pair-contiguous connectivity & coefficients */
        double*tcg=new double[N*2];pack_pairs(tcg,s.cg,N);coeff_gradekin=nr1(tcg,N*2);delete[]tcg;
        double*tcl=new double[N*2];pack_pairs(tcl,s.cl,N);c_lin_e=nr1(tcl,N*2);delete[]tcl;
        int*tci=new int[N*2];pack_pairs_i(tci,s.ci,N);cell_idx=nri(tci,(size_t)N*2);delete[]tci;
        int*tvi=new int[N*2];pack_pairs_i(tvi,s.vi,N);vert_idx=nri(tvi,(size_t)N*2);delete[]tvi;
        out=na<double>(se);_Pragma("omp parallel for schedule(static)")for(size_t i=0;i<se;i++)out[i]=0;}
    void nfree(){nf_(grpA,szA);nf_(grpD,szD);nf_(z_kin_hor_e,se);nf_(ddqz_z_full_e,se);nf_(f_e,(size_t)N);
        nf_(coeff_gradekin,(size_t)N*2);nf_(c_lin_e,(size_t)N*2);nf_(zeta,(size_t)Nv*nlev);
        nf_(cell_idx,(size_t)N*2);nf_(vert_idx,(size_t)N*2);nf_(out,se);}};

struct AoSoA_D{int N,Nc,Nv,nlev;size_t se,szA,szD;
    double*grpA,*grpD,*z_kin_hor_e,*ddqz_z_full_e,*f_e,*coeff_gradekin,*c_lin_e,*zeta;int*cell_idx,*vert_idx;double*out;
    void build(const Stg&s,bool nf_){N=s.N;Nc=s.Nc;Nv=s.Nv;nlev=s.nlev;se=s.se;int np=nlev+1;
        szA=szAao(N,np);szD=szDao(Nc,nlev);
        double*tA=new double[szA];memset(tA,0,szA*8);double*tD=new double[szD];memset(tD,0,szD*8);
        if(!nf_){
            for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++){tA[IAao(A_VT,je,jk,N)]=s.vt[I2(je,jk,N)];tA[IAao(A_VN_IE,je,jk,N)]=s.vn[I2(je,jk,N)];}
            for(int je=0;je<N;je++)tA[IAao(A_VN_IE,je,nlev,N)]=s.vn[I2(je,nlev,N)];
            for(int jk=0;jk<nlev;jk++)for(int ic=0;ic<Nc;ic++){tD[IDao(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)];tD[IDao(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}
            grpA=nraA(tA,N,np);grpD=nraD(tD,Nc,nlev);
            z_kin_hor_e=nr2(s.ek,N,nlev);ddqz_z_full_e=nr2(s.dq,N,nlev);zeta=nr2(s.zt,Nv,nlev);
        }else{
            for(int je=0;je<N;je++)for(int jk=0;jk<nlev;jk++){tA[IAaonf(A_VT,je,jk,np)]=s.vt[I2(je,jk,N)];tA[IAaonf(A_VN_IE,je,jk,np)]=s.vn[I2(je,jk,N)];}
            for(int je=0;je<N;je++)tA[IAaonf(A_VN_IE,je,nlev,np)]=s.vn[I2(je,nlev,N)];
            for(int ic=0;ic<Nc;ic++)for(int jk=0;jk<nlev;jk++){tD[IDaonf(D_EKINH,ic,jk,nlev)]=s.eh[I2(ic,jk,Nc)];tD[IDaonf(D_WCON,ic,jk,nlev)]=s.wc[I2(ic,jk,Nc)];}
            grpA=nraAnf(tA,N,np);grpD=nraDnf(tD,Nc,nlev);
            double*t;
            t=new double[se];for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++)t[I2nf(je,jk,nlev)]=s.ek[I2(je,jk,N)];z_kin_hor_e=nr2nf(t,N,nlev);delete[]t;
            t=new double[se];for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++)t[I2nf(je,jk,nlev)]=s.dq[I2(je,jk,N)];ddqz_z_full_e=nr2nf(t,N,nlev);delete[]t;
            t=new double[(size_t)Nv*nlev];for(int jk=0;jk<nlev;jk++)for(int iv=0;iv<Nv;iv++)t[I2nf(iv,jk,nlev)]=s.zt[I2(iv,jk,Nv)];zeta=nr2nf(t,Nv,nlev);delete[]t;}
        delete[]tA;delete[]tD;
        f_e=nr1(s.fe,N);
        double*tcg=new double[N*2];pack_pairs(tcg,s.cg,N);coeff_gradekin=nr1(tcg,N*2);delete[]tcg;
        double*tcl=new double[N*2];pack_pairs(tcl,s.cl,N);c_lin_e=nr1(tcl,N*2);delete[]tcl;
        int*tci=new int[N*2];pack_pairs_i(tci,s.ci,N);cell_idx=nri(tci,(size_t)N*2);delete[]tci;
        int*tvi=new int[N*2];pack_pairs_i(tvi,s.vi,N);vert_idx=nri(tvi,(size_t)N*2);delete[]tvi;
        out=na<double>(se);_Pragma("omp parallel for schedule(static)")for(size_t i=0;i<se;i++)out[i]=0;}
    void nfree(){nf_(grpA,szA);nf_(grpD,szD);nf_(z_kin_hor_e,se);nf_(ddqz_z_full_e,se);nf_(f_e,(size_t)N);
        nf_(coeff_gradekin,(size_t)N*2);nf_(c_lin_e,(size_t)N*2);nf_(zeta,(size_t)Nv*nlev);
        nf_(cell_idx,(size_t)N*2);nf_(vert_idx,(size_t)N*2);nf_(out,se);}};

/* ─── Pointer extraction ──────────────────────────────────────── */
#define PS(d) const int N=(d).N,N_c=(d).Nc,N_v=(d).Nv,nlev=(d).nlev;double*__restrict__ out=(d).out;\
    const double*__restrict__ vt=(d).vt,*__restrict__ vn_ie=(d).vn_ie,*__restrict__ z_kin_hor_e=(d).z_kin_hor_e;\
    const double*__restrict__ ddqz_z_full_e=(d).ddqz_z_full_e,*__restrict__ coeff_gradekin=(d).coeff_gradekin;\
    const double*__restrict__ c_lin_e=(d).c_lin_e,*__restrict__ f_e=(d).f_e,*__restrict__ z_ekinh=(d).z_ekinh;\
    const double*__restrict__ z_w_con_c_full=(d).z_w_con_c_full,*__restrict__ zeta=(d).zeta;\
    const int*__restrict__ cell_idx=(d).cell_idx,*__restrict__ vert_idx=(d).vert_idx;
#define PA(d) const int N=(d).N,N_c=(d).Nc,N_v=(d).Nv,nlev=(d).nlev;double*__restrict__ out=(d).out;\
    const double*__restrict__ aos_e3d=(d).aos_e3d,*__restrict__ vn_ie=(d).vn_ie;\
    const double*__restrict__ z_kin_hor_e=(d).z_kin_hor_e,*__restrict__ f_e=(d).f_e;\
    const double*__restrict__ aos_e2d=(d).aos_e2d;const int*__restrict__ aos_conn=(d).aos_conn;\
    const double*__restrict__ aos_cell=(d).aos_cell,*__restrict__ zeta=(d).zeta;
#define PG(d) const int N=(d).N,N_c=(d).Nc,N_v=(d).Nv,nlev=(d).nlev;double*__restrict__ out=(d).out;\
    const double*__restrict__ grpA=(d).grpA,*__restrict__ grpD=(d).grpD;\
    const double*__restrict__ z_kin_hor_e=(d).z_kin_hor_e,*__restrict__ ddqz_z_full_e=(d).ddqz_z_full_e;\
    const double*__restrict__ coeff_gradekin=(d).coeff_gradekin,*__restrict__ c_lin_e=(d).c_lin_e;\
    const double*__restrict__ f_e=(d).f_e,*__restrict__ zeta=(d).zeta;\
    const int*__restrict__ cell_idx=(d).cell_idx,*__restrict__ vert_idx=(d).vert_idx;
#define PO PG /* AoSoA uses same member names as Grp */

/* ─── Kernels ─────────────────────────────────────────────────── */
#define JK_O(B) for(int jk=0;jk<nlev;jk++)for(int je=0;je<N;je++)B
#define JE_O(B) for(int je=0;je<N;je++)for(int jk=0;jk<nlev;jk++)B

static void ks(SoA_D&d){PS(d)_Pragma("omp parallel for schedule(static)")JK_O(SOA_BODY())}
static void ks_c(SoA_D&d){PS(d)_Pragma("omp parallel for collapse(2) schedule(static)")JK_O(SOA_BODY())}
static void ksn(SoA_D&d){PS(d)_Pragma("omp parallel for schedule(static)")JE_O(SOA_NF_BODY())}
static void ksn_c(SoA_D&d){PS(d)_Pragma("omp parallel for collapse(2) schedule(static)")JE_O(SOA_NF_BODY())}
static void ka(AoS_D&d){PA(d)_Pragma("omp parallel for schedule(static)")JK_O(AOS_BODY())}
static void ka_c(AoS_D&d){PA(d)_Pragma("omp parallel for collapse(2) schedule(static)")JK_O(AOS_BODY())}
static void kan(AoS_D&d){PA(d)_Pragma("omp parallel for schedule(static)")JE_O(AOS_NF_BODY())}
static void kan_c(AoS_D&d){PA(d)_Pragma("omp parallel for collapse(2) schedule(static)")JE_O(AOS_NF_BODY())}
static void kg(Grp_D&d){PG(d)_Pragma("omp parallel for schedule(static)")JK_O(GRP_BODY())}
static void kg_c(Grp_D&d){PG(d)_Pragma("omp parallel for collapse(2) schedule(static)")JK_O(GRP_BODY())}
static void kgn(Grp_D&d){PG(d)_Pragma("omp parallel for schedule(static)")JE_O(GRP_NF_BODY())}
static void kgn_c(Grp_D&d){PG(d)_Pragma("omp parallel for collapse(2) schedule(static)")JE_O(GRP_NF_BODY())}
static void ko(AoSoA_D&d){PO(d)_Pragma("omp parallel for schedule(static)")JK_O(AOSOA_BODY())}
static void ko_c(AoSoA_D&d){PO(d)_Pragma("omp parallel for collapse(2) schedule(static)")JK_O(AOSOA_BODY())}
static void kon(AoSoA_D&d){PO(d)_Pragma("omp parallel for schedule(static)")JE_O(AOSOA_NF_BODY())}
static void kon_c(AoSoA_D&d){PO(d)_Pragma("omp parallel for collapse(2) schedule(static)")JE_O(AOSOA_NF_BODY())}

static void kref(SoA_D&d,double*ro){const int N=d.N,N_c=d.Nc,N_v=d.Nv,nlev=d.nlev;double*__restrict__ out=ro;
    const double*__restrict__ vt=d.vt,*__restrict__ vn_ie=d.vn_ie,*__restrict__ z_kin_hor_e=d.z_kin_hor_e;
    const double*__restrict__ ddqz_z_full_e=d.ddqz_z_full_e,*__restrict__ coeff_gradekin=d.coeff_gradekin;
    const double*__restrict__ c_lin_e=d.c_lin_e,*__restrict__ f_e=d.f_e,*__restrict__ z_ekinh=d.z_ekinh;
    const double*__restrict__ z_w_con_c_full=d.z_w_con_c_full,*__restrict__ zeta=d.zeta;
    const int*__restrict__ cell_idx=d.cell_idx,*__restrict__ vert_idx=d.vert_idx;
    JK_O(SOA_BODY())}

static void transpose_nf(const double*nfp,double*jf,int N,int nlev){
    for(int je=0;je<N;je++)for(int jk=0;jk<nlev;jk++)jf[I2(je,jk,N)]=nfp[I2nf(je,jk,nlev)];}

/* ─── Benchmark runner ────────────────────────────────────────── */
static void trun(FILE*csv,const char*lay,int col,const char*dn,int nlev,int N,
    void(*fn)(),double*outp,const double*ref,size_t se,bool isnf){
    fn();int nfail;double mr;
    if(!isnf)verify(outp,ref,se,nfail,mr);
    else{double*t=new double[se];transpose_nf(outp,t,N,nlev);verify(t,ref,se,nfail,mr);delete[]t;}
    printf("  %-12s c=%d %-12s %s %.2e\n",lay,col,dn,nfail?"FAIL":"OK",mr);
    for(int w=0;w<WARMUP;w++){flush();fn();}
    for(int r=0;r<NRUNS;r++){flush();auto t0=std::chrono::high_resolution_clock::now();fn();
        auto t1=std::chrono::high_resolution_clock::now();double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
        fprintf(csv,"%s,%d,%d,%s,%d,%d,%.9f\n",lay,nlev,N,dn,col,r,ms);}fflush(csv);}

static SoA_D*_s,*_sn;static AoS_D*_a,*_an;static Grp_D*_g,*_gn;static AoSoA_D*_o,*_on;
static void ws(){ks(*_s);}static void wsc(){ks_c(*_s);}static void wsn(){ksn(*_sn);}static void wsnc(){ksn_c(*_sn);}
static void wa(){ka(*_a);}static void wac(){ka_c(*_a);}static void wan(){kan(*_an);}static void wanc(){kan_c(*_an);}
static void wg(){kg(*_g);}static void wgc(){kg_c(*_g);}static void wgn(){kgn(*_gn);}static void wgnc(){kgn_c(*_gn);}
static void wo(){ko(*_o);}static void woc(){ko_c(*_o);}static void won(){kon(*_on);}static void wonc(){kon_c(*_on);}

int main(int argc,char*argv[]){
    std::vector<int>nlevs;if(argc>=2)nlevs.push_back(atoi(argv[1]));
    else for(int i=0;i<DEFAULT_N_NLEVS;i++)nlevs.push_back(DEFAULT_NLEVS[i]);
    const int N=NPROMA,Nc=N,Nv=N;
    printf("ddt_vn CPU sweep  N=%d thr=%d NRUNS=%d VW=%d NA=%d ND=%d NE3=%d NE2=%d\n",N,omp_get_max_threads(),NRUNS,VW,NA,ND,NE3,NE2);
    FILE*csv=fopen("ddt_vn_apc_pc_cpu_sweep.csv","w");
    fprintf(csv,"layout,nlev,nproma,cell_dist,omp_collapse,run_id,time_ms\n");
    std::mt19937 rng(42);flush_init();srand((unsigned)time(NULL));flush();printf("Ready\n\n");

    /* ── Load ICON exact connectivity ─────────────────────────────── */
    int icon_step = 9;
    std::string gp = icon_global_path(icon_step);
    std::string pp = icon_patch_path(icon_step);
    int icon_nproma = icon_read_nproma(gp.c_str());
    if (icon_nproma <= 0)
        fprintf(stderr, "WARNING: could not read nproma from '%s'\n", gp.c_str());
    printf("Loading ICON: %s  (nproma=%d)\n", pp.c_str(), icon_nproma);
    IconEdgeData ied;
    bool have_exact = (icon_nproma > 0) &&
        icon_load_patch(pp.c_str(), icon_nproma, ied);
    if (have_exact)
        printf("ICON: nproma=%d  n_edges=%d (valid=%d)  n_cells=%d  n_verts=%d\n",
               ied.nproma, ied.n_edges, ied.n_edges_valid, ied.n_cells, ied.n_verts);
    else
        printf("ICON data not available — skipping 'exact' distribution\n");

    /* Prepare remapped ICON connectivity for the benchmark's [0,N) range */
    int *icon_ci = nullptr, *icon_vi = nullptr;
    if (have_exact) {
        icon_ci = new int[N * 2];
        icon_vi = new int[N * 2];
        for (int i = 0; i < N; i++) {
            int src = i % ied.n_edges;  /* wrap if ICON has fewer edges */
            icon_ci[i * 2 + 0] = ied.cell_idx[src * 2 + 0] % Nc;
            icon_ci[i * 2 + 1] = ied.cell_idx[src * 2 + 1] % Nc;
            icon_vi[i * 2 + 0] = ied.vert_idx[src * 2 + 0] % Nv;
            icon_vi[i * 2 + 1] = ied.vert_idx[src * 2 + 1] % Nv;
        }
        ied.free_all();  /* raw ICON data no longer needed */
    }

    int n_dists = have_exact ? 5 : 4;

    for(int nlev:nlevs){size_t se=(size_t)N*nlev;
        for(int di=0;di<n_dists;di++){printf("=== nlev=%d dist=%s ===\n",nlev,dist_names[di]);
            int*cl=new int[N*2],*vl=new int[N*2];
            if (di == EXACT) {
                memcpy(cl, icon_ci, N * 2 * sizeof(int));
                memcpy(vl, icon_vi, N * 2 * sizeof(int));
            } else {
                gen_conn(cl,N,Nc,(CellDist)di,rng);
                gen_conn(vl,N,Nv,(CellDist)di,rng);
            }
            Stg stg;stg.alloc(N,Nc,Nv,nlev);stg.fill();stg.set_conn(cl,vl);

            SoA_D sj,sn; sj.build(stg,false); sn.build(stg,true);
            AoS_D aj,an; aj.build(stg,false); an.build(stg,true);
            Grp_D gj,gn; gj.build(stg,false); gn.build(stg,true);
            AoSoA_D oj,on_; oj.build(stg,false); on_.build(stg,true);

            double*ref=na<double>(se);_Pragma("omp parallel for schedule(static)")for(size_t i=0;i<se;i++)ref[i]=0;
            kref(sj,ref);

            _s=&sj;_sn=&sn;_a=&aj;_an=&an;_g=&gj;_gn=&gn;_o=&oj;_on=&on_;
            const char*dn=dist_names[di];
            trun(csv,"soa",0,dn,nlev,N,ws,sj.out,ref,se,false);
            trun(csv,"soa",1,dn,nlev,N,wsc,sj.out,ref,se,false);
            trun(csv,"soa_nf",0,dn,nlev,N,wsn,sn.out,ref,se,true);
            trun(csv,"soa_nf",1,dn,nlev,N,wsnc,sn.out,ref,se,true);
            trun(csv,"aos",0,dn,nlev,N,wa,aj.out,ref,se,false);
            trun(csv,"aos",1,dn,nlev,N,wac,aj.out,ref,se,false);
            trun(csv,"aos_nf",0,dn,nlev,N,wan,an.out,ref,se,true);
            trun(csv,"aos_nf",1,dn,nlev,N,wanc,an.out,ref,se,true);
            trun(csv,"grp",0,dn,nlev,N,wg,gj.out,ref,se,false);
            trun(csv,"grp",1,dn,nlev,N,wgc,gj.out,ref,se,false);
            trun(csv,"grp_nf",0,dn,nlev,N,wgn,gn.out,ref,se,true);
            trun(csv,"grp_nf",1,dn,nlev,N,wgnc,gn.out,ref,se,true);
            trun(csv,"aosoa16",0,dn,nlev,N,wo,oj.out,ref,se,false);
            trun(csv,"aosoa16",1,dn,nlev,N,woc,oj.out,ref,se,false);
            trun(csv,"aosoa16_nf",0,dn,nlev,N,won,on_.out,ref,se,true);
            trun(csv,"aosoa16_nf",1,dn,nlev,N,wonc,on_.out,ref,se,true);

            nf_(ref,se);sj.nfree();sn.nfree();aj.nfree();an.nfree();
            gj.nfree();gn.nfree();oj.nfree();on_.nfree();
            stg.free_all();delete[]cl;delete[]vl;}}
    if (icon_ci) { delete[] icon_ci; delete[] icon_vi; }
    flush_destroy();fclose(csv);
    printf("Written: ddt_vn_apc_pc_cpu_sweep.csv\n");
    printf("  %zu nlevs × %d dists × 8 layouts × 2 collapse × %d reps = %zu rows\n",
           nlevs.size(),n_dists,NRUNS,nlevs.size()*n_dists*8*2*NRUNS);return 0;}