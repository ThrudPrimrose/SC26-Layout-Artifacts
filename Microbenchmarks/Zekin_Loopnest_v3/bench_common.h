/*
 * bench_common.h -- shared definitions for z_v_grad_w stencil benchmark
 *
 * Arrays have 3 distinct horizontal dimensions matching ICON topology:
 *   N_e  = number of edges   (direct arrays: vn_ie, z_vt_ie, out, inv_dual, etc.)
 *   N_c  = number of cells   (indirect target: w; cell_idx values in [0,N_c))
 *   N_v  = number of vertices (indirect target: z_w_v; vert_idx values in [0,N_v))
 *
 * For ICON R02B05: N_e=122880, N_c=81920, N_v=81920
 */
#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <random>

#ifndef __CUDACC__
#include <sys/mman.h>
#include <omp.h>
#endif

static constexpr int NLEVS[] = {90};
static constexpr int N_NLEVS = 1;
static constexpr int WARMUP = 5;
static constexpr int NRUNS = 100;
static constexpr int BLOCK_SIZES[] = {8, 16, 32, 64, 128};
static constexpr int N_BLOCK_SIZES = 5;
static constexpr int NUMA_DOMAINS = 4;

#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#else
#ifdef __HIP_PLATFORM_AMD__
#define HD __host__ __device__ __forceinline__
#else
#define HD inline
#endif
#endif

/* ================================================================ */
/*  2D index helpers -- unblocked                                    */
/*  N is the horizontal stride of the specific array being indexed   */
/* ================================================================ */
template <int V> HD int IC(int i, int jk, int N, int nlev) {
  if constexpr (V <= 2) return i + jk * N;
  else                  return jk + i * nlev;
}
template <int V> HD int IN(int je, int n, int N_e) {
  if constexpr (V == 1 || V == 3) return je + n * N_e;
  else                            return n + je * 2;
}

/* ================================================================ */
/*  2D index helpers -- blocked                                      */
/* ================================================================ */
HD int IC_blocked(int i, int jk, int B, int nlev) {
  return (i % B) + jk * B + (i / B) * B * nlev;
}
HD int IN_blocked(int je, int n, int B) {
  return (je % B) + n * B + (je / B) * B * 2;
}

static inline double l1_working_set_bytes(int B) { return (double)B * (8*8+4*4); }
static inline double l1_ratio(int B, int L1_bytes) { return l1_working_set_bytes(B)/(double)L1_bytes; }

static inline uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
static void fill_xor(double *arr, size_t n, unsigned seed) {
  for (size_t i = 0; i < n; i++) {
    uint64_t h = splitmix64((uint64_t)seed * 2654435761ULL + i);
    arr[i] = (double)(int64_t)(h & 0xFFFFF) / 100000.0 - 5.0;
  }
}

/* ================================================================ */
/*  distributions                                                    */
/* ================================================================ */
enum CellDist { UNIFORM = 0, NORMAL1 = 1 };
static const char *dist_name[] = {"uniform", "normal_var1", "exact"};

static void gen_idx_logical(int *L, int N_e, int idx_max, CellDist dist,
                            std::mt19937 &rng) {
  switch (dist) {
  case UNIFORM: {
    std::uniform_int_distribution<int> u(0, idx_max - 1);
    for (int i = 0; i < N_e; i++) { L[i*2]=u(rng); L[i*2+1]=u(rng); }
    break;
  }
  case NORMAL1: {
    std::normal_distribution<double> nd(0, 1);
    for (int i = 0; i < N_e; i++) {
      int base = (int)((long long)i * idx_max / N_e);
      L[i*2]   = ((base+1+(int)std::round(nd(rng))) % idx_max + idx_max) % idx_max;
      L[i*2+1] = ((base-1+(int)std::round(nd(rng))) % idx_max + idx_max) % idx_max;
    }
    break;
  }
  }
}

enum SchedKind { SCHED_JK_OUTER=0, SCHED_JE_OUTER=1, SCHED_COLLAPSE2=2,
                 SCHED_NUMA4=3, SCHED_COLLAPSE2_JE=4 };

/* ================================================================ */
/*  stencil body macros                                              */
/*  Expects N_e, N_c, N_v, nlev in scope.                           */
/* ================================================================ */
#define STENCIL_BODY(V)                                                        \
  int c2d = IC<V>(je, jk, N_e, nlev);                                          \
  int ci0 = cell_idx[IN<V>(je, 0, N_e)];                                       \
  int ci1 = cell_idx[IN<V>(je, 1, N_e)];                                       \
  int vi0 = vert_idx[IN<V>(je, 0, N_e)];                                       \
  int vi1 = vert_idx[IN<V>(je, 1, N_e)];                                       \
  out[c2d] =                                                                   \
      vn_ie[c2d] * inv_dual[je] *                                              \
          (w[IC<V>(ci0, jk, N_c, nlev)] - w[IC<V>(ci1, jk, N_c, nlev)]) +      \
      z_vt_ie[c2d] * inv_primal[je] * tangent[je] *                            \
          (z_w_v[IC<V>(vi0, jk, N_v, nlev)] - z_w_v[IC<V>(vi1, jk, N_v, nlev)]);

/* Blocked: uses B for indexing, N_e for edge loop, N_c/N_v implicit via IC_blocked */
#define STENCIL_BODY_BLOCKED(B_)                                               \
  int c2d = IC_blocked(je, jk, B_, nlev);                                      \
  int ci0 = cell_idx[IN_blocked(je, 0, B_)];                                   \
  int ci1 = cell_idx[IN_blocked(je, 1, B_)];                                   \
  int vi0 = vert_idx[IN_blocked(je, 0, B_)];                                   \
  int vi1 = vert_idx[IN_blocked(je, 1, B_)];                                   \
  out[c2d] = vn_ie[c2d] * inv_dual[je] *                                       \
                 (w[IC_blocked(ci0, jk, B_, nlev)] -                           \
                  w[IC_blocked(ci1, jk, B_, nlev)]) +                          \
             z_vt_ie[c2d] * inv_primal[je] * tangent[je] *                     \
                 (z_w_v[IC_blocked(vi0, jk, B_, nlev)] -                       \
                  z_w_v[IC_blocked(vi1, jk, B_, nlev)]);

/* ================================================================ */
/*  layout transforms -- unblocked                                   */
/* ================================================================ */
template <int V> static void layout_idx(int *dst, const int *logical, int N_e) {
  for (int je = 0; je < N_e; je++) {
    dst[IN<V>(je,0,N_e)] = logical[je*2+0];
    dst[IN<V>(je,1,N_e)] = logical[je*2+1];
  }
}
template <int V>
static void layout_2d(double *dst, const double *src, int N, int nlev) {
#pragma omp parallel for schedule(static) collapse(2)
  for (int jk = 0; jk < nlev; jk++)
    for (int i = 0; i < N; i++)
      dst[IC<V>(i,jk,N,nlev)] = src[i+jk*N];
}
static void rearrange_idx(int V, int *dst, const int *logical, int N_e) {
  switch(V){case 1:layout_idx<1>(dst,logical,N_e);break;case 2:layout_idx<2>(dst,logical,N_e);break;
            case 3:layout_idx<3>(dst,logical,N_e);break;case 4:layout_idx<4>(dst,logical,N_e);break;}
}
static void rearrange_2d(int V, double *dst, const double *src, int N, int nlev) {
  switch(V){case 1:layout_2d<1>(dst,src,N,nlev);break;case 2:layout_2d<2>(dst,src,N,nlev);break;
            case 3:layout_2d<3>(dst,src,N,nlev);break;case 4:layout_2d<4>(dst,src,N,nlev);break;}
}

static void layout_2d_blocked(double *dst, const double *src, int N, int nlev, int B) {
#pragma omp parallel for schedule(static) collapse(2)
  for (int jk = 0; jk < nlev; jk++)
    for (int i = 0; i < N; i++)
      dst[IC_blocked(i,jk,B,nlev)] = src[i+jk*N];
}
static void layout_idx_blocked(int *dst, const int *logical, int N_e, int B) {
  for (int je = 0; je < N_e; je++) {
    dst[IN_blocked(je,0,B)] = logical[je*2+0];
    dst[IN_blocked(je,1,B)] = logical[je*2+1];
  }
}

/* ================================================================ */
/*  NUMA allocation + redistribution (CPU only)                      */
/* ================================================================ */
#ifndef __CUDACC__

template <typename T> static T *numa_alloc_unfaulted(size_t count) {
  size_t bytes = count * sizeof(T);
  void *p = mmap(nullptr, bytes, PROT_READ|PROT_WRITE,
                 MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0);
  if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
  madvise(p, bytes, MADV_HUGEPAGE);
  return (T *)p;
}
template <typename T> static void numa_dealloc(T *p, size_t count) {
  if (p) munmap(p, count * sizeof(T));
}

template <int V>
static double *redistribute_2d(const double *src, int N, int nlev, SchedKind target) {
  size_t total = (size_t)N * nlev;
  double *dst = numa_alloc_unfaulted<double>(total);
  switch (target) {
  case SCHED_COLLAPSE2:
#pragma omp parallel for collapse(2) schedule(static)
    for (int jk = 0; jk < nlev; jk++)
      for (int i = 0; i < N; i++) { int x=IC<V>(i,jk,N,nlev); dst[x]=src[x]; }
    break;
  case SCHED_JE_OUTER:
#pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++)
      for (int jk = 0; jk < nlev; jk++) { int x=IC<V>(i,jk,N,nlev); dst[x]=src[x]; }
    break;
  case SCHED_JK_OUTER:
#pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++)
      for (int i = 0; i < N; i++) { int x=IC<V>(i,jk,N,nlev); dst[x]=src[x]; }
    break;
  case SCHED_NUMA4: {
    constexpr int ND = NUMA_DOMAINS;
    int tpd = std::max(1, omp_get_max_threads() / ND);
    int chunk = N / ND;
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int d = std::min(tid / tpd, ND - 1);
      int i0 = d * chunk;
      int i1 = (d == ND-1) ? N : i0 + chunk;
      int ltid = tid - d * tpd;
      int ln = (d == ND-1) ? (omp_get_num_threads() - d*tpd) : tpd;
      int my0 = i0 + (int)((long long)(i1-i0) * ltid / ln);
      int my1 = i0 + (int)((long long)(i1-i0) * (ltid+1) / ln);
      for (int i = my0; i < my1; i++)
        for (int jk = 0; jk < nlev; jk++) { int x=IC<V>(i,jk,N,nlev); dst[x]=src[x]; }
    }
    break;
  }
  case SCHED_COLLAPSE2_JE:
#pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < N; i++)
      for (int jk = 0; jk < nlev; jk++) { int x=IC<V>(i,jk,N,nlev); dst[x]=src[x]; }
    break;
  }
  return dst;
}

static double *redistribute_1d(const double *src, int N, SchedKind sched) {
  double *dst = numa_alloc_unfaulted<double>(N);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++) dst[i] = src[i];
  return dst;
}

template <int V>
static int *redistribute_idx(const int *src, int N_e, SchedKind sched) {
  int *dst = numa_alloc_unfaulted<int>(N_e * 2);
#pragma omp parallel for schedule(static)
  for (int je = 0; je < N_e; je++) {
    int i0=IN<V>(je,0,N_e), i1=IN<V>(je,1,N_e);
    dst[i0]=src[i0]; dst[i1]=src[i1];
  }
  return dst;
}

static double *redistribute_2d_blocked(const double *src, int N, int nlev,
                                       int B, SchedKind target) {
  size_t total = (size_t)N * nlev;
  double *dst = numa_alloc_unfaulted<double>(total);
  int nblocks = N / B;
#pragma omp parallel for schedule(static)
  for (int jb = 0; jb < nblocks; jb++)
    for (int jk = 0; jk < nlev; jk++)
      for (int jl = 0; jl < B; jl++) { int x=IC_blocked(jb*B+jl,jk,B,nlev); dst[x]=src[x]; }
  return dst;
}
static int *redistribute_idx_blocked(const int *src, int N_e, int B, SchedKind) {
  int *dst = numa_alloc_unfaulted<int>(N_e * 2);
  int nblocks = N_e / B;
#pragma omp parallel for schedule(static)
  for (int jb = 0; jb < nblocks; jb++)
    for (int jl = 0; jl < B; jl++) {
      int je = jb*B+jl;
      int i0=IN_blocked(je,0,B), i1=IN_blocked(je,1,B);
      dst[i0]=src[i0]; dst[i1]=src[i1];
    }
  return dst;
}

#endif /* !__CUDACC__ */

/* ================================================================ */
/*  BenchData -- separate N_e, N_c, N_v allocation                   */
/* ================================================================ */
struct BenchData {
  int N_e=0, N_c=0, N_v=0, nlev=0;
  size_t sz_e=0, sz_c=0, sz_v=0;
  int cur_V=0, cur_B=0;
  /* source (row-major, pre-transform) */
  double *src_vn_ie=nullptr, *src_w=nullptr, *src_z_vt_ie=nullptr, *src_z_w_v=nullptr;
  /* 1D edge arrays */
  double *inv_dual=nullptr, *inv_primal=nullptr, *tangent_o=nullptr;
  /* layout-transformed */
  double *h_vn_ie=nullptr, *h_w=nullptr, *h_z_vt_ie=nullptr, *h_z_w_v=nullptr;
  double *h_out=nullptr;
  int *h_cidx=nullptr, *h_vidx=nullptr;
#ifndef __CUDACC__
  bool numa_owned = false;
#endif

  void alloc(int Ne, int Nc, int Nv, int nl) {
    N_e=Ne; N_c=Nc; N_v=Nv; nlev=nl;
    sz_e=(size_t)N_e*nlev; sz_c=(size_t)N_c*nlev; sz_v=(size_t)N_v*nlev;
    src_vn_ie=new double[sz_e]; src_w=new double[sz_c];
    src_z_vt_ie=new double[sz_e]; src_z_w_v=new double[sz_v];
    inv_dual=new double[N_e]; inv_primal=new double[N_e]; tangent_o=new double[N_e];
    h_vn_ie=new double[sz_e]; h_w=new double[sz_c];
    h_z_vt_ie=new double[sz_e]; h_z_w_v=new double[sz_v];
    h_out=new double[sz_e];
    h_cidx=new int[N_e*2]; h_vidx=new int[N_e*2];
  }
  void fill(int nl) {
    fill_xor(src_vn_ie,sz_e,100+nl); fill_xor(src_w,sz_c,200+nl);
    fill_xor(src_z_vt_ie,sz_e,300+nl); fill_xor(src_z_w_v,sz_v,400+nl);
    fill_xor(inv_dual,N_e,500); fill_xor(inv_primal,N_e,600); fill_xor(tangent_o,N_e,700);
  }
  void set_variant(int V, const int *cell_logical, const int *vert_logical,
                   SchedKind sched = SCHED_JK_OUTER) {
    cur_V=V; cur_B=0;
    rearrange_2d(V, h_vn_ie, src_vn_ie, N_e, nlev);
    rearrange_2d(V, h_w, src_w, N_c, nlev);
    rearrange_2d(V, h_z_vt_ie, src_z_vt_ie, N_e, nlev);
    rearrange_2d(V, h_z_w_v, src_z_w_v, N_v, nlev);
    rearrange_idx(V, h_cidx, cell_logical, N_e);
    rearrange_idx(V, h_vidx, vert_logical, N_e);
#ifndef __CUDACC__
    _numa_redistribute(V, sched);
#else
    (void)sched;
#endif
  }
  void set_variant_blocked(int B, const int *cell_logical, const int *vert_logical,
                           SchedKind sched = SCHED_JK_OUTER) {
    cur_V=0; cur_B=B;
    layout_2d_blocked(h_vn_ie, src_vn_ie, N_e, nlev, B);
    layout_2d_blocked(h_w, src_w, N_c, nlev, B);
    layout_2d_blocked(h_z_vt_ie, src_z_vt_ie, N_e, nlev, B);
    layout_2d_blocked(h_z_w_v, src_z_w_v, N_v, nlev, B);
    layout_idx_blocked(h_cidx, cell_logical, N_e, B);
    layout_idx_blocked(h_vidx, vert_logical, N_e, B);
#ifndef __CUDACC__
    _numa_redistribute_blocked(B, sched);
#else
    (void)sched;
#endif
  }

#ifndef __CUDACC__
  void change_schedule(SchedKind sched) {
    if (cur_B > 0) _numa_redistribute_blocked(cur_B, sched);
    else           _numa_redistribute(cur_V, sched);
  }
private:
  void _free_numa() {
    if (!numa_owned) return;
    numa_dealloc(h_vn_ie,sz_e); numa_dealloc(h_w,sz_c);
    numa_dealloc(h_z_vt_ie,sz_e); numa_dealloc(h_z_w_v,sz_v);
    numa_dealloc(h_out,sz_e);
    numa_dealloc(h_cidx,(size_t)N_e*2); numa_dealloc(h_vidx,(size_t)N_e*2);
    numa_dealloc(inv_dual,(size_t)N_e); numa_dealloc(inv_primal,(size_t)N_e);
    numa_dealloc(tangent_o,(size_t)N_e);
    numa_owned = false;
  }
  void _numa_redistribute(int V, SchedKind sched) {
    double *nv,*nw,*nvt,*nzw; int *nc,*nvi;
    /* edge arrays redistriibuted at N_e, cell at N_c, vert at N_v */
    switch(V){
    case 1:
      nv=redistribute_2d<1>(h_vn_ie,N_e,nlev,sched); nw=redistribute_2d<1>(h_w,N_c,nlev,sched);
      nvt=redistribute_2d<1>(h_z_vt_ie,N_e,nlev,sched); nzw=redistribute_2d<1>(h_z_w_v,N_v,nlev,sched);
      nc=redistribute_idx<1>(h_cidx,N_e,sched); nvi=redistribute_idx<1>(h_vidx,N_e,sched); break;
    case 2:
      nv=redistribute_2d<2>(h_vn_ie,N_e,nlev,sched); nw=redistribute_2d<2>(h_w,N_c,nlev,sched);
      nvt=redistribute_2d<2>(h_z_vt_ie,N_e,nlev,sched); nzw=redistribute_2d<2>(h_z_w_v,N_v,nlev,sched);
      nc=redistribute_idx<2>(h_cidx,N_e,sched); nvi=redistribute_idx<2>(h_vidx,N_e,sched); break;
    case 3:
      nv=redistribute_2d<3>(h_vn_ie,N_e,nlev,sched); nw=redistribute_2d<3>(h_w,N_c,nlev,sched);
      nvt=redistribute_2d<3>(h_z_vt_ie,N_e,nlev,sched); nzw=redistribute_2d<3>(h_z_w_v,N_v,nlev,sched);
      nc=redistribute_idx<3>(h_cidx,N_e,sched); nvi=redistribute_idx<3>(h_vidx,N_e,sched); break;
    default:
      nv=redistribute_2d<4>(h_vn_ie,N_e,nlev,sched); nw=redistribute_2d<4>(h_w,N_c,nlev,sched);
      nvt=redistribute_2d<4>(h_z_vt_ie,N_e,nlev,sched); nzw=redistribute_2d<4>(h_z_w_v,N_v,nlev,sched);
      nc=redistribute_idx<4>(h_cidx,N_e,sched); nvi=redistribute_idx<4>(h_vidx,N_e,sched); break;
    }
    double *nout = numa_alloc_unfaulted<double>(sz_e);
#pragma omp parallel for schedule(static)
    for (size_t i=0;i<sz_e;i++) nout[i]=0.0;
    double *nid=redistribute_1d(inv_dual,N_e,sched);
    double *nip=redistribute_1d(inv_primal,N_e,sched);
    double *ntg=redistribute_1d(tangent_o,N_e,sched);
    if (numa_owned) _free_numa();
    else { delete[]h_vn_ie;delete[]h_w;delete[]h_z_vt_ie;delete[]h_z_w_v;
           delete[]h_out;delete[]h_cidx;delete[]h_vidx;
           delete[]inv_dual;delete[]inv_primal;delete[]tangent_o; }
    h_vn_ie=nv;h_w=nw;h_z_vt_ie=nvt;h_z_w_v=nzw;h_out=nout;
    h_cidx=nc;h_vidx=nvi;inv_dual=nid;inv_primal=nip;tangent_o=ntg;
    numa_owned=true;
  }
  void _numa_redistribute_blocked(int B, SchedKind sched) {
    double *nv=redistribute_2d_blocked(h_vn_ie,N_e,nlev,B,sched);
    double *nw=redistribute_2d_blocked(h_w,N_c,nlev,B,sched);
    double *nvt=redistribute_2d_blocked(h_z_vt_ie,N_e,nlev,B,sched);
    double *nzw=redistribute_2d_blocked(h_z_w_v,N_v,nlev,B,sched);
    int *nc=redistribute_idx_blocked(h_cidx,N_e,B,sched);
    int *nvi=redistribute_idx_blocked(h_vidx,N_e,B,sched);
    double *nout=numa_alloc_unfaulted<double>(sz_e);
#pragma omp parallel for schedule(static)
    for (size_t i=0;i<sz_e;i++) nout[i]=0.0;
    double *nid=redistribute_1d(inv_dual,N_e,sched);
    double *nip=redistribute_1d(inv_primal,N_e,sched);
    double *ntg=redistribute_1d(tangent_o,N_e,sched);
    if (numa_owned) _free_numa();
    else { delete[]h_vn_ie;delete[]h_w;delete[]h_z_vt_ie;delete[]h_z_w_v;
           delete[]h_out;delete[]h_cidx;delete[]h_vidx;
           delete[]inv_dual;delete[]inv_primal;delete[]tangent_o; }
    h_vn_ie=nv;h_w=nw;h_z_vt_ie=nvt;h_z_w_v=nzw;h_out=nout;
    h_cidx=nc;h_vidx=nvi;inv_dual=nid;inv_primal=nip;tangent_o=ntg;
    numa_owned=true;
  }
public:
#endif

  void free_all() {
#ifndef __CUDACC__
    if (numa_owned) { _free_numa(); } else {
#endif
      delete[]h_vn_ie;delete[]h_w;delete[]h_z_vt_ie;delete[]h_z_w_v;
      delete[]h_out;delete[]h_cidx;delete[]h_vidx;
      delete[]inv_dual;delete[]inv_primal;delete[]tangent_o;
#ifndef __CUDACC__
    }
#endif
    delete[]src_vn_ie;delete[]src_w;delete[]src_z_vt_ie;delete[]src_z_w_v;
  }
};