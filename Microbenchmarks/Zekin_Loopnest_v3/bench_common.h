/*
 * bench_common.h -- shared definitions for z_v_grad_w stencil benchmark
 *
 * 4 unblocked layout variants (V1-V4) + blocked layout (parameterized by B).
 * 4+1 cell_idx distributions: uniform, normal(var=1), normal(var=4),
 *                              sequential, exact (from ICON data)
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
#endif

/* ================================================================ */
/*  constants                                                        */
/* ================================================================ */
static constexpr int NPROMA = 163840;
static constexpr int NLEVS[] = {90};
static constexpr int N_NLEVS = 1;
static constexpr int WARMUP = 5;
static constexpr int NRUNS = 100;

static constexpr int BLOCK_SIZES[] = {8, 16, 32, 64, 128};
static constexpr int N_BLOCK_SIZES = 5;

/* ================================================================ */
/*  host/device qualifier                                            */
/* ================================================================ */
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
/*  2D index helpers -- unblocked variants                           */
/* ================================================================ */
template <int V> HD int IC(int je, int jk, int N, int nlev) {
  if constexpr (V <= 2)
    return je + jk * N;
  else
    return jk + je * nlev;
}
template <int V> HD int IN(int je, int n, int N) {
  if constexpr (V == 1 || V == 3)
    return je + n * N;
  else
    return n + je * 2;
}

/* ================================================================ */
/*  2D index helpers -- blocked layout                               */
/*  2D: (je%B) + jk*B + (je/B)*B*nlev                               */
/*  Idx: (je%B) + n*B + (je/B)*B*2                                  */
/* ================================================================ */
HD int IC_blocked(int je, int jk, int B, int nlev) {
  return (je % B) + jk * B + (je / B) * B * nlev;
}
HD int IN_blocked(int je, int n, int B) {
  return (je % B) + n * B + (je / B) * B * 2;
}

/* ================================================================ */
/*  L1 working-set ratio                                             */
/* ================================================================ */
static inline double l1_working_set_bytes(int B) {
  return (double)B * (8 * 8 + 4 * 4); /* ~80B bytes */
}
static inline double l1_ratio(int B, int L1_bytes) {
  return l1_working_set_bytes(B) / (double)L1_bytes;
}

/* ================================================================ */
/*  splitmix64 + fill                                                */
/* ================================================================ */
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
enum CellDist { UNIFORM = 0, NORMAL1 = 1, NORMAL4 = 2, SEQUENTIAL = 3 };
static const char *dist_name[] = {"uniform", "normal_var1", "normal_var4",
                                  "sequential", "exact"};

static void gen_permutation(int *p, int N, std::mt19937 &rng) {
  std::iota(p, p + N, 0);
  std::shuffle(p, p + N, rng);
}

static void gen_cell_idx_logical(int *L, int N, CellDist dist,
                                 std::mt19937 &rng) {
  switch (dist) {
  case UNIFORM: {
    std::uniform_int_distribution<int> u(0, N - 1);
    for (int i = 0; i < N; i++) {
      L[i * 2] = u(rng);
      L[i * 2 + 1] = u(rng);
    }
    break;
  }
  case NORMAL1: {
    std::normal_distribution<double> nd(0, 1);
    for (int i = 0; i < N; i++) {
      L[i * 2] = ((i + 1 + (int)std::round(nd(rng))) % N + N) % N;
      L[i * 2 + 1] = ((i - 1 + (int)std::round(nd(rng))) % N + N) % N;
    }
    break;
  }
  case NORMAL4: {
    std::normal_distribution<double> nd(0, 2);
    for (int i = 0; i < N; i++) {
      L[i * 2] = ((i + 1 + (int)std::round(nd(rng))) % N + N) % N;
      L[i * 2 + 1] = ((i - 1 + (int)std::round(nd(rng))) % N + N) % N;
    }
    break;
  }
  case SEQUENTIAL:
    for (int i = 0; i < N; i++) {
      L[i * 2] = (i + 1) % N;
      L[i * 2 + 1] = (i + 1) % N;
    }
    break;
  }
}

/* ================================================================ */
/*  NUMA schedule kinds                                              */
/* ================================================================ */
enum SchedKind { SCHED_JK_OUTER = 0, SCHED_JE_OUTER = 1, SCHED_COLLAPSE2 = 2 };

/* ================================================================ */
/*  stencil body macros                                              */
/* ================================================================ */
#define STENCIL_BODY(V)                                                        \
  int c2d = IC<V>(je, jk, N, nlev);                                            \
  int ci0 = cell_idx[IN<V>(je, 0, N)];                                         \
  int ci1 = cell_idx[IN<V>(je, 1, N)];                                         \
  int vi0 = vert_idx[IN<V>(je, 0, N)];                                         \
  int vi1 = vert_idx[IN<V>(je, 1, N)];                                         \
  out[c2d] =                                                                   \
      vn_ie[c2d] * inv_dual[je] *                                              \
          (w[IC<V>(ci0, jk, N, nlev)] - w[IC<V>(ci1, jk, N, nlev)]) +          \
      z_vt_ie[c2d] * inv_primal[je] * tangent[je] *                            \
          (z_w_v[IC<V>(vi0, jk, N, nlev)] - z_w_v[IC<V>(vi1, jk, N, nlev)]);

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
template <int V> static void layout_idx(int *dst, const int *logical, int N) {
  for (int je = 0; je < N; je++) {
    dst[IN<V>(je, 0, N)] = logical[je * 2 + 0];
    dst[IN<V>(je, 1, N)] = logical[je * 2 + 1];
  }
}
template <int V>
static void layout_2d(double *dst, const double *src, int N, int nlev) {
#pragma omp parallel for schedule(static) collapse(2)
  for (int jk = 0; jk < nlev; jk++)
    for (int je = 0; je < N; je++)
      dst[IC<V>(je, jk, N, nlev)] = src[je + jk * N];
}
static void rearrange_idx(int V, int *dst, const int *logical, int N) {
  switch (V) {
  case 1:
    layout_idx<1>(dst, logical, N);
    break;
  case 2:
    layout_idx<2>(dst, logical, N);
    break;
  case 3:
    layout_idx<3>(dst, logical, N);
    break;
  case 4:
    layout_idx<4>(dst, logical, N);
    break;
  }
}
static void rearrange_2d(int V, double *dst, const double *src, int N,
                         int nlev) {
  switch (V) {
  case 1:
    layout_2d<1>(dst, src, N, nlev);
    break;
  case 2:
    layout_2d<2>(dst, src, N, nlev);
    break;
  case 3:
    layout_2d<3>(dst, src, N, nlev);
    break;
  case 4:
    layout_2d<4>(dst, src, N, nlev);
    break;
  }
}

/* ================================================================ */
/*  layout transforms -- blocked                                     */
/* ================================================================ */
static void layout_2d_blocked(double *dst, const double *src, int N, int nlev,
                              int B) {
#pragma omp parallel for schedule(static) collapse(2)
  for (int jk = 0; jk < nlev; jk++)
    for (int je = 0; je < N; je++)
      dst[IC_blocked(je, jk, B, nlev)] = src[je + jk * N];
}
static void layout_idx_blocked(int *dst, const int *logical, int N, int B) {
  for (int je = 0; je < N; je++) {
    dst[IN_blocked(je, 0, B)] = logical[je * 2 + 0];
    dst[IN_blocked(je, 1, B)] = logical[je * 2 + 1];
  }
}

/* ================================================================ */
/*  NUMA allocation + redistribution (CPU only)                      */
/* ================================================================ */
#ifndef __CUDACC__

template <typename T> static T *numa_alloc_unfaulted(size_t count) {
  size_t bytes = count * sizeof(T);
  void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                 MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
  if (p == MAP_FAILED) {
    perror("mmap");
    std::abort();
  }
  madvise(p, bytes, MADV_HUGEPAGE);
  return (T *)p;
}
template <typename T> static void numa_dealloc(T *p, size_t count) {
  if (p)
    munmap(p, count * sizeof(T));
}

template <int V>
static double *redistribute_2d(const double *src, int N, int nlev,
                               SchedKind target) {
  size_t total = (size_t)N * nlev;
  double *dst = numa_alloc_unfaulted<double>(total);
  switch (target) {
  case SCHED_COLLAPSE2:
#pragma omp parallel for collapse(2) schedule(static)
    for (int jk = 0; jk < nlev; jk++)
      for (int je = 0; je < N; je++) {
        int i = IC<V>(je, jk, N, nlev);
        dst[i] = src[i];
      }
    break;
  case SCHED_JE_OUTER:
#pragma omp parallel for schedule(static)
    for (int je = 0; je < N; je++)
      for (int jk = 0; jk < nlev; jk++) {
        int i = IC<V>(je, jk, N, nlev);
        dst[i] = src[i];
      }
    break;
  case SCHED_JK_OUTER:
#pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++)
      for (int je = 0; je < N; je++) {
        int i = IC<V>(je, jk, N, nlev);
        dst[i] = src[i];
      }
    break;
  }
  return dst;
}
static double *redistribute_1d(const double *src, int N, SchedKind) {
  double *dst = numa_alloc_unfaulted<double>(N);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++)
    dst[i] = src[i];
  return dst;
}
template <int V>
static int *redistribute_idx(const int *src, int N, SchedKind) {
  int *dst = numa_alloc_unfaulted<int>(N * 2);
#pragma omp parallel for schedule(static)
  for (int je = 0; je < N; je++) {
    int i0 = IN<V>(je, 0, N), i1 = IN<V>(je, 1, N);
    dst[i0] = src[i0];
    dst[i1] = src[i1];
  }
  return dst;
}
static double *redistribute_2d_blocked(const double *src, int N, int nlev,
                                       int B, SchedKind target) {
  size_t total = (size_t)N * nlev;
  double *dst = numa_alloc_unfaulted<double>(total);
  int nblocks = N / B;
  switch (target) {
  case SCHED_COLLAPSE2:
#pragma omp parallel for collapse(2) schedule(static)
    for (int jb = 0; jb < nblocks; jb++)
      for (int jk = 0; jk < nlev; jk++)
        for (int jl = 0; jl < B; jl++) {
          int i = IC_blocked(jb * B + jl, jk, B, nlev);
          dst[i] = src[i];
        }
    break;
  default:
#pragma omp parallel for schedule(static)
    for (int jb = 0; jb < nblocks; jb++)
      for (int jk = 0; jk < nlev; jk++)
        for (int jl = 0; jl < B; jl++) {
          int i = IC_blocked(jb * B + jl, jk, B, nlev);
          dst[i] = src[i];
        }
    break;
  }
  return dst;
}
static int *redistribute_idx_blocked(const int *src, int N, int B, SchedKind) {
  int *dst = numa_alloc_unfaulted<int>(N * 2);
  int nblocks = N / B;
#pragma omp parallel for schedule(static)
  for (int jb = 0; jb < nblocks; jb++)
    for (int jl = 0; jl < B; jl++) {
      int je = jb * B + jl;
      int i0 = IN_blocked(je, 0, B), i1 = IN_blocked(je, 1, B);
      dst[i0] = src[i0];
      dst[i1] = src[i1];
    }
  return dst;
}

#endif /* !__CUDACC__ */

/* ================================================================ */
/*  BenchData                                                        */
/* ================================================================ */
struct BenchData {
  int N = 0, nlev = 0;
  size_t sz2d = 0;
  int cur_V = 0, cur_B = 0;
  double *src_vn_ie = nullptr, *src_w = nullptr, *src_z_vt_ie = nullptr,
         *src_z_w_v = nullptr;
  double *inv_dual = nullptr, *inv_primal = nullptr, *tangent_o = nullptr;
  double *h_vn_ie = nullptr, *h_w = nullptr, *h_z_vt_ie = nullptr,
         *h_z_w_v = nullptr, *h_out = nullptr;
  int *h_cidx = nullptr, *h_vidx = nullptr;
#ifndef __CUDACC__
  bool numa_owned = false;
#endif

  void alloc(int N_, int nlev_) {
    N = N_;
    nlev = nlev_;
    sz2d = (size_t)N * nlev;
    src_vn_ie = new double[sz2d];
    src_w = new double[sz2d];
    src_z_vt_ie = new double[sz2d];
    src_z_w_v = new double[sz2d];
    inv_dual = new double[N];
    inv_primal = new double[N];
    tangent_o = new double[N];
    h_vn_ie = new double[sz2d];
    h_w = new double[sz2d];
    h_z_vt_ie = new double[sz2d];
    h_z_w_v = new double[sz2d];
    h_out = new double[sz2d];
    h_cidx = new int[N * 2];
    h_vidx = new int[N * 2];
  }
  void fill(int nlev_) {
    fill_xor(src_vn_ie, sz2d, 100 + nlev_);
    fill_xor(src_w, sz2d, 200 + nlev_);
    fill_xor(src_z_vt_ie, sz2d, 300 + nlev_);
    fill_xor(src_z_w_v, sz2d, 400 + nlev_);
    fill_xor(inv_dual, N, 500);
    fill_xor(inv_primal, N, 600);
    fill_xor(tangent_o, N, 700);
  }
  void set_variant(int V, const int *cell_logical, const int *vert_logical,
                   SchedKind sched = SCHED_JK_OUTER) {
    cur_V = V;
    cur_B = 0;
    rearrange_2d(V, h_vn_ie, src_vn_ie, N, nlev);
    rearrange_2d(V, h_w, src_w, N, nlev);
    rearrange_2d(V, h_z_vt_ie, src_z_vt_ie, N, nlev);
    rearrange_2d(V, h_z_w_v, src_z_w_v, N, nlev);
    rearrange_idx(V, h_cidx, cell_logical, N);
    rearrange_idx(V, h_vidx, vert_logical, N);
#ifndef __CUDACC__
    _numa_redistribute(V, sched);
#else
    (void)sched;
#endif
  }
  void set_variant_blocked(int B, const int *cell_logical,
                           const int *vert_logical,
                           SchedKind sched = SCHED_JK_OUTER) {
    cur_V = 0;
    cur_B = B;
    layout_2d_blocked(h_vn_ie, src_vn_ie, N, nlev, B);
    layout_2d_blocked(h_w, src_w, N, nlev, B);
    layout_2d_blocked(h_z_vt_ie, src_z_vt_ie, N, nlev, B);
    layout_2d_blocked(h_z_w_v, src_z_w_v, N, nlev, B);
    layout_idx_blocked(h_cidx, cell_logical, N, B);
    layout_idx_blocked(h_vidx, vert_logical, N, B);
#ifndef __CUDACC__
    _numa_redistribute_blocked(B, sched);
#else
    (void)sched;
#endif
  }

#ifndef __CUDACC__
  void change_schedule(SchedKind sched) {
    if (cur_B > 0)
      _numa_redistribute_blocked(cur_B, sched);
    else
      _numa_redistribute(cur_V, sched);
  }

private:
  void _free_numa() {
    if (!numa_owned)
      return;
    numa_dealloc(h_vn_ie, sz2d);
    numa_dealloc(h_w, sz2d);
    numa_dealloc(h_z_vt_ie, sz2d);
    numa_dealloc(h_z_w_v, sz2d);
    numa_dealloc(h_out, sz2d);
    numa_dealloc(h_cidx, (size_t)N * 2);
    numa_dealloc(h_vidx, (size_t)N * 2);
    numa_dealloc(inv_dual, (size_t)N);
    numa_dealloc(inv_primal, (size_t)N);
    numa_dealloc(tangent_o, (size_t)N);
    numa_owned = false;
  }
  void _numa_redistribute(int V, SchedKind sched) {
    double *nv, *nw, *nvt, *nzw;
    int *nc, *nvi;
    switch (V) {
    case 1:
      nv = redistribute_2d<1>(h_vn_ie, N, nlev, sched);
      nw = redistribute_2d<1>(h_w, N, nlev, sched);
      nvt = redistribute_2d<1>(h_z_vt_ie, N, nlev, sched);
      nzw = redistribute_2d<1>(h_z_w_v, N, nlev, sched);
      nc = redistribute_idx<1>(h_cidx, N, sched);
      nvi = redistribute_idx<1>(h_vidx, N, sched);
      break;
    case 2:
      nv = redistribute_2d<2>(h_vn_ie, N, nlev, sched);
      nw = redistribute_2d<2>(h_w, N, nlev, sched);
      nvt = redistribute_2d<2>(h_z_vt_ie, N, nlev, sched);
      nzw = redistribute_2d<2>(h_z_w_v, N, nlev, sched);
      nc = redistribute_idx<2>(h_cidx, N, sched);
      nvi = redistribute_idx<2>(h_vidx, N, sched);
      break;
    case 3:
      nv = redistribute_2d<3>(h_vn_ie, N, nlev, sched);
      nw = redistribute_2d<3>(h_w, N, nlev, sched);
      nvt = redistribute_2d<3>(h_z_vt_ie, N, nlev, sched);
      nzw = redistribute_2d<3>(h_z_w_v, N, nlev, sched);
      nc = redistribute_idx<3>(h_cidx, N, sched);
      nvi = redistribute_idx<3>(h_vidx, N, sched);
      break;
    default:
      nv = redistribute_2d<4>(h_vn_ie, N, nlev, sched);
      nw = redistribute_2d<4>(h_w, N, nlev, sched);
      nvt = redistribute_2d<4>(h_z_vt_ie, N, nlev, sched);
      nzw = redistribute_2d<4>(h_z_w_v, N, nlev, sched);
      nc = redistribute_idx<4>(h_cidx, N, sched);
      nvi = redistribute_idx<4>(h_vidx, N, sched);
      break;
    }
    double *nout = numa_alloc_unfaulted<double>(sz2d);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < sz2d; i++)
      nout[i] = 0.0;
    double *nid = redistribute_1d(inv_dual, N, sched);
    double *nip = redistribute_1d(inv_primal, N, sched);
    double *ntg = redistribute_1d(tangent_o, N, sched);
    if (numa_owned)
      _free_numa();
    else {
      delete[] h_vn_ie;
      delete[] h_w;
      delete[] h_z_vt_ie;
      delete[] h_z_w_v;
      delete[] h_out;
      delete[] h_cidx;
      delete[] h_vidx;
      delete[] inv_dual;
      delete[] inv_primal;
      delete[] tangent_o;
    }
    h_vn_ie = nv;
    h_w = nw;
    h_z_vt_ie = nvt;
    h_z_w_v = nzw;
    h_out = nout;
    h_cidx = nc;
    h_vidx = nvi;
    inv_dual = nid;
    inv_primal = nip;
    tangent_o = ntg;
    numa_owned = true;
  }
  void _numa_redistribute_blocked(int B, SchedKind sched) {
    double *nv = redistribute_2d_blocked(h_vn_ie, N, nlev, B, sched);
    double *nw = redistribute_2d_blocked(h_w, N, nlev, B, sched);
    double *nvt = redistribute_2d_blocked(h_z_vt_ie, N, nlev, B, sched);
    double *nzw = redistribute_2d_blocked(h_z_w_v, N, nlev, B, sched);
    int *nc = redistribute_idx_blocked(h_cidx, N, B, sched);
    int *nvi = redistribute_idx_blocked(h_vidx, N, B, sched);
    double *nout = numa_alloc_unfaulted<double>(sz2d);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < sz2d; i++)
      nout[i] = 0.0;
    double *nid = redistribute_1d(inv_dual, N, sched);
    double *nip = redistribute_1d(inv_primal, N, sched);
    double *ntg = redistribute_1d(tangent_o, N, sched);
    if (numa_owned)
      _free_numa();
    else {
      delete[] h_vn_ie;
      delete[] h_w;
      delete[] h_z_vt_ie;
      delete[] h_z_w_v;
      delete[] h_out;
      delete[] h_cidx;
      delete[] h_vidx;
      delete[] inv_dual;
      delete[] inv_primal;
      delete[] tangent_o;
    }
    h_vn_ie = nv;
    h_w = nw;
    h_z_vt_ie = nvt;
    h_z_w_v = nzw;
    h_out = nout;
    h_cidx = nc;
    h_vidx = nvi;
    inv_dual = nid;
    inv_primal = nip;
    tangent_o = ntg;
    numa_owned = true;
  }

public:
#endif

  void free_all() {
#ifndef __CUDACC__
    if (numa_owned) {
      _free_numa();
    } else {
#endif
      delete[] h_vn_ie;
      delete[] h_w;
      delete[] h_z_vt_ie;
      delete[] h_z_w_v;
      delete[] h_out;
      delete[] h_cidx;
      delete[] h_vidx;
      delete[] inv_dual;
      delete[] inv_primal;
      delete[] tangent_o;
#ifndef __CUDACC__
    }
#endif
    delete[] src_vn_ie;
    delete[] src_w;
    delete[] src_z_vt_ie;
    delete[] src_z_w_v;
  }
};

struct VertData {
  int *logical = nullptr;
  void init(int N, std::mt19937 &rng) {
    logical = new int[N * 2];
    int *tmp = new int[N];
    gen_permutation(tmp, N, rng);
    for (int i = 0; i < N; i++)
      logical[i * 2] = tmp[i];
    gen_permutation(tmp, N, rng);
    for (int i = 0; i < N; i++)
      logical[i * 2 + 1] = tmp[i];
    delete[] tmp;
  }
  void free_all() { delete[] logical; }
};