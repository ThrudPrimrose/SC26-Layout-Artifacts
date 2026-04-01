/*
 * cost_metrics.cpp -- Cost metrics for z_v_grad_w stencil
 *
 * 4 metrics: mu, mu*Delta, mu*Delta_numa, Delta_max
 * Supports unblocked (V1-V4) + blocked (B=8..128) + exact ICON dist
 *
 * Compile: g++ -O3 -std=c++17 -fopenmp -o cost_metrics cost_metrics.cpp
 * Run:     ./cost_metrics [N] [nlev] [beta] [alpha] [gamma] [P_NUMA] [L1_bytes]
 */
#include "icon_data_loader.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <random>
#include <unistd.h>
#include <vector>

static int BETA = 1;
static double ALPHA = 0.012, GAMMA = 1.8;
static int P_NUMA = 4, BLOCK_BYTES_G = 64, L1_BYTES = 32768;
#pragma omp threadprivate(BLOCK_BYTES_G)

/* ---- Index functions ---- */
inline int IC_v(int V, int je, int jk, int N, int nl) {
  return (V <= 2) ? je + jk * N : jk + je * nl;
}
inline int IN_v(int V, int je, int n, int N) {
  return (V == 1 || V == 3) ? je + n * N : n + je * 2;
}
inline int IC_b(int je, int jk, int B, int nl) {
  return (je % B) + jk * B + (je / B) * B * nl;
}
inline int IN_b(int je, int n, int B) {
  return (je % B) + n * B + (je / B) * B * 2;
}

enum ArrID {
  A_OUT = 0,
  A_VN_IE,
  A_W,
  A_Z_VT_IE,
  A_Z_W_V,
  A_INV_DUAL,
  A_INV_PRIMAL,
  A_TANGENT,
  A_CELL_IDX,
  A_VERT_IDX,
  NUM_ARR
};
static const int eb[] = {8, 8, 8, 8, 8, 8, 8, 8, 4, 4};
inline int64_t blk(int ei, int a) {
  return (int64_t)ei * eb[a] / BLOCK_BYTES_G;
}

enum Schedule { SCHED_OMP_FOR = 0, SCHED_OMP_COLLAPSE2 = 1 };
static const char *sched_name[] = {"omp_for", "omp_collapse2"};
enum LoopOrder { KLON_FIRST = 0, KLEV_FIRST = 1 };
inline LoopOrder iter_order(Schedule s, int V) {
  if (s == SCHED_OMP_COLLAPSE2)
    return KLON_FIRST;
  return (V <= 2) ? KLON_FIRST : KLEV_FIRST;
}

/* ---- NUMA ---- */
static Schedule G_SCHED;
#pragma omp threadprivate(G_SCHED)
static int NU_HOME = 0;
#pragma omp threadprivate(NU_HOME)

inline int numa_2d(Schedule s, int V, int e, int N, int nl) {
  int64_t tot = (int64_t)N * nl;
  if (s == SCHED_OMP_COLLAPSE2) {
    int64_t lin;
    if (V <= 2)
      lin = e;
    else {
      int jk = e % nl, je = e / nl;
      lin = (int64_t)jk * N + je;
    }
    int64_t ch = (tot + P_NUMA - 1) / P_NUMA;
    int d = (int)(lin / ch);
    return (d < P_NUMA) ? d : P_NUMA - 1;
  }
  if (V <= 2) {
    int jk = e / N;
    int ch = (nl + P_NUMA - 1) / P_NUMA;
    int d = jk / ch;
    return (d < P_NUMA) ? d : P_NUMA - 1;
  } else {
    int je = e / nl;
    int ch = (N + P_NUMA - 1) / P_NUMA;
    int d = je / ch;
    return (d < P_NUMA) ? d : P_NUMA - 1;
  }
}
inline int numa_2d_b(int e, int B, int N, int nl) {
  int jb = e / (B * nl);
  int nb = N / B;
  int ch = (nb + P_NUMA - 1) / P_NUMA;
  int d = jb / ch;
  return (d < P_NUMA) ? d : P_NUMA - 1;
}
inline int numa_1d(int e, int N) {
  int ch = (N + P_NUMA - 1) / P_NUMA;
  if (!ch)
    return 0;
  int d = e / ch;
  return (d < P_NUMA) ? d : P_NUMA - 1;
}
inline int nu_blk_v(int a, int V, int64_t ba, int N, int nl) {
  int e = (int)(ba * BLOCK_BYTES_G / eb[a]);
  switch (a) {
  case A_OUT:
  case A_VN_IE:
  case A_W:
  case A_Z_VT_IE:
  case A_Z_W_V:
    return numa_2d(G_SCHED, V, e, N, nl);
  case A_INV_DUAL:
  case A_INV_PRIMAL:
  case A_TANGENT:
    return numa_1d(e, N);
  case A_CELL_IDX:
  case A_VERT_IDX: {
    int je = (V == 1 || V == 3) ? e % N : e / 2;
    return numa_1d(je, N);
  }
  default:
    return 0;
  }
}
inline int nu_blk_b(int a, int B, int64_t ba, int N, int nl) {
  int e = (int)(ba * BLOCK_BYTES_G / eb[a]);
  switch (a) {
  case A_OUT:
  case A_VN_IE:
  case A_W:
  case A_Z_VT_IE:
  case A_Z_W_V:
    return numa_2d_b(e, B, N, nl);
  case A_INV_DUAL:
  case A_INV_PRIMAL:
  case A_TANGENT:
    return numa_1d(e, N);
  case A_CELL_IDX:
  case A_VERT_IDX: {
    int jb = e / (B * 2);
    int nb = N / B;
    int ch = (nb + P_NUMA - 1) / P_NUMA;
    int d = jb / ch;
    return (d < P_NUMA) ? d : P_NUMA - 1;
  }
  default:
    return 0;
  }
}
inline double wn_v(int a, int V, int64_t ba, int64_t rd, int N, int nl) {
  if (nu_blk_v(a, V, ba, N, nl) != NU_HOME)
    return GAMMA;
  return (rd < BETA) ? ALPHA : 1.0;
}
inline double wn_b(int a, int B, int64_t ba, int64_t rd, int N, int nl) {
  if (nu_blk_b(a, B, ba, N, nl) != NU_HOME)
    return GAMMA;
  return (rd < BETA) ? ALPHA : 1.0;
}

/* ---- Refs ---- */
static constexpr int NR = 14;
struct Ref {
  int a;
  int64_t b;
};
inline void refs_v(int V, int je, int jk, int N, int nl, const int *ci,
                   const int *vi, Ref r[NR]) {
  int c = IC_v(V, je, jk, N, nl);
  int c0 = ci[IN_v(V, je, 0, N)], c1 = ci[IN_v(V, je, 1, N)],
      v0 = vi[IN_v(V, je, 0, N)], v1 = vi[IN_v(V, je, 1, N)];
  r[0] = {A_OUT, blk(c, A_OUT)};
  r[1] = {A_VN_IE, blk(c, A_VN_IE)};
  r[2] = {A_INV_DUAL, blk(je, A_INV_DUAL)};
  r[3] = {A_W, blk(IC_v(V, c0, jk, N, nl), A_W)};
  r[4] = {A_W, blk(IC_v(V, c1, jk, N, nl), A_W)};
  r[5] = {A_Z_VT_IE, blk(c, A_Z_VT_IE)};
  r[6] = {A_INV_PRIMAL, blk(je, A_INV_PRIMAL)};
  r[7] = {A_TANGENT, blk(je, A_TANGENT)};
  r[8] = {A_Z_W_V, blk(IC_v(V, v0, jk, N, nl), A_Z_W_V)};
  r[9] = {A_Z_W_V, blk(IC_v(V, v1, jk, N, nl), A_Z_W_V)};
  r[10] = {A_CELL_IDX, blk(IN_v(V, je, 0, N), A_CELL_IDX)};
  r[11] = {A_CELL_IDX, blk(IN_v(V, je, 1, N), A_CELL_IDX)};
  r[12] = {A_VERT_IDX, blk(IN_v(V, je, 0, N), A_VERT_IDX)};
  r[13] = {A_VERT_IDX, blk(IN_v(V, je, 1, N), A_VERT_IDX)};
}
inline void refs_b(int B, int je, int jk, int N, int nl, const int *ci,
                   const int *vi, Ref r[NR]) {
  int c = IC_b(je, jk, B, nl);
  int c0 = ci[IN_b(je, 0, B)], c1 = ci[IN_b(je, 1, B)], v0 = vi[IN_b(je, 0, B)],
      v1 = vi[IN_b(je, 1, B)];
  r[0] = {A_OUT, blk(c, A_OUT)};
  r[1] = {A_VN_IE, blk(c, A_VN_IE)};
  r[2] = {A_INV_DUAL, blk(je, A_INV_DUAL)};
  r[3] = {A_W, blk(IC_b(c0, jk, B, nl), A_W)};
  r[4] = {A_W, blk(IC_b(c1, jk, B, nl), A_W)};
  r[5] = {A_Z_VT_IE, blk(c, A_Z_VT_IE)};
  r[6] = {A_INV_PRIMAL, blk(je, A_INV_PRIMAL)};
  r[7] = {A_TANGENT, blk(je, A_TANGENT)};
  r[8] = {A_Z_W_V, blk(IC_b(v0, jk, B, nl), A_Z_W_V)};
  r[9] = {A_Z_W_V, blk(IC_b(v1, jk, B, nl), A_Z_W_V)};
  r[10] = {A_CELL_IDX, blk(IN_b(je, 0, B), A_CELL_IDX)};
  r[11] = {A_CELL_IDX, blk(IN_b(je, 1, B), A_CELL_IDX)};
  r[12] = {A_VERT_IDX, blk(IN_b(je, 0, B), A_VERT_IDX)};
  r[13] = {A_VERT_IDX, blk(IN_b(je, 1, B), A_VERT_IDX)};
}

struct BS {
  std::vector<int64_t> a[NUM_ARR];
  void clear() {
    for (int i = 0; i < NUM_ARR; i++)
      a[i].clear();
  }
  void add(int ar, int64_t b) { a[ar].push_back(b); }
  void fin() {
    for (int i = 0; i < NUM_ARR; i++) {
      auto &v = a[i];
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
    }
  }
  int tot() const {
    int n = 0;
    for (int i = 0; i < NUM_ARR; i++)
      n += (int)a[i].size();
    return n;
  }
  bool has(int ar, int64_t b) const {
    return std::binary_search(a[ar].begin(), a[ar].end(), b);
  }
};

struct Metrics {
  double mu, mu_d, mu_dn, dmax;
  int64_t T;
};

static inline void step(int W, int N, int nl, LoopOrder lo, int je0, int jk0,
                        const int *ci, const int *vi, BS &prev, BS &curr,
                        int64_t &T, double &smu, double &smd, double &smdn,
                        double &sdm, int V, int B) {
  curr.clear();
  for (int w = 0; w < W; w++) {
    int je = (lo == KLON_FIRST) ? je0 + w : je0;
    int jk = (lo == KLON_FIRST) ? jk0 : jk0 + w;
    Ref r[NR];
    if (V > 0)
      refs_v(V, je, jk, N, nl, ci, vi, r);
    else
      refs_b(B, je, jk, N, nl, ci, vi, r);
    for (int i = 0; i < NR; i++)
      curr.add(r[i].a, r[i].b);
  }
  curr.fin();
  if (T == 0) {
    int nb = curr.tot();
    smu += nb;
    smd += nb;
    smdn += nb;
    sdm += 1.0;
  } else {
    int nc = 0;
    double sm = 0, sx = 0, sn = 0;
    for (int a = 0; a < NUM_ARR; a++)
      for (int64_t b : curr.a[a]) {
        if (prev.has(a, b))
          continue;
        nc++;
        int64_t bm = INT64_MAX, bx = 0;
        for (int64_t bp : prev.a[a]) {
          int64_t d = std::abs(b - bp);
          if (d < bm)
            bm = d;
          if (d > bx)
            bx = d;
        }
        if (bm == INT64_MAX) {
          sm += 1;
          sx += 1;
          double wn = (V > 0) ? wn_v(a, V, b, INT64_MAX, N, nl)
                              : wn_b(a, B, b, INT64_MAX, N, nl);
          sn += wn;
        } else {
          double wn =
              (V > 0) ? wn_v(a, V, b, bm, N, nl) : wn_b(a, B, b, bm, N, nl);
          sm += (double)bm;
          sx += (double)bx;
          sn += wn * (double)bm;
        }
      }
    smu += nc;
    smd += sm;
    smdn += sn;
    if (nc > 0)
      sdm += sx / nc;
  }
  T++;
  std::swap(prev, curr);
}

Metrics compute_slice(int V, int B, int W, int N, int nl, Schedule sc,
                      const int *ci, const int *vi, int64_t lo, int64_t hi) {
  G_SCHED = sc;
  LoopOrder lp = (V > 0) ? iter_order(sc, V) : KLON_FIRST;
  BS prev, curr;
  int64_t T = 0;
  double smu = 0, smd = 0, smdn = 0, sdm = 0;
  if (B > 0) {
    for (int64_t jb = lo; jb < hi; jb++)
      for (int jk = 0; jk < nl; jk++)
        step(B, N, nl, KLON_FIRST, (int)jb * B, jk, ci, vi, prev, curr, T, smu,
             smd, smdn, sdm, 0, B);
  } else if (sc == SCHED_OMP_COLLAPSE2) {
    int64_t ln = lo;
    while (ln < hi) {
      int jk = (int)(ln / N);
      int js = (int)(ln % N);
      int64_t re = std::min((int64_t)(jk + 1) * N, hi);
      int je_end = (int)(re - (int64_t)jk * N);
      for (int j = js; j + W <= je_end; j += W)
        step(W, N, nl, KLON_FIRST, j, jk, ci, vi, prev, curr, T, smu, smd, smdn,
             sdm, V, 0);
      ln = re;
    }
  } else {
    int inn = (lp == KLON_FIRST) ? N : nl;
    for (int64_t o = lo; o < hi; o++)
      for (int i = 0; i + W <= inn; i += W) {
        int je = (lp == KLON_FIRST) ? i : (int)o;
        int jk = (lp == KLON_FIRST) ? (int)o : i;
        step(W, N, nl, lp, je, jk, ci, vi, prev, curr, T, smu, smd, smdn, sdm,
             V, 0);
      }
  }
  if (!T)
    return {};
  double d = (double)T;
  return {smu / d, smd / d, smdn / d, sdm / d, T};
}

Metrics compute_metrics(int V, int B, int W, int N, int nl, Schedule sc,
                        const int *ci, const int *vi) {
  LoopOrder lp = (V > 0) ? iter_order(sc, V) : KLON_FIRST;
  int on;
  if (B > 0)
    on = N / B;
  else
    on = (lp == KLON_FIRST) ? nl : N;
  Metrics ac = {};
  int nd = 0;
  auto add = [&](Metrics &m) {
    if (!m.T)
      return;
    ac.mu += m.mu;
    ac.mu_d += m.mu_d;
    ac.mu_dn += m.mu_dn;
    ac.dmax += m.dmax;
    ac.T += m.T;
    nd++;
  };
  if (B == 0 && sc == SCHED_OMP_COLLAPSE2) {
    int64_t tot = (int64_t)on * ((lp == KLON_FIRST) ? N : nl);
    int64_t ch = (tot + P_NUMA - 1) / P_NUMA;
    for (int d = 0; d < P_NUMA; d++) {
      int64_t l = d * ch, h = std::min((d + 1) * ch, tot);
      if (l >= h)
        continue;
      NU_HOME = d;
      auto m = compute_slice(V, B, W, N, nl, sc, ci, vi, l, h);
      add(m);
    }
  } else {
    int ch = (on + P_NUMA - 1) / P_NUMA;
    for (int d = 0; d < P_NUMA; d++) {
      int l = d * ch, h = std::min((d + 1) * ch, on);
      if (l >= h)
        continue;
      NU_HOME = d;
      auto m = compute_slice(V, B, W, N, nl, sc, ci, vi, l, h);
      add(m);
    }
  }
  if (!nd)
    return {};
  double dn = (double)nd;
  return {ac.mu / dn, ac.mu_d / dn, ac.mu_dn / dn, ac.dmax / dn, ac.T};
}

/* ---- Distributions ---- */
enum CellDist { UNIFORM = 0, NORMAL1 = 1, NORMAL4 = 2, SEQUENTIAL = 3 };
static const char *dname[] = {"uniform", "normal_var1", "normal_var4",
                              "sequential", "exact"};
static void gen_ci(int *d, int V, int N, CellDist dist, std::mt19937 &rng) {
  std::vector<int> L(N * 2);
  switch (dist) {
  case UNIFORM: {
    std::uniform_int_distribution<int> u(0, N - 1);
    for (int i = 0; i < N; i++) {
      L[2 * i] = u(rng);
      L[2 * i + 1] = u(rng);
    }
    break;
  }
  case NORMAL1: {
    std::normal_distribution<double> nd(0, 1);
    for (int i = 0; i < N; i++) {
      L[2 * i] = ((i + 1 + (int)std::round(nd(rng))) % N + N) % N;
      L[2 * i + 1] = ((i - 1 + (int)std::round(nd(rng))) % N + N) % N;
    }
    break;
  }
  case NORMAL4: {
    std::normal_distribution<double> nd(0, 2);
    for (int i = 0; i < N; i++) {
      L[2 * i] = ((i + 1 + (int)std::round(nd(rng))) % N + N) % N;
      L[2 * i + 1] = ((i - 1 + (int)std::round(nd(rng))) % N + N) % N;
    }
    break;
  }
  case SEQUENTIAL:
    for (int i = 0; i < N; i++) {
      L[2 * i] = (i + 1) % N;
      L[2 * i + 1] = (i + 1) % N;
    }
    break;
  }
  for (int je = 0; je < N; je++) {
    d[IN_v(V, je, 0, N)] = L[2 * je];
    d[IN_v(V, je, 1, N)] = L[2 * je + 1];
  }
}
static void gen_ci_b(int *d, int N, int B, CellDist dist, std::mt19937 &rng) {
  std::vector<int> L(N * 2);
  switch (dist) {
  case UNIFORM: {
    std::uniform_int_distribution<int> u(0, N - 1);
    for (int i = 0; i < N; i++) {
      L[2 * i] = u(rng);
      L[2 * i + 1] = u(rng);
    }
    break;
  }
  case NORMAL1: {
    std::normal_distribution<double> nd(0, 1);
    for (int i = 0; i < N; i++) {
      L[2 * i] = ((i + 1 + (int)std::round(nd(rng))) % N + N) % N;
      L[2 * i + 1] = ((i - 1 + (int)std::round(nd(rng))) % N + N) % N;
    }
    break;
  }
  case NORMAL4: {
    std::normal_distribution<double> nd(0, 2);
    for (int i = 0; i < N; i++) {
      L[2 * i] = ((i + 1 + (int)std::round(nd(rng))) % N + N) % N;
      L[2 * i + 1] = ((i - 1 + (int)std::round(nd(rng))) % N + N) % N;
    }
    break;
  }
  case SEQUENTIAL:
    for (int i = 0; i < N; i++) {
      L[2 * i] = (i + 1) % N;
      L[2 * i + 1] = (i + 1) % N;
    }
    break;
  }
  for (int je = 0; je < N; je++) {
    d[IN_b(je, 0, B)] = L[2 * je];
    d[IN_b(je, 1, B)] = L[2 * je + 1];
  }
}
static void gen_vi(int *d, int V, int N, std::mt19937 &rng) {
  std::vector<int> p(N);
  std::iota(p.begin(), p.end(), 0);
  std::shuffle(p.begin(), p.end(), rng);
  for (int je = 0; je < N; je++)
    d[IN_v(V, je, 0, N)] = p[je];
  std::iota(p.begin(), p.end(), 0);
  std::shuffle(p.begin(), p.end(), rng);
  for (int je = 0; je < N; je++)
    d[IN_v(V, je, 1, N)] = p[je];
}
static void gen_vi_b(int *d, int N, int B, std::mt19937 &rng) {
  std::vector<int> p(N);
  std::iota(p.begin(), p.end(), 0);
  std::shuffle(p.begin(), p.end(), rng);
  for (int je = 0; je < N; je++)
    d[IN_b(je, 0, B)] = p[je];
  std::iota(p.begin(), p.end(), 0);
  std::shuffle(p.begin(), p.end(), rng);
  for (int je = 0; je < N; je++)
    d[IN_b(je, 1, B)] = p[je];
}

struct Tgt {
  const char *n, *c;
  int bb, w;
};
static const Tgt tgts[] = {{"CPU_scalar", "cpu_scalar", 64, 1},
                           {"CPU_AVX512", "cpu_avx512", 64, 8},
                           {"GPU_scalar", "gpu_scalar", 128, 1},
                           {"GPU_Warp32", "gpu_warp32", 128, 32}};
static constexpr int NT = 4;
static const int BSZ[] = {8, 16, 32, 64, 128};
static constexpr int NB = 5;

struct Res {
  const char *tgt, *sch, *dist;
  int V, B, bb, w;
  double l1r;
  Metrics m;
};

int main(int argc, char **argv) {
  int N = (argc > 1) ? atoi(argv[1]) : 81920,
      nl = (argc > 2) ? atoi(argv[2]) : 96;
  BETA = (argc > 3) ? atoi(argv[3]) : 1;
  ALPHA = (argc > 4) ? atof(argv[4]) : 0.012;
  GAMMA = (argc > 5) ? atof(argv[5]) : 1.8;
  P_NUMA = (argc > 6) ? atoi(argv[6]) : 4;
  L1_BYTES = (argc > 7) ? atoi(argv[7]) : 32768;

  fprintf(stderr, "N=%d nl=%d beta=%d alpha=%.4f gamma=%.3f P=%d L1=%d\n", N,
          nl, BETA, ALPHA, GAMMA, P_NUMA, L1_BYTES);
  std::string pp = icon_patch_path(9);
  IconEdgeData ied;
  bool hex = icon_load_patch(pp.c_str(), ied);
  if (hex)
    fprintf(stderr, "ICON Ne=%d\n", ied.n_edges_valid);

  FILE *csv = fopen("metrics.csv", "w");
  fprintf(csv, "nlev,variant,blocking,cell_dist,schedule,target,block_bytes,"
               "vector_width,"
               "mu,mu_delta,mu_delta_numa,delta_max,l1_ratio,"
               "beta,alpha,gamma,P_NUMA,T\n");

  int ndists = hex ? 5 : 4;
  Schedule scs[] = {SCHED_OMP_FOR, SCHED_OMP_COLLAPSE2};
  std::vector<Res> res;

/* Unblocked */
#pragma omp parallel
  {
    std::vector<Res> loc;
#pragma omp for schedule(dynamic, 1) collapse(4)
    for (int V = 1; V <= 4; V++)
      for (int si = 0; si < 2; si++)
        for (int di = 0; di < ndists; di++)
          for (int ti = 0; ti < NT; ti++) {
            Schedule sc = scs[si];
            LoopOrder lp = iter_order(sc, V);
            int Wv = tgts[ti].w;
            int cN = (di == 4) ? ied.n_edges_valid : N;
            int inn = (lp == KLON_FIRST) ? cN : nl;
            if (Wv > inn)
              continue;
            std::mt19937 rng(42 + di);
            std::vector<int> ci(cN * 2), vi(cN * 2);
            if (di < 4) {
              gen_ci(ci.data(), V, cN, (CellDist)di, rng);
              gen_vi(vi.data(), V, cN, rng);
            } else {
              for (int je = 0; je < cN; je++) {
                ci[IN_v(V, je, 0, cN)] = ied.cell_idx[je * 2];
                ci[IN_v(V, je, 1, cN)] = ied.cell_idx[je * 2 + 1];
                vi[IN_v(V, je, 0, cN)] = ied.vert_idx[je * 2];
                vi[IN_v(V, je, 1, cN)] = ied.vert_idx[je * 2 + 1];
              }
            }
            BLOCK_BYTES_G = tgts[ti].bb;
            Metrics m =
                compute_metrics(V, 0, Wv, cN, nl, sc, ci.data(), vi.data());
            loc.push_back({tgts[ti].c, sched_name[si], dname[di], V, 0,
                           tgts[ti].bb, Wv, 0, m});
          }
#pragma omp critical
    res.insert(res.end(), loc.begin(), loc.end());
  }

/* Blocked */
#pragma omp parallel
  {
    std::vector<Res> loc;
#pragma omp for schedule(dynamic, 1) collapse(4)
    for (int bi = 0; bi < NB; bi++)
      for (int si = 0; si < 2; si++)
        for (int di = 0; di < ndists; di++)
          for (int ti = 0; ti < NT; ti++) {
            int B = BSZ[bi];
            Schedule sc = scs[si];
            int Wv = tgts[ti].w;
            int cN = (di == 4) ? ied.n_edges_valid : N;
            if (cN % B != 0)
              continue;
            std::mt19937 rng(42 + di);
            std::vector<int> ci(cN * 2), vi(cN * 2);
            if (di < 4) {
              gen_ci_b(ci.data(), cN, B, (CellDist)di, rng);
              gen_vi_b(vi.data(), cN, B, rng);
            } else {
              for (int je = 0; je < cN; je++) {
                ci[IN_b(je, 0, B)] = ied.cell_idx[je * 2];
                ci[IN_b(je, 1, B)] = ied.cell_idx[je * 2 + 1];
                vi[IN_b(je, 0, B)] = ied.vert_idx[je * 2];
                vi[IN_b(je, 1, B)] = ied.vert_idx[je * 2 + 1];
              }
            }
            BLOCK_BYTES_G = tgts[ti].bb;
            Metrics m =
                compute_metrics(0, B, Wv, cN, nl, sc, ci.data(), vi.data());
            double ws = ((double)B * (8 * 8 + 4 * 4)) / (double)L1_BYTES;
            loc.push_back({tgts[ti].c, sched_name[si], dname[di], 0, B,
                           tgts[ti].bb, Wv, ws, m});
          }
#pragma omp critical
    res.insert(res.end(), loc.begin(), loc.end());
  }

  /* CSV + table */
  for (auto &r : res) {
    if (!r.m.T)
      continue;
    fprintf(csv,
            "%d,%d,%d,%s,%s,%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.4f,%.4f,%d,%"
            "ld\n",
            nl, r.V, r.B, r.dist, r.sch, r.tgt, r.bb, r.w, r.m.mu, r.m.mu_d,
            r.m.mu_dn, r.m.dmax, r.l1r, BETA, ALPHA, GAMMA, P_NUMA,
            (long)r.m.T);
  }
  fclose(csv);

  auto fmt = [](double v, char *b) {
    if (v < 100)
      snprintf(b, 16, "%8.3f", v);
    else if (v < 10000)
      snprintf(b, 16, "%8.1f", v);
    else
      snprintf(b, 16, "%8.0f", v);
  };
  printf("\n  N=%d nl=%d beta=%d alpha=%.4f gamma=%.3f P=%d L1=%d\n\n", N, nl,
         BETA, ALPHA, GAMMA, P_NUMA, L1_BYTES);
  printf("  %-4s %4s %-11s %-8s %-11s %8s %8s %8s %8s %8s\n", "V", "B",
         "Target", "Sched", "Dist", "mu", "mu*D", "mu*Dn", "D_max", "L1%");
  printf("  %-4s %4s %-11s %-8s %-11s %8s %8s %8s %8s %8s\n", "----", "----",
         "-----------", "--------", "-----------", "--------", "--------",
         "--------", "--------", "--------");
  for (auto &r : res) {
    if (!r.m.T)
      continue;
    char b[5][16];
    fmt(r.m.mu, b[0]);
    fmt(r.m.mu_d, b[1]);
    fmt(r.m.mu_dn, b[2]);
    fmt(r.m.dmax, b[3]);
    if (r.B > 0)
      snprintf(b[4], 16, "%7.1f%%", r.l1r * 100);
    else
      snprintf(b[4], 16, "%8s", "n/a");
    char vl[8], bl[8];
    if (r.V > 0)
      snprintf(vl, 8, "V%d", r.V);
    else
      snprintf(vl, 8, "-");
    if (r.B > 0)
      snprintf(bl, 8, "%d", r.B);
    else
      snprintf(bl, 8, "-");
    printf("  %-4s %4s %-11s %-8s %-11s %s %s %s %s %s\n", vl, bl, r.tgt, r.sch,
           r.dist, b[0], b[1], b[2], b[3], b[4]);
  }
  printf("\n");
  if (hex)
    ied.free_all();
  return 0;
}