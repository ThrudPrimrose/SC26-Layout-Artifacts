/*
 * bench_common.h -- shared definitions for z_v_grad_w stencil benchmark
 *
 * 4 layout variants:
 *   V1: compute(je,jk), idx(je,n)  -- nothing permuted
 *   V2: compute(je,jk), idx(n,je)  -- only idx permuted
 *   V3: compute(jk,je), idx(je,n)  -- only compute permuted, loops switched
 *   V4: compute(jk,je), idx(n,je)  -- both permuted, loops switched
 *
 * 4 cell_idx distributions: uniform, normal(var=1), normal(var=4), sequential
 * vertex_idx: always a permutation of 0..N-1
 */
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <chrono>
#include <random>
#include <algorithm>
#include <numeric>
#include <cstdint>

/* ================================================================ */
/*  constants                                                        */
/* ================================================================ */
static constexpr int NPROMA   = 81920;
static constexpr int NLEVS[]  = {96};
static constexpr int N_NLEVS  = 1;
static constexpr int WARMUP   = 5;
static constexpr int NRUNS    = 100;

/* ================================================================ */
/*  2D index helpers                                                 */
/* ================================================================ */

#ifdef __HIP_PLATFORM_AMD__
#define HD __host__ __device__ __forceinline__
#else
#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#else
#define HD inline
#endif
#endif

/* compute-array:  V1,V2 -> (je,jk),   V3,V4 -> (jk,je) */
template<int V>
HD int IC(int je, int jk, int N, int nlev) {
    if constexpr (V <= 2) return je + jk * N;
    else                  return jk + je * nlev;
}

/* neighbor-table (2 neighbors per edge):
   V1,V3 -> (je,n)    V2,V4 -> (n,je) */
template<int V>
HD int IN(int je, int n, int N) {
    if constexpr (V == 1 || V == 3) return je + n * N;
    else                            return n + je * 2;
}

/* ================================================================ */
/*  fast xor fill (splitmix64, trivially parallel)                   */
/* ================================================================ */
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x  = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x  = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

static void fill_xor(double* a, size_t n, uint64_t seed) {
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) {
        uint64_t h = splitmix64(seed + i);
        a[i] = 2.0 * ((double)(h >> 11) / (double)(1ULL << 53)) - 1.0;
    }
}

/* ================================================================ */
/*  index-generation distributions                                   */
/* ================================================================ */
enum CellDist { UNIFORM = 0, NORMAL1 = 1, NORMAL4 = 2, SEQUENTIAL = 3 };
static const char* dist_name[] = {
    "uniform", "normal_var1", "normal_var4", "sequential"};

static void gen_permutation(int* arr, int N, std::mt19937& rng) {
    std::iota(arr, arr + N, 0);
    std::shuffle(arr, arr + N, rng);
}

/* logical layout: idx[je*2+n], n in {0,1} */
static void gen_cell_idx_logical(int* idx, int N,
                                 CellDist dist, std::mt19937& rng)
{
    switch (dist) {
    case UNIFORM: {
        std::uniform_int_distribution<int> ud(0, N - 1);
        for (int i = 0; i < N; i++) {
            idx[i*2+0] = ud(rng);
            idx[i*2+1] = ud(rng);
        }
        break;
    }
    case NORMAL1: {
        std::normal_distribution<double> nd(0.0, 1.0);
        for (int i = 0; i < N; i++) {
            int v0 = i + 1 + (int)std::round(nd(rng));
            int v1 = i - 1 + (int)std::round(nd(rng));
            idx[i*2+0] = ((v0 % N) + N) % N;
            idx[i*2+1] = ((v1 % N) + N) % N;
        }
        break;
    }
    case NORMAL4: {
        std::normal_distribution<double> nd(0.0, 2.0);
        for (int i = 0; i < N; i++) {
            int v0 = i + 1 + (int)std::round(nd(rng));
            int v1 = i - 1 + (int)std::round(nd(rng));
            idx[i*2+0] = ((v0 % N) + N) % N;
            idx[i*2+1] = ((v1 % N) + N) % N;
        }
        break;
    }
    case SEQUENTIAL:
        for (int i = 0; i < N; i++) {
            idx[i*2+0] = (i + 1) % N;
            idx[i*2+1] = (i + 1) % N;
        }
        break;
    }
}

/* scatter logical (je,n) into variant-specific flat layout */
template<int V>
static void layout_idx(int* dst, const int* logical, int N) {
    for (int je = 0; je < N; je++) {
        dst[IN<V>(je, 0, N)] = logical[je*2+0];
        dst[IN<V>(je, 1, N)] = logical[je*2+1];
    }
}

/* scatter canonical (je,jk) source into variant layout */
template<int V>
static void layout_2d(double* dst, const double* src, int N, int nlev) {
    #pragma omp parallel for schedule(static) collapse(2)
    for (int jk = 0; jk < nlev; jk++)
        for (int je = 0; je < N; je++)
            dst[IC<V>(je, jk, N, nlev)] = src[je + jk * N];
}

/* runtime dispatchers */
static void rearrange_idx(int V, int* dst, const int* logical, int N) {
    switch(V) {
        case 1: layout_idx<1>(dst, logical, N); break;
        case 2: layout_idx<2>(dst, logical, N); break;
        case 3: layout_idx<3>(dst, logical, N); break;
        case 4: layout_idx<4>(dst, logical, N); break;
    }
}
static void rearrange_2d(int V, double* dst, const double* src,
                          int N, int nlev) {
    switch(V) {
        case 1: layout_2d<1>(dst, src, N, nlev); break;
        case 2: layout_2d<2>(dst, src, N, nlev); break;
        case 3: layout_2d<3>(dst, src, N, nlev); break;
        case 4: layout_2d<4>(dst, src, N, nlev); break;
    }
}

/* ================================================================ */
/*  stencil body macro (used by CPU kernels)                         */
/* ================================================================ */
#define STENCIL_BODY(V)                                                    \
    int c2d = IC<V>(je, jk, N, nlev);                                     \
    int ci0 = cell_idx[IN<V>(je, 0, N)];                                  \
    int ci1 = cell_idx[IN<V>(je, 1, N)];                                  \
    int vi0 = vert_idx[IN<V>(je, 0, N)];                                  \
    int vi1 = vert_idx[IN<V>(je, 1, N)];                                  \
    out[c2d] =                                                             \
        vn_ie[c2d] * inv_dual[je] *                                        \
            (w[IC<V>(ci0, jk, N, nlev)] - w[IC<V>(ci1, jk, N, nlev)])     \
      + z_vt_ie[c2d] * inv_primal[je] * tangent[je] *                     \
            (z_w_v[IC<V>(vi0, jk, N, nlev)] -                             \
             z_w_v[IC<V>(vi1, jk, N, nlev)]);

/* ================================================================ */
/*  data generation helpers                                          */
/* ================================================================ */
struct BenchData {
    int N, nlev;
    size_t sz2d;
    /* canonical source */
    double *src_vn_ie, *src_w, *src_z_vt_ie, *src_z_w_v;
    double *inv_dual, *inv_primal, *tangent_o;
    /* variant-layout host buffers */
    double *h_vn_ie, *h_w, *h_z_vt_ie, *h_z_w_v, *h_out;
    int    *h_cidx, *h_vidx;

    void alloc(int N_, int nlev_) {
        N = N_; nlev = nlev_;
        sz2d = (size_t)N * nlev;
        src_vn_ie   = new double[sz2d];
        src_w       = new double[sz2d];
        src_z_vt_ie = new double[sz2d];
        src_z_w_v   = new double[sz2d];
        inv_dual    = new double[N];
        inv_primal  = new double[N];
        tangent_o   = new double[N];
        h_vn_ie   = new double[sz2d];
        h_w       = new double[sz2d];
        h_z_vt_ie = new double[sz2d];
        h_z_w_v   = new double[sz2d];
        h_out     = new double[sz2d];
        h_cidx    = new int[N * 2];
        h_vidx    = new int[N * 2];
    }

    void fill(int nlev_) {
        fill_xor(src_vn_ie,   sz2d, 100 + nlev_);
        fill_xor(src_w,       sz2d, 200 + nlev_);
        fill_xor(src_z_vt_ie, sz2d, 300 + nlev_);
        fill_xor(src_z_w_v,   sz2d, 400 + nlev_);
        fill_xor(inv_dual,    N,    500);
        fill_xor(inv_primal,  N,    600);
        fill_xor(tangent_o,   N,    700);
    }

    void set_variant(int V, const int* cell_logical,
                     const int* vert_logical) {
        rearrange_2d(V, h_vn_ie,   src_vn_ie,   N, nlev);
        rearrange_2d(V, h_w,       src_w,       N, nlev);
        rearrange_2d(V, h_z_vt_ie, src_z_vt_ie, N, nlev);
        rearrange_2d(V, h_z_w_v,   src_z_w_v,   N, nlev);
        rearrange_idx(V, h_cidx,   cell_logical, N);
        rearrange_idx(V, h_vidx,   vert_logical, N);
    }

    void free_all() {
        delete[] src_vn_ie;  delete[] src_w;
        delete[] src_z_vt_ie; delete[] src_z_w_v;
        delete[] inv_dual;   delete[] inv_primal; delete[] tangent_o;
        delete[] h_vn_ie;    delete[] h_w;
        delete[] h_z_vt_ie;  delete[] h_z_w_v;    delete[] h_out;
        delete[] h_cidx;     delete[] h_vidx;
    }
};

/* vertex logical indices -- generated once */
struct VertData {
    int* logical;
    void init(int N, std::mt19937& rng) {
        logical = new int[N * 2];
        int* tmp = new int[N];
        gen_permutation(tmp, N, rng);
        for (int i = 0; i < N; i++) logical[i*2+0] = tmp[i];
        gen_permutation(tmp, N, rng);
        for (int i = 0; i < N; i++) logical[i*2+1] = tmp[i];
        delete[] tmp;
    }
    void free_all() { delete[] logical; }
};
