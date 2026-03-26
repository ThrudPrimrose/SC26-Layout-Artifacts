/*
 * bench_common.h -- NUMA-aware shared definitions for z_v_grad_w benchmark
 *
 * All arrays allocated via mmap (unfaulted) + MADV_HUGEPAGE.
 * Page placement is driven entirely by first-touch under OpenMP
 * schedule(static), matching the compute kernel's access pattern.
 *
 * Before switching to a kernel with a different parallelization
 * strategy, call redistribute_for_*() to copy arrays into fresh
 * buffers whose pages are faulted by the correct threads.
 *
 * 4 layout variants:
 *   V1: compute(je,jk), idx(je,n)
 *   V2: compute(je,jk), idx(n,je)
 *   V3: compute(jk,je), idx(je,n)
 *   V4: compute(jk,je), idx(n,je)
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
#include <sys/mman.h>
#include <omp.h>

/* ================================================================ */
/*  constants                                                        */
/* ================================================================ */
static constexpr int NPROMA   = 81920;
static constexpr int NLEVS[]  = {96};
static constexpr int N_NLEVS  = 1;
static constexpr int WARMUP   = 5;
static constexpr int NRUNS    = 100;

/* ================================================================ */
/*  NUMA allocation primitives                                       */
/* ================================================================ */

/* Reserve virtual pages without faulting.  First parallel touch
   decides NUMA placement.  Hugepages cut TLB pressure. */
template<typename T>
static T* numa_alloc_unfaulted(size_t count)
{
    size_t bytes = count * sizeof(T);
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); exit(1); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return static_cast<T*>(p);
}

template<typename T>
static void numa_dealloc(T* p, size_t count)
{
    if (p) munmap(p, count * sizeof(T));
}

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
/*  Schedule types -- describes how a compute kernel iterates        */
/* ================================================================ */
enum SchedKind {
    SCHED_JE_OUTER,       /* #pragma omp for schedule(static) over je  */
    SCHED_JK_OUTER,       /* #pragma omp for schedule(static) over jk  */
    SCHED_COLLAPSE2,      /* #pragma omp for collapse(2) schedule(static)
                             outer jk, inner je  (default collapse order) */
};

/* ================================================================ */
/*  Redistribution kernels                                           */
/*                                                                   */
/*  Copy src[] into a fresh mmap'd buffer.  The parallel-for that    */
/*  does the copy uses the TARGET schedule, so first-touch places    */
/*  every page on the NUMA node that will access it in the compute   */
/*  kernel.                                                          */
/*                                                                   */
/*  For 2D arrays the template parameter V controls the index        */
/*  function IC<V>(je,jk,N,nlev) so we touch the same addresses     */
/*  the stencil will.                                                */
/* ================================================================ */

/* -- 2D compute arrays (size N * nlev) -- */
template<int V>
static double* redistribute_2d(const double* src, int N, int nlev,
                                SchedKind target)
{
    size_t total = (size_t)N * nlev;
    double* dst = numa_alloc_unfaulted<double>(total);

    switch (target) {

    case SCHED_COLLAPSE2:
        #pragma omp parallel for collapse(2) schedule(static)
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) {
                int idx = IC<V>(je, jk, N, nlev);
                dst[idx] = src[idx];
            }
        break;

    case SCHED_JE_OUTER:
        #pragma omp parallel for schedule(static)
        for (int je = 0; je < N; je++)
            for (int jk = 0; jk < nlev; jk++) {
                int idx = IC<V>(je, jk, N, nlev);
                dst[idx] = src[idx];
            }
        break;

    case SCHED_JK_OUTER:
        #pragma omp parallel for schedule(static)
        for (int jk = 0; jk < nlev; jk++)
            for (int je = 0; je < N; je++) {
                int idx = IC<V>(je, jk, N, nlev);
                dst[idx] = src[idx];
            }
        break;
    }

    return dst;
}

/* Runtime dispatcher (no template V) */
static double* redistribute_2d_rt(int V, const double* src, int N, int nlev,
                                   SchedKind target)
{
    switch (V) {
        case 1: return redistribute_2d<1>(src, N, nlev, target);
        case 2: return redistribute_2d<2>(src, N, nlev, target);
        case 3: return redistribute_2d<3>(src, N, nlev, target);
        case 4: return redistribute_2d<4>(src, N, nlev, target);
    }
    return nullptr;
}

/* -- 1D arrays indexed by je only (size N) --
   Distribute pages to match the je-iteration of the target schedule.
   For collapse(2) over (jk,je), the first N iterations (jk=0) cover
   all je values, so a flat schedule(static) over je gives the same
   first-touch mapping.  For je-outer, same thing.
   For jk-outer, every thread touches all je, so we just spread evenly. */
static double* redistribute_1d(const double* src, int N, SchedKind /*target*/)
{
    double* dst = numa_alloc_unfaulted<double>(N);
    #pragma omp parallel for schedule(static)
    for (int je = 0; je < N; je++)
        dst[je] = src[je];
    return dst;
}

/* -- Index arrays (size N*2) -- */
template<int V>
static int* redistribute_idx(const int* src, int N, SchedKind /*target*/)
{
    int* dst = numa_alloc_unfaulted<int>(N * 2);
    /* Touch in je order, same as the stencil's je-loop */
    #pragma omp parallel for schedule(static)
    for (int je = 0; je < N; je++) {
        dst[IN<V>(je, 0, N)] = src[IN<V>(je, 0, N)];
        dst[IN<V>(je, 1, N)] = src[IN<V>(je, 1, N)];
    }
    return dst;
}

static int* redistribute_idx_rt(int V, const int* src, int N,
                                 SchedKind target)
{
    switch (V) {
        case 1: return redistribute_idx<1>(src, N, target);
        case 2: return redistribute_idx<2>(src, N, target);
        case 3: return redistribute_idx<3>(src, N, target);
        case 4: return redistribute_idx<4>(src, N, target);
    }
    return nullptr;
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

/* First-touch aware fill: pages land on the right NUMA node */
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
/*  data generation helpers -- NUMA-aware                            */
/* ================================================================ */
struct BenchData {
    int N, nlev;
    size_t sz2d;

    /* canonical source (staging, allocated once, layout-agnostic) */
    double *src_vn_ie, *src_w, *src_z_vt_ie, *src_z_w_v;
    double *inv_dual, *inv_primal, *tangent_o;

    /* variant-layout host buffers -- NUMA-placed via first-touch */
    double *h_vn_ie, *h_w, *h_z_vt_ie, *h_z_w_v, *h_out;
    int    *h_cidx, *h_vidx;

    /* track current schedule so we know when to redistribute */
    SchedKind current_sched;
    int       current_variant;

    void alloc(int N_, int nlev_) {
        N = N_; nlev = nlev_;
        sz2d = (size_t)N * nlev;

        /* Canonical sources: small staging area, schedule doesn't matter.
           Use mmap + first-touch with flat schedule(static). */
        src_vn_ie   = numa_alloc_unfaulted<double>(sz2d);
        src_w       = numa_alloc_unfaulted<double>(sz2d);
        src_z_vt_ie = numa_alloc_unfaulted<double>(sz2d);
        src_z_w_v   = numa_alloc_unfaulted<double>(sz2d);
        inv_dual    = numa_alloc_unfaulted<double>(N);
        inv_primal  = numa_alloc_unfaulted<double>(N);
        tangent_o   = numa_alloc_unfaulted<double>(N);

        /* Variant buffers start null, allocated on first set_variant */
        h_vn_ie = h_w = h_z_vt_ie = h_z_w_v = h_out = nullptr;
        h_cidx = h_vidx = nullptr;
        current_sched   = SCHED_COLLAPSE2;
        current_variant = 0;
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

    /* Free only the variant (NUMA-placed) buffers */
    void free_variant() {
        numa_dealloc(h_vn_ie,   sz2d);
        numa_dealloc(h_w,       sz2d);
        numa_dealloc(h_z_vt_ie, sz2d);
        numa_dealloc(h_z_w_v,   sz2d);
        numa_dealloc(h_out,     sz2d);
        numa_dealloc(h_cidx,    (size_t)N * 2);
        numa_dealloc(h_vidx,    (size_t)N * 2);
        h_vn_ie = h_w = h_z_vt_ie = h_z_w_v = h_out = nullptr;
        h_cidx = h_vidx = nullptr;
    }

    /* Set layout variant + initial NUMA placement for a given schedule.
       Allocates fresh mmap buffers and first-touch populates them
       using the target schedule's loop structure. */
    void set_variant(int V, const int* cell_logical,
                     const int* vert_logical,
                     SchedKind sched = SCHED_COLLAPSE2) {
        free_variant();

        /* Allocate staging (non-NUMA, just for rearranging) */
        double* tmp = numa_alloc_unfaulted<double>(sz2d);
        int*    itmp = numa_alloc_unfaulted<int>((size_t)N * 2);

        /* -- 2D arrays: rearrange into layout V, then redistribute -- */
        auto do_2d = [&](double* src) -> double* {
            rearrange_2d(V, tmp, src, N, nlev);
            double* placed = redistribute_2d_rt(V, tmp, N, nlev, sched);
            return placed;
        };

        h_vn_ie   = do_2d(src_vn_ie);
        h_w       = do_2d(src_w);
        h_z_vt_ie = do_2d(src_z_vt_ie);
        h_z_w_v   = do_2d(src_z_w_v);

        /* output array: just allocate + first-touch with target schedule */
        h_out = redistribute_2d_rt(V, tmp /* content doesn't matter */,
                                   N, nlev, sched);

        /* -- index arrays -- */
        rearrange_idx(V, itmp, cell_logical, N);
        h_cidx = redistribute_idx_rt(V, itmp, N, sched);

        rearrange_idx(V, itmp, vert_logical, N);
        h_vidx = redistribute_idx_rt(V, itmp, N, sched);

        /* -- 1D arrays: redistribute for je-local access -- */
        double* old_dual    = inv_dual;
        double* old_primal  = inv_primal;
        double* old_tangent = tangent_o;
        inv_dual   = redistribute_1d(old_dual,    N, sched);
        inv_primal = redistribute_1d(old_primal,  N, sched);
        tangent_o  = redistribute_1d(old_tangent, N, sched);
        numa_dealloc(old_dual,    (size_t)N);
        numa_dealloc(old_primal,  (size_t)N);
        numa_dealloc(old_tangent, (size_t)N);

        numa_dealloc(tmp,  sz2d);
        numa_dealloc(itmp, (size_t)N * 2);

        current_sched   = sched;
        current_variant = V;
    }

    /* Redistribute all variant buffers in-place for a new schedule.
       Call this BEFORE running a kernel with a different schedule.
       The layout variant stays the same; only NUMA placement changes. */
    void change_schedule(SchedKind new_sched) {
        if (new_sched == current_sched) return;

        int V = current_variant;
        printf("  [NUMA] redistributing pages: %s -> %s\n",
               new_sched == SCHED_COLLAPSE2 ? "collapse2" :
               new_sched == SCHED_JE_OUTER  ? "je_outer"  : "jk_outer",
               new_sched == SCHED_COLLAPSE2 ? "collapse2" :
               new_sched == SCHED_JE_OUTER  ? "je_outer"  : "jk_outer");

        auto swap_2d = [&](double*& buf) {
            double* old = buf;
            buf = redistribute_2d_rt(V, old, N, nlev, new_sched);
            numa_dealloc(old, sz2d);
        };

        swap_2d(h_vn_ie);
        swap_2d(h_w);
        swap_2d(h_z_vt_ie);
        swap_2d(h_z_w_v);
        swap_2d(h_out);

        auto swap_idx = [&](int*& buf) {
            int* old = buf;
            buf = redistribute_idx_rt(V, old, N, new_sched);
            numa_dealloc(old, (size_t)N * 2);
        };

        swap_idx(h_cidx);
        swap_idx(h_vidx);

        /* 1D arrays */
        auto swap_1d = [&](double*& buf) {
            double* old = buf;
            buf = redistribute_1d(old, N, new_sched);
            numa_dealloc(old, (size_t)N);
        };

        swap_1d(inv_dual);
        swap_1d(inv_primal);
        swap_1d(tangent_o);

        current_sched = new_sched;
    }

    void free_all() {
        free_variant();
        numa_dealloc(src_vn_ie,   sz2d);
        numa_dealloc(src_w,       sz2d);
        numa_dealloc(src_z_vt_ie, sz2d);
        numa_dealloc(src_z_w_v,   sz2d);
        /* 1D arrays may have been swapped, only free if non-null */
        numa_dealloc(inv_dual,    (size_t)N);
        numa_dealloc(inv_primal,  (size_t)N);
        numa_dealloc(tangent_o,   (size_t)N);
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