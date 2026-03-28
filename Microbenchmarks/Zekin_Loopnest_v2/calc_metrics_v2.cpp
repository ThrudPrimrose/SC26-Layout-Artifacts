/*
 * cost_metrics.cpp -- NUMA-aware cost-model metrics for z_v_grad_w stencil
 *
 * Implements the metrics from Section 3 (Cost Model):
 *   μ(π,Φ)  = average new-block count   (Eq. 4)
 *   Δ(π,Φ)  = NUMA-aware average block distance (Eq. 7 + NUMA weighting)
 *
 * The block distance uses a piecewise locality weight:
 *   ω(b,b*) = α   if ρ < β        (within prefetcher radius)
 *           = γ   if ν(b) ≠ ν(b*)  (cross-NUMA)
 *           = 1   otherwise         (same NUMA, beyond prefetcher)
 *
 * NUMA domain assignment ν(b) models first-touch under schedule(static):
 *   V1/V2 (klon-first layout): outer jk loop → NUMA slices along jk
 *   V3/V4 (klev-first layout): outer je loop → NUMA slices along je
 *   1D / index arrays: schedule(static) over je
 *
 * Compile:  g++ -O3 -std=c++17 -o cost_metrics cost_metrics.cpp
 * Run:      ./cost_metrics [N] [nlev] [β] [α] [γ] [P_NUMA]
 *           defaults: 8192  96    8   0.16 2.0  4
 * Output:   CSV to results_full.csv, progress on stderr
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <climits>
#include <cassert>
#include <cstring>
#include <unistd.h>

/* ============================================================
 *  NUMA-aware weighting parameters (set from argv)
 * ============================================================ */

static int    BETA    = 8;     // prefetcher radius in cache lines
static double ALPHA   = 0.16;  // page-locality discount
static double GAMMA   = 2.0;   // NUMA penalty factor
static int    P_NUMA  = 4;     // number of NUMA domains

/* ============================================================
 *  Layout index functions  (same semantics as bench_common.h)
 * ============================================================ */

inline int IC(int V, int je, int jk, int N, int nlev) {
    return (V <= 2) ? je + jk * N : jk + je * nlev;
}

inline int IN(int V, int je, int n, int N) {
    return (V == 1 || V == 3) ? je + n * N : n + je * 2;
}

/* ============================================================
 *  Array identifiers and metadata
 * ============================================================ */

enum ArrID {
    A_OUT = 0, A_VN_IE, A_W, A_Z_VT_IE, A_Z_W_V,
    A_INV_DUAL, A_INV_PRIMAL, A_TANGENT,
    A_CELL_IDX, A_VERT_IDX,
    NUM_ARR
};

static const char* arr_name[] = {
    "out","vn_ie","w","z_vt_ie","z_w_v",
    "inv_dual","inv_primal","tangent","cell_idx","vert_idx"
};

static const int ebytes[] = { 8,8,8,8,8, 8,8,8, 4,4 };

/* ============================================================
 *  Block address helper
 * ============================================================ */

static int BLOCK_BYTES = 64;

inline int64_t blk(int elem_idx, int arr) {
    return (int64_t)elem_idx * ebytes[arr] / BLOCK_BYTES;
}

/* ============================================================
 *  Schedule (OpenMP parallelization strategy)
 *
 *  Determines both the iteration order and NUMA first-touch.
 *
 *  OMP_FOR:
 *    V1/V2: outer jk parallel schedule(static), inner je seq
 *           → first-touch by jk slice
 *           → iteration order: outer jk, inner je  (KLON_FIRST)
 *    V3/V4: outer je parallel schedule(static), inner jk seq
 *           → first-touch by je slice
 *           → iteration order: outer je, inner jk  (KLEV_FIRST)
 *
 *  OMP_COLLAPSE2:
 *    All V: collapse(2) schedule(static) over (jk, je)
 *           → first-touch by linearized (jk*N + je) slice
 *           → iteration order: outer jk, inner je  (always)
 * ============================================================ */

enum Schedule { SCHED_OMP_FOR = 0, SCHED_OMP_COLLAPSE2 = 1 };
static const char* sched_name[] = {"omp_for", "omp_collapse2"};

enum LoopOrder { KLON_FIRST = 0, KLEV_FIRST = 1 };
static const char* loop_name[] = {"klon_first", "klev_first"};

/* Iteration order implied by (schedule, variant) */
inline LoopOrder iteration_order(Schedule sched, int V) {
    if (sched == SCHED_OMP_COLLAPSE2)
        return KLON_FIRST;  /* collapse(2) always linearizes as jk*N+je */
    return (V <= 2) ? KLON_FIRST : KLEV_FIRST;
}

/* ============================================================
 *  NUMA domain assignment  ν(elem_idx)
 *
 *  Depends on schedule, variant, and array type.
 * ============================================================ */

/* 2D arrays: NUMA by the parallelized dimension's slice */
inline int numa_domain_2d(Schedule sched, int V, int elem_idx,
                           int N, int nlev) {
    int64_t total = (int64_t)N * nlev;
    int64_t chunk;

    if (sched == SCHED_OMP_COLLAPSE2) {
        /* First-touch by linearized iteration: lin = jk*N + je.
         * For V1/V2: elem = je + jk*N = lin  (identity).
         * For V3/V4: elem = jk + je*nlev → jk = elem%nlev, je = elem/nlev
         *            → lin = (elem%nlev)*N + (elem/nlev). */
        int64_t lin;
        if (V <= 2) {
            lin = elem_idx;
        } else {
            int jk = elem_idx % nlev;
            int je = elem_idx / nlev;
            lin = (int64_t)jk * N + je;
        }
        chunk = (total + P_NUMA - 1) / P_NUMA;
        int d = (int)(lin / chunk);
        return (d < P_NUMA) ? d : P_NUMA - 1;
    }

    /* SCHED_OMP_FOR: parallel over outer loop */
    if (V <= 2) {
        /* outer jk schedule(static) → NUMA by jk */
        int jk = elem_idx / N;
        int chunk_jk = (nlev + P_NUMA - 1) / P_NUMA;
        int d = jk / chunk_jk;
        return (d < P_NUMA) ? d : P_NUMA - 1;
    } else {
        /* outer je schedule(static) → NUMA by je */
        int je = elem_idx / nlev;
        int chunk_je = (N + P_NUMA - 1) / P_NUMA;
        int d = je / chunk_je;
        return (d < P_NUMA) ? d : P_NUMA - 1;
    }
}

/* 1D arrays: always schedule(static) over je */
inline int numa_domain_1d(int elem_idx, int N) {
    int chunk = (N + P_NUMA - 1) / P_NUMA;
    if (chunk == 0) return 0;
    int d = elem_idx / chunk;
    return (d < P_NUMA) ? d : P_NUMA - 1;
}

/* Index arrays: schedule(static) over je, recover je from layout */
inline int numa_domain_idx(int V, int elem_idx, int N) {
    int je;
    if (V == 1 || V == 3) {
        je = elem_idx % N;     // IN = je + n*N
    } else {
        je = elem_idx / 2;     // IN = n + je*2
    }
    return numa_domain_1d(je, N);
}

/* Dispatch: 2D arrays use schedule-aware, 1D/idx are schedule-independent */
static Schedule G_SCHED;  /* current schedule, set before compute_metrics */

inline int numa_domain(int arr_id, int V, int elem_idx, int N, int nlev) {
    switch (arr_id) {
        case A_OUT: case A_VN_IE: case A_W:
        case A_Z_VT_IE: case A_Z_W_V:
            return numa_domain_2d(G_SCHED, V, elem_idx, N, nlev);
        case A_INV_DUAL: case A_INV_PRIMAL: case A_TANGENT:
            return numa_domain_1d(elem_idx, N);
        case A_CELL_IDX: case A_VERT_IDX:
            return numa_domain_idx(V, elem_idx, N);
        default: return 0;
    }
}

/* ============================================================
 *  Block-to-NUMA: given a block address for an array, recover
 *  the representative element index and return its NUMA domain.
 *
 *  block_addr = floor(elem_idx * elem_bytes / BLOCK_BYTES)
 *  → representative elem = block_addr * BLOCK_BYTES / elem_bytes
 * ============================================================ */

inline int block_numa(int arr_id, int V, int64_t block_addr,
                      int N, int nlev) {
    int elem = (int)(block_addr * BLOCK_BYTES / ebytes[arr_id]);
    return numa_domain(arr_id, V, elem, N, nlev);
}

/* ============================================================
 *  Reference descriptors
 * ============================================================ */

static constexpr int NREFS = 14;

struct Ref { int arr; int64_t block; int elem; };

inline void collect_refs(int V, int je, int jk, int N, int nlev,
                         const int* cidx, const int* vidx,
                         Ref refs[NREFS])
{
    int c2d = IC(V, je, jk, N, nlev);
    int ci0 = cidx[IN(V, je, 0, N)];
    int ci1 = cidx[IN(V, je, 1, N)];
    int vi0 = vidx[IN(V, je, 0, N)];
    int vi1 = vidx[IN(V, je, 1, N)];

    int e3  = IC(V, ci0, jk, N, nlev);
    int e4  = IC(V, ci1, jk, N, nlev);
    int e8  = IC(V, vi0, jk, N, nlev);
    int e9  = IC(V, vi1, jk, N, nlev);
    int e10 = IN(V, je, 0, N);
    int e11 = IN(V, je, 1, N);

    refs[ 0] = { A_OUT,        blk(c2d, A_OUT),        c2d };
    refs[ 1] = { A_VN_IE,      blk(c2d, A_VN_IE),      c2d };
    refs[ 2] = { A_INV_DUAL,   blk(je,  A_INV_DUAL),   je  };
    refs[ 3] = { A_W,          blk(e3,  A_W),           e3  };
    refs[ 4] = { A_W,          blk(e4,  A_W),           e4  };
    refs[ 5] = { A_Z_VT_IE,    blk(c2d, A_Z_VT_IE),    c2d };
    refs[ 6] = { A_INV_PRIMAL, blk(je,  A_INV_PRIMAL),  je  };
    refs[ 7] = { A_TANGENT,    blk(je,  A_TANGENT),      je  };
    refs[ 8] = { A_Z_W_V,      blk(e8,  A_Z_W_V),       e8  };
    refs[ 9] = { A_Z_W_V,      blk(e9,  A_Z_W_V),       e9  };
    refs[10] = { A_CELL_IDX,   blk(e10, A_CELL_IDX),    e10 };
    refs[11] = { A_CELL_IDX,   blk(e11, A_CELL_IDX),    e11 };
    refs[12] = { A_VERT_IDX,   blk(e10, A_VERT_IDX),    e10 };
    refs[13] = { A_VERT_IDX,   blk(e11, A_VERT_IDX),    e11 };
}

static const char* ref_label[] = {
    "out[c2d]","vn_ie[c2d]","inv_dual[je]",
    "w[ci0,jk]","w[ci1,jk]",
    "z_vt_ie[c2d]","inv_primal[je]","tangent[je]",
    "z_w_v[vi0,jk]","z_w_v[vi1,jk]",
    "cell_idx[je,0]","cell_idx[je,1]",
    "vert_idx[je,0]","vert_idx[je,1]"
};

/* ============================================================
 *  Per-array block set  (B_t)
 * ============================================================ */

struct BlockSet {
    std::vector<int64_t> a[NUM_ARR];

    void clear() { for (int i = 0; i < NUM_ARR; i++) a[i].clear(); }
    void add(int arr, int64_t b) { a[arr].push_back(b); }

    void finalize() {
        for (int i = 0; i < NUM_ARR; i++) {
            auto& v = a[i];
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        }
    }

    int total() const {
        int n = 0;
        for (int i = 0; i < NUM_ARR; i++) n += (int)a[i].size();
        return n;
    }

    bool contains(int arr, int64_t b) const {
        return std::binary_search(a[arr].begin(), a[arr].end(), b);
    }
};

/* ============================================================
 *  Locality weight  ω(b, b*)
 *
 *  ω = α   if |b - b*| < β         (within prefetcher radius)
 *  ω = γ   if ν(b) ≠ ν(b*)         (cross-NUMA domain)
 *  ω = 1   otherwise                (same NUMA, beyond prefetcher)
 *
 *  Priority: prefetcher wins (top-to-bottom).
 * ============================================================ */

inline double locality_weight(int64_t b, int64_t b_star,
                              int arr_id, int V, int N, int nlev) {
    int64_t dist = std::abs(b - b_star);
    if (dist < BETA)
        return ALPHA;
    int nu_b  = block_numa(arr_id, V, b,      N, nlev);
    int nu_bs = block_numa(arr_id, V, b_star,  N, nlev);
    if (nu_b != nu_bs)
        return GAMMA;
    return 1.0;
}

/* ============================================================
 *  Core metric computation
 * ============================================================ */

struct MetricResult {
    double mu_total, delta_total, sigma_total;
    int64_t T;
};

MetricResult compute_metrics(int V, int W, int N, int nlev,
                             Schedule sched,
                             const int* cidx, const int* vidx)
{
    MetricResult res = {};
    BlockSet prev, curr;
    Ref refs[NREFS];

    G_SCHED = sched;  /* set global for numa_domain dispatch */
    LoopOrder loop = iteration_order(sched, V);

    int64_t T = 0;
    double sum_new = 0, sum_rho = 0, sum_sigma = 0;

    int outer_n = (loop == KLON_FIRST) ? nlev : N;
    int inner_n = (loop == KLON_FIRST) ? N    : nlev;

    for (int outer = 0; outer < outer_n; outer++) {
        for (int inner0 = 0; inner0 + W <= inner_n; inner0 += W) {

            curr.clear();
            for (int w = 0; w < W; w++) {
                int je, jk;
                if (loop == KLON_FIRST) {
                    jk = outer;
                    je = inner0 + w;
                } else {
                    je = outer;
                    jk = inner0 + w;
                }
                collect_refs(V, je, jk, N, nlev, cidx, vidx, refs);
                for (int r = 0; r < NREFS; r++)
                    curr.add(refs[r].arr, refs[r].block);
            }
            curr.finalize();

            if (T == 0) {
                int nb = curr.total();
                sum_new += nb;
                sum_rho += 1.0;
                sum_sigma += nb;  /* σ_1 = |B_1| · 1 */
            } else {
                int    new_total = 0;
                double dist_total = 0.0;

                for (int a = 0; a < NUM_ARR; a++) {
                    for (int64_t b : curr.a[a]) {
                        if (!prev.contains(a, b)) {
                            new_total++;

                            int64_t best_raw = INT64_MAX;
                            int64_t b_star = b;
                            for (int64_t bp : prev.a[a]) {
                                int64_t d = std::abs(b - bp);
                                if (d < best_raw) {
                                    best_raw = d;
                                    b_star = bp;
                                }
                            }

                            if (best_raw == INT64_MAX) {
                                dist_total += 1.0;
                            } else {
                                double wt = locality_weight(
                                    b, b_star, a, V, N, nlev);
                                dist_total += wt * (double)best_raw;
                            }
                        }
                    }
                }

                sum_new += new_total;
                sum_rho += (new_total > 0) ? dist_total / new_total : 0.0;
                sum_sigma += dist_total;  /* σ_t = Σ ρ_w(b), NOT divided by |N_t| */
            }

            T++;
            std::swap(prev, curr);
        }
    }

    res.T = T;
    res.mu_total    = sum_new / T;
    res.delta_total = sum_rho / T;
    res.sigma_total = sum_sigma / T;
    return res;
}

/* ============================================================
 *  Index array generation  (mirrors bench_common.h)
 * ============================================================ */

enum CellDist { UNIFORM=0, NORMAL1=1, NORMAL4=2, SEQUENTIAL=3 };
static const char* dist_name[] = {"uniform","normal1","normal4","sequential"};

static void gen_cell_idx(int* dst, int V, int N,
                         CellDist dist, std::mt19937& rng)
{
    std::vector<int> logical(N * 2);
    switch (dist) {
    case UNIFORM: {
        std::uniform_int_distribution<int> ud(0, N - 1);
        for (int i = 0; i < N; i++) {
            logical[i*2+0] = ud(rng);
            logical[i*2+1] = ud(rng);
        }
        break;
    }
    case NORMAL1: {
        std::normal_distribution<double> nd(0.0, 1.0);
        for (int i = 0; i < N; i++) {
            int v0 = i + 1 + (int)std::round(nd(rng));
            int v1 = i - 1 + (int)std::round(nd(rng));
            logical[i*2+0] = ((v0 % N) + N) % N;
            logical[i*2+1] = ((v1 % N) + N) % N;
        }
        break;
    }
    case NORMAL4: {
        std::normal_distribution<double> nd(0.0, 2.0);
        for (int i = 0; i < N; i++) {
            int v0 = i + 1 + (int)std::round(nd(rng));
            int v1 = i - 1 + (int)std::round(nd(rng));
            logical[i*2+0] = ((v0 % N) + N) % N;
            logical[i*2+1] = ((v1 % N) + N) % N;
        }
        break;
    }
    case SEQUENTIAL:
        for (int i = 0; i < N; i++) {
            logical[i*2+0] = (i + 1) % N;
            logical[i*2+1] = (i + 1) % N;
        }
        break;
    }
    for (int je = 0; je < N; je++) {
        dst[IN(V, je, 0, N)] = logical[je*2+0];
        dst[IN(V, je, 1, N)] = logical[je*2+1];
    }
}

static void gen_vert_idx(int* dst, int V, int N, std::mt19937& rng) {
    std::vector<int> perm(N);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    for (int je = 0; je < N; je++)
        dst[IN(V, je, 0, N)] = perm[je];
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    for (int je = 0; je < N; je++)
        dst[IN(V, je, 1, N)] = perm[je];
}

/* ============================================================
 *  Named target configurations
 * ============================================================ */

struct Target {
    const char* name;
    int block_bytes;
    int vec_width;
};

static const Target targets[] = {
    { "CPU_scalar",    64,  1 },   // B=8 doubles, W=1
    { "CPU_NEON",      64,  2 },   // B=8, W=2 (Neoverse V2 NEON)
    { "CPU_AVX512",    64,  8 },   // B=8, W=8 (Zen 4 AVX-512)
    { "GPU_scalar",   128,  1 },   // B=16 doubles, W=1
    { "GPU_halfw",    128, 16 },   // B=16, W=16 (128B sector)
    { "GPU_warp32",   128, 32 },   // B=16, W=32 (NVIDIA warp / AMD wave32)
};
static constexpr int N_TARGETS = sizeof(targets) / sizeof(targets[0]);

/* ============================================================
 *  Main
 * ============================================================ */

int main(int argc, char** argv) {
    int N    = (argc > 1) ? atoi(argv[1]) : 8192;
    int nlev = (argc > 2) ? atoi(argv[2]) : 96;
    BETA     = (argc > 3) ? atoi(argv[3]) : (int)(sysconf(_SC_PAGESIZE) / 64);
    ALPHA    = (argc > 4) ? atof(argv[4]) : 0.16;
    GAMMA    = (argc > 5) ? atof(argv[5]) : 2.0;
    P_NUMA   = (argc > 6) ? atoi(argv[6]) : 4;

    fprintf(stderr, "NUMA-aware cost metrics for z_v_grad_w\n");
    fprintf(stderr, "  N=%d  nlev=%d  beta=%d  alpha=%.3f  gamma=%.3f  P_NUMA=%d\n",
            N, nlev, BETA, ALPHA, GAMMA, P_NUMA);
    fprintf(stderr, "  Page size: %ld B  (beta = %d CLs)\n\n",
            sysconf(_SC_PAGESIZE), BETA);

    std::mt19937 rng(42);
    std::vector<int> cidx(N * 2), vidx(N * 2);

    FILE* csv = fopen("results_full.csv", "w");
    if (!csv) { perror("fopen results_full.csv"); return 1; }

    fprintf(csv, "target,block_bytes,vec_width,variant,schedule,loop_order,dist,"
                 "T,mu,delta,sigma,beta,alpha,gamma,P_NUMA\n");

    Schedule scheds[] = { SCHED_OMP_FOR, SCHED_OMP_COLLAPSE2 };
    int n_scheds = 2;

    for (int V = 1; V <= 4; V++) {
        for (int si = 0; si < n_scheds; si++) {
            Schedule sched = scheds[si];
            LoopOrder loop = iteration_order(sched, V);

            for (int d = 0; d < 4; d++) {
                CellDist dist = (CellDist)d;
                rng.seed(42 + d);
                gen_cell_idx(cidx.data(), V, N, dist, rng);
                gen_vert_idx(vidx.data(), V, N, rng);

                for (int ti = 0; ti < N_TARGETS; ti++) {
                    const Target& tgt = targets[ti];
                    BLOCK_BYTES = tgt.block_bytes;
                    int W = tgt.vec_width;

                    int inner_n = (loop == KLON_FIRST) ? N : nlev;
                    if (W > inner_n) continue;

                    fprintf(stderr, "  %-12s V%d  %-14s  %-12s  %-10s ...",
                            tgt.name, V, sched_name[sched],
                            loop_name[loop], dist_name[d]);

                    MetricResult m = compute_metrics(V, W, N, nlev, sched,
                                                     cidx.data(), vidx.data());

                    fprintf(stderr, "  mu=%.3f  delta=%.3f  sigma=%.3f  (T=%ld)\n",
                            m.mu_total, m.delta_total, m.sigma_total, (long)m.T);

                    fprintf(csv, "%s,%d,%d,V%d,%s,%s,%s,"
                                 "%ld,%.6f,%.6f,%.6f,%d,%.4f,%.4f,%d\n",
                            tgt.name, tgt.block_bytes, W,
                            V, sched_name[sched], loop_name[loop],
                            dist_name[d],
                            (long)m.T, m.mu_total, m.delta_total,
                            m.sigma_total,
                            BETA, ALPHA, GAMMA, P_NUMA);
                }
            }
        }
    }

    fclose(csv);
    fprintf(stderr, "\nResults written to results_full.csv\n");
    return 0;
}