/*
 * cost_metrics.cpp -- Cost-model metrics for the z_v_grad_w stencil
 *
 * Implements the metrics from Section 3 (Cost Model):
 *   μ(π,Φ)  = average new-block count   (Eq. 4)
 *   Δ(π,Φ)  = average block distance    (Eq. 7)
 *
 * Schedule π:  for jk in [0,nlev)  for je in [0,N)  { STENCIL_BODY }
 *
 * Vectorized variants treat W consecutive je elements as one "iteration",
 * so B_t collects blocks across all W lanes.
 *
 * Additionally reports per-reference element strides (|addr(je+1)-addr(je)|)
 * and per-array aggregated μ / Δ for detailed analysis.
 *
 * Compile:  g++ -O3 -std=c++17 -o cost_metrics cost_metrics.cpp
 * Run:      ./cost_metrics [N] [nlev]      (defaults: 8192, 96)
 * Output:   CSV on stdout, progress on stderr
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

/* ============================================================
 *  Layout index functions  (same semantics as bench_common.h)
 * ============================================================ */

// 2D compute-array index: (je, jk) -> flat
//   V1,V2: row-of-je  = je + jk*N
//   V3,V4: row-of-jk  = jk + je*nlev
inline int IC(int V, int je, int jk, int N, int nlev) {
    return (V <= 2) ? je + jk * N : jk + je * nlev;
}

// Neighbor-table index: (je, n) -> flat,  n ∈ {0,1}
//   V1,V3: (je,n) = je + n*N
//   V2,V4: (n,je) = n  + je*2
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

// Element size in bytes (double=8, int=4)
static const int ebytes[] = { 8,8,8,8,8, 8,8,8, 4,4 };

/* ============================================================
 *  Block address helper
 *  block_addr = floor(element_index * elem_bytes / B)
 * ============================================================ */

static int BLOCK_BYTES = 64;  // cache line size, set at runtime

// Return block address for element index of array
inline int64_t blk(int elem_idx, int arr) {
    return (int64_t)elem_idx * ebytes[arr] / BLOCK_BYTES;
}

/* ============================================================
 *  Reference descriptors for a single scalar point (je, jk)
 *
 *  14 memory references in the stencil body:
 *    0  out[c2d]                   (write, direct)
 *    1  vn_ie[c2d]                 (read,  direct)
 *    2  inv_dual[je]               (read,  direct)
 *    3  w[IC(ci0,jk)]              (read,  indirect via cell_idx)
 *    4  w[IC(ci1,jk)]              (read,  indirect via cell_idx)
 *    5  z_vt_ie[c2d]               (read,  direct)
 *    6  inv_primal[je]             (read,  direct)
 *    7  tangent[je]                (read,  direct)
 *    8  z_w_v[IC(vi0,jk)]          (read,  indirect via vert_idx)
 *    9  z_w_v[IC(vi1,jk)]          (read,  indirect via vert_idx)
 *   10  cell_idx[IN(je,0)]         (read,  direct)
 *   11  cell_idx[IN(je,1)]         (read,  direct)
 *   12  vert_idx[IN(je,0)]         (read,  direct)
 *   13  vert_idx[IN(je,1)]         (read,  direct)
 * ============================================================ */

static constexpr int NREFS = 14;

struct Ref { int arr; int64_t block; };

// Fills refs[0..13] for one scalar point
inline void collect_refs(int V, int je, int jk, int N, int nlev,
                         const int* cidx, const int* vidx,
                         Ref refs[NREFS])
{
    int c2d = IC(V, je, jk, N, nlev);
    int ci0 = cidx[IN(V, je, 0, N)];
    int ci1 = cidx[IN(V, je, 1, N)];
    int vi0 = vidx[IN(V, je, 0, N)];
    int vi1 = vidx[IN(V, je, 1, N)];

    refs[ 0] = { A_OUT,        blk(c2d, A_OUT)       };
    refs[ 1] = { A_VN_IE,      blk(c2d, A_VN_IE)     };
    refs[ 2] = { A_INV_DUAL,   blk(je,  A_INV_DUAL)  };
    refs[ 3] = { A_W,          blk(IC(V, ci0, jk, N, nlev), A_W) };
    refs[ 4] = { A_W,          blk(IC(V, ci1, jk, N, nlev), A_W) };
    refs[ 5] = { A_Z_VT_IE,    blk(c2d, A_Z_VT_IE)   };
    refs[ 6] = { A_INV_PRIMAL, blk(je,  A_INV_PRIMAL) };
    refs[ 7] = { A_TANGENT,    blk(je,  A_TANGENT)    };
    refs[ 8] = { A_Z_W_V,      blk(IC(V, vi0, jk, N, nlev), A_Z_W_V) };
    refs[ 9] = { A_Z_W_V,      blk(IC(V, vi1, jk, N, nlev), A_Z_W_V) };
    refs[10] = { A_CELL_IDX,   blk(IN(V, je, 0, N), A_CELL_IDX) };
    refs[11] = { A_CELL_IDX,   blk(IN(V, je, 1, N), A_CELL_IDX) };
    refs[12] = { A_VERT_IDX,   blk(IN(V, je, 0, N), A_VERT_IDX) };
    refs[13] = { A_VERT_IDX,   blk(IN(V, je, 1, N), A_VERT_IDX) };
}

// Also collect raw element indices (for stride computation)
inline void collect_elem_indices(int V, int je, int jk, int N, int nlev,
                                 const int* cidx, const int* vidx,
                                 int elems[NREFS])
{
    int c2d = IC(V, je, jk, N, nlev);
    int ci0 = cidx[IN(V, je, 0, N)];
    int ci1 = cidx[IN(V, je, 1, N)];
    int vi0 = vidx[IN(V, je, 0, N)];
    int vi1 = vidx[IN(V, je, 1, N)];

    elems[ 0] = c2d;
    elems[ 1] = c2d;
    elems[ 2] = je;
    elems[ 3] = IC(V, ci0, jk, N, nlev);
    elems[ 4] = IC(V, ci1, jk, N, nlev);
    elems[ 5] = c2d;
    elems[ 6] = je;
    elems[ 7] = je;
    elems[ 8] = IC(V, vi0, jk, N, nlev);
    elems[ 9] = IC(V, vi1, jk, N, nlev);
    elems[10] = IN(V, je, 0, N);
    elems[11] = IN(V, je, 1, N);
    elems[12] = IN(V, je, 0, N);
    elems[13] = IN(V, je, 1, N);
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
 *  Per-array block set  (used as B_t)
 *
 *  Stored as sorted, deduplicated vectors per array.
 *  Sizes are tiny (~2-4 entries per array per scalar iteration)
 *  so linear operations are fast.
 * ============================================================ */

struct BlockSet {
    std::vector<int64_t> a[NUM_ARR];   // blocks for each array

    void clear() { for (int i = 0; i < NUM_ARR; i++) a[i].clear(); }

    void add(int arr, int64_t b) { a[arr].push_back(b); }

    // Sort + deduplicate each array's block list
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
 *  Core metric computation
 *
 *  For a given (variant V, vector width W), iterates through
 *  the full loop nest and accumulates:
 *    mu_total    = average |N_t|  across all vector-iterations
 *    delta_total = average rho_bar_t
 *    mu_arr[a]   = per-array average new-block count
 *    delta_arr[a]= per-array average block distance
 * ============================================================ */

struct MetricResult {
    double mu_total, delta_total;
    double mu_arr[NUM_ARR];
    double delta_arr[NUM_ARR];
    int64_t T;   // number of vector-iterations
};

MetricResult compute_metrics(int V, int W, int N, int nlev,
                             const int* cidx, const int* vidx)
{
    MetricResult res = {};
    BlockSet prev, curr;
    Ref refs[NREFS];

    int64_t T = 0;
    double sum_new = 0, sum_rho = 0;
    double sum_new_a[NUM_ARR] = {}, sum_rho_a[NUM_ARR] = {};

    for (int jk = 0; jk < nlev; jk++) {
        for (int je0 = 0; je0 + W <= N; je0 += W) {

            /* --- Build B_t: collect blocks across all W lanes --- */
            curr.clear();
            for (int w = 0; w < W; w++) {
                collect_refs(V, je0 + w, jk, N, nlev, cidx, vidx, refs);
                for (int r = 0; r < NREFS; r++)
                    curr.add(refs[r].arr, refs[r].block);
            }
            curr.finalize();

            /* --- First iteration: N_1 = B_1, rho_1(b) = 1 --- */
            if (T == 0) {
                int nb = curr.total();
                sum_new += nb;
                sum_rho += 1.0;  // rho_bar_1 = 1 by convention
                for (int a = 0; a < NUM_ARR; a++) {
                    int na = (int)curr.a[a].size();
                    sum_new_a[a] += na;
                    sum_rho_a[a] += (na > 0) ? 1.0 : 0.0;
                }
            } else {
                /* --- N_t = B_t \ B_{t-1} and distance rho_t --- */
                int    new_total = 0;
                double dist_total = 0.0;

                for (int a = 0; a < NUM_ARR; a++) {
                    int    new_a = 0;
                    double dist_a = 0.0;

                    for (int64_t b : curr.a[a]) {
                        if (!prev.contains(a, b)) {
                            new_a++;
                            // rho_t(b) = min distance to any block in B_{t-1}[a]
                            int64_t min_d = INT64_MAX;
                            for (int64_t bp : prev.a[a])
                                min_d = std::min(min_d, std::abs(b - bp));
                            // If no previous blocks for this array, use 1
                            dist_a += (min_d == INT64_MAX) ? 1.0 : (double)min_d;
                        }
                    }

                    sum_new_a[a] += new_a;
                    // Per-array rho_bar: mean distance of new blocks
                    sum_rho_a[a] += (new_a > 0) ? dist_a / new_a : 0.0;

                    new_total += new_a;
                    dist_total += dist_a;
                }

                sum_new += new_total;
                // rho_bar_t = mean distance over all new blocks (Eq. 6)
                sum_rho += (new_total > 0) ? dist_total / new_total : 0.0;
            }

            T++;
            std::swap(prev, curr);
        }
    }

    /* --- Normalize: mu = (1/T) sum |N_t|,  delta = (1/T) sum rho_bar_t --- */
    res.T = T;
    res.mu_total    = sum_new / T;
    res.delta_total = sum_rho / T;
    for (int a = 0; a < NUM_ARR; a++) {
        res.mu_arr[a]    = sum_new_a[a] / T;
        res.delta_arr[a] = sum_rho_a[a] / T;
    }
    return res;
}

/* ============================================================
 *  Per-reference stride computation
 *
 *  For each of the 14 references, compute the average absolute
 *  element-index change between consecutive je iterations
 *  (within the same jk row).  This captures spatial locality
 *  at the element level, independent of blocking.
 * ============================================================ */

struct StrideResult {
    double avg_stride[NREFS];  // average |elem(je+1) - elem(je)| per ref
};

StrideResult compute_strides(int V, int N, int nlev,
                             const int* cidx, const int* vidx)
{
    StrideResult res = {};
    double sum[NREFS] = {};
    int64_t count = 0;

    int prev_elem[NREFS], curr_elem[NREFS];

    for (int jk = 0; jk < nlev; jk++) {
        // First je in this row: just record, no stride yet
        collect_elem_indices(V, 0, jk, N, nlev, cidx, vidx, prev_elem);

        for (int je = 1; je < N; je++) {
            collect_elem_indices(V, je, jk, N, nlev, cidx, vidx, curr_elem);
            for (int r = 0; r < NREFS; r++)
                sum[r] += std::abs(curr_elem[r] - prev_elem[r]);
            memcpy(prev_elem, curr_elem, sizeof(prev_elem));
            count++;
        }
    }

    for (int r = 0; r < NREFS; r++)
        res.avg_stride[r] = sum[r] / count;

    return res;
}

/* ============================================================
 *  Index array generation  (mirrors bench_common.h logic)
 * ============================================================ */

enum CellDist { UNIFORM=0, NORMAL1=1, NORMAL4=2, SEQUENTIAL=3 };
static const char* dist_name[] = {"uniform","normal1","normal4","sequential"};

// Generate cell_idx in logical (je*2+n) format, then scatter into variant layout
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

    // Scatter logical[je*2+n] -> dst[IN(V, je, n, N)]
    for (int je = 0; je < N; je++) {
        dst[IN(V, je, 0, N)] = logical[je*2+0];
        dst[IN(V, je, 1, N)] = logical[je*2+1];
    }
}

// Generate vert_idx as two independent permutations
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
 *  Main
 * ============================================================ */

int main(int argc, char** argv) {
    int N    = (argc > 1) ? atoi(argv[1]) : 8192;
    int nlev = (argc > 2) ? atoi(argv[2]) : 96;

    int block_sizes[] = {64, 128};
    int n_bsizes = 2;

    int vec_widths[] = {1, 2, 8, 32, 64};
    int n_widths = 5;

    std::mt19937 rng(42);

    // Allocate index arrays (big enough for any variant)
    std::vector<int> cidx(N * 2), vidx(N * 2);

    /* --- CSV header --- */
    printf("block_bytes,variant,dist,vector_width,T,"
           "mu_total,delta_total");
    for (int a = 0; a < NUM_ARR; a++)
        printf(",mu_%s,delta_%s", arr_name[a], arr_name[a]);
    printf("\n");

    /* --- Open CSV output --- */
    FILE* csv = fopen("results_full.csv", "w");
    if (!csv) { perror("fopen results_full.csv"); return 1; }

    /* --- CSV header --- */
    fprintf(csv, "block_bytes,variant,dist,vec_width,T,"
                 "mu_total,delta_total");
    for (int a = 0; a < NUM_ARR; a++)
        fprintf(csv, ",mu_%s,delta_%s", arr_name[a], arr_name[a]);
    fprintf(csv, "\n");

    /* --- Main loop over configurations --- */
    for (int bi = 0; bi < n_bsizes; bi++) {
      BLOCK_BYTES = block_sizes[bi];
      fprintf(stderr, "=== Block size B=%d bytes ===\n", BLOCK_BYTES);

      for (int V = 1; V <= 4; V++) {
        for (int d = 0; d < 4; d++) {
            CellDist dist = (CellDist)d;

            rng.seed(42 + d);
            gen_cell_idx(cidx.data(), V, N, dist, rng);
            gen_vert_idx(vidx.data(), V, N, rng);

            for (int wi = 0; wi < n_widths; wi++) {
                int W = vec_widths[wi];

                fprintf(stderr, "  B=%-3d V%d  dist=%-10s  W=%-3d ...",
                        BLOCK_BYTES, V, dist_name[d], W);

                MetricResult m = compute_metrics(V, W, N, nlev,
                                                 cidx.data(), vidx.data());

                fprintf(stderr, "  mu=%.3f  delta=%.3f  (T=%ld)\n",
                        m.mu_total, m.delta_total, (long)m.T);

                /* CSV row */
                fprintf(csv, "%d,V%d,%s,%d,%ld,%.6f,%.6f",
                        BLOCK_BYTES, V, dist_name[d], W, (long)m.T,
                        m.mu_total, m.delta_total);
                for (int a = 0; a < NUM_ARR; a++)
                    fprintf(csv, ",%.6f,%.6f", m.mu_arr[a], m.delta_arr[a]);
                fprintf(csv, "\n");
            }
        }
      }
    } /* block sizes */

    fclose(csv);

    /* --- Per-reference strides (independent of block size) --- */
    fprintf(stderr, "\n=== Per-reference average element strides ===\n");
    fprintf(stderr, "%-4s %-12s", "V", "dist");
    for (int r = 0; r < NREFS; r++)
        fprintf(stderr, "  %16s", ref_label[r]);
    fprintf(stderr, "\n");

    for (int V = 1; V <= 4; V++) {
        for (int d = 0; d < 4; d++) {
            CellDist dist = (CellDist)d;
            rng.seed(42 + d);
            gen_cell_idx(cidx.data(), V, N, dist, rng);
            gen_vert_idx(vidx.data(), V, N, rng);

            StrideResult s = compute_strides(V, N, nlev,
                                             cidx.data(), vidx.data());

            fprintf(stderr, "V%-3d %-12s", V, dist_name[d]);
            for (int r = 0; r < NREFS; r++)
                fprintf(stderr, "  %16.2f", s.avg_stride[r]);
            fprintf(stderr, "\n");
        }
    }

    return 0;
}