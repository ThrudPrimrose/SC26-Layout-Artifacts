/*
 * calc_metrics.cpp -- Cost model metrics for z_v_grad_w stencil
 *
 * Each "step" processes VW consecutive inner-loop iterations (vector tile).
 * The block set for a step is the union of blocks touched by all VW lanes.
 *
 *   mu    = avg new-block count per vector step
 *   delta = avg block distance of new blocks to previous step's set
 *   sigma = avg element stride (lane-0 to lane-0) across refs / 16
 *
 * Configurations:
 *   CPU:  B=64,  VW=1 (scalar), VW=8 (AVX-512 doubles)
 *   GPU:  B=128, VW=1 (scalar), VW=32 (warp), VW=64 (AMD wavefront)
 *
 * Compile: g++ -O3 -std=c++17 -march=native calc_metrics.cpp -o calc_metrics
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <chrono>
#include <random>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <vector>

/* ================================================================ */
/*  constants                                                        */
/* ================================================================ */
static constexpr int NPROMA = 81920;
static constexpr int NLEVS[] = {90, 96};
static constexpr int N_NLEVS = 2;
static constexpr int N_REFS  = 14;

/* ================================================================ */
/*  index helpers                                                    */
/* ================================================================ */
static inline int ic(int V, int je, int jk, int N, int nlev) {
    return (V <= 2) ? (je + jk * N) : (jk + je * nlev);
}

static inline int in_idx(int V, int je, int n, int N) {
    return (V == 1 || V == 3) ? (je + n * N) : (n + je * 2);
}

/* ================================================================ */
/*  index generation                                                 */
/* ================================================================ */
enum CellDist { UNIFORM = 0, NORMAL1 = 1, NORMAL4 = 2, SEQUENTIAL = 3 };
static const char* dist_name[] = {"uniform", "normal_var1", "normal_var4", "sequential"};

static void gen_permutation(int* arr, int N, std::mt19937& rng) {
    std::iota(arr, arr + N, 0);
    std::shuffle(arr, arr + N, rng);
}

static void gen_cell_idx_logical(int* idx, int N, CellDist dist, std::mt19937& rng) {
    switch (dist) {
    case UNIFORM: {
        std::uniform_int_distribution<int> ud(0, N - 1);
        for (int i = 0; i < N; i++) { idx[i*2] = ud(rng); idx[i*2+1] = ud(rng); }
        break;
    }
    case NORMAL1: {
        std::normal_distribution<double> nd(0.0, 1.0);
        for (int i = 0; i < N; i++) {
            idx[i*2]   = (((i+1+(int)std::round(nd(rng))) % N) + N) % N;
            idx[i*2+1] = (((i-1+(int)std::round(nd(rng))) % N) + N) % N;
        }
        break;
    }
    case NORMAL4: {
        std::normal_distribution<double> nd(0.0, 2.0);
        for (int i = 0; i < N; i++) {
            idx[i*2]   = (((i+1+(int)std::round(nd(rng))) % N) + N) % N;
            idx[i*2+1] = (((i-1+(int)std::round(nd(rng))) % N) + N) % N;
        }
        break;
    }
    case SEQUENTIAL:
        for (int i = 0; i < N; i++) { idx[i*2] = (i+1)%N; idx[i*2+1] = (i+1)%N; }
        break;
    }
}

static void layout_idx(int V, int* dst, const int* logical, int N) {
    for (int je = 0; je < N; je++) {
        dst[in_idx(V, je, 0, N)] = logical[je*2];
        dst[in_idx(V, je, 1, N)] = logical[je*2+1];
    }
}

/* ================================================================ */
/*  array layout                                                     */
/* ================================================================ */
enum ArrID {
    A_OUT=0, A_VNIE, A_W, A_ZVTIE, A_ZWV,
    A_INVDUAL, A_INVPRIMAL, A_TANGENT,
    A_CIDX, A_VIDX,
    N_ARRAYS
};

static constexpr int ref_arr[N_REFS] = {
    A_OUT, A_VNIE, A_INVDUAL,
    A_W, A_W,
    A_ZVTIE, A_INVPRIMAL, A_TANGENT,
    A_ZWV, A_ZWV,
    A_CIDX, A_CIDX,
    A_VIDX, A_VIDX,
};

struct ArrInfo {
    int64_t base;
    int     elem_size;
};

static void compute_bases(ArrInfo* info, int N, int nlev, int B) {
    int64_t sizes[N_ARRAYS] = {
        (int64_t)N*nlev*8, (int64_t)N*nlev*8, (int64_t)N*nlev*8,
        (int64_t)N*nlev*8, (int64_t)N*nlev*8,
        (int64_t)N*8, (int64_t)N*8, (int64_t)N*8,
        (int64_t)N*2*4, (int64_t)N*2*4,
    };
    int esz[N_ARRAYS] = {8,8,8,8,8,8,8,8,4,4};
    int64_t off = 0;
    for (int i = 0; i < N_ARRAYS; i++) {
        info[i].base = off;
        info[i].elem_size = esz[i];
        off += sizes[i];
        off = ((off + B - 1) / B) * B;
    }
}

static inline int64_t baddr(const ArrInfo& a, int elem, int B) {
    return (a.base + (int64_t)elem * a.elem_size) / B;
}
static inline int64_t byteaddr(const ArrInfo& a, int elem) {
    return a.base + (int64_t)elem * a.elem_size;
}

/* ================================================================ */
/*  block set -- capacity for VW*14 refs (max ~900)                  */
/* ================================================================ */
struct BlockSet {
    std::vector<int64_t> d;
    BlockSet() { d.reserve(1024); }
    void clear() { d.clear(); }
    void add(int64_t v) { d.push_back(v); }  /* no dedup yet */
    void finish() {
        std::sort(d.begin(), d.end());
        d.erase(std::unique(d.begin(), d.end()), d.end());
    }
    int size() const { return (int)d.size(); }

    bool contains(int64_t v) const {
        auto it = std::lower_bound(d.begin(), d.end(), v);
        return it != d.end() && *it == v;
    }
    int64_t nearest(int64_t v) const {
        auto it = std::lower_bound(d.begin(), d.end(), v);
        int64_t best = INT64_MAX;
        if (it != d.end())
            best = std::min(best, std::abs(*it - v));
        if (it != d.begin()) {
            --it;
            best = std::min(best, std::abs(*it - v));
        }
        return best;
    }
};

/* ================================================================ */
/*  metrics                                                          */
/* ================================================================ */
struct Metrics {
    double mu;
    double delta;
    double sigma;
};

enum LoopOrder { KLON_FIRST = 0, KLEV_FIRST = 1 };
static const char* loop_name[] = {"klon_first", "klev_first"};

/* ================================================================ */
/*  compute metrics with vector-width tiling                         */
/*                                                                   */
/*  VW = vector width (number of inner-loop iters per step)          */
/*  One "step" = VW consecutive inner iterations at fixed outer      */
/*  Block set = union of all blocks across VW lanes x 14 refs        */
/*  sigma = stride of lane-0 address between consecutive steps       */
/* ================================================================ */
static Metrics compute(
    int V, int N, int nlev, int B, int VW, int loop_order,
    const int* cidx, const int* vidx)
{
    ArrInfo arr[N_ARRAYS];
    compute_bases(arr, N, nlev, B);

    int outer_n = (loop_order == KLON_FIRST) ? nlev : N;
    int inner_n = (loop_order == KLON_FIRST) ? N    : nlev;

    /* number of vector steps */
    int inner_steps = (inner_n + VW - 1) / VW;
    int64_t T_vec = (int64_t)outer_n * inner_steps;

    double sum_new   = 0.0;
    double sum_delta = 0.0;
    double sum_sigma = 0.0;

    BlockSet prev, cur;

    /* lane-0 byte addresses from previous step (for sigma) */
    int64_t prev_byte[N_REFS];
    for (int i = 0; i < N_REFS; i++) prev_byte[i] = -1;

    int64_t step = 0;
    for (int outer = 0; outer < outer_n; outer++) {
        for (int inner_base = 0; inner_base < inner_n; inner_base += VW, step++) {

            int vw_end = std::min(inner_base + VW, inner_n);

            cur.clear();
            int64_t lane0_byte[N_REFS];
            bool lane0_set = false;

            /* collect blocks from all VW lanes */
            for (int lane = inner_base; lane < vw_end; lane++) {
                int je, jk;
                if (loop_order == KLON_FIRST) {
                    jk = outer; je = lane;
                } else {
                    je = outer; jk = lane;
                }

                int ci0 = cidx[in_idx(V, je, 0, N)];
                int ci1 = cidx[in_idx(V, je, 1, N)];
                int vi0 = vidx[in_idx(V, je, 0, N)];
                int vi1 = vidx[in_idx(V, je, 1, N)];

                int elem[N_REFS] = {
                    ic(V, je,  jk, N, nlev),   /* 0  out */
                    ic(V, je,  jk, N, nlev),   /* 1  vn_ie */
                    je,                         /* 2  inv_dual */
                    ic(V, ci0, jk, N, nlev),   /* 3  w[ci0] */
                    ic(V, ci1, jk, N, nlev),   /* 4  w[ci1] */
                    ic(V, je,  jk, N, nlev),   /* 5  z_vt_ie */
                    je,                         /* 6  inv_primal */
                    je,                         /* 7  tangent */
                    ic(V, vi0, jk, N, nlev),   /* 8  z_w_v[vi0] */
                    ic(V, vi1, jk, N, nlev),   /* 9  z_w_v[vi1] */
                    in_idx(V, je, 0, N),       /* 10 cell_idx[0] */
                    in_idx(V, je, 1, N),       /* 11 cell_idx[1] */
                    in_idx(V, je, 0, N),       /* 12 vert_idx[0] */
                    in_idx(V, je, 1, N),       /* 13 vert_idx[1] */
                };

                for (int r = 0; r < N_REFS; r++)
                    cur.add(baddr(arr[ref_arr[r]], elem[r], B));

                /* record lane 0 for sigma */
                if (!lane0_set) {
                    for (int r = 0; r < N_REFS; r++)
                        lane0_byte[r] = byteaddr(arr[ref_arr[r]], elem[r]);
                    lane0_set = true;
                }
            }
            cur.finish();

            if (step == 0) {
                sum_new   += cur.size();
                sum_delta += 1.0 * cur.size();
            } else {
                /* new blocks = cur \ prev */
                int n_new = 0;
                double step_dist = 0.0;
                for (int i = 0; i < cur.size(); i++) {
                    if (!prev.contains(cur.d[i])) {
                        n_new++;
                        step_dist += (double)prev.nearest(cur.d[i]);
                    }
                }
                sum_new += n_new;
                if (n_new > 0)
                    sum_delta += step_dist / n_new;

                /* sigma: lane-0 stride between consecutive steps */
                double step_stride = 0.0;
                int n_valid = 0;
                for (int r = 0; r < N_REFS; r++) {
                    if (prev_byte[r] >= 0) {
                        step_stride += (double)std::abs(lane0_byte[r] - prev_byte[r]);
                        n_valid++;
                    }
                }
                if (n_valid > 0)
                    sum_sigma += (step_stride / n_valid) / 16.0;
            }

            std::swap(prev, cur);
            for (int r = 0; r < N_REFS; r++) prev_byte[r] = lane0_byte[r];
        }
    }

    return {sum_new / T_vec, sum_delta / T_vec, sum_sigma / T_vec};
}

/* ================================================================ */
/*  target configs                                                   */
/* ================================================================ */
struct Target {
    int         B;       /* cache line bytes */
    int         VW;      /* vector width (inner-loop tile) */
    const char* label;
};

static constexpr Target TARGETS[] = {
    { 64,  1, "cpu_scalar"},
    { 64,  8, "cpu_avx512"},
    { 64,  2, "cpu_neon"  },
    {128,  1, "gpu_scalar"},
    {128, 32, "gpu_warp32"},
    {128, 64, "gpu_wave64"},
};
static constexpr int N_TARGETS = sizeof(TARGETS) / sizeof(TARGETS[0]);

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main() {
    const int N = NPROMA;
    std::mt19937 rng(42);

    int* vert_logical = new int[N * 2];
    {
        int* tmp = new int[N];
        gen_permutation(tmp, N, rng);
        for (int i = 0; i < N; i++) vert_logical[i*2] = tmp[i];
        gen_permutation(tmp, N, rng);
        for (int i = 0; i < N; i++) vert_logical[i*2+1] = tmp[i];
        delete[] tmp;
    }

    int* cell_logical = new int[N * 2];
    int* cidx = new int[N * 2];
    int* vidx = new int[N * 2];

    FILE* fcsv = fopen("metrics.csv", "w");
    fprintf(fcsv,
        "nlev,variant,cell_dist,loop_order,target,block_bytes,vector_width,"
        "mu,delta,sigma\n");

    for (int nlev_i = 0; nlev_i < N_NLEVS; nlev_i++) {
        int nlev = NLEVS[nlev_i];

        for (int di = 0; di < 4; di++) {
            CellDist dist = (CellDist)di;
            gen_cell_idx_logical(cell_logical, N, dist, rng);

            for (int V = 1; V <= 4; V++) {
                layout_idx(V, cidx, cell_logical, N);
                layout_idx(V, vidx, vert_logical, N);

                for (int lo = 0; lo < 2; lo++) {
                    for (int ti = 0; ti < N_TARGETS; ti++) {
                        auto& tgt = TARGETS[ti];

                        auto t0 = std::chrono::steady_clock::now();
                        Metrics m = compute(V, N, nlev,
                                            tgt.B, tgt.VW, lo,
                                            cidx, vidx);
                        auto t1 = std::chrono::steady_clock::now();
                        double sec = std::chrono::duration<double>(t1-t0).count();

                        fprintf(fcsv,
                            "%d,%d,%s,%s,%s,%d,%d,%.6f,%.6f,%.6f\n",
                            nlev, V, dist_name[di], loop_name[lo],
                            tgt.label, tgt.B, tgt.VW,
                            m.mu, m.delta, m.sigma);
                        fflush(fcsv);

                        printf("nlev=%d V=%d dist=%-12s loop=%-11s "
                               "%-12s B=%3d VW=%2d  "
                               "mu=%.3f delta=%.3f sigma=%.3f  (%.1fs)\n",
                               nlev, V, dist_name[di], loop_name[lo],
                               tgt.label, tgt.B, tgt.VW,
                               m.mu, m.delta, m.sigma, sec);
                    }
                }
            }
        }
    }

    fclose(fcsv);
    delete[] vert_logical;
    delete[] cell_logical;
    delete[] cidx;
    delete[] vidx;

    printf("\nResults written to metrics.csv\n");
    return 0;
}