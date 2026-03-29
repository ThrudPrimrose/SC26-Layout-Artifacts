/*
 * amx_gemm_cost.cpp -- Cost model for GEMM with atomic access shapes
 *
 * Atomic access shape: the minimal contiguous region the hardware
 * consumes per memory instruction. The cost model counts cache
 * blocks touched per atomic load, not per element.
 *
 *   Target          Atomic shape          Size     Cache lines (best/worst)
 *   ──────────────  ────────────────────  ───────  ────────────────────────
 *   Scalar          [1] element           1-8 B    1 / 1
 *   AVX-512 fp64    [8] contiguous fp64   64 B     1 / 1 (contiguous), 8 (gather)
 *   AVX-512 int8    [64] contiguous int8  64 B     1 / 1 (contiguous), 64 (gather)
 *   AMX A-tile row  [64] contiguous int8  64 B     1 per row, 16 rows
 *   AMX B-tile grp  [16][4] VNNI int8     64 B     1 per group, 16 groups
 *   GPU warp fp32   [32] coalesced fp32   128 B    1 sector / 32 sectors
 *
 * The cost model computes:
 *   mu    = avg new cache blocks per step
 *   delta = avg min distance of new blocks to previous set
 *   lines_per_load = avg cache lines touched per atomic load (utilization)
 *
 * Compile: g++ -O3 -std=c++17 -o amx_gemm_cost amx_gemm_cost.cpp
 * Run:     ./amx_gemm_cost [M] [N] [K] [block_bytes]
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <climits>
#include <functional>

static int BLOCK_BYTES = 64;

static constexpr int TILE_M = 16;
static constexpr int TILE_N = 16;
static constexpr int TILE_K = 64;

inline int64_t blk(int64_t byte_offset) {
    return byte_offset / BLOCK_BYTES;
}

/* ═══════════════════════════════════════════════════════════
 *  Layout address functions
 * ═══════════════════════════════════════════════════════════ */

/* Row-major A[M][K] int8 */
inline int64_t addr_A_rm(int m, int k, int K) {
    return (int64_t)m * K + k;
}

/* Row-major B[K][N] int8 */
inline int64_t addr_B_rm(int k, int n, int N) {
    return (int64_t)k * N + n;
}

/* Blocked B: tiles [K/64][N/16], row-major within tile [64][16] */
inline int64_t addr_B_blocked(int k, int n, int K, int N) {
    int64_t tile_idx = (int64_t)(k/TILE_K) * (N/TILE_N) + (n/TILE_N);
    return tile_idx * (TILE_K * TILE_N) + (k%TILE_K) * TILE_N + (n%TILE_N);
}

/* VNNI B: tiles [K/64][N/16], VNNI within tile [16][16][4] */
inline int64_t addr_B_vnni(int k, int n, int K, int N) {
    int64_t tile_idx = (int64_t)(k/TILE_K) * (N/TILE_N) + (n/TILE_N);
    int ki = k % TILE_K, ni = n % TILE_N;
    return tile_idx * (TILE_K * TILE_N) + (ki/4) * (TILE_N*4) + ni*4 + (ki%4);
}

/* Row-major C[M][N] int32 */
inline int64_t addr_C_rm(int m, int n, int N) {
    return ((int64_t)m * N + n) * 4;
}

/* ═══════════════════════════════════════════════════════════
 *  Atomic access shape: describes one hardware memory op
 * ═══════════════════════════════════════════════════════════ */

struct AtomicLoad {
    std::vector<int64_t> byte_addrs;  /* addresses of all elements in this load */
};

/*
 * Collect atomic loads for one GEMM step (one TDPBUSD).
 *
 * AMX issues:
 *   - 16 row-loads for A tile, each reading 64 contiguous int8
 *   - 16 group-loads for B tile, each reading 64 bytes
 *     (VNNI: contiguous 64B; row-major: 4 rows × 16 cols strided)
 *   - 16 row-loads for C tile, each reading 16 contiguous int32 = 64B
 *
 * For scalar/AVX targets, each "load" is a single vector register fill.
 */

enum BLayout { B_ROWMAJOR=0, B_BLOCKED=1, B_VNNI=2 };
static const char* blayout_name[] = {"row-major", "blocked", "VNNI"};

static void collect_amx_atomic_loads(
    int mt, int nt, int kt, int K, int N,
    BLayout b_layout,
    int64_t base_A, int64_t base_B, int64_t base_C,
    std::vector<AtomicLoad>& loads)
{
    loads.clear();

    /* A tile: 16 row-loads, each 64 contiguous int8 */
    for (int m = 0; m < TILE_M; m++) {
        AtomicLoad ld;
        for (int k = 0; k < TILE_K; k++)
            ld.byte_addrs.push_back(base_A + addr_A_rm(
                mt*TILE_M + m, kt*TILE_K + k, K));
        loads.push_back(std::move(ld));
    }

    /* B tile: depends on layout */
    if (b_layout == B_VNNI) {
        /* VNNI: 16 groups of 64 contiguous bytes each */
        for (int g = 0; g < 16; g++) {
            AtomicLoad ld;
            for (int ni = 0; ni < TILE_N; ni++)
                for (int ki = 0; ki < 4; ki++)
                    ld.byte_addrs.push_back(base_B + addr_B_vnni(
                        kt*TILE_K + g*4 + ki, nt*TILE_N + ni, K, N));
            loads.push_back(std::move(ld));
        }
    } else if (b_layout == B_BLOCKED) {
        /* Blocked: 16 groups, each 64 contiguous bytes (row-major within tile) */
        for (int g = 0; g < 16; g++) {
            AtomicLoad ld;
            for (int ki = g*4; ki < g*4+4; ki++)
                for (int ni = 0; ni < TILE_N; ni++)
                    ld.byte_addrs.push_back(base_B + addr_B_blocked(
                        kt*TILE_K + ki, nt*TILE_N + ni, K, N));
            loads.push_back(std::move(ld));
        }
    } else {
        /* Row-major: 64 row-loads of 16 bytes each (strided!) */
        for (int k = 0; k < TILE_K; k++) {
            AtomicLoad ld;
            for (int n = 0; n < TILE_N; n++)
                ld.byte_addrs.push_back(base_B + addr_B_rm(
                    kt*TILE_K + k, nt*TILE_N + n, N));
            loads.push_back(std::move(ld));
        }
    }

    /* C tile: 16 row-loads of 16 int32 = 64 bytes each */
    for (int m = 0; m < TILE_M; m++) {
        AtomicLoad ld;
        for (int n = 0; n < TILE_N; n++)
            ld.byte_addrs.push_back(base_C + addr_C_rm(
                mt*TILE_M + m, nt*TILE_N + n, N));
        loads.push_back(std::move(ld));
    }
}

/*
 * For scalar CPU: each element is its own atomic load.
 * For AVX-512 int8: each vector of 64 contiguous elements is one atomic load.
 * We parametrize by vector_width (1 = scalar, 64 = AVX-512 int8, etc.)
 */
static void collect_scalar_atomic_loads(
    int mt, int nt, int kt, int K, int N,
    BLayout b_layout, int vector_width,
    int64_t base_A, int64_t base_B, int64_t base_C,
    std::vector<AtomicLoad>& loads)
{
    loads.clear();

    /* A tile elements */
    for (int m = 0; m < TILE_M; m++) {
        for (int k = 0; k < TILE_K; k += vector_width) {
            AtomicLoad ld;
            for (int v = 0; v < vector_width && (k+v) < TILE_K; v++)
                ld.byte_addrs.push_back(base_A + addr_A_rm(
                    mt*TILE_M + m, kt*TILE_K + k + v, K));
            loads.push_back(std::move(ld));
        }
    }

    /* B tile elements */
    for (int k = 0; k < TILE_K; k++) {
        for (int n = 0; n < TILE_N; n += vector_width) {
            AtomicLoad ld;
            for (int v = 0; v < vector_width && (n+v) < TILE_N; v++) {
                int64_t addr;
                if (b_layout == B_VNNI)
                    addr = base_B + addr_B_vnni(kt*TILE_K+k, nt*TILE_N+n+v, K, N);
                else if (b_layout == B_BLOCKED)
                    addr = base_B + addr_B_blocked(kt*TILE_K+k, nt*TILE_N+n+v, K, N);
                else
                    addr = base_B + addr_B_rm(kt*TILE_K+k, nt*TILE_N+n+v, N);
                ld.byte_addrs.push_back(addr);
            }
            loads.push_back(std::move(ld));
        }
    }

    /* C tile elements (int32) */
    for (int m = 0; m < TILE_M; m++) {
        for (int n = 0; n < TILE_N; n += vector_width) {
            AtomicLoad ld;
            for (int v = 0; v < vector_width && (n+v) < TILE_N; v++)
                ld.byte_addrs.push_back(base_C + addr_C_rm(
                    mt*TILE_M + m, nt*TILE_N + n + v, N));
            loads.push_back(std::move(ld));
        }
    }
}

/* ═══════════════════════════════════════════════════════════
 *  Cost model: compute mu, delta, lines_per_load from
 *  atomic loads
 * ═══════════════════════════════════════════════════════════ */

struct CostMetrics {
    double mu;              /* avg new cache blocks per step */
    double delta;           /* avg min distance of new blocks */
    double mu_delta;        /* composite */
    double lines_per_load;  /* avg cache lines per atomic load */
    double total_loads;     /* avg atomic loads per step */
    int64_t T;
};

typedef void (*LoadCollector)(int mt, int nt, int kt, int K, int N,
                               BLayout b_layout, int64_t base_A,
                               int64_t base_B, int64_t base_C,
                               std::vector<AtomicLoad>& loads);

/* Generic cost computation: works with any load collector */
CostMetrics compute_cost_generic(
    int M, int N, int K,
    BLayout b_layout,
    int64_t base_A, int64_t base_B, int64_t base_C,
    /* Function that collects atomic loads for one step */
    std::function<void(int,int,int, std::vector<AtomicLoad>&)> collect)
{
    int mt_max = M / TILE_M;
    int nt_max = N / TILE_N;
    int kt_max = K / TILE_K;

    std::set<int64_t> prev_blocks, curr_blocks;
    std::vector<AtomicLoad> loads;

    int64_t T = 0;
    double s_mu = 0, s_delta = 0;
    double s_lines_per_load = 0, s_total_loads = 0;

    for (int mt = 0; mt < mt_max; mt++) {
        for (int nt = 0; nt < nt_max; nt++) {
            for (int kt = 0; kt < kt_max; kt++) {
                curr_blocks.clear();
                collect(mt, nt, kt, loads);

                /* Per-load statistics */
                double step_lines = 0;
                for (auto& ld : loads) {
                    std::set<int64_t> load_blocks;
                    for (int64_t addr : ld.byte_addrs) {
                        int64_t b = blk(addr);
                        load_blocks.insert(b);
                        curr_blocks.insert(b);
                    }
                    step_lines += load_blocks.size();
                }
                s_total_loads += loads.size();
                s_lines_per_load += step_lines;

                /* mu and delta */
                if (T == 0) {
                    s_mu += (double)curr_blocks.size();
                    s_delta += 1.0;
                } else {
                    int n_new = 0; double sum_dist = 0;
                    for (int64_t b : curr_blocks) {
                        if (prev_blocks.count(b)) continue;
                        n_new++;
                        int64_t best = INT64_MAX;
                        for (int64_t bp : prev_blocks)
                            best = std::min(best, std::abs(b - bp));
                        if (best == INT64_MAX) best = 1;
                        sum_dist += (double)best;
                    }
                    s_mu += n_new;
                    if (n_new > 0) s_delta += sum_dist / n_new;
                }

                prev_blocks = curr_blocks;
                T++;
            }
        }
    }

    if (T == 0) return {};
    double dT = (double)T;
    double mu = s_mu / dT, delta = s_delta / dT;
    return {mu, delta, mu*delta,
            s_lines_per_load / s_total_loads,
            s_total_loads / dT, T};
}

/* ═══════════════════════════════════════════════════════════ */

int main(int argc, char** argv) {
    int M = (argc > 1) ? atoi(argv[1]) : 256;
    int N = (argc > 2) ? atoi(argv[2]) : 256;
    int K = (argc > 3) ? atoi(argv[3]) : 256;
    BLOCK_BYTES = (argc > 4) ? atoi(argv[4]) : 64;

    M = (M/TILE_M)*TILE_M; N = (N/TILE_N)*TILE_N; K = (K/TILE_K)*TILE_K;

    printf("AMX GEMM Cost Model with Atomic Access Shapes\n");
    printf("  M=%d  N=%d  K=%d  block=%dB\n", M, N, K, BLOCK_BYTES);
    printf("  Tiles: %dx%dx%d = %d TDPBUSD ops\n\n",
           M/TILE_M, N/TILE_N, K/TILE_K,
           (M/TILE_M)*(N/TILE_N)*(K/TILE_K));

    int64_t base_A = 0;
    int64_t base_B = (int64_t)M * K;
    int64_t base_C = base_B + (int64_t)K * N;

    struct Scenario {
        const char* name;
        const char* target;
        const char* atomic_shape;
        BLayout b_layout;
        bool use_amx_loads;
        int vec_width;  /* for scalar/simd */
    };

    Scenario scenarios[] = {
        /* AMX scenarios */
        {"AMX + B VNNI",          "AMX",      "[16]x[64B] rows",  B_VNNI,     true,  0},
        {"AMX + B blocked",       "AMX",      "[16]x[64B] rows",  B_BLOCKED,  true,  0},
        {"AMX + B row-major",     "AMX",      "[64]x[16B] rows",  B_ROWMAJOR, true,  0},
        /* Scalar CPU scenarios */
        {"Scalar + B row-major",  "Scalar",   "[1] element",      B_ROWMAJOR, false, 1},
        {"Scalar + B blocked",    "Scalar",   "[1] element",      B_BLOCKED,  false, 1},
        {"Scalar + B VNNI",       "Scalar",   "[1] element",      B_VNNI,     false, 1},
        /* AVX-512 scenarios */
        {"AVX512 + B row-major",  "AVX-512",  "[64] int8 vec",    B_ROWMAJOR, false, 16},
        {"AVX512 + B blocked",    "AVX-512",  "[64] int8 vec",    B_BLOCKED,  false, 16},
    };
    int n_scenarios = sizeof(scenarios) / sizeof(scenarios[0]);

    printf("  %-28s %-10s %-18s  %6s  %6s  %8s  %6s  %6s\n",
           "Scenario", "Target", "Atomic shape",
           "mu", "Delta", "mu*Delta", "L/load", "loads");
    printf("  %-28s %-10s %-18s  %6s  %6s  %8s  %6s  %6s\n",
           "----------------------------", "----------", "------------------",
           "------", "------", "--------", "------", "------");

    double best_mu_delta = 1e18;

    for (int s = 0; s < n_scenarios; s++) {
        auto& sc = scenarios[s];

        CostMetrics cm;

        if (sc.use_amx_loads) {
            cm = compute_cost_generic(M, N, K, sc.b_layout, base_A, base_B, base_C,
                [&](int mt, int nt, int kt, std::vector<AtomicLoad>& loads) {
                    collect_amx_atomic_loads(mt, nt, kt, K, N,
                        sc.b_layout, base_A, base_B, base_C, loads);
                });
        } else {
            int vw = sc.vec_width;
            BLayout bl = sc.b_layout;
            cm = compute_cost_generic(M, N, K, sc.b_layout, base_A, base_B, base_C,
                [&, vw, bl](int mt, int nt, int kt, std::vector<AtomicLoad>& loads) {
                    collect_scalar_atomic_loads(mt, nt, kt, K, N,
                        bl, vw, base_A, base_B, base_C, loads);
                });
        }

        if (cm.mu_delta < best_mu_delta) best_mu_delta = cm.mu_delta;

        printf("  %-28s %-10s %-18s  %6.1f  %6.1f  %8.1f  %6.2f  %6.0f",
               sc.name, sc.target, sc.atomic_shape,
               cm.mu, cm.delta, cm.mu_delta,
               cm.lines_per_load, cm.total_loads);
        if (s > 0 && best_mu_delta > 0)
            printf("  (%.2fx)", cm.mu_delta / best_mu_delta);
        printf("\n");
    }

    printf("\n  Key:\n");
    printf("    mu       = avg new cache blocks per TDPBUSD step\n");
    printf("    Delta    = avg min distance of new blocks to previous set\n");
    printf("    mu*Delta = composite memory cost (lower = better)\n");
    printf("    L/load   = avg cache lines touched per atomic load (1.0 = perfect)\n");
    printf("    loads    = avg atomic loads per step\n");
    printf("\n");

    printf("  Observations:\n");
    printf("    - mu*Delta is identical for VNNI and blocked B because both\n");
    printf("      pack each 64x16 tile into 16 contiguous cache lines.\n");
    printf("      The VNNI shuffle permutes within cache lines, invisible to mu*Delta.\n");
    printf("    - L/load reveals utilization: AMX+row-major B issues 64 loads of\n");
    printf("      16 bytes each (0.25 lines/load effective), while AMX+VNNI issues\n");
    printf("      16 loads of 64 bytes each (1.0 lines/load).\n");
    printf("    - Scalar sees the same mu*Delta as AMX for matching B layouts\n");
    printf("      because cache blocks touched are layout-dependent, not target-dependent.\n");
    printf("    - The atomic shape affects L/load (utilization) and loads (count),\n");
    printf("      but not mu*Delta (which depends only on cache-block access pattern).\n");

    return 0;
}