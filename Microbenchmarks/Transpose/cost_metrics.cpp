/*  transpose_metrics.cpp — Compute μ (avg new-block count) and Δ (avg block distance)
 *  for the matrix transpose kernel: out[c*N+r] = in[r*N+c]
 *  under row-major and blocked(SB) layouts, with various schedules.
 *
 *  Compile: g++ -O3 -fopenmp -o transpose_metrics transpose_metrics.cpp
 *  Usage:   ./transpose_metrics <N> <B> [SB_list]
 *           N  = matrix dimension
 *           B  = cache-line size in floats (e.g. 16 for 64-byte lines of fp32)
 *           SB_list = comma-separated block sizes (default: 8,16,32,64,128)
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <omp.h>

/* ── Layout address functions ─────────────────────────────────────── */

// Row-major: element (r,c) in NxN matrix → linear address
static inline long rm_addr(int r, int c, int N) {
    return (long)r * N + c;
}

// Blocked(SB): element (r,c) in NxN matrix → linear address
// Block (br,bc) at linear offset (br*NB+bc)*SB*SB, within-block row-major
static inline long blk_addr(int r, int c, int N, int SB) {
    int NB = N / SB;
    int br = r / SB, bc = c / SB;
    int lr = r % SB, lc = c % SB;
    return (long)(br * NB + bc) * SB * SB + lr * SB + lc;
}

/* ── Metric computation ───────────────────────────────────────────── */

struct Metrics {
    double mu;
    double delta;
};

// Given a schedule (list of (r,c) pairs), layout functions for in/out,
// and block size B, compute μ and Δ.
// We parallelise by splitting the iteration sequence into chunks per thread,
// then combining. Since μ and Δ are defined sequentially (step t depends on t-1),
// we compute serially but the outer sweep over configurations is parallel.

static Metrics compute_metrics(
    const std::vector<std::pair<int,int>>& schedule,
    int N, int B, bool blocked_in, bool blocked_out, int SB)
{
    long T = (long)schedule.size();
    if (T == 0) return {0.0, 0.0};

    // Block address for input and output at iteration t
    auto in_block = [&](int r, int c) -> long {
        long addr = blocked_in ? blk_addr(r, c, N, SB) : rm_addr(r, c, N);
        return addr / B;
    };
    auto out_block = [&](int r, int c) -> long {
        // transpose: out[c*N+r] = in[r*N+c]
        // out element is at position (c,r) in the output matrix
        long addr = blocked_out ? blk_addr(c, r, N, SB) : rm_addr(c, r, N);
        return addr / B;
    };

    double sum_N = 0;       // sum of |N_t|
    double sum_rho_bar = 0; // sum of rho_bar_t

    // Previous step's block set (at most 2 blocks: one in, one out)
    long prev_in = in_block(schedule[0].first, schedule[0].second);
    long prev_out = out_block(schedule[0].first, schedule[0].second);

    // t=1 (index 0): cold misses
    int n1 = (prev_in == prev_out) ? 1 : 2;
    sum_N += n1;
    sum_rho_bar += 1.0; // convention: rho_1 = 1

    for (long t = 1; t < T; t++) {
        int r = schedule[t].first, c = schedule[t].second;
        long bi = in_block(r, c);
        long bo = out_block(r, c);

        // B_{t-1} = {prev_in, prev_out}
        // B_t = {bi, bo}  (may overlap)
        // N_t = B_t \ B_{t-1}

        // Collect new blocks and their distances
        int n_new = 0;
        double rho_sum = 0.0;

        auto check_new = [&](long b) {
            if (b != prev_in && b != prev_out) {
                // New block — compute distance to nearest in B_{t-1}
                long d1 = std::abs(b - prev_in);
                long d2 = std::abs(b - prev_out);
                long d = std::min(d1, d2);
                rho_sum += (double)d;
                n_new++;
            }
        };

        check_new(bi);
        if (bo != bi) check_new(bo);

        sum_N += n_new;
        if (n_new > 0)
            sum_rho_bar += rho_sum / n_new;
        // else rho_bar_t = 0 when N_t is empty

        prev_in = bi;
        prev_out = bo;
    }

    return { sum_N / T, sum_rho_bar / T };
}

/* ── Schedule generators ──────────────────────────────────────────── */

// Naive: row-major iteration
static std::vector<std::pair<int,int>> sched_naive(int N) {
    std::vector<std::pair<int,int>> s;
    s.reserve((long)N * N);
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            s.push_back({r, c});
    return s;
}

// Tiled: TB x TB tiles, row-major tile order, row-major within tile
static std::vector<std::pair<int,int>> sched_tiled(int N, int TB) {
    std::vector<std::pair<int,int>> s;
    s.reserve((long)N * N);
    int NT = (N + TB - 1) / TB;
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++)
            for (int lr = 0; lr < TB && tr * TB + lr < N; lr++)
                for (int lc = 0; lc < TB && tc * TB + lc < N; lc++)
                    s.push_back({tr * TB + lr, tc * TB + lc});
    return s;
}

// Block-aligned: iterate over SB x SB blocks, row-major within each block
static std::vector<std::pair<int,int>> sched_block_aligned(int N, int SB) {
    std::vector<std::pair<int,int>> s;
    s.reserve((long)N * N);
    int NB = N / SB;
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++)
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    s.push_back({br * SB + lr, bc * SB + lc});
    return s;
}

/* ── Main ─────────────────────────────────────────────────────────── */

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <N> <B_floats> [SB_list]\n"
                        "  N        = matrix dimension\n"
                        "  B_floats = block size in elements (e.g. 16 for 64B lines, fp32)\n"
                        "  SB_list  = comma-separated blocking factors (default 8,16,32,64,128)\n",
                argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int B = atoi(argv[2]);

    std::vector<int> SBs;
    if (argc > 3) {
        std::istringstream ss(argv[3]);
        std::string tok;
        while (std::getline(ss, tok, ','))
            SBs.push_back(atoi(tok.c_str()));
    } else {
        SBs = {8, 16, 32, 64, 128};
    }

    // Filter SBs that divide N
    std::vector<int> valid_SBs;
    for (int sb : SBs)
        if (N % sb == 0) valid_SBs.push_back(sb);

    printf("N=%d  B=%d elements (%d bytes for fp32)\n\n", N, B, B * 4);
    printf("%-30s  %10s  %10s  %12s\n", "Configuration", "mu", "delta", "mu*delta");
    printf("%-30s  %10s  %10s  %12s\n",
           "------------------------------", "----------", "----------", "------------");

    // Collect all configurations to evaluate
    struct Config {
        std::string name;
        int sched_type; // 0=naive, 1=tiled(SB), 2=block_aligned(SB)
        int sched_param;
        bool blk_in, blk_out;
        int sb;
    };
    std::vector<Config> configs;

    // Row-major, naive schedule
    configs.push_back({"rm / naive", 0, 0, false, false, 0});

    // Row-major, tiled schedules
    for (int sb : valid_SBs) {
        char name[64];
        snprintf(name, sizeof(name), "rm / tiled(%d)", sb);
        configs.push_back({name, 1, sb, false, false, 0});
    }

    // Row-major data, block-aligned schedule
    for (int sb : valid_SBs) {
        char name[64];
        snprintf(name, sizeof(name), "rm / blk_sched(%d)", sb);
        configs.push_back({name, 2, sb, false, false, sb});
    }

    // Blocked layout, block-aligned schedule
    for (int sb : valid_SBs) {
        char name[64];
        snprintf(name, sizeof(name), "blk(%d) / blk_sched(%d)", sb, sb);
        configs.push_back({name, 2, sb, true, true, sb});
    }

    // Blocked layout, naive schedule (to show layout effect alone)
    for (int sb : valid_SBs) {
        char name[64];
        snprintf(name, sizeof(name), "blk(%d) / naive", sb);
        configs.push_back({name, 0, 0, true, true, sb});
    }

    int nc = (int)configs.size();
    std::vector<Metrics> results(nc);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nc; i++) {
        auto& cfg = configs[i];

        std::vector<std::pair<int,int>> sched;
        switch (cfg.sched_type) {
            case 0: sched = sched_naive(N); break;
            case 1: sched = sched_tiled(N, cfg.sched_param); break;
            case 2: sched = sched_block_aligned(N, cfg.sched_param); break;
        }

        results[i] = compute_metrics(sched, N, B, cfg.blk_in, cfg.blk_out,
                                     cfg.sb > 0 ? cfg.sb : 1);
    }

    for (int i = 0; i < nc; i++) {
        printf("%-30s  %10.4f  %10.4f  %12.4f\n",
               configs[i].name.c_str(),
               results[i].mu, results[i].delta,
               results[i].mu * results[i].delta);
    }

    return 0;
}