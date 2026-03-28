/*
 * cost_metrics_conj.cpp -- Cost metrics for complex conjugate kernel
 *
 * Layouts: AoS, SoA, AoSoA-K
 * Modes:   in-place, out-of-place
 *
 * Metrics per config:
 *   mu         = avg new blocks per step
 *   W          = avg total distinct blocks per step (working set)
 *   delta_min  = avg block distance (nearest in prev set)
 *   delta_max  = avg block distance (farthest in prev set)
 *   sigma      = avg element-address stride
 *   delta_numa = avg NUMA-weighted min distance
 *   W_delta    = W * delta_min  (composite)
 *   span       = avg (max_block - min_block) per step
 *
 * Compile: g++ -O3 -std=c++17 -o cost_metrics_conj cost_metrics_conj.cpp
 * Run:     ./cost_metrics_conj [N] [block_bytes] [P_NUMA]
 */
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <climits>
#include <cstring>
#include <set>

static int BLOCK_BYTES = 64;
static int P_NUMA = 4;
static double ALPHA = 0.125;
static double GAMMA = 8.0;
static int BETA = 4;

/* ---- Block address from byte offset ---- */
inline int64_t blk(int64_t byte_offset) {
    return byte_offset / BLOCK_BYTES;
}

/* ---- NUMA domain from byte offset (first-touch static schedule) ---- */
inline int numa_domain(int64_t byte_offset, int64_t total_bytes) {
    int64_t chunk = (total_bytes + P_NUMA - 1) / P_NUMA;
    int d = (int)(byte_offset / chunk);
    return (d < P_NUMA) ? d : P_NUMA - 1;
}

/* ---- Layout: compute byte offsets of refs for element i ---- */
struct LayoutConfig {
    const char* name;
    const char* csv_name;
    int K;  /* AoSoA group size, 0=AoS, -1=SoA */
    bool out_of_place;
    const char* mode;
};

static void get_refs(const LayoutConfig& L, int64_t N, int64_t i,
                     std::vector<int64_t>& byte_offsets) {
    byte_offsets.clear();

    if (L.K == 0) {
        /* AoS: element i at [i*16, i*16+16) */
        int64_t base_in = i * 16;
        byte_offsets.push_back(base_in);  /* real+imag together */
        if (L.out_of_place) {
            int64_t base_out = N * 16 + i * 16;
            byte_offsets.push_back(base_out);
        }
    } else if (L.K < 0) {
        /* SoA: real[i] at i*8, imag[i] at N*8 + i*8 */
        int64_t arr_size = N * 8;
        if (L.out_of_place) {
            byte_offsets.push_back(i * 8);                  /* real_in */
            byte_offsets.push_back(arr_size + i * 8);       /* imag_in */
            byte_offsets.push_back(2*arr_size + i * 8);     /* real_out */
            byte_offsets.push_back(3*arr_size + i * 8);     /* imag_out */
        } else {
            byte_offsets.push_back(arr_size + i * 8);       /* imag (RMW) */
            /* real not touched for conjugate */
        }
    } else {
        /* AoSoA-K: group g = i/K, pos p = i%K
         * real[i] at g*(2*K*8) + p*8
         * imag[i] at g*(2*K*8) + K*8 + p*8 */
        int K = L.K;
        int64_t g = i / K, p = i % K;
        int64_t group_bytes = 2 * K * 8;
        if (L.out_of_place) {
            int64_t in_base = g * group_bytes;
            int64_t out_base = N / K * group_bytes + g * group_bytes;
            byte_offsets.push_back(in_base + p * 8);            /* real_in */
            byte_offsets.push_back(in_base + K * 8 + p * 8);   /* imag_in */
            byte_offsets.push_back(out_base + p * 8);           /* real_out */
            byte_offsets.push_back(out_base + K * 8 + p * 8);  /* imag_out */
        } else {
            int64_t base = g * group_bytes;
            byte_offsets.push_back(base + K * 8 + p * 8);      /* imag (RMW) */
        }
    }
}

static int64_t total_bytes(const LayoutConfig& L, int64_t N) {
    if (L.K == 0) return (L.out_of_place ? 2 : 1) * N * 16;
    if (L.K < 0)  return (L.out_of_place ? 4 : 1) * N * 8;
    return (L.out_of_place ? 2 : 1) * N * 2 * L.K * 8 / L.K; /* = 2*N*8 per side */
}

/* ---- Metrics struct ---- */
struct Metrics {
    double mu;          /* avg new blocks per step */
    double W;           /* avg total distinct blocks per step */
    double delta_min;   /* avg distance to nearest in prev set */
    double delta_max;   /* avg distance to farthest in prev set */
    double sigma;       /* avg element-address stride (bytes) */
    double delta_numa;  /* avg NUMA-weighted min distance */
    double W_delta;     /* W * delta_min */
    double span;        /* avg (max_block - min_block) per step */
    int64_t T;
    int n_refs;         /* refs per step */
};

/* ---- Compute metrics for one layout config ---- */
Metrics compute_metrics(const LayoutConfig& L, int64_t N, int vec_width) {
    int64_t tot_bytes = total_bytes(L, N);
    int64_t max_steps = std::min(N, (int64_t)200000);  /* cap for speed */

    std::vector<int64_t> curr_offsets, prev_offsets;
    std::set<int64_t> prev_blocks, curr_blocks_set;

    int64_t T = 0;
    double s_mu = 0, s_W = 0;
    double s_dmin = 0, s_dmax = 0;
    double s_sigma = 0, s_dnuma = 0;
    double s_span = 0;

    /* Per-reference tracking for sigma */
    std::vector<int64_t> prev_ref_blocks;
    bool have_prev = false;

    for (int64_t i = 0; i < max_steps; i += vec_width) {
        /* Collect all refs for this W-wide step */
        curr_blocks_set.clear();
        curr_offsets.clear();

        for (int w = 0; w < vec_width && (i+w) < N; w++) {
            std::vector<int64_t> offsets;
            get_refs(L, N, i + w, offsets);
            for (int64_t off : offsets) {
                curr_offsets.push_back(off);
                curr_blocks_set.insert(blk(off));
            }
        }

        std::vector<int64_t> curr_bvec(curr_blocks_set.begin(), curr_blocks_set.end());
        int W_t = (int)curr_bvec.size();
        s_W += W_t;

        /* Span */
        int64_t bmin = curr_bvec.front();
        int64_t bmax = curr_bvec.back();
        s_span += (double)(bmax - bmin);

        if (T == 0) {
            s_mu += W_t;
            s_dmin += 1.0;
            s_dmax += 1.0;
            s_dnuma += 1.0;
            s_sigma += W_t;
        } else {
            /* New blocks */
            int n_new = 0;
            double sum_dmin = 0, sum_dmax = 0, sum_dnuma = 0;

            for (int64_t b : curr_bvec) {
                if (prev_blocks.count(b)) continue;
                n_new++;

                /* Min and max distance to prev set */
                int64_t best_min = INT64_MAX, best_max = 0;
                for (int64_t bp : prev_blocks) {
                    int64_t d = std::abs(b - bp);
                    if (d < best_min) best_min = d;
                    if (d > best_max) best_max = d;
                }
                if (best_min == INT64_MAX) best_min = 1;

                sum_dmin += (double)best_min;
                sum_dmax += (double)best_max;

                /* NUMA-weighted min distance */
                int nu_b = numa_domain(b * BLOCK_BYTES, tot_bytes);
                int nu_home = numa_domain(i * 8, tot_bytes);  /* home = where we are */
                double w_numa;
                if (nu_b != nu_home) w_numa = GAMMA;
                else if (best_min < BETA) w_numa = ALPHA;
                else w_numa = 1.0;
                sum_dnuma += w_numa * (double)best_min;
            }

            s_mu += n_new;
            if (n_new > 0) {
                s_dmin += sum_dmin / n_new;
                s_dmax += sum_dmax / n_new;
                s_dnuma += sum_dnuma / n_new;
            }

            /* Per-reference sigma (element stride) */
            std::vector<int64_t> offsets0;
            get_refs(L, N, i, offsets0);
            if (have_prev && offsets0.size() == prev_ref_blocks.size()) {
                double stride_sum = 0;
                for (size_t r = 0; r < offsets0.size(); r++) {
                    stride_sum += std::abs(blk(offsets0[r]) - prev_ref_blocks[r]);
                }
                s_sigma += stride_sum;
            }
        }

        /* Save state for next step */
        prev_blocks.clear();
        prev_blocks.insert(curr_bvec.begin(), curr_bvec.end());

        {
            std::vector<int64_t> offsets0;
            get_refs(L, N, i, offsets0);
            prev_ref_blocks.clear();
            for (int64_t off : offsets0) prev_ref_blocks.push_back(blk(off));
            have_prev = true;
        }

        T++;
    }

    if (T == 0) return {};
    double dT = (double)T;
    double mu = s_mu / dT;
    double Wavg = s_W / dT;
    double dmin = s_dmin / dT;
    double dmax = s_dmax / dT;
    double sig = s_sigma / dT;
    double dnuma = s_dnuma / dT;

    /* Count refs per step */
    std::vector<int64_t> tmp;
    get_refs(L, N, 0, tmp);
    int n_refs = (int)tmp.size();

    return {mu, Wavg, dmin, dmax, sig, dnuma, Wavg * dmin, s_span / dT, T, n_refs};
}

int main(int argc, char** argv) {
    int64_t N = (argc > 1) ? atol(argv[1]) : 1 << 24;  /* 16M complex numbers */
    BLOCK_BYTES = (argc > 2) ? atoi(argv[2]) : 64;
    P_NUMA      = (argc > 3) ? atoi(argv[3]) : 4;
    ALPHA       = (argc > 4) ? atof(argv[4]) : 0.125;
    GAMMA       = (argc > 5) ? atof(argv[5]) : 8.0;
    BETA        = (argc > 6) ? atoi(argv[6]) : 4;

    fprintf(stderr, "Conjugate cost metrics: N=%ld  block=%d  P=%d  alpha=%.4f  gamma=%.1f  beta=%d\n\n",
            (long)N, BLOCK_BYTES, P_NUMA, ALPHA, GAMMA, BETA);

    LayoutConfig configs[] = {
        {"AoS",       "AoS",       0,  false, "in-place"},
        {"AoS",       "AoS",       0,  true,  "out-of-place"},
        {"SoA",       "SoA",      -1,  false, "in-place"},
        {"SoA",       "SoA",      -1,  true,  "out-of-place"},
        {"AoSoA-8",   "AoSoA-8",   8,  false, "in-place"},
        {"AoSoA-8",   "AoSoA-8",   8,  true,  "out-of-place"},
        {"AoSoA-16",  "AoSoA-16", 16,  false, "in-place"},
        {"AoSoA-16",  "AoSoA-16", 16,  true,  "out-of-place"},
    };
    int n_configs = sizeof(configs) / sizeof(configs[0]);

    /* CSV */
    FILE* csv = fopen("metrics_conj.csv", "w");
    fprintf(csv, "layout,mode,n_refs,mu,W,delta_min,delta_max,sigma,delta_numa,W_delta,span,block_bytes,N\n");

    /* Table header */
    printf("  %-10s %-12s %4s %7s %5s %9s %9s %9s %9s %9s %9s\n",
           "Layout", "Mode", "refs", "mu", "W", "D_min", "D_max", "sigma", "D_numa", "W*D_min", "span");
    printf("  %-10s %-12s %4s %7s %5s %9s %9s %9s %9s %9s %9s\n",
           "----------", "------------", "----", "-------", "-----",
           "---------", "---------", "---------", "---------", "---------", "---------");

    for (int c = 0; c < n_configs; c++) {
        const auto& L = configs[c];
        Metrics m = compute_metrics(L, N, 1);

        auto fmt = [](double v, char* buf) {
            if (v < 0.01)       snprintf(buf, 16, "%9.6f", v);
            else if (v < 100)   snprintf(buf, 16, "%9.4f", v);
            else if (v < 100000) snprintf(buf, 16, "%9.1f", v);
            else                 snprintf(buf, 16, "%9.0f", v);
        };

        char b[9][16];
        fmt(m.mu, b[0]); fmt(m.W, b[1]); fmt(m.delta_min, b[2]);
        fmt(m.delta_max, b[3]); fmt(m.sigma, b[4]); fmt(m.delta_numa, b[5]);
        fmt(m.W_delta, b[6]); fmt(m.span, b[7]);

        printf("  %-10s %-12s %4d %s %s %s %s %s %s %s %s\n",
               L.name, L.mode, m.n_refs,
               b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);

        fprintf(csv, "%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%ld\n",
                L.csv_name, L.mode, m.n_refs,
                m.mu, m.W, m.delta_min, m.delta_max, m.sigma, m.delta_numa,
                m.W_delta, m.span, BLOCK_BYTES, (long)N);
    }

    fclose(csv);
    fprintf(stderr, "\nWritten: metrics_conj.csv\n");
    return 0;
}