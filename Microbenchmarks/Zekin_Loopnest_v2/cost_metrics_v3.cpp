/*
 * cost_metrics.cpp -- Three-tier cost metrics for z_v_grad_w stencil
 *
 * Metrics (Section 3):
 *   μ(π,Φ)      = avg new-block count per step           (Eq. 4)
 *   Δ(π,Φ)      = avg block distance, UMA                (Eq. 7)
 *   σ(π,Φ,ν)    = avg weighted traffic, NUMA              (Eq. 10)
 *
 * Weight functions:
 *   ω_U(b)  = α              if ρ_t(b) < β                (UMA)
 *           = 1              otherwise
 *   ω(b)    = γ              if ν(b) ≠ ν₀                 (NUMA: remote priority)
 *           = α              if ν(b) = ν₀ and ρ_t(b) < β
 *           = 1              otherwise
 *
 * ν : block → NUMA domain     (placement function, from first-touch)
 * ν₀ = 0                      (home domain of evaluated iteration slice)
 *
 * Compile:  g++ -O3 -std=c++17 -fopenmp -o cost_metrics cost_metrics.cpp
 * Run:      ./cost_metrics [N] [nlev] [β] [α] [γ] [P_NUMA]
 * Output:   results_full.csv
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
#include <cstring>
#include <unistd.h>

/* ---- Parameters (from argv) ---- */
static int    BETA       = 8;
static double ALPHA      = 0.16;
static double GAMMA      = 2.0;
static int    P_NUMA     = 4;
static int    BLOCK_BYTES = 64;
#pragma omp threadprivate(BLOCK_BYTES)

/* ---- Layout index functions ---- */
inline int IC(int V, int je, int jk, int N, int nlev) {
    return (V <= 2) ? je + jk * N : jk + je * nlev;
}
inline int IN(int V, int je, int n, int N) {
    return (V == 1 || V == 3) ? je + n * N : n + je * 2;
}

/* ---- Arrays ---- */
enum ArrID {
    A_OUT=0, A_VN_IE, A_W, A_Z_VT_IE, A_Z_W_V,
    A_INV_DUAL, A_INV_PRIMAL, A_TANGENT,
    A_CELL_IDX, A_VERT_IDX, NUM_ARR
};
static const int ebytes[] = { 8,8,8,8,8, 8,8,8, 4,4 };

inline int64_t blk(int elem_idx, int arr) {
    return (int64_t)elem_idx * ebytes[arr] / BLOCK_BYTES;
}

/* ---- Schedule / loop order ---- */
enum Schedule  { SCHED_OMP_FOR=0, SCHED_OMP_COLLAPSE2=1 };
static const char* sched_name[] = {"omp_for", "omp_collapse2"};

enum LoopOrder { KLON_FIRST=0, KLEV_FIRST=1 };
static const char* loop_name[] = {"klon_first", "klev_first"};

inline LoopOrder iteration_order(Schedule sched, int V) {
    if (sched == SCHED_OMP_COLLAPSE2) return KLON_FIRST;
    return (V <= 2) ? KLON_FIRST : KLEV_FIRST;
}

/* ============================================================
 *  Placement function  ν : block_addr → NUMA domain
 *
 *  External input to the cost model.  Encodes first-touch under
 *  schedule(static) without introducing threads into the model.
 * ============================================================ */

static Schedule G_SCHED;
#pragma omp threadprivate(G_SCHED)

inline int numa_2d(Schedule sched, int V, int elem, int N, int nlev) {
    int64_t total = (int64_t)N * nlev;
    if (sched == SCHED_OMP_COLLAPSE2) {
        int64_t lin;
        if (V <= 2) lin = elem;
        else { int jk = elem % nlev, je = elem / nlev;
               lin = (int64_t)jk * N + je; }
        int64_t chunk = (total + P_NUMA - 1) / P_NUMA;
        int d = (int)(lin / chunk);
        return (d < P_NUMA) ? d : P_NUMA - 1;
    }
    if (V <= 2) {
        int jk = elem / N;
        int chunk = (nlev + P_NUMA - 1) / P_NUMA;
        int d = jk / chunk;
        return (d < P_NUMA) ? d : P_NUMA - 1;
    } else {
        int je = elem / nlev;
        int chunk = (N + P_NUMA - 1) / P_NUMA;
        int d = je / chunk;
        return (d < P_NUMA) ? d : P_NUMA - 1;
    }
}

inline int numa_1d(int elem, int N) {
    int chunk = (N + P_NUMA - 1) / P_NUMA;
    if (chunk == 0) return 0;
    int d = elem / chunk;
    return (d < P_NUMA) ? d : P_NUMA - 1;
}

inline int numa_idx(int V, int elem, int N) {
    int je = (V == 1 || V == 3) ? elem % N : elem / 2;
    return numa_1d(je, N);
}

inline int nu_elem(int arr, int V, int elem, int N, int nlev) {
    switch (arr) {
    case A_OUT: case A_VN_IE: case A_W: case A_Z_VT_IE: case A_Z_W_V:
        return numa_2d(G_SCHED, V, elem, N, nlev);
    case A_INV_DUAL: case A_INV_PRIMAL: case A_TANGENT:
        return numa_1d(elem, N);
    case A_CELL_IDX: case A_VERT_IDX:
        return numa_idx(V, elem, N);
    default: return 0;
    }
}

inline int nu_block(int arr, int V, int64_t baddr, int N, int nlev) {
    int elem = (int)(baddr * BLOCK_BYTES / ebytes[arr]);
    return nu_elem(arr, V, elem, N, nlev);
}

/* ============================================================
 *  Weight functions
 * ============================================================ */

static int NU_HOME = 0;
#pragma omp threadprivate(NU_HOME)

/* ω_U(b): UMA — prefetcher discount only (Eq. 6) */
inline double w_uma(int64_t raw_dist) {
    return (raw_dist < BETA) ? ALPHA : 1.0;
}

/* ω(b): NUMA — remote priority (Eq. 9) */
inline double w_numa(int arr, int V, int64_t baddr,
                     int64_t raw_dist, int N, int nlev) {
    if (nu_block(arr, V, baddr, N, nlev) != NU_HOME) return GAMMA;
    return (raw_dist < BETA) ? ALPHA : 1.0;
}

/* ---- References ---- */
static constexpr int NREFS = 14;
struct Ref { int arr; int64_t block; };

inline void collect_refs(int V, int je, int jk, int N, int nlev,
                         const int* cidx, const int* vidx, Ref refs[NREFS]) {
    int c2d = IC(V, je, jk, N, nlev);
    int ci0 = cidx[IN(V,je,0,N)], ci1 = cidx[IN(V,je,1,N)];
    int vi0 = vidx[IN(V,je,0,N)], vi1 = vidx[IN(V,je,1,N)];
    refs[ 0] = { A_OUT,        blk(c2d, A_OUT) };
    refs[ 1] = { A_VN_IE,      blk(c2d, A_VN_IE) };
    refs[ 2] = { A_INV_DUAL,   blk(je,  A_INV_DUAL) };
    refs[ 3] = { A_W,          blk(IC(V,ci0,jk,N,nlev), A_W) };
    refs[ 4] = { A_W,          blk(IC(V,ci1,jk,N,nlev), A_W) };
    refs[ 5] = { A_Z_VT_IE,    blk(c2d, A_Z_VT_IE) };
    refs[ 6] = { A_INV_PRIMAL, blk(je,  A_INV_PRIMAL) };
    refs[ 7] = { A_TANGENT,    blk(je,  A_TANGENT) };
    refs[ 8] = { A_Z_W_V,      blk(IC(V,vi0,jk,N,nlev), A_Z_W_V) };
    refs[ 9] = { A_Z_W_V,      blk(IC(V,vi1,jk,N,nlev), A_Z_W_V) };
    refs[10] = { A_CELL_IDX,   blk(IN(V,je,0,N), A_CELL_IDX) };
    refs[11] = { A_CELL_IDX,   blk(IN(V,je,1,N), A_CELL_IDX) };
    refs[12] = { A_VERT_IDX,   blk(IN(V,je,0,N), A_VERT_IDX) };
    refs[13] = { A_VERT_IDX,   blk(IN(V,je,1,N), A_VERT_IDX) };
}

/* ---- Per-array block set B_t ---- */
struct BlockSet {
    std::vector<int64_t> a[NUM_ARR];
    void clear() { for (int i=0;i<NUM_ARR;i++) a[i].clear(); }
    void add(int arr, int64_t b) { a[arr].push_back(b); }
    void finalize() {
        for (int i=0;i<NUM_ARR;i++){
            auto& v=a[i]; std::sort(v.begin(),v.end());
            v.erase(std::unique(v.begin(),v.end()),v.end());
        }
    }
    int total() const {
        int n=0; for (int i=0;i<NUM_ARR;i++) n+=(int)a[i].size(); return n;
    }
    bool has(int arr, int64_t b) const {
        return std::binary_search(a[arr].begin(),a[arr].end(),b);
    }
};

/* ============================================================
 *  Core: compute μ, Δ, σ for one NUMA domain's iteration slice
 *
 *  Two modes:
 *    omp_for:      outer ∈ [range_lo, range_hi), inner full
 *    collapse(2):  linearized lin ∈ [range_lo, range_hi),
 *                  jk = lin / N,  je = lin % N
 *                  may cut mid-row at boundaries
 *
 *  NU_HOME must be set before calling.
 *
 *  Δ: (1/T) Σ_t [ (1/|N_t|) Σ_{b∈N_t} ω_U(b)·ρ ]   (Eq. 7)
 *  σ: (1/T) Σ_t [            Σ_{b∈N_t} ω(b)·ρ   ]   (Eq. 10)
 * ============================================================ */

struct Metrics { double mu, delta, sigma, delta_max, mu_delta_max; int64_t T; };

/* Process one W-wide step — shared by both modes.
 * For KLON_FIRST: je varies [je0, je0+W), jk fixed.
 * For KLEV_FIRST: jk varies [jk0, jk0+W), je fixed. */
static inline void process_step(
    int V, int W, int N, int nlev,
    LoopOrder loop, int je0, int jk0,
    const int* cidx, const int* vidx,
    BlockSet& prev, BlockSet& curr,
    int64_t& T, double& s_mu, double& s_delta, double& s_sigma,
    double& s_delta_max, double& s_mu_delta_max)
{
    Ref refs[NREFS];
    curr.clear();
    for (int w = 0; w < W; w++) {
        int je = (loop == KLON_FIRST) ? je0 + w : je0;
        int jk = (loop == KLON_FIRST) ? jk0     : jk0 + w;
        collect_refs(V, je, jk, N, nlev, cidx, vidx, refs);
        for (int r = 0; r < NREFS; r++) curr.add(refs[r].arr, refs[r].block);
    }
    curr.finalize();

    if (T == 0) {
        int nb = curr.total();
        s_mu          += nb;
        s_delta       += 1.0;
        s_sigma       += nb;
        s_delta_max   += 1.0;
        s_mu_delta_max += nb * 1.0;
    } else {
        int    nc = 0;
        double d_uma = 0, d_numa = 0;
        double max_uma = 0;                /* max ω_U·ρ this step */
        for (int a = 0; a < NUM_ARR; a++) {
          for (int64_t b : curr.a[a]) {
            if (prev.has(a, b)) continue;
            nc++;
            int64_t best = INT64_MAX;
            for (int64_t bp : prev.a[a]) {
                int64_t d = std::abs(b - bp);
                if (d < best) best = d;
            }
            double cost_uma, cost_numa;
            if (best == INT64_MAX) {
                cost_uma  = 1.0;
                cost_numa = w_numa(a, V, b, INT64_MAX, N, nlev);
            } else {
                cost_uma  = w_uma(best)                    * (double)best;
                cost_numa = w_numa(a, V, b, best, N, nlev) * (double)best;
            }
            d_uma  += cost_uma;
            d_numa += cost_numa;
            if (cost_uma > max_uma) max_uma = cost_uma;
          }
        }
        s_mu    += nc;
        s_delta += (nc > 0) ? d_uma / nc : 0.0;
        s_sigma += d_numa;
        s_delta_max    += max_uma;              /* Δ_max: max distance this step */
        s_mu_delta_max += (double)nc * max_uma; /* μ·Δ_max: blocks × max distance */
    }
    T++;
    std::swap(prev, curr);
}

Metrics compute_slice(int V, int W, int N, int nlev,
                      Schedule sched,
                      const int* cidx, const int* vidx,
                      int64_t range_lo, int64_t range_hi) {
    G_SCHED = sched;
    LoopOrder loop = iteration_order(sched, V);
    BlockSet prev, curr;

    int64_t T = 0;
    double s_mu = 0, s_delta = 0, s_sigma = 0;
    double s_dmax = 0, s_mudmax = 0;

    if (sched == SCHED_OMP_COLLAPSE2) {
        int64_t lin = range_lo;
        while (lin < range_hi) {
            int jk = (int)(lin / N);
            int je_start = (int)(lin % N);
            int64_t row_end = std::min((int64_t)(jk + 1) * N, range_hi);
            int je_end = (int)(row_end - (int64_t)jk * N);

            for (int je0 = je_start; je0 + W <= je_end; je0 += W) {
                process_step(V, W, N, nlev, KLON_FIRST, je0, jk, cidx, vidx,
                             prev, curr, T, s_mu, s_delta, s_sigma,
                             s_dmax, s_mudmax);
            }
            lin = row_end;
        }
    } else {
        int inner_n = (loop == KLON_FIRST) ? N : nlev;
        for (int64_t outer = range_lo; outer < range_hi; outer++) {
            for (int inner0 = 0; inner0 + W <= inner_n; inner0 += W) {
                int je = (loop == KLON_FIRST) ? inner0    : (int)outer;
                int jk = (loop == KLON_FIRST) ? (int)outer : inner0;
                process_step(V, W, N, nlev, loop, je, jk, cidx, vidx,
                             prev, curr, T, s_mu, s_delta, s_sigma,
                             s_dmax, s_mudmax);
            }
        }
    }

    if (T == 0) return {0, 0, 0, 0, 0, 0};
    return { s_mu/T, s_delta/T, s_sigma/T, s_dmax/T, s_mudmax/T, T };
}

/* ============================================================
 *  Wrapper: evaluate all P_NUMA domains and average
 *
 *  omp_for:      each domain d owns outer ∈ [d·chunk, (d+1)·chunk)
 *  collapse(2):  each domain d owns lin ∈ [d·chunk, (d+1)·chunk)
 *                where total = outer_n × inner_n
 *
 *  σ differs per domain (ν₀ = d).  We average all three.
 * ============================================================ */

Metrics compute_metrics(int V, int W, int N, int nlev,
                        Schedule sched,
                        const int* cidx, const int* vidx) {
    LoopOrder loop = iteration_order(sched, V);
    int outer_n = (loop == KLON_FIRST) ? nlev : N;
    int inner_n = (loop == KLON_FIRST) ? N    : nlev;

    double tot_mu = 0, tot_delta = 0, tot_sigma = 0;
    double tot_dmax = 0, tot_mudmax = 0;
    int64_t tot_T = 0;
    int n_domains = 0;

    if (sched == SCHED_OMP_COLLAPSE2) {
        int64_t total = (int64_t)outer_n * inner_n;
        int64_t chunk = (total + P_NUMA - 1) / P_NUMA;

        for (int d = 0; d < P_NUMA; d++) {
            int64_t lo = d * chunk;
            int64_t hi = std::min((d + 1) * chunk, total);
            if (lo >= hi) continue;

            NU_HOME = d;
            Metrics m = compute_slice(V, W, N, nlev, sched, cidx, vidx, lo, hi);
            if (m.T == 0) continue;

            tot_mu += m.mu; tot_delta += m.delta; tot_sigma += m.sigma;
            tot_dmax += m.delta_max; tot_mudmax += m.mu_delta_max;
            tot_T += m.T; n_domains++;
        }
    } else {
        int chunk = (outer_n + P_NUMA - 1) / P_NUMA;

        for (int d = 0; d < P_NUMA; d++) {
            int lo = d * chunk;
            int hi = std::min((d + 1) * chunk, outer_n);
            if (lo >= hi) continue;

            NU_HOME = d;
            Metrics m = compute_slice(V, W, N, nlev, sched, cidx, vidx, lo, hi);
            if (m.T == 0) continue;

            tot_mu += m.mu; tot_delta += m.delta; tot_sigma += m.sigma;
            tot_dmax += m.delta_max; tot_mudmax += m.mu_delta_max;
            tot_T += m.T; n_domains++;
        }
    }

    if (n_domains == 0) return {0, 0, 0, 0, 0, 0};
    return { tot_mu / n_domains, tot_delta / n_domains,
             tot_sigma / n_domains,
             tot_dmax / n_domains, tot_mudmax / n_domains,
             tot_T };
}

/* ---- Index generation ---- */
enum CellDist { UNIFORM=0, NORMAL1=1, NORMAL4=2, SEQUENTIAL=3 };
static const char* dist_name[] = {"uniform","normal1","normal4","sequential"};

static void gen_cell_idx(int* dst, int V, int N,
                         CellDist dist, std::mt19937& rng) {
    std::vector<int> L(N*2);
    switch (dist) {
    case UNIFORM: {
        std::uniform_int_distribution<int> u(0,N-1);
        for(int i=0;i<N;i++){L[2*i]=u(rng);L[2*i+1]=u(rng);} break;}
    case NORMAL1: {
        std::normal_distribution<double> nd(0,1);
        for(int i=0;i<N;i++){
            int v0=i+1+(int)std::round(nd(rng)),v1=i-1+(int)std::round(nd(rng));
            L[2*i]=((v0%N)+N)%N;L[2*i+1]=((v1%N)+N)%N;} break;}
    case NORMAL4: {
        std::normal_distribution<double> nd(0,2);
        for(int i=0;i<N;i++){
            int v0=i+1+(int)std::round(nd(rng)),v1=i-1+(int)std::round(nd(rng));
            L[2*i]=((v0%N)+N)%N;L[2*i+1]=((v1%N)+N)%N;} break;}
    case SEQUENTIAL:
        for(int i=0;i<N;i++){L[2*i]=(i+1)%N;L[2*i+1]=(i+1)%N;} break;
    }
    for(int je=0;je<N;je++){
        dst[IN(V,je,0,N)]=L[2*je]; dst[IN(V,je,1,N)]=L[2*je+1];
    }
}

static void gen_vert_idx(int* dst, int V, int N, std::mt19937& rng) {
    std::vector<int> p(N);
    std::iota(p.begin(),p.end(),0); std::shuffle(p.begin(),p.end(),rng);
    for(int je=0;je<N;je++) dst[IN(V,je,0,N)]=p[je];
    std::iota(p.begin(),p.end(),0); std::shuffle(p.begin(),p.end(),rng);
    for(int je=0;je<N;je++) dst[IN(V,je,1,N)]=p[je];
}

/* ---- Targets ---- */
struct Target { const char* name; int block_bytes; int vec_width; };
static const Target targets[] = {
    {"CPU_scalar",  64,  1}, {"CPU_NEON",   64,  2}, {"CPU_AVX512", 64,  8},
    {"GPU_scalar", 128,  1}, {"GPU_halfw", 128, 16}, {"GPU_warp32",128, 32},
};
static constexpr int N_TGT = sizeof(targets)/sizeof(targets[0]);

/* ---- Main ---- */
int main(int argc, char** argv) {
    int N    = (argc>1)?atoi(argv[1]):81920;
    int nlev = (argc>2)?atoi(argv[2]):96;
    BETA     = (argc>3)?atoi(argv[3]):(int)(sysconf(_SC_PAGESIZE)/64);
    ALPHA    = (argc>4)?atof(argv[4]):0.16;
    GAMMA    = (argc>5)?atof(argv[5]):2.0;
    P_NUMA   = (argc>6)?atoi(argv[6]):4;

    fprintf(stderr,"Cost metrics: N=%d nlev=%d β=%d α=%.3f γ=%.3f P=%d\n\n",
            N,nlev,BETA,ALPHA,GAMMA,P_NUMA);

    FILE* csv = fopen("results_full.csv","w");
    fprintf(csv,"target,block_bytes,vec_width,variant,schedule,loop_order,"
                "dist,T,mu,delta,sigma,delta_max,mu_delta_max,"
                "beta,alpha,gamma,P_NUMA\n");
    fflush(csv);

    /* Flatten: 4 variants × 2 schedules × 4 dists × N_TGT targets */
    int n_jobs = 4 * 2 * 4 * N_TGT;

    #pragma omp parallel for schedule(dynamic,1)
    for (int job = 0; job < n_jobs; job++) {
        int rem = job;
        int ti  = rem % N_TGT;  rem /= N_TGT;
        int d   = rem % 4;      rem /= 4;
        int si  = rem % 2;      rem /= 2;
        int V   = rem + 1;

        Schedule sched = (si == 0) ? SCHED_OMP_FOR : SCHED_OMP_COLLAPSE2;
        LoopOrder loop = iteration_order(sched, V);

        int local_block_bytes = targets[ti].block_bytes;
        int W = targets[ti].vec_width;
        int inner = (loop == KLON_FIRST) ? N : nlev;
        if (W > inner) continue;

        std::mt19937 rng(42 + d);
        std::vector<int> cidx(N*2), vidx(N*2);
        gen_cell_idx(cidx.data(), V, N, (CellDist)d, rng);
        gen_vert_idx(vidx.data(), V, N, rng);

        BLOCK_BYTES = local_block_bytes;

        auto m = compute_metrics(V, W, N, nlev, sched,
                                  cidx.data(), vidx.data());

        flockfile(stderr);
        fprintf(stderr,"  %-12s V%d %-14s %-12s %-10s  "
                "μ=%.2f Δ=%.1f σ=%.1f Δ_max=%.1f μΔ_max=%.1f\n",
                targets[ti].name, V, sched_name[si],
                loop_name[loop], dist_name[d],
                m.mu, m.delta, m.sigma, m.delta_max, m.mu_delta_max);
        funlockfile(stderr);

        flockfile(csv);
        fprintf(csv,"%s,%d,%d,V%d,%s,%s,%s,%ld,"
                "%.6f,%.6f,%.6f,%.6f,%.6f,"
                "%d,%.4f,%.4f,%d\n",
                targets[ti].name, targets[ti].block_bytes, W,
                V, sched_name[si], loop_name[loop], dist_name[d],
                (long)m.T,
                m.mu, m.delta, m.sigma, m.delta_max, m.mu_delta_max,
                BETA, ALPHA, GAMMA, P_NUMA);
        fflush(csv);
        funlockfile(csv);
    }

    fclose(csv);
    fprintf(stderr,"\nWritten: results_full.csv\n");
    return 0;
}