// bench_goto_vec.cpp
// ═══════════════════════════════════════════════════════════════════════════
// GOTO GEMM — no-pack microkernel + CANNON over NUMA domains
// ═══════════════════════════════════════════════════════════════════════════
//
// Key design decisions (see DECISION markers throughout):
//
//   1. NO PACKING.  The microkernel reads A and B directly from domain-local
//      arrays with their natural strides.  Packing was copying ~280 MB per
//      GEMM call for zero vectorization benefit (A is scalar-broadcast,
//      B rows are already contiguous in row-major layout).
//
//   2. CANNON over 2×2 NUMA grid.  C, A, B are block-distributed across 4
//      NUMA domains.  CANNON's shift-and-multiply structure ensures all
//      domains work on local data simultaneously, with global barriers
//      between steps so pages migrate together.
//
//   3. MR=8 register blocking.  8 output rows × 32 columns = 16 ZMM
//      accumulators.  Each B row loaded once, broadcast-multiplied by 8
//      different A scalars.  8× reduction in B traffic from L1.
//
//   4. fp32 accumulation.  fp16 inputs, fp32 accumulator.  Required for
//      numerical correctness: fp16 has 3.3 decimal digits, K=8192
//      additions would overflow.  This matches AMX _tdpfp16ps semantics.
//
// Build:
//   icpx -O3 -g -qopenmp -march=native -ffast-math -std=c++17 \
//       -I${MKLROOT}/include -o bench_goto_vec bench_goto_vec.cpp \
//       -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread \
//       -lmkl_core -liomp5 -lpthread -lm
//
#include <mkl.h>
#include <omp.h>
#include <sched.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

using IN_T  = _Float16;
using OUT_T = float;

// ═══════════════════════════════════════════════════════════════════════════
// §1  CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════

// DECISION: BM=BN=32 matches AVX-512 sweet spot (32 fp32 = 2 ZMM vectors).
// BK=32 gives 32 k-steps per rank-update, enough to amortize loop overhead.
static constexpr int BM = 32;
static constexpr int BN = 32;
static constexpr int BK = 32;

// DECISION: MR=8.  8 rows × 2 ZMM = 16 accumulators out of 32 ZMM.
// MR=4: only 8 accumulators, wastes register file.
// MR=16: needs 32 accumulators, guaranteed register spills.
static constexpr int MR = 8;

// DECISION: KC controls the rank-update depth.  MC controls IC-parallel
// granularity.  NC is irrelevant when N_local ≤ NC (JC loop runs once).
// Best empirical config from sweep: MC=256, KC=512.
static constexpr int DEFAULT_MC = 256;
static constexpr int DEFAULT_KC = 512;

static int  g_warmup  = 3;
static int  g_iters   = 20;
static int  g_threads = 96;
static int  g_MC      = DEFAULT_MC;
static int  g_KC      = DEFAULT_KC;
static const char* g_outfile = "bench_goto_vec.csv";

static inline double now_ns() {
    auto t = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::nano>(t.time_since_epoch()).count();
}

// ═══════════════════════════════════════════════════════════════════════════
// §2  NUMA UTILITIES
// ═══════════════════════════════════════════════════════════════════════════

static int get_num_numa_nodes() {
    int n = 0; char path[256];
    for (int i = 0; i < 64; ++i) {
        snprintf(path, sizeof(path), "/sys/devices/system/node/node%d", i);
        if (access(path, F_OK) == 0) ++n; else break;
    }
    return n > 0 ? n : 1;
}

static std::vector<int> cpus_on_node(int node) {
    std::vector<int> cpus; char path[256];
    snprintf(path, sizeof(path), "/sys/devices/system/node/node%d/cpulist", node);
    FILE* f = fopen(path, "r"); if (!f) return cpus;
    char buf[4096];
    if (!fgets(buf, sizeof(buf), f)) { fclose(f); return cpus; }
    fclose(f);
    char* end = buf + strlen(buf);
    while (end > buf && (end[-1]=='\n'||end[-1]=='\r'||end[-1]==' ')) *--end = '\0';
    char* p = buf;
    while (*p) {
        char* next = p;
        int lo = (int)strtol(p, &next, 10);
        if (next == p) break;
        p = next; int hi = lo;
        if (*p == '-') { ++p; hi = (int)strtol(p, &p, 10); }
        for (int c = lo; c <= hi; ++c) cpus.push_back(c);
        if (*p == ',') ++p;
    }
    std::sort(cpus.begin(), cpus.end());
    return cpus;
}

static void* numa_alloc(size_t bytes) {
    void* p = mmap(nullptr, bytes, PROT_READ|PROT_WRITE,
                   MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return p;
}
static void numa_free(void* p, size_t bytes) { if (p && bytes) munmap(p, bytes); }
template<typename T> T* typed_alloc(size_t n) { return (T*)numa_alloc(n*sizeof(T)); }
template<typename T> void typed_free(T* p, size_t n) { numa_free(p, n*sizeof(T)); }

static void pin_to_cpu(int cpu) {
    cpu_set_t s; CPU_ZERO(&s); CPU_SET(cpu, &s);
    sched_setaffinity(0, sizeof(s), &s);
}

// ═══════════════════════════════════════════════════════════════════════════
// §3  NO-PACK MICROKERNEL
// ═══════════════════════════════════════════════════════════════════════════
//
// DECISION: No packing.  The microkernel reads A and B with their natural
// strides (lda, ldb).  This eliminates ~280 MB of redundant copies per
// GEMM call.  It works because:
//
//   A access: scalar A[row * lda + k], one element per k-step, broadcast.
//     Scalar loads are stride-insensitive.  The hardware does not care
//     whether lda=32 (packed) or lda=4096 (natural).
//
//   B access: B[k * ldb + j] for j=0..BN-1, consecutive in memory.
//     Row-major B means BN consecutive fp16 values are contiguous.
//     The j-loop vectorizes as: load 16 fp16 → vcvtph2ps → vfmadd231ps.
//     Stride ldb between k-steps is handled by the hardware prefetcher.
//
// ── No-pack rank-BK update for MR=8 rows ────────────────────────────────
// C[i0:i0+MR, j0:j0+BN] += A[i0:i0+MR, k0:k0+BK] × B[k0:k0+BK, j0:j0+BN]
//
static inline void microkernel_mr8(
    const IN_T* __restrict__ A, int lda,
    const IN_T* __restrict__ B, int ldb,
    int i0, int j0, int k0,
    float* __restrict__ d0, float* __restrict__ d1,
    float* __restrict__ d2, float* __restrict__ d3,
    float* __restrict__ d4, float* __restrict__ d5,
    float* __restrict__ d6, float* __restrict__ d7)
{
    for (int k = k0; k < k0 + BK; ++k) {
        // Scalar loads from A: one per row, broadcast.
        // lda-strided, but the access is one element — no cache line waste.
        const float a0 = (float)A[(i0+0)*lda + k];
        const float a1 = (float)A[(i0+1)*lda + k];
        const float a2 = (float)A[(i0+2)*lda + k];
        const float a3 = (float)A[(i0+3)*lda + k];
        const float a4 = (float)A[(i0+4)*lda + k];
        const float a5 = (float)A[(i0+5)*lda + k];
        const float a6 = (float)A[(i0+6)*lda + k];
        const float a7 = (float)A[(i0+7)*lda + k];

        // Vector load from B: BN=32 consecutive fp16, contiguous in memory.
        const IN_T* __restrict__ brow = &B[k * ldb + j0];
        #pragma omp simd
        for (int j = 0; j < BN; ++j) {
            const float bv = (float)brow[j];
            d0[j] += a0 * bv;  d1[j] += a1 * bv;
            d2[j] += a2 * bv;  d3[j] += a3 * bv;
            d4[j] += a4 * bv;  d5[j] += a5 * bv;
            d6[j] += a6 * bv;  d7[j] += a7 * bv;
        }
    }
}

// ── Full tile: loop over KC in BK steps, store to C ─────────────────────
static inline void tile_gemm(
    const IN_T* __restrict__ A, int lda,
    const IN_T* __restrict__ B, int ldb,
    OUT_T* __restrict__ C, int ldc,
    int i0, int j0, int k0, int kc,
    bool accum)
{
    for (int ri = 0; ri < BM; ri += MR) {
        alignas(64) float d0[BN]={},d1[BN]={},d2[BN]={},d3[BN]={};
        alignas(64) float d4[BN]={},d5[BN]={},d6[BN]={},d7[BN]={};

        // Accumulate over KC in steps of BK.
        for (int kk = k0; kk < k0 + kc; kk += BK) {
            microkernel_mr8(A, lda, B, ldb, i0+ri, j0, kk,
                            d0,d1,d2,d3,d4,d5,d6,d7);
        }

        // Store to C.
        float* rows[MR] = {d0,d1,d2,d3,d4,d5,d6,d7};
        for (int r = 0; r < MR; ++r) {
            OUT_T* __restrict__ ci = &C[(i0+ri+r)*ldc + j0];
            const float* __restrict__ di = rows[r];
            if (accum) {
                #pragma omp simd
                for (int j = 0; j < BN; ++j) ci[j] += di[j];
            } else {
                #pragma omp simd
                for (int j = 0; j < BN; ++j) ci[j] = di[j];
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §4  GOTO 5-LOOP (no packing, direct array access)
// ═══════════════════════════════════════════════════════════════════════════
//
// DECISION: Eliminated JC/NC loop entirely.  With 2×2 grid, N_local=4096.
// The entire N_local is processed without column partitioning.  The loop
// nest reduces to: PC → IC(parallel) → JR → IR.
//
// DECISION: No Bp/Ap buffers allocated.  The microkernel reads A_local
// and B_local directly.  Workspace is zero.

static void goto_gemm_nopack(
    const IN_T* __restrict__ A, int lda,  // [M_local × K_dim]
    const IN_T* __restrict__ B, int ldb,  // [K_dim × N_local]
    OUT_T* __restrict__ C, int ldc,       // [M_local × N_local]
    int M, int N, int K,
    int nthreads, int KC)
{
    const int tn = N / BN;   // tiles in N direction

    for (int pc = 0; pc < K; pc += KC) {
        const int kc = std::min(KC, K - pc);
        const bool accum = (pc > 0);

        #pragma omp parallel num_threads(nthreads) proc_bind(close)
        {
            int tid = omp_get_thread_num();
            int nth = omp_get_num_threads();

            // IC partition: divide M rows across threads.
            int ic_total = M / BM;  // number of BM-row blocks
            int ic_base = ic_total / nth, ic_extra = ic_total % nth;
            int ic_start = tid * ic_base + std::min(tid, ic_extra);
            int ic_count = ic_base + (tid < ic_extra ? 1 : 0);

            for (int ic_idx = ic_start; ic_idx < ic_start + ic_count; ++ic_idx) {
                int i0 = ic_idx * BM;
                // JR loop: sweep across N in BN-wide tiles.
                for (int jr = 0; jr < tn; ++jr) {
                    int j0 = jr * BN;
                    tile_gemm(A, lda, B, ldb, C, ldc,
                              i0, j0, pc, kc, accum);
                }
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §5  NUMA GRID
// ═══════════════════════════════════════════════════════════════════════════

struct NumaDomain {
    int node_id;
    std::vector<int> cpus;
    int cores;
};

struct NumaGrid {
    int P;                           // number of NUMA domains
    int total_threads;
    std::vector<NumaDomain> domains; // P domains

    void init(int p, int total_cores) {
        P = p; total_threads = total_cores;
        domains.resize(p);
        int num_nodes = get_num_numa_nodes();
        int cores_per = total_cores / p;
        for (int d = 0; d < p; ++d) {
            domains[d].node_id = d % num_nodes;
            domains[d].cores = cores_per;
            auto all = cpus_on_node(domains[d].node_id);
            int off = (d / num_nodes) * cores_per;
            for (int i = 0; i < cores_per && (off+i) < (int)all.size(); ++i)
                domains[d].cpus.push_back(all[off + i]);
        }
    }
};

// ═══════════════════════════════════════════════════════════════════════════
// §6  CANNON OVER NUMA DOMAINS
// ═══════════════════════════════════════════════════════════════════════════
//
// DECISION: CANNON (not SUMMA) for the 2×2 case.
//
//   CANNON requires a square P×P grid and K divisible by P.
//   For P=2, K=8192: K_blk = 4096.  Exactly 2 steps.
//
//   Advantages over SUMMA:
//     • Constant memory per domain: one M_loc × K_blk A-block and one
//       K_blk × N_loc B-block.  No growing panel buffers.
//     • Nearest-neighbor shifts (not broadcasts).
//     • All domains shift simultaneously → pages migrate together.
//
//   Algorithm:
//     1. Partition C into P×P blocks.  Domain (i,j) owns C(i,j).
//     2. Initial skew: domain (i,j) gets A(i, (j+i)%P) and B((i+j)%P, j).
//     3. For step = 0 to P-1:
//          a. C_local += A_block × B_block   (local GOTO, all domains)
//          b. Shift A blocks LEFT by 1 in each row (circular).
//          c. Shift B blocks UP by 1 in each column (circular).
//          d. Global barrier (all threads sync before next step).
//
// DECISION: All buffers are allocated and first-touched inside the flat
// parallel region so pages land on the correct NUMA node.  The "shift"
// is a memcpy into a double-buffer, with the destination first-touched
// by the receiving domain's thread → pages migrate to the receiver.
//
static void cannon_gemm(const IN_T* A, const IN_T* B, OUT_T* C,
                        int M, int N, int K,
                        NumaGrid& grid, int KC) {
    const int nd = grid.P;                       // total domains (e.g. 4)
    const int Pc = (int)sqrt((double)nd);        // grid side (e.g. 2)
    const int Ml = M / Pc, Nl = N / Pc, Kl = K / Pc;
    const int cores_per = grid.domains[0].cores;

    std::vector<IN_T*>  A_blk(nd), B_blk(nd), A_tmp(nd), B_tmp(nd);
    std::vector<OUT_T*> C_loc(nd);
    for (int d = 0; d < nd; ++d) {
        A_blk[d] = typed_alloc<IN_T>((size_t)Ml * Kl);
        B_blk[d] = typed_alloc<IN_T>((size_t)Kl * Nl);
        A_tmp[d] = typed_alloc<IN_T>((size_t)Ml * Kl);
        B_tmp[d] = typed_alloc<IN_T>((size_t)Kl * Nl);
        C_loc[d] = typed_alloc<OUT_T>((size_t)Ml * Nl);
    }

    #pragma omp parallel num_threads(nd * cores_per) proc_bind(close)
    {
        const int gtid = omp_get_thread_num();
        const int d    = gtid / cores_per;
        const int ltid = gtid % cores_per;
        const int pi = d / Pc, pj = d % Pc;
        auto& dom = grid.domains[d];

        // Pin thread to its domain's CPU.
        if (ltid < (int)dom.cpus.size()) pin_to_cpu(dom.cpus[ltid]);

        // ── First-touch: zero C_local, copy initial skewed A/B blocks ──
        // CANNON initial skew: domain (pi,pj) gets:
        //   A_block from column-block (pj+pi) % Pc
        //   B_block from row-block    (pi+pj) % Pc
        {
            int a_col = (pj + pi) % Pc;
            int b_row = (pi + pj) % Pc;

            // Parallel first-touch copy of A_block.
            size_t a_elts = (size_t)Ml * Kl;
            size_t a_lo = ltid * a_elts / cores_per;
            size_t a_hi = (ltid+1) * a_elts / cores_per;
            for (size_t idx = a_lo; idx < a_hi; ++idx) {
                int r = (int)(idx / Kl), c = (int)(idx % Kl);
                A_blk[d][idx] = A[(pi*Ml + r) * K + a_col*Kl + c];
            }

            // Parallel first-touch copy of B_block.
            size_t b_elts = (size_t)Kl * Nl;
            size_t b_lo = ltid * b_elts / cores_per;
            size_t b_hi = (ltid+1) * b_elts / cores_per;
            for (size_t idx = b_lo; idx < b_hi; ++idx) {
                int r = (int)(idx / Nl), c = (int)(idx % Nl);
                B_blk[d][idx] = B[(b_row*Kl + r) * N + pj*Nl + c];
            }

            // Zero C_local.
            size_t c_elts = (size_t)Ml * Nl;
            size_t c_lo = ltid * c_elts / cores_per;
            size_t c_hi = (ltid+1) * c_elts / cores_per;
            memset(C_loc[d] + c_lo, 0, (c_hi - c_lo) * sizeof(OUT_T));

            // First-touch the double-buffers too (for shift step).
            memset(A_tmp[d] + a_lo, 0, (a_hi - a_lo) * sizeof(IN_T));
            memset(B_tmp[d] + b_lo, 0, (b_hi - b_lo) * sizeof(IN_T));
        }

        #pragma omp barrier  // all domains ready

        // ── P steps of local-GEMM + shift ──────────────────────────
        for (int step = 0; step < Pc; ++step) {
            // (a) Local GOTO GEMM: C_local += A_block × B_block
            //     All cores_per threads in this domain collaborate via
            //     IC-parallel partitioning (no packing needed).
            {
                int tn = Nl / BN;
                int ic_total = Ml / BM;
                int ic_base = ic_total / cores_per;
                int ic_extra = ic_total % cores_per;
                int ic_start = ltid * ic_base + std::min(ltid, ic_extra);
                int ic_count = ic_base + (ltid < ic_extra ? 1 : 0);

                for (int pc = 0; pc < Kl; pc += KC) {
                    int kc = std::min(KC, Kl - pc);
                    bool accum = (step > 0 || pc > 0);

                    for (int ic = ic_start; ic < ic_start + ic_count; ++ic) {
                        int i0 = ic * BM;
                        for (int jr = 0; jr < tn; ++jr) {
                            int j0 = jr * BN;
                            tile_gemm(A_blk[d], Kl,   // lda = Kl
                                      B_blk[d], Nl,   // ldb = Nl
                                      C_loc[d], Nl,   // ldc = Nl
                                      i0, j0, pc, kc, accum);
                        }
                    }
                }
            }

            #pragma omp barrier  // all domains finish local GEMM

            // (b) Shift A left, B up (except after last step).
            if (step < Pc - 1) {
                // A: domain (pi, pj) sends to (pi, (pj-1+P)%P)
                int dst_a = pi * Pc + (pj - 1 + Pc) % Pc;
                {
                    size_t a_elts = (size_t)Ml * Kl;
                    size_t lo = ltid * a_elts / cores_per;
                    size_t hi = (ltid+1) * a_elts / cores_per;
                    // Copy into RECEIVER's tmp buffer → first-touch on receiver's node.
                    memcpy(A_tmp[dst_a] + lo, A_blk[d] + lo, (hi-lo)*sizeof(IN_T));
                }

                // B: domain (pi, pj) sends to ((pi-1+P)%P, pj)
                int dst_b = ((pi - 1 + Pc) % Pc) * Pc + pj;
                {
                    size_t b_elts = (size_t)Kl * Nl;
                    size_t lo = ltid * b_elts / cores_per;
                    size_t hi = (ltid+1) * b_elts / cores_per;
                    memcpy(B_tmp[dst_b] + lo, B_blk[d] + lo, (hi-lo)*sizeof(IN_T));
                }

                #pragma omp barrier  // all shifts complete

                // Swap: A_blk[d] ↔ A_tmp[d], B_blk[d] ↔ B_tmp[d]
                // Only one thread per domain needs to do the pointer swap.
                #pragma omp single
                {
                    std::swap(A_blk, A_tmp);
                    std::swap(B_blk, B_tmp);
                }
                // implicit barrier after omp single
            }
        }
    } // end flat parallel

    // ── Gather C blocks back to global layout ───────────────────────
    for (int pi = 0; pi < Pc; ++pi)
        for (int pj = 0; pj < Pc; ++pj) {
            int d = pi * Pc + pj;
            for (int i = 0; i < Ml; ++i)
                memcpy(&C[(pi*Ml+i)*N + pj*Nl],
                       &C_loc[d][i*Nl], Nl * sizeof(OUT_T));
        }

    for (int d = 0; d < nd; ++d) {
        typed_free(A_blk[d], (size_t)Ml*Kl); typed_free(B_blk[d], (size_t)Kl*Nl);
        typed_free(A_tmp[d], (size_t)Ml*Kl); typed_free(B_tmp[d], (size_t)Kl*Nl);
        typed_free(C_loc[d], (size_t)Ml*Nl);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §6b  GOTO INNER LOOP (no OMP parallel — called by threads in a region)
// ═══════════════════════════════════════════════════════════════════════════
//
// Same loop structure as goto_gemm_nopack, but called by individual
// threads inside an existing parallel region.  Takes (ltid, nth) instead
// of creating its own #pragma omp parallel.

static void goto_inner(
    const IN_T* __restrict__ A, int lda,
    const IN_T* __restrict__ B, int ldb,
    OUT_T* __restrict__ C, int ldc,
    int M, int N, int K,
    int ltid, int nth, int KC)
{
    const int tn = N / BN;
    const int ic_total = M / BM;
    int ic_base = ic_total / nth, ic_extra = ic_total % nth;
    int ic_start = ltid * ic_base + std::min(ltid, ic_extra);
    int ic_count = ic_base + (ltid < ic_extra ? 1 : 0);

    for (int pc = 0; pc < K; pc += KC) {
        int kc = std::min(KC, K - pc);
        bool accum = (pc > 0);
        for (int ic = ic_start; ic < ic_start + ic_count; ++ic) {
            int i0 = ic * BM;
            for (int jr = 0; jr < tn; ++jr)
                tile_gemm(A, lda, B, ldb, C, ldc,
                          i0, jr * BN, pc, kc, accum);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §6c  OUTER-PRODUCT SWEEP OVER P NUMA DOMAINS
// ═══════════════════════════════════════════════════════════════════════════
//
// DECISION: 1D distribution — works for ANY number of NUMA domains.
//
//   Distribute:
//     A by rows:     domain d owns A[d·M/P : (d+1)·M/P, :]   → (M/P) × K
//     B by columns:  domain d owns B[:, d·N/P : (d+1)·N/P]   → K × (N/P)
//     C by rows:     domain d owns C[d·M/P : (d+1)·M/P, :]   → (M/P) × N
//
//   Outer product sweep (P steps):
//     step s: each domain d computes
//       C_d[:, src_cols] += A_d × B_src     where src = (d + s) % P
//
//   After P steps, every column block of C has been computed.
//   Total FLOPs = P × P × 2·(M/P)·(N/P)·K = 2·M·N·K.  Correct.
//
// DECISION: Global barrier between steps.  Not strictly necessary (B is
//   read-only after scatter), but keeps all domain threads streaming the
//   same B_src → better L3 cache utilization within each CCD.
//
// DECISION: A_d points directly into global A (contiguous row strip →
//   no copy needed, NUMA-local via first-touch init).  Only B is copied
//   into per-domain NUMA-local buffers.  C_d is a per-domain buffer
//   because the outer-product writes to non-contiguous column slices.
//

static void outerproduct_gemm(const IN_T* A, const IN_T* B, OUT_T* C,
                               int M, int N, int K,
                               NumaGrid& grid, int KC) {
    const int P = grid.P;
    const int Ml = M / P;          // rows per domain
    const int Nl = N / P;          // columns per domain
    const int cores_per = grid.domains[0].cores;

    // Per-domain buffers.
    //   B_loc[d]: K × Nl, NUMA-local column strip.
    //   C_loc[d]: Ml × N, full row strip (writes to different col slices per step).
    // A is NOT copied — domain d reads A + d*Ml*K directly (contiguous rows,
    // already NUMA-local from the first-touch init in main).
    std::vector<IN_T*>  B_loc(P);
    std::vector<OUT_T*> C_loc(P);
    for (int d = 0; d < P; ++d) {
        B_loc[d] = typed_alloc<IN_T>((size_t)K * Nl);
        C_loc[d] = typed_alloc<OUT_T>((size_t)Ml * N);
    }

    #pragma omp parallel num_threads(P * cores_per) proc_bind(close)
    {
        const int gtid = omp_get_thread_num();
        const int d    = gtid / cores_per;
        const int ltid = gtid % cores_per;
        auto& dom = grid.domains[d];

        if (ltid < (int)dom.cpus.size()) pin_to_cpu(dom.cpus[ltid]);

        // ── First-touch scatter ─────────────────────────────────────
        // B columns: non-contiguous in source (stride N), contiguous in dest (stride Nl).
        {
            int col_lo = ltid * Nl / cores_per;
            int col_hi = (ltid + 1) * Nl / cores_per;
            int ncols = col_hi - col_lo;
            if (ncols > 0)
                for (int k = 0; k < K; ++k)
                    memcpy(&B_loc[d][(size_t)k * Nl + col_lo],
                           &B[(size_t)k * N + d * Nl + col_lo],
                           ncols * sizeof(IN_T));
        }

        // C: zero entire Ml × N buffer.
        {
            size_t elts = (size_t)Ml * N;
            size_t lo = ltid * elts / cores_per;
            size_t hi = (ltid + 1) * elts / cores_per;
            memset(C_loc[d] + lo, 0, (hi - lo) * sizeof(OUT_T));
        }

        #pragma omp barrier  // scatter complete

        // ── P-step outer product sweep ──────────────────────────────
        for (int step = 0; step < P; ++step) {
            int src = (d + step) % P;

            // Sub-GEMM: C_loc[d][:, src*Nl..(src+1)*Nl] += A_d × B_loc[src]
            //
            // A_d:       read from global A, rows d*Ml..(d+1)*Ml.
            //            lda = K (global row stride).  NUMA-local (first-touch).
            // B_loc[src]: K × Nl, stride Nl.  At step 0, src=d → local.
            //            At other steps, src ≠ d → remote NUMA read.
            //            The hardware prefetcher hides the latency because
            //            GEMM is compute-bound (AI ≈ 1365 FLOP/byte).
            // C_loc[d]:  Ml × N, stride N.  Column offset = src * Nl.
            //            Always local.
            goto_inner(A + (size_t)d * Ml * K, K,      // A_d, lda=K
                       B_loc[src], Nl,                  // B_src, ldb=Nl
                       C_loc[d] + src * Nl, N,          // C_d + col offset, ldc=N
                       Ml, Nl, K,
                       ltid, cores_per, KC);

            #pragma omp barrier  // all domains finish step before next
        }
    } // end parallel

    // ── Gather C ────────────────────────────────────────────────────
    for (int d = 0; d < P; ++d)
        for (int i = 0; i < Ml; ++i)
            memcpy(&C[(d * Ml + i) * N], &C_loc[d][i * N], N * sizeof(OUT_T));

    for (int d = 0; d < P; ++d) {
        typed_free(B_loc[d], (size_t)K * Nl);
        typed_free(C_loc[d], (size_t)Ml * N);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §7  STANDALONE GOTO (no NUMA distribution, for comparison)
// ═══════════════════════════════════════════════════════════════════════════

static void goto_gemm_flat(const IN_T* A, const IN_T* B, OUT_T* C,
                           int M, int N, int K, int nthreads, int KC) {
    // Zero C with NUMA first-touch.
    #pragma omp parallel num_threads(nthreads) proc_bind(close)
    {
        int t = omp_get_thread_num(), n = omp_get_num_threads();
        size_t lo = (size_t)t * M * N / n, hi = (size_t)(t+1) * M * N / n;
        memset(C + lo, 0, (hi - lo) * sizeof(OUT_T));
    }
    goto_gemm_nopack(A, K, B, N, C, N, M, N, K, nthreads, KC);
}

// ═══════════════════════════════════════════════════════════════════════════
// §8  MKL BASELINE
// ═══════════════════════════════════════════════════════════════════════════

static float *g_A32 = nullptr, *g_B32 = nullptr;

static void mkl_sgemm(const IN_T* A, const IN_T* B, OUT_T* C,
                      int M, int N, int K) {
    mkl_set_dynamic(0); mkl_set_num_threads(g_threads);
    #pragma omp parallel for num_threads(g_threads) schedule(static)
    for (size_t i = 0; i < (size_t)M*K; ++i) g_A32[i] = (float)A[i];
    #pragma omp parallel for num_threads(g_threads) schedule(static)
    for (size_t i = 0; i < (size_t)K*N; ++i) g_B32[i] = (float)B[i];
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0f, g_A32, K, g_B32, N, 0.0f, C, N);
}

// ═══════════════════════════════════════════════════════════════════════════
// §9  VERIFICATION + BENCHMARK
// ═══════════════════════════════════════════════════════════════════════════

static bool verify(const OUT_T* C, const OUT_T* ref, int M, int N,
                   const char* label) {
    int bad = 0; double max_abs = 0, sum_abs = 0;
    size_t total = (size_t)M * N;
    for (size_t i = 0; i < total; ++i) {
        double d = fabs((double)C[i] - (double)ref[i]);
        double den = fabs((double)ref[i]);
        double rel = (den > 1e-6) ? d / den : d;
        sum_abs += d; if (d > max_abs) max_abs = d;
        if (rel > 0.05 && d > 0.1f) ++bad;
    }
    printf("    [%s] %-20s max_abs=%.6f mean=%.6f\n",
           bad ? "FAIL" : "PASS", label, max_abs, sum_abs / total);
    fflush(stdout);
    return bad == 0;
}

using KernFn = std::function<void(const IN_T*, const IN_T*, OUT_T*, int, int, int)>;

static void bench_one(FILE* f, const char* name, KernFn fn,
                      const IN_T* A, const IN_T* B, OUT_T* C,
                      int M, int N, int K) {
    double flops = 2.0 * (double)M * N * K;
    for (int w = 0; w < g_warmup; ++w) {
        double t0 = now_ns(); fn(A,B,C,M,N,K); double t1 = now_ns();
        printf("    %-20s warmup %d/%d  %.2fs  %.0f GF/s\n",
               name, w+1, g_warmup, (t1-t0)*1e-9, flops/(t1-t0));
        fflush(stdout);
    }
    std::vector<double> times(g_iters);
    for (int i = 0; i < g_iters; ++i) {
        double t0 = now_ns(); fn(A,B,C,M,N,K); double t1 = now_ns();
        times[i] = t1 - t0;
        fprintf(f, "%s,%d,%d,%d,%d,%d,%.2f,%.6f\n",
                name, M, N, K, g_threads, i, times[i], flops/times[i]);
        fflush(f);
    }
    std::sort(times.begin(), times.end());
    double med = times[g_iters/2];
    printf("    %-20s  med=%.2fs  %.0f GF/s (peak %.0f)\n",
           name, med*1e-9, flops/med, flops/times[0]);
    fflush(stdout);
}

// ═══════════════════════════════════════════════════════════════════════════
// §10  MAIN
// ═══════════════════════════════════════════════════════════════════════════

struct MNK { int M, N, K; };

static std::vector<MNK> parse_args(int argc, char** argv) {
    std::vector<MNK> sizes; std::vector<int> nums;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i],"-w") && i+1<argc) g_warmup = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-i") && i+1<argc) g_iters = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-t") && i+1<argc) g_threads = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-mc") && i+1<argc) g_MC = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-kc") && i+1<argc) g_KC = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-o") && i+1<argc) g_outfile = argv[++i];
        else nums.push_back(atoi(argv[i]));
    }
    for (size_t i = 0; i+2 < nums.size(); i += 3)
        sizes.push_back({nums[i], nums[i+1], nums[i+2]});
    return sizes;
}

int main(int argc, char** argv) {
    setvbuf(stdout, nullptr, _IOLBF, 0);
    if (argc < 4) {
        printf("Usage: %s M N K [-t T] [-mc MC] [-kc KC] [-w W] [-i I] [-o F]\n", argv[0]);
        return 1;
    }
    auto sizes = parse_args(argc, argv);

    // DECISION: P = number of NUMA nodes.  Works for any P (1, 2, 3, 4, 6, ...).
    // CANNON additionally requires P to be a perfect square.
    int nn = get_num_numa_nodes();
    int P = nn;
    bool can_cannon = false;
    {   int sq = (int)sqrt((double)P);
        can_cannon = (sq * sq == P && P > 1); }

    printf("════════════════════════════════════════════════════════\n");
    printf("  GOTO GEMM — no-pack, MR=%d, P=%d NUMA domains\n", MR, P);
    printf("  Threads=%d  Cores/dom=%d  KC=%d\n",
           g_threads, g_threads / P, g_KC);
    printf("  Outer loops: flat, outerproduct_%d%s\n",
           P, can_cannon ? ", cannon" : "");
    printf("════════════════════════════════════════════════════════\n\n");

    FILE* fout = fopen(g_outfile, "w");
    if (!fout) { perror("fopen"); return 1; }
    fprintf(fout, "kernel,M,N,K,threads,run,time_ns,gflops\n");

    NumaGrid grid;
    grid.init(P, g_threads);

    // For CANNON we need a separate square grid.
    NumaGrid cannon_grid;
    int Pc = (int)sqrt((double)P);
    if (can_cannon) cannon_grid.init(Pc * Pc, g_threads);

    for (auto& sz : sizes) {
        int M = sz.M, N = sz.N, K = sz.K;
        printf("═══ %d × %d × %d ═══\n", M, N, K);
        if (M % (P*BM) || N % (P*BN) || K % BK) {
            printf("  [SKIP] divisibility\n"); continue;
        }

        size_t szA = (size_t)M*K, szB = (size_t)K*N, szC = (size_t)M*N;
        IN_T*  A    = typed_alloc<IN_T>(szA);
        IN_T*  B    = typed_alloc<IN_T>(szB);
        OUT_T* Cref = typed_alloc<OUT_T>(szC);
        OUT_T* Cchk = typed_alloc<OUT_T>(szC);
        g_A32       = typed_alloc<float>(szA);
        g_B32       = typed_alloc<float>(szB);

        // NUMA first-touch init: distribute pages across all nodes.
        #pragma omp parallel num_threads(g_threads) proc_bind(close)
        {
            int t = omp_get_thread_num(), n = omp_get_num_threads();
            unsigned seed = 42 + t;
            size_t lo, hi;
            lo = t*szA/n; hi = (t+1)*szA/n;
            for (size_t i = lo; i < hi; ++i)
                A[i] = (IN_T)((float)rand_r(&seed)/RAND_MAX - 0.5f);
            lo = t*szB/n; hi = (t+1)*szB/n;
            for (size_t i = lo; i < hi; ++i)
                B[i] = (IN_T)((float)rand_r(&seed)/RAND_MAX - 0.5f);
            lo = t*szA/n; hi = (t+1)*szA/n;
            memset(g_A32+lo, 0, (hi-lo)*sizeof(float));
            lo = t*szB/n; hi = (t+1)*szB/n;
            memset(g_B32+lo, 0, (hi-lo)*sizeof(float));
            lo = t*szC/n; hi = (t+1)*szC/n;
            memset(Cref+lo, 0, (hi-lo)*sizeof(OUT_T));
            memset(Cchk+lo, 0, (hi-lo)*sizeof(OUT_T));
        }

        // Reference
        printf("  MKL reference...\n"); fflush(stdout);
        mkl_sgemm(A, B, Cref, M, N, K);

        printf("  Kernels:\n"); fflush(stdout);

        // MKL benchmark
        bench_one(fout, "mkl_sgemm", mkl_sgemm, A, B, Cchk, M, N, K);
        verify(Cchk, Cref, M, N, "mkl_sgemm");

        // Flat GOTO (no NUMA distribution, for comparison)
        auto flat_fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
            goto_gemm_flat(a, b, c, m, n, k, g_threads, g_KC);
        };
        bench_one(fout, "goto_flat", flat_fn, A, B, Cchk, M, N, K);
        verify(Cchk, Cref, M, N, "goto_flat");

        // Outer-product sweep over P NUMA domains (works for ANY P)
        {
            char lbl[64]; snprintf(lbl, sizeof(lbl), "outerproduct_%d", P);
            auto op_fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
                outerproduct_gemm(a, b, c, m, n, k, grid, g_KC);
            };
            bench_one(fout, lbl, op_fn, A, B, Cchk, M, N, K);
            verify(Cchk, Cref, M, N, lbl);
        }

        // CANNON (only if P is a perfect square and K divisible by sqrt(P))
        if (can_cannon && K % (Pc * BK) == 0 && M % (Pc * BM) == 0 && N % (Pc * BN) == 0) {
            char lbl[64]; snprintf(lbl, sizeof(lbl), "cannon_%dx%d", Pc, Pc);
            auto cn_fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
                cannon_gemm(a, b, c, m, n, k, cannon_grid, g_KC);
            };
            bench_one(fout, lbl, cn_fn, A, B, Cchk, M, N, K);
            verify(Cchk, Cref, M, N, lbl);
        }

        typed_free(A, szA); typed_free(B, szB);
        typed_free(Cref, szC); typed_free(Cchk, szC);
        typed_free(g_A32, szA); typed_free(g_B32, szB);
        g_A32 = g_B32 = nullptr;
        printf("\n"); fflush(stdout);
    }

    fclose(fout);
    printf("Done. Results: %s\n", g_outfile);
    return 0;
}