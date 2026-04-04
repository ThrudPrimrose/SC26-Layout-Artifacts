// bench_goto_vec.cpp
// ═══════════════════════════════════════════════════════════════════════════
// GOTO GEMM with auto-vectorized microkernel + SUMMA / Cannon over NUMA
// ═══════════════════════════════════════════════════════════════════════════
//
// Implements the GOTO 5-loop nest with a pure C microkernel that the
// compiler auto-vectorizes via #pragma omp simd and -O3 -march=native.
// No inline intrinsics, no AMX, no vendor tile extensions.
// Portable: works on AVX-512, AVX2, SVE — whatever -march=native implies.
//
// Target: AMD Zen4 (MI300A CPU side, EPYC 9004 "Genoa")
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  THREE PRIMITIVES (matching the AMX tile.hpp interface)                 │
// │                                                                        │
// │  Peak (init)   = NO-OP.  No tile config register, no special state.    │
// │                                                                        │
// │  MMAD (mma)    = Plain C loop with #pragma omp simd.                   │
// │                  Computes D[32×32] += A[32×32] · B[32×32] using       │
// │                  ikj loop order: broadcast A[i,k], stream B[k,0..31]. │
// │                  The compiler auto-vectorizes the inner j-loop with    │
// │                  -O3 -march=native → vfmadd231ps on AVX-512 targets.  │
// │                  No inline intrinsics anywhere in this file.            │
// │                                                                        │
// │  Unpack (store)= NO-OP.  D is a plain fp32[32][32] array on the       │
// │                  stack.  Stored to strided C via #pragma omp simd loop.│
// │                                                                        │
// │  B packing:     Plain row-major tiles (BK×BN), NOT vnni-interleaved.  │
// │                  The SIMD j-loop reads B as BN consecutive fp16 values │
// │                  which the compiler converts to fp32 before the FMA.   │
// └─────────────────────────────────────────────────────────────────────────┘
//
// Build (AMD Zen4 / MI300A with oneMKL):
//   g++ -O3 -fopenmp -march=native -mtune=native -ffast-math -std=c++17 \
//       -fno-vect-cost-model                                               \
//       -I${MKLROOT}/include -o bench_goto_vec bench_goto_vec.cpp           \
//       -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread        \
//       -lmkl_core -lgomp -lpthread -lm
//
// -fno-vect-cost-model ensures the compiler vectorizes all #pragma omp simd
// loops unconditionally, rather than second-guessing profitability.
//
// Run (MI300A: 4 NUMA nodes × 24 cores = 96 threads):
//   ./bench_goto_vec 8192 8192 8192 -t 96 -p 2x2
//
#include <mkl.h>          // oneMKL — cblas_sgemm baseline
#include <omp.h>
#include <sched.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

// ═══════════════════════════════════════════════════════════════════════════
// §1  COMPILE-TIME CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════
//
// We keep the same BM=BN=BK=32 tile dimensions as the AMX version so that
// blocking parameters, packing layouts, and GOTO loop structure are
// directly comparable.  The only difference is how the 32×32×32 tile
// multiply-accumulate is implemented (vector FMA vs. AMX tdp).
//
// NOTE: _Float16 requires GCC 12+ with -march=znver4 / -march=native on Zen4.
//       Will fail if compiled without proper arch flags.
//
using IN_T  = _Float16;
using OUT_T = float;

static constexpr int BM = 32;            // tile rows   (= MR in GOTO)
static constexpr int BN = 32;            // tile cols   (= NR in GOTO)
static constexpr int BK = 32;            // tile depth
static constexpr int A_TSZ = BM * BK;    // packed-A tile: 32×32 = 1024 fp16 elements
static constexpr int B_TSZ = BK * BN;    // packed-B tile: 32×32 = 1024 fp16 elements
                                          // ^^^ plain row-major, NOT vnni

static constexpr int DEFAULT_MC = 1024;
static constexpr int DEFAULT_KC = 256;
static constexpr int DEFAULT_NC = 4096;

// ═══════════════════════════════════════════════════════════════════════════
// §2  GLOBALS
// ═══════════════════════════════════════════════════════════════════════════

static int  g_warmup  = 5;
static int  g_iters   = 50;
static int  g_threads = 32;
static int  g_PX      = 2;
static int  g_PY      = 2;
static int  g_MC      = DEFAULT_MC;
static int  g_KC      = DEFAULT_KC;
static int  g_NC      = DEFAULT_NC;
static const char* g_outfile = "bench_goto_vec.csv";

static inline double now_ns() {
    auto t = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::nano>(t.time_since_epoch()).count();
}

// ═══════════════════════════════════════════════════════════════════════════
// §3  NUMA UTILITIES  (identical to bench_goto_numa.cpp)
// ═══════════════════════════════════════════════════════════════════════════

static int get_numa_node() {
    unsigned cpu, node;
    syscall(__NR_getcpu, &cpu, &node, nullptr);
    return (int)node;
}

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
    char* p = buf;
    while (*p) {
        int lo = (int)strtol(p, &p, 10), hi = lo;
        if (*p == '-') { ++p; hi = (int)strtol(p, &p, 10); }
        for (int c = lo; c <= hi; ++c) cpus.push_back(c);
        if (*p == ',') ++p;
    }
    std::sort(cpus.begin(), cpus.end());
    return cpus;
}

static void* numa_alloc(size_t bytes) {
    void* p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (p == MAP_FAILED) {
        fprintf(stderr, "mmap failed for %zu bytes: ", bytes);
        perror("");
        std::abort();
    }
    madvise(p, bytes, MADV_HUGEPAGE);
    return p;
}
static void numa_free(void* p, size_t bytes) { if (p && bytes) munmap(p, bytes); }

template <typename T> static T* typed_alloc(size_t n) {
    return static_cast<T*>(numa_alloc(n * sizeof(T)));
}
template <typename T> static void typed_free(T* p, size_t n) {
    numa_free(p, n * sizeof(T));
}

static void pin_thread_to_cpus(const std::vector<int>& cpus) {
    cpu_set_t set; CPU_ZERO(&set);
    for (int c : cpus) CPU_SET(c, &set);
    sched_setaffinity(0, sizeof(set), &set);
}

// ═══════════════════════════════════════════════════════════════════════════
// §3b  OPENMP REFERENCE GEMM (scalar ground truth)
// ═══════════════════════════════════════════════════════════════════════════

static void ref_gemm_omp(const IN_T* __restrict__ A,
                         const IN_T* __restrict__ B,
                         OUT_T* __restrict__ C,
                         int M, int N, int K) {
    memset(C, 0, (size_t)M * N * sizeof(OUT_T));
    constexpr int BLK = 64;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int bi = 0; bi < M; bi += BLK) {
        for (int bj = 0; bj < N; bj += BLK) {
            for (int bk = 0; bk < K; bk += BLK) {
                int mi = std::min(bi + BLK, M);
                int nj = std::min(bj + BLK, N);
                int pk = std::min(bk + BLK, K);
                for (int i = bi; i < mi; ++i)
                    for (int k = bk; k < pk; ++k) {
                        float a_val = (float)A[i * K + k];
                        for (int j = bj; j < nj; ++j)
                            C[i * N + j] += a_val * (float)B[k * N + j];
                    }
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §4  PACKING ROUTINES (plain row-major — no VNNI)
// ═══════════════════════════════════════════════════════════════════════════
//
// A packing is identical to the AMX version: each BM×BK tile is copied
// contiguously in row-major order.
//
// B packing is DIFFERENT from AMX: we use plain row-major tiles (BK×BN),
// NOT vnni-interleaved.  The auto-vectorized microkernel reads B as
// BN consecutive fp16 values in the N dimension, which the compiler
// converts to fp32 before the FMA.  Plain row-major within the tile.
//

static void pack_A_block(const IN_T* __restrict__ A, IN_T* __restrict__ dst,
                         int row0, int col0, int mc, int kc, int lda) {
    const int tm_ = mc / BM, tk_ = kc / BK;
    for (int tm = 0; tm < tm_; ++tm)
        for (int tk = 0; tk < tk_; ++tk) {
            IN_T* d = dst + (size_t)(tm * tk_ + tk) * A_TSZ;
            int gr = row0 + tm * BM, gc = col0 + tk * BK;
            for (int li = 0; li < BM; ++li)
                for (int lk = 0; lk < BK; ++lk)
                    d[li * BK + lk] = A[(gr + li) * lda + (gc + lk)];
        }
}

// ── B packing: plain row-major tiles (NOT vnni) ─────────────────────────
//
// For tile (tn, tk), the packed layout is:
//   dst[(tn * tk_ + tk) * B_TSZ + k * BN + n]
//     = B[(row0 + tk*BK + k) * ldb + (col0 + tn*BN + n)]
//
// This is a simple contiguous copy of each BK×BN sub-matrix of B.
// The auto-vectorized microkernel reads B row-by-row:
//   for j in 0..BN-1:  D[i][j] += A_val * (float)B[k*BN + j]
// The compiler vectorizes this inner j-loop with #pragma omp simd.
//
static void pack_B_panel(const IN_T* __restrict__ B, IN_T* __restrict__ dst,
                         int row0, int col0, int kc, int nc, int ldb) {
    const int tk_ = kc / BK, tn_ = nc / BN;
    for (int tn = 0; tn < tn_; ++tn)
        for (int tk = 0; tk < tk_; ++tk) {
            IN_T* d = dst + (size_t)(tn * tk_ + tk) * B_TSZ;
            int gr = row0 + tk * BK, gc = col0 + tn * BN;
            for (int k = 0; k < BK; ++k)
                for (int n = 0; n < BN; ++n)
                    d[k * BN + n] = B[(gr + k) * ldb + (gc + n)];
        }
}

// Parallel packing helpers
static void pack_A_block_par(const IN_T* A, IN_T* dst,
                             int row0, int col0, int mc, int kc, int lda,
                             int tid, int nth) {
    const int tm_ = mc / BM, tk_ = kc / BK, total = tm_ * tk_;
    int base = total / nth, extra = total % nth;
    int start = tid * base + std::min(tid, extra);
    int count = base + (tid < extra ? 1 : 0);
    for (int idx = start; idx < start + count; ++idx) {
        int tm = idx / tk_, tk = idx % tk_;
        IN_T* d = dst + (size_t)(tm * tk_ + tk) * A_TSZ;
        int gr = row0 + tm * BM, gc = col0 + tk * BK;
        for (int li = 0; li < BM; ++li)
            for (int lk = 0; lk < BK; ++lk)
                d[li * BK + lk] = A[(gr + li) * lda + (gc + lk)];
    }
}

static void pack_B_panel_par(const IN_T* B, IN_T* dst,
                             int row0, int col0, int kc, int nc, int ldb,
                             int tid, int nth) {
    const int tk_ = kc / BK, tn_ = nc / BN, total = tn_ * tk_;
    int base = total / nth, extra = total % nth;
    int start = tid * base + std::min(tid, extra);
    int count = base + (tid < extra ? 1 : 0);
    for (int idx = start; idx < start + count; ++idx) {
        int tn = idx / tk_, tk = idx % tk_;
        IN_T* d = dst + (size_t)(tn * tk_ + tk) * B_TSZ;
        int gr = row0 + tk * BK, gc = col0 + tn * BN;
        for (int k = 0; k < BK; ++k)
            for (int n = 0; n < BN; ++n)
                d[k * BN + n] = B[(gr + k) * ldb + (gc + n)];
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §5  AUTO-VECTORIZED MICROKERNEL
// ═══════════════════════════════════════════════════════════════════════════
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  No inline intrinsics.  Pure C with #pragma omp simd.                  │
// │                                                                        │
// │  The compiler (-O3 -march=native -ffast-math) will:                    │
// │    1. Convert the fp16 → fp32 cast to vcvtph2ps (F16C)               │
// │    2. Broadcast the A scalar to vbroadcastss                          │
// │    3. Fuse multiply-add into vfmadd231ps (FMA3)                       │
// │    4. Vectorize the j-loop with BN=32 / SIMD_width iterations         │
// │                                                                        │
// │  On Zen4 (AVX-512, 16 floats/vector): j-loop = 2 iterations           │
// │  On Zen2/3 (AVX2, 8 floats/vector): j-loop = 4 iterations            │
// │                                                                        │
// │  This is portable: the same source works on AVX2, AVX-512, SVE, etc.  │
// │  The compiler picks the right instructions for -march=native.          │
// └─────────────────────────────────────────────────────────────────────────┘

// ── Peak: no-op ─────────────────────────────────────────────────────────
//
// AMX needs _tile_loadconfig to program tile register dimensions.
// Vector units have no such configuration — they're always ready.
static inline void vec_init() { /* no-op */ }

// ── MMAD: auto-vectorized FMA microkernel ───────────────────────────────
//
// Computes D[BM×BN] += A_tile[BM×BK] × B_tile[BK×BN], fp16→fp32.
//
// No inline intrinsics.  The compiler auto-vectorizes the inner j-loop
// with -O3 -march=native, producing AVX-512 vfmadd231ps instructions
// (or whatever the target supports).  The #pragma omp simd ensures
// vectorization even when the compiler's cost model is conservative.
//
// The loop structure is ikj (broadcast A[i,k], stream B[k,0..BN-1]):
//   - The j-loop is the SIMD dimension (BN=32 = 2 AVX-512 vectors).
//   - The i-loop iterates over output rows.
//   - The k-loop accumulates the rank-BK partial product.
//
// This handles ONE BK=32 K-tile step.  The caller loops over kc_tiles.
//
static inline void vec_mma_rank32(const IN_T* __restrict__ ap,
                                  const IN_T* __restrict__ bp,
                                  OUT_T* __restrict__ D)
{
    for (int i = 0; i < BM; ++i) {
        OUT_T* __restrict__ di = D + i * BN;
        const IN_T* __restrict__ ai = ap + i * BK;
        for (int k = 0; k < BK; ++k) {
            const float a_val = (float)ai[k];
            const IN_T* __restrict__ bk = bp + k * BN;
            #pragma omp simd
            for (int j = 0; j < BN; ++j) {
                di[j] += a_val * (float)bk[j];
            }
        }
    }
}

// ── Full microkernel: accumulate over kc_tiles K-tile steps ─────────────
//
// Overwrites C_ptr with the accumulated result.
static inline void vec_microkernel(const IN_T* __restrict__ ap,
                                   const IN_T* __restrict__ bp,
                                   OUT_T* __restrict__ C_ptr,
                                   int kc_tiles, int ldc,
                                   int a_stride, int b_stride)
{
    // D accumulator: BM×BN = 32×32 = 1024 fp32 = 4 KB on the stack.
    alignas(64) OUT_T D[BM * BN] = {};   // zero-init

    for (int tk = 0; tk < kc_tiles; ++tk) {
        vec_mma_rank32(ap, bp, D);
        ap += a_stride;
        bp += b_stride;
    }

    // ── Unpack: no-op.  D is already row-major fp32. ────────────────
    // Store D to C (D has stride BN=32, C has stride ldc).
    for (int i = 0; i < BM; ++i) {
        const OUT_T* __restrict__ di = D + i * BN;
        OUT_T* __restrict__ ci = C_ptr + i * ldc;
        #pragma omp simd
        for (int j = 0; j < BN; ++j)
            ci[j] = di[j];
    }
}

// Accumulate variant: adds MMA result to existing C.
static inline void vec_microkernel_accum(const IN_T* __restrict__ ap,
                                         const IN_T* __restrict__ bp,
                                         OUT_T* __restrict__ C_ptr,
                                         int kc_tiles, int ldc,
                                         int a_stride, int b_stride)
{
    alignas(64) OUT_T D[BM * BN] = {};

    for (int tk = 0; tk < kc_tiles; ++tk) {
        vec_mma_rank32(ap, bp, D);
        ap += a_stride;
        bp += b_stride;
    }

    // Add D to existing C.
    for (int i = 0; i < BM; ++i) {
        const OUT_T* __restrict__ di = D + i * BN;
        OUT_T* __restrict__ ci = C_ptr + i * ldc;
        #pragma omp simd
        for (int j = 0; j < BN; ++j)
            ci[j] += di[j];
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §6  GOTO GEMM — 5-LOOP NEST (vector version)
// ═══════════════════════════════════════════════════════════════════════════
//
// Identical loop structure to bench_goto_numa.cpp §6.
// Only the microkernel call is replaced: vec_microkernel instead of
// amx_microkernel.  The B packing stride changes because B_TSZ = BK*BN
// (plain) instead of (BK/2)*(2*BN) (vnni) — but both equal 1024.
//

struct GotoParams { int MC, KC, NC; };

static void goto_gemm(const IN_T* __restrict__ A,
                      const IN_T* __restrict__ B,
                      OUT_T* __restrict__ C,
                      int M, int N, int K,
                      int nthreads,
                      const GotoParams& p,
                      bool accum = false)
{
    const int MC = p.MC, KC = p.KC, NC = p.NC;
    const size_t Bp_sz = (size_t)KC * NC;
    const size_t Ap_sz = (size_t)MC * KC;
    IN_T* Bp      = typed_alloc<IN_T>(Bp_sz);
    IN_T* Ap_pool = typed_alloc<IN_T>(Ap_sz * nthreads);

    if (!accum) memset(C, 0, (size_t)M * N * sizeof(OUT_T));

    // Loop 5 (JC)
    for (int jc = 0; jc < N; jc += NC) {
        const int nc = std::min(NC, N - jc);

        // Loop 4 (PC)
        for (int pc = 0; pc < K; pc += KC) {
            const int kc = std::min(KC, K - pc);
            const int kc_tiles = kc / BK;

            // Pack B̃ (plain row-major tiles, NOT vnni)
            #pragma omp parallel num_threads(nthreads) proc_bind(close)
            {
                pack_B_panel_par(B, Bp, pc, jc, kc, nc, N,
                                 omp_get_thread_num(), omp_get_num_threads());
            }

            const int tn_ = nc / BN;
            const int tk_ = kc / BK;

            // Loop 3 (IC) — parallelised
            #pragma omp parallel num_threads(nthreads) proc_bind(close)
            {
                // Peak: no-op (no tile config needed for vector units)
                vec_init();

                int tid = omp_get_thread_num();
                int nth = omp_get_num_threads();
                IN_T* Ap = Ap_pool + (size_t)tid * Ap_sz;

                int ic_total = (M + MC - 1) / MC;
                int ic_base = ic_total / nth, ic_extra = ic_total % nth;
                int ic_start = tid * ic_base + std::min(tid, ic_extra);
                int ic_count = ic_base + (tid < ic_extra ? 1 : 0);

                for (int ic_idx = ic_start; ic_idx < ic_start + ic_count; ++ic_idx) {
                    int ic = ic_idx * MC;
                    int mc = std::min(MC, M - ic);
                    int tm_ = mc / BM;

                    // Pack Ã (private per thread)
                    pack_A_block(A, Ap, ic, pc, mc, kc, K);

                    // Loop 2 (JR)
                    for (int jr = 0; jr < tn_; ++jr) {
                        const IN_T* bp_base = Bp + (size_t)jr * tk_ * B_TSZ;

                        // Loop 1 (IR)
                        for (int ir = 0; ir < tm_; ++ir) {
                            const IN_T* ap_base = Ap + (size_t)ir * tk_ * A_TSZ;
                            OUT_T* C_ptr = C + (ic + ir * BM) * N + (jc + jr * BN);

                            if (accum || pc > 0) {
                                vec_microkernel_accum(ap_base, bp_base, C_ptr,
                                                      kc_tiles, N, A_TSZ, B_TSZ);
                            } else {
                                vec_microkernel(ap_base, bp_base, C_ptr,
                                                kc_tiles, N, A_TSZ, B_TSZ);
                            }
                        }
                    }
                }
            } // end IC parallel
        } // end PC
    } // end JC

    typed_free(Bp, Bp_sz);
    typed_free(Ap_pool, Ap_sz * nthreads);
}

// ═══════════════════════════════════════════════════════════════════════════
// §6b  PER-DOMAIN SPINNING BARRIER
// ═══════════════════════════════════════════════════════════════════════════
//
// Within a single flat OMP parallel region (all 96 threads), threads are
// partitioned into domains (e.g. 4 × 24).  We need barriers that only
// synchronize the threads within one domain, not all 96.
//
// OMP's #pragma omp barrier always syncs the entire team, so we use a
// simple spinning atomic barrier per domain.  On NUMA this is cheap:
// all threads in a domain are on the same CCD sharing L3, so the atomic
// stays in the local cache coherence domain.
//
struct DomainBarrier {
    alignas(64) std::atomic<int> count{0};
    alignas(64) std::atomic<int> generation{0};
    int team_size = 0;

    void init(int n) { count = 0; generation = 0; team_size = n; }

    void wait() {
        int gen = generation.load(std::memory_order_relaxed);
        if (count.fetch_add(1, std::memory_order_acq_rel) == team_size - 1) {
            // Last thread to arrive: reset counter, advance generation.
            count.store(0, std::memory_order_relaxed);
            generation.fetch_add(1, std::memory_order_release);
        } else {
            // Spin until the last thread advances the generation.
            while (generation.load(std::memory_order_acquire) == gen) {
                // Yield hint — compiler will emit PAUSE on x86, YIELD on ARM
                #if defined(__x86_64__) || defined(__i386__)
                    __builtin_ia32_pause();
                #elif defined(__aarch64__)
                    asm volatile("yield");
                #endif
            }
        }
    }
};

// ═══════════════════════════════════════════════════════════════════════════
// §7-11  NUMA GRID + STATIC / SUMMA / CANNON
// ═══════════════════════════════════════════════════════════════════════════
//
// All NUMA variants use the same approach:
//   1. ONE flat OMP parallel region with all threads (e.g. 96).
//   2. Threads are partitioned into domains: thread gtid belongs to
//      domain d = gtid / cores_per, with local ID ltid = gtid % cores_per.
//   3. Within each domain, threads collaborate on the GOTO loop body:
//      - Collaborative B̃ packing (each thread packs its share of tiles).
//      - Per-domain barrier (DomainBarrier, spinning atomic).
//      - IC-loop: each thread takes its share of MC-blocks.
//      - Per-domain barrier (before the next PC step overwrites B̃).
//   4. All domains run simultaneously → all cores active at all times.

struct NumaDomain {
    int node_id; std::vector<int> cpus; int cores;
    IN_T* A_local; size_t A_bytes;
    IN_T* B_local; size_t B_bytes;
    OUT_T* C_local; size_t C_bytes;
    int M_local, N_local;
};

struct NumaGrid {
    int PX, PY, total_threads;
    std::vector<NumaDomain> domains;

    void init(int px, int py, int total_cores) {
        PX = px; PY = py; total_threads = total_cores;
        int nd = px * py; domains.resize(nd);
        int num_nodes = get_num_numa_nodes();
        int cores_per = total_cores / nd;
        for (int d = 0; d < nd; ++d) {
            domains[d].node_id = d % num_nodes;
            domains[d].cores = cores_per;
            auto all = cpus_on_node(domains[d].node_id);
            int off = (d / num_nodes) * cores_per;
            for (int i = 0; i < cores_per && (off + i) < (int)all.size(); ++i)
                domains[d].cpus.push_back(all[off + i]);
            while ((int)domains[d].cpus.size() < cores_per && !all.empty())
                domains[d].cpus.push_back(all[domains[d].cpus.size() % all.size()]);
        }
    }
    void alloc_buffers(int M, int N, int K) {
        int Mp = M / PX, Np = N / PY;
        for (int px = 0; px < PX; ++px) for (int py = 0; py < PY; ++py) {
            int d = px * PY + py; auto& dom = domains[d];
            dom.M_local = Mp; dom.N_local = Np;
            dom.A_bytes = (size_t)Mp * K * sizeof(IN_T);
            dom.B_bytes = (size_t)K * Np * sizeof(IN_T);
            dom.C_bytes = (size_t)Mp * Np * sizeof(OUT_T);
            dom.A_local = static_cast<IN_T*>(numa_alloc(dom.A_bytes));
            dom.B_local = static_cast<IN_T*>(numa_alloc(dom.B_bytes));
            dom.C_local = static_cast<OUT_T*>(numa_alloc(dom.C_bytes));
        }
    }
    void free_buffers() {
        for (auto& dom : domains) {
            numa_free(dom.A_local, dom.A_bytes);
            numa_free(dom.B_local, dom.B_bytes);
            numa_free(dom.C_local, dom.C_bytes);
        }
    }
};

static void scatter_A_rows(const IN_T* A, IN_T* dst, int row0, int nrows, int K) {
    memcpy(dst, A + (size_t)row0 * K, (size_t)nrows * K * sizeof(IN_T));
}
static void scatter_B_cols(const IN_T* B, IN_T* dst, int col0, int ncols, int K, int N) {
    for (int k = 0; k < K; ++k)
        memcpy(dst + (size_t)k * ncols, B + (size_t)k * N + col0, ncols * sizeof(IN_T));
}
static void gather_C_block(OUT_T* C, const OUT_T* src, int r0, int c0, int nr, int nc, int N) {
    for (int i = 0; i < nr; ++i)
        memcpy(C + (size_t)(r0 + i) * N + c0, src + (size_t)i * nc, nc * sizeof(OUT_T));
}

static void numa_static_gemm(const IN_T* A, const IN_T* B, OUT_T* C,
                             int M, int N, int K, NumaGrid& grid, const GotoParams& gp) {
    const int PX = grid.PX, PY = grid.PY;
    const int Ml = M / PX, Nl = N / PY, nd = PX * PY;
    const int cores_per = grid.domains[0].cores;
    const int MC = gp.MC, KC = gp.KC, NC = gp.NC;
    const size_t Bp_sz = (size_t)KC * NC;
    const size_t Ap_sz = (size_t)MC * KC;

    // Per-domain workspace: shared B̃ + per-thread Ã pool
    std::vector<IN_T*>  Bp_bufs(nd);
    std::vector<IN_T*>  Ap_pools(nd);
    std::vector<DomainBarrier> barriers(nd);
    for (int d = 0; d < nd; ++d) {
        Bp_bufs[d]  = typed_alloc<IN_T>(Bp_sz);
        Ap_pools[d] = typed_alloc<IN_T>(Ap_sz * cores_per);
        barriers[d].init(cores_per);
    }

    // ── One flat parallel region: all threads, all domains simultaneously ──
    #pragma omp parallel num_threads(nd * cores_per) proc_bind(close)
    {
        const int gtid = omp_get_thread_num();
        const int d    = gtid / cores_per;         // which NUMA domain
        const int ltid = gtid % cores_per;         // local thread ID within domain
        auto& dom = grid.domains[d];

        // Pin this thread to its domain's CPUs for NUMA locality.
        if (ltid < (int)dom.cpus.size())
            pin_thread_to_cpus({dom.cpus[ltid]});

        vec_init();  // Peak: no-op

        // ── Phase 1: Scatter A rows / B cols to NUMA-local buffers ──────
        // Each thread copies its share of the domain's data (first-touch).
        {
            int px = d / PY, py = d % PY;
            // A rows: contiguous, partition rows across domain threads.
            size_t a_elts = (size_t)Ml * K;
            size_t a_lo = ltid * a_elts / cores_per;
            size_t a_hi = (ltid + 1) * a_elts / cores_per;
            memcpy(dom.A_local + a_lo, A + (size_t)px * Ml * K + a_lo,
                   (a_hi - a_lo) * sizeof(IN_T));

            // B cols: not contiguous in source, but contiguous in dest.
            size_t b_elts = (size_t)K * Nl;
            size_t b_lo = ltid * b_elts / cores_per;
            size_t b_hi = (ltid + 1) * b_elts / cores_per;
            // First-touch the dest pages, then copy row-by-row.
            memset(dom.B_local + b_lo, 0, (b_hi - b_lo) * sizeof(IN_T));
            // Full scatter (each thread does all K rows for its column range).
            int col_lo = ltid * Nl / cores_per;
            int col_hi = (ltid + 1) * Nl / cores_per;
            int ncols = col_hi - col_lo;
            if (ncols > 0) {
                for (int k = 0; k < K; ++k)
                    memcpy(dom.B_local + (size_t)k * Nl + col_lo,
                           B + (size_t)k * N + py * Nl + col_lo,
                           ncols * sizeof(IN_T));
            }

            // Zero local C
            size_t c_elts = (size_t)Ml * Nl;
            size_t c_lo = ltid * c_elts / cores_per;
            size_t c_hi = (ltid + 1) * c_elts / cores_per;
            memset(dom.C_local + c_lo, 0, (c_hi - c_lo) * sizeof(OUT_T));
        }
        barriers[d].wait();

        // ── Phase 2: GOTO 5-loop inline, per domain ─────────────────────
        // All cores_per threads in domain d collaborate on this domain's GEMM.
        IN_T* Bp = Bp_bufs[d];
        IN_T* Ap = Ap_pools[d] + (size_t)ltid * Ap_sz;

        // JC loop
        for (int jc = 0; jc < Nl; jc += NC) {
            const int nc = std::min(NC, Nl - jc);

            // PC loop
            for (int pc = 0; pc < K; pc += KC) {
                const int kc = std::min(KC, K - pc);
                const int kc_tiles = kc / BK;
                const int tn_ = nc / BN;
                const int tk_ = kc / BK;

                // Collaborative B̃ packing within domain.
                pack_B_panel_par(dom.B_local, Bp, pc, jc, kc, nc, Nl,
                                 ltid, cores_per);
                barriers[d].wait();  // all domain threads must see packed B̃

                // IC loop: partitioned across domain threads.
                int ic_total = (Ml + MC - 1) / MC;
                int ic_base = ic_total / cores_per, ic_extra = ic_total % cores_per;
                int ic_start = ltid * ic_base + std::min(ltid, ic_extra);
                int ic_count = ic_base + (ltid < ic_extra ? 1 : 0);

                for (int ic_idx = ic_start; ic_idx < ic_start + ic_count; ++ic_idx) {
                    int ic = ic_idx * MC;
                    int mc = std::min(MC, Ml - ic);
                    int tm_ = mc / BM;

                    // Private Ã packing (NUMA-local, in this thread's L2).
                    pack_A_block(dom.A_local, Ap, ic, pc, mc, kc, K);

                    // JR loop (over B̃ micro-panels)
                    for (int jr = 0; jr < tn_; ++jr) {
                        const IN_T* bp_base = Bp + (size_t)jr * tk_ * B_TSZ;
                        // IR loop (over Ã micro-panels)
                        for (int ir = 0; ir < tm_; ++ir) {
                            const IN_T* ap_base = Ap + (size_t)ir * tk_ * A_TSZ;
                            OUT_T* C_ptr = dom.C_local
                                         + (ic + ir * BM) * Nl + (jc + jr * BN);
                            if (pc > 0)
                                vec_microkernel_accum(ap_base, bp_base, C_ptr,
                                                      kc_tiles, Nl, A_TSZ, B_TSZ);
                            else
                                vec_microkernel(ap_base, bp_base, C_ptr,
                                                kc_tiles, Nl, A_TSZ, B_TSZ);
                        }
                    }
                }
                barriers[d].wait();  // sync before next PC step overwrites B̃
            } // end PC
        } // end JC
    } // end flat parallel

    // Gather C blocks back to global layout.
    for (int px = 0; px < PX; ++px)
        for (int py = 0; py < PY; ++py) {
            int d = px * PY + py;
            gather_C_block(C, grid.domains[d].C_local,
                           px * Ml, py * Nl, Ml, Nl, N);
        }

    for (int d = 0; d < nd; ++d) {
        typed_free(Bp_bufs[d], Bp_sz);
        typed_free(Ap_pools[d], Ap_sz * cores_per);
    }
}

static void numa_summa_gemm(const IN_T* A, const IN_T* B, OUT_T* C,
                            int M, int N, int K, NumaGrid& grid,
                            const GotoParams& gp, int KC_summa = 0) {
    const int PX = grid.PX, PY = grid.PY;
    const int Ml = M / PX, Nl = N / PY, nd = PX * PY;
    const int cores_per = grid.domains[0].cores;
    const int MC = gp.MC, KC = gp.KC, NC = gp.NC;

    if (KC_summa <= 0) KC_summa = K;
    KC_summa = ((KC_summa + BK - 1) / BK) * BK;
    if (KC_summa > K) KC_summa = K;

    // Per-domain: SUMMA panel buffers + GOTO workspace
    const size_t Bp_sz = (size_t)KC * NC;
    const size_t Ap_sz = (size_t)MC * KC;
    size_t panel_A_sz = (size_t)Ml * KC_summa;
    size_t panel_B_sz = (size_t)KC_summa * Nl;

    std::vector<IN_T*>  panel_A(nd), panel_B(nd);
    std::vector<IN_T*>  Bp_bufs(nd), Ap_pools(nd);
    std::vector<DomainBarrier> barriers(nd);

    for (int d = 0; d < nd; ++d) {
        panel_A[d]  = typed_alloc<IN_T>(panel_A_sz);
        panel_B[d]  = typed_alloc<IN_T>(panel_B_sz);
        Bp_bufs[d]  = typed_alloc<IN_T>(Bp_sz);
        Ap_pools[d] = typed_alloc<IN_T>(Ap_sz * cores_per);
        barriers[d].init(cores_per);
        memset(grid.domains[d].C_local, 0, grid.domains[d].C_bytes);
    }
    memset(C, 0, (size_t)M * N * sizeof(OUT_T));

    // ── One flat parallel region for the entire SUMMA outer loop ─────────
    #pragma omp parallel num_threads(nd * cores_per) proc_bind(close)
    {
        const int gtid = omp_get_thread_num();
        const int d    = gtid / cores_per;
        const int ltid = gtid % cores_per;
        const int px = d / PY, py = d % PY;
        auto& dom = grid.domains[d];

        if (ltid < (int)dom.cpus.size())
            pin_thread_to_cpus({dom.cpus[ltid]});
        vec_init();

        IN_T* Bp_local = Bp_bufs[d];
        IN_T* Ap_local = Ap_pools[d] + (size_t)ltid * Ap_sz;

        // ── SUMMA outer K-panel loop ────────────────────────────────
        for (int kp = 0; kp < K; kp += KC_summa) {
            int kc = std::min(KC_summa, K - kp);

            // Phase 1: Copy panels to NUMA-local buffers (parallel within domain).
            // A_panel: rows [px*Ml..(px+1)*Ml), cols [kp..kp+kc)
            {
                int row_lo = ltid * Ml / cores_per;
                int row_hi = (ltid + 1) * Ml / cores_per;
                for (int i = row_lo; i < row_hi; ++i)
                    memcpy(panel_A[d] + (size_t)i * kc,
                           A + (size_t)(px * Ml + i) * K + kp,
                           kc * sizeof(IN_T));
            }
            // B_panel: rows [kp..kp+kc), cols [py*Nl..(py+1)*Nl)
            {
                int k_lo = ltid * kc / cores_per;
                int k_hi = (ltid + 1) * kc / cores_per;
                for (int k = k_lo; k < k_hi; ++k)
                    memcpy(panel_B[d] + (size_t)k * Nl,
                           B + (size_t)(kp + k) * N + py * Nl,
                           Nl * sizeof(IN_T));
            }
            barriers[d].wait();

            // Phase 2: Local GOTO on the panel (inline, all domain threads).
            // C_local += panel_A[Ml × kc] · panel_B[kc × Nl]
            bool accum = (kp > 0);
            int local_KC = std::min(KC, kc);
            local_KC = ((local_KC + BK - 1) / BK) * BK;

            for (int jc = 0; jc < Nl; jc += NC) {
                int nc = std::min(NC, Nl - jc);
                for (int pc = 0; pc < kc; pc += local_KC) {
                    int lkc = std::min(local_KC, kc - pc);
                    int kc_tiles = lkc / BK;
                    int tn_ = nc / BN, tk_ = lkc / BK;

                    // Collaborative B̃ packing within domain
                    pack_B_panel_par(panel_B[d], Bp_local, pc, jc, lkc, nc, Nl,
                                     ltid, cores_per);
                    barriers[d].wait();

                    // IC loop: partitioned across domain threads
                    int ic_total = (Ml + MC - 1) / MC;
                    int ic_base = ic_total / cores_per, ic_extra = ic_total % cores_per;
                    int ic_start = ltid * ic_base + std::min(ltid, ic_extra);
                    int ic_count = ic_base + (ltid < ic_extra ? 1 : 0);

                    for (int ic_idx = ic_start; ic_idx < ic_start + ic_count; ++ic_idx) {
                        int ic = ic_idx * MC;
                        int mc = std::min(MC, Ml - ic);
                        int tm_ = mc / BM;

                        pack_A_block(panel_A[d], Ap_local, ic, pc, mc, lkc, kc);

                        for (int jr = 0; jr < tn_; ++jr) {
                            const IN_T* bp_base = Bp_local + (size_t)jr * tk_ * B_TSZ;
                            for (int ir = 0; ir < tm_; ++ir) {
                                const IN_T* ap_base = Ap_local + (size_t)ir * tk_ * A_TSZ;
                                OUT_T* C_ptr = dom.C_local
                                             + (ic + ir * BM) * Nl + (jc + jr * BN);
                                if (accum || pc > 0)
                                    vec_microkernel_accum(ap_base, bp_base, C_ptr,
                                                          kc_tiles, Nl, A_TSZ, B_TSZ);
                                else
                                    vec_microkernel(ap_base, bp_base, C_ptr,
                                                    kc_tiles, Nl, A_TSZ, B_TSZ);
                            }
                        }
                    }
                    barriers[d].wait();
                } // end PC (local)
            } // end JC (local)
        } // end SUMMA K-panel loop
    } // end flat parallel

    // Gather C
    for (int px = 0; px < PX; ++px)
        for (int py = 0; py < PY; ++py) {
            int d = px * PY + py;
            gather_C_block(C, grid.domains[d].C_local,
                           px * Ml, py * Nl, Ml, Nl, N);
        }

    for (int d = 0; d < nd; ++d) {
        typed_free(panel_A[d], panel_A_sz);
        typed_free(panel_B[d], panel_B_sz);
        typed_free(Bp_bufs[d], Bp_sz);
        typed_free(Ap_pools[d], Ap_sz * cores_per);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §12  MKL BASELINE
// ═══════════════════════════════════════════════════════════════════════════
//
// oneMKL sgemm after upcasting fp16 → fp32.  On Zen4 MKL will use its
// AVX-512 codepath (not AMX, since there is no AMX on AMD).  This gives
// the vendor-optimised baseline for the same vector ISA we use.
//
// Allocates temporary fp32 buffers on demand (512 MB for 8192²).  These
// are freed after each call so they don't bloat the persistent footprint.

static void mkl_sgemm_f32(const IN_T* A, const IN_T* B, OUT_T* C, int M, int N, int K) {
    mkl_set_dynamic(0); mkl_set_num_threads(g_threads);
    float* A32 = typed_alloc<float>((size_t)M * K);
    float* B32 = typed_alloc<float>((size_t)K * N);
    #pragma omp parallel for num_threads(g_threads) schedule(static)
    for (size_t i = 0; i < (size_t)M*K; ++i) A32[i] = (float)A[i];
    #pragma omp parallel for num_threads(g_threads) schedule(static)
    for (size_t i = 0; i < (size_t)K*N; ++i) B32[i] = (float)B[i];
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0f, A32, K, B32, N, 0.0f, C, N);
    typed_free(A32, (size_t)M * K);
    typed_free(B32, (size_t)K * N);
}

// ═══════════════════════════════════════════════════════════════════════════
// §13-14  VERIFICATION + BENCHMARK
// ═══════════════════════════════════════════════════════════════════════════

static bool verify(const OUT_T* C, const OUT_T* ref, int M, int N,
                   const char* label, float tol = 0.02f) {
    int bad = 0; double max_abs = 0, max_rel = 0, sum_abs = 0;
    int first_bad[3] = {-1, -1, -1}; size_t total = (size_t)M * N;
    for (size_t i = 0; i < total; ++i) {
        double d = fabs((double)C[i] - (double)ref[i]);
        double den = fabs((double)ref[i]);
        double rel = (den > 1e-6) ? d / den : d;
        sum_abs += d; if (d > max_abs) max_abs = d; if (rel > max_rel) max_rel = rel;
        if (rel > tol && d > 0.05f) { if (bad < 3) first_bad[bad] = (int)i; ++bad; }
    }
    if (bad == 0)
        printf("    [PASS] %-24s max_abs=%.6f max_rel=%.4f%% mean_abs=%.6f\n",
               label, max_abs, max_rel*100, sum_abs/total);
    else {
        printf("    [FAIL] %-24s %d mismatches max_abs=%.6f max_rel=%.4f%%\n",
               label, bad, max_abs, max_rel*100);
        for (int k = 0; k < 3 && first_bad[k] >= 0; ++k) {
            int idx = first_bad[k];
            printf("           idx=%d (%d,%d) got=%.6f ref=%.6f\n",
                   idx, idx/N, idx%N, C[idx], ref[idx]);
        }
    }
    return bad == 0;
}

using KernFn = std::function<void(const IN_T*, const IN_T*, OUT_T*, int, int, int)>;

static void bench_one(FILE* f, const char* name, KernFn fn,
                      const IN_T* A, const IN_T* B, OUT_T* C,
                      int M, int N, int K, int nthreads) {
    for (int w = 0; w < g_warmup; ++w) fn(A, B, C, M, N, K);
    double flops = 2.0 * (double)M * N * K;
    std::vector<double> times(g_iters); double total = 0;
    for (int i = 0; i < g_iters; ++i) {
        double t0 = now_ns(); fn(A, B, C, M, N, K); double t1 = now_ns();
        times[i] = t1 - t0; total += times[i];
    }
    double avg = total / g_iters;
    std::sort(times.begin(), times.end());
    double med = (g_iters % 2) ? times[g_iters/2] : (times[g_iters/2-1]+times[g_iters/2])/2;
    for (int i = 0; i < g_iters; ++i)
        fprintf(f, "%s,%d,%d,%d,%d,%d,%.2f,%.6f\n", name, M, N, K, nthreads, i, times[i], flops/times[i]);
    printf("    %-28s t=%2d avg=%.0f ns  med=%.0f ns  %.4f GF/s (peak %.4f)\n",
           name, nthreads, avg, med, flops/avg, flops/times[0]);
}

// ═══════════════════════════════════════════════════════════════════════════
// §15-16  ARG PARSING + MAIN
// ═══════════════════════════════════════════════════════════════════════════

struct MNK { int M, N, K; };

static void usage(const char* p) {
    printf("Usage: %s M N K [...] [-w W] [-i I] [-t T] [-p PXxPY] [-mc MC] [-kc KC] [-nc NC] [-o F]\n", p);
}

static std::vector<MNK> parse_args(int argc, char** argv) {
    std::vector<MNK> sizes; std::vector<int> nums;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-w") && i+1<argc) g_warmup = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-i") && i+1<argc) g_iters = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t") && i+1<argc) g_threads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-p") && i+1<argc) { ++i; sscanf(argv[i], "%dx%d", &g_PX, &g_PY); }
        else if (!strcmp(argv[i], "-mc") && i+1<argc) g_MC = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-kc") && i+1<argc) g_KC = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nc") && i+1<argc) g_NC = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1<argc) g_outfile = argv[++i];
        else if (!strcmp(argv[i], "-h")) { usage(argv[0]); exit(0); }
        else nums.push_back(atoi(argv[i]));
    }
    if (nums.size() % 3 || nums.empty()) { fprintf(stderr, "Need M N K triples\n"); exit(1); }
    for (size_t i = 0; i < nums.size(); i += 3) sizes.push_back({nums[i], nums[i+1], nums[i+2]});
    return sizes;
}

int main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 1; }
    auto sizes = parse_args(argc, argv);
    GotoParams gp{g_MC, g_KC, g_NC};
    int nd = g_PX * g_PY;

    printf("════════════════════════════════════════════════════════\n");
    printf("  GOTO + NUMA GEMM Benchmark  (Auto-Vectorized Version)\n");
    printf("════════════════════════════════════════════════════════\n");
    printf("  Microkernel   : auto-vectorized C (#pragma omp simd)\n");
    printf("  Threads       : %d\n", g_threads);
    printf("  NUMA nodes    : %d (detected)\n", get_num_numa_nodes());
    printf("  Domain grid   : %d x %d = %d\n", g_PX, g_PY, nd);
    printf("  Cores/domain  : %d\n", g_threads / nd);
    printf("  GOTO blocking : MC=%d KC=%d NC=%d\n", gp.MC, gp.KC, gp.NC);
    printf("  Tile size     : %dx%dx%d (same as AMX version)\n", BM, BN, BK);
    printf("  B packing     : plain row-major (NOT vnni)\n");
    printf("  Peak (init)   : no-op\n");
    printf("  Unpack (store): no-op (D is row-major fp32 array)\n");
    printf("════════════════════════════════════════════════════════\n\n");

    FILE* fout = fopen(g_outfile, "w");
    if (!fout) { perror("fopen"); return 1; }
    fprintf(fout, "kernel,M,N,K,threads,run,time_ns,gflops\n");

    NumaGrid grid;
    grid.init(g_PX, g_PY, g_threads);

    for (auto& sz : sizes) {
        int M = sz.M, N = sz.N, K = sz.K;
        printf("═══ M=%d  N=%d  K=%d ═══════════════════════════════\n", M, N, K);
        if (M % (g_PX * BM) || N % (g_PY * BN) || K % BK) {
            printf("  [SKIP] divisibility\n"); continue;
        }

        size_t szA = (size_t)M*K, szB = (size_t)K*N, szC = (size_t)M*N;
        double mem_mb = (szA * sizeof(IN_T) + szB * sizeof(IN_T)
                        + szC * sizeof(OUT_T) * 2) / (1024.0 * 1024.0);
        printf("  Global alloc: A=%.0f MB  B=%.0f MB  C×2=%.0f MB  total=%.0f MB\n",
               szA * sizeof(IN_T) / 1e6, szB * sizeof(IN_T) / 1e6,
               szC * sizeof(OUT_T) * 2 / 1e6, mem_mb);

        IN_T*  A    = typed_alloc<IN_T>(szA);
        IN_T*  B    = typed_alloc<IN_T>(szB);
        OUT_T* Cref = typed_alloc<OUT_T>(szC);
        OUT_T* Cchk = typed_alloc<OUT_T>(szC);

        #pragma omp parallel num_threads(g_threads) proc_bind(close)
        {
            int t = omp_get_thread_num(), n_ = omp_get_num_threads();
            unsigned seed = 42 + t; size_t lo, hi;
            lo = t*szA/n_; hi = (t+1)*szA/n_;
            for (size_t i = lo; i < hi; ++i) A[i] = (IN_T)((float)rand_r(&seed)/RAND_MAX - 0.5f);
            lo = t*szB/n_; hi = (t+1)*szB/n_;
            for (size_t i = lo; i < hi; ++i) B[i] = (IN_T)((float)rand_r(&seed)/RAND_MAX - 0.5f);
        }

        grid.alloc_buffers(M, N, K);

        // ── Reference: MKL sgemm (fast, vendor-optimised) ───────────
        // The scalar OpenMP ref is unusable at 8192³ (takes minutes).
        // We use MKL as the reference and verify all other kernels against it.
        printf("  Reference (MKL sgemm)...\n");
        mkl_sgemm_f32(A, B, Cref, M, N, K);

        printf("  Kernels:\n");

        // 0. MKL baseline (same as reference — benchmark timing only)
        bench_one(fout, "mkl_sgemm_f32", mkl_sgemm_f32, A, B, Cchk, M, N, K, g_threads);
        verify(Cchk, Cref, M, N, "mkl_sgemm_f32");

        // 2. GOTO GEMM (vector microkernel)
        auto goto_fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
            goto_gemm(a, b, c, m, n, k, g_threads, gp);
        };
        bench_one(fout, "goto_vec", goto_fn, A, B, Cchk, M, N, K, g_threads);
        verify(Cchk, Cref, M, N, "goto_vec");

        // 3. NUMA-Static
        {
            char lbl[64]; snprintf(lbl, sizeof(lbl), "numa_static_%dx%d", g_PX, g_PY);
            auto fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
                numa_static_gemm(a, b, c, m, n, k, grid, gp);
            };
            bench_one(fout, lbl, fn, A, B, Cchk, M, N, K, g_threads);
            verify(Cchk, Cref, M, N, lbl);
        }

        // 4. NUMA-SUMMA
        {
            char lbl[64]; snprintf(lbl, sizeof(lbl), "numa_summa_%dx%d_kc%d", g_PX, g_PY, gp.KC);
            auto fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
                numa_summa_gemm(a, b, c, m, n, k, grid, gp, gp.KC);
            };
            bench_one(fout, lbl, fn, A, B, Cchk, M, N, K, g_threads);
            verify(Cchk, Cref, M, N, lbl);
        }

        grid.free_buffers();
        typed_free(A, szA); typed_free(B, szB);
        typed_free(Cref, szC); typed_free(Cchk, szC);
        printf("\n");
    }

    fclose(fout);
    printf("All runs written to %s\n", g_outfile);
    return 0;
}