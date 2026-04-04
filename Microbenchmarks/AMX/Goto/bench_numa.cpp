// bench_goto_numa.cpp
// ═══════════════════════════════════════════════════════════════════════════
// GOTO GEMM with AMX microkernel  +  SUMMA / Cannon over NUMA domains
// ═══════════════════════════════════════════════════════════════════════════
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  BACKGROUND: THE GOTO ALGORITHM                                        │
// │  (Goto & van de Geijn, "Anatomy of High-Performance Matrix             │
// │   Multiplication", ACM TOMS 2008)                                      │
// │                                                                        │
// │  The key insight is to decompose C = A·B into a hierarchy of blocked   │
// │  operations that map onto the CPU cache hierarchy:                     │
// │                                                                        │
// │    C[M×N] += A[M×K] · B[K×N]                                          │
// │                                                                        │
// │  Five nested loops, from outermost to innermost:                       │
// │                                                                        │
// │    Loop 5 (JC): partition N into panels of width NC                    │
// │      Loop 4 (PC): partition K into panels of depth KC (rank-KC update) │
// │        → Pack B̃ : KC×NC panel of B into contiguous L3-resident buffer │
// │        Loop 3 (IC): partition M into blocks of height MC               │
// │          → Pack Ã : MC×KC block of A into contiguous L2-resident buf   │
// │          Loop 2 (JR): sweep B̃ in micro-panels of width NR             │
// │            Loop 1 (IR): sweep Ã in micro-panels of height MR           │
// │              → Microkernel: MR×NR tile of C, rank-KC update            │
// │                                                                        │
// │  Cache mapping (Sapphire Rapids fp16):                                 │
// │    • KC × NR × 2B ≈ L1 capacity  (streaming B̃ micro-panel)           │
// │    • MC × KC × 2B ≈ ½ L2         (Ã block resident in L2)            │
// │    • KC × NC × 2B ≤ L3           (B̃ panel shared in L3)              │
// │                                                                        │
// │  Data reuse logic:                                                     │
// │    The JR loop iterates over B̃ micro-panels while Ã stays in L2.     │
// │    The IR loop iterates over Ã micro-panels while one B̃ micro-panel  │
// │    stays in L1.  This maximises data reuse at every cache level.       │
// │                                                                        │
// │  Packing serves three purposes:                                        │
// │    1. Eliminates TLB misses (packed buffer is contiguous)              │
// │    2. Aligns data for SIMD / AMX tile loads                            │
// │    3. Reorders data to match the microkernel's access pattern          │
// │                                                                        │
// │  Threading (BLIS IC-parallel strategy):                                │
// │    We parallelise the IC loop (over M blocks).  All IC threads share  │
// │    the same packed B̃ (in L3), but each packs its own Ã (in L2).      │
// │                                                                        │
// │  For multi-socket / multi-NUMA, we layer SUMMA or Cannon on top:       │
// │    • Each NUMA domain runs an independent GOTO GEMM on local data.     │
// │    • Data is pre-distributed (Static) or iteratively communicated      │
// │      (SUMMA panels / Cannon shifts).                                   │
// │    • This eliminates cross-NUMA traffic in the hot inner loops.        │
// │    • Analogy: JC-parallel in BLIS replicates B̃ per socket; our       │
// │      SUMMA extends this to a full 2D decomposition of C.               │
// └─────────────────────────────────────────────────────────────────────────┘
//
// Build:
//   g++ -O3 -fopenmp -march=sapphirerapids -mamx-tile -mamx-bf16 -mamx-int8 \
//       -I. -I${MKLROOT}/include -o bench_goto_numa bench_goto_numa.cpp       \
//       -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread           \
//       -lmkl_core -lgomp -lpthread -lm
//
// Run examples:
//   ./bench_goto_numa 4096 4096 4096 -t 32 -p 2x2        # 4 NUMA domains
//   ./bench_goto_numa 4096 4096 4096 -t 24 -p 3x1        # 3 NUMA domains
//   ./bench_goto_numa 4096 4096 4096 -t 48 -p 3x2        # 6 NUMA domains
//   ./bench_goto_numa 4096 4096 4096 -t 32 -p 2x2 -mc 512 -kc 128
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
#include <numeric>
#include <string>
#include <vector>

#include "core/amx/tile.hpp"

using namespace core_ir;

// ═══════════════════════════════════════════════════════════════════════════
// §1  COMPILE-TIME CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════
//
// AMX hardware tile dimensions.  These are fixed by the tile.hpp abstraction
// which maps the raw AMX 16-row tiles into a 32×32 logical tile:
//
//   FragmentA = 32×32 fp16  → 2 AMX tile registers  (2 × 16×32)
//   FragmentB = 32×32 fp16  → 2 AMX tile registers  (vnni packed)
//   FragmentD = 32×32 fp32  → 4 AMX tile registers  (4 × 16×16)
//                              ─────────────────────────────────
//                              Total = 8 AMX tile registers (all of them)
//
// This means one 32×32 output tile saturates the register file.
// In GOTO terminology:  MR = BM = 32,  NR = BN = 32.
// There is no room for register blocking beyond 1×1 output tiles
// without spilling tile registers to the stack.
//
static constexpr int BM = 32;          // tile rows   (= MR in GOTO)
static constexpr int BN = 32;          // tile cols   (= NR in GOTO)
static constexpr int BK = 32;          // tile depth  (K-step per single AMX MMA)
static constexpr int A_TSZ = BM * BK;  // elements per packed-A tile (32×32 = 1024 fp16)
static constexpr int B_TSZ = (BK / 2) * (2 * BN);  // elements per vnni-packed B tile

// ── GOTO blocking parameter defaults ────────────────────────────────────
//
// These control how the three outer loops (JC, PC, IC) tile the problem
// to fit different levels of the cache hierarchy on Sapphire Rapids.
//
//   MC (IC block height):
//     Ã occupies MC × KC × sizeof(fp16) bytes.
//     Target: fit in ≈ ½ L2.
//     SPR L2 = 2 MB → 1024 × 256 × 2B = 512 KB.  Good headroom.
//
//   KC (PC rank-update depth):
//     One NR-wide micro-panel of B̃ streams from L1 during the microkernel.
//     KC × NR × sizeof(fp16) should fit comfortably in L1.
//     SPR L1 = 48 KB → 256 × 32 × 2B = 16 KB.  Plenty of room.
//     KC also affects TLB: MC × KC × 2B / 4KB = 128 pages for defaults.
//
//   NC (JC panel width):
//     B̃ occupies KC × NC × sizeof(fp16) bytes, shared in L3.
//     SPR L3 share ≈ 2 MB/core → 256 × 4096 × 2B = 2 MB.
//
using IN_T  = _Float16;
using OUT_T = float;
using GEMM  = TileOps<IN_T, IN_T, OUT_T, OUT_T,
                      MajorAxis::ROW, MajorAxis::ROW,
                      MajorAxis::ROW, MajorAxis::ROW,
                      BM, BN, BK>;

// ═══════════════════════════════════════════════════════════════════════
// §2  GLOBALS & PARAMETERS
// ═══════════════════════════════════════════════════════════════════════

static int  g_warmup  = 5;
static int  g_iters   = 50;
static int  g_threads = 32;
static int  g_PX      = 2;       // NUMA grid rows
static int  g_PY      = 2;       // NUMA grid cols
static int  g_MC      = DEFAULT_MC;
static int  g_KC      = DEFAULT_KC;
static int  g_NC      = DEFAULT_NC;
static const char* g_outfile = "bench_goto_numa.csv";

static inline double now_ns() {
    auto t = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::nano>(t.time_since_epoch()).count();
}

// ═══════════════════════════════════════════════════════════════════════════
// §3  NUMA UTILITIES
// ═══════════════════════════════════════════════════════════════════════════
//
// All buffers use anonymous mmap + MADV_HUGEPAGE.  This gives us:
//   • 2 MB transparent huge pages (THP), reducing TLB pressure on the
//     large packed A/B buffers (hundreds of KB to tens of MB).
//   • Pages are MAP_PRIVATE | MAP_ANONYMOUS → *migrateable* by the kernel:
//     move_pages() can relocate them between NUMA nodes, and Linux's
//     automatic NUMA balancing can migrate frequently-accessed remote pages.
//   • NOT locked (no MAP_LOCKED / mlock), so the kernel can also swap or
//     compact them.  For GEMM workloads the working set fits in DRAM.
//   • First-touch policy places pages on the NUMA node of the first
//     writing thread.  We pin threads before touching → deterministic
//     NUMA placement without explicit mbind().
//
static int get_numa_node() {
    unsigned cpu, node;
    syscall(__NR_getcpu, &cpu, &node, nullptr);
    return (int)node;
}

static int get_num_numa_nodes() {
    // count directories /sys/devices/system/node/node*
    int n = 0;
    char path[256];
    for (int i = 0; i < 64; ++i) {
        snprintf(path, sizeof(path), "/sys/devices/system/node/node%d", i);
        if (access(path, F_OK) == 0) ++n;
        else break;
    }
    return n > 0 ? n : 1;
}

// return sorted list of online CPU ids on `node`
static std::vector<int> cpus_on_node(int node) {
    std::vector<int> cpus;
    char path[256];
    snprintf(path, sizeof(path),
             "/sys/devices/system/node/node%d/cpulist", node);
    FILE* f = fopen(path, "r");
    if (!f) return cpus;
    char buf[4096];
    if (!fgets(buf, sizeof(buf), f)) { fclose(f); return cpus; }
    fclose(f);
    // parse comma-separated ranges: "0-7,16-23"
    char* p = buf;
    while (*p) {
        int lo = (int)strtol(p, &p, 10);
        int hi = lo;
        if (*p == '-') { ++p; hi = (int)strtol(p, &p, 10); }
        for (int c = lo; c <= hi; ++c) cpus.push_back(c);
        if (*p == ',') ++p;
    }
    std::sort(cpus.begin(), cpus.end());
    return cpus;
}

// hugepage-backed anonymous mmap
static void* numa_alloc(size_t bytes) {
    void* p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return p;
}

static void numa_free(void* p, size_t bytes) {
    if (p && bytes) munmap(p, bytes);
}

template <typename T>
static T* typed_alloc(size_t n) {
    return static_cast<T*>(numa_alloc(n * sizeof(T)));
}

template <typename T>
static void typed_free(T* p, size_t n) {
    numa_free(p, n * sizeof(T));
}

// first-touch: have thread `tid` (pinned to a specific NUMA node) write zeros
static void first_touch_zero(void* p, size_t bytes, int nthreads) {
    char* base = (char*)p;
    #pragma omp parallel num_threads(nthreads) proc_bind(close)
    {
        int t = omp_get_thread_num(), n = omp_get_num_threads();
        size_t lo = t * bytes / n;
        size_t hi = (t + 1) * bytes / n;
        memset(base + lo, 0, hi - lo);
    }
}

// bind calling thread to a specific CPU set
static void pin_thread_to_cpus(const std::vector<int>& cpus) {
    cpu_set_t set;
    CPU_ZERO(&set);
    for (int c : cpus) CPU_SET(c, &set);
    sched_setaffinity(0, sizeof(set), &set);
}

// ═══════════════════════════════════════════════════════════════════════════
// §3b  OPENMP REFERENCE GEMM (scalar, trivially correct)
// ═══════════════════════════════════════════════════════════════════════════
//
// A simple 2D block-tiled scalar GEMM: fp16 inputs, fp32 output, row-major.
// No packing, no SIMD, no AMX — just nested loops with OpenMP collapse(2).
// This is the numerical ground truth: easy to inspect, impossible to get
// wrong, and independent of any library (MKL could silently use AMX or
// different rounding).
//
// The 64×64 blocking is for cache-friendliness, not correctness.
// The ikj loop order within each block maximises A-row reuse.
//
static void ref_gemm_omp(const IN_T* __restrict__ A,
                         const IN_T* __restrict__ B,
                         OUT_T* __restrict__ C,
                         int M, int N, int K)
{
    memset(C, 0, (size_t)M * N * sizeof(OUT_T));
    constexpr int BLK = 64;

    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int bi = 0; bi < M; bi += BLK) {
        for (int bj = 0; bj < N; bj += BLK) {
            // Iterate over K blocks outside the collapse to avoid races on C.
            // Each (bi,bj) block is owned by exactly one thread → safe.
            for (int bk = 0; bk < K; bk += BLK) {
                int mi = std::min(bi + BLK, M);
                int nj = std::min(bj + BLK, N);
                int pk = std::min(bk + BLK, K);
                // ikj order: load A[i,k] once, broadcast across j
                for (int i = bi; i < mi; ++i) {
                    for (int k = bk; k < pk; ++k) {
                        float a_val = (float)A[i * K + k];
                        for (int j = bj; j < nj; ++j) {
                            C[i * N + j] += a_val * (float)B[k * N + j];
                        }
                    }
                }
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §4  PACKING ROUTINES (AMX format)
// ═══════════════════════════════════════════════════════════════════════════
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  WHY PACKING?                                                          │
// │                                                                        │
// │  In the original matrix, elements of a BM×BK tile of A are spread     │
// │  across BM rows with stride K (the full matrix width).  For K=4096,   │
// │  consecutive rows of a tile are 8 KB apart (fp16).  Loading a 32-row  │
// │  tile touches 32 separate cache lines / pages — bad for TLB and       │
// │  prefetcher.                                                           │
// │                                                                        │
// │  Packing copies each tile into a contiguous buffer:                    │
// │    • All BM×BK = 1024 elements of one A tile occupy 2 KB contiguous.  │
// │    • AMX _tile_loadd loads the tile with stride BK*2B = 64B = 1 cache │
// │      line → perfect streaming, no wasted prefetch bandwidth.           │
// │    • The full MC×KC packed Ã block = (MC/BM)*(KC/BK) tiles × 2 KB    │
// │      e.g. 32 × 8 × 2 KB = 512 KB → fits in L2.                       │
// │                                                                        │
// │  B packing is more complex: AMX uses "VNNI format" where pairs of     │
// │  consecutive K-values are interleaved with N-values:                   │
// │    packed[(k/2)*(2*BN) + 2*n + 0] = B[k,   n]                         │
// │    packed[(k/2)*(2*BN) + 2*n + 1] = B[k+1, n]                         │
// │  This layout lets the AMX hardware load two fp16 values from adjacent  │
// │  K indices and multiply them in a single fused operation.              │
// │                                                                        │
// │  Packing cost is amortised:                                            │
// │    Pack Ã: O(MC·KC) work, reused across NC/NR = NC/32 JR iterations. │
// │    Pack B̃: O(KC·NC) work, reused across MC/MR = MC/32 IR iterations  │
// │             and across all M/MC IC blocks.                             │
// └─────────────────────────────────────────────────────────────────────────┘
//

// ── Pack Ã block ────────────────────────────────────────────────────────
//
// Copies A[row0..row0+mc, col0..col0+kc] into contiguous tile layout.
//
// Packed layout for tile index (tm, tk):
//   dst[(tm * tiles_k + tk) * A_TSZ + li * BK + lk]
//     = A[(row0 + tm*BM + li) * lda + (col0 + tk*BK + lk)]
//
// Tiles stored in row-major tile order; within each tile, row-major with
// stride BK.  This is exactly what FragmentA::load() expects.
static void pack_A_block(const IN_T* __restrict__ A, IN_T* __restrict__ dst,
                         int row0, int col0, int mc, int kc, int lda) {
    const int tm_ = mc / BM;
    const int tk_ = kc / BK;
    for (int tm = 0; tm < tm_; ++tm) {
        for (int tk = 0; tk < tk_; ++tk) {
            IN_T* __restrict__ d = dst + (size_t)(tm * tk_ + tk) * A_TSZ;
            const int gr = row0 + tm * BM;
            const int gc = col0 + tk * BK;
            for (int li = 0; li < BM; ++li) {
                for (int lk = 0; lk < BK; ++lk) {
                    d[li * BK + lk] = A[(gr + li) * lda + (gc + lk)];
                }
            }
        }
    }
}

// ── Pack B̃ panel (VNNI format) ──────────────────────────────────────────
//
// Copies B[row0..row0+kc, col0..col0+nc] into VNNI-interleaved tiles.
// B is row-major (K rows × N cols).
//
// Tile order: [tn][tk] (N-tiles outer, K-tiles inner) so that all K-tiles
// for the same NR-wide micro-panel are contiguous.  The JR loop indexes
// B̃ by advancing tn, jumping tk_ * B_TSZ elements to the next micro-panel.
//
// VNNI interleaving within each tile:
//   for k ∈ [0, BK) step 2, n ∈ [0, BN):
//     dst[(k/2)*(2*BN) + 2*n + 0] = B[row0+tk*BK+k,   col0+tn*BN+n]
//     dst[(k/2)*(2*BN) + 2*n + 1] = B[row0+tk*BK+k+1, col0+tn*BN+n]
static void pack_B_panel(const IN_T* __restrict__ B, IN_T* __restrict__ dst,
                         int row0, int col0, int kc, int nc, int ldb) {
    const int tk_ = kc / BK;
    const int tn_ = nc / BN;
    for (int tn = 0; tn < tn_; ++tn) {
        for (int tk = 0; tk < tk_; ++tk) {
            IN_T* __restrict__ d = dst + (size_t)(tn * tk_ + tk) * B_TSZ;
            const int gr = row0 + tk * BK;
            const int gc = col0 + tn * BN;
            for (int k = 0; k < BK; k += 2) {
                for (int n = 0; n < BN; ++n) {
                    d[(k / 2) * (2 * BN) + 2 * n + 0] = B[(gr + k)     * ldb + (gc + n)];
                    d[(k / 2) * (2 * BN) + 2 * n + 1] = B[(gr + k + 1) * ldb + (gc + n)];
                }
            }
        }
    }
}

// ── Parallel packing helpers ────────────────────────────────────────────
//
// In the GOTO algorithm, B̃ packing happens once per (JC, PC) iteration
// and is shared by all IC threads.  We parallelise it across all threads
// so each touches a disjoint subset of tiles.
//
// NUMA note: first-touch on the packed buffer places pages on the node of
// the writing thread.  For B̃ (shared in L3), this is acceptable when all
// IC threads are on the same socket.  For multi-socket, the NUMA-SUMMA
// wrapper ensures each domain packs its own B̃ into NUMA-local memory
// (analogous to BLIS's JC-parallel B̃ replication per socket).
static void pack_A_block_par(const IN_T* A, IN_T* dst,
                             int row0, int col0, int mc, int kc, int lda,
                             int tid, int nth) {
    const int tm_ = mc / BM;
    const int tk_ = kc / BK;
    const int total = tm_ * tk_;
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
    const int tk_ = kc / BK;
    const int tn_ = nc / BN;
    const int total = tn_ * tk_;
    int base = total / nth, extra = total % nth;
    int start = tid * base + std::min(tid, extra);
    int count = base + (tid < extra ? 1 : 0);
    for (int idx = start; idx < start + count; ++idx) {
        int tn = idx / tk_, tk = idx % tk_;
        IN_T* d = dst + (size_t)(tn * tk_ + tk) * B_TSZ;
        int gr = row0 + tk * BK, gc = col0 + tn * BN;
        for (int k = 0; k < BK; k += 2)
            for (int n = 0; n < BN; ++n) {
                d[(k / 2) * (2 * BN) + 2 * n + 0] = B[(gr + k)     * ldb + (gc + n)];
                d[(k / 2) * (2 * BN) + 2 * n + 1] = B[(gr + k + 1) * ldb + (gc + n)];
            }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §5  AMX MICROKERNEL
// ═══════════════════════════════════════════════════════════════════════════
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  THE MICROKERNEL IN GOTO / BLIS                                        │
// │                                                                        │
// │  The microkernel is the innermost computation unit.  It computes:      │
// │                                                                        │
// │    C_tile[MR × NR] += Ã_micropanel[MR × KC] · B̃_micropanel[KC × NR] │
// │                                                                        │
// │  In our case MR = NR = 32 (one AMX logical tile = 32×32).             │
// │  The KC dimension is processed in BK=32 steps; each step is one AMX   │
// │  MMA instruction (tdpfp16ps), performing a rank-32 update of the      │
// │  fp32 accumulator.                                                     │
// │                                                                        │
// │  AMX register budget per invocation:                                   │
// │    tmm0-tmm3 : FragmentD accumulator  (4 tiles = 32×32 fp32)          │
// │    tmm4-tmm5 : FragmentA operand      (2 tiles = 32×32 fp16)          │
// │    tmm6-tmm7 : FragmentB operand      (2 tiles = vnni-packed fp16)    │
// │    ─────────────────────────────────────────────────────────────       │
// │    Total: 8/8 tile registers.  Fully saturated.                        │
// │                                                                        │
// │  Data flow per inner iteration (one K-tile step):                      │
// │    1. _tile_loadd(tmm4-5, Ã + tk*A_TSZ)  → load 32×32 fp16 A tile    │
// │    2. _tile_loadd(tmm6-7, B̃ + tk*B_TSZ)  → load 32×32 vnni B tile   │
// │    3. _tdpfp16ps(tmm0-3, tmm4-5, tmm6-7) → rank-32 FMA into D       │
// │    4. Advance ap += A_TSZ, bp += B_TSZ                                │
// │                                                                        │
// │  After KC/BK iterations, store D to C:                                │
// │    5. _tile_stored(tmm0-3, C_ptr, ldc*sizeof(float))                  │
// └─────────────────────────────────────────────────────────────────────────┘
//

// Write-new variant: D starts at zero, result overwrites C.
static inline void amx_microkernel(const IN_T* __restrict__ ap,
                                   const IN_T* __restrict__ bp,
                                   OUT_T* __restrict__ C_ptr,
                                   int kc_tiles, int ldc,
                                   int a_stride,  // stride between A k-tiles
                                   int b_stride)  // stride between B k-tiles
{
    // Zero-initialise the 32×32 fp32 accumulator (lives in 4 AMX tile regs).
    typename GEMM::FragmentD df{};

    // Inner loop: KC/BK rank-32 updates.
    // Each iteration loads one A tile and one B tile, then does a fused
    // multiply-accumulate (MMA) adding rank-32 outer product to df.
    for (int tk = 0; tk < kc_tiles; ++tk) {
        typename GEMM::FragmentA af;
        typename GEMM::FragmentB bf;
        af.template load(ap);        // _tile_loadd from contiguous packed Ã
        bf.template load(bp);        // _tile_loadd from contiguous packed B̃ (vnni)
        GEMM::MMA::mma(af, bf, df);  // _tdpfp16ps: D += A × B (rank-32 update)
        ap += a_stride;              // advance to next K-tile of Ã (= A_TSZ elts)
        bp += b_stride;              // advance to next K-tile of B̃ (= B_TSZ elts)
    }

    // Store the 32×32 fp32 result to C (row-major, strided by ldc = N).
    df.template unload_unpack<MajorAxis::ROW>(C_ptr, ldc);
}

// Accumulate variant: result is *added* to existing C (not overwritten).
//
// Needed in two situations:
//   1. GOTO PC loop (pc > 0): C tile already has partial sums from
//      earlier rank-KC updates within the same GOTO call.
//   2. SUMMA K-panel iteration: C_local is accumulated across multiple
//      SUMMA steps, each contributing a rank-KC_summa partial product.
//
// Implementation: accumulate MMA result into zero-initialised D, unpack
// into a stack temporary (BM×BN contiguous), then element-wise add to C.
// The extra add pass costs O(MR×NR) = O(1024) FLOPs vs. O(MR×NR×KC/BK)
// ≈ O(8192) FLOPs in the MMA loop, so overhead is ~12%.
static inline void amx_microkernel_accum(const IN_T* __restrict__ ap,
                                         const IN_T* __restrict__ bp,
                                         OUT_T* __restrict__ C_ptr,
                                         int kc_tiles, int ldc,
                                         int a_stride, int b_stride)
{
    // Load existing C into FragmentD
    // Since FragmentD doesn't have a direct "load from strided" we do it manually:
    // accumulate MMA result into zeroed D, then add to C.
    typename GEMM::FragmentD df{};
    for (int tk = 0; tk < kc_tiles; ++tk) {
        typename GEMM::FragmentA af;
        typename GEMM::FragmentB bf;
        af.template load(ap);
        bf.template load(bp);
        GEMM::MMA::mma(af, bf, df);
        ap += a_stride;
        bp += b_stride;
    }
    // Unpack into temp, add to C
    alignas(64) OUT_T tmp[BM * BN];
    df.template unload_unpack<MajorAxis::ROW>(tmp, BN);
    for (int i = 0; i < BM; ++i)
        for (int j = 0; j < BN; ++j)
            C_ptr[i * ldc + j] += tmp[i * BN + j];
}

// ═══════════════════════════════════════════════════════════════════════════
// §6  GOTO GEMM — 5-LOOP NEST
// ═══════════════════════════════════════════════════════════════════════════
//
// C[M×N] += A[M×K] · B[K×N],  all row-major, fp16 inputs, fp32 output.
//
// ┌───────────────────────────────────────────────────────────────────────┐
// │  DATA FLOW DIAGRAM                                                   │
// │                                                                      │
// │  Global A [M×K]          Global B [K×N]           Global C [M×N]    │
// │       │                       │                        ↑             │
// │  ┌────▼──────┐           ┌────▼──────┐                 │             │
// │  │ Pack Ã    │           │ Pack B̃    │  ← 1 per (JC,PC), shared    │
// │  │ MC×KC     │           │ KC×NC     │    across all IC threads     │
// │  │ per IC    │           │ in L3     │                              │
// │  │ thread    │           └─────┬─────┘                              │
// │  │ in L2     │                 │                                     │
// │  └─────┬─────┘                 │                                     │
// │        │         ┌─────────────▼─────────────┐                      │
// │   Ã[ir_tile]     │  B̃[jr_tile]               │                      │
// │   MR×KC          │  KC×NR                    │                      │
// │        └────┬────┘                           │                      │
// │             ▼                                │                      │
// │     ┌───────────────┐                        │                      │
// │     │ AMX Microkernel│  32×32 tile of C      │                      │
// │     │ rank-KC update │  accumulated in regs  │                      │
// │     └───────┬───────┘                        │                      │
// │             ▼                                │                      │
// │       C[ic+ir, jc+jr] ──────────────────────-┘                      │
// └───────────────────────────────────────────────────────────────────────┘
//
// Threading strategy:
//   • B̃ packing: all threads collaborate (parallel over tiles).
//   • IC loop: statically partitioned across threads.
//     Each thread packs its own Ã (private → NUMA-local via first-touch).
//   • JR and IR loops: sequential within each thread.
//     JR streams through B̃ micro-panels (L3); IR streams Ã micro-panels (L2).
//   • K dimension (PC loop) is NOT parallelised: all threads share the
//     same rank-KC update and write to disjoint C blocks → no reduction.
//
struct GotoParams {
    int MC, KC, NC;
};

static void goto_gemm(const IN_T* __restrict__ A,
                      const IN_T* __restrict__ B,
                      OUT_T* __restrict__ C,
                      int M, int N, int K,
                      int nthreads,
                      const GotoParams& p,
                      bool accum = false)  // if true, C += A*B; else C = A*B
{
    const int MC = p.MC, KC = p.KC, NC = p.NC;

    // ── Workspace allocation ────────────────────────────────────────────
    // B̃ (shared): one KC×NC panel, read by all IC threads from L3.
    //   Allocated once, reused across all (JC, PC) iterations.
    // Ã (private): one MC×KC block per thread, lives in each thread's L2.
    //   Pool of nthreads buffers; each thread indexes Ap_pool[tid * Ap_sz].
    const size_t Bp_sz = (size_t)KC * NC;           // max fp16 elements
    const size_t Ap_sz = (size_t)MC * KC;           // max per thread
    IN_T* Bp      = typed_alloc<IN_T>(Bp_sz);             // shared B̃
    IN_T* Ap_pool = typed_alloc<IN_T>(Ap_sz * nthreads);  // per-thread Ã

    // Zero C if not accumulating
    if (!accum) memset(C, 0, (size_t)M * N * sizeof(OUT_T));

    // ════════════════════════════════════════════════════════════════════
    // LOOP 5 (JC): partition columns of B and C into panels of width NC.
    //
    // Purpose: limit the packed B̃ panel to fit in L3.
    // KC × NC × 2B bytes → must be ≤ L3 capacity.
    // For N ≤ NC this loop runs once (common case for medium matrices).
    // ════════════════════════════════════════════════════════════════════
    for (int jc = 0; jc < N; jc += NC) {
        const int nc = std::min(NC, N - jc);

        // ════════════════════════════════════════════════════════════════
        // LOOP 4 (PC): rank-KC update.
        //
        // Partition K into panels of depth KC.  At each step:
        //   C[:, jc:jc+nc] += A[:, pc:pc+kc] · B[pc:pc+kc, jc:jc+nc]
        //
        // First step (pc=0, not accumulating) writes fresh to C.
        // Subsequent steps add to existing C values (rank-KC partial sums).
        // The K dimension is NEVER parallelised in standard GOTO: doing so
        // would require a reduction across threads on the same C elements.
        // ════════════════════════════════════════════════════════════════
        for (int pc = 0; pc < K; pc += KC) {
            const int kc = std::min(KC, K - pc);
            const int kc_tiles = kc / BK;

            // ── Pack B̃[pc:pc+kc, jc:jc+nc] ─────────────────────────────
            // All IC threads collaborate on packing B̃.  The packed panel
            // is KC×NC in VNNI format, shared by all IC threads via L3.
            //
            // NUMA note: first-touch places B̃ pages on the node of whichever
            // thread writes them.  On a single socket this is fine (shared L3).
            // On multi-socket, cross-NUMA B̃ reads are the dominant bottleneck
            // (see OpenBLAS issue #611).  The NUMA-SUMMA wrapper solves this
            // by giving each NUMA domain its own B̃ copy — analogous to BLIS's
            // JC-parallel strategy that replicates B̃ per socket.
            #pragma omp parallel num_threads(nthreads) proc_bind(close)
            {
                int tid = omp_get_thread_num();
                int nth = omp_get_num_threads();
                pack_B_panel_par(B, Bp, pc, jc, kc, nc, N, tid, nth);
            }

            const int tn_ = nc / BN;   // tiles in N for this panel
            const int tk_ = kc / BK;

            // ════════════════════════════════════════════════════════════
            // LOOP 3 (IC): partition rows of A into blocks of height MC.
            // This is the parallelised loop (BLIS "IC-parallel" strategy).
            //
            // Each thread:
            //   1. Packs Ã[ic:ic+mc, pc:pc+kc] into its private buffer.
            //      Ã fits in L2 (MC×KC×2B ≈ 512 KB for defaults).
            //      Because each thread writes its own Ã, first-touch places
            //      these pages on the thread's NUMA node → Ã is NUMA-local.
            //   2. Runs the JR × IR microkernel loops over its IC block.
            //      All threads read from the *same* shared B̃ in L3, but
            //      write to disjoint rows of C → no synchronisation needed.
            //
            // GOTO insight for NUMA:
            //   Ã is always NUMA-local (private per thread).
            //   B̃ is shared (L3 resident, same socket is fine).
            //   C writes are partitioned by rows → no false sharing.
            // ════════════════════════════════════════════════════════════
            #pragma omp parallel num_threads(nthreads) proc_bind(close)
            {
                GEMM::init();
                int tid = omp_get_thread_num();
                int nth = omp_get_num_threads();

                // per-thread Ã buffer
                IN_T* Ap = Ap_pool + (size_t)tid * Ap_sz;

                // static partition of IC range
                int ic_total = (M + MC - 1) / MC;
                int ic_base  = ic_total / nth, ic_extra = ic_total % nth;
                int ic_start = tid * ic_base + std::min(tid, ic_extra);
                int ic_count = ic_base + (tid < ic_extra ? 1 : 0);

                for (int ic_idx = ic_start; ic_idx < ic_start + ic_count; ++ic_idx) {
                    int ic = ic_idx * MC;
                    int mc = std::min(MC, M - ic);
                    int tm_ = mc / BM;

                    // Pack Ã[ic:ic+mc, pc:pc+kc]
                    pack_A_block(A, Ap, ic, pc, mc, kc, K);

                    // ════════════════════════════════════════════════════
                    // LOOP 2 (JR): iterate over NR-wide micro-panels of B̃.
                    //
                    // Each step processes one BN=32-wide column strip of B̃.
                    // While JR iterates, the packed Ã block stays in L2
                    // (being reused across all NC/NR column strips of B̃).
                    // The current B̃ micro-panel (KC × NR = 256×32 = 16 KB)
                    // streams from L3 → L1.
                    //
                    // This is where Ã reuse happens: one MC×KC pack of Ã
                    // is amortised over NC/NR = e.g. 4096/32 = 128 JR steps.
                    // ════════════════════════════════════════════════════
                    for (int jr = 0; jr < tn_; ++jr) {
                        // B̃ micro-panel pointer: tile index jr, all K-tiles.
                        // Layout is [tn][tk][B_TSZ], so jr advances by tk_*B_TSZ.
                        const IN_T* bp_base = Bp + (size_t)jr * tk_ * B_TSZ;

                        // ════════════════════════════════════════════════
                        // LOOP 1 (IR): iterate over MR-tall micro-panels of Ã.
                        //
                        // Each step processes one BM=32-tall row strip of Ã.
                        // While IR iterates, the current B̃ micro-panel
                        // (KC × NR) stays in L1, being reused across all
                        // MC/MR row strips.
                        //
                        // The Ã micro-panel (MR × KC = 32×256 = 16 KB)
                        // streams from L2.
                        //
                        // At the innermost level, the AMX microkernel:
                        //   - Loads one 32×32 A tile and one 32×32 B tile
                        //     per K-step (BK=32).
                        //   - Accumulates KC/BK = 8 rank-32 updates in the
                        //     8 AMX tile registers.
                        //   - Stores the 32×32 fp32 result to C.
                        // ════════════════════════════════════════════════
                        for (int ir = 0; ir < tm_; ++ir) {
                            // Ã micro-panel pointer: tile row ir, all K-tiles.
                            const IN_T* ap_base = Ap + (size_t)ir * tk_ * A_TSZ;

                            OUT_T* C_ptr = C + (ic + ir * BM) * N + (jc + jr * BN);

                            // Dispatch: overwrite vs. accumulate.
                            //
                            // On the first PC step (pc==0) with accum=false,
                            // this is the first time we touch this C tile →
                            // we can directly store the MMA result (overwrite).
                            //
                            // On subsequent PC steps (pc > 0) or when called
                            // with accum=true (e.g. from SUMMA K-panel
                            // iteration), the C tile already has partial sums
                            // from earlier rank-KC updates → we must *add*
                            // the new MMA result to the existing C values.
                            if (accum || pc > 0) {
                                amx_microkernel_accum(ap_base, bp_base, C_ptr,
                                                      kc_tiles, N, A_TSZ, B_TSZ);
                            } else {
                                amx_microkernel(ap_base, bp_base, C_ptr,
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
// §7  NUMA DOMAIN CONFIGURATION
// ═══════════════════════════════════════════════════════════════════════════
//
// We model the system as a PX × PY grid of "NUMA domains" (logical units).
// Each domain has a set of pinned CPU cores and per-domain matrix buffers.
//
// ┌──────────────────────────────────────────────────────────────────────┐
// │  Mapping to hardware NUMA nodes:                                    │
// │                                                                     │
// │    domain d  →  NUMA node (d % num_nodes)                          │
// │                                                                     │
// │  Examples:                                                          │
// │    • 2-socket (2 HW nodes), 2×2 grid:                              │
// │        domain 0,2 → node 0;  domain 1,3 → node 1                  │
// │    • SNC4 mode (4 HW nodes), 2×2 grid:                             │
// │        domain 0 → node 0, domain 1 → node 1, ...  (1:1 mapping)   │
// │    • 8-CCD AMD EPYC (8 HW nodes), 4×2 grid:                       │
// │        round-robin across all 8 CCDs                               │
// │                                                                     │
// │  C is partitioned in 2D: domain (px,py) owns                       │
// │    C[px*M/PX : (px+1)*M/PX,  py*N/PY : (py+1)*N/PY]             │
// │  This is the same decomposition used by SUMMA and Cannon.          │
// └──────────────────────────────────────────────────────────────────────┘
//

struct NumaDomain {
    int node_id;                // Linux NUMA node
    std::vector<int> cpus;      // CPUs pinned to this domain
    int cores;                  // number of cores for local GOTO

    // per-domain buffers (hugepage-backed, first-touch on this domain)
    IN_T*  A_local;    size_t A_bytes;
    IN_T*  B_local;    size_t B_bytes;
    OUT_T* C_local;    size_t C_bytes;

    // local matrix dimensions
    int M_local, N_local;
};

struct NumaGrid {
    int PX, PY;                          // domain grid
    int total_threads;
    std::vector<NumaDomain> domains;     // PX*PY domains

    void init(int px, int py, int total_cores) {
        PX = px; PY = py;
        total_threads = total_cores;
        int nd = px * py;
        domains.resize(nd);

        int num_nodes = get_num_numa_nodes();
        int cores_per_domain = total_cores / nd;

        // map domains to NUMA nodes round-robin
        for (int d = 0; d < nd; ++d) {
            int node = d % num_nodes;
            domains[d].node_id = node;
            domains[d].cores   = cores_per_domain;

            // assign CPUs
            auto all_cpus = cpus_on_node(node);
            // take a slice of CPUs for this domain
            int offset = (d / num_nodes) * cores_per_domain;
            for (int i = 0; i < cores_per_domain && (offset + i) < (int)all_cpus.size(); ++i) {
                domains[d].cpus.push_back(all_cpus[offset + i]);
            }
            // fallback: if not enough CPUs, wrap around
            while ((int)domains[d].cpus.size() < cores_per_domain && !all_cpus.empty()) {
                domains[d].cpus.push_back(all_cpus[domains[d].cpus.size() % all_cpus.size()]);
            }
        }
    }

    void alloc_buffers(int M, int N, int K) {
        int M_per = M / PX, N_per = N / PY;
        for (int px = 0; px < PX; ++px) {
            for (int py = 0; py < PY; ++py) {
                int d = px * PY + py;
                auto& dom = domains[d];
                dom.M_local = M_per;
                dom.N_local = N_per;
                dom.A_bytes = (size_t)M_per * K * sizeof(IN_T);
                dom.B_bytes = (size_t)K * N_per * sizeof(IN_T);
                dom.C_bytes = (size_t)M_per * N_per * sizeof(OUT_T);
                dom.A_local = static_cast<IN_T*>(numa_alloc(dom.A_bytes));
                dom.B_local = static_cast<IN_T*>(numa_alloc(dom.B_bytes));
                dom.C_local = static_cast<OUT_T*>(numa_alloc(dom.C_bytes));
            }
        }
    }

    void free_buffers() {
        for (auto& dom : domains) {
            numa_free(dom.A_local, dom.A_bytes);
            numa_free(dom.B_local, dom.B_bytes);
            numa_free(dom.C_local, dom.C_bytes);
            dom.A_local = nullptr;
            dom.B_local = nullptr;
            dom.C_local = nullptr;
        }
    }
};

// ═══════════════════════════════════════════════════════════════════════════
// §8  SCATTER / GATHER between global and per-domain buffers
// ═══════════════════════════════════════════════════════════════════════════
//
// Move data between the global row-major matrices and per-domain buffers.
//
// NUMA locality rule: the *destination* buffer is already allocated on the
// target NUMA node (via hugepage mmap + first-touch by a pinned thread).
// The memcpy writes to local DRAM, even if it reads from remote DRAM.
// This remote-read cost is paid once at distribution time; all subsequent
// accesses during the local GOTO GEMM are NUMA-local.
//
// For the SUMMA variant, these copies happen at each K-panel step.
// The cost is bounded by the panel size (M_loc × KC_summa + KC_summa × N_loc),
// which is small relative to the O(M_loc × N_loc × KC_summa) compute.

// Copy rows [row0..row0+nrows) of global A[M×K] into domain A_local[nrows×K]
static void scatter_A_rows(const IN_T* A, IN_T* dst,
                           int row0, int nrows, int K) {
    memcpy(dst, A + (size_t)row0 * K, (size_t)nrows * K * sizeof(IN_T));
}

// Copy columns [col0..col0+ncols) of global B[K×N] into domain B_local[K×ncols]
static void scatter_B_cols(const IN_T* B, IN_T* dst,
                           int col0, int ncols, int K, int N) {
    for (int k = 0; k < K; ++k) {
        memcpy(dst + (size_t)k * ncols,
               B + (size_t)k * N + col0,
               ncols * sizeof(IN_T));
    }
}

// Gather domain C_local[nrows×ncols] back into global C[M×N]
static void gather_C_block(OUT_T* C, const OUT_T* src,
                           int row0, int col0, int nrows, int ncols, int N) {
    for (int i = 0; i < nrows; ++i) {
        memcpy(C + (size_t)(row0 + i) * N + col0,
               src + (size_t)i * ncols,
               ncols * sizeof(OUT_T));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §9  NUMA-STATIC: pre-distribute, local GOTO, gather
// ═══════════════════════════════════════════════════════════════════════════
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  ANALOGY TO BLIS's JC-parallel strategy                                │
// │                                                                        │
// │  In BLIS, JC-parallel assigns each socket a different NC-wide column   │
// │  panel of B and the corresponding columns of C.  Each socket packs    │
// │  its own B̃ → NUMA-local.  IC-parallel within the socket handles A.   │
// │                                                                        │
// │  Our NUMA-Static extends this to a 2D decomposition:                   │
// │    • PX partitions along M (rows of A and C)                           │
// │    • PY partitions along N (columns of B and C)                        │
// │    • Domain (px,py) owns C[px_rows, py_cols] and stores locally:      │
// │        A_local = A[px_rows, :]     (full K columns — row strip)       │
// │        B_local = B[:, py_cols]     (full K rows — column strip)        │
// │    • Each domain runs GOTO independently → zero cross-NUMA traffic.   │
// │                                                                        │
// │  Memory per domain: M_loc × K + K × N_loc + M_loc × N_loc            │
// │  For M=N=4096, K=4096, 2×2 grid, fp16:                               │
// │    A_local = 2048 × 4096 × 2B = 16 MB                                │
// │    B_local = 4096 × 2048 × 2B = 16 MB                                │
// │    C_local = 2048 × 2048 × 4B = 16 MB                                │
// │  Total: 48 MB per domain.  Fits in DRAM, pages are NUMA-local.        │
// │                                                                        │
// │  The SUMMA variant (§10) reduces A/B memory to KC_summa depth.        │
// └─────────────────────────────────────────────────────────────────────────┘
//
static void numa_static_gemm(const IN_T* A, const IN_T* B, OUT_T* C,
                             int M, int N, int K,
                             NumaGrid& grid, const GotoParams& gp)
{
    int PX = grid.PX, PY = grid.PY;
    int M_loc = M / PX, N_loc = N / PY;
    int nd = PX * PY;

    // Scatter A rows and B cols to per-domain buffers
    // Use domain threads for first-touch locality
    #pragma omp parallel num_threads(nd) proc_bind(spread)
    {
        int d = omp_get_thread_num();
        int px = d / PY, py = d % PY;
        auto& dom = grid.domains[d];

        // pin to domain CPUs
        if (!dom.cpus.empty()) pin_thread_to_cpus(dom.cpus);

        scatter_A_rows(A, dom.A_local, px * M_loc, M_loc, K);
        scatter_B_cols(B, dom.B_local, py * N_loc, N_loc, K, N);
        memset(dom.C_local, 0, dom.C_bytes);
    }

    // Run GOTO GEMM on each domain in parallel
    // Each domain uses its own thread team
    int cores_per = grid.domains[0].cores;

    #pragma omp parallel num_threads(nd * cores_per) proc_bind(close)
    {
        int gtid = omp_get_thread_num();
        int d = gtid / cores_per;
        auto& dom = grid.domains[d];

        // Only first thread per domain runs the GOTO
        // The GOTO itself will spawn nested parallelism
        if (gtid % cores_per == 0) {
            // Adjust OMP nesting for local GOTO
            goto_gemm(dom.A_local, dom.B_local, dom.C_local,
                      dom.M_local, dom.N_local, K,
                      dom.cores, gp);
        }
    }

    // Gather C blocks back
    for (int px = 0; px < PX; ++px) {
        for (int py = 0; py < PY; ++py) {
            int d = px * PY + py;
            gather_C_block(C, grid.domains[d].C_local,
                           px * M_loc, py * N_loc, M_loc, N_loc, N);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §10  NUMA-SUMMA: K-panel iteration with local GOTO
// ═══════════════════════════════════════════════════════════════════════════
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  THE SUMMA ALGORITHM (van de Geijn & Watts, 1995)                      │
// │                                                                        │
// │  SUMMA distributes C = A·B across a PX×PY processor grid.             │
// │  Each processor (px,py) owns a block of C, and the K dimension is     │
// │  processed in panels of width KC_summa:                                │
// │                                                                        │
// │    for kp = 0 to K step KC_summa:                                      │
// │      1. BROADCAST A_panel[px_rows, kp:kp+kc] along processor row px   │
// │      2. BROADCAST B_panel[kp:kp+kc, py_cols] along processor col py   │
// │      3. C_local += A_panel · B_panel  (local GEMM)                    │
// │                                                                        │
// │  In distributed MPI, "broadcast" = MPI_Bcast along row/column comm.   │
// │  In our shared-memory NUMA context, "broadcast" = memcpy from the     │
// │  global matrix to a NUMA-local buffer.  The destination is on the     │
// │  receiving domain's NUMA node; the source may be remote but that cost │
// │  is paid once per SUMMA step, not in the hot inner loops.             │
// │                                                                        │
// │  CORRESPONDENCE: SUMMA ↔ GOTO                                         │
// │                                                                        │
// │  SUMMA's outer K-panel loop is structurally identical to GOTO's PC    │
// │  loop.  When KC_summa == KC, the two loops collapse: SUMMA step =     │
// │  GOTO PC step, and the local GOTO's PC loop runs exactly once.        │
// │  When KC_summa > KC, the local GOTO's PC loop does KC_summa/KC steps, │
// │  all on NUMA-local data.                                               │
// │                                                                        │
// │  ADVANTAGE OVER STATIC: reduced per-domain memory.                    │
// │    Static:  M_loc × K   of A per domain (full K depth)               │
// │    SUMMA:   M_loc × KC_summa of A per domain (one panel)             │
// │                                                                        │
// │  COMMUNICATION / COMPUTE RATIO:                                        │
// │    Comm per step: (M_loc × kc + kc × N_loc) × 2B  bytes             │
// │    Comp per step: 2 × M_loc × N_loc × kc  FLOPs                     │
// │    Ratio ≈ O(1 / min(M_loc, N_loc))  → favourable for large M, N.   │
// └─────────────────────────────────────────────────────────────────────────┘
//
static void numa_summa_gemm(const IN_T* A, const IN_T* B, OUT_T* C,
                            int M, int N, int K,
                            NumaGrid& grid, const GotoParams& gp,
                            int KC_summa = 0)  // 0 = use full K (static equiv.)
{
    int PX = grid.PX, PY = grid.PY;
    int M_loc = M / PX, N_loc = N / PY;
    int nd = PX * PY;

    if (KC_summa <= 0) KC_summa = K;  // full-K = one SUMMA step = static
    // round to multiple of BK
    KC_summa = ((KC_summa + BK - 1) / BK) * BK;
    if (KC_summa > K) KC_summa = K;

    // Allocate SUMMA panel buffers per domain
    size_t A_panel_sz = (size_t)M_loc * KC_summa;
    size_t B_panel_sz = (size_t)KC_summa * N_loc;
    std::vector<IN_T*> A_panels(nd), B_panels(nd);
    for (int d = 0; d < nd; ++d) {
        A_panels[d] = typed_alloc<IN_T>(A_panel_sz);
        B_panels[d] = typed_alloc<IN_T>(B_panel_sz);
    }

    // Zero C
    memset(C, 0, (size_t)M * N * sizeof(OUT_T));

    // Zero local C
    for (int d = 0; d < nd; ++d) {
        memset(grid.domains[d].C_local, 0, grid.domains[d].C_bytes);
    }

    // SUMMA outer loop over K panels
    for (int kp = 0; kp < K; kp += KC_summa) {
        int kc = std::min(KC_summa, K - kp);

        // Phase 1: Copy panels to NUMA-local buffers (parallel per domain)
        #pragma omp parallel num_threads(nd) proc_bind(spread)
        {
            int d = omp_get_thread_num();
            int px = d / PY, py = d % PY;
            auto& dom = grid.domains[d];
            if (!dom.cpus.empty()) pin_thread_to_cpus(dom.cpus);

            // Copy A panel: rows [px*M_loc .. (px+1)*M_loc), cols [kp..kp+kc)
            for (int i = 0; i < M_loc; ++i) {
                memcpy(A_panels[d] + (size_t)i * kc,
                       A + (size_t)(px * M_loc + i) * K + kp,
                       kc * sizeof(IN_T));
            }
            // Copy B panel: rows [kp..kp+kc), cols [py*N_loc .. (py+1)*N_loc)
            for (int k = 0; k < kc; ++k) {
                memcpy(B_panels[d] + (size_t)k * N_loc,
                       B + (size_t)(kp + k) * N + py * N_loc,
                       N_loc * sizeof(IN_T));
            }
        }

        // Phase 2: Local GOTO GEMM on each domain (C_local += A_panel · B_panel)
        // Run sequentially across domains, each using its own thread team
        for (int d = 0; d < nd; ++d) {
            auto& dom = grid.domains[d];
            GotoParams local_gp = gp;
            // Adjust KC if kc < gp.KC
            if (kc < local_gp.KC) local_gp.KC = ((kc + BK - 1) / BK) * BK;

            goto_gemm(A_panels[d], B_panels[d], dom.C_local,
                      dom.M_local, dom.N_local, kc,
                      dom.cores, local_gp,
                      /*accum=*/ (kp > 0));
        }
    }

    // Gather C
    for (int px = 0; px < PX; ++px) {
        for (int py = 0; py < PY; ++py) {
            int d = px * PY + py;
            gather_C_block(C, grid.domains[d].C_local,
                           px * M_loc, py * N_loc, M_loc, N_loc, N);
        }
    }

    for (int d = 0; d < nd; ++d) {
        typed_free(A_panels[d], A_panel_sz);
        typed_free(B_panels[d], B_panel_sz);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// §11  NUMA-CANNON: shift-rotate on square grid
// ═══════════════════════════════════════════════════════════════════════════
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │  CANNON'S ALGORITHM (1969)                                             │
// │                                                                        │
// │  Requires a square P×P processor grid and K divisible by P.           │
// │  A, B, C are partitioned into P×P blocks of size                      │
// │  (M/P)×(K/P), (K/P)×(N/P), (M/P)×(N/P) respectively.               │
// │                                                                        │
// │  Algorithm:                                                            │
// │    1. Initial skew (pre-rotation):                                     │
// │         Domain (i,j) receives A_block(i, (j+i) mod P)                │
// │                           and B_block((i+j) mod P, j)                │
// │    2. For step = 0 to P-1:                                            │
// │         a. C_local += A_block · B_block  (local GOTO GEMM)           │
// │         b. Shift A blocks LEFT by 1 within each row of the grid      │
// │         c. Shift B blocks UP by 1 within each column of the grid     │
// │                                                                        │
// │  After P steps, each domain has accumulated all P partial products,   │
// │  i.e. C(i,j) = Σ_{k=0}^{P-1} A(i,k) × B(k,j).                     │
// │                                                                        │
// │  WHY THE INITIAL SKEW?                                                 │
// │  Without it, domain (i,j) starts with A(i,j) and B(i,j).  After P    │
// │  left-shifts of A and P up-shifts of B, domain (i,j) sees blocks:    │
// │    A(i,j), A(i,j-1), A(i,j-2), ...                                   │
// │    B(i,j), B(i-1,j), B(i-2,j), ...                                   │
// │  which gives Σ_s A(i, j-s) × B(i-s, j) — only correct when the      │
// │  skew aligns the initial blocks so the k-indices match at each step.  │
// │                                                                        │
// │  vs. SUMMA:                                                            │
// │    • Cannon: nearest-neighbour shifts (lower latency per step)        │
// │    • SUMMA:  broadcasts (higher bandwidth utilisation)                │
// │    • Cannon: constant memory (one A + one B block per domain)         │
// │    • SUMMA:  panel can be wider (KC_summa) for better compute ratio   │
// │    • Cannon: requires square grid; SUMMA works on any PX×PY           │
// │                                                                        │
// │  In our shared-memory NUMA implementation, "shift" = memcpy between   │
// │  per-domain buffers.  The memcpy destination is always NUMA-local     │
// │  (first-touched by the receiving domain's thread).                    │
// │  We double-buffer (A_blk + A_tmp) so shifts don't overwrite in-use   │
// │  data, then swap the pointer vectors after each shift.                │
// └─────────────────────────────────────────────────────────────────────────┘
//
static void numa_cannon_gemm(const IN_T* A, const IN_T* B, OUT_T* C,
                             int M, int N, int K,
                             NumaGrid& grid, const GotoParams& gp)
{
    int P = grid.PX;
    if (grid.PX != grid.PY) {
        fprintf(stderr, "Cannon requires square grid (PX==PY), got %dx%d\n",
                grid.PX, grid.PY);
        return;
    }

    int M_loc = M / P, N_loc = N / P, K_loc = K / P;
    int nd = P * P;

    // Allocate per-domain A/B blocks (M_loc × K_loc each)
    size_t A_blk_sz = (size_t)M_loc * K_loc;
    size_t B_blk_sz = (size_t)K_loc * N_loc;
    std::vector<IN_T*> A_blk(nd), B_blk(nd), A_tmp(nd), B_tmp(nd);
    for (int d = 0; d < nd; ++d) {
        A_blk[d] = typed_alloc<IN_T>(A_blk_sz);
        B_blk[d] = typed_alloc<IN_T>(B_blk_sz);
        A_tmp[d] = typed_alloc<IN_T>(A_blk_sz);
        B_tmp[d] = typed_alloc<IN_T>(B_blk_sz);
    }

    // Zero local C
    for (int d = 0; d < nd; ++d)
        memset(grid.domains[d].C_local, 0, grid.domains[d].C_bytes);

    // ── Step 1: Initial skew (pre-rotation) ────────────────────────────
    //
    // Domain (i,j) receives:
    //   A block at column index (j+i) mod P  →  A[i, (j+i)%P]
    //   B block at row    index (i+j) mod P  →  B[(i+j)%P, j]
    //
    // After the skew, domain (i,j) holds blocks whose K-indices differ
    // by exactly (i+j) from the "natural" position.  The P subsequent
    // shift-and-multiply steps then sweep each domain through all P
    // K-blocks in the correct order.
    //
    // Example for P=2:
    //   Dom(0,0): A(0,0), B(0,0)  →  step 0: C(0,0) += A(0,0)·B(0,0)
    //   after shift: A(0,1), B(1,0) →  step 1: C(0,0) += A(0,1)·B(1,0)  ✓
    //   Dom(0,1): A(0,1), B(1,1)  →  step 0: C(0,1) += A(0,1)·B(1,1)
    //   after shift: A(0,0), B(0,1) →  step 1: C(0,1) += A(0,0)·B(0,1)  ✓
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < P; ++j) {
            int d = i * P + j;
            int a_col_blk = (j + i) % P;
            int b_row_blk = (i + j) % P;

            // Copy A block [i*M_loc .. (i+1)*M_loc, a_col_blk*K_loc .. ]
            for (int r = 0; r < M_loc; ++r) {
                memcpy(A_blk[d] + (size_t)r * K_loc,
                       A + (size_t)(i * M_loc + r) * K + a_col_blk * K_loc,
                       K_loc * sizeof(IN_T));
            }
            // Copy B block [b_row_blk*K_loc .. , j*N_loc .. ]
            for (int r = 0; r < K_loc; ++r) {
                memcpy(B_blk[d] + (size_t)r * N_loc,
                       B + (size_t)(b_row_blk * K_loc + r) * N + j * N_loc,
                       N_loc * sizeof(IN_T));
            }
        }
    }

    // ── Step 2: P iterations of local-GEMM + shift ────────────────────
    //
    // At each step:
    //   2a. Each domain computes C_local += A_blk × B_blk via local GOTO.
    //       All data is NUMA-local → full local bandwidth.
    //   2b. Circular left-shift of A within each row: (i,j) → (i, (j-1)%P).
    //   2c. Circular up-shift of B within each column: (i,j) → ((i-1)%P, j).
    //
    // Shifts use double-buffering: write into A_tmp/B_tmp, then swap vectors.
    // The memcpy destination is on the *receiving* domain's NUMA node,
    // so after the shift the data is again NUMA-local for the next GEMM.
    //
    // After P steps, domain (i,j) has seen all P A-column-blocks and
    // B-row-blocks → C(i,j) = Σ_{k} A(i,k) × B(k,j).  Complete.
    GotoParams local_gp = gp;
    if (K_loc < local_gp.KC) local_gp.KC = ((K_loc + BK - 1) / BK) * BK;

    for (int step = 0; step < P; ++step) {
        // Local GOTO GEMM on each domain: C_local += A_blk · B_blk
        for (int d = 0; d < nd; ++d) {
            auto& dom = grid.domains[d];
            goto_gemm(A_blk[d], B_blk[d], dom.C_local,
                      M_loc, N_loc, K_loc,
                      dom.cores, local_gp,
                      /*accum=*/ true);
        }

        if (step < P - 1) {
            // Shift A blocks left within each row: domain (i,j) sends to (i, (j-1+P)%P)
            for (int i = 0; i < P; ++i) {
                for (int j = 0; j < P; ++j) {
                    int src = i * P + j;
                    int dst_j = (j - 1 + P) % P;
                    int dst_d = i * P + dst_j;
                    memcpy(A_tmp[dst_d], A_blk[src], A_blk_sz * sizeof(IN_T));
                }
            }
            std::swap(A_blk, A_tmp);

            // Shift B blocks up within each column: domain (i,j) sends to ((i-1+P)%P, j)
            for (int i = 0; i < P; ++i) {
                for (int j = 0; j < P; ++j) {
                    int src = i * P + j;
                    int dst_i = (i - 1 + P) % P;
                    int dst_d = dst_i * P + j;
                    memcpy(B_tmp[dst_d], B_blk[src], B_blk_sz * sizeof(IN_T));
                }
            }
            std::swap(B_blk, B_tmp);
        }
    }

    // Gather C
    memset(C, 0, (size_t)M * N * sizeof(OUT_T));
    for (int i = 0; i < P; ++i)
        for (int j = 0; j < P; ++j)
            gather_C_block(C, grid.domains[i * P + j].C_local,
                           i * M_loc, j * N_loc, M_loc, N_loc, N);

    for (int d = 0; d < nd; ++d) {
        typed_free(A_blk[d], A_blk_sz); typed_free(B_blk[d], B_blk_sz);
        typed_free(A_tmp[d], A_blk_sz); typed_free(B_tmp[d], B_blk_sz);
    }
}

// ═══════════════════════════════════════════════════════════════════════
// §12  MKL BASELINES
// ═══════════════════════════════════════════════════════════════════════

static float* g_A32_buf = nullptr;
static float* g_B32_buf = nullptr;

static void mkl_sgemm_f32(const IN_T* A, const IN_T* B, OUT_T* C,
                           int M, int N, int K) {
    mkl_set_dynamic(0);
    mkl_set_num_threads(g_threads);
    #pragma omp parallel for num_threads(g_threads) schedule(static)
    for (size_t i = 0; i < (size_t)M * K; ++i) g_A32_buf[i] = (float)A[i];
    #pragma omp parallel for num_threads(g_threads) schedule(static)
    for (size_t i = 0; i < (size_t)K * N; ++i) g_B32_buf[i] = (float)B[i];
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0f, g_A32_buf, K, g_B32_buf, N, 0.0f, C, N);
}

static void mkl_f16_gemm(const IN_T* A, const IN_T* B, OUT_T* C,
                          int M, int N, int K) {
    mkl_set_dynamic(0);
    mkl_set_num_threads(g_threads);
    cblas_gemm_f16f16f32(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                         M, N, K, 1.0f,
                         (const MKL_F16*)A, K,
                         (const MKL_F16*)B, N,
                         0.0f, C, N);
}

// ═══════════════════════════════════════════════════════════════════════
// §13  VERIFICATION
// ═══════════════════════════════════════════════════════════════════════
//
// Compare kernel output C against the OpenMP scalar reference.
// Reports: max absolute error, max relative error, mean absolute error,
// number of mismatches, and the first 3 mismatched locations.
//
// Tolerance is generous (2% relative OR 0.05 absolute) because fp16
// inputs have only ~3 decimal digits of precision, and different
// accumulation orders (scalar ikj vs. packed tiled AMX) produce
// legitimately different rounding.
//
static bool verify(const OUT_T* C, const OUT_T* ref, int M, int N,
                   const char* label, float tol = 0.02f) {
    int bad = 0;
    double max_abs = 0, max_rel = 0, sum_abs = 0;
    int first_bad[3] = {-1, -1, -1};
    size_t total = (size_t)M * N;

    for (size_t i = 0; i < total; ++i) {
        double d = fabs((double)C[i] - (double)ref[i]);
        double den = fabs((double)ref[i]);
        double rel = (den > 1e-6) ? d / den : d;
        sum_abs += d;
        if (d > max_abs) max_abs = d;
        if (rel > max_rel) max_rel = rel;
        if (rel > tol && d > 0.05f) {
            if (bad < 3) first_bad[bad] = (int)i;
            ++bad;
        }
    }

    double mean_abs = sum_abs / total;

    if (bad == 0) {
        printf("    [PASS] %-24s max_abs=%.6f max_rel=%.4f%% mean_abs=%.6f\n",
               label, max_abs, max_rel * 100.0, mean_abs);
    } else {
        printf("    [FAIL] %-24s %d mismatches  max_abs=%.6f max_rel=%.4f%%\n",
               label, bad, max_abs, max_rel * 100.0);
        for (int k = 0; k < 3 && first_bad[k] >= 0; ++k) {
            int idx = first_bad[k];
            int row = idx / N, col = idx % N;
            printf("           idx=%d (%d,%d) got=%.6f ref=%.6f diff=%.6f\n",
                   idx, row, col, C[idx], ref[idx],
                   fabsf(C[idx] - ref[idx]));
        }
    }
    return bad == 0;
}

// ═══════════════════════════════════════════════════════════════════════
// §14  BENCHMARK HARNESS
// ═══════════════════════════════════════════════════════════════════════

using KernFn = std::function<void(const IN_T*, const IN_T*, OUT_T*, int, int, int)>;

static void bench_one(FILE* f, const char* name, KernFn fn,
                      const IN_T* A, const IN_T* B, OUT_T* C,
                      int M, int N, int K, int nthreads) {
    for (int w = 0; w < g_warmup; ++w) fn(A, B, C, M, N, K);

    double flops = 2.0 * (double)M * N * K;
    std::vector<double> times(g_iters);
    double total = 0;
    for (int i = 0; i < g_iters; ++i) {
        double t0 = now_ns();
        fn(A, B, C, M, N, K);
        double t1 = now_ns();
        times[i] = t1 - t0;
        total += times[i];
    }
    double avg = total / g_iters;
    std::sort(times.begin(), times.end());
    double med = (g_iters % 2) ? times[g_iters / 2]
                                : (times[g_iters / 2 - 1] + times[g_iters / 2]) / 2.0;

    for (int i = 0; i < g_iters; ++i)
        fprintf(f, "%s,%d,%d,%d,%d,%d,%.2f,%.6f\n",
                name, M, N, K, nthreads, i, times[i], flops / times[i]);

    printf("    %-28s t=%2d avg=%.0f ns  med=%.0f ns  %.4f GF/s (peak %.4f)\n",
           name, nthreads, avg, med, flops / avg, flops / times[0]);
}

// ═══════════════════════════════════════════════════════════════════════
// §15  ARG PARSING
// ═══════════════════════════════════════════════════════════════════════

struct MNK { int M, N, K; };

static void usage(const char* p) {
    printf("Usage: %s M N K [M2 N2 K2 ...] [options]\n", p);
    printf("  -w W         warmup iterations (default %d)\n", g_warmup);
    printf("  -i I         benchmark iterations (default %d)\n", g_iters);
    printf("  -t T         total threads (default %d)\n", g_threads);
    printf("  -p PXxPY     NUMA domain grid (default %dx%d)\n", g_PX, g_PY);
    printf("  -mc MC       GOTO MC blocksize (default %d)\n", g_MC);
    printf("  -kc KC       GOTO KC blocksize (default %d)\n", g_KC);
    printf("  -nc NC       GOTO NC blocksize (default %d)\n", g_NC);
    printf("  -o FILE      CSV output file (default %s)\n", g_outfile);
}

static std::vector<MNK> parse_args(int argc, char** argv) {
    std::vector<MNK> sizes;
    std::vector<int> nums;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-w") && i + 1 < argc) g_warmup = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-i") && i + 1 < argc) g_iters = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t") && i + 1 < argc) g_threads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-p") && i + 1 < argc) {
            ++i;
            if (sscanf(argv[i], "%dx%d", &g_PX, &g_PY) != 2) {
                fprintf(stderr, "Bad -p format, expected PXxPY (e.g. 2x2)\n");
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-mc") && i + 1 < argc) g_MC = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-kc") && i + 1 < argc) g_KC = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nc") && i + 1 < argc) g_NC = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o")  && i + 1 < argc) g_outfile = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(argv[0]); exit(0); }
        else nums.push_back(atoi(argv[i]));
    }
    if (nums.size() % 3 || nums.empty()) {
        fprintf(stderr, "Need M N K triples\n"); usage(argv[0]); exit(1);
    }
    for (size_t i = 0; i < nums.size(); i += 3)
        sizes.push_back({nums[i], nums[i + 1], nums[i + 2]});
    return sizes;
}

// ═══════════════════════════════════════════════════════════════════════
// §16  MAIN
// ═══════════════════════════════════════════════════════════════════════

int main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 1; }

    auto sizes = parse_args(argc, argv);
    GotoParams gp{g_MC, g_KC, g_NC};

    int nd = g_PX * g_PY;
    int num_nodes = get_num_numa_nodes();

    printf("════════════════════════════════════════════════════════\n");
    printf("  GOTO + NUMA GEMM Benchmark\n");
    printf("════════════════════════════════════════════════════════\n");
    printf("  Threads       : %d\n", g_threads);
    printf("  NUMA nodes    : %d (detected)\n", num_nodes);
    printf("  Domain grid   : %d x %d = %d domains\n", g_PX, g_PY, nd);
    printf("  Cores/domain  : %d\n", g_threads / nd);
    printf("  GOTO blocking : MC=%d KC=%d NC=%d\n", gp.MC, gp.KC, gp.NC);
    printf("  AMX tile      : %dx%dx%d\n", BM, BN, BK);
    printf("  Warmup/Iters  : %d / %d\n", g_warmup, g_iters);
    printf("  Output        : %s\n", g_outfile);
    printf("════════════════════════════════════════════════════════\n\n");

    // Print NUMA topology
    for (int n = 0; n < num_nodes; ++n) {
        auto cpus = cpus_on_node(n);
        printf("  Node %d: %zu CPUs [%d..%d]\n", n, cpus.size(),
               cpus.empty() ? -1 : cpus.front(),
               cpus.empty() ? -1 : cpus.back());
    }
    printf("\n");

    FILE* fout = fopen(g_outfile, "w");
    if (!fout) { perror("fopen"); return 1; }
    fprintf(fout, "kernel,M,N,K,threads,run,time_ns,gflops\n");

    // Setup NUMA grid
    NumaGrid grid;
    grid.init(g_PX, g_PY, g_threads);

    for (auto& dom : grid.domains)
        printf("  Domain %d: node=%d cores=%d cpus=[%d..%d]\n",
               (int)(&dom - &grid.domains[0]),
               dom.node_id, dom.cores,
               dom.cpus.empty() ? -1 : dom.cpus.front(),
               dom.cpus.empty() ? -1 : dom.cpus.back());
    printf("\n");

    // Enable nested parallelism for NUMA-domain-level GOTO
    omp_set_nested(1);
    omp_set_max_active_levels(3);

    for (auto& sz : sizes) {
        int M = sz.M, N = sz.N, K = sz.K;
        printf("═══ M=%d  N=%d  K=%d ═══════════════════════════════\n", M, N, K);

        // Check divisibility
        if (M % (g_PX * BM) || N % (g_PY * BN) || K % BK) {
            printf("  [SKIP] M must be div by PX*%d=%d, N by PY*%d=%d, K by %d\n",
                   BM, g_PX * BM, BN, g_PY * BN, BK);
            continue;
        }

        size_t szA = (size_t)M * K, szB = (size_t)K * N, szC = (size_t)M * N;

        // Allocate matrices (hugepage-backed)
        IN_T*  A    = typed_alloc<IN_T>(szA);
        IN_T*  B    = typed_alloc<IN_T>(szB);
        OUT_T* Cref = typed_alloc<OUT_T>(szC);
        OUT_T* Cchk = typed_alloc<OUT_T>(szC);

        g_A32_buf = typed_alloc<float>(szA);
        g_B32_buf = typed_alloc<float>(szB);

        // First-touch init with random data
        #pragma omp parallel num_threads(g_threads) proc_bind(close)
        {
            int t = omp_get_thread_num(), n_ = omp_get_num_threads();
            size_t lo, hi;
            unsigned seed = 42 + t;
            lo = t * szA / n_; hi = (t + 1) * szA / n_;
            for (size_t i = lo; i < hi; ++i)
                A[i] = (IN_T)((float)rand_r(&seed) / RAND_MAX - 0.5f);
            lo = t * szB / n_; hi = (t + 1) * szB / n_;
            for (size_t i = lo; i < hi; ++i)
                B[i] = (IN_T)((float)rand_r(&seed) / RAND_MAX - 0.5f);
        }

        // Allocate per-domain buffers
        grid.alloc_buffers(M, N, K);

        // ── Reference: OpenMP scalar GEMM (ground truth) ────────────
        //
        // We use the trivially-correct scalar GEMM as the primary reference.
        // This is independent of MKL/AMX/any optimised path.
        // We also compute MKL sgemm and cross-validate the two references
        // to catch any systematic issues.
        OUT_T* Cmkl = typed_alloc<OUT_T>(szC);  // second reference for cross-check

        printf("  Reference (OpenMP scalar GEMM)...\n");
        ref_gemm_omp(A, B, Cref, M, N, K);

        printf("  Cross-check (MKL sgemm fp32)...\n");
        mkl_sgemm_f32(A, B, Cmkl, M, N, K);

        // Cross-validate: OpenMP ref vs MKL.  If these disagree, something
        // is wrong with the test setup (bad data, overflow, etc.).
        {
            int cross_bad = 0;
            double max_diff = 0;
            for (size_t i = 0; i < szC; ++i) {
                double d = fabs((double)Cref[i] - (double)Cmkl[i]);
                double den = fabs((double)Cref[i]);
                double rel = (den > 1e-6) ? d / den : d;
                if (d > max_diff) max_diff = d;
                if (rel > 0.02 && d > 0.05) ++cross_bad;
            }
            if (cross_bad == 0)
                printf("    [PASS] ref_omp vs mkl_sgemm cross-check (max_diff=%.6f)\n", max_diff);
            else
                printf("    [WARN] ref_omp vs mkl_sgemm: %d mismatches (max_diff=%.6f)\n"
                       "           fp16 accumulation order may differ — using ref_omp as ground truth\n",
                       cross_bad, max_diff);
        }

        // ── Verify & Benchmark each kernel ──────────────────────────
        //
        // Each kernel is:
        //   1. Run once and verified against Cref (the OpenMP scalar result).
        //   2. Benchmarked for g_iters iterations with timing.
        //
        // Verification tolerance: 2% relative OR 0.05 absolute.
        // fp16 inputs have ~3 decimal digits of precision, so small
        // differences in accumulation order produce ~0.01–0.1 absolute error.

        printf("  Kernels:\n");

        // 0. OpenMP scalar baseline (benchmark only — it IS the reference)
        bench_one(fout, "ref_omp", ref_gemm_omp, A, B, Cchk, M, N, K, g_threads);
        printf("    [REF ] ref_omp (this is the ground truth)\n");

        // 1. MKL sgemm (fp32 upcast)
        bench_one(fout, "mkl_sgemm_f32", mkl_sgemm_f32, A, B, Cchk, M, N, K, g_threads);
        verify(Cchk, Cref, M, N, "mkl_sgemm_f32");

        // 2. MKL f16f16f32
        bench_one(fout, "mkl_f16f16f32", mkl_f16_gemm, A, B, Cchk, M, N, K, g_threads);
        verify(Cchk, Cref, M, N, "mkl_f16f16f32");

        // 3. GOTO GEMM (shared-memory, all threads)
        auto goto_fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
            goto_gemm(a, b, c, m, n, k, g_threads, gp);
        };
        bench_one(fout, "goto_amx", goto_fn, A, B, Cchk, M, N, K, g_threads);
        verify(Cchk, Cref, M, N, "goto_amx");

        // 4. NUMA-Static
        auto numa_static_fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
            numa_static_gemm(a, b, c, m, n, k, grid, gp);
        };
        {
            char label[64];
            snprintf(label, sizeof(label), "numa_static_%dx%d", g_PX, g_PY);
            bench_one(fout, label, numa_static_fn, A, B, Cchk, M, N, K, g_threads);
            verify(Cchk, Cref, M, N, label);
        }

        // 5. NUMA-SUMMA (with KC_summa = KC, i.e. fine-grained panels)
        auto numa_summa_fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
            numa_summa_gemm(a, b, c, m, n, k, grid, gp, gp.KC);
        };
        {
            char label[64];
            snprintf(label, sizeof(label), "numa_summa_%dx%d_kc%d", g_PX, g_PY, gp.KC);
            bench_one(fout, label, numa_summa_fn, A, B, Cchk, M, N, K, g_threads);
            verify(Cchk, Cref, M, N, label);
        }

        // 5b. NUMA-SUMMA with larger KC_summa panels
        int KC_summa_large = std::min(K, gp.KC * 4);
        if (KC_summa_large != gp.KC) {
            auto numa_summa_fn2 = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
                numa_summa_gemm(a, b, c, m, n, k, grid, gp, KC_summa_large);
            };
            char label[64];
            snprintf(label, sizeof(label), "numa_summa_%dx%d_kc%d", g_PX, g_PY, KC_summa_large);
            bench_one(fout, label, numa_summa_fn2, A, B, Cchk, M, N, K, g_threads);
            verify(Cchk, Cref, M, N, label);
        }

        // 6. NUMA-Cannon (only if square grid)
        if (g_PX == g_PY && K % (g_PX * BK) == 0) {
            auto numa_cannon_fn = [&](const IN_T* a, const IN_T* b, OUT_T* c, int m, int n, int k) {
                numa_cannon_gemm(a, b, c, m, n, k, grid, gp);
            };
            char label[64];
            snprintf(label, sizeof(label), "numa_cannon_%dx%d", g_PX, g_PY);
            bench_one(fout, label, numa_cannon_fn, A, B, Cchk, M, N, K, g_threads);
            verify(Cchk, Cref, M, N, label);
        }

        // Cleanup
        grid.free_buffers();
        typed_free(A, szA);
        typed_free(B, szB);
        typed_free(Cref, szC);
        typed_free(Cmkl, szC);
        typed_free(Cchk, szC);
        typed_free(g_A32_buf, szA);
        typed_free(g_B32_buf, szB);
        g_A32_buf = nullptr;
        g_B32_buf = nullptr;
        printf("\n");
    }

    fclose(fout);
    printf("All runs written to %s\n", g_outfile);
    return 0;
}