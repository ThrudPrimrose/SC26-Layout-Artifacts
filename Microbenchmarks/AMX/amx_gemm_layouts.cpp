// bench_amx_gemm_layouts.cpp
// g++ -O3 -fopenmp -march=sapphirerapids -mamx-tile -mamx-bf16 -mamx-int8 \
//     -I. -o bench2 bench_amx_gemm_layouts.cpp -lopenblas
//
// Thread dispatch:
//   Static  <BX,BY>:         BX*BY threads, 2D grid, each owns tiles_m/BX x tiles_n/BY tiles
//   Hypertile <BX,BY,TX,TY>: BX*BY threads, each computes TX*TY tiles per step,
//                             hypertile = (BX*TX)x(BY*TY), swept over C
#include <cblas.h>
#include <omp.h>
#include <sys/mman.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "core/amx/tile.hpp"

using namespace core_ir;

static int g_warmup  = 5;
static int g_iters   = 100;

#if !defined(_NTHREADS)
    #define _NTHREADS 32
#endif

static int g_threads = _NTHREADS;
static const char* g_outfile = "bench_amx_layouts.csv";

static constexpr int T = 32;
static constexpr int TSZ = T * T;

// ── helpers ─────────────────────────────────────────────────────────
static inline double now_ns() {
    const auto t = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::nano>(t.time_since_epoch()).count();
}

template <typename X>
static X* __restrict__ mmap_alloc(const size_t n) {
    void* const p = mmap(nullptr, n * sizeof(X), PROT_READ | PROT_WRITE,
                         MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    return static_cast<X*>(p);
}

template <typename X>
static void mmap_free(X* __restrict__ p, const size_t n) {
    munmap(p, n * sizeof(X));
}

static void rand_fill_f16(_Float16* __restrict__ p, const size_t n) {
    for (size_t i = 0; i < n; ++i) {
        p[i] = (_Float16)((float)rand() / (float)RAND_MAX - 0.5f);
    }
}

// ════════════════════════════════════════════════════════════════════
// FLAT reference
// ════════════════════════════════════════════════════════════════════

static void ref_flat(const _Float16* __restrict__ A, const _Float16* __restrict__ B,
                     float* __restrict__ C, const int M, const int N, const int K) {
    constexpr int BLK = 64;
    memset(C, 0, (size_t)M * N * sizeof(float));
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int bi = 0; bi < M; bi += BLK) {
        for (int bj = 0; bj < N; bj += BLK) {
            for (int bk = 0; bk < K; bk += BLK) {
                const int mi = (bi + BLK < M) ? bi + BLK : M;
                const int nj = (bj + BLK < N) ? bj + BLK : N;
                const int pk = (bk + BLK < K) ? bk + BLK : K;
                for (int i = bi; i < mi; ++i) {
                    for (int k = bk; k < pk; ++k) {
                        const float av = (float)A[i * K + k];
                        for (int j = bj; j < nj; ++j) {
                            C[i * N + j] += av * (float)B[k + K * j];
                        }
                    }
                }
            }
        }
    }
}

static float* __restrict__ g_A32;
static float* __restrict__ g_B32;

static void blas_flat(const _Float16* __restrict__ A, const _Float16* __restrict__ B,
                      float* __restrict__ C, const int M, const int N, const int K) {
    for (size_t i = 0; i < (size_t)M * K; ++i) { g_A32[i] = (float)A[i]; }
    for (size_t i = 0; i < (size_t)K * N; ++i) { g_B32[i] = (float)B[i]; }
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                M, N, K, 1.0f, g_A32, K, g_B32, K, 0.0f, C, N);
}

// ════════════════════════════════════════════════════════════════════
// TILED copyin/copyout
// ════════════════════════════════════════════════════════════════════

static void copyin_A(const _Float16* __restrict__ f, _Float16* __restrict__ t,
                     const int M, const int K) {
    const int tm_ = M / T, tk_ = K / T;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tm = 0; tm < tm_; ++tm) { for (int tk = 0; tk < tk_; ++tk) {
        _Float16* __restrict__ d = t + (size_t)(tm * tk_ + tk) * TSZ;
        for (int li = 0; li < T; ++li) { for (int lk = 0; lk < T; ++lk) {
            d[li * T + lk] = f[(tm * T + li) * K + (tk * T + lk)];
        }}
    }}
}

static void copyin_B(const _Float16* __restrict__ f, _Float16* __restrict__ t,
                     const int K, const int N) {
    const int tk_ = K / T, tn_ = N / T;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tn = 0; tn < tn_; ++tn) { for (int tk = 0; tk < tk_; ++tk) {
        _Float16* __restrict__ d = t + (size_t)(tn * tk_ + tk) * TSZ;
        for (int lk = 0; lk < T; ++lk) { for (int ln = 0; ln < T; ++ln) {
            d[lk * T + ln] = f[(tk * T + lk) + K * (tn * T + ln)];
        }}
    }}
}

static void copyin_C(const float* __restrict__ f, float* __restrict__ t,
                     const int M, const int N) {
    const int tm_ = M / T, tn_ = N / T;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tm = 0; tm < tm_; ++tm) { for (int tn = 0; tn < tn_; ++tn) {
        float* __restrict__ d = t + (size_t)(tm * tn_ + tn) * TSZ;
        for (int li = 0; li < T; ++li) { for (int ln = 0; ln < T; ++ln) {
            d[li * T + ln] = f[(tm * T + li) * N + (tn * T + ln)];
        }}
    }}
}

static void zero_C_tiled(float* __restrict__ t, const int M, const int N) {
    memset(t, 0, (size_t)M * N * sizeof(float));
}

static void copyout_A(_Float16* __restrict__ f, const _Float16* __restrict__ t,
                      const int M, const int K) {
    const int tm_ = M / T, tk_ = K / T;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tm = 0; tm < tm_; ++tm) { for (int tk = 0; tk < tk_; ++tk) {
        const _Float16* __restrict__ s = t + (size_t)(tm * tk_ + tk) * TSZ;
        for (int li = 0; li < T; ++li) { for (int lk = 0; lk < T; ++lk) {
            f[(tm * T + li) * K + (tk * T + lk)] = s[li * T + lk];
        }}
    }}
}

static void copyout_B(_Float16* __restrict__ f, const _Float16* __restrict__ t,
                      const int K, const int N) {
    const int tk_ = K / T, tn_ = N / T;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tn = 0; tn < tn_; ++tn) { for (int tk = 0; tk < tk_; ++tk) {
        const _Float16* __restrict__ s = t + (size_t)(tn * tk_ + tk) * TSZ;
        for (int lk = 0; lk < T; ++lk) { for (int ln = 0; ln < T; ++ln) {
            f[(tk * T + lk) + K * (tn * T + ln)] = s[lk * T + ln];
        }}
    }}
}

static void copyout_C(float* __restrict__ f, const float* __restrict__ t,
                      const int M, const int N) {
    const int tm_ = M / T, tn_ = N / T;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tm = 0; tm < tm_; ++tm) { for (int tn = 0; tn < tn_; ++tn) {
        const float* __restrict__ s = t + (size_t)(tm * tn_ + tn) * TSZ;
        for (int li = 0; li < T; ++li) { for (int ln = 0; ln < T; ++ln) {
            f[(tm * T + li) * N + (tn * T + ln)] = s[li * T + ln];
        }}
    }}
}

static void ref_tiled(const _Float16* __restrict__ At, const _Float16* __restrict__ Bt,
                      float* __restrict__ Ct, const int M, const int N, const int K) {
    const int tm_ = M / T, tn_ = N / T, tk_ = K / T;
    memset(Ct, 0, (size_t)M * N * sizeof(float));
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tm = 0; tm < tm_; ++tm) { for (int tn = 0; tn < tn_; ++tn) {
        float* __restrict__ ct = Ct + (size_t)(tm * tn_ + tn) * TSZ;
        for (int tk = 0; tk < tk_; ++tk) {
            const _Float16* __restrict__ at = At + (size_t)(tm * tk_ + tk) * TSZ;
            const _Float16* __restrict__ bt = Bt + (size_t)(tn * tk_ + tk) * TSZ;
            for (int li = 0; li < T; ++li) { for (int lk = 0; lk < T; ++lk) {
                const float av = (float)at[li * T + lk];
                for (int ln = 0; ln < T; ++ln) {
                    ct[li * T + ln] += av * (float)bt[lk * T + ln];
                }
            }}
        }
    }}
}

// ════════════════════════════════════════════════════════════════════
// Shared: packing routines
// ════════════════════════════════════════════════════════════════════

using IN_T  = _Float16;
using OUT_T = float;
constexpr int BM = 32, BN = 32, BK = 32;

using GEMM = TileOps<IN_T, IN_T, OUT_T, OUT_T,
                     MajorAxis::ROW, MajorAxis::ROW,
                     MajorAxis::ROW, MajorAxis::ROW, BM, BN, BK>;
static constexpr int A_TSZ = BM * BK;
static constexpr int B_TSZ = (BK / 2) * (2 * BN);

static void pack_A(const IN_T* __restrict__ A, IN_T* __restrict__ ap,
                   const int M, const int K) {
    const int tm_ = M / BM, tk_ = K / BK;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tm = 0; tm < tm_; ++tm) { for (int tk = 0; tk < tk_; ++tk) {
        IN_T* __restrict__ d = ap + (size_t)(tm * tk_ + tk) * A_TSZ;
        for (int li = 0; li < BM; ++li) { for (int lk = 0; lk < BK; ++lk) {
            d[li * BK + lk] = A[(tm * BM + li) * K + (tk * BK + lk)];
        }}
    }}
}

static void pack_B(const IN_T* __restrict__ B, IN_T* __restrict__ bp,
                   const int K, const int N) {
    const int tk_ = K / BK, tn_ = N / BN;
    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int tn = 0; tn < tn_; ++tn) { for (int tk = 0; tk < tk_; ++tk) {
        IN_T* __restrict__ d = bp + (size_t)(tn * tk_ + tk) * B_TSZ;
        for (int k = 0; k < BK; k += 2) { for (int n = 0; n < BN; ++n) {
            const int gk = tk * BK + k, gn = tn * BN + n;
            d[(k / 2) * (2 * BN) + 2 * n + 0] = B[gk + K * gn];
            d[(k / 2) * (2 * BN) + 2 * n + 1] = B[gk + 1 + K * gn];
        }}
    }}
}

static void repack_B_tiles(const IN_T* __restrict__ bt, IN_T* __restrict__ ba,
                           const int K, const int N) {
    const int tk_ = K / BK, tn_ = N / BN, tot = tn_ * tk_;
    #pragma omp parallel for num_threads(g_threads) schedule(static)
    for (int t = 0; t < tot; ++t) {
        const int tn = t / tk_, tk = t % tk_;
        const IN_T* __restrict__ s = bt + (size_t)(tn * tk_ + tk) * TSZ;
        IN_T* __restrict__ d = ba + (size_t)(tn * tk_ + tk) * B_TSZ;
        for (int k = 0; k < BK; k += 2) { for (int n = 0; n < BN; ++n) {
            d[(k / 2) * (2 * BN) + 2 * n + 0] = s[(k + 0) * BN + n];
            d[(k / 2) * (2 * BN) + 2 * n + 1] = s[(k + 1) * BN + n];
        }}
    }
}

// ════════════════════════════════════════════════════════════════════
// MMA body: single tile (tm,tn) over all K tiles
// ════════════════════════════════════════════════════════════════════

// Packed A/B -> flat C
static inline void mma_tile_flat(const IN_T* __restrict__ apk, const IN_T* __restrict__ bpk,
                                 OUT_T* __restrict__ C,
                                 const int tm, const int tn,
                                 const int tk_, const int N) {
    typename GEMM::FragmentD df{};
    for (int tk = 0; tk < tk_; ++tk) {
        const IN_T* ap = apk + (size_t)(tm * tk_ + tk) * A_TSZ;
        const IN_T* bp = bpk + (size_t)(tn * tk_ + tk) * B_TSZ;
        typename GEMM::FragmentA af; af.template load(ap);
        typename GEMM::FragmentB bf; bf.template load(bp);
        GEMM::MMA::mma(af, bf, df);
    }
    OUT_T* p = C + tm * BM * N + tn * BN;
    df.template unload_unpack<MajorAxis::ROW>(p, N);
}

// Packed A / repacked B -> tiled C
static inline void mma_tile_tiled(const IN_T* __restrict__ at, const IN_T* __restrict__ ba,
                                  OUT_T* __restrict__ ct,
                                  const int tm, const int tn,
                                  const int tk_, const int tn_) {
    typename GEMM::FragmentD df{};
    for (int tk = 0; tk < tk_; ++tk) {
        const IN_T* ap = at + (size_t)(tm * tk_ + tk) * A_TSZ;
        const IN_T* bp = ba + (size_t)(tn * tk_ + tk) * B_TSZ;
        typename GEMM::FragmentA af; af.template load(ap);
        typename GEMM::FragmentB bf; bf.template load(bp);
        GEMM::MMA::mma(af, bf, df);
    }
    OUT_T* p = ct + (size_t)(tm * tn_ + tn) * TSZ;
    df.template unload_unpack<MajorAxis::ROW>(p, BN);
}

// ════════════════════════════════════════════════════════════════════
// TX×TY register-blocked body
// ════════════════════════════════════════════════════════════════════

template <int TX, int TY>
static inline void mma_block_flat(const IN_T* __restrict__ apk, const IN_T* __restrict__ bpk,
                                  OUT_T* __restrict__ C,
                                  const int tm0, const int tn0,
                                  const int tk_, const int N) {
    typename GEMM::FragmentD df[TX][TY]{};
    for (int tk = 0; tk < tk_; ++tk) {
        typename GEMM::FragmentA af[TX];
        for (int tx = 0; tx < TX; ++tx) {
            const IN_T* ap = apk + (size_t)((tm0 + tx) * tk_ + tk) * A_TSZ;
            af[tx].template load(ap);
        }
        typename GEMM::FragmentB bf[TY];
        for (int ty = 0; ty < TY; ++ty) {
            const IN_T* bp = bpk + (size_t)((tn0 + ty) * tk_ + tk) * B_TSZ;
            bf[ty].template load(bp);
        }
        for (int tx = 0; tx < TX; ++tx) {
            for (int ty = 0; ty < TY; ++ty) {
                GEMM::MMA::mma(af[tx], bf[ty], df[tx][ty]);
            }
        }
    }
    for (int tx = 0; tx < TX; ++tx) {
        for (int ty = 0; ty < TY; ++ty) {
            OUT_T* p = C + (tm0 + tx) * BM * N + (tn0 + ty) * BN;
            df[tx][ty].template unload_unpack<MajorAxis::ROW>(p, N);
        }
    }
}

template <int TX, int TY>
static inline void mma_block_tiled(const IN_T* __restrict__ at, const IN_T* __restrict__ ba,
                                   OUT_T* __restrict__ ct,
                                   const int tm0, const int tn0,
                                   const int tk_, const int tn_) {
    typename GEMM::FragmentD df[TX][TY]{};
    for (int tk = 0; tk < tk_; ++tk) {
        typename GEMM::FragmentA af[TX];
        for (int tx = 0; tx < TX; ++tx) {
            const IN_T* ap = at + (size_t)((tm0 + tx) * tk_ + tk) * A_TSZ;
            af[tx].template load(ap);
        }
        typename GEMM::FragmentB bf[TY];
        for (int ty = 0; ty < TY; ++ty) {
            const IN_T* bp = ba + (size_t)((tn0 + ty) * tk_ + tk) * B_TSZ;
            bf[ty].template load(bp);
        }
        for (int tx = 0; tx < TX; ++tx) {
            for (int ty = 0; ty < TY; ++ty) {
                GEMM::MMA::mma(af[tx], bf[ty], df[tx][ty]);
            }
        }
    }
    for (int tx = 0; tx < TX; ++tx) {
        for (int ty = 0; ty < TY; ++ty) {
            OUT_T* p = ct + (size_t)((tm0 + tx) * tn_ + (tn0 + ty)) * TSZ;
            df[tx][ty].template unload_unpack<MajorAxis::ROW>(p, BN);
        }
    }
}

// On-the-fly: flat A(row-major) + B(col-major) -> stack pack -> MMA -> flat C
template <int TX, int TY>
static inline void mma_block_flat_otf(const IN_T* __restrict__ A, const IN_T* __restrict__ B,
                                      OUT_T* __restrict__ C,
                                      const int tm0, const int tn0,
                                      const int tk_, const int K, const int N) {
    alignas(64) IN_T ab[TX][A_TSZ];
    alignas(64) IN_T bb[TY][B_TSZ];
    typename GEMM::FragmentD df[TX][TY]{};
    for (int tk = 0; tk < tk_; ++tk) {
        typename GEMM::FragmentA af[TX];
        for (int tx = 0; tx < TX; ++tx) {
            for (int li = 0; li < BM; ++li) { for (int lk = 0; lk < BK; ++lk) {
                ab[tx][li * BK + lk] = A[((tm0 + tx) * BM + li) * K + (tk * BK + lk)];
            }}
            const IN_T* ap = ab[tx]; af[tx].template load(ap);
        }
        typename GEMM::FragmentB bf[TY];
        for (int ty = 0; ty < TY; ++ty) {
            for (int k = 0; k < BK; k += 2) { for (int n = 0; n < BN; ++n) {
                const int gk = tk * BK + k, gn = (tn0 + ty) * BN + n;
                bb[ty][(k / 2) * (2 * BN) + 2 * n + 0] = B[gk + K * gn];
                bb[ty][(k / 2) * (2 * BN) + 2 * n + 1] = B[gk + 1 + K * gn];
            }}
            const IN_T* bp = bb[ty]; bf[ty].template load(bp);
        }
        for (int tx = 0; tx < TX; ++tx) {
            for (int ty = 0; ty < TY; ++ty) {
                GEMM::MMA::mma(af[tx], bf[ty], df[tx][ty]);
            }
        }
    }
    for (int tx = 0; tx < TX; ++tx) {
        for (int ty = 0; ty < TY; ++ty) {
            OUT_T* p = C + (tm0 + tx) * BM * N + (tn0 + ty) * BN;
            df[tx][ty].template unload_unpack<MajorAxis::ROW>(p, N);
        }
    }
}

// On-the-fly: tiled A + raw tiled B -> stack repack B -> MMA -> tiled C
template <int TX, int TY>
static inline void mma_block_tiled_otf(const IN_T* __restrict__ at, const IN_T* __restrict__ bt,
                                       OUT_T* __restrict__ ct,
                                       const int tm0, const int tn0,
                                       const int tk_, const int tn_) {
    alignas(64) IN_T bb[TY][B_TSZ];
    typename GEMM::FragmentD df[TX][TY]{};
    for (int tk = 0; tk < tk_; ++tk) {
        typename GEMM::FragmentA af[TX];
        for (int tx = 0; tx < TX; ++tx) {
            const IN_T* ap = at + (size_t)((tm0 + tx) * tk_ + tk) * A_TSZ;
            af[tx].template load(ap);
        }
        typename GEMM::FragmentB bf[TY];
        for (int ty = 0; ty < TY; ++ty) {
            const IN_T* s = bt + (size_t)((tn0 + ty) * tk_ + tk) * TSZ;
            for (int k = 0; k < BK; k += 2) { for (int n = 0; n < BN; ++n) {
                bb[ty][(k / 2) * (2 * BN) + 2 * n + 0] = s[(k + 0) * BN + n];
                bb[ty][(k / 2) * (2 * BN) + 2 * n + 1] = s[(k + 1) * BN + n];
            }}
            const IN_T* bp = bb[ty]; bf[ty].template load(bp);
        }
        for (int tx = 0; tx < TX; ++tx) {
            for (int ty = 0; ty < TY; ++ty) {
                GEMM::MMA::mma(af[tx], bf[ty], df[tx][ty]);
            }
        }
    }
    for (int tx = 0; tx < TX; ++tx) {
        for (int ty = 0; ty < TY; ++ty) {
            OUT_T* p = ct + (size_t)((tm0 + tx) * tn_ + (tn0 + ty)) * TSZ;
            df[tx][ty].template unload_unpack<MajorAxis::ROW>(p, BN);
        }
    }
}

// ════════════════════════════════════════════════════════════════════
// Static kernels <BX,BY>: BX*BY cores, each owns tiles_m/BX x tiles_n/BY
// ════════════════════════════════════════════════════════════════════

// flat: pack + compute -> flat C
template <int BX, int BY>
static void flat_pk_S(const IN_T* __restrict__ A, const IN_T* __restrict__ B,
                      OUT_T* __restrict__ C, const int M, const int N, const int K,
                      IN_T* __restrict__ apk, IN_T* __restrict__ bpk) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    pack_A(A, apk, M, K); pack_B(B, bpk, K, N);
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        const int tm_per = tm_ / BX, tn_per = tn_ / BY;
        const int tm0 = bx * tm_per, tn0 = by * tn_per;
        for (int tm = tm0; tm < tm0 + tm_per; ++tm) {
            for (int tn = tn0; tn < tn0 + tn_per; ++tn) {
                mma_tile_flat(apk, bpk, C, tm, tn, tk_, N);
            }
        }
    }
}

// flat: compute only -> flat C
template <int BX, int BY>
static void flat_co_S(const IN_T* __restrict__ apk, const IN_T* __restrict__ bpk,
                      OUT_T* __restrict__ C, const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        const int tm_per = tm_ / BX, tn_per = tn_ / BY;
        const int tm0 = bx * tm_per, tn0 = by * tn_per;
        for (int tm = tm0; tm < tm0 + tm_per; ++tm) {
            for (int tn = tn0; tn < tn0 + tn_per; ++tn) {
                mma_tile_flat(apk, bpk, C, tm, tn, tk_, N);
            }
        }
    }
}

// flat: on-the-fly -> flat C
template <int BX, int BY>
static void flat_otf_S(const IN_T* __restrict__ A, const IN_T* __restrict__ B,
                       OUT_T* __restrict__ C, const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        const int tm_per = tm_ / BX, tn_per = tn_ / BY;
        const int tm0 = bx * tm_per, tn0 = by * tn_per;
        for (int tm = tm0; tm < tm0 + tm_per; ++tm) {
            for (int tn = tn0; tn < tn0 + tn_per; ++tn) {
                mma_block_flat_otf<1, 1>(A, B, C, tm, tn, tk_, K, N);
            }
        }
    }
}

// tiled: repack + compute -> tiled C
template <int BX, int BY>
static void tiled_rk_S(const IN_T* __restrict__ at, const IN_T* __restrict__ bt,
                       IN_T* __restrict__ ba, OUT_T* __restrict__ ct,
                       const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    repack_B_tiles(bt, ba, K, N);
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        const int tm_per = tm_ / BX, tn_per = tn_ / BY;
        const int tm0 = bx * tm_per, tn0 = by * tn_per;
        for (int tm = tm0; tm < tm0 + tm_per; ++tm) {
            for (int tn = tn0; tn < tn0 + tn_per; ++tn) {
                mma_tile_tiled(at, ba, ct, tm, tn, tk_, tn_);
            }
        }
    }
}

// tiled: compute only -> tiled C
template <int BX, int BY>
static void tiled_co_S(const IN_T* __restrict__ at, const IN_T* __restrict__ ba,
                       OUT_T* __restrict__ ct, const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        const int tm_per = tm_ / BX, tn_per = tn_ / BY;
        const int tm0 = bx * tm_per, tn0 = by * tn_per;
        for (int tm = tm0; tm < tm0 + tm_per; ++tm) {
            for (int tn = tn0; tn < tn0 + tn_per; ++tn) {
                mma_tile_tiled(at, ba, ct, tm, tn, tk_, tn_);
            }
        }
    }
}

// tiled: on-the-fly -> tiled C
template <int BX, int BY>
static void tiled_otf_S(const IN_T* __restrict__ at, const IN_T* __restrict__ bt,
                        OUT_T* __restrict__ ct, const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        const int tm_per = tm_ / BX, tn_per = tn_ / BY;
        const int tm0 = bx * tm_per, tn0 = by * tn_per;
        for (int tm = tm0; tm < tm0 + tm_per; ++tm) {
            for (int tn = tn0; tn < tn0 + tn_per; ++tn) {
                mma_block_tiled_otf<1, 1>(at, bt, ct, tm, tn, tk_, tn_);
            }
        }
    }
}

// ════════════════════════════════════════════════════════════════════
// Hypertile kernels <BX,BY,TX,TY>:
//   BX*BY cores, each computes TX*TY tiles per step
//   Hypertile = (BX*TX) x (BY*TY) tiles, swept over C
// ════════════════════════════════════════════════════════════════════

template <int BX, int BY, int TX, int TY>
static void flat_pk_H(const IN_T* __restrict__ A, const IN_T* __restrict__ B,
                      OUT_T* __restrict__ C, const int M, const int N, const int K,
                      IN_T* __restrict__ apk, IN_T* __restrict__ bpk) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    pack_A(A, apk, M, K); pack_B(B, bpk, K, N);
    constexpr int HM = BX * TX, HN = BY * TY;
    const int nhm = tm_ / HM, nhn = tn_ / HN;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        for (int hm = 0; hm < nhm; ++hm) {
            for (int hn = 0; hn < nhn; ++hn) {
                const int tm0 = hm * HM + bx * TX;
                const int tn0 = hn * HN + by * TY;
                mma_block_flat<TX, TY>(apk, bpk, C, tm0, tn0, tk_, N);
            }
        }
    }
}

template <int BX, int BY, int TX, int TY>
static void flat_co_H(const IN_T* __restrict__ apk, const IN_T* __restrict__ bpk,
                      OUT_T* __restrict__ C, const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    constexpr int HM = BX * TX, HN = BY * TY;
    const int nhm = tm_ / HM, nhn = tn_ / HN;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        for (int hm = 0; hm < nhm; ++hm) {
            for (int hn = 0; hn < nhn; ++hn) {
                const int tm0 = hm * HM + bx * TX;
                const int tn0 = hn * HN + by * TY;
                mma_block_flat<TX, TY>(apk, bpk, C, tm0, tn0, tk_, N);
            }
        }
    }
}

template <int BX, int BY, int TX, int TY>
static void flat_otf_H(const IN_T* __restrict__ A, const IN_T* __restrict__ B,
                       OUT_T* __restrict__ C, const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    constexpr int HM = BX * TX, HN = BY * TY;
    const int nhm = tm_ / HM, nhn = tn_ / HN;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        for (int hm = 0; hm < nhm; ++hm) {
            for (int hn = 0; hn < nhn; ++hn) {
                const int tm0 = hm * HM + bx * TX;
                const int tn0 = hn * HN + by * TY;
                mma_block_flat_otf<TX, TY>(A, B, C, tm0, tn0, tk_, K, N);
            }
        }
    }
}

template <int BX, int BY, int TX, int TY>
static void tiled_rk_H(const IN_T* __restrict__ at, const IN_T* __restrict__ bt,
                       IN_T* __restrict__ ba, OUT_T* __restrict__ ct,
                       const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    repack_B_tiles(bt, ba, K, N);
    constexpr int HM = BX * TX, HN = BY * TY;
    const int nhm = tm_ / HM, nhn = tn_ / HN;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        for (int hm = 0; hm < nhm; ++hm) {
            for (int hn = 0; hn < nhn; ++hn) {
                const int tm0 = hm * HM + bx * TX;
                const int tn0 = hn * HN + by * TY;
                mma_block_tiled<TX, TY>(at, ba, ct, tm0, tn0, tk_, tn_);
            }
        }
    }
}

template <int BX, int BY, int TX, int TY>
static void tiled_co_H(const IN_T* __restrict__ at, const IN_T* __restrict__ ba,
                       OUT_T* __restrict__ ct, const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    constexpr int HM = BX * TX, HN = BY * TY;
    const int nhm = tm_ / HM, nhn = tn_ / HN;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        for (int hm = 0; hm < nhm; ++hm) {
            for (int hn = 0; hn < nhn; ++hn) {
                const int tm0 = hm * HM + bx * TX;
                const int tn0 = hn * HN + by * TY;
                mma_block_tiled<TX, TY>(at, ba, ct, tm0, tn0, tk_, tn_);
            }
        }
    }
}

template <int BX, int BY, int TX, int TY>
static void tiled_otf_H(const IN_T* __restrict__ at, const IN_T* __restrict__ bt,
                        OUT_T* __restrict__ ct, const int M, const int N, const int K) {
    const int tm_ = M / BM, tn_ = N / BN, tk_ = K / BK;
    constexpr int HM = BX * TX, HN = BY * TY;
    const int nhm = tm_ / HM, nhn = tn_ / HN;
    #pragma omp parallel num_threads(BX * BY) proc_bind(close)
    {
        GEMM::init();
        const int tid = omp_get_thread_num();
        const int bx = tid / BY, by = tid % BY;
        for (int hm = 0; hm < nhm; ++hm) {
            for (int hn = 0; hn < nhn; ++hn) {
                const int tm0 = hm * HM + bx * TX;
                const int tn0 = hn * HN + by * TY;
                mma_block_tiled_otf<TX, TY>(at, bt, ct, tm0, tn0, tk_, tn_);
            }
        }
    }
}

// ════════════════════════════════════════════════════════════════════
// Verification + benchmark
// ════════════════════════════════════════════════════════════════════

static bool check_f32(const float* __restrict__ C, const float* __restrict__ R,
                      const int M, const int N) {
    for (int i = 0; i < M * N; ++i) {
        const float d = fabsf(C[i] - R[i]), den = fabsf(R[i]);
        if ((den > 1e-6f ? d / den : d) > 0.01f && d > 0.01f) { return false; }
    }
    return true;
}

static bool verify_rt_f16(const _Float16* o, const _Float16* r, size_t n, const char* l) {
    int bad = 0;
    for (size_t i = 0; i < n; ++i) { if (o[i] != r[i]) { if (!bad) { printf("  [FAIL] %s rt idx=%zu\n", l, i); } ++bad; } }
    printf(bad ? "  [FAIL] %s rt: %d\n" : "  [PASS] %s rt\n", l, bad); return !bad;
}

static bool verify_rt_f32(const float* o, const float* r, size_t n, const char* l) {
    int bad = 0;
    for (size_t i = 0; i < n; ++i) { if (o[i] != r[i]) { if (!bad) { printf("  [FAIL] %s rt idx=%zu\n", l, i); } ++bad; } }
    printf(bad ? "  [FAIL] %s rt: %d\n" : "  [PASS] %s rt\n", l, bad); return !bad;
}

// ── kernel table ────────────────────────────────────────────────────
using KernelFn = void(*)(int, int, int);
struct KernelEntry { std::string name; KernelFn fn; bool tiled_out; int req_m; int req_n; int nthreads; };

// Globals
static _Float16* __restrict__ g_A;     static _Float16* __restrict__ g_B;
static float* __restrict__ g_C;        static _Float16* __restrict__ g_At;
static _Float16* __restrict__ g_Bt;    static float* __restrict__ g_Ct;
static _Float16* __restrict__ g_apk;   static _Float16* __restrict__ g_bpk;
static _Float16* __restrict__ g_b_amx;
static float* __restrict__ g_Cref;     static float* __restrict__ g_Cchk;
static int g_M, g_N, g_K;

// ── wrappers (thunks to typed functions via globals) ────────────────
static void w_ref_flat(int M, int N, int K)  { ref_flat(g_A, g_B, g_C, M, N, K); }
static void w_blas_flat(int M, int N, int K) { blas_flat(g_A, g_B, g_C, M, N, K); }
static void w_ref_tiled(int M, int N, int K) { ref_tiled(g_At, g_Bt, g_Ct, M, N, K); }

// Default flat: pack + compute, flat in/out
static void w_default_flat(int M, int N, int K) {
    flat_pk_S<8, 4>(g_A, g_B, g_C, M, N, K, g_apk, g_bpk);
}
// Default tiled: copyin + repack + compute + copyout, flat in/out
static void w_default_tiled(int M, int N, int K) {
    copyin_A(g_A, g_At, M, K); copyin_B(g_B, g_Bt, K, N);
    zero_C_tiled(g_Ct, M, N);
    tiled_rk_S<8, 4>(g_At, g_Bt, g_b_amx, g_Ct, M, N, K);
    copyout_C(g_C, g_Ct, M, N);
}

// Static wrappers
#define W_FLAT_PK_S(BX, BY)  static void w_flat_pk_S_##BX##x##BY(int M,int N,int K) { flat_pk_S<BX,BY>(g_A,g_B,g_C,M,N,K,g_apk,g_bpk); }
#define W_FLAT_CO_S(BX, BY)  static void w_flat_co_S_##BX##x##BY(int M,int N,int K) { flat_co_S<BX,BY>(g_apk,g_bpk,g_C,M,N,K); }
#define W_FLAT_OTF_S(BX, BY) static void w_flat_otf_S_##BX##x##BY(int M,int N,int K) { flat_otf_S<BX,BY>(g_A,g_B,g_C,M,N,K); }
#define W_TILED_RK_S(BX, BY) static void w_tiled_rk_S_##BX##x##BY(int M,int N,int K) { tiled_rk_S<BX,BY>(g_At,g_Bt,g_b_amx,g_Ct,M,N,K); }
#define W_TILED_CO_S(BX, BY) static void w_tiled_co_S_##BX##x##BY(int M,int N,int K) { tiled_co_S<BX,BY>(g_At,g_b_amx,g_Ct,M,N,K); }
#define W_TILED_OTF_S(BX, BY) static void w_tiled_otf_S_##BX##x##BY(int M,int N,int K) { tiled_otf_S<BX,BY>(g_At,g_Bt,g_Ct,M,N,K); }

#define ALL_S_WRAPPERS(BX, BY) \
    W_FLAT_PK_S(BX,BY) W_FLAT_CO_S(BX,BY) W_FLAT_OTF_S(BX,BY) \
    W_TILED_RK_S(BX,BY) W_TILED_CO_S(BX,BY) W_TILED_OTF_S(BX,BY)

// Hypertile wrappers
#define W_FLAT_PK_H(BX,BY,TX,TY)  static void w_flat_pk_H_##BX##x##BY##_##TX##x##TY(int M,int N,int K) { flat_pk_H<BX,BY,TX,TY>(g_A,g_B,g_C,M,N,K,g_apk,g_bpk); }
#define W_FLAT_CO_H(BX,BY,TX,TY)  static void w_flat_co_H_##BX##x##BY##_##TX##x##TY(int M,int N,int K) { flat_co_H<BX,BY,TX,TY>(g_apk,g_bpk,g_C,M,N,K); }
#define W_FLAT_OTF_H(BX,BY,TX,TY) static void w_flat_otf_H_##BX##x##BY##_##TX##x##TY(int M,int N,int K) { flat_otf_H<BX,BY,TX,TY>(g_A,g_B,g_C,M,N,K); }
#define W_TILED_RK_H(BX,BY,TX,TY) static void w_tiled_rk_H_##BX##x##BY##_##TX##x##TY(int M,int N,int K) { tiled_rk_H<BX,BY,TX,TY>(g_At,g_Bt,g_b_amx,g_Ct,M,N,K); }
#define W_TILED_CO_H(BX,BY,TX,TY) static void w_tiled_co_H_##BX##x##BY##_##TX##x##TY(int M,int N,int K) { tiled_co_H<BX,BY,TX,TY>(g_At,g_b_amx,g_Ct,M,N,K); }
#define W_TILED_OTF_H(BX,BY,TX,TY) static void w_tiled_otf_H_##BX##x##BY##_##TX##x##TY(int M,int N,int K) { tiled_otf_H<BX,BY,TX,TY>(g_At,g_Bt,g_Ct,M,N,K); }

#define ALL_H_WRAPPERS(BX,BY,TX,TY) \
    W_FLAT_PK_H(BX,BY,TX,TY) W_FLAT_CO_H(BX,BY,TX,TY) W_FLAT_OTF_H(BX,BY,TX,TY) \
    W_TILED_RK_H(BX,BY,TX,TY) W_TILED_CO_H(BX,BY,TX,TY) W_TILED_OTF_H(BX,BY,TX,TY)

// Instantiate all BX*BY static + hypertile(1x1) for thread counts 1,2,4,8,16,32
#define ALL_CONFIGS(BX, BY) ALL_S_WRAPPERS(BX, BY) ALL_H_WRAPPERS(BX, BY, 1, 1)

// 1 thread
ALL_CONFIGS(1, 1)
// 2 threads
ALL_CONFIGS(2, 1) ALL_CONFIGS(1, 2)
// 4 threads
ALL_CONFIGS(4, 1) ALL_CONFIGS(2, 2) ALL_CONFIGS(1, 4)
// 8 threads
ALL_CONFIGS(8, 1) ALL_CONFIGS(4, 2) ALL_CONFIGS(2, 4) ALL_CONFIGS(1, 8)
// 16 threads
ALL_CONFIGS(16, 1) ALL_CONFIGS(8, 2) ALL_CONFIGS(4, 4) ALL_CONFIGS(2, 8) ALL_CONFIGS(1, 16)
// 32 threads
ALL_CONFIGS(32, 1) ALL_CONFIGS(16, 2) ALL_CONFIGS(8, 4) ALL_CONFIGS(4, 8) ALL_CONFIGS(2, 16) ALL_CONFIGS(1, 32)

// ── registration macros ─────────────────────────────────────────────
#define REG_S(V, BX, BY) \
    V.push_back({"flat_pk_S_"  #BX "x" #BY, w_flat_pk_S_##BX##x##BY, false, BX*T, BY*T, BX*BY}); \
    V.push_back({"flat_co_S_"  #BX "x" #BY, w_flat_co_S_##BX##x##BY, false, BX*T, BY*T, BX*BY}); \
    V.push_back({"flat_otf_S_" #BX "x" #BY, w_flat_otf_S_##BX##x##BY, false, BX*T, BY*T, BX*BY}); \
    V.push_back({"tiled_rk_S_"  #BX "x" #BY, w_tiled_rk_S_##BX##x##BY, true, BX*T, BY*T, BX*BY}); \
    V.push_back({"tiled_co_S_"  #BX "x" #BY, w_tiled_co_S_##BX##x##BY, true, BX*T, BY*T, BX*BY}); \
    V.push_back({"tiled_otf_S_" #BX "x" #BY, w_tiled_otf_S_##BX##x##BY, true, BX*T, BY*T, BX*BY});

#define REG_H(V, BX, BY, TX, TY) \
    V.push_back({"flat_pk_H"  #BX "x" #BY "_" #TX "x" #TY, w_flat_pk_H_##BX##x##BY##_##TX##x##TY, false, BX*TX*T, BY*TY*T, BX*BY}); \
    V.push_back({"flat_co_H"  #BX "x" #BY "_" #TX "x" #TY, w_flat_co_H_##BX##x##BY##_##TX##x##TY, false, BX*TX*T, BY*TY*T, BX*BY}); \
    V.push_back({"flat_otf_H" #BX "x" #BY "_" #TX "x" #TY, w_flat_otf_H_##BX##x##BY##_##TX##x##TY, false, BX*TX*T, BY*TY*T, BX*BY}); \
    V.push_back({"tiled_rk_H"  #BX "x" #BY "_" #TX "x" #TY, w_tiled_rk_H_##BX##x##BY##_##TX##x##TY, true, BX*TX*T, BY*TY*T, BX*BY}); \
    V.push_back({"tiled_co_H"  #BX "x" #BY "_" #TX "x" #TY, w_tiled_co_H_##BX##x##BY##_##TX##x##TY, true, BX*TX*T, BY*TY*T, BX*BY}); \
    V.push_back({"tiled_otf_H" #BX "x" #BY "_" #TX "x" #TY, w_tiled_otf_H_##BX##x##BY##_##TX##x##TY, true, BX*TX*T, BY*TY*T, BX*BY});

#define REG_CONFIGS(V, BX, BY) REG_S(V, BX, BY) REG_H(V, BX, BY, 1, 1)

static std::vector<KernelEntry> build_kernels() {
    std::vector<KernelEntry> v;
    v.push_back({"ref_flat",      w_ref_flat,      false, T, T, 0});
    v.push_back({"blas_flat",     w_blas_flat,     false, T, T, 0});
    v.push_back({"default_flat",  w_default_flat,  false, 8*T, 4*T, 32});
    v.push_back({"default_tiled", w_default_tiled, false, 8*T, 4*T, 32});
    v.push_back({"ref_tiled",     w_ref_tiled,     true,  T, T, 0});

    // 1 thread
    REG_CONFIGS(v, 1, 1)
    // 2 threads
    REG_CONFIGS(v, 2, 1) REG_CONFIGS(v, 1, 2)
    // 4 threads
    REG_CONFIGS(v, 4, 1) REG_CONFIGS(v, 2, 2) REG_CONFIGS(v, 1, 4)
    // 8 threads
    REG_CONFIGS(v, 8, 1) REG_CONFIGS(v, 4, 2) REG_CONFIGS(v, 2, 4) REG_CONFIGS(v, 1, 8)
    // 16 threads
    REG_CONFIGS(v, 16, 1) REG_CONFIGS(v, 8, 2) REG_CONFIGS(v, 4, 4) REG_CONFIGS(v, 2, 8) REG_CONFIGS(v, 1, 16)
    // 32 threads
    REG_CONFIGS(v, 32, 1) REG_CONFIGS(v, 16, 2) REG_CONFIGS(v, 8, 4) REG_CONFIGS(v, 4, 8) REG_CONFIGS(v, 2, 16) REG_CONFIGS(v, 1, 32)

    return v;
}

// ── per-kernel correctness + benchmark ──────────────────────────────
static bool check_kernel(const KernelEntry& ke) {
    const int M = g_M, N = g_N, K = g_K;
    if (ke.tiled_out) { zero_C_tiled(g_Ct, M, N); }
    else { memset(g_C, 0, (size_t)M * N * sizeof(float)); }

    ke.fn(M, N, K);

    const float* result;
    if (ke.tiled_out) { copyout_C(g_Cchk, g_Ct, M, N); result = g_Cchk; }
    else { result = g_C; }

    if (check_f32(result, g_Cref, M, N)) { return true; }

    for (int i = 0; i < M * N; ++i) {
        const float d = fabsf(result[i] - g_Cref[i]), den = fabsf(g_Cref[i]);
        if ((den > 1e-6f ? d / den : d) > 0.01f && d > 0.01f) {
            printf("    [FAIL] %-36s idx=%d got=%.4f ref=%.4f\n",
                   ke.name.c_str(), i, result[i], g_Cref[i]);
            break;
        }
    }
    return false;
}

static void bench_kernel(FILE* __restrict__ f, const KernelEntry& ke,
                        const std::vector<int>& thread_filter) {
    const int M = g_M, N = g_N, K = g_K;
    const int nt = ke.nthreads ? ke.nthreads : g_threads;

    // Skip if thread count not in filter (empty filter = run all)
    if (!thread_filter.empty()) {
        bool found = false;
        for (int tf : thread_filter) {
            if (tf == nt) { found = true; break; }
        }
        if (!found) { return; }
    }

    // Skip if matrix dimensions not divisible by this kernel's tile/thread config
    if (M % ke.req_m != 0 || N % ke.req_n != 0) {
        printf("    %-36s [SKIP] t=%2d M=%d%%req_m=%d N=%d%%req_n=%d\n",
               ke.name.c_str(), nt, M, ke.req_m, N, ke.req_n);
        for (int i = 0; i < g_iters; ++i) {
            fprintf(f, "%s,%d,%d,%d,%d,%d,-2.00,-2.000000\n", ke.name.c_str(), M, N, K, nt, i);
        }
        return;
    }

    if (!check_kernel(ke)) {
        for (int i = 0; i < g_iters; ++i) {
            fprintf(f, "%s,%d,%d,%d,%d,%d,-1.00,-1.000000\n", ke.name.c_str(), M, N, K, nt, i);
        }
        return;
    }
    for (int w = 0; w < g_warmup; ++w) { ke.fn(M, N, K); }
    const double flops = 2.0 * (double)M * N * K;
    double total = 0; std::vector<double> times(g_iters);
    for (int i = 0; i < g_iters; ++i) {
        const double t0 = now_ns(); ke.fn(M, N, K); const double t1 = now_ns();
        times[i] = t1 - t0; total += times[i];
    }
    const double avg = total / g_iters;
    std::vector<double> s(times); std::sort(s.begin(), s.end());
    const double med = (g_iters % 2) ? s[g_iters / 2] : (s[g_iters / 2 - 1] + s[g_iters / 2]) / 2.0;
    for (int i = 0; i < g_iters; ++i) {
        fprintf(f, "%s,%d,%d,%d,%d,%d,%.2f,%.6f\n", ke.name.c_str(), M, N, K, nt, i, times[i], flops / times[i]);
    }
    printf("    %-36s [PASS] t=%2d avg=%.0f med=%.0f %.4f GF/s\n", ke.name.c_str(), nt, avg, med, flops / avg);
}

// ── arg parsing ─────────────────────────────────────────────────────
struct MNK { int M, N, K; };
static std::vector<int> g_thread_filter; // empty = run all

static void usage(const char* p) {
    printf("Usage: %s M N K [...] [-w W] [-i I] [-t T1,T2,...] [-o F]\n", p);
    printf("  -t: comma-separated thread counts to benchmark (default: all)\n");
    printf("      e.g. -t 1,4,32 or -t 32\n");
}

static std::vector<int> parse_thread_list(const char* s) {
    std::vector<int> v;
    const char* p = s;
    while (*p) {
        v.push_back(atoi(p));
        while (*p && *p != ',') { ++p; }
        if (*p == ',') { ++p; }
    }
    return v;
}

static std::vector<MNK> parse_args(int argc, char** argv) {
    std::vector<MNK> sz; std::vector<int> nums;
    for (int i = 1; i < argc; ++i) {
        if      (!strcmp(argv[i], "-w") && i+1<argc) { g_warmup  = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-i") && i+1<argc) { g_iters   = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-t") && i+1<argc) { g_thread_filter = parse_thread_list(argv[++i]); }
        else if (!strcmp(argv[i], "-o") && i+1<argc) { g_outfile = argv[++i]; }
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(argv[0]); exit(0); }
        else { nums.push_back(atoi(argv[i])); }
    }
    if (nums.size() % 3 || nums.empty()) { fprintf(stderr, "Need M N K triples\n"); exit(1); }
    for (size_t i = 0; i < nums.size(); i += 3) { sz.push_back({nums[i], nums[i+1], nums[i+2]}); }
    // Set g_threads to max of filter for base kernels
    if (!g_thread_filter.empty()) {
        g_threads = *std::max_element(g_thread_filter.begin(), g_thread_filter.end());
    }
    return sz;
}

// ── main ────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 1; }
    const auto sizes = parse_args(argc, argv);
    const auto kernels = build_kernels();
    printf("Config: warmup=%d iters=%d threads=%d output=%s kernels=%zu",
           g_warmup, g_iters, g_threads, g_outfile, kernels.size());
    if (!g_thread_filter.empty()) {
        printf(" filter=[");
        for (size_t i = 0; i < g_thread_filter.size(); ++i) {
            printf("%s%d", i ? "," : "", g_thread_filter[i]);
        }
        printf("]");
    }
    printf("\n\n");

    FILE* fout = fopen(g_outfile, "w");
    if (!fout) { perror("fopen"); return 1; }
    fprintf(fout, "kernel,M,N,K,threads,run,time_ns,gflops\n");

    for (size_t c = 0; c < sizes.size(); ++c) {
        const int M = sizes[c].M, N = sizes[c].N, K = sizes[c].K;
        if (M % T || N % T || K % T) { printf("SKIP %dx%dx%d\n\n", M, N, K); continue; }
        printf("=== M=%d N=%d K=%d ===\n", M, N, K);
        g_M = M; g_N = N; g_K = K;
        const size_t szA = (size_t)M*K, szB = (size_t)K*N, szC = (size_t)M*N;

        g_A = mmap_alloc<_Float16>(szA); g_B = mmap_alloc<_Float16>(szB);
        g_C = mmap_alloc<float>(szC); g_A32 = mmap_alloc<float>(szA); g_B32 = mmap_alloc<float>(szB);
        g_Cref = mmap_alloc<float>(szC); g_Cchk = mmap_alloc<float>(szC);

        #pragma omp parallel num_threads(g_threads) proc_bind(close)
        { int t=omp_get_thread_num(),n=omp_get_num_threads();
          size_t l=t*szA/n,h=(t+1)*szA/n; rand_fill_f16(g_A+l,h-l);
          l=t*szB/n; h=(t+1)*szB/n; rand_fill_f16(g_B+l,h-l); }

        g_At = mmap_alloc<_Float16>(szA); g_Bt = mmap_alloc<_Float16>(szB); g_Ct = mmap_alloc<float>(szC);
        copyin_A(g_A, g_At, M, K); copyin_B(g_B, g_Bt, K, N); zero_C_tiled(g_Ct, M, N);
        g_apk = mmap_alloc<_Float16>(szA); g_bpk = mmap_alloc<_Float16>(szB);
        g_b_amx = mmap_alloc<_Float16>((size_t)(N/32)*(K/32)*B_TSZ);
        pack_A(g_A, g_apk, M, K); pack_B(g_B, g_bpk, K, N);
        repack_B_tiles(g_Bt, g_b_amx, K, N);

        printf("  Roundtrips:\n");
        { auto t=mmap_alloc<_Float16>(szA); copyout_A(t,g_At,M,K); verify_rt_f16(g_A,t,szA,"A"); mmap_free(t,szA); }
        { auto t=mmap_alloc<_Float16>(szB); copyout_B(t,g_Bt,K,N); verify_rt_f16(g_B,t,szB,"B"); mmap_free(t,szB); }
        { auto o=mmap_alloc<float>(szC); auto ti=mmap_alloc<float>(szC); auto b=mmap_alloc<float>(szC);
          for(size_t i=0;i<szC;++i){o[i]=(float)rand()/(float)RAND_MAX-0.5f;}
          copyin_C(o,ti,M,N); copyout_C(b,ti,M,N); verify_rt_f32(o,b,szC,"C");
          mmap_free(o,szC); mmap_free(ti,szC); mmap_free(b,szC); }

        ref_flat(g_A, g_B, g_Cref, M, N, K);

        printf("  Kernels (%d warmup, %d iters, %zu total):\n", g_warmup, g_iters, kernels.size());
        for (const auto& ke : kernels) { bench_kernel(fout, ke, g_thread_filter); }

        mmap_free(g_A,szA); mmap_free(g_B,szB); mmap_free(g_C,szC);
        mmap_free(g_A32,szA); mmap_free(g_B32,szB);
        mmap_free(g_Cref,szC); mmap_free(g_Cchk,szC);
        mmap_free(g_At,szA); mmap_free(g_Bt,szB); mmap_free(g_Ct,szC);
        mmap_free(g_apk,szA); mmap_free(g_bpk,szB);
        mmap_free(g_b_amx,(size_t)(N/32)*(K/32)*B_TSZ);
        printf("\n");
    }
    fclose(fout);
    printf("All runs written to %s\n", g_outfile);
    return 0;
}