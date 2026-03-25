// bench_amx_gemm.cpp
// Usage: ./bench M1 N1 K1 [M2 N2 K2 ...] [-w warmup] [-i iters] [-t threads] [-o file.csv]
// g++ -O3 -fopenmp -march=sapphirerapids -mamx-tile -mamx-bf16 -mamx-int8
//     -I. -I${MKLROOT}/include -o bench bench_amx_gemm.cpp
//     -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm
#include <mkl.h>
#include <omp.h>
#include <sys/mman.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "core/amx/tile.hpp"

using namespace core_ir;

// ── defaults ────────────────────────────────────────────────────────
static int g_warmup  = 5;
static int g_iters   = 100;
#if !defined(_NTHREADS)
    #define _NTHREADS 32
#endif

static int g_threads = _NTHREADS;
static const char* g_outfile = "bench_amx_gemm.csv";

// ── timing ──────────────────────────────────────────────────────────
static inline double now_ns() {
    auto t = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::nano>(t.time_since_epoch()).count();
}

// ── NUMA-safe mmap alloc ────────────────────────────────────────────
template <typename T>
static T* mmap_alloc(size_t n) {
    void* p = mmap(nullptr, n * sizeof(T), PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    return static_cast<T*>(p);
}

template <typename T>
static void mmap_free(T* p, size_t n) {
    munmap(p, n * sizeof(T));
}

// ── random init ─────────────────────────────────────────────────────
static void rand_fill_f16(_Float16* p, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        p[i] = (_Float16)((float)rand() / (float)RAND_MAX - 0.5f);
    }
}

// ── reference: OpenMP 2D block-tiled scalar GEMM (fp16 in, fp32 out) ──
static void ref_gemm_omp(const _Float16* A, const _Float16* B, float* C,
                         int M, int N, int K, bool b_row_major) {
    constexpr int BLK = 64;
    memset(C, 0, (size_t)M * N * sizeof(float));

    #pragma omp parallel for collapse(2) num_threads(g_threads) schedule(static)
    for (int bi = 0; bi < M; bi += BLK) {
        for (int bj = 0; bj < N; bj += BLK) {
            for (int bk = 0; bk < K; bk += BLK) {
                int mi = (bi + BLK < M) ? bi + BLK : M;
                int nj = (bj + BLK < N) ? bj + BLK : N;
                int pk = (bk + BLK < K) ? bk + BLK : K;
                for (int i = bi; i < mi; ++i) {
                    for (int k = bk; k < pk; ++k) {
                        float a_val = (float)A[i * K + k];
                        for (int j = bj; j < nj; ++j) {
                            float b_val = b_row_major
                                ? (float)B[k * N + j]
                                : (float)B[j * K + k];
                            C[i * N + j] += a_val * b_val;
                        }
                    }
                }
            }
        }
    }
}

// ── MKL sgemm (upcast fp16->fp32) ──────────────────────────────────
static void mkl_sgemm_wrap(const _Float16* A16, const _Float16* B16, float* C,
                           int M, int N, int K, bool b_row_major,
                           float* A32_buf, float* B32_buf) {
    mkl_set_num_threads(g_threads);
    for (size_t i = 0; i < (size_t)M * K; ++i) {
        A32_buf[i] = (float)A16[i];
    }
    for (size_t i = 0; i < (size_t)K * N; ++i) {
        B32_buf[i] = (float)B16[i];
    }
    if (b_row_major) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    M, N, K, 1.0f, A32_buf, K, B32_buf, N, 0.0f, C, N);
    } else {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    M, N, K, 1.0f, A32_buf, K, B32_buf, K, 0.0f, C, N);
    }
}

// ── MKL f16 GEMM (native fp16 in, fp32 out) ────────────────────────
static void mkl_f16_wrap(const _Float16* A, const _Float16* B, float* C,
                         int M, int N, int K, bool b_row_major) {
    mkl_set_num_threads(g_threads);
    if (b_row_major) {
        cblas_gemm_f16f16f32(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                             M, N, K, 1.0f,
                             (const MKL_F16*)A, K,
                             (const MKL_F16*)B, N,
                             0.0f, C, N);
    } else {
        cblas_gemm_f16f16f32(CblasRowMajor, CblasNoTrans, CblasTrans,
                             M, N, K, 1.0f,
                             (const MKL_F16*)A, K,
                             (const MKL_F16*)B, K,
                             0.0f, C, N);
    }
}

// ── AMX GEMM (hypertile + static variants) ─────────────────────────
template <typename IN_T, typename OUT_T, bool R2, int BM, int BN, int BK>
class amx_gemm_omp {
public:
    using GEMM = TileOps<IN_T, IN_T, OUT_T, OUT_T, MajorAxis::ROW,
                         R2 ? MajorAxis::ROW : MajorAxis::COLUMN,
                         MajorAxis::ROW, MajorAxis::ROW, BM, BN, BK>;

    static constexpr int HM = 8;
    static constexpr int HN = 4;
    static constexpr int a_tile_sz = BM * BK;
    static constexpr int b_tile_sz = (BK / 2) * (2 * BN);

    static void pack_parallel(const IN_T* src_a, const IN_T* src_b,
                              IN_T* a_pk, IN_T* b_pk,
                              int M, int N, int K) {
        int tiles_m = M / BM, tiles_n = N / BN, tiles_k = K / BK;
        typename GEMM::TileLayoutA la(src_a, M, K);
        typename GEMM::TileLayoutB lb(src_b, N, K);

        #pragma omp parallel num_threads(g_threads) proc_bind(close)
        {
            int tid = omp_get_thread_num();
            for (int tm = tid; tm < tiles_m; tm += g_threads) {
                for (int tk = 0; tk < tiles_k; ++tk) {
                    IN_T* dst = a_pk + (size_t)tm * tiles_k * a_tile_sz
                                     + (size_t)tk * a_tile_sz;
                    la.pack(tm, tk, dst);
                }
            }
            if constexpr (R2) {
                for (int tn = tid; tn < tiles_n; tn += g_threads) {
                    for (int tk = 0; tk < tiles_k; ++tk) {
                        IN_T* dst = b_pk + (size_t)tn * tiles_k * b_tile_sz
                                         + (size_t)tk * b_tile_sz;
                        lb.pack(tn, tk, dst);
                    }
                }
            } else {
                for (int tk = 0; tk < tiles_k; ++tk) {
                    for (int tn = tid; tn < tiles_n; tn += g_threads) {
                        IN_T* dst = b_pk + (size_t)tn * tiles_k * b_tile_sz
                                         + (size_t)tk * b_tile_sz;
                        lb.pack(tn, tk, dst);
                    }
                }
            }
        }
    }

    static void compute_hypertile(const IN_T* a_pk, const IN_T* b_pk,
                                  OUT_T* d, int M, int N, int K) {
        int tiles_m = M / BM, tiles_n = N / BN, tiles_k = K / BK;
        int ht_n = tiles_n / HN;
        int ht_total = (tiles_m / HM) * ht_n;

        #pragma omp parallel num_threads(g_threads) proc_bind(close)
        {
            GEMM::init();
            int tid = omp_get_thread_num();
            for (int ht = tid; ht < ht_total; ht += g_threads) {
                int tm0 = (ht / ht_n) * HM;
                int tn0 = (ht % ht_n) * HN;
                for (int ltm = 0; ltm < HM; ++ltm) {
                    for (int ltn = 0; ltn < HN; ++ltn) {
                        int t_m = tm0 + ltm, t_n = tn0 + ltn;
                        OUT_T* dt = d + t_m * BM * N + t_n * BN;
                        typename GEMM::FragmentD df{};
                        const IN_T* ap = a_pk + (size_t)t_m * tiles_k * a_tile_sz;
                        const IN_T* bp = b_pk + (size_t)t_n * tiles_k * b_tile_sz;
                        for (int tk = 0; tk < tiles_k; ++tk) {
                            typename GEMM::FragmentA af;
                            typename GEMM::FragmentB bf;
                            af.template load(ap);
                            bf.template load(bp);
                            GEMM::MMA::mma(af, bf, df);
                            ap += a_tile_sz;
                            bp += b_tile_sz;
                        }
                        df.template unload_unpack<MajorAxis::ROW>(dt, N);
                    }
                }
            }
        }
    }

    static void compute_static(const IN_T* a_pk, const IN_T* b_pk,
                               OUT_T* d, int M, int N, int K) {
        int tiles_m = M / BM, tiles_n = N / BN, tiles_k = K / BK;
        int total = tiles_m * tiles_n;

        #pragma omp parallel num_threads(g_threads) proc_bind(close)
        {
            GEMM::init();
            int tid = omp_get_thread_num(), nth = omp_get_num_threads();
            int base = total / nth, extra = total % nth;
            int start = tid * base + (tid < extra ? tid : extra);
            int count = base + (tid < extra ? 1 : 0);
            for (int idx = start; idx < start + count; ++idx) {
                int t_m = idx / tiles_n, t_n = idx % tiles_n;
                OUT_T* dt = d + t_m * BM * N + t_n * BN;
                typename GEMM::FragmentD df{};
                const IN_T* ap = a_pk + (size_t)t_m * tiles_k * a_tile_sz;
                const IN_T* bp = b_pk + (size_t)t_n * tiles_k * b_tile_sz;
                for (int tk = 0; tk < tiles_k; ++tk) {
                    typename GEMM::FragmentA af;
                    typename GEMM::FragmentB bf;
                    af.template load(ap);
                    bf.template load(bp);
                    GEMM::MMA::mma(af, bf, df);
                    ap += a_tile_sz;
                    bp += b_tile_sz;
                }
                df.template unload_unpack<MajorAxis::ROW>(dt, N);
            }
        }
    }
};

// ── verification ────────────────────────────────────────────────────
static bool verify(const float* C, const float* C_ref, int M, int N,
                   const char* label) {
    int mismatches = 0;
    for (int i = 0; i < M * N; ++i) {
        float diff = fabsf(C[i] - C_ref[i]);
        float denom = fabsf(C_ref[i]);
        float rel = (denom > 1e-6f) ? diff / denom : diff;
        if (rel > 0.01f && diff > 0.01f) {
            if (mismatches == 0) {
                printf("  [FAIL] %s: idx=%d got=%.4f ref=%.4f diff=%.6f\n",
                       label, i, C[i], C_ref[i], diff);
            }
            ++mismatches;
        }
    }
    if (mismatches == 0) {
        printf("  [PASS] %s\n", label);
    } else {
        printf("  [FAIL] %s: %d mismatches\n", label, mismatches);
    }
    return mismatches == 0;
}

// ── benchmark + write ALL runs ──────────────────────────────────────
using KernelFn = void(*)(const _Float16*, const _Float16*, float*, int, int, int);

static void bench_and_write(FILE* f, const char* name, KernelFn fn,
                            const _Float16* A, const _Float16* B, float* C,
                            int M, int N, int K) {
    for (int w = 0; w < g_warmup; ++w) {
        fn(A, B, C, M, N, K);
    }

    const double flops = 2.0 * (double)M * N * K;
    double total = 0;
    std::vector<double> times(g_iters);

    for (int i = 0; i < g_iters; ++i) {
        double t0 = now_ns();
        fn(A, B, C, M, N, K);
        double t1 = now_ns();
        times[i] = t1 - t0;
        total += times[i];
    }

    double avg = total / g_iters;
    std::vector<double> sorted(times);
    std::sort(sorted.begin(), sorted.end());
    double median = (g_iters % 2 == 0)
        ? (sorted[g_iters / 2 - 1] + sorted[g_iters / 2]) / 2.0
        : sorted[g_iters / 2];

    for (int i = 0; i < g_iters; ++i) {
        fprintf(f, "%s,%d,%d,%d,%d,%d,%.2f,%.6f\n",
                name, M, N, K, g_threads, i, times[i], flops / times[i]);
    }

    printf("    %-18s avg=%.0f ns  med=%.0f ns  min=%.0f  max=%.0f  %.4f GF/s\n",
           name, avg, median, sorted[0], sorted[g_iters - 1], flops / avg);
}

// ── kernel wrappers ─────────────────────────────────────────────────
static float *g_A32, *g_B32;
static _Float16 *g_apk, *g_bpk;

static void wrap_mkl_sgemm(const _Float16* A, const _Float16* B, float* C,
                            int M, int N, int K) {
    mkl_sgemm_wrap(A, B, C, M, N, K, true, g_A32, g_B32);
}

static void wrap_mkl_f16(const _Float16* A, const _Float16* B, float* C,
                          int M, int N, int K) {
    mkl_f16_wrap(A, B, C, M, N, K, true);
}

static void wrap_ref(const _Float16* A, const _Float16* B, float* C,
                     int M, int N, int K) {
    ref_gemm_omp(A, B, C, M, N, K, true);
}

template <bool Hypertile>
static void wrap_amx(const _Float16*, const _Float16*, float* C,
                     int M, int N, int K) {
    if constexpr (Hypertile) {
        amx_gemm_omp<_Float16, float, true, 32, 32, 32>::compute_hypertile(
            g_apk, g_bpk, C, M, N, K);
    } else {
        amx_gemm_omp<_Float16, float, true, 32, 32, 32>::compute_static(
            g_apk, g_bpk, C, M, N, K);
    }
}

// ── arg parsing ─────────────────────────────────────────────────────
struct MNK { int M, N, K; };

static void usage(const char* prog) {
    printf("Usage: %s M1 N1 K1 [M2 N2 K2 ...] [-w warmup] [-i iters] [-t threads] [-o file.csv]\n", prog);
    printf("  Defaults: -w %d -i %d -t %d -o %s\n", g_warmup, g_iters, g_threads, g_outfile);
}

static std::vector<MNK> parse_args(int argc, char** argv) {
    std::vector<MNK> sizes;
    std::vector<int> nums;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-w") == 0 && i + 1 < argc) {
            g_warmup = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            g_iters = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            g_threads = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            g_outfile = argv[++i];
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            exit(0);
        } else {
            nums.push_back(atoi(argv[i]));
        }
    }

    if (nums.size() % 3 != 0 || nums.empty()) {
        fprintf(stderr, "Error: need M N K triples as positional args\n");
        usage(argv[0]);
        exit(1);
    }
    for (size_t i = 0; i < nums.size(); i += 3) {
        sizes.push_back({nums[i], nums[i + 1], nums[i + 2]});
    }
    return sizes;
}

// ── main ────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }

    std::vector<MNK> sizes = parse_args(argc, argv);

    FILE* fout = fopen(g_outfile, "w");
    if (!fout) {
        perror("fopen");
        return 1;
    }

    fprintf(fout, "kernel,M,N,K,threads,run,time_ns,gflops\n");

    printf("Config: warmup=%d iters=%d threads=%d output=%s\n\n",
           g_warmup, g_iters, g_threads, g_outfile);

    for (size_t c = 0; c < sizes.size(); ++c) {
        int M = sizes[c].M, N = sizes[c].N, K = sizes[c].K;
        printf("=== M=%d N=%d K=%d ===\n", M, N, K);

        size_t szA = (size_t)M * K, szB = (size_t)K * N, szC = (size_t)M * N;
        _Float16* A = mmap_alloc<_Float16>(szA);
        _Float16* B = mmap_alloc<_Float16>(szB);
        float* C_ref      = mmap_alloc<float>(szC);
        float* C_amx_ht   = mmap_alloc<float>(szC);
        float* C_amx_st   = mmap_alloc<float>(szC);
        float* C_mkl_s    = mmap_alloc<float>(szC);
        float* C_mkl_f16  = mmap_alloc<float>(szC);
        g_A32 = mmap_alloc<float>(szA);
        g_B32 = mmap_alloc<float>(szB);

        // parallel first-touch init
        #pragma omp parallel num_threads(g_threads) proc_bind(close)
        {
            int tid = omp_get_thread_num(), nth = omp_get_num_threads();
            size_t lo, hi;
            lo = tid * szA / nth; hi = (tid + 1) * szA / nth;
            rand_fill_f16(A + lo, hi - lo);
            lo = tid * szB / nth; hi = (tid + 1) * szB / nth;
            rand_fill_f16(B + lo, hi - lo);
            lo = tid * szC / nth; hi = (tid + 1) * szC / nth;
            memset(C_ref + lo, 0, (hi - lo) * sizeof(float));
            memset(C_amx_ht + lo, 0, (hi - lo) * sizeof(float));
            memset(C_amx_st + lo, 0, (hi - lo) * sizeof(float));
            memset(C_mkl_s + lo, 0, (hi - lo) * sizeof(float));
            memset(C_mkl_f16 + lo, 0, (hi - lo) * sizeof(float));
        }

        // pack for AMX
        g_apk = mmap_alloc<_Float16>(szA);
        g_bpk = mmap_alloc<_Float16>(szB);
        amx_gemm_omp<_Float16, float, true, 32, 32, 32>::pack_parallel(
            A, B, g_apk, g_bpk, M, N, K);

        // correctness
        printf("  Correctness:\n");
        ref_gemm_omp(A, B, C_ref, M, N, K, true);
        wrap_amx<true>(A, B, C_amx_ht, M, N, K);
        verify(C_amx_ht, C_ref, M, N, "amx_hypertile");
        wrap_amx<false>(A, B, C_amx_st, M, N, K);
        verify(C_amx_st, C_ref, M, N, "amx_static");
        mkl_sgemm_wrap(A, B, C_mkl_s, M, N, K, true, g_A32, g_B32);
        verify(C_mkl_s, C_ref, M, N, "mkl_sgemm");
        mkl_f16_wrap(A, B, C_mkl_f16, M, N, K, true);
        verify(C_mkl_f16, C_ref, M, N, "mkl_f16");

        // benchmark
        printf("  Benchmark (%d warmup, %d iters):\n", g_warmup, g_iters);
        bench_and_write(fout, "ref_omp",       wrap_ref,          A, B, C_ref,      M, N, K);
        bench_and_write(fout, "mkl_sgemm",     wrap_mkl_sgemm,    A, B, C_mkl_s,    M, N, K);
        bench_and_write(fout, "mkl_f16",       wrap_mkl_f16,      A, B, C_mkl_f16,  M, N, K);
        bench_and_write(fout, "amx_hypertile", wrap_amx<true>,    A, B, C_amx_ht,   M, N, K);
        bench_and_write(fout, "amx_static",    wrap_amx<false>,   A, B, C_amx_st,   M, N, K);

        // cleanup
        mmap_free(A, szA); mmap_free(B, szB);
        mmap_free(C_ref, szC); mmap_free(C_amx_ht, szC);
        mmap_free(C_amx_st, szC); mmap_free(C_mkl_s, szC);
        mmap_free(C_mkl_f16, szC);
        mmap_free(g_A32, szA); mmap_free(g_B32, szB);
        mmap_free(g_apk, szA); mmap_free(g_bpk, szB);
        printf("\n");
    }

    fclose(fout);
    printf("All %d runs per kernel written to %s\n", g_iters, g_outfile);
    return 0;
}