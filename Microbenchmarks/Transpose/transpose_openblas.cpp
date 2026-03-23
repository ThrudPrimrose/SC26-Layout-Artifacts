#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <cblas.h>

// ══════════════════════════════════════════════════════════════════════
//  V0: cblas_somatcopy, row-major
//      Single call: transpose entire N x N matrix.
// ══════════════════════════════════════════════════════════════════════
static void tr_blas_rm(const float* __restrict__ in, float* __restrict__ out, int N) {
    cblas_somatcopy(CblasRowMajor, CblasTrans,
                    N, N, 1.0f,
                    in,  N,   // src,  lda
                    out, N);  // dst,  ldb
}

// ══════════════════════════════════════════════════════════════════════
//  V1: cblas_somatcopy per SB x SB block, blocked layout
//      Source block (br,bc) is contiguous SB*SB. Dest block (bc,br) is
//      contiguous SB*SB. Each sub-block is row-major, so somatcopy can
//      transpose it with lda=SB, ldb=SB.
// ══════════════════════════════════════════════════════════════════════
static void tr_blas_blk(const float* __restrict__ in, float* __restrict__ out,
                         int N, int SB) {
    int NB = N / SB;
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            cblas_somatcopy(CblasRowMajor, CblasTrans,
                            SB, SB, 1.0f,
                            src, SB,
                            dst, SB);
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
//  V2: cblas_somatcopy per block, blocked layout, OMP parallel
// ══════════════════════════════════════════════════════════════════════
static void tr_blas_blk_omp(const float* __restrict__ in, float* __restrict__ out,
                              int N, int SB) {
    int NB = N / SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            cblas_somatcopy(CblasRowMajor, CblasTrans,
                            SB, SB, 1.0f,
                            src, SB,
                            dst, SB);
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
//  V3: cblas_somatcopy per tile, row-major layout, OMP parallel
//      Tile the N x N matrix into SB x SB tiles.  Each tile is a
//      sub-matrix of the full row-major array, so lda = ldb = N.
//      src tile (tr,tc) at row tr*SB, col tc*SB  -->
//      dst tile at row tc*SB, col tr*SB (transposed position).
//
//      somatcopy(Trans, rows=SB, cols=SB, alpha=1,
//                &in[tr*SB * N + tc*SB], lda=N,
//                &out[tc*SB * N + tr*SB], ldb=N)
//
//      This transposes the SB x SB sub-matrix and writes it into the
//      transposed position in the output, with full-matrix strides.
// ══════════════════════════════════════════════════════════════════════
static void tr_blas_rm_omp(const float* __restrict__ in, float* __restrict__ out,
                            int N, int SB) {
    int NT = (N + SB - 1) / SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * SB, c0 = tc * SB;
            int tile_h = std::min(SB, N - r0);
            int tile_w = std::min(SB, N - c0);
            cblas_somatcopy(CblasRowMajor, CblasTrans,
                            tile_h, tile_w, 1.0f,
                            in  + r0 * N + c0, N,
                            out + c0 * N + r0, N);
        }
    }
}

// ══════════════════════════════════════════════════════════════════════

static const char* V_NAMES[] = {
    "openblas",           // 0
    "openblas_blk",       // 1
    "openblas_blk_omp",   // 2
    "openblas_rm_omp",    // 3
};
static const int N_VARIANTS = 4;

static bool is_blocked(int var) { return var == 1 || var == 2; }

// ── Verification ──
static void ref_transpose(const float* in, float* ref, int N) {
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            ref[c * N + r] = in[r * N + c];
}

static float verify(const float* out, const float* ref, int N, int SB, bool blocked_out) {
    float maxerr = 0.0f;
    if (blocked_out) {
        int NB = N / SB;
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++) {
                float got = out[(r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB];
                float exp = ref[r * N + c];
                float err = fabsf(got - exp);
                if (err > maxerr) maxerr = err;
            }
    } else {
        for (size_t i = 0; i < (size_t)N*N; i++) {
            float err = fabsf(out[i] - ref[i]);
            if (err > maxerr) maxerr = err;
        }
    }
    return maxerr;
}

// ══════════════════════════════════════════════════════════════════════

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr,
            "Usage: %s <N> <variant> <csv> [SB=32] [WARMUP=3] [REPS=20] [THREADS=0]\n"
            "\n"
            "  variant:\n"
            "    0 = openblas           (somatcopy, row-major, single call)\n"
            "    1 = openblas_blk       (somatcopy per SB x SB block, serial)\n"
            "    2 = openblas_blk_omp   (somatcopy per SB x SB block, OMP)\n"
            "    3 = openblas_rm_omp    (somatcopy per SB x SB tile, row-major, OMP)\n",
            argv[0]);
        return 1;
    }

    int N   = atoi(argv[1]);
    int VAR = atoi(argv[2]);
    const char* csv = argv[3];
    int SB      = (argc > 4) ? atoi(argv[4]) : 32;
    int WARMUP  = (argc > 5) ? atoi(argv[5]) : 3;
    int REPS    = (argc > 6) ? atoi(argv[6]) : 20;
    int THREADS = (argc > 7) ? atoi(argv[7]) : 0;

    if (VAR < 0 || VAR >= N_VARIANTS) {
        fprintf(stderr, "Unknown variant %d\n", VAR);
        return 1;
    }
    if (is_blocked(VAR) && N % SB != 0) {
        fprintf(stderr, "N=%d not divisible by SB=%d for blocked variant\n", N, SB);
        return 1;
    }
    if (THREADS > 0) {
        omp_set_num_threads(THREADS);
    }
    // For OMP-parallel variants, pin OpenBLAS to 1 thread per call
    // to avoid oversubscription inside the parallel region.
    // For V0/V1 (serial), let OpenBLAS use all threads (though
    // somatcopy may still be single-threaded internally).
    if (VAR == 2 || VAR == 3) {
        openblas_set_num_threads(1);
    } else if (THREADS > 0) {
        openblas_set_num_threads(THREADS);
    }
    int nthreads;
    #pragma omp parallel
    {
        #pragma omp single
        nthreads = omp_get_num_threads();
    }

    size_t elems = (size_t)N * N;
    size_t bytes = elems * sizeof(float);

    float* h_row = (float*)aligned_alloc(64, bytes);
    float* h_ref = (float*)aligned_alloc(64, bytes);
    float* h_in  = (float*)aligned_alloc(64, bytes);
    float* h_out = (float*)aligned_alloc(64, bytes);

    for (size_t i = 0; i < elems; i++) h_row[i] = (float)i / (float)N;
    ref_transpose(h_row, h_ref, N);

    if (is_blocked(VAR)) {
        int NB = N / SB;
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                h_in[(r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB] = h_row[r * N + c];
    } else {
        memcpy(h_in, h_row, bytes);
    }

    auto launch = [&]() {
        switch (VAR) {
            case 0: tr_blas_rm     (h_in, h_out, N);      break;
            case 1: tr_blas_blk    (h_in, h_out, N, SB);  break;
            case 2: tr_blas_blk_omp(h_in, h_out, N, SB);  break;
            case 3: tr_blas_rm_omp (h_in, h_out, N, SB);  break;
        }
    };

    // Warmup
    for (int i = 0; i < WARMUP; i++) launch();

    // Verify
    memset(h_out, 0, bytes);
    launch();
    float maxerr = verify(h_out, h_ref, N, SB, is_blocked(VAR));
    bool pass = (maxerr == 0.0f);

    // Timed runs
    double* times = (double*)malloc(REPS * sizeof(double));
    for (int i = 0; i < REPS; i++) {
        double t0 = omp_get_wtime();
        launch();
        double t1 = omp_get_wtime();
        times[i] = t1 - t0;
    }

    double total = 0;
    for (int i = 0; i < REPS; i++) total += times[i];
    double mean_s = total / REPS;

    std::sort(times, times + REPS);
    double med_s = times[REPS / 2];
    double p5_s  = times[(int)(REPS * 0.05)];
    double p95_s = times[(int)(REPS * 0.95)];

    double bpi   = 2.0 * N * (double)N * sizeof(float);
    double mean_gbps = bpi / mean_s / 1e9;
    double med_gbps  = bpi / med_s  / 1e9;

    double cksum = 0;
    for (size_t i = 0; i < elems; i++) cksum += h_out[i];

    // Use TB=0, MT=0 placeholders to match cpu_transpose CSV format
    int TB = 0, MT = 0;

    printf("%s N=%d SB=%d threads=%d | "
           "mean %.4f ms (%.1f GB/s)  med %.4f ms (%.1f GB/s)  "
           "p5 %.4f ms  p95 %.4f ms  maxerr=%.1e  %s  cksum=%.6e\n",
           V_NAMES[VAR], N, SB, nthreads,
           mean_s*1e3, mean_gbps, med_s*1e3, med_gbps,
           p5_s*1e3, p95_s*1e3, maxerr, pass ? "PASS" : "FAIL", cksum);

    // CSV: same format as cpu_transpose.cpp
    // variant, N, TB, SB, MT, nthreads, rep, time_s, gbs, cksum, PASS/FAIL
    FILE* f = fopen(csv, "a");
    if (f) {
        for (int i = 0; i < REPS; i++) {
            double gbs = bpi / times[i] / 1e9;
            fprintf(f, "%s,%d,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e,%s\n",
                    V_NAMES[VAR], N, TB, SB, MT, nthreads, i,
                    times[i], gbs, cksum, pass ? "PASS" : "FAIL");
        }
        fclose(f);
    }

    free(times); free(h_row); free(h_ref); free(h_in); free(h_out);
    return pass ? 0 : 1;
}