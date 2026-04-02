/*  transpose_hptt.cpp — HPTT-based matrix transpose benchmark.
 *
 *  Variants:
 *    0 = hptt             full N×N, HPTT-managed threading, ESTIMATE
 *    1 = hptt_blk         blocked layout, serial per-block
 *    2 = hptt_blk_omp     blocked layout, OMP parallel per-block
 *    3 = hptt_rm_omp      row-major tiled, OMP parallel per-tile
 *    4 = hptt_patient     full N×N, HPTT-managed threading, PATIENT auto-tune
 *
 *  Build:
 *    g++ -O3 -march=native -fopenmp -I<hptt>/include -L<hptt>/lib \
 *        -o transpose_hptt transpose_hptt.cpp -lhptt
 *
 *  CSV format matches transpose_cpu.cpp:
 *    variant, N, TB, SB, MT, nthreads, rep, time_s, gbs, cksum, status
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <memory>
#include <omp.h>
#include <hptt.h>

/* ═══════════════════════════════════════════════════════════════════════
 *  Helpers
 * ═══════════════════════════════════════════════════════════════════════ */

/* Transpose full N×N row-major matrix using a pre-created plan. */
static std::shared_ptr<hptt::Transpose<float>>
make_full_plan(const float *A, float *B, int N, int numThreads,
               hptt::SelectionMethod method) {
    int perm[2] = {1, 0};
    int size[2] = {N, N};
    float alpha = 1.0f, beta = 0.0f;
    return hptt::create_plan(perm, 2, alpha, A, size, nullptr,
                             beta, B, nullptr,
                             method, numThreads,
                             nullptr, true /* useRowMajor */);
}

/* Transpose a single SB×SB contiguous sub-block (stride = SB). */
static void transpose_block(const float *src, float *dst, int SB) {
    int perm[2] = {1, 0};
    int size[2] = {SB, SB};
    float alpha = 1.0f, beta = 0.0f;
    auto plan = hptt::create_plan(perm, 2, alpha, src, size, nullptr,
                                  beta, dst, nullptr,
                                  hptt::ESTIMATE, 1,
                                  nullptr, true);
    plan->execute();
}

/* Transpose an SB×SB tile of a row-major N×N matrix.
 * src at (r0,c0), dst at (c0,r0) in the output matrix. */
static void transpose_tile(const float *in, float *out,
                           int N, int r0, int c0, int tile_h, int tile_w) {
    int perm[2] = {1, 0};
    int size[2] = {tile_h, tile_w};
    int outerA[2] = {N, N};
    int outerB[2] = {N, N};
    float alpha = 1.0f, beta = 0.0f;
    auto plan = hptt::create_plan(perm, 2, alpha,
                                  in + r0 * N + c0, size, outerA,
                                  beta,
                                  out + c0 * N + r0, nullptr, outerB,
                                  hptt::ESTIMATE, 1,
                                  nullptr, true);
    plan->execute();
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V0: Full matrix, HPTT threading, ESTIMATE
 * ═══════════════════════════════════════════════════════════════════════ */
static void tr_hptt(const float *in, float *out, int N,
                    std::shared_ptr<hptt::Transpose<float>> &plan) {
    plan->execute();
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V1: Blocked layout, serial per-block
 * ═══════════════════════════════════════════════════════════════════════ */
static void tr_hptt_blk(const float *in, float *out, int N, int SB) {
    int NB = N / SB;
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            const float *src = in + (br * NB + bc) * SB * SB;
            float *dst = out + (bc * NB + br) * SB * SB;
            transpose_block(src, dst, SB);
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V2: Blocked layout, OMP parallel per-block
 * ═══════════════════════════════════════════════════════════════════════ */
static void tr_hptt_blk_omp(const float *in, float *out, int N, int SB) {
    int NB = N / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            const float *src = in + (br * NB + bc) * SB * SB;
            float *dst = out + (bc * NB + br) * SB * SB;
            transpose_block(src, dst, SB);
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V3: Row-major tiled, OMP parallel per-tile
 * ═══════════════════════════════════════════════════════════════════════ */
static void tr_hptt_rm_omp(const float *in, float *out, int N, int SB) {
    int NT = (N + SB - 1) / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * SB, c0 = tc * SB;
            int tile_h = std::min(SB, N - r0);
            int tile_w = std::min(SB, N - c0);
            transpose_tile(in, out, N, r0, c0, tile_h, tile_w);
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V4: Full matrix, HPTT threading, PATIENT auto-tune
 * ═══════════════════════════════════════════════════════════════════════ */
/* (uses the same tr_hptt function, just with a PATIENT plan) */

/* ═══════════════════════════════════════════════════════════════════════ */

static const char *V_NAMES[] = {
    "hptt",          /* 0 */
    "hptt_blk",      /* 1 */
    "hptt_blk_omp",  /* 2 */
    "hptt_rm_omp",   /* 3 */
    "hptt_patient",  /* 4 */
};
static const int N_VARIANTS = 5;

static bool is_blocked(int var) { return var == 1 || var == 2; }

/* ── Verification ──────────────────────────────────────────────────── */

static void ref_transpose(const float *in, float *ref, int N) {
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            ref[c * N + r] = in[r * N + c];
}

static float verify(const float *out, const float *ref, int N, int SB, bool blocked_out) {
    float maxerr = 0.0f;
    if (blocked_out) {
        int NB = N / SB;
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++) {
                float got = out[(r / SB * NB + c / SB) * SB * SB + (r % SB) * SB + c % SB];
                float exp = ref[r * N + c];
                float err = fabsf(got - exp);
                if (err > maxerr)
                    maxerr = err;
            }
    } else {
        for (size_t i = 0; i < (size_t)N * N; i++) {
            float err = fabsf(out[i] - ref[i]);
            if (err > maxerr)
                maxerr = err;
        }
    }
    return maxerr;
}

/* ═══════════════════════════════════════════════════════════════════════
 *  main
 * ═══════════════════════════════════════════════════════════════════════ */

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <N> <variant> <csv> [SB=32] [WARMUP=3] [REPS=20] [THREADS=0]\n"
                "\n"
                "  variant:\n"
                "    0 = hptt           (full N×N, ESTIMATE, HPTT threading)\n"
                "    1 = hptt_blk       (blocked, serial per-block)\n"
                "    2 = hptt_blk_omp   (blocked, OMP per-block)\n"
                "    3 = hptt_rm_omp    (row-major tiled, OMP per-tile)\n"
                "    4 = hptt_patient   (full N×N, PATIENT auto-tune)\n",
                argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int VAR = atoi(argv[2]);
    const char *csv = argv[3];
    int SB = (argc > 4) ? atoi(argv[4]) : 32;
    int WARMUP = (argc > 5) ? atoi(argv[5]) : 3;
    int REPS = (argc > 6) ? atoi(argv[6]) : 20;
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
    int nthreads;
#pragma omp parallel
    {
#pragma omp single
        nthreads = omp_get_num_threads();
    }

    /* For V0/V4: HPTT manages its own threading.
     * For V2/V3: OMP parallel, each tile uses 1 HPTT thread. */
    int hptt_threads = (VAR == 0 || VAR == 4) ? nthreads : 1;

    size_t elems = (size_t)N * N;
    size_t bytes = elems * sizeof(float);

    float *h_row = (float *)aligned_alloc(64, bytes);
    float *h_ref = (float *)aligned_alloc(64, bytes);
    float *h_in = (float *)aligned_alloc(64, bytes);
    float *h_out = (float *)aligned_alloc(64, bytes);

    /* Parallel first-touch */
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < elems; i++)
        h_row[i] = (float)i / (float)N;

    ref_transpose(h_row, h_ref, N);

    if (is_blocked(VAR)) {
        int NB = N / SB;
#pragma omp parallel for schedule(static)
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                h_in[(r / SB * NB + c / SB) * SB * SB + (r % SB) * SB + c % SB] =
                    h_row[r * N + c];
    } else {
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < elems; i++)
            h_in[i] = h_row[i];
    }

    /* Parallel first-touch output */
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < elems; i++)
        h_out[i] = 0.0f;

    /* Pre-create plan for V0/V4 (outside timed region) */
    std::shared_ptr<hptt::Transpose<float>> full_plan;
    if (VAR == 0) {
        full_plan = make_full_plan(h_in, h_out, N, hptt_threads, hptt::ESTIMATE);
    } else if (VAR == 4) {
        printf("  [INFO] HPTT PATIENT auto-tuning (may take a while)...\n");
        full_plan = make_full_plan(h_in, h_out, N, hptt_threads, hptt::PATIENT);
        printf("  [INFO] Auto-tuning done.\n");
    }

    auto launch = [&]() {
        switch (VAR) {
        case 0:
        case 4:
            full_plan->execute();
            break;
        case 1:
            tr_hptt_blk(h_in, h_out, N, SB);
            break;
        case 2:
            tr_hptt_blk_omp(h_in, h_out, N, SB);
            break;
        case 3:
            tr_hptt_rm_omp(h_in, h_out, N, SB);
            break;
        }
    };

    /* Warmup */
    for (int i = 0; i < WARMUP; i++)
        launch();

    /* Verify */
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < elems; i++)
        h_out[i] = 0.0f;
    launch();
    float maxerr = verify(h_out, h_ref, N, SB, is_blocked(VAR));
    bool pass = (maxerr == 0.0f);

    /* Timed runs */
    double *times = (double *)malloc(REPS * sizeof(double));
    for (int i = 0; i < REPS; i++) {
        double t0 = omp_get_wtime();
        launch();
        double t1 = omp_get_wtime();
        times[i] = t1 - t0;
    }

    std::sort(times, times + REPS);
    double bpi = 2.0 * N * (double)N * sizeof(float);
    double med_s = times[REPS / 2];
    double med_gbps = bpi / med_s / 1e9;

    double cksum = 0;
    for (size_t i = 0; i < elems; i++)
        cksum += h_out[i];

    int TB = 0, MT = 0;

    printf("%s N=%d SB=%d threads=%d hptt_threads=%d | "
           "med %.4f ms (%.1f GB/s)  p5 %.4f ms  p95 %.4f ms  "
           "maxerr=%.1e  %s  cksum=%.6e\n",
           V_NAMES[VAR], N, SB, nthreads, hptt_threads,
           med_s * 1e3, med_gbps,
           times[(int)(REPS * 0.05)] * 1e3,
           times[(int)(REPS * 0.95)] * 1e3,
           maxerr, pass ? "PASS" : "FAIL", cksum);

    /* CSV: same 11-column format as transpose_cpu.cpp */
    FILE *f = fopen(csv, "a");
    if (f) {
        for (int i = 0; i < REPS; i++) {
            double gbs = bpi / times[i] / 1e9;
            fprintf(f, "%s,%d,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e,%s\n",
                    V_NAMES[VAR], N, TB, SB, MT, nthreads, i,
                    times[i], gbs, cksum, pass ? "PASS" : "FAIL");
        }
        fclose(f);
    }

    free(times);
    free(h_row);
    free(h_ref);
    free(h_in);
    free(h_out);
    return pass ? 0 : 1;
}