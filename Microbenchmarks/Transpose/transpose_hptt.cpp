/*  transpose_hptt.cpp — HPTT transpose benchmark with NUMA-aware allocation.
 *
 *  All variants use HPTT's internal threading (numThreads = OMP thread count).
 *  Plans are pre-created outside the timed region.
 *
 *  Variants:
 *    0 = hptt             full N×N 2D transpose, ESTIMATE
 *    1 = hptt_blk         blocked layout as 4D tensor [NB,NB,SB,SB], perm {1,0,3,2}
 *    2 = hptt_rm_tiled    row-major tiled, serial loop over per-tile 2D plans
 *    3 = hptt_patient     full N×N 2D transpose, PATIENT auto-tune
 *
 *  Build:
 *    g++ -O3 -march=native -fopenmp -I<hptt>/include -L<hptt>/lib \
 *        -o transpose_hptt transpose_hptt.cpp -lhptt
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <memory>
#include <vector>
#include <omp.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <hptt.h>

using Plan = std::shared_ptr<hptt::Transpose<float>>;

/* ═══════════════════════════════════════════════════════════════════════
 *  NUMA helpers (no libnuma — direct syscalls)
 * ═══════════════════════════════════════════════════════════════════════ */
#ifndef MPOL_BIND
#define MPOL_BIND 2
#endif

static long sys_mbind(void *a, unsigned long l, int m,
                      const unsigned long *nm, unsigned long mx, unsigned f) {
    return syscall(SYS_mbind, a, l, m, nm, mx, f);
}

static int detect_numa_nodes() {
    int n = 0;
    for (int i = 0; i < 128; i++) {
        char p[128];
        snprintf(p, 128, "/sys/devices/system/node/node%d", i);
        if (access(p, F_OK) == 0)
            n = i + 1;
        else if (n)
            break;
    }
    return n > 0 ? n : 1;
}

static void bind_pages(void *a, size_t l, int node) {
    if (!l || node < 0)
        return;
    unsigned long m[4] = {};
    m[node / 64] |= 1UL << (node % 64);
    sys_mbind(a, l, MPOL_BIND, m, 257, 0);
}

static size_t g_pagesz = 0;
static size_t pagesz() {
    if (!g_pagesz)
        g_pagesz = (size_t)sysconf(_SC_PAGESIZE);
    return g_pagesz;
}

template <typename T>
static T *numa_alloc(size_t count) {
    size_t bytes = count * sizeof(T);
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) {
        perror("mmap");
        std::abort();
    }
    madvise(p, bytes, MADV_NOHUGEPAGE);
    return (T *)p;
}

template <typename T>
static void numa_dealloc(T *p, size_t count) {
    if (p)
        munmap(p, count * sizeof(T));
}

static void bind_by_rows(void *base, int N, int D) {
    size_t ps = pagesz();
    size_t row_bytes = (size_t)N * sizeof(float);
    size_t total = (size_t)N * row_bytes;
    int rows_per = (N + D - 1) / D;

    for (int d = 0; d < D; d++) {
        size_t lo = (size_t)d * rows_per * row_bytes;
        size_t hi = std::min((size_t)(d + 1) * rows_per * row_bytes, total);
        lo = (lo / ps) * ps;
        hi = ((hi + ps - 1) / ps) * ps;
        if (hi > ((total + ps - 1) / ps) * ps)
            hi = ((total + ps - 1) / ps) * ps;
        if (hi > lo)
            bind_pages((char *)base + lo, hi - lo, d);
    }
}

static void bind_by_cols(void *base, int N, int D) {
    size_t ps = pagesz();
    size_t elem_sz = sizeof(float);
    size_t cols_per_domain = (size_t)N / D;
    size_t row_bytes = (size_t)N * elem_sz;
    size_t total_bytes = (size_t)N * row_bytes;
    size_t total_pages = (total_bytes + ps - 1) / ps;

    for (size_t pg = 0; pg < total_pages; pg++) {
        size_t byte_off = pg * ps;
        size_t elem0 = byte_off / elem_sz;
        size_t col0 = elem0 % (size_t)N;
        int domain = (int)(col0 / cols_per_domain);
        if (domain >= D)
            domain = D - 1;
        size_t len = ps;
        if (byte_off + len > total_bytes)
            len = total_bytes - byte_off;
        bind_pages((char *)base + byte_off, len, domain);
    }
}

static void bind_contiguous(void *base, size_t total_bytes, int D) {
    size_t ps = pagesz();
    size_t per = total_bytes / D;
    for (int d = 0; d < D; d++) {
        size_t lo = (size_t)d * per;
        size_t hi = (d == D - 1) ? total_bytes : (size_t)(d + 1) * per;
        lo = (lo / ps) * ps;
        hi = ((hi + ps - 1) / ps) * ps;
        if (hi > ((total_bytes + ps - 1) / ps) * ps)
            hi = ((total_bytes + ps - 1) / ps) * ps;
        if (hi > lo)
            bind_pages((char *)base + lo, hi - lo, d);
    }
}

static void parallel_init(float *p, size_t n, float val) {
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++)
        p[i] = val;
}

/* ═══════════════════════════════════════════════════════════════════════ */

static const char *V_NAMES[] = {
    "hptt", "hptt_blk", "hptt_rm_tiled", "hptt_patient",
};
static const int N_VARIANTS = 4;
static bool is_blocked(int var) { return var == 1; }

/* ═══════════════════════════════════════════════════════════════════════
 *  Plan factories — all use HPTT-internal threading
 * ═══════════════════════════════════════════════════════════════════════ */

/* 2D full transpose: perm {1, 0} on [N, N] */
static Plan make_full_plan(const float *A, float *B, int N,
                           int numThreads, hptt::SelectionMethod method) {
    int perm[2] = {1, 0};
    int size[2] = {N, N};
    float alpha = 1.0f, beta = 0.0f;
    return hptt::create_plan(perm, 2, alpha, A, size, nullptr,
                             beta, B, nullptr,
                             method, numThreads, nullptr, true);
}

/* 4D blocked transpose: perm {1, 0, 3, 2} on [NB, NB, SB, SB]
 *
 * Input  layout:  A[br][bc][lr][lc]  — block (br,bc), local element (lr,lc)
 * Output layout:  B[bc][br][lc][lr]  — swap block indices AND local indices
 *
 * Single HPTT call transposes the entire blocked matrix. */
static Plan make_blocked_plan(const float *src, float *dst,
                              int N, int SB, int numThreads,
                              hptt::SelectionMethod method) {
    int NB = N / SB;
    int perm[4] = {1, 0, 3, 2};
    int size[4] = {NB, NB, SB, SB};
    float alpha = 1.0f, beta = 0.0f;
    return hptt::create_plan(perm, 4, alpha, src, size, nullptr,
                             beta, dst, nullptr,
                             method, numThreads, nullptr, true);
}

/* 2D tile plan: transpose a th×tw sub-matrix at (r0,c0) within N×N row-major */
static Plan make_tile_plan(const float *in, float *out, int N,
                           int r0, int c0, int th, int tw, int numThreads) {
    int perm[2] = {1, 0};
    int size[2] = {th, tw};
    int outerA[2] = {N, N};
    int outerB[2] = {N, N};
    float alpha = 1.0f, beta = 0.0f;
    return hptt::create_plan(perm, 2, alpha,
                             in + r0 * N + c0, size, outerA,
                             beta,
                             out + c0 * N + r0, outerB,
                             hptt::ESTIMATE, numThreads, nullptr, true);
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Verification
 * ═══════════════════════════════════════════════════════════════════════ */

static void ref_transpose(const float *in, float *ref, int N) {
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            ref[c * N + r] = in[r * N + c];
}

static float verify(const float *out, const float *ref, int N, int SB, bool blk) {
    float mx = 0.0f;
    if (blk) {
        int NB = N / SB;
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++) {
                float e = fabsf(out[(r / SB * NB + c / SB) * SB * SB +
                                    (r % SB) * SB + c % SB] -
                                ref[r * N + c]);
                if (e > mx)
                    mx = e;
            }
    } else {
        for (size_t i = 0; i < (size_t)N * N; i++) {
            float e = fabsf(out[i] - ref[i]);
            if (e > mx)
                mx = e;
        }
    }
    return mx;
}

/* ═══════════════════════════════════════════════════════════════════════
 *  main
 * ═══════════════════════════════════════════════════════════════════════ */

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <N> <variant> <csv> [SB=32] [WARMUP=3] [REPS=20] [THREADS=0]\n"
                "\n"
                "  0 = hptt           full N×N 2D, ESTIMATE\n"
                "  1 = hptt_blk       blocked 4D [NB,NB,SB,SB], perm {1,0,3,2}\n"
                "  2 = hptt_rm_tiled  row-major tiled, per-tile 2D plans\n"
                "  3 = hptt_patient   full N×N 2D, PATIENT auto-tune\n",
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
        fprintf(stderr, "N=%d not div by SB=%d\n", N, SB);
        return 1;
    }
    if (THREADS > 0)
        omp_set_num_threads(THREADS);

    int nthreads;
#pragma omp parallel
    {
#pragma omp single
        nthreads = omp_get_num_threads();
    }

    int D = detect_numa_nodes();
    int hptt_threads = nthreads;

    size_t elems = (size_t)N * N;
    size_t bytes = elems * sizeof(float);

    printf("  %s N=%d SB=%d thr=%d hptt_thr=%d numa_nodes=%d\n",
           V_NAMES[VAR], N, SB, nthreads, hptt_threads, D);

    /* ── Allocate with NUMA binding ───────────────────────────────── */
    float *h_in, *h_out, *h_row, *h_ref;
    bool use_numa = (D >= 2);

    if (use_numa) {
        h_in  = numa_alloc<float>(elems);
        h_out = numa_alloc<float>(elems);
        h_row = numa_alloc<float>(elems);
        h_ref = (float *)aligned_alloc(64, bytes);

        if (is_blocked(VAR)) {
            bind_contiguous(h_in, bytes, D);
            bind_contiguous(h_out, bytes, D);
        } else {
            bind_by_rows(h_in, N, D);
            bind_by_cols(h_out, N, D);
        }
        bind_by_rows(h_row, N, D);

        printf("  NUMA: %s\n", is_blocked(VAR)
                   ? "blocked contiguous split"
                   : "input rows->domains, output cols->domains");
    } else {
        h_in  = (float *)aligned_alloc(64, bytes);
        h_out = (float *)aligned_alloc(64, bytes);
        h_row = (float *)aligned_alloc(64, bytes);
        h_ref = (float *)aligned_alloc(64, bytes);
    }

    /* ── Init ─────────────────────────────────────────────────────── */
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < elems; i++)
        h_row[i] = (float)i / (float)N;

    ref_transpose(h_row, h_ref, N);

    if (is_blocked(VAR)) {
        int NB = N / SB;
#pragma omp parallel for schedule(static)
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                h_in[(r / SB * NB + c / SB) * SB * SB +
                      (r % SB) * SB + c % SB] = h_row[r * N + c];
    } else {
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < elems; i++)
            h_in[i] = h_row[i];
    }

    parallel_init(h_out, elems, 0.0f);

    /* ── Pre-create plans (all use hptt_threads) ──────────────────── */
    Plan main_plan;
    std::vector<Plan> tile_plans;
    int NT = 0;

    printf("  Creating plans...\n");
    double pt0 = omp_get_wtime();

    switch (VAR) {
    case 0:
        main_plan = make_full_plan(h_in, h_out, N, hptt_threads, hptt::ESTIMATE);
        break;
    case 1:
        main_plan = make_blocked_plan(h_in, h_out, N, SB, hptt_threads, hptt::ESTIMATE);
        break;
    case 2: {
        NT = (N + SB - 1) / SB;
        tile_plans.resize(NT * NT);
        for (int tr = 0; tr < NT; tr++)
            for (int tc = 0; tc < NT; tc++) {
                int r0 = tr * SB, c0 = tc * SB;
                int th = std::min(SB, N - r0), tw = std::min(SB, N - c0);
                tile_plans[tr * NT + tc] =
                    make_tile_plan(h_in, h_out, N, r0, c0, th, tw, hptt_threads);
            }
        break;
    }
    case 3:
        printf("  PATIENT auto-tuning...\n");
        main_plan = make_full_plan(h_in, h_out, N, hptt_threads, hptt::PATIENT);
        break;
    }

    printf("  Plans: %.3f s (%d plans)\n",
           omp_get_wtime() - pt0,
           VAR == 2 ? (int)tile_plans.size() : 1);

    /* ── Launch (HPTT parallelizes internally) ────────────────────── */
    auto launch = [&]() {
        switch (VAR) {
        case 0:
        case 1:
        case 3:
            main_plan->execute();
            break;
        case 2:
            for (int i = 0; i < NT * NT; i++)
                tile_plans[i]->execute();
            break;
        }
    };

    for (int i = 0; i < WARMUP; i++)
        launch();

    /* Verify */
    parallel_init(h_out, elems, 0.0f);
    launch();
    float maxerr = verify(h_out, h_ref, N, SB, is_blocked(VAR));
    bool pass = (maxerr <= 1e-5f);
    if (!pass)
        fprintf(stderr, "  [%s] VERIFY FAIL maxerr=%.6e\n", V_NAMES[VAR], maxerr);

    /* Timed */
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

    double cksum = 0;
    for (size_t i = 0; i < elems; i++)
        cksum += h_out[i];

    int TB = 0, MT = 0;
    printf("%s N=%d SB=%d thr=%d hptt_thr=%d | "
           "med %.4f ms (%.1f GB/s)  p5 %.4f  p95 %.4f  "
           "maxerr=%.1e  %s  cksum=%.6e\n",
           V_NAMES[VAR], N, SB, nthreads, hptt_threads,
           med_s * 1e3, bpi / med_s / 1e9,
           times[(int)(REPS * 0.05)] * 1e3,
           times[(int)(REPS * 0.95)] * 1e3,
           maxerr, pass ? "PASS" : "FAIL", cksum);

    FILE *f = fopen(csv, "a");
    if (f) {
        for (int i = 0; i < REPS; i++)
            fprintf(f, "%s,%d,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e,%s\n",
                    V_NAMES[VAR], N, TB, SB, MT, nthreads, i,
                    times[i], bpi / times[i] / 1e9, cksum,
                    pass ? "PASS" : "FAIL");
        fclose(f);
    }

    free(times);
    if (use_numa) {
        numa_dealloc(h_in, elems);
        numa_dealloc(h_out, elems);
        numa_dealloc(h_row, elems);
        free(h_ref);
    } else {
        free(h_in);
        free(h_out);
        free(h_row);
        free(h_ref);
    }
    return pass ? 0 : 1;
}