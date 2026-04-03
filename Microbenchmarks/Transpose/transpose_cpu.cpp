/*  transpose_cpu.cpp — Fully-templated CPU transpose with NUMA.
 *  TB, SB, MT as template params; dispatch macro instantiates all swept combos.
 *  Compile:  g++ -O3 -march=native -fopenmp -o transpose_cpu transpose_cpu.cpp
 */
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <vector>

/* ── NUMA (no libnuma) ────────────────────────────────────────────── */
#ifndef MPOL_BIND
#define MPOL_BIND 2
#endif


static long sys_mbind(void *a, unsigned long l, int m, const unsigned long *nm, unsigned long mx,
                      unsigned f) {
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
    return n ? n : 1;
}
static void bind_pages(void *a, size_t l, int node) {
    if (!l || node < 0)
        return;
    unsigned long m[4] = {};
    m[node / 64] |= 1UL << (node % 64);
    sys_mbind(a, l, MPOL_BIND, m, 257, 0);
}
template <typename T> static T *numa_alloc(size_t c) {
    size_t b = c * sizeof(T);
    void *p = mmap(nullptr, b, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE,
                   -1, 0);
    if (p == MAP_FAILED) {
        perror("mmap");
        std::abort();
    }
#if NOHUGEPAGE
    madvise(p, b, MADV_NOHUGEPAGE);
#else
    madvise(p, b, MADV_HUGEPAGE);
#endif
    return (T *)p;
}
template <typename T> static void numa_dealloc(T *p, size_t c) {
    if (p)
        munmap(p, c * sizeof(T));
}

enum NumaPolicy { NP_NONE = 0, NP_CONTIG = 1, NP_CYCLIC = 2, NP_TRANSPOSE = 3 };
static const char *NP_SUFFIX[] = {"", "_nd", "_nc", "_nt"};
static const int NP_COUNT = 4;
static size_t pagesz() {
    static size_t p = sysconf(_SC_PAGESIZE);
    return p;
}

static void bind_contiguous(void *base, size_t tot, int D) {
    size_t ps = pagesz();
    for (int d = 0; d < D; d++) {
        size_t lo = ((size_t)d * tot / D / ps) * ps;
        size_t hi = (((size_t)(d + 1) * tot / D + ps - 1) / ps) * ps;
        if (hi > ((tot + ps - 1) / ps) * ps)
            hi = ((tot + ps - 1) / ps) * ps;
        if (hi > lo)
            bind_pages((char *)base + lo, hi - lo, d);
    }
}
static void bind_by_map(void *base, int NB2, size_t bb, const std::vector<int> &nm) {
    size_t ps = pagesz();
    int cn = nm[0];
    size_t rs = 0;
    for (int b = 1; b <= NB2; b++) {
        int n = (b < NB2) ? nm[b] : -1;
        if (n != cn) {
            size_t re = (size_t)b * bb;
            bind_pages((char *)base + (rs / ps) * ps, ((re + ps - 1) / ps) * ps - (rs / ps) * ps,
                       cn);
            cn = n;
            rs = re;
        }
    }
}
/* Bind columns [d*N/D, (d+1)*N/D) to NUMA node d.
 * In row-major, each page within a row is bound to the domain owning its column range. */
static void bind_by_cols(void *base, int N, int D) {
    size_t ps = pagesz();
    size_t esz = sizeof(float);
    size_t cols_per_d = (size_t)N / D;
    size_t total_bytes = (size_t)N * N * esz;
    size_t total_pages = (total_bytes + ps - 1) / ps;
    size_t elems_per_page = ps / esz;

    for (size_t pg = 0; pg < total_pages; pg++) {
        size_t byte_off = pg * ps;
        size_t elem0 = byte_off / esz;
        size_t col0 = elem0 % (size_t)N;
        int domain = std::min((int)(col0 / cols_per_d), D - 1);
        size_t len = std::min(ps, total_bytes - byte_off);
        bind_pages((char *)base + byte_off, len, domain);
    }
}

static void apply_numa(void *ptr, size_t tot, int N, int SB, bool blk, bool is_out, NumaPolicy pol,
                       int D) {
    if (D <= 1 || pol == NP_NONE)
        return;
    if (!blk) {
        /* Row-major: input by rows, output by columns */
        if (is_out)
            bind_by_cols(ptr, N, D);
        else
            bind_contiguous(ptr, tot, D);
        return;
    }
    int NB = N / SB, NB2 = NB * NB;
    size_t bb = (size_t)SB * SB * sizeof(float);
    int per = (NB2 + D - 1) / D;
    std::vector<int> nm(NB2);
    if (pol == NP_CONTIG)
        for (int b = 0; b < NB2; b++)
            nm[b] = std::min(b / per, D - 1);
    else if (pol == NP_CYCLIC)
        for (int b = 0; b < NB2; b++)
            nm[b] = b % D;
    else if (pol == NP_TRANSPOSE) {
        if (!is_out)
            for (int b = 0; b < NB2; b++)
                nm[b] = std::min(b / per, D - 1);
        else
            for (int ob = 0; ob < NB2; ob++) {
                int bc = ob / NB, br = ob % NB;
                nm[ob] = std::min((br * NB + bc) / per, D - 1);
            }
    }
    bind_by_map(ptr, NB2, bb, nm);
}
/* Block-row NUMA binding for row-major data with blocked schedule.
 * Input: binds rows [br*SB, (br+1)*SB) as a unit to a NUMA domain.
 * Output: binds columns to NUMA domains (transpose target locality). */
static void apply_numa_rm_blk(void *ptr, size_t /*tot*/, int N, int SB, bool is_out, NumaPolicy pol,
                              int D) {
    if (D <= 1 || pol == NP_NONE)
        return;

    /* Output: column-based binding (same as row-major output) */
    if (is_out) {
        bind_by_cols(ptr, N, D);
        return;
    }

    /* Input: block-row-based binding */
    int NB = N / SB;
    size_t row_bytes = (size_t)N * sizeof(float);
    size_t blk_row_bytes = (size_t)SB * row_bytes;
    int per = (NB + D - 1) / D;
    size_t ps = pagesz();

    for (int br = 0; br < NB; br++) {
        int domain;
        if (pol == NP_CONTIG)
            domain = std::min(br / per, D - 1);
        else if (pol == NP_CYCLIC)
            domain = br % D;
        else
            domain = std::min(br / per, D - 1);

        size_t lo = (size_t)br * blk_row_bytes;
        size_t hi = lo + blk_row_bytes;
        size_t a_lo = (lo / ps) * ps;
        size_t a_hi = ((hi + ps - 1) / ps) * ps;
        bind_pages((char *)ptr + a_lo, a_hi - a_lo, domain);
    }
}

static void parallel_zero(float *p, size_t n) {
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++)
        p[i] = 0.0f;
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Templated kernels — all tile dims are constexpr
 * ═══════════════════════════════════════════════════════════════════════ */

/* V0/V1: naive (no template params needed) */
static void tr_naive(const float *__restrict__ i, float *__restrict__ o, int N) {
#pragma omp parallel for schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            o[c * N + r] = i[r * N + c];
}
static void tr_naive_c2(const float *__restrict__ i, float *__restrict__ o, int N) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            o[c * N + r] = i[r * N + c];
}

/* V2/V3: tiled<TB> */
template <int TB>
static void tr_tiled(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NT = (N + TB - 1) / TB;
#pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
#pragma unroll
            for (int lr = 0; lr < TB; lr++) {
                int r = r0 + lr;
                if (r >= N)
                    break;
#pragma unroll
                for (int lc = 0; lc < TB; lc++) {
                    int c = c0 + lc;
                    if (c >= N)
                        break;
                    out[c * N + r] = in[r * N + c];
                }
            }
        }
}
template <int TB>
static void tr_tiled_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NT = (N + TB - 1) / TB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
#pragma unroll
            for (int lr = 0; lr < TB; lr++) {
                int r = r0 + lr;
                if (r >= N)
                    break;
#pragma unroll
                for (int lc = 0; lc < TB; lc++) {
                    int c = c0 + lc;
                    if (c >= N)
                        break;
                    out[c * N + r] = in[r * N + c];
                }
            }
        }
}

/* V4/V5: blk_naive<SB> */
template <int SB>
static void tr_blk_naive(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            int si = (r / SB * NB + c / SB) * SB * SB + (r % SB) * SB + c % SB;
            int di = (c / SB * NB + r / SB) * SB * SB + (c % SB) * SB + r % SB;
            out[di] = in[si];
        }
}
template <int SB>
static void tr_blk_naive_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            int si = (r / SB * NB + c / SB) * SB * SB + (r % SB) * SB + c % SB;
            int di = (c / SB * NB + r / SB) * SB * SB + (c % SB) * SB + r % SB;
            out[di] = in[si];
        }
}

/* V6/V7: blk_tiled<SB,TB> */
template <int SB, int TB>
static void tr_blk_tiled(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB, NT = (N + TB - 1) / TB;
#pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
            for (int r = r0; r < std::min(r0 + TB, N); r++)
                for (int c = c0; c < std::min(c0 + TB, N); c++) {
                    int si = (r / SB * NB + c / SB) * SB * SB + (r % SB) * SB + c % SB;
                    int di = (c / SB * NB + r / SB) * SB * SB + (c % SB) * SB + r % SB;
                    out[di] = in[si];
                }
        }
}
template <int SB, int TB>
static void tr_blk_tiled_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB, NT = (N + TB - 1) / TB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
            for (int r = r0; r < std::min(r0 + TB, N); r++)
                for (int c = c0; c < std::min(c0 + TB, N); c++) {
                    int si = (r / SB * NB + c / SB) * SB * SB + (r % SB) * SB + c % SB;
                    int di = (c / SB * NB + r / SB) * SB * SB + (c % SB) * SB + r % SB;
                    out[di] = in[si];
                }
        }
}

/* V8/V9: blk_aligned<SB> */
template <int SB>
static void tr_blk_aligned(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            const float *s = in + (br * NB + bc) * SB * SB;
            float *d = out + (bc * NB + br) * SB * SB;
#pragma unroll
            for (int lr = 0; lr < SB; lr++)
#pragma unroll
                for (int lc = 0; lc < SB; lc++)
                    d[lc * SB + lr] = s[lr * SB + lc];
        }
}
template <int SB>
static void tr_blk_aligned_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            const float *s = in + (br * NB + bc) * SB * SB;
            float *d = out + (bc * NB + br) * SB * SB;
#pragma unroll
            for (int lr = 0; lr < SB; lr++)
#pragma unroll
                for (int lc = 0; lc < SB; lc++)
                    d[lc * SB + lr] = s[lr * SB + lc];
        }
}

/* V10/V11: blk_aligned_mt<SB,MT> */
template <int SB, int MT>
static void tr_blk_aligned_mt(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            const float *s = in + (br * NB + bc) * SB * SB;
            float *d = out + (bc * NB + br) * SB * SB;
#pragma unroll
            for (int lr0 = 0; lr0 < SB; lr0 += MT)
#pragma unroll
                for (int lc0 = 0; lc0 < SB; lc0 += MT) {
#pragma unroll
                    for (int lr = lr0; lr < lr0 + MT; lr++)
#pragma unroll
                        for (int lc = lc0; lc < lc0 + MT; lc++)
                            d[lc * SB + lr] = s[lr * SB + lc];
                }
        }
}
template <int SB, int MT>
static void tr_blk_aligned_mt_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            const float *s = in + (br * NB + bc) * SB * SB;
            float *d = out + (bc * NB + br) * SB * SB;
#pragma unroll
            for (int lr0 = 0; lr0 < SB; lr0 += MT)
#pragma unroll
                for (int lc0 = 0; lc0 < SB; lc0 += MT) {
#pragma unroll
                    for (int lr = lr0; lr < lr0 + MT; lr++)
#pragma unroll
                        for (int lc = lc0; lc < lc0 + MT; lc++)
                            d[lc * SB + lr] = s[lr * SB + lc];
                }
        }
}

/* V12/V13: locbuf<TB> */
template <int TB>
static void tr_locbuf(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NT = (N + TB - 1) / TB;
#pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        float buf[TB][TB];
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB, bh = std::min(TB, N - r0), bw = std::min(TB, N - c0);
            for (int lr = 0; lr < bh; lr++) {
                const float *row = in + (r0 + lr) * N + c0;
#pragma unroll
                for (int lc = 0; lc < TB; lc++)
                    if (lc < bw)
                        buf[lr][lc] = row[lc];
            }
            for (int lc = 0; lc < bw; lc++) {
                float *dr = out + (c0 + lc) * N + r0;
#pragma unroll
                for (int lr = 0; lr < TB; lr++)
                    if (lr < bh)
                        dr[lr] = buf[lr][lc];
            }
        }
    }
}
template <int TB>
static void tr_locbuf_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NT = (N + TB - 1) / TB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            float buf[TB][TB];
            int r0 = tr * TB, c0 = tc * TB, bh = std::min(TB, N - r0), bw = std::min(TB, N - c0);
            for (int lr = 0; lr < bh; lr++) {
                const float *row = in + (r0 + lr) * N + c0;
#pragma unroll
                for (int lc = 0; lc < TB; lc++)
                    if (lc < bw)
                        buf[lr][lc] = row[lc];
            }
            for (int lc = 0; lc < bw; lc++) {
                float *dr = out + (c0 + lc) * N + r0;
#pragma unroll
                for (int lr = 0; lr < TB; lr++)
                    if (lr < bh)
                        dr[lr] = buf[lr][lc];
            }
        }
}

/* V14/V15: locbuf_blk<SB> */
template <int SB>
static void tr_locbuf_blk(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++) {
        float buf[SB][SB];
        for (int bc = 0; bc < NB; bc++) {
            const float *s = in + (br * NB + bc) * SB * SB;
            float *d = out + (bc * NB + br) * SB * SB;
#pragma unroll
            for (int lr = 0; lr < SB; lr++)
#pragma unroll
                for (int lc = 0; lc < SB; lc++)
                    buf[lr][lc] = s[lr * SB + lc];
#pragma unroll
            for (int lc = 0; lc < SB; lc++) {
                float *dr = d + lc * SB;
#pragma unroll
                for (int lr = 0; lr < SB; lr++)
                    dr[lr] = buf[lr][lc];
            }
        }
    }
}
template <int SB>
static void tr_locbuf_blk_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            float buf[SB][SB];
            const float *s = in + (br * NB + bc) * SB * SB;
            float *d = out + (bc * NB + br) * SB * SB;
#pragma unroll
            for (int lr = 0; lr < SB; lr++)
#pragma unroll
                for (int lc = 0; lc < SB; lc++)
                    buf[lr][lc] = s[lr * SB + lc];
#pragma unroll
            for (int lc = 0; lc < SB; lc++) {
                float *dr = d + lc * SB;
#pragma unroll
                for (int lr = 0; lr < SB; lr++)
                    dr[lr] = buf[lr][lc];
            }
        }
}

/* V16/V17: locbuf_blk_mt<SB,MT> */
template <int SB, int MT>
static void tr_locbuf_blk_mt(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++) {
        float buf[SB][SB];
        for (int bc = 0; bc < NB; bc++) {
            const float *s = in + (br * NB + bc) * SB * SB;
            float *d = out + (bc * NB + br) * SB * SB;
#pragma unroll
            for (int lr = 0; lr < SB; lr++)
#pragma unroll
                for (int lc = 0; lc < SB; lc++)
                    buf[lr][lc] = s[lr * SB + lc];
#pragma unroll
            for (int lc0 = 0; lc0 < SB; lc0 += MT)
#pragma unroll
                for (int lr0 = 0; lr0 < SB; lr0 += MT) {
#pragma unroll
                    for (int lc = lc0; lc < lc0 + MT; lc++) {
                        float *dr = d + lc * SB + lr0;
#pragma unroll
                        for (int lr = lr0; lr < lr0 + MT; lr++)
                            dr[lr - lr0] = buf[lr][lc];
                    }
                }
        }
    }
}
template <int SB, int MT>
static void tr_locbuf_blk_mt_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            float buf[SB][SB];
            const float *s = in + (br * NB + bc) * SB * SB;
            float *d = out + (bc * NB + br) * SB * SB;
#pragma unroll
            for (int lr = 0; lr < SB; lr++)
#pragma unroll
                for (int lc = 0; lc < SB; lc++)
                    buf[lr][lc] = s[lr * SB + lc];
#pragma unroll
            for (int lc0 = 0; lc0 < SB; lc0 += MT)
#pragma unroll
                for (int lr0 = 0; lr0 < SB; lr0 += MT) {
#pragma unroll
                    for (int lc = lc0; lc < lc0 + MT; lc++) {
                        float *dr = d + lc * SB + lr0;
#pragma unroll
                        for (int lr = lr0; lr < lr0 + MT; lr++)
                            dr[lr - lr0] = buf[lr][lc];
                    }
                }
        }
}

/* V18/V19: locbuf_2buf<TB> */
template <int TB>
static void tr_locbuf_2buf(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NT = (N + TB - 1) / TB;
#pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        float brd[TB][TB], bwr[TB][TB];
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB, bh = std::min(TB, N - r0), bw = std::min(TB, N - c0);
            for (int lr = 0; lr < bh; lr++) {
                const float *row = in + (r0 + lr) * N + c0;
                for (int lc = 0; lc < bw; lc++)
                    brd[lr][lc] = row[lc];
            }
#pragma unroll
            for (int lr = 0; lr < TB; lr++)
#pragma unroll
                for (int lc = 0; lc < TB; lc++)
                    bwr[lc][lr] = brd[lr][lc];
            for (int lc = 0; lc < bw; lc++) {
                float *dr = out + (c0 + lc) * N + r0;
                for (int lr = 0; lr < bh; lr++)
                    dr[lr] = bwr[lc][lr];
            }
        }
    }
}
template <int TB>
static void tr_locbuf_2buf_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NT = (N + TB - 1) / TB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            float brd[TB][TB], bwr[TB][TB];
            int r0 = tr * TB, c0 = tc * TB, bh = std::min(TB, N - r0), bw = std::min(TB, N - c0);
            for (int lr = 0; lr < bh; lr++) {
                const float *row = in + (r0 + lr) * N + c0;
                for (int lc = 0; lc < bw; lc++)
                    brd[lr][lc] = row[lc];
            }
#pragma unroll
            for (int lr = 0; lr < TB; lr++)
#pragma unroll
                for (int lc = 0; lc < TB; lc++)
                    bwr[lc][lr] = brd[lr][lc];
            for (int lc = 0; lc < bw; lc++) {
                float *dr = out + (c0 + lc) * N + r0;
                for (int lr = 0; lr < bh; lr++)
                    dr[lr] = bwr[lc][lr];
            }
        }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  V20-V23: Blocked schedule on ROW-MAJOR data.
 *  Same iteration order as blk_aligned variants, but arrays stay
 *  row-major:  in[(br*SB+lr)*N + bc*SB+lc]  instead of blocked indexing.
 *  Combined with blocked NUMA mapping this gives NUMA-local block access
 *  without changing the data layout.
 * ═══════════════════════════════════════════════════════════════════════ */
template <int SB>
static void tr_rm_blk(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
#pragma unroll
            for (int lr = 0; lr < SB; lr++)
#pragma unroll
                for (int lc = 0; lc < SB; lc++)
                    out[(bc * SB + lc) * N + br * SB + lr] = in[(br * SB + lr) * N + bc * SB + lc];
        }
}
template <int SB>
static void tr_rm_blk_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
#pragma unroll
            for (int lr = 0; lr < SB; lr++)
#pragma unroll
                for (int lc = 0; lc < SB; lc++)
                    out[(bc * SB + lc) * N + br * SB + lr] = in[(br * SB + lr) * N + bc * SB + lc];
        }
}
template <int SB, int MT>
static void tr_rm_blk_mt(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
#pragma unroll
            for (int lr0 = 0; lr0 < SB; lr0 += MT)
#pragma unroll
                for (int lc0 = 0; lc0 < SB; lc0 += MT) {
#pragma unroll
                    for (int lr = lr0; lr < lr0 + MT; lr++)
#pragma unroll
                        for (int lc = lc0; lc < lc0 + MT; lc++)
                            out[(bc * SB + lc) * N + br * SB + lr] =
                                in[(br * SB + lr) * N + bc * SB + lc];
                }
        }
}
template <int SB, int MT>
static void tr_rm_blk_mt_c2(const float *__restrict__ in, float *__restrict__ out, int N) {
    int NB = N / SB;
#pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
#pragma unroll
            for (int lr0 = 0; lr0 < SB; lr0 += MT)
#pragma unroll
                for (int lc0 = 0; lc0 < SB; lc0 += MT) {
#pragma unroll
                    for (int lr = lr0; lr < lr0 + MT; lr++)
#pragma unroll
                        for (int lc = lc0; lc < lc0 + MT; lc++)
                            out[(bc * SB + lc) * N + br * SB + lr] =
                                in[(br * SB + lr) * N + bc * SB + lc];
                }
        }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Dispatch — runtime (TB,SB,MT) → template instantiation
 * ═══════════════════════════════════════════════════════════════════════ */
static const char *V_NAMES[] = {
    "naive",         "naive_c2",         "tiled",          "tiled_c2",
    "blk_naive",     "blk_naive_c2",     "blk_tiled",      "blk_tiled_c2",
    "blk_aligned",   "blk_aligned_c2",   "blk_aligned_mt", "blk_aligned_mt_c2",
    "locbuf",        "locbuf_c2",        "locbuf_blk",     "locbuf_blk_c2",
    "locbuf_blk_mt", "locbuf_blk_mt_c2", "locbuf_2buf",    "locbuf_2buf_c2",
    "rm_blk",        "rm_blk_c2",        "rm_blk_mt",      "rm_blk_mt_c2"};
static const int N_VARIANTS = 24;
/* Data is stored in blocked layout */
static bool is_blocked(int v) {
    return (v >= 4 && v <= 11) || (v >= 14 && v <= 17);
}
/* Row-major data but blocked schedule — allows blocked NUMA policies */
static bool is_rm_blk_sched(int v) {
    return v >= 20 && v <= 23;
}

#define D_TB(V, TB_, i, o, N)                                                                      \
    switch (V) {                                                                                   \
    case 2:                                                                                        \
        tr_tiled<TB_>(i, o, N);                                                                    \
        return;                                                                                    \
    case 3:                                                                                        \
        tr_tiled_c2<TB_>(i, o, N);                                                                 \
        return;                                                                                    \
    case 12:                                                                                       \
        tr_locbuf<TB_>(i, o, N);                                                                   \
        return;                                                                                    \
    case 13:                                                                                       \
        tr_locbuf_c2<TB_>(i, o, N);                                                                \
        return;                                                                                    \
    case 18:                                                                                       \
        tr_locbuf_2buf<TB_>(i, o, N);                                                              \
        return;                                                                                    \
    case 19:                                                                                       \
        tr_locbuf_2buf_c2<TB_>(i, o, N);                                                           \
        return;                                                                                    \
    default:                                                                                       \
        break;                                                                                     \
    }
#define D_SB(V, SB_, i, o, N)                                                                      \
    switch (V) {                                                                                   \
    case 4:                                                                                        \
        tr_blk_naive<SB_>(i, o, N);                                                                \
        return;                                                                                    \
    case 5:                                                                                        \
        tr_blk_naive_c2<SB_>(i, o, N);                                                             \
        return;                                                                                    \
    case 8:                                                                                        \
        tr_blk_aligned<SB_>(i, o, N);                                                              \
        return;                                                                                    \
    case 9:                                                                                        \
        tr_blk_aligned_c2<SB_>(i, o, N);                                                           \
        return;                                                                                    \
    case 14:                                                                                       \
        tr_locbuf_blk<SB_>(i, o, N);                                                               \
        return;                                                                                    \
    case 15:                                                                                       \
        tr_locbuf_blk_c2<SB_>(i, o, N);                                                            \
        return;                                                                                    \
    case 20:                                                                                       \
        tr_rm_blk<SB_>(i, o, N);                                                                   \
        return;                                                                                    \
    case 21:                                                                                       \
        tr_rm_blk_c2<SB_>(i, o, N);                                                                \
        return;                                                                                    \
    default:                                                                                       \
        break;                                                                                     \
    }
#define D_SB_TB(V, SB_, TB_, i, o, N)                                                              \
    switch (V) {                                                                                   \
    case 6:                                                                                        \
        tr_blk_tiled<SB_, TB_>(i, o, N);                                                           \
        return;                                                                                    \
    case 7:                                                                                        \
        tr_blk_tiled_c2<SB_, TB_>(i, o, N);                                                        \
        return;                                                                                    \
    default:                                                                                       \
        break;                                                                                     \
    }
#define D_SB_MT(V, SB_, MT_, i, o, N)                                                              \
    switch (V) {                                                                                   \
    case 10:                                                                                       \
        tr_blk_aligned_mt<SB_, MT_>(i, o, N);                                                      \
        return;                                                                                    \
    case 11:                                                                                       \
        tr_blk_aligned_mt_c2<SB_, MT_>(i, o, N);                                                   \
        return;                                                                                    \
    case 16:                                                                                       \
        tr_locbuf_blk_mt<SB_, MT_>(i, o, N);                                                       \
        return;                                                                                    \
    case 17:                                                                                       \
        tr_locbuf_blk_mt_c2<SB_, MT_>(i, o, N);                                                    \
        return;                                                                                    \
    case 22:                                                                                       \
        tr_rm_blk_mt<SB_, MT_>(i, o, N);                                                           \
        return;                                                                                    \
    case 23:                                                                                       \
        tr_rm_blk_mt_c2<SB_, MT_>(i, o, N);                                                        \
        return;                                                                                    \
    default:                                                                                       \
        break;                                                                                     \
    }

static void dispatch(int V, int tb, int sb, int mt, const float *in, float *out, int N) {
    if (V == 0) {
        tr_naive(in, out, N);
        return;
    }
    if (V == 1) {
        tr_naive_c2(in, out, N);
        return;
    }
    /* TB-only */
    if (V == 2 || V == 3 || V == 12 || V == 13 || V == 18 || V == 19) {
        if (tb == 16) {
            D_TB(V, 16, in, out, N)
        }
        if (tb == 32) {
            D_TB(V, 32, in, out, N)
        }
        if (tb == 64) {
            D_TB(V, 64, in, out, N)
        }
        if (tb == 128) {
            D_TB(V, 128, in, out, N)
        }
        fprintf(stderr, "No TB=%d instantiation\n", tb);
        exit(1);
    }
    /* SB-only */
    if (V == 4 || V == 5 || V == 8 || V == 9 || V == 14 || V == 15 || V == 20 || V == 21) {
#define T(S)                                                                                       \
    if (sb == S) {                                                                                 \
        D_SB(V, S, in, out, N)                                                                     \
    }
        T(8)
        T(16)
        T(32) T(64) T(128) T(256) T(512)
#undef T
            fprintf(stderr, "No SB=%d instantiation\n", sb);
        exit(1);
    }
    /* SB+TB */
    if (V == 6 || V == 7) {
#define T(S, B)                                                                                    \
    if (sb == S && tb == B) {                                                                      \
        D_SB_TB(V, S, B, in, out, N)                                                               \
    }
        T(8, 64)
        T(16, 64)
        T(32, 64) T(64, 64) T(128, 64) T(256, 64) T(512, 64) T(32, 32) T(64, 32) T(128, 128)
#undef T
            fprintf(stderr, "No SB=%d TB=%d instantiation\n", sb, tb);
        exit(1);
    }
    /* SB+MT */
    if (V == 10 || V == 11 || V == 16 || V == 17 || V == 22 || V == 23) {
#define T(S, M)                                                                                    \
    if (sb == S && mt == M) {                                                                      \
        D_SB_MT(V, S, M, in, out, N)                                                               \
    }
        T(8, 4)
        T(8, 8)
        T(16, 4) T(16, 8) T(16, 16) T(32, 4) T(32, 8) T(32, 16) T(32, 32) T(64, 4) T(64, 8)
            T(64, 16) T(64, 32) T(128, 4) T(128, 8) T(128, 16) T(128, 32) T(256, 4) T(256, 8)
                T(256, 16) T(256, 32) T(512, 4) T(512, 8) T(512, 16) T(512, 32)
#undef T
                    fprintf(stderr, "No SB=%d MT=%d instantiation\n", sb, mt);
        exit(1);
    }
    fprintf(stderr, "Unhandled variant %d\n", V);
    exit(1);
}

/* ── Verification ───────────────────────────────────────────────────── */
static void ref_transpose(const float *in, float *ref, int N) {
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            ref[c * N + r] = in[r * N + c];
}
static float verify(const float *out, const float *ref, int N, int SB, bool blk) {
    float mx = 0;
    if (blk) {
        int NB = N / SB;
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++) {
                float e = fabsf(out[(r / SB * NB + c / SB) * SB * SB + (r % SB) * SB + c % SB] -
                                ref[r * N + c]);
                if (e > mx)
                    mx = e;
            }
    } else
        for (size_t i = 0; i < (size_t)N * N; i++) {
            float e = fabsf(out[i] - ref[i]);
            if (e > mx)
                mx = e;
        }
    return mx;
}

/* ═══════════════════════════════════════════════════════════════════════
 *  main
 * ═══════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <N> <variant> <csv> [TB=64] [SB=32] [MT=8] "
                "[WARMUP=3] [REPS=20] [THREADS=0] [NUMA=0]\n",
                argv[0]);
        return 1;
    }
    int N = atoi(argv[1]), VAR = atoi(argv[2]);
    const char *csv = argv[3];
    int TB = argc > 4 ? atoi(argv[4]) : 64, SB = argc > 5 ? atoi(argv[5]) : 32,
        MT = argc > 6 ? atoi(argv[6]) : 8;
    int WARMUP = argc > 7 ? atoi(argv[7]) : 3, REPS = argc > 8 ? atoi(argv[8]) : 20;
    int THREADS = argc > 9 ? atoi(argv[9]) : 0, NPOL = argc > 10 ? atoi(argv[10]) : 0;
    if (VAR < 0 || VAR >= N_VARIANTS) {
        fprintf(stderr, "bad variant\n");
        return 1;
    }
    bool blocked = is_blocked(VAR);
    bool rm_blk = is_rm_blk_sched(VAR);
    if ((blocked || rm_blk) && N % SB) {
        fprintf(stderr, "N%%SB!=0\n");
        return 1;
    }
    if (NPOL < 0 || NPOL >= NP_COUNT) {
        fprintf(stderr, "bad NUMA\n");
        return 1;
    }
    if ((NPOL == NP_CYCLIC || NPOL == NP_TRANSPOSE) && !blocked && !rm_blk) {
        fprintf(stderr, "cyclic/transpose NUMA only for blocked or rm_blk variants\n");
        return 1;
    }
    NumaPolicy npol = (NumaPolicy)NPOL;
    if (THREADS > 0)
        omp_set_num_threads(THREADS);
    int nthreads;
#pragma omp parallel
    {
#pragma omp single
        nthreads = omp_get_num_threads();
    }
    int D = detect_numa_nodes();
    char vname[128];
    snprintf(vname, sizeof(vname), "%s%s", V_NAMES[VAR], NP_SUFFIX[NPOL]);
    size_t elems = (size_t)N * N, bytes = elems * sizeof(float);
    float *h_row = (float *)aligned_alloc(64, bytes);
    float *h_ref = (float *)aligned_alloc(64, bytes);
    bool use_numa = (npol != NP_NONE);
    float *h_in = use_numa ? numa_alloc<float>(elems) : (float *)aligned_alloc(64, bytes);
    float *h_out = use_numa ? numa_alloc<float>(elems) : (float *)aligned_alloc(64, bytes);
    if (use_numa) {
        if (rm_blk) {
            /* Row-major data with block-row NUMA binding */
            apply_numa_rm_blk(h_in, bytes, N, SB, false, npol, D);
            apply_numa_rm_blk(h_out, bytes, N, SB, true, npol, D);
        } else {
            apply_numa(h_in, bytes, N, SB, blocked, false, npol, D);
            apply_numa(h_out, bytes, N, SB, blocked, true, npol, D);
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < elems; i++)
        h_row[i] = (float)i / (float)N;
    ref_transpose(h_row, h_ref, N);
    if (blocked) {
        int NB = N / SB;
#pragma omp parallel for schedule(static)
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                h_in[(r / SB * NB + c / SB) * SB * SB + (r % SB) * SB + c % SB] = h_row[r * N + c];
    } else {
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < elems; i++)
            h_in[i] = h_row[i];
    }
    parallel_zero(h_out, elems);
    printf("  %s N=%d TB=%d SB=%d MT=%d thr=%d numa=%d nodes=%d\n", vname, N, TB, SB, MT, nthreads,
           NPOL, D);
    for (int i = 0; i < WARMUP; i++)
        dispatch(VAR, TB, SB, MT, h_in, h_out, N);
    parallel_zero(h_out, elems);
    dispatch(VAR, TB, SB, MT, h_in, h_out, N);
    float maxerr = verify(h_out, h_ref, N, SB, blocked);
    bool pass = (maxerr == 0.0f);
    double *times = (double *)malloc(REPS * sizeof(double));
    for (int i = 0; i < REPS; i++) {
        double t0 = omp_get_wtime();
        dispatch(VAR, TB, SB, MT, h_in, h_out, N);
        times[i] = omp_get_wtime() - t0;
    }
    std::sort(times, times + REPS);
    double bpi = 2.0 * N * (double)N * sizeof(float), med_s = times[REPS / 2];
    double cksum = 0;
    for (size_t i = 0; i < elems; i++)
        cksum += h_out[i];
    printf("%s N=%d TB=%d SB=%d MT=%d thr=%d | med %.4f ms (%.1f GB/s) "
           "p5 %.4f p95 %.4f maxerr=%.1e %s cksum=%.6e\n",
           vname, N, TB, SB, MT, nthreads, med_s * 1e3, bpi / med_s / 1e9,
           times[(int)(REPS * .05)] * 1e3, times[(int)(REPS * .95)] * 1e3, maxerr,
           pass ? "PASS" : "FAIL", cksum);
    FILE *f = fopen(csv, "a");
    if (f) {
        for (int i = 0; i < REPS; i++)
            fprintf(f, "%s,%d,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e,%s\n", vname, N, TB, SB, MT, nthreads,
                    i, times[i], bpi / times[i] / 1e9, cksum, pass ? "PASS" : "FAIL");
        fclose(f);
    }
    free(times);
    free(h_row);
    free(h_ref);
    if (use_numa) {
        numa_dealloc(h_in, elems);
        numa_dealloc(h_out, elems);
    } else {
        free(h_in);
        free(h_out);
    }
    return pass ? 0 : 1;
}