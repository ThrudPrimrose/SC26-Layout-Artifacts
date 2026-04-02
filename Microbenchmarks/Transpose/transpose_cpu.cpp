/*  transpose_cpu.cpp — CPU matrix transpose sweep with NUMA-aware allocation.
 *
 *  NUMA policies (last CLI arg):
 *    0 = none       (aligned_alloc, backward compat)
 *    1 = contiguous (_nd)  rows or block-ranges split across NUMA domains
 *    2 = cyclic     (_nc)  block_id % D   (blocked variants only)
 *    3 = transpose  (_nt)  in(i,j) & out(j,i) on same domain (blocked only)
 *
 *  Compile:  g++ -O3 -march=native -fopenmp -o transpose_cpu transpose_cpu.cpp
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>

/* ═══════════════════════════════════════════════════════════════════════
 *  NUMA helpers  (no libnuma — direct syscalls)
 * ═══════════════════════════════════════════════════════════════════════ */
#ifndef MPOL_BIND
#define MPOL_BIND 2
#endif

static long sys_mbind(void *addr, unsigned long len, int mode,
                      const unsigned long *nodemask, unsigned long maxnode,
                      unsigned flags) {
    return syscall(SYS_mbind, addr, len, mode, nodemask, maxnode, flags);
}

static int detect_numa_nodes() {
    int n = 0;
    for (int i = 0; i < 128; i++) {
        char path[128];
        snprintf(path, sizeof(path), "/sys/devices/system/node/node%d", i);
        if (access(path, F_OK) == 0) n = i + 1;
        else if (n > 0) break;
    }
    return n > 0 ? n : 1;
}

static void bind_pages(void *addr, size_t len, int node) {
    if (len == 0 || node < 0) return;
    unsigned long mask[4] = {};
    mask[node / 64] |= 1UL << (node % 64);
    sys_mbind(addr, len, MPOL_BIND, mask, 257, 0);
}

template <typename T>
static T *numa_alloc(size_t count) {
    size_t bytes = count * sizeof(T);
    void *p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); std::abort(); }
    madvise(p, bytes, MADV_NOHUGEPAGE);
    return static_cast<T *>(p);
}

template <typename T>
static void numa_dealloc(T *p, size_t count) {
    if (p) munmap(p, count * sizeof(T));
}

/* ── NUMA binding ─────────────────────────────────────────────────── */

enum NumaPolicy { NP_NONE = 0, NP_CONTIG = 1, NP_CYCLIC = 2, NP_TRANSPOSE = 3 };
static const char *NP_SUFFIX[] = { "", "_nd", "_nc", "_nt" };
static const int   NP_COUNT    = 4;

static size_t g_page_size = 0;
static size_t pagesz() {
    if (!g_page_size) g_page_size = (size_t)sysconf(_SC_PAGESIZE);
    return g_page_size;
}

/* Bind contiguous byte range split evenly across D domains. */
static void bind_contiguous(void *base, size_t total_bytes, int D) {
    size_t ps = pagesz();
    for (int d = 0; d < D; d++) {
        size_t lo = (size_t)d       * total_bytes / D;
        size_t hi = (size_t)(d + 1) * total_bytes / D;
        // page-align
        lo = (lo / ps) * ps;
        hi = ((hi + ps - 1) / ps) * ps;
        if (hi > total_bytes) hi = ((total_bytes + ps - 1) / ps) * ps;
        if (hi > lo)
            bind_pages(static_cast<char*>(base) + lo, hi - lo, d);
    }
}

/* Bind blocked array where node_map[bid] gives the NUMA node for block bid.
 * Coalesces consecutive same-node blocks into single mbind calls. */
static void bind_by_map(void *base, int NB2, size_t blk_bytes,
                        const std::vector<int> &node_map) {
    size_t ps = pagesz();
    int cur_node = node_map[0];
    size_t run_start = 0;

    for (int bid = 1; bid <= NB2; bid++) {
        int node = (bid < NB2) ? node_map[bid] : -1;
        if (node != cur_node) {
            size_t run_end = (size_t)bid * blk_bytes;
            size_t a_start = (run_start / ps) * ps;
            size_t a_end   = ((run_end + ps - 1) / ps) * ps;
            bind_pages(static_cast<char*>(base) + a_start, a_end - a_start, cur_node);
            cur_node = node;
            run_start = run_end;
        }
    }
}

/* Apply NUMA policy to one array (input or output). */
static void apply_numa(void *ptr, size_t total_bytes,
                       int N, int SB, bool is_blocked, bool is_output,
                       NumaPolicy pol, int D) {
    if (D <= 1 || pol == NP_NONE) return;

    if (!is_blocked) {
        /* Row-major: all NUMA policies degenerate to contiguous split. */
        bind_contiguous(ptr, total_bytes, D);
        return;
    }

    int NB  = N / SB;
    int NB2 = NB * NB;
    size_t blk_bytes = (size_t)SB * SB * sizeof(float);
    int per = (NB2 + D - 1) / D;          // blocks per domain (ceil)

    std::vector<int> nmap(NB2);

    switch (pol) {
    case NP_CONTIG:
        for (int b = 0; b < NB2; b++)
            nmap[b] = std::min(b / per, D - 1);
        break;

    case NP_CYCLIC:
        for (int b = 0; b < NB2; b++)
            nmap[b] = b % D;
        break;

    case NP_TRANSPOSE:
        if (!is_output) {
            /* Input: contiguous split of blocks. */
            for (int b = 0; b < NB2; b++)
                nmap[b] = std::min(b / per, D - 1);
        } else {
            /* Output: block (bc, br) in storage gets same domain as
             * input block (br, bc). Storage linear id = bc*NB + br,
             * input linear id = br*NB + bc. */
            for (int out_bid = 0; out_bid < NB2; out_bid++) {
                int bc = out_bid / NB;
                int br = out_bid % NB;
                int in_bid = br * NB + bc;
                nmap[out_bid] = std::min(in_bid / per, D - 1);
            }
        }
        break;

    default: break;
    }

    bind_by_map(ptr, NB2, blk_bytes, nmap);
}

/* Parallel first-touch helpers. */
static void parallel_zero(float *p, size_t n) {
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) p[i] = 0.0f;
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Row-major variants (V0–V3)
 * ═══════════════════════════════════════════════════════════════════════ */
static void tr_naive(const float* __restrict__ in, float* __restrict__ out, int N) {
    #pragma omp parallel for schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            out[c * N + r] = in[r * N + c];
}
static void tr_naive_c2(const float* __restrict__ in, float* __restrict__ out, int N) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            out[c * N + r] = in[r * N + c];
}
static void tr_tiled(const float* __restrict__ in, float* __restrict__ out, int N, int TB) {
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr*TB, c0 = tc*TB;
            int re = std::min(r0+TB, N), ce = std::min(c0+TB, N);
            for (int r = r0; r < re; r++)
                for (int c = c0; c < ce; c++)
                    out[c*N+r] = in[r*N+c];
        }
}
static void tr_tiled_c2(const float* __restrict__ in, float* __restrict__ out, int N, int TB) {
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr*TB, c0 = tc*TB;
            int re = std::min(r0+TB, N), ce = std::min(c0+TB, N);
            for (int r = r0; r < re; r++)
                for (int c = c0; c < ce; c++)
                    out[c*N+r] = in[r*N+c];
        }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Blocked-layout variants (V4–V11)
 * ═══════════════════════════════════════════════════════════════════════ */
static void tr_blk_naive(const float* __restrict__ in, float* __restrict__ out, int N, int SB) {
    int NB = N/SB;
    #pragma omp parallel for schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            int si = (r/SB*NB + c/SB)*SB*SB + (r%SB)*SB + c%SB;
            int di = (c/SB*NB + r/SB)*SB*SB + (c%SB)*SB + r%SB;
            out[di] = in[si];
        }
}
static void tr_blk_naive_c2(const float* __restrict__ in, float* __restrict__ out, int N, int SB) {
    int NB = N/SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            int si = (r/SB*NB + c/SB)*SB*SB + (r%SB)*SB + c%SB;
            int di = (c/SB*NB + r/SB)*SB*SB + (c%SB)*SB + r%SB;
            out[di] = in[si];
        }
}
static void tr_blk_tiled(const float* __restrict__ in, float* __restrict__ out,
                          int N, int SB, int TB) {
    int NB = N/SB, NT = (N+TB-1)/TB;
    #pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr*TB, c0 = tc*TB;
            int re = std::min(r0+TB, N), ce = std::min(c0+TB, N);
            for (int r = r0; r < re; r++)
                for (int c = c0; c < ce; c++) {
                    int si = (r/SB*NB + c/SB)*SB*SB + (r%SB)*SB + c%SB;
                    int di = (c/SB*NB + r/SB)*SB*SB + (c%SB)*SB + r%SB;
                    out[di] = in[si];
                }
        }
}
static void tr_blk_tiled_c2(const float* __restrict__ in, float* __restrict__ out,
                              int N, int SB, int TB) {
    int NB = N/SB, NT = (N+TB-1)/TB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr*TB, c0 = tc*TB;
            int re = std::min(r0+TB, N), ce = std::min(c0+TB, N);
            for (int r = r0; r < re; r++)
                for (int c = c0; c < ce; c++) {
                    int si = (r/SB*NB + c/SB)*SB*SB + (r%SB)*SB + c%SB;
                    int di = (c/SB*NB + r/SB)*SB*SB + (c%SB)*SB + r%SB;
                    out[di] = in[si];
                }
        }
}
static void tr_blk_aligned(const float* __restrict__ in, float* __restrict__ out, int N, int SB) {
    int NB = N/SB;
    #pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br*NB+bc)*SB*SB;
            float*       dst = out + (bc*NB+br)*SB*SB;
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    dst[lc*SB+lr] = src[lr*SB+lc];
        }
}
static void tr_blk_aligned_c2(const float* __restrict__ in, float* __restrict__ out, int N, int SB) {
    int NB = N/SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br*NB+bc)*SB*SB;
            float*       dst = out + (bc*NB+br)*SB*SB;
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    dst[lc*SB+lr] = src[lr*SB+lc];
        }
}
static void tr_blk_aligned_mt(const float* __restrict__ in, float* __restrict__ out,
                                int N, int SB, int MT) {
    int NB = N/SB;
    #pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br*NB+bc)*SB*SB;
            float*       dst = out + (bc*NB+br)*SB*SB;
            for (int lr0 = 0; lr0 < SB; lr0 += MT)
                for (int lc0 = 0; lc0 < SB; lc0 += MT) {
                    int lre = std::min(lr0+MT, SB), lce = std::min(lc0+MT, SB);
                    for (int lr = lr0; lr < lre; lr++)
                        for (int lc = lc0; lc < lce; lc++)
                            dst[lc*SB+lr] = src[lr*SB+lc];
                }
        }
}
static void tr_blk_aligned_mt_c2(const float* __restrict__ in, float* __restrict__ out,
                                   int N, int SB, int MT) {
    int NB = N/SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br*NB+bc)*SB*SB;
            float*       dst = out + (bc*NB+br)*SB*SB;
            for (int lr0 = 0; lr0 < SB; lr0 += MT)
                for (int lc0 = 0; lc0 < SB; lc0 += MT) {
                    int lre = std::min(lr0+MT, SB), lce = std::min(lc0+MT, SB);
                    for (int lr = lr0; lr < lre; lr++)
                        for (int lc = lc0; lc < lce; lc++)
                            dst[lc*SB+lr] = src[lr*SB+lc];
                }
        }
}

/* ═══════════════════════════════════════════════════════════════════════
 *  Local-buffer variants (V12–V19)
 * ═══════════════════════════════════════════════════════════════════════ */
static constexpr int MAX_TB = 128;

static void tr_locbuf(const float* __restrict__ in, float* __restrict__ out, int N, int TB) {
    int NT = (N+TB-1)/TB;
    #pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        float buf[MAX_TB][MAX_TB];
        for (int tc = 0; tc < NT; tc++) {
            int r0=tr*TB, c0=tc*TB, bh=std::min(TB,N-r0), bw=std::min(TB,N-c0);
            for (int lr=0; lr<bh; lr++) { const float* row=in+(r0+lr)*N+c0; for(int lc=0;lc<bw;lc++) buf[lr][lc]=row[lc]; }
            for (int lc=0; lc<bw; lc++) { float* dr=out+(c0+lc)*N+r0; for(int lr=0;lr<bh;lr++) dr[lr]=buf[lr][lc]; }
        }
    }
}
static void tr_locbuf_c2(const float* __restrict__ in, float* __restrict__ out, int N, int TB) {
    int NT = (N+TB-1)/TB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            float buf[MAX_TB][MAX_TB];
            int r0=tr*TB, c0=tc*TB, bh=std::min(TB,N-r0), bw=std::min(TB,N-c0);
            for (int lr=0; lr<bh; lr++) { const float* row=in+(r0+lr)*N+c0; for(int lc=0;lc<bw;lc++) buf[lr][lc]=row[lc]; }
            for (int lc=0; lc<bw; lc++) { float* dr=out+(c0+lc)*N+r0; for(int lr=0;lr<bh;lr++) dr[lr]=buf[lr][lc]; }
        }
}
static void tr_locbuf_blk(const float* __restrict__ in, float* __restrict__ out, int N, int SB) {
    int NB = N/SB;
    #pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++) {
        float buf[MAX_TB][MAX_TB];
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in+(br*NB+bc)*SB*SB; float* dst = out+(bc*NB+br)*SB*SB;
            for (int lr=0;lr<SB;lr++) for(int lc=0;lc<SB;lc++) buf[lr][lc]=src[lr*SB+lc];
            for (int lc=0;lc<SB;lc++) { float* dr=dst+lc*SB; for(int lr=0;lr<SB;lr++) dr[lr]=buf[lr][lc]; }
        }
    }
}
static void tr_locbuf_blk_c2(const float* __restrict__ in, float* __restrict__ out, int N, int SB) {
    int NB = N/SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            float buf[MAX_TB][MAX_TB];
            const float* src = in+(br*NB+bc)*SB*SB; float* dst = out+(bc*NB+br)*SB*SB;
            for (int lr=0;lr<SB;lr++) for(int lc=0;lc<SB;lc++) buf[lr][lc]=src[lr*SB+lc];
            for (int lc=0;lc<SB;lc++) { float* dr=dst+lc*SB; for(int lr=0;lr<SB;lr++) dr[lr]=buf[lr][lc]; }
        }
}
static void tr_locbuf_blk_mt(const float* __restrict__ in, float* __restrict__ out,
                               int N, int SB, int MT) {
    int NB = N/SB;
    #pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++) {
        float buf[MAX_TB][MAX_TB];
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in+(br*NB+bc)*SB*SB; float* dst = out+(bc*NB+br)*SB*SB;
            for (int lr=0;lr<SB;lr++) for(int lc=0;lc<SB;lc++) buf[lr][lc]=src[lr*SB+lc];
            for (int lc0=0; lc0<SB; lc0+=MT) for (int lr0=0; lr0<SB; lr0+=MT) {
                int lce=std::min(lc0+MT,SB), lre=std::min(lr0+MT,SB);
                for(int lc=lc0;lc<lce;lc++) { float*dr=dst+lc*SB+lr0; for(int lr=lr0;lr<lre;lr++) dr[lr-lr0]=buf[lr][lc]; }
            }
        }
    }
}
static void tr_locbuf_blk_mt_c2(const float* __restrict__ in, float* __restrict__ out,
                                  int N, int SB, int MT) {
    int NB = N/SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++)
        for (int bc = 0; bc < NB; bc++) {
            float buf[MAX_TB][MAX_TB];
            const float* src = in+(br*NB+bc)*SB*SB; float* dst = out+(bc*NB+br)*SB*SB;
            for (int lr=0;lr<SB;lr++) for(int lc=0;lc<SB;lc++) buf[lr][lc]=src[lr*SB+lc];
            for (int lc0=0; lc0<SB; lc0+=MT) for (int lr0=0; lr0<SB; lr0+=MT) {
                int lce=std::min(lc0+MT,SB), lre=std::min(lr0+MT,SB);
                for(int lc=lc0;lc<lce;lc++) { float*dr=dst+lc*SB+lr0; for(int lr=lr0;lr<lre;lr++) dr[lr-lr0]=buf[lr][lc]; }
            }
        }
}
static void tr_locbuf_2buf(const float* __restrict__ in, float* __restrict__ out, int N, int TB) {
    int NT = (N+TB-1)/TB;
    #pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        float brd[MAX_TB][MAX_TB], bwr[MAX_TB][MAX_TB];
        for (int tc = 0; tc < NT; tc++) {
            int r0=tr*TB, c0=tc*TB, bh=std::min(TB,N-r0), bw=std::min(TB,N-c0);
            for(int lr=0;lr<bh;lr++){const float*row=in+(r0+lr)*N+c0; for(int lc=0;lc<bw;lc++) brd[lr][lc]=row[lc];}
            for(int lr=0;lr<bh;lr++) for(int lc=0;lc<bw;lc++) bwr[lc][lr]=brd[lr][lc];
            for(int lc=0;lc<bw;lc++){float*dr=out+(c0+lc)*N+r0; for(int lr=0;lr<bh;lr++) dr[lr]=bwr[lc][lr];}
        }
    }
}
static void tr_locbuf_2buf_c2(const float* __restrict__ in, float* __restrict__ out, int N, int TB) {
    int NT = (N+TB-1)/TB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++)
        for (int tc = 0; tc < NT; tc++) {
            float brd[MAX_TB][MAX_TB], bwr[MAX_TB][MAX_TB];
            int r0=tr*TB, c0=tc*TB, bh=std::min(TB,N-r0), bw=std::min(TB,N-c0);
            for(int lr=0;lr<bh;lr++){const float*row=in+(r0+lr)*N+c0; for(int lc=0;lc<bw;lc++) brd[lr][lc]=row[lc];}
            for(int lr=0;lr<bh;lr++) for(int lc=0;lc<bw;lc++) bwr[lc][lr]=brd[lr][lc];
            for(int lc=0;lc<bw;lc++){float*dr=out+(c0+lc)*N+r0; for(int lr=0;lr<bh;lr++) dr[lr]=bwr[lc][lr];}
        }
}

/* ═══════════════════════════════════════════════════════════════════════ */

static const char* V_NAMES[] = {
    "naive","naive_c2","tiled","tiled_c2",
    "blk_naive","blk_naive_c2","blk_tiled","blk_tiled_c2",
    "blk_aligned","blk_aligned_c2","blk_aligned_mt","blk_aligned_mt_c2",
    "locbuf","locbuf_c2","locbuf_blk","locbuf_blk_c2",
    "locbuf_blk_mt","locbuf_blk_mt_c2","locbuf_2buf","locbuf_2buf_c2",
};
static const int N_VARIANTS = 20;
static bool is_blocked(int v) { return (v>=4 && v<=11) || (v>=14 && v<=17); }

/* ── Verification ─────────────────────────────────────────────────── */

static void ref_transpose(const float* in, float* ref, int N) {
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            ref[c*N+r] = in[r*N+c];
}

static float verify(const float* out, const float* ref, int N, int SB, bool blk) {
    float mx = 0;
    if (blk) {
        int NB = N/SB;
        for (int r=0;r<N;r++) for(int c=0;c<N;c++) {
            float e = fabsf(out[(r/SB*NB+c/SB)*SB*SB+(r%SB)*SB+c%SB] - ref[r*N+c]);
            if (e > mx) mx = e;
        }
    } else {
        for (size_t i=0; i<(size_t)N*N; i++) { float e=fabsf(out[i]-ref[i]); if(e>mx) mx=e; }
    }
    return mx;
}

/* ═══════════════════════════════════════════════════════════════════════
 *  main
 * ═══════════════════════════════════════════════════════════════════════ */

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr,
            "Usage: %s <N> <variant> <csv> [TB=64] [SB=32] [MT=8] "
            "[WARMUP=3] [REPS=20] [THREADS=0] [NUMA=0]\n"
            "  NUMA: 0=none  1=contig(_nd)  2=cyclic(_nc)  3=transpose(_nt)\n", argv[0]);
        return 1;
    }

    int N=atoi(argv[1]), VAR=atoi(argv[2]); const char* csv=argv[3];
    int TB     = argc> 4 ? atoi(argv[ 4]) : 64;
    int SB     = argc> 5 ? atoi(argv[ 5]) : 32;
    int MT     = argc> 6 ? atoi(argv[ 6]) : 8;
    int WARMUP = argc> 7 ? atoi(argv[ 7]) : 3;
    int REPS   = argc> 8 ? atoi(argv[ 8]) : 20;
    int THREADS= argc> 9 ? atoi(argv[ 9]) : 0;
    int NPOL   = argc>10 ? atoi(argv[10]) : 0;

    if (VAR<0||VAR>=N_VARIANTS) { fprintf(stderr,"bad variant\n"); return 1; }
    bool blocked = is_blocked(VAR);
    if (blocked && N%SB) { fprintf(stderr,"N%%SB!=0\n"); return 1; }
    if (VAR>=12 && (blocked?SB:TB)>MAX_TB) { fprintf(stderr,"tile>MAX_TB\n"); return 1; }
    if (NPOL<0||NPOL>=NP_COUNT) { fprintf(stderr,"bad NUMA policy\n"); return 1; }
    if ((NPOL==NP_CYCLIC||NPOL==NP_TRANSPOSE) && !blocked) {
        fprintf(stderr,"NUMA cyclic/transpose only for blocked variants\n"); return 1;
    }
    NumaPolicy npol = (NumaPolicy)NPOL;

    if (THREADS>0) omp_set_num_threads(THREADS);
    int nthreads;
    #pragma omp parallel
    { #pragma omp single nthreads = omp_get_num_threads(); }

    int D = detect_numa_nodes();
    char vname[128];
    snprintf(vname, sizeof(vname), "%s%s", V_NAMES[VAR], NP_SUFFIX[NPOL]);

    size_t elems = (size_t)N*N, bytes = elems*sizeof(float);

    /* ── Allocate & bind ────────────────────────────────────────────── */
    /* h_row, h_ref: temporaries, no NUMA needed. */
    float* h_row = (float*)aligned_alloc(64, bytes);
    float* h_ref = (float*)aligned_alloc(64, bytes);

    /* h_in, h_out: benchmark arrays — NUMA-aware when policy != NONE. */
    float *h_in, *h_out;
    bool use_numa = (npol != NP_NONE);
    if (use_numa) {
        h_in  = numa_alloc<float>(elems);
        h_out = numa_alloc<float>(elems);
        /* Bind pages BEFORE first-touch. */
        apply_numa(h_in,  bytes, N, SB, blocked, false, npol, D);
        apply_numa(h_out, bytes, N, SB, blocked, true,  npol, D);
    } else {
        h_in  = (float*)aligned_alloc(64, bytes);
        h_out = (float*)aligned_alloc(64, bytes);
    }

    /* ── Init source (parallel first-touch) ─────────────────────────── */
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < elems; i++) h_row[i] = (float)i / (float)N;

    ref_transpose(h_row, h_ref, N);

    /* ── Prepare h_in (parallel — correct first-touch) ──────────────── */
    if (blocked) {
        int NB = N/SB;
        #pragma omp parallel for schedule(static)
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                h_in[(r/SB*NB+c/SB)*SB*SB + (r%SB)*SB + c%SB] = h_row[r*N+c];
    } else {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < elems; i++) h_in[i] = h_row[i];
    }

    /* Parallel zero of output — faults pages on bound NUMA nodes. */
    parallel_zero(h_out, elems);

    printf("  %s N=%d TB=%d SB=%d MT=%d thr=%d numa=%d(%s) nodes=%d\n",
           vname, N, TB, SB, MT, nthreads, NPOL, NP_SUFFIX[NPOL]+(*NP_SUFFIX[NPOL]?0:0), D);

    /* ── Dispatch ───────────────────────────────────────────────────── */
    auto launch = [&]() {
        switch (VAR) {
            case  0: tr_naive(h_in,h_out,N); break;
            case  1: tr_naive_c2(h_in,h_out,N); break;
            case  2: tr_tiled(h_in,h_out,N,TB); break;
            case  3: tr_tiled_c2(h_in,h_out,N,TB); break;
            case  4: tr_blk_naive(h_in,h_out,N,SB); break;
            case  5: tr_blk_naive_c2(h_in,h_out,N,SB); break;
            case  6: tr_blk_tiled(h_in,h_out,N,SB,TB); break;
            case  7: tr_blk_tiled_c2(h_in,h_out,N,SB,TB); break;
            case  8: tr_blk_aligned(h_in,h_out,N,SB); break;
            case  9: tr_blk_aligned_c2(h_in,h_out,N,SB); break;
            case 10: tr_blk_aligned_mt(h_in,h_out,N,SB,MT); break;
            case 11: tr_blk_aligned_mt_c2(h_in,h_out,N,SB,MT); break;
            case 12: tr_locbuf(h_in,h_out,N,TB); break;
            case 13: tr_locbuf_c2(h_in,h_out,N,TB); break;
            case 14: tr_locbuf_blk(h_in,h_out,N,SB); break;
            case 15: tr_locbuf_blk_c2(h_in,h_out,N,SB); break;
            case 16: tr_locbuf_blk_mt(h_in,h_out,N,SB,MT); break;
            case 17: tr_locbuf_blk_mt_c2(h_in,h_out,N,SB,MT); break;
            case 18: tr_locbuf_2buf(h_in,h_out,N,TB); break;
            case 19: tr_locbuf_2buf_c2(h_in,h_out,N,TB); break;
        }
    };

    /* ── Warmup ─────────────────────────────────────────────────────── */
    for (int i = 0; i < WARMUP; i++) launch();

    /* ── Verify (parallel zero, not serial memset) ──────────────────── */
    parallel_zero(h_out, elems);
    launch();
    float maxerr = verify(h_out, h_ref, N, SB, blocked);
    bool pass = (maxerr == 0.0f);

    /* ── Timed runs ─────────────────────────────────────────────────── */
    double* times = (double*)malloc(REPS * sizeof(double));
    for (int i = 0; i < REPS; i++) {
        double t0 = omp_get_wtime();
        launch();
        double t1 = omp_get_wtime();
        times[i] = t1 - t0;
    }

    std::sort(times, times+REPS);
    double bpi = 2.0*N*(double)N*sizeof(float);
    double med_s = times[REPS/2];
    double cksum = 0;
    for (size_t i = 0; i < elems; i++) cksum += h_out[i];

    printf("%s N=%d TB=%d SB=%d MT=%d thr=%d | med %.4f ms (%.1f GB/s)  "
           "p5 %.4f ms  p95 %.4f ms  maxerr=%.1e  %s  cksum=%.6e\n",
           vname, N, TB, SB, MT, nthreads,
           med_s*1e3, bpi/med_s/1e9,
           times[(int)(REPS*0.05)]*1e3, times[(int)(REPS*0.95)]*1e3,
           maxerr, pass?"PASS":"FAIL", cksum);

    /* ── CSV ────────────────────────────────────────────────────────── */
    FILE* f = fopen(csv, "a");
    if (f) {
        for (int i = 0; i < REPS; i++)
            fprintf(f, "%s,%d,%d,%d,%d,%d,%d,%.9f,%.3f,%.6e,%s\n",
                    vname, N, TB, SB, MT, nthreads, i,
                    times[i], bpi/times[i]/1e9, cksum, pass?"PASS":"FAIL");
        fclose(f);
    }

    free(times); free(h_row); free(h_ref);
    if (use_numa) { numa_dealloc(h_in,elems); numa_dealloc(h_out,elems); }
    else { free(h_in); free(h_out); }
    return pass ? 0 : 1;
}