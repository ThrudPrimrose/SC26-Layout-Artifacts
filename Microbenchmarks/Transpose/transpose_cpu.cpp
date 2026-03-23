#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <omp.h>

// ── Blocked-layout index helpers ──
static inline int blk_idx(int r, int c, int N, int SB) {
    int NB = N / SB;
    return (r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB;
}

// ══════════════════════════════════════════════════════════════════════
//  Row-major variants
// ══════════════════════════════════════════════════════════════════════

// ── V0: Naive, outer-loop parallel ──
static void tr_naive(const float* __restrict__ in, float* __restrict__ out, int N) {
    #pragma omp parallel for schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            out[c * N + r] = in[r * N + c];
}

// ── V1: Naive, collapse(2) ──
static void tr_naive_c2(const float* __restrict__ in, float* __restrict__ out, int N) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            out[c * N + r] = in[r * N + c];
}

// ── V2: Tiled, outer-loop parallel ──
static void tr_tiled(const float* __restrict__ in, float* __restrict__ out, int N, int TB) {
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
            int rend = std::min(r0 + TB, N), cend = std::min(c0 + TB, N);
            for (int r = r0; r < rend; r++)
                for (int c = c0; c < cend; c++)
                    out[c * N + r] = in[r * N + c];
        }
    }
}

// ── V3: Tiled, collapse(2) ──
static void tr_tiled_c2(const float* __restrict__ in, float* __restrict__ out, int N, int TB) {
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
            int rend = std::min(r0 + TB, N), cend = std::min(c0 + TB, N);
            for (int r = r0; r < rend; r++)
                for (int c = c0; c < cend; c++)
                    out[c * N + r] = in[r * N + c];
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
//  Blocked-layout variants (SB-blocked storage in both src and dst)
// ══════════════════════════════════════════════════════════════════════

// ── V4: Blocked naive, outer-loop parallel ──
static void tr_blk_naive(const float* __restrict__ in, float* __restrict__ out, int N, int SB) {
    int NB = N / SB;
    #pragma omp parallel for schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            int si = (r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB;
            int di = (c/SB * NB + r/SB) * SB*SB + (c%SB)*SB + r%SB;
            out[di] = in[si];
        }
}

// ── V5: Blocked naive, collapse(2) ──
static void tr_blk_naive_c2(const float* __restrict__ in, float* __restrict__ out, int N, int SB) {
    int NB = N / SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            int si = (r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB;
            int di = (c/SB * NB + r/SB) * SB*SB + (c%SB)*SB + r%SB;
            out[di] = in[si];
        }
}

// ── V6: Blocked tiled, outer-loop parallel ──
//    Outer loops tile over storage blocks; inner loops over elements within tile.
static void tr_blk_tiled(const float* __restrict__ in, float* __restrict__ out,
                          int N, int SB, int TB) {
    int NB = N / SB;
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
            int rend = std::min(r0 + TB, N), cend = std::min(c0 + TB, N);
            for (int r = r0; r < rend; r++)
                for (int c = c0; c < cend; c++) {
                    int si = (r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB;
                    int di = (c/SB * NB + r/SB) * SB*SB + (c%SB)*SB + r%SB;
                    out[di] = in[si];
                }
        }
    }
}

// ── V7: Blocked tiled, collapse(2) ──
static void tr_blk_tiled_c2(const float* __restrict__ in, float* __restrict__ out,
                              int N, int SB, int TB) {
    int NB = N / SB;
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
            int rend = std::min(r0 + TB, N), cend = std::min(c0 + TB, N);
            for (int r = r0; r < rend; r++)
                for (int c = c0; c < cend; c++) {
                    int si = (r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB;
                    int di = (c/SB * NB + r/SB) * SB*SB + (c%SB)*SB + r%SB;
                    out[di] = in[si];
                }
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
//  Block-aware tiled variants: tile size = SB so inner loops walk
//  contiguously within a single storage block.
// ══════════════════════════════════════════════════════════════════════

// ── V8: Block-aligned tiled, outer-loop parallel ──
//    Iterate over (block_row, block_col) pairs; inner loop walks the
//    SB x SB sub-block which is contiguous in memory.
static void tr_blk_aligned(const float* __restrict__ in, float* __restrict__ out,
                            int N, int SB) {
    int NB = N / SB;
    #pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            // src is row-major SB x SB block at (br, bc)
            // dst is row-major SB x SB block at (bc, br)
            // Transpose the sub-block:
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    dst[lc * SB + lr] = src[lr * SB + lc];
        }
    }
}

// ── V9: Block-aligned tiled, collapse(2) ──
static void tr_blk_aligned_c2(const float* __restrict__ in, float* __restrict__ out,
                                int N, int SB) {
    int NB = N / SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    dst[lc * SB + lr] = src[lr * SB + lc];
        }
    }
}

// ── V10: Block-aligned + micro-tiled inner loop, outer-loop parallel ──
//    The SB x SB sub-block transpose is itself cache-tiled with micro-tile MT.
static void tr_blk_aligned_mt(const float* __restrict__ in, float* __restrict__ out,
                                int N, int SB, int MT) {
    int NB = N / SB;
    #pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            for (int lr0 = 0; lr0 < SB; lr0 += MT)
                for (int lc0 = 0; lc0 < SB; lc0 += MT) {
                    int lrend = std::min(lr0 + MT, SB);
                    int lcend = std::min(lc0 + MT, SB);
                    for (int lr = lr0; lr < lrend; lr++)
                        for (int lc = lc0; lc < lcend; lc++)
                            dst[lc * SB + lr] = src[lr * SB + lc];
                }
        }
    }
}

// ── V11: Block-aligned + micro-tiled inner loop, collapse(2) ──
static void tr_blk_aligned_mt_c2(const float* __restrict__ in, float* __restrict__ out,
                                   int N, int SB, int MT) {
    int NB = N / SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            for (int lr0 = 0; lr0 < SB; lr0 += MT)
                for (int lc0 = 0; lc0 < SB; lc0 += MT) {
                    int lrend = std::min(lr0 + MT, SB);
                    int lcend = std::min(lc0 + MT, SB);
                    for (int lr = lr0; lr < lrend; lr++)
                        for (int lc = lc0; lc < lcend; lc++)
                            dst[lc * SB + lr] = src[lr * SB + lc];
                }
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
//  Local-buffer variants (CPU analogue of GPU shared memory)
//
//  Each thread owns a stack-allocated TB x TB (or SB x SB) buffer that
//  sits in L1 / registers.  Two-phase approach:
//    Phase 1 (LOAD):  sequential reads from global -> fill buf row-major
//    Phase 2 (STORE): sequential writes to global  <- read buf column-major
//  The strided access only hits the tiny local buffer, never DRAM.
// ══════════════════════════════════════════════════════════════════════

static constexpr int MAX_TB = 128;   // max tile side; 128x128 = 64 KB per thread

// ── V12: Local-buf, row-major, outer-loop parallel ──
//   Load:  buf[lr][lc] = in[(r0+lr)*N + c0+lc]        (sequential read per row)
//   Store: out[(c0+lc)*N + r0+lr] = buf[lr][lc]        (sequential write per row)
//          iterate lc outer, lr inner => writes stride-1
static void tr_locbuf(const float* __restrict__ in, float* __restrict__ out,
                       int N, int TB) {
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        float buf[MAX_TB][MAX_TB];
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
            int rend = std::min(r0 + TB, N), cend = std::min(c0 + TB, N);
            int bh = rend - r0, bw = cend - c0;
            // Phase 1: sequential reads from in (row-major source)
            for (int lr = 0; lr < bh; lr++) {
                const float* row = in + (r0 + lr) * N + c0;
                for (int lc = 0; lc < bw; lc++)
                    buf[lr][lc] = row[lc];
            }
            // Phase 2: sequential writes to out
            // out[(c0+lc)*N + r0+lr] — vary lr in inner loop => stride-1 writes
            for (int lc = 0; lc < bw; lc++) {
                float* dst_row = out + (c0 + lc) * N + r0;
                for (int lr = 0; lr < bh; lr++)
                    dst_row[lr] = buf[lr][lc];   // column read of buf (in L1)
            }
        }
    }
}

// ── V13: Local-buf, row-major, collapse(2) ──
static void tr_locbuf_c2(const float* __restrict__ in, float* __restrict__ out,
                           int N, int TB) {
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        for (int tc = 0; tc < NT; tc++) {
            float buf[MAX_TB][MAX_TB];
            int r0 = tr * TB, c0 = tc * TB;
            int rend = std::min(r0 + TB, N), cend = std::min(c0 + TB, N);
            int bh = rend - r0, bw = cend - c0;
            for (int lr = 0; lr < bh; lr++) {
                const float* row = in + (r0 + lr) * N + c0;
                for (int lc = 0; lc < bw; lc++)
                    buf[lr][lc] = row[lc];
            }
            for (int lc = 0; lc < bw; lc++) {
                float* dst_row = out + (c0 + lc) * N + r0;
                for (int lr = 0; lr < bh; lr++)
                    dst_row[lr] = buf[lr][lc];
            }
        }
    }
}

// ── V14: Local-buf, blocked layout, outer-loop parallel ──
//   Source block (br,bc) is contiguous SB*SB floats => sequential load.
//   Dest block (bc,br) is contiguous SB*SB floats => sequential store.
//   Transpose happens purely inside the local buffer.
static void tr_locbuf_blk(const float* __restrict__ in, float* __restrict__ out,
                            int N, int SB) {
    int NB = N / SB;
    #pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++) {
        float buf[MAX_TB][MAX_TB];
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            // Phase 1: load contiguous source block into buf (sequential read)
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    buf[lr][lc] = src[lr * SB + lc];
            // Phase 2: write transposed into contiguous dest block
            // dst[lc*SB+lr] — vary lr in inner loop => stride-1 writes
            for (int lc = 0; lc < SB; lc++) {
                float* dst_row = dst + lc * SB;
                for (int lr = 0; lr < SB; lr++)
                    dst_row[lr] = buf[lr][lc];   // column read of buf (in L1)
            }
        }
    }
}

// ── V15: Local-buf, blocked layout, collapse(2) ──
static void tr_locbuf_blk_c2(const float* __restrict__ in, float* __restrict__ out,
                               int N, int SB) {
    int NB = N / SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            float buf[MAX_TB][MAX_TB];
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    buf[lr][lc] = src[lr * SB + lc];
            for (int lc = 0; lc < SB; lc++) {
                float* dst_row = dst + lc * SB;
                for (int lr = 0; lr < SB; lr++)
                    dst_row[lr] = buf[lr][lc];
            }
        }
    }
}

// ── V16: Local-buf, blocked layout + micro-tiled buffer access, outer-loop parallel ──
//   The SB x SB buffer transpose is itself tiled with MT to stay in registers.
//   Load:  sequential read of whole SB*SB source block into buf
//   Store: micro-tiled write — MT x MT sub-tiles keep the write pointer hot
static void tr_locbuf_blk_mt(const float* __restrict__ in, float* __restrict__ out,
                               int N, int SB, int MT) {
    int NB = N / SB;
    #pragma omp parallel for schedule(static)
    for (int br = 0; br < NB; br++) {
        float buf[MAX_TB][MAX_TB];
        for (int bc = 0; bc < NB; bc++) {
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            // Phase 1: sequential load
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    buf[lr][lc] = src[lr * SB + lc];
            // Phase 2: micro-tiled transposed write
            for (int lc0 = 0; lc0 < SB; lc0 += MT) {
                int lcend = std::min(lc0 + MT, SB);
                for (int lr0 = 0; lr0 < SB; lr0 += MT) {
                    int lrend = std::min(lr0 + MT, SB);
                    for (int lc = lc0; lc < lcend; lc++) {
                        float* dst_row = dst + lc * SB + lr0;
                        for (int lr = lr0; lr < lrend; lr++)
                            dst_row[lr - lr0] = buf[lr][lc];
                    }
                }
            }
        }
    }
}

// ── V17: Local-buf, blocked layout + micro-tiled, collapse(2) ──
static void tr_locbuf_blk_mt_c2(const float* __restrict__ in, float* __restrict__ out,
                                  int N, int SB, int MT) {
    int NB = N / SB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int br = 0; br < NB; br++) {
        for (int bc = 0; bc < NB; bc++) {
            float buf[MAX_TB][MAX_TB];
            const float* src = in  + (br * NB + bc) * SB * SB;
            float*       dst = out + (bc * NB + br) * SB * SB;
            for (int lr = 0; lr < SB; lr++)
                for (int lc = 0; lc < SB; lc++)
                    buf[lr][lc] = src[lr * SB + lc];
            for (int lc0 = 0; lc0 < SB; lc0 += MT) {
                int lcend = std::min(lc0 + MT, SB);
                for (int lr0 = 0; lr0 < SB; lr0 += MT) {
                    int lrend = std::min(lr0 + MT, SB);
                    for (int lc = lc0; lc < lcend; lc++) {
                        float* dst_row = dst + lc * SB + lr0;
                        for (int lr = lr0; lr < lrend; lr++)
                            dst_row[lr - lr0] = buf[lr][lc];
                    }
                }
            }
        }
    }
}

// ── V18: Local-buf, row-major, read-optimized then write-optimized ──
//   Two buffers: load into buf_rd with sequential reads, then copy-transpose
//   into buf_wr so the write phase is a simple sequential memcpy-style loop.
//   Cost: 2x buffer memory, but both DRAM read and write are perfectly sequential.
static void tr_locbuf_2buf(const float* __restrict__ in, float* __restrict__ out,
                            int N, int TB) {
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        float buf_rd[MAX_TB][MAX_TB];   // source tile  (row = source row)
        float buf_wr[MAX_TB][MAX_TB];   // transposed   (row = source col = dest row)
        for (int tc = 0; tc < NT; tc++) {
            int r0 = tr * TB, c0 = tc * TB;
            int rend = std::min(r0 + TB, N), cend = std::min(c0 + TB, N);
            int bh = rend - r0, bw = cend - c0;
            // Phase 1: sequential read from in
            for (int lr = 0; lr < bh; lr++) {
                const float* row = in + (r0 + lr) * N + c0;
                for (int lc = 0; lc < bw; lc++)
                    buf_rd[lr][lc] = row[lc];
            }
            // Phase 2: local transpose (all in L1)
            for (int lr = 0; lr < bh; lr++)
                for (int lc = 0; lc < bw; lc++)
                    buf_wr[lc][lr] = buf_rd[lr][lc];
            // Phase 3: sequential write from buf_wr
            for (int lc = 0; lc < bw; lc++) {
                float* dst_row = out + (c0 + lc) * N + r0;
                for (int lr = 0; lr < bh; lr++)
                    dst_row[lr] = buf_wr[lc][lr];   // row read of buf_wr (sequential)
            }
        }
    }
}

// ── V19: Local-buf, row-major, 2-buffer, collapse(2) ──
static void tr_locbuf_2buf_c2(const float* __restrict__ in, float* __restrict__ out,
                                int N, int TB) {
    int NT = (N + TB - 1) / TB;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int tr = 0; tr < NT; tr++) {
        for (int tc = 0; tc < NT; tc++) {
            float buf_rd[MAX_TB][MAX_TB];
            float buf_wr[MAX_TB][MAX_TB];
            int r0 = tr * TB, c0 = tc * TB;
            int rend = std::min(r0 + TB, N), cend = std::min(c0 + TB, N);
            int bh = rend - r0, bw = cend - c0;
            for (int lr = 0; lr < bh; lr++) {
                const float* row = in + (r0 + lr) * N + c0;
                for (int lc = 0; lc < bw; lc++)
                    buf_rd[lr][lc] = row[lc];
            }
            for (int lr = 0; lr < bh; lr++)
                for (int lc = 0; lc < bw; lc++)
                    buf_wr[lc][lr] = buf_rd[lr][lc];
            for (int lc = 0; lc < bw; lc++) {
                float* dst_row = out + (c0 + lc) * N + r0;
                for (int lr = 0; lr < bh; lr++)
                    dst_row[lr] = buf_wr[lc][lr];
            }
        }
    }
}

// ══════════════════════════════════════════════════════════════════════

static const char* V_NAMES[] = {
    "naive",              // 0
    "naive_c2",           // 1
    "tiled",              // 2
    "tiled_c2",           // 3
    "blk_naive",          // 4
    "blk_naive_c2",       // 5
    "blk_tiled",          // 6
    "blk_tiled_c2",       // 7
    "blk_aligned",        // 8
    "blk_aligned_c2",     // 9
    "blk_aligned_mt",     // 10
    "blk_aligned_mt_c2",  // 11
    "locbuf",             // 12
    "locbuf_c2",          // 13
    "locbuf_blk",         // 14
    "locbuf_blk_c2",      // 15
    "locbuf_blk_mt",      // 16
    "locbuf_blk_mt_c2",   // 17
    "locbuf_2buf",        // 18
    "locbuf_2buf_c2",     // 19
};
static const int N_VARIANTS = 20;

static bool is_blocked(int var) { return (var >= 4 && var <= 11) || (var >= 14 && var <= 17); }

// ══════════════════════════════════════════════════════════════════════
//  Verification
// ══════════════════════════════════════════════════════════════════════

// Compute reference transpose into ref[] (row-major -> row-major).
static void ref_transpose(const float* in, float* ref, int N) {
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            ref[c * N + r] = in[r * N + c];
}

// Convert blocked -> row-major for comparison.
static void blocked_to_row(const float* blk, float* row, int N, int SB) {
    int NB = N / SB;
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++)
            row[r * N + c] = blk[(r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB];
}

// Returns max absolute error; -1 on dimension mismatch.
static float verify(const float* out, const float* ref, int N, int SB, bool blocked_out) {
    float maxerr = 0.0f;
    if (blocked_out) {
        // out is in blocked layout; convert and compare
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
            "Usage: %s <N> <variant> <csv> [TB=64] [SB=32] [MT=8] [WARMUP=3] [REPS=20] [THREADS=0]\n"
            "\n"
            "  variant:\n"
            "    Row-major:\n"
            "      0 = naive              (omp parallel for)\n"
            "      1 = naive_c2           (collapse(2))\n"
            "      2 = tiled              (omp parallel for)\n"
            "      3 = tiled_c2           (collapse(2))\n"
            "    Blocked layout (requires N %% SB == 0):\n"
            "      4 = blk_naive          (omp parallel for)\n"
            "      5 = blk_naive_c2       (collapse(2))\n"
            "      6 = blk_tiled          (omp parallel for, tile TB)\n"
            "      7 = blk_tiled_c2       (collapse(2), tile TB)\n"
            "      8 = blk_aligned        (omp parallel for, tile=SB)\n"
            "      9 = blk_aligned_c2     (collapse(2), tile=SB)\n"
            "     10 = blk_aligned_mt     (omp parallel for, micro-tile MT)\n"
            "     11 = blk_aligned_mt_c2  (collapse(2), micro-tile MT)\n"
            "    Local-buffer (CPU 'shared memory'):\n"
            "     12 = locbuf             (row-major, omp parallel for)\n"
            "     13 = locbuf_c2          (row-major, collapse(2))\n"
            "     14 = locbuf_blk         (blocked, omp parallel for)\n"
            "     15 = locbuf_blk_c2      (blocked, collapse(2))\n"
            "     16 = locbuf_blk_mt      (blocked + micro-tile, omp parallel for)\n"
            "     17 = locbuf_blk_mt_c2   (blocked + micro-tile, collapse(2))\n"
            "     18 = locbuf_2buf        (row-major, 2-buf, omp parallel for)\n"
            "     19 = locbuf_2buf_c2     (row-major, 2-buf, collapse(2))\n"
            "\n"
            "  TB:      cache tile size for row-major tiled variants (default 64)\n"
            "  SB:      storage block size for blocked variants (default 32)\n"
            "  MT:      micro-tile size for blk_aligned_mt variants (default 8)\n"
            "  THREADS: 0 = use OMP_NUM_THREADS (default)\n",
            argv[0]);
        return 1;
    }

    int N   = atoi(argv[1]);
    int VAR = atoi(argv[2]);
    const char* csv = argv[3];
    int TB      = (argc > 4) ? atoi(argv[4]) : 64;
    int SB      = (argc > 5) ? atoi(argv[5]) : 32;
    int MT      = (argc > 6) ? atoi(argv[6]) : 8;
    int WARMUP  = (argc > 7) ? atoi(argv[7]) : 3;
    int REPS    = (argc > 8) ? atoi(argv[8]) : 20;
    int THREADS = (argc > 9) ? atoi(argv[9]) : 0;

    if (VAR < 0 || VAR >= N_VARIANTS) {
        fprintf(stderr, "Unknown variant %d\n", VAR);
        return 1;
    }
    if (is_blocked(VAR) && N % SB != 0) {
        fprintf(stderr, "N=%d not divisible by SB=%d for blocked variant\n", N, SB);
        return 1;
    }
    if (VAR >= 12) {  // locbuf variants
        int tile = is_blocked(VAR) ? SB : TB;
        if (tile > MAX_TB) {
            fprintf(stderr, "Tile size %d exceeds MAX_TB=%d for locbuf variant\n", tile, MAX_TB);
            return 1;
        }
    }
    if (THREADS > 0) omp_set_num_threads(THREADS);
    int nthreads;
    #pragma omp parallel
    {
        #pragma omp single
        nthreads = omp_get_num_threads();
    }

    size_t elems = (size_t)N * N;
    size_t bytes = elems * sizeof(float);

    float* h_row = (float*)aligned_alloc(64, bytes);  // row-major source
    float* h_ref = (float*)aligned_alloc(64, bytes);  // reference transpose (row-major)
    float* h_in  = (float*)aligned_alloc(64, bytes);  // actual input  (may be blocked)
    float* h_out = (float*)aligned_alloc(64, bytes);  // actual output

    // Fill row-major source
    for (size_t i = 0; i < elems; i++) h_row[i] = (float)i / (float)N;

    // Reference transpose (always row-major)
    ref_transpose(h_row, h_ref, N);

    // Prepare input
    if (is_blocked(VAR)) {
        int NB = N / SB;
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                h_in[(r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB] = h_row[r * N + c];
    } else {
        memcpy(h_in, h_row, bytes);
    }

    // ── Dispatch lambda ──
    auto launch = [&]() {
        switch (VAR) {
            case  0: tr_naive      (h_in, h_out, N);         break;
            case  1: tr_naive_c2   (h_in, h_out, N);         break;
            case  2: tr_tiled      (h_in, h_out, N, TB);     break;
            case  3: tr_tiled_c2   (h_in, h_out, N, TB);     break;
            case  4: tr_blk_naive     (h_in, h_out, N, SB);        break;
            case  5: tr_blk_naive_c2  (h_in, h_out, N, SB);        break;
            case  6: tr_blk_tiled     (h_in, h_out, N, SB, TB);    break;
            case  7: tr_blk_tiled_c2  (h_in, h_out, N, SB, TB);    break;
            case  8: tr_blk_aligned   (h_in, h_out, N, SB);        break;
            case  9: tr_blk_aligned_c2(h_in, h_out, N, SB);        break;
            case 10: tr_blk_aligned_mt   (h_in, h_out, N, SB, MT); break;
            case 11: tr_blk_aligned_mt_c2(h_in, h_out, N, SB, MT); break;
            case 12: tr_locbuf           (h_in, h_out, N, TB);     break;
            case 13: tr_locbuf_c2        (h_in, h_out, N, TB);     break;
            case 14: tr_locbuf_blk       (h_in, h_out, N, SB);    break;
            case 15: tr_locbuf_blk_c2    (h_in, h_out, N, SB);    break;
            case 16: tr_locbuf_blk_mt    (h_in, h_out, N, SB, MT); break;
            case 17: tr_locbuf_blk_mt_c2 (h_in, h_out, N, SB, MT); break;
            case 18: tr_locbuf_2buf      (h_in, h_out, N, TB);     break;
            case 19: tr_locbuf_2buf_c2   (h_in, h_out, N, TB);     break;
        }
    };

    // ── Warmup ──
    for (int i = 0; i < WARMUP; i++) launch();

    // ── Verify ──
    memset(h_out, 0, bytes);
    launch();
    float maxerr = verify(h_out, h_ref, N, SB, is_blocked(VAR));
    bool pass = (maxerr == 0.0f);

    // ── Timed runs ──
    double* times = (double*)malloc(REPS * sizeof(double));
    for (int i = 0; i < REPS; i++) {
        double t0 = omp_get_wtime();
        launch();
        double t1 = omp_get_wtime();
        times[i] = t1 - t0;
    }

    // Stats
    double total = 0;
    for (int i = 0; i < REPS; i++) total += times[i];
    double mean_s = total / REPS;

    // Sort for median / percentiles
    std::sort(times, times + REPS);
    double med_s = times[REPS / 2];
    double p5_s  = times[(int)(REPS * 0.05)];
    double p95_s = times[(int)(REPS * 0.95)];

    double bpi   = 2.0 * N * (double)N * sizeof(float);
    double mean_gbps = bpi / mean_s / 1e9;
    double med_gbps  = bpi / med_s  / 1e9;

    double cksum = 0;
    for (size_t i = 0; i < elems; i++) cksum += h_out[i];

    printf("%s N=%d TB=%d SB=%d MT=%d threads=%d | "
           "mean %.4f ms (%.1f GB/s)  med %.4f ms (%.1f GB/s)  "
           "p5 %.4f ms  p95 %.4f ms  maxerr=%.1e  %s  cksum=%.6e\n",
           V_NAMES[VAR], N, TB, SB, MT, nthreads,
           mean_s*1e3, mean_gbps, med_s*1e3, med_gbps,
           p5_s*1e3, p95_s*1e3, maxerr, pass ? "PASS" : "FAIL", cksum);

    // ── CSV ──
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