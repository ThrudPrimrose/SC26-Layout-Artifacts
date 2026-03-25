#pragma once

#include <cstring>
#include <cstdio>

#include "amx_tile_context.hpp"
// #include "svmopa.hpp"

#ifndef __forceinline
#define __forceinline inline __attribute__((always_inline))
#endif

namespace core_ir {

//============================================================================
// TYPE CONVERSION HELPERS - All type conversion intrinsics in one place
//============================================================================
namespace FragConvert {}  // namespace FragConvert

// ============================================================================
// TensorFragment - AMX tile register abstraction
// ============================================================================
template <Tensor T, typename IN_T, typename OUT_T, int TILE_M, int TILE_N,
          typename AMX = DefaultAMXContext>
class TensorFragment {
  // AMX uses 8 tile registers:
  // Tiles 0-1: A fragments (16 rows x 32 cols each for FP16)
  // Tiles 2-3: B fragments (16 rows x 32 cols each for FP16, interleaved)
  // Tiles 4-7: D/C fragments (16 rows x 16 cols each for FP32)
  static constexpr int TILE_M_2 = TILE_M / 2;
  static constexpr int TILE_N_2 = TILE_N / 2;

 public:
  template <bool INIT_ZERO = true>
  __forceinline TensorFragment() {
    if constexpr (INIT_ZERO && T == Tensor::D) {
      // Zero the accumulator tiles
      _tile_zero(4);
      _tile_zero(5);
      _tile_zero(6);
      _tile_zero(7);
    }
  }

  // Load D fragment from memory (for accumulation)
  template <MajorAxis OUT_AXIS = MajorAxis::ROW>
  __forceinline void load(const void* __restrict__ pmem, int stride = 0) {
    if constexpr (T == Tensor::D) {
      static_assert(OUT_AXIS == MajorAxis::ROW,
                    "Only ROW major output is supported for now");
      const OUT_T* __restrict__ mem = reinterpret_cast<const OUT_T*>(pmem);
      _tile_loadd(4, mem, stride * sizeof(OUT_T));
      _tile_loadd(5, mem + TILE_N_2, stride * sizeof(OUT_T));
      _tile_loadd(6, mem + TILE_M_2 * stride, stride * sizeof(OUT_T));
      _tile_loadd(7, mem + TILE_M_2 * stride + TILE_N_2,
                  stride * sizeof(OUT_T));
    } else if constexpr (T == Tensor::A) {
      const IN_T* mem = reinterpret_cast<const IN_T*>(pmem);
      _tile_loadd(0, mem, TILE_N * sizeof(IN_T));
      _tile_loadd(1, mem + (TILE_M / 2) * TILE_N, TILE_N * sizeof(IN_T));
    } else if constexpr (T == Tensor::B) {
      const IN_T* mem = reinterpret_cast<const IN_T*>(pmem);
      _tile_loadd(2, mem, 2 * TILE_N * sizeof(IN_T));
      _tile_loadd(3, mem + TILE_N, 2 * TILE_N * sizeof(IN_T));
    }
  }

  // Unload D fragment to memory
  template <MajorAxis OUT_AXIS>
  __forceinline void unload_unpack(OUT_T* mem, int stride) {
    static_assert(T == Tensor::D,
                  "Fragment unload only supported for D tensor");
    static_assert(OUT_AXIS == MajorAxis::ROW,
                  "Only ROW major output is supported for now");

    _tile_stored(4, mem, stride * sizeof(OUT_T));
    _tile_stored(5, mem + TILE_N_2, stride * sizeof(OUT_T));
    _tile_stored(6, mem + TILE_M_2 * stride, stride * sizeof(OUT_T));
    _tile_stored(7, mem + TILE_M_2 * stride + TILE_N_2, stride * sizeof(OUT_T));
  }
};

// ============================================================================
// PUBLIC API: TENSOR LAYOUT - Unified tensor abstraction
//
// Template Parameters:
//   T       - Element type (float, half, etc.)
//   AMXTC   - AMXTileContext<>
//   TENSOR  - Tensor role (A, B, or D)
//   TILE_X  - Tile rows (M for A, K for B)
//   TILE_Y  - Tile cols (K for A, N for B)
// ============================================================================

template <typename T, typename AMXTC, Tensor TENSOR, int TILE_MN, int TILE_K>
class TensorLayout {
  // Private configuration from AMXTC
  static constexpr MajorAxis AXIS = AMXTC::major_axis;
  static constexpr bool kmajor =
      (TENSOR == Tensor::A && AXIS == MajorAxis::ROW) ||
      (TENSOR == Tensor::B && AXIS == MajorAxis::COLUMN);

  const int _dim_mn, _dim_k;
  const T* __restrict__ _ptr;

 public:
  TensorLayout(const T* __restrict__ ptr, const int dim_mn, const int dim_k)
      : _ptr(ptr), _dim_mn(dim_mn), _dim_k(dim_k) {}

  __forceinline void print(__m512i& reg) {
    alignas(64) T vals[32];
    _mm512_storeu_si512(vals, reg);
    for (int i = 0; i < 32; i++) {
      printf("%.4f ", (float)vals[i]);
    }
    printf("\n");
  }

  __forceinline int pack_bytes() {
    return (((_dim_mn + TILE_MN - 1) / TILE_MN) * TILE_MN) *
           (((_dim_k + TILE_K - 1) / TILE_K) * TILE_K) * sizeof(T);
  }

  __forceinline int pack_elements() {
    return (((_dim_mn + TILE_MN - 1) / TILE_MN) * TILE_MN) *
           (((_dim_k + TILE_K - 1) / TILE_K) * TILE_K);
  }

  // Pack tile (tmn, tk) into contiguous buffer smem
  __forceinline void pack(int tmn, int tk, T* __restrict__& smem) {
    if constexpr (std::is_same_v<T, _Float16> ||
                  std::is_same_v<T, __bfloat16>) {
      if constexpr (TENSOR == Tensor::A) {
        if constexpr (AXIS == MajorAxis::ROW) {
          // A is row-major: tile at (tmn * TILE_MN, tk * TILE_K)
          const T* __restrict__ tile_ptr = _ptr + tmn * TILE_MN * _dim_k + tk * TILE_K;
#pragma unroll
          for (int i = 0; i < TILE_MN; i++) {
#pragma unroll
            for (int j = 0; j < TILE_K; j += 64 / sizeof(T)) {
              __m512i reg =
                  _mm512_loadu_si512((__m512i*)&tile_ptr[i * _dim_k + j]);
              _mm512_storeu_si512((__m512i*)&smem[i * TILE_K + j], reg);
            }
          }
        } else {
          // transpose A from column-major to row-major
          const T* __restrict__ tile_ptr = _ptr + tmn * TILE_MN + tk * TILE_K * _dim_mn;
          __mmask32 k1 = 0x0000FFFF;
          constexpr int REGS = 16;
#pragma unroll
          for (int i = 0; i < TILE_MN; i += 16) {
#pragma unroll
            for (int j = 0; j < TILE_K; j += 16) {
              __m512i rows[REGS];
#pragma unroll
              for (int r = 0; r < REGS; r++) {
                rows[r] = _mm512_maskz_loadu_epi16(
                    k1, &tile_ptr[(j + r) * _dim_mn + i]);
              }
              // Stage 1: epi16 interleave
#pragma unroll
              __m512i transposed[REGS];
              for (int r = 0; r < REGS; r += 2) {
                transposed[r] = _mm512_unpacklo_epi16(rows[r], rows[r + 1]);
                transposed[r + 1] = _mm512_unpackhi_epi16(rows[r], rows[r + 1]);
              }
              // Stage 2: epi32 interleave
#pragma unroll
              for (int r = 0; r < REGS; r += 4) {
                for (int l = 0; l < 2; ++l) {
                  rows[r + l] = _mm512_unpacklo_epi32(transposed[r + l],
                                                      transposed[r + l + 2]);
                  rows[r + l + 2] = _mm512_unpackhi_epi32(
                      transposed[r + l], transposed[r + l + 2]);
                }
              }
              // Stage 3: epi64 interleave
#pragma unroll
              for (int r = 0; r < REGS; r += 8) {
                for (int l = 0; l < 4; ++l) {
                  transposed[r + l] =
                      _mm512_unpacklo_epi64(rows[r + l], rows[r + l + 4]);
                  transposed[r + l + 4] =
                      _mm512_unpackhi_epi64(rows[r + l], rows[r + l + 4]);
                }
              }
              // Step 4: Shuffle
#pragma unroll
              for (int r = 0; r < REGS / 2; r++) {
                rows[r] = _mm512_shuffle_i32x4(transposed[r], transposed[r + 8],
                                               0x44);
              }
#pragma unroll
              for (int r = 0; r < REGS / 2; r++) {
                transposed[r] = _mm512_shuffle_i32x4(rows[r], rows[r], 0x28);
                transposed[r + 8] =
                    _mm512_shuffle_i32x4(rows[r], rows[r], 0x9D);
              }
              constexpr int perm[16] = {0, 4,  2,  6,  1, 5,  3,  7,
                                        8, 12, 10, 14, 9, 13, 11, 15};
#pragma unroll
              for (int r = 0; r < REGS; r++) {
                __m256i ymm = _mm512_castsi512_si256(transposed[r]);
                _mm256_storeu_si256((__m256i*)&smem[(i + perm[r]) * TILE_K + j],
                                    ymm);
              }
            }
          }
        }
      } else {
        // Tensor B packing
        if constexpr (AXIS == MajorAxis::ROW) {
          // B is row-major (K x N): tile at (tk * TILE_K, tmn * TILE_MN)
          const T* __restrict__ tile_ptr = _ptr + tk * TILE_K * _dim_mn + tmn * TILE_MN;
          // repack B with interleave for AMX VNNI format
#pragma unroll
          for (int i = 0; i < TILE_K; i += 4) {
#pragma unroll
            for (int j = 0; j < TILE_MN; j += 32) {
              __m512i row1 = _mm512_loadu_epi16(&tile_ptr[i * _dim_mn + j]);
              __m512i row2 =
                  _mm512_loadu_epi16(&tile_ptr[(i + 1) * _dim_mn + j]);
              __m512i row3 =
                  _mm512_loadu_epi16(&tile_ptr[(i + 2) * _dim_mn + j]);
              __m512i row4 =
                  _mm512_loadu_epi16(&tile_ptr[(i + 3) * _dim_mn + j]);
              __m512i interleaved1 = _mm512_unpacklo_epi16(row1, row2);
              __m512i interleaved2 = _mm512_unpackhi_epi16(row1, row2);
              __m512i interleaved3 = _mm512_unpacklo_epi16(row3, row4);
              __m512i interleaved4 = _mm512_unpackhi_epi16(row3, row4);
              __m512i shuff1 =
                  _mm512_shuffle_i32x4(interleaved1, interleaved2, 0x44);
              __m512i shuff2 =
                  _mm512_shuffle_i32x4(interleaved1, interleaved2, 0xEE);
              __m512i shuff3 =
                  _mm512_shuffle_i32x4(interleaved3, interleaved4, 0x44);
              __m512i shuff4 =
                  _mm512_shuffle_i32x4(interleaved3, interleaved4, 0xEE);
              __m512i out1 = _mm512_shuffle_i32x4(shuff1, shuff1, 0xD8);
              __m512i out2 = _mm512_shuffle_i32x4(shuff2, shuff2, 0xD8);
              __m512i out3 = _mm512_shuffle_i32x4(shuff3, shuff3, 0xD8);
              __m512i out4 = _mm512_shuffle_i32x4(shuff4, shuff4, 0xD8);
              _mm512_storeu_si512(
                  (__m512i*)&smem[(i / 2) * (2 * TILE_MN) + 2 * j], out1);
              _mm512_storeu_si512(
                  (__m512i*)&smem[(i / 2) * (2 * TILE_MN) + 2 * j + 32], out2);
              _mm512_storeu_si512(
                  (__m512i*)&smem[(i / 2 + 1) * (2 * TILE_MN) + 2 * j], out3);
              _mm512_storeu_si512(
                  (__m512i*)&smem[(i / 2 + 1) * (2 * TILE_MN) + 2 * j + 32],
                  out4);
            }
          }
        } else {
          // transpose B from column-major (N x K) to row-major with
          // interleaving
          const T* __restrict__ tile_ptr = _ptr + tmn * TILE_MN + tk * TILE_K * _dim_mn;
          __mmask32 k1 = 0x0000FFFF;
          constexpr int REGS = 16;
#pragma unroll
          for (int i = 0; i < TILE_K; i += 16) {
#pragma unroll
            for (int j = 0; j < TILE_MN; j += 16) {
              __m512i transposed[REGS];
              __m512i rows[REGS];
              for (int r = 0; r < REGS; r++) {
                rows[r] = _mm512_maskz_loadu_epi16(
                    k1, &tile_ptr[(i + r) * _dim_mn + j]);
              }
              // Stage 1: epi32 interleave
#pragma unroll
              for (int r = 0; r < REGS; r += 2) {
                transposed[r] = _mm512_unpacklo_epi32(rows[r], rows[r + 1]);
                transposed[r + 1] = _mm512_unpackhi_epi32(rows[r], rows[r + 1]);
              }
              // Stage 2: epi64 interleave
#pragma unroll
              for (int r = 0; r < REGS; r += 4) {
                for (int k = 0; k < 2; ++k) {
                  rows[r + k] = _mm512_unpacklo_epi64(transposed[r + k],
                                                      transposed[r + k + 2]);
                  rows[r + k + 2] = _mm512_unpackhi_epi64(
                      transposed[r + k], transposed[r + k + 2]);
                }
              }
              // Stage 3: Shuffle to combine 128-bit lanes
#pragma unroll
              for (int r = 0; r < REGS / 4; r++) {
                transposed[r] =
                    _mm512_shuffle_i32x4(rows[r], rows[r + 4], 0x44);
                transposed[r + 8] =
                    _mm512_shuffle_i32x4(rows[r + 8], rows[r + 12], 0x44);
              }
#pragma unroll
              for (int r = 0; r < REGS / 4; r++) {
                rows[r] =
                    _mm512_shuffle_i32x4(transposed[r], transposed[r], 0x28);
                rows[r + 4] =
                    _mm512_shuffle_i32x4(transposed[r], transposed[r], 0x9D);
                rows[r + 8] = _mm512_shuffle_i32x4(transposed[r + 8],
                                                   transposed[r + 8], 0x28);
                rows[r + 12] = _mm512_shuffle_i32x4(transposed[r + 8],
                                                    transposed[r + 8], 0x9D);
              }
              // Permutation to correct column ordering
              constexpr int perm[16] = {0, 2, 1, 3, 4, 6, 5, 7};
#pragma unroll
              for (int r = 0; r < REGS / 2; r++) {
                __m256i ymm_1 = _mm512_castsi512_si256(rows[r]);
                __m256i ymm_2 = _mm512_castsi512_si256(rows[r + 8]);
                _mm256_storeu_si256(
                    (__m256i*)&smem[(j / 2 + perm[r]) * (2 * TILE_MN) +
                                    (2 * i)],
                    ymm_1);
                _mm256_storeu_si256(
                    (__m256i*)&smem[(j / 2 + perm[r]) * (2 * TILE_MN) +
                                    (2 * i + 16)],
                    ymm_2);
              }
            }
          }
        }
      }
    } else if constexpr (std::is_same_v<T, uint8_t> ||
                         std::is_same_v<T, int8_t>) {
      static_assert(sizeof(T) == 2, "Not supported for now");
    } else {
      static_assert(sizeof(T) == 0, "Unsupported type for packing");
    }
  }
};

template <typename T_A, typename T_B, typename T_D, typename AMXTC_A,
          typename AMXTC_B, typename AMXTC_D, int TILE_M, int TILE_N,
          int TILE_K>
class TensorMMA {
  static constexpr MajorAxis AXIS_D = AMXTC_D::major_axis;
  static_assert(AXIS_D == MajorAxis::ROW,
                "Only ROW major output is supported for now");

  // AMX tile dimensions
  static constexpr int TILE_M_2 = TILE_M / 2;  // 16 for 32x32 tiles
  static constexpr int TILE_N_2 = TILE_N / 2;  // 16 for 32x32 tiles
  static constexpr int TILE_K_2 = TILE_K / 2;  // For VNNI B layout

 public:
  // MMA API matching SME: takes packed A/B pointers and D fragment
  static __forceinline void mma(
      TensorFragment<Tensor::A, T_A, T_D, TILE_M, TILE_K, AMXTC_A>& frag_a,
      TensorFragment<Tensor::B, T_B, T_D, TILE_K, TILE_N, AMXTC_B>& frag_b,
      TensorFragment<Tensor::D, T_A, T_D, TILE_M, TILE_N, AMXTC_D>& frag_d) {
    // (void)frag_d;  // frag_d state is implicit in tile registers
    if constexpr (std::is_same_v<T_A, _Float16> ||
                  std::is_same_v<T_A, __bfloat16>) {
      // // Load A tiles (contiguously packed: row-major TILE_M x TILE_K)
      // // Tiles 0,1: each 16 rows x 32 cols of FP16
      // _tile_loadd(0, ptr_a, TILE_K * sizeof(T_A));
      // _tile_loadd(1, ptr_a + TILE_M_2 * TILE_K, TILE_K * sizeof(T_A));

      // // Load B tiles (VNNI packed: (TILE_K/2) x (2*TILE_N))
      // // Tiles 2,3: each (TILE_K/2) rows x 32 cols of FP16-pairs
      // _tile_loadd(2, ptr_b, 2 * TILE_N * sizeof(T_B));
      // _tile_loadd(3, ptr_b + TILE_N_2 * 2, 2 * TILE_N * sizeof(T_B));

      // Compute: D += A * B
      if constexpr (std::is_same_v<T_A, _Float16>) {
        _tile_dpfp16ps(4, 0, 2);
        _tile_dpfp16ps(5, 0, 3);
        _tile_dpfp16ps(6, 1, 2);
        _tile_dpfp16ps(7, 1, 3);
      } else {
        _tile_dpbf16ps(4, 0, 2);
        _tile_dpbf16ps(5, 0, 3);
        _tile_dpbf16ps(6, 1, 2);
        _tile_dpbf16ps(7, 1, 3);
      }
    } else {
      static_assert(sizeof(T_A) == 0, "Unsupported type for MMA");
    }
  }
};

}  // namespace core_ir
