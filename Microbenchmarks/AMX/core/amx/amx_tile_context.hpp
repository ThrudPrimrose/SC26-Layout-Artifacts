#pragma once

#include <immintrin.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <type_traits>

#ifndef __forceinline
#define __forceinline inline __attribute__((always_inline))
#endif

namespace core_ir {

// ============================================================================
// Intel AMX-SPECIFIC ENUMS - Standalone at namespace level
// ============================================================================

// Major axis for memory layout
enum class MajorAxis { ROW = 0, COLUMN = 1 };

// Tensor role in matrix operations
enum class Tensor { A = 0, B = 1, D = 4 };

template <MajorAxis AXIS>
struct AMXTileContext {
  static constexpr MajorAxis major_axis = AXIS;
};

// ============================================================================
// Default AMX Parameters
// ============================================================================

#define ARCH_GET_XCOMP_PERM 0x1022
#define ARCH_REQ_XCOMP_PERM 0x1023
#define XFEATURE_XTILECFG 17
#define XFEATURE_XTILEDATA 18

// Define tile config data structure
typedef struct __tile_config {
  uint8_t palette_id;
  uint8_t start_row;
  uint8_t reserved_0[14];
  uint16_t colsb[16];
  uint8_t rows[16];
} __tilecfg;

// Initialize tile config
template <typename T, int TILE_M, int TILE_N, int TILE_K>
static __forceinline void init_tile_config(__tilecfg* __restrict__ tileinfo) {
  int i;
  tileinfo->palette_id = 1;
  tileinfo->start_row = 0;

  static_assert(TILE_M <= 16 && TILE_N <= 16 && TILE_K <= 32,
                "Tile dimensions must fit within AMX limits");

  if constexpr (std::is_same_v<T, _Float16> || std::is_same_v<T, __bfloat16>) {
    // A tiles: 2 tiles of t_m * t_k
    for (i = 0; i < 2; i++) {
      tileinfo->colsb[i] = TILE_K * sizeof(T);  // upto 64 bytes per row
      tileinfo->rows[i] = TILE_M;               // upto 16 rows
    }
    // B tiles: 2 tiles of t_k * t_n rearranged as (t_k / 2) * (2 * t_n)
    for (i = 2; i < 4; i++) {
      tileinfo->colsb[i] = (2 * TILE_N) * sizeof(T);
      tileinfo->rows[i] = TILE_K / 2;
    }
    // C tiles: 4 tiles of t_m * t_n
    for (i = 4; i < 8; i++) {
      tileinfo->colsb[i] = TILE_N * sizeof(float);  // upto 64 bytes per row
      tileinfo->rows[i] = TILE_M;                   // upto 16 rows
    }

  } else if constexpr (std::is_same_v<T, uint8_t> ||
                       std::is_same_v<T, int8_t>) {
    static_assert(sizeof(T) == 2, "Not supported for now");
  } else {
    static_assert(sizeof(T) == 0, "Unsupported type for tile config");
  }

  _tile_loadconfig(tileinfo);
}

// Set_tiledata_use() - Invoke syscall to set ARCH_SET_STATE_USE
static bool set_tiledata_use() {
  if (syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA)) {
    // printf("\n Fail to do XFEATURE_XTILEDATA \n\n");
    return false;
  } else {
    // printf("\n TILE DATA USE SET - OK \n\n");
    return true;
  }
  return true;
}

// ============================================================================
// Default context alias for convenience
// ============================================================================
using DefaultAMXContext = AMXTileContext<MajorAxis::ROW>;

enum class StrideMajor { K, MN };  // For layout packing

}  // namespace core_ir
