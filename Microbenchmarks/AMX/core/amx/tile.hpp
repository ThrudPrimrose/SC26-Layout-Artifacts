#pragma once

#include "amx_tile_context.hpp"
#include "tensor.hpp"

namespace core_ir {

template <typename T_A, typename T_B, typename T_C, typename T_D,
          MajorAxis AXIS_A, MajorAxis AXIS_B, MajorAxis AXIS_C,
          MajorAxis AXIS_D, int TILE_M, int TILE_N, int TILE_K>
struct TileOps {
  using AMXTC_A = AMXTileContext<AXIS_A>;
  using AMXTC_B = AMXTileContext<AXIS_B>;
  using AMXTC_C = AMXTileContext<AXIS_C>;
  using AMXTC_D = AMXTileContext<AXIS_D>;

  using TileLayoutA = TensorLayout<T_A, AMXTC_A, Tensor::A, TILE_M, TILE_K>;
  using TileLayoutB = TensorLayout<T_B, AMXTC_B, Tensor::B, TILE_N, TILE_K>;
  using FragmentA =
      TensorFragment<Tensor::A, T_A, T_C, TILE_M, TILE_K, AMXTC_A>;
  using FragmentB =
      TensorFragment<Tensor::B, T_B, T_C, TILE_K, TILE_N, AMXTC_B>;
  using FragmentD =
      TensorFragment<Tensor::D, T_A, T_D, TILE_M, TILE_N, AMXTC_D>;
  using MMA = TensorMMA<T_A, T_B, T_D, AMXTC_A, AMXTC_B, AMXTC_D, TILE_M,
                        TILE_N, TILE_K>;

  static void init() {
    // Request permission to linux kernel to run AMX
    if (!set_tiledata_use()) exit(-1);

    // Load tile configuration
    static __tilecfg tile_data = {0};
    init_tile_config<T_A, std::min(16, TILE_M / 2), std::min(16, TILE_N / 2),
                     std::min(32, TILE_K)>(&tile_data);
  }
};

}  // namespace core_ir