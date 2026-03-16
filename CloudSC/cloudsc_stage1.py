from typing import Dict

import numpy as np
import dace
from pathlib import Path
import copy
from dace.transformation.dataflow import MapUnroll
from dace.transformation.interstate.loop_unroll import LoopUnroll
from dace.transformation.passes import ScalarToSymbolPromotion
from dace.transformation.layout.split_array import LoopRegion, SplitArray
from dace.transformation.passes.analysis import loop_analysis
from dace.transformation.passes.clean_access_node_to_scalar_slice_to_tasklet_pattern import CleanAccessNodeToScalarSliceToTaskletPattern



def unroll(sdfg: dace.SDFG, symbol_map: Dict[str, int]):
    """Unroll maps and loops whose iteration range matches a split-dimension extent.
    After ``replace_dict``, extents that depended only on split symbols become
    concrete integers. We unroll one at a time and rescan, because each
    unroll mutates the graph and invalidates node references.
    """
    array_dim_map = dict()
    for array_name, desc in sdfg.arrays.items():
        array_dim_map[array_name] = copy.deepcopy(desc.shape)

    sdfg.replace_dict(symbol_map)
    potential_ranges = set(symbol_map.values())

    while True:
        target = None
        target_extent = None

        for n, g in sdfg.all_nodes_recursive():
            if isinstance(n, dace.nodes.MapEntry):
                has_split_dim = False
                for _, r in zip(n.map.params, n.map.range):
                    extent = ((r[1] + 1) - r[0]) // r[2]
                    if extent.free_symbols:
                        continue
                    try:
                        val = int(extent)
                    except (TypeError, ValueError):
                        continue
                    if val in potential_ranges:
                        has_split_dim = True
                        target_extent = val
                        break
                if has_split_dim:
                    target = ("map", n, g)
                    break

            elif isinstance(n, LoopRegion):
                beg = loop_analysis.get_init_assignment(n)
                end = loop_analysis.get_loop_end(n)
                step = loop_analysis.get_loop_stride(n)
                if beg is None or end is None or step is None:
                    continue
                extent = ((end + 1) - beg) // step
                if extent.free_symbols:
                    continue
                try:
                    val = int(extent)
                except (TypeError, ValueError):
                    continue
                if val in potential_ranges:
                    target = ("loop", n, g)
                    target_extent = val
                    break

        if target is None:
            break

        kind, n, g = target
        if kind == "map":
            MapUnroll().apply_to(sdfg=g.sdfg, options={}, map_entry=n)
        else:
            LoopUnroll().apply_to(sdfg=n.sdfg, options={"inline_iterations": True}, loop=n)


if __name__ == '__main__':
    NAME_ORDER = ["ncldql", "ncldqi", "ncldqr", "ncldqs", "ncldqv"]
    NAME_MAP = {i: NAME_ORDER[i] for i in range(5)}
    # Split along the first dimension (size = nclv-1) into separate 2D arrays
    symbol_map = {
        "nclv": 5,
        "ncldql": 1,
        "ncldqi": 2,
        "ncldqr": 3,
        "ncldqs": 4,
        "ncldqv": 5,
    }
    name_map = {"nclv": NAME_ORDER}

    sdfg = dace.SDFG.from_file("cloudsc_cleaned.sdfgz")
    sdfg.validate()
    unroll(sdfg, symbol_map=symbol_map)
    sdfg.name += "_unrolled"
    sdfg.validate()
    sdfg.save("cloudsc_unrolled.sdfgz", compress=True)