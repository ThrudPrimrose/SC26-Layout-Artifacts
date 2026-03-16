"""lowprec.py — Precision lowering for CloudSC CPU pipeline.

SDFG-level transforms to lower float64 arrays/scalars to a target dtype.
Ported from velocity/utils/boundary_cast.py and compile_gpu_stage8.py.

Strategy (same as velocity):
  - Lower EVERYTHING we can to the target dtype
  - Non-transient arrays (function parameters / ABI) stay double externally
  - Boundary cast Maps at read/write sites convert double ↔ lowered dtype
  - Transient arrays are lowered directly (no cast needed)
  - Scalars: top-level params stay double, nested transients are lowered
"""

from __future__ import annotations

import re

import dace
from dace import dtypes, nodes


LOWPREC_MAP = {
    "fp64": dace.float64,
    "fp32": dace.float32,
    "fp16": dace.float16,
    "f32": dace.float32,
    "f64": dace.float64,
    "f16": dace.float16,
    # TODO: BFP modes
    # "bfp8": dace.float16,
    # "bfp16": dace.float16,
    # "bfp32": dace.float32,
}


# --- Core dtype propagation (from velocity/utils/boundary_cast.py) ---


def _propagate_dtype(
    sdfg: dace.SDFG,
    array_name: str,
    new_dtype: dace.typeclass,
    visited: set | None = None,
):
    """Recursively change an array's dtype through all NestedSDFGs.

    When an array's dtype is changed at the top level, NestedSDFGs that
    access it via connectors still have the old dtype in their internal
    array descriptors.  DaCe codegen uses the internal descriptors for
    parameter types, so they must match the allocation dtype.

    This function walks the SDFG hierarchy and updates:
    - The array descriptor dtype at each level
    - NestedSDFG connector-mapped internal arrays (recursively)
    """
    if visited is None:
        visited = set()
    key = (id(sdfg), array_name)
    if key in visited:
        return
    visited.add(key)

    if array_name not in sdfg.arrays:
        return
    sdfg.arrays[array_name].dtype = new_dtype

    for state in sdfg.states():
        for node in state.nodes():
            if not isinstance(node, nodes.NestedSDFG):
                continue
            for edge in state.in_edges(node):
                if edge.data.data == array_name and edge.dst_conn:
                    _propagate_dtype(node.sdfg, edge.dst_conn, new_dtype, visited)
            for edge in state.out_edges(node):
                if edge.data.data == array_name and edge.src_conn:
                    _propagate_dtype(node.sdfg, edge.src_conn, new_dtype, visited)


# --- Boundary cast for CPU non-transient arrays ---
#
# Same strategy as velocity's H2D/D2H boundary cast, treating the lowered
# computation as a "device":
#
#   Entry state ("H2D"): AN(param, double) → [cast Map] → AN(param_lowered, float)
#   ... all computation states use param_lowered ...
#   Exit state ("D2H"):  AN(param_lowered, float) → [cast Map] → AN(param, double)
#
# Non-transient arrays keep double in the function signature (ABI).
# All internal references are renamed to the lowered transient.


def inject_cpu_boundary_cast(
    sdfg: dace.SDFG,
    array_names: list[str],
    external_dtype: dace.typeclass,
):
    """Insert boundary cast states for non-transient CPU arrays.

    For each non-transient array:
    1. Create a lowered transient ``{name}_lowered`` with external_dtype
    2. Rename ALL internal references (AccessNodes + memlets) to the lowered name
    3. Propagate the lowered dtype into NestedSDFGs
    4. After all renames, add entry/exit cast states at the SDFG boundary

    The original non-transient array stays double (ABI preserved).
    All computation runs in the lowered dtype.
    """
    lowered_pairs = []  # (orig_name, lowered_name) for entry/exit state creation

    for name in array_names:
        if name not in sdfg.arrays:
            print(f"  cpu_boundary_cast: {name} not in SDFG, skipping.")
            continue

        arr = sdfg.arrays[name]
        orig_dtype = arr.dtype
        lowered_name = f"{name}_lowered"

        # Create lowered transient with same shape
        sdfg.add_array(
            lowered_name,
            shape=arr.shape,
            dtype=external_dtype,
            transient=True,
            storage=dtypes.StorageType.CPU_Heap,
        )

        print(
            f"  cpu_boundary_cast: {name} ({orig_dtype}) → "
            f"{lowered_name} ({external_dtype})"
        )

        # Rename ALL AccessNodes and memlets in ALL states
        _rename_array_references(sdfg, name, lowered_name)

        # Propagate lowered dtype into NestedSDFGs
        _propagate_dtype(sdfg, lowered_name, external_dtype)

        lowered_pairs.append((name, lowered_name))

    if not lowered_pairs:
        return

    # Add entry state: cast double→float for all inputs ("H2D")
    _add_entry_cast_state(sdfg, lowered_pairs, external_dtype)

    # Add exit state: cast float→double for all outputs ("D2H")
    _add_exit_cast_state(sdfg, lowered_pairs)

    # Re-index internal CFG list after adding new states
    sdfg.reset_cfg_list()


def _rename_array_references(sdfg: dace.SDFG, old_name: str, new_name: str):
    """Rename all AccessNodes and memlets from old_name to new_name in all states."""
    for state in sdfg.states():
        for node in state.nodes():
            if isinstance(node, nodes.AccessNode) and node.data == old_name:
                node.data = new_name
        for edge in state.edges():
            if edge.data.data == old_name:
                edge.data.data = new_name


def _add_entry_cast_state(
    sdfg: dace.SDFG,
    pairs: list[tuple[str, str]],
    external_dtype: dace.typeclass,
):
    """Add a new entry state with cast Maps: orig(double) → lowered(float).

    Inserts before the current start state, like an H2D copy state.
    """
    # Find the current source states (no incoming interstate edges)
    source_states = sdfg.source_nodes()

    cast_state = sdfg.add_state("boundary_cast_entry", is_start_state=True)

    for orig_name, lowered_name in pairs:
        orig_arr = sdfg.arrays[orig_name]

        src_an = cast_state.add_access(orig_name)
        dst_an = cast_state.add_access(lowered_name)

        _add_cast_map(
            sdfg, cast_state,
            src_an, orig_name,
            dst_an, lowered_name,
            orig_arr.shape,
            external_dtype,
            f"h2d_{orig_name}",
        )

    # Connect to previous source states
    for src_state in source_states:
        sdfg.add_edge(cast_state, src_state, dace.InterstateEdge())


def _add_exit_cast_state(
    sdfg: dace.SDFG,
    pairs: list[tuple[str, str]],
):
    """Add a new exit state with cast Maps: lowered(float) → orig(double).

    Appends after the current sink states, like a D2H copy state.
    """
    # Find the current sink states (no outgoing interstate edges)
    sink_states = sdfg.sink_nodes()

    cast_state = sdfg.add_state("boundary_cast_exit")

    for orig_name, lowered_name in pairs:
        orig_arr = sdfg.arrays[orig_name]

        src_an = cast_state.add_access(lowered_name)
        dst_an = cast_state.add_access(orig_name)

        _add_cast_map(
            sdfg, cast_state,
            src_an, lowered_name,
            dst_an, orig_name,
            orig_arr.shape,
            orig_arr.dtype,  # cast back to double
            f"d2h_{orig_name}",
        )

    # Connect from previous sink states
    for sink_state in sink_states:
        sdfg.add_edge(sink_state, cast_state, dace.InterstateEdge())


def _add_cast_map(
    sdfg: dace.SDFG,
    state: dace.SDFGState,
    src_an: nodes.AccessNode,
    src_name: str,
    dst_an: nodes.AccessNode,
    dst_name: str,
    shape: tuple,
    target_dtype: dace.typeclass,
    label: str,
):
    """Add a Map that casts every element from src to dst.

    Ported from velocity/utils/boundary_cast.py _add_cast_map.
    """
    src_desc = sdfg.arrays[src_name]
    dst_desc = sdfg.arrays[dst_name]

    map_ranges = {f"__i{i}": f"0:{s}" for i, s in enumerate(shape)}
    indices = ", ".join(f"__i{i}" for i in range(len(shape)))

    map_entry, map_exit = state.add_map(label, map_ranges)

    tasklet = state.add_tasklet(
        name=label,
        inputs={"_in"},
        outputs={"_out"},
        code=f"_out = static_cast<{target_dtype.ctype}>(_in);",
        language=dtypes.Language.CPP,
    )

    # Map entry connectors
    in_conn = f"IN_{src_name}"
    out_conn = f"OUT_{src_name}"
    map_entry.add_in_connector(in_conn)
    map_entry.add_out_connector(out_conn)

    # Map exit connectors
    exit_in_conn = f"IN_{dst_name}"
    exit_out_conn = f"OUT_{dst_name}"
    map_exit.add_in_connector(exit_in_conn)
    map_exit.add_out_connector(exit_out_conn)

    # Wire: src_an → map_entry
    state.add_edge(
        src_an, None, map_entry, in_conn,
        dace.Memlet.from_array(src_name, src_desc),
    )

    # Wire: map_entry → tasklet
    state.add_edge(
        map_entry, out_conn, tasklet, "_in",
        dace.Memlet(f"{src_name}[{indices}]"),
    )

    # Wire: tasklet → map_exit
    state.add_edge(
        tasklet, "_out", map_exit, exit_in_conn,
        dace.Memlet(f"{dst_name}[{indices}]"),
    )

    # Wire: map_exit → dst_an
    state.add_edge(
        map_exit, exit_out_conn, dst_an, None,
        dace.Memlet.from_array(dst_name, dst_desc),
    )


# --- Direct lowering for transient arrays ---


def lower_transient_arrays(
    sdfg: dace.SDFG,
    external_dtype: dace.typeclass,
    exclude: set[str] | None = None,
):
    """Lower transient float64 arrays directly (no cast needed — internal only)."""
    exclude = exclude or set()

    candidates = [
        name
        for name, arr in sdfg.arrays.items()
        if isinstance(arr, dace.data.Array)
        and arr.dtype == dace.float64
        and arr.transient
        and arr.total_size != 1
        and name not in exclude
    ]

    if candidates:
        print(f"Lowering {len(candidates)} transient arrays to {external_dtype}")
        for name in candidates:
            _propagate_dtype(sdfg, name, external_dtype)

    return candidates


# --- Top-level transient scalar lowering ---


def lower_top_level_scalars(
    sdfg: dace.SDFG,
    external_dtype: dace.typeclass,
    exclude: set[str] | None = None,
):
    """Lower transient double scalars at the root SDFG level.

    These are intermediates like zvqx_index, zvqx_slice that DaCe creates
    for array element indexing.  They must match the dtype of the arrays
    they feed into, otherwise CopyND emits mismatched types.

    Non-transient scalars are function parameters — they stay double (ABI).
    """
    exclude = exclude or set()
    count = 0

    for name, arr in list(sdfg.arrays.items()):
        if (
            isinstance(arr, dace.data.Scalar)
            and arr.dtype == dace.float64
            and arr.transient
            and name not in exclude
            and not any(name.startswith(ex + "_") for ex in exclude)
        ):
            arr.dtype = external_dtype
            count += 1

    if count:
        print(f"Lowered {count} top-level transient scalars to {external_dtype}")


# --- Scalar lowering (from compile_gpu_stage8.py) ---


def lower_nested_scalars(
    sdfg: dace.SDFG,
    external_dtype: dace.typeclass,
    exclude: set[str] | None = None,
):
    """Lower transient double scalars inside nested SDFGs.

    Top-level non-transient scalars are function parameters — they stay double
    (ABI). Only transient scalars in nested SDFGs (kernel intermediates) are
    lowered. Comparison results get int32 instead (logically boolean).

    Ported from compile_gpu_stage8.py.
    """
    exclude = exclude or set()

    # First pass: find scalars that are comparison outputs
    comparison_scalars: set[str] = set()
    for nsdfg in sdfg.all_sdfgs_recursive():
        if nsdfg is sdfg:
            continue
        for state in nsdfg.states():
            for node in state.nodes():
                if not isinstance(node, nodes.Tasklet):
                    continue
                code = node.code.as_string
                if not any(
                    op in code
                    for op in (" > ", " < ", " >= ", " <= ", " == ", " != ")
                ):
                    continue
                for e in state.out_edges(node):
                    if (
                        isinstance(e.dst, nodes.AccessNode)
                        and e.dst.data in nsdfg.arrays
                        and isinstance(nsdfg.arrays[e.dst.data], dace.data.Scalar)
                    ):
                        comparison_scalars.add(e.dst.data)

    count = 0
    for nsdfg in sdfg.all_sdfgs_recursive():
        if nsdfg is sdfg:
            continue
        for name, arr in list(nsdfg.arrays.items()):
            if (
                isinstance(arr, dace.data.Scalar)
                and arr.dtype == dace.float64
                and arr.transient
                and name not in exclude
            ):
                if name in comparison_scalars:
                    arr.dtype = dace.int32
                else:
                    arr.dtype = external_dtype
                count += 1

    if count:
        print(f"Lowered {count} nested transient scalars")


# --- Floatify tasklet code (from compile_gpu_stage8.py) ---


def floatify_tasklets(sdfg: dace.SDFG):
    """Replace double literals and math functions in tasklets with
    float equivalents to prevent FP64 promotion.

    Ported from compile_gpu_stage8.py — handles both CPP and Python tasklets.
    """
    constants = ["0.5", "0.85", "1.0", "0.0", "0.05", "0.65", "1.15"]

    for node, _state in sdfg.all_nodes_recursive():
        if not isinstance(node, nodes.Tasklet):
            continue

        code = node.code.as_string
        new_code = code

        if node.language == dace.Language.CPP:
            for val in constants:
                new_code = re.sub(
                    rf"(?<![0-9fF]){re.escape(val)}(?![0-9fF])",
                    f"{val}f",
                    new_code,
                )
                new_code = re.sub(
                    rf"(?<![0-9fF])\-{re.escape(val)}(?![0-9fF])",
                    f"-{val}f",
                    new_code,
                )
            # abs() -> fabsf()
            new_code = re.sub(r"\babs\(", "fabsf(", new_code)
            # ipow(x, 2) -> (x * x)
            new_code = re.sub(
                r"dace::math::ipow\(([^,]+),\s*2\)",
                r"((\1) * (\1))",
                new_code,
            )
        else:
            # Python tasklets: float() wrapper
            for val in constants:
                new_code = re.sub(
                    rf"(?<![0-9fF]){re.escape(val)}(?![0-9fF])",
                    f"float({val})",
                    new_code,
                )
                new_code = re.sub(
                    rf"(?<![0-9fF])\-{re.escape(val)}(?![0-9fF])",
                    f"float(-{val})",
                    new_code,
                )
            # x ** 2 or x ** 2.0 -> (x * x)
            new_code = re.sub(
                r"(\w+)\s*\*\*\s*2(?:\.0)?\b",
                r"((\1) * (\1))",
                new_code,
            )

        if new_code != code:
            node.code = dace.properties.CodeBlock(new_code, node.language)


# --- Top-level entry point ---


def apply_lowprec(sdfg: dace.SDFG, lowprec: str):
    """Apply all precision lowering transforms to the SDFG.

    Strategy (mirrors velocity):
      1. Non-transient arrays: boundary cast (keep double at ABI, cast to lowered)
      2. Transient arrays: lower dtype directly
      3. Nested transient scalars: lower dtype directly
      4. Floatify tasklet literals

    Args:
        sdfg: The SDFG to transform.
        lowprec: Precision mode string (fp64, fp32, fp16, etc.)
    """
    if lowprec not in LOWPREC_MAP:
        raise ValueError(
            f"Unknown lowprec value: {lowprec!r}. "
            f"Valid options: {', '.join(LOWPREC_MAP)}"
        )

    external_dtype = LOWPREC_MAP[lowprec]

    # TODO: BFP mode — separate path like velocity's inject_bfp_packing
    # if lowprec.startswith("bfp"):
    #     ...

    print(f"Applying precision lowering: fp64 → {external_dtype}")

    # 0. Remove the return block — boundary_cast_exit will be the new sink
    from dace.sdfg.state import ReturnBlock
    return_blocks = [s for s in sdfg.nodes() if isinstance(s, ReturnBlock)]
    assert len(return_blocks) == 1, (
        f"Expected exactly 1 ReturnBlock, found {len(return_blocks)}"
    )
    rb = return_blocks[0]
    for edge in list(sdfg.in_edges(rb)):
        sdfg.remove_edge(edge)
    sdfg.remove_node(rb)
    print(f"  Removed ReturnBlock: {rb.label}")

    # Exclusion list: arrays/scalars that must stay fp64.
    # Everything NOT in this set gets lowered.
    _LOWERING_EXCLUDE = {
        # === MUST STAY FP64 ===
        # high-sensitivity inputs (SNR < 10 dB)
        "pt", "pa", "tendency_tmp_a",
        # accumulation chain (SNR -51..52 dB)
        "zsolqa", "psum_solqa", "za", "zlcust",
        # moisture (tiny values, SNR 18 dB)
        "zqxfg",
        # temperature (large magnitude, small differences, SNR 75 dB)
        "ztp1",
        # condensation / evaporation (SNR 51..70 dB)
        "zlcond1", "zlcond2",
        # input params (SNR 62..105 dB)
        "psupsat", "pq", "pap", "paph",
        # thermodynamic (SNR ~105 dB)
        # "zqlhs", "zdqsliqdt", "zdqsicedt",
        #
        # === SAFE: verified lowerable to fp32 ===
        # (commented out — kept for reference)
        #
        # -- flux outputs --
        # "pfsqlf", "pfsqrf", "pfsqif", "pfsqsf",
        # "pfcqlng", "pfcqrng", "pfcqnng", "pfcqsng",
        # "pfhpsn", "pfplsn",
        # -- high-sensitivity inputs (safe subgroups) --
        # "tendency_tmp_q", "tendency_tmp_t", "tendency_tmp_cld",
        # "pclv", "plude",
        # "phrlw", "pvfl", "pvfi",
        # -- accumulation chain --
        # "zgdph_r", "zalfaw", "zqx", "zexplicit",
        # "zqxn2d", "zqx0", "zqxn", "zqxnm1",
        # "zlneg", "zfoealfa",
        # "zconvsrce", "zconvsink",
        # -- moisture --
        # "zvqx",
        # -- temperature / pressure --
        # "ztold", "zdp",
        # -- solution arrays --
        # "zsolqb", "zratio",
        # "zsinksum", "zfallsink", "zfallsrce",
        # -- phase fractions --
        # "zliqfrac", "zicefrac", "zli", "zlfinalsum",
        # -- flux arrays --
        # "zfluxq", "zpfplsx",
        # -- thermodynamic (dqs, saturation, foeew) --
        # "zqold", "zdtgdp", "zrdtgdp",
        # "zdqs", "zdqsmixdt",
        # "zqsmix", "zqsliq", "zqsice",
        # "zfoeewmt", "zfoeew", "zfoeeliqt",
        # -- condensation / evaporation --
        # "zrainaut", "zsnowaut",
        # -- cloud fractions / cover --
        # "zliqcld", "zicecld", "zlicld",
        # "zcovpclr", "zcovptot", "zcovpmax",
        # "zaorig", "zanewm1", "zda", "zacust",
        # -- saturation / supersaturation --
        # "zsupsat", "zcorqsliq", "zcorqsice", "zcorqsmix",
        # -- ice nucleation --
        # "zfokoop", "zicenuclei", "zicetot",
        # -- precipitation / fluxes --
        # "zqpretot", "zldefr", "zldifdt", "zmf", "zrho",
        # "zsolab", "zsolac", "zmeltmax", "zfrzmax",
        # "zevaplimice", "zevaplimmix", "zcldtopdist",
        # "zrainacc", "zraincld", "zsnowcld", "zsnowrime",
        # -- pressure + other --
        # "zgdp", "zpsupsatsrce",
        # -- non-transient input arrays --
        # "pvfa", "pdyna", "pdynl", "pdyni",
        # "phrsw", "pvervel", "plsm", "plu", "psnde",
        # "pmfu", "pmfd",
        # "plcrit_aer", "picrit_aer",
        # "pre_ice", "pccn", "pnice",
        # "pcovptot", "prainfrac_toprfz",
        # "pfsqltur", "pfsqitur", "pfplsl", "pfhpsl",
        # -- output arrays --
        # "tendency_loc_t", "tendency_loc_q",
        # "tendency_loc_a", "tendency_loc_cld",
    }

    # --- Scan: find everything that CAN be lowered ---

    # All non-transient fp64 arrays (candidates for boundary cast)
    all_non_transient = [
        name
        for name, arr in sdfg.arrays.items()
        if isinstance(arr, dace.data.Array)
        and arr.dtype == dace.float64
        and not arr.transient
        and arr.total_size != 1
    ]

    # All transient fp64 arrays
    all_transient_arrays = [
        name
        for name, arr in sdfg.arrays.items()
        if isinstance(arr, dace.data.Array)
        and arr.dtype == dace.float64
        and arr.transient
        and arr.total_size != 1
    ]

    # All top-level transient fp64 scalars
    all_top_scalars = [
        name
        for name, arr in sdfg.arrays.items()
        if isinstance(arr, dace.data.Scalar)
        and arr.dtype == dace.float64
        and arr.transient
    ]

    # All nested transient fp64 scalars
    all_nested_scalars = []
    for nsdfg in sdfg.all_sdfgs_recursive():
        if nsdfg is sdfg:
            continue
        for name, arr in nsdfg.arrays.items():
            if (isinstance(arr, dace.data.Scalar)
                and arr.dtype == dace.float64
                and arr.transient):
                all_nested_scalars.append((nsdfg, name, arr))

    def is_excluded(name):
        return (name in _LOWERING_EXCLUDE
                or any(name.startswith(ex + "_") for ex in _LOWERING_EXCLUDE))

    # --- Apply: lower everything not excluded ---

    # 1. Non-transient arrays: boundary cast
    bc_targets = [n for n in all_non_transient if not is_excluded(n)]
    if bc_targets:
        print(f"Boundary cast: {len(bc_targets)} non-transient arrays")
        inject_cpu_boundary_cast(sdfg, bc_targets, external_dtype)

    # 2. Transient arrays: lower directly
    ta_targets = [n for n in all_transient_arrays if not is_excluded(n)]
    if ta_targets:
        print(f"Lowering {len(ta_targets)} transient arrays to {external_dtype}")
        for name in ta_targets:
            _propagate_dtype(sdfg, name, external_dtype)

    # Build set of all lowered array names (for scalar matching)
    _lowered_arrays = set(bc_targets) | set(ta_targets)

    def is_derived_from_lowered(name):
        """A scalar should only be lowered if its parent array was lowered."""
        for arr_name in _lowered_arrays:
            if name == arr_name or name.startswith(arr_name + "_"):
                return True
        return False

    # 3. Top-level transient scalars
    # WARN: using is_derived_from_lowered, not is_excluded — standalone scalars
    # (ztw1, zqadj, etc.) stay double. Lowering them causes severe errors.
    ts_count = 0
    for name in all_top_scalars:
        if is_derived_from_lowered(name):
            sdfg.arrays[name].dtype = external_dtype
            ts_count += 1
    if ts_count:
        print(f"Lowered {ts_count} top-level transient scalars to {external_dtype}")

    # 4. Nested transient scalars
    # Find comparison outputs first (should become int32, not float)
    comparison_scalars: set[str] = set()
    for nsdfg in sdfg.all_sdfgs_recursive():
        if nsdfg is sdfg:
            continue
        for state in nsdfg.states():
            for node in state.nodes():
                if not isinstance(node, nodes.Tasklet):
                    continue
                code = node.code.as_string
                if any(op in code for op in (" > ", " < ", " >= ", " <= ", " == ", " != ")):
                    for e in state.out_edges(node):
                        if (isinstance(e.dst, nodes.AccessNode)
                            and e.dst.data in nsdfg.arrays
                            and isinstance(nsdfg.arrays[e.dst.data], dace.data.Scalar)):
                            comparison_scalars.add(e.dst.data)

    # WARN: same as above — standalone nested scalars stay double.
    ns_count = 0
    for nsdfg, name, arr in all_nested_scalars:
        if is_derived_from_lowered(name):
            if name in comparison_scalars:
                arr.dtype = dace.int32
            else:
                arr.dtype = external_dtype
            ns_count += 1
    if ns_count:
        print(f"Lowered {ns_count} nested transient scalars to {external_dtype}")

    # 5. Floatify tasklet literals
    floatify_tasklets(sdfg)

    sdfg.validate()
