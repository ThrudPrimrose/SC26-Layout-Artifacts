import dace
from dace.sdfg.state import LoopRegion, ReturnBlock
from dace.transformation.interstate.loop_to_map import LoopToMap
from dace.transformation.passes.analysis import loop_analysis
from dace.transformation.passes.clean_access_node_to_scalar_slice_to_tasklet_pattern import CleanAccessNodeToScalarSliceToTaskletPattern
from dace.transformation.passes.ssa_loop_iterators import SSALoopIterators
import copy
from dace.sdfg.state import LoopRegion
import dace.sdfg.utils as sdutil
import dace.sdfg.construction_utils as cutil
from dace.transformation.interstate.loop_to_map import LoopToMap
from dace.transformation.passes.split_tasklets import SplitTasklets
from dace.transformation.passes.analysis import loop_analysis
from dace.transformation.passes.scalar_fission import ScalarFission
from dace.transformation.passes import RemoveUnusedSymbols, analysis as ap

sdfg1 = dace.SDFG.from_file("cloudsc_split.sdfgz")
sdfg2 = dace.SDFG.from_file("cloudsc_unrolled.sdfgz")


for name, sdfg in [("cloudsc_unrolled", sdfg2),]: #[("cloudsc_split", sdfg1), ("cloudsc_unrolled", sdfg2)]:
    for n in sdfg.nodes():
        if isinstance(n, ReturnBlock):
            ies = sdfg.in_edges(n)
            oes = sdfg.out_edges(n)
            sdfg.remove_node(n)
            s = sdfg.add_state(n.label + "_hull")
            for e in ies:
                sdfg.add_edge(e.src, e.src_conn, s, None, copy.deepcopy(e.data))
            for e in oes:
                sdfg.add_edge(s, None, e.dst, e.dst_conn, copy.deepcopy(e.data))
    
    # Breaks
    #pipeline_results = dict()
    #pipeline_results[ap.ControlFlowBlockReachability.__name__] = ap.ControlFlowBlockReachability().apply_pass(sdfg, pipeline_results)
    #pipeline_results[ap.FindAccessNodes.__name__] = ap.FindAccessNodes().apply_pass(sdfg, pipeline_results)
    #pipeline_results[ap.AccessSets.__name__] = ap.AccessSets().apply_pass(sdfg, pipeline_results)
    #pipeline_results[ap.ScalarWriteShadowScopes.__name__] = ap.ScalarWriteShadowScopes().apply_pass(sdfg, pipeline_results)
    #ScalarFission().apply_pass(sdfg, pipeline_results)
    #sdfg.validate()

    #SSALoopIterators().apply_pass(sdfg, None)
    #sdfg.validate()
    #CleanAccessNodeToScalarSliceToTaskletPattern(permissive=True).apply_pass(sdfg, None)
    #sdfg.validate()

    sdfg.simplify()
    sdfg.name += "_simplified"
    sdfg.save(f"{sdfg.name}.sdfgz", compress=True)
    exit(0)


    # Collect loops before
    loops_before = set()
    for node, _ in sdfg.all_nodes_recursive():
        if isinstance(node, LoopRegion):
            loop_range = f"for ({loop_analysis.get_init_assignment(node)}; {loop_analysis.get_loop_end(node)}; {loop_analysis.get_update_assignment(node)})"
            loops_before.add((id(node), node.label, node.loop_variable, loop_range))

    sdfg.apply_transformations_repeated(LoopToMap)

    # Collect loops after
    loops_after = set()
    for node, _ in sdfg.all_nodes_recursive():
        if isinstance(node, LoopRegion):
            loop_range = f"for ({loop_analysis.get_init_assignment(node)}; {loop_analysis.get_loop_end(node)}; {loop_analysis.get_update_assignment(node)})"
            loops_after.add((node.label, node.loop_variable, loop_range))

    # Remaining loops = not converted
    labels_before = {(label, var, rng) for _, label, var, rng in loops_before}
    converted = labels_before - loops_after
    remaining = loops_after  # these still exist as loops

    print(f"\n=== {name} ===")
    print(f"Loops before: {len(loops_before)}")
    print(f"Converted to maps: {len(converted)}")
    print(f"Remaining loops ({len(remaining)}):")
    for label, var, rng in sorted(remaining, key=lambda x: x[0]):
        print(f"  {label}  (var={var}, range={rng})")

    sdfg.name += "_l2m"
    sdfg.save(f"{sdfg.name}.sdfgz", compress=True)