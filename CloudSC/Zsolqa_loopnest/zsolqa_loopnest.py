#!/usr/bin/env python3
"""DaCe benchmark: 3D zsolqa[nclv,nclv,klon] vs split zsolqa_i_j[klon].
Run with --variant=3d|split|3d-gpu|split-gpu. Logs per-iteration timings to CSV."""
import argparse
import csv
import dace
import numpy as np
import time

from dace.sdfg.state import LoopRegion

"""
export LD_LIBRARY_PATH=$(spack location -i cuda@12.9)/lib64:$LD_LIBRARY_PATH
spack load cuda@12.9
pyenv activate dace_py_12
spack load gcc/76jw6nu
"""

VALID_VARIANTS = ['3d', 'split', '3d-gpu', 'split-gpu']

def parse_variant(variant):
    """Return (base, is_gpu) from variant string."""
    if variant.endswith('-gpu'):
        return variant[:-4], True
    return variant, False

# ---------- symbols ----------
klev = dace.symbol('klev', dtype=dace.int32)
klon = dace.symbol('klon', dtype=dace.int32)
nclv = dace.symbol('nclv', dtype=dace.int32)
ncldql = dace.symbol('ncldql', dtype=dace.int32)
ncldqi = dace.symbol('ncldqi', dtype=dace.int32)
ncldqv = dace.symbol('ncldqv', dtype=dace.int32)

# ================================================================
# Variant 1: original 3D zsolqa[nclv, nclv, klon]
# ================================================================
@dace.program
def condense_3d(
    za: dace.float64[klev, klon],
    zdqs: dace.float64[klon],
    zqsmix: dace.float64[klev, klon],
    zqv: dace.float64[klev, klon],
    ztp1: dace.float64[klev, klon],
    zsolqa: dace.float64[nclv, nclv, klon],
    zqxfg: dace.float64[nclv, klon],
    retv: dace.float64, rtice: dace.float64, rtwat: dace.float64,
    rtwat_rtice_r: dace.float64, r5alvcp: dace.float64, r4les: dace.float64,
    r5alscp: dace.float64, r4ies: dace.float64, rthomo: dace.float64,
    rlmin: dace.float64,
):
    for jk in range(klev):
        for jl in range(klon):
            if za[jk, jl] > 1e-14:
                if zdqs[jl] <= -rlmin:
                    lc = max(-zdqs[jl], 0.0)
                    af = min(1.0, ((max(rtice, min(rtwat, ztp1[jk, jl])) - rtice) * rtwat_rtice_r) ** 2)
                    zcor = 1.0 / (1.0 - retv * zqsmix[jk, jl])
                    cdm_full = (zqv[jk, jl] - zqsmix[jk, jl]) / (1.0 + zcor * zqsmix[jk, jl] * (
                        af * r5alvcp / (ztp1[jk, jl] - r4les) ** 2
                        + (1.0 - af) * r5alscp / (ztp1[jk, jl] - r4ies) ** 2))
                    cdm_part = (zqv[jk, jl] - za[jk, jl] * zqsmix[jk, jl]) / za[jk, jl]
                    if za[jk, jl] > 0.99:
                        cdm = cdm_full
                    else:
                        cdm = cdm_part
                    lc = za[jk, jl] * max(min(lc, cdm), 0.0)
                    if lc >= rlmin:
                        if ztp1[jk, jl] > rthomo:
                            zsolqa[ncldqv - 1, ncldql - 1, jl] = zsolqa[ncldqv - 1, ncldql - 1, jl] + lc
                            zsolqa[ncldql - 1, ncldqv - 1, jl] = zsolqa[ncldql - 1, ncldqv - 1, jl] - lc
                            zqxfg[ncldql - 1, jl] = zqxfg[ncldql - 1, jl] + lc
                        else:
                            zsolqa[ncldqv - 1, ncldqi - 1, jl] = zsolqa[ncldqv - 1, ncldqi - 1, jl] + lc
                            zsolqa[ncldqi - 1, ncldqv - 1, jl] = zsolqa[ncldqi - 1, ncldqv - 1, jl] - lc
                            zqxfg[ncldqi - 1, jl] = zqxfg[ncldqi - 1, jl] + lc


# ================================================================
# Variant 2: split zsolqa_{src}_{dst}[klon]
# ================================================================
@dace.program
def condense_split(
    za: dace.float64[klev, klon],
    zdqs: dace.float64[klon],
    zqsmix: dace.float64[klev, klon],
    zqv: dace.float64[klev, klon],
    ztp1: dace.float64[klev, klon],
    zsolqa_ncldqv_ncldql: dace.float64[klon],
    zsolqa_ncldql_ncldqv: dace.float64[klon],
    zsolqa_ncldqv_ncldqi: dace.float64[klon],
    zsolqa_ncldqi_ncldqv: dace.float64[klon],
    zqxfg_ncldql: dace.float64[klon],
    zqxfg_ncldqi: dace.float64[klon],
    retv: dace.float64, rtice: dace.float64, rtwat: dace.float64,
    rtwat_rtice_r: dace.float64, r5alvcp: dace.float64, r4les: dace.float64,
    r5alscp: dace.float64, r4ies: dace.float64, rthomo: dace.float64,
    rlmin: dace.float64,
):
    for jk in range(klev):
        for jl in range(klon):
            if za[jk, jl] > 1e-14:
                if zdqs[jl] <= -rlmin:
                    lc = max(-zdqs[jl], 0.0)
                    af = min(1.0, ((max(rtice, min(rtwat, ztp1[jk, jl])) - rtice) * rtwat_rtice_r) ** 2)
                    zcor = 1.0 / (1.0 - retv * zqsmix[jk, jl])
                    cdm_full = (zqv[jk, jl] - zqsmix[jk, jl]) / (1.0 + zcor * zqsmix[jk, jl] * (
                        af * r5alvcp / (ztp1[jk, jl] - r4les) ** 2
                        + (1.0 - af) * r5alscp / (ztp1[jk, jl] - r4ies) ** 2))
                    cdm_part = (zqv[jk, jl] - za[jk, jl] * zqsmix[jk, jl]) / za[jk, jl]
                    if za[jk, jl] > 0.99:
                        cdm = cdm_full
                    else:
                        cdm = cdm_part
                    lc = za[jk, jl] * max(min(lc, cdm), 0.0)
                    if lc >= rlmin:
                        if ztp1[jk, jl] > rthomo:
                            zsolqa_ncldqv_ncldql[jl] = zsolqa_ncldqv_ncldql[jl] + lc
                            zsolqa_ncldql_ncldqv[jl] = zsolqa_ncldql_ncldqv[jl] - lc
                            zqxfg_ncldql[jl] = zqxfg_ncldql[jl] + lc
                        else:
                            zsolqa_ncldqv_ncldqi[jl] = zsolqa_ncldqv_ncldqi[jl] + lc
                            zsolqa_ncldqi_ncldqv[jl] = zsolqa_ncldqi_ncldqv[jl] - lc
                            zqxfg_ncldqi[jl] = zqxfg_ncldqi[jl] + lc


# ================================================================
# Helpers
# ================================================================
def xfill(shape, lo=0.0, hi=1.0):
    a = np.empty(shape, np.float64)
    n = a.size
    i = np.arange(1, n + 1, dtype=np.uint64)
    i ^= i << np.uint64(13)
    i ^= i >> np.uint64(7)
    i ^= i << np.uint64(17)
    a.ravel()[:] = lo + (hi - lo) * (i % np.uint64(10000)).astype(np.float64) / 10000.0
    return a


_TRASH = np.empty(128 * 1024 * 1024, dtype=np.float64)

def flush_cache():
    n = _TRASH.shape[0]
    _TRASH[1:n-1] = 0.5 * (_TRASH[0:n-2] + _TRASH[2:n])
    return _TRASH[n >> 1]


def build_sdfg(variant, gpu=False):
    """Build and optimize the SDFG for the given variant."""
    base, _ = parse_variant(variant) if not gpu else (variant, True)
    if base == '3d':
        sdfg = condense_3d.to_sdfg()
    else:
        sdfg = condense_split.to_sdfg()

    from dace.transformation.passes.clean_access_node_to_scalar_slice_to_tasklet_pattern import CleanAccessNodeToScalarSliceToTaskletPattern
    CleanAccessNodeToScalarSliceToTaskletPattern(permissive=True).apply_pass(sdfg, {})

    symbol_map = {
        "nclv": 5, "ncldql": 1, "ncldqi": 2,
        "ncldqr": 3, "ncldqs": 4, "ncldqv": 5,
    }
    sdfg.replace_dict(symbol_map)

    from dace.transformation.interstate.loop_to_map import LoopToMap
    sdfg.apply_transformations_repeated([LoopToMap])
    sdfg.simplify()

    from dace.transformation.passes.fusion_inline import InlineSDFGs
    from dace.transformation.dataflow.map_collapse import MapCollapse
    InlineSDFGs().apply_pass(sdfg, {})
    sdfg.apply_transformations_repeated(MapCollapse)

    if gpu:
        from dace.transformation.auto.auto_optimize import auto_optimize, apply_gpu_storage
        #Interchanges klev (loop) and klon (map), in original klev is above
        auto_optimize(sdfg, device=dace.dtypes.DeviceType.GPU, use_gpu_storage=True)
        #apply_gpu_storage(sdfg)
        #sdfg.apply_gpu_transformations()
        dace.config.Config.set('compiler', 'cuda', 'default_block_size', value="256,1,1")

    sdfg.save(f"condense_{variant}{'_gpu' if gpu else ''}.sdfg")
    return sdfg


def run_variant(variant, reps, csv_path):
    """Run a single variant for `reps` iterations, log to CSV."""
    base, is_gpu = parse_variant(variant)
    KLEV, KLON, NCLV = 96, 20480, 5

    cst = dict(
        retv=0.6077, rtice=250.16, rtwat=273.16,
        rtwat_rtice_r=1.0 / (273.16 - 250.16),
        r5alvcp=5.4697e6, r4les=35.86,
        r5alscp=6.3147e6, r4ies=7.66,
        rthomo=250.0, rlmin=1e-12,
    )
    sym = dict(klev=KLEV, klon=KLON, nclv=NCLV, ncldql=1, ncldqi=2, ncldqv=5)

    # Choose array library: cupy for GPU, numpy for CPU
    if is_gpu:
        import cupy as xp
        print("Using cupy for GPU array allocation", flush=True)
    else:
        xp = np

    # Generate data on CPU then transfer (xfill uses numpy internals)
    za_cpu     = xfill((KLEV, KLON), 0.0, 1.0)
    zdqs_cpu   = xfill((KLON,), -0.01, 0.01)
    zqsmix_cpu = xfill((KLEV, KLON), 0.001, 0.02)
    zqv_cpu    = xfill((KLEV, KLON), 0.001, 0.02)
    ztp1_cpu   = xfill((KLEV, KLON), 200.0, 300.0)

    za     = xp.asarray(za_cpu)
    zdqs   = xp.asarray(zdqs_cpu)
    zqsmix = xp.asarray(zqsmix_cpu)
    zqv    = xp.asarray(zqv_cpu)
    ztp1   = xp.asarray(ztp1_cpu)

    print(f"Building {variant} SDFG...", flush=True)
    sdfg = build_sdfg(base, gpu=is_gpu)
    if not is_gpu:
        sdfg.instrument = dace.dtypes.InstrumentationType.Timer
    else:
        sdfg.instrument = dace.dtypes.InstrumentationType.Timer
        #for state in sdfg.all_states():
        #    for node in state.nodes():
        #        if isinstance(node, dace.nodes.MapEntry) and node.map.schedule == dace.ScheduleType.GPU_Device:
        #            node.map.instrument = dace.dtypes.InstrumentationType.GPU_Events

    print(f"Compiling {variant}...", flush=True)
    csdfg = sdfg.compile()
    print("Done.\n", flush=True)

    if base == '3d':
        sq3 = xp.zeros((NCLV, NCLV, KLON), dtype=xp.float64)
        fg3 = xp.zeros((NCLV, KLON), dtype=xp.float64)
        def call():
            sq3[:] = 0; fg3[:] = 0
            csdfg(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
                  zsolqa=sq3, zqxfg=fg3, **cst, **sym)
    else:
        svl = xp.zeros(KLON, dtype=xp.float64); slv = xp.zeros(KLON, dtype=xp.float64)
        svi = xp.zeros(KLON, dtype=xp.float64); siv = xp.zeros(KLON, dtype=xp.float64)
        fl  = xp.zeros(KLON, dtype=xp.float64); fi  = xp.zeros(KLON, dtype=xp.float64)
        def call():
            svl[:] = 0; slv[:] = 0; svi[:] = 0; siv[:] = 0; fl[:] = 0; fi[:] = 0
            csdfg(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
                  zsolqa_ncldqv_ncldql=svl, zsolqa_ncldql_ncldqv=slv,
                  zsolqa_ncldqv_ncldqi=svi, zsolqa_ncldqi_ncldqv=siv,
                  zqxfg_ncldql=fl, zqxfg_ncldqi=fi, **cst, **sym)

    # Warmup
    print("Warmup (3 reps)...", flush=True)
    for _ in range(3):
        if not is_gpu:
            flush_cache()
        call()
    if is_gpu:
        import cupy
        cupy.cuda.Device().synchronize()


    # Timed reps
    print(f"\n{'':=<60}")
    print(f"  {variant.upper()} variant  ({reps} iterations)")
    print(f"{'':=<60}")
    print(f"  {'Iter':>4s}  {'Time (ms)':>10s}  {'Running med (ms)':>16s}")
    print(f"  {'----':>4s}  {'----------':>10s}  {'----------------':>16s}")

    times = []
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['variant', 'rep', 'time_ms'])
        for rep in range(reps):
            if not is_gpu:
                flush_cache()
                flush_cache()
            call()
            report = sdfg.get_latest_report()
            elapsed_ms = report.events[-1].duration * 1e-3
            times.append(elapsed_ms)
            writer.writerow([variant, rep + 1, f"{elapsed_ms:.6f}"])
            running_med = np.median(times)
            print(f"  {rep+1:4d}  {elapsed_ms:10.3f}  {running_med:16.3f}")

    med = np.median(times)
    print(f"\n  Summary: median={med:.3f} ms  min={np.min(times):.3f} ms  "
          f"max={np.max(times):.3f} ms  std={np.std(times):.3f} ms")
    print(f"  Results written to: {csv_path}\n")


def main():
    parser = argparse.ArgumentParser(description="Condensation kernel benchmark")
    parser.add_argument('--variant', required=True, choices=VALID_VARIANTS,
                        help="Which variant to run: '3d', 'split', '3d-gpu', or 'split-gpu'")
    parser.add_argument('--reps', type=int, default=100,
                        help="Number of timed repetitions (default: 100)")
    parser.add_argument('--csv', type=str, default=None,
                        help="Output CSV path (default: bench_{variant}.csv)")
    args = parser.parse_args()

    csv_path = args.csv or f"bench_{args.variant}.csv"
    run_variant(args.variant, args.reps, csv_path)


if __name__ == "__main__":
    main()