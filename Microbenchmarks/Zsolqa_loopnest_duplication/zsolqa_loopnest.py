#!/usr/bin/env python3
"""DaCe benchmark: 3D zsolqa[nclv,nclv,klon] vs split zsolqa_i_j[klon],
plus jk-expanded "reduce" variants.

Variants:
  3d           - original 3D zsolqa[nclv, nclv, klon]
  split        - split zsolqa_{src}_{dst}[klon]
  3d-reduce    - kernel writes zsolqa_buf[nclv, nclv, klev, klon],
                 then reduces over klev into zsolqa[nclv, nclv, klon]
  split-reduce - kernel writes zsolqa_buf_{src}_{dst}[klev, klon],
                 then reduces over klev into zsolqa_{src}_{dst}[klon]
  (all with optional -gpu suffix)

Run with --variant=<n>. Logs per-iteration timings to CSV.
"""
import argparse
import csv
import dace
import numpy as np
import time

from dace.sdfg.state import LoopRegion

VALID_VARIANTS = [
    '3d', 'split', '3d-reduce', 'split-reduce',
    '3d-gpu', 'split-gpu', '3d-reduce-gpu', 'split-reduce-gpu',
]

def parse_variant(variant):
    """Return (base, is_gpu) from variant string."""
    is_gpu = variant.endswith('-gpu')
    base = variant[:-4] if is_gpu else variant
    return base, is_gpu

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
# Variant 3: 3D-reduce
#   Kernel writes to expanded buffers: zsolqa_buf[nclv, nclv, klev, klon]
#                                      zqxfg_buf[nclv, klev, klon]
#   Caller reduces over klev into:     zsolqa[nclv, nclv, klon]
#                                      zqxfg[nclv, klon]
# ================================================================
@dace.program
def condense_3d_reduce(
    za: dace.float64[klev, klon],
    zdqs: dace.float64[klon],
    zqsmix: dace.float64[klev, klon],
    zqv: dace.float64[klev, klon],
    ztp1: dace.float64[klev, klon],
    zsolqa_buf: dace.float64[nclv, nclv, klev, klon],
    zqxfg_buf: dace.float64[nclv, klev, klon],
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
                            zsolqa_buf[ncldqv - 1, ncldql - 1, jk, jl] = lc
                            zsolqa_buf[ncldql - 1, ncldqv - 1, jk, jl] = -lc
                            zqxfg_buf[ncldql - 1, jk, jl] = lc
                        else:
                            zsolqa_buf[ncldqv - 1, ncldqi - 1, jk, jl] = lc
                            zsolqa_buf[ncldqi - 1, ncldqv - 1, jk, jl] = -lc
                            zqxfg_buf[ncldqi - 1, jk, jl] = lc
    for jl in range(klon):
        sum_qv_ql = 0.0
        sum_ql_qv = 0.0
        sum_qv_qi = 0.0
        sum_qi_qv = 0.0
        sum_fgl = 0.0
        sum_fgi = 0.0
        for jk in range(klev):
            sum_qv_ql += zsolqa_buf[ncldqv - 1, ncldql - 1, jk, jl]
            sum_ql_qv += zsolqa_buf[ncldql - 1, ncldqv - 1, jk, jl]
            sum_qv_qi += zsolqa_buf[ncldqv - 1, ncldqi - 1, jk, jl]
            sum_qi_qv += zsolqa_buf[ncldqi - 1, ncldqv - 1, jk, jl]
            sum_fgl += zqxfg_buf[ncldql - 1, jk, jl]
            sum_fgi += zqxfg_buf[ncldqi - 1, jk, jl]
        zsolqa[ncldqv - 1, ncldql - 1, jl] = sum_qv_ql
        zsolqa[ncldql - 1, ncldqv - 1, jl] = sum_ql_qv
        zsolqa[ncldqv - 1, ncldqi - 1, jl] = sum_qv_qi
        zsolqa[ncldqi - 1, ncldqv - 1, jl] = sum_qi_qv
        zqxfg[ncldql - 1, jl] = sum_fgl
        zqxfg[ncldqi - 1, jl] = sum_fgi

# ================================================================
# Variant 4: split-reduce
#   Kernel writes to expanded buffers: zsolqa_*_buf[klev, klon]
#                                      zqxfg_*_buf[klev, klon]
#   Caller reduces over klev into:     zsolqa_*[klon]
#                                      zqxfg_*[klon]
# ================================================================
@dace.program
def condense_split_reduce(
    za: dace.float64[klev, klon],
    zdqs: dace.float64[klon],
    zqsmix: dace.float64[klev, klon],
    zqv: dace.float64[klev, klon],
    ztp1: dace.float64[klev, klon],
    zsolqa_ncldqv_ncldql_buf: dace.float64[klev, klon],
    zsolqa_ncldql_ncldqv_buf: dace.float64[klev, klon],
    zsolqa_ncldqv_ncldqi_buf: dace.float64[klev, klon],
    zsolqa_ncldqi_ncldqv_buf: dace.float64[klev, klon],
    zqxfg_ncldql_buf: dace.float64[klev, klon],
    zqxfg_ncldqi_buf: dace.float64[klev, klon],
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
                            zsolqa_ncldqv_ncldql_buf[jk, jl] = lc
                            zsolqa_ncldql_ncldqv_buf[jk, jl] = - lc
                            zqxfg_ncldql_buf[jk, jl] = lc
                        else:
                            zsolqa_ncldqv_ncldqi_buf[jk, jl] = lc
                            zsolqa_ncldqi_ncldqv_buf[jk, jl] = - lc
                            zqxfg_ncldqi_buf[jk, jl] = lc

    for jl in range(klon):
        sum_ncldqv_ncldql = 0.0
        sum_ncldql_ncldqv = 0.0
        sum_ncldqv_ncldqi = 0.0
        sum_ncldqi_ncldqv = 0.0
        sum_qxfg_ncldql = 0.0
        sum_qxfg_ncldqi = 0.0
        for jk in range(klev):
            sum_ncldqv_ncldql += zsolqa_ncldqv_ncldql_buf[jk, jl]
            sum_ncldql_ncldqv += zsolqa_ncldql_ncldqv_buf[jk, jl]
            sum_ncldqv_ncldqi += zsolqa_ncldqv_ncldqi_buf[jk, jl]
            sum_ncldqi_ncldqv += zsolqa_ncldqi_ncldqv_buf[jk, jl]
            sum_qxfg_ncldql += zqxfg_ncldql_buf[jk, jl]
            sum_qxfg_ncldqi += zqxfg_ncldqi_buf[jk, jl]
        zsolqa_ncldqv_ncldql[jl] = sum_ncldqv_ncldql
        zsolqa_ncldql_ncldqv[jl] = sum_ncldql_ncldqv
        zsolqa_ncldqv_ncldqi[jl] = sum_ncldqv_ncldqi
        zsolqa_ncldqi_ncldqv[jl] = sum_ncldqi_ncldqv
        zqxfg_ncldql[jl] = sum_qxfg_ncldql
        zqxfg_ncldqi[jl] = sum_qxfg_ncldqi

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


FLUSH_SIZE = 128 * 1024 * 1024  # 128M doubles = 1 GiB

_TRASH_CPU = None
_TRASH_GPU = None

def _ensure_trash_cpu():
    global _TRASH_CPU
    if _TRASH_CPU is None:
        _TRASH_CPU = np.empty(FLUSH_SIZE, dtype=np.float64)
        _TRASH_CPU[:] = np.arange(FLUSH_SIZE, dtype=np.float64) * 1e-15
    return _TRASH_CPU

def _ensure_trash_gpu():
    global _TRASH_GPU
    if _TRASH_GPU is None:
        import cupy
        _TRASH_GPU = cupy.empty(FLUSH_SIZE, dtype=cupy.float64)
        _TRASH_GPU[:] = cupy.arange(FLUSH_SIZE, dtype=cupy.float64) * 1e-15
    return _TRASH_GPU

def flush_cache_cpu():
    t = _ensure_trash_cpu()
    n = t.shape[0]
    t[1:n-1] = 0.5 * (t[0:n-2] + t[2:n])
    h = float(np.sum(t[::4096]))
    return h

def flush_cache_gpu():
    import cupy
    t = _ensure_trash_gpu()
    n = t.shape[0]
    t[1:n-1] = 0.5 * (t[0:n-2] + t[2:n])
    cupy.cuda.Device().synchronize()
    h = float(cupy.sum(t[::4096]))
    return h


def build_sdfg(base, gpu=False):
    """Build and optimize the SDFG for the given variant."""
    if base == '3d':
        sdfg = condense_3d.to_sdfg()
    elif base == 'split':
        sdfg = condense_split.to_sdfg()
    elif base == '3d-reduce':
        sdfg = condense_3d_reduce.to_sdfg()
    elif base == 'split-reduce':
        sdfg = condense_split_reduce.to_sdfg()
    else:
        raise ValueError(f"Unknown base variant: {base}")

    if gpu is True:
        sdfg.name += "_gpu"
    else:
        sdfg.name += "_cpu"

    from dace.transformation.passes.clean_access_node_to_scalar_slice_to_tasklet_pattern import (
        CleanAccessNodeToScalarSliceToTaskletPattern,
    )
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
    InlineSDFGs().apply_pass(sdfg, {})
    sdfg.apply_transformations_repeated(MapCollapse)

    if gpu:
        from dace.transformation.auto.auto_optimize import auto_optimize
        auto_optimize(sdfg, device=dace.dtypes.DeviceType.GPU, use_gpu_storage=True)
        dace.config.Config.set('compiler', 'cuda', 'default_block_size', value="256,1,1")

    tag = f"{base}{'_gpu' if gpu else ''}"
    sdfg.save(f"condense_{tag}.sdfg")
    return sdfg


def run_variant(variant, reps, csv_path):
    """Run a single variant for `reps` iterations, log to CSV."""
    base, is_gpu = parse_variant(variant)
    KLEV, KLON, NCLV = 96, 20480*16, 5

    cst = dict(
        retv=0.6077, rtice=250.16, rtwat=273.16,
        rtwat_rtice_r=1.0 / (273.16 - 250.16),
        r5alvcp=5.4697e6, r4les=35.86,
        r5alscp=6.3147e6, r4ies=7.66,
        rthomo=250.0, rlmin=1e-12,
    )
    sym = dict(klev=KLEV, klon=KLON, nclv=NCLV, ncldql=1, ncldqi=2, ncldqv=5)

    # Choose array library
    if is_gpu:
        import cupy as xp
        print("Using cupy for GPU array allocation", flush=True)
    else:
        xp = np

    # Choose flush function
    flush_fn = flush_cache_gpu if is_gpu else flush_cache_cpu

    # Allocate flush-cache buffer and print hash
    h = flush_fn()
    print(f"Flush-cache hash: {h:.10e}", flush=True)

    # Generate data on CPU then transfer
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

    # ---- Allocate ALL buffers up front (passed from outside the kernel) ----

    if base == '3d':
        zsolqa = xp.zeros((NCLV, NCLV, KLON), dtype=xp.float64)
        zqxfg  = xp.zeros((NCLV, KLON), dtype=xp.float64)

        def call():
            zsolqa[:] = 0; zqxfg[:] = 0
            csdfg(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
                  zsolqa=zsolqa, zqxfg=zqxfg, **cst, **sym)

    elif base == 'split':
        svl = xp.zeros(KLON, dtype=xp.float64)
        slv = xp.zeros(KLON, dtype=xp.float64)
        svi = xp.zeros(KLON, dtype=xp.float64)
        siv = xp.zeros(KLON, dtype=xp.float64)
        fl  = xp.zeros(KLON, dtype=xp.float64)
        fi  = xp.zeros(KLON, dtype=xp.float64)

        def call():
            svl[:] = 0; slv[:] = 0; svi[:] = 0; siv[:] = 0
            fl[:] = 0; fi[:] = 0
            csdfg(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
                  zsolqa_ncldqv_ncldql=svl, zsolqa_ncldql_ncldqv=slv,
                  zsolqa_ncldqv_ncldqi=svi, zsolqa_ncldqi_ncldqv=siv,
                  zqxfg_ncldql=fl, zqxfg_ncldqi=fi, **cst, **sym)

    elif base == '3d-reduce':
        # Intermediate jk-expanded buffers (kernel writes here)
        zsolqa_buf = xp.zeros((NCLV, NCLV, KLEV, KLON), dtype=xp.float64)
        zqxfg_buf  = xp.zeros((NCLV, KLEV, KLON), dtype=xp.float64)
        # Final reduced outputs (same shape as '3d')
        zsolqa = xp.zeros((NCLV, NCLV, KLON), dtype=xp.float64)
        zqxfg  = xp.zeros((NCLV, KLON), dtype=xp.float64)

        def call():
            zsolqa_buf[:] = 0; zqxfg_buf[:] = 0
            zsolqa[:] = 0; zqxfg[:] = 0
            csdfg(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
                  zsolqa_buf=zsolqa_buf, zqxfg_buf=zqxfg_buf,
                  zsolqa=zsolqa, zqxfg=zqxfg, **cst, **sym)

    elif base == 'split-reduce':
        # Intermediate jk-expanded buffers (kernel writes here)
        svl_buf = xp.zeros((KLEV, KLON), dtype=xp.float64)
        slv_buf = xp.zeros((KLEV, KLON), dtype=xp.float64)
        svi_buf = xp.zeros((KLEV, KLON), dtype=xp.float64)
        siv_buf = xp.zeros((KLEV, KLON), dtype=xp.float64)
        fl_buf  = xp.zeros((KLEV, KLON), dtype=xp.float64)
        fi_buf  = xp.zeros((KLEV, KLON), dtype=xp.float64)
        # Final reduced outputs (same shape as 'split')
        svl = xp.zeros(KLON, dtype=xp.float64)
        slv = xp.zeros(KLON, dtype=xp.float64)
        svi = xp.zeros(KLON, dtype=xp.float64)
        siv = xp.zeros(KLON, dtype=xp.float64)
        fl  = xp.zeros(KLON, dtype=xp.float64)
        fi  = xp.zeros(KLON, dtype=xp.float64)

        def call():
            svl_buf[:] = 0; slv_buf[:] = 0; svi_buf[:] = 0; siv_buf[:] = 0
            fl_buf[:] = 0; fi_buf[:] = 0
            svl[:] = 0; slv[:] = 0; svi[:] = 0; siv[:] = 0
            fl[:] = 0; fi[:] = 0
            csdfg(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
                  zsolqa_ncldqv_ncldql_buf=svl_buf,
                  zsolqa_ncldql_ncldqv_buf=slv_buf,
                  zsolqa_ncldqv_ncldqi_buf=svi_buf,
                  zsolqa_ncldqi_ncldqv_buf=siv_buf,
                  zqxfg_ncldql_buf=fl_buf,
                  zqxfg_ncldqi_buf=fi_buf,
                  zsolqa_ncldqv_ncldql=svl,
                  zsolqa_ncldql_ncldqv=slv,
                  zsolqa_ncldqv_ncldqi=svi,
                  zsolqa_ncldqi_ncldqv=siv,
                  zqxfg_ncldql=fl,
                  zqxfg_ncldqi=fi, **cst, **sym)

    else:
        raise ValueError(f"Unknown base variant: {base}")

    # ---- Build & compile ----
    print(f"Building {variant} SDFG...", flush=True)
    sdfg = build_sdfg(base, gpu=is_gpu)
    states = {s for s in sdfg.all_states()}
    assert len(states) == 1, "Expected exactly one state in the SDFG"
    state = next(iter(states))
    
    if is_gpu:
        state.instrument = dace.dtypes.InstrumentationType.Timer
    else:
        state.instrument = dace.dtypes.InstrumentationType.GPU_Events

    print(f"Compiling {variant}...", flush=True)
    csdfg = sdfg.compile()
    print("Done.\n", flush=True)

    # ---- Warmup ----
    print("Warmup (3 reps)...", flush=True)
    for _ in range(3):
        flush_fn()
        call()
    if is_gpu:
        import cupy
        cupy.cuda.Device().synchronize()

    # ---- Timed reps ----
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
            flush_fn()
            call()
            flush_fn()
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
                        help="Which variant to run: " + ", ".join(VALID_VARIANTS))
    parser.add_argument('--reps', type=int, default=100,
                        help="Number of timed repetitions (default: 100)")
    parser.add_argument('--csv', type=str, default=None,
                        help="Output CSV path (default: bench_{variant}.csv)")
    args = parser.parse_args()

    csv_path = args.csv or f"bench_{args.variant}.csv"
    run_variant(args.variant, args.reps, csv_path)


if __name__ == "__main__":
    main()