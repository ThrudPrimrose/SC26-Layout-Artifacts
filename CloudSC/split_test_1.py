#!/usr/bin/env python3
"""DaCe benchmark: 3D zsolqa[nclv,nclv,klon] vs split zsolqa_i_j[klon].
All dimensions symbolic. Extracted condensation loop nest from CLOUDSC."""
import dace
import numpy as np
import time

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
# Benchmark harness
# ================================================================
def xfill(shape, lo=0.0, hi=1.0):
    """Fast quasi-random fill via xorshift on flat index."""
    a = np.empty(shape, np.float64)
    n = a.size
    i = np.arange(1, n + 1, dtype=np.uint64)
    i ^= i << np.uint64(13)
    i ^= i >> np.uint64(7)
    i ^= i << np.uint64(17)
    a.ravel()[:] = lo + (hi - lo) * (i % np.uint64(10000)).astype(np.float64) / 10000.0
    return a


# ================================================================
# Benchmark harness
# ================================================================



# ~32 MB trash buffer — larger than typical L3 per core
_TRASH = np.empty(4 * 1024 * 1024, dtype=np.float64)

def flush_cache():
    """3-point Jacobi sweep on _TRASH to evict all cache lines."""
    n = _TRASH.shape[0]
    _TRASH[1:n-1] = 0.5 * (_TRASH[0:n-2] + _TRASH[2:n])
    return _TRASH[n >> 1]  # force materialization


def run():
    KLEV, KLON, NCLV = 96, 81920, 5
    QL, QI, QV = 1, 2, 5

    cst = dict(
        retv=0.6077, rtice=250.16, rtwat=273.16,
        rtwat_rtice_r=1.0 / (273.16 - 250.16),
        r5alvcp=5.4697e6, r4les=35.86,
        r5alscp=6.3147e6, r4ies=7.66,
        rthomo=250.0, rlmin=1e-12,
    )
    sym = dict(klev=KLEV, klon=KLON, nclv=NCLV, ncldql=QL, ncldqi=QI, ncldqv=QV)

    za     = xfill((KLEV, KLON), 0.0, 1.0)
    zdqs   = xfill((KLON,), -0.01, 0.01)
    zqsmix = xfill((KLEV, KLON), 0.001, 0.02)
    zqv    = xfill((KLEV, KLON), 0.001, 0.02)
    ztp1   = xfill((KLEV, KLON), 200.0, 300.0)

    rng = np.random.RandomState(42)

    print("Compiling 3D variant...", flush=True)
    csdfg_3d = condense_3d.compile(**sym)
    print("Compiling split variant...", flush=True)
    csdfg_sp = condense_split.compile(**sym)
    print("Done.\n", flush=True)

    REPS = 100

    # --- warmup both ---
    sq3 = np.zeros((NCLV, NCLV, KLON)); fg3 = np.zeros((NCLV, KLON))
    csdfg_3d(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
             zsolqa=sq3, zqxfg=fg3, **cst, **sym)
    svl = np.zeros(KLON); slv = np.zeros(KLON)
    svi = np.zeros(KLON); siv = np.zeros(KLON)
    fl  = np.zeros(KLON); fi  = np.zeros(KLON)
    csdfg_sp(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
             zsolqa_ncldqv_ncldql=svl, zsolqa_ncldql_ncldqv=slv,
             zsolqa_ncldqv_ncldqi=svi, zsolqa_ncldqi_ncldqv=siv,
             zqxfg_ncldql=fl, zqxfg_ncldqi=fi, **cst, **sym)

    # --- 3D variant ---
    t3 = []
    for _ in range(REPS):
        sq3[:] = 0; fg3[:] = 0
        trash_val = flush_cache()
        t0 = time.perf_counter()
        csdfg_3d(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
                 zsolqa=sq3, zqxfg=fg3, **cst, **sym)
        t3.append(time.perf_counter() - t0)
        print(f"\r  3D  rep {_+1:3d}/{REPS}  trash={trash_val:.4e}  "
              f"sq3[{rng.randint(KLON)}]={sq3[QV-1, QL-1, rng.randint(KLON)]:.6e}",
              end="", flush=True)
    print()
    ref_vl = sq3[QV - 1, QL - 1].copy()
    ref_iv = sq3[QI - 1, QV - 1].copy()

    # --- Split variant ---
    ts = []
    for _ in range(REPS):
        svl[:] = 0; slv[:] = 0; svi[:] = 0; siv[:] = 0; fl[:] = 0; fi[:] = 0
        trash_val = flush_cache()
        t0 = time.perf_counter()
        csdfg_sp(za=za, zdqs=zdqs, zqsmix=zqsmix, zqv=zqv, ztp1=ztp1,
                 zsolqa_ncldqv_ncldql=svl, zsolqa_ncldql_ncldqv=slv,
                 zsolqa_ncldqv_ncldqi=svi, zsolqa_ncldqi_ncldqv=siv,
                 zqxfg_ncldql=fl, zqxfg_ncldqi=fi, **cst, **sym)
        ts.append(time.perf_counter() - t0)
        print(f"\r  Split rep {_+1:3d}/{REPS}  trash={trash_val:.4e}  "
              f"svl[{rng.randint(KLON)}]={svl[rng.randint(KLON)]:.6e}",
              end="", flush=True)
    print()

    # --- Verify ---
    assert np.allclose(ref_vl, svl, atol=1e-12), f"MISMATCH vl: {np.max(np.abs(ref_vl - svl))}"
    assert np.allclose(ref_iv, siv, atol=1e-12), f"MISMATCH iv: {np.max(np.abs(ref_iv - siv))}"
    print("Verification PASSED\n")

    # --- Report ---
    m3 = np.median(t3); ms = np.median(ts)
    print(f"KLEV={KLEV}  KLON={KLON}  NCLV={NCLV}  REPS={REPS}")
    print(f"{'Variant':<22s} {'median(s)':>10s} {'min(s)':>10s} {'max(s)':>10s} {'vs 3D':>8s}")
    print("-" * 62)
    for nm, tt in [("3D zsolqa[5,5,klon]", t3), ("Split zsolqa_i_j", ts)]:
        print(f"{nm:<22s} {np.median(tt):10.4f} {np.min(tt):10.4f} {np.max(tt):10.4f} {m3/np.median(tt):8.3f}x")
    print(f"\n3D footprint : {NCLV*NCLV*KLON*8/1e6:.1f} MB")
    print(f"Split touched: {4*KLON*8/1e6:.1f} MB (4 slices)")


if __name__ == "__main__":
    run()