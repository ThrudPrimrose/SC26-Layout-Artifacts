#!/usr/bin/env python3
"""
N-Body Simulation: AoS / SoA / AoSoA — local vs buffer accel variants
======================================================================
Based on LRZ COW 2022, Jonathan Coles.

Each layout has TWO accel variants:
  _local:   scalar axi, ayi, azi accumulators.  j is sequential.  O(N) memory.
  _buf:     fx[N,N], fy[N,N], fz[N,N] passed from outside.  i,j both parallel.
            Separate map phase (fill buffer) and reduce phase (sum over j).

All kernels are plain functions.  No classes.  No semicolons.
"""

import numpy as np
import time
import argparse
import math

G      = 1.0
RSOFT2 = 0.01 * 0.01


# ═══════════════════════════════════════════════════════════════════════════
#  Initial conditions
# ═══════════════════════════════════════════════════════════════════════════

def init_circular(N, seed=42):
    rng = np.random.default_rng(seed)
    m = np.ones(N, dtype=np.float64) / N
    r = np.zeros((N, 3), dtype=np.float64)
    v = np.zeros((N, 3), dtype=np.float64)
    for i in range(N):
        radius = 0.1 + 0.9 * rng.random()
        theta  = 2.0 * np.pi * rng.random()
        r[i, 0] = radius * np.cos(theta)
        r[i, 1] = radius * np.sin(theta)
        r[i, 2] = 0.05 * (rng.random() - 0.5)
        vc = np.sqrt(G / radius) * 0.5
        v[i, 0] = -vc * np.sin(theta)
        v[i, 1] =  vc * np.cos(theta)
    return m, r, v


def init_random(N, seed=42):
    rng = np.random.default_rng(seed)
    m = np.ones(N, dtype=np.float64) / N
    r = (rng.random((N, 3)) - 0.5).astype(np.float64)
    v = (0.1 * (rng.random((N, 3)) - 0.5)).astype(np.float64)
    return m, r, v


def compute_energy(N, m, r, v, rsoft2):
    KE = 0.0
    for i in range(N):
        KE += 0.5 * m[i] * (v[i, 0]**2 + v[i, 1]**2 + v[i, 2]**2)
    PE = 0.0
    for i in range(N - 1):
        for j in range(i + 1, N):
            dx = r[i, 0] - r[j, 0]
            dy = r[i, 1] - r[j, 1]
            dz = r[i, 2] - r[j, 2]
            PE -= G * m[i] * m[j] / np.sqrt(dx*dx + dy*dy + dz*dz + rsoft2)
    return KE, PE, KE + PE


# ═══════════════════════════════════════════════════════════════════════════
#  Conversion helpers
# ═══════════════════════════════════════════════════════════════════════════

def aos_to_soa(r):
    return r[:, 0].copy(), r[:, 1].copy(), r[:, 2].copy()

def soa_to_aos(rx, ry, rz):
    return np.column_stack([rx, ry, rz])

def aos_to_aosoa(arr, VL):
    N = arr.shape[0]
    NB = (N + VL - 1) // VL
    if arr.ndim == 2:
        out = np.zeros((NB, 3, VL), dtype=arr.dtype)
        for i in range(N):
            out[i // VL, :, i % VL] = arr[i, :]
    else:
        out = np.zeros((NB, VL), dtype=arr.dtype)
        for i in range(N):
            out[i // VL, i % VL] = arr[i]
    return out

def aosoa_to_aos(arr, N, VL):
    out = np.zeros((N, 3), dtype=arr.dtype)
    for i in range(N):
        out[i, :] = arr[i // VL, :, i % VL]
    return out


# ═══════════════════════════════════════════════════════════════════════════
#  1.  AoS — r[N, 3], v[N, 3], a[N, 3], m[N]
# ═══════════════════════════════════════════════════════════════════════════

def aos_kick(N, v, a, dt):
    for i in range(N):
        v[i, 0] += a[i, 0] * dt
        v[i, 1] += a[i, 1] * dt
        v[i, 2] += a[i, 2] * dt


def aos_drift(N, r, v, dt):
    for i in range(N):
        r[i, 0] += v[i, 0] * dt
        r[i, 1] += v[i, 1] * dt
        r[i, 2] += v[i, 2] * dt


def aos_accel_local(N, m, r, a, rsoft2):
    """Full N×N with local scalar accumulators.
    i: parallel map.  j: sequential local reduction.  O(N) memory."""
    for i in range(N):
        axi = 0.0
        ayi = 0.0
        azi = 0.0
        rxi = r[i, 0]
        ryi = r[i, 1]
        rzi = r[i, 2]
        for j in range(N):
            dx = rxi - r[j, 0]
            dy = ryi - r[j, 1]
            dz = rzi - r[j, 2]
            ir  = 1.0 / np.sqrt(dx*dx + dy*dy + dz*dz + rsoft2)
            ir3 = ir * ir * ir
            axi -= m[j] * dx * ir3
            ayi -= m[j] * dy * ir3
            azi -= m[j] * dz * ir3
        a[i, 0] = axi
        a[i, 1] = ayi
        a[i, 2] = azi


def aos_accel_buf(N, m, r, a, rsoft2, fx, fy, fz):
    """Full N×N with external buffer fx[N,N], fy[N,N], fz[N,N].
    Phase 1 — MAP(i, j):   fill fx[i,j] independently.  Both i,j parallel.
    Phase 2 — REDUCE(j):   a[i] = sum_j fx[i,j].  Parallel over i."""
    # Phase 1: map — each (i,j) writes to its own slot, no conflicts
    for i in range(N):
        rxi = r[i, 0]
        ryi = r[i, 1]
        rzi = r[i, 2]
        for j in range(N):
            dx = rxi - r[j, 0]
            dy = ryi - r[j, 1]
            dz = rzi - r[j, 2]
            ir  = 1.0 / np.sqrt(dx*dx + dy*dy + dz*dz + rsoft2)
            ir3 = ir * ir * ir
            fx[i, j] = -m[j] * dx * ir3
            fy[i, j] = -m[j] * dy * ir3
            fz[i, j] = -m[j] * dz * ir3

    # Phase 2: reduce — sum over j
    for i in range(N):
        axi = 0.0
        ayi = 0.0
        azi = 0.0
        for j in range(N):
            axi += fx[i, j]
            ayi += fy[i, j]
            azi += fz[i, j]
        a[i, 0] = axi
        a[i, 1] = ayi
        a[i, 2] = azi


def aos_run(N, m, r, v, a, dt, nsteps, rsoft2, accel_fn, **kw):
    half_dt = 0.5 * dt
    accel_fn(N, m, r, a, rsoft2, **kw)
    for _ in range(nsteps):
        aos_kick(N, v, a, half_dt)
        aos_drift(N, r, v, dt)
        accel_fn(N, m, r, a, rsoft2, **kw)
        aos_kick(N, v, a, half_dt)


# ═══════════════════════════════════════════════════════════════════════════
#  2.  SoA — rx[N], ry[N], rz[N], ...
# ═══════════════════════════════════════════════════════════════════════════

def soa_kick(N, vx, vy, vz, ax, ay, az, dt):
    for i in range(N):
        vx[i] += ax[i] * dt
        vy[i] += ay[i] * dt
        vz[i] += az[i] * dt


def soa_drift(N, rx, ry, rz, vx, vy, vz, dt):
    for i in range(N):
        rx[i] += vx[i] * dt
        ry[i] += vy[i] * dt
        rz[i] += vz[i] * dt


def soa_accel_local(N, m, rx, ry, rz, ax, ay, az, rsoft2):
    """Full N×N with local scalar accumulators.  O(N) memory."""
    for i in range(N):
        axi = 0.0
        ayi = 0.0
        azi = 0.0
        rxi = rx[i]
        ryi = ry[i]
        rzi = rz[i]
        for j in range(N):
            dx = rxi - rx[j]
            dy = ryi - ry[j]
            dz = rzi - rz[j]
            ir  = 1.0 / np.sqrt(dx*dx + dy*dy + dz*dz + rsoft2)
            ir3 = ir * ir * ir
            axi -= m[j] * dx * ir3
            ayi -= m[j] * dy * ir3
            azi -= m[j] * dz * ir3
        ax[i] = axi
        ay[i] = ayi
        az[i] = azi


def soa_accel_buf(N, m, rx, ry, rz, ax, ay, az, rsoft2, fx, fy, fz):
    """Full N×N with external buffer fx[N,N], fy[N,N], fz[N,N].
    Phase 1 — MAP(i, j):   both parallel.
    Phase 2 — REDUCE(j):   parallel over i."""
    # Phase 1: map
    for i in range(N):
        rxi = rx[i]
        ryi = ry[i]
        rzi = rz[i]
        for j in range(N):
            dx = rxi - rx[j]
            dy = ryi - ry[j]
            dz = rzi - rz[j]
            ir  = 1.0 / np.sqrt(dx*dx + dy*dy + dz*dz + rsoft2)
            ir3 = ir * ir * ir
            fx[i, j] = -m[j] * dx * ir3
            fy[i, j] = -m[j] * dy * ir3
            fz[i, j] = -m[j] * dz * ir3

    # Phase 2: reduce
    for i in range(N):
        axi = 0.0
        ayi = 0.0
        azi = 0.0
        for j in range(N):
            axi += fx[i, j]
            ayi += fy[i, j]
            azi += fz[i, j]
        ax[i] = axi
        ay[i] = ayi
        az[i] = azi


def soa_run(N, m, rx, ry, rz, vx, vy, vz, ax, ay, az,
            dt, nsteps, rsoft2, accel_fn, **kw):
    half_dt = 0.5 * dt
    accel_fn(N, m, rx, ry, rz, ax, ay, az, rsoft2, **kw)
    for _ in range(nsteps):
        soa_kick(N, vx, vy, vz, ax, ay, az, half_dt)
        soa_drift(N, rx, ry, rz, vx, vy, vz, dt)
        accel_fn(N, m, rx, ry, rz, ax, ay, az, rsoft2, **kw)
        soa_kick(N, vx, vy, vz, ax, ay, az, half_dt)


# ═══════════════════════════════════════════════════════════════════════════
#  3.  AoSoA — r[NB, 3, VL], m[NB, VL]
# ═══════════════════════════════════════════════════════════════════════════

def aosoa_kick(N, VL, v, a, dt):
    for i in range(N):
        bi = i // VL
        li = i % VL
        v[bi, 0, li] += a[bi, 0, li] * dt
        v[bi, 1, li] += a[bi, 1, li] * dt
        v[bi, 2, li] += a[bi, 2, li] * dt


def aosoa_drift(N, VL, r, v, dt):
    for i in range(N):
        bi = i // VL
        li = i % VL
        r[bi, 0, li] += v[bi, 0, li] * dt
        r[bi, 1, li] += v[bi, 1, li] * dt
        r[bi, 2, li] += v[bi, 2, li] * dt


def aosoa_accel_local(N, VL, m, r, a, rsoft2):
    """Full N×N with local scalar accumulators.  O(N) memory."""
    for i in range(N):
        bi_i = i // VL
        li_i = i % VL
        axi = 0.0
        ayi = 0.0
        azi = 0.0
        rxi = r[bi_i, 0, li_i]
        ryi = r[bi_i, 1, li_i]
        rzi = r[bi_i, 2, li_i]
        for j in range(N):
            bi_j = j // VL
            li_j = j % VL
            dx = rxi - r[bi_j, 0, li_j]
            dy = ryi - r[bi_j, 1, li_j]
            dz = rzi - r[bi_j, 2, li_j]
            ir  = 1.0 / np.sqrt(dx*dx + dy*dy + dz*dz + rsoft2)
            ir3 = ir * ir * ir
            axi -= m[bi_j, li_j] * dx * ir3
            ayi -= m[bi_j, li_j] * dy * ir3
            azi -= m[bi_j, li_j] * dz * ir3
        a[bi_i, 0, li_i] = axi
        a[bi_i, 1, li_i] = ayi
        a[bi_i, 2, li_i] = azi


def aosoa_accel_buf(N, VL, m, r, a, rsoft2, fx, fy, fz):
    """Full N×N with external buffer fx[N,N], fy[N,N], fz[N,N].
    Phase 1 — MAP(i, j):   both parallel.  Flat indexing into buffer.
    Phase 2 — REDUCE(j):   parallel over i."""
    # Phase 1: map
    for i in range(N):
        bi_i = i // VL
        li_i = i % VL
        rxi = r[bi_i, 0, li_i]
        ryi = r[bi_i, 1, li_i]
        rzi = r[bi_i, 2, li_i]
        for j in range(N):
            bi_j = j // VL
            li_j = j % VL
            dx = rxi - r[bi_j, 0, li_j]
            dy = ryi - r[bi_j, 1, li_j]
            dz = rzi - r[bi_j, 2, li_j]
            ir  = 1.0 / np.sqrt(dx*dx + dy*dy + dz*dz + rsoft2)
            ir3 = ir * ir * ir
            fx[i, j] = -m[bi_j, li_j] * dx * ir3
            fy[i, j] = -m[bi_j, li_j] * dy * ir3
            fz[i, j] = -m[bi_j, li_j] * dz * ir3

    # Phase 2: reduce
    for i in range(N):
        bi_i = i // VL
        li_i = i % VL
        axi = 0.0
        ayi = 0.0
        azi = 0.0
        for j in range(N):
            axi += fx[i, j]
            ayi += fy[i, j]
            azi += fz[i, j]
        a[bi_i, 0, li_i] = axi
        a[bi_i, 1, li_i] = ayi
        a[bi_i, 2, li_i] = azi


def aosoa_run(N, VL, m, r, v, a, dt, nsteps, rsoft2, accel_fn, **kw):
    half_dt = 0.5 * dt
    accel_fn(N, VL, m, r, a, rsoft2, **kw)
    for _ in range(nsteps):
        aosoa_kick(N, VL, v, a, half_dt)
        aosoa_drift(N, VL, r, v, dt)
        accel_fn(N, VL, m, r, a, rsoft2, **kw)
        aosoa_kick(N, VL, v, a, half_dt)


# ═══════════════════════════════════════════════════════════════════════════
#  4.  DaCe — SoA, symbolic n, both local and buffer variants
# ═══════════════════════════════════════════════════════════════════════════

def build_dace_kernels():
    import dace

    n = dace.symbol('n', dtype=dace.int64)

    @dace.program
    def dace_kick(
        vx: dace.float64[n],
        vy: dace.float64[n],
        vz: dace.float64[n],
        ax: dace.float64[n],
        ay: dace.float64[n],
        az: dace.float64[n],
        dt: dace.float64,
    ):
        for i in dace.map[0:n]:
            vx[i] = vx[i] + ax[i] * dt
            vy[i] = vy[i] + ay[i] * dt
            vz[i] = vz[i] + az[i] * dt

    @dace.program
    def dace_drift(
        rx: dace.float64[n],
        ry: dace.float64[n],
        rz: dace.float64[n],
        vx: dace.float64[n],
        vy: dace.float64[n],
        vz: dace.float64[n],
        dt: dace.float64,
    ):
        for i in dace.map[0:n]:
            rx[i] = rx[i] + vx[i] * dt
            ry[i] = ry[i] + vy[i] * dt
            rz[i] = rz[i] + vz[i] * dt

    @dace.program
    def dace_accel_local(
        m:      dace.float64[n],
        rx:     dace.float64[n],
        ry:     dace.float64[n],
        rz:     dace.float64[n],
        ax:     dace.float64[n],
        ay:     dace.float64[n],
        az:     dace.float64[n],
        rsoft2: dace.float64,
    ):
        for i in dace.map[0:n]:
            axi = 0.0
            ayi = 0.0
            azi = 0.0
            rxi = rx[i]
            ryi = ry[i]
            rzi = rz[i]
            for j in range(n):
                dx  = rxi - rx[j]
                dy  = ryi - ry[j]
                dz  = rzi - rz[j]
                dsq = dx*dx + dy*dy + dz*dz + rsoft2
                inv_dist  = 1.0 / np.sqrt(dsq)
                inv_dist3 = inv_dist * inv_dist * inv_dist
                axi = axi - m[j] * dx * inv_dist3
                ayi = ayi - m[j] * dy * inv_dist3
                azi = azi - m[j] * dz * inv_dist3
            ax[i] = axi
            ay[i] = ayi
            az[i] = azi

    @dace.program
    def dace_accel_buf(
        m:      dace.float64[n],
        rx:     dace.float64[n],
        ry:     dace.float64[n],
        rz:     dace.float64[n],
        ax:     dace.float64[n],
        ay:     dace.float64[n],
        az:     dace.float64[n],
        rsoft2: dace.float64,
        fx:     dace.float64[n, n],
        fy:     dace.float64[n, n],
        fz:     dace.float64[n, n],
    ):
        # Phase 1: 2D map — both i and j parallel, no conflicts
        for i, j in dace.map[0:n, 0:n]:
            dx  = rx[i] - rx[j]
            dy  = ry[i] - ry[j]
            dz  = rz[i] - rz[j]
            dsq = dx*dx + dy*dy + dz*dz + rsoft2
            inv_dist  = 1.0 / np.sqrt(dsq)
            inv_dist3 = inv_dist * inv_dist * inv_dist
            fx[i, j] = -m[j] * dx * inv_dist3
            fy[i, j] = -m[j] * dy * inv_dist3
            fz[i, j] = -m[j] * dz * inv_dist3

        # Phase 2: reduce over j
        for i in dace.map[0:n]:
            axi = 0.0
            ayi = 0.0
            azi = 0.0
            for j in range(n):
                axi = axi + fx[i, j]
                ayi = ayi + fy[i, j]
                azi = azi + fz[i, j]
            ax[i] = axi
            ay[i] = ayi
            az[i] = azi

    return dace_kick, dace_drift, dace_accel_local, dace_accel_buf


# ═══════════════════════════════════════════════════════════════════════════
#  Driver
# ═══════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N",     type=int,   default=128)
    ap.add_argument("--steps", type=int,   default=5)
    ap.add_argument("--dt",    type=float, default=0.001)
    ap.add_argument("--vl",    type=int,   default=8)
    ap.add_argument("--seed",  type=int,   default=42)
    ap.add_argument("--ic",    choices=["circular", "random"], default="circular")
    ap.add_argument("--dace",  action="store_true")
    args = ap.parse_args()

    N = args.N
    nsteps = args.steps
    dt = args.dt
    VL = args.vl
    rsoft2 = RSOFT2

    if N % VL != 0:
        N = ((N + VL - 1) // VL) * VL
        print(f"  Padded N to {N} (multiple of VL={VL})")

    print("=" * 70)
    print(f"  N-body  N={N}  steps={nsteps}  dt={dt}  VL={VL}")
    print("=" * 70)

    m0, r0, v0 = (init_circular if args.ic == "circular" else init_random)(N, args.seed)
    _, _, E0 = compute_energy(N, m0, r0, v0, rsoft2)
    print(f"  E0 = {E0:.10f}")
    print()

    # ---- allocate N×N buffers once, reused by all _buf variants -----------
    fx = np.zeros((N, N), dtype=np.float64)
    fy = np.zeros((N, N), dtype=np.float64)
    fz = np.zeros((N, N), dtype=np.float64)

    timings = {}

    # ── AoS local ──────────────────────────────────────────────────────────
    m = m0.copy()
    r = r0.copy()
    v = v0.copy()
    a = np.zeros_like(r)
    t0 = time.perf_counter()
    aos_run(N, m, r, v, a, dt, nsteps, rsoft2, aos_accel_local)
    timings["AoS local"] = time.perf_counter() - t0
    r_ref = r.copy()
    _, _, E = compute_energy(N, m, r, v, rsoft2)
    print(f"  AoS  local : {timings['AoS local']:7.4f}s  dE/E0={abs((E-E0)/E0):.2e}")

    # ── AoS buf ────────────────────────────────────────────────────────────
    m = m0.copy()
    r = r0.copy()
    v = v0.copy()
    a = np.zeros_like(r)
    t0 = time.perf_counter()
    aos_run(N, m, r, v, a, dt, nsteps, rsoft2, aos_accel_buf, fx=fx, fy=fy, fz=fz)
    timings["AoS buf"] = time.perf_counter() - t0
    err = np.max(np.abs(r_ref - r))
    print(f"  AoS  buf   : {timings['AoS buf']:7.4f}s  max|Δr|={err:.2e}")

    # ── SoA local ──────────────────────────────────────────────────────────
    m = m0.copy()
    rx, ry, rz = aos_to_soa(r0)
    vx, vy, vz = aos_to_soa(v0)
    ax = np.zeros(N)
    ay = np.zeros(N)
    az = np.zeros(N)
    t0 = time.perf_counter()
    soa_run(N, m, rx, ry, rz, vx, vy, vz, ax, ay, az,
            dt, nsteps, rsoft2, soa_accel_local)
    timings["SoA local"] = time.perf_counter() - t0
    err = np.max(np.abs(r_ref - soa_to_aos(rx, ry, rz)))
    print(f"  SoA  local : {timings['SoA local']:7.4f}s  max|Δr|={err:.2e}")

    # ── SoA buf ────────────────────────────────────────────────────────────
    m = m0.copy()
    rx, ry, rz = aos_to_soa(r0)
    vx, vy, vz = aos_to_soa(v0)
    ax = np.zeros(N)
    ay = np.zeros(N)
    az = np.zeros(N)
    t0 = time.perf_counter()
    soa_run(N, m, rx, ry, rz, vx, vy, vz, ax, ay, az,
            dt, nsteps, rsoft2, soa_accel_buf, fx=fx, fy=fy, fz=fz)
    timings["SoA buf"] = time.perf_counter() - t0
    err = np.max(np.abs(r_ref - soa_to_aos(rx, ry, rz)))
    print(f"  SoA  buf   : {timings['SoA buf']:7.4f}s  max|Δr|={err:.2e}")

    # ── AoSoA local ────────────────────────────────────────────────────────
    m_a = aos_to_aosoa(m0, VL)
    r_a = aos_to_aosoa(r0, VL)
    v_a = aos_to_aosoa(v0, VL)
    a_a = np.zeros_like(r_a)
    t0 = time.perf_counter()
    aosoa_run(N, VL, m_a, r_a, v_a, a_a, dt, nsteps, rsoft2, aosoa_accel_local)
    timings["AoSoA local"] = time.perf_counter() - t0
    err = np.max(np.abs(r_ref - aosoa_to_aos(r_a, N, VL)))
    print(f"  AoSoA local: {timings['AoSoA local']:7.4f}s  max|Δr|={err:.2e}")

    # ── AoSoA buf ──────────────────────────────────────────────────────────
    m_a = aos_to_aosoa(m0, VL)
    r_a = aos_to_aosoa(r0, VL)
    v_a = aos_to_aosoa(v0, VL)
    a_a = np.zeros_like(r_a)
    t0 = time.perf_counter()
    aosoa_run(N, VL, m_a, r_a, v_a, a_a, dt, nsteps, rsoft2,
              aosoa_accel_buf, fx=fx, fy=fy, fz=fz)
    timings["AoSoA buf"] = time.perf_counter() - t0
    err = np.max(np.abs(r_ref - aosoa_to_aos(r_a, N, VL)))
    print(f"  AoSoA buf  : {timings['AoSoA buf']:7.4f}s  max|Δr|={err:.2e}")

    # ── DaCe ───────────────────────────────────────────────────────────────
    if args.dace:
        try:
            print("  DaCe: compiling...", end=" ", flush=True)
            dk, dd, da_local, da_buf = build_dace_kernels()
            dk_c = dk.compile()
            dd_c = dd.compile()
            dal_c = da_local.compile()
            dab_c = da_buf.compile()
            print("done")

            rs = np.float64(rsoft2)
            half_dt = np.float64(0.5 * dt)
            dt_ = np.float64(dt)

            for label, accel_c, needs_buf in [("DaCe local", dal_c, False),
                                               ("DaCe buf",   dab_c, True)]:
                m = m0.copy()
                rx, ry, rz = aos_to_soa(r0)
                vx, vy, vz = aos_to_soa(v0)
                ax = np.zeros(N)
                ay = np.zeros(N)
                az = np.zeros(N)

                akw = dict(m=m, rx=rx, ry=ry, rz=rz,
                           ax=ax, ay=ay, az=az, rsoft2=rs, n=N)
                if needs_buf:
                    akw.update(fx=fx, fy=fy, fz=fz)

                t0 = time.perf_counter()
                accel_c(**akw)
                for _ in range(nsteps):
                    dk_c(vx=vx, vy=vy, vz=vz, ax=ax, ay=ay, az=az,
                         dt=half_dt, n=N)
                    dd_c(rx=rx, ry=ry, rz=rz, vx=vx, vy=vy, vz=vz,
                         dt=dt_, n=N)
                    accel_c(**akw)
                    dk_c(vx=vx, vy=vy, vz=vz, ax=ax, ay=ay, az=az,
                         dt=half_dt, n=N)
                elapsed = time.perf_counter() - t0

                err = np.max(np.abs(r_ref - soa_to_aos(rx, ry, rz)))
                print(f"  {label:13s}: {elapsed:7.4f}s  max|Δr|={err:.2e}")
                timings[label] = elapsed

        except ImportError:
            print("  DaCe: skipped (pip install dace)")
        except Exception as e:
            print(f"  DaCe: error: {e}")

    # ── summary ────────────────────────────────────────────────────────────
    print()
    print("-" * 70)
    t_base = timings["AoS local"]
    for k, t in timings.items():
        print(f"  {k:14s} {t:7.4f}s  {t_base/t:5.2f}× vs AoS local")

    print()
    print(f"  Buffer memory: 3 × {N}×{N} × 8 = "
          f"{3 * N * N * 8 / 1024:.1f} KiB")
    print(f"  Layout shapes:")
    print(f"    AoS    r[{N}, 3]")
    print(f"    SoA    rx[{N}]")
    print(f"    AoSoA  r[{N//VL}, 3, {VL}]")
    print(f"    buf    fx[{N}, {N}]  fy[{N}, {N}]  fz[{N}, {N}]")


if __name__ == "__main__":
    main()