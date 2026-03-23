#!/usr/bin/env python3
"""CPU matrix transpose sweep: custom kernels + OpenBLAS somatcopy, with STREAM-like roofline."""
import subprocess, sys, os, csv, platform, re
from collections import defaultdict
import numpy as np
from pathlib import Path

# ── Paths ──
BINARY_KERN = "./transpose_cpu"
BINARY_BLAS = "./transpose_openblas"
CSV_RAW     = "transpose_cpu_raw.csv"
CSV_AGG     = "transpose_cpu_results.csv"

# ── Defaults (overridable via env) ──
N       = int(os.environ.get("CPU_TR_N", 16384))
REPS    = int(os.environ.get("CPU_TR_REPS", 100))
WARMUP  = int(os.environ.get("CPU_TR_WARMUP", 5))
THREADS = int(os.environ.get("OMP_NUM_THREADS", 64))

# ── Variant metadata ──
KERN_NAMES = {
    0:  "naive",             1:  "naive_c2",
    2:  "tiled",             3:  "tiled_c2",
    4:  "blk_naive",         5:  "blk_naive_c2",
    6:  "blk_tiled",         7:  "blk_tiled_c2",
    8:  "blk_aligned",       9:  "blk_aligned_c2",
    10: "blk_aligned_mt",    11: "blk_aligned_mt_c2",
    12: "locbuf",            13: "locbuf_c2",
    14: "locbuf_blk",        15: "locbuf_blk_c2",
    16: "locbuf_blk_mt",     17: "locbuf_blk_mt_c2",
    18: "locbuf_2buf",       19: "locbuf_2buf_c2",
}
BLAS_NAMES = {0: "openblas", 1: "openblas_blk", 2: "openblas_blk_omp"}
LIB_NAMES  = set(BLAS_NAMES.values())

# ── Sweep configs ──
# (variant, TB, SB, MT)
# Row-major variants: sweep TB
RM_VARIANTS = [0, 1, 2, 3]
TB_VALS     = [16, 32, 64, 128]

# Blocked variants (no locbuf): sweep SB
BLK_VARIANTS = [4, 5, 6, 7, 8, 9, 10, 11]
SB_VALS      = [8, 16, 32, 64]
MT_VALS      = [4, 8, 16]

# Locbuf row-major: sweep TB
LOCBUF_RM_VARIANTS = [12, 13, 18, 19]

# Locbuf blocked: sweep SB, MT
LOCBUF_BLK_VARIANTS = [14, 15, 16, 17]

# OpenBLAS blocked: sweep SB
BLAS_BLK_SB_VALS = [8, 16, 32, 64]


# ══════════════════════════════════════════════════════════════════════
#  System info + STREAM-like bandwidth measurement
# ══════════════════════════════════════════════════════════════════════

STREAM_SRC = r'''
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <cstring>

int main(int argc, char** argv) {
    int threads = (argc > 1) ? atoi(argv[1]) : 0;
    if (threads > 0) omp_set_num_threads(threads);

    const long long SZ = 64LL * 1024 * 1024;  // 64M doubles = 512 MB
    const size_t bytes = SZ * sizeof(double);
    double* a = (double*)aligned_alloc(64, bytes);
    double* b = (double*)aligned_alloc(64, bytes);

    // Init (parallel, first-touch)
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < SZ; i++) { a[i] = 1.0; b[i] = 0.0; }

    // Warmup
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < SZ; i++) b[i] = a[i];

    const int REPS = 100;
    double best = 1e30;
    for (int r = 0; r < REPS; r++) {
        double t0 = omp_get_wtime();
        #pragma omp parallel for schedule(static)
        for (long long i = 0; i < SZ; i++) b[i] = a[i];
        double t1 = omp_get_wtime();
        best = std::min(best, t1 - t0);
    }
    // bytes moved: read a + write b = 2 * bytes
    printf("%.1f\n", 2.0 * bytes / best / 1e9);
    free(a); free(b);
    return 0;
}
'''

def get_cpu_info():
    """Return a dict with CPU name and core count."""
    info = {"name": "unknown", "cores": 64}
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    info["name"] = line.split(":", 1)[1].strip()
                    break
    except FileNotFoundError:
        pass
    return info


def measure_bandwidth(threads):
    """Compile and run a STREAM-copy benchmark, return empirical GB/s."""
    cache = Path(".cache")
    cache.mkdir(exist_ok=True)
    src  = cache / "stream_copy.cpp"
    bbin = cache / "stream_copy"

    src.write_text(STREAM_SRC)
    cmd = f"g++ -O3 -march=native -fopenmp -o {bbin} {src}"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [ERROR] STREAM compile failed:\n{r.stderr.strip()}")
        sys.exit(1)

    args = [str(bbin)]
    if threads > 0:
        args.append(str(threads))
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [ERROR] STREAM binary failed:\n{r.stderr.strip()}")
        sys.exit(1)

    return float(r.stdout.strip())


# ══════════════════════════════════════════════════════════════════════
#  Compile
# ══════════════════════════════════════════════════════════════════════

def compile_kernels(force=False):
    if Path(BINARY_KERN).exists() and not force:
        print(f"  {BINARY_KERN} exists, skipping (use --compile to force)")
        return True
    cmd = f"g++ -O3 -march=native -fopenmp -o {BINARY_KERN} transpose_cpu.cpp"
    print(f"  Compiling kernels: {cmd}")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [ERROR] kernel compile failed:\n{r.stderr.strip()}")
        return False
    return True


def compile_openblas(force=False):
    if Path(BINARY_BLAS).exists() and not force:
        print(f"  {BINARY_BLAS} exists, skipping")
        return True

    # Try common OpenBLAS locations
    blas_flags = ""
    for prefix in ["/usr", "/usr/local", "/opt/OpenBLAS"]:
        inc = Path(prefix) / "include" / "cblas.h"
        if not inc.exists():
            inc = Path(prefix) / "include" / "openblas" / "cblas.h"
        if inc.exists():
            inc_dir = str(inc.parent)
            lib_dir = str(Path(prefix) / "lib")
            blas_flags = f"-I{inc_dir} -L{lib_dir}"
            break

    cmd = (f"g++ -O3 -march=native -fopenmp {blas_flags} "
           f"-o {BINARY_BLAS} transpose_openblas.cpp -lopenblas")
    print(f"  Compiling OpenBLAS: {cmd}")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [WARN] OpenBLAS compile failed: {r.stderr.strip()[:200]}")
        print(f"  OpenBLAS benchmark will be skipped.")
        return False
    return True


# ══════════════════════════════════════════════════════════════════════
#  Sweep
# ══════════════════════════════════════════════════════════════════════

def run_one(binary, args, label=""):
    """Run a single benchmark, print output, return True on success."""
    try:
        r = subprocess.run([binary] + [str(a) for a in args],
                           capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            print(f"  FAIL {label}: {r.stderr.strip()[:120]}")
            return False
        print(f"  {r.stdout.strip()}")
        return True
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT {label}")
        return False


def sweep_kernels():
    """Sweep all transpose_cpu variants."""
    # transpose_cpu args: N variant csv TB SB MT WARMUP REPS THREADS

    # Row-major variants: only TB matters
    for var in RM_VARIANTS:
        tbs = TB_VALS if var in (2, 3) else [64]  # naive doesn't use TB
        for tb in tbs:
            label = f"{KERN_NAMES[var]} TB={tb}"
            run_one(BINARY_KERN, [N, var, CSV_RAW, tb, 32, 8, WARMUP, REPS, THREADS], label)

    # Blocked variants (non-locbuf)
    for var in BLK_VARIANTS:
        for sb in SB_VALS:
            if N % sb != 0:
                continue
            mts = MT_VALS if var in (10, 11) else [8]  # only _mt variants use MT
            for mt in mts:
                label = f"{KERN_NAMES[var]} SB={sb} MT={mt}"
                run_one(BINARY_KERN, [N, var, CSV_RAW, 64, sb, mt, WARMUP, REPS, THREADS], label)

    # Locbuf row-major
    for var in LOCBUF_RM_VARIANTS:
        for tb in TB_VALS:
            label = f"{KERN_NAMES[var]} TB={tb}"
            run_one(BINARY_KERN, [N, var, CSV_RAW, tb, 32, 8, WARMUP, REPS, THREADS], label)

    # Locbuf blocked
    for var in LOCBUF_BLK_VARIANTS:
        for sb in SB_VALS:
            if N % sb != 0:
                continue
            mts = MT_VALS if var in (16, 17) else [8]
            for mt in mts:
                label = f"{KERN_NAMES[var]} SB={sb} MT={mt}"
                run_one(BINARY_KERN, [N, var, CSV_RAW, 64, sb, mt, WARMUP, REPS, THREADS], label)


def sweep_openblas():
    """Sweep OpenBLAS somatcopy variants."""
    if not Path(BINARY_BLAS).exists():
        return
    print(f"\n  -- Library: OpenBLAS --")
    # transpose_openblas args: N variant csv SB WARMUP REPS THREADS

    # V0: row-major
    run_one(BINARY_BLAS, [N, 0, CSV_RAW, 32, WARMUP, REPS, THREADS], "openblas row-major")

    # V1, V2: blocked per SB
    for var in (1, 2):
        for sb in BLAS_BLK_SB_VALS:
            if N % sb != 0:
                continue
            label = f"{BLAS_NAMES[var]} SB={sb}"
            run_one(BINARY_BLAS, [N, var, CSV_RAW, sb, WARMUP, REPS, THREADS], label)


# ══════════════════════════════════════════════════════════════════════
#  Aggregate
# ══════════════════════════════════════════════════════════════════════

def aggregate(emp_bw):
    """Read CSV_RAW, compute per-config stats, write CSV_AGG, return rows."""
    if not os.path.exists(CSV_RAW):
        print(f"  [ERROR] {CSV_RAW} not found")
        sys.exit(1)

    # CSV columns: variant, N, TB, SB, MT, nthreads, rep, time_s, gbs, cksum, status
    groups = defaultdict(list)
    with open(CSV_RAW) as f:
        for line in f:
            p = line.strip().split(",")
            if len(p) < 11:
                continue
            key = (p[0], p[1], p[2], p[3], p[4], p[5])  # variant, N, TB, SB, MT, threads
            groups[key].append((float(p[7]), float(p[8]), p[10]))  # time_s, gbs, status

    if not groups:
        print(f"  [ERROR] {CSV_RAW} empty")
        sys.exit(1)

    rows = []
    for key, iters in groups.items():
        gbs = np.array([g for _, g, _ in iters])
        ts  = np.array([t for t, _, _ in iters])
        status = iters[0][2]
        med_gbs = float(np.median(gbs))
        rows.append(dict(
            variant=key[0], N=key[1], TB=key[2], SB=key[3], MT=key[4], threads=key[5],
            med_time_ms=f"{np.median(ts)*1000:.4f}",
            std_time_ms=f"{np.std(ts)*1000:.4f}",
            p5_time_ms=f"{np.percentile(ts, 5)*1000:.4f}",
            p95_time_ms=f"{np.percentile(ts, 95)*1000:.4f}",
            min_gbs=f"{np.min(gbs):.1f}",
            med_gbs=f"{med_gbs:.1f}",
            max_gbs=f"{np.max(gbs):.1f}",
            pct_peak_bw=f"{100*med_gbs/emp_bw:.1f}" if emp_bw > 0 else "—",
            status=status,
        ))

    fields = ["variant", "N", "TB", "SB", "MT", "threads",
              "med_time_ms", "std_time_ms", "p5_time_ms", "p95_time_ms",
              "min_gbs", "med_gbs", "max_gbs", "pct_peak_bw", "status"]
    with open(CSV_AGG, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: -float(r["med_gbs"])))
    return rows


# ══════════════════════════════════════════════════════════════════════
#  Report
# ══════════════════════════════════════════════════════════════════════

def report(cpu_info, emp_bw, rows):
    bpi = 2.0 * N * N * 4  # bytes per invocation (read + write, float32)

    print(f"\n{'='*110}")
    print(f" CPU: {cpu_info['name']}   Cores: {cpu_info['cores']}   Emp BW: {emp_bw:.1f} GB/s")
    print(f" Transpose: N={N}  bytes/call = {bpi/1e9:.3f} GB (read+write, float32)")
    print(f"{'='*110}")

    kern_rows = [r for r in rows if r["variant"] not in LIB_NAMES]
    lib_rows  = [r for r in rows if r["variant"] in LIB_NAMES]

    # Best library row-major as speedup reference
    lib_rm  = [r for r in lib_rows if "blk" not in r["variant"]]
    ref_gbs = float(max(lib_rm, key=lambda x: float(x["med_gbs"]))["med_gbs"]) if lib_rm else None

    # Header
    hdr = (f" {'variant':<22} {'TB':>3} {'SB':>3} {'MT':>2} {'thr':>3}"
           f"  {'medGB/s':>8} {'minGB/s':>8} {'maxGB/s':>8} {'%pkBW':>6}"
           f"  {'medMs':>8} {'p5ms':>8} {'p95ms':>8}"
           f"  {'vs lib':>7}  {'stat':>4}")
    print(f"\n{hdr}")
    print(f" {'-'*len(hdr)}")

    for r in sorted(rows, key=lambda x: -float(x["med_gbs"])):
        med = float(r["med_gbs"])
        vs  = f"{100*med/ref_gbs:.0f}%" if ref_gbs else "--"
        is_lib = r["variant"] in LIB_NAMES
        tag = " *" if is_lib else "  "
        print(f" {r['variant']:<22} {r['TB']:>3} {r['SB']:>3} {r['MT']:>2} {r['threads']:>3}"
              f"  {med:8.1f} {float(r['min_gbs']):8.1f} {float(r['max_gbs']):8.1f}"
              f"  {r['pct_peak_bw']:>5}%"
              f"  {r['med_time_ms']:>8} {r['p5_time_ms']:>8} {r['p95_time_ms']:>8}"
              f"  {vs:>6}  {r['status']}{tag}")

    # Summary
    if lib_rows:
        best_lib = max(lib_rows, key=lambda x: float(x["med_gbs"]))
        print(f"\n Best library : {best_lib['variant']} SB={best_lib['SB']}"
              f" -> {best_lib['med_gbs']} GB/s ({best_lib['pct_peak_bw']}% peak)")
    if kern_rows:
        best_kern = max(kern_rows, key=lambda x: float(x["med_gbs"]))
        bk = float(best_kern["med_gbs"])
        extra = f"  ({100*bk/ref_gbs:.0f}% of best library)" if ref_gbs else ""
        print(f" Best kernel  : {best_kern['variant']} TB={best_kern['TB']} SB={best_kern['SB']}"
              f" MT={best_kern['MT']}"
              f" -> {best_kern['med_gbs']} GB/s ({best_kern['pct_peak_bw']}% peak){extra}")

    # Per-category bests
    categories = {
        "Row-major naive":      lambda r: r["variant"] in ("naive", "naive_c2"),
        "Row-major tiled":      lambda r: r["variant"] in ("tiled", "tiled_c2"),
        "Blocked aligned":      lambda r: r["variant"].startswith("blk_aligned"),
        "Locbuf row-major":     lambda r: r["variant"] in ("locbuf", "locbuf_c2", "locbuf_2buf", "locbuf_2buf_c2"),
        "Locbuf blocked":       lambda r: r["variant"].startswith("locbuf_blk"),
        "OpenBLAS":             lambda r: r["variant"].startswith("openblas"),
    }
    print(f"\n Per-category bests:")
    for cat, pred in categories.items():
        subset = [r for r in rows if pred(r)]
        if subset:
            best = max(subset, key=lambda x: float(x["med_gbs"]))
            print(f"   {cat:<22} {best['variant']:<22} SB={best['SB']} TB={best['TB']} MT={best['MT']}"
                  f"  -> {best['med_gbs']} GB/s ({best['pct_peak_bw']}% peak)")


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    force = "--compile" in sys.argv

    print("-- System info --")
    cpu_info = get_cpu_info()
    print(f"  CPU: {cpu_info['name']}")
    print(f"  Cores: {cpu_info['cores']}")

    print("\n-- STREAM copy bandwidth --")
    emp_bw = measure_bandwidth(THREADS)
    print(f"  Empirical BW: {emp_bw:.1f} GB/s")

    print("\n-- Compile --")
    if not compile_kernels(force):
        sys.exit(1)
    has_blas = compile_openblas(force)

    print(f"\n-- Sweep: N={N} reps={REPS} warmup={WARMUP} threads={THREADS or 'auto'} --")
    open(CSV_RAW, "w").close()  # truncate
    sweep_kernels()
    if has_blas:
        sweep_openblas()

    print(f"\n-- Aggregate --")
    rows = aggregate(emp_bw)
    if rows:
        report(cpu_info, emp_bw, rows)

    print(f"\nPer-iteration: {CSV_RAW}\nAggregated:    {CSV_AGG}")