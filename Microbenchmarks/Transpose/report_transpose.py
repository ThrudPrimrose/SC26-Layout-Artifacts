#!/usr/bin/env python3
"""
Read transpose_cpu_raw.csv and transpose_raw.csv (GPU),
aggregate per-config stats, and print summary tables.

No benchmarks are run — this is a reporting-only script.

Usage:
    python report_transpose.py [--cpu FILE] [--gpu FILE] [--bw-cpu GB/s] [--bw-gpu GB/s]

Defaults:
    --cpu  transpose_cpu_raw.csv
    --gpu  transpose_raw.csv
    --bw-cpu / --bw-gpu   0  (peak BW for %-of-peak; 0 = skip column)
"""
import argparse, csv, os, sys
from collections import defaultdict
import numpy as np

# ── Library variant names (used to tag rows in the table) ──
LIB_NAMES = {"openblas", "openblas_blk", "openblas_blk_omp",
             "hptt", "hptt_blk", "hptt_blk_omp", "hptt_rm_omp", "hptt_patient",
             "cublas", "cutensor", "cutlass"}

# NUMA suffixes appended by transpose_cpu
NUMA_SUFFIXES = ("_nd", "_nc", "_nt")

def strip_numa(name):
    for sfx in NUMA_SUFFIXES:
        if name.endswith(sfx):
            return name[:-len(sfx)]
    return name


# ══════════════════════════════════════════════════════════════════════
#  Aggregate
# ══════════════════════════════════════════════════════════════════════

def load_and_aggregate(path, peak_bw=0.0):
    """
    Auto-detect CSV format and aggregate.

    CPU raw (11 cols):
        variant, N, TB, SB, MT, nthreads, rep, time_s, gbs, cksum, status

    GPU raw (12 cols):
        variant, N, BX, BY, TX, TY, SB, PAD, rep, time_s, gbs, cksum
    """
    if not os.path.exists(path):
        print(f"  [SKIP] {path} not found")
        return []

    groups = defaultdict(list)
    fmt = None  # will be set on first valid row

    with open(path) as f:
        for line in f:
            p = line.strip().split(",")
            ncols = len(p)

            # Auto-detect format from first parseable row
            if fmt is None:
                if ncols == 12:
                    # GPU: try parsing p[9] as time
                    try:
                        float(p[9])
                        fmt = "gpu"
                    except ValueError:
                        continue
                elif ncols == 11:
                    try:
                        float(p[7])
                        fmt = "cpu"
                    except ValueError:
                        continue
                else:
                    continue

            if fmt == "gpu" and ncols == 12:
                try:
                    time_s = float(p[9])
                    gbs_v  = float(p[10])
                except ValueError:
                    continue
                # key: variant, N, BX, BY, TX, TY, SB, PAD
                key = tuple(p[:8])
                groups[key].append((time_s, gbs_v, "PASS"))

            elif fmt == "cpu" and ncols == 11:
                try:
                    time_s = float(p[7])
                    gbs_v  = float(p[8])
                except ValueError:
                    continue
                # key: variant, N, TB, SB, MT, nthreads
                key = tuple(p[:6])
                groups[key].append((time_s, gbs_v, p[10]))

            # else: skip malformed lines

    if not groups:
        print(f"  [WARN] {path} has no valid data rows")
        return []

    first_key = list(groups.keys())[0]
    N = first_key[1]
    print(f"  {path}: {len(groups)} configs, {fmt.upper()} format, N={N}")

    rows = []
    for key, iters in groups.items():
        gbs = np.array([g for _, g, _ in iters])
        ts  = np.array([t for t, _, _ in iters])
        status = iters[0][2]
        med_gbs = float(np.median(gbs))

        if fmt == "gpu":
            # variant, N, BX, BY, TX, TY, SB, PAD
            label = key[0]
            tb_str = f"{key[2]}x{key[3]}"  # BX x BY as "TB"
            sb_str = key[6]                  # SB
            mt_str = key[7]                  # PAD (repurpose MT column)
            thr_str = f"{key[4]}x{key[5]}"  # TX x TY as "thr"
        else:
            # variant, N, TB, SB, MT, nthreads
            label = key[0]
            tb_str = key[2]
            sb_str = key[3]
            mt_str = key[4]
            thr_str = key[5]

        rows.append(dict(
            variant=label, N=N, TB=tb_str, SB=sb_str, MT=mt_str,
            threads=thr_str,
            reps=len(iters),
            med_time_ms=np.median(ts) * 1000,
            std_time_ms=np.std(ts) * 1000,
            p5_time_ms=np.percentile(ts, 5) * 1000,
            p95_time_ms=np.percentile(ts, 95) * 1000,
            min_gbs=float(np.min(gbs)),
            med_gbs=med_gbs,
            max_gbs=float(np.max(gbs)),
            pct_peak=100 * med_gbs / peak_bw if peak_bw > 0 else None,
            status=status,
            _fmt=fmt,
        ))

    rows.sort(key=lambda r: -r["med_gbs"])
    return rows


# ══════════════════════════════════════════════════════════════════════
#  Print table
# ══════════════════════════════════════════════════════════════════════

def print_table(title, rows, peak_bw=0.0):
    if not rows:
        return

    N = rows[0]["N"]
    bpi = 2.0 * int(N) ** 2 * 4  # read + write, float32

    # reference: best library row-major
    lib_rm = [r for r in rows
              if strip_numa(r["variant"]) in LIB_NAMES and "blk" not in r["variant"]]
    ref_gbs = max((float(r["med_gbs"]) for r in lib_rm), default=None)

    print(f"\n{'=' * 120}")
    print(f"  {title}")
    bw_str = f"  Peak BW: {peak_bw:.1f} GB/s" if peak_bw > 0 else ""
    print(f"  N={N}  bytes/call={bpi / 1e9:.3f} GB (R+W fp32)  "
          f"reps={rows[0]['reps']}{bw_str}")
    print(f"{'=' * 120}")

    show_peak = peak_bw > 0
    is_gpu = rows[0].get("_fmt") == "gpu"

    # Header
    if is_gpu:
        hdr = (f" {'variant':<24} {'BXxBY':>5} {'SB':>4} {'PAD':>3} {'TXxTY':>5}"
               f"  {'medGB/s':>8} {'minGB/s':>8} {'maxGB/s':>8}")
    else:
        hdr = (f" {'variant':<24} {'TB':>3} {'SB':>4} {'MT':>2} {'thr':>3}"
               f"  {'medGB/s':>8} {'minGB/s':>8} {'maxGB/s':>8}")
    if show_peak:
        hdr += f" {'%pkBW':>6}"
    hdr += (f"  {'medMs':>8} {'p5ms':>8} {'p95ms':>8}")
    if ref_gbs:
        hdr += f"  {'vs lib':>7}"
    hdr += f"  {'stat':>4}"
    print(f"\n{hdr}")
    print(f" {'-' * (len(hdr) - 1)}")

    for r in rows:
        med = r["med_gbs"]
        if is_gpu:
            line = (f" {r['variant']:<24} {r['TB']:>5} {r['SB']:>4} {r['MT']:>3}"
                    f" {r['threads']:>5}"
                    f"  {med:8.1f} {r['min_gbs']:8.1f} {r['max_gbs']:8.1f}")
        else:
            line = (f" {r['variant']:<24} {r['TB']:>3} {r['SB']:>4} {r['MT']:>2}"
                    f" {r['threads']:>3}"
                    f"  {med:8.1f} {r['min_gbs']:8.1f} {r['max_gbs']:8.1f}")
        if show_peak:
            line += f" {r['pct_peak']:5.1f}%"
        line += (f"  {r['med_time_ms']:8.4f} {r['p5_time_ms']:8.4f}"
                 f" {r['p95_time_ms']:8.4f}")
        if ref_gbs:
            vs = f"{100 * med / ref_gbs:.0f}%"
            line += f"  {vs:>6}"
        is_lib = strip_numa(r["variant"]) in LIB_NAMES
        line += f"  {r['status']}"
        if is_lib:
            line += " *"
        print(line)

    # ── Summary ──
    kern_rows = [r for r in rows if strip_numa(r["variant"]) not in LIB_NAMES]
    lib_rows  = [r for r in rows if strip_numa(r["variant"]) in LIB_NAMES]

    if lib_rows:
        best = max(lib_rows, key=lambda x: x["med_gbs"])
        pk = f" ({best['pct_peak']:.1f}% peak)" if show_peak else ""
        print(f"\n  Best library : {best['variant']} SB={best['SB']}"
              f" -> {best['med_gbs']:.1f} GB/s{pk}")
    if kern_rows:
        best = max(kern_rows, key=lambda x: x["med_gbs"])
        pk = f" ({best['pct_peak']:.1f}% peak)" if show_peak else ""
        extra = ""
        if ref_gbs:
            extra = f"  ({100 * best['med_gbs'] / ref_gbs:.0f}% of best library)"
        print(f"  Best kernel  : {best['variant']} TB={best['TB']} SB={best['SB']}"
              f" MT={best['MT']} -> {best['med_gbs']:.1f} GB/s{pk}{extra}")

    # ── Per-category bests ──
    if is_gpu:
        categories = {
            "Naive":            lambda r: strip_numa(r["variant"]) == "naive",
            "Blocked":          lambda r: strip_numa(r["variant"]) == "blocked",
            "Shared mem":       lambda r: strip_numa(r["variant"]) == "smem",
            "Smem blocked":     lambda r: strip_numa(r["variant"]) == "smem_blk",
            "Smem padded":      lambda r: strip_numa(r["variant"]) == "smem_pad",
            "Smem swizzle":     lambda r: strip_numa(r["variant"]) == "smem_swiz",
            "Blk swizzle":      lambda r: strip_numa(r["variant"]) == "blk_swiz",
            "cuBLAS":           lambda r: r["variant"].startswith("cublas"),
            "cuTENSOR":         lambda r: r["variant"].startswith("cutensor"),
        }
    else:
        categories = {
            "Row-major naive":  lambda r: strip_numa(r["variant"]) in ("naive", "naive_c2"),
            "Row-major tiled":  lambda r: strip_numa(r["variant"]) in ("tiled", "tiled_c2"),
            "Blocked aligned":  lambda r: r["variant"].startswith("blk_aligned"),
            "Locbuf row-major": lambda r: strip_numa(r["variant"]) in (
                "locbuf", "locbuf_c2", "locbuf_2buf", "locbuf_2buf_c2"),
            "Locbuf blocked":   lambda r: r["variant"].startswith("locbuf_blk"),
            "RM blk schedule":  lambda r: strip_numa(r["variant"]).startswith("rm_blk"),
            "HPTT":             lambda r: r["variant"].startswith("hptt"),
            "OpenBLAS":         lambda r: r["variant"].startswith("openblas"),
        }
    any_cat = False
    for cat, pred in categories.items():
        subset = [r for r in rows if pred(r)]
        if subset:
            if not any_cat:
                print(f"\n  Per-category bests:")
                any_cat = True
            best = max(subset, key=lambda x: x["med_gbs"])
            pk = f" ({best['pct_peak']:.1f}% peak)" if show_peak else ""
            print(f"    {cat:<22} {best['variant']:<24} SB={best['SB']}"
                  f" TB={best['TB']} MT={best['MT']}"
                  f"  -> {best['med_gbs']:.1f} GB/s{pk}")
    print()


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="Report transpose benchmark results from raw CSVs")
    ap.add_argument("--cpu", default="transpose_cpu_raw.csv",
                    help="CPU raw CSV (default: transpose_cpu_raw.csv)")
    ap.add_argument("--gpu", default="transpose_raw.csv",
                    help="GPU raw CSV (default: transpose_raw.csv)")
    ap.add_argument("--bw-cpu", type=float, default=0,
                    help="CPU peak BW in GB/s for %%-of-peak (0=skip)")
    ap.add_argument("--bw-gpu", type=float, default=0,
                    help="GPU peak BW in GB/s for %%-of-peak (0=skip)")
    ap.add_argument("--top", type=int, default=0,
                    help="Show only top N rows per table (0=all)")
    args = ap.parse_args()

    cpu_rows = load_and_aggregate(args.cpu, args.bw_cpu)
    gpu_rows = load_and_aggregate(args.gpu, args.bw_gpu)

    if not cpu_rows and not gpu_rows:
        print("No data found. Check file paths.")
        sys.exit(1)

    if cpu_rows:
        display = cpu_rows[:args.top] if args.top else cpu_rows
        print_table("CPU Transpose", display, args.bw_cpu)

    if gpu_rows:
        display = gpu_rows[:args.top] if args.top else gpu_rows
        print_table("GPU Transpose", display, args.bw_gpu)


if __name__ == "__main__":
    main()