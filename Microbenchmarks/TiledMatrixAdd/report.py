#!/usr/bin/env python3
"""
report_addition.py
Aggregate per-config stats from layout-conflict addition benchmark CSVs
and print detailed summary tables.

Usage:
    python report_addition.py
    python report_addition.py --cpu cpu_raw.csv --gpu gpu_raw.csv --bw-cpu 1700 --bw-gpu 3720
    python report_addition.py --top 20

CSV formats expected:
    CPU: variant,M,N,tile,nthreads,rep,time_s,bw_gbs,checksum,status  (10 cols)
    GPU: "kernel",M,N,tile_rows,tile_cols,rep,time_ms,bw_gbs,checksum,status  (10 cols, kernel may be quoted)
"""
import argparse, os, sys, re
from collections import defaultdict
import numpy as np

# ── Category sets ──
CPU_NAIVE   = {"row_major", "col_major"}
CPU_TILED   = {"tiled"}
CPU_CONTROL = {"all_rowmajor"}

GPU_NAIVE   = {"direct", "direct_T"}
GPU_TILED   = {"tiled+smem"}
GPU_CONTROL = {"all_rowmajor"}


def classify_gpu_kernel(name):
    name = name.strip().strip('"')
    if name.startswith("tiled+smem") or name.startswith("tiled_smem"):
        return "tiled+smem"
    elif name.startswith("direct_T"):
        return "direct_T"
    elif name.startswith("direct"):
        return "direct"
    elif name.startswith("all_rowmajor"):
        return "all_rowmajor"
    return "unknown"


def extract_gpu_params(name):
    """Extract BX, BY, TX, TY, tile from kernel name string."""
    name = name.strip().strip('"')
    params = {}
    for m in re.finditer(r'(\w+)=(\S+)', name):
        params[m.group(1)] = m.group(2)
    return params


# ══════════════════════════════════════════════════════════════════════
#  Loading & Aggregation
# ══════════════════════════════════════════════════════════════════════

def load_cpu(path, peak_bw=0.0):
    """Load CPU CSV, aggregate per (variant, tile, nthreads)."""
    if not os.path.exists(path):
        print(f"  [SKIP] {path} not found")
        return []

    groups = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("variant"):
                continue
            p = line.split(",")
            if len(p) < 10:
                continue
            try:
                key = (p[0], int(p[1]), int(p[2]), int(p[3]), int(p[4]))
                groups[key].append((float(p[6]), float(p[7]), p[9]))
            except (ValueError, IndexError):
                continue

    if not groups:
        print(f"  [WARN] {path}: no valid data")
        return []

    first = list(groups.keys())[0]
    print(f"  {path}: {len(groups)} configs, CPU format, M={first[1]} N={first[2]}")

    rows = []
    for key, iters in groups.items():
        variant, M, N, tile, nthreads = key
        ts  = np.array([t for t, _, _ in iters])
        gbs = np.array([g for _, g, _ in iters])
        status = iters[0][2]
        med_gbs = float(np.median(gbs))

        rows.append(dict(
            variant=variant, M=M, N=N, tile=tile, nthreads=nthreads,
            config_str=f"T={tile}" if tile > 0 else "-",
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
            _fmt="cpu",
        ))

    rows.sort(key=lambda r: -r["med_gbs"])
    return rows


def load_gpu(path, peak_bw=0.0):
    """Load GPU CSV, aggregate per unique kernel name."""
    if not os.path.exists(path):
        print(f"  [SKIP] {path} not found")
        return []

    groups = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("kernel"):
                continue
            # Handle quoted kernel name
            if line.startswith('"'):
                end_q = line.index('"', 1)
                kernel = line[1:end_q]
                rest = line[end_q+1:].lstrip(",").split(",")
            else:
                parts = line.split(",")
                kernel = parts[0]
                rest = parts[1:]
            if len(rest) < 9:
                continue
            try:
                key = (kernel, int(rest[0]), int(rest[1]),
                       int(rest[2]), int(rest[3]))
                groups[key].append((float(rest[5]) * 1e-3,  # ms -> s
                                    float(rest[6]), rest[8]))
            except (ValueError, IndexError):
                continue

    if not groups:
        print(f"  [WARN] {path}: no valid data")
        return []

    first = list(groups.keys())[0]
    print(f"  {path}: {len(groups)} configs, GPU format, M={first[1]} N={first[2]}")

    rows = []
    for key, iters in groups.items():
        kernel, M, N, tr, tc = key
        ts  = np.array([t for t, _, _ in iters])
        gbs = np.array([g for _, g, _ in iters])
        status = iters[0][2]
        med_gbs = float(np.median(gbs))

        cat = classify_gpu_kernel(kernel)
        params = extract_gpu_params(kernel)
        bx = params.get("BX", "?")
        by = params.get("BY", "?")
        tx = params.get("TX", "?")
        ty = params.get("TY", "?")

        rows.append(dict(
            variant=cat, kernel=kernel,
            M=M, N=N, tile_rows=tr, tile_cols=tc,
            BX=bx, BY=by, TX=tx, TY=ty,
            config_str=f"{bx}x{by} TX={tx} TY={ty} tile={tr}x{tc}",
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
            _fmt="gpu",
        ))

    rows.sort(key=lambda r: -r["med_gbs"])
    return rows


# ══════════════════════════════════════════════════════════════════════
#  Table printing
# ══════════════════════════════════════════════════════════════════════

def print_table(title, rows, peak_bw=0.0):
    if not rows:
        return

    M, N = rows[0]["M"], rows[0]["N"]
    data_bytes = M * N * 8 * 3  # read A + read B + write C, doubles
    is_gpu = rows[0]["_fmt"] == "gpu"
    show_peak = peak_bw > 0

    bw_str = f"  Peak BW: {peak_bw:.1f} GB/s" if show_peak else ""
    print(f"\n{'=' * 120}")
    print(f"  {title}")
    print(f"  M={M}  N={N}  bytes/call={data_bytes / 1e9:.3f} GB (3× fp64 arrays)"
          f"  reps={rows[0]['reps']}{bw_str}")
    print(f"  Layouts: A row-major, B col-major, C row-major")
    print(f"{'=' * 120}")

    if is_gpu:
        hdr = (f" {'variant':<14} {'BXxBY':>6} {'TXxTY':>6} {'tile':>8}"
               f"  {'medGB/s':>8} {'minGB/s':>8} {'maxGB/s':>8}")
    else:
        hdr = (f" {'variant':<14} {'tile':>4} {'thr':>4}"
               f"  {'medGB/s':>8} {'minGB/s':>8} {'maxGB/s':>8}")
    if show_peak:
        hdr += f" {'%peak':>6}"
    hdr += f"  {'medMs':>9} {'p5ms':>9} {'p95ms':>9}  {'stat':>4}"
    print(f"\n{hdr}")
    print(f" {'-' * (len(hdr) - 1)}")

    for r in rows:
        med = r["med_gbs"]
        if is_gpu:
            line = (f" {r['variant']:<14} {r['BX']+' x'+r['BY']:>6}"
                    f" {r['TX']+' x'+r['TY']:>6}"
                    f" {str(r['tile_rows'])+'x'+str(r['tile_cols']):>8}"
                    f"  {med:8.1f} {r['min_gbs']:8.1f} {r['max_gbs']:8.1f}")
        else:
            tile_s = str(r["tile"]) if r["tile"] > 0 else "-"
            line = (f" {r['variant']:<14} {tile_s:>4} {r['nthreads']:>4}"
                    f"  {med:8.1f} {r['min_gbs']:8.1f} {r['max_gbs']:8.1f}")
        if show_peak:
            line += f" {r['pct_peak']:5.1f}%"
        line += (f"  {r['med_time_ms']:9.4f} {r['p5_time_ms']:9.4f}"
                 f" {r['p95_time_ms']:9.4f}  {r['status']}")
        print(line)

    # ── Per-category bests ──
    if is_gpu:
        categories = [
            ("direct (row-coal.)",   GPU_NAIVE,   lambda r: r["variant"] == "direct"),
            ("direct_T (col-coal.)", GPU_NAIVE,   lambda r: r["variant"] == "direct_T"),
            ("tiled+smem",           GPU_TILED,   lambda r: r["variant"] == "tiled+smem"),
            ("all_rowmajor (ctrl)",  GPU_CONTROL, lambda r: r["variant"] == "all_rowmajor"),
        ]
    else:
        categories = [
            ("row_major (B stride)", CPU_NAIVE,   lambda r: r["variant"] == "row_major"),
            ("col_major (A,C str.)", CPU_NAIVE,   lambda r: r["variant"] == "col_major"),
            ("tiled",                CPU_TILED,   lambda r: r["variant"] == "tiled"),
            ("all_rowmajor (ctrl)",  CPU_CONTROL, lambda r: r["variant"] == "all_rowmajor"),
        ]

    print(f"\n  Per-category bests:")
    ctrl_gbs = None
    for cat_name, _, pred in categories:
        subset = [r for r in rows if pred(r)]
        if not subset:
            continue
        best = max(subset, key=lambda x: x["med_gbs"])
        pk = f" ({best['pct_peak']:.1f}% peak)" if show_peak else ""
        print(f"    {cat_name:<26} {best['config_str']:<30} -> {best['med_gbs']:.1f} GB/s{pk}")
        if "ctrl" in cat_name or "control" in cat_name:
            ctrl_gbs = best["med_gbs"]

    # ── Model predictions (CPU only) ──
    if not is_gpu:
        B_eff = 64.0 / 8  # cache line / sizeof(double)
        print(f"\n  Model prediction (B_eff = {B_eff:.0f} doubles/line):")
        print(f"    row_major:  cost = 1 + {B_eff:.0f} + 1 = {2+B_eff:.0f}"
              f"  -> predicted {(2+B_eff)/3:.1f}x vs ideal")
        print(f"    col_major:  cost = {B_eff:.0f} + 1 + {B_eff:.0f} = {1+2*B_eff:.0f}"
              f"  -> predicted {(1+2*B_eff)/3:.1f}x vs ideal")

        # Actual ratios
        tiled_best = [r for r in rows if r["variant"] == "tiled"]
        naive_best = [r for r in rows if r["variant"] == "row_major"]
        col_best   = [r for r in rows if r["variant"] == "col_major"]
        if tiled_best and naive_best:
            tb = max(tiled_best, key=lambda x: x["med_gbs"])["med_gbs"]
            nb = max(naive_best, key=lambda x: x["med_gbs"])["med_gbs"]
            print(f"    row_major actual:  {tb/nb:.1f}x slowdown vs best tiled")
        if tiled_best and col_best:
            tb = max(tiled_best, key=lambda x: x["med_gbs"])["med_gbs"]
            cb = max(col_best,   key=lambda x: x["med_gbs"])["med_gbs"]
            print(f"    col_major actual:  {tb/cb:.1f}x slowdown vs best tiled")

    # ── Tile-size sensitivity (CPU) ──
    if not is_gpu:
        tiled_rows = [r for r in rows if r["variant"] == "tiled"]
        if tiled_rows:
            print(f"\n  Tile-size sensitivity:")
            for r in sorted(tiled_rows, key=lambda x: x["tile"]):
                pk = f" ({r['pct_peak']:.1f}%)" if show_peak else ""
                print(f"    T={r['tile']:<4}  {r['med_gbs']:7.1f} GB/s{pk}")

    print()


# ══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="Report layout-conflict addition benchmark results")
    ap.add_argument("--cpu", default="addition_cpu_raw.csv")
    ap.add_argument("--gpu", default="addition_gpu_raw.csv")
    ap.add_argument("--bw-cpu", type=float, default=0,
                    help="CPU STREAM peak in GB/s (0=skip)")
    ap.add_argument("--bw-gpu", type=float, default=0,
                    help="GPU STREAM peak in GB/s (0=skip)")
    ap.add_argument("--top", type=int, default=0,
                    help="Show only top N rows (0=all)")
    args = ap.parse_args()

    cpu_rows = load_cpu(args.cpu, args.bw_cpu)
    gpu_rows = load_gpu(args.gpu, args.bw_gpu)

    if not cpu_rows and not gpu_rows:
        print("No data found. Check file paths.")
        sys.exit(1)

    if cpu_rows:
        display = cpu_rows[:args.top] if args.top else cpu_rows
        print_table("CPU Addition: Layout Conflict", display, args.bw_cpu)

    if gpu_rows:
        display = gpu_rows[:args.top] if args.top else gpu_rows
        print_table("GPU Addition: Layout Conflict", display, args.bw_gpu)


if __name__ == "__main__":
    main()