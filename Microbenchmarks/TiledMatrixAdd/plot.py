#!/usr/bin/env python3
"""
plot_addition_violins.py
2×2 violin bandwidth plots for layout-conflict matrix-addition benchmark.

Orange  = Best naive schedule (no tiling) — exposes layout penalty
Blue    = Best tiled / shared-memory schedule — recovers locality
Green   = All row-major control — layout-free peak

Each violin shows the GB/s distribution across all 100 repetitions for
the best-performing configuration within that category.

Usage:
    python plot_addition_violins.py
    python plot_addition_violins.py --add-peak
    python plot_addition_violins.py --cpu-amd cpu_amd.csv --gpu-nv gpu_nv.csv
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import pandas as pd, numpy as np, argparse, sys, re
from collections import defaultdict

# ══════════════════════════════════════════════════════════════════════
#  Config
# ══════════════════════════════════════════════════════════════════════

STREAM_PEAK_GBS = {
    "MI300A CPU":  1160.35,
    "Grace CPU":   1700.62,
    "MI300A GPU":  3457.5,
    "H200 GPU":    3720.48,
}

SUBPLOT_W = 5.5
SUBPLOT_H = 3.8

OUT_STEM = "addition_violins"

# ── Category classification ──
# "naive":   no tiling, directly exposes layout conflict
# "tiled":   tiled + local transpose (stack buf on CPU, smem on GPU)
# "control": all row-major, no layout conflict (peak reference)

CPU_NAIVE    = {"row_major", "col_major"}
CPU_TILED    = {"tiled"}
CPU_CONTROL  = {"all_rowmajor"}

GPU_NAIVE    = {"direct", "direct_T"}
GPU_TILED    = {"tiled+smem"}
GPU_CONTROL  = {"all_rowmajor"}

VCOL = {"naive": "#e67e22", "tiled": "#2980b9", "ctrl": "#27ae60"}
VLAB = {"naive": "Naive schedule", "tiled": "Tiled + local transpose", "ctrl": "All row-major (control)"}

GRID = [
    [("MI300A CPU", "addition_cpu_raw.csv",  "cpu"),
     ("MI300A GPU", "addition_gpu_raw.csv",  "gpu")],
    [("Grace CPU",  "addition_cpu_raw.csv",  "cpu"),
     ("H200 GPU",   "addition_gpu_raw.csv",  "gpu")],
]


# ══════════════════════════════════════════════════════════════════════
#  CSV Parsing
# ══════════════════════════════════════════════════════════════════════

def classify_gpu_kernel(name):
    """Extract category and config key from GPU kernel name.

    Kernel names look like:
        direct        BX=32  BY=8   TX=1  TY=1   tile=8x32
        direct_T      BX=32  BY=8   TX=1  TY=1   tile=8x32
        tiled+smem    BX=16  BY=16  TX=2  TY=2   tile=32x32
        all_rowmajor  BX=32  BY=8   TX=4  TY=4   (control)
    """
    name = name.strip().strip('"')
    if name.startswith("tiled+smem") or name.startswith("tiled_smem"):
        return "tiled+smem", name
    elif name.startswith("direct_T"):
        return "direct_T", name
    elif name.startswith("direct"):
        return "direct", name
    elif name.startswith("all_rowmajor"):
        return "all_rowmajor", name
    return "unknown", name


def parse_cpu_csv(path):
    """Parse CPU CSV: variant,M,N,tile,nthreads,rep,time_s,bw_gbs,checksum,status"""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("variant"):
                continue
            p = line.split(",")
            if len(p) < 10:
                continue
            try:
                rows.append(dict(
                    variant=p[0],
                    M=int(p[1]), N=int(p[2]),
                    tile=int(p[3]),
                    nthreads=int(p[4]),
                    rep=int(p[5]),
                    time_s=float(p[6]),
                    gbs=float(p[7]),
                    checksum=p[8],
                    status=p[9],
                ))
            except (ValueError, IndexError):
                continue
    return rows


def parse_gpu_csv(path):
    """Parse GPU CSV: kernel,M,N,tile_rows,tile_cols,rep,time_ms,bw_gbs,checksum,status
    kernel field may be quoted."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("kernel"):
                continue
            # Handle quoted kernel name
            if line.startswith('"'):
                end_quote = line.index('"', 1)
                kernel = line[1:end_quote]
                rest = line[end_quote+1:].lstrip(",").split(",")
            else:
                parts = line.split(",")
                kernel = parts[0]
                rest = parts[1:]
            if len(rest) < 9:
                continue
            try:
                cat, cfg_key = classify_gpu_kernel(kernel)
                rows.append(dict(
                    kernel=kernel,
                    category=cat,
                    config_key=cfg_key,
                    M=int(rest[0]), N=int(rest[1]),
                    tile_rows=int(rest[2]), tile_cols=int(rest[3]),
                    rep=int(rest[4]),
                    time_ms=float(rest[5]),
                    gbs=float(rest[6]),
                    checksum=rest[7],
                    status=rest[8],
                ))
            except (ValueError, IndexError):
                continue
    return rows


# ══════════════════════════════════════════════════════════════════════
#  Data helpers
# ══════════════════════════════════════════════════════════════════════

def remove_outliers(vals, k=3.0):
    if len(vals) < 4:
        return vals
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    clean = vals[(vals >= lo) & (vals <= hi)]
    return clean if len(clean) > 2 else vals


def best_cpu_config(rows, variant_set):
    """Best (variant, tile) combo by median GB/s among variant_set."""
    groups = defaultdict(list)
    for r in rows:
        if r["variant"] in variant_set:
            groups[(r["variant"], r["tile"])].append(r["gbs"])
    if not groups:
        return None
    best_key = max(groups, key=lambda k: np.median(groups[k]))
    return best_key[0], best_key[1], np.array(groups[best_key])


def best_gpu_config(rows, category_set):
    """Best config_key by median GB/s among category_set."""
    groups = defaultdict(list)
    for r in rows:
        if r["category"] in category_set:
            groups[(r["category"], r["config_key"])].append(r["gbs"])
    if not groups:
        return None
    best_key = max(groups, key=lambda k: np.median(groups[k]))
    gbs_arr = np.array(groups[best_key])
    # Extract a short label
    label = best_key[1]  # full kernel name
    return best_key[0], label, gbs_arr


# ══════════════════════════════════════════════════════════════════════
#  Plotting
# ══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--add-peak", action="store_true",
                        help="Draw STREAM peak as dashed red line")
    parser.add_argument("--cpu-amd", default=None)
    parser.add_argument("--cpu-nv", default=None)
    parser.add_argument("--gpu-amd", default=None)
    parser.add_argument("--gpu-nv", default=None)
    args = parser.parse_args()

    file_map = {
        "MI300A CPU": args.cpu_amd or "addition_cpu_raw.csv",
        "MI300A GPU": args.gpu_amd or "addition_gpu_raw.csv",
        "Grace CPU":  args.cpu_nv  or args.cpu_amd or "addition_cpu_raw.csv",
        "H200 GPU":   args.gpu_nv  or args.gpu_amd or "addition_gpu_raw.csv",
    }

    # Parse
    parsed = {}
    for label, path in file_map.items():
        if path not in parsed:
            try:
                fmt = "gpu" if "GPU" in label else "cpu"
                if fmt == "cpu":
                    parsed[path] = ("cpu", parse_cpu_csv(path))
                else:
                    parsed[path] = ("gpu", parse_gpu_csv(path))
            except FileNotFoundError:
                print(f"  [WARN] {path} not found, skipping {label}")
                parsed[path] = (None, [])

    # ── Figure ────────────────────────────────────────────────────────
    plt.rcParams.update({
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
        "xtick.labelsize": 11, "ytick.labelsize": 12, "legend.fontsize": 12,
    })

    nrows = len(GRID)
    ncols = max(len(row) for row in GRID)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(SUBPLOT_W * ncols, SUBPLOT_H * nrows),
                             squeeze=False)

    fig.suptitle(r"$C = A + B$: Layout Conflict (A row-major, B col-major, C row-major)",
                 fontsize=15, y=0.98)
    fig.text(0.5, 0.945,
             "% annotations relative to STREAM peak bandwidth",
             ha="center", va="top", fontsize=12, color="dimgray")

    vkeys = ["naive", "tiled", "ctrl"]
    vlabels_short = ["Naive\nschedule", "Tiled + local\ntranspose", "All row-major\n(control)"]

    for ri, row_data in enumerate(GRID):
        for ci, (label, csv_default, fmt) in enumerate(row_data):
            ax = axes[ri, ci]
            path = file_map[label]
            detected_fmt, all_rows = parsed.get(path, (None, []))

            if not all_rows:
                ax.set_title(label)
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=14, color="gray")
                continue

            # Find best per category
            results = {}
            if detected_fmt == "cpu":
                for vk, vset in [("naive", CPU_NAIVE), ("tiled", CPU_TILED),
                                 ("ctrl", CPU_CONTROL)]:
                    res = best_cpu_config(all_rows, vset)
                    if res:
                        results[vk] = res
            else:
                for vk, cset in [("naive", GPU_NAIVE), ("tiled", GPU_TILED),
                                 ("ctrl", GPU_CONTROL)]:
                    res = best_gpu_config(all_rows, cset)
                    if res:
                        results[vk] = res

            if not results:
                ax.set_title(label)
                ax.text(0.5, 0.5, "no matching data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=12, color="gray")
                continue

            # ── Build violins ─────────────────────────────────────────
            positions, data_all, col_all = [], [], []
            medians_info = []
            xticks, xlabels = [], []
            pos = 0

            for vi, vk in enumerate(vkeys):
                if vk not in results:
                    pos += 1
                    continue
                vname, cfg, gbs_arr = results[vk]
                gbs_arr = remove_outliers(gbs_arr)
                if len(gbs_arr) == 0:
                    pos += 1
                    continue

                data_all.append(gbs_arr)
                positions.append(pos)
                col_all.append(VCOL[vk])
                medians_info.append((pos, np.median(gbs_arr), vk, vname, cfg))
                xticks.append(pos)
                # Short label with variant info
                if detected_fmt == "cpu":
                    detail = f"({vname}" + (f" T={cfg}" if cfg else "") + ")"
                else:
                    # Shorten GPU config
                    m = re.search(r'tile=(\S+)', str(cfg))
                    tile_str = m.group(1) if m else ""
                    detail = f"({vname}" + (f" {tile_str}" if tile_str else "") + ")"
                xlabels.append(f"{vlabels_short[vi]}\n{detail}")
                pos += 1

            # ── Y-axis ────────────────────────────────────────────────
            all_flat = np.concatenate(data_all) if data_all else np.array([0.0])
            max_val = float(np.max(all_flat))

            locator = MaxNLocator(nbins=5, min_n_ticks=5)
            ticks = locator.tick_values(0.0, max_val * 1.14)
            ticks = ticks[ticks >= 0]
            if len(ticks) > 6:
                ticks = ticks[:6]
            top_lim = ticks[-1] * 1.06 if len(ticks) else max_val * 1.15

            # ── Draw ──────────────────────────────────────────────────
            if data_all:
                parts = ax.violinplot(data_all, positions=positions,
                                      showmeans=True, showmedians=True,
                                      showextrema=False, widths=0.85)
                for i, body in enumerate(parts["bodies"]):
                    body.set_facecolor(col_all[i])
                    body.set_edgecolor("black")
                    body.set_alpha(0.75)
                parts["cmeans"].set_color("black")
                parts["cmedians"].set_color("white")

            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabels, fontsize=9)
            ax.set_yticks(ticks)
            ax.set_ylim(bottom=0, top=top_lim)
            if ci == 0:
                ax.set_ylabel("Bandwidth [GB/s]")
            ax.set_title(label)
            ax.grid(axis="y", alpha=0.3)

            # ── STREAM peak ───────────────────────────────────────────
            peak_gbs = STREAM_PEAK_GBS.get(label)
            if peak_gbs:
                ax.text(0.03, 0.97, f"{peak_gbs:.0f} GB/s STREAM Peak",
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=10, color="dimgray")
                if args.add_peak and peak_gbs <= top_lim * 1.1:
                    ax.axhline(y=peak_gbs, color="red", linestyle="--",
                               linewidth=1.5, alpha=0.7)

            # ── % of peak annotations ─────────────────────────────────
            if peak_gbs:
                ylo, yhi = ax.get_ylim()
                offset = 0.045 * (yhi - ylo)
                for p, med, vk, vname, cfg in medians_info:
                    pct = 100.0 * med / peak_gbs
                    ax.text(p, med - offset, f"{pct:.0f}%",
                            ha="center", va="top",
                            fontsize=12, color=VCOL[vk], fontweight="bold")

    # ── Legend ────────────────────────────────────────────────────────
    handles = [Patch(facecolor=VCOL[k], edgecolor="black", label=VLAB[k])
               for k in vkeys]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.01), ncol=3,
               framealpha=0.9, columnspacing=1.0)

    fig.tight_layout(rect=[0, 0.05, 1, 0.999])
    sfx = "_w_stream_peak" if args.add_peak else ""
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT_STEM}{sfx}.{ext}", dpi=180, bbox_inches="tight")

    # ── Summary table ─────────────────────────────────────────────────
    print(f"\n{'Category':<28} {'Platform':<14} {'Variant':<24} "
          f"{'Config':<20} {'Med GB/s':>9} {'%Peak':>6}")
    print("-" * 105)
    for ri, row_data in enumerate(GRID):
        for ci, (label, _, fmt) in enumerate(row_data):
            path = file_map[label]
            detected_fmt, all_rows = parsed.get(path, (None, []))
            if not all_rows:
                continue
            if detected_fmt == "cpu":
                cats = [("naive", CPU_NAIVE), ("tiled", CPU_TILED), ("ctrl", CPU_CONTROL)]
                find = best_cpu_config
            else:
                cats = [("naive", GPU_NAIVE), ("tiled", GPU_TILED), ("ctrl", GPU_CONTROL)]
                find = best_gpu_config
            for vk, vset in cats:
                res = find(all_rows, vset)
                if not res:
                    continue
                vname, cfg, gbs_arr = res
                med = float(np.median(gbs_arr))
                peak = STREAM_PEAK_GBS.get(label, 0)
                pct = f"{100*med/peak:.1f}%" if peak else "-"
                cfg_s = str(cfg) if not isinstance(cfg, int) else f"T={cfg}"
                print(f"{VLAB[vk]:<28} {label:<14} {vname:<24} "
                      f"{cfg_s:<20} {med:9.1f} {pct:>6}")

    print(f"\nSaved: {OUT_STEM}{sfx}.png, {OUT_STEM}{sfx}.pdf")


if __name__ == "__main__":
    main()