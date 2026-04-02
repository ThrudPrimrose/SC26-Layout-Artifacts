#!/usr/bin/env python3
"""
plot_transpose_violins.py
2×2 violin bandwidth plots for CPU/GPU transpose benchmarks.

Orange  = Row-major naive (best config)
Green   = Library (OpenBLAS / cuTENSOR / hipTensor), best config
Blue    = Best blocked kernel (best variant × config)

Each violin shows the GB/s distribution across repetitions for the
best-performing configuration within that category.

Usage:
    python plot_transpose_violins.py                         # defaults
    python plot_transpose_violins.py --add-peak              # add STREAM line
    python plot_transpose_violins.py --cpu-file my_cpu.csv   # custom paths
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import pandas as pd, numpy as np, argparse, sys
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

OUT_STEM = "transpose_violins"

# ── Variant classification ──
# Row-major schedule: any custom kernel that operates on row-major data
#   (includes smem/swizzle — they optimise the schedule, data stays row-major).
# Library (row-major): library transpose on row-major data only.
# Blocked schedule: custom kernels whose data access is blocked/tiled.
# No library variants in either custom category.

ROW_MAJOR_VARIANTS = {
    "cpu": {"naive", "naive_c2", "tiled", "tiled_c2",
            "locbuf", "locbuf_c2", "locbuf_2buf", "locbuf_2buf_c2"},
    "gpu": {"naive", "smem", "smem_pad", "smem_swiz"},
}

LIBRARY_VARIANTS = {
    "cpu": {"hptt", "hptt_patient", "hptt_rm_omp", "openblas"},
    "gpu": {"cutensor", "hiptensor", "cublas", "hipblas"},
}

BLOCKED_VARIANTS = {
    "cpu": {"blk_naive", "blk_naive_c2", "blk_tiled", "blk_tiled_c2",
            "blk_aligned", "blk_aligned_c2",
            "blk_aligned_mt", "blk_aligned_mt_c2",
            "locbuf_blk", "locbuf_blk_c2",
            "locbuf_blk_mt", "locbuf_blk_mt_c2",
            "rm_blk", "rm_blk_c2", "rm_blk_mt", "rm_blk_mt_c2",
            "hptt_blk", "hptt_blk_omp"},
    "gpu": {"blocked", "smem_blk", "blk_swiz"},
}

VCOL = {"rm": "#e67e22", "lib": "#27ae60", "blk": "#2980b9"}
VLAB = {"rm": "Row-major schedule", "lib": "Library (row-major)", "blk": "Blocked schedule"}

# NUMA suffixes appended by transpose_cpu — strip before classifying
NUMA_SUFFIXES = ("_nd", "_nc", "_nt")

def strip_numa(name):
    """Strip NUMA policy suffix so 'blk_aligned_mt_nd' matches 'blk_aligned_mt'."""
    for sfx in NUMA_SUFFIXES:
        if name.endswith(sfx):
            return name[:-len(sfx)]
    return name

GRID = [
    # (label,        csv_default,             fmt)
    [("MI300A CPU",  "transpose_cpu_raw.csv",  "cpu"),
     ("MI300A GPU",  "transpose_raw.csv",      "gpu")],
    [("Grace CPU",   "transpose_cpu_raw.csv",  "cpu"),
     ("H200 GPU",    "transpose_raw.csv",      "gpu")],
]


# ══════════════════════════════════════════════════════════════════════
#  CSV Parsing  (auto-detect CPU 11-col vs GPU 12-col)
# ══════════════════════════════════════════════════════════════════════

def parse_raw_csv(path):
    """
    Return list of dicts with keys: variant, config_key, time_s, gbs.
    config_key is a hashable tuple identifying the parameter combination.
    """
    rows = []
    fmt = None
    with open(path) as f:
        for line in f:
            p = line.strip().split(",")
            ncols = len(p)

            if fmt is None:
                if ncols == 12:
                    try: float(p[9]); fmt = "gpu"
                    except ValueError: continue
                elif ncols == 11:
                    try: float(p[7]); fmt = "cpu"
                    except ValueError: continue
                else:
                    continue

            try:
                if fmt == "gpu" and ncols == 12:
                    rows.append(dict(
                        variant=p[0],
                        config_key=(p[2], p[3], p[4], p[5], p[6], p[7]),  # BX,BY,TX,TY,SB,PAD
                        time_s=float(p[9]),
                        gbs=float(p[10]),
                    ))
                elif fmt == "cpu" and ncols == 11:
                    rows.append(dict(
                        variant=p[0],
                        config_key=(p[2], p[3], p[4], p[5]),  # TB,SB,MT,nthreads
                        time_s=float(p[7]),
                        gbs=float(p[8]),
                    ))
            except (ValueError, IndexError):
                continue

    return rows, fmt


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


def best_config_gbs(rows, variant_set):
    """
    Among all (variant, config_key) combos where variant ∈ variant_set,
    pick the one with the highest median GB/s.
    Return (best_variant, best_config, gbs_array) or None.
    """
    # Group by (variant, config_key)
    groups = defaultdict(list)
    for r in rows:
        if strip_numa(r["variant"]) in variant_set:
            groups[(r["variant"], r["config_key"])].append(r["gbs"])

    if not groups:
        return None

    # Find best by median
    best_key = max(groups, key=lambda k: np.median(groups[k]))
    return best_key[0], best_key[1], np.array(groups[best_key])


# ══════════════════════════════════════════════════════════════════════
#  Plotting
# ══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--add-peak", action="store_true",
                        help="Draw STREAM peak as dashed red line")
    parser.add_argument("--cpu-amd", default=None,
                        help="CPU raw CSV for MI300A (default: transpose_cpu_raw.csv)")
    parser.add_argument("--cpu-nv", default=None,
                        help="CPU raw CSV for Grace (default: same as --cpu-amd)")
    parser.add_argument("--gpu-amd", default=None,
                        help="GPU raw CSV for MI300A (default: transpose_raw.csv)")
    parser.add_argument("--gpu-nv", default=None,
                        help="GPU raw CSV for H200 (default: same as --gpu-amd)")
    args = parser.parse_args()

    # Resolve file paths per subplot
    file_map = {
        "MI300A CPU": args.cpu_amd or "transpose_cpu_raw.csv",
        "MI300A GPU": args.gpu_amd or "transpose_raw.csv",
        "Grace CPU":  args.cpu_nv  or args.cpu_amd or "transpose_cpu_raw.csv",
        "H200 GPU":   args.gpu_nv  or args.gpu_amd or "transpose_raw.csv",
    }

    # Parse all unique files
    parsed_cache = {}
    for label, path in file_map.items():
        if path not in parsed_cache:
            try:
                parsed_cache[path] = parse_raw_csv(path)
            except FileNotFoundError:
                print(f"  [WARN] {path} not found, skipping {label}")
                parsed_cache[path] = ([], None)

    # ── Figure setup ──────────────────────────────────────────────────
    plt.rcParams.update({
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
    })

    nrows = len(GRID)
    ncols = max(len(row) for row in GRID)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(SUBPLOT_W * ncols, SUBPLOT_H * nrows),
                             squeeze=False)

    fig.suptitle("Matrix Transpose: Row-major vs Library vs Blocked Schedule",
                 fontsize=15, y=0.98)
    fig.text(0.5, 0.945,
             "% annotations relative to STREAM peak bandwidth",
             ha="center", va="top", fontsize=12, color="dimgray")

    violin_keys = ["rm", "lib", "blk"]
    violin_labels = ["Row-major\nschedule", "Library\n(row-major)", "Blocked\nschedule"]

    for ri, row_data in enumerate(GRID):
        for ci, (label, csv_default, fmt) in enumerate(row_data):
            ax = axes[ri, ci]
            path = file_map[label]
            all_rows, detected_fmt = parsed_cache[path]

            if not all_rows:
                ax.set_title(label)
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=14, color="gray")
                continue

            actual_fmt = detected_fmt or fmt
            rm_set  = ROW_MAJOR_VARIANTS.get(actual_fmt, set())
            lib_set = LIBRARY_VARIANTS.get(actual_fmt, set())
            blk_set = BLOCKED_VARIANTS.get(actual_fmt, set())

            # Get best config per category
            results = {}
            for vk, vset in [("rm", rm_set), ("lib", lib_set), ("blk", blk_set)]:
                res = best_config_gbs(all_rows, vset)
                if res is not None:
                    results[vk] = res  # (variant, config, gbs_array)

            if not results:
                ax.set_title(label)
                ax.text(0.5, 0.5, "no matching variants",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=12, color="gray")
                continue

            # ── Build violin data ─────────────────────────────────────
            positions, data_all, col_all = [], [], []
            medians_info = []
            xticks, xlabels = [], []
            pos = 0

            for vi, vk in enumerate(violin_keys):
                if vk not in results:
                    pos += 1
                    continue
                var_name, cfg, gbs_arr = results[vk]
                gbs_arr = remove_outliers(gbs_arr)
                if len(gbs_arr) == 0:
                    pos += 1
                    continue

                data_all.append(gbs_arr)
                positions.append(pos)
                col_all.append(VCOL[vk])
                medians_info.append((pos, np.median(gbs_arr), vk, var_name, cfg))
                xticks.append(pos)
                # Build x-label: category + variant name
                xlabels.append(f"{violin_labels[vi]}\n({var_name})")
                pos += 1

            # ── Y-axis ticks ──────────────────────────────────────────
            all_flat = np.concatenate(data_all) if data_all else np.array([0.0])
            max_val = float(np.max(all_flat))

            locator = MaxNLocator(nbins=5, min_n_ticks=5)
            ticks = locator.tick_values(0.0, max_val * 1.14)
            ticks = ticks[ticks >= 0]
            if len(ticks) > 6:
                ticks = ticks[:6]
            top_lim = ticks[-1] * 1.06 if len(ticks) else max_val * 1.15

            # ── Draw violins ─────────────────────────────────────────
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

            # ── Axes ─────────────────────────────────────────────────
            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabels, fontsize=10)
            ax.set_yticks(ticks)
            ax.set_ylim(bottom=0, top=top_lim)
            if ci == 0:
                ax.set_ylabel("Bandwidth [GB/s]")
            ax.set_title(label)
            ax.grid(axis="y", alpha=0.3)

            # ── STREAM peak line ─────────────────────────────────────
            peak_gbs = STREAM_PEAK_GBS.get(label)
            if peak_gbs:
                ax.text(0.03, 0.97,
                        f"{peak_gbs:.0f} GB/s STREAM Peak",
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=10, color="dimgray")

                if args.add_peak and peak_gbs <= top_lim * 1.1:
                    ax.axhline(y=peak_gbs, color="red", linestyle="--",
                               linewidth=1.5, alpha=0.7)

            # ── % of peak annotations ────────────────────────────────
            if peak_gbs:
                ylo, yhi = ax.get_ylim()
                offset = 0.045 * (yhi - ylo)
                for p, med, vk, vname, cfg in medians_info:
                    pct = 100.0 * med / peak_gbs
                    ax.text(p, med - offset,
                            f"{pct:.0f}%",
                            ha="center", va="top",
                            fontsize=12, color=VCOL[vk], fontweight="bold")

    # ── Legend ────────────────────────────────────────────────────────
    handles = [Patch(facecolor=VCOL[k], edgecolor="black", label=VLAB[k])
               for k in violin_keys]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.01), ncol=3,
               framealpha=0.9, columnspacing=1.0)

    fig.tight_layout(rect=[0, 0.05, 1, 0.999])
    sfx = "_w_stream_peak" if args.add_peak else ""
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT_STEM}{sfx}.{ext}", dpi=180, bbox_inches="tight")

    # ── Summary table ─────────────────────────────────────────────────
    print(f"\n{'Category':<16} {'Platform':<14} {'Variant':<20} {'Config':<28} "
          f"{'Med GB/s':>9} {'%Peak':>6}")
    print("-" * 95)
    for ri, row_data in enumerate(GRID):
        for ci, (label, _, fmt) in enumerate(row_data):
            path = file_map[label]
            all_rows, detected_fmt = parsed_cache[path]
            if not all_rows:
                continue
            actual_fmt = detected_fmt or fmt
            for vk, vset in [("rm", ROW_MAJOR_VARIANTS.get(actual_fmt, set())),
                             ("lib", LIBRARY_VARIANTS.get(actual_fmt, set())),
                             ("blk", BLOCKED_VARIANTS.get(actual_fmt, set()))]:
                res = best_config_gbs(all_rows, vset)
                if res is None:
                    continue
                var_name, cfg, gbs_arr = res
                med = float(np.median(gbs_arr))
                peak = STREAM_PEAK_GBS.get(label, 0)
                pct = f"{100 * med / peak:.1f}%" if peak else "—"
                cfg_str = ",".join(str(c) for c in cfg)
                print(f"{VLAB[vk]:<16} {label:<14} {var_name:<20} {cfg_str:<28} "
                      f"{med:9.1f} {pct:>6}")

    print(f"\nSaved: {OUT_STEM}{sfx}.png, {OUT_STEM}{sfx}.pdf")


if __name__ == "__main__":
    main()