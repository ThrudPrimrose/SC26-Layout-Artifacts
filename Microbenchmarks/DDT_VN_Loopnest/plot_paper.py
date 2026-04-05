#!/usr/bin/env python3
"""
plot_ddt_vn_violins.py
2×2 violin-plot grid for the ddt_vn_apc_pc stencil layout comparison.

Five violins per connectivity distribution:
    Orange = SoA baseline, best of {soa, soa_nf} × {collapse 0/1}
    Green  = Full AoS,     best of {aos, aos_nf} × {collapse 0/1}
    Blue   = Best Grouped / AoSoA
    Purple = Best Blocked SoA (blk8..blk128)
    Red    = Best Tiled (tiled_nf_{64..512}, tiled_hf_{64..512})
"""

import re
import sys
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np

# ═══════════════════════════════════════════════════════════════════════════
#  Defaults
# ═══════════════════════════════════════════════════════════════════════════

GPU_BEVERIN_CSV = "ddt_vn_apc_pc_gpu_sweep.csv"
GPU_DAINT_CSV   = "ddt_vn_apc_pc_gpu_sweep.csv"
CPU_BEVERIN_CSV = "ddt_vn_apc_pc_cpu_sweep.csv"
CPU_DAINT_CSV   = "ddt_vn_apc_pc_cpu_sweep.csv"

NPROMA = 81920

STREAM_PEAK = {
    "MI300A Zen CPU":  1228   * 1e-3,
    "Grace CPU":       1700.62 * 1e-3,
    "MI300A GPU":      4294   * 1e-3,
    "GH200 GPU":       3780   * 1e-3,
}

SUBPLOT_W = 9.5
SUBPLOT_H = 3.8

DISTS = ["uniform", "normal_var1", "normal_var4", "sequential", "exact"]
DIST_LABEL = {
    "uniform":     "Uniform",
    "normal_var1": r"Normal ($\sigma^2\!=\!1$)",
    "normal_var4": r"Normal ($\sigma^2\!=\!4$)",
    "sequential":  "Sequential",
    "exact":       "R02B05",
}

# Five violin categories
VCOL = {
    "soa":   "#e67e22",
    "aos":   "#27ae60",
    "best":  "#2980b9",
    "blk":   "#8e44ad",
    "tiled": "#c0392b",
}
VLAB = {
    "soa":   "SoA Baseline",
    "aos":   "Full AoS",
    "best":  "Best Grouped / AoSoA",
    "blk":   "Best Blocked SoA",
    "tiled": "Best Tiled (NF/HF)",
}
VIOLIN_KEYS = ["soa", "aos", "best", "blk", "tiled"]

# Layout-family membership
SOA_FAMILY   = {"soa"}
AOS_FAMILY   = {"aos"}
GRP_FAMILY   = {"grp", "aosoa16", "aosoa32", "aosoa64"}
BLK_FAMILY   = {"blk8", "blk16", "blk32", "blk64", "blk128"}
TILED_FAMILY = {
    "tiled_nf_64",  "tiled_nf_128",  "tiled_nf_256",  "tiled_nf_512",
    "tiled_hf_64",  "tiled_hf_128",  "tiled_hf_256",  "tiled_hf_512",
}


def _base_layout(layout):
    return re.sub(r"_nf$", "", layout)


def _in_family(layout, family):
    base = _base_layout(layout)
    return base in family or layout in family


# ═══════════════════════════════════════════════════════════════════════════
#  Byte-count model
# ═══════════════════════════════════════════════════════════════════════════

def compute_bytes(N, nlev):
    N_c = N_v = N
    b  = 4 * N * nlev * 8
    b += N * (nlev + 1) * 8
    b += 2 * N_c * nlev * 8
    b += N_v * nlev * 8
    b += 2 * N * 2 * 8
    b += N * 8
    b += 2 * N * 2 * 4
    return b


# ═══════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════

def remove_outliers(vals, k=3.0):
    if len(vals) < 4:
        return vals
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    clean = vals[(vals >= lo) & (vals <= hi)]
    return clean if len(clean) > 2 else vals


def _best_group_bw(sub, group_col, bw_bytes):
    if sub.empty:
        return np.array([])
    sub = sub.copy()
    sub["bandwidth"] = bw_bytes / (sub["time_ms"].values * 1e-3) / 1e12
    medians = sub.groupby(group_col)["bandwidth"].median()
    if medians.empty:
        return np.array([])
    best_key = medians.idxmax()
    return sub.loc[sub[group_col] == best_key, "bandwidth"].values


def _best_group_key(sub, group_col, bw_bytes):
    if sub.empty:
        return None, 0.0
    sub = sub.copy()
    sub["bandwidth"] = bw_bytes / (sub["time_ms"].values * 1e-3) / 1e12
    medians = sub.groupby(group_col)["bandwidth"].median()
    if medians.empty:
        return None, 0.0
    best_key = medians.idxmax()
    return best_key, medians[best_key]


# ═══════════════════════════════════════════════════════════════════════════
#  CPU data loaders
# ═══════════════════════════════════════════════════════════════════════════

def load_cpu(path):
    df = pd.read_csv(path)
    df["cfg"] = df["layout"] + "_c" + df["omp_collapse"].astype(str)
    return df


def _cpu_family(df, family, dist, nlev):
    mask = (df["layout"].apply(lambda l: _in_family(l, family))
            & (df["cell_dist"] == dist)
            & (df["nlev"] == nlev))
    return df[mask]


def cpu_soa(df, dist, nlev, bw_bytes):
    return _best_group_bw(_cpu_family(df, SOA_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_aos(df, dist, nlev, bw_bytes):
    return _best_group_bw(_cpu_family(df, AOS_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_best(df, dist, nlev, bw_bytes):
    return _best_group_bw(_cpu_family(df, GRP_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_blk(df, dist, nlev, bw_bytes):
    return _best_group_bw(_cpu_family(df, BLK_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_tiled(df, dist, nlev, bw_bytes):
    return _best_group_bw(_cpu_family(df, TILED_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_soa_key(df, dist, nlev, bw_bytes):
    return _best_group_key(_cpu_family(df, SOA_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_aos_key(df, dist, nlev, bw_bytes):
    return _best_group_key(_cpu_family(df, AOS_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_best_key(df, dist, nlev, bw_bytes):
    return _best_group_key(_cpu_family(df, GRP_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_blk_key(df, dist, nlev, bw_bytes):
    return _best_group_key(_cpu_family(df, BLK_FAMILY, dist, nlev), "cfg", bw_bytes)

def cpu_tiled_key(df, dist, nlev, bw_bytes):
    return _best_group_key(_cpu_family(df, TILED_FAMILY, dist, nlev), "cfg", bw_bytes)


# ═══════════════════════════════════════════════════════════════════════════
#  GPU data loaders
# ═══════════════════════════════════════════════════════════════════════════

def load_gpu(path):
    df = pd.read_csv(path)

    def _layout_full(row):
        if row["layout"] == "aosoa":
            return f"aosoa{int(row['vec_width'])}"
        return row["layout"]

    df["layout_full"] = df.apply(_layout_full, axis=1)
    df["cfg"] = (df["layout_full"] + "_"
                 + df["block_x"].astype(str) + "x" + df["block_y"].astype(str)
                 + "_" + df["coarsen_x"].astype(str) + "x" + df["coarsen_y"].astype(str))
    return df


def _gpu_family(df, family, dist):
    mask = (df["layout_full"].apply(lambda l: _in_family(l, family))
            & (df["cell_dist"] == dist))
    return df[mask]


def gpu_soa(df, dist, bw_bytes):
    return _best_group_bw(_gpu_family(df, SOA_FAMILY, dist), "cfg", bw_bytes)

def gpu_aos(df, dist, bw_bytes):
    return _best_group_bw(_gpu_family(df, AOS_FAMILY, dist), "cfg", bw_bytes)

def gpu_best(df, dist, bw_bytes):
    return _best_group_bw(_gpu_family(df, GRP_FAMILY, dist), "cfg", bw_bytes)

def gpu_blk(df, dist, bw_bytes):
    return _best_group_bw(_gpu_family(df, BLK_FAMILY, dist), "cfg", bw_bytes)

def gpu_tiled(df, dist, bw_bytes):
    return _best_group_bw(_gpu_family(df, TILED_FAMILY, dist), "cfg", bw_bytes)

def gpu_soa_key(df, dist, bw_bytes):
    return _best_group_key(_gpu_family(df, SOA_FAMILY, dist), "cfg", bw_bytes)

def gpu_aos_key(df, dist, bw_bytes):
    return _best_group_key(_gpu_family(df, AOS_FAMILY, dist), "cfg", bw_bytes)

def gpu_best_key(df, dist, bw_bytes):
    return _best_group_key(_gpu_family(df, GRP_FAMILY, dist), "cfg", bw_bytes)

def gpu_blk_key(df, dist, bw_bytes):
    return _best_group_key(_gpu_family(df, BLK_FAMILY, dist), "cfg", bw_bytes)

def gpu_tiled_key(df, dist, bw_bytes):
    return _best_group_key(_gpu_family(df, TILED_FAMILY, dist), "cfg", bw_bytes)


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="2×2 violin plot for ddt_vn_apc_pc layout comparison")
    parser.add_argument("--add-peak", action="store_true")
    parser.add_argument("--nlev", type=int, default=90)
    parser.add_argument("--gpu-beverin", default=GPU_BEVERIN_CSV)
    parser.add_argument("--gpu-daint",   default=GPU_DAINT_CSV)
    parser.add_argument("--cpu-beverin", default=CPU_BEVERIN_CSV)
    parser.add_argument("--cpu-daint",   default=CPU_DAINT_CSV)
    args = parser.parse_args()

    nlev = args.nlev
    bw_bytes = compute_bytes(NPROMA, nlev)
    out_stem = f"ddt_vn_violins_nlev{nlev}"

    GRID = [
        [("MI300A Zen CPU", "MI300A Zen CPU", args.cpu_beverin, "cpu"),
         ("MI300A GPU",     "MI300A GPU",     args.gpu_beverin, "gpu")],
        [("Grace CPU",      "Grace CPU",      args.cpu_daint,   "cpu"),
         ("GH200 GPU",      "GH200 GPU",      args.gpu_daint,   "gpu")],
    ]

    loaded = {}
    for row in GRID:
        for _, _, csv_path, kind in row:
            if csv_path in loaded:
                continue
            try:
                loaded[csv_path] = (
                    load_cpu(csv_path) if kind == "cpu"
                    else load_gpu(csv_path)
                )
            except FileNotFoundError:
                print(f"  [WARN] {csv_path} not found — skipping")
                loaded[csv_path] = None

    if all(v is None for v in loaded.values()):
        print("No CSV files found. Nothing to plot.")
        sys.exit(1)

    plt.rcParams.update({
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
        "xtick.labelsize": 9, "ytick.labelsize": 12, "legend.fontsize": 10,
    })

    nrows = len(GRID)
    ncols = max(len(row) for row in GRID)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(SUBPLOT_W * ncols, SUBPLOT_H * nrows),
        squeeze=False,
    )
    fig.suptitle(
        "ddt_vn_apc_pc  Stencil — Data Layout Comparison",
        fontsize=15, y=0.975,
    )
    fig.text(
        0.5, 0.935,
        "% annotations relative to STREAM Triad peak bandwidth",
        ha="center", va="top", fontsize=12, color="dimgray",
    )

    n_violins = len(VIOLIN_KEYS)   # 5
    group_spacing = n_violins + 1  # 5 violins + 1 gap

    cpu_dispatch = {
        "soa": cpu_soa, "aos": cpu_aos, "best": cpu_best,
        "blk": cpu_blk, "tiled": cpu_tiled,
    }
    gpu_dispatch = {
        "soa": gpu_soa, "aos": gpu_aos, "best": gpu_best,
        "blk": gpu_blk, "tiled": gpu_tiled,
    }

    for ri, row_data in enumerate(GRID):
        for ci, (label, stream_key, csv_path, kind) in enumerate(row_data):
            ax = axes[ri, ci]
            df = loaded.get(csv_path)
            if df is None:
                ax.set_title(f"{label}\n(no data)")
                continue

            positions, data_all, col_all = [], [], []
            medians_for_pct = []
            xticks, xlabels = [], []

            for di, dist in enumerate(DISTS):
                for vi, vk in enumerate(VIOLIN_KEYS):
                    if kind == "cpu":
                        vals = cpu_dispatch[vk](df, dist, nlev, bw_bytes)
                    else:
                        vals = gpu_dispatch[vk](df, dist, bw_bytes)

                    vals = remove_outliers(vals)
                    pos = di * group_spacing + vi

                    if len(vals) > 0:
                        data_all.append(vals)
                        positions.append(pos)
                        col_all.append(VCOL[vk])
                        medians_for_pct.append((pos, np.median(vals), vk))

                xticks.append(di * group_spacing + (n_violins - 1) / 2.0)
                xlabels.append(DIST_LABEL[dist])

            # Y-axis
            all_flat = np.concatenate(data_all) if data_all else np.array([0.0])
            max_val = float(np.max(all_flat))
            locator = MaxNLocator(nbins=5, min_n_ticks=5)
            ticks = locator.tick_values(0.0, max_val * 1.14)
            ticks = ticks[ticks >= 0]
            if len(ticks) > 6:
                ticks = ticks[:6]
            top_lim = ticks[-1] * 1.06 if len(ticks) else max_val * 1.15

            # Draw violins
            if data_all:
                parts = ax.violinplot(
                    data_all, positions=positions,
                    showmeans=True, showmedians=True,
                    showextrema=False, widths=0.85,
                )
                for i, body in enumerate(parts["bodies"]):
                    body.set_facecolor(col_all[i])
                    body.set_edgecolor("black")
                    body.set_alpha(0.75)
                parts["cmeans"].set_color("black")
                parts["cmedians"].set_color("white")

            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabels)
            ax.set_yticks(ticks)
            ax.set_ylim(bottom=0, top=top_lim)
            if ci == 0:
                ax.set_ylabel("Bandwidth [TB/s]")
            ax.set_title(label)
            ax.grid(axis="y", alpha=0.3)

            # Optional STREAM peak line
            if args.add_peak and stream_key in STREAM_PEAK:
                peak = STREAM_PEAK[stream_key]
                ax.axhline(y=peak, color="red", linestyle="--",
                           linewidth=1.5, alpha=0.7)
                xlo, xhi = ax.get_xlim()
                ax.text(xhi - 0.02 * (xhi - xlo), peak * 0.98,
                        f"STREAM Peak: {peak:.2f} TB/s",
                        ha="right", va="top", fontsize=9, color="red")

            # STREAM peak corner label
            if stream_key in STREAM_PEAK:
                peak = STREAM_PEAK[stream_key]
                ax.text(0.03, 0.97, f"{peak:.2f} TB/s STREAM Peak",
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=10, color="dimgray")

            # Per-violin % of peak
            if stream_key in STREAM_PEAK:
                peak = STREAM_PEAK[stream_key]
                yhi = ax.get_ylim()[1]
                offset_down = 0.04 * yhi

                offsets = {
                    "soa":   (-0.25, "right"),
                    "aos":   (-0.10, "right"),
                    "best":  ( 0.00, "center"),
                    "blk":   ( 0.10, "left"),
                    "tiled": ( 0.25, "left"),
                }
                for pos, med_bw, vk in medians_for_pct:
                    pct = 100.0 * med_bw / peak
                    xoff, ha = offsets[vk]
                    ax.text(
                        pos + xoff, med_bw - offset_down,
                        f"{pct:.0f}%", ha=ha, va="top",
                        fontsize=9, color=VCOL[vk], fontweight="bold",
                    )

    # Legend
    handles = [
        Patch(facecolor=VCOL[k], edgecolor="black", label=VLAB[k])
        for k in VIOLIN_KEYS
    ]
    fig.legend(
        handles=handles, loc="lower center",
        bbox_to_anchor=(0.5, -0.01), ncol=5,
        framealpha=0.9, columnspacing=0.6,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.999])

    sfx = "_w_stream_peak" if args.add_peak else ""
    fig.savefig(f"{out_stem}{sfx}.png", dpi=180, bbox_inches="tight")
    fig.savefig(f"{out_stem}{sfx}.pdf", dpi=180, bbox_inches="tight")

    # Summary table
    hdr = (f"{'Platform':<20} {'Violin':<10} {'Distribution':<16} "
           f"{'Best Config':<36} {'Median TB/s':<13} {'% Peak':<8}")
    print(f"\n{hdr}")
    print("─" * len(hdr))

    cpu_key_dispatch = {
        "soa": cpu_soa_key, "aos": cpu_aos_key, "best": cpu_best_key,
        "blk": cpu_blk_key, "tiled": cpu_tiled_key,
    }
    gpu_key_dispatch = {
        "soa": gpu_soa_key, "aos": gpu_aos_key, "best": gpu_best_key,
        "blk": gpu_blk_key, "tiled": gpu_tiled_key,
    }
    violin_names = {
        "soa": "orange", "aos": "green", "best": "blue",
        "blk": "purple", "tiled": "red",
    }

    for row_data in GRID:
        for label_disp, stream_key, csv_path, kind in row_data:
            df = loaded.get(csv_path)
            if df is None:
                continue
            peak = STREAM_PEAK.get(stream_key, 1.0)

            for dist in DISTS:
                for vk in VIOLIN_KEYS:
                    if kind == "cpu":
                        best_cfg, med = cpu_key_dispatch[vk](
                            df, dist, nlev, bw_bytes)
                    else:
                        best_cfg, med = gpu_key_dispatch[vk](
                            df, dist, bw_bytes)

                    if best_cfg is None:
                        continue
                    pct = 100.0 * med / peak
                    print(
                        f"{label_disp:<20} {violin_names[vk]:<10} "
                        f"{dist:<16} {str(best_cfg):<36} "
                        f"{med:<13.4f} {pct:.1f}%"
                    )

    print(f"\nSaved: {out_stem}{sfx}.png, {out_stem}{sfx}.pdf")
    print(f"Byte model: {bw_bytes:,} bytes  (N={NPROMA}, nlev={nlev})")


if __name__ == "__main__":
    main()