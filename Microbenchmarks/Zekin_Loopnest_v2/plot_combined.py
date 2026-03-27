#!/usr/bin/env python3
"""
plot_best_schedule.py

Top 2x2:  violins of bandwidth per platform.
Bottom 1x2: scatter of cost metrics (mu, delta).
            x = 4 combos: V1/Uniform, V1/Normal, V4/Uniform, V4/Normal
            Bracket annotations show klon-first vs klev-first grouping.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyArrowPatch
from matplotlib.lines import Line2D
import pandas as pd, numpy as np

# ---- CONFIG ----
GPU_AMD_CSV = "z_v_grad_w_gpu_beverin.csv"
GPU_NV_CSV  = "z_v_grad_w_gpu_daint.csv"
CPU_AMD_CSV = "z_v_grad_w_cpu_beverin.csv"
CPU_NV_CSV  = "z_v_grad_w_cpu_daint.csv"
COST_CSV    = "metrics.csv"

STREAM_PEAK = {
    "MI300A Zen Cores":   811.35*1e-3,
    "Grace Neoverse v2": 1607.62*1e-3,
    "MI300A":   3457.5*1e-3,
    "H200":    3720.48*1e-3,
}

NPROMA   = 81920
NLEV     = 96
VARIANTS = [1, 4]
DISTS    = ["uniform", "normal_var1"]
OUT_STEM = f"best_sched_nlev{NLEV}"

BYTES = (2 * NPROMA * 4 +
         2 * NPROMA * 4 +
         NLEV * NPROMA * 8 +
         NLEV * NPROMA * 8 +
         NPROMA * 8 +
         NLEV * NPROMA * 8 +
         NLEV * NPROMA * 8 +
         NPROMA * 8 +
         NPROMA * 8 +
         NLEV * NPROMA * 8)

GRID = [
    [("MI300A Zen Cores", CPU_AMD_CSV, "parallelization", None,        "cpu_scalar"),
     ("MI300A",           GPU_AMD_CSV, "config_label",    "1x1_32x16", "gpu_scalar")],
    [("Grace Neoverse v2", CPU_NV_CSV, "parallelization", None,        "cpu_scalar"),
     ("H200",              GPU_NV_CSV, "config_label",    "1x1_32x16", "gpu_scalar")],
]

LOOP_FOR_V = {1: "klon_first", 4: "klev_first"}
VCOL = {1: "#e67e22", 4: "#2980b9"}
VLAB = {1: "Klon-first Layout", 4: "Klev-first Layout"}

DIST_LABEL = {
    "uniform": "Uniform",
    "normal_var1": r"Normal ($\sigma^2 = 1$)",
}

PLAT_STYLE = {
    "MI300A Zen Cores":  {"color": "#e67e22", "marker": "o"},
    "Grace Neoverse v2": {"color": "#2ecc71", "marker": "s"},
    "MI300A":            {"color": "#e74c3c", "marker": "D"},
    "H200":              {"color": "#3498db", "marker": "^"},
}

SCHED_SHORT = {
    "omp_for":           "for",
    "omp_collapse2":     "col2",
    "opt_for":           "opt-for",
    "opt_collapse_tile": "opt-col",
}


COMPUTE_VARIANTS = [
    ("BlockWidth:8, ComputeWidth:1",  "cpu_scalar", 64, 1,   "#2c3e50", "o"),
    ("BlockWidth:8, ComputeWidth:8",  "cpu_avx512", 64, 8,   "#e67e22", "s"),
    ("BlockWidth:16, ComputeWidth:1", "gpu_scalar", 128, 1,  "#2ecc71", "D"),
    ("BlockWidth:16, ComputeWidth:16","gpu_warp32", 128, 32, "#3498db", "^"),
    ("BlockWidth:16, ComputeWidth:32","gpu_wave64", 128, 64, "#e74c3c", "P"),
]

COST_COMBOS = [
    (1, "uniform",    "Uniform"),
    (1, "normal_var1", r"Normal ($\sigma^2 = 1$)"),
    (4, "uniform",    "Uniform"),
    (4, "normal_var1", r"Normal ($\sigma^2 = 1$)"),
]

BRACKET_GROUPS = [
    (0, 1, "Klon-First"),
    (2, 3, "Klev-First"),
]

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--add-peak', action='store_true',
                    help='Add stream peak lines and set ylim bottom=0')
args = parser.parse_args()

# ---- helpers ----
def remove_outliers(vals, k=3.0):
    if len(vals) < 4:
        return vals
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    clean = vals[(vals >= lo) & (vals <= hi)]
    return clean if len(clean) > 2 else vals

def compute_bandwidth(df):
    return BYTES / (df["time_ms"] * 1e-3) / 1e12

def pick_best_schedule(df_full, V, dist, nlev):
    sub = df_full[(df_full["variant"] == V) &
                  (df_full["cell_dist"] == dist) &
                  (df_full["nlev"] == nlev)].copy()
    if sub.empty:
        return None, sub
    sub["bandwidth"] = compute_bandwidth(sub)
    medians = sub.groupby("parallelization")["bandwidth"].median()
    if medians.empty:
        return None, sub
    best_par = medians.idxmax()
    return best_par, sub[sub["parallelization"] == best_par]

def get_best_data(df_full, label, fcol, fval, V, dist, nlev):
    if fcol == "parallelization":
        best_par, best_df = pick_best_schedule(df_full, V, dist, nlev)
        return best_par, best_df["bandwidth"].values if not best_df.empty else np.array([])
    else:
        sub = df_full[(df_full[fcol] == fval) &
                      (df_full["variant"] == V) &
                      (df_full["cell_dist"] == dist) &
                      (df_full["nlev"] == nlev)].copy()
        if sub.empty:
            return fval, np.array([])
        sub["bandwidth"] = compute_bandwidth(sub)
        return fval, sub["bandwidth"].values

def draw_bracket(ax, x_start, x_end, label, y_frac=-0.16, height=0.06):
    """Draw a horizontal curly brace below the axis with a centered label."""
    inv = ax.transData + ax.transAxes.inverted()
    x0 = inv.transform((x_start, 0))[0]
    x1 = inv.transform((x_end, 0))[0]
    xm = (x0 + x1) / 2.0
    w = x1 - x0

    # Parametric curly brace: left half then right half
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches

    # Control points for curly brace (horizontal, opening downward)
    y_top = y_frac
    y_bot = y_frac - height
    y_mid = (y_top + y_bot) / 2.0
    q = w / 4.0  # quarter width

    verts = [
        (x0, y_top),
        (x0, y_mid),           # left arm down
        (x0, y_bot),
        (xm - 0.01, y_bot),   # to center tip
        (xm, y_bot - height * 0.5),  # tip
        (xm + 0.01, y_bot),   # from center tip
        (x1, y_bot),
        (x1, y_mid),          # right arm up
        (x1, y_top),
    ]
    codes = [
        mpath.Path.MOVETO,
        mpath.Path.CURVE3,
        mpath.Path.CURVE3,
        mpath.Path.CURVE3,
        mpath.Path.CURVE3,
        mpath.Path.CURVE3,
        mpath.Path.CURVE3,
        mpath.Path.CURVE3,
        mpath.Path.CURVE3,
    ]
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor="none", edgecolor="black",
                                lw=1.3, transform=ax.transAxes, clip_on=False)
    ax.add_patch(patch)

    ax.text(xm, y_bot - height * 0.5 - 0.04, label,
            transform=ax.transAxes, ha="center", va="top",
            fontsize=10)

# ---- main ----
def main():
    raw = {}
    have_runtime = True
    for row in GRID:
        for label, csv, fcol, fval, _ in row:
            if label not in raw:
                try:
                    raw[label] = (pd.read_csv(csv), fcol, fval)
                except FileNotFoundError:
                    print(f"  [WARN] {csv} not found, skipping runtime plots")
                    have_runtime = False

    cost = pd.read_csv(COST_CSV)

    plt.rcParams.update({
        "font.size": 14,
        "axes.titlesize": 15,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
    })

    if have_runtime:
        fig = plt.figure(figsize=(9, 10))
        gs_top = fig.add_gridspec(2, 2, top=0.92, bottom=0.51,
                                   hspace=0.35, wspace=0.35,
                                   left=0.10, right=0.96)
        gs_bot = fig.add_gridspec(1, 2, top=0.40, bottom=0.24,
                                   wspace=0.35,
                                   left=0.10, right=0.96)
        fig.suptitle("Semi-Structured Weather Stencil Pattern from ICON",
                     fontsize=15, )
    else:
        fig = plt.figure(figsize=(9, 4.5))
        gs_bot = fig.add_gridspec(1, 2, top=0.80, bottom=0.30,
                                   wspace=0.35,
                                   left=0.10, right=0.96)
        fig.suptitle("Cost Metrics", fontsize=15, )

    # =========================================================
    #  TOP 2x2: violin plots
    # =========================================================
    if have_runtime:
        for ri, row_data in enumerate(GRID):
            for ci, (label, csv, fcol, fval, _) in enumerate(row_data):
                ax = fig.add_subplot(gs_top[ri, ci])
                df_full = raw[label][0]
                positions, data_all, col_all = [], [], []
                xticks, xlabels = [], []

                for di, dist in enumerate(DISTS):
                    for vi, V in enumerate(VARIANTS):
                        best_par, vals = get_best_data(df_full, label, fcol,
                                                        fval, V, dist, NLEV)
                        vals = remove_outliers(vals)
                        pos = di * 3 + vi
                        if len(vals) > 0:
                            data_all.append(vals)
                            positions.append(pos)
                            col_all.append(VCOL[V])
                    xticks.append(di * 3 + 0.5)
                    xlabels.append(DIST_LABEL[dist])

                if data_all:
                    parts = ax.violinplot(data_all, positions=positions,
                                          showmeans=True, showmedians=True,
                                          showextrema=False, widths=0.8)
                    for i, body in enumerate(parts["bodies"]):
                        body.set_facecolor(col_all[i])
                        body.set_edgecolor("black")
                        body.set_alpha(0.75)
                    parts["cmeans"].set_color("black")
                    parts["cmedians"].set_color("white")

                ax.set_xticks(xticks)
                ax.set_xticklabels(xlabels)
                ax.set_ylabel("Bandwidth [TB/s]")
                from matplotlib.ticker import MaxNLocator
                ax.yaxis.set_major_locator(MaxNLocator(nbins=6, min_n_ticks=5))
                ax.margins(y=0.05)
                ax.set_title(label, )
                ax.grid(axis="y", alpha=0.3)
                # --- NEW: add peak line and ylim(bottom=0) ---
                print(f"Label: {label}, add_peak: {args.add_peak}")
                print(STREAM_PEAK)
                if args.add_peak and label in STREAM_PEAK:
                    peak = STREAM_PEAK[label]
                    ax.axhline(y=peak, color='red', linestyle='--',
                               linewidth=1.5, alpha=0.7,
                               label=f'Peak {peak:.1f} TB/s')
                    ax.set_ylim(bottom=0)
                    x_min, x_max = ax.get_xlim()
                    x_text = x_max - 0.02 * (x_max - x_min)  # small offset from right edge
                    y_text = peak * 0.98
                    ax.text(x_text, y_text, f'STREAM Peak: {peak:.1f} TB/s',
                            ha='right', va='top', fontsize=9, color='red')
        # Violin legend below top grid
        v_handles = [
            Patch(facecolor=VCOL[1], edgecolor="black", label=VLAB[1]),
            Patch(facecolor=VCOL[4], edgecolor="black", label=VLAB[4]),
        ]
        fig.legend(handles=v_handles,
                   loc='lower center', bbox_to_anchor=(0.5, 0.425),
                   ncol=2, framealpha=0.9)


    # =========================================================
    #  BOTTOM 1x2: cost metric scatter
    # =========================================================
    n_combos = len(COST_COMBOS)
    x_positions = np.arange(n_combos)
    x_labels = [lbl for _, _, lbl in COST_COMBOS]
    sep_x = 1.5

    n_cv = len(COMPUTE_VARIANTS)
    offsets = np.linspace(-0.25, 0.25, n_cv)

    metrics_to_plot = [
        ("delta", r"Avg. Block Distance $\Delta$"),
        ("mu",    r"Avg. New-Block Count $\mu$"),
    ]

    # Explicit x-axis mapping
    x_map = {
        ("klon_first", "uniform"): 0,
        ("klon_first", "normal_var1"): 1,
        ("klev_first", "uniform"): 2,
        ("klev_first", "normal_var1"): 3,
    }

    xticks = [0, 1, 2, 3]
    xlabels = [
        "Uniform",
        "Normal\n($\\sigma^2 = 1$)",
        "Uniform",
        "Normal\n($\\sigma^2 = 1$)",
    ]

    # offsets for compute variants (visual separation)
    n_cv = len(COMPUTE_VARIANTS)
    offsets = np.linspace(-0.25, 0.25, n_cv)

    metrics_to_plot = [
        ("delta", r"Avg. Block Distance $\Delta$"),
        ("mu",    r"Avg. New-Block Count $\mu$"),
    ]

    for mi, (mcol, mlab) in enumerate(metrics_to_plot):
        ax = fig.add_subplot(gs_bot[0, mi])

        for cvi, (cv_label, cv_target, cv_block, cv_width, cv_color, cv_marker) in enumerate(COMPUTE_VARIANTS):
            for lay, var in [("klon_first", 1), ("klev_first", 4)]:
                xs, ys = [], []
                V = var
                for (_, dist, _) in COST_COMBOS:   # ignore V in the tuple
                    key = (lay, dist)
                    if key not in x_map:
                        continue
                    row = cost[
                        (cost["variant"] == V) &
                        (cost["cell_dist"] == dist) &
                        (cost["target"] == cv_target) &
                        (cost["nlev"] == NLEV) &
                        (cost["loop_order"] == lay) &
                        (cost["block_bytes"] == cv_block) &
                        (cost["vector_width"] == cv_width)
                    ]
                    # now row should be exactly one entry
                    if len(row) == 0:
                        continue
                    print(V, key, cv_label, row[mcol].values[0])

                    assert len(row) == 1, f"Non-unique row for {key}, {cv_target}"
                    x = x_map[key] + offsets[cvi]
                    y = row[mcol].values[0]
                    print(x, y)
                    xs.append(x)
                    ys.append(y)
                    

                if len(xs) > 0:
                    ax.scatter(
                        xs, ys,
                        c=[cv_color],
                        marker=cv_marker,
                        s=65,
                        alpha=0.85,
                        edgecolors="black",
                        linewidths=0.4,
                        label=cv_label,
                        zorder=3
                    )

        # separator
        ax.axvline(x=1.5, color="gray", linestyle="--",
                linewidth=0.8, alpha=0.6)

        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)

        ax.set_title(mlab, fontsize=13)
        ax.set_ylabel("Metric Value")
        ax.set_yscale("log")

        ax.margins(y=0.05)
        ax.grid(axis="y", alpha=0.3, which="both")
        ax.set_xlim(-0.5, 3.5)

        # brackets
        for start, end, blabel in BRACKET_GROUPS:
            draw_bracket(ax, start-0.28, end+0.28, blabel, y_frac=-0.3, height=0.04)

   # Single compute-variant legend below bottom subplots
    cv_handles = [
        Line2D([0], [0], marker=m, color="none",
               markerfacecolor=c, markeredgecolor="black",
               markersize=9, label=lab)
        for lab, _, _, _, c, m in COMPUTE_VARIANTS
    ]
    fig.legend(handles=cv_handles,
               loc='lower center', bbox_to_anchor=(0.52, 0.06),
               ncol=2, )

    fig.tight_layout()
    if args.add_peak:
        PPref = "_w_stream_peak"
    else:
        PPref = ""
    fig.savefig(f"{OUT_STEM}{PPref}.png", dpi=180)
    fig.savefig(f"{OUT_STEM}{PPref}.pdf", dpi=180)

    if have_runtime:
        print("\nBest schedule selection:")
        print(f"{'Platform':<22} {'Variant':<8} {'Distribution':<20} {'Best Schedule':<22} {'Median TB/s':<12}")
        print("-" * 84)
        for row_data in GRID:
            for label, csv, fcol, fval, _ in row_data:
                if fcol != "parallelization":
                    continue
                df_full = raw[label][0]
                for V in VARIANTS:
                    for dist in DISTS:
                        best_par, best_df = pick_best_schedule(df_full, V, dist, NLEV)
                        if best_par and not best_df.empty:
                            med = best_df["bandwidth"].median()
                            print(f"{label:<22} V{V:<7} {dist:<20} {best_par:<22} {med:.4f}")
        print()

    print(f"Saved: {OUT_STEM}.png, {OUT_STEM}.pdf")

if __name__ == "__main__":
    main()