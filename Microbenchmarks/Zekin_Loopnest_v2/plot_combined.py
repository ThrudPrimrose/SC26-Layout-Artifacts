#!/usr/bin/env python3
"""
plot_best_schedule.py

Top 2x2:  violins of bandwidth per platform.
Bottom 1x2: scatter of cost metrics (mu, delta) with schedule dimension.
            Filled markers = omp_for, open markers = omp_collapse2.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyArrowPatch
from matplotlib.lines import Line2D
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import pandas as pd, numpy as np

# ---- CONFIG ----
GPU_AMD_CSV = "z_v_grad_w_gpu_beverin.csv"
GPU_NV_CSV  = "z_v_grad_w_gpu_daint.csv"
CPU_AMD_CSV = "z_v_grad_w_cpu_beverin.csv"
CPU_NV_CSV  = "z_v_grad_w_cpu_daint.csv"
COST_CSV    = "results_full.csv"

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

VCOL = {1: "#e67e22", 4: "#2980b9"}
VLAB = {1: "Klon-first Layout", 4: "Klev-first Layout"}

DIST_LABEL = {
    "uniform": "Uniform",
    "normal_var1": r"Normal ($\sigma^2 = 1$)",
}

COMPUTE_VARIANTS = [
    ("BlockWidth:8, ComputeWidth:1",  "CPU_scalar", 64, 1,   "#2c3e50", "o"),
    ("BlockWidth:8, ComputeWidth:8",  "CPU_AVX512", 64, 8,   "#e67e22", "s"),
    ("BlockWidth:16, ComputeWidth:1", "GPU_scalar", 128, 1,  "#2ecc71", "D"),
    ("BlockWidth:16, ComputeWidth:16","GPU_halfw",  128, 16, "#3498db", "^"),
    ("BlockWidth:16, ComputeWidth:32","GPU_warp32", 128, 32, "#e74c3c", "P"),
]

# Schedule visual encoding
SCHED_STYLE = {
    "omp_for":       {"fillstyle": "full",  "label": "omp_for"},
    "omp_collapse2": {"fillstyle": "none",  "label": "collapse2"},
}

# x-axis: (layout, dist) → position
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

LOOP_FOR_V = {1: "klon_first", 4: "klev_first"}

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
                  (df_full["cell_dist"] == dist)].copy()
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
    inv = ax.transData + ax.transAxes.inverted()
    x0 = inv.transform((x_start, 0))[0]
    x1 = inv.transform((x_end, 0))[0]
    xm = (x0 + x1) / 2.0

    verts = [
        (x0, y_frac), (x0, y_frac - height/2),
        (x0, y_frac - height), (xm - 0.01, y_frac - height),
        (xm, y_frac - height * 1.5),
        (xm + 0.01, y_frac - height), (x1, y_frac - height),
        (x1, y_frac - height/2), (x1, y_frac),
    ]
    codes = [mpath.Path.MOVETO] + [mpath.Path.CURVE3] * 8
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor="none", edgecolor="black",
                                lw=1.3, transform=ax.transAxes, clip_on=False)
    ax.add_patch(patch)
    ax.text(xm, y_frac - height * 1.5 - 0.04, label,
            transform=ax.transAxes, ha="center", va="top", fontsize=10)


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
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
    })

    if have_runtime:
        fig = plt.figure(figsize=(9, 10))
        gs_top = fig.add_gridspec(2, 2, top=0.92, bottom=0.51,
                                   hspace=0.35, wspace=0.35,
                                   left=0.10, right=0.96)
        gs_bot = fig.add_gridspec(1, 2, top=0.40, bottom=0.24,
                                   wspace=0.35, left=0.10, right=0.96)
        fig.suptitle("Semi-Structured Weather Stencil Pattern from ICON",
                     fontsize=15)
    else:
        fig = plt.figure(figsize=(9, 4.5))
        gs_bot = fig.add_gridspec(1, 2, top=0.80, bottom=0.30,
                                   wspace=0.35, left=0.10, right=0.96)
        fig.suptitle("Cost Metrics", fontsize=15)

    # =========================================================
    #  TOP 2x2: violin plots (unchanged)
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
                ax.set_title(label)
                ax.grid(axis="y", alpha=0.3)

                if args.add_peak and label in STREAM_PEAK:
                    peak = STREAM_PEAK[label]
                    ax.axhline(y=peak, color='red', linestyle='--',
                               linewidth=1.5, alpha=0.7)
                    ax.set_ylim(bottom=0)
                    x_min, x_max = ax.get_xlim()
                    ax.text(x_max - 0.02 * (x_max - x_min), peak * 0.98,
                            f'STREAM Peak: {peak:.1f} TB/s',
                            ha='right', va='top', fontsize=9, color='red')

        v_handles = [
            Patch(facecolor=VCOL[1], edgecolor="black", label=VLAB[1]),
            Patch(facecolor=VCOL[4], edgecolor="black", label=VLAB[4]),
        ]
        fig.legend(handles=v_handles,
                   loc='lower center', bbox_to_anchor=(0.5, 0.425),
                   ncol=2, framealpha=0.9)

    # =========================================================
    #  BOTTOM 1x2: cost metric scatter with schedule dimension
    #
    #  x-axis:  4 positions (klon/uni, klon/norm, klev/uni, klev/norm)
    #  markers: compute variants (shape + color)
    #  fill:    omp_for = filled, omp_collapse2 = open
    # =========================================================

    # Map (layout, dist) → x position
    x_map = {
        ("klon_first", "uniform"): 0,
        ("klon_first", "normal1"): 1,
        ("klev_first", "uniform"): 2,
        ("klev_first", "normal1"): 3,
    }

    xticks = [0, 1, 2, 3]
    xlabels = [
        "Uniform", r"Normal" + "\n" + r"($\sigma^2 = 1$)",
        "Uniform", r"Normal" + "\n" + r"($\sigma^2 = 1$)",
    ]

    n_cv = len(COMPUTE_VARIANTS)
    schedules = ["omp_for", "omp_collapse2"]
    # Offsets: compute variants spread, schedules stacked within
    cv_offsets = np.linspace(-0.30, 0.30, n_cv)
    sched_offsets = {"omp_for": -0.04, "omp_collapse2": 0.04}

    metrics_to_plot = [
        ("delta", r"Avg. Block Distance $\Delta$"),
        ("mu",    r"Avg. New-Block Count $\mu$"),
    ]

    for mi, (mcol, mlab) in enumerate(metrics_to_plot):
        ax = fig.add_subplot(gs_bot[0, mi])

        seen_labels = set()

        for cvi, (cv_label, cv_target, cv_block, cv_width, cv_color, cv_marker) in enumerate(COMPUTE_VARIANTS):
            for sched in schedules:
                sty = SCHED_STYLE[sched]

                for lay, var in [("klon_first", 1), ("klev_first", 4)]:
                    xs, ys = [], []

                    for dist_csv in ["uniform", "normal1"]:
                        key = (lay, dist_csv)
                        if key not in x_map:
                            continue

                        row = cost[
                            (cost["variant"] == f"V{var}") &
                            (cost["dist"] == dist_csv) &
                            (cost["target"] == cv_target) &
                            (cost["schedule"] == sched) &
                            (cost["block_bytes"] == cv_block) &
                            (cost["vec_width"] == cv_width)
                        ]
                        if len(row) == 0:
                            continue
                        if len(row) > 1:
                            row = row.iloc[:1]

                        x = x_map[key] + cv_offsets[cvi] + sched_offsets[sched]
                        y = row[mcol].values[0]
                        if np.isfinite(y) and y > 0:
                            xs.append(x)
                            ys.append(y)

                    if len(xs) > 0:
                        label_key = (cv_label, sched)
                        lbl = None
                        if label_key not in seen_labels:
                            seen_labels.add(label_key)
                            # Only add to legend once per combo

                        fc = cv_color if sty["fillstyle"] == "full" else "none"
                        ax.scatter(
                            xs, ys,
                            facecolors=[fc] * len(xs),
                            edgecolors=[cv_color] * len(xs),
                            marker=cv_marker,
                            s=65,
                            alpha=0.85,
                            linewidths=1.2 if sty["fillstyle"] == "none" else 0.4,
                            zorder=3,
                        )

        # Separator between layouts
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

        for start, end, blabel in BRACKET_GROUPS:
            draw_bracket(ax, start - 0.28, end + 0.28, blabel,
                         y_frac=-0.3, height=0.04)

    # ---- Legends ----
    # Compute variant legend (shape + color)
    cv_handles = [
        Line2D([0], [0], marker=m, color="none",
               markerfacecolor=c, markeredgecolor="black",
               markersize=9, label=lab)
        for lab, _, _, _, c, m in COMPUTE_VARIANTS
    ]
    # Schedule legend (filled vs open)
    sched_handles = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="gray", markeredgecolor="black",
               markersize=8, label="omp_for (filled)"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="none", markeredgecolor="gray",
               markeredgewidth=1.5,
               markersize=8, label="collapse2 (open)"),
    ]

    all_handles = cv_handles + sched_handles
    fig.legend(handles=all_handles,
               loc='lower center', bbox_to_anchor=(0.52, 0.04),
               ncol=3, fontsize=10)

    fig.tight_layout()
    sfx = "_w_stream_peak" if args.add_peak else ""
    fig.savefig(f"{OUT_STEM}{sfx}.png", dpi=180, bbox_inches='tight')
    fig.savefig(f"{OUT_STEM}{sfx}.pdf", dpi=180, bbox_inches='tight')

    # ---- Summary table ----
    if have_runtime:
        print("\nBest schedule selection:")
        fmt = f"{'Platform':<22} {'Variant':<8} {'Distribution':<20} {'Best Schedule':<22} {'Median TB/s':<12}"
        print(fmt)
        print("-" * len(fmt))
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

    print(f"Saved: {OUT_STEM}{sfx}.png, {OUT_STEM}{sfx}.pdf")

if __name__ == "__main__":
    main()