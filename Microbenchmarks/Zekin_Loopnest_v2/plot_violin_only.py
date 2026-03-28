#!/usr/bin/env python3
"""
plot_violins_only.py
Only the 2x2 violin bandwidth plots. Annotates median as % of STREAM peak.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd, numpy as np, argparse

# ---- CONFIG ----
GPU_AMD_CSV = "z_v_grad_w_gpu_beverin.csv"
GPU_NV_CSV  = "z_v_grad_w_gpu_daint.csv"
CPU_AMD_CSV = "z_v_grad_w_cpu_beverin.csv"
CPU_NV_CSV  = "z_v_grad_w_cpu_daint.csv"

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
OUT_STEM = f"violins_nlev{NLEV}"

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
DIST_LABEL = {"uniform": "Uniform", "normal_var1": r"Normal ($\sigma^2 = 1$)"}

parser = argparse.ArgumentParser()
parser.add_argument('--add-peak', action='store_true')
args = parser.parse_args()

def remove_outliers(vals, k=3.0):
    if len(vals) < 4: return vals
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
    if sub.empty: return None, sub
    sub["bandwidth"] = compute_bandwidth(sub)
    medians = sub.groupby("parallelization")["bandwidth"].median()
    if medians.empty: return None, sub
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
        if sub.empty: return fval, np.array([])
        sub["bandwidth"] = compute_bandwidth(sub)
        return fval, sub["bandwidth"].values

def main():
    raw = {}
    for row in GRID:
        for label, csv, fcol, fval, _ in row:
            if label not in raw:
                try:
                    raw[label] = (pd.read_csv(csv), fcol, fval)
                except FileNotFoundError:
                    print(f"  [WARN] {csv} not found"); return

    plt.rcParams.update({
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
    })

    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    fig.suptitle("Semi-Structured Weather Stencil Pattern from ICON", fontsize=15)

    for ri, row_data in enumerate(GRID):
        for ci, (label, csv, fcol, fval, _) in enumerate(row_data):
            ax = axes[ri, ci]
            df_full = raw[label][0]
            positions, data_all, col_all = [], [], []
            medians_for_pct = []   # (pos, median_bw, variant)
            xticks, xlabels = [], []

            for di, dist in enumerate(DISTS):
                for vi, V in enumerate(VARIANTS):
                    best_par, vals = get_best_data(df_full, label, fcol, fval, V, dist, NLEV)
                    vals = remove_outliers(vals)
                    pos = di * 3 + vi
                    if len(vals) > 0:
                        data_all.append(vals)
                        positions.append(pos)
                        col_all.append(VCOL[V])
                        medians_for_pct.append((pos, np.median(vals), V))
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
                ax.axhline(y=peak, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
                ax.set_ylim(bottom=0)
                x_min, x_max = ax.get_xlim()
                ax.text(x_max - 0.02 * (x_max - x_min), peak * 0.98,
                        f'STREAM Peak: {peak:.1f} TB/s',
                        ha='right', va='top', fontsize=9, color='red')

            # Annotate each violin with "X% of Peak"
            if label in STREAM_PEAK:
                # Increase y lim a bit
                ylo, yhi = ax.get_ylim()
                ax.set_ylim(bottom=0,top=yhi * 1.05)
                # Plot line
                peak = STREAM_PEAK[label]
                ylo, yhi = ax.get_ylim()
                # Increase space a bit more
                ax.set_ylim(bottom=0,top=yhi * 1.025)
                ax.axhline(y=yhi, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
                pct = 100.0 * yhi / peak
                pct_rounded = 5 * round(pct / 5)
                ax.text(1.4, yhi - 0.07 * (yhi * 1.2 - ylo), f'{pct_rounded:.0f}% of STREAM Peak',
                        ha='right', va='bottom', fontsize=9, color='gray')
                

    v_handles = [
        Patch(facecolor=VCOL[1], edgecolor="black", label=VLAB[1]),
        Patch(facecolor=VCOL[4], edgecolor="black", label=VLAB[4]),
    ]
    fig.legend(handles=v_handles, loc='lower center',
               bbox_to_anchor=(0.5, -0.02), ncol=2, framealpha=0.9)

    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    sfx = "_w_stream_peak" if args.add_peak else ""
    fig.savefig(f"{OUT_STEM}{sfx}.png", dpi=180, bbox_inches='tight')
    fig.savefig(f"{OUT_STEM}{sfx}.pdf", dpi=180, bbox_inches='tight')

    print("\nBest schedule selection:")
    print(f"{'Platform':<22} {'Variant':<8} {'Distribution':<20} {'Best Schedule':<22} {'Median TB/s':<12} {'% Peak':<8}")
    print("-" * 92)
    for row_data in GRID:
        for label, csv, fcol, fval, _ in row_data:
            if fcol != "parallelization": continue
            df_full = raw[label][0]
            for V in VARIANTS:
                for dist in DISTS:
                    best_par, best_df = pick_best_schedule(df_full, V, dist, NLEV)
                    if best_par and not best_df.empty:
                        med = best_df["bandwidth"].median()
                        pct = 100.0 * med / STREAM_PEAK.get(label, 1.0)
                        print(f"{label:<22} V{V:<7} {dist:<20} {best_par:<22} {med:.4f}      {pct:.1f}%")

    print(f"\nSaved: {OUT_STEM}{sfx}.png, {OUT_STEM}{sfx}.pdf")

if __name__ == "__main__":
    main()