#!/usr/bin/env python3
"""
plot_zaxpy_sweep.py
Violin bandwidth plots for indirect zaxpy benchmark.

Two PNGs: one for "small" (nat20 original sizes), one for "1gb" (tiled).
Columns: AMD (left), NVIDIA (right)
Rows:    GPU only (1 row) or CPU + GPU (2×2) when --cpu-*-csv given.

X-axis groups: Uniform, Normal, QE
Colours = 4 kernel variants:
    aos_scatter (orange), aos_sorted (blue),
    soa_scatter (red),    soa_sorted (green)

Y-axis: effective bandwidth [TB/s], with % of STREAM peak annotated.

Usage:
    python plot_zaxpy_sweep.py \
        --gpu-amd-small-csv amd/zaxpy_sweep_small.csv \
        --gpu-amd-1gb-csv   amd/zaxpy_sweep_1gb.csv   \
        --gpu-nv-small-csv  nv/zaxpy_sweep_small.csv   \
        --gpu-nv-1gb-csv    nv/zaxpy_sweep_1gb.csv
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import pandas as pd, numpy as np, argparse, sys, re

# ── Colour / label scheme ────────────────────────────────────────────────────
VCOL = {
    "aos_scatter": "#e67e22",
    "aos_sorted":  "#2980b9",
    "soa_scatter": "#c0392b",
    "soa_sorted":  "#27ae60",
}
VLAB = {
    "aos_scatter": "AoS  Scatter (original)",
    "aos_sorted":  "AoS  Sorted",
    "soa_scatter": "SoA  Scatter",
    "soa_sorted":  "SoA  Sorted",
}
VARIANT_ORDER = ["aos_scatter", "aos_sorted", "soa_scatter", "soa_sorted"]

# Bytes per element
BYTES_PER_ELEM = {
    "aos_scatter": 52, "aos_sorted": 56,
    "soa_scatter": 52, "soa_sorted": 56,
}

# Distribution keys and how to match experiment names in the CSV
#   small CSV: "qe", "uniform_s0", "normal_s2", ...
#   1gb CSV:   "qe_tiled", "uniform_tiled_s0", "normal_tiled_s1", ...
DIST_ORDER = ["uniform", "normal", "qe"]
DIST_LABEL = {"uniform": "Uniform", "normal": "Normal", "qe": "QE"}

STREAM_PEAK = {
    "MI300A Zen CPU":  1228    * 1e-3,
    "Grace CPU":       1700.62 * 1e-3,
    "MI300A GPU":      4294    * 1e-3,
    "GH200 GPU":       3780    * 1e-3,
}

SUBPLOT_W = 5.5
SUBPLOT_H = 3.5

# ── Helpers ──────────────────────────────────────────────────────────────────

def remove_outliers(vals, k=3.0):
    if len(vals) < 4:
        return vals
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    c = vals[(vals >= lo) & (vals <= hi)]
    return c if len(c) > 2 else vals


def compute_bandwidth(time_ms, nx, variant):
    return nx * BYTES_PER_ELEM[variant] / (time_ms * 1e-3) / 1e12


def best_config_bw(df, variant):
    """
    Best (tpb, coarsen) by highest median BW.
    nx is read per-row from the CSV.
    """
    sub = df[df["variant"] == variant].copy()
    if sub.empty:
        return np.array([])
    sub["bw"] = compute_bandwidth(sub["time_ms"].values,
                                  sub["nx"].values,
                                  variant)
    medians = sub.groupby(["tpb", "coarsen"])["bw"].median()
    if medians.empty:
        return np.array([])
    best = medians.idxmax()
    vals = sub[(sub["tpb"] == best[0]) & (sub["coarsen"] == best[1])]["bw"].values
    return vals


def select_dist(df, dist_key, tiled):
    """
    Filter CSV rows for a given distribution.
      small CSV experiments: "qe", "uniform_s0" .. "uniform_s4", "normal_s0" ..
      1gb CSV experiments:   "qe_tiled", "uniform_tiled_s0" .., "normal_tiled_s0" ..
    For uniform/normal, pool all samples together.
    """
    if dist_key == "qe":
        pat = "qe_tiled" if tiled else "qe"
        return df[df["experiment"] == pat]
    elif dist_key == "uniform":
        prefix = "uniform_tiled_s" if tiled else "uniform_s"
        return df[df["experiment"].str.startswith(prefix)]
    else:  # normal
        prefix = "normal_tiled_s" if tiled else "normal_s"
        return df[df["experiment"].str.startswith(prefix)]


# ── Subplot painter ──────────────────────────────────────────────────────────

def paint_subplot(ax, df, peak_label, tiled, add_peak=False):
    nvars = len(VARIANT_ORDER)
    group_spacing = nvars + 1.5

    positions, data_all, col_all = [], [], []
    medians_annot = []
    xticks, xlabels = [], []

    for di, dist_key in enumerate(DIST_ORDER):
        df_slice = select_dist(df, dist_key, tiled)

        for vi, vk in enumerate(VARIANT_ORDER):
            vals = best_config_bw(df_slice, vk)
            vals = remove_outliers(vals)
            pos = di * group_spacing + vi
            if len(vals) > 0:
                data_all.append(vals)
                positions.append(pos)
                col_all.append(VCOL[vk])
                medians_annot.append((pos, np.median(vals), vk))

        xticks.append(di * group_spacing + (nvars - 1) / 2.0)
        xlabels.append(DIST_LABEL[dist_key])

    # Y-axis
    all_flat = np.concatenate(data_all) if data_all else np.array([0.0])
    mx = float(np.max(all_flat))
    loc = MaxNLocator(nbins=5, min_n_ticks=5)
    ticks = loc.tick_values(0.0, mx * 1.14)
    ticks = ticks[ticks >= 0]
    if len(ticks) > 6: ticks = ticks[:6]
    top = ticks[-1] * 1.06 if len(ticks) else mx * 1.15

    # Violins
    if data_all:
        parts = ax.violinplot(data_all, positions=positions,
                              showmeans=True, showmedians=True,
                              showextrema=False, widths=0.9)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(col_all[i])
            body.set_edgecolor("black")
            body.set_alpha(0.75)
        parts["cmeans"].set_color("black")
        parts["cmedians"].set_color("white")

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticks(ticks)
    ax.set_ylim(bottom=0, top=top)
    ax.grid(axis="y", alpha=0.3)

    # STREAM peak corner label
    if peak_label in STREAM_PEAK:
        peak = STREAM_PEAK[peak_label]
        ax.text(0.03, 0.97, f"{peak:.2f} TB/s STREAM Peak",
                transform=ax.transAxes, ha='left', va='top',
                fontsize=10, color='dimgray')

    # Optional horizontal line
    if add_peak and peak_label in STREAM_PEAK:
        peak = STREAM_PEAK[peak_label]
        ax.axhline(y=peak, color='red', ls='--', lw=1.5, alpha=0.7)
        xmin, xmax = ax.get_xlim()
        ax.text(xmax - 0.02*(xmax-xmin), peak*0.98,
                f'STREAM Peak: {peak:.2f} TB/s',
                ha='right', va='top', fontsize=9, color='red')

    # % annotations
    if peak_label in STREAM_PEAK:
        peak = STREAM_PEAK[peak_label]
        yr = top
        off = 0.045 * yr
        for pos, med, vk in medians_annot:
            pct = 100.0 * med / peak
            ax.text(pos, med - off, f'{pct:.0f}%',
                    ha='center', va='top',
                    fontsize=9, color=VCOL[vk], fontweight='bold')


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()

    # GPU (required)
    ap.add_argument("--gpu-amd-small-csv", required=True)
    ap.add_argument("--gpu-amd-1gb-csv",   required=True)
    ap.add_argument("--gpu-nv-small-csv",  required=True)
    ap.add_argument("--gpu-nv-1gb-csv",    required=True)

    # CPU (optional)
    ap.add_argument("--cpu-amd-small-csv", default=None)
    ap.add_argument("--cpu-amd-1gb-csv",   default=None)
    ap.add_argument("--cpu-nv-small-csv",  default=None)
    ap.add_argument("--cpu-nv-1gb-csv",    default=None)

    ap.add_argument("--add-peak", action="store_true")
    args = ap.parse_args()

    # Load
    def load(p):
        try: return pd.read_csv(p)
        except FileNotFoundError: print(f"WARN: {p} not found"); return None

    gpu_amd_sm = load(args.gpu_amd_small_csv)
    gpu_amd_1g = load(args.gpu_amd_1gb_csv)
    gpu_nv_sm  = load(args.gpu_nv_small_csv)
    gpu_nv_1g  = load(args.gpu_nv_1gb_csv)
    if any(x is None for x in [gpu_amd_sm, gpu_amd_1g, gpu_nv_sm, gpu_nv_1g]):
        print("ERROR: GPU CSVs required."); return

    have_cpu = all(v is not None for v in [
        args.cpu_amd_small_csv, args.cpu_amd_1gb_csv,
        args.cpu_nv_small_csv,  args.cpu_nv_1gb_csv])
    cpu_amd_sm = cpu_amd_1g = cpu_nv_sm = cpu_nv_1g = None
    if have_cpu:
        cpu_amd_sm = load(args.cpu_amd_small_csv)
        cpu_amd_1g = load(args.cpu_amd_1gb_csv)
        cpu_nv_sm  = load(args.cpu_nv_small_csv)
        cpu_nv_1g  = load(args.cpu_nv_1gb_csv)
        if any(x is None for x in [cpu_amd_sm, cpu_amd_1g, cpu_nv_sm, cpu_nv_1g]):
            print("WARN: CPU CSVs incomplete, GPU-only."); have_cpu = False

    plt.rcParams.update({
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
    })

    # ── One figure per scale ─────────────────────────────────────────────
    scales = [
        ("small", False, {
            "gpu_amd": gpu_amd_sm, "gpu_nv": gpu_nv_sm,
            "cpu_amd": cpu_amd_sm, "cpu_nv": cpu_nv_sm,
        }),
        ("1gb", True, {
            "gpu_amd": gpu_amd_1g, "gpu_nv": gpu_nv_1g,
            "cpu_amd": cpu_amd_1g, "cpu_nv": cpu_nv_1g,
        }),
    ]

    for scale_name, tiled, dfs in scales:

        if have_cpu:
            GRID = [
                [("MI300A Zen CPU", dfs["cpu_amd"], "MI300A Zen CPU"),
                 ("Grace CPU",      dfs["cpu_nv"],  "Grace CPU")],
                [("MI300A GPU",     dfs["gpu_amd"], "MI300A GPU"),
                 ("GH200 GPU",      dfs["gpu_nv"],  "GH200 GPU")],
            ]
        else:
            GRID = [
                [("MI300A GPU", dfs["gpu_amd"], "MI300A GPU"),
                 ("GH200 GPU",  dfs["gpu_nv"],  "GH200 GPU")],
            ]

        nrows = len(GRID)
        ncols = 2

        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(SUBPLOT_W*ncols, SUBPLOT_H*nrows),
                                 squeeze=False)

        title_scale = "nat20 original" if not tiled else "~1 GiB tiled"
        fig.suptitle(f"Indirect ZAXPY \u2014 {title_scale}",
                     fontsize=15, y=0.975)
        fig.text(0.5, 0.935,
                 "% annotations relative to STREAM peak bandwidth",
                 ha='center', va='top', fontsize=12, color='dimgray')

        for ri, row in enumerate(GRID):
            for ci, (label, df, peak_label) in enumerate(row):
                ax = axes[ri, ci]
                ax.set_title(label)
                if ci == 0:
                    ax.set_ylabel("Bandwidth [TB/s]")
                paint_subplot(ax, df, peak_label, tiled,
                              add_peak=args.add_peak)

        handles = [Patch(facecolor=VCOL[v], edgecolor="black", label=VLAB[v])
                   for v in VARIANT_ORDER]
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.5, -0.008), ncol=len(VARIANT_ORDER),
                   framealpha=0.9, columnspacing=1.0)

        fig.tight_layout(rect=[0, 0.04, 1, 0.93])
        sfx = "_w_stream_peak" if args.add_peak else ""
        stem = f"zaxpy_violins_{scale_name}{sfx}"
        fig.savefig(f"{stem}.png", dpi=180, bbox_inches='tight')
        fig.savefig(f"{stem}.pdf", dpi=180, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {stem}.png / .pdf")

    # ── Summary table ────────────────────────────────────────────────────
    print(f"\n{'Scale':<8} {'Platform':<22} {'Dist':<10} {'Variant':<16} "
          f"{'Best (tpb,cf)':<16} {'Median TB/s':<12} {'% Peak'}")
    print("-" * 96)
    for scale_name, tiled, dfs in scales:
        grid_flat = [("MI300A GPU", dfs["gpu_amd"], "MI300A GPU"),
                     ("GH200 GPU",  dfs["gpu_nv"],  "GH200 GPU")]
        if have_cpu:
            grid_flat += [("MI300A Zen CPU", dfs["cpu_amd"], "MI300A Zen CPU"),
                          ("Grace CPU",      dfs["cpu_nv"],  "Grace CPU")]
        for label, df, pk in grid_flat:
            peak = STREAM_PEAK.get(pk, 1.0)
            for dk in DIST_ORDER:
                sl = select_dist(df, dk, tiled)
                for vk in VARIANT_ORDER:
                    sub = sl[sl["variant"] == vk].copy()
                    if sub.empty: continue
                    sub["bw"] = compute_bandwidth(sub["time_ms"].values,
                                                  sub["nx"].values, vk)
                    meds = sub.groupby(["tpb","coarsen"])["bw"].median()
                    if meds.empty: continue
                    bt, bc = meds.idxmax()
                    m = meds.max()
                    p = 100*m/peak
                    print(f"{scale_name:<8} {label:<22} {dk:<10} {vk:<16} "
                          f"({bt},{bc}){'':>8} {m:.4f}      {p:.1f}%")


if __name__ == "__main__":
    main()