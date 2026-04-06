#!/usr/bin/env python3
"""
plot_zaxpy_sweep.py
Violin plots for indirect zaxpy benchmark.

One PNG per QE problem size (nat05, nat10, nat20).
X-axis groups: Uniform, Normal, QE  (from 3 separate CSV files)
Colours = 4 kernel variants:
    aos_scatter (orange)  — original interleaved layout, random write order
    aos_sorted  (blue)    — interleaved layout, sorted write order
    soa_scatter (red)     — split re/im layout, random write order
    soa_sorted  (green)   — split re/im layout, sorted write order

For each (variant, distribution) we pick the (tpb, coarsen) config with the
highest median runtime and show the violin for that config's 100 repetitions.

Usage:
    python plot_zaxpy_sweep.py [--qe-csv FILE] [--uniform-csv FILE] [--normal-csv FILE]
                               [--ny-mult INT]

python plot.py --qe-csv zaxpy_sweep_qe.csv \
                           --uniform-csv zaxpy_sweep_uniform.csv \
                           --normal-csv zaxpy_sweep_normal.csv \
                           --ny-mult 2 --sample 0
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import pandas as pd, numpy as np, argparse, sys

# ── Colour / label scheme ────────────────────────────────────────────────────
VCOL = {
    "aos_scatter": "#e67e22",   # orange
    "aos_sorted":  "#2980b9",   # blue
    "soa_scatter": "#c0392b",   # red
    "soa_sorted":  "#27ae60",   # green
}
VLAB = {
    "aos_scatter": "AoS  Scatter (original)",
    "aos_sorted":  "AoS  Sorted",
    "soa_scatter": "SoA  Scatter",
    "soa_sorted":  "SoA  Sorted",
}
VARIANT_ORDER = ["aos_scatter", "aos_sorted", "soa_scatter", "soa_sorted"]

DIST_LABEL = {
    "uniform": "Uniform",
    "normal":  "Normal",
    "qe":      "QE",
}
DIST_ORDER = ["uniform", "normal", "qe"]

QE_CASES = {
    "qe_nat05": {"ny": 64001,  "nx": 27609,  "title": "BaTiO₃  5-atom"},
    "qe_nat10": {"ny": 120001, "nx": 55191,  "title": "BaTiO₃ 10-atom"},
    "qe_nat20": {"ny": 225001, "nx": 110273, "title": "BaTiO₃ 20-atom"},
}

SUBPLOT_W = 5.5
SUBPLOT_H = 4.0

# ── Helpers ──────────────────────────────────────────────────────────────────

def remove_outliers(vals, k=3.0):
    if len(vals) < 4:
        return vals
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    clean = vals[(vals >= lo) & (vals <= hi)]
    return clean if len(clean) > 2 else vals


def best_config_data(df, variant):
    """
    For a given variant, find the (tpb, coarsen) with the lowest median
    runtime and return the time_ms array for that config.
    """
    sub = df[df["variant"] == variant].copy()
    if sub.empty:
        return np.array([])
    medians = sub.groupby(["tpb", "coarsen"])["time_ms"].median()
    if medians.empty:
        return np.array([])
    best_tpb, best_cf = medians.idxmin()
    vals = sub[(sub["tpb"] == best_tpb) & (sub["coarsen"] == best_cf)]["time_ms"].values
    return vals


def filter_synthetic(df, ny_mult, sample=0):
    """
    From a uniform / normal CSV, keep only rows matching a specific
    ny_multiple and sample index (e.g. 'uniform_mult2_s0').
    Falls back to first available sample if requested one is missing.
    """
    # experiment column looks like  uniform_mult2_s0  or  normal_mult4_s3
    candidates = df[df["experiment"].str.contains(f"mult{ny_mult}_")]
    if candidates.empty:
        return candidates
    # pick one sample
    tag = f"mult{ny_mult}_s{sample}"
    sub = candidates[candidates["experiment"].str.contains(tag)]
    if sub.empty:
        # fall back to first available sample
        first_exp = candidates["experiment"].iloc[0]
        sub = candidates[candidates["experiment"] == first_exp]
    return sub


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qe-csv",      default="zaxpy_sweep_qe.csv")
    ap.add_argument("--uniform-csv", default="zaxpy_sweep_uniform.csv")
    ap.add_argument("--normal-csv",  default="zaxpy_sweep_normal.csv")
    ap.add_argument("--ny-mult",     type=int, default=2,
                    help="Which ny/nx multiple to use from synthetic CSVs (default 2)")
    ap.add_argument("--sample",      type=int, default=0,
                    help="Which sample index to use from synthetic CSVs (default 0)")
    args = ap.parse_args()

    try:
        df_qe  = pd.read_csv(args.qe_csv)
        df_uni = pd.read_csv(args.uniform_csv)
        df_nor = pd.read_csv(args.normal_csv)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr); return

    plt.rcParams.update({
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11,
    })

    ndists = len(DIST_ORDER)
    nvars  = len(VARIANT_ORDER)
    group_spacing = nvars + 1.5          # positions per distribution group

    for qe_name, qe_info in QE_CASES.items():

        fig, ax = plt.subplots(1, 1,
                               figsize=(SUBPLOT_W * 1.6, SUBPLOT_H))
        fig.suptitle(f"Indirect ZAXPY — {qe_info['title']}   "
                     f"(ny={qe_info['ny']:,}  nx={qe_info['nx']:,})",
                     fontsize=14, y=0.98)

        positions, data_all, col_all = [], [], []
        medians_annot = []          # (pos, median, variant_key)
        xticks, xlabels = [], []

        for di, dist_key in enumerate(DIST_ORDER):
            # ── Select the right dataframe slice ──────────────────────────
            if dist_key == "qe":
                df_slice = df_qe[df_qe["experiment"] == qe_name]
            elif dist_key == "uniform":
                df_slice = filter_synthetic(df_uni, args.ny_mult, args.sample)
            else:  # normal
                df_slice = filter_synthetic(df_nor, args.ny_mult, args.sample)

            for vi, vk in enumerate(VARIANT_ORDER):
                vals = best_config_data(df_slice, vk)
                vals = remove_outliers(vals)
                pos = di * group_spacing + vi
                if len(vals) > 0:
                    data_all.append(vals)
                    positions.append(pos)
                    col_all.append(VCOL[vk])
                    medians_annot.append((pos, np.median(vals), vk))

            xticks.append(di * group_spacing + (nvars - 1) / 2.0)
            xlabels.append(DIST_LABEL[dist_key])

        # ── Y-axis limits (6 nice ticks) ─────────────────────────────────
        all_flat = np.concatenate(data_all) if data_all else np.array([0.0])
        max_val = float(np.max(all_flat))
        locator = MaxNLocator(nbins=5, min_n_ticks=5)
        ticks = locator.tick_values(0.0, max_val * 1.14)
        ticks = ticks[ticks >= 0]
        if len(ticks) > 6:
            ticks = ticks[:6]
        top_lim = ticks[-1] * 1.06 if len(ticks) else max_val * 1.15

        # ── Violins ──────────────────────────────────────────────────────
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

        # ── Median annotations (below each violin) ──────────────────────
        ylo, yhi = 0, top_lim
        y_range = yhi - ylo
        for pos, med, vk in medians_annot:
            ax.text(pos, med - 0.04 * y_range,
                    f'{med:.3f}',
                    ha='center', va='top',
                    fontsize=9, color=VCOL[vk], fontweight='bold')

        # ── Decoration ───────────────────────────────────────────────────
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)
        ax.set_yticks(ticks)
        ax.set_ylim(bottom=0, top=top_lim)
        ax.set_ylabel("Kernel time [ms]")
        ax.grid(axis="y", alpha=0.3)

        # ── Legend ────────────────────────────────────────────────────────
        handles = [Patch(facecolor=VCOL[vk], edgecolor="black", label=VLAB[vk])
                   for vk in VARIANT_ORDER]
        ax.legend(handles=handles, loc='upper right', framealpha=0.9,
                  fontsize=10)

        fig.tight_layout(rect=[0, 0, 1, 0.96])
        out = f"zaxpy_violins_{qe_name}.png"
        fig.savefig(out, dpi=180, bbox_inches='tight')
        fig.savefig(out.replace(".png", ".pdf"), dpi=180, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {out}")

    # ── Summary table ────────────────────────────────────────────────────────
    print(f"\n{'QE Case':<12} {'Dist':<10} {'Variant':<16} {'Best (tpb,cf)':<16} {'Median ms':<12}")
    print("-" * 66)
    for qe_name in QE_CASES:
        for dist_key in DIST_ORDER:
            if dist_key == "qe":
                df_s = df_qe[df_qe["experiment"] == qe_name]
            elif dist_key == "uniform":
                df_s = filter_synthetic(df_uni, args.ny_mult, args.sample)
            else:
                df_s = filter_synthetic(df_nor, args.ny_mult, args.sample)
            for vk in VARIANT_ORDER:
                sub = df_s[df_s["variant"] == vk]
                if sub.empty:
                    continue
                meds = sub.groupby(["tpb", "coarsen"])["time_ms"].median()
                if meds.empty:
                    continue
                best_tpb, best_cf = meds.idxmin()
                med = meds.min()
                print(f"{qe_name:<12} {dist_key:<10} {vk:<16} ({best_tpb},{best_cf}){'':>8} {med:.6f}")


if __name__ == "__main__":
    main()