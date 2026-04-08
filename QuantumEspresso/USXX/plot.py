#!/usr/bin/env python3
"""
plot_addusxx_sweep.py
2 violins per panel: Baseline (AoS) vs Transformed (SoA).

GPU CSV: variant,tblock,coarsen,rep,time_ms
CPU CSV: variant,blocksize,ngms,rhoc_size,nthreads,rep,time_ms
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter
import pandas as pd
import numpy as np
import argparse, warnings
from scipy.stats import bootstrap

# ── Variant groups ───────────────────────────────────────────────────────────

GPU_AOS_BASELINE = ["gpu_baseline_aos", "gpu_baseline_aos_u"]
GPU_SOA_ALL = [
    "gpu_baseline_soa", "gpu_baseline_soa_u",
    "gpu_eigts_t_soa", "gpu_eigts_t_soa_u",
    "gpu_shared_bec_soa", "gpu_shared_bec_soa_u",
    "gpu_optimized_soa", "gpu_optimized_soa_u",
]

CPU_AOS_BASELINE = ["baseline_aos"]
CPU_SOA_ALL      = ["baseline_soa", "eigts_t_soa", "sorted_soa"]

GROUP_COLORS = {"baseline": "#e67e22", "transformed": "#2980b9"}

MIN_TIME_MS = 0.1

# ── Helpers ──────────────────────────────────────────────────────────────────

def remove_outliers(vals, k=3.0):
    if len(vals) < 4:
        return vals
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    c = vals[(vals >= lo) & (vals <= hi)]
    return c if len(c) > 2 else vals

def bootstrap_ci(arr, confidence_level=0.95, n_resamples=10000):
    """Bootstrap CI of the median using scipy BCa method, fallback to percentile."""
    if len(arr) < 3:
        med = np.median(arr)
        return med, med
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            res = bootstrap((arr,), np.median,
                            confidence_level=confidence_level,
                            n_resamples=n_resamples, method='BCa')
            lo, hi = res.confidence_interval.low, res.confidence_interval.high
            if np.isnan(lo) or np.isnan(hi):
                raise ValueError("BCa returned nan")
            return lo, hi
        except Exception:
            try:
                res = bootstrap((arr,), np.median,
                                confidence_level=confidence_level,
                                n_resamples=n_resamples, method='percentile')
                lo, hi = res.confidence_interval.low, res.confidence_interval.high
                if np.isnan(lo) or np.isnan(hi):
                    raise ValueError("percentile returned nan")
                return lo, hi
            except Exception:
                med = np.median(arr)
                return med, med

def filter_valid(df):
    return df[df["time_ms"] >= MIN_TIME_MS].copy()


def best_gpu_group(df, variant_list):
    sub = df[df["variant"].isin(variant_list)]
    if sub.empty:
        return np.array([]), ""
    medians = sub.groupby(["variant", "tblock", "coarsen"])["time_ms"].median()
    best_key = medians.idxmin()
    vals = sub[
        (sub["variant"] == best_key[0]) &
        (sub["tblock"] == best_key[1]) &
        (sub["coarsen"] == best_key[2])
    ]["time_ms"].values
    return vals, ""


def best_cpu_group(df, variant_list):
    sub = df[df["variant"].isin(variant_list)]
    if sub.empty:
        return np.array([]), ""
    medians = sub.groupby(["variant", "blocksize", "nthreads"])["time_ms"].median()
    best_key = medians.idxmin()
    vals = sub[
        (sub["variant"] == best_key[0]) &
        (sub["blocksize"] == best_key[1]) &
        (sub["nthreads"] == best_key[2])
    ]["time_ms"].values
    return vals, ""


# ── Subplot ──────────────────────────────────────────────────────────────────

def paint_subplot(ax, df, platform_label, is_gpu):
    ax.set_box_aspect(2.8 / 4.8)
    df = filter_valid(df)
    get_best = best_gpu_group if is_gpu else best_cpu_group
    aos_variants = GPU_AOS_BASELINE if is_gpu else CPU_AOS_BASELINE
    soa_variants = GPU_SOA_ALL if is_gpu else CPU_SOA_ALL

    entries = [
        ("baseline",    aos_variants, "Baseline"),
        ("transformed", soa_variants, "Transformed"),
    ]

    positions, data_all, col_all, tick_labels = [], [], [], []
    for i, (gkey, variants, xlabel) in enumerate(entries):
        vals, _ = get_best(df, variants)
        vals = remove_outliers(vals)
        if len(vals) > 0:
            positions.append(i)
            data_all.append(vals)
            col_all.append(GROUP_COLORS[gkey])
            tick_labels.append(xlabel)

    if not data_all:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
        return

    parts = ax.violinplot(data_all, positions=positions,
                          showmeans=True, showmedians=True,
                          showextrema=False, widths=0.6)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(col_all[i])
        body.set_edgecolor("black")
        body.set_alpha(0.75)
        body.set_zorder(2)
    parts["cmeans"].set_color("black")
    parts["cmeans"].set_zorder(3)
    parts["cmedians"].set_color("white")
    parts["cmedians"].set_zorder(3)

    # Bootstrap 95% CI of median
    ci_color = "#FF69B4"  # light pink
    global_max = max(np.max(d) for d in data_all)
    nice_top = np.ceil(global_max * 1.25 * 4) / 4
    for i, arr in enumerate(data_all):
        ci_lo, ci_hi = bootstrap_ci(arr)
        ci_mid = (ci_lo + ci_hi) / 2
        min_height = 0.01 * nice_top if nice_top > 0 else 0.001
        if (ci_hi - ci_lo) < min_height:
            ci_lo = ci_mid - min_height / 2
            ci_hi = ci_mid + min_height / 2
        ax.vlines(positions[i], ci_lo, ci_hi,
                  color="black", lw=3.5, zorder=10)
        ax.vlines(positions[i], ci_lo, ci_hi,
                  color=ci_color, lw=1.8, zorder=11)
        cap_w = 0.08
        for y in (ci_lo, ci_hi):
            ax.hlines(y, positions[i] - cap_w, positions[i] + cap_w,
                      color="black", lw=3, zorder=10)
            ax.hlines(y, positions[i] - cap_w, positions[i] + cap_w,
                      color=ci_color, lw=1.5, zorder=11)

    # Median annotations — italic, above violin
    medians = []
    for i, pos in enumerate(positions):
        med = np.median(data_all[i])
        medians.append(med)
        ymax = np.max(data_all[i])
        ax.text(pos, ymax * 1.04, f"med. {med:.3f} ms",
                ha="center", va="bottom", fontsize=9,
                fontstyle="italic", fontweight="normal",
                color=col_all[i])

    # Speedup below transformed violin
    if len(medians) == 2 and medians[1] > 0:
        speedup = medians[0] / medians[1]
        ymin = np.min(data_all[1])
        ax.text(positions[1] + 0.025, ymin * 0.95, f"{speedup:.2f}×",
                ha="center", va="top", fontsize=10,
                fontweight="bold", color=GROUP_COLORS["transformed"])

    # ylim + ticks + grid
    ax.set_ylim(bottom=0, top=nice_top)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=5))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis='y', which='minor', length=3)
    ax.grid(axis="y", which="major", alpha=0.25, zorder=0)
    ax.grid(axis="y", which="minor", alpha=0.12, ls=':', zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", which="minor", bottom=False)

    ax.set_xticks(positions)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.set_title(platform_label, fontsize=11)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu-amd-csv", default="results/beverin/addusxx_gpu_sweep.csv")
    ap.add_argument("--gpu-nv-csv", default="results/daint/addusxx_gpu_sweep.csv")
    ap.add_argument("--cpu-amd-csv", default="results/beverin/addusxx_cpu_sweep.csv")
    ap.add_argument("--cpu-nv-csv", default="results/daint/addusxx_cpu_sweep.csv")
    args = ap.parse_args()

    def load(p):
        if p is None:
            return None
        try:
            return pd.read_csv(p)
        except FileNotFoundError:
            print(f"WARN: {p} not found")
            return None

    gpu_amd = load(args.gpu_amd_csv)
    gpu_nv  = load(args.gpu_nv_csv)
    cpu_amd = load(args.cpu_amd_csv)
    cpu_nv  = load(args.cpu_nv_csv)

    have_gpu = gpu_amd is not None or gpu_nv is not None
    have_cpu = cpu_amd is not None or cpu_nv is not None

    if not have_gpu and not have_cpu:
        print("ERROR: No CSV files provided.")
        return

    GRID = []
    if have_cpu:
        GRID.append([
            ("MI300A Zen CPU", cpu_amd, False),
            ("GH200 Grace CPU", cpu_nv, False),
        ])
    if have_gpu:
        GRID.append([
            ("MI300A GPU", gpu_amd, True),
            ("GH200 Hopper GPU", gpu_nv, True),
        ])

    nrows = len(GRID)
    ncols = 2

    # No STREAM subtitle → try -0.4
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.6 * ncols, 2.8 * nrows),
                             squeeze=False)

    fig.suptitle("`addusxx_g` Kernel: Original vs Transformed Layout",
                 fontsize=15, y=0.89 if nrows > 1 else 0.98)

    for ri, row in enumerate(GRID):
        for ci, (label, df, is_gpu) in enumerate(row):
            ax = axes[ri, ci]
            if ci == 0:
                ax.set_ylabel("Kernel time [ms]", fontsize=11)
            if df is not None:
                paint_subplot(ax, df, label, is_gpu)
            else:
                ax.set_title(label, fontsize=11)
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=11, color="gray")

    handles = [
        Patch(facecolor=GROUP_COLORS["baseline"], edgecolor="black", label="Baseline (AoS)"),
        Patch(facecolor=GROUP_COLORS["transformed"], edgecolor="black", label="Transformed (SoA+Permute+Shuffle)"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=2,
               framealpha=0.9, columnspacing=1.5, fontsize=10)

    fig.tight_layout(rect=[0, 0.06, 1, 0.89 if nrows > 1 else 0.92], h_pad=4.0)
    fig.savefig("addusxx_sweep.png", dpi=200, bbox_inches="tight")
    fig.savefig("addusxx_sweep.pdf", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved addusxx_sweep.png / .pdf")

    # Summary
    print(f"\n{'Platform':<28} {'Group':<14} {'Median ms':<12}")
    print("-" * 54)
    for row in GRID:
        for label, df, is_gpu in row:
            if df is None:
                continue
            df_f = filter_valid(df)
            get_best = best_gpu_group if is_gpu else best_cpu_group
            aos_v = GPU_AOS_BASELINE if is_gpu else CPU_AOS_BASELINE
            soa_v = GPU_SOA_ALL if is_gpu else CPU_SOA_ALL
            for gname, variants in [("Original", aos_v), ("Transformed", soa_v)]:
                vals, _ = get_best(df_f, variants)
                if len(vals) == 0:
                    continue
                med = np.median(vals)
                print(f"{label:<28} {gname:<14} {med:.4f}")


if __name__ == "__main__":
    main()