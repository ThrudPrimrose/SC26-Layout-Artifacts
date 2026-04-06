#!/usr/bin/env python3
"""
plot_addusxx_sweep.py
Violin plots for addusxx kernel sweep: 2×2 grid (CPU/GPU × AMD/NV).

GPU CSV format: variant,tblock,coarsen,time_ms,correct
CPU CSV format: variant,ngms,rhoc_size,nthreads,rep,time_ms

Usage:
    python plot_addusxx_sweep.py \
        --gpu-amd-csv amd/addusxx_gpu_sweep.csv \
        --gpu-nv-csv  nv/addusxx_gpu_sweep.csv  \
        --cpu-amd-csv amd/addusxx_cpu_sweep.csv  \
        --cpu-nv-csv  nv/addusxx_cpu_sweep.csv
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import numpy as np
import argparse

# ── Variant definitions ──────────────────────────────────────────────────────

# Unified variant keys (GPU prefix "gpu_" stripped, CPU "sorted" = GPU "optimized")
VARIANT_ORDER = [
    "baseline_aos", "eigts_t_aos", "shared_bec_aos", "optimized_aos",
    "baseline_soa", "eigts_t_soa", "shared_bec_soa", "optimized_soa",
]

VARIANT_ORDER_CPU = [
    "baseline_aos", "eigts_t_aos", "sorted_aos",
    "baseline_soa", "eigts_t_soa", "sorted_soa",
]

VARIANT_LABELS = {
    "baseline_aos":   "AoS Baseline",
    "eigts_t_aos":    "AoS Eigts-T",
    "shared_bec_aos": "AoS Shared-Bec",
    "optimized_aos":  "AoS Sorted",
    "sorted_aos":     "AoS Sorted",
    "baseline_soa":   "SoA Baseline",
    "eigts_t_soa":    "SoA Eigts-T",
    "shared_bec_soa": "SoA Shared-Bec",
    "optimized_soa":  "SoA Sorted",
    "sorted_soa":     "SoA Sorted",
}

VARIANT_COLORS = {
    "baseline_aos":   "#f39c12",  # light orange
    "eigts_t_aos":    "#e67e22",  # orange
    "shared_bec_aos": "#d35400",  # dark orange
    "optimized_aos":  "#a04000",  # darker orange
    "sorted_aos":     "#a04000",
    "baseline_soa":   "#5dade2",  # light blue
    "eigts_t_soa":    "#2e86c1",  # blue
    "shared_bec_soa": "#1a5276",  # dark blue
    "optimized_soa":  "#0b3d5c",  # darker blue
    "sorted_soa":     "#0b3d5c",
}

# AoS group and SoA group for x-axis
GROUPS = {
    "gpu": {
        "AoS": ["baseline_aos", "eigts_t_aos", "shared_bec_aos", "optimized_aos"],
        "SoA": ["baseline_soa", "eigts_t_soa", "shared_bec_soa", "optimized_soa"],
    },
    "cpu": {
        "AoS": ["baseline_aos", "eigts_t_aos", "sorted_aos"],
        "SoA": ["baseline_soa", "eigts_t_soa", "sorted_soa"],
    },
}

# ── Helpers ──────────────────────────────────────────────────────────────────

def remove_outliers(vals, k=3.0):
    if len(vals) < 4:
        return vals
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    c = vals[(vals >= lo) & (vals <= hi)]
    return c if len(c) > 2 else vals


def best_gpu_times(df, variant):
    """Best (tblock, coarsen) by lowest median time for a GPU variant."""
    gpu_name = f"gpu_{variant}" if not variant.startswith("gpu_") else variant
    sub = df[df["variant"] == gpu_name]
    if sub.empty:
        return np.array([])
    medians = sub.groupby(["tblock", "coarsen"])["time_ms"].median()
    if medians.empty:
        return np.array([])
    best = medians.idxmin()
    vals = sub[(sub["tblock"] == best[0]) & (sub["coarsen"] == best[1])]["time_ms"].values
    return vals


def best_cpu_times(df, variant):
    """Best nthreads by lowest median time for a CPU variant."""
    sub = df[df["variant"] == variant]
    if sub.empty:
        return np.array([])
    medians = sub.groupby("nthreads")["time_ms"].median()
    if medians.empty:
        return np.array([])
    best_nt = medians.idxmin()
    vals = sub[sub["nthreads"] == best_nt]["time_ms"].values
    return vals


# ── Subplot painter ──────────────────────────────────────────────────────────

def paint_subplot(ax, df, platform_label, is_gpu):
    groups = GROUPS["gpu"] if is_gpu else GROUPS["cpu"]
    get_times = best_gpu_times if is_gpu else best_cpu_times

    positions, data_all, col_all = [], [], []
    xticks, xlabels = [], []
    group_offset = 0

    for gname, variants in groups.items():
        nvars = len(variants)
        for vi, vk in enumerate(variants):
            vals = get_times(df, vk)
            vals = remove_outliers(vals)
            pos = group_offset + vi
            if len(vals) > 0:
                data_all.append(vals)
                positions.append(pos)
                col_all.append(VARIANT_COLORS[vk])

        xticks.append(group_offset + (nvars - 1) / 2.0)
        xlabels.append(gname)
        group_offset += nvars + 1.5

    if not data_all:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
        return

    # Violins
    parts = ax.violinplot(data_all, positions=positions,
                          showmeans=True, showmedians=True,
                          showextrema=False, widths=0.85)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(col_all[i])
        body.set_edgecolor("black")
        body.set_alpha(0.75)
    parts["cmeans"].set_color("black")
    parts["cmedians"].set_color("white")

    # Median annotation
    for i, pos in enumerate(positions):
        med = np.median(data_all[i])
        ax.text(pos, med * 0.92, f"{med:.3f}",
                ha="center", va="top", fontsize=8, fontweight="bold",
                color=col_all[i])

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title(platform_label)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu-amd-csv", default=None)
    ap.add_argument("--gpu-nv-csv", default=None)
    ap.add_argument("--cpu-amd-csv", default=None)
    ap.add_argument("--cpu-nv-csv", default=None)
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
    gpu_nv = load(args.gpu_nv_csv)
    cpu_amd = load(args.cpu_amd_csv)
    cpu_nv = load(args.cpu_nv_csv)

    have_gpu = gpu_amd is not None or gpu_nv is not None
    have_cpu = cpu_amd is not None or cpu_nv is not None

    if not have_gpu and not have_cpu:
        print("ERROR: No CSV files provided.")
        return

    plt.rcParams.update({
        "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13,
        "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10,
    })

    # ── Build 2×2 grid ───────────────────────────────────────────────────
    # Row 0: CPU (AMD, NV)
    # Row 1: GPU (AMD, NV)
    GRID = []
    if have_cpu:
        GRID.append([
            ("MI300A Zen4 CPU", cpu_amd, False),
            ("Grace Neoverse V2 CPU", cpu_nv, False),
        ])
    if have_gpu:
        GRID.append([
            ("MI300A GPU", gpu_amd, True),
            ("GH200 GPU", gpu_nv, True),
        ])

    nrows = len(GRID)
    ncols = 2
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(6.0 * ncols, 3.8 * nrows),
                             squeeze=False)

    fig.suptitle("addusxx_g Kernel — Layout × Optimization Sweep",
                 fontsize=15, y=0.98)

    for ri, row in enumerate(GRID):
        for ci, (label, df, is_gpu) in enumerate(row):
            ax = axes[ri, ci]
            if ci == 0:
                ax.set_ylabel("Kernel time [ms]")
            if df is not None:
                paint_subplot(ax, df, label, is_gpu)
            else:
                ax.set_title(label)
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=14, color="gray")

    # ── Legend ────────────────────────────────────────────────────────────
    # Collect all variants that appear
    all_variants = list(VARIANT_ORDER)
    if have_cpu:
        for v in VARIANT_ORDER_CPU:
            if v not in all_variants:
                all_variants.append(v)
    # Deduplicate sorted_aos/optimized_aos (same color/label)
    seen_labels = set()
    handles = []
    for v in all_variants:
        lab = VARIANT_LABELS[v]
        if lab in seen_labels:
            continue
        seen_labels.add(lab)
        handles.append(Patch(facecolor=VARIANT_COLORS[v], edgecolor="black",
                             label=lab))

    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.01), ncol=4,
               framealpha=0.9, columnspacing=1.0)

    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig("addusxx_sweep.png", dpi=180, bbox_inches="tight")
    fig.savefig("addusxx_sweep.pdf", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("Saved addusxx_sweep.png / .pdf")

    # ── Summary table ────────────────────────────────────────────────────
    print(f"\n{'Platform':<28} {'Variant':<22} {'Best Config':<18} {'Median ms':<12}")
    print("-" * 80)

    for label, df, is_gpu in [(l, d, g) for row in GRID for l, d, g in row]:
        if df is None:
            continue
        groups = GROUPS["gpu"] if is_gpu else GROUPS["cpu"]
        get_fn = best_gpu_times if is_gpu else best_cpu_times
        for gname, variants in groups.items():
            for vk in variants:
                vals = get_fn(df, vk)
                if len(vals) == 0:
                    continue
                med = np.median(vals)

                # Find best config
                if is_gpu:
                    gpu_name = f"gpu_{vk}"
                    sub = df[df["variant"] == gpu_name]
                    meds = sub.groupby(["tblock", "coarsen"])["time_ms"].median()
                    bt, bc = meds.idxmin()
                    cfg = f"tblock={bt},coarsen={bc}"
                else:
                    sub = df[df["variant"] == vk]
                    meds = sub.groupby("nthreads")["time_ms"].median()
                    bnt = meds.idxmin()
                    cfg = f"nthreads={bnt}"

                print(f"{label:<28} {VARIANT_LABELS[vk]:<22} {cfg:<18} {med:.4f}")


if __name__ == "__main__":
    main()