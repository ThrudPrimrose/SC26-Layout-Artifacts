#!/usr/bin/env python3
"""
Plot elementwise_16x16 benchmark results.

Expected CSV columns:
  kernel_name,layout,use_smem,M,N,run,time_ms,bw_gbps,checksum

Produces:
  1. violins_by_layout.pdf   -- 4 panels (one per layout), 12 violins each
  2. violins_smem_compare.pdf -- direct vs smem paired, one panel per base kernel
  3. heatmap_median.pdf       -- heatmap: kernel x layout, color = median time
  4. summary_bw.pdf           -- bandwidth violin, all configs
"""

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ---------- style ----------
sns.set_theme(style="whitegrid", font_scale=0.9)
DIRECT_COLOR = "#534AB7"  # purple
SMEM_COLOR   = "#1D9E75"  # teal
PALETTE = {0: DIRECT_COLOR, 1: SMEM_COLOR}

# nice short labels for kernels (strip layout prefix)
def short_kernel(name: str) -> str:
    """Ar_Br_direct -> direct, Ac_Br_smem_2x2 -> smem_2x2"""
    parts = name.split("_", 2)  # ['Ar', 'Br', 'direct'] or ['Ar', 'Br', 'smem_flat_row']
    return parts[2] if len(parts) >= 3 else name

# base kernel (strip smem_ prefix to pair direct<->smem)
def base_kernel(short: str) -> str:
    """direct -> direct, smem -> direct, smem_2x2 -> direct_2x2, flat_row -> flat_row, smem_flat_row -> flat_row"""
    if short == "smem":
        return "direct"
    if short.startswith("smem_flat"):
        return short.replace("smem_flat", "flat")
    if short.startswith("smem_"):
        return "direct_" + short[5:]
    return short

BASE_ORDER = ["direct", "flat_row", "flat_col", "direct_2x2", "direct_4x1", "direct_1x4"]
LAYOUT_ORDER = ["Ar_Br", "Ar_Bc", "Ac_Br", "Ac_Bc"]

# ---------- load ----------
def load(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["short"] = df["kernel_name"].apply(short_kernel)
    df["base"] = df["short"].apply(base_kernel)
    df["mem"] = df["use_smem"].map({0: "direct", 1: "smem"})
    return df


# ---------- plot 1: violins by layout ----------
def plot_violins_by_layout(df: pd.DataFrame, out: str):
    """One subplot per layout, 12 violins (kernels) each, colored direct/smem."""
    layouts = [l for l in LAYOUT_ORDER if l in df["layout"].unique()]
    n = len(layouts)
    fig, axes = plt.subplots(n, 1, figsize=(14, 3.5 * n), sharex=False)
    if n == 1:
        axes = [axes]

    # consistent kernel order across panels
    kernel_order = []
    for b in BASE_ORDER:
        kernel_order.append(b)
        # find smem counterpart
        smem_name = "smem" if b == "direct" else b.replace("direct_", "smem_").replace("flat_", "smem_flat_")
        kernel_order.append(smem_name)

    for ax, layout in zip(axes, layouts):
        sub = df[df["layout"] == layout].copy()
        # keep only kernels present
        present = [k for k in kernel_order if k in sub["short"].unique()]
        sub["short"] = pd.Categorical(sub["short"], categories=present, ordered=True)

        sns.violinplot(
            data=sub, x="short", y="time_ms", hue="use_smem",
            palette=PALETTE, inner="box", linewidth=0.6,
            density_norm="width", cut=0, ax=ax, dodge=False,
            legend=False,
        )
        ax.set_title(f"Layout: {layout}", fontweight="bold", fontsize=12)
        ax.set_xlabel("")
        ax.set_ylabel("Time (ms)")
        ax.tick_params(axis="x", rotation=35)

    # shared legend
    handles = [
        mpatches.Patch(color=DIRECT_COLOR, label="direct (no smem)"),
        mpatches.Patch(color=SMEM_COLOR, label="shared memory"),
    ]
    fig.legend(handles=handles, loc="upper right", frameon=True, fontsize=10)
    fig.suptitle("Kernel runtime distribution by layout", fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=150)
    print(f"  Saved {out}")
    plt.close(fig)


# ---------- plot 2: direct vs smem paired ----------
def plot_smem_compare(df: pd.DataFrame, out: str):
    """One subplot per base kernel (6), each shows 4 layouts x direct/smem violins."""
    bases = [b for b in BASE_ORDER if b in df["base"].unique()]
    n = len(bases)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for idx, (ax, bk) in enumerate(zip(axes, bases)):
        sub = df[df["base"] == bk].copy()
        sub["layout"] = pd.Categorical(sub["layout"], categories=LAYOUT_ORDER, ordered=True)

        sns.violinplot(
            data=sub, x="layout", y="time_ms", hue="mem",
            palette={"direct": DIRECT_COLOR, "smem": SMEM_COLOR},
            inner="box", linewidth=0.6, density_norm="width",
            cut=0, ax=ax, dodge=True, legend=(idx == 0),
        )
        ax.set_title(bk, fontweight="bold", fontsize=11)
        ax.set_xlabel("")
        ax.set_ylabel("Time (ms)" if idx % 3 == 0 else "")

    # hide unused axes
    for j in range(len(bases), len(axes)):
        axes[j].set_visible(False)

    if axes[0].get_legend():
        axes[0].get_legend().set_title("")
    fig.suptitle("Direct vs shared memory by base kernel", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=150)
    print(f"  Saved {out}")
    plt.close(fig)


# ---------- plot 3: heatmap of median times ----------
def plot_heatmap(df: pd.DataFrame, out: str):
    """Heatmap: rows = kernel variants, cols = layouts, color = median time_ms."""
    kernel_order = []
    for b in BASE_ORDER:
        kernel_order.append(b)
        smem_name = "smem" if b == "direct" else b.replace("direct_", "smem_").replace("flat_", "smem_flat_")
        kernel_order.append(smem_name)

    pivot = df.groupby(["short", "layout"])["time_ms"].median().unstack("layout")
    # reorder
    present_k = [k for k in kernel_order if k in pivot.index]
    present_l = [l for l in LAYOUT_ORDER if l in pivot.columns]
    pivot = pivot.loc[present_k, present_l]

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(
        pivot, annot=True, fmt=".3f", cmap="YlOrRd",
        linewidths=0.5, linecolor="white", ax=ax,
        cbar_kws={"label": "Median time (ms)"},
    )
    ax.set_title("Median kernel time (ms)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Kernel variant")
    ax.set_xlabel("Layout")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=150)
    print(f"  Saved {out}")
    plt.close(fig)


# ---------- plot 4: bandwidth violins ----------
def plot_bandwidth(df: pd.DataFrame, out: str):
    """Violin plot of effective bandwidth, all configs, grouped by kernel, colored by smem."""
    kernel_order = []
    for b in BASE_ORDER:
        kernel_order.append(b)
        smem_name = "smem" if b == "direct" else b.replace("direct_", "smem_").replace("flat_", "smem_flat_")
        kernel_order.append(smem_name)
    present = [k for k in kernel_order if k in df["short"].unique()]
    df_copy = df.copy()
    df_copy["short"] = pd.Categorical(df_copy["short"], categories=present, ordered=True)

    fig, ax = plt.subplots(figsize=(14, 5))
    sns.violinplot(
        data=df_copy, x="short", y="bw_gbps", hue="use_smem",
        palette=PALETTE, inner="box", linewidth=0.6,
        density_norm="width", cut=0, ax=ax, dodge=False,
        legend=False,
    )
    handles = [
        mpatches.Patch(color=DIRECT_COLOR, label="direct"),
        mpatches.Patch(color=SMEM_COLOR, label="smem"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=True)
    ax.set_xlabel("")
    ax.set_ylabel("Effective bandwidth (GB/s)")
    ax.set_title("Bandwidth distribution across all layouts", fontsize=13, fontweight="bold")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=150)
    print(f"  Saved {out}")
    plt.close(fig)


# ---------- plot 5: speedup heatmap (smem time / direct time) ----------
def plot_speedup_heatmap(df: pd.DataFrame, out: str):
    """Heatmap of smem speedup over direct: median(direct) / median(smem).
       >1 means smem is faster."""
    medians = df.groupby(["base", "layout", "mem"])["time_ms"].median().reset_index()

    direct = medians[medians["mem"] == "direct"].set_index(["base", "layout"])["time_ms"]
    smem   = medians[medians["mem"] == "smem"].set_index(["base", "layout"])["time_ms"]

    speedup = (direct / smem).unstack("layout")
    present_b = [b for b in BASE_ORDER if b in speedup.index]
    present_l = [l for l in LAYOUT_ORDER if l in speedup.columns]
    speedup = speedup.loc[present_b, present_l]

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.heatmap(
        speedup, annot=True, fmt=".2f",
        cmap="RdYlGn", center=1.0,
        linewidths=0.5, linecolor="white", ax=ax,
        cbar_kws={"label": "Speedup (direct / smem)"},
        vmin=0.7, vmax=1.3,
    )
    ax.set_title("Shared memory speedup (>1 = smem faster)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Base kernel")
    ax.set_xlabel("Layout")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=150)
    print(f"  Saved {out}")
    plt.close(fig)


# ---------- main ----------
def main():
    parser = argparse.ArgumentParser(description="Plot elementwise_16x16 benchmark results")
    parser.add_argument("csv", help="Path to results CSV")
    parser.add_argument("--prefix", default="bench", help="Output filename prefix")
    args = parser.parse_args()

    df = load(args.csv)
    print(f"Loaded {len(df)} rows: {df['kernel_name'].nunique()} configs, "
          f"{df['run'].nunique()} runs each")
    print(f"Layouts: {sorted(df['layout'].unique())}")
    print(f"Kernels: {sorted(df['short'].unique())}")
    print()

    plot_violins_by_layout(df, f"{args.prefix}_violins_by_layout.pdf")
    plot_smem_compare(df, f"{args.prefix}_smem_compare.pdf")
    plot_heatmap(df, f"{args.prefix}_heatmap_median.pdf")
    plot_bandwidth(df, f"{args.prefix}_bandwidth.pdf")
    plot_speedup_heatmap(df, f"{args.prefix}_speedup_heatmap.pdf")

    print("\nDone.")


if __name__ == "__main__":
    main()