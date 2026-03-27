#!/usr/bin/env python3
"""
Plot cost model metrics from metrics.csv.
Both loop orders combined per chart, hatched bars for klev_first.
Labels: "V{n}:layout | dist | loop_order"
Saves PDF + PNG.

Usage:
    python plot_metrics.py
    python plot_metrics.py --input metrics.csv -o plots/
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import pandas as pd
except ImportError:
    print("pip install pandas matplotlib"); sys.exit(1)

VARIANT_COLORS = {1: "#4c72b0", 2: "#dd8452", 3: "#55a868", 4: "#c44e52"}
VARIANT_LABELS = {1: "V1:none", 2: "V2:idx", 3: "V3:compute", 4: "V4:both"}
LOOP_SHORT = {"klon_first": "klon1st", "klev_first": "klev1st"}

METRICS = [
    ("mu",    r"$\mu$ : avg new-block count per vector step"),
    ("delta", r"$\Delta$ : avg block distance"),
    ("sigma", r"$\sigma$ : avg element stride / 16"),
]


def _save(fig, out_dir, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"{stem}.{ext}"), dpi=180)
    print(f"  saved {stem}.{{pdf,png}}")


def _make_label(row):
    vl = VARIANT_LABELS[row["variant"]]
    lo = LOOP_SHORT.get(row["loop_order"], row["loop_order"])
    return f"{vl}  |  {row['cell_dist']}  |  {lo}"


def _hatches(series):
    return ["///" if lo == "klev_first" else "" for lo in series]


def plot_metric(df, metric, title, out_dir):
    """One figure per (metric, nlev, target). Both loops combined."""
    for (nlev, target), g in df.groupby(["nlev", "target"]):
        g = g.copy()
        g["label"] = g.apply(_make_label, axis=1)
        g = g.sort_values(metric, ascending=True).reset_index(drop=True)

        vw = g["vector_width"].iloc[0]
        B  = g["block_bytes"].iloc[0]

        n = len(g)
        fig, ax = plt.subplots(figsize=(9, max(4, 0.3 * n)),
                               constrained_layout=True)

        colors = [VARIANT_COLORS[v] for v in g["variant"]]
        hatches = _hatches(g["loop_order"])
        y = np.arange(n)

        bars = ax.barh(y, g[metric].values, color=colors,
                       edgecolor="black", linewidth=0.3)
        for bar, h in zip(bars, hatches):
            bar.set_hatch(h)

        xmax = g[metric].max()
        for i, val in enumerate(g[metric].values):
            ax.text(val + xmax * 0.01, i, f"{val:.2f}",
                    va="center", ha="left", fontsize=6)

        ax.set_yticks(y)
        ax.set_yticklabels(g["label"].values, fontsize=7)
        ax.set_xlabel(title, fontsize=10)
        ax.set_title(f"{title}\nnlev={nlev}  |  {target}  (B={B}, VW={vw})",
                     fontsize=11, fontweight="bold")
        ax.set_xlim(0, xmax * 1.15)
        ax.grid(axis="x", alpha=0.3)
        ax.invert_yaxis()

        from matplotlib.patches import Patch
        handles = [Patch(facecolor=VARIANT_COLORS[v], edgecolor="black",
                         label=VARIANT_LABELS[v]) for v in sorted(VARIANT_COLORS)]
        handles.append(Patch(facecolor="white", edgecolor="black",
                             hatch="", label="klon_first"))
        handles.append(Patch(facecolor="white", edgecolor="black",
                             hatch="///", label="klev_first"))
        ax.legend(handles=handles, loc="lower right", fontsize=6,
                  ncol=3, framealpha=0.9)

        _save(fig, out_dir, f"metric_{metric}_nlev{nlev}_{target}")
        plt.close(fig)


def plot_combined(df, out_dir):
    """One figure per (nlev, target): 3 metric subplots, both loops."""
    for (nlev, target), g in df.groupby(["nlev", "target"]):
        g = g.copy()
        vw = g["vector_width"].iloc[0]
        B  = g["block_bytes"].iloc[0]

        fig, axes = plt.subplots(1, 3, figsize=(22, max(5, 0.3 * len(g))),
                                 constrained_layout=True)
        fig.suptitle(
            f"Cost Model Metrics  |  nlev={nlev}  |  {target}  (B={B}, VW={vw})",
            fontsize=13, fontweight="bold")

        for col, (metric, title) in enumerate(METRICS):
            ax = axes[col]
            gs = g.copy()
            gs["label"] = gs.apply(_make_label, axis=1)
            gs = gs.sort_values(metric, ascending=True).reset_index(drop=True)

            n = len(gs)
            colors = [VARIANT_COLORS[v] for v in gs["variant"]]
            hatches = _hatches(gs["loop_order"])
            y = np.arange(n)

            bars = ax.barh(y, gs[metric].values, color=colors,
                           edgecolor="black", linewidth=0.3)
            for bar, h in zip(bars, hatches):
                bar.set_hatch(h)

            xmax = gs[metric].max()
            for i, val in enumerate(gs[metric].values):
                ax.text(val + xmax * 0.01, i, f"{val:.1f}",
                        va="center", ha="left", fontsize=5)

            ax.set_yticks(y)
            ax.set_yticklabels(gs["label"].values, fontsize=6)
            ax.set_xlabel(metric, fontsize=9)
            ax.set_title(title, fontsize=9)
            ax.set_xlim(0, xmax * 1.18)
            ax.grid(axis="x", alpha=0.3)
            ax.invert_yaxis()

        from matplotlib.patches import Patch
        handles = [Patch(facecolor=VARIANT_COLORS[v], edgecolor="black",
                         label=VARIANT_LABELS[v]) for v in sorted(VARIANT_COLORS)]
        handles.append(Patch(facecolor="white", edgecolor="black",
                             hatch="", label="klon_first"))
        handles.append(Patch(facecolor="white", edgecolor="black",
                             hatch="///", label="klev_first"))
        axes[0].legend(handles=handles, loc="lower right", fontsize=5,
                       ncol=2, framealpha=0.9)

        _save(fig, out_dir, f"metrics_combined_nlev{nlev}_{target}")
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot cost model metrics")
    parser.add_argument("--input", default="metrics.csv", help="CSV path")
    parser.add_argument("-o", "--outdir", default="plots", help="Output dir")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Not found: {args.input}"); sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.input)
    print(f"Read {len(df)} rows from {args.input}")

    for metric, title in METRICS:
        plot_metric(df, metric, title, args.outdir)

    plot_combined(df, args.outdir)
    print(f"\nAll plots in {args.outdir}/")


if __name__ == "__main__":
    main()