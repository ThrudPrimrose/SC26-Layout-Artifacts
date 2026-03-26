#!/usr/bin/env python3
"""
Plot z_v_grad_w benchmark results.
Each layout variant (V1-V4) is a separate violin within each group.

Usage:
    python plot_bench.py                          # looks for csvs in cwd
    python plot_bench.py --cpu cpu.csv --gpu gpu.csv
    python plot_bench.py --cpu cpu.csv            # cpu only
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

try:
    import pandas as pd
except ImportError:
    print("pip install pandas matplotlib"); sys.exit(1)

# ── style ──────────────────────────────────────────────────────────
VARIANT_LABELS = {1: "V1: none", 2: "V2: idx", 3: "V3: compute", 4: "V4: both"}
VARIANT_COLORS = {1: "#4c72b0", 2: "#dd8452", 3: "#55a868", 4: "#c44e52"}


def _violin(ax, data, positions, colors, width=0.7):
    """Draw violin plots manually so we control color per violin."""
    parts = ax.violinplot(
        data, positions=positions, showmeans=True,
        showmedians=True, showextrema=False, widths=width,
    )
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(colors[i])
        body.set_edgecolor("black")
        body.set_alpha(0.75)
    parts["cmeans"].set_color("black")
    parts["cmedians"].set_color("white")
    parts["cmedians"].set_linewidth(1.5)
    return parts


def plot_cpu(df, out_dir):
    """One figure per nlev.  Subplots = dist.  Groups = parallelization.  Violins = V1-V4."""
    for nlev, gn in df.groupby("nlev"):
        dists = sorted(gn["cell_dist"].unique())
        n_dist = len(dists)

        fig, axes = plt.subplots(n_dist, 1, figsize=(8, 4 * n_dist),
                                 squeeze=False, constrained_layout=True)
        fig.suptitle(f"CPU  |  nlev = {nlev}  |  nproma = {gn['nproma'].iloc[0]}",
                     fontsize=14, fontweight="bold")

        for row, dist in enumerate(dists):
            ax = axes[row, 0]
            gd = gn[gn["cell_dist"] == dist]
            pars = sorted(gd["parallelization"].unique())

            tick_positions = []
            tick_labels = []
            group_offset = 0
            spacing = 5  # gap between groups

            for pi, par in enumerate(pars):
                gp = gd[gd["parallelization"] == par]
                variants = sorted(gp["variant"].unique())
                data_list = []
                pos_list = []
                col_list = []
                for vi, v in enumerate(variants):
                    times = gp[gp["variant"] == v]["time_s"].values * 1e3
                    data_list.append(times)
                    pos_list.append(group_offset + vi)
                    col_list.append(VARIANT_COLORS[v])

                _violin(ax, data_list, pos_list, col_list, width=0.8)
                center = group_offset + (len(variants) - 1) / 2.0
                tick_positions.append(center)
                tick_labels.append(par)
                group_offset += len(variants) + spacing

            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, fontsize=10)
            ax.set_ylabel("time [ms]")
            ax.set_title(f"dist = {dist}", fontsize=11)
            ax.grid(axis="y", alpha=0.3)

        # legend
        from matplotlib.patches import Patch
        handles = [Patch(facecolor=VARIANT_COLORS[v], edgecolor="black",
                         label=VARIANT_LABELS[v]) for v in sorted(VARIANT_COLORS)]
        axes[0, 0].legend(handles=handles, loc="upper right", fontsize=8,
                          ncol=2, framealpha=0.9)

        path = os.path.join(out_dir, f"cpu_nlev{nlev}.pdf")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"  saved {path}")


def plot_gpu(df, out_dir):
    """One figure per (nlev, dist).  X-axis = config_label.  Violins = V1-V4."""
    for (nlev, dist), gnd in df.groupby(["nlev", "cell_dist"]):
        configs = sorted(gnd["config_label"].unique(),
                         key=lambda c: gnd[gnd["config_label"] == c]["time_s"].median())

        n_cfg = len(configs)
        fig_w = max(10, 0.9 * n_cfg)
        fig, ax = plt.subplots(1, 1, figsize=(fig_w, 5), constrained_layout=True)
        fig.suptitle(
            f"GPU  |  nlev = {nlev}  |  dist = {dist}  |  nproma = {gnd['nproma'].iloc[0]}",
            fontsize=13, fontweight="bold")

        tick_positions = []
        tick_labels = []
        group_offset = 0
        spacing = 2

        for ci, cfg in enumerate(configs):
            gc = gnd[gnd["config_label"] == cfg]
            variants = sorted(gc["variant"].unique())
            data_list = []
            pos_list = []
            col_list = []
            for vi, v in enumerate(variants):
                times = gc[gc["variant"] == v]["time_s"].values * 1e3
                if len(times) == 0:
                    continue
                data_list.append(times)
                pos_list.append(group_offset + vi)
                col_list.append(VARIANT_COLORS[v])

            if data_list:
                _violin(ax, data_list, pos_list, col_list, width=0.8)
            center = group_offset + (len(variants) - 1) / 2.0
            tick_positions.append(center)
            tick_labels.append(cfg)
            group_offset += len(variants) + spacing

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=7, rotation=45, ha="right")
        ax.set_ylabel("time [ms]")
        ax.grid(axis="y", alpha=0.3)

        from matplotlib.patches import Patch
        handles = [Patch(facecolor=VARIANT_COLORS[v], edgecolor="black",
                         label=VARIANT_LABELS[v]) for v in sorted(VARIANT_COLORS)]
        ax.legend(handles=handles, loc="upper right", fontsize=8,
                  ncol=2, framealpha=0.9)

        path = os.path.join(out_dir, f"gpu_nlev{nlev}_{dist}.pdf")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"  saved {path}")


def plot_gpu_summary(df, out_dir):
    """One figure per nlev: bar chart of median time per (config, variant), averaged over dists."""
    for nlev, gn in df.groupby("nlev"):
        med = gn.groupby(["config_label", "variant"])["time_s"].median().reset_index()
        med["time_ms"] = med["time_s"] * 1e3

        configs = sorted(med["config_label"].unique(),
                         key=lambda c: med[med["config_label"] == c]["time_ms"].mean())
        variants = sorted(med["variant"].unique())
        n_v = len(variants)
        x = np.arange(len(configs))
        bar_w = 0.8 / n_v

        fig, ax = plt.subplots(1, 1, figsize=(max(10, 0.7 * len(configs)), 5),
                               constrained_layout=True)
        fig.suptitle(f"GPU summary (median across dists)  |  nlev = {nlev}",
                     fontsize=13, fontweight="bold")

        for vi, v in enumerate(variants):
            vals = []
            for cfg in configs:
                row = med[(med["config_label"] == cfg) & (med["variant"] == v)]
                vals.append(row["time_ms"].values[0] if len(row) else 0)
            ax.bar(x + vi * bar_w, vals, bar_w, color=VARIANT_COLORS[v],
                   edgecolor="black", linewidth=0.3, label=VARIANT_LABELS[v])

        ax.set_xticks(x + bar_w * (n_v - 1) / 2)
        ax.set_xticklabels(configs, fontsize=7, rotation=45, ha="right")
        ax.set_ylabel("median time [ms]")
        ax.legend(fontsize=8, ncol=2)
        ax.grid(axis="y", alpha=0.3)

        path = os.path.join(out_dir, f"gpu_summary_nlev{nlev}.pdf")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"  saved {path}")


def main():
    parser = argparse.ArgumentParser(description="Plot z_v_grad_w benchmarks")
    parser.add_argument("--cpu", default="z_v_grad_w_cpu.csv",
                        help="CPU CSV path (default: z_v_grad_w_cpu.csv)")
    parser.add_argument("--gpu", default="z_v_grad_w_gpu.csv",
                        help="GPU CSV path (default: z_v_grad_w_gpu.csv)")
    parser.add_argument("-o", "--outdir", default="plots",
                        help="Output directory (default: plots/)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    have_cpu = os.path.isfile(args.cpu)
    have_gpu = os.path.isfile(args.gpu)

    if not have_cpu and not have_gpu:
        print(f"No CSV files found ({args.cpu}, {args.gpu}). Run benchmarks first.")
        sys.exit(1)

    if have_cpu:
        print(f"Reading {args.cpu} ...")
        df_cpu = pd.read_csv(args.cpu)
        print(f"  {len(df_cpu)} rows, variants {sorted(df_cpu['variant'].unique())}")
        plot_cpu(df_cpu, args.outdir)

    if have_gpu:
        print(f"Reading {args.gpu} ...")
        df_gpu = pd.read_csv(args.gpu)
        print(f"  {len(df_gpu)} rows, configs {df_gpu['config_label'].nunique()}")
        plot_gpu(df_gpu, args.outdir)
        plot_gpu_summary(df_gpu, args.outdir)

    print(f"\nAll plots in {args.outdir}/")


if __name__ == "__main__":
    main()