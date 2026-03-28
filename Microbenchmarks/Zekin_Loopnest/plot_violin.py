#!/usr/bin/env python3
"""
plot_violin.py — Violin plot comparing selected kernel variants from benchmark CSV.

Usage:
    python plot_violin.py timings.csv "A1-01" "B-02"
    python plot_violin.py timings.csv "A1-09" "B-05" --nlev 90
    python plot_violin.py timings.csv "A1-01" "A1-09" "B-02" "B-05" --nlev 90 -o compare.pdf

Reads the CSV written by benchmark.cu and draws a violin plot of per-iteration
timings for each variant you name (matched by prefix in the variant column).
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'serif'

def main():
    parser = argparse.ArgumentParser(description="Violin plot of kernel timings")
    parser.add_argument("csv", help="Path to timings.csv")
    parser.add_argument("variants", nargs="+", help="Variant prefixes to compare (e.g. 'A1-01' 'B-02')")
    parser.add_argument("--nlev", type=int, default=None, help="Filter by nlev (default: all)")
    parser.add_argument("-o", "--output", default=None, help="Save plot to file (e.g. plot.pdf)")
    parser.add_argument("--title", default=None, help="Plot title")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    # Strip whitespace from variant names (C printf may pad)
    df["variant"] = df["variant"].str.strip()

    if args.nlev is not None:
        df = df[df["nlev"] == args.nlev]
        if df.empty:
            print(f"No data for nlev={args.nlev}. Available: {pd.read_csv(args.csv)['nlev'].unique()}")
            sys.exit(1)

    # Match each requested prefix against the variant column
    matched = []
    for prefix in args.variants:
        mask = df["variant"].str.startswith(prefix)
        if not mask.any():
            print(f"Warning: no rows match prefix '{prefix}'. Available variants:")
            for v in sorted(df["variant"].unique()):
                print(f"  {v}")
            sys.exit(1)
        # Take the first matching variant name (exact unique match)
        vname = df.loc[mask, "variant"].iloc[0]
        matched.append(vname)

    # Deduplicate while preserving order
    seen = set()
    unique_matched = []
    for v in matched:
        if v not in seen:
            seen.add(v)
            unique_matched.append(v)
    matched = unique_matched

    # Collect data per variant
    data = []
    labels = []
    for vname in matched:
        times = df.loc[df["variant"] == vname, "time_ms"].values
        data.append(times)
        # Shorten label: take the variant ID (first token)
        short = vname.split()[0] if " " in vname else vname
        labels.append(f"{short}\n({len(times)} samples)")

    # Plot
    fig, ax = plt.subplots(figsize=(max(3, 1.8 * len(data)), 5))

    parts = ax.violinplot(data, positions=range(len(data)), showmeans=True, showmedians=True, showextrema=False)

    # Color A variants blue, B variants orange
    for i, (body, vname) in enumerate(zip(parts["bodies"], matched)):
        if vname.startswith("A"):
            body.set_facecolor("#4C72B0")
        else:
            body.set_facecolor("#DD8452")
        body.set_alpha(0.7)

    parts["cmeans"].set_color("black")
    parts["cmedians"].set_color("red")

    # Overlay individual points (jittered)
    import numpy as np
    for i, d in enumerate(data):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(d))
        ax.scatter(np.full(len(d), i) + jitter, d, s=8, alpha=0.3, color="gray", zorder=3)

    ax.set_xticks(range(len(data)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Time (ms)")
    ax.grid(axis="y", alpha=0.3)

    title = args.title
    if title is None:
        nlev_str = f"nlev={args.nlev}" if args.nlev else "all nlev"
        title = f"Kernel Timing Distribution ({nlev_str})"
    ax.set_title(title)

    # Legend
    from matplotlib.patches import Patch
    legend_items = [Patch(facecolor="#4C72B0", alpha=0.7, label="Layout A [nproma,nlev,nblks]"),
                    Patch(facecolor="#DD8452", alpha=0.7, label="Layout B [nlev,nproma,nblks]")]
    ax.legend(handles=legend_items, loc="upper right", fontsize=8)

    # Stats annotation
    for i, d in enumerate(data):
        med = np.median(d)
        ax.annotate(f"{med:.3f}", xy=(i, med), fontsize=7, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points", color="red")

    plt.tight_layout()
    if args.output:
        fig.savefig(args.output, dpi=200, bbox_inches="tight")
        print(f"Saved to {args.output}")
    else:
        plt.show()

if __name__ == "__main__":
    main()