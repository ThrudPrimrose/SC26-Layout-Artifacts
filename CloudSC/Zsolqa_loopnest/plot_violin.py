#!/usr/bin/env python3
"""Read benchmark CSVs and produce a violin + trace plot."""
import argparse
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_csv(path):
    """Read a benchmark CSV, return (variant_name, list_of_time_ms)."""
    times = []
    variant = None
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            variant = row['variant']
            times.append(float(row['time_ms']))
    return variant, times


def make_plot(datasets, output):
    """
    datasets: list of (label, times_ms_array)
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                             gridspec_kw={'width_ratios': [3, 1]})
    colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

    # ---- Left: violin + box ----
    ax = axes[0]
    data = [np.array(t) for _, t in datasets]
    labels = [l for l, _ in datasets]
    positions = list(range(1, len(data) + 1))

    parts = ax.violinplot(data, positions=positions,
                          showmeans=False, showmedians=False, showextrema=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i % len(colors)])
        pc.set_edgecolor('black')
        pc.set_linewidth(0.8)
        pc.set_alpha(0.7)

    bp = ax.boxplot(data, positions=positions, widths=0.15, patch_artist=True,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.4))
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[i % len(colors)])
        patch.set_alpha(0.9)
    for element in ['whiskers', 'caps']:
        for line in bp[element]:
            line.set_color('black')
            line.set_linewidth(1.0)
    for line in bp['medians']:
        line.set_color('white')
        line.set_linewidth(2.0)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("Time (ms)", fontsize=12)
    ax.set_title("Condensation kernel — isolated process runs", fontsize=13, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    # Annotate medians + speedup
    medians = [np.median(d) for d in data]
    for i, (d, med) in enumerate(zip(data, medians)):
        suffix = ""
        if i > 0 and medians[0] > 0:
            ratio = medians[0] / med
            suffix = f"  ({ratio:.2f}x vs {labels[0]})"
        ax.annotate(f"med={med:.2f} ms{suffix}",
                    xy=(positions[i], med),
                    xytext=(15, 12 + i * 8), textcoords='offset points',
                    fontsize=9, color=colors[i % len(colors)], fontweight='bold',
                    arrowprops=dict(arrowstyle='->', color=colors[i % len(colors)], lw=1.2))

    # ---- Right: per-iteration trace ----
    ax2 = axes[1]
    markers = ['o', 's', '^', 'D']
    for i, (label, times) in enumerate(datasets):
        iters = np.arange(1, len(times) + 1)
        ax2.plot(iters, times,
                 marker=markers[i % len(markers)], markersize=2.5,
                 linewidth=0.8, alpha=0.7, color=colors[i % len(colors)],
                 label=label)
    ax2.set_xlabel("Iteration", fontsize=11)
    ax2.set_ylabel("Time (ms)", fontsize=11)
    ax2.set_title("Per-iteration trace", fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to: {output}")

    # Print summary table
    print(f"\n{'Variant':<28s} {'median':>8s} {'min':>8s} {'max':>8s} {'std':>8s} {'vs first':>10s}")
    print("-" * 72)
    for i, (label, times) in enumerate(datasets):
        arr = np.array(times)
        ratio_str = "1.000x" if i == 0 else f"{medians[0]/np.median(arr):.3f}x"
        print(f"{label:<28s} {np.median(arr):8.3f} {np.min(arr):8.3f} "
              f"{np.max(arr):8.3f} {np.std(arr):8.3f} {ratio_str:>10s}")


def main():
    parser = argparse.ArgumentParser(description="Plot benchmark results")
    parser.add_argument('--csv', nargs='+', required=True,
                        help="One or more CSV files from condense_bench.py")
    parser.add_argument('--output', type=str, default='condense_violin.png',
                        help="Output image path (default: condense_violin.png)")
    args = parser.parse_args()

    datasets = []
    label_map = {'3d': '3D zsolqa[5,5,klon]', 'split': 'Split zsolqa_i_j'}
    for path in args.csv:
        variant, times = read_csv(path)
        label = label_map.get(variant, variant)
        datasets.append((label, times))

    make_plot(datasets, args.output)


if __name__ == "__main__":
    main()