#!/usr/bin/env python3
"""Read benchmark CSVs and produce a violin + trace plot with statistical tests.

Follows Hoefler & Belli (SC'15) recommendations:
  - Nonparametric tests (data is typically non-normal on HPC systems)
  - Bootstrap CIs of the median
  - Mann-Whitney U for pairwise comparison against baseline
  - Kruskal-Wallis omnibus test when >2 variants
  - Effect size (rank-biserial r) for practical significance
"""
import argparse
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def bootstrap_median_ci(data, n_boot=10_000, alpha=0.05, rng=None):
    """Return (median, ci_lo, ci_hi) via percentile bootstrap."""
    if rng is None:
        rng = np.random.default_rng(42)
    data = np.asarray(data)
    boot_medians = np.array([
        np.median(rng.choice(data, size=len(data), replace=True))
        for _ in range(n_boot)
    ])
    lo = np.percentile(boot_medians, 100 * alpha / 2)
    hi = np.percentile(boot_medians, 100 * (1 - alpha / 2))
    return np.median(data), lo, hi


def rank_biserial_r(u_stat, n1, n2):
    """Rank-biserial correlation as effect size for Mann-Whitney U.
    r = 1 - 2U/(n1*n2).  |r| ~ 0.1 small, 0.3 medium, 0.5 large."""
    return 1.0 - (2.0 * u_stat) / (n1 * n2)


def pairwise_tests(datasets):
    """Run Mann-Whitney U (two-sided) of each variant vs the first (baseline).
    Returns list of dicts with test results."""
    baseline_label, baseline_times = datasets[0]
    base = np.asarray(baseline_times)
    results = []
    for i, (label, times) in enumerate(datasets):
        arr = np.asarray(times)
        med, ci_lo, ci_hi = bootstrap_median_ci(arr)
        entry = dict(
            label=label, n=len(arr),
            median=med, ci_lo=ci_lo, ci_hi=ci_hi,
            min=np.min(arr), max=np.max(arr), std=np.std(arr),
        )
        if i == 0:
            entry.update(U=None, p=None, r_effect=None, significant=None)
        else:
            u_stat, p_val = stats.mannwhitneyu(base, arr, alternative='two-sided')
            r = rank_biserial_r(u_stat, len(base), len(arr))
            entry.update(U=u_stat, p=p_val, r_effect=r,
                         significant=(p_val < 0.05))
        results.append(entry)
    return results


def omnibus_test(datasets):
    """Kruskal-Wallis H-test across all variants."""
    arrays = [np.asarray(t) for _, t in datasets]
    if len(arrays) < 2:
        return None, None
    h_stat, p_val = stats.kruskal(*arrays)
    return h_stat, p_val


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_plot(datasets, output, test_results, kw_h, kw_p):
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
    ax.set_title("Condensation kernel \u2014 isolated process runs", fontsize=13, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    # Annotate medians + speedup + significance star
    medians = [r['median'] for r in test_results]
    for i, r in enumerate(test_results):
        suffix = ""
        if i > 0 and medians[0] > 0:
            ratio = medians[0] / r['median']
            star = "*" if r['significant'] else ""
            suffix = f"  ({ratio:.2f}x vs {labels[0]}){star}"
        ax.annotate(f"med={r['median']:.2f} ms{suffix}",
                    xy=(positions[i], r['median']),
                    xytext=(15, 12 + i * 8), textcoords='offset points',
                    fontsize=9, color=colors[i % len(colors)], fontweight='bold',
                    arrowprops=dict(arrowstyle='->', color=colors[i % len(colors)], lw=1.2))

        # Draw bootstrap CI as a horizontal error bar
        ax.plot([positions[i] - 0.25, positions[i] - 0.25],
                [r['ci_lo'], r['ci_hi']],
                color=colors[i % len(colors)], linewidth=2.5, alpha=0.8,
                solid_capstyle='round')

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


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_report(test_results, kw_h, kw_p):
    medians = [r['median'] for r in test_results]
    base_med = medians[0]

    hdr = (f"{'Variant':<28s} {'n':>5s} {'median':>8s} "
           f"{'95% CI':>16s} {'min':>8s} {'max':>8s} {'std':>8s} "
           f"{'vs base':>8s} {'U':>10s} {'p':>10s} {'r':>7s} {'sig':>5s}")
    print(hdr)
    print("-" * len(hdr))

    for i, r in enumerate(test_results):
        ratio_str = "1.000x" if i == 0 else f"{base_med / r['median']:.3f}x"
        ci_str = f"[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]"
        if i == 0:
            u_str = p_str = r_str = sig_str = ""
        else:
            u_str = f"{r['U']:.0f}"
            p_str = f"{r['p']:.2e}" if r['p'] >= 1e-10 else f"{r['p']:.1e}"
            r_str = f"{r['r_effect']:.3f}"
            sig_str = "***" if r['p'] < 0.001 else ("**" if r['p'] < 0.01 else ("*" if r['p'] < 0.05 else "n.s."))
        print(f"{r['label']:<28s} {r['n']:5d} {r['median']:8.3f} "
              f"{ci_str:>16s} {r['min']:8.3f} {r['max']:8.3f} {r['std']:8.3f} "
              f"{ratio_str:>8s} {u_str:>10s} {p_str:>10s} {r_str:>7s} {sig_str:>5s}")

    print()
    if kw_h is not None:
        print(f"Kruskal-Wallis H = {kw_h:.3f},  p = {kw_p:.2e}"
              f"  {'=> at least one variant differs (p<0.05)' if kw_p < 0.05 else '=> no significant difference across variants'}")

    # CI overlap check
    print()
    base = test_results[0]
    for r in test_results[1:]:
        overlap = base['ci_lo'] <= r['ci_hi'] and r['ci_lo'] <= base['ci_hi']
        print(f"  Median CI overlap  {base['label']:<20s} vs {r['label']:<20s}: "
              f"{'YES (overlapping)' if overlap else 'NO  (separated)'}"
              f"  -- but see Mann-Whitney p={r['p']:.2e}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Plot benchmark results with statistical tests")
    parser.add_argument('--csv', nargs='+', required=True,
                        help="One or more CSV files from condense_bench.py")
    parser.add_argument('--output', type=str, default='condense_violin.pdf',
                        help="Output image path")
    parser.add_argument('--n-boot', type=int, default=10_000,
                        help="Number of bootstrap resamples (default: 10000)")
    args = parser.parse_args()

    datasets = []
    label_map = {'3d': '3D zsolqa[5,5,klon]', 'split': 'Split zsolqa_i_j'}
    for path in args.csv:
        variant, times = read_csv(path)
        label = label_map.get(variant, variant)
        datasets.append((label, times))

    # Statistical tests
    test_results = pairwise_tests(datasets)
    kw_h, kw_p = omnibus_test(datasets)

    # Plot
    make_plot(datasets, args.output, test_results, kw_h, kw_p)

    # Text report
    print_report(test_results, kw_h, kw_p)


if __name__ == "__main__":
    main()