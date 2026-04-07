#!/usr/bin/env python3
"""
plot_zaxpy_sweep.py
Violin bandwidth plots for indirect scatter-accumulate benchmark.

    y[σ(i)] += x[i],  i = 0 … nx−1

Two violins per distribution group:
    Orange = Baseline (original index order — sequential reads, scattered writes)
    Blue   = Best write-optimized layout (sorted or NUMA-partitioned)

GPU picks best of {aos_sorted, soa_sorted}.
CPU picks best of {aos_sorted, soa_sorted, aos_partitioned, soa_partitioned}.

Usage:
    python plot_zaxpy_sweep.py --amd-dir results/beverin --nv-dir results/daint --with-cpu --add-peak
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import pandas as pd, numpy as np, argparse, sys, os

# ── Two-violin scheme ────────────────────────────────────────────────────────
VCOL = {"baseline": "#e67e22", "optimized": "#2980b9"}
VLAB = {
    "baseline":  "Original Order",
    "optimized": "Sorted (Shuffle)",
}
VIOLIN_KEYS = ["baseline", "optimized"]

# Which CSV variants feed each violin
BASELINE_CANDIDATES  = ["aos_scatter", "soa_scatter"]
GPU_OPT_CANDIDATES   = ["aos_sorted", "soa_sorted"]
CPU_OPT_CANDIDATES   = ["aos_sorted", "soa_sorted",
                         "aos_partitioned", "soa_partitioned"]

# Bytes per element (minimum traffic model)
BYTES_PER_ELEM = {
    "aos_scatter": 52, "aos_sorted": 56,
    "soa_scatter": 52, "soa_sorted": 56,
    "aos_partitioned": 56, "soa_partitioned": 56,
}

DIST_ORDER = ["uniform", "normal", "qe"]
DIST_LABEL = {
    "uniform": "Uniform",
    "normal":  r"Normal",
    "qe":      "Exact",
}

STREAM_PEAK = {
    "MI300A Zen CPU":  1228    * 1e-3,
    "Grace CPU":       1700.62 * 1e-3,
    "MI300A GPU":      4294    * 1e-3,
    "GH200 GPU":       3780    * 1e-3,
}

GPU_W, GPU_H = 4.2*(1/1.7), 4.5*(1/1.7)
CPU_W, CPU_H = 2.8, 1.6

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


def best_config_bw(df, variant, max_threads_only=False):
    """Best (tpb, coarsen) for a single variant by highest median BW."""
    sub = df[df["variant"] == variant].copy()
    if sub.empty:
        return np.array([]), -1.0, variant
    if max_threads_only and "tpb" in sub.columns:
        sub = sub[sub["tpb"] == sub["tpb"].max()]
    sub["bw"] = compute_bandwidth(sub["time_ms"].values,
                                  sub["nx"].values, variant)
    medians = sub.groupby(["tpb", "coarsen"])["bw"].median()
    if medians.empty:
        return np.array([]), -1.0, variant
    best = medians.idxmax()
    vals = sub[(sub["tpb"] == best[0]) & (sub["coarsen"] == best[1])]["bw"].values
    return vals, np.median(vals), variant


def best_across_variants(df, candidate_list, max_threads_only=False):
    """Pick the variant+config with highest median BW from a candidate list."""
    best_vals, best_med, best_vk = np.array([]), -1.0, ""
    for vk in candidate_list:
        vals, med, _ = best_config_bw(df, vk, max_threads_only)
        if med > best_med:
            best_vals, best_med, best_vk = vals, med, vk
    return best_vals, best_vk


def select_dist(df, dist_key, tiled):
    if dist_key == "qe":
        pat = "qe_tiled" if tiled else "qe"
        return df[df["experiment"] == pat]
    elif dist_key == "uniform":
        prefix = "uniform_tiled_s" if tiled else "uniform_s"
        return df[df["experiment"].str.startswith(prefix)]
    else:
        prefix = "normal_tiled_s" if tiled else "normal_s"
        return df[df["experiment"].str.startswith(prefix)]


def load_dir(d, cpu=False):
    sfx = "_cpu" if cpu else ""
    out = {}
    for key, base in [("small", "zaxpy_sweep_small"),
                      ("1gb",   "zaxpy_sweep_1gb")]:
        p = os.path.join(d, f"{base}{sfx}.csv")
        try:
            out[key] = pd.read_csv(p)
            print(f"  Loaded {p}  ({len(out[key])} rows)")
        except FileNotFoundError:
            print(f"  [WARN] {p} not found")
            out[key] = None
    return out


# ── Subplot painter ──────────────────────────────────────────────────────────

def paint_subplot(ax, df, peak_label, tiled, is_cpu, add_peak=False):
    opt_candidates = CPU_OPT_CANDIDATES if is_cpu else GPU_OPT_CANDIDATES

    group_spacing = 3   # 2 violins + 1 gap
    positions, data_all, col_all = [], [], []
    medians_annot = []
    xticks, xlabels = [], []

    for di, dist_key in enumerate(DIST_ORDER):
        df_slice = select_dist(df, dist_key, tiled)

        for vi, vk in enumerate(VIOLIN_KEYS):
            if vk == "baseline":
                vals, chosen = best_across_variants(df_slice, BASELINE_CANDIDATES, is_cpu)
            else:
                vals, chosen = best_across_variants(df_slice, opt_candidates, is_cpu)

            vals = remove_outliers(vals)
            pos = di * group_spacing + vi
            if len(vals) > 0:
                data_all.append(vals)
                positions.append(pos)
                col_all.append(VCOL[vk])
                medians_annot.append((pos, np.median(vals), np.min(vals), vk))

        xticks.append(di * group_spacing + 0.5)
        xlabels.append(DIST_LABEL[dist_key])

    # Y-axis
    all_flat = np.concatenate(data_all) if data_all else np.array([0.0])
    mx = float(np.max(all_flat))
    loc = MaxNLocator(nbins=5, min_n_ticks=5)
    ticks = loc.tick_values(0.0, mx * 1.14)
    ticks = ticks[ticks >= 0]
    if len(ticks) > 6:
        ticks = ticks[:6]
    top = ticks[-1] * 1.06 if len(ticks) else mx * 1.15

    # Violins
    if data_all:
        parts = ax.violinplot(data_all, positions=positions,
                              showmeans=True, showmedians=True,
                              showextrema=False, widths=1.2)
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
        ax.text(0.03, 0.97, f"{peak:.2f} TB/s STREAM",
                transform=ax.transAxes, ha='left', va='top',
                fontsize=8, color='dimgray')

    # Horizontal peak line
    if add_peak and peak_label in STREAM_PEAK:
        peak = STREAM_PEAK[peak_label]
        ax.axhline(y=peak, color='dimgray', ls='--', lw=0.8, alpha=0.7)

    # % annotations
    if peak_label in STREAM_PEAK:
        peak = STREAM_PEAK[peak_label]
        yrange = top
        off = 0.045 * yrange
        for pos, med, vmin, vk in medians_annot:
            pct = 100.0 * med / peak
            color = VCOL[vk]
            if vk == "baseline":
                xoff, ha = 0.3, 'right'
            else:
                xoff, ha = -0.3, 'left'
            ax.text(pos + xoff, vmin - off, f'{pct:.0f}%',
                    ha=ha, va='top',
                    fontsize=8.5, color=color, fontweight='bold')

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amd-dir", default="results/beverin")
    ap.add_argument("--nv-dir",  default="results/daint")
    ap.add_argument("--with-cpu", action="store_true")
    ap.add_argument("--add-peak", action="store_true")
    args = ap.parse_args()

    print("Loading GPU AMD:")
    gpu_amd = load_dir(args.amd_dir, cpu=False)
    print("Loading GPU NV:")
    gpu_nv  = load_dir(args.nv_dir, cpu=False)

    have_cpu = args.with_cpu
    cpu_amd, cpu_nv = {}, {}
    if have_cpu:
        print("Loading CPU AMD:")
        cpu_amd = load_dir(args.amd_dir, cpu=True)
        print("Loading CPU NV:")
        cpu_nv  = load_dir(args.nv_dir, cpu=True)
        if any(cpu_amd.get(k) is None for k in ["small","1gb"]) or \
           any(cpu_nv.get(k) is None for k in ["small","1gb"]):
            print("WARN: CPU CSVs incomplete, GPU-only.")
            have_cpu = False

    plt.rcParams.update({
        "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 9,
    })

    scales = [
        ("small", False, gpu_amd.get("small"), gpu_nv.get("small"),
                         cpu_amd.get("small"), cpu_nv.get("small")),
        ("1gb",   True,  gpu_amd.get("1gb"),   gpu_nv.get("1gb"),
                         cpu_amd.get("1gb"),   cpu_nv.get("1gb")),
    ]

    for scale_name, tiled, ga, gn, ca, cn in scales:
        if ga is None or gn is None:
            print(f"Skipping {scale_name}: missing GPU data"); continue

        # Build grid: (label, df, peak_key, is_cpu, subplot_w, subplot_h)
        rows = []
        if have_cpu and ca is not None and cn is not None:
            rows.append([
                ("MI300A Zen CPU", ca, "MI300A Zen CPU", True, CPU_W, CPU_H),
                ("Grace CPU",      cn, "Grace CPU",      True, CPU_W, CPU_H),
            ])
        rows.append([
            ("MI300A GPU", ga, "MI300A GPU", False, GPU_W, GPU_H),
            ("GH200 GPU",  gn, "GH200 GPU", False, GPU_W, GPU_H),
        ])

        nrows = len(rows)
        ncols = 2
        col_w = max(rows[r][c][4] for r in range(nrows) for c in range(ncols))
        row_heights = [rows[r][0][5] for r in range(nrows)]
        fig_w = col_w * ncols
        fig_h = sum(row_heights)

        fig, axes = plt.subplots(
            nrows, ncols, figsize=(fig_w, fig_h), squeeze=False,
        )

        title_scale = "nat20 original" if not tiled else "\u22481 GiB tiled"
        fig.suptitle(
            f"Indirect Scatter-Accumulate",
            fontsize=10, y=0.998,
        )
        fig.text(0.5, 0.95,
                 r"Original: $y[\sigma(i)]\ {+}{=}\ x[i]$"
                 "\n"
                 r"Write-coalesced: $y[\sigma(\pi(i))]\ {+}{=}\ x[\pi(i)]$,"
                 r"$\;\pi = \mathrm{argsort}(\sigma)$",
                 ha='center', va='top', fontsize=9.5, color='dimgray')

        for ri, row_data in enumerate(rows):
            for ci, (label, df, peak_label, is_cpu, sw, sh) in enumerate(row_data):
                ax = axes[ri, ci]
                ax.set_title(label, fontsize=9, pad=3)
                if ci == 0:
                    ax.set_ylabel("BW [TB/s]", fontsize=8)
                paint_subplot(ax, df, peak_label, tiled, is_cpu,
                              add_peak=args.add_peak)

        # Legend
        handles = [Patch(facecolor=VCOL[v], edgecolor="black", label=VLAB[v])
                   for v in VIOLIN_KEYS]
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.5, -0.01), ncol=2,
                   framealpha=0.9, columnspacing=1.0, fontsize=8)

        fig.subplots_adjust(left=0.04, right=0.99,
                            top=0.98, bottom=0.12, hspace=0.05, wspace=0.05)
        fig.tight_layout(rect=[0, 0.05, 1, 0.95])
        sfx = "_w_stream_peak" if args.add_peak else ""
        stem = f"zaxpy_violins_{scale_name}{sfx}"
        fig.savefig(f"{stem}.png", dpi=200, bbox_inches='tight')
        fig.savefig(f"{stem}.pdf", dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {stem}.png / .pdf")

    # ── Summary table ────────────────────────────────────────────────────
    print(f"\n{'Scale':<8} {'Platform':<22} {'Dist':<10} {'Violin':<12} "
          f"{'Chosen variant':<22} {'Median TB/s':<12} {'% Peak'}")
    print("-" * 96)
    for scale_name, tiled, ga, gn, ca, cn in scales:
        entries = []
        if ga is not None:
            entries.append(("MI300A GPU", ga, "MI300A GPU", False))
        if gn is not None:
            entries.append(("GH200 GPU",  gn, "GH200 GPU", False))
        if have_cpu and ca is not None:
            entries.append(("MI300A Zen CPU", ca, "MI300A Zen CPU", True))
        if have_cpu and cn is not None:
            entries.append(("Grace CPU",      cn, "Grace CPU",      True))

        for label, df, pk, is_cpu in entries:
            peak = STREAM_PEAK.get(pk, 1.0)
            opt_cands = CPU_OPT_CANDIDATES if is_cpu else GPU_OPT_CANDIDATES
            for dk in DIST_ORDER:
                sl = select_dist(df, dk, tiled)
                # Baseline
                vals, chosen = best_across_variants(sl, BASELINE_CANDIDATES, is_cpu)
                if len(vals) > 0:
                    med = np.median(vals)
                    pct = 100*med/peak
                    print(f"{scale_name:<8} {label:<22} {dk:<10} {'baseline':<12} "
                          f"{chosen:<22} {med:.4f}      {pct:.1f}%")
                # Optimized
                vals, chosen = best_across_variants(sl, opt_cands, is_cpu)
                if len(vals) > 0:
                    med = np.median(vals)
                    pct = 100*med/peak
                    print(f"{scale_name:<8} {label:<22} {dk:<10} {'optimized':<12} "
                          f"{chosen:<22} {med:.4f}      {pct:.1f}%")


if __name__ == "__main__":
    main()