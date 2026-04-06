#!/usr/bin/env python3
"""
plot_zaxpy_sweep.py
Violin bandwidth plots for indirect zaxpy benchmark.

Two PNGs: "small" (nat20 original) and "1gb" (tiled).
Columns: AMD (left), NVIDIA (right).
Rows: GPU only (1 row) or CPU + GPU (2 rows) when --cpu-amd-dir given.

Usage:
    python plot_zaxpy_sweep.py --amd-dir results/beverin --nv-dir results/daint

    # With CPU row:
python plot_zaxpy_sweep.py \
  --amd-dir results/beverin --nv-dir results/daint \
  --cpu-amd-dir results/beverin_cpu --cpu-nv-dir results/daint_cpu \
  --add-peak

python plot_zaxpy_sweep.py \
  --amd-dir . --nv-dir . \
  --cpu-amd-dir . --cpu-nv-dir . \
  --add-peak
  
Each directory should contain zaxpy_sweep_small.csv and zaxpy_sweep_1gb.csv.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import pandas as pd, numpy as np, argparse, sys, os

# ── Colour / label scheme ────────────────────────────────────────────────────
VCOL = {
    "aos_scatter":      "#e67e22",
    "aos_sorted":       "#2980b9",
    "soa_scatter":      "#c0392b",
    "soa_sorted":       "#27ae60",
    "aos_partitioned":  "#8e44ad",   # purple
    "soa_partitioned":  "#16a085",   # teal
}
VLAB = {
    "aos_scatter":      "AoS Scatter (original)",
    "aos_sorted":       "AoS Sorted",
    "soa_scatter":      "SoA Scatter",
    "soa_sorted":       "SoA Sorted",
    "aos_partitioned":  "AoS NUMA-part.",
    "soa_partitioned":  "SoA NUMA-part.",
}

# GPU has 4 variants; CPU has 6 (with partitioned)
GPU_VARIANTS = ["aos_scatter", "aos_sorted", "soa_scatter", "soa_sorted"]
CPU_VARIANTS = ["aos_scatter", "aos_sorted", "soa_scatter", "soa_sorted",
                "aos_partitioned", "soa_partitioned"]

# Bytes per element (minimum traffic model)
BYTES_PER_ELEM = {
    "aos_scatter": 52, "aos_sorted": 56,
    "soa_scatter": 52, "soa_sorted": 56,
    "aos_partitioned": 56, "soa_partitioned": 56,
}

DIST_ORDER = ["uniform", "normal", "qe"]
DIST_LABEL = {"uniform": "Uniform", "normal": "Normal", "qe": "QE"}

STREAM_PEAK = {
    "MI300A Zen CPU":  1228    * 1e-3,   # TB/s
    "Grace CPU":       1700.62 * 1e-3,
    "MI300A GPU":      4294    * 1e-3,
    "GH200 GPU":       3780    * 1e-3,
}

# Per-subplot dimensions (inches)
GPU_W, GPU_H = 3.6, 1.6
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


def best_config_bw(df, variant):
    sub = df[df["variant"] == variant].copy()
    if sub.empty:
        return np.array([])
    sub["bw"] = compute_bandwidth(sub["time_ms"].values,
                                  sub["nx"].values, variant)
    medians = sub.groupby(["tpb", "coarsen"])["bw"].median()
    if medians.empty:
        return np.array([])
    best = medians.idxmax()
    vals = sub[(sub["tpb"] == best[0]) & (sub["coarsen"] == best[1])]["bw"].values
    return vals


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


def load_dir(d):
    """Load small and 1gb CSVs from a directory."""
    out = {}
    for key, fn in [("small", "zaxpy_sweep_small.csv"),
                    ("1gb",   "zaxpy_sweep_1gb.csv")]:
        p = os.path.join(d, fn)
        try:
            out[key] = pd.read_csv(p)
            print(f"  Loaded {p}  ({len(out[key])} rows)")
        except FileNotFoundError:
            print(f"  [WARN] {p} not found")
            out[key] = None
    return out


# ── Subplot painter ──────────────────────────────────────────────────────────

def paint_subplot(ax, df, peak_label, tiled, variant_list, add_peak=False):
    nvars = len(variant_list)
    group_spacing = nvars + 1.5

    positions, data_all, col_all = [], [], []
    medians_annot = []
    xticks, xlabels = [], []

    for di, dist_key in enumerate(DIST_ORDER):
        df_slice = select_dist(df, dist_key, tiled)

        for vi, vk in enumerate(variant_list):
            vals = best_config_bw(df_slice, vk)
            vals = remove_outliers(vals)
            pos = di * group_spacing + vi
            if len(vals) > 0:
                data_all.append(vals)
                positions.append(pos)
                col_all.append(VCOL[vk])
                medians_annot.append((pos, np.median(vals), vk))

        xticks.append(di * group_spacing + (nvars - 1) / 2.0)
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
                              showextrema=False, widths=0.9)
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
        ax.axhline(y=peak, color='red', ls='--', lw=1.2, alpha=0.7)

    # % annotations
    if peak_label in STREAM_PEAK:
        peak = STREAM_PEAK[peak_label]
        off = 0.045 * top
        for pos, med, vk in medians_annot:
            pct = 100.0 * med / peak
            ax.text(pos, med - off, f'{pct:.0f}%',
                    ha='center', va='top',
                    fontsize=7, color=VCOL[vk], fontweight='bold')


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amd-dir", default="results/beverin",
                    help="Dir with GPU AMD CSVs (default: results/beverin)")
    ap.add_argument("--nv-dir",  default="results/daint",
                    help="Dir with GPU NV CSVs (default: results/daint)")
    ap.add_argument("--cpu-amd-dir", default=None,
                    help="Dir with CPU AMD CSVs (enables 2×2 grid)")
    ap.add_argument("--cpu-nv-dir",  default=None,
                    help="Dir with CPU NV CSVs (enables 2×2 grid)")
    ap.add_argument("--add-peak", action="store_true",
                    help="Draw horizontal STREAM peak line")
    args = ap.parse_args()

    # ── Load ─────────────────────────────────────────────────────────────
    print("Loading GPU AMD:")
    gpu_amd = load_dir(args.amd_dir)
    print("Loading GPU NV:")
    gpu_nv  = load_dir(args.nv_dir)

    have_cpu = args.cpu_amd_dir is not None and args.cpu_nv_dir is not None
    cpu_amd, cpu_nv = {}, {}
    if have_cpu:
        print("Loading CPU AMD:")
        cpu_amd = load_dir(args.cpu_amd_dir)
        print("Loading CPU NV:")
        cpu_nv  = load_dir(args.cpu_nv_dir)
        if any(cpu_amd.get(k) is None for k in ["small","1gb"]) or \
           any(cpu_nv.get(k) is None for k in ["small","1gb"]):
            print("WARN: CPU CSVs incomplete, GPU-only.")
            have_cpu = False

    plt.rcParams.update({
        "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    })

    # ── One figure per scale ─────────────────────────────────────────────
    scales = [
        ("small", False, gpu_amd.get("small"), gpu_nv.get("small"),
                         cpu_amd.get("small"), cpu_nv.get("small")),
        ("1gb",   True,  gpu_amd.get("1gb"),   gpu_nv.get("1gb"),
                         cpu_amd.get("1gb"),   cpu_nv.get("1gb")),
    ]

    for scale_name, tiled, ga, gn, ca, cn in scales:
        if ga is None or gn is None:
            print(f"Skipping {scale_name}: missing GPU data"); continue

        # Build grid: list of (label, df, peak_key, variant_list, subplot_w, subplot_h)
        rows = []
        if have_cpu and ca is not None and cn is not None:
            rows.append([
                ("MI300A Zen CPU", ca, "MI300A Zen CPU", CPU_VARIANTS, CPU_W, CPU_H),
                ("Grace CPU",      cn, "Grace CPU",      CPU_VARIANTS, CPU_W, CPU_H),
            ])
        rows.append([
            ("MI300A GPU", ga, "MI300A GPU", GPU_VARIANTS, GPU_W, GPU_H),
            ("GH200 GPU",  gn, "GH200 GPU", GPU_VARIANTS, GPU_W, GPU_H),
        ])

        nrows = len(rows)
        ncols = 2

        # Compute figure size from per-row heights and per-col widths
        # Use max width across columns for each col (they're equal here)
        col_w = max(rows[0][0][4], rows[0][1][4])  # both same per row
        row_heights = [rows[r][0][5] for r in range(nrows)]
        fig_w = col_w * ncols
        fig_h = sum(row_heights)

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(fig_w, fig_h),
            squeeze=False,
            gridspec_kw={"height_ratios": row_heights},
        )

        title_scale = "nat20 original" if not tiled else "\u22481 GiB tiled"
        fig.suptitle(f"Indirect ZAXPY \u2014 {title_scale}",
                     fontsize=11, y=0.995)

        for ri, row_data in enumerate(rows):
            for ci, (label, df, peak_label, vlist, sw, sh) in enumerate(row_data):
                ax = axes[ri, ci]
                ax.set_title(label, fontsize=9, pad=3)
                if ci == 0:
                    ax.set_ylabel("BW [TB/s]", fontsize=8)
                paint_subplot(ax, df, peak_label, tiled, vlist,
                              add_peak=args.add_peak)

        # ── Legend ────────────────────────────────────────────────────────
        # Collect all variants used
        all_vk = GPU_VARIANTS if not have_cpu else CPU_VARIANTS
        handles = [Patch(facecolor=VCOL[v], edgecolor="black", label=VLAB[v])
                   for v in all_vk]
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.5, -0.01),
                   ncol=3 if len(all_vk) > 4 else len(all_vk),
                   framealpha=0.9, columnspacing=0.8, fontsize=7)

        fig.tight_layout(rect=[0, 0.06, 1, 0.97])
        sfx = "_w_stream_peak" if args.add_peak else ""
        stem = f"zaxpy_violins_{scale_name}{sfx}"
        fig.savefig(f"{stem}.png", dpi=200, bbox_inches='tight')
        fig.savefig(f"{stem}.pdf", dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {stem}.png / .pdf")

    # ── Summary table ────────────────────────────────────────────────────
    print(f"\n{'Scale':<8} {'Platform':<22} {'Dist':<10} {'Variant':<20} "
          f"{'Best (tpb,cf)':<16} {'Median TB/s':<12} {'% Peak'}")
    print("-" * 100)
    for scale_name, tiled, ga, gn, ca, cn in scales:
        entries = []
        if ga is not None:
            entries.append(("MI300A GPU", ga, "MI300A GPU", GPU_VARIANTS))
        if gn is not None:
            entries.append(("GH200 GPU",  gn, "GH200 GPU", GPU_VARIANTS))
        if have_cpu and ca is not None:
            entries.append(("MI300A Zen CPU", ca, "MI300A Zen CPU", CPU_VARIANTS))
        if have_cpu and cn is not None:
            entries.append(("Grace CPU",      cn, "Grace CPU",      CPU_VARIANTS))

        for label, df, pk, vlist in entries:
            peak = STREAM_PEAK.get(pk, 1.0)
            for dk in DIST_ORDER:
                sl = select_dist(df, dk, tiled)
                for vk in vlist:
                    sub = sl[sl["variant"] == vk].copy()
                    if sub.empty: continue
                    sub["bw"] = compute_bandwidth(sub["time_ms"].values,
                                                  sub["nx"].values, vk)
                    meds = sub.groupby(["tpb","coarsen"])["bw"].median()
                    if meds.empty: continue
                    bt, bc = meds.idxmax()
                    m = meds.max()
                    p = 100*m/peak
                    print(f"{scale_name:<8} {label:<22} {dk:<10} {vk:<20} "
                          f"({bt},{bc}){'':>8} {m:.4f}      {p:.1f}%")


if __name__ == "__main__":
    main()