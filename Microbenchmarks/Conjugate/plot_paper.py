#!/usr/bin/env python3
"""
plot_conjugate.py – Line+scatter plots for the conjugate microbenchmark.

2×2 grid:  rows = {AMD MI300A, NVIDIA GH200},  cols = {CPU, GPU}
X-axis:    P (number of complex pairs) = 3, 6, 9, 12, 15, 18, 21
Lines:     AoS, SoA, AoSoA-8, AoSoA-16, AoSoA-32
Markers:   median ± 95% CI,  shaded band = CI

Generates two figures: one for in-place, one for out-of-place.

Usage:
    python plot_conjugate.py                         # defaults
    python plot_conjugate.py --bev-dir results/beverin --daint-dir results/daint
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
import numpy as np
import argparse, os, sys
from collections import defaultdict

# ══════════════════════════════════════════════════════════════════════
#  Constants
# ══════════════════════════════════════════════════════════════════════

STREAM_PEAK = {
    "MI300A Zen CPU":  1160.35,   "GH200 Grace CPU":  1946.62,
    "MI300A GPU":  4294,      "GH200 Hopper GPU":  3780,
}

LAYOUTS = ["AoS", "SoA", "AoSoA-16"] # , "AoSoA-16", "AoSoA-32"

STYLE = {
    "AoS":      dict(color="#e74c3c", marker="o",  ls="-",  lw=1.0),
    "SoA":      dict(color="#2ecc71", marker="s",  ls="-",  lw=1.0),
    #"AoSoA-8":  dict(color="#3498db", marker="^",  ls="--", lw=1.6),
    "AoSoA-16": dict(color="#9b59b6", marker="D",  ls="--", lw=1.0),
    #"AoSoA-32": dict(color="#e67e22", marker="v",  ls="--", lw=1.6),
}

# ══════════════════════════════════════════════════════════════════════
#  CSV parsing
# ══════════════════════════════════════════════════════════════════════

def parse_csv(path):
    """Return {(P, layout): [gbps, …]}."""
    groups = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("P,"): continue
            parts = line.split(",")
            if len(parts) < 5: continue
            try:
                P      = int(parts[0])
                layout = parts[1].strip()
                gbps   = float(parts[4])
                groups[(P, layout)].append(gbps)
            except (ValueError, IndexError):
                continue
    return groups

# ══════════════════════════════════════════════════════════════════════
#  Statistics
# ══════════════════════════════════════════════════════════════════════

def stats(vals):
    """Return (median, p5, p95) — 5th/95th percentile spread of the data."""
    a = np.array(vals)
    if len(a) < 3:
        m = np.median(a)
        return m, m, m
    med = np.median(a)
    p5  = np.percentile(a, 5)
    p95 = np.percentile(a, 95)
    return med, p5, p95

# ══════════════════════════════════════════════════════════════════════
#  Panel drawing
# ══════════════════════════════════════════════════════════════════════

def draw_panel(ax, groups, title, peak=None, add_peak=True):
    """Draw line+scatter+CI for one subplot."""

    ps_all = sorted({p for (p, _) in groups})
    if not ps_all:
        ax.set_title(title, fontsize=11)
        return

    ymax = 0

    for lay in LAYOUTS:
        ps, meds, lo, hi = [], [], [], []
        for p in ps_all:
            key = (p, lay)
            if key not in groups: continue
            m, cl, ch = stats(groups[key])
            ps.append(p)
            meds.append(m)
            lo.append(cl)
            hi.append(ch)
        if not ps: continue

        ps   = np.array(ps)
        meds = np.array(meds)
        lo   = np.array(lo)
        hi   = np.array(hi)
        ymax = max(ymax, float(np.max(hi)))

        sty = STYLE[lay]

        # Shaded P5–P95 band
        ax.fill_between(ps, lo, hi, color=sty["color"], alpha=0.22, edgecolor="none")

        # Line (median)
        ax.plot(ps, meds,
                color=sty["color"], ls=sty["ls"], lw=sty["lw"],
                marker=sty["marker"], markersize=6, markeredgecolor="white",
                markeredgewidth=0.8, label=lay, zorder=3)

        # Error bars (P5–P95)
        yerr = np.array([meds - lo, hi - meds])
        ax.errorbar(ps, meds, yerr=yerr, fmt="none",
                    ecolor=sty["color"], elinewidth=1.4, capsize=4,
                    capthick=1.2, alpha=0.8, zorder=2)

    # STREAM peak line
    if add_peak and peak:
        ax.axhline(y=peak, color="dimgray", ls="--", lw=1, alpha=0.5, zorder=1)
        ax.text(ps_all[0], peak * 1.02,
                f"STREAM {peak:.0f} GB/s",
                fontsize=7.5, color="dimgray", va="bottom")
        if peak > ymax:
            ymax = peak

    # Axes — force exactly 5 major y-ticks
    ax.set_xticks(ps_all)
    ax.set_xlim(ps_all[0] - 1, ps_all[-1] + 1)
    top = ymax * 1.12
    ax.set_ylim(0, top)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=5, integer=True))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.tick_params(axis="y", which="minor", length=3)
    ax.grid(axis="y", alpha=0.25)
    ax.grid(axis="y", which="minor", alpha=0.12, ls=":")
    ax.grid(axis="x", alpha=0.12, ls=":")
    ax.set_title(title, fontsize=11)

# ══════════════════════════════════════════════════════════════════════
#  Grid assembly  (one figure per mode: IP / OOP)
# ══════════════════════════════════════════════════════════════════════

def build_figure(grid, mode_label, out_stem):
    """
    grid: dict  (row_key, col_key) → (title, groups, peak)
    mode_label: "In-Place" or "Out-of-Place"
    """
    if not grid:
        print(f"  [{mode_label}] no data"); return

    active_rows = sorted({r for r, _ in grid}, key=lambda x: ["amd", "nv"].index(x))
    active_cols = sorted({c for _, c in grid}, key=lambda x: ["cpu", "gpu"].index(x))
    nrows = len(active_rows)
    ncols = len(active_cols)

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.6 * ncols, 2.8 * nrows),
                             squeeze=False)

    for ri, rk in enumerate(active_rows):
        for ci, ck in enumerate(active_cols):
            ax = axes[ri, ci]
            if (rk, ck) not in grid:
                ax.set_visible(False); continue
            title, groups, peak = grid[(rk, ck)]
            draw_panel(ax, groups, title, peak)
            if ci == 0:
                ax.set_ylabel("Bandwidth [GB/s]", fontsize=10)
            if ri == nrows - 1:
                ax.set_xlabel("Number of complex arrays (P)", fontsize=10)

    # Legend from bottom-right panel (most likely to have all layouts)
    handles, labels = axes[-1, -1].get_legend_handles_labels()
    if not handles:
        for ax in axes.flat:
            handles, labels = ax.get_legend_handles_labels()
            if handles: break
    if handles:
        fig.legend(handles, labels,
                   loc="lower center", ncol=len(labels),
                   fontsize=9, frameon=True, fancybox=True,
                   bbox_to_anchor=(0.5, -0.01))

    titles = {
        "In-Place":      r"In-Place Conjugate: $A = \bar{A}$ for P Complex Arrays",
        "Out-of-Place":  r"Out-of-Place Conjugate: $B = \bar{A}$ P Complex Arrays",
    }
    fig.suptitle(
        f"{titles[mode_label]}"
        if mode_label == "In-Place" else
        f"{titles[mode_label]}",
        fontsize=14, y=0.89 if nrows > 1 else 0.99)
    #fig.text(0.5, 0.85 if nrows > 1 else 0.95,
    #         "Shaded band / error bars = P5–P95 spread;  dashed line = STREAM peak",
    #         ha="center", va="top", fontsize=10, color="dimgray")
    fig.tight_layout(rect=[0, 0.03, 1, 0.90 if nrows > 1 else 0.92])

    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=200, bbox_inches="tight")
    print(f"  Saved {out_stem}.pdf / .png")
    plt.close(fig)

# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bev-dir",   default="results/beverin",
                    help="Directory with Beverin CSVs")
    ap.add_argument("--daint-dir", default="results/daint",
                    help="Directory with Daint CSVs")
    ap.add_argument("--out-prefix", default="conjugate", help="Output filename prefix")
    args = ap.parse_args()

    # ── File mapping ──
    # Each entry: (row_key, col_key, label, peak_key, filename)
    slots = [
        ("amd", "cpu", "MI300A Zen CPU",  "MI300A Zen CPU",  args.bev_dir,   "results_cpu"),
        ("amd", "gpu", "MI300A GPU",  "MI300A GPU",  args.bev_dir,   "results_gpu"),
        ("nv",  "cpu", "GH200 Grace CPU",   "GH200 Grace CPU",   args.daint_dir, "results_cpu"),
        ("nv",  "gpu", "GH200 Hopper GPU",   "GH200 Hopper GPU",  args.daint_dir, "results_gpu"),
    ]

    for mode, mode_label in [("oop", "Out-of-Place"), ("inplace", "In-Place")]:
        grid = {}
        for rk, ck, label, pk, dirpath, stem in slots:
            fname = f"{stem}_{mode}.csv"
            path  = os.path.join(dirpath, fname)
            if not os.path.isfile(path):
                # try alternate name
                path2 = os.path.join(dirpath, f"{stem}{mode}.csv")
                if os.path.isfile(path2):
                    path = path2
                else:
                    continue
            groups = parse_csv(path)
            if not groups: continue
            peak = STREAM_PEAK.get(pk, None)
            grid[(rk, ck)] = (label, groups, peak)

        # ── Sanity check: no data should exceed STREAM peak ──
        for (rk, ck), (label, groups, peak) in grid.items():
            if peak is None:
                continue
            all_vals = np.concatenate([np.array(v) for v in groups.values()])
            data_max = float(np.max(all_vals))
            if data_max > peak:
                raise ValueError(
                    f"[{mode_label}] {label}: measured BW ({data_max:.1f} GB/s) "
                    f"exceeds STREAM peak ({peak:.1f} GB/s). "
                    f"Check your BW formula or update STREAM_PEAK."
                )

        build_figure(grid, mode_label, f"{args.out_prefix}_{mode}")

    # ── Print summary table ──
    for mode in ["oop", "inplace"]:
        print(f"\n{'='*60}")
        print(f"  {mode.upper()} summary (medians)")
        print(f"{'='*60}")
        for rk, ck, label, pk, dirpath, stem in slots:
            fname = f"{stem}_{mode}.csv"
            path  = os.path.join(dirpath, fname)
            if not os.path.isfile(path): continue
            groups = parse_csv(path)
            if not groups: continue
            peak = STREAM_PEAK.get(pk, 0)
            print(f"\n  {label}:")
            ps = sorted({p for p, _ in groups})
            header = f"    {'P':>3}  " + "  ".join(f"{l:>10}" for l in LAYOUTS)
            print(header)
            for p in ps:
                vals = []
                for lay in LAYOUTS:
                    key = (p, lay)
                    if key in groups:
                        m, _, _ = stats(groups[key])
                        pct = f"({100*m/peak:.0f}%)" if peak else ""
                        vals.append(f"{m:7.1f}{pct:>5}")
                    else:
                        vals.append(f"{'—':>12}")
                print(f"    {p:3d}  " + "  ".join(vals))

if __name__ == "__main__":
    main()