#!/usr/bin/env python3
"""
plot_addition_paper.py
Up to 2×2 violin grid:  rows = {AMD, NVIDIA},  cols = {CPU, GPU}

Usage:
    # Single platform (1×2)
    python plot_addition_paper.py --amd-cpu cpu.csv --amd-gpu gpu.csv --add-peak

    # Both platforms (2×2)
    python plot_addition_paper.py \
        --amd-cpu amd_cpu.csv  --amd-gpu amd_gpu.csv \
        --nv-cpu  nv_cpu.csv   --nv-gpu  nv_gpu.csv  --add-peak

    # Any subset works (e.g. CPUs only → 2×1)
    python plot_addition_paper.py --amd-cpu amd.csv --nv-cpu nv.csv --add-peak
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter
import numpy as np, argparse, re
from collections import defaultdict
from scipy.stats import bootstrap

# ══════════════════════════════════════════════════════════════════════
#  Constants
# ══════════════════════════════════════════════════════════════════════

STREAM_PEAK = {
    "MI300A Zen CPU": 1161*1e-3,  "GH200 Grace CPU": 1806.62*1e-3,
    "MI300A GPU":       4294*1e-3,     "GH200 Hopper GPU": 3780*1e-3,
}

VCOL = {"naive": "#e67e22", "tiled": "#27ae60", "perm": "#2980b9", "blk": "#9b59b6"}

XLABEL_CPU = {"naive": "Naive 1D", "tiled": "2D Tiled", "perm": "Permuted", "blk": "Blocked"}
XLABEL_GPU = {"naive": "Naive 1D", "tiled": "2D Tiled", "perm": "Permuted", "blk": "Blocked"}

OUT_STEM = "addition_paper"

# ══════════════════════════════════════════════════════════════════════
#  CSV parsing
# ══════════════════════════════════════════════════════════════════════

def parse_cpu(path):
    rows = []
    with open(path) as f:
        for l in f:
            l = l.strip()
            if not l or l.startswith("variant"): continue
            p = l.split(",")
            if len(p) < 10: continue
            try: rows.append(dict(variant=p[0], tile=int(p[3]), gbs=float(p[7])*1e-3))
            except: continue
    return rows

def parse_gpu(path):
    rows = []
    with open(path) as f:
        for l in f:
            l = l.strip()
            if not l or l.startswith("kernel"): continue
            if l.startswith('"'):
                eq = l.index('"', 1); kernel = l[1:eq]
                rest = l[eq+1:].lstrip(",").split(",")
            else:
                parts = l.split(","); kernel = parts[0]; rest = parts[1:]
            if len(rest) < 9: continue
            try: rows.append(dict(kernel=kernel.strip(), gbs=float(rest[6])*1e-3))
            except: continue
    return rows

# ══════════════════════════════════════════════════════════════════════
#  Grouping
# ══════════════════════════════════════════════════════════════════════

def remove_outliers(v, k=3.0):
    if len(v) < 4: return v
    q1, q3 = np.percentile(v, [25, 75]); iqr = q3 - q1
    c = v[(v >= q1 - k * iqr) & (v <= q3 + k * iqr)]
    return c if len(c) > 2 else v

import warnings

def bootstrap_ci(arr, confidence_level=0.95, n_resamples=10000):
    """Bootstrap CI of the median using scipy BCa method, fallback to percentile."""
    if len(arr) < 3:
        med = np.median(arr)
        return med, med
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            res = bootstrap((arr,), np.median,
                            confidence_level=confidence_level,
                            n_resamples=n_resamples, method='BCa')
            lo, hi = res.confidence_interval.low, res.confidence_interval.high
            if np.isnan(lo) or np.isnan(hi):
                raise ValueError("BCa returned nan")
            return lo, hi
        except Exception:
            try:
                res = bootstrap((arr,), np.median,
                                confidence_level=confidence_level,
                                n_resamples=n_resamples, method='percentile')
                lo, hi = res.confidence_interval.low, res.confidence_interval.high
                if np.isnan(lo) or np.isnan(hi):
                    raise ValueError("percentile returned nan")
                return lo, hi
            except Exception:
                med = np.median(arr)
                return med, med

def best_of(groups):
    if not groups: return None
    bk = max(groups, key=lambda k: np.median(groups[k]))
    return np.array(groups[bk])

def cpu_cats(rows):
    out = {}
    g = [r["gbs"] for r in rows if r["variant"] == "row_major"]
    if g: out["naive"] = np.array(g)
    gs = defaultdict(list)
    for r in rows:
        if r["variant"] == "tiled": gs[r["tile"]].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["tiled"] = v
    g = [r["gbs"] for r in rows if r["variant"] == "all_rowmajor"]
    if g: out["perm"] = np.array(g)
    gs = defaultdict(list)
    for r in rows:
        if r["variant"].startswith("blk_"): gs[(r["variant"], r["tile"])].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["blk"] = v
    return out

def gpu_cats(rows):
    out = {}
    gs = defaultdict(list)
    for r in rows:
        k = r["kernel"]
        if k.startswith("direct ") and not k.startswith("direct_T"):
            m = re.search(r'BY=(\d+)', k)
            if m and m.group(1) == "1": gs[k].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["naive"] = v
    gs = defaultdict(list)
    for r in rows:
        k = r["kernel"]
        if (k.startswith("direct ") or k.startswith("direct_T")) and \
           not k.startswith("tiled") and not k.startswith("all_"):
            gs[k].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["tiled"] = v
    gs = defaultdict(list)
    for r in rows:
        if r["kernel"].startswith("all_rowmajor"): gs[r["kernel"]].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["perm"] = v
    gs = defaultdict(list)
    for r in rows:
        k = r["kernel"]
        if k.startswith("tiled+smem") or k.startswith("tiled_smem"):
            gs[k].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["blk"] = v
    return out

# ══════════════════════════════════════════════════════════════════════
#  Panel drawing
# ══════════════════════════════════════════════════════════════════════

def draw_panel(ax, cats, title, peak=None, add_peak=False, xlabels_map=None, gpu=False):
    order = ["naive", "tiled", "perm", "blk"]
    present = [k for k in order if k in cats]
    if not present:
        ax.set_title(title); return

    positions, data, colors, xlabels = [], [], [], []
    medians = []
    pos = 0
    sep_x = None

    for vk in present:
        if vk == "perm" and sep_x is None and pos > 0:
            sep_x = pos - 0.5
            pos += 0.4
        arr = remove_outliers(cats[vk])
        if len(arr) == 0: pos += 1; continue
        data.append(arr)
        positions.append(pos)
        colors.append(VCOL[vk])
        medians.append((pos, float(np.median(arr)), vk, float(np.min(arr))))
        xlabels.append((xlabels_map or XLABEL_CPU)[vk])
        pos += 1

    # y-axis
    ymax = float(np.max(np.concatenate(data))) if data else 1
    if add_peak and peak and peak > ymax:
        ymax = peak
    if ymax > 3.5 and ymax < 4.9:
        ymax += 0.15
    ax.set_ylim(0, ymax)
    loc = MaxNLocator(nbins=5, min_n_ticks=5)
    ticks = loc.tick_values(0, ymax)
    ticks = ticks[ticks >= 0]
    if len(ticks) > 7: ticks = ticks[:7]
    top = ticks[-1] * 1.01 if len(ticks) else ymax

    # violins
    if data:
        parts = ax.violinplot(data, positions=positions,
                              showmeans=True, showmedians=True,
                              showextrema=False, widths=0.9)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(colors[i])
            body.set_edgecolor("black")
            body.set_alpha(0.75)
            body.set_zorder(2)
        parts["cmeans"].set_color("black")
        parts["cmeans"].set_zorder(3)
        parts["cmedians"].set_color("white")
        parts["cmedians"].set_zorder(3)

        # Bootstrap 95% CI of median
        ci_color = "#FF69B4"  # soft hot pink
        for i, arr in enumerate(data):
            ci_lo, ci_hi = bootstrap_ci(arr)
            print(f"  CI [{title}] pos={positions[i]}: lo={ci_lo:.4f} hi={ci_hi:.4f} span={ci_hi-ci_lo:.6f} top={top:.4f}")
            # Ensure minimum visual height so CI is always visible
            ci_mid = (ci_lo + ci_hi) / 2
            min_height = 0.01 * top if top > 0 else 0.001
            if (ci_hi - ci_lo) < min_height:
                ci_lo = ci_mid - min_height / 2
                ci_hi = ci_mid + min_height / 2
            ax.vlines(positions[i], ci_lo, ci_hi,
                      color="black", lw=3.5, zorder=10)
            ax.vlines(positions[i], ci_lo, ci_hi,
                      color=ci_color, lw=1.8, zorder=11)
            cap_w = 0.08
            for y in (ci_lo, ci_hi):
                ax.hlines(y, positions[i] - cap_w, positions[i] + cap_w,
                          color="black", lw=3, zorder=10)
                ax.hlines(y, positions[i] - cap_w, positions[i] + cap_w,
                          color=ci_color, lw=1.5, zorder=11)

    # separator
    if sep_x is not None:
        ax.axvline(x=sep_x, color="gray", ls="--", lw=1.5, alpha=0.6)
        ax.text(sep_x - 0.1, top * 0.158, "Schedule\nOnly",
                ha="right", va="top", fontsize=9, color="gray", fontweight="bold")
        ax.text(sep_x + 0.1, top * 0.158, "With Layout\nTransformations",
                ha="left", va="top", fontsize=9, color="gray", fontweight="bold")

    ax.set_xticks(positions)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.tick_params(axis='y', which='minor', length=3)
    ax.set_ylim(0, top)
    ax.set_title(title, fontsize=11)
    ax.grid(axis="y", alpha=0.25, zorder=0)
    ax.grid(axis="y", which='minor', alpha=0.12, ls=':', zorder=0)
    ax.set_axisbelow(True)  # force grid behind all data

    # STREAM peak
    if peak:
        ax.text(0.03, 0.97, f"STREAM {peak:.2f} TB/s",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=9, color="dimgray", zorder=1)
        ax.axhline(y=peak, color="dimgray", ls="--", lw=1, alpha=0.5, zorder=1)

    # % annotations
    if peak:
        off = 0.04 * top
        for p, med, vk, vmin in medians:
            pct = 100.0 * med / peak
            if pct < 10:
                ax.text(p, med + off, f"{pct:.0f}%",
                        ha="center", va="bottom", fontsize=10,
                        color=VCOL[vk], fontweight="bold")
            else:
                ax.text(p, vmin - off, f"{pct:.0f}%",
                        ha="center", va="top", fontsize=10,
                        color=VCOL[vk], fontweight="bold")

        if peak > 3000:
            if top < peak*1.1:
                ax.set_ylim(0, peak * 1.1)

    # ── Gap arrow: best schedule → best layout ──
    med_dict = {vk: (p, m) for p, m, vk, _ in medians}
    min_dict = {vk: vmin for _, _, vk, vmin in medians}

    # Pick target: if blk >= 10% of peak above perm, arrow goes to blk
    layout_target = None
    if "perm" in med_dict:
        layout_target = "perm"
        if "blk" in med_dict and peak:
            _, y_p = med_dict["perm"]
            _, y_b = med_dict["blk"]
            if (y_b - y_p) >= 0.10 * peak:
                layout_target = "blk"

    # Arrow 1: tiled → best layout target
    if "tiled" in med_dict and layout_target and layout_target in med_dict:
        p_tiled, y_tiled = med_dict["tiled"]
        p_target, y_target = med_dict[layout_target]
        arrow_x = p_tiled + 0.82
        ax.plot([p_tiled, arrow_x], [y_tiled, y_tiled],
                color="#555555", ls=":", lw=1.5, alpha=0.5)
        ax.plot([p_target, arrow_x], [y_target, y_target],
                color="#555555", ls=":", lw=1.5, alpha=0.5)
        style = "|-|, widthA=0.15, widthB=0.15" if abs(y_target - y_tiled) < 0.05 * ymax else "<->"
        ax.annotate("", xy=(arrow_x, y_target), xytext=(arrow_x, y_tiled),
                    arrowprops=dict(arrowstyle=style, color="#555555",
                                    lw=1.5, shrinkA=0, shrinkB=0))
        mid_y = (y_tiled + y_target) / 2
        y_off = 0.21 * ymax if gpu else 0.12 * ymax
        text_y = mid_y - y_off
        perm_min = min_dict.get("perm", text_y)
        print(perm_min, text_y)
        if perm_min < text_y * 1.3:
            text_y = perm_min - 0.25 * ymax
        if y_target - y_tiled > 0.15 * peak:
            ax.text(arrow_x + 0.08, text_y - 0.05,
                    "Gap between\nbest schedule\nand best layout",
                    fontsize=9, color="#555555", va="center", ha="left",
                    style="italic")

    # Arrow 2: perm → blk, only if blk wasn't already Arrow 1 target
    if "perm" in med_dict and "blk" in med_dict and layout_target != "blk":
        p_perm, y_perm = med_dict["perm"]
        p_blk,  y_blk  = med_dict["blk"]
        arrow_x = p_blk + 0.42
        mid_y = (y_perm + y_blk) / 2


# ══════════════════════════════════════════════════════════════════════
#  Grid assembly
# ══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    # AMD
    ap.add_argument("--amd-cpu", default="results/beverin/madd_beverin_cpu.csv", help="CSV for AMD CPU")
    ap.add_argument("--amd-gpu", default="results/beverin/madd_beverin_gpu.csv", help="CSV for AMD GPU")
    ap.add_argument("--amd-cpu-label", default="MI300A Zen CPU")
    ap.add_argument("--amd-gpu-label", default="MI300A GPU")
    # NVIDIA
    ap.add_argument("--nv-cpu",  default="results/daint/madd_daint_cpu.csv", help="CSV for NVIDIA CPU")
    ap.add_argument("--nv-gpu",  default="results/daint/madd_daint_gpu.csv", help="CSV for NVIDIA GPU")
    ap.add_argument("--nv-cpu-label", default="GH200 Grace CPU")
    ap.add_argument("--nv-gpu-label", default="GH200 Hopper GPU")
    # Legacy single-platform (maps to AMD slots)
    ap.add_argument("--cpu", default=None, help="(legacy) same as --amd-cpu")
    ap.add_argument("--gpu", default=None, help="(legacy) same as --amd-gpu")
    ap.add_argument("--cpu-label", default=None)
    ap.add_argument("--gpu-label", default=None)
    # BW correction factor (e.g. 0.75 to fix 4x→3x)
    ap.add_argument("--bw-scale", type=float, default=None,
                    help="Multiply all BW values by this factor (e.g. 0.75)")
    ap.add_argument("--add-peak", action="store_true")
    args = ap.parse_args()

    # Legacy fallback
    if args.cpu and not args.amd_cpu:
        args.amd_cpu = args.cpu
        if args.cpu_label: args.amd_cpu_label = args.cpu_label
    if args.gpu and not args.amd_gpu:
        args.amd_gpu = args.gpu
        if args.gpu_label: args.amd_gpu_label = args.gpu_label

    # ── Build the panel grid ──
    grid = {}

    def try_add(csv_path, parser, grouper, row, col, label, is_gpu):
        if csv_path is None: return
        try:
            cats = grouper(parser(csv_path))
            if args.bw_scale:
                cats = {k: v * args.bw_scale for k, v in cats.items()}
            peak = STREAM_PEAK.get(label, 0)
            grid[(row, col)] = (label, cats, peak, is_gpu)
        except FileNotFoundError:
            print(f"[WARN] {csv_path} not found")

    try_add(args.amd_cpu, parse_cpu, cpu_cats, "amd", "cpu", args.amd_cpu_label, False)
    try_add(args.amd_gpu, parse_gpu, gpu_cats, "amd", "gpu", args.amd_gpu_label, True)
    try_add(args.nv_cpu,  parse_cpu, cpu_cats, "nv",  "cpu", args.nv_cpu_label,  False)
    try_add(args.nv_gpu,  parse_gpu, gpu_cats, "nv",  "gpu", args.nv_gpu_label,  True)

    if not grid:
        print("No data."); return

    # ── Determine active rows / cols ──
    active_rows = sorted({r for r, c in grid}, key=["amd", "nv"].index)
    active_cols = sorted({c for r, c in grid}, key=["cpu", "gpu"].index)
    nrows = len(active_rows)
    ncols = len(active_cols)

    print(f"Grid: {nrows}x{ncols}  rows={active_rows}  cols={active_cols}")

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.6 * ncols, 2.8 * nrows),
                             squeeze=False)

    for ri, row_key in enumerate(active_rows):
        for ci, col_key in enumerate(active_cols):
            ax = axes[ri, ci]
            ax.set_box_aspect(2.8 / 4.8)
            if (row_key, col_key) not in grid:
                ax.set_visible(False)
                continue
            title, cats, peak, is_gpu = grid[(row_key, col_key)]
            xmap = XLABEL_GPU if is_gpu else XLABEL_CPU
            draw_panel(ax, cats, title, peak, args.add_peak, xmap, is_gpu)
            if ci == 0:
                ax.set_ylabel("Bandwidth [TB/s]", fontsize=11)

    fig.suptitle("Matrix Addition (C += A + B) with Suboptimal Layouts",
                 fontsize=15, y=0.89 if nrows > 1 else 0.98)
    fig.text(0.5, 0.85 if nrows > 1 else 0.95,
            "% annotations relative to STREAM peak bandwidth",
            ha='center', va='top', fontsize=12, color='dimgray')
    fig.tight_layout(rect=[0, 0, 1, 0.89 if nrows > 1 else 0.92])

    tags = sorted({c for _, c in grid})
    if any(r == "cpu" for r, _ in grid): tags.append("cpu")
    if any(r == "gpu" for r, _ in grid): tags.append("gpu")
    sfx = "_" + "_".join(sorted(set(tags)))
    if args.add_peak: sfx += "_w_peak"
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT_STEM}{sfx}.{ext}", dpi=200, bbox_inches="tight")
    print(f"Saved {OUT_STEM}{sfx}.pdf/png")

    for (row_key, col_key) in sorted(grid.keys()):
        title, cats, peak, is_gpu = grid[(row_key, col_key)]
        xmap = XLABEL_GPU if is_gpu else XLABEL_CPU
        print(f"\n  {title}:")
        for vk in ["naive", "tiled", "perm", "blk"]:
            if vk not in cats: continue
            med = float(np.median(cats[vk]))
            pk = f"  ({100 * med / peak:.2f}%)" if peak else ""
            print(f"    {xmap[vk]:<24} {med:8.2f} TB/s{pk}")

if __name__ == "__main__":
    main()