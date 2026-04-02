#!/usr/bin/env python3
"""
plot_addition_paper.py
1x2 violin: Schedule Only | W. Layout Transformation

Usage:
    python plot_addition_paper.py --cpu cpu.csv --gpu gpu.csv
    python plot_addition_paper.py --cpu cpu.csv --gpu gpu.csv --add-peak
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np, argparse, re
from collections import defaultdict

STREAM_PEAK = {
    "MI300A Zen Cores": 1160.35, "Grace CPU": 1700.62,
    "MI300A": 4271,  "GH200":  3720.48,
}

VCOL = {"naive": "#e67e22", "tiled": "#27ae60", "perm": "#2980b9", "blk": "#9b59b6"}

XLABEL_CPU = {
    "naive": "Naive 1D",
    "tiled": "2D Tiled",
    "perm":  "Permuted",
    "blk":   "Blocked",
}
XLABEL_GPU = {
    "naive": "Naive 1D",
    "tiled": "2D Tiled",
    "perm":  "Permuted",
    "blk":   "Blocked",
}

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
            try: rows.append(dict(variant=p[0], tile=int(p[3]), gbs=float(p[7])))
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
            try: rows.append(dict(kernel=kernel.strip(), gbs=float(rest[6])))
            except: continue
    return rows

# ══════════════════════════════════════════════════════════════════════
#  Grouping
# ══════════════════════════════════════════════════════════════════════

def remove_outliers(v, k=3.0):
    if len(v) < 4: return v
    q1, q3 = np.percentile(v, [25, 75]); iqr = q3 - q1
    c = v[(v >= q1-k*iqr) & (v <= q3+k*iqr)]
    return c if len(c) > 2 else v

def best_of(groups):
    if not groups: return None
    bk = max(groups, key=lambda k: np.median(groups[k]))
    return np.array(groups[bk])

def cpu_cats(rows):
    out = {}
    # naive: row_major
    g = [r["gbs"] for r in rows if r["variant"] == "row_major"]
    if g: out["naive"] = np.array(g)
    # tiled: best T
    gs = defaultdict(list)
    for r in rows:
        if r["variant"] == "tiled": gs[r["tile"]].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["tiled"] = v
    # perm: all_rowmajor
    g = [r["gbs"] for r in rows if r["variant"] == "all_rowmajor"]
    if g: out["perm"] = np.array(g)
    # blk: best of any blk_*
    gs = defaultdict(list)
    for r in rows:
        if r["variant"].startswith("blk_"): gs[(r["variant"], r["tile"])].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["blk"] = v
    return out

def gpu_cats(rows):
    out = {}
    # naive: direct BY=1 (1D schedule, worst coalescing)
    gs = defaultdict(list)
    for r in rows:
        k = r["kernel"]
        if k.startswith("direct ") and not k.startswith("direct_T"):
            m = re.search(r'BY=(\d+)', k)
            if m and m.group(1) == "1": gs[k].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["naive"] = v
    # tiled: best direct (any BY, no smem) — best schedule-only
    gs = defaultdict(list)
    for r in rows:
        k = r["kernel"]
        if (k.startswith("direct ") or k.startswith("direct_T")) and \
           not k.startswith("tiled") and not k.startswith("all_"):
            gs[k].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["tiled"] = v
    # perm: all_rowmajor (layout fix — permute B to row-major)
    gs = defaultdict(list)
    for r in rows:
        if r["kernel"].startswith("all_rowmajor"): gs[r["kernel"]].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["perm"] = v
    # blk: tiled+smem (smem as local layout transformation)
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
        medians.append((pos, float(np.median(arr)), vk))
        xlabels.append((xlabels_map or XLABEL_CPU)[vk])
        pos += 1

    # y-axis
    ymax = float(np.max(np.concatenate(data))) if data else 1
    print(f"{title} max: {ymax:.1f} GB/s", end="")
    if add_peak:
        print(f", STREAM peak: {peak:.1f} GB/s", end="")
    print()
    if add_peak and peak and peak > ymax:
        ymax = peak  # extend axis to show STREAM line
    loc = MaxNLocator(nbins=5, min_n_ticks=5)
    ticks = loc.tick_values(0, ymax * 1.15 if not add_peak else ymax)
    ticks = ticks[ticks >= 0]
    if len(ticks) > 7: ticks = ticks[:7]
    top = ticks[-1] * 1.01 if len(ticks) else ymax * 1.2

    # violins
    if data:
        parts = ax.violinplot(data, positions=positions,
                              showmeans=True, showmedians=True,
                              showextrema=False, widths=0.9)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(colors[i])
            body.set_edgecolor("black")
            body.set_alpha(0.75)
        parts["cmeans"].set_color("black")
        parts["cmedians"].set_color("white")

    # separator
    if sep_x is not None:
        ax.axvline(x=sep_x, color="gray", ls="--", lw=1.5, alpha=0.6)
        ax.text(sep_x - 0.3, top * 0.128, "Schedule\nOnly",
                ha="right", va="top", fontsize=9, color="gray",
                fontweight="bold")
        ax.text(sep_x + 0.9, top * 0.128, "With Layout\nTransformations",
                ha="left", va="top", fontsize=9, color="gray",
                fontweight="bold")

    ax.set_xticks(positions)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_yticks(ticks)
    ax.set_ylim(0, top)
    ax.set_title(title, fontsize=13)
    ax.grid(axis="y", alpha=0.25)

    # STREAM
    if peak:
        ax.text(0.03, 0.97, f"STREAM peak {peak:.0f} GB/s",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=9, color="dimgray")
        if add_peak:
            ax.axhline(y=peak, color="dimgray", ls="--", lw=1, alpha=0.5)

    # % annotations
    if peak:
        off = 0.04 * top
        for p, med, vk in medians:
            pct = 100.0 * med / peak
            if pct < 10:
                ax.text(p, med + off, f"{pct:.0f}%",
                        ha="center", va="bottom", fontsize=11,
                        color=VCOL[vk], fontweight="bold")
            else:
                ax.text(p, med - off, f"{pct:.0f}%",
                        ha="center", va="top", fontsize=11,
                        color=VCOL[vk], fontweight="bold")

    # ── Gap arrows ──
    med_dict = {vk: (p, m) for p, m, vk in medians}

    # Arrow 1: tiled → perm (schedule-only gap)
    if "tiled" in med_dict and "perm" in med_dict:
        p_tiled, y_tiled = med_dict["tiled"]
        p_perm,  y_perm  = med_dict["perm"]
        # Place arrow to the right of the last schedule violin
        arrow_x = med_dict["tiled"][0] + 0.82
        ax.plot([p_tiled, arrow_x], [y_tiled, y_tiled],
                color="#555555", ls=":", lw=1, alpha=0.5)
        ax.plot([p_perm, arrow_x], [y_perm, y_perm],
                color="#555555", ls=":", lw=1, alpha=0.5)
        ax.annotate("", xy=(arrow_x, y_perm), xytext=(arrow_x, y_tiled),
                    arrowprops=dict(arrowstyle="<->", color="#555555",
                                   lw=1.5, shrinkA=2, shrinkB=2))
        mid_y = (y_tiled + y_perm) / 2

        if gpu:
            ax.text(arrow_x + 0.08, mid_y - 0.21*ymax,
                    "Gap between\nbest schedule\nand best layout",
                    fontsize=9.25, color="#555555", va="center", ha="left",
                    style="italic")
        else:
            ax.text(arrow_x + 0.08, mid_y - 0.12*ymax,
                    "Gap between\nbest schedule\nand best layout",
                    fontsize=9.25, color="#555555", va="center", ha="left",
                    style="italic")

    # Arrow 2: perm → blk (better layout)
    if "perm" in med_dict and "blk" in med_dict:
        _, y_perm = med_dict["perm"]
        _, y_blk  = med_dict["blk"]
        arrow_x = med_dict["blk"][0] + 0.42
        mid_y = (y_perm + y_blk) / 2

# ══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", default=None)
    ap.add_argument("--gpu", default=None)
    ap.add_argument("--cpu-label", default="MI300A CPU")
    ap.add_argument("--gpu-label", default="MI300A GPU")
    ap.add_argument("--add-peak", action="store_true")
    args = ap.parse_args()

    panels = []
    if args.cpu:
        try:
            cats = cpu_cats(parse_cpu(args.cpu))
            panels.append((args.cpu_label, cats, STREAM_PEAK.get(args.cpu_label, 0), False))
        except FileNotFoundError:
            print(f"[WARN] {args.cpu} not found")
    if args.gpu:
        try:
            cats = gpu_cats(parse_gpu(args.gpu))
            panels.append((args.gpu_label, cats, STREAM_PEAK.get(args.gpu_label, 0), True))
        except FileNotFoundError:
            print(f"[WARN] {args.gpu} not found")

    if not panels:
        print("No data."); return

    ncols = len(panels)
    fig, axes = plt.subplots(1, ncols, figsize=(4.2 * ncols, 3.5), squeeze=False)

    for ci, (title, cats, peak, is_gpu) in enumerate(panels):
        ax = axes[0, ci]
        xmap = XLABEL_GPU if is_gpu else XLABEL_CPU
        assert args.add_peak
        draw_panel(ax, cats, title, peak, args.add_peak, xmap, is_gpu)
        if ci == 0: ax.set_ylabel("Bandwidth [GB/s]", fontsize=12)

    fig.suptitle("Addition ($C[:] += A[:] + B[:]$) with Suboptimal Layouts",
                fontsize=14, y=0.88)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    sfx = "_w_peak" if args.add_peak else ""
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT_STEM}{sfx}.{ext}", dpi=200, bbox_inches="tight")

    # summary
    for title, cats, peak, is_gpu in panels:
        xmap = XLABEL_GPU if is_gpu else XLABEL_CPU
        print(f"\n  {title}:")
        for vk in ["naive", "tiled", "perm", "blk"]:
            if vk not in cats: continue
            med = float(np.median(cats[vk]))
            pk = f"  ({100*med/peak:.0f}%)" if peak else ""
            print(f"    {xmap[vk].replace(chr(10),' '):<24} {med:8.1f} GB/s{pk}")

    # Build filename with platform tags
    tags = []
    if args.cpu: tags.append("cpu")
    if args.gpu: tags.append("gpu")
    sfx = "_" + "_".join(tags)
    if args.add_peak: sfx += "_w_peak"
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT_STEM}{sfx}.{ext}", dpi=200, bbox_inches="tight")

if __name__ == "__main__":
    main()