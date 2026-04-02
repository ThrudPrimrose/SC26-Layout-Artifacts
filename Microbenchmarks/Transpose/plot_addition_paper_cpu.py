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
    "MI300A CPU": 1160.35, "Grace CPU": 1700.62,
    "MI300A GPU": 3457.5,  "H200 GPU":  3720.48,
}

VCOL = {"naive": "#e67e22", "tiled": "#27ae60", "perm": "#27ae60", "blk": "#9b59b6"}

XLABEL = {
    "naive": "Naive 1D\nSchedule",
    "tiled": "Tiled\nSched.",
    "perm":  "Permuted\nAll Row-Major",
    "blk":   "Blocked\nLayout",
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
    # naive: direct BY=1
    gs = defaultdict(list)
    for r in rows:
        k = r["kernel"]
        if k.startswith("direct ") and not k.startswith("direct_T"):
            m = re.search(r'BY=(\d+)', k)
            if m and m.group(1) == "1": gs[k].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["naive"] = v
    # tiled+smem
    gs = defaultdict(list)
    for r in rows:
        if r["kernel"].startswith("tiled+smem") or r["kernel"].startswith("tiled_smem"):
            gs[r["kernel"]].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["tiled"] = v
    # perm: all_rowmajor
    gs = defaultdict(list)
    for r in rows:
        if r["kernel"].startswith("all_rowmajor"): gs[r["kernel"]].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["perm"] = v
    # blk
    gs = defaultdict(list)
    for r in rows:
        if r["kernel"].startswith("blk_") or r["kernel"].startswith("blocked"):
            gs[r["kernel"]].append(r["gbs"])
    v = best_of(gs)
    if v is not None: out["blk"] = v
    return out

# ══════════════════════════════════════════════════════════════════════
#  Panel drawing
# ══════════════════════════════════════════════════════════════════════

def draw_panel(ax, cats, title, peak=None, add_peak=False):
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
        xlabels.append(XLABEL[vk])
        pos += 1

    # y-axis
    ymax = float(np.max(np.concatenate(data))) if data else 1
    loc = MaxNLocator(nbins=5, min_n_ticks=5)
    ticks = loc.tick_values(0, ymax * 1.15)
    ticks = ticks[ticks >= 0]
    if len(ticks) > 7: ticks = ticks[:7]
    top = ticks[-1] * 1.06 if len(ticks) else ymax * 1.2

    # violins
    if data:
        parts = ax.violinplot(data, positions=positions,
                              showmeans=True, showmedians=True,
                              showextrema=False, widths=0.7)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(colors[i])
            body.set_edgecolor("black")
            body.set_alpha(0.75)
        parts["cmeans"].set_color("black")
        parts["cmedians"].set_color("white")

    # separator
    if sep_x is not None:
        ax.axvline(x=sep_x, color="gray", ls="--", lw=1.5, alpha=0.6)
        ax.text(sep_x - 0.12, top * 0.97, "Schedule\nOnly",
                ha="right", va="top", fontsize=9, color="gray",
                fontweight="bold")
        ax.text(sep_x + 0.12, top * 0.97, "W. Layout\nTransformation",
                ha="left", va="top", fontsize=9, color="gray",
                fontweight="bold")

    ax.set_xticks(positions)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_yticks(ticks)
    ax.set_ylim(0, top)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.grid(axis="y", alpha=0.25)

    # STREAM
    if peak:
        ax.text(0.03, 0.97, f"STREAM {peak:.0f} GB/s",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=9, color="dimgray")
        if add_peak and peak <= top * 1.1:
            ax.axhline(y=peak, color="red", ls="--", lw=1.2, alpha=0.5)

    # % annotations
    if peak:
        off = 0.04 * top
        for p, med, vk in medians:
            pct = 100.0 * med / peak
            ax.text(p, med - off, f"{pct:.0f}%",
                    ha="center", va="top", fontsize=11,
                    color=VCOL[vk], fontweight="bold")

    # ── Gap arrows ──
    med_dict = {vk: (p, m) for p, m, vk in medians}

    # Arrow 1: tiled → perm (schedule-only gap)
    if "tiled" in med_dict and "perm" in med_dict:
        _, y_tiled = med_dict["tiled"]
        _, y_perm  = med_dict["perm"]
        # Place arrow to the right of the last schedule violin
        arrow_x = med_dict["tiled"][0] + 0.42
        ax.annotate("", xy=(arrow_x, y_perm), xytext=(arrow_x, y_tiled),
                    arrowprops=dict(arrowstyle="<->", color="#555555",
                                   lw=1.5, shrinkA=2, shrinkB=2))
        mid_y = (y_tiled + y_perm) / 2
        ax.text(arrow_x + 0.08, mid_y,
                "Gap: schedule\nonly vs. default\nlayout",
                fontsize=7.5, color="#555555", va="center", ha="left",
                style="italic")

    # Arrow 2: perm → blk (better layout)
    if "perm" in med_dict and "blk" in med_dict:
        _, y_perm = med_dict["perm"]
        _, y_blk  = med_dict["blk"]
        arrow_x = med_dict["blk"][0] + 0.42
        ax.annotate("", xy=(arrow_x, y_blk), xytext=(arrow_x, y_perm),
                    arrowprops=dict(arrowstyle="<->", color="#8e44ad",
                                   lw=1.5, shrinkA=2, shrinkB=2))
        mid_y = (y_perm + y_blk) / 2
        ax.text(arrow_x + 0.08, mid_y,
                "Better layout\nfor CPUs",
                fontsize=7.5, color="#8e44ad", va="center", ha="left",
                fontweight="bold")

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
            panels.append((args.cpu_label, cats, STREAM_PEAK.get(args.cpu_label, 0)))
        except FileNotFoundError:
            print(f"[WARN] {args.cpu} not found")
    if args.gpu:
        try:
            cats = gpu_cats(parse_gpu(args.gpu))
            panels.append((args.gpu_label, cats, STREAM_PEAK.get(args.gpu_label, 0)))
        except FileNotFoundError:
            print(f"[WARN] {args.gpu} not found")

    if not panels:
        print("No data."); return

    ncols = len(panels)
    fig, axes = plt.subplots(1, ncols, figsize=(5.5 * ncols, 4.5), squeeze=False)

    for ci, (title, cats, peak) in enumerate(panels):
        ax = axes[0, ci]
        draw_panel(ax, cats, title, peak, args.add_peak)
        if ci == 0: ax.set_ylabel("Bandwidth [GB/s]", fontsize=12)

    fig.suptitle("$C$ += $A + B$:  A,C row-major, B col-major  ($M = N = 16384$, fp64)",
                 fontsize=12, y=0.99)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    sfx = "_w_peak" if args.add_peak else ""
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT_STEM}{sfx}.{ext}", dpi=200, bbox_inches="tight")

    # summary
    for title, cats, peak in panels:
        print(f"\n  {title}:")
        for vk in ["naive", "tiled", "perm", "blk"]:
            if vk not in cats: continue
            med = float(np.median(cats[vk]))
            pk = f"  ({100*med/peak:.0f}%)" if peak else ""
            print(f"    {XLABEL[vk].replace(chr(10),' '):<24} {med:8.1f} GB/s{pk}")

    print(f"\nSaved: {OUT_STEM}{sfx}.pdf, {OUT_STEM}{sfx}.png")

if __name__ == "__main__":
    main()