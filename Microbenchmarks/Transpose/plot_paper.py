#!/usr/bin/env python3
"""
plot_transpose_paper.py
2×2 violin grid:  rows = {AMD, NVIDIA},  cols = {CPU, GPU}

Four violins per panel:
  1. Library (row-major)     — HPTT on CPU, cuTENSOR/hipTensor on GPU
  2. Ours (row-major) — best non-library row-major variant
  3. Library (blocked)       — HPTT(blk) on CPU, cuTENSOR/hipTensor(blk) on GPU
  4. Ours (blocked)   — best non-library blocked variant

Usage:
    python plot_transpose_paper.py --add-peak
    python plot_transpose_paper.py --add-peak --amd-cpu X.csv --amd-gpu Y.csv ...
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
import numpy as np, argparse, csv as csvmod
from collections import defaultdict

# ══════════════════════════════════════════════════════════════════════
#  Constants
# ══════════════════════════════════════════════════════════════════════

STREAM_PEAK = {
    "MI300A Zen CPU":  1160.35*1e-3,   "GH200 Grace CPU":  1806.62*1e-3,
    "MI300A GPU":  4294*1e-3,      "GH200 Hopper GPU":  3780*1e-3,
}
VCOL = {
    "lib_rm":  "#e67e22",   # orange
    "kern_rm": "#27ae60",   # green
    "lib_blk": "#2980b9",   # blue
    "kern_blk":"#9b59b6",   # purple
}

OUT_STEM = "transpose_paper"

# ── Library name sets ──
CPU_LIB = {"hptt", "hptt_patient", "hptt_blk", "hptt_patient_blk"}
GPU_LIB = {"cutensor", "cutensor_blk", "hiptensor", "hiptensor_blk"}

# ══════════════════════════════════════════════════════════════════════
#  Helpers
# ══════════════════════════════════════════════════════════════════════

def is_blocked_cpu(v):
    """Blocked-layout variant on CPU. 'rm_blk*' is row-major schedule, not blocked layout."""
    return ("blk" in v) and not v.startswith("rm_blk")

def is_blocked_gpu(v):
    """Blocked-layout variant on GPU."""
    return ("blk" in v) or (v == "blocked")

def remove_outliers(v, k=3.0):
    if len(v) < 4: return v
    q1, q3 = np.percentile(v, [25, 75]); iqr = q3 - q1
    c = v[(v >= q1 - k * iqr) & (v <= q3 + k * iqr)]
    return c if len(c) > 2 else v

def best_of(groups):
    """Return (best_key, GBS_array) of the config with highest median."""
    if not groups: return None, None
    bk = max(groups, key=lambda k: np.median(groups[k]))
    return bk, np.array(groups[bk])

# ══════════════════════════════════════════════════════════════════════
#  CSV parsing  (raw per-iteration CSVs)
# ══════════════════════════════════════════════════════════════════════

def parse_cpu_raw(path):
    """CPU raw CSV: variant,N,TB,SB,MT,threads,rep,time_s,gbs,cksum,status"""
    groups = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.strip().split(",")
            if len(p) < 11 or p[0] == "variant":
                continue
            try:
                key = (p[0], p[1], p[2], p[3], p[4], p[5])
                groups[key].append(float(p[8])*1e-3)
            except (ValueError, IndexError):
                continue
    return groups

def parse_gpu_raw(path):
    """GPU raw CSV: variant,N,BX,BY,TX,TY,SB,PAD,rep,time_s,gbs,cksum"""
    groups = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.strip().split(",")
            if len(p) < 12 or p[0] == "variant":
                continue
            try:
                key = tuple(p[:8])
                groups[key].append(float(p[10])*1e-3)
            except (ValueError, IndexError):
                continue
    return groups

# ══════════════════════════════════════════════════════════════════════
#  Categorisation — now returns best key alongside data
# ══════════════════════════════════════════════════════════════════════

def categorise_cpu(groups, max_sb=None):
    """Split CPU groups into {lib_rm, kern_rm, lib_blk, kern_blk}.
    Returns (cats_data, cats_bestkey).  If max_sb is set, blocked non-library
    configs with SB > max_sb are excluded."""
    buckets = {c: {} for c in ("lib_rm", "kern_rm", "lib_blk", "kern_blk")}
    for key, gbs_list in groups.items():
        variant = key[0]
        is_lib = variant in CPU_LIB
        is_blk = is_blocked_cpu(variant)
        if is_lib and not is_blk:     cat = "lib_rm"
        elif is_lib and is_blk:       cat = "lib_blk"
        elif not is_lib and is_blk:   cat = "kern_blk"
        else:                         cat = "kern_rm"
        # Filter by max_sb for non-library blocked configs
        if max_sb is not None and cat == "kern_blk":
            try:
                sb = int(key[3])  # CPU key: (variant, N, TB, SB, MT, threads)
                if sb > max_sb:
                    continue
            except (ValueError, IndexError):
                pass
        buckets[cat][key] = np.array(gbs_list)
    cats_data, cats_key = {}, {}
    for c, b in buckets.items():
        bk, bd = best_of(b)
        if bd is not None:
            cats_data[c] = bd
            cats_key[c] = bk
    return cats_data, cats_key

def categorise_gpu(groups, max_sb=None):
    """Split GPU groups into {lib_rm, kern_rm, lib_blk, kern_blk}.
    Returns (cats_data, cats_bestkey).  If max_sb is set, blocked non-library
    configs with SB > max_sb are excluded."""
    buckets = {c: {} for c in ("lib_rm", "kern_rm", "lib_blk", "kern_blk")}
    for key, gbs_list in groups.items():
        variant = key[0]
        is_lib = variant in GPU_LIB
        is_blk = is_blocked_gpu(variant)
        if is_lib and not is_blk:     cat = "lib_rm"
        elif is_lib and is_blk:       cat = "lib_blk"
        elif not is_lib and is_blk:   cat = "kern_blk"
        else:                         cat = "kern_rm"
        # Filter by max_sb for non-library blocked configs
        if max_sb is not None and cat == "kern_blk":
            try:
                sb = int(key[6])  # GPU key: (variant, N, BX, BY, TX, TY, SB, PAD)
                if sb > max_sb:
                    continue
            except (ValueError, IndexError):
                pass
        buckets[cat][key] = np.array(gbs_list)
    cats_data, cats_key = {}, {}
    for c, b in buckets.items():
        bk, bd = best_of(b)
        if bd is not None:
            cats_data[c] = bd
            cats_key[c] = bk
    return cats_data, cats_key

def format_cpu_key(key):
    """(variant, N, TB, SB, MT, threads) → readable string."""
    return f"{key[0]} N={key[1]} TB={key[2]} SB={key[3]} MT={key[4]} thr={key[5]}"

def format_gpu_key(key):
    """(variant, N, BX, BY, TX, TY, SB, PAD) → readable string."""
    return f"{key[0]} N={key[1]} BX={key[2]} BY={key[3]} TX={key[4]} TY={key[5]} SB={key[6]} PAD={key[7]}"

# ══════════════════════════════════════════════════════════════════════
#  Metrics CSV loading
# ══════════════════════════════════════════════════════════════════════

def load_metrics(path):
    """Load transpose_metrics.csv → dict[(layout, schedule, SB)] → {mu, delta, mu_delta}."""
    metrics = {}
    try:
        with open(path) as f:
            reader = csvmod.DictReader(f)
            for row in reader:
                key = (row["layout"], row["schedule"], int(row["SB"]))
                metrics[key] = {
                    "mu": float(row["mu"]),
                    "delta": float(row["delta"]),
                    "mu_delta": float(row["mu_delta"]),
                }
    except FileNotFoundError:
        pass
    return metrics

# ══════════════════════════════════════════════════════════════════════
#  Panel drawing
# ══════════════════════════════════════════════════════════════════════

def draw_panel(ax, cats, title, peak, add_peak, xlabels_map):
    order = ["lib_rm", "kern_rm", "lib_blk", "kern_blk"]
    present = [k for k in order if k in cats]
    if not present:
        ax.set_title(title); return

    positions, data, colors, xlabels = [], [], [], []
    medians = []
    pos = 0
    sep_x = None

    for vk in present:
        # separator between row-major and blocked
        if vk == "lib_blk" and sep_x is None and pos > 0:
            sep_x = pos - 0.5
            pos += 0.4
        arr = remove_outliers(cats[vk])
        if len(arr) == 0:
            pos += 1; continue
        data.append(arr)
        positions.append(pos)
        colors.append(VCOL[vk])
        medians.append((pos, float(np.median(arr)), vk, float(np.min(arr))))
        xlabels.append(xlabels_map[vk])
        pos += 1

    # y-axis
    ymax = float(np.max(np.concatenate(data))) if data else 1
    if add_peak and peak and peak > ymax:
        ymax = peak
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
        parts["cmeans"].set_color("black")
        parts["cmedians"].set_color("white")

    # separator line
    if sep_x is not None:
        ax.axvline(x=sep_x, color="gray", ls="--", lw=1.5, alpha=0.6)
        ax.text(sep_x - 0.1, top * 0.08, "Row-Major",
                ha="right", va="top", fontsize=8, color="gray", fontweight="bold")
        ax.text(sep_x + 0.1, top * 0.08, "Blocked",
                ha="left", va="top", fontsize=8, color="gray", fontweight="bold")

    ax.set_xticks(positions)
    ax.set_xticklabels(xlabels, fontsize=7.5)
    ax.set_yticks(ticks)
    if peak > 3.76 and peak < 3.9:
        top *= 1.033
    ax.set_ylim(0, top)
    from matplotlib.ticker import FormatStrFormatter
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.tick_params(axis='y', which='minor', length=3)
    ax.set_title(title, fontsize=11)
    ax.grid(axis="y", alpha=0.25)
    ax.grid(axis="y", which='minor', alpha=0.12, ls=':')

    # STREAM peak line + label
    if peak:
        ax.axhline(y=peak, color="dimgray", ls="--", lw=1, alpha=0.5)
        ax.text(0.03, 0.97, f"STREAM {peak:.2f} TB/s",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=8, color="dimgray")
        if peak > ymax and top < peak * 1.1:
            ax.set_ylim(0, peak * 1.1)

    # % annotations
    if peak:
        off = 0.04 * top
        for p, med, vk, vmin in medians:
            pct = 100.0 * med / peak
            if pct < 14:
                ax.text(p, med + off, f"{pct:.0f}%",
                        ha="center", va="bottom", fontsize=10,
                        color=VCOL[vk], fontweight="bold")
            else:
                ax.text(p, vmin - off, f"{pct:.0f}%",
                        ha="center", va="top", fontsize=10,
                        color=VCOL[vk], fontweight="bold")

# ══════════════════════════════════════════════════════════════════════
#  LaTeX + text table generation
# ══════════════════════════════════════════════════════════════════════

def generate_tables(grid):
    """Print text + LaTeX table combining bandwidth results with cost-model metrics."""
    # Collect rows: (platform, layout_label, variant_label, BW, %peak, mu, delta, mu_delta)
    rows = []
    for (rk, ck) in sorted(grid.keys()):
        title, cats, cats_key, peak, xmap, fmt_key, metrics = grid[(rk, ck)]
        is_cpu = (ck == "cpu")
        for cat in ["lib_rm", "kern_rm", "lib_blk", "kern_blk"]:
            if cat not in cats:
                continue
            med = float(np.median(cats[cat]))
            pct = 100.0 * med / peak if peak else 0
            best_key = cats_key.get(cat)
            best_str = fmt_key(best_key) if best_key else "?"

            # Determine layout/schedule for metrics lookup
            is_blk = cat.endswith("_blk")
            # Try to extract SB from best_key
            sb = 0
            if best_key:
                if is_cpu:
                    # CPU key: (variant, N, TB, SB, MT, threads)
                    try: sb = int(best_key[3])
                    except: pass
                else:
                    # GPU key: (variant, N, BX, BY, TX, TY, SB, PAD)
                    try: sb = int(best_key[6])
                    except: pass

            mu_val, delta_val, mudelta_val = None, None, None
            if metrics:
                # For row-major "Ours": look up best tiled schedule
                # For blocked "Ours": look up blocked layout with blk_aligned schedule
                if not is_blk:
                    # row-major naive or best tiled
                    mk = ("row_major", "naive", 0)
                    if mk in metrics:
                        mu_val = metrics[mk]["mu"]
                        delta_val = metrics[mk]["delta"]
                        mudelta_val = metrics[mk]["mu_delta"]
                    # Also try tiled with matching SB
                    for sched in ["tiled", "blk_aligned"]:
                        mk2 = ("row_major", sched, sb)
                        if mk2 in metrics:
                            mu_val = metrics[mk2]["mu"]
                            delta_val = metrics[mk2]["delta"]
                            mudelta_val = metrics[mk2]["mu_delta"]
                else:
                    mk = ("blocked", "blk_aligned", sb)
                    if mk in metrics:
                        mu_val = metrics[mk]["mu"]
                        delta_val = metrics[mk]["delta"]
                        mudelta_val = metrics[mk]["mu_delta"]

            rows.append({
                "platform": title,
                "cat": cat,
                "label": xmap[cat].replace("\n", " "),
                "best_config": best_str,
                "bw": med,
                "pct": pct,
                "mu": mu_val,
                "delta": delta_val,
                "mu_delta": mudelta_val,
            })

    # ── Text table ──
    print("\n" + "=" * 120)
    print("  Summary Table: Bandwidth + Cost Model Metrics")
    print("=" * 120)
    hdr = f"{'Platform':<22} {'Category':<18} {'BW (TB/s)':>10} {'%Peak':>7} {'mu':>8} {'delta':>10} {'mu*delta':>12}  Best Config"
    print(hdr)
    print("-" * 120)
    for r in rows:
        mu_s = f"{r['mu']:.4f}" if r['mu'] is not None else "—"
        d_s  = f"{r['delta']:.2f}" if r['delta'] is not None else "—"
        md_s = f"{r['mu_delta']:.2f}" if r['mu_delta'] is not None else "—"
        print(f"{r['platform']:<22} {r['label']:<18} {r['bw']:10.2f} {r['pct']:6.1f}% {mu_s:>8} {d_s:>10} {md_s:>12}  {r['best_config']}")

    # ── LaTeX table ──
    print("\n" + "=" * 80)
    print("  LaTeX Table")
    print("=" * 80)
    print(r"""\begin{table}[!htbp]
\centering\small
\renewcommand{\arraystretch}{1.1}
\begin{tabular}{@{}llrrrrrr@{}}
\toprule
\textbf{Platform} & \textbf{Layout} & \textbf{BW (TB/s)} & \textbf{\%Peak} & $\mu$ & $\Delta$ & $\mu \cdot \Delta$ \\
\midrule""")

    prev_platform = None
    for r in rows:
        if r["platform"] != prev_platform:
            if prev_platform is not None:
                print(r"\midrule")
            prev_platform = r["platform"]

        mu_s  = f"{r['mu']:.4f}" if r['mu'] is not None else "---"
        d_s   = f"{r['delta']:.2f}" if r['delta'] is not None else "---"
        md_s  = f"{r['mu_delta']:.2f}" if r['mu_delta'] is not None else "---"
        pct_s = f"{r['pct']:.1f}\\%"
        label = r["label"].replace("(blocked)", "(blk.)")
        plat  = r["platform"].replace("_", r"\_")

        print(f"{plat} & {label} & {r['bw']:.2f} & {pct_s} & {mu_s} & {d_s} & {md_s} \\\\")

    print(r"""\bottomrule
\end{tabular}
\caption{Matrix transpose: achieved bandwidth and cost-model metrics ($\mu$, $\Delta$, $\mu \cdot \Delta$) for row-major and blocked layouts on each platform.}
\label{tab:transpose-metrics}
\end{table}""")

# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amd-cpu", default="results/beverin/transpose_cpu_raw.csv")
    ap.add_argument("--amd-gpu", default="results/beverin/transpose_raw.csv")
    ap.add_argument("--nv-cpu",  default="results/daint/transpose_cpu_raw.csv")
    ap.add_argument("--nv-gpu",  default="results/daint/transpose_raw.csv")
    ap.add_argument("--amd-cpu-label", default="MI300A Zen CPU")
    ap.add_argument("--amd-gpu-label", default="MI300A GPU")
    ap.add_argument("--nv-cpu-label",  default="GH200 Grace CPU")
    ap.add_argument("--nv-gpu-label",  default="GH200 Hopper GPU")
    ap.add_argument("--metrics-csv", default=None,
                    help="Single metrics CSV (used for both CPU and GPU if separate not given)")
    ap.add_argument("--metrics-cpu", default="transpose_metrics_cpu.csv",
                    help="CPU metrics CSV (B=8, 32-byte lines)")
    ap.add_argument("--metrics-gpu", default="transpose_metrics_gpu.csv",
                    help="GPU metrics CSV (B=16, 64-byte lines)")
    ap.add_argument("--add-peak", action="store_true", default=True, help="Enable STREAM peak line")
    ap.add_argument("--amd-max-sb", type=int, default=128,
                    help="Max SB for AMD non-library blocked configs (default 128)")
    ap.add_argument("--nv-max-sb", type=int, default=None,
                    help="Max SB for NVIDIA non-library blocked configs (default: no limit)")
    args = ap.parse_args()

    if args.metrics_csv:
        metrics_cpu = load_metrics(args.metrics_csv)
        metrics_gpu = metrics_cpu
        print(f"Loaded {len(metrics_cpu)} metric entries from {args.metrics_csv} (shared)")
    else:
        metrics_cpu = load_metrics(args.metrics_cpu)
        metrics_gpu = load_metrics(args.metrics_gpu)
        if metrics_cpu:
            print(f"Loaded {len(metrics_cpu)} CPU metric entries from {args.metrics_cpu}")
        else:
            print(f"[WARN] No CPU metrics loaded from {args.metrics_cpu}")
        if metrics_gpu:
            print(f"Loaded {len(metrics_gpu)} GPU metric entries from {args.metrics_gpu}")
        else:
            print(f"[WARN] No GPU metrics loaded from {args.metrics_gpu}")

    # ── Build 2×2 grid ──
    # (row, col) → (title, cats_data, cats_key, peak, xlabels_map, fmt_key_fn, metrics_dict)
    grid = {}

    def try_cpu(csv, row, label, max_sb=None):
        try:
            cats, keys = categorise_cpu(parse_cpu_raw(csv), max_sb=max_sb)
            peak = STREAM_PEAK.get(label, 0)
            xmap = {"lib_rm": "HPTT", "kern_rm": "Ours",
                    "lib_blk": "HPTT\n(blocked)", "kern_blk": "Ours\n(blocked)"}
            grid[(row, "cpu")] = (label, cats, keys, peak, xmap, format_cpu_key, metrics_cpu)
            print(f"  {label}: {len(cats)} categories" +
                  (f" (max_sb={max_sb})" if max_sb else ""))
            for cat in ["lib_rm", "kern_rm", "lib_blk", "kern_blk"]:
                if cat in keys:
                    print(f"    {cat:12s} best: {format_cpu_key(keys[cat])}")
        except FileNotFoundError:
            print(f"  [WARN] {csv} not found")

    def try_gpu(csv, row, label, amd, max_sb=None):
        try:
            cats, keys = categorise_gpu(parse_gpu_raw(csv), max_sb=max_sb)
            peak = STREAM_PEAK.get(label, 0)
            lib_name = "hipTensor" if amd else "cuTENSOR"
            xmap = {"lib_rm": lib_name, "kern_rm": "Ours",
                    "lib_blk": f"{lib_name}\n(blocked)", "kern_blk": "Ours\n(blocked)"}
            grid[(row, "gpu")] = (label, cats, keys, peak, xmap, format_gpu_key, metrics_gpu)
            print(f"  {label}: {len(cats)} categories" +
                  (f" (max_sb={max_sb})" if max_sb else ""))
            for cat in ["lib_rm", "kern_rm", "lib_blk", "kern_blk"]:
                if cat in keys:
                    print(f"    {cat:12s} best: {format_gpu_key(keys[cat])}")
        except FileNotFoundError:
            print(f"  [WARN] {csv} not found")

    print("── Loading data ──")
    try_cpu(args.amd_cpu, "amd", args.amd_cpu_label, max_sb=args.amd_max_sb)
    try_gpu(args.amd_gpu, "amd", args.amd_gpu_label, amd=True, max_sb=args.amd_max_sb)
    try_cpu(args.nv_cpu,  "nv",  args.nv_cpu_label,  max_sb=args.nv_max_sb)
    try_gpu(args.nv_gpu,  "nv",  args.nv_gpu_label,  amd=False, max_sb=args.nv_max_sb)

    if not grid:
        print("No data."); return

    # ── Plot 2×2 ──
    rows_order = ["amd", "nv"]
    cols_order = ["cpu", "gpu"]

    fig, axes = plt.subplots(2, 2, figsize=(3.6 * 2, 2.8 * 2), squeeze=False)

    for ri, rk in enumerate(rows_order):
        for ci, ck in enumerate(cols_order):
            ax = axes[ri, ci]
            if (rk, ck) not in grid:
                ax.set_visible(False); continue
            title, cats, keys, peak, xmap, fmt, mdict = grid[(rk, ck)]
            draw_panel(ax, cats, title, peak, args.add_peak, xmap)
            if ci == 0:
                ax.set_ylabel("Bandwidth [TB/s]", fontsize=11)

    fig.suptitle("Matrix Transpose: Row-Major vs Blocked Layout",
                 fontsize=15, y=0.89)
    fig.text(0.5, 0.85,
             "% annotations relative to STREAM peak bandwidth",
             ha="center", va="top", fontsize=12, color="dimgray")
    fig.tight_layout(rect=[0, 0, 1, 0.89])

    sfx = "" if args.add_peak else "_no_peak"
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT_STEM}{sfx}.{ext}", dpi=200, bbox_inches="tight")
    print(f"\nSaved {OUT_STEM}{sfx}.pdf/png")

    # ── Summary with best config ──
    for (rk, ck) in sorted(grid.keys()):
        title, cats, keys, peak, xmap, fmt, mdict = grid[(rk, ck)]
        print(f"\n  {title}:")
        for cat in ["lib_rm", "kern_rm", "lib_blk", "kern_blk"]:
            if cat not in cats: continue
            med = float(np.median(cats[cat]))
            pk = f"  ({100 * med / peak:.2f}%)" if peak else ""
            best_str = fmt(keys[cat]) if cat in keys else "?"
            print(f"    {xmap[cat]:<24} {med:8.2f} TB/s{pk}")
            print(f"      → config: {best_str}")

    # ── Generate combined table ──
    generate_tables(grid)

if __name__ == "__main__":
    main()