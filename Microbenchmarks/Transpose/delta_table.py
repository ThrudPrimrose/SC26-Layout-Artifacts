#!/usr/bin/env python3
"""
gen_delta_table.py — Generate LaTeX table of Δ vs SB for the transpose,
plus a mapping of which platform chose which block size.

Usage:
    python gen_delta_table.py [--cpu FILE] [--gpu FILE]
"""
import argparse, csv

# ══════════════════════════════════════════════════════════════════════
#  Known best SB per platform (from benchmark sweep)
# ══════════════════════════════════════════════════════════════════════
BEST_SB = {
    "MI300A Zen4 CPU":   16,   # locbuf_blk_c2_nd SB=16
    "Grace CPU":          8,   # locbuf_blk_mt_c2_nd SB=8
    "MI300A GPU":       256,   # smem_blk_swiz SB=256
    "Hopper GPU":        32,   # smem_pad_blk SB=32
}

def load_metrics(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            key = (row["layout"], row["schedule"], int(row["SB"]))
            out[key] = {"mu": float(row["mu"]), "delta": float(row["delta"])}
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", default="transpose_metrics_cpu.csv")
    ap.add_argument("--gpu", default="transpose_metrics_gpu.csv")
    ap.add_argument("--N", type=int, default=16384)
    ap.add_argument("--max-sb", type=int, default=128,
                    help="Maximum block size to include (default 128)")
    args = ap.parse_args()

    N = args.N
    mcpu = load_metrics(args.cpu)
    mgpu = load_metrics(args.gpu)

    # Collect SBs, filtered by max-sb
    sbs = sorted({k[2] for k in mcpu if k[0] == "blocked" and k[1] == "blk_aligned" and k[2] <= args.max_sb}
               & {k[2] for k in mgpu if k[0] == "blocked" and k[1] == "blk_aligned" and k[2] <= args.max_sb})

    # Row-major baselines
    rm_cpu = mcpu.get(("row_major", "naive", 0), {})
    rm_gpu = mgpu.get(("row_major", "naive", 0), {})

    # Blocked entries
    rows = []
    for sb in sbs:
        dc = mcpu.get(("blocked", "blk_aligned", sb), {})
        dg = mgpu.get(("blocked", "blk_aligned", sb), {})
        rows.append({"sb": sb, "delta_cpu": dc.get("delta"), "delta_gpu": dg.get("delta")})

    # ── Text table ──
    print(f"N = {N}")
    print(f"Row-major: Δ_CPU = {rm_cpu.get('delta', 0):.1f},  Δ_GPU = {rm_gpu.get('delta', 0):.1f}\n")
    print(f"{'SB':>4}  {'Δ (64B)':>8}  {'Δ (128B)':>9}  Selected by")
    print("-" * 60)

    # Filter BEST_SB: if a platform's chosen SB exceeds max-sb,
    # pick the best SB within the allowed range (lowest Δ on its cache line)
    filtered_best = {}
    for plat, sb in BEST_SB.items():
        if sb <= args.max_sb:
            filtered_best[plat] = sb
        else:
            # Reselect: CPU platforms use 64B metrics, GPU use 128B
            is_gpu = "GPU" in plat
            m = mgpu if is_gpu else mcpu
            best_sb, best_delta = None, float('inf')
            for s in sbs:
                k = ("blocked", "blk_aligned", s)
                if k in m and m[k]["delta"] < best_delta:
                    best_delta = m[k]["delta"]
                    best_sb = s
            if best_sb:
                filtered_best[plat] = best_sb
                print(f"  [NOTE] {plat}: SB={sb} > max-sb={args.max_sb}, "
                      f"reselected SB={best_sb} (Δ={best_delta:.1f})")

    sb_to_plats = {}
    for plat, sb in filtered_best.items():
        sb_to_plats.setdefault(sb, []).append(plat)

    for r in rows:
        dc = f"{r['delta_cpu']:.1f}" if r["delta_cpu"] else "—"
        dg = f"{r['delta_gpu']:.1f}" if r["delta_gpu"] else "—"
        plats = ", ".join(sb_to_plats.get(r["sb"], []))
        print(f"{r['sb']:>4}  {dc:>8}  {dg:>9}  {plats}")
    dc = f"{rm_cpu.get('delta', 0):.1f}"
    dg = f"{rm_gpu.get('delta', 0):.1f}"
    print(f"{'—':>4}  {dc:>8}  {dg:>9}  (row-major)")

    # ── LaTeX ──
    min_dc = min((r["delta_cpu"] for r in rows if r["delta_cpu"] is not None), default=None)
    min_dg = min((r["delta_gpu"] for r in rows if r["delta_gpu"] is not None), default=None)

    print("\n\n% ─── LaTeX table ───")
    print(r"""\begin{table}[!htbp]
\centering\small
\renewcommand{\arraystretch}{1.15}
\begin{tabular}{@{}r rr l@{}}
\toprule
$SB$ & $\Delta_{\text{64\,B}}$ & $\Delta_{\text{128\,B}}$ & Selected by \\
\midrule""")

    for r in rows:
        sb = r["sb"]

        def fmt(v, vmin):
            if v is None: return "---"
            s = f"{v:.1f}"
            return f"\\textbf{{{s}}}" if v == vmin else s

        dc = fmt(r["delta_cpu"], min_dc)
        dg = fmt(r["delta_gpu"], min_dg)
        plats = sb_to_plats.get(sb, [])
        plat_str = ", ".join(plats) if plats else ""

        print(f"{sb} & {dc} & {dg} & {plat_str} \\\\")

    dc_rm = f"{rm_cpu.get('delta', 0):.1f}"
    dg_rm = f"{rm_gpu.get('delta', 0):.1f}"
    print(r"\midrule")
    print(f"--- & {dc_rm} & {dg_rm} & (row-major) \\\\")

    print(r"""\bottomrule
\end{tabular}
\caption{Average block distance~$\Delta$ for the $""" + str(N) + r"""\times""" + str(N) + r"""$~\texttt{fp32} matrix transpose under blocked layout with block-aligned schedule.  $\Delta_{\text{64\,B}}$ and $\Delta_{\text{128\,B}}$ correspond to 64-byte (CPU) and 128-byte (GPU) cache lines, respectively.  The last column lists which platform's best configuration selected each block size.  Bold marks the minimum~$\Delta$ per column.}
\label{tab:delta-blocking}
\end{table}""")

if __name__ == "__main__":
    main()