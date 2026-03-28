#!/usr/bin/env python3
"""
analyze_metrics.py

Rank-correlates cost-model metrics against measured runtime bandwidth.

Required:
    --target   e.g. cpu_scalar, gpu_scalar, cpu_avx512, gpu_warp32
    --runtime  Runtime CSV

Output (in order):
    1. Table sorted by measured BW, each metric with its rank
    2. Per-metric rank correlation vs BW (Spearman ρ, R², Kendall τ, Pearson r)
    3. Best predictor summary
    4. Scatter plots: metric (x) vs BW (y), ordered by BW on y-axis

Usage:
    python analyze_metrics.py --target cpu_scalar --runtime z_v_grad_w_cpu_beverin.csv
    python analyze_metrics.py --target gpu_scalar --runtime z_v_grad_w_gpu_beverin.csv --gpu-config 1x1_32x16
"""
import argparse
import pandas as pd
import numpy as np
import sys
import os
from scipy.stats import spearmanr, kendalltau, pearsonr

parser = argparse.ArgumentParser()
parser.add_argument("--csv", default="metrics.csv", help="Cost metrics CSV")
parser.add_argument("--target", required=True, help="e.g. cpu_scalar, gpu_scalar")
parser.add_argument("--runtime", required=True, help="Runtime CSV")
parser.add_argument("--gpu-config", default="1x1_32x16", help="GPU config_label filter")
parser.add_argument("--nlev", type=int, default=96)
parser.add_argument("--no-plot", action="store_true")
args = parser.parse_args()

NPROMA = 81920
BYTES  = (2*NPROMA*4 + 2*NPROMA*4 +
          args.nlev*NPROMA*8 + args.nlev*NPROMA*8 + NPROMA*8 +
          args.nlev*NPROMA*8 + args.nlev*NPROMA*8 +
          NPROMA*8 + NPROMA*8 + args.nlev*NPROMA*8)

ALL_METRICS = [
    ("mu",            "μ"),
    ("delta_raw",     "Δ_raw"),
    ("delta_uma",     "Δ_uma"),
    ("delta_numa",    "Δ_numa"),
    ("mu_delta_raw",  "μ·Δ_raw"),
    ("mu_delta_uma",  "μ·Δ_uma"),
    ("mu_delta_numa", "μ·Δ_numa"),
    ("sigma",         "σ"),
    ("cost_uma",      "ω_uma"),
    ("cost_numa",     "ω_numa"),
]

# ═══════════════════════════════════════════════════════════════
#  Load cost metrics
# ═══════════════════════════════════════════════════════════════

cost = pd.read_csv(args.csv)
cost = cost[cost["target"] == args.target].copy()
if cost.empty:
    print(f"ERROR: target={args.target} not in {args.csv}")
    print(f"Available: {sorted(pd.read_csv(args.csv)['target'].unique())}")
    sys.exit(1)

metrics = [(c, s) for c, s in ALL_METRICS if c in cost.columns]
mcols = [c for c, _ in metrics]
msyms = [s for _, s in metrics]

print(f"Cost metrics: {args.csv} → {len(cost)} rows for {args.target}")
print(f"Metrics: {msyms}")
print(f"Variants: {sorted(cost['variant'].unique())}")
print(f"Dists: {sorted(cost['cell_dist'].unique())}")
print(f"Loops: {sorted(cost['loop_order'].unique())}")

# ═══════════════════════════════════════════════════════════════
#  Load runtime → median BW per (variant, cell_dist)
# ═══════════════════════════════════════════════════════════════

rt = pd.read_csv(args.runtime)
print(f"\nRuntime: {args.runtime} → {len(rt)} rows")

# GPU filter
if "config_label" in rt.columns and args.gpu_config:
    before = len(rt)
    rt = rt[rt["config_label"] == args.gpu_config]
    print(f"  GPU filter config_label={args.gpu_config}: {before} → {len(rt)} rows")

# nlev filter
if "nlev" in rt.columns:
    rt = rt[rt["nlev"] == args.nlev]

rt["bw_tbs"] = BYTES / (rt["time_ms"] * 1e-3) / 1e12

# Group by (variant, cell_dist, parallelization) → median BW
group_cols = ["variant", "cell_dist"]
if "parallelization" in rt.columns:
    group_cols.append("parallelization")

bw_agg = (rt.groupby(group_cols)["bw_tbs"]
            .agg(bw_median="median", bw_std="std", bw_n="count")
            .reset_index())

# Best schedule per (variant, cell_dist)
if "parallelization" in bw_agg.columns:
    idx = bw_agg.groupby(["variant", "cell_dist"])["bw_median"].idxmax()
    bw_best = bw_agg.loc[idx].reset_index(drop=True)
else:
    bw_best = bw_agg

print(f"  Best-schedule configs: {len(bw_best)}")

# ═══════════════════════════════════════════════════════════════
#  Join: runtime BW ↔ cost metrics
# ═══════════════════════════════════════════════════════════════

rows = []
for _, rb in bw_best.iterrows():
    var = int(rb["variant"])
    dist = rb["cell_dist"]

    mask = (cost["variant"] == var) & (cost["cell_dist"] == dist)
    crows = cost[mask]

    if crows.empty:
        print(f"  [WARN] No metric for V{var} {dist}")
        continue

    # Pick expected loop_order
    if len(crows) > 1:
        expected = "klon_first" if var <= 2 else "klev_first"
        filtered = crows[crows["loop_order"] == expected]
        if not filtered.empty:
            crows = filtered

    c = crows.iloc[0]
    row = {
        "variant": var, "dist": dist,
        "loop_order": c["loop_order"],
        "bw_median": rb["bw_median"],
        "bw_std": rb["bw_std"],
        "bw_n": int(rb["bw_n"]),
    }
    if "parallelization" in rb.index:
        row["schedule"] = rb["parallelization"]
    for mc, _ in metrics:
        row[mc] = c[mc]
    rows.append(row)

if not rows:
    print("\nERROR: No rows joined. Check variant/dist names match between CSVs.")
    sys.exit(1)

J = pd.DataFrame(rows)

# ═══════════════════════════════════════════════════════════════
#  Sort by BW, assign ranks
# ═══════════════════════════════════════════════════════════════

J = J.sort_values("bw_median", ascending=False).reset_index(drop=True)
J["bw_rank"] = J["bw_median"].rank(ascending=False).astype(int)
for mc, _ in metrics:
    J[f"{mc}_rank"] = J[mc].rank(ascending=True).astype(int)

n = len(J)
has_sched = "schedule" in J.columns

# ═══════════════════════════════════════════════════════════════
#  1) TABLE: sorted by BW (best first), metrics with ranks
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*200}")
print(f"  RANKED TABLE: {args.target}  (n={n}, sorted by measured BW)")
print(f"  Higher BW = better.  Metric rank 1 = lowest metric value = predicted best.")
print(f"  '=' after rank means metric rank matches BW rank exactly.")
print(f"{'='*200}\n")

# Header
hdr = f"  {'#':>2} {'V':>2} {'Dist':<14}"
if has_sched: hdr += f" {'Schedule':<18}"
hdr += f" {'Loop':<12} {'BW TB/s':>9} {'±':>6}"
for _, s in metrics:
    hdr += f"  {s:>8}(rk)"
print(hdr)
print(f"  {'-' * (len(hdr) + 5)}")

for _, r in J.iterrows():
    line = f"  {r['bw_rank']:>2} V{int(r['variant']):>1} {r['dist']:<14}"
    if has_sched: line += f" {r['schedule']:<18}"
    line += f" {r['loop_order']:<12} {r['bw_median']:>9.4f} {r['bw_std']:>6.4f}"
    for mc, _ in metrics:
        v = r[mc]
        rk = r[f"{mc}_rank"]
        match = "=" if rk == r["bw_rank"] else " "
        if abs(v) < 100:
            line += f"  {v:>7.2f}({rk:>2}{match})"
        elif abs(v) < 100000:
            line += f"  {v:>7.0f}({rk:>2}{match})"
        else:
            line += f"  {v:>7.0f}({rk:>2}{match})"
    print(line)

# ═══════════════════════════════════════════════════════════════
#  2) RANK CORRELATION QUALITY: per metric vs BW
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*120}")
print(f"  METRIC QUALITY: rank correlation vs measured BW  (n={n})")
print(f"  Positive ρ = higher metric → lower BW (metric correctly predicts cost)")
print(f"{'='*120}\n")

print(f"  {'Metric':<14} {'Spearman ρ':>10} {'ρ²':>8} {'Kendall τ':>10} "
      f"{'Pearson r':>10} {'r²':>8} {'p(ρ)':>12} {'p(τ)':>12} {'exact':>7}")
print(f"  {'-'*100}")

corr_results = []
bw = J["bw_median"].values

for mc, ms in metrics:
    vals = J[mc].values
    ok = np.isfinite(bw) & np.isfinite(vals)
    nok = ok.sum()
    if nok < 3:
        print(f"  {ms:<14} n<3")
        continue

    rho, p_rho  = spearmanr(vals[ok], -bw[ok])
    tau, p_tau  = kendalltau(vals[ok], -bw[ok])
    pr, p_pr    = pearsonr(vals[ok], bw[ok])

    rho2 = rho**2
    pr2  = pr**2
    n_exact = int(np.sum(J[f"{mc}_rank"].values[ok] == J["bw_rank"].values[ok]))

    sig = "***" if p_rho < 0.001 else "**" if p_rho < 0.01 else "*" if p_rho < 0.05 else ""

    print(f"  {ms:<14} {rho:>+10.4f} {rho2:>8.4f} {tau:>+10.4f} "
          f"{pr:>+10.4f} {pr2:>8.4f} {p_rho:>12.2e} {p_tau:>12.2e} {n_exact:>3}/{nok} {sig}")

    corr_results.append({
        "metric": ms, "col": mc,
        "rho": rho, "rho2": rho2, "tau": tau,
        "pr": pr, "pr2": pr2,
        "p_rho": p_rho, "p_tau": p_tau,
        "n_exact": n_exact, "n": nok,
    })

# Best / worst
if corr_results:
    by_rho = sorted(corr_results, key=lambda x: abs(x["rho"]), reverse=True)
    print(f"\n  Best predictor (|ρ|):  {by_rho[0]['metric']:<14} ρ={by_rho[0]['rho']:+.4f}  ρ²={by_rho[0]['rho2']:.4f}  exact={by_rho[0]['n_exact']}/{by_rho[0]['n']}")
    print(f"  Worst predictor (|ρ|): {by_rho[-1]['metric']:<14} ρ={by_rho[-1]['rho']:+.4f}  ρ²={by_rho[-1]['rho2']:.4f}  exact={by_rho[-1]['n_exact']}/{by_rho[-1]['n']}")

    by_exact = sorted(corr_results, key=lambda x: x["n_exact"], reverse=True)
    if by_exact[0]["metric"] != by_rho[0]["metric"]:
        print(f"  Most exact matches:    {by_exact[0]['metric']:<14} exact={by_exact[0]['n_exact']}/{by_exact[0]['n']}")

# ═══════════════════════════════════════════════════════════════
#  3) PLOTS: metric vs BW, sorted by BW on y-axis
# ═══════════════════════════════════════════════════════════════

if not args.no_plot:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        LOOP_COLOR = {"klon_first": "#e67e22", "klev_first": "#2980b9"}
        LOOP_MARKER = {"klon_first": "o", "klev_first": "s"}

        nm = len(metrics)
        ncols = min(4, nm)
        nrows = (nm + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.5*ncols, 4*nrows), squeeze=False)
        fig.suptitle(f"Metric vs Measured BW — {args.target}", fontsize=14)

        for mi, (mc, ms) in enumerate(metrics):
            ax = axes[mi // ncols, mi % ncols]

            for lo in ["klon_first", "klev_first"]:
                sub = J[J["loop_order"] == lo]
                if sub.empty: continue
                ax.scatter(sub[mc], sub["bw_median"],
                           c=LOOP_COLOR.get(lo, "gray"),
                           marker=LOOP_MARKER.get(lo, "o"),
                           s=60, alpha=0.85, edgecolors="black", linewidths=0.4,
                           label=lo.replace("_", "-"), zorder=3)

                for _, r in sub.iterrows():
                    label = f"V{int(r['variant'])}/{r['dist'][:4]}"
                    ax.annotate(label, (r[mc], r["bw_median"]),
                                fontsize=5, alpha=0.7, xytext=(4, 4),
                                textcoords="offset points")

            # Trend line (log-log if range is large)
            ok = np.isfinite(J[mc]) & np.isfinite(J["bw_median"])
            xv, yv = J.loc[ok, mc].values, J.loc[ok, "bw_median"].values
            use_log = len(xv) > 0 and xv.max() / max(xv.min(), 1e-12) > 10

            if ok.sum() >= 3:
                rho_i = [x for x in corr_results if x["col"] == mc]
                if rho_i:
                    r = rho_i[0]
                    ax.text(0.04, 0.04,
                            f"ρ={r['rho']:+.3f}  ρ²={r['rho2']:.3f}\n"
                            f"τ={r['tau']:+.3f}  r={r['pr']:+.3f}\n"
                            f"exact={r['n_exact']}/{r['n']}",
                            transform=ax.transAxes, fontsize=7, va="bottom",
                            family="monospace",
                            bbox=dict(fc="lightyellow", alpha=0.9, pad=2))

            ax.set_xlabel(ms, fontsize=11)
            ax.set_ylabel("BW [TB/s]" if mi % ncols == 0 else "")
            ax.set_title(ms, fontsize=11)
            ax.grid(alpha=0.2)
            if use_log:
                ax.set_xscale("log")
            if mi == 0:
                ax.legend(fontsize=7, loc="upper right")

        for mi in range(nm, nrows * ncols):
            axes[mi // ncols, mi % ncols].set_visible(False)

        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fname = f"metric_vs_bw_{args.target}"
        fig.savefig(f"{fname}.png", dpi=180, bbox_inches="tight")
        fig.savefig(f"{fname}.pdf", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"\n  Scatter saved: {fname}.png/pdf")

        # ── Bar chart: ρ² per metric ──
        if corr_results:
            fig2, ax2 = plt.subplots(figsize=(max(6, 0.8*nm), 4))
            names = [r["metric"] for r in corr_results]
            rho2s = [r["rho2"] for r in corr_results]
            colors = ["#2ecc71" if r["p_rho"] < 0.05 else "#e74c3c" for r in corr_results]

            bars = ax2.bar(range(len(names)), rho2s, color=colors,
                           edgecolor="black", linewidth=0.5)
            ax2.set_xticks(range(len(names)))
            ax2.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
            ax2.set_ylabel("Spearman ρ²")
            ax2.set_title(f"Metric Quality: ρ² vs BW — {args.target}")
            ax2.set_ylim(0, 1.05)
            ax2.axhline(y=1.0, color="gray", linestyle=":", linewidth=0.5)
            ax2.grid(axis="y", alpha=0.3)

            # Annotate bars
            for i, (b, r) in enumerate(zip(bars, corr_results)):
                ax2.text(i, b.get_height() + 0.02,
                         f"ρ={r['rho']:+.2f}\n{r['n_exact']}/{r['n']}",
                         ha="center", va="bottom", fontsize=7)

            # Legend: green = significant, red = not
            from matplotlib.patches import Patch
            ax2.legend(handles=[
                Patch(facecolor="#2ecc71", edgecolor="black", label="p < 0.05"),
                Patch(facecolor="#e74c3c", edgecolor="black", label="p ≥ 0.05"),
            ], fontsize=8, loc="upper right")

            fig2.tight_layout()
            fname2 = f"metric_quality_{args.target}"
            fig2.savefig(f"{fname2}.png", dpi=180, bbox_inches="tight")
            fig2.savefig(f"{fname2}.pdf", dpi=180, bbox_inches="tight")
            plt.close(fig2)
            print(f"  Quality bar saved: {fname2}.png/pdf")

    except ImportError:
        print("\n  [WARN] matplotlib not available")

print("\nDone.")