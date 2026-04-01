#!/usr/bin/env python3
"""
analyze_metrics.py

Rank-correlates cost-model metrics against measured runtime bandwidth.

Required: --target, --runtime
Usage:
    python analyze_metrics.py --target cpu_scalar --runtime z_v_grad_w_cpu.csv
    python analyze_metrics.py --target gpu_scalar --runtime z_v_grad_w_gpu.csv --gpu-config 1x1_32x16
    python analyze_metrics.py --target cpu_scalar --runtime z_v_grad_w_cpu.csv --blocked
"""
import argparse, sys, warnings
import pandas as pd, numpy as np
from scipy.stats import spearmanr, kendalltau, pearsonr
warnings.filterwarnings("ignore", message="An input array is constant")

parser = argparse.ArgumentParser()
parser.add_argument("--csv", default="metrics.csv")
parser.add_argument("--target", required=True)
parser.add_argument("--runtime", required=True)
parser.add_argument("--gpu-config", default="1x1_32x16")
parser.add_argument("--nlev", type=int, default=90)
parser.add_argument("--blocked", action="store_true",
                    help="Analyze blocked configs (variant=0, blocking>0)")
parser.add_argument("--no-plot", action="store_true")
args = parser.parse_args()

def loop_order_for(variant):
    """V1/V2 -> klon_first, V3/V4 -> klev_first, blocked -> klon_first"""
    if variant <= 0:
        return "klon_first"
    return "klon_first" if variant <= 2 else "klev_first"

# ── Byte count for bandwidth calculation ──
# 2 idx arrays (int, N*2 each) + 5 2D arrays (double, N*nlev each) + 3 1D arrays (double, N each)
def compute_bytes(N, nlev):
    return (2*N*2*4 + 5*N*nlev*8 + 3*N*8)

ALL_METRICS = [
    ("mu",              "mu",       "cost"),
    ("mu_delta",        "muD",      "cost"),
    ("mu_delta_numa",   "muDn",     "cost"),
    ("delta_max",       "D_max",    "cost"),
    ("delta_min_numa",  "Dn_avg",   "cost"),
]

# ══════════════════════════════════════════════════════════════
#  Load cost metrics
# ══════════════════════════════════════════════════════════════
cost = pd.read_csv(args.csv)
cost = cost[cost["target"] == args.target].copy()
if cost.empty:
    print(f"ERROR: target={args.target} not in {args.csv}")
    print(f"Available: {sorted(pd.read_csv(args.csv)['target'].unique())}")
    sys.exit(1)

# Derive loop_order from variant
cost["loop_order"] = cost["variant"].apply(loop_order_for)

# Filter: unblocked (variant>0, blocking==0) vs blocked (variant==0, blocking>0)
if "blocking" in cost.columns:
    if args.blocked:
        cost = cost[cost["blocking"] > 0].copy()
        print(f"Mode: BLOCKED (B values: {sorted(cost['blocking'].unique())})")
    else:
        cost = cost[cost["blocking"] == 0].copy()
        print(f"Mode: UNBLOCKED (variants: {sorted(cost['variant'].unique())})")

metrics = [(c, s, d) for c, s, d in ALL_METRICS if c in cost.columns]
print(f"Target: {args.target}  ({len(cost)} cost rows)")
print(f"Metrics: {[s for _,s,_ in metrics]}")

# ══════════════════════════════════════════════════════════════
#  Load runtime
# ══════════════════════════════════════════════════════════════
rt = pd.read_csv(args.runtime)

# Determine N from runtime or default
if "nproma" in rt.columns:
    N_vals = rt["nproma"].unique()
    N = int(N_vals[0]) if len(N_vals) == 1 else int(rt["nproma"].mode().iloc[0])
else:
    N = 81920

BYTES = compute_bytes(N, args.nlev)
rt["bw_tbs"] = BYTES / (rt["time_ms"] * 1e-3) / 1e12

# GPU: filter by config
if "config_label" in rt.columns and args.gpu_config:
    rt = rt[rt["config_label"] == args.gpu_config]

# Filter nlev
if "nlev" in rt.columns:
    rt = rt[rt["nlev"] == args.nlev]

# Filter blocked vs unblocked in runtime
if "blocking" in rt.columns:
    if args.blocked:
        rt = rt[rt["blocking"] > 0].copy()
    else:
        rt = rt[rt["blocking"] == 0].copy()

# ── Group and pick best schedule ──
group_cols = []
if args.blocked and "blocking" in rt.columns:
    group_cols = ["blocking", "cell_dist"]
else:
    group_cols = ["variant", "cell_dist"]

if "parallelization" in rt.columns:
    group_cols.append("parallelization")

bw_agg = rt.groupby(group_cols)["bw_tbs"].agg(
    bw_median="median", bw_std="std", bw_n="count").reset_index()

# Pick best schedule per (variant/blocking, dist)
key_cols = [c for c in group_cols if c != "parallelization"]
if "parallelization" in bw_agg.columns and len(key_cols) > 0:
    idx = bw_agg.groupby(key_cols)["bw_median"].idxmax()
    bw_best = bw_agg.loc[idx].reset_index(drop=True)
else:
    bw_best = bw_agg

print(f"Runtime: {len(bw_best)} configs after best-schedule selection")

# ══════════════════════════════════════════════════════════════
#  Join cost metrics with runtime
# ══════════════════════════════════════════════════════════════
rows = []
for _, rb in bw_best.iterrows():
    dist = rb["cell_dist"]

    if args.blocked:
        blk = int(rb["blocking"])
        mask = (cost["blocking"] == blk) & (cost["cell_dist"] == dist)
        crows = cost[mask]
        label_key = f"B{blk}"
        variant_val = 0
    else:
        var = int(rb["variant"])
        mask = (cost["variant"] == var) & (cost["cell_dist"] == dist)
        crows = cost[mask]
        label_key = f"V{var}"
        variant_val = var

    if crows.empty:
        print(f"  [WARN] No metric for {label_key} {dist}")
        continue

    # If multiple schedule rows, pick the one matching the runtime's best schedule
    if len(crows) > 1 and "parallelization" in rb.index:
        rt_sched = rb["parallelization"]
        # Map runtime schedule names to cost schedule names
        sched_map = {
            "omp_for": "omp_for",
            "omp_collapse2": "omp_collapse2",
            "blocked_omp_for": "omp_for",
            "blocked_collapse2": "omp_collapse2",
        }
        cost_sched = sched_map.get(rt_sched, rt_sched)
        f = crows[crows["schedule"] == cost_sched]
        if not f.empty:
            crows = f

    if len(crows) > 1:
        # Last resort: pick first
        crows = crows.head(1)

    c = crows.iloc[0]
    lo = loop_order_for(variant_val)

    row = {
        "variant": variant_val,
        "dist": dist,
        "loop_order": lo,
        "bw_median": rb["bw_median"],
        "bw_std": rb["bw_std"],
    }
    if args.blocked:
        row["blocking"] = int(rb["blocking"])
    if "parallelization" in rb.index:
        row["schedule"] = rb["parallelization"]
    for mc, _, _ in metrics:
        row[mc] = c[mc]
    rows.append(row)

if not rows:
    print("ERROR: No rows joined")
    sys.exit(1)

J = pd.DataFrame(rows).sort_values("bw_median", ascending=False).reset_index(drop=True)
J["bw_rank"] = J["bw_median"].rank(ascending=False).astype(int)
n = len(J)
has_sched = "schedule" in J.columns
has_block = "blocking" in J.columns

# Compute ranks
for mc, _, direction in metrics:
    if direction == "benefit":
        J[f"{mc}_rank"] = J[mc].rank(ascending=False).astype(int)
    else:
        J[f"{mc}_rank"] = J[mc].rank(ascending=True).astype(int)

# ══════════════════════════════════════════════════════════════
#  1) RANKED TABLE
# ══════════════════════════════════════════════════════════════

mode_str = "BLOCKED" if args.blocked else "UNBLOCKED"
print(f"\n{'='*200}")
print(f"  RANKED TABLE: {args.target} [{mode_str}]  (n={n})")
print(f"{'='*200}\n")

hdr = f"  {'#':>2}"
if has_block:
    hdr += f" {'B':>4}"
else:
    hdr += f" {'V':>2}"
hdr += f" {'Dist':<14}"
if has_sched: hdr += f" {'Schedule':<18}"
hdr += f" {'Loop':<12} {'BW':>9} {'+-':>6}"
for _, s, _ in metrics: hdr += f"  {s:>8}(rk)"
print(hdr)
print(f"  {'-'*(len(hdr)+5)}")

for _, r in J.iterrows():
    if has_block:
        line = f"  {r['bw_rank']:>2} {int(r['blocking']):>4} {r['dist']:<14}"
    else:
        line = f"  {r['bw_rank']:>2} V{int(r['variant']):>1} {r['dist']:<14}"
    if has_sched: line += f" {r['schedule']:<18}"
    line += f" {r['loop_order']:<12} {r['bw_median']:>9.4f} {r['bw_std']:>6.4f}"
    for mc, _, _ in metrics:
        v, rk = r[mc], r[f"{mc}_rank"]
        m = "=" if rk == r["bw_rank"] else " "
        if abs(v) < 100:   line += f"  {v:>7.2f}({rk:>2}{m})"
        else:               line += f"  {v:>7.0f}({rk:>2}{m})"
    print(line)

# ══════════════════════════════════════════════════════════════
#  2) RANK CORRELATION
# ══════════════════════════════════════════════════════════════

print(f"\n{'='*130}")
print(f"  METRIC QUALITY: rank correlation vs BW  (n={n})")
print(f"{'='*130}\n")

print(f"  {'Metric':<12} {'dir':>4} {'Spearman':>10} {'rho^2':>8} {'Kendall':>10} "
      f"{'Pearson':>10} {'r^2':>8} {'p(rho)':>12} {'p(tau)':>12} {'exact':>7}")
print(f"  {'-'*110}")

bw = J["bw_median"].values
corr_results = []

for mc, ms, direction in metrics:
    vals = J[mc].values
    ok = np.isfinite(bw) & np.isfinite(vals)
    nok = ok.sum()
    if nok < 3: print(f"  {ms:<12} {direction:>4}  n<3"); continue
    if np.std(vals[ok]) < 1e-12: print(f"  {ms:<12} {direction:>4}  constant"); continue

    if direction == "benefit":
        rho, p_rho = spearmanr(vals[ok], bw[ok])
        tau, p_tau = kendalltau(vals[ok], bw[ok])
    else:
        rho, p_rho = spearmanr(vals[ok], -bw[ok])
        tau, p_tau = kendalltau(vals[ok], -bw[ok])

    pr, _ = pearsonr(vals[ok], bw[ok])
    rho2, pr2 = rho**2, pr**2
    n_exact = int(np.sum(J[f"{mc}_rank"].values[ok] == J["bw_rank"].values[ok]))
    sig = "***" if p_rho < 0.001 else "**" if p_rho < 0.01 else "*" if p_rho < 0.05 else ""

    print(f"  {ms:<12} {direction:>4} {rho:>+10.4f} {rho2:>8.4f} {tau:>+10.4f} "
          f"{pr:>+10.4f} {pr2:>8.4f} {p_rho:>12.2e} {p_tau:>12.2e} {n_exact:>3}/{nok} {sig}")
    corr_results.append({"metric": ms, "col": mc, "dir": direction,
                         "rho": rho, "rho2": rho2, "tau": tau,
                         "pr": pr, "pr2": pr2, "p_rho": p_rho,
                         "n_exact": n_exact, "n": nok})

# ── Summary ──
if corr_results:
    by_rho2 = sorted(corr_results, key=lambda x: x["rho2"], reverse=True)
    by_pr2  = sorted(corr_results, key=lambda x: x["pr2"],  reverse=True)
    by_exact = sorted(corr_results, key=lambda x: x["n_exact"], reverse=True)

    print(f"\n  === SUMMARY ===")
    print(f"  Best Spearman rho^2:  {by_rho2[0]['metric']:<12} rho={by_rho2[0]['rho']:+.4f}  rho^2={by_rho2[0]['rho2']:.4f}")
    print(f"  Best Pearson  r^2:    {by_pr2[0]['metric']:<12}    r={by_pr2[0]['pr']:+.4f}    r^2={by_pr2[0]['pr2']:.4f}")
    print(f"  Most exact matches:   {by_exact[0]['metric']:<12} exact={by_exact[0]['n_exact']}/{by_exact[0]['n']}")

    print(f"\n  Top-5 by Spearman rho^2:")
    print(f"  {'#':>3} {'Metric':<12} {'dir':>7} {'rho':>8} {'rho^2':>8} {'Pearson r^2':>12} {'exact':>7} {'p(rho)':>12}")
    print(f"  {'-'*75}")
    for i, r in enumerate(by_rho2[:5], 1):
        sig = "***" if r["p_rho"] < 0.001 else "**" if r["p_rho"] < 0.01 else "*" if r["p_rho"] < 0.05 else ""
        print(f"  {i:>3} {r['metric']:<12} {r['dir']:>7} {r['rho']:>+8.4f} {r['rho2']:>8.4f} {r['pr2']:>12.4f} "
              f"{r['n_exact']:>3}/{r['n']} {r['p_rho']:>12.2e} {sig}")

# ══════════════════════════════════════════════════════════════
#  3) PLOTS
# ══════════════════════════════════════════════════════════════

if not args.no_plot:
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        from matplotlib.patches import Patch

        LOOP_COLOR = {"klon_first": "#e67e22", "klev_first": "#2980b9"}
        LOOP_MARKER = {"klon_first": "o", "klev_first": "s"}

        nm = len(metrics)
        ncols = min(4, nm)
        nrows = (nm + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.5*ncols, 4*nrows), squeeze=False)
        suffix = "blocked" if args.blocked else "unblocked"
        fig.suptitle(f"Metric vs BW - {args.target} [{suffix}]", fontsize=14)

        for mi, (mc, ms, _) in enumerate(metrics):
            ax = axes[mi // ncols, mi % ncols]
            for lo in ["klon_first", "klev_first"]:
                sub = J[J["loop_order"] == lo]
                if sub.empty: continue
                ax.scatter(sub[mc], sub["bw_median"],
                           c=LOOP_COLOR.get(lo, "gray"), marker=LOOP_MARKER.get(lo, "o"),
                           s=55, alpha=0.85, edgecolors="black", linewidths=0.4,
                           label=lo.replace("_", "-"))
                for _, r in sub.iterrows():
                    if has_block:
                        lbl = f"B{int(r['blocking'])}"
                    else:
                        lbl = f"V{int(r['variant'])}"
                    ax.annotate(lbl, (r[mc], r["bw_median"]),
                                fontsize=5, alpha=0.7, xytext=(3, 3), textcoords="offset points")

            ri = [x for x in corr_results if x["col"] == mc]
            if ri:
                r = ri[0]
                ax.text(0.04, 0.04,
                        f"rho={r['rho']:+.3f} rho2={r['rho2']:.3f}\n"
                        f"  r={r['pr']:+.3f}   r2={r['pr2']:.3f}\n"
                        f"exact={r['n_exact']}/{r['n']}",
                        transform=ax.transAxes, fontsize=7, va="bottom", family="monospace",
                        bbox=dict(fc="lightyellow", alpha=0.9, pad=2))

            ax.set_xlabel(ms); ax.set_ylabel("BW [TB/s]" if mi % ncols == 0 else "")
            ax.set_title(ms, fontsize=11); ax.grid(alpha=0.2)
            vals = J[mc].dropna()
            if len(vals) > 0 and vals.max() / max(vals.min(), 1e-12) > 10:
                ax.set_xscale("log")
            if mi == 0: ax.legend(fontsize=7)

        for mi in range(nm, nrows * ncols):
            axes[mi // ncols, mi % ncols].set_visible(False)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fname = f"metric_vs_bw_{args.target}_{suffix}"
        fig.savefig(f"{fname}.png", dpi=180, bbox_inches="tight")
        fig.savefig(f"{fname}.pdf", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"\n  Scatter saved: {fname}.png/pdf")

        # ── Quality bar chart ──
        if corr_results:
            fig2, (ax_rho, ax_pr) = plt.subplots(2, 1, figsize=(max(10, 0.6*nm), 7))
            names = [r["metric"] for r in corr_results]
            rho2s = [r["rho2"] for r in corr_results]
            pr2s  = [r["pr2"] for r in corr_results]

            def color_bar(vals):
                return ["#2ecc71" if v > 0.5 else "#f39c12" if v > 0.2 else "#e74c3c" for v in vals]

            bars = ax_rho.bar(range(len(names)), rho2s, color=color_bar(rho2s),
                              edgecolor="black", linewidth=0.5)
            ax_rho.set_xticks(range(len(names)))
            ax_rho.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
            ax_rho.set_ylabel("Spearman rho^2"); ax_rho.set_ylim(0, 1.1)
            ax_rho.set_title(f"Rank Correlation (Spearman) - {args.target} [{suffix}]")
            ax_rho.axhline(y=1.0, color="gray", linestyle=":", linewidth=0.5)
            ax_rho.grid(axis="y", alpha=0.3)
            for i, (b, r) in enumerate(zip(bars, corr_results)):
                ax_rho.text(i, b.get_height() + 0.02,
                            f"{r['rho']:+.2f}", ha="center", va="bottom", fontsize=6)

            bars2 = ax_pr.bar(range(len(names)), pr2s, color=color_bar(pr2s),
                              edgecolor="black", linewidth=0.5)
            ax_pr.set_xticks(range(len(names)))
            ax_pr.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
            ax_pr.set_ylabel("Pearson r^2"); ax_pr.set_ylim(0, 1.1)
            ax_pr.set_title(f"Linear Correlation (Pearson) - {args.target} [{suffix}]")
            ax_pr.axhline(y=1.0, color="gray", linestyle=":", linewidth=0.5)
            ax_pr.grid(axis="y", alpha=0.3)
            for i, (b, r) in enumerate(zip(bars2, corr_results)):
                ax_pr.text(i, b.get_height() + 0.02,
                            f"{r['pr']:+.2f}", ha="center", va="bottom", fontsize=6)

            ax_pr.legend(handles=[
                Patch(facecolor="#2ecc71", edgecolor="black", label="r^2 > 0.5"),
                Patch(facecolor="#f39c12", edgecolor="black", label="0.2 < r^2 < 0.5"),
                Patch(facecolor="#e74c3c", edgecolor="black", label="r^2 < 0.2"),
            ], fontsize=7, loc="upper right")

            fig2.tight_layout()
            fname2 = f"metric_quality_{args.target}_{suffix}"
            fig2.savefig(f"{fname2}.png", dpi=180, bbox_inches="tight")
            fig2.savefig(f"{fname2}.pdf", dpi=180, bbox_inches="tight")
            plt.close(fig2)
            print(f"  Quality saved: {fname2}.png/pdf")

    except ImportError:
        print("\n  [WARN] matplotlib not available")

print("\nDone.")