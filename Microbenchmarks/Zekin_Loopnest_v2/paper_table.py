#!/usr/bin/env python3
"""
paper_table.py – Compact ranked table for SC26 paper.

Joins cost-model metrics with measured bandwidth, picks the best schedule
per (layout-group, distribution), and emits a terminal + LaTeX table.

Two metric modes (auto-selected or overridden):
  --metric muDnuma   →  μ·Δ_NUMA   (default for cpu_* targets)
  --metric mu        →  μ           (default for gpu_* targets)

Usage:
  python paper_table.py --target cpu_scalar --runtime z_v_grad_w_cpu_daint.csv
  python paper_table.py --target gpu_scalar --runtime z_v_grad_w_gpu_beverin.csv --gpu-config 1x1_32x16 --metric mu
"""

import argparse, sys, textwrap, warnings
import pandas as pd, numpy as np
from scipy.stats import spearmanr, pearsonr

warnings.filterwarnings("ignore", message="An input array is constant")

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser()
p.add_argument("--csv",        default="metrics.csv")
p.add_argument("--target",     required=True)
p.add_argument("--runtime",    required=True)
p.add_argument("--gpu-config", default="1x1_32x16")
p.add_argument("--nlev",       type=int, default=96)
p.add_argument("--metric",     choices=["muDnuma", "mu"], default=None,
               help="muDnuma (default CPU) or mu (default GPU)")
p.add_argument("--out",        default=None, help="LaTeX output file (default: table_<target>.tex)")
args = p.parse_args()

NPROMA = 81920
BYTES = (2*NPROMA*4 + 2*NPROMA*4 +
         args.nlev*NPROMA*8 + args.nlev*NPROMA*8 + NPROMA*8 +
         args.nlev*NPROMA*8 + args.nlev*NPROMA*8 +
         NPROMA*8 + NPROMA*8 + args.nlev*NPROMA*8)

# ── Resolve metric ───────────────────────────────────────────────────
if args.metric is None:
    args.metric = "mu" if args.target.startswith("gpu") else "muDnuma"

METRIC_COL   = {"muDnuma": "mu_delta_numa", "mu": "mu"}[args.metric]
METRIC_LABEL = {"muDnuma": r"$\mu\!\cdot\!\Delta_{\mathrm{NUMA}}$",
                "mu":      r"$\mu$"}[args.metric]
METRIC_SHORT = {"muDnuma": "μ·Δ_NUMA", "mu": "μ"}[args.metric]
METRIC_DIR   = "cost"  # lower is better for both

LAYOUT_GROUP = {1: "klon-first", 2: "klon-first", 3: "klev-first", 4: "klev-first"}

# ── Load cost metrics ────────────────────────────────────────────────
cost = pd.read_csv(args.csv)
cost = cost[cost["target"] == args.target].copy()
if cost.empty:
    avail = sorted(pd.read_csv(args.csv)["target"].unique())
    print(f"ERROR: target={args.target} not in {args.csv}\nAvailable: {avail}")
    sys.exit(1)
if METRIC_COL not in cost.columns:
    print(f"ERROR: column '{METRIC_COL}' not in {args.csv}")
    sys.exit(1)

# ── Load runtime ─────────────────────────────────────────────────────
rt = pd.read_csv(args.runtime)
rt["bw_tbs"] = BYTES / (rt["time_ms"] * 1e-3) / 1e12
if "config_label" in rt.columns and args.gpu_config:
    rt = rt[rt["config_label"] == args.gpu_config]
if "nlev" in rt.columns:
    rt = rt[rt["nlev"] == args.nlev]

group_cols = ["variant", "cell_dist"]
if "parallelization" in rt.columns:
    group_cols.append("parallelization")

bw_agg = rt.groupby(group_cols)["bw_tbs"].agg(
    bw_median="median", bw_std="std", bw_n="count"
).reset_index()

# Best schedule per (variant, dist)
if "parallelization" in bw_agg.columns:
    idx = bw_agg.groupby(["variant", "cell_dist"])["bw_median"].idxmax()
    bw_best = bw_agg.loc[idx].reset_index(drop=True)
else:
    bw_best = bw_agg.copy()

# ── Join ─────────────────────────────────────────────────────────────
rows = []
for _, rb in bw_best.iterrows():
    var, dist = int(rb["variant"]), rb["cell_dist"]
    mask = (cost["variant"] == var) & (cost["cell_dist"] == dist)
    crows = cost[mask]
    if crows.empty:
        continue
    if len(crows) > 1:
        expected = "klon_first" if var <= 2 else "klev_first"
        f = crows[crows["loop_order"] == expected]
        if not f.empty:
            crows = f
    c = crows.iloc[0]
    row = {
        "variant":    var,
        "dist":       dist,
        "layout":     LAYOUT_GROUP[var],
        "loop_order": c["loop_order"],
        "bw":         rb["bw_median"],
        "bw_std":     rb["bw_std"],
        "metric_val": c[METRIC_COL],
    }
    if "parallelization" in rb.index:
        row["schedule"] = rb["parallelization"]
    rows.append(row)

if not rows:
    print("ERROR: no rows after join")
    sys.exit(1)

J = pd.DataFrame(rows)

# ── Assign schedule labels (S1, S2, …) ──────────────────────────────
sched_map = {}
if "schedule" in J.columns:
    for s in sorted(J["schedule"].unique()):
        sched_map[s] = f"S{len(sched_map)+1}"
    J["sched_label"] = J["schedule"].map(sched_map)
else:
    J["sched_label"] = "S1"
    sched_map["default"] = "S1"

# ── Rank ─────────────────────────────────────────────────────────────
J = J.sort_values("bw", ascending=False).reset_index(drop=True)
J["bw_rank"]     = J["bw"].rank(ascending=False).astype(int)
J["metric_rank"] = J["metric_val"].rank(ascending=True).astype(int)  # cost: lower=better=rank 1

n = len(J)

# ── Correlation summary ──────────────────────────────────────────────
vals, bws = J["metric_val"].values, J["bw"].values
rho, p_rho = spearmanr(vals, -bws)  # cost metric vs -BW
pr, _      = pearsonr(vals, bws)
n_exact    = int(np.sum(J["metric_rank"].values == J["bw_rank"].values))

# ── Pretty distribution names ────────────────────────────────────────
DIST_PRETTY = {
    "sequential":  "Seq",
    "normal_var1": "N(1)",
    "normal_var4": "N(4)",
    "uniform":     "Unif",
}

# ═════════════════════════════════════════════════════════════════════
#  TERMINAL TABLE
# ═════════════════════════════════════════════════════════════════════

print(f"\n  Target: {args.target}   Metric: {METRIC_SHORT}   n={n}")
print(f"  Spearman ρ={rho:+.4f}  (ρ²={rho**2:.4f}, p={p_rho:.2e})")
print(f"  Pearson  r={pr:+.4f}   (r²={pr**2:.4f})")
print(f"  Exact rank matches: {n_exact}/{n}\n")

if sched_map:
    print("  Schedule legend:")
    for orig, label in sched_map.items():
        print(f"    {label} = {orig}")
    print()

hdr = f"  {'#':>2}  {'Layout':<11} {'Dist':<6} {'Sched':<5}  {METRIC_SHORT:>10} {'(rk)':>4}   {'BW [TB/s]':>10} {'(rk)':>4}"
sep = f"  {'-'*len(hdr)}"
print(hdr)
print(sep)

for _, r in J.iterrows():
    match = "◆" if r["metric_rank"] == r["bw_rank"] else " "
    print(f"  {r['bw_rank']:>2}  {r['layout']:<11} "
          f"{DIST_PRETTY.get(r['dist'], r['dist']):<6} "
          f"{r['sched_label']:<5}  "
          f"{r['metric_val']:>10.1f} ({r['metric_rank']:>2})   "
          f"{r['bw']:>10.4f} ({r['bw_rank']:>2}) {match}")

print(sep)
print(f"  ◆ = exact rank match\n")

# ═════════════════════════════════════════════════════════════════════
#  LATEX TABLE
# ═════════════════════════════════════════════════════════════════════

out_path = args.out or f"table_{args.target}.tex"

# Group rows by layout for midrule placement
prev_layout = None
body_lines = []
for _, r in J.iterrows():
    if prev_layout is not None and r["layout"] != prev_layout:
        body_lines.append(r"    \midrule")
    prev_layout = r["layout"]

    rk_m = int(r["metric_rank"])
    rk_b = int(r["bw_rank"])

    # Bold if metric rank == bw rank (exact prediction)
    m_str = f"{r['metric_val']:.1f}"
    b_str = f"{r['bw']:.4f}"
    if rk_m == rk_b:
        m_str = r"\textbf{" + m_str + "}"
        b_str = r"\textbf{" + b_str + "}"

    body_lines.append(
        f"    {int(r['bw_rank'])} & {r['layout']} & "
        f"{DIST_PRETTY.get(r['dist'], r['dist'])} & "
        f"{r['sched_label']} & "
        f"{m_str} & {rk_m} & "
        f"{b_str} & {rk_b} \\\\"
    )

body = "\n".join(body_lines)

caption_metric = METRIC_LABEL
tex = textwrap.dedent(rf"""
% Auto-generated by paper_table.py
% Target: {args.target}  |  Metric: {METRIC_SHORT}
% Spearman rho={rho:+.4f} (rho^2={rho**2:.4f})  Pearson r={pr:+.4f} (r^2={pr**2:.4f})
% Exact rank matches: {n_exact}/{n}
\begin{{table}}[t]
  \centering
  \caption{{Layout ranking by {caption_metric} vs.\ measured bandwidth
           ({args.target.replace("_", r"\_")}).
           Spearman $\rho^2\!=\!{rho**2:.2f}$, Pearson $r^2\!=\!{pr**2:.2f}$.
           Bold = exact rank match.}}
  \label{{tab:ranking_{args.target}_{args.metric}}}
  \small
  \begin{{tabular}}{{rl l l r@{{\;}}c r@{{\;}}c}}
    \toprule
    \# & Layout & Dist & Sched
       & \multicolumn{{2}}{{c}}{{{caption_metric}}}
       & \multicolumn{{2}}{{c}}{{BW [TB/s]}} \\
    \midrule
{body}
    \bottomrule
  \end{{tabular}}
\end{{table}}
""").strip() + "\n"

with open(out_path, "w") as f:
    f.write(tex)
print(f"  LaTeX written to {out_path}")

# ── Schedule legend as LaTeX footnote suggestion ─────────────────────
if len(sched_map) > 1:
    legend = ", ".join(f"{v}\\,=\\,\\texttt{{{k.replace('_', '-')}}}" for k, v in sched_map.items())
    print(f"  Schedule legend for caption/footnote:\n    {legend}")

print("\nDone.")