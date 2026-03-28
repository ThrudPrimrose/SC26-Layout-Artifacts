#!/usr/bin/env python3
"""
gen_metrics_table.py

- LaTeX table with mu/delta/sigma for both nlev=90 and nlev=96
- Inter-metric correlation (Pearson, Spearman, Kendall) matrix
- Per-compute-variant rank analysis: do metrics agree on ordering?
"""
import pandas as pd
import numpy as np
from scipy import stats
import argparse, sys, itertools

parser = argparse.ArgumentParser()
parser.add_argument('--metrics', default='metrics.csv')
parser.add_argument('--out-table', default='metrics_table.tex')
parser.add_argument('--out-corr', default='metrics_corr.tex')
args = parser.parse_args()

cost = pd.read_csv(args.metrics)

# ---- Debug ----
print("=== metrics.csv summary ===")
for c in cost.columns:
    uvals = cost[c].unique()
    if len(uvals) < 20:
        print(f"  {c}: {sorted(uvals.tolist())}")
    else:
        print(f"  {c}: {len(uvals)} unique, range [{cost[c].min():.4f}, {cost[c].max():.4f}]")
print()

COMPUTE_VARIANTS = [
    ("cpu_scalar", 64,   1, "CPU Scalar"),
    ("cpu_avx512", 64,   8, "CPU AVX-512"),
    ("gpu_scalar", 128,  1, "GPU Scalar"),
    ("gpu_warp32", 128, 32, "GPU Warp32"),
    ("gpu_wave64", 128, 64, "GPU Wave64"),
]

LAYOUT_COMBOS = [
    ("klon_first", 1, "uniform",     "V1 / Uniform"),
    ("klon_first", 1, "normal_var1", "V1 / Normal"),
    ("klev_first", 4, "uniform",     "V4 / Uniform"),
    ("klev_first", 4, "normal_var1", "V4 / Normal"),
]

NLEVS = sorted(cost["nlev"].unique())
metric_cols = [c for c in ["mu", "delta", "sigma"] if c in cost.columns]
print(f"Metrics: {metric_cols},  nlev values: {NLEVS}\n")

# ---- Build full data frame ----
rows = []
for nlev in NLEVS:
    for cv_target, cv_block, cv_width, cv_label in COMPUTE_VARIANTS:
        for lay, var, dist, combo_label in LAYOUT_COMBOS:
            mask = (
                (cost["variant"] == var) &
                (cost["cell_dist"] == dist) &
                (cost["target"] == cv_target) &
                (cost["nlev"] == nlev) &
                (cost["loop_order"] == lay) &
                (cost["block_bytes"] == cv_block) &
                (cost["vector_width"] == cv_width)
            )
            row = cost[mask]
            if len(row) != 1:
                continue
            r = row.iloc[0]
            entry = {
                "nlev": nlev,
                "compute": cv_label,
                "layout": combo_label,
                "loop_order": lay,
                "variant": var,
                "dist": dist,
            }
            for m in metric_cols:
                entry[m] = r[m]
            rows.append(entry)

df = pd.DataFrame(rows)
print(f"Matched {len(df)} rows total ({len(df)//len(NLEVS)} per nlev).\n")

if len(df) == 0:
    print("ERROR: No data matched."); sys.exit(1)

# ═══════════════════════════════════════════════════════════════
#  TABLE 1: Metrics per compute variant × layout, both nlev
# ═══════════════════════════════════════════════════════════════

metric_tex = {"mu": r"$\bar\mu$", "delta": r"$\bar\Delta$", "sigma": r"$\bar\sigma$"}

# Build column headers: for each metric, one sub-column per nlev
# Format: Compute | Layout | mu_90 | mu_96 | delta_90 | delta_96 | sigma_90 | sigma_96
sub_headers = []
for m in metric_cols:
    for nlev in NLEVS:
        sub_headers.append((m, nlev))

n_data_cols = len(sub_headers)
col_spec = "ll" + "r" * n_data_cols

lines = []
lines.append(r"\begin{table}[t]")
lines.append(r"\centering")
lines.append(r"\caption{Cost model metrics for \texttt{z\_v\_grad\_w} across two problem sizes.}")
lines.append(r"\label{tab:cost-metrics}")
lines.append(r"\small")
lines.append(f"\\begin{{tabular}}{{{col_spec}}}")
lines.append(r"\toprule")

# Two-level header
# Top: metric names spanning 2 cols each
mc_line = "& "
for mi, m in enumerate(metric_cols):
    mc_line += rf" & \multicolumn{{2}}{{c}}{{{metric_tex[m]}}}"
mc_line += r" \\"
lines.append(mc_line)

# Bottom: nlev values
nlev_line = "Compute & Layout"
for m in metric_cols:
    for nlev in NLEVS:
        nlev_line += f" & $n_{{\\ell}}$={nlev}"
nlev_line += r" \\"
# cmidrules
cmi = []
col_idx = 3  # first data col
for m in metric_cols:
    cmi.append(f"\\cmidrule(lr){{{col_idx}-{col_idx+1}}}")
    col_idx += 2
lines.append(" ".join(cmi))
lines.append(nlev_line)
lines.append(r"\midrule")

prev_compute = None
for cv_target, cv_block, cv_width, cv_label in COMPUTE_VARIANTS:
    if prev_compute is not None:
        lines.append(r"\addlinespace")
    prev_compute = cv_label

    for lay, var, dist, combo_label in LAYOUT_COMBOS:
        vals = []
        for m in metric_cols:
            for nlev in NLEVS:
                sub = df[(df["compute"] == cv_label) &
                         (df["layout"] == combo_label) &
                         (df["nlev"] == nlev)]
                if len(sub) == 1:
                    v = sub.iloc[0][m]
                    if v < 100:
                        vals.append(f"{v:.2f}")
                    elif v < 10000:
                        vals.append(f"{v:.1f}")
                    else:
                        vals.append(f"{v:.0f}")
                else:
                    vals.append("--")
        vals_str = " & ".join(vals)
        lines.append(f"  {cv_label} & {combo_label} & {vals_str} \\\\")

lines.append(r"\bottomrule")
lines.append(r"\end{tabular}")
lines.append(r"\end{table}")

tex_table = "\n".join(lines)
with open(args.out_table, 'w') as f:
    f.write(tex_table)
print(f"Wrote: {args.out_table}")
print("\n--- Table Preview ---")
print(tex_table)

# ═══════════════════════════════════════════════════════════════
#  CORRELATION ANALYSIS
# ═══════════════════════════════════════════════════════════════

print("\n\n" + "=" * 70)
print("  CORRELATION ANALYSIS")
print("=" * 70)

# --- A) Inter-metric correlations (across all rows) ---
print("\n--- A) Inter-Metric Correlations (all rows pooled) ---\n")

corr_types = [
    ("Pearson r",  lambda x, y: stats.pearsonr(x, y)),
    ("Spearman ρ", lambda x, y: stats.spearmanr(x, y)),
    ("Kendall τ",  lambda x, y: stats.kendalltau(x, y)),
]

pairs = list(itertools.combinations(metric_cols, 2))

for cname, cfunc in corr_types:
    print(f"  {cname}:")
    for m1, m2 in pairs:
        vals1 = df[m1].values
        vals2 = df[m2].values
        r, p = cfunc(vals1, vals2)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"    {metric_tex[m1]:>15s} vs {metric_tex[m2]:>15s}:  {r:+.4f}  (p={p:.2e}) {sig}")
    print()

# --- B) Per-compute-variant: rank correlation of layouts ---
print("--- B) Per-Compute-Variant Layout Rankings ---\n")
print("  For each compute variant, rank the 4 layout combos by each metric.")
print("  Then compare rankings across metrics via Spearman and Kendall.\n")

for nlev in NLEVS:
    print(f"  nlev = {nlev}:")
    sub_nlev = df[df["nlev"] == nlev]

    for cv_target, cv_block, cv_width, cv_label in COMPUTE_VARIANTS:
        sub = sub_nlev[sub_nlev["compute"] == cv_label].copy()
        if len(sub) < 3:
            continue

        # Rank by each metric (lower = better for all three)
        rankings = {}
        for m in metric_cols:
            rankings[m] = sub[m].rank().values

        print(f"\n    {cv_label}:")
        # Print the ranking table
        print(f"      {'Layout':<20s}", end="")
        for m in metric_cols:
            print(f"  {m:>8s} (rank)", end="")
        print()
        for idx, (_, r) in enumerate(sub.iterrows()):
            print(f"      {r['layout']:<20s}", end="")
            for m in metric_cols:
                rv = rankings[m][idx]
                print(f"  {r[m]:>8.1f} ({rv:.0f})  ", end="")
            print()

        # Pairwise rank correlations
        if len(metric_cols) >= 2:
            print(f"      Rank correlations:")
            for m1, m2 in pairs:
                rho, p_rho = stats.spearmanr(rankings[m1], rankings[m2])
                tau, p_tau = stats.kendalltau(rankings[m1], rankings[m2])
                print(f"        {m1} vs {m2}:  Spearman={rho:+.3f} (p={p_rho:.3f})  Kendall={tau:+.3f} (p={p_tau:.3f})")

    print()

# --- C) Cross-nlev stability: does ranking change between nlev=90 and nlev=96? ---
if len(NLEVS) >= 2:
    print("--- C) Cross-nlev Ranking Stability ---\n")
    print("  Spearman ρ between metric rankings at nlev=90 vs nlev=96:\n")

    for cv_target, cv_block, cv_width, cv_label in COMPUTE_VARIANTS:
        print(f"  {cv_label}:")
        for m in metric_cols:
            vals_per_nlev = []
            for nlev in NLEVS[:2]:
                sub = df[(df["nlev"] == nlev) & (df["compute"] == cv_label)]
                sub = sub.sort_values("layout")  # align rows
                vals_per_nlev.append(sub[m].values)

            if len(vals_per_nlev[0]) == len(vals_per_nlev[1]) and len(vals_per_nlev[0]) >= 3:
                rho, p = stats.spearmanr(vals_per_nlev[0], vals_per_nlev[1])
                print(f"    {m:>8s}:  ρ={rho:+.4f}  (p={p:.3f})")
            else:
                print(f"    {m:>8s}:  insufficient data")
        print()

# ═══════════════════════════════════════════════════════════════
#  LATEX: Correlation summary table
# ═══════════════════════════════════════════════════════════════

clines = []
clines.append(r"\begin{table}[t]")
clines.append(r"\centering")
clines.append(r"\caption{Inter-metric correlations (pooled across all compute variants and layouts).}")
clines.append(r"\label{tab:metric-corr}")
clines.append(r"\small")
clines.append(r"\begin{tabular}{lrrr}")
clines.append(r"\toprule")
clines.append(r"Metric Pair & Pearson $r$ & Spearman $\rho$ & Kendall $\tau$ \\")
clines.append(r"\midrule")

for m1, m2 in pairs:
    v1, v2 = df[m1].values, df[m2].values
    pr, _ = stats.pearsonr(v1, v2)
    sr, _ = stats.spearmanr(v1, v2)
    kt, _ = stats.kendalltau(v1, v2)
    label = f"{metric_tex[m1]} vs {metric_tex[m2]}"
    clines.append(f"  {label} & {pr:+.3f} & {sr:+.3f} & {kt:+.3f} \\\\")

clines.append(r"\bottomrule")
clines.append(r"\end{tabular}")
clines.append(r"\end{table}")

tex_corr = "\n".join(clines)
with open(args.out_corr, 'w') as f:
    f.write(tex_corr)
print(f"\nWrote: {args.out_corr}")
print("\n--- Correlation Table Preview ---")
print(tex_corr)