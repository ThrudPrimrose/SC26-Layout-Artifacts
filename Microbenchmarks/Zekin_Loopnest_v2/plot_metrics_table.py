#!/usr/bin/env python3
"""
gen_metrics_table.py

- LaTeX table: metrics + median runtime per platform/layout combo
- Correlation (Spearman, Kendall, Pearson R²) of each metric vs runtime
- Tests all compute models per platform, reports which fits best
"""
import pandas as pd
import numpy as np
from scipy import stats
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--metrics', default='metrics.csv')
parser.add_argument('--cpu-amd', default='z_v_grad_w_cpu_beverin.csv')
parser.add_argument('--cpu-nv',  default='z_v_grad_w_cpu_daint.csv')
parser.add_argument('--gpu-amd', default='z_v_grad_w_gpu_beverin.csv')
parser.add_argument('--gpu-nv',  default='z_v_grad_w_gpu_daint.csv')
parser.add_argument('--nlev', type=int, default=96)
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

# ---- Constants (from main plot script) ----
NPROMA = 81920
NLEV = args.nlev
BYTES = (2*NPROMA*4 + 2*NPROMA*4 +
         NLEV*NPROMA*8 + NLEV*NPROMA*8 + NPROMA*8 +
         NLEV*NPROMA*8 + NLEV*NPROMA*8 + NPROMA*8 +
         NPROMA*8 + NLEV*NPROMA*8)

STREAM_PEAK = {
    "MI300A Zen Cores":   811.35e-3,
    "Grace Neoverse v2": 1607.62e-3,
    "MI300A":   3457.5e-3,
    "H200":    3720.48e-3,
}

# Platform definitions: (label, csv_path, filter_col, filter_val)
PLATFORMS = [
    ("MI300A Zen Cores",  args.cpu_amd, "parallelization", None),
    ("Grace Neoverse v2", args.cpu_nv,  "parallelization", None),
    ("MI300A",            args.gpu_amd, "config_label",    "1x1_32x16"),
    ("H200",              args.gpu_nv,  "config_label",    "1x1_32x16"),
]

# Natural platform -> compute model mapping
PLATFORM_COMPUTE = {
    "MI300A Zen Cores":  [("cpu_scalar", 64, 1), ("cpu_avx512", 64, 8)],
    "Grace Neoverse v2": [("cpu_scalar", 64, 1)],
    "MI300A":            [("gpu_scalar", 128, 1), ("gpu_wave64", 128, 64)],
    "H200":              [("gpu_scalar", 128, 1), ("gpu_warp32", 128, 32)],
}

LAYOUT_COMBOS = [
    ("klon_first", 1, "uniform",     "V1/Uni"),
    ("klon_first", 1, "normal_var1", "V1/Norm"),
    ("klev_first", 4, "uniform",     "V4/Uni"),
    ("klev_first", 4, "normal_var1", "V4/Norm"),
]

metric_cols = [c for c in ["mu", "delta", "sigma"] if c in cost.columns]
metric_tex = {"mu": r"$\bar\mu$", "delta": r"$\bar\Delta$", "sigma": r"$\bar\sigma$"}

# ---- Load runtimes and compute median bandwidth ----
def compute_bandwidth(df):
    return BYTES / (df["time_ms"] * 1e-3) / 1e12

def pick_best_schedule(df_full, V, dist, nlev):
    sub = df_full[(df_full["variant"] == V) &
                  (df_full["cell_dist"] == dist) &
                  (df_full["nlev"] == nlev)].copy()
    if sub.empty:
        return None, np.nan
    sub["bandwidth"] = compute_bandwidth(sub)
    medians = sub.groupby("parallelization")["bandwidth"].median()
    if medians.empty:
        return None, np.nan
    best_par = medians.idxmax()
    return best_par, medians[best_par]

def get_median_bw(df_full, fcol, fval, V, dist, nlev):
    if fcol == "parallelization":
        _, med = pick_best_schedule(df_full, V, dist, nlev)
        return med
    else:
        sub = df_full[(df_full[fcol] == fval) &
                      (df_full["variant"] == V) &
                      (df_full["cell_dist"] == dist) &
                      (df_full["nlev"] == nlev)].copy()
        if sub.empty:
            return np.nan
        sub["bandwidth"] = compute_bandwidth(sub)
        return sub["bandwidth"].median()

# Load runtime CSVs
rt_data = {}
for plat_label, csv_path, fcol, fval in PLATFORMS:
    try:
        rt_data[plat_label] = (pd.read_csv(csv_path), fcol, fval)
    except FileNotFoundError:
        print(f"  [WARN] {csv_path} not found, skipping {plat_label}")

# Build joined table: platform × layout → (median_bw, metrics per compute model)
rows = []
for plat_label, csv_path, fcol, fval in PLATFORMS:
    if plat_label not in rt_data:
        continue
    df_rt, rt_fcol, rt_fval = rt_data[plat_label]
    peak = STREAM_PEAK.get(plat_label, np.nan)

    for lay, var, dist, combo_label in LAYOUT_COMBOS:
        med_bw = get_median_bw(df_rt, rt_fcol, rt_fval, var, dist, NLEV)

        for cv_target, cv_block, cv_width in PLATFORM_COMPUTE[plat_label]:
            mask = (
                (cost["variant"] == var) &
                (cost["cell_dist"] == dist) &
                (cost["target"] == cv_target) &
                (cost["nlev"] == NLEV) &
                (cost["loop_order"] == lay) &
                (cost["block_bytes"] == cv_block) &
                (cost["vector_width"] == cv_width)
            )
            mrow = cost[mask]
            if len(mrow) != 1:
                continue
            r = mrow.iloc[0]
            entry = {
                "platform": plat_label,
                "layout": combo_label,
                "variant": var,
                "dist": dist,
                "compute_model": cv_target,
                "block_bytes": cv_block,
                "vector_width": cv_width,
                "median_bw": med_bw,
                "pct_peak": 100.0 * med_bw / peak if not np.isnan(med_bw) else np.nan,
            }
            for m in metric_cols:
                entry[m] = r[m]
            rows.append(entry)

df = pd.DataFrame(rows)
print(f"Joined {len(df)} rows (platform × layout × compute model).\n")

if len(df) == 0:
    print("ERROR: No data matched."); sys.exit(1)

# ═══════════════════════════════════════════════════════════════
#  TABLE 1: Metrics + Runtime per platform
# ═══════════════════════════════════════════════════════════════

# Pick the "natural" compute model per platform for the table
NATURAL = {
    "MI300A Zen Cores":  "cpu_scalar",
    "Grace Neoverse v2": "cpu_scalar",
    "MI300A":            "gpu_wave64",
    "H200":              "gpu_warp32",
}

df_nat = df[df.apply(lambda r: r["compute_model"] == NATURAL[r["platform"]], axis=1)].copy()

metric_headers = " & ".join([metric_tex[m] for m in metric_cols])
col_spec = "llr" + "r" * len(metric_cols) + "r"

lines = []
lines.append(r"\begin{table}[t]")
lines.append(r"\centering")
lines.append(r"\caption{Cost metrics and measured bandwidth for \texttt{z\_v\_grad\_w}"
             f" ($n_{{\\text{{lev}}}}={NLEV}$).}}")
lines.append(r"\label{tab:cost-metrics-runtime}")
lines.append(r"\small")
lines.append(f"\\begin{{tabular}}{{{col_spec}}}")
lines.append(r"\toprule")
lines.append(f"Platform & Layout & BW [TB/s] & {metric_headers} & \\% Peak \\\\")
lines.append(r"\midrule")

prev_plat = None
for _, r in df_nat.sort_values(["platform", "layout"]).iterrows():
    if prev_plat is not None and r["platform"] != prev_plat:
        lines.append(r"\addlinespace")
    prev_plat = r["platform"]

    vals = []
    for m in metric_cols:
        v = r[m]
        if v < 100:
            vals.append(f"{v:.2f}")
        elif v < 10000:
            vals.append(f"{v:.1f}")
        else:
            vals.append(f"{v:.0f}")
    vals_str = " & ".join(vals)
    bw_str = f"{r['median_bw']:.3f}" if not np.isnan(r['median_bw']) else "--"
    pct_str = f"{r['pct_peak']:.1f}" if not np.isnan(r['pct_peak']) else "--"
    lines.append(f"  {r['platform']} & {r['layout']} & {bw_str} & {vals_str} & {pct_str} \\\\")

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
#  CORRELATION: Each metric vs runtime
# ═══════════════════════════════════════════════════════════════

print("\n\n" + "=" * 70)
print("  METRIC vs RUNTIME CORRELATION")
print("=" * 70)

def safe_corr(func, x, y):
    """Return (corr, pval) or (nan, nan) if insufficient data."""
    if len(x) < 3:
        return np.nan, np.nan
    try:
        return func(x, y)
    except Exception:
        return np.nan, np.nan

def pearson_r2(x, y):
    r, p = stats.pearsonr(x, y)
    return r**2, p

corr_funcs = [
    ("Spearman ρ",  stats.spearmanr),
    ("Kendall τ",   stats.kendalltau),
    ("Pearson R²",  pearson_r2),
]

# ---- A) Per-platform, per-compute-model ----
print("\n--- A) Per-Platform Correlation (metric vs median BW) ---")
print("  Higher metric = worse layout → expect NEGATIVE correlation with BW\n")

corr_rows = []
for plat_label in rt_data:
    print(f"  {plat_label}:")
    for cv_target, cv_block, cv_width in PLATFORM_COMPUTE[plat_label]:
        sub = df[(df["platform"] == plat_label) &
                 (df["compute_model"] == cv_target) &
                 (df["block_bytes"] == cv_block) &
                 (df["vector_width"] == cv_width)].copy()
        sub = sub.dropna(subset=["median_bw"])
        if len(sub) < 3:
            print(f"    {cv_target} (w={cv_width}): insufficient data ({len(sub)} points)")
            continue

        print(f"    {cv_target} (block={cv_block}, w={cv_width}), n={len(sub)}:")
        bw = sub["median_bw"].values

        for m in metric_cols:
            mv = sub[m].values
            for cname, cfunc in corr_funcs:
                r, p = safe_corr(cfunc, mv, bw)
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
                is_r2 = "R²" in cname
                val_str = f"{r:.4f}" if not is_r2 else f"{r:.4f}"
                print(f"      {metric_tex[m]:>15s} vs BW  {cname}: {val_str:>8s}  (p={p:.3f}) {sig}")
                corr_rows.append({
                    "platform": plat_label,
                    "compute_model": cv_target,
                    "vector_width": cv_width,
                    "metric": m,
                    "corr_type": cname,
                    "value": r,
                    "p_value": p,
                })
        print()

# ---- B) Pooled across platforms (using natural compute model) ----
print("--- B) Pooled Correlation (natural compute model per platform) ---\n")

bw_all = df_nat["median_bw"].dropna().values
for m in metric_cols:
    mv = df_nat.loc[df_nat["median_bw"].notna(), m].values
    print(f"  {metric_tex[m]:>15s} vs BW (n={len(mv)}):")
    for cname, cfunc in corr_funcs:
        r, p = safe_corr(cfunc, mv, bw_all)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"    {cname}: {r:+.4f}  (p={p:.3f}) {sig}")
    print()

# ---- C) Rank agreement: do metrics rank layouts the same as runtime? ----
print("--- C) Rank Agreement per Platform ---\n")
print("  Layout ranking by runtime vs by each metric (1=best BW / lowest metric)\n")

for plat_label in rt_data:
    sub = df_nat[df_nat["platform"] == plat_label].sort_values("median_bw", ascending=False).copy()
    if len(sub) < 2:
        continue
    sub["bw_rank"] = sub["median_bw"].rank(ascending=False)  # 1 = highest BW = best
    print(f"  {plat_label}:")
    header = f"    {'Layout':<12s} {'BW':>8s} {'BW_rk':>6s}"
    for m in metric_cols:
        sub[f"{m}_rank"] = sub[m].rank()  # 1 = lowest metric = best
        header += f"  {m:>8s} {m+'_rk':>6s}"
    print(header)
    for _, r in sub.iterrows():
        line = f"    {r['layout']:<12s} {r['median_bw']:>8.4f} {r['bw_rank']:>6.0f}"
        for m in metric_cols:
            line += f"  {r[m]:>8.2f} {r[f'{m}_rank']:>6.0f}"
        print(line)
    print()

# ═══════════════════════════════════════════════════════════════
#  LATEX: Correlation summary table
# ═══════════════════════════════════════════════════════════════

clines = []
clines.append(r"\begin{table}[t]")
clines.append(r"\centering")
clines.append(r"\caption{Metric--runtime correlations per platform (natural compute model).}")
clines.append(r"\label{tab:metric-runtime-corr}")
clines.append(r"\small")
clines.append(r"\begin{tabular}{ll" + "r" * len(metric_cols) + "}")
clines.append(r"\toprule")
metric_hdr = " & ".join([metric_tex[m] for m in metric_cols])
clines.append(f"Platform & Correlation & {metric_hdr} \\\\")
clines.append(r"\midrule")

prev_plat = None
for plat_label in [p for p, *_ in PLATFORMS if p in rt_data]:
    if prev_plat is not None:
        clines.append(r"\addlinespace")
    prev_plat = plat_label

    cv_target = NATURAL[plat_label]
    sub = df[(df["platform"] == plat_label) & (df["compute_model"] == cv_target)].dropna(subset=["median_bw"])
    if len(sub) < 3:
        continue
    bw = sub["median_bw"].values

    for cname, cfunc in corr_funcs:
        vals = []
        for m in metric_cols:
            r, p = safe_corr(cfunc, sub[m].values, bw)
            sig = "^{***}" if p < 0.001 else "^{**}" if p < 0.01 else "^{*}" if p < 0.05 else ""
            if np.isnan(r):
                vals.append("--")
            else:
                vals.append(f"${r:+.3f}{sig}$")
        vals_str = " & ".join(vals)
        clines.append(f"  {plat_label} & {cname} & {vals_str} \\\\")

clines.append(r"\addlinespace")
clines.append(r"\midrule")
# Pooled row
bw_all = df_nat["median_bw"].dropna().values
for cname, cfunc in corr_funcs:
    vals = []
    for m in metric_cols:
        mv = df_nat.loc[df_nat["median_bw"].notna(), m].values
        r, p = safe_corr(cfunc, mv, bw_all)
        sig = "^{***}" if p < 0.001 else "^{**}" if p < 0.01 else "^{*}" if p < 0.05 else ""
        vals.append(f"${r:+.3f}{sig}$" if not np.isnan(r) else "--")
    vals_str = " & ".join(vals)
    clines.append(f"  Pooled & {cname} & {vals_str} \\\\")

clines.append(r"\bottomrule")
clines.append(r"\end{tabular}")
clines.append(r"\end{table}")

tex_corr = "\n".join(clines)
with open(args.out_corr, 'w') as f:
    f.write(tex_corr)
print(f"\nWrote: {args.out_corr}")
print("\n--- Correlation Table Preview ---")
print(tex_corr)