#!/usr/bin/env python3
"""
paper_table.py – One compact table per invocation.

Single --target determines the metric (cpu_* → μ·Δ_NUMA, gpu_* → μ).
Multiple --runtime label:file entries share that metric column.
Only V1 (klon-first) and V4 (klev-first) reported.

Usage:
  # Grace CPU
  python paper_table.py --csv metrics_daint.csv --target cpu_scalar \
    --runtime grace:z_v_grad_w_cpu_daint.csv --drop-dist sequential \
    --out table_grace.tex

  # AMD CPU
  python paper_table.py --csv metrics_beverin.csv --target cpu_scalar \
    --runtime amd:z_v_grad_w_cpu_beverin.csv --drop-dist sequential \
    --out table_amd.tex

  # GPU (two machines, same metric)
  python paper_table.py --csv metrics_beverin.csv --target gpu_scalar \
    --runtime nvidia:z_v_grad_w_gpu_daint.csv \
    --runtime amd:z_v_grad_w_gpu_beverin.csv \
    --drop-dist sequential --gpu-config 1x1_32x16 --out table_gpu.tex
"""

import argparse, sys, textwrap, warnings
import pandas as pd, numpy as np
from scipy.stats import spearmanr, pearsonr

warnings.filterwarnings("ignore", message="An input array is constant")

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser()
p.add_argument("--csv",        required=True, help="Cost-model CSV")
p.add_argument("--target",     required=True, help="Cost-model target row filter")
p.add_argument("--runtime",    action="append", required=True,
               help="label:file  (repeatable)")
p.add_argument("--gpu-config", default="1x1_32x16")
p.add_argument("--nlev",       type=int, default=96)
p.add_argument("--drop-dist",  action="append", default=[],
               help="Drop distribution(s), e.g. --drop-dist sequential")
p.add_argument("--out",        default=None)
args = p.parse_args()

NPROMA = 81920
BYTES = (2*NPROMA*4 + 2*NPROMA*4 +
         args.nlev*NPROMA*8 + args.nlev*NPROMA*8 + NPROMA*8 +
         args.nlev*NPROMA*8 + args.nlev*NPROMA*8 +
         NPROMA*8 + NPROMA*8 + args.nlev*NPROMA*8)

KEEP_VARIANTS = [1, 4]
LAYOUT_NAME   = {1: "klon-first", 4: "klev-first"}
DIST_PRETTY   = {"sequential": "Seq", "normal_var1": "N(1)",
                 "normal_var4": "N(4)", "uniform": "Unif"}
DIST_ORDER    = {"sequential": 0, "normal_var1": 1, "normal_var4": 2, "uniform": 3}

# ── Metric selection ─────────────────────────────────────────────────
if args.target.startswith("gpu"):
    MET_COL, MET_LATEX, MET_SHORT = "mu", r"$\mu$", "μ"
else:
    MET_COL, MET_LATEX, MET_SHORT = (
        "mu_delta_numa",
        r"$\mu\!\cdot\!\Delta_{\mathrm{NUMA}}$",
        "μ·Δ_NUMA")

out_path = args.out or f"table_{args.target}.tex"

# ── Parse --runtime entries ──────────────────────────────────────────
machines = []
for spec in args.runtime:
    parts = spec.split(":", 1)
    if len(parts) != 2:
        sys.exit(f"ERROR: --runtime expects label:file, got '{spec}'")
    machines.append({"label": parts[0], "file": parts[1]})

# ── Load cost model ──────────────────────────────────────────────────
cost = pd.read_csv(args.csv)
cost = cost[cost["target"] == args.target].copy()
if cost.empty:
    avail = sorted(pd.read_csv(args.csv)["target"].unique())
    sys.exit(f"ERROR: target={args.target} not in {args.csv}\nAvailable: {avail}")
if MET_COL not in cost.columns:
    sys.exit(f"ERROR: column '{MET_COL}' not in {args.csv}")

# ── Load runtime ─────────────────────────────────────────────────────
def load_bw(filepath):
    rt = pd.read_csv(filepath)
    rt["bw_tbs"] = BYTES / (rt["time_ms"] * 1e-3) / 1e12
    if "config_label" in rt.columns and args.target.startswith("gpu"):
        rt = rt[rt["config_label"] == args.gpu_config]
    if "nlev" in rt.columns:
        rt = rt[rt["nlev"] == args.nlev]
    rt = rt[rt["variant"].isin(KEEP_VARIANTS)]
    for d in args.drop_dist:
        rt = rt[rt["cell_dist"] != d]

    gcols = ["variant", "cell_dist"]
    has_sched = "parallelization" in rt.columns
    if has_sched:
        gcols.append("parallelization")

    agg = rt.groupby(gcols)["bw_tbs"].agg(
        bw_median="median", bw_std="std"
    ).reset_index()

    if not has_sched:
        agg["parallelization"] = ""

    # Keep only the best schedule per (variant, dist)
    idx = agg.groupby(["variant", "cell_dist"])["bw_median"].idxmax()
    return agg.loc[idx].reset_index(drop=True)

# ── Join cost + runtime ──────────────────────────────────────────────
canonical = set()       # (var, dist, sched)
all_scheds = set()
bw_data = {}            # label → {key: bw}
met_data = {}           # (var, dist) → metric_val  (same for all machines)

for m in machines:
    bw_all = load_bw(m["file"])
    bw_data[m["label"]] = {}

    for _, rb in bw_all.iterrows():
        var, dist, sched = int(rb["variant"]), rb["cell_dist"], rb["parallelization"]
        key = (var, dist, sched)
        canonical.add(key)
        all_scheds.add(sched)
        bw_data[m["label"]][key] = rb["bw_median"]

        if (var, dist) not in met_data:
            mask = (cost["variant"] == var) & (cost["cell_dist"] == dist)
            crows = cost[mask]
            if crows.empty:
                continue
            if len(crows) > 1:
                expected = "klon_first" if var <= 2 else "klev_first"
                f = crows[crows["loop_order"] == expected]
                if not f.empty:
                    crows = f
            met_data[(var, dist)] = crows.iloc[0][MET_COL]

# ── Schedule labels ──────────────────────────────────────────────────
sched_map = {s: f"S{i+1}" for i, s in enumerate(sorted(all_scheds - {""}))}
if "" in all_scheds:
    sched_map[""] = ""
has_sched = any(s != "" for s in sched_map)

# ── Assemble table ───────────────────────────────────────────────────
rows = []
for var, dist, sched in sorted(canonical,
        key=lambda x: (x[0], DIST_ORDER.get(x[1], 99), x[2])):
    row = {"layout": LAYOUT_NAME[var], "dist": dist,
           "sched": sched_map.get(sched, sched),
           "metric": met_data.get((var, dist), np.nan)}
    for m in machines:
        row[f"bw_{m['label']}"] = bw_data[m["label"]].get((var, dist, sched), np.nan)
    rows.append(row)

T = pd.DataFrame(rows)
n = len(T)

# Ranks
T["rk_met"] = T["metric"].rank(ascending=True, method="min").astype("Int64")
for m in machines:
    T[f"rk_{m['label']}"] = T[f"bw_{m['label']}"].rank(
        ascending=False, method="min").astype("Int64")

# ── Correlations ─────────────────────────────────────────────────────
print(f"\n  Target: {args.target}   Metric: {MET_SHORT}   n={n}")
if args.drop_dist:
    print(f"  Dropped: {', '.join(args.drop_dist)}")
print()

corr_info = []
for m in machines:
    lbl = m["label"]
    mv = T["metric"].values.astype(float)
    bv = T[f"bw_{lbl}"].values.astype(float)
    ok = np.isfinite(mv) & np.isfinite(bv)
    if ok.sum() < 3:
        continue
    rho, p_rho = spearmanr(mv[ok], -bv[ok])
    pr, _ = pearsonr(mv[ok], bv[ok])
    n_ex = int(np.sum(T["rk_met"].values[ok] == T[f"rk_{lbl}"].values[ok]))
    sig = "***" if p_rho < 0.001 else "**" if p_rho < 0.01 else "*" if p_rho < 0.05 else ""
    print(f"  {lbl:>8}  ρ={rho:+.4f} (ρ²={rho**2:.4f})  "
          f"r={pr:+.4f} (r²={pr**2:.4f})  exact={n_ex}/{ok.sum()} {sig}")
    corr_info.append({"label": lbl, "rho": rho, "rho2": rho**2, "r2": pr**2})

# ═════════════════════════════════════════════════════════════════════
#  TERMINAL TABLE
# ═════════════════════════════════════════════════════════════════════
print()
if has_sched:
    for orig, label in sorted(sched_map.items(), key=lambda x: x[1]):
        if orig:
            print(f"  {label} = {orig}")
    print()

hdr_parts = [f"{'Layout':<11}", f"{'Dist':<5}"]
if has_sched:
    hdr_parts.append(f"{'S':<3}")
hdr_parts.append(f"{MET_SHORT:>10}")
for m in machines:
    hdr_parts.append(f"{'BW_'+m['label']:>10}")
hdr = "  ".join(hdr_parts)
print(f"  {hdr}")
print(f"  {'─' * len(hdr)}")

# Sort by first machine BW descending
first = machines[0]["label"]
T = T.sort_values(f"bw_{first}", ascending=False, na_position="last").reset_index(drop=True)

prev_layout = None
for _, r in T.iterrows():
    if prev_layout is not None and r["layout"] != prev_layout:
        print(f"  {'─' * len(hdr)}")
    prev_layout = r["layout"]

    parts = [f"{r['layout']:<11}", f"{DIST_PRETTY.get(r['dist'], r['dist']):<5}"]
    if has_sched:
        parts.append(f"{r['sched']:<3}")
    v = r["metric"]
    parts.append(f"{v:>10.1f}" if pd.notna(v) else f"{'—':>10}")
    for m in machines:
        bv = r[f"bw_{m['label']}"]
        parts.append(f"{bv:>10.4f}" if pd.notna(bv) else f"{'—':>10}")
    print(f"  {'  '.join(parts)}")

print()

# ═════════════════════════════════════════════════════════════════════
#  LATEX TABLE
# ═════════════════════════════════════════════════════════════════════

n_mach = len(machines)
colspec = "l l"
if has_sched:
    colspec += " l"
colspec += " r" + " r" * n_mach

# Header
tex_parts = ["Layout & Dist"]
if has_sched:
    tex_parts.append("& S")
tex_parts.append(f"& {MET_LATEX}")
for m in machines:
    if n_mach == 1:
        tex_parts.append(r"& BW [TB/s]")
    else:
        tex_parts.append(rf"& \textsc{{{m['label']}}}")
tex_header = " ".join(tex_parts) + r" \\"

# Body
body_lines = []
prev_layout = None
for _, r in T.iterrows():
    if prev_layout is not None and r["layout"] != prev_layout:
        body_lines.append(r"    \midrule")
    prev_layout = r["layout"]

    dist_tex = DIST_PRETTY.get(r["dist"], r["dist"])
    parts = [f"    {r['layout']}", dist_tex]
    if has_sched:
        parts.append(r["sched"])
    v = r["metric"]
    parts.append(f"{v:.1f}" if pd.notna(v) else "---")
    for m in machines:
        bv = r[f"bw_{m['label']}"]
        parts.append(f"{bv:.3f}" if pd.notna(bv) else "---")
    body_lines.append(" & ".join(parts) + r" \\")

body = "\n".join(body_lines)

# Caption
cap_parts = []
for ci in corr_info:
    if n_mach == 1:
        cap_parts.append(
            rf"$\rho^2\!=\!{ci['rho2']:.2f}$, $r^2\!=\!{ci['r2']:.2f}$")
    else:
        cap_parts.append(
            rf"\textsc{{{ci['label']}}}: $\rho^2\!=\!{ci['rho2']:.2f}$, $r^2\!=\!{ci['r2']:.2f}$")
cap_corr = "; ".join(cap_parts)

sched_note = ""
if has_sched and len(sched_map) > 1:
    legend = ", ".join(f"{v}\\,=\\,\\texttt{{{k.replace('_', '-')}}}"
                       for k, v in sorted(sched_map.items(), key=lambda x: x[1]) if k)
    sched_note = f"\n           {legend}."

tex = textwrap.dedent(rf"""
% Auto-generated by paper_table.py  target={args.target}
\begin{{table}}[t]
  \centering
  \caption{{{MET_LATEX} vs.\ measured bandwidth ({args.target.replace("_", r"\_")}).
           {cap_corr}.{sched_note}}}
  \label{{tab:{args.target}}}
  \small
  \begin{{tabular}}{{{colspec}}}
    \toprule
    {tex_header}
    \midrule
{body}
    \bottomrule
  \end{{tabular}}
\end{{table}}
""").strip() + "\n"

with open(out_path, "w") as f:
    f.write(tex)
print(f"  LaTeX → {out_path}")
print("Done.")