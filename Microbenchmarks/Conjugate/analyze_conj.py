#!/usr/bin/env python3
"""
analyze_conj.py
Joins conjugate cost metrics with runtime. Computes rank correlation.
Usage:
    python analyze_conj.py --runtime conj_bench.csv --device CPU --mode out-of-place
"""
import argparse, sys, warnings
import pandas as pd, numpy as np
from scipy.stats import spearmanr, kendalltau, pearsonr
warnings.filterwarnings("ignore", message="An input array is constant")

parser = argparse.ArgumentParser()
parser.add_argument("--metrics", default="metrics_conj.csv")
parser.add_argument("--runtime", required=True)
parser.add_argument("--device", default=None)
parser.add_argument("--mode", default=None)
parser.add_argument("--no-plot", action="store_true")
args = parser.parse_args()

cost = pd.read_csv(args.metrics)
print(f"Cost metrics: {len(cost)} rows, layouts={sorted(cost['layout'].unique())}, modes={sorted(cost['mode'].unique())}")

ALL_METRICS = [("mu","mu"),("W","W"),("delta_min","D_min"),("delta_max","D_max"),
               ("sigma","sigma"),("delta_numa","D_numa"),("W_delta","W*D_min"),("span","span")]
metrics = [(c,s) for c,s in ALL_METRICS if c in cost.columns]

# Load runtime — just layout/run/ms/gbps
rt = pd.read_csv(args.runtime)
bw_agg = rt.groupby("layout")["gbps"].agg(
    bw_median="median", bw_std="std", bw_n="count").reset_index()

# Filter cost metrics by mode (if specified)
if args.mode and "mode" in cost.columns:
    cost = cost[cost["mode"] == args.mode]

bw_col = "gbps" if "gbps" in rt.columns else "bw_tbs"
gcols = [c for c in ["layout","mode","device"] if c in rt.columns]
bw_agg = rt.groupby(gcols)[bw_col].agg(bw_median="median",bw_std="std",bw_n="count").reset_index()

print(f"\nRuntime ({len(bw_agg)} configs):")
for _,r in bw_agg.iterrows():
    tag = " / ".join(str(r[c]) for c in gcols)
    print(f"  {tag:<35s} median={r['bw_median']:>8.2f}  std={r['bw_std']:>6.2f}")

rows = []
for _, rb in bw_agg.iterrows():
    layout = rb["layout"]
    crows = cost[cost["layout"] == layout]
    if crows.empty:
        print(f"  [WARN] No metric for {layout}")
        continue
    c = crows.iloc[0]
    row = {"layout":layout,"bw_median":rb["bw_median"],"bw_std":rb["bw_std"],"n_refs":c["n_refs"]}
    if "mode" in c.index: row["mode"]=c["mode"]
    for mc,_ in metrics: row[mc]=c[mc]
    rows.append(row)

if not rows: print("\nERROR: No rows joined"); sys.exit(1)
J = pd.DataFrame(rows).sort_values("bw_median",ascending=False).reset_index(drop=True)
J["bw_rank"] = J["bw_median"].rank(ascending=False).astype(int)
n = len(J)

for mc,_ in metrics:
    J[f"{mc}_rank"] = J[mc].rank(ascending=True).astype(int)

# ── Table ──
print(f"\n{'='*140}\n  RANKED TABLE (n={n})\n{'='*140}\n")
hdr = f"  {'#':>2} {'Layout':<12}"
if "mode" in J.columns: hdr += f" {'Mode':<14}"
hdr += f" {'refs':>4} {'BW':>9} {'+-':>6}"
for _,s in metrics: hdr += f"  {s:>9}(rk)"
print(hdr)
print(f"  {'-'*(len(hdr)+5)}")
for _,r in J.iterrows():
    line = f"  {r['bw_rank']:>2} {r['layout']:<12}"
    if "mode" in J.columns: line += f" {r['mode']:<14}"
    line += f" {int(r['n_refs']):>4} {r['bw_median']:>9.2f} {r['bw_std']:>6.2f}"
    for mc,_ in metrics:
        v,rk = r[mc], r[f"{mc}_rank"]
        m = "=" if rk==r["bw_rank"] else " "
        if abs(v)<1: line+=f"  {v:>8.4f}({rk:>2}{m})"
        elif abs(v)<1000: line+=f"  {v:>8.2f}({rk:>2}{m})"
        else: line+=f"  {v:>8.0f}({rk:>2}{m})"
    print(line)

# ── Correlation ──
print(f"\n{'='*100}\n  METRIC QUALITY (n={n})\n{'='*100}\n")
print(f"  {'Metric':<12} {'Spearman':>10} {'r^2':>8} {'Kendall':>10} {'Pearson':>10} {'r^2':>8} {'p(rho)':>12} {'exact':>7}")
print(f"  {'-'*85}")

bw = J["bw_median"].values
corr_results = []
for mc,ms in metrics:
    vals = J[mc].values
    ok = np.isfinite(bw) & np.isfinite(vals)
    nok = ok.sum()
    if nok < 3: print(f"  {ms:<12} n<3"); continue
    # Check for constant input
    if np.std(vals[ok]) < 1e-12:
        print(f"  {ms:<12} {'constant':>10}")
        continue
    rho,p_rho = spearmanr(vals[ok],-bw[ok])
    tau,p_tau = kendalltau(vals[ok],-bw[ok])
    pr,_ = pearsonr(vals[ok],bw[ok])
    rho2,pr2 = rho**2,pr**2
    n_exact = int(np.sum(J[f"{mc}_rank"].values[ok]==J["bw_rank"].values[ok]))
    sig = "***" if p_rho<0.001 else "**" if p_rho<0.01 else "*" if p_rho<0.05 else ""
    print(f"  {ms:<12} {rho:>+10.4f} {rho2:>8.4f} {tau:>+10.4f} {pr:>+10.4f} {pr2:>8.4f} {p_rho:>12.2e} {n_exact:>3}/{nok} {sig}")
    corr_results.append({"metric":ms,"col":mc,"rho":rho,"rho2":rho2,"tau":tau,"pr":pr,"n_exact":n_exact,"n":nok})

if corr_results:
    best = max(corr_results, key=lambda x: x["rho2"])
    print(f"\n  >>> Best: {best['metric']:<12} rho={best['rho']:+.4f}  rho^2={best['rho2']:.4f}  exact={best['n_exact']}/{best['n']}")

# ── Plots ──
if not args.no_plot and len(J)>=3:
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        COLORS={"AoS":"#e74c3c","SoA":"#3498db","AoSoA-8":"#2ecc71","AoSoA-16":"#f39c12"}
        MARKERS={"AoS":"o","SoA":"s","AoSoA-8":"D","AoSoA-16":"^"}
        nm=len(metrics); ncols=min(4,nm); nrows=(nm+ncols-1)//ncols
        dev=args.device or "all"; mod=args.mode or "all"
        fig,axes=plt.subplots(nrows,ncols,figsize=(4.5*ncols,4*nrows),squeeze=False)
        fig.suptitle(f"Conjugate: Metric vs BW ({dev}/{mod})",fontsize=14)
        for mi,(mc,ms) in enumerate(metrics):
            ax=axes[mi//ncols,mi%ncols]
            for layout in J["layout"].unique():
                sub=J[J["layout"]==layout]
                ax.scatter(sub[mc],sub["bw_median"],c=COLORS.get(layout,"gray"),
                           marker=MARKERS.get(layout,"o"),s=70,alpha=0.85,
                           edgecolors="black",linewidths=0.4,label=layout,zorder=3)
            ri=[x for x in corr_results if x["col"]==mc]
            if ri:
                r=ri[0]
                ax.text(0.04,0.04,f"rho={r['rho']:+.3f}\nrho2={r['rho2']:.3f}",
                        transform=ax.transAxes,fontsize=7,va="bottom",family="monospace",
                        bbox=dict(fc="lightyellow",alpha=0.9,pad=2))
            ax.set_xlabel(ms); ax.set_ylabel("BW [GB/s]" if mi%ncols==0 else "")
            ax.set_title(ms,fontsize=11); ax.grid(alpha=0.2)
            vals=J[mc].dropna()
            if len(vals)>0 and vals.max()/max(vals.min(),1e-12)>10: ax.set_xscale("log")
            if mi==0: ax.legend(fontsize=7)
        for mi in range(nm,nrows*ncols): axes[mi//ncols,mi%ncols].set_visible(False)
        fig.tight_layout(rect=[0,0,1,0.94])
        fname=f"conj_metric_vs_bw_{dev}_{mod}"
        fig.savefig(f"{fname}.png",dpi=180,bbox_inches="tight")
        plt.close(fig); print(f"\n  Saved: {fname}.png")
    except ImportError: print("\n  [WARN] matplotlib not available")
print("\nDone.")