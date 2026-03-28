#!/usr/bin/env python3
"""
join_bw_metrics.py

Joins cost-model metrics with measured bandwidth from runtime CSVs.
Computes rank correlation for all 5 metrics: μ, Δ, σ, Δ_max, μ·Δ_max.

Usage:
    python join_bw_metrics.py --platform beverin
    python join_bw_metrics.py --platform beverin --platform daint
"""
import argparse
import pandas as pd
import numpy as np
import sys

NPROMA = 81920
NLEV   = 96
BYTES  = (2*NPROMA*4 + 2*NPROMA*4 +
          NLEV*NPROMA*8 + NLEV*NPROMA*8 + NPROMA*8 +
          NLEV*NPROMA*8 + NLEV*NPROMA*8 +
          NPROMA*8 + NPROMA*8 + NLEV*NPROMA*8)

DIST_TO_COST = {
    "uniform": "uniform", "normal_var1": "normal1",
    "normal_var4": "normal4", "sequential": "sequential",
}
SCHED_TO_COST = {
    "omp_for": "omp_for", "omp_collapse2": "omp_collapse2",
    "opt_for": "omp_for", "opt_collapse_tile": "omp_collapse2",
}
RUNTIME_TO_COST_TARGET = {"cpu": "CPU_scalar", "gpu": "GPU_scalar"}

# All metrics to evaluate (col_name, display_symbol)
ALL_METRICS = [
    ("mu",           "μ"),
    ("delta",        "Δ"),
    ("sigma",        "σ"),
    ("delta_max",    "Δ_max"),
    ("mu_delta_max", "μ·Δ_max"),
]

PLATFORM_CFG = {
    "beverin": {
        "label": "MI300A Zen 4",
        "cpu_csv": "z_v_grad_w_cpu_beverin.csv",
        "gpu_csv": "z_v_grad_w_gpu_beverin.csv",
        "cost_csv": "results_full_beverin.csv",
        "gpu_filter_col": "config_label",
        "gpu_filter_val": "1x1_32x16",
    },
    "daint": {
        "label": "Grace / H200",
        "cpu_csv": "z_v_grad_w_cpu_daint.csv",
        "gpu_csv": "z_v_grad_w_gpu_daint.csv",
        "cost_csv": "results_full_daint.csv",
        "gpu_filter_col": "config_label",
        "gpu_filter_val": "1x1_32x16",
    },
}


def load_safe(path):
    try:
        df = pd.read_csv(path)
        return df if not df.empty else None
    except Exception:
        return None


def compute_bw(df):
    df = df.copy()
    df["bw_tbs"] = BYTES / (df["time_ms"] * 1e-3) / 1e12
    return df


def extract_runtime_bw(df, hw_type, gpu_filter_col=None, gpu_filter_val=None):
    if df is None:
        return pd.DataFrame()
    df = df.copy()
    if hw_type == "gpu" and gpu_filter_col and gpu_filter_val:
        df = df[df[gpu_filter_col] == gpu_filter_val]
    if "time_ms" not in df.columns:
        return pd.DataFrame()
    df = compute_bw(df)

    group_cols = (["variant", "cell_dist", "parallelization"] if hw_type == "cpu"
                  else (["variant", "cell_dist", "parallelization"]
                        if "parallelization" in df.columns
                        else ["variant", "cell_dist"]))

    rows = []
    for keys, grp in df.groupby(group_cols):
        if hw_type == "cpu":
            var, dist, sched = keys
        elif len(keys) == 3:
            var, dist, sched = keys
        else:
            var, dist = keys; sched = "gpu_kernel"

        rows.append({
            "variant_num": int(var), "variant": f"V{int(var)}",
            "dist_cost": DIST_TO_COST.get(dist, dist),
            "schedule_rt": sched,
            "schedule_cost": SCHED_TO_COST.get(sched, sched),
            "bw_median": grp["bw_tbs"].median(),
            "bw_std": grp["bw_tbs"].std(),
            "hw_type": hw_type,
        })
    return pd.DataFrame(rows)


def best_schedule_per_config(rt_df):
    if rt_df.empty:
        return rt_df
    idx = rt_df.groupby(["variant", "dist_cost", "hw_type"])["bw_median"].idxmax()
    return rt_df.loc[idx].reset_index(drop=True)


def join_with_cost(rt_best, cost_df, hw_type):
    if rt_best.empty or cost_df is None or cost_df.empty:
        return pd.DataFrame()

    cost_target = RUNTIME_TO_COST_TARGET[hw_type]
    metric_cols = [m[0] for m in ALL_METRICS]

    rows = []
    for _, r in rt_best[rt_best["hw_type"] == hw_type].iterrows():
        mask = ((cost_df["variant"] == r["variant"]) &
                (cost_df["dist"] == r["dist_cost"]) &
                (cost_df["schedule"] == r["schedule_cost"]) &
                (cost_df["target"] == cost_target))
        crows = cost_df[mask]
        if crows.empty:
            mask2 = ((cost_df["variant"] == r["variant"]) &
                     (cost_df["dist"] == r["dist_cost"]) &
                     (cost_df["target"] == cost_target))
            crows = cost_df[mask2]

        row = {
            "variant": r["variant"], "dist": r["dist_cost"],
            "schedule": r["schedule_rt"],
            "loop_order": crows.iloc[0]["loop_order"] if not crows.empty else "?",
            "hw_type": hw_type, "cost_target": cost_target,
            "bw_median": r["bw_median"], "bw_std": r["bw_std"],
        }
        for mc in metric_cols:
            row[mc] = crows.iloc[0][mc] if (not crows.empty and mc in crows.columns) else np.nan
        rows.append(row)

    return pd.DataFrame(rows)


def fmt(v, w=10, prec=1):
    """Format a float, or '---' if nan."""
    if np.isfinite(v):
        return f"{v:>{w}.{prec}f}"
    return f"{'---':>{w}}"


def print_table(joined, platform_label):
    if joined.empty:
        return
    print(f"\n{'='*155}")
    print(f"  {platform_label}")
    print(f"{'='*155}")

    for hw in ["cpu", "gpu"]:
        sub = joined[joined["hw_type"] == hw].copy()
        if sub.empty:
            continue

        # Primary ranking metric: σ for CPU, μ for GPU
        primary = "sigma" if hw == "cpu" else "mu"
        primary_sym = "σ" if hw == "cpu" else "μ"

        sub["bw_rank"] = sub["bw_median"].rank(ascending=False).astype(int)
        sub["metric_rank"] = sub[primary].rank(ascending=True).astype(int)
        sub = sub.sort_values("bw_rank")

        hw_label = "CPU" if hw == "cpu" else "GPU"
        target = sub.iloc[0]["cost_target"]
        print(f"\n  --- {hw_label} (metric target: {target}, primary ranking: {primary_sym}) ---")

        print(f"  {'#bw':>3} {'#'+primary_sym:>5} {'Variant':<8} {'Dist':<14} "
              f"{'Schedule':<18} {'Loop':<14} {'BW':>9} {'±':>6} "
              f"{'μ':>8} {'Δ':>10} {'σ':>12} {'Δ_max':>10} {'μ·Δ_max':>12}")
        print(f"  {'-'*145}")

        for _, r in sub.iterrows():
            match = "✓" if r["bw_rank"] == r["metric_rank"] else " "
            print(f"  {r['bw_rank']:>3} {r['metric_rank']:>5}{match} "
                  f"{r['variant']:<8} {r['dist']:<14} {r['schedule']:<18} "
                  f"{r['loop_order']:<14} {r['bw_median']:>9.4f} {r['bw_std']:>6.4f} "
                  f"{fmt(r['mu'], 8, 2)} {fmt(r['delta'], 10)} {fmt(r['sigma'], 12)} "
                  f"{fmt(r['delta_max'], 10)} {fmt(r['mu_delta_max'], 12)}")


def compute_correlations(combined):
    """Compute Spearman ρ, Rank R², Kendall τ for every metric."""
    try:
        from scipy.stats import spearmanr, kendalltau
        have_scipy = True
    except ImportError:
        have_scipy = False

    print(f"\n{'='*90}")
    print("  Rank correlation: lower metric ↔ higher BW")
    print(f"{'='*90}")

    for (plat, hw), grp in combined.groupby(["platform", "hw_type"]):
        print(f"\n  {plat} {hw.upper()} (n={len(grp)}):")

        metrics = ALL_METRICS if hw == "cpu" else [("mu", "μ")]

        bw = grp["bw_median"].values

        print(f"    {'Metric':<12} {'ρ':>8} {'R²':>8} {'τ':>8} "
              f"{'p(ρ)':>12} {'p(τ)':>12} {'exact':>8}")
        print(f"    {'-'*68}")

        for col, sym in metrics:
            if col not in grp.columns:
                print(f"    {sym:<12} {'N/A':>8}")
                continue
            vals = grp[col].values
            ok = np.isfinite(bw) & np.isfinite(vals)
            nok = ok.sum()
            if nok < 3:
                print(f"    {sym:<12} {'n<3':>8}")
                continue

            m_rank = pd.Series(vals[ok]).rank(ascending=True).values
            bw_r = pd.Series(bw[ok]).rank(ascending=False).values

            if have_scipy:
                rho, p_rho = spearmanr(m_rank, bw_r)
                tau, p_tau = kendalltau(m_rank, bw_r)
            else:
                d6 = 6 * np.sum((bw_r - m_rank)**2)
                rho = 1 - d6 / (nok * (nok**2 - 1))
                p_rho = tau = p_tau = float('nan')

            r2 = rho**2
            n_exact = int(np.sum(m_rank == bw_r))

            print(f"    {sym:<12} {rho:>+8.4f} {r2:>8.4f} {tau:>+8.4f} "
                  f"{p_rho:>12.2e} {p_tau:>12.2e} {n_exact:>4}/{nok}")


def plot_correlations(combined, output_stem="bw_vs_metrics"):
    """Scatter plots: BW vs each metric, one subplot per metric."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not available, skipping plots")
        return

    LAYOUT_COLOR = {"klon_first": "#e67e22", "klev_first": "#2980b9"}
    LAYOUT_MARKER = {"klon_first": "o", "klev_first": "s"}

    for (plat, hw), grp in combined.groupby(["platform", "hw_type"]):
        metrics = ALL_METRICS if hw == "cpu" else [("mu", "μ")]
        n_met = len(metrics)
        fig, axes = plt.subplots(1, n_met, figsize=(4.5 * n_met, 4), squeeze=False)
        fig.suptitle(f"{plat} — {hw.upper()}: BW vs Cost Metrics", fontsize=14)

        for mi, (col, sym) in enumerate(metrics):
            ax = axes[0, mi]
            if col not in grp.columns:
                ax.set_visible(False)
                continue

            for lo in ["klon_first", "klev_first"]:
                sub = grp[grp["loop_order"] == lo]
                if sub.empty:
                    continue
                ax.scatter(sub[col], sub["bw_median"],
                           c=LAYOUT_COLOR.get(lo, "gray"),
                           marker=LAYOUT_MARKER.get(lo, "o"),
                           s=50, alpha=0.8, edgecolors="black", linewidths=0.4,
                           label=lo.replace("_", "-"))

            # Compute Spearman for annotation
            ok = np.isfinite(grp[col]) & np.isfinite(grp["bw_median"])
            if ok.sum() >= 3:
                try:
                    from scipy.stats import spearmanr
                    rho, _ = spearmanr(grp.loc[ok, col], -grp.loc[ok, "bw_median"])
                    ax.text(0.05, 0.05, f"ρ={rho:+.3f}\nR²={rho**2:.3f}",
                            transform=ax.transAxes, fontsize=9,
                            verticalalignment="bottom",
                            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
                except ImportError:
                    pass

            ax.set_xlabel(sym, fontsize=12)
            ax.set_ylabel("BW [TB/s]" if mi == 0 else "", fontsize=12)
            ax.set_title(sym, fontsize=13)
            ax.grid(alpha=0.3)

            # Log scale if range > 10×
            vals = grp[col].dropna()
            if len(vals) > 0 and vals.max() / max(vals.min(), 1e-12) > 10:
                ax.set_xscale("log")

            if mi == 0:
                ax.legend(fontsize=9, loc="upper right")

        plat_short = plat.replace(" ", "_").replace("/", "_")
        fname = f"{output_stem}_{plat_short}_{hw}"
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        fig.savefig(f"{fname}.png", dpi=180, bbox_inches="tight")
        fig.savefig(f"{fname}.pdf", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"  Plot saved: {fname}.png/.pdf")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--platform", action="append", default=None,
                   choices=["beverin", "daint"])
    p.add_argument("-o", "--output", default="bw_vs_metrics.csv")
    p.add_argument("--no-plot", action="store_true", help="Skip plots")
    args = p.parse_args()

    if args.platform is None:
        args.platform = ["beverin"]

    all_joined = []
    for plat in args.platform:
        cfg = PLATFORM_CFG[plat]
        print(f"\nProcessing {cfg['label']} ({plat})...")

        cost_df = load_safe(cfg["cost_csv"])
        if cost_df is None:
            cost_df = load_safe("results_full.csv")
        if cost_df is None:
            print(f"  [ERROR] No cost CSV found for {plat}")
            continue

        cpu_df = load_safe(cfg["cpu_csv"])
        cpu_rt = extract_runtime_bw(cpu_df, "cpu")
        cpu_best = best_schedule_per_config(cpu_rt)
        cpu_joined = join_with_cost(cpu_best, cost_df, "cpu")

        gpu_df = load_safe(cfg["gpu_csv"])
        gpu_rt = extract_runtime_bw(gpu_df, "gpu",
                                     cfg["gpu_filter_col"], cfg["gpu_filter_val"])
        gpu_best = best_schedule_per_config(gpu_rt)
        gpu_joined = join_with_cost(gpu_best, cost_df, "gpu")

        joined = pd.concat([cpu_joined, gpu_joined], ignore_index=True)
        if not joined.empty:
            joined["platform"] = cfg["label"]

        print_table(joined, cfg["label"])
        all_joined.append(joined)

    combined = pd.concat(all_joined, ignore_index=True)
    if not combined.empty:
        combined.to_csv(args.output, index=False)
        print(f"\nCombined CSV written to: {args.output}  ({len(combined)} rows)")
        compute_correlations(combined)
        if not args.no_plot:
            plot_correlations(combined)
    else:
        print("\nNo data to write.")


if __name__ == "__main__":
    main()