#!/usr/bin/env python3
"""
join_bw_metrics.py

Joins cost-model metrics (mu, delta) with measured bandwidth from runtime CSVs.
Produces a combined CSV and prints a readable table.

Runtime CSVs expected columns:
  CPU: variant, cell_dist, parallelization, nlev, time_ms
  GPU: variant, cell_dist, config_label,    nlev, time_ms

Cost CSV expected columns:
  target, block_bytes, vec_width, variant, schedule, loop_order, dist,
  T, mu, delta, beta, alpha, gamma, P_NUMA

Usage:
    python join_bw_metrics.py                              # defaults
    python join_bw_metrics.py --platform beverin           # MI300A
    python join_bw_metrics.py --platform daint             # Grace
    python join_bw_metrics.py --platform beverin --platform daint  # both
"""
import argparse
import pandas as pd
import numpy as np
import sys

# ---- Transfer volume for bandwidth calculation ----
NPROMA = 81920
NLEV   = 96
BYTES  = (2*NPROMA*4 + 2*NPROMA*4 +
          NLEV*NPROMA*8 + NLEV*NPROMA*8 + NPROMA*8 +
          NLEV*NPROMA*8 + NLEV*NPROMA*8 +
          NPROMA*8 + NPROMA*8 + NLEV*NPROMA*8)

# ---- Naming maps: runtime CSV → cost CSV ----
DIST_TO_COST = {
    "uniform":     "uniform",
    "normal_var1": "normal1",
    "normal_var4": "normal4",
    "sequential":  "sequential",
}

SCHED_TO_COST = {
    "omp_for":           "omp_for",
    "omp_collapse2":     "omp_collapse2",
    "opt_for":           "omp_for",
    "opt_collapse_tile": "omp_collapse2",
}

# ---- Which cost-model target to match for each runtime type ----
# CPU runtime → CPU_scalar metric (B=64, W=1)
# GPU runtime → GPU_warp32 metric (B=128, W=32)
RUNTIME_TO_COST_TARGET = {
    "cpu": "CPU_AVX512",
    "gpu": "GPU_warp32",
}

# ---- Platform configs ----
PLATFORM_CFG = {
    "beverin": {
        "label":    "MI300A Zen 4",
        "cpu_csv":  "z_v_grad_w_cpu_beverin.csv",
        "gpu_csv":  "z_v_grad_w_gpu_beverin.csv",
        "cost_csv": "results_full_beverin.csv",
        "gpu_filter_col": "config_label",
        "gpu_filter_val": "1x1_32x16",
    },
    "daint": {
        "label":    "Grace / H200",
        "cpu_csv":  "z_v_grad_w_cpu_daint.csv",
        "gpu_csv":  "z_v_grad_w_gpu_daint.csv",
        "cost_csv": "results_full_daint.csv",
        "gpu_filter_col": "config_label",
        "gpu_filter_val": "1x1_32x16",
    },
}


def load_safe(path):
    try:
        df = pd.read_csv(path)
        if df.empty:
            raise ValueError
        return df
    except Exception:
        return None


def compute_bw(df):
    """Add bandwidth column in TB/s."""
    df = df.copy()
    df["bw_tbs"] = BYTES / (df["time_ms"] * 1e-3) / 1e12
    return df


def extract_runtime_bw(df, hw_type, gpu_filter_col=None, gpu_filter_val=None):
    """
    For each (variant, cell_dist, parallelization), compute median bandwidth.
    Returns DataFrame with columns: variant, dist_cost, schedule_cost, bw_median, bw_std, hw_type
    """
    if df is None:
        return pd.DataFrame()

    df = df.copy()

    # GPU: filter to specific config
    if hw_type == "gpu" and gpu_filter_col and gpu_filter_val:
        df = df[df[gpu_filter_col] == gpu_filter_val]

    if "time_ms" not in df.columns:
        print(f"  [WARN] no time_ms column in {hw_type} CSV", file=sys.stderr)
        return pd.DataFrame()

    df = compute_bw(df)

    # Group by (variant, cell_dist, parallelization)
    if hw_type == "cpu":
        group_cols = ["variant", "cell_dist", "parallelization"]
    else:
        group_cols = ["variant", "cell_dist"]
        if "parallelization" in df.columns:
            group_cols.append("parallelization")

    rows = []
    for keys, grp in df.groupby(group_cols):
        if hw_type == "cpu":
            var, dist, sched = keys
        else:
            if len(keys) == 3:
                var, dist, sched = keys
            else:
                var, dist = keys
                sched = "gpu_kernel"

        dist_cost = DIST_TO_COST.get(dist, dist)
        sched_cost = SCHED_TO_COST.get(sched, sched)

        rows.append({
            "variant_num": int(var),
            "variant":     f"V{int(var)}",
            "dist_cost":   dist_cost,
            "schedule_rt": sched,
            "schedule_cost": sched_cost,
            "bw_median":   grp["bw_tbs"].median(),
            "bw_std":      grp["bw_tbs"].std(),
            "n_samples":   len(grp),
            "hw_type":     hw_type,
        })

    return pd.DataFrame(rows)


def best_schedule_per_config(rt_df):
    """For each (variant, dist, hw_type), keep only the schedule with highest median BW."""
    if rt_df.empty:
        return rt_df
    idx = rt_df.groupby(["variant", "dist_cost", "hw_type"])["bw_median"].idxmax()
    return rt_df.loc[idx].reset_index(drop=True)


def join_with_cost(rt_best, cost_df, hw_type):
    """
    Join runtime BW with cost-model metrics.
    Matches on (variant, dist, schedule) using the appropriate cost target.
    """
    if rt_best.empty or cost_df is None or cost_df.empty:
        return pd.DataFrame()

    cost_target = RUNTIME_TO_COST_TARGET[hw_type]

    rows = []
    for _, r in rt_best[rt_best["hw_type"] == hw_type].iterrows():
        mask = (
            (cost_df["variant"] == r["variant"]) &
            (cost_df["dist"] == r["dist_cost"]) &
            (cost_df["schedule"] == r["schedule_cost"]) &
            (cost_df["target"] == cost_target)
        )
        crows = cost_df[mask]
        if crows.empty:
            # Try without schedule match (GPU kernels don't map 1:1)
            mask2 = (
                (cost_df["variant"] == r["variant"]) &
                (cost_df["dist"] == r["dist_cost"]) &
                (cost_df["target"] == cost_target)
            )
            crows = cost_df[mask2]

        if crows.empty:
            mu, delta = np.nan, np.nan
            loop_order = "?"
        else:
            c = crows.iloc[0]
            mu = c["mu"]
            delta = c["delta"]
            loop_order = c["loop_order"]

        rows.append({
            "variant":      r["variant"],
            "dist":         r["dist_cost"],
            "schedule":     r["schedule_rt"],
            "loop_order":   loop_order,
            "hw_type":      hw_type,
            "cost_target":  cost_target,
            "bw_median":    r["bw_median"],
            "bw_std":       r["bw_std"],
            "mu":           mu,
            "delta":        delta,
        })

    return pd.DataFrame(rows)


def print_table(joined, platform_label):
    """Pretty-print the joined table."""
    if joined.empty:
        return

    print(f"\n{'='*100}")
    print(f"  {platform_label}")
    print(f"{'='*100}")

    for hw in ["cpu", "gpu"]:
        sub = joined[joined["hw_type"] == hw].copy()
        if sub.empty:
            continue

        # CPU ranks by Δ (lower = better), GPU ranks by μ (lower = better)
        # Also compute μ·Δ as combined CPU metric
        if hw == "cpu":
            sub["mu_delta"] = sub["mu"] * sub["delta"]
            metric_col = "mu_delta"
            metric_sym = "μ·Δ"
        else:
            metric_col = "mu"
            metric_sym = "μ"

        sub["bw_rank"] = sub["bw_median"].rank(ascending=False).astype(int)
        sub["metric_rank"] = sub[metric_col].rank(ascending=True).astype(int)
        sub = sub.sort_values("bw_rank")

        hw_label = "CPU" if hw == "cpu" else "GPU"
        target = sub.iloc[0]["cost_target"]
        print(f"\n  --- {hw_label} (metric target: {target}, ranking by: {metric_sym}) ---")

        if hw == "cpu":
            print(f"  {'#bw':>3} {'#'+metric_sym:>5} {'Variant':<8} {'Distribution':<14} {'Schedule':<18} "
                  f"{'Loop Order':<14} {'BW [TB/s]':>10} {'±':>6} "
                  f"{'μ':>8} {'Δ':>10} {'μ·Δ':>12}")
            print(f"  {'-'*115}")
        else:
            print(f"  {'#bw':>3} {'#'+metric_sym:>3} {'Variant':<8} {'Distribution':<14} {'Schedule':<18} "
                  f"{'Loop Order':<14} {'BW [TB/s]':>10} {'±':>6} "
                  f"{'μ':>8} {'Δ':>10}")
            print(f"  {'-'*100}")

        for _, r in sub.iterrows():
            bw_s = f"{r['bw_median']:.4f}"
            std_s = f"{r['bw_std']:.4f}" if np.isfinite(r['bw_std']) else "---"
            mu_s = f"{r['mu']:.2f}" if np.isfinite(r['mu']) else "---"
            delta_s = f"{r['delta']:.1f}" if np.isfinite(r['delta']) else "---"
            match = "✓" if r["bw_rank"] == r["metric_rank"] else " "

            if hw == "cpu":
                mud = r["mu"] * r["delta"]
                mud_s = f"{mud:.0f}" if np.isfinite(mud) else "---"
                print(f"  {r['bw_rank']:>3} {r['metric_rank']:>5}{match} {r['variant']:<8} {r['dist']:<14} {r['schedule']:<18} "
                      f"{r['loop_order']:<14} {bw_s:>10} {std_s:>6} "
                      f"{mu_s:>8} {delta_s:>10} {mud_s:>12}")
            else:
                print(f"  {r['bw_rank']:>3} {r['metric_rank']:>3}{match} {r['variant']:<8} {r['dist']:<14} {r['schedule']:<18} "
                      f"{r['loop_order']:<14} {bw_s:>10} {std_s:>6} "
                      f"{mu_s:>8} {delta_s:>10}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--platform", action="append", default=None,
                   choices=["beverin", "daint"],
                   help="Platform(s) to process (can repeat)")
    p.add_argument("-o", "--output", default="bw_vs_metrics.csv",
                   help="Output CSV path")
    args = p.parse_args()

    if args.platform is None:
        args.platform = ["beverin"]

    all_joined = []

    for plat in args.platform:
        cfg = PLATFORM_CFG[plat]
        print(f"\nProcessing {cfg['label']} ({plat})...")

        cost_df = load_safe(cfg["cost_csv"])
        if cost_df is None:
            # Fallback to generic results_full.csv
            cost_df = load_safe("results_full.csv")
        if cost_df is None:
            print(f"  [ERROR] No cost CSV found for {plat}")
            continue

        # --- CPU ---
        cpu_df = load_safe(cfg["cpu_csv"])
        cpu_rt = extract_runtime_bw(cpu_df, "cpu")
        cpu_best = best_schedule_per_config(cpu_rt)
        cpu_joined = join_with_cost(cpu_best, cost_df, "cpu")

        # --- GPU ---
        gpu_df = load_safe(cfg["gpu_csv"])
        gpu_rt = extract_runtime_bw(gpu_df, "gpu",
                                     cfg["gpu_filter_col"],
                                     cfg["gpu_filter_val"])
        gpu_best = best_schedule_per_config(gpu_rt)
        gpu_joined = join_with_cost(gpu_best, cost_df, "gpu")

        joined = pd.concat([cpu_joined, gpu_joined], ignore_index=True)
        if not joined.empty:
            joined["platform"] = cfg["label"]

        print_table(joined, cfg["label"])
        all_joined.append(joined)

    # --- Output combined CSV ---
    combined = pd.concat(all_joined, ignore_index=True)
    if not combined.empty:
        # Add ranks per (platform, hw_type) group
        # CPU ranks by μ·Δ, GPU ranks by μ
        for (plat, hw), idx in combined.groupby(["platform", "hw_type"]).groups.items():
            if hw == "cpu":
                combined.loc[idx, "mu_delta"] = combined.loc[idx, "mu"] * combined.loc[idx, "delta"]
                metric_col = "mu_delta"
                metric_name = "mu_delta"
            else:
                metric_col = "mu"
                metric_name = "mu"
            combined.loc[idx, "bw_rank"] = combined.loc[idx, "bw_median"].rank(ascending=False).astype(int)
            combined.loc[idx, "metric_rank"] = combined.loc[idx, metric_col].rank(ascending=True).astype(int)
            combined.loc[idx, "ranking_metric"] = metric_name

        combined.to_csv(args.output, index=False)
        print(f"\nCombined CSV written to: {args.output}")
        print(f"  Rows: {len(combined)}")

        # --- Correlation summary ---
        print(f"\n{'='*60}")
        print("  Rank correlation: lower metric ↔ higher BW")
        print(f"  (CPU uses Δ, GPU uses μ)")
        print(f"{'='*60}")
        for (plat, hw), grp in combined.groupby(["platform", "hw_type"]):
            if hw == "cpu":
                metric_col = "delta"
                metric_sym = "μ·Δ"
            else:
                metric_col = "mu"
                metric_sym = "μ"
            mask = np.isfinite(grp["bw_median"]) & np.isfinite(grp[metric_col])
            n = mask.sum()
            if n < 3:
                print(f"  {plat} {hw}: too few points ({n})")
                continue

            bw_rank = grp.loc[mask, "bw_rank"].values
            metric_rank = grp.loc[mask, "metric_rank"].values
            n_exact = int(np.sum(bw_rank == metric_rank))

            try:
                from scipy.stats import spearmanr, kendalltau
                rho, p_rho = spearmanr(metric_rank, bw_rank)
                tau, p_tau = kendalltau(metric_rank, bw_rank)
            except ImportError:
                d = 6 * np.sum((bw_rank - metric_rank)**2)
                rho = 1 - d / (n * (n**2 - 1))
                p_rho = float('nan')
                tau = float('nan')
                p_tau = float('nan')

            r2 = rho ** 2

            print(f"  {plat} {hw.upper()} — ranking by {metric_sym} (n={n}):")
            print(f"    Spearman ρ  = {rho:+.4f}  (p = {p_rho:.2e})")
            print(f"    Rank R²     = {r2:.4f}   ({r2*100:.1f}% of BW rank variance explained by {metric_sym} rank)")
            try:
                print(f"    Kendall τ   = {tau:+.4f}  (p = {p_tau:.2e})")
            except:
                pass
            print(f"    Exact match = {n_exact}/{n}  ({100*n_exact/n:.0f}%)")
    else:
        print("\nNo data to write.")


if __name__ == "__main__":
    main()