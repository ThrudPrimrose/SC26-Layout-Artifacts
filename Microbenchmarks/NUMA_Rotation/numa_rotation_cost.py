#!/usr/bin/env python3
"""
numa_rotation_cost.py

Computes the three delta metrics for the NUMA rotation benchmark pattern:
  Δ_raw  : unweighted mean block distance
  Δ_uma  : with α discount for nearby blocks (prefetcher locality)
  Δ_numa : with α discount + γ penalty for cross-NUMA blocks (σ in the paper)

Also reads the runtime CSV (from numa_rotation_bench) and computes
correlations between each metric variant and measured bandwidth.

Outputs:
  - Console table: metrics per phase
  - LaTeX table: metrics + runtime + correlation
  - CSV: joined metrics + runtime

Usage:
    python numa_rotation_cost.py --runtime numa_rotation.csv
    python numa_rotation_cost.py --N 16777216 --P 4 --beta 64 --alpha 0.4 --gamma 2.0
"""
import numpy as np
import pandas as pd
from scipy import stats
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--N', type=int, default=1<<24, help='Elements per array per domain')
parser.add_argument('--P', type=int, default=4, help='Number of NUMA domains')
parser.add_argument('--B', type=int, default=8, help='Elements per cache block (64B / 8B)')
parser.add_argument('--beta', type=int, default=4, help='Page span in blocks (page_size / block_size)')
parser.add_argument('--alpha', type=float, default=0.3, help='Intra-page locality discount')
parser.add_argument('--gamma', type=float, default=2.0, help='Cross-NUMA penalty')
parser.add_argument('--tpd', type=int, default=0, help='Threads per domain (0=auto from P)')
parser.add_argument('--runtime', default=None, help='CSV from numa_rotation_bench')
parser.add_argument('--out-csv', default='numa_rotation_metrics.csv')
parser.add_argument('--out-tex', default='numa_rotation_table.tex')
args = parser.parse_args()

N       = args.N
P_N     = args.P
B       = args.B
BETA    = args.beta
ALPHA   = args.alpha
GAMMA   = args.gamma
TPD     = args.tpd if args.tpd > 0 else P_N  # threads per domain

# Total elements across all domains
N_total = N * P_N
N_blocks_per_domain = N // B

# ═══════════════════════════════════════════════════════════════
#  NUMA domain mapping
# ═══════════════════════════════════════════════════════════════

def block_to_domain(block_addr):
    """
    Under first-touch with schedule(static), domain d owns blocks
    [d * N_blocks_per_domain, (d+1) * N_blocks_per_domain).
    """
    return int(block_addr // N_blocks_per_domain) % P_N


# ═══════════════════════════════════════════════════════════════
#  Cost model: simulate the STREAM triad access pattern
#  for a single thread processing a contiguous chunk
# ═══════════════════════════════════════════════════════════════

def compute_metrics_for_phase(phase, thread_domain=0):
    """
    Simulate one thread from `thread_domain` accessing arrays
    owned by domain (thread_domain + phase) % P_N.

    STREAM triad: A[i] = B[i] + scalar * C[i]
    3 arrays, streaming access pattern.

    Returns dict with mu, delta_raw, delta_uma, delta_numa (sigma).
    """
    target_domain = (thread_domain + phase) % P_N

    # The thread processes a chunk of size N / TPD elements
    chunk_size = N // TPD
    # Starting element index within the target domain's array
    local_tid = 0  # first thread in the domain
    start_elem = local_tid * chunk_size

    # Block addresses: target domain's arrays start at target_domain * N
    # 3 arrays: offsets 0, N_total, 2*N_total in the address space
    arr_base = [target_domain * N, N_total + target_domain * N, 2 * N_total + target_domain * N]

    # For STREAM triad, at each step we access 3 consecutive elements
    # (one per array), all at the same index i.
    # The "step" in our model = one iteration of the inner loop.
    # With vector width W, each step processes W elements.
    W = 1  # scalar for now

    # Track block sets
    prev_blocks = set()
    sum_mu = 0.0
    sum_raw_dist = 0.0
    sum_uma_dist = 0.0
    sum_numa_dist = 0.0
    T = 0

    for i in range(start_elem, start_elem + chunk_size, W):
        curr_blocks = set()
        for a_off in arr_base:
            elem_addr = a_off + i
            blk = elem_addr // B
            curr_blocks.add(blk)

        # New blocks
        new_blocks = curr_blocks - prev_blocks
        n_new = len(new_blocks)
        sum_mu += n_new

        if n_new > 0 and len(prev_blocks) > 0:
            raw_total = 0.0
            uma_total = 0.0
            numa_total = 0.0

            prev_list = sorted(prev_blocks)
            for b in new_blocks:
                # Nearest block in prev set
                raw_dist = min(abs(b - bp) for bp in prev_list)

                # Raw
                raw_total += raw_dist

                # UMA: α discount if within page
                w_uma = ALPHA if raw_dist < BETA else 1.0
                uma_total += w_uma * raw_dist

                # NUMA: γ if cross-domain, else α/1.0
                # Domain of the new block
                # Map block back to domain: which domain's array does this block belong to?
                # Block addresses: domain d's arrays span
                #   [d*N/B, (d+1)*N/B) for array 0
                #   [N_total/B + d*N/B, ...] for array 1, etc.
                global_elem = b * B
                arr_idx = global_elem // N_total  # which of the 3 arrays
                within_arr = global_elem % N_total
                blk_domain = within_arr // N  # which NUMA domain owns this

                is_remote = (blk_domain != thread_domain)

                if is_remote:
                    w_numa = GAMMA
                else:
                    w_numa = ALPHA if raw_dist < BETA else 1.0
                numa_total += w_numa * raw_dist

            sum_raw_dist += raw_total / n_new
            sum_uma_dist += uma_total / n_new
            sum_numa_dist += numa_total / n_new
        elif n_new > 0:
            # First step: distance = 1 by convention
            sum_raw_dist += 1.0
            sum_uma_dist += ALPHA  # within page
            if phase == 0:
                sum_numa_dist += ALPHA  # local
            else:
                sum_numa_dist += GAMMA  # remote

        T += 1
        prev_blocks = curr_blocks

        # Early termination for large N (sample first 10K steps, pattern is periodic)
        if T >= 10000 and T >= chunk_size // (B * 10):
            break

    return {
        "phase": phase,
        "target_domain": target_domain,
        "T": T,
        "mu": sum_mu / T if T > 0 else 0,
        "delta_raw": sum_raw_dist / T if T > 0 else 0,
        "delta_uma": sum_uma_dist / T if T > 0 else 0,
        "delta_numa": sum_numa_dist / T if T > 0 else 0,
        "is_local": phase == 0,
    }


# ═══════════════════════════════════════════════════════════════
#  Compute metrics for all phases
# ═══════════════════════════════════════════════════════════════

print(f"=== NUMA Rotation Cost Model ===")
print(f"  N={N} elems/domain, P_N={P_N} domains, B={B} elems/block")
print(f"  β={BETA} blocks/page, α={ALPHA}, γ={GAMMA}")
print(f"  N_blocks/domain={N_blocks_per_domain}\n")

rows = []
for phase in range(P_N):
    m = compute_metrics_for_phase(phase, thread_domain=0)
    rows.append(m)

df_metrics = pd.DataFrame(rows)

# ── Console table ──
print(f"{'Phase':>6s} {'Target':>7s} {'Local':>6s} {'μ':>8s} "
      f"{'Δ_raw':>10s} {'Δ_uma':>10s} {'Δ_numa':>10s}")
print("-" * 65)
for _, r in df_metrics.iterrows():
    print(f"{int(r['phase']):>6d} {int(r['target_domain']):>7d} "
          f"{'yes' if r['is_local'] else 'no':>6s} "
          f"{r['mu']:>8.3f} {r['delta_raw']:>10.3f} "
          f"{r['delta_uma']:>10.3f} {r['delta_numa']:>10.3f}")

# ═══════════════════════════════════════════════════════════════
#  Join with runtime if available
# ═══════════════════════════════════════════════════════════════

if args.runtime:
    print(f"\n=== Loading runtime from {args.runtime} ===")
    rt = pd.read_csv(args.runtime)

    # Median BW per phase
    rt_med = rt.groupby("phase").agg(
        bw_median=("bw_gbs", "median"),
        bw_std=("bw_gbs", "std"),
        bw_min=("bw_gbs", "min"),
        bw_max=("bw_gbs", "max"),
    ).reset_index()

    df_joined = df_metrics.merge(rt_med, on="phase")
    df_joined["gamma_measured"] = df_joined["bw_median"].iloc[0] / df_joined["bw_median"]

    print(f"\n{'Phase':>6s} {'BW med':>10s} {'BW std':>8s} "
          f"{'γ_meas':>8s} {'Δ_raw':>10s} {'Δ_uma':>10s} {'Δ_numa':>10s}")
    print("-" * 70)
    for _, r in df_joined.iterrows():
        print(f"{int(r['phase']):>6d} {r['bw_median']:>10.2f} {r['bw_std']:>8.2f} "
              f"{r['gamma_measured']:>8.3f} {r['delta_raw']:>10.3f} "
              f"{r['delta_uma']:>10.3f} {r['delta_numa']:>10.3f}")

    # ── Correlations: each delta variant vs BW ──
    print("\n=== Correlation: Metric vs Measured BW ===")
    print("  (Higher metric → expect lower BW → negative correlation)\n")

    delta_cols = ["delta_raw", "delta_uma", "delta_numa", "mu"]
    delta_labels = ["Δ_raw", "Δ_UMA", "Δ_NUMA (σ)", "μ"]

    def safe_corr(func, x, y):
        if len(x) < 3: return np.nan, np.nan
        try: return func(x, y)
        except: return np.nan, np.nan

    def pearson_r2(x, y):
        r, p = stats.pearsonr(x, y)
        return r**2, p

    corr_funcs = [
        ("Spearman ρ",  stats.spearmanr),
        ("Kendall τ",   stats.kendalltau),
        ("Pearson R²",  pearson_r2),
    ]

    bw = df_joined["bw_median"].values
    print(f"  {'Metric':<15s}", end="")
    for cname, _ in corr_funcs:
        print(f"  {cname:>12s}", end="")
    print()
    print("  " + "-" * 55)

    corr_results = []
    for dc, dl in zip(delta_cols, delta_labels):
        mv = df_joined[dc].values
        print(f"  {dl:<15s}", end="")
        row_corr = {"metric": dl}
        for cname, cfunc in corr_funcs:
            r, p = safe_corr(cfunc, mv, bw)
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"  {r:>+10.4f}{sig:2s}", end="")
            row_corr[cname] = r
            row_corr[f"{cname}_p"] = p
        print()
        corr_results.append(row_corr)

    # ── Save joined CSV ──
    df_joined.to_csv(args.out_csv, index=False)
    print(f"\nWrote: {args.out_csv}")

    # ═══════════════════════════════════════════════════════════
    #  LaTeX table
    # ═══════════════════════════════════════════════════════════

    lines = []
    lines.append(r"\begin{table}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{NUMA rotation benchmark: cost metrics vs measured bandwidth."
                 f" $N={N}$, $P_{{\\mathcal{{N}}}}={P_N}$,"
                 f" $\\beta={BETA}$, $\\alpha={ALPHA}$, $\\gamma={GAMMA}$.}}")
    lines.append(r"\label{tab:numa-rotation}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{@{}crrrrrr@{}}")
    lines.append(r"\toprule")
    lines.append(r"Phase & $\bar\mu$ & $\Delta_{\text{raw}}$ & "
                 r"$\Delta_{\text{UMA}}$ & $\Delta_{\text{NUMA}}$ & "
                 r"BW [GB/s] & $\hat\gamma$ \\")
    lines.append(r"\midrule")

    for _, r in df_joined.iterrows():
        phase = int(r["phase"])
        label = "local" if phase == 0 else f"+{phase}"
        lines.append(
            f"  {label} & {r['mu']:.3f} & {r['delta_raw']:.2f} & "
            f"{r['delta_uma']:.2f} & {r['delta_numa']:.2f} & "
            f"{r['bw_median']:.1f} & {r['gamma_measured']:.3f} \\\\"
        )

    lines.append(r"\midrule")
    lines.append(r"\multicolumn{7}{@{}l}{\small Correlation with BW (Spearman $\rho$):} \\")
    lines.append(r"\addlinespace")

    # Correlation row
    corr_vals = []
    for dc, dl in zip(delta_cols[:4], delta_labels[:4]):
        mv = df_joined[dc].values
        r, _ = safe_corr(stats.spearmanr, mv, bw)
        corr_vals.append(f"{r:+.3f}" if not np.isnan(r) else "--")

    lines.append(f"  $\\rho$ & {corr_vals[3]} & {corr_vals[0]} & "
                 f"{corr_vals[1]} & {corr_vals[2]} & -- & -- \\\\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    tex = "\n".join(lines)
    with open(args.out_tex, 'w') as f:
        f.write(tex)
    print(f"Wrote: {args.out_tex}")
    print("\n--- LaTeX Preview ---")
    print(tex)

else:
    # No runtime: just save metrics
    df_metrics.to_csv(args.out_csv, index=False)
    print(f"\nWrote: {args.out_csv}")

    # Minimal LaTeX table (metrics only)
    lines = []
    lines.append(r"\begin{table}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{Predicted cost metrics for NUMA rotation pattern."
                 f" $\\beta={BETA}$, $\\alpha={ALPHA}$, $\\gamma={GAMMA}$.}}")
    lines.append(r"\label{tab:numa-rotation-predicted}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{@{}crrrr@{}}")
    lines.append(r"\toprule")
    lines.append(r"Phase & $\bar\mu$ & $\Delta_{\text{raw}}$ & "
                 r"$\Delta_{\text{UMA}}$ & $\Delta_{\text{NUMA}}$ \\")
    lines.append(r"\midrule")

    for _, r in df_metrics.iterrows():
        phase = int(r["phase"])
        label = "local" if phase == 0 else f"+{phase}"
        lines.append(
            f"  {label} & {r['mu']:.3f} & {r['delta_raw']:.2f} & "
            f"{r['delta_uma']:.2f} & {r['delta_numa']:.2f} \\\\"
        )

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    tex = "\n".join(lines)
    with open(args.out_tex, 'w') as f:
        f.write(tex)
    print(f"Wrote: {args.out_tex}")
    print("\n--- LaTeX Preview ---")
    print(tex)