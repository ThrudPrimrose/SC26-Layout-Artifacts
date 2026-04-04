#!/usr/bin/env python3
"""
plot_sweep.py — Analyse and plot GOTO + NUMA GEMM sweep results.

Usage:
    python3 plot_sweep.py results/sweep_*/

Reads all .csv files in the given directory and produces:
    1. thread_scaling.pdf       — GF/s vs threads for each kernel
    2. numa_grid.pdf            — GF/s for different PXxPY grids
    3. blocking_heatmap.pdf     — MC vs KC heatmap (best kernel)
    4. amx_vs_vec.pdf           — AMX vs Vector microkernel comparison
    5. numa_distribution.pdf    — Static vs SUMMA vs Cannon bar chart
    6. violin_all.pdf           — Violin plots of all kernels
    7. summary.csv              — Aggregated statistics table
"""

import sys
import os
import glob
import numpy as np
import pandas as pd

# ── Plotting setup ───────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[WARN] matplotlib not found — skipping plots, writing summary only")

try:
    import seaborn as sns
    HAS_SNS = True
except ImportError:
    HAS_SNS = False

COLORS = {
    "ref_omp":        "#999999",
    "mkl_sgemm_f32":  "#1f77b4",
    "mkl_f16f16f32":  "#2ca02c",
    "goto_amx":       "#d62728",
    "goto_vec":       "#ff7f0e",
}

def get_color(name):
    for k, v in COLORS.items():
        if k in name:
            return v
    return None


# ── Data loading ─────────────────────────────────────────────────────────

def load_all(directory):
    """Load and concatenate all CSV files in the directory."""
    csvs = sorted(glob.glob(os.path.join(directory, "*.csv")))
    if not csvs:
        print(f"No CSV files found in {directory}")
        sys.exit(1)

    frames = []
    for path in csvs:
        basename = os.path.splitext(os.path.basename(path))[0]
        if basename == "all_results":
            continue  # skip the merged master
        try:
            df = pd.read_csv(path)
            df["source"] = basename
            frames.append(df)
        except Exception as e:
            print(f"  [WARN] skipping {path}: {e}")

    if not frames:
        print("No valid CSV data found")
        sys.exit(1)

    data = pd.concat(frames, ignore_index=True)

    # Filter out failed runs (time_ns < 0)
    data = data[data["time_ns"] > 0].copy()

    # Parse sweep metadata from source filename
    data["sweep"] = data["source"].apply(classify_sweep)

    print(f"Loaded {len(data)} rows from {len(frames)} files")
    return data


def classify_sweep(source):
    """Classify a source filename into sweep category."""
    if source.startswith("thread_"):
        return "thread_scaling"
    elif source.startswith("numa_grid_"):
        return "numa_grid"
    elif source.startswith("block_"):
        return "blocking"
    elif source.startswith("amx_") or source.startswith("vec_"):
        return "amx_vs_vec"
    elif source.startswith("numa_dist_"):
        return "numa_dist"
    else:
        return "other"


# ── Statistics ───────────────────────────────────────────────────────────

def summarise(data):
    """Compute per-kernel statistics."""
    stats = (data.groupby(["source", "kernel", "M", "N", "K", "threads"])
             .agg(
                 median_ns=("time_ns", "median"),
                 mean_ns=("time_ns", "mean"),
                 min_ns=("time_ns", "min"),
                 max_ns=("time_ns", "max"),
                 std_ns=("time_ns", "std"),
                 median_gflops=("gflops", "median"),
                 mean_gflops=("gflops", "mean"),
                 max_gflops=("gflops", "max"),
                 count=("gflops", "count"),
             )
             .reset_index())
    return stats


# ── Plot 1: Thread scaling ──────────────────────────────────────────────

def plot_thread_scaling(data, outdir):
    """GF/s vs thread count for each kernel."""
    df = data[data["sweep"] == "thread_scaling"].copy()
    if df.empty:
        print("  [SKIP] thread_scaling: no data")
        return

    stats = (df.groupby(["kernel", "threads"])
             .agg(med=("gflops", "median"), lo=("gflops", lambda x: x.quantile(0.1)),
                  hi=("gflops", lambda x: x.quantile(0.9)))
             .reset_index())

    fig, ax = plt.subplots(figsize=(10, 6))
    for kern in stats["kernel"].unique():
        kd = stats[stats["kernel"] == kern].sort_values("threads")
        color = get_color(kern)
        ax.plot(kd["threads"], kd["med"], "o-", label=kern, color=color, markersize=5)
        ax.fill_between(kd["threads"], kd["lo"], kd["hi"], alpha=0.15, color=color)

    ax.set_xlabel("Threads")
    ax.set_ylabel("GF/s (median)")
    ax.set_title("Thread Scaling — 8192×8192×8192 GEMM")
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "thread_scaling.pdf"), dpi=150)
    plt.close(fig)
    print("  → thread_scaling.pdf")


# ── Plot 2: NUMA grid comparison ────────────────────────────────────────

def plot_numa_grid(data, outdir):
    """Bar chart: GF/s for different NUMA grid configs."""
    df = data[data["sweep"] == "numa_grid"].copy()
    if df.empty:
        print("  [SKIP] numa_grid: no data")
        return

    # Extract grid from source name
    df["grid"] = df["source"].str.extract(r"numa_grid_(\d+x\d+)")

    # Only keep interesting kernels
    interesting = ["goto_amx", "goto_vec"]
    numa_kernels = df["kernel"].unique()
    keep = [k for k in numa_kernels if any(x in k for x in
            ["goto_amx", "goto_vec", "numa_static", "numa_summa", "numa_cannon", "mkl"])]
    df = df[df["kernel"].isin(keep)]

    stats = (df.groupby(["kernel", "grid"])
             .agg(med=("gflops", "median"))
             .reset_index())

    if stats.empty:
        print("  [SKIP] numa_grid: no matching kernels")
        return

    grids = sorted(stats["grid"].dropna().unique(),
                   key=lambda g: tuple(int(x) for x in g.split("x")))
    kernels = sorted(stats["kernel"].unique())

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(grids))
    width = 0.8 / max(len(kernels), 1)

    for i, kern in enumerate(kernels):
        kd = stats[stats["kernel"] == kern]
        vals = [kd[kd["grid"] == g]["med"].values[0] if len(kd[kd["grid"] == g]) > 0 else 0
                for g in grids]
        color = get_color(kern)
        ax.bar(x + i * width - 0.4 + width/2, vals, width,
               label=kern, color=color, alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(grids, rotation=45, ha="right")
    ax.set_xlabel("NUMA Grid (PX×PY)")
    ax.set_ylabel("GF/s (median)")
    ax.set_title("NUMA Grid Configuration — 8192×8192×8192")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "numa_grid.pdf"), dpi=150)
    plt.close(fig)
    print("  → numa_grid.pdf")


# ── Plot 3: Blocking parameter heatmap ──────────────────────────────────

def plot_blocking(data, outdir):
    """Heatmap of MC vs KC for the GOTO kernel."""
    df = data[data["sweep"] == "blocking"].copy()
    if df.empty:
        print("  [SKIP] blocking: no data")
        return

    # Extract MC, KC from source name
    df["MC"] = df["source"].str.extract(r"mc(\d+)").astype(float)
    df["KC"] = df["source"].str.extract(r"kc(\d+)").astype(float)

    # Focus on the goto_amx kernel
    goto_df = df[df["kernel"] == "goto_amx"]
    if goto_df.empty:
        goto_df = df[df["kernel"].str.contains("goto")]
    if goto_df.empty:
        print("  [SKIP] blocking heatmap: no goto kernel data")
        return

    pivot = (goto_df.groupby(["MC", "KC"])
             .agg(med=("gflops", "median"))
             .reset_index()
             .pivot(index="MC", columns="KC", values="med"))

    fig, ax = plt.subplots(figsize=(8, 6))
    if HAS_SNS:
        sns.heatmap(pivot, annot=True, fmt=".0f", cmap="YlOrRd", ax=ax,
                    cbar_kws={"label": "GF/s"})
    else:
        im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", origin="lower")
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns.astype(int))
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index.astype(int))
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = pivot.values[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8)
        plt.colorbar(im, ax=ax, label="GF/s")

    ax.set_xlabel("KC (rank-update depth)")
    ax.set_ylabel("MC (IC block height)")
    ax.set_title("GOTO Blocking Parameters — goto_amx 8192³")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "blocking_heatmap.pdf"), dpi=150)
    plt.close(fig)
    print("  → blocking_heatmap.pdf")


# ── Plot 4: AMX vs Vector ───────────────────────────────────────────────

def plot_amx_vs_vec(data, outdir):
    """Side-by-side AMX vs Vector comparison."""
    df = data[data["sweep"] == "amx_vs_vec"].copy()
    if df.empty:
        print("  [SKIP] amx_vs_vec: no data")
        return

    # Tag AMX vs VEC from source
    df["backend"] = df["source"].apply(lambda s: "AMX" if s.startswith("amx") else "VEC")

    # Keep only the GOTO kernels
    goto_df = df[df["kernel"].str.contains("goto")]
    if goto_df.empty:
        print("  [SKIP] amx_vs_vec: no goto kernels")
        return

    stats = (goto_df.groupby(["backend", "threads"])
             .agg(med=("gflops", "median"))
             .reset_index())

    fig, ax = plt.subplots(figsize=(8, 5))
    for backend in ["AMX", "VEC"]:
        bd = stats[stats["backend"] == backend].sort_values("threads")
        color = "#d62728" if backend == "AMX" else "#ff7f0e"
        ax.plot(bd["threads"], bd["med"], "o-", label=f"GOTO ({backend})",
                color=color, markersize=7, linewidth=2)

    ax.set_xlabel("Threads")
    ax.set_ylabel("GF/s (median)")
    ax.set_title("AMX vs AVX-512 Microkernel — 8192³ GEMM")
    ax.set_xscale("log", base=2)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "amx_vs_vec.pdf"), dpi=150)
    plt.close(fig)
    print("  → amx_vs_vec.pdf")


# ── Plot 5: NUMA distribution comparison ────────────────────────────────

def plot_numa_dist(data, outdir):
    """Bar chart: Static vs SUMMA vs Cannon."""
    df = data[data["sweep"] == "numa_dist"].copy()
    if df.empty:
        print("  [SKIP] numa_dist: no data")
        return

    # Keep NUMA kernels + baselines
    keep = df["kernel"].unique()
    stats = (df.groupby(["kernel", "source"])
             .agg(med=("gflops", "median"))
             .reset_index())

    fig, ax = plt.subplots(figsize=(12, 6))
    sources = stats["source"].unique()

    for src_idx, src in enumerate(sources):
        sd = stats[stats["source"] == src].sort_values("med", ascending=True)
        y = np.arange(len(sd))
        grid_label = src.replace("numa_dist_", "").replace("_", " ")
        ax.barh(y + src_idx * 0.3, sd["med"], height=0.25,
                label=grid_label, alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(sd["kernel"], fontsize=7)

    ax.set_xlabel("GF/s (median)")
    ax.set_title("NUMA Distribution — Static vs SUMMA vs Cannon")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "numa_distribution.pdf"), dpi=150)
    plt.close(fig)
    print("  → numa_distribution.pdf")


# ── Plot 6: Violin plots ────────────────────────────────────────────────

def plot_violins(data, outdir):
    """Violin plots showing distribution of run times for all kernels."""
    # Use the largest thread count sweep
    df = data.copy()
    max_t = df["threads"].max()
    df = df[df["threads"] == max_t]

    if df.empty or not HAS_SNS:
        print("  [SKIP] violins: no data or seaborn missing")
        return

    # Top 10 kernels by median GF/s
    top = (df.groupby("kernel")["gflops"].median()
           .nlargest(10).index.tolist())
    df = df[df["kernel"].isin(top)]

    fig, ax = plt.subplots(figsize=(14, 6))
    order = (df.groupby("kernel")["gflops"].median()
             .sort_values(ascending=False).index.tolist())

    sns.violinplot(data=df, x="kernel", y="gflops", order=order,
                   inner="quartile", cut=0, ax=ax, scale="width")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("GF/s")
    ax.set_title(f"Runtime Distribution — Top 10 Kernels (t={max_t})")
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "violin_all.pdf"), dpi=150)
    plt.close(fig)
    print("  → violin_all.pdf")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <results_directory>")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Not a directory: {directory}")
        sys.exit(1)

    print(f"Loading data from {directory}")
    data = load_all(directory)

    # Summary CSV
    stats = summarise(data)
    summary_path = os.path.join(directory, "summary.csv")
    stats.to_csv(summary_path, index=False, float_format="%.2f")
    print(f"\n  → summary.csv ({len(stats)} rows)")

    # Print top kernels
    print("\n  Top kernels by median GF/s:")
    top = (stats.sort_values("median_gflops", ascending=False)
           .head(20)[["source", "kernel", "threads", "median_gflops", "max_gflops"]])
    print(top.to_string(index=False))
    print()

    if not HAS_MPL:
        print("Install matplotlib + seaborn for plots:")
        print("  pip install matplotlib seaborn")
        return

    print("Generating plots...")
    plot_thread_scaling(data, directory)
    plot_numa_grid(data, directory)
    plot_blocking(data, directory)
    plot_amx_vs_vec(data, directory)
    plot_numa_dist(data, directory)
    plot_violins(data, directory)
    print("\nDone.")


if __name__ == "__main__":
    main()