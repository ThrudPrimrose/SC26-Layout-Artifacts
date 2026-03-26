#!/usr/bin/env python3
"""
2x2 violin grid (runtime) + 1x3 scatter (metric vs runtime).
Top:    rows=dist, cols=variant. Each cell: 2 violins (AMD, NV).
Bottom: scatter per metric. x=metric, y=median runtime.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

# ---- CONFIG ----
GPU_AMD_CSV = "z_v_grad_w_gpu_beverin.csv"
GPU_NV_CSV  = "z_v_grad_w_gpu_daint.csv"
METRICS_CSV = "metrics.csv"

NLEV = 96

VARIANTS = [1, 4]
DISTS    = ["uniform", "normal_var4"]

AMD_CONFIG = "1x1_128x1"
NV_CONFIG  = "1x1_32x8"

# variant -> (loop_order, metric_target)
METRIC_AMD = {1: ("klon_first", "gpu_wave64"), 4: ("klev_first", "gpu_wave64")}
METRIC_NV  = {1: ("klon_first", "gpu_warp32"), 4: ("klev_first", "gpu_warp32")}

OUT_STEM = "combined"
# ---- END CONFIG ----

VLAB = {1: "V1:none", 2: "V2:idx", 3: "V3:compute", 4: "V4:both"}
VCOL = {1: "#4c72b0", 2: "#dd8452", 3: "#55a868", 4: "#c44e52"}
PLATFORMS = [
    ("AMD", GPU_AMD_CSV, AMD_CONFIG, METRIC_AMD, "o"),
    ("NV",  GPU_NV_CSV,  NV_CONFIG,  METRIC_NV,  "s"),
]
PCOL = {"AMD": "#e74c3c", "NV": "#3498db"}


def main():
    # load runtime
    rt = {}
    for name, csv, cfg, _, _ in PLATFORMS:
        df = pd.read_csv(csv)
        rt[name] = df[(df["config_label"] == cfg) & (df["nlev"] == NLEV)]

    met = pd.read_csv(METRICS_CSV)
    met = met[met["nlev"] == NLEV]

    # figure with subfigures
    fig = plt.figure(figsize=(9, 11))
    sf_top, sf_bot = fig.subfigures(2, 1, height_ratios=[1.1, 1])

    sf_top.suptitle(
        f"GPU Runtime  (nlev={NLEV})  |  AMD:{AMD_CONFIG}  NV:{NV_CONFIG}",
        fontsize=11, fontweight="bold")

    # ---- top 2x2 violins ----
    axes_v = sf_top.subplots(len(DISTS), len(VARIANTS))

    for ri, dist in enumerate(DISTS):
        for ci, V in enumerate(VARIANTS):
            ax = axes_v[ri][ci]
            data, positions, colors = [], [], []

            for pi, (pname, _, _, _, _) in enumerate(PLATFORMS):
                vals = rt[pname][(rt[pname]["variant"] == V) &
                                 (rt[pname]["cell_dist"] == dist)]["time_ms"].values
                if len(vals) == 0:
                    continue
                data.append(vals)
                positions.append(pi)
                colors.append(PCOL[pname])

            if data:
                parts = ax.violinplot(data, positions=positions,
                                      showmeans=True, showmedians=True,
                                      showextrema=False, widths=0.6)
                for i, body in enumerate(parts["bodies"]):
                    body.set_facecolor(colors[i])
                    body.set_edgecolor("black")
                    body.set_alpha(0.7)
                parts["cmeans"].set_color("black")
                parts["cmedians"].set_color("white")
                parts["cmedians"].set_linewidth(1.5)

            pnames = [p[0] for p in PLATFORMS]
            ax.set_xticks(range(len(pnames)))
            ax.set_xticklabels(pnames, fontsize=9)
            ax.set_ylabel("time [ms]", fontsize=8)
            ax.set_title(f"{VLAB[V]}  |  {dist}", fontsize=9, fontweight="bold")
            ax.grid(axis="y", alpha=0.3)

    # ---- bottom 1x3 scatter ----
    sf_bot.suptitle("Cost Model Metrics vs Median Runtime", fontsize=11,
                    fontweight="bold")
    axes_s = sf_bot.subplots(1, 3)

    # collect points
    points = []
    for pname, _, cfg, mmap, marker in PLATFORMS:
        for V in VARIANTS:
            lo, tgt = mmap[V]
            for dist in DISTS:
                vals = rt[pname][(rt[pname]["variant"] == V) &
                                 (rt[pname]["cell_dist"] == dist)]["time_ms"]
                if len(vals) == 0:
                    continue
                mrow = met[(met["variant"] == V) &
                           (met["cell_dist"] == dist) &
                           (met["loop_order"] == lo) &
                           (met["target"] == tgt)]
                if len(mrow) == 0:
                    continue
                points.append(dict(
                    plat=pname, V=V, dist=dist, marker=marker,
                    med=vals.median(),
                    mu=mrow.iloc[0]["mu"],
                    delta=mrow.iloc[0]["delta"],
                    sigma=mrow.iloc[0]["sigma"]))

    mkeys = [("mu", r"$\mu$ (new blocks)"),
             ("delta", r"$\Delta$ (block distance)"),
             ("sigma", r"$\sigma$ (stride/16)")]

    dist_fill = {"uniform": True, "normal_var4": False}

    for mi, (mk, mlab) in enumerate(mkeys):
        ax = axes_s[mi]
        for p in points:
            fc = VCOL[p["V"]] if dist_fill[p["dist"]] else "none"
            ec = VCOL[p["V"]]
            ax.scatter(p[mk], p["med"], marker=p["marker"], s=90,
                       facecolors=fc, edgecolors=ec, linewidths=1.5, zorder=5)
        ax.set_xlabel(mlab, fontsize=9)
        if mi == 0:
            ax.set_ylabel("median time [ms]", fontsize=9)
        ax.grid(alpha=0.3)

    # legend
    handles = []
    for V in VARIANTS:
        handles.append(Line2D([0], [0], marker="o", color="w", markersize=8,
                              markerfacecolor=VCOL[V], markeredgecolor=VCOL[V],
                              label=VLAB[V]))
    for pname, _, _, _, mk in PLATFORMS:
        handles.append(Line2D([0], [0], marker=mk, color="w", markersize=8,
                              markerfacecolor=PCOL[pname],
                              markeredgecolor="black", label=pname))
    handles.append(Line2D([0], [0], marker="o", color="w", markersize=8,
                          markerfacecolor="gray", markeredgecolor="black",
                          label="uniform"))
    handles.append(Line2D([0], [0], marker="o", color="w", markersize=8,
                          markerfacecolor="none", markeredgecolor="black",
                          label="normal_var4"))
    sf_bot.legend(handles=handles, loc="lower center", ncol=6, fontsize=7,
                  bbox_to_anchor=(0.5, -0.02))

    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT_STEM}_nlev{NLEV}.{ext}", dpi=180)
    print(f"Saved {OUT_STEM}_nlev{NLEV}.{{pdf,png}}")


if __name__ == "__main__":
    main()