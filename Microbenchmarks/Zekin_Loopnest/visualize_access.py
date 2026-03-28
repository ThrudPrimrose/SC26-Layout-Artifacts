#!/usr/bin/env python3
"""
visualize_access.py — Visualize gather access patterns for ker_A vs ker_B
on a small 2D grid (single block).

Shows which elements of kin[] each warp touches, colored by warp ID,
illustrating why layout A gets L2 reuse and layout B gets L1 coalescing.

Usage:
    python visualize_access.py
    python visualize_access.py --np 16 --nl 8 --neighbors local
    python visualize_access.py --np 16 --nl 8 --neighbors random --seed 42 -o access.pdf
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import cm

def make_neighbors(NP, mode, seed=42):
    """Generate neighbor index arrays: 3 neighbors per cell."""
    rng = np.random.default_rng(seed)
    ix = np.zeros((3, NP), dtype=int)
    if mode == "local":
        ix[0] = np.arange(NP)                          # self
        ix[1] = np.clip(np.arange(NP) - 1, 0, NP - 1)  # left
        ix[2] = np.clip(np.arange(NP) + 1, 0, NP - 1)  # right
    elif mode == "random":
        for n in range(3):
            ix[n] = rng.integers(0, NP, size=NP)
    elif mode == "clustered":
        # neighbors within ±4 of self
        for n in range(3):
            offsets = rng.integers(-4, 5, size=NP)
            ix[n] = np.clip(np.arange(NP) + offsets, 0, NP - 1)
    return ix


def addr_layout_A(jc, jk, NP, NL):
    """Memory address for kin[jc, jk] in layout A: [nproma, nlev]"""
    return jc + jk * NP


def addr_layout_B(jc, jk, NP, NL):
    """Memory address for kin[jc, jk] in layout B: [nlev, nproma]"""
    return jk + jc * NL


def get_warp_accesses(NP, NL, ix, layout, warp_size=32):
    """
    Simulate kernel launch and return per-warp access info.

    Layout A: threads map (threadIdx.x → jc), grid.y → jk
      warp covers 32 consecutive jc at one jk
    Layout B: threads map (threadIdx.x → jk), grid.y → jc
      warp covers 32 consecutive jk at one jc

    Returns list of dicts: {warp_id, threads: [(jc, jk)], accesses: [(jc_neighbor, jk, addr)]}
    """
    warps = []
    warp_id = 0

    if layout == "A":
        # ker_A: blockDim.x covers jc, grid.y covers jk
        for jk in range(NL):
            for jc_start in range(0, NP, warp_size):
                threads = []
                accesses = []
                for t in range(min(warp_size, NP - jc_start)):
                    jc = jc_start + t
                    threads.append((jc, jk))
                    for n in range(3):
                        neighbor_jc = ix[n, jc]
                        addr = addr_layout_A(neighbor_jc, jk, NP, NL)
                        accesses.append((neighbor_jc, jk, addr))
                warps.append({"warp_id": warp_id, "threads": threads, "accesses": accesses})
                warp_id += 1
    else:  # layout B
        # ker_B: blockDim.x covers jk, grid.y covers jc
        for jc in range(NP):
            for jk_start in range(0, NL, warp_size):
                threads = []
                accesses = []
                for t in range(min(warp_size, NL - jk_start)):
                    jk = jk_start + t
                    threads.append((jc, jk))
                    for n in range(3):
                        neighbor_jc = ix[n, jc]
                        addr = addr_layout_B(neighbor_jc, jk, NP, NL)
                        accesses.append((neighbor_jc, jk, addr))
                warps.append({"warp_id": warp_id, "threads": threads, "accesses": accesses})
                warp_id += 1

    return warps


def visualize(NP, NL, ix, layout, warps, ax, title, max_warps_to_show=6,
              cacheline_bytes=128, elem_bytes=4):
    """
    Draw a heatmap of the kin[] array (NP × NL logical grid).
    Color each cell by which warp accesses it.
    Overlay cache line boundaries.
    """
    elems_per_cl = cacheline_bytes // elem_bytes  # 32 floats per 128B cacheline

    # Build access map: for each (jc, jk) → set of warp_ids that read it
    access_map = {}  # (jc, jk) → set of warp_ids
    for w in warps[:max_warps_to_show]:
        wid = w["warp_id"]
        for (njc, njk, addr) in w["accesses"]:
            key = (njc, njk)
            if key not in access_map:
                access_map[key] = set()
            access_map[key].add(wid)

    # Create color grid: -1 = not accessed, warp_id for single, -2 for multi
    grid = np.full((NL, NP), -1, dtype=int)
    multi_grid = np.zeros((NL, NP), dtype=bool)

    for (jc, jk), wids in access_map.items():
        if 0 <= jc < NP and 0 <= jk < NL:
            if len(wids) == 1:
                grid[jk, jc] = list(wids)[0] % max_warps_to_show
            else:
                grid[jk, jc] = -2
                multi_grid[jk, jc] = True

    # Color map: distinct colors for warps, gray for not accessed, white/hatched for shared
    warp_colors = plt.cm.Set2(np.linspace(0, 1, max_warps_to_show))
    colors = [(0.92, 0.92, 0.92, 1.0)]  # -1: not accessed (light gray)
    colors.append((1.0, 0.4, 0.4, 1.0))  # -2: multi-warp (red)
    for i in range(max_warps_to_show):
        colors.append(tuple(warp_colors[i]))
    cmap = ListedColormap(colors)
    # Map grid values: -1→0, -2→1, 0→2, 1→3, ...
    display = np.where(grid == -1, 0, np.where(grid == -2, 1, grid + 2))
    bounds = np.arange(-0.5, max_warps_to_show + 2.5, 1)
    norm = BoundaryNorm(bounds, cmap.N)

    ax.imshow(display, cmap=cmap, norm=norm, origin="upper", aspect="auto",
              extent=[-0.5, NP - 0.5, NL - 0.5, -0.5])

    # Draw cache line boundaries
    if layout == "A":
        # Layout A: address = jc + jk*NP, so cache lines are horizontal bands in jc
        for jk in range(NL):
            for cl_start in range(0, NP, elems_per_cl):
                cl_end = min(cl_start + elems_per_cl, NP)
                rect = plt.Rectangle((cl_start - 0.5, jk - 0.5), cl_end - cl_start, 1,
                                     linewidth=0.8, edgecolor="black", facecolor="none",
                                     linestyle="-")
                ax.add_patch(rect)
    else:
        # Layout B: address = jk + jc*NL, so cache lines are vertical bands in jk
        for jc in range(NP):
            for cl_start in range(0, NL, elems_per_cl):
                cl_end = min(cl_start + elems_per_cl, NL)
                rect = plt.Rectangle((jc - 0.5, cl_start - 0.5), 1, cl_end - cl_start,
                                     linewidth=0.8, edgecolor="black", facecolor="none",
                                     linestyle="-")
                ax.add_patch(rect)

    # Mark the thread positions of the shown warps
    for w in warps[:max_warps_to_show]:
        wid = w["warp_id"] % max_warps_to_show
        for (jc, jk) in w["threads"]:
            ax.plot(jc, jk, marker="x", color="black", markersize=4, markeredgewidth=0.8)

    ax.set_xlabel("jc (nproma)")
    ax.set_ylabel("jk (nlev)")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlim(-0.5, NP - 0.5)
    ax.set_ylim(NL - 0.5, -0.5)

    # Tick marks
    ax.set_xticks(range(0, NP, max(1, NP // 8)))
    ax.set_yticks(range(0, NL, max(1, NL // 8)))

    return warp_colors


def count_cachelines(warps, NP, NL, layout, cacheline_bytes=128, elem_bytes=4):
    """Count unique cache lines touched by each warp and total."""
    elems_per_cl = cacheline_bytes // elem_bytes
    per_warp = []
    all_cls = set()
    for w in warps:
        cls = set()
        for (njc, njk, addr) in w["accesses"]:
            cls.add(addr // elems_per_cl)
        per_warp.append(len(cls))
        all_cls.update(cls)
    return per_warp, len(all_cls)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--np", type=int, default=16, help="nproma (columns)")
    parser.add_argument("--nl", type=int, default=8, help="nlev (levels)")
    parser.add_argument("--neighbors", choices=["local", "random", "clustered"], default="local")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--warp-size", type=int, default=8,
                        help="Simulated warp size (use 8 for visualization, 32 for real)")
    parser.add_argument("--max-warps", type=int, default=6, help="Max warps to color")
    parser.add_argument("-o", "--output", default=None)
    args = parser.parse_args()

    NP, NL = args.np, args.nl
    ix = make_neighbors(NP, args.neighbors, args.seed)

    warps_A = get_warp_accesses(NP, NL, ix, "A", args.warp_size)
    warps_B = get_warp_accesses(NP, NL, ix, "B", args.warp_size)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    wc = visualize(NP, NL, ix, "A", warps_A, axes[0],
                   f"Layout A [nproma, nlev]\nkin[jc + jk·NP]   (warp → horizontal jc band)",
                   args.max_warps)
    visualize(NP, NL, ix, "B", warps_B, axes[1],
              f"Layout B [nlev, nproma]\nkin[jk + jc·NL]   (warp → vertical jk band)",
              args.max_warps)

    # Legend
    legend_items = [
        mpatches.Patch(color=(0.92, 0.92, 0.92), label="Not accessed"),
        mpatches.Patch(color=(1.0, 0.4, 0.4), label="Multi-warp overlap"),
        plt.Line2D([0], [0], marker="x", color="black", linestyle="", markersize=5,
                   label="Thread position"),
    ]
    for i in range(min(args.max_warps, 6)):
        legend_items.append(mpatches.Patch(color=wc[i], label=f"Warp {i}"))
    fig.legend(handles=legend_items, loc="lower center", ncol=5, fontsize=8,
              bbox_to_anchor=(0.5, -0.02))

    # Cache line stats
    cl_A, total_A = count_cachelines(warps_A, NP, NL, "A")
    cl_B, total_B = count_cachelines(warps_B, NP, NL, "B")
    stats = (f"Neighbors: {args.neighbors} | Warp size: {args.warp_size}\n"
             f"Layout A: avg {np.mean(cl_A):.1f} cache lines/warp, {total_A} total unique\n"
             f"Layout B: avg {np.mean(cl_B):.1f} cache lines/warp, {total_B} total unique")
    fig.text(0.5, 0.01, stats, ha="center", fontsize=8, family="monospace",
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)

    if args.output:
        fig.savefig(args.output, dpi=200, bbox_inches="tight")
        print(f"Saved to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()