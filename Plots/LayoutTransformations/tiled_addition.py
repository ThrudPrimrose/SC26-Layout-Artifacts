import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

import layout_stages as ls

rainbow_colors = {
    'red':    '#e8908a', 'orange': '#eab070',
    'yellow': '#e8d870', 'green':  '#82c482',
    'blue':   '#80a8d8', 'purple': '#a080d8',
    'teal': '#82c4b4',  'pink':   '#e88ab8',
}

def tiled_addition(name, tile_row_colors, grid):
    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.8))
    nr, nc = 4, 4

    for r in range(nr):
        lc, rc = tile_row_colors[r]
        ax.add_patch(plt.Rectangle((0, (nr-1-r)*ls.cell), 2*ls.cell, ls.cell,
            facecolor=rainbow_colors[lc], edgecolor='none'))
        ax.add_patch(plt.Rectangle((2*ls.cell, (nr-1-r)*ls.cell), 2*ls.cell, ls.cell,
            facecolor=rainbow_colors[rc], edgecolor='none'))

    ls.dashed_grid(ax, nr, nc)

    for r in range(nr):
        for c in range(nc):
            lc, rc = tile_row_colors[r]
            color = rainbow_colors[lc if c < 2 else rc]
            x0, y0 = ls.gx(c) - ls.cell/2, ls.gy(r, nr) - ls.cell/2
            ax.add_patch(plt.Rectangle((x0, y0), ls.cell, ls.cell,
                facecolor='#ececec', edgecolor='none'))
            ax.add_patch(mpatches.FancyBboxPatch(
                (x0+ls.off, y0+ls.off), ls.inner, ls.inner,
                boxstyle="round,pad=0.02,rounding_size=0.06",
                facecolor=color, edgecolor='none'))
            ax.text(ls.gx(c), ls.gy(r, nr), str(grid[r][c]),
                    ha='center', va='center', fontsize=9, color='#222')

    ls.finish(ax, nc, nr, -0.15, name)

if __name__ == "__main__":
    ## 1 - original
    tile_row_colors = [
        ('pink', 'red'), ('orange', 'yellow'),
        ('green', 'teal'), ('blue', 'purple'),
    ]
    grid = [[0,2,8,10],[1,3,9,11],[4,6,12,14],[5,7,13,15]]
    tiled_addition("tiled_add_1", tile_row_colors, grid)

    ## 2 - after permute
    tile_row_colors = [
        ('pink', 'red'), ('orange', 'yellow'),
        ('green', 'teal'), ('blue', 'purple'),
    ]
    grid = [[0,1,4,5],[2,3,6,7],[8,9,12,13],[10,11,14,15]]
    tiled_addition("tiled_add_2", tile_row_colors, grid)

    # 3 - after blocking
    tile_row_colors = [
        ('pink', 'orange'), ('red', 'yellow'),
        ('green', 'blue'), ('teal', 'purple'),
    ]
    grid = [[0,1,4,5],[2,3,6,7],[8,9,12,13],[10,11,14,15]]
    tiled_addition("tiled_add_3", tile_row_colors, grid)