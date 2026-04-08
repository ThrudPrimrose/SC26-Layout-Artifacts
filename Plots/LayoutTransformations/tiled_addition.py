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

tile_row_colors = [
    ('red', 'pink'), ('orange', 'yellow'),
    ('green', 'teal'), ('blue', 'purple'),
]

arrow_colors = {
    'dark gray':'#222222',
    'gray':'#888888',
    'black': '#000000'
}

def dashed_rows(ax, nr, nc, x_off=0):
    for c in [0, nc]:
        ax.plot([x_off+c*ls.cell]*2, [0, nr*ls.cell], ls='--', color='#666', lw=1.8, dashes=(3,3))
    for r in range(nr+1):
        ax.plot([x_off, x_off+nc*ls.cell], [r*ls.cell]*2, ls='--', color='#666', lw=1.8, dashes=(3,3))

def dashed_blocks(ax, nr, nc, x_off=0):
    for c in [0, nc//2, nc]:
        ax.plot([x_off+c*ls.cell]*2, [0, nr*ls.cell], ls='--', color='#666', lw=1.8, dashes=(3,3))
    for r in [0, nr//2, nr]:
        ax.plot([x_off, x_off+nc*ls.cell], [r*ls.cell]*2, ls='--', color='#666', lw=1.8, dashes=(3,3))

def box_point(r, c, nr, dx, dy):
    return ls.gx(c) + dx * ls.inner/2, ls.gy(r, nr) + dy * ls.inner/2

def straight_arrow(ax, x0, y0, x1, y1, color, lw=1.2):
    lw=2.0
    head_width=0.3
    head_length=0.52
    alpha=1.0

    if color == arrow_colors['black']:
        lw = 2.2
        head_width=0.3
        head_length=0.8

    ax.annotate('',
        xy=(x1, y1), xytext=(x0, y0),
        arrowprops=dict(
            arrowstyle=f'->,head_width={head_width},head_length={head_length}',
            color=color, lw=lw, connectionstyle='arc3,rad=0.0',
            shrinkA=0.8, shrinkB=0,
            alpha=alpha
        ),
        zorder=4,
    )

def tiled_addition(name, tile_row_colors, grid, arrow_specs, order='row-major'):
    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.8))
    nr, nc = 4, 4

    for r in range(nr):
        lc, rc = tile_row_colors[r]
        ax.add_patch(plt.Rectangle((0, (nr-1-r)*ls.cell), 2*ls.cell, ls.cell,
            facecolor=rainbow_colors[lc], edgecolor='none'))
        ax.add_patch(plt.Rectangle((2*ls.cell, (nr-1-r)*ls.cell), 2*ls.cell, ls.cell,
            facecolor=rainbow_colors[rc], edgecolor='none'))

    for r in range(nr):
        for c in range(nc):
            lc, rc = tile_row_colors[r]
            color = rainbow_colors[lc if c < 2 else rc]
            x0, y0 = ls.gx(c) - ls.cell/2, ls.gy(r, nr) - ls.cell/2
            ax.add_patch(plt.Rectangle((x0, y0), ls.cell, ls.cell,
                facecolor='#ececec', edgecolor='none', zorder=1))
            ax.add_patch(mpatches.FancyBboxPatch(
                (x0+ls.off, y0+ls.off), ls.inner, ls.inner,
                boxstyle="round,pad=0.02,rounding_size=0.06",
                facecolor=color, edgecolor='none', zorder=4))
            ax.text(ls.gx(c), ls.gy(r, nr), str(grid[r][c]),
                    ha='center', va='center', fontsize=12, color='#222', zorder=5)

    # ls.dashed_grid(ax, nr, nc)

    if (order=='row-major'):
        dashed_rows(ax, nr, nc)
    elif (order=='blocked'):
        dashed_blocks(ax, nr, nc)
    else:
        raise ValueError()

    for spec in arrow_specs:
        (r0, c0), (r1, c1), direction, color = spec
        if direction == "bottom edge to top edge":
            dx0, dy0 = 0, -1
            dx1, dy1 = 0, +1
        elif direction == "top right corner to bottom left corner":
            dx0, dy0 = +1, +1
            dx1, dy1 = -1, -1
        elif direction == "bottom left corner to top right corner":
            dx0, dy0 = -1, -1
            dx1, dy1 = +1, +1
        else:
            raise ValueError(f"direction unrecognised: {direction}")

        x0, y0 = box_point(r0, c0, nr, dx0, dy0)
        x1, y1 = box_point(r1, c1, nr, dx1, dy1)
        straight_arrow(ax, x0, y0, x1, y1, color)

    ls.finish(ax, nc, nr, -0.15, name)

def memory_blocks(name):
    fig, ax = plt.subplots(1, 1, figsize=(12.8, 1.2))
    nr, nc = 2, 8
    ribbon_colors = [
        
        ['green', 'green', 'teal', 'teal', 'blue', 'blue', 'purple', 'purple'],
        ['pink', 'pink', 'red', 'red', 'orange', 'orange', 'yellow', 'yellow'],
    ]

    for r in range(nr):
        for c in range(nc):
            color_name = ribbon_colors[r][c]
    
            x0 = c * ls.cell
            y0 = r * ls.cell
            ax.add_patch(plt.Rectangle((x0, y0), ls.cell, ls.cell,
                facecolor='#ececec', edgecolor='none'))
            ax.add_patch(mpatches.FancyBboxPatch(
                (x0+ls.off, y0+ls.off), ls.inner, ls.inner,
                boxstyle="round,pad=0.02,rounding_size=0.06",
                facecolor=rainbow_colors[color_name], edgecolor='none', zorder=4))

    for c in range(nc+1):
        lw = 2.0 if c % 2 == 0 else 1.0
        col = '#999' if c % 2 == 0 else '#999'
        ax.plot([c*ls.cell, c*ls.cell], [0, (nr)*ls.cell], ls='--', color=col, lw=lw, dashes=(2,2), zorder=2)
    ax.plot([0, nc*ls.cell], [0, 0],         ls='--', color='#999', lw=2.0, dashes=(2,2), zorder=2)
    ax.plot([0, nc*ls.cell], [ls.cell, ls.cell], ls='--', color='#999', lw=2.0, dashes=(2,2), zorder=2)
    ax.plot([0, nc*ls.cell], [nr*ls.cell, nr*ls.cell], ls='--', color='#999', lw=2.0, dashes=(2,2), zorder=2)

    ls.finish(ax, nc, nr, -0.15, name)

if __name__ == "__main__":
    tile_row_colors = {
        0: [
            ('pink', 'red'), ('orange', 'yellow'),
            ('green', 'teal'), ('blue', 'purple'),
        ],
        1: [
            ('pink', 'red'), ('orange', 'yellow'),
            ('green', 'teal'), ('blue', 'purple'),
        ],
        2: [
            ('pink', 'orange'), ('red', 'yellow'),
            ('green', 'blue'), ('teal', 'purple'),
        ]
    }

    grids = {
        0: [[0,2,8,10],[1,3,9,11],[4,6,12,14],[5,7,13,15]],
        1: [[0,1,4,5],[2,3,6,7],[8,9,12,13],[10,11,14,15]],
        2: [[0,1,4,5],[2,3,6,7],[8,9,12,13],[10,11,14,15]]
    }

    arrow_specs = {
        0: [
            ((0, 0), (1, 0), "bottom edge to top edge", arrow_colors['dark gray']),
            ((1, 0), (0, 1), "top right corner to bottom left corner", arrow_colors['dark gray']),
            ((0, 1), (1, 1), "bottom edge to top edge", arrow_colors['dark gray']),
            ((1, 1), (2, 0), "bottom left corner to top right corner", arrow_colors['dark gray']),
            ((2, 0), (3, 0), "bottom edge to top edge", arrow_colors['dark gray']),
            ((3, 0), (2, 1), "top right corner to bottom left corner", arrow_colors['dark gray']),
            ((2, 1), (3, 1), "bottom edge to top edge", arrow_colors['dark gray']),
            ((3, 1), (0, 2), "top right corner to bottom left corner", arrow_colors['black']),
            ((0, 2), (1, 2), "bottom edge to top edge", arrow_colors['dark gray']),
            ((1, 2), (0, 3), "top right corner to bottom left corner", arrow_colors['dark gray']),
            ((0, 3), (1, 3), "bottom edge to top edge", arrow_colors['dark gray']),
            ((1, 3), (2, 2), "bottom left corner to top right corner", arrow_colors['dark gray']),
            ((2, 2), (3, 2), "bottom edge to top edge", arrow_colors['dark gray']),
            ((3, 2), (2, 3), "top right corner to bottom left corner", arrow_colors['dark gray']),
            ((2, 3), (3, 3), "bottom edge to top edge", arrow_colors['dark gray']),
        ],
        1: [
            ((0, 1), (1, 0), "bottom left corner to top right corner", arrow_colors['dark gray']),
            ((1, 1), (0, 2), "top right corner to bottom left corner", arrow_colors['gray']),
            ((0, 3), (1, 2), "bottom left corner to top right corner", arrow_colors['dark gray']),
            ((1, 3), (2, 0), "bottom left corner to top right corner", arrow_colors['gray']),
            ((2, 1), (3, 0), "bottom left corner to top right corner", arrow_colors['dark gray']),
            ((3, 1), (2, 2), "top right corner to bottom left corner", arrow_colors['gray']),
            ((2, 3), (3, 2), "bottom left corner to top right corner", arrow_colors['dark gray']),
        ],
        2: [
            ((0, 1), (1, 0), "bottom left corner to top right corner", arrow_colors['gray']),
            ((1, 1), (0, 2), "top right corner to bottom left corner", arrow_colors['gray']),
            ((0, 3), (1, 2), "bottom left corner to top right corner", arrow_colors['gray']),
            ((1, 3), (2, 0), "bottom left corner to top right corner", arrow_colors['gray']),
            ((2, 1), (3, 0), "bottom left corner to top right corner", arrow_colors['gray']),
            ((3, 1), (2, 2), "top right corner to bottom left corner", arrow_colors['gray']),
            ((2, 3), (3, 2), "bottom left corner to top right corner", arrow_colors['gray']),
        ]
    }

    tiled_addition("tiled_add_0", tile_row_colors[0], grids[0], arrow_specs[0], 'row-major')
    tiled_addition("tiled_add_1", tile_row_colors[1], grids[1], arrow_specs[1], 'row-major')
    tiled_addition("tiled_add_2", tile_row_colors[2], grids[2], arrow_specs[2], 'blocked')

    ## UNCOMMENT BELOW FOR MEMORY BLOCKS
    # memory_blocks('./memory_blocks2')