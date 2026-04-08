import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

colors = {
    'red':    '#e8908a',
    'green':  '#82c482',
    'orange': '#eab070',
    'blue':   '#80a8d8',
}
bg = {
    'red':    '#f8cecc',
    'green':  '#d5e8d4',
    'orange': '#ffe6cc',
    'blue':   '#dae8fc',
}
gray_fill = '#e0e0e0'
group_names = ['red', 'green', 'orange', 'blue']

def is_pad(eid, npg):
    return npg == 4 and (eid % npg == 3)

def elem_color(eid, npg):
    if is_pad(eid, npg): return gray_fill
    return colors[group_names[eid // npg]]

def elem_tc(eid, npg):
    return '#999' if is_pad(eid, npg) else '#222'

def elem_label(eid, npg):
    return '*' if is_pad(eid, npg) else str(eid)

cell = 0.8
inner = 0.66
off = (cell - inner) / 2

def draw_cell(ax, cx, cy, eid, npg):
    x0, y0 = cx - cell/2, cy - cell/2
    ax.add_patch(plt.Rectangle((x0, y0), cell, cell,
        facecolor='#ececec', edgecolor='none'))
    ax.add_patch(mpatches.FancyBboxPatch(
        (x0+off, y0+off), inner, inner,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        facecolor=elem_color(eid, npg), edgecolor='none'))
    ax.text(cx, cy, elem_label(eid, npg), ha='center', va='center',
            fontsize=9, color=elem_tc(eid, npg))

def gx(c, x_off=0): return x_off + c * cell + cell/2
def gy(r, nr): return (nr-1-r) * cell + cell/2

def dashed_grid(ax, nr, nc, x_off=0):
    for c in range(nc+1):
        ax.plot([x_off+c*cell]*2, [0, nr*cell], ls='--', color='#777', lw=2, dashes=(3,3))
    for r in range(nr+1):
        ax.plot([x_off, x_off+nc*cell], [r*cell]*2, ls='--', color='#777', lw=2, dashes=(3,3))

def dashed_rows(ax, nr, nc, x_off=0):
    for c in [0, nc]:
        ax.plot([x_off+c*cell]*2, [0, nr*cell], ls='--', color='#777', lw=2, dashes=(3,3))
    for r in range(nr+1):
        ax.plot([x_off, x_off+nc*cell], [r*cell]*2, ls='--', color='#777', lw=2, dashes=(3,3))

def dashed_blocks(ax, nr, nc, x_off=0):
    for c in [0, nc//2, nc]:
        ax.plot([x_off+c*cell]*2, [0, nr*cell], ls='--', color='#777', lw=2, dashes=(3,3))
    for r in [0, nr//2, nr]:
        ax.plot([x_off, x_off+nc*cell], [r*cell]*2, ls='--', color='#777', lw=2, dashes=(3,3))

def draw_cache_outlines_4x3(ax, x_off=0):
    """Draw L/Z/L dashed outlines for cache lines (B=4) on a 4×3 grid.
    Row-major: CL0=addr0-3, CL1=addr4-7, CL2=addr8-11.
    CL0: (0,0)(0,1)(0,2)(1,0) → L
    CL1: (1,1)(1,2)(2,0)(2,1) → Z
    CL2: (2,2)(3,0)(3,1)(3,2) → L inverted"""
    c = cell
    cl_colors = ['#777', '#777', '#777']

    # CL0: L shape
    poly0 = [(0,4*c),(3*c,4*c),(3*c,3*c),(1*c,3*c),(1*c,2*c),(0,2*c)]
    # CL1: Z shape
    poly1 = [(1*c,3*c),(3*c,3*c),(3*c,2*c),(2*c,2*c),(2*c,1*c),(0,1*c),(0,2*c),(1*c,2*c)]
    # CL2: inverted L
    poly2 = [(2*c,2*c),(3*c,2*c),(3*c,0),(0,0),(0,1*c),(2*c,1*c)]

    for poly, color in zip([poly0, poly1, poly2], cl_colors):
        shifted = [(x + x_off, y) for x, y in poly]
        p = ax.add_patch(plt.Polygon(shifted, closed=True, fill=False,
            edgecolor=color, linewidth=1.5, linestyle='--', zorder=5))
        p.set_linestyle((0, (4, 3)))  # explicit dash pattern

def bg_cols(ax, nc, nr, names, x_off=0):
    for c, n in enumerate(names):
        ax.add_patch(plt.Rectangle((x_off+c*cell, 0), cell, nr*cell,
            facecolor=bg[n], edgecolor='none'))

def bg_rows(ax, nc, nr, names, x_off=0):
    for r, n in enumerate(names):
        ax.add_patch(plt.Rectangle((x_off, (nr-1-r)*cell), nc*cell, cell,
            facecolor=bg[n], edgecolor='none'))

def bg_tiles(ax, nr, tiles, x_off=0):
    ts = 2*cell
    for (tr, tc), n in tiles:
        ax.add_patch(plt.Rectangle((x_off+tc*ts, (1-tr)*ts), ts, ts,
            facecolor=bg[n], edgecolor='none'))

def ribbon(ax, x0, ry, order, nc, npg):
    n = len(order)
    rh, ew, B = 0.5, nc*cell/n, 2
    ax.text(x0 + nc*cell/2, ry+rh+0.22, 'Linearized memory',
            ha='center', va='center', fontsize=6.5, fontstyle='italic', color='#555')
    for i, eid in enumerate(order):
        x = x0 + i*ew
        ax.add_patch(mpatches.FancyBboxPatch(
            (x+0.008, ry+0.025), ew-0.016, rh-0.05,
            boxstyle="round,pad=0.008,rounding_size=0.025",
            facecolor=elem_color(eid, npg), edgecolor='#bbb', linewidth=0.3))
        ax.text(x+ew/2, ry+rh/2, elem_label(eid, npg), ha='center', va='center',
                fontsize=4.5, color=elem_tc(eid, npg))
    for b in range(0, n+1, B):
        x = x0 + b*ew
        ax.plot([x,x], [ry-0.02, ry+rh+0.02], ls='--', color='#aaa', lw=0.8)
    for b in range(n//B):
        x = x0 + (b*B + B/2)*ew
        ax.text(x, ry-0.15, f'B$_{b}$', ha='center', va='center', fontsize=5, color='#444')

def draw_cells(ax, grid, nr, nc, npg, x_off=0):
    for r in range(nr):
        for c in range(nc):
            draw_cell(ax, gx(c, x_off), gy(r, nr), grid[r][c], npg)

def finish(ax, nc, nr, ry, name):
    ax.set_xlim(-0.15, nc*cell+0.15)
    ax.set_ylim(ry-0.45, nr*cell+0.15)
    ax.set_aspect('equal'); ax.axis('off')
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    plt.savefig(f'{name}.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(f'{name}.png', bbox_inches='tight', dpi=300)
    plt.close(); print(f"Done: {name}")


# ── Stage 1: Original row-major 4×3 ──────────────────────────────────
# Padded numbering: row r, col c → eid = r*4 + c  (skipping padding slots)
# Grid: 0,1,2 / 4,5,6 / 8,9,10 / 12,13,14
# Colors: 0-3=red, 4-7=green, 8-11=orange, 12-15=blue (npg=4)
# Cache outlines: L/Z/L showing misalignment
def stage1(name, no_ribbon=True):
    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.8))
    nr, nc, npg = 4, 3, 4
    # Light neutral background
    ax.add_patch(plt.Rectangle((0, 0), nc*cell, nr*cell,
        facecolor='#f5f5f5', edgecolor='none'))
    grid = [[0,1,2],[4,5,6],[8,9,10],[12,13,14]]
    draw_cells(ax, grid, nr, nc, npg)
    draw_cache_outlines_4x3(ax)
    ry = -0.9 if not no_ribbon else -0.15
    if not no_ribbon:
        ro = [0,1,2,4,5,6,8,9,10,12,13,14]
        ribbon(ax, 0, -0.9, ro, nc, npg)
    finish(ax, nc, nr, ry, name)


# ── Stage 2: After pad → 4×4 ─────────────────────────────────────────
# Padding added: elements 3,7,11,15 (shown as *)
# Cache lines now aligned with rows
def stage2(name, no_ribbon=True):
    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.8))
    nr, nc, npg = 4, 4, 4
    bg_rows(ax, nc, nr, ['red','green','orange','blue'])
    dashed_grid(ax, nr, nc)
    draw_cells(ax, [[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]], nr, nc, npg)
    ry = -0.9 if not no_ribbon else -0.15
    if not no_ribbon:
        ribbon(ax, 0, -0.9, list(range(16)), nc, npg)
    finish(ax, nc, nr, ry, name)


# ── Stage 3: After permute (transpose) ───────────────────────────────
def stage3(name, no_ribbon=True):
    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.8))
    nr, nc, npg = 4, 4, 4
    bg_cols(ax, nc, nr, ['red','green','orange','blue'])
    dashed_grid(ax, nr, nc)
    draw_cells(ax, [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]], nr, nc, npg)
    ry = -0.9 if not no_ribbon else -0.15
    if not no_ribbon:
        ribbon(ax, 0, -0.9, [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15], nc, npg)
    finish(ax, nc, nr, ry, name)


# ── Stage 4: After block (2×2 tiles) ─────────────────────────────────
def stage4(name, no_ribbon=True):
    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.8))
    nr, nc, npg = 4, 4, 4
    bg_tiles(ax, nr, [((0,0),'red'),((0,1),'green'),((1,0),'orange'),((1,1),'blue')])
    dashed_grid(ax, nr, nc)
    grid = [[0,1,4,5],[2,3,6,7],[8,9,12,13],[10,11,14,15]]
    draw_cells(ax, grid, nr, nc, npg)
    ry = -0.9 if not no_ribbon else -0.15
    if not no_ribbon:
        ro = [grid[r][c] for r in range(nr) for c in range(nc)]
        ribbon(ax, 0, -0.9, ro, nc, npg)
    finish(ax, nc, nr, ry, name)


# ── Stage 5: After shuffle ───────────────────────────────────────────
def stage5(name, no_ribbon=True):
    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.8))
    nr, nc, npg = 4, 4, 4
    ax.add_patch(plt.Rectangle((0,0), nc*cell, nr*cell,
        facecolor='#f9f9f9', edgecolor='none'))
    dashed_grid(ax, nr, nc)
    grid = [[0,4,8,12],[5,1,13,9],[10,14,2,6],[15,11,7,3]]
    draw_cells(ax, grid, nr, nc, npg)
    ry = -0.9 if not no_ribbon else -0.15
    if not no_ribbon:
        ro = [grid[r][c] for r in range(nr) for c in range(nc)]
        ribbon(ax, 0, -0.9, ro, nc, npg)
    finish(ax, nc, nr, ry, name)


# ── Stage 6: Combined figure ─────────────────────────────────────────
def stage6(name):
    nr = 4
    npg_list = [4, 4, 4, 4, 4]
    nc_list  = [3, 4, 4, 4, 4]
    grids = [
        [[0,1,2],[4,5,6],[8,9,10],[12,13,14]],
        [[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]],
        [[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]],
        [[0,1,4,5],[2,3,6,7],[8,9,12,13],[10,11,14,15]],
        [[0,4,8,12],[5,1,13,9],[10,14,2,6],[15,11,7,3]],
    ]
    arrow_labels = ['Pad', 'Permute', 'Block', 'Shuffle']
    arrow_gap = 1.45

    x_offsets = [0.0]
    for i in range(4):
        x_offsets.append(x_offsets[i] + nc_list[i] * cell + arrow_gap)

    total_w = x_offsets[-1] + nc_list[-1] * cell
    fig, ax = plt.subplots(1, 1, figsize=(total_w * 0.55, nr * cell * 0.55 + 0.3))

    for idx in range(5):
        xo = x_offsets[idx]
        nc = nc_list[idx]
        npg = npg_list[idx]
        grid = grids[idx]

        if idx == 0:
            ax.add_patch(plt.Rectangle((xo, 0), nc*cell, nr*cell,
                facecolor='#f5f5f5', edgecolor='none'))
            draw_cache_outlines_4x3(ax, x_off=xo)
        elif idx == 1:
            bg_rows(ax, nc, nr, ['red','green','orange','blue'], x_off=xo)
            dashed_rows(ax, nr, nc, x_off=xo)
        elif idx == 2:
            bg_cols(ax, nc, nr, ['red','green','orange','blue'], x_off=xo)
            dashed_rows(ax, nr, nc, x_off=xo)
        elif idx == 3:
            bg_tiles(ax, nr, [((0,0),'red'),((0,1),'green'),
                               ((1,0),'orange'),((1,1),'blue')], x_off=xo)
            dashed_blocks(ax, nr, nc, x_off=xo)
        elif idx == 4:
            ax.add_patch(plt.Rectangle((xo, 0), nc*cell, nr*cell,
                facecolor='#f9f9f9', edgecolor='none'))
            dashed_grid(ax, nr, nc, x_off=xo)

        draw_cells(ax, grid, nr, nc, npg, x_off=xo)

    for i in range(4):
        x_start = x_offsets[i] + nc_list[i] * cell + 0.12
        x_end   = x_offsets[i+1] - 0.12
        y_mid   = nr * cell / 2

        ax.annotate('',
            xy=(x_end, y_mid), xytext=(x_start, y_mid),
            arrowprops=dict(
                arrowstyle='->,head_width=0.35,head_length=0.2',
                color='#555', lw=2.5))

        mx = (x_start + x_end) / 2
        ax.text(mx, y_mid + 0.35, arrow_labels[i],
                ha='center', va='center', fontsize=11, fontstyle='italic',
                color='#333',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                          edgecolor='#999', linewidth=0.8))

    ax.set_xlim(-0.15, total_w + 0.15)
    ax.set_ylim(-0.15, nr*cell + 0.15)
    ax.set_aspect('equal'); ax.axis('off')
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(f'{name}.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(f'{name}.png', bbox_inches='tight', dpi=300)
    plt.close(); print(f"Done: {name}")

if __name__ == "__main__":
    stage1('layout_1_original')
    stage2('layout_2_pad')
    stage3('layout_3_permute')
    stage4('layout_4_block')
    stage5('layout_5_shuffle')
    stage6('layout_6_combined')