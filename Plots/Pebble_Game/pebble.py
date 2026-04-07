import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np

# ===========================================================================
# Shared style — all constants in one place
# ===========================================================================
color_A      = '#F0A050'
color_B      = '#70C070'
color_C      = '#B8A0D0'
color_A_fill = '#FDF0E0'
color_B_fill = '#E4F5E4'
color_C_fill = '#EDE4F5'

# Dominator-set backgrounds
bg_V1       = '#EDE5D8';  bg_V1_edge = '#B8A88C'
bg_V2       = '#F0D8D8';  bg_V2_edge = '#C49898'
bg_cV1      = '#E8DDD0';  bg_cV1_edge = '#A8977C'
bg_cV2      = '#F2D0D5';  bg_cV2_edge = '#C08890'

# Consistent sizes
NODE_FS   = 9        # node label font size (all plots)
LABEL_FS  = 13        # V1/V2 label font size
ARROW_LW  = 1.0      # arrow line width (all plots)
NODE_LW   = 0.8      # node circle line width
BORDER_LW = 0.9      # background border line width
ARROW_MS  = 8        # arrow mutation scale

FIGSIZE_A = (6, 1.8)   # plot A
FIGSIZE_BC = (3.5, 1.8)  # plots B and C — same size

def draw_node(ax, cx, cy, r, label, face_color, edge_color, alpha=0.60):
    circle = mpatches.Circle(
        (cx, cy), r,
        facecolor=face_color, edgecolor=edge_color,
        linewidth=NODE_LW, alpha=alpha, zorder=3
    )
    ax.add_patch(circle)
    ax.text(cx, cy, label, ha='center', va='center',
            fontsize=NODE_FS, fontweight='normal', color='#333', zorder=4,
            linespacing=1.2)

def draw_arrow(ax, x0, y0, x1, y1, color='#555'):
    arrow = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle='->', mutation_scale=ARROW_MS,
        color=color, lw=ARROW_LW, alpha=0.7,
        shrinkA=0, shrinkB=0, zorder=2
    )
    ax.add_patch(arrow)

def arrow_pts(src_xy, dst_xy, r_src, r_dst):
    dx = dst_xy[0] - src_xy[0]
    dy = dst_xy[1] - src_xy[1]
    d = np.hypot(dx, dy)
    if d < 1e-9:
        return src_xy, dst_xy
    ux, uy = dx / d, dy / d
    return ((src_xy[0] + ux * r_src, src_xy[1] + uy * r_src),
            (dst_xy[0] - ux * r_dst, dst_xy[1] - uy * r_dst))

def draw_bg_box(ax, x0, y0, w, h, facecolor, edgecolor, rounding=0.10):
    bg = mpatches.FancyBboxPatch(
        (x0, y0), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={rounding}",
        facecolor=facecolor, edgecolor='none', alpha=0.15, zorder=0
    )
    ax.add_patch(bg)
    border = mpatches.FancyBboxPatch(
        (x0, y0), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={rounding}",
        facecolor='none', edgecolor=edgecolor,
        linewidth=BORDER_LW, linestyle='--', alpha=0.55, zorder=1
    )
    ax.add_patch(border)

# ===========================================================================
# PLOT A — Unit-granular CDAG
# ===========================================================================
fig_a, ax_a = plt.subplots(1, 1, figsize=FIGSIZE_A)

groups_a = [
    {'A': 'A[0, 0]', 'B': 'B[0, 0]', 'C': 'C[0, 0]'},
    {'A': 'A[0, 1]', 'B': 'B[0, 1]', 'C': 'C[0, 1]'},
    {'A': 'A[1, 0]', 'B': 'B[1, 0]', 'C': 'C[1, 0]'},
    {'A': 'A[1, 1]', 'B': 'B[1, 1]', 'C': 'C[1, 1]'},
]

n_a = len(groups_a)
gw_a, gap_a = 1.2, 0.28
tw_a = n_a * gw_a + (n_a - 1) * gap_a
r_a = 0.28
top_a, bot_a = 1.15, 0.42
dx_a = 0.38

for i, g in enumerate(groups_a):
    gcx = i * (gw_a + gap_a) + gw_a / 2
    axy = (gcx - dx_a, top_a)
    bxy = (gcx + dx_a, top_a)
    cxy = (gcx, bot_a)
    draw_node(ax_a, *axy, r_a, g['A'], color_A_fill, color_A)
    draw_node(ax_a, *bxy, r_a, g['B'], color_B_fill, color_B)
    draw_node(ax_a, *cxy, r_a, g['C'], color_C_fill, color_C)
    p0, p1 = arrow_pts(axy, cxy, r_a, r_a)
    draw_arrow(ax_a, *p0, *p1)
    p0, p1 = arrow_pts(bxy, cxy, r_a, r_a)
    draw_arrow(ax_a, *p0, *p1)

ax_a.set_xlim(-0.15, tw_a + 0.15)
ax_a.set_ylim(0.0, 1.58)
ax_a.set_aspect('equal')
ax_a.axis('off')
fig_a.tight_layout()
fig_a.savefig('plot_a_unit_granular.png', bbox_inches='tight', dpi=200)
fig_a.savefig('plot_a_unit_granular.pdf', bbox_inches='tight', dpi=200)
print("Done: pplot_a")

# ===========================================================================
# PLOT B — Row-major blocked CDAG with V1/V2 dominator sets
# ===========================================================================
fig_b, ax_b = plt.subplots(1, 1, figsize=FIGSIZE_BC)

groups_b = [
    {'A': 'A[0, 0]\nA[0, 1]', 'B': 'B[0, 0]\nB[0, 1]', 'C': 'C[0, 0]\nC[0, 1]'},
    {'A': 'A[1, 0]\nA[1, 1]', 'B': 'B[1, 0]\nB[1, 1]', 'C': 'C[1, 0]\nC[1, 1]'},
]

n_b = len(groups_b)
gw_b, gap_b = 1.5, 0.35
tw_b = n_b * gw_b + (n_b - 1) * gap_b
r_b = 0.32
top_b, bot_b = 1.35, 0.46
dx_b = 0.45
pad = 0.08

group_positions_b = []
for i, g in enumerate(groups_b):
    gcx = i * (gw_b + gap_b) + gw_b / 2
    group_positions_b.append({
        'A': (gcx - dx_b, top_b),
        'B': (gcx + dx_b, top_b),
        'C': (gcx, bot_b),
    })

# V1 background — first group inputs
g0 = group_positions_b[0]
v1_xs = [g0['A'][0], g0['B'][0]]
v1_x0 = min(v1_xs) - r_b - pad
v1_y0 = top_b - r_b - pad
v1_w  = max(v1_xs) - min(v1_xs) + 2 * r_b + 2 * pad
v1_h  = 2 * r_b + 2 * pad
draw_bg_box(ax_b, v1_x0, v1_y0, v1_w, v1_h, bg_V1, bg_V1_edge)
ax_b.text(v1_x0 + v1_w - 0.06, v1_y0 - 0.06, '$V_1$',
          ha='right', va='top', fontsize=LABEL_FS, color=bg_V1_edge)

# V2 background — first group output
v2_x0 = g0['C'][0] - r_b - pad
v2_y0 = bot_b - r_b - pad
v2_w  = 2 * r_b + 2 * pad
v2_h  = 2 * r_b + 2 * pad
draw_bg_box(ax_b, v2_x0, v2_y0, v2_w, v2_h, bg_V2, bg_V2_edge)
ax_b.text(v2_x0 + v2_w + 0.06, v2_y0 + v2_h / 2, '$V_2$',
          ha='left', va='center', fontsize=LABEL_FS, color=bg_V2_edge)

for i, g in enumerate(groups_b):
    pos = group_positions_b[i]
    axy, bxy, cxy = pos['A'], pos['B'], pos['C']
    draw_node(ax_b, *axy, r_b, g['A'], color_A_fill, color_A)
    draw_node(ax_b, *bxy, r_b, g['B'], color_B_fill, color_B)
    draw_node(ax_b, *cxy, r_b, g['C'], color_C_fill, color_C)
    p0, p1 = arrow_pts(axy, cxy, r_b, r_b)
    draw_arrow(ax_b, *p0, *p1)
    p0, p1 = arrow_pts(bxy, cxy, r_b, r_b)
    draw_arrow(ax_b, *p0, *p1)

all_xs_b = [pos[k][0] for pos in group_positions_b for k in ('A','B','C')]

# Precompute Plot C's x-extent so B and C get identical coordinate ranges
# (Plot C is wider due to 4 input nodes; match it here with an invisible artist)
_in_spacing_c = 1.05
_in_total_c = 3 * _in_spacing_c
_c_xs = [(k - 1.5) * _in_spacing_c + _in_total_c / 2 + 0.3 for k in range(4)]
_c_out_xs = [(k - 0.5) * 1.8 + _in_total_c / 2 + 0.3 for k in range(2)]
_all_c_x = _c_xs + _c_out_xs
shared_xlim = (min(_all_c_x) - r_b - pad - 0.15, max(_all_c_x) + r_b + pad + 0.15)

ax_b.set_xlim(*shared_xlim)
ax_b.set_ylim(v2_y0 - 0.10, v1_y0 + v1_h + 0.10)
ax_b.set_aspect('equal')
ax_b.axis('off')
fig_b.tight_layout()
fig_b.savefig('plot_b_row_major_blocked.png', bbox_inches='tight', dpi=200)
fig_b.savefig('plot_b_row_major_blocked.pdf', bbox_inches='tight', dpi=200)
print("Done: pplot_b")

# ===========================================================================
# PLOT C — Conflicting layout: A col-major, B col-major, C row-major
# ===========================================================================
fig_c, ax_c = plt.subplots(1, 1, figsize=FIGSIZE_BC)

r_c = 0.32
top_c = 1.35
bot_c = 0.46

input_nodes = [
    {'label': 'A[0, 0]\nA[1, 0]', 'fill': color_A_fill, 'edge': color_A, 'id': 'A0'},
    {'label': 'B[0, 0]\nB[1, 0]', 'fill': color_B_fill, 'edge': color_B, 'id': 'B0'},
    {'label': 'B[0, 1]\nB[1, 1]', 'fill': color_B_fill, 'edge': color_B, 'id': 'B1'},
    {'label': 'A[0, 1]\nA[1, 1]', 'fill': color_A_fill, 'edge': color_A, 'id': 'A1'},
]
output_nodes = [
    {'label': 'C[0, 0]\nC[0, 1]', 'fill': color_C_fill, 'edge': color_C, 'id': 'C0'},
    {'label': 'C[1, 0]\nC[1, 1]', 'fill': color_C_fill, 'edge': color_C, 'id': 'C1'},
]

edges = [
    ('A0', 'C0'), ('A0', 'C1'),
    ('B0', 'C0'), ('B0', 'C1'),
    ('B1', 'C0'), ('B1', 'C1'),
    ('A1', 'C0'), ('A1', 'C1'),
]

in_spacing = 1.05
in_total = (len(input_nodes) - 1) * in_spacing
input_positions = {}
for k, nd in enumerate(input_nodes):
    x = (k - (len(input_nodes) - 1) / 2) * in_spacing + in_total / 2 + 0.3
    input_positions[nd['id']] = (x, top_c)

out_spacing = 1.8
output_positions = {}
for k, nd in enumerate(output_nodes):
    x = (k - (len(output_nodes) - 1) / 2) * out_spacing + in_total / 2 + 0.3
    output_positions[nd['id']] = (x, bot_c)

all_pos = {**input_positions, **output_positions}
all_input_x  = [input_positions[nd['id']][0] for nd in input_nodes]
all_output_x = [output_positions[nd['id']][0] for nd in output_nodes]

# V1 background
v1_x0 = min(all_input_x) - r_c - pad
v1_y0 = top_c - r_c - pad
v1_w  = max(all_input_x) - min(all_input_x) + 2 * r_c + 2 * pad
v1_h  = 2 * r_c + 2 * pad
draw_bg_box(ax_c, v1_x0, v1_y0, v1_w, v1_h, bg_cV1, bg_cV1_edge, rounding=0.12)
ax_c.text(v1_x0 + v1_w - 0.06, v1_y0 - 0.06, '$V_1$',
          ha='right', va='top', fontsize=LABEL_FS, color=bg_cV1_edge)

# V2 background
v2_x0 = min(all_output_x) - r_c - pad
v2_y0 = bot_c - r_c - pad
v2_w  = 2 * r_c + 2 * pad
v2_h  = 2 * r_c + 2 * pad
draw_bg_box(ax_c, v2_x0, v2_y0, v2_w, v2_h, bg_cV2, bg_cV2_edge, rounding=0.12)
ax_c.text(v2_x0 + v2_w + 0.06, v2_y0 + v2_h / 2, '$V_2$',
          ha='left', va='center', fontsize=LABEL_FS, color=bg_cV2_edge)

for nd in input_nodes:
    px, py = input_positions[nd['id']]
    draw_node(ax_c, px, py, r_c, nd['label'], nd['fill'], nd['edge'])

for nd in output_nodes:
    px, py = output_positions[nd['id']]
    draw_node(ax_c, px, py, r_c, nd['label'], nd['fill'], nd['edge'])

for src_id, dst_id in edges:
    sxy = all_pos[src_id]
    dxy = all_pos[dst_id]
    p0, p1 = arrow_pts(sxy, dxy, r_c, r_c)
    draw_arrow(ax_c, *p0, *p1)

all_x = all_input_x + all_output_x
lx = min(all_x) - r_c - pad
rx = max(all_x) + r_c + pad
ax_c.set_xlim(lx - 0.15, rx + 0.15)
# ylim matches Plot B (same top/bot/r/pad)
ax_c.set_ylim(v2_y0 - 0.10, v1_y0 + v1_h + 0.10)
ax_c.set_aspect('equal')
ax_c.axis('off')
fig_c.tight_layout()
fig_c.savefig('plot_c_conflicting_layout.png', bbox_inches='tight', dpi=200)
fig_c.savefig('plot_c_conflicting_layout.pdf', bbox_inches='tight', dpi=200)
print("Done: pplot_c")