import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(6.0, 2.0))

# --- Color scheme (no red/blue for pebble game) ---
color_A = '#F0A050'   # amber/orange for A
color_B = '#70C070'   # green for B
color_C = '#B8A0D0'   # lavender for C

color_A_fill = '#FDF0E0'  # light amber fill
color_B_fill = '#E4F5E4'  # light green fill
color_C_fill = '#EDE4F5'  # light lavender fill

# --- Computation groups: C[i,j] = f(A[i,j], B[i,j]) ---
groups = [
    {'A': 'A[0, 0]', 'B': 'B[0, 0]', 'C': 'C[0, 0]'},
    {'A': 'A[0, 1]', 'B': 'B[0, 1]', 'C': 'C[0, 1]'},
    {'A': 'A[1, 0]', 'B': 'B[1, 0]', 'C': 'C[1, 0]'},
    {'A': 'A[1, 1]', 'B': 'B[1, 1]', 'C': 'C[1, 1]'},
]

n_groups = len(groups)
group_w = 1.2        # width per group (tighter)
gap = 0.28           # gap between groups
total_w = n_groups * group_w + (n_groups - 1) * gap

# Node geometry
node_r = 0.28        # circle radius
top_y = 1.35         # y-center of A, B nodes
bot_y = 0.45         # y-center of C node
dx = 0.38            # horizontal offset of A,B from group center

def draw_node(ax, cx, cy, r, label, face_color, edge_color, alpha=0.60):
    """Draw a circle node with label."""
    circle = mpatches.Circle(
        (cx, cy), r,
        facecolor=face_color, edgecolor=edge_color,
        linewidth=0.8, alpha=alpha, zorder=3
    )
    ax.add_patch(circle)
    ax.text(cx, cy, label, ha='center', va='center',
            fontsize=7, fontweight='normal', color='#333', zorder=4)

def draw_arrow(ax, x0, y0, x1, y1, color='#555'):
    """Draw an arrow from (x0,y0) to (x1,y1), shrunk to avoid node overlap."""
    arrow = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle='->', mutation_scale=10,
        color=color, lw=0.8, alpha=0.7,
        shrinkA=0, shrinkB=0,
        zorder=2
    )
    ax.add_patch(arrow)

for i, g in enumerate(groups):
    # Group center x
    gcx = i * (group_w + gap) + group_w / 2

    # --- Node positions ---
    ax_pos = (gcx - dx, top_y)
    bx_pos = (gcx + dx, top_y)
    cx_pos = (gcx, bot_y)

    # --- Draw nodes ---
    draw_node(ax, *ax_pos, node_r, g['A'], color_A_fill, color_A)
    draw_node(ax, *bx_pos, node_r, g['B'], color_B_fill, color_B)
    draw_node(ax, *cx_pos, node_r, g['C'], color_C_fill, color_C)

    # --- Arrows: edge exits from bottom of A/B circles, enters top of C circle ---
    # A -> C
    a_bot = (ax_pos[0] + 0.08, ax_pos[1] - node_r)
    c_tl  = (cx_pos[0] - 0.12, cx_pos[1] + node_r)
    draw_arrow(ax, a_bot[0], a_bot[1], c_tl[0], c_tl[1])

    # B -> C
    b_bot = (bx_pos[0] - 0.08, bx_pos[1] - node_r)
    c_tr  = (cx_pos[0] + 0.12, cx_pos[1] + node_r)
    draw_arrow(ax, b_bot[0], b_bot[1], c_tr[0], c_tr[1])

# --- Axes ---
ax.set_xlim(-0.2, total_w + 0.2)
ax.set_ylim(-0.05, 1.85)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig('plot_a_unit_granular.png', bbox_inches='tight', dpi=200)
plt.savefig('plot_a_unit_granular.pdf', bbox_inches='tight', dpi=200)
print("Done: plot_a_unit_granular.png / .pdf")