import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(5.6, 3.6))

cell = 1.0
inner = 0.82
offset_c = (cell - inner) / 2
white = np.array([1.0, 1.0, 1.0])

# Colors
col_i_hex = '#e8706a'
col_f_hex = '#70c070'
col_i_rgb = np.array([int(col_i_hex[i:i+2], 16) / 255 for i in (1, 3, 5)])
col_f_rgb = np.array([int(col_f_hex[i:i+2], 16) / 255 for i in (1, 3, 5)])

col_i = col_i_hex
col_i_light = '#f4b8b4'
col_f = '#70c070'
col_f_light = '#b8e0b8'

opacities = [0.95, 0.55]

def elem_color_rgb(base_rgb, idx):
    a = opacities[idx % len(opacities)]
    return a * base_rgb + (1 - a) * white


def draw_cell(ax, cx, cy, label, facecolor, subscript):
    x0 = cx - cell / 2
    y0 = cy - cell / 2
    bg = mpatches.FancyBboxPatch(
        (x0, y0), cell, cell,
        boxstyle="square,pad=0",
        facecolor='#ececec', edgecolor='none', linewidth=0
    )
    ax.add_patch(bg)
    rect = mpatches.FancyBboxPatch(
        (x0 + offset_c, y0 + offset_c), inner, inner,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        facecolor=facecolor, edgecolor='none', linewidth=0
    )
    ax.add_patch(rect)
    ax.text(cx, cy, f'{label}$_{subscript}$',
            ha='center', va='center', fontsize=12, color='#222')


def draw_outline(ax, x0, y0, w, h):
    rect = mpatches.FancyBboxPatch(
        (x0, y0), w, h,
        boxstyle="square,pad=0",
        facecolor='none', edgecolor='#777', linewidth=2.0,
        linestyle='--'
    )
    ax.add_patch(rect)


# ── Cell positions ────────────────────────────────────────────────────
# Top-left: i array
i_y = 3.0
i_x0 = 1.0
draw_cell(ax, i_x0, i_y, 'i', col_i, '0')
draw_cell(ax, i_x0 + cell, i_y, 'i', col_i_light, '1')
draw_outline(ax, i_x0 - cell/2, i_y - cell/2, 2*cell, cell)

# Bottom-left: f array
f_y = 1.4
f_x0 = 1.0
draw_cell(ax, f_x0, f_y, 'f', col_f, '0')
draw_cell(ax, f_x0 + cell, f_y, 'f', col_f_light, '1')
draw_outline(ax, f_x0 - cell/2, f_y - cell/2, 2*cell, cell)

# Right: packed array [i₀, f₀, i₁, f₁]
pack_y = 2.2
pack_x0 = 4.8
packed_colors = [col_i, col_f, col_i_light, col_f_light]
packed_labels = [('i', '0'), ('f', '0'), ('i', '1'), ('f', '1')]
for k in range(4):
    draw_cell(ax, pack_x0 + k * cell, pack_y,
              packed_labels[k][0], packed_colors[k], packed_labels[k][1])
draw_outline(ax, pack_x0 - cell/2, pack_y - cell/2, 4*cell, cell)

# ── Arrows ────────────────────────────────────────────────────────────
arrow_kw = dict(
    arrowstyle='->,head_width=0.3,head_length=0.18',
    color='#555555', lw=2.5,
    connectionstyle='arc3,rad=0'
)

# pack
a1x = i_x0 + cell + cell/2 + 0.15
a1y = i_y - 0.15
a2x = pack_x0 - cell/2 - 0.15
a2y = pack_y + 0.15
ax.annotate('', xy=(a2x, a2y), xytext=(a1x, a1y), arrowprops=arrow_kw)
mx, my = (a1x+a2x)/2, (a1y+a2y)/2
ax.text(mx+0.1 , my + 0.28, 'Zip', ha='center', va='center',
        fontsize=13, fontstyle='italic', color='#333',
        bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                  edgecolor='#999999', linewidth=0.8))

# unpack
b1x = pack_x0 - cell/2 - 0.15
b1y = pack_y - 0.15
b2x = f_x0 + cell + cell/2 + 0.15
b2y = f_y + 0.15
ax.annotate('', xy=(b2x, b2y), xytext=(b1x, b1y), arrowprops=arrow_kw)
mx2, my2 = (b1x+b2x)/2, (b1y+b2y)/2
ax.text(mx2+0.1 , my2 - 0.28, 'Unzip', ha='center', va='center',
        fontsize=13, fontstyle='italic', color='#333',
        bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                  edgecolor='#999999', linewidth=0.8))

# ── Memory ribbons (reference style) ─────────────────────────────────
"""
n_elems = 4  # total elements to show per ribbon (same length)
B = 2
ribbon_h = 0.5

# Both ribbons same total width = 4 * elem_w
# SoA ribbon: i₀,i₁,f₀,f₁  under the left side
# AoS ribbon: i₀,f₀,i₁,f₁  under the right side

# SoA ribbon (left)
soa_ribbon_y = -0.3
soa_x0 = i_x0 - cell / 2  # align with left edge of i-array
soa_total_w = 2 * cell      # same width as the 2-cell arrays
soa_elem_w = soa_total_w / n_elems

soa_colors = [
    elem_color_rgb(col_i_rgb, 0), elem_color_rgb(col_i_rgb, 1),
    elem_color_rgb(col_f_rgb, 0), elem_color_rgb(col_f_rgb, 1),
]
soa_labels = ['i$_0$', 'i$_1$', 'f$_0$', 'f$_1$']

ax.text(soa_x0 + soa_total_w / 2, soa_ribbon_y + ribbon_h + 0.25,
        'SoA layout', ha='center', va='center',
        fontsize=7, fontstyle='italic', color='#555')

for addr in range(n_elems):
    x = soa_x0 + addr * soa_elem_w
    rect = mpatches.FancyBboxPatch(
        (x + 0.01, soa_ribbon_y + 0.03), soa_elem_w - 0.02, ribbon_h - 0.06,
        boxstyle="round,pad=0.01,rounding_size=0.03",
        facecolor=soa_colors[addr], edgecolor='#bbb', linewidth=0.3
    )
    ax.add_patch(rect)
    ax.text(x + soa_elem_w / 2, soa_ribbon_y + ribbon_h / 2,
            soa_labels[addr], ha='center', va='center',
            fontsize=7, color='#222')

for b in range(0, n_elems + 1, B):
    x = soa_x0 + b * soa_elem_w
    ax.plot([x, x], [soa_ribbon_y - 0.02, soa_ribbon_y + ribbon_h + 0.02],
            ls='--', color='#aaaaaa', lw=1.0)

for b in range(n_elems // B):
    x = soa_x0 + (b * B + 1) * soa_elem_w
    ax.text(x, soa_ribbon_y - 0.18, f'B$_{b}$',
            ha='center', va='center', fontsize=6, color='#444')

# AoS ribbon (right) — half-width per element, same total width as SoA
aos_ribbon_y = -0.3
aos_x0 = (pack_x0 - cell / 2) + (4 * cell - soa_total_w) / 2  # center under packed array
aos_total_w = soa_total_w     # same total width as SoA ribbon
aos_elem_w = aos_total_w / n_elems  # half of SoA elem_w since same width, same count

aos_colors = [
    elem_color_rgb(col_i_rgb, 0), elem_color_rgb(col_f_rgb, 0),
    elem_color_rgb(col_i_rgb, 1), elem_color_rgb(col_f_rgb, 1),
]
aos_labels = ['i$_0$', 'f$_0$', 'i$_1$', 'f$_1$']

ax.text(aos_x0 + aos_total_w / 2, aos_ribbon_y + ribbon_h + 0.25,
        'AoS layout', ha='center', va='center',
        fontsize=7, fontstyle='italic', color='#555')

for addr in range(n_elems):
    x = aos_x0 + addr * aos_elem_w
    rect = mpatches.FancyBboxPatch(
        (x + 0.01, aos_ribbon_y + 0.03), aos_elem_w - 0.02, ribbon_h - 0.06,
        boxstyle="round,pad=0.01,rounding_size=0.03",
        facecolor=aos_colors[addr], edgecolor='#bbb', linewidth=0.3
    )
    ax.add_patch(rect)
    ax.text(x + aos_elem_w / 2, aos_ribbon_y + ribbon_h / 2,
            aos_labels[addr], ha='center', va='center',
            fontsize=6, color='#222')

# Block separators (every B elements) + middle dashed line
for b in range(0, n_elems + 1, B):
    x = aos_x0 + b * aos_elem_w
    ax.plot([x, x], [aos_ribbon_y - 0.02, aos_ribbon_y + ribbon_h + 0.02],
            ls='--', color='#aaaaaa', lw=1.0)

# Block labels
for b in range(n_elems // B):
    x = aos_x0 + (b * B + B / 2) * aos_elem_w
    ax.text(x, aos_ribbon_y - 0.18, f'B$_{b}$',
            ha='center', va='center', fontsize=6, color='#444')

# ── Limits ────────────────────────────────────────────────────────────
"""
ax.set_xlim(-0.1, pack_x0 + 3.5*cell + 0.3)
# soa_ribbon_y - 0.55
ax.set_ylim(-0.1, i_y + cell/2 + 0.25)
ax.set_aspect('equal')
ax.axis('off')

plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
plt.savefig('pack_unpack.pdf', bbox_inches='tight', dpi=300)
plt.savefig('pack_unpack.png', bbox_inches='tight', dpi=300)
print("Done: pack_unpack")