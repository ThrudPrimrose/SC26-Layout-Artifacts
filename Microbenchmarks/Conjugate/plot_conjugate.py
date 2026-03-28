#!/usr/bin/env python3
"""2x2 grids for runtime and bandwidth. Usage: python plot_conj.py [--vl 16] [--N 27]"""
import argparse, os
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch

p = argparse.ArgumentParser()
p.add_argument('--vl-gpu', type=int, default=16)
p.add_argument('--vl-cpu', type=int, default=8)
p.add_argument('--N', type=int, default=27, help='log2 of number of complex elements')
args = p.parse_args()
VL_GPU, VL_CPU, N = args.vl_gpu, args.vl_cpu, 1 << args.N

# in-place:  read im + write im = 2*N doubles
# out-of-place: read re+im + write re+im = 4*N doubles
BYTES_IP  = 2 * N * 8
BYTES_OOP = 4 * N * 8

plt.rcParams.update({'font.size':10, 'figure.facecolor':'white'})
C_IP, C_OOP = '#2563eb', '#f59e0b'
LAYS_CPU = ['AoS', 'SoA', f'AoSoA-{VL_CPU}']
LAYS_GPU = ['AoS', 'SoA', f'AoSoA-{VL_GPU}']
WARMUP = 5

LEGEND = [Patch(facecolor=C_IP, alpha=0.7, label='In-Place'),
          Patch(facecolor=C_OOP, alpha=0.7, label='Out-of-Place')]

def read(path):
    if not os.path.exists(path): return None
    df = pd.read_csv(path); df.columns = [c.strip().lower() for c in df.columns]
    return df[df['run'] >= WARMUP]

def get_ms(df, lay):
    s = df[df['layout']==lay]['ms']
    return s.values if not s.empty else np.array([])

def violin(ax, pos, vals, color):
    if len(vals) == 0: return
    vp = ax.violinplot([vals], positions=[pos], widths=0.35, showmedians=True, showextrema=False)
    for body in vp['bodies']:
        body.set_facecolor(color); body.set_alpha(0.7); body.set_edgecolor(color)
    vp['cmedians'].set_color('black')

def grid_style(ax):
    ax.yaxis.grid(True, alpha=0.3); ax.xaxis.grid(False)
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.grid(True, which='minor', alpha=0.15, linewidth=0.5)


P = {
    (0,0): ('results/daint/results_cpu.csv',       'results/daint/results_cpu_oop.csv',
            'Grace CPU (Neoverse V2)'),
    (0,1): ('results/daint/results_gpu_inplace.csv','results/daint/results_gpu_oop.csv',
            'Hopper GPU'),
    (1,0): ('results/beverin/results_cpu.csv',      'results/beverin/results_cpu_oop.csv',
            'AMD EPYC (MI300A)'),
    (1,1): ('results/beverin/results_gpu_inplace.csv','results/beverin/results_gpu_oop.csv',
            'AMD MI300A GPU'),
}

os.makedirs('results', exist_ok=True)

def panel(ax, ip_path, oop_path, title, metric, gpu=False):
    lays = LAYS_GPU if gpu else LAYS_CPU
    ip, oop = read(ip_path), read(oop_path)
    if ip is None and oop is None:
        ax.text(.5,.5,'no data',ha='center',va='center',transform=ax.transAxes,color='grey')
        ax.set_title(title); return
    x, w = np.arange(len(lays)), 0.2
    for off, df, c, nbytes in [(-w, ip, C_IP, BYTES_IP), (w, oop, C_OOP, BYTES_OOP)]:
        if df is None: continue
        for i, lay in enumerate(lays):
            ms = get_ms(df, lay)
            if len(ms) == 0: continue
            vals = (nbytes / (ms * 1e9)) if metric == 'gbps' else ms
            violin(ax, i+off, vals, c)
    ax.set_xticks(x); ax.set_xticklabels(lays)
    ax.set_ylabel('TB/s' if metric == 'gbps' else 'ms')
    ax.set_title(title); ax.set_ylim(bottom=0); grid_style(ax)

for metric, suffix, title in [('gbps', 'bandwidth', 'Complex Conjugate Bandwidth'),
                               ('ms',   'runtime',   'Complex Conjugate Runtime')]:
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5), sharey="col")
    fig.suptitle(title, fontsize=13, y=.9)
    for (r,c),(ip,oop,t) in P.items(): panel(axes[r][c], ip, oop, t, metric, c==1)
    fig.legend(handles=LEGEND, loc='lower center', ncol=2, fontsize=11,
               frameon=False, bbox_to_anchor=(0.5, 0.0))
    plt.tight_layout(rect=[0, 0.045, 1, 0.93])
    for ax_row in axes:
        for ax in ax_row:
            ax.autoscale(axis='y')
            ax.set_ylim(bottom=0)
            ax.set_ylim(top=ax.get_ylim()[1] * 1.2)
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    for ext in ['png', 'pdf']:
        path = f'results/conj_{suffix}.{ext}'
        fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote results/conj_{suffix}.png / .pdf')

print(f'\nN = {N} ({args.N} bits), sizeof(complex) = 16 B')
print(f'  in-place bytes:      {BYTES_IP/1e9:.3f} GB  (read+write im only)')
print(f'  out-of-place bytes:  {BYTES_OOP/1e9:.3f} GB  (read+write re+im)')