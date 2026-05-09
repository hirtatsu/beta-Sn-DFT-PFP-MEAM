"""Cij comparison figure (fig3) and MAPE ranking (fig4) including Ko 2018 MEAM
(NIST IPR id 2018--Ko-W-S-Kim-D-H-Kwon-Y-J-Lee-M--Sn) as 'MEAM [3]'.

Based on make_figures_v2.py with these changes:
  - 'MEAM [3]' bar added to fig3 with values from the LAMMPS lmp_kokkos run
  - 'MEAM [3]' entry added to fig4 with MAPE = 28.24% (rank 3 of 8)
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
OUTDIR = os.path.join(PROJECT_ROOT, 'figures')
os.makedirs(OUTDIR, exist_ok=True)

methods = ['DFT/PBE', 'PFP/PBE', 'PFP/PBE+D3',
           'PFP/r$^2$SCAN', 'PFP/r$^2$SCAN+D3',
           'MEAM [1]', 'MEAM [2]', 'MEAM [3]']
cij_labels = [r'$C_{11}$', r'$C_{12}$', r'$C_{13}$',
              r'$C_{33}$', r'$C_{44}$', r'$C_{66}$']
cij = {
    'DFT/PBE':          [ 89.7, 17.4, 31.5,  91.8, 17.9, 17.6],
    'PFP/PBE':          [114.0, 41.5, 41.2, 104.6, 29.5, 31.9],
    'PFP/PBE+D3':       [ 98.5, 36.2, 36.5, 121.1, 22.9, 16.2],
    'PFP/r$^2$SCAN':    [118.0, 15.2, 40.7, 129.7, 29.7, 24.8],
    'PFP/r$^2$SCAN+D3': [123.9, 20.7, 43.0, 132.1, 34.3, 27.2],
    'MEAM [1]':         [108.2, 63.8, 24.4, 139.6,  1.5, 22.4],
    'MEAM [2]':         [131.6, 48.7, 16.6, 157.1,  2.4, 28.0],
    'MEAM [3]':         [ 89.75, 46.67, 36.92,  93.73,  7.93, 10.63],
}
exp = [73.4, 59.1, 35.8, 90.7, 22.0, 24.0]

palette = plt.get_cmap('tab10').colors
method_colors = {
    'DFT/PBE':          palette[0],
    'PFP/PBE':          palette[1],
    'PFP/PBE+D3':       palette[2],
    'PFP/r$^2$SCAN':    palette[3],
    'PFP/r$^2$SCAN+D3': palette[4],
    'MEAM [1]':         palette[7],
    'MEAM [2]':         palette[8],
    'MEAM [3]':         palette[9],
}

plt.rcParams.update({
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 11,
    'legend.fontsize': 9, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'sans',
    'mathtext.it': 'sans:italic',
    'mathtext.bf': 'sans:bold',
    'axes.linewidth': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 3.0, 'ytick.major.size': 3.0,
})


def savetiff(fig, name, tight=True):
    path = os.path.join(OUTDIR, f'{name}.tif')
    kwargs = {'bbox_inches': 'tight'} if tight else {}
    fig.savefig(path, format='tiff', dpi=600,
                pil_kwargs={'compression': 'tiff_lzw'},
                **kwargs)
    fig.savefig(os.path.join(OUTDIR, f'{name}_preview.png'), dpi=150,
                **kwargs)
    print(f'wrote {path}')


fig1, ax = plt.subplots(figsize=(8.0, 4.6))
n_methods = len(methods)
n_comps = len(cij_labels)
width = 0.10
x = np.arange(n_comps)
for i, m in enumerate(methods):
    offset = (i - (n_methods - 1) / 2) * width
    ax.bar(x + offset, cij[m], width, label=m,
           color=method_colors[m], edgecolor='black', linewidth=0.4)

half = n_methods / 2 * width
for i, e in enumerate(exp):
    ax.hlines(e, x[i] - half - width / 2, x[i] + half + width / 2,
              colors='red', linewidth=1.8, linestyles=(0, (4, 1.5)), zorder=5)
exp_handle = plt.Line2D([], [], color='red', linewidth=1.8,
                        linestyle=(0, (4, 1.5)), label='Experiment [4]')

ax.set_xticks(x)
ax.set_xticklabels(cij_labels)
ax.set_ylabel('Elastic constant (GPa)')
ax.set_ylim(0, 170)
ax.grid(axis='y', linestyle=':', linewidth=0.5, alpha=0.6)
ax.set_axisbelow(True)

bar_handles = [mpatches.Patch(facecolor=method_colors[m], edgecolor='black',
                              linewidth=0.4, label=m) for m in methods]
handles = bar_handles + [exp_handle]

ncol_legend = 5
nrow_legend = (len(handles) + ncol_legend - 1) // ncol_legend
reordered = []
for col in range(ncol_legend):
    for row in range(nrow_legend):
        idx = row * ncol_legend + col
        if idx < len(handles):
            reordered.append(handles[idx])

ax.legend(handles=reordered, loc='lower center',
          bbox_to_anchor=(0.5, 1.02), ncol=ncol_legend,
          frameon=False, columnspacing=1.2, handlelength=1.6,
          handletextpad=0.5)
fig1.tight_layout()
savetiff(fig1, 'fig3_cij_comparison')
plt.close(fig1)


# ============================================================
# Fig 4: MAPE ranking with MEAM [3] (Ko 2018) added
# ============================================================
# Compute MEAM[3] MAPE explicitly so it stays in sync with cij values above
def mape(theory, exp_):
    return float(np.mean([abs(t - e) / e * 100 for t, e in zip(theory, exp_)]))

mape_meam3 = mape(cij['MEAM [3]'], exp)
print(f'MEAM [3] (Ko 2018) MAPE = {mape_meam3:.2f}%')

mape_ranking = sorted([
    ('PFP/PBE+D3', 24.2), ('DFT/PBE', 25.2), ('PFP/PBE', 30.5),
    ('PFP/r$^2$SCAN', 38.3), ('MEAM [1]', 40.2),
    ('PFP/r$^2$SCAN+D3', 44.8), ('MEAM [2]', 54.9),
    ('MEAM [3]', round(mape_meam3, 1)),
], key=lambda x: x[1])

fig2, ax = plt.subplots(figsize=(6.4, 4.2))
names = [m for m, _ in mape_ranking]
values = [v for _, v in mape_ranking]
ypos = np.arange(len(names))[::-1]
colors = [method_colors[m] for m in names]
bars = ax.barh(ypos, values, color=colors, edgecolor='black', linewidth=0.5)
for b, v in zip(bars, values):
    ax.text(v + 0.6, b.get_y() + b.get_height() / 2, f'{v:.1f}%',
            va='center', fontsize=9)
ax.set_yticks(ypos)
ax.set_yticklabels(names)
ax.set_xlabel('MAPE vs experiment (%)  — average over six $C_{ij}$')
ax.set_xlim(0, max(values) * 1.18)
ax.grid(axis='x', linestyle=':', linewidth=0.5, alpha=0.6)
ax.set_axisbelow(True)
# highlight the best (lowest MAPE) — bar at top of plot
bars[0].set_edgecolor('black')
bars[0].set_linewidth(1.5)
fig2.tight_layout()
savetiff(fig2, 'fig4_mape_ranking')
plt.close(fig2)
