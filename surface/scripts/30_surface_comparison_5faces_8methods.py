"""5-face × 8-method comparison figure: DFT + 4 PFP modes + 3 MEAM (incl. Ko 2018).
Faces: 100, 101, 111, 001, 110 (user-specified subset).
Based on 23_surface_comparison_full.py with Ko 2018 added and FIG path corrected.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

FIG = os.path.expanduser('~/project_dir/surface_gb/figures')

FACES = ['(100)', '(101)', '(111)', '(001)', '(110)']
DATA = {
    'DFT/PBE (BFGS)':    [0.4922, 0.5056, 0.5500, 0.5505, 0.5269],
    'PFP/PBE':           [0.3855, 0.4127, 0.4171, 0.4601, 0.4870],
    'PFP/PBE+D3':        [0.6422, 0.8325, 0.6851, 0.7344, 0.7375],
    'PFP/R²SCAN':        [0.7846, 0.4766, 0.5388, 0.5275, 0.5993],
    'PFP/R²SCAN+D3':     [0.6148, 0.7073, 0.6901, 0.6804, 0.7306],
    'MEAM/Ravelo 1997':  [0.7350, 0.8183, 0.9093, 0.8968, 0.9542],
    'MEAM/Etesami 2018': [1.3763, 1.3547, 1.3847, 1.2997, 1.5874],
    'MEAM/Ko 2018':      [0.3454, 0.3672, 0.4728, 0.4853, 0.4943],
}
EXP_AVG = 0.670

colors = {
    'DFT/PBE (BFGS)':     '#000000',
    'PFP/PBE':            '#2066a8',
    'PFP/PBE+D3':         '#6fb7c8',
    'PFP/R²SCAN':         '#d97c2b',
    'PFP/R²SCAN+D3':      '#f2b16d',
    'MEAM/Ravelo 1997':   '#6a3d9a',
    'MEAM/Etesami 2018':  '#b15928',
    'MEAM/Ko 2018':       '#17becf',  # cyan; matches MEAM[3] in fig3/fig4
}

dft = np.array(DATA['DFT/PBE (BFGS)'])
mae = {k: float(np.mean(np.abs(np.array(v) - dft)))
       for k, v in DATA.items() if k != 'DFT/PBE (BFGS)'}

x = np.arange(len(FACES))
methods = list(DATA.keys())
n = len(methods)
w = 0.85 / n

fig, ax = plt.subplots(figsize=(11.5, 6.4))
for i, m in enumerate(methods):
    offset = (i - (n - 1) / 2) * w
    ax.bar(x + offset, DATA[m], w, label=m, color=colors[m],
           edgecolor='black' if m == 'DFT/PBE (BFGS)' else None,
           linewidth=1.4 if m == 'DFT/PBE (BFGS)' else 0.2)
    for xi, v in zip(x + offset, DATA[m]):
        ax.text(xi, v + 0.02, f'{v:.2f}', ha='center', fontsize=5.8, rotation=90)

ax.axhline(EXP_AVG, linestyle='--', color='red', lw=1.2, alpha=0.7,
           label=f'exp. avg ~{EXP_AVG:.2f} J/m² [Tyson–Miller 1977]')
ax.set_xticks(x); ax.set_xticklabels(FACES)
ax.set_xlabel('surface')
ax.set_ylabel('surface energy γ (J/m²)')
ax.set_title('β-Sn surface energies — 5 faces × 8 methods (with Ko 2018 MEAM added)')
ax.set_ylim(0, 1.75)
ax.legend(loc='upper left', fontsize=8.5, ncol=2, framealpha=0.95)
ax.grid(axis='y', alpha=0.3)

# face ordering per method
orderings = {}
for m, v in DATA.items():
    idxs = np.argsort(v)
    orderings[m] = '<'.join(FACES[i].replace('(','').replace(')','') for i in idxs)

txt_lines = ['face ordering (low → high γ):']
for m, o in orderings.items():
    txt_lines.append(f'  {m:22s}: {o}')
txt_lines.append('')
txt_lines.append('MAE vs DFT/PBE (5 faces, J/m²):')
for m, v in sorted(mae.items(), key=lambda x: x[1]):
    txt_lines.append(f'  {m:22s}: {v:>6.3f}')
txt = '\n'.join(txt_lines)
ax.text(0.995, 0.99, txt, transform=ax.transAxes,
        ha='right', va='top', fontsize=7,
        family='monospace',
        bbox=dict(facecolor='white', edgecolor='gray', alpha=0.95, boxstyle='round,pad=0.5'))

plt.tight_layout()
outp_png = os.path.join(FIG, 'surface_energy_comparison_5faces_8methods.png')
outp_tif = os.path.join(FIG, 'surface_energy_comparison_5faces_8methods.tif')
plt.savefig(outp_png, dpi=150, bbox_inches='tight')
plt.savefig(outp_tif, dpi=300, bbox_inches='tight')
print(f'saved -> {outp_png}')
print(f'saved -> {outp_tif}')

print('\nMAE vs DFT/PBE (5 faces, J/m²):')
for m, v in sorted(mae.items(), key=lambda x: x[1]):
    print(f'  {m:22s}: {v:>6.3f}')

print('\nFace ordering:')
for m, o in orderings.items():
    print(f'  {m:22s}: {o}')
