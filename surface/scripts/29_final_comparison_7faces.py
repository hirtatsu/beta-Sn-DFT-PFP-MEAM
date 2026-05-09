"""Final bar chart + MAE for 7 faces × 7 methods."""
import os
import json
import numpy as np
import matplotlib.pyplot as plt

FIG = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/figures')

FACES = [(1,0,0), (1,0,1), (1,1,0), (0,0,1), (1,1,1), (2,1,1), (1,1,2)]
FACE_LABELS = ['(100)', '(101)', '(110)', '(001)', '(111)', '(211)', '(112)']

DATA = {
    'DFT/PBE (BFGS)':    [492.2, 505.6, 526.9, 550.5, 550.0, 538.4, 358.9],
    'PFP/PBE':           [385.5, 412.7, 487.0, 460.1, 417.1, 424.0, 442.5],
    'PFP/PBE+D3':        [642.2, 832.5, 737.5, 734.4, 685.1, 704.1, 668.2],
    'PFP/R²SCAN':        [784.6, 476.6, 599.3, 527.5, 538.8, 529.2, 483.7],
    'PFP/R²SCAN+D3':     [614.8, 707.3, 730.6, 680.4, 690.1, 680.9, 633.9],
    'MEAM/Ravelo 1997':  [735.0, 818.3, 954.2, 896.8, 909.3, 941.4, 814.9],
    'MEAM/Etesami 2018': [1376.3, 1354.7, 1587.4, 1299.7, 1384.7, 1224.8, 1303.3],
}
EXP_AVG = 670

dft = np.array(DATA['DFT/PBE (BFGS)'])
mae = {m: float(np.mean(np.abs(np.array(v) - dft))) for m, v in DATA.items() if m != 'DFT/PBE (BFGS)'}

# ordering per method
orderings = {}
for m, v in DATA.items():
    idxs = np.argsort(v)
    orderings[m] = '<'.join(FACE_LABELS[i].replace('(','').replace(')','') for i in idxs)

colors = {
    'DFT/PBE (BFGS)':     '#000000',
    'PFP/PBE':            '#2066a8',
    'PFP/PBE+D3':         '#6fb7c8',
    'PFP/R²SCAN':         '#d97c2b',
    'PFP/R²SCAN+D3':      '#f2b16d',
    'MEAM/Ravelo 1997':   '#6a3d9a',
    'MEAM/Etesami 2018':  '#b15928',
}

x = np.arange(len(FACE_LABELS))
methods = list(DATA.keys())
n = len(methods)
w = 0.85 / n

fig, ax = plt.subplots(figsize=(14.2, 7.0))
for i, m in enumerate(methods):
    offset = (i - (n - 1) / 2) * w
    ax.bar(x + offset, DATA[m], w, label=m, color=colors[m],
           edgecolor='black' if m == 'DFT/PBE (BFGS)' else None,
           linewidth=1.3 if m == 'DFT/PBE (BFGS)' else 0.2)
    for xi, v in zip(x + offset, DATA[m]):
        ax.text(xi, v + 20, f'{v:.0f}', ha='center', fontsize=5.5, rotation=90)

ax.axhline(EXP_AVG, linestyle='--', color='red', lw=1.2, alpha=0.7,
           label=f'exp. avg ~{EXP_AVG} mJ/m² [Tyson–Miller]')
ax.set_xticks(x); ax.set_xticklabels(FACE_LABELS)
ax.set_xlabel('surface')
ax.set_ylabel('surface energy γ (mJ/m²)')
ax.set_title('β-Sn surface energies — 7 faces × 7 methods')
ax.set_ylim(0, 1750)
ax.legend(loc='upper left', fontsize=8.5, ncol=2, framealpha=0.95)
ax.grid(axis='y', alpha=0.3)

txt_lines = ['face ordering (low → high γ):']
for m, o in orderings.items():
    txt_lines.append(f'  {m:22s}: {o}')
txt_lines.append('')
txt_lines.append('MAE vs DFT/PBE (7 faces, mJ/m²):')
for m, v in sorted(mae.items(), key=lambda x: x[1]):
    txt_lines.append(f'  {m:22s}: {v:>6.0f}')
txt = '\n'.join(txt_lines)
ax.text(0.995, 0.99, txt, transform=ax.transAxes,
        ha='right', va='top', fontsize=6.8,
        family='monospace',
        bbox=dict(facecolor='white', edgecolor='gray', alpha=0.95, boxstyle='round,pad=0.5'))

plt.tight_layout()
outp_png = os.path.join(FIG, 'surface_energy_comparison_7faces.png')
outp_tif = os.path.join(FIG, 'surface_energy_comparison_7faces.tif')
plt.savefig(outp_png, dpi=150, bbox_inches='tight')
plt.savefig(outp_tif, dpi=300, bbox_inches='tight')
print(f'saved -> {outp_png}')

print('\nMAE vs DFT/PBE (7 faces):')
for m, v in sorted(mae.items(), key=lambda x: x[1]):
    print(f'  {m:22s}: {v:>6.1f}')

print('\nFace ordering:')
for m, o in orderings.items():
    print(f'  {m:22s}: {o}')
