"""Final 7-method comparison figure: DFT + 4 PFP modes + 2 MEAM."""
import os
import numpy as np
import matplotlib.pyplot as plt

FIG = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/figures')

FACES = ['(100)', '(101)', '(110)', '(001)']
DATA = {
    'DFT/PBE (BFGS)':    [492.2, 505.6, 526.9, 550.5],
    'PFP/PBE':           [385.5, 412.7, 487.0, 460.1],
    'PFP/PBE+D3':        [642.2, 832.5, 737.5, 734.4],
    'PFP/R²SCAN':        [784.6, 476.6, 599.3, 527.5],
    'PFP/R²SCAN+D3':     [614.8, 707.3, 730.6, 680.4],
    'MEAM/Ravelo 1997':  [735.0, 818.3, 954.2, 896.8],
    'MEAM/Etesami 2018': [1376.3, 1354.7, 1587.4, 1299.7],
}
EXP_AVG = 670

colors = {
    'DFT/PBE (BFGS)':     '#000000',
    'PFP/PBE':            '#2066a8',
    'PFP/PBE+D3':         '#6fb7c8',
    'PFP/R²SCAN':         '#d97c2b',
    'PFP/R²SCAN+D3':      '#f2b16d',
    'MEAM/Ravelo 1997':   '#6a3d9a',
    'MEAM/Etesami 2018':  '#b15928',
}

# MAE vs DFT
dft = np.array(DATA['DFT/PBE (BFGS)'])
mae = {}
for k, v in DATA.items():
    if k == 'DFT/PBE (BFGS)':
        continue
    mae[k] = float(np.mean(np.abs(np.array(v) - dft)))

# Bar chart
x = np.arange(len(FACES))
methods = list(DATA.keys())
n = len(methods)
w = 0.85 / n

fig, ax = plt.subplots(figsize=(11.5, 6.2))
for i, m in enumerate(methods):
    offset = (i - (n - 1) / 2) * w
    ax.bar(x + offset, DATA[m], w, label=m, color=colors[m],
           edgecolor='black' if m == 'DFT/PBE (BFGS)' else None,
           linewidth=1.4 if m == 'DFT/PBE (BFGS)' else 0.2)
    for xi, v in zip(x + offset, DATA[m]):
        ax.text(xi, v + 20, f'{v:.0f}', ha='center', fontsize=5.8, rotation=90)

ax.axhline(EXP_AVG, linestyle='--', color='red', lw=1.2, alpha=0.7,
           label=f'exp. avg ~{EXP_AVG} mJ/m² [Tyson–Miller 1977]')
ax.set_xticks(x); ax.set_xticklabels(FACES)
ax.set_xlabel('surface (ordered by DFT/PBE γ: low → high)')
ax.set_ylabel('surface energy γ (mJ/m²)')
ax.set_title('β-Sn surface energies — DFT, 4 PFP modes, 2 MEAM\n'
             '(mode-self-consistent bulk lattice, full atom relax)')
ax.set_ylim(0, 1750)
ax.legend(loc='upper left', fontsize=8.5, ncol=2, framealpha=0.95)
ax.grid(axis='y', alpha=0.3)

# Ordering + MAE box
txt = "face ordering (low → high γ):\n"
txt += "• DFT/PBE              : (100)<(101)<(110)<(001)\n"
txt += "• PFP/PBE              : (100)<(101)<(110)<(001)  ✓ match\n"
txt += "• PFP/PBE+D3           : (100)<(001)≈(110)<(101)  ✗\n"
txt += "• PFP/R²SCAN           : (101)<(001)<(110)<(100)  ✗ reversed\n"
txt += "• PFP/R²SCAN+D3        : (100)<(001)<(101)<(110)  ✗\n"
txt += "• MEAM/Ravelo 1997     : (100)<(101)<(001)<(110)  partial\n"
txt += "• MEAM/Etesami 2018    : (001)<(101)<(100)<(110)  ✗ reversed\n"
txt += "\nMAE vs DFT/PBE (mJ/m²):\n"
for k, v in mae.items():
    txt += f"  {k:20s} : {v:>6.0f}\n"
ax.text(0.985, 0.97, txt, transform=ax.transAxes,
        ha='right', va='top', fontsize=7,
        family='monospace',
        bbox=dict(facecolor='white', edgecolor='gray', alpha=0.95, boxstyle='round,pad=0.5'))

plt.tight_layout()
outp_png = os.path.join(FIG, 'surface_energy_comparison_full.png')
outp_tif = os.path.join(FIG, 'surface_energy_comparison_full.tif')
plt.savefig(outp_png, dpi=150, bbox_inches='tight')
plt.savefig(outp_tif, dpi=300, bbox_inches='tight')
print(f'saved -> {outp_png}')

# print summary
print('\nMAE vs DFT/PBE (mJ/m²):')
for k, v in sorted(mae.items(), key=lambda x: x[1]):
    print(f"  {k:20s} : {v:>6.1f}")
