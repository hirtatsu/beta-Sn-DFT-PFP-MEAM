"""Regenerate surface energy comparison figure with 4 PFP modes + DFT BFGS."""
import os
import numpy as np
import matplotlib.pyplot as plt

FIG = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/figures')

# Face order: lowest-γ first per DFT order
FACES = ['(100)', '(101)', '(110)', '(001)']
DATA = {
    'DFT/PBE (BFGS)':    [492.2, 505.6, 526.9, 550.5],
    'PFP/PBE':           [385.5, 412.7, 487.0, 460.1],
    'PFP/PBE+D3':        [642.2, 832.5, 737.5, 734.4],
    'PFP/R²SCAN':        [784.6, 476.6, 599.3, 527.5],
    'PFP/R²SCAN+D3':     [614.8, 707.3, 730.6, 680.4],
}
EXP_AVG = 670  # Tyson-Miller 1977, mJ/m²

colors = {
    'DFT/PBE (BFGS)': '#000000',
    'PFP/PBE':        '#2066a8',
    'PFP/PBE+D3':     '#6fb7c8',
    'PFP/R²SCAN':     '#d97c2b',
    'PFP/R²SCAN+D3':  '#f2b16d',
}

x = np.arange(len(FACES))
methods = list(DATA.keys())
n_methods = len(methods)
w = 0.85 / n_methods

fig, ax = plt.subplots(figsize=(9.6, 5.4))
for i, method in enumerate(methods):
    offset = (i - (n_methods - 1) / 2) * w
    bar = ax.bar(x + offset, DATA[method], w, label=method,
                 color=colors[method],
                 edgecolor='black' if method == 'DFT/PBE (BFGS)' else None,
                 linewidth=1.2 if method == 'DFT/PBE (BFGS)' else 0.3)
    # annotate
    for xi, v in zip(x + offset, DATA[method]):
        ax.text(xi, v + 12, f'{v:.0f}', ha='center', fontsize=6.5, rotation=0)

ax.axhline(EXP_AVG, linestyle='--', color='red', lw=1.2, alpha=0.7,
           label=f'exp. avg ~{EXP_AVG} mJ/m² [Tyson–Miller 1977]')
ax.set_xticks(x)
ax.set_xticklabels(FACES)
ax.set_xlabel('surface (ordered by DFT/PBE γ: low → high)')
ax.set_ylabel('surface energy γ (mJ/m²)')
ax.set_title('β-Sn surface energies — 4 PFP modes vs. DFT/PBE (BFGS)\n'
             'conditions matched: mode-self-consistent bulk lattice, '
             'all atoms free, fmax~1.5×10⁻² eV/Å')
ax.set_ylim(0, 900)
ax.legend(loc='upper left', fontsize=8.5, ncol=2, framealpha=0.95)
ax.grid(axis='y', alpha=0.3)

# inset: ordering summary text
order_text = (
    "face ordering (low → high γ):\n"
    "• DFT/PBE    : (100) < (101) < (110) < (001)\n"
    "• PFP/PBE    : (100) < (101) < (110) < (001)  ✓ match\n"
    "• PFP/PBE+D3 : (100) < (001) ≈ (110) < (101)  ✗\n"
    "• PFP/R²SCAN : (101) < (001) < (110) < (100)  ✗ reversed\n"
    "• PFP/R²SCAN+D3 : (100) < (001) < (101) < (110)  ✗"
)
ax.text(0.98, 0.97, order_text, transform=ax.transAxes,
        ha='right', va='top', fontsize=7.5,
        bbox=dict(facecolor='white', edgecolor='gray', alpha=0.9, boxstyle='round,pad=0.4'))

plt.tight_layout()
outp_png = os.path.join(FIG, 'surface_energy_comparison_4modes.png')
outp_tif = os.path.join(FIG, 'surface_energy_comparison_4modes.tif')
plt.savefig(outp_png, dpi=150, bbox_inches='tight')
plt.savefig(outp_tif, dpi=300, bbox_inches='tight')
print(f'saved -> {outp_png}')
print(f'saved -> {outp_tif}')
