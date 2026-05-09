"""Generate paper figures for β-Sn surface energies — 5 low-index faces × 8 methods.

Reads ./data/surface_energies_5faces.json (8-method canonical table) and writes
into the project-root figures/ directory:
  fig5_surface_energies_8methods.{png,tif}           — manuscript Fig. 5 (bar chart)
  supp_wulff_8methods.{png,tif}                      — supplementary (Wulff 4×2 overview)
  supp_wulff_<METHOD>.{png,tif}                      — one per method (8 files)
  supp_parity_vs_DFT_8methods.{png,tif}              — supplementary (parity plot)

Dependencies: numpy, matplotlib, ase, wulffpack.
Run with: python surface/make_figures.py
"""
import os, json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3D projection)
from matplotlib.patches import Patch
from ase import Atoms
from wulffpack import SingleCrystal

HERE = os.path.dirname(os.path.abspath(__file__))                  # surface/
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..'))            # beta-sn-method-comparison/
FIG = os.path.join(PROJECT_ROOT, 'figures')
DATA = os.path.join(HERE, 'data')
os.makedirs(FIG, exist_ok=True)

with open(os.path.join(DATA, 'surface_energies_5faces.json')) as f:
    payload = json.load(f)

FACE_KEYS = payload['faces_order']                      # e.g. ['100','101','110','111','001']
FACE_HKL = [tuple(h) for h in payload['faces_miller']]
FACES = [f'({k})' for k in FACE_KEYS]
units = payload.get('units', 'J/m²')
unit_scale = 1.0 if units.startswith('mJ') else 1000.0  # convert to mJ/m² internally
EXP_AVG = 670.0  # mJ/m² (Tyson–Miller 1977 polycrystalline average)

METHOD_ORDER = [
    'DFT/PBE (BFGS)',
    'PFP/PBE', 'PFP/PBE+D3', 'PFP/R²SCAN', 'PFP/R²SCAN+D3',
    'MEAM/Ravelo 1997', 'MEAM/Etesami 2018', 'MEAM/Ko 2018',
]
gamma = {m: np.array([payload['gamma_table'][m][k] * unit_scale for k in FACE_KEYS])
         for m in METHOD_ORDER}
dft = gamma['DFT/PBE (BFGS)']
methods = list(METHOD_ORDER)
n = len(methods)

mae = {m: float(np.mean(np.abs(gamma[m] - dft))) for m in methods if m != 'DFT/PBE (BFGS)'}
orderings = {m: '<'.join(FACE_KEYS[i] for i in np.argsort(gamma[m])) for m in methods}

colors = {
    'DFT/PBE (BFGS)':    '#000000',
    'PFP/PBE':           '#2066a8',
    'PFP/PBE+D3':        '#6fb7c8',
    'PFP/R²SCAN':        '#d97c2b',
    'PFP/R²SCAN+D3':     '#f2b16d',
    'MEAM/Ravelo 1997':  '#6a3d9a',
    'MEAM/Etesami 2018': '#b15928',
    'MEAM/Ko 2018':      '#33a02c',
}
markers = {
    'PFP/PBE':           'o',
    'PFP/PBE+D3':        's',
    'PFP/R²SCAN':        '^',
    'PFP/R²SCAN+D3':     'D',
    'MEAM/Ravelo 1997':  'P',
    'MEAM/Etesami 2018': 'X',
    'MEAM/Ko 2018':      'v',
}

# ============== Fig 1: Bar chart (8 methods × 5 faces) ==============
x = np.arange(len(FACES))
w = 0.85 / n
fig, ax = plt.subplots(figsize=(12.5, 6.2))
for i, m in enumerate(methods):
    offset = (i - (n - 1) / 2) * w
    is_ref = (m == 'DFT/PBE (BFGS)')
    ax.bar(x + offset, gamma[m], w, label=m, color=colors[m],
           edgecolor='black' if is_ref else 'none',
           linewidth=1.3 if is_ref else 0.0)
    for xi, v in zip(x + offset, gamma[m]):
        ax.text(xi, v + 18, f'{v:.0f}', ha='center', fontsize=6, rotation=90)

ax.axhline(EXP_AVG, linestyle='--', color='red', lw=1.2, alpha=0.7,
           label=f'exp. avg ~{EXP_AVG:.0f} mJ/m² [Tyson–Miller 1977]')
ax.set_xticks(x); ax.set_xticklabels(FACES)
ax.set_xlabel('surface')
ax.set_ylabel('surface energy γ (mJ/m²)')
ax.set_title('β-Sn surface energies — 5 low-index faces × 8 methods')
ax.set_ylim(0, 1750)
ax.legend(loc='upper left', fontsize=8.5, ncol=2, framealpha=0.95)
ax.grid(axis='y', alpha=0.3)

txt = 'face ordering (low → high γ):\n'
for m in methods:
    txt += f'  {m:22s}: {orderings[m]}\n'
txt += '\nMAE vs DFT/PBE (mJ/m²):\n'
for m, v in sorted(mae.items(), key=lambda kv: kv[1]):
    txt += f'  {m:22s}: {v:>6.0f}\n'
ax.text(0.995, 0.99, txt, transform=ax.transAxes,
        ha='right', va='top', fontsize=6.8, family='monospace',
        bbox=dict(facecolor='white', edgecolor='gray', alpha=0.95, boxstyle='round,pad=0.5'))

plt.tight_layout()
for ext, dpi in [('png', 150), ('tif', 300)]:
    plt.savefig(os.path.join(FIG, f'fig5_surface_energies_8methods.{ext}'),
                dpi=dpi, bbox_inches='tight')
print('wrote fig5_surface_energies_8methods.{png,tif}')
plt.close(fig)

# ============== Fig 2: Wulff shapes (8 methods, 4x2 grid) ==============
A_DFT, C_DFT = 5.970, 3.218  # use a single common reference lattice for shape comparability
BETA_SN = Atoms('Sn4',
                scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5],
                                  [0, 0.5, 0.25], [0.5, 0, 0.75]],
                cell=[[A_DFT, 0, 0], [0, A_DFT, 0], [0, 0, C_DFT]], pbc=True)

WULFF_COLORS = {
    (1, 0, 0): '#4C72B0', (1, 0, 1): '#8172B3', (1, 1, 0): '#C44E52',
    (1, 1, 1): '#CCB974', (0, 0, 1): '#55A868',
}


def make_wulff(gam_arr):
    se = {hkl: g for hkl, g in zip(FACE_HKL, gam_arr)}
    return SingleCrystal(surface_energies=se, primitive_structure=BETA_SN, natoms=2000)


cols = 4
rows = (n + cols - 1) // cols
fig = plt.figure(figsize=(4.2 * cols, 4.2 * rows))
for i, m in enumerate(methods, 1):
    ax = fig.add_subplot(rows, cols, i, projection='3d')
    sc = make_wulff(gamma[m])
    sc.make_plot(ax, alpha=0.85, linewidth=0.4, colors=WULFF_COLORS)
    gam_txt = ', '.join(f'({k}):{g:.0f}' for k, g in zip(FACE_KEYS, gamma[m]))
    ax.set_title(f'{m}\n{gam_txt}', fontsize=7)
    ax.view_init(elev=22, azim=35)
    ax.set_box_aspect([1, 1, 1])

handles = [Patch(color=c, label=f'({k[0]}{k[1]}{k[2]})') for k, c in WULFF_COLORS.items()]
fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=9,
           bbox_to_anchor=(0.5, -0.01))
plt.suptitle('β-Sn Wulff shapes — 5 low-index faces × 8 methods  '
             f'(common primitive lattice a={A_DFT:.3f}, c={C_DFT:.3f} Å)',
             fontsize=11, y=1.00)
plt.tight_layout(rect=[0, 0.03, 1, 0.98])
for ext, dpi in [('png', 150), ('tif', 300)]:
    plt.savefig(os.path.join(FIG, f'supp_wulff_8methods.{ext}'),
                dpi=dpi, bbox_inches='tight')
print('wrote supp_wulff_8methods.{png,tif}')
plt.close(fig)

# -------- Per-method Wulff figures (one PNG/TIF per calculation method) --------
def slugify(method_name):
    """Convert method label to a filesystem-safe slug.

    'DFT/PBE (BFGS)'    -> 'DFT_PBE'
    'PFP/PBE+D3'        -> 'PFP_PBE_D3'
    'PFP/R²SCAN'        -> 'PFP_R2SCAN'
    'PFP/R²SCAN+D3'     -> 'PFP_R2SCAN_D3'
    'MEAM/Ravelo 1997'  -> 'MEAM_Ravelo1997'
    'MEAM/Etesami 2018' -> 'MEAM_Etesami2018'
    'MEAM/Ko 2018'      -> 'MEAM_Ko2018'
    """
    s = method_name.replace('²', '2').replace('+', '_')
    s = s.split(' (')[0]                       # drop trailing parenthetical
    s = s.replace('/', '_').replace(' ', '')   # path separators and spaces
    return s


for m in methods:
    figm = plt.figure(figsize=(7.0, 7.0))
    ax = figm.add_subplot(111, projection='3d')
    sc = make_wulff(gamma[m])
    sc.make_plot(ax, alpha=0.88, linewidth=0.5, colors=WULFF_COLORS)
    ax.view_init(elev=22, azim=35)
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()  # bare Wulff shape — no title, no axes, no annotations

    plt.tight_layout(pad=0)
    slug = slugify(m)
    for ext, dpi in [('png', 150), ('tif', 300)]:
        plt.savefig(os.path.join(FIG, f'supp_wulff_{slug}.{ext}'),
                    dpi=dpi, bbox_inches='tight', pad_inches=0)
    print(f'wrote supp_wulff_{slug}.{{png,tif}}')
    plt.close(figm)

# ============== Fig 3: Parity plot (γ_method vs γ_DFT, 7 non-DFT methods) ==============
fig, ax = plt.subplots(figsize=(7.2, 7.2))
lim = [200, 1700]
ax.plot(lim, lim, 'k--', lw=0.9, alpha=0.6, label='y = x (perfect agreement)')

for m in methods:
    if m == 'DFT/PBE (BFGS)':
        continue
    ax.scatter(dft, gamma[m], s=85, color=colors[m],
               marker=markers[m], edgecolor='black', linewidth=0.4,
               label=f'{m}  (MAE = {mae[m]:.0f} mJ/m²)', zorder=3)

# annotate face labels next to PFP/PBE points as a visual anchor
for k, gd, gp in zip(FACE_KEYS, dft, gamma['PFP/PBE']):
    ax.annotate(f'({k})', xy=(gd, gp), xytext=(6, -10),
                textcoords='offset points', fontsize=7,
                color=colors['PFP/PBE'], alpha=0.85)

ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel(r'$\gamma_{\rm DFT/PBE}$ (mJ/m²)')
ax.set_ylabel(r'$\gamma_{\rm method}$ (mJ/m²)')
ax.set_title('β-Sn surface energies — methods vs DFT/PBE (5 low-index faces)')
ax.legend(loc='upper left', fontsize=8.5, framealpha=0.95)
ax.grid(alpha=0.3)
ax.set_aspect('equal')
plt.tight_layout()
for ext, dpi in [('png', 150), ('tif', 300)]:
    plt.savefig(os.path.join(FIG, f'supp_parity_vs_DFT_8methods.{ext}'),
                dpi=dpi, bbox_inches='tight')
print('wrote supp_parity_vs_DFT_8methods.{png,tif}')
plt.close(fig)

print('\n=== Summary ===')
hdr = f"{'method':>22} " + ' '.join(f'{f:>6}' for f in FACES) + '   MAE   MAPE'
print(hdr)
for m in methods:
    row = f"{m:>22} " + ' '.join(f'{g:>6.1f}' for g in gamma[m])
    if m in mae:
        mape = float(np.mean(np.abs(gamma[m] - dft) / dft) * 100)
        row += f"  {mae[m]:>5.0f}  {mape:>5.1f}%"
    print(row)
