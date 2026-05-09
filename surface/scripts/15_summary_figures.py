"""
15_summary_figures.py — Build summary figures + table for β-Sn surface & GB energies.

Outputs:
  figures/surface_energy_comparison.png/tif  — γ_surf PFP vs DFT
  figures/gb_energy_bar.png/tif              — γ_GB all grain boundaries
  figures/gb_vs_misorientation.png/tif       — γ_GB vs θ for [001] tilt series
  summary_table.md                           — human-readable
  summary_all.json                           — machine-readable
"""
import json
import os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb')
FIG = os.path.join(BASE, 'figures')
os.makedirs(FIG, exist_ok=True)

# ---------- data ---------- #
# Surface energies (mJ/m²)
# γ_PFP: PFP/PBE+D3 at PFP-relaxed geometry
# γ_DFT_SP: OpenMX/PBE single-point on PFP-relaxed slabs (intel, Δk=0.15, SCF 1e-7)
SURF = [
    {'face': '(001)', 'L': 8,  'pfp': 736.2, 'dft_sp': 579.6},
    {'face': '(110)', 'L': 8,  'pfp': 737.5, 'dft_sp': 532.5},
    {'face': '(100)', 'L': 12, 'pfp': 834.5, 'dft_sp': 619.7},
    {'face': '(101)', 'L': 12, 'pfp': 833.1, 'dft_sp': 624.1},
]

# GB energies (mJ/m²), PFP/PBE+D3
# Σ tilt: 4 RBT positions → minimum
# Twins: 10×10 RBT γ-surface → physical minimum (excluding bulk-recovery γ<20)
GB = [
    # [001] symmetric tilts, misorientation angle θ (degrees)
    {'name': 'Σ5(210)[001]',  'type': 'tilt_001',   'theta': 36.87, 'gamma': 181.2, 'Σ': 5,  'gb_plane': '(210)'},
    {'name': 'Σ17(410)[001]', 'type': 'tilt_001',   'theta': 28.07, 'gamma': 300.4, 'Σ': 17, 'gb_plane': '(410)'},
    {'name': 'Σ13(320)[001]', 'type': 'tilt_001',   'theta': 22.62, 'gamma': 373.7, 'Σ': 13, 'gb_plane': '(320)'},
    {'name': 'Σ5(310)[001]',  'type': 'tilt_001',   'theta': 53.13, 'gamma': 369.7, 'Σ': 5,  'gb_plane': '(310)'},
    # compound twins
    {'name': '(101) twin',    'type': 'twin',       'theta': None,  'gamma': 126.4, 'note': 'compound, γ-surface physical min'},
    {'name': '(301) twin',    'type': 'twin',       'theta': None,  'gamma': 21.6,  'note': 'compound, γ-surface physical min — β-Sn "tin cry" twin'},
]

# Literature / experimental references (mJ/m²)
LIT_SURF = {'avg_experimental': 670, 'src': 'Tyson–Miller 1977, semi-empirical avg for β-Sn'}

# ---------- figure 1: surface energies ---------- #
def fig_surface():
    labels = [s['face'] for s in SURF]
    gpfp = [s['pfp'] for s in SURF]
    gdft = [s['dft_sp'] for s in SURF]
    x = np.arange(len(labels))
    w = 0.35

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.bar(x - w/2, gpfp, w, label='PFP/PBE+D3 (relaxed)', color='#3366cc')
    ax.bar(x + w/2, gdft, w, label='OpenMX/PBE SP (on PFP geometry)', color='#dc3912')
    ax.axhline(LIT_SURF['avg_experimental'], linestyle='--', color='k', lw=1,
               label=f"exp. avg ~{LIT_SURF['avg_experimental']} mJ/m² [Tyson–Miller 1977]")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel('surface energy γ (mJ/m²)')
    ax.set_title('β-Sn surface energies — PFP vs intel DFT(SP)')
    ax.legend(loc='upper left', fontsize=9)
    ax.set_ylim(0, max(gpfp) * 1.18)
    for xi, v in zip(x - w/2, gpfp):
        ax.text(xi, v + 8, f'{v:.0f}', ha='center', fontsize=8)
    for xi, v in zip(x + w/2, gdft):
        ax.text(xi, v + 8, f'{v:.0f}', ha='center', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'surface_energy_comparison.png'), dpi=150)
    plt.savefig(os.path.join(FIG, 'surface_energy_comparison.tif'), dpi=300)
    print('wrote surface_energy_comparison.{png,tif}')
    plt.close(fig)


# ---------- figure 2: GB energy bar chart ---------- #
def fig_gb_bar():
    names = [g['name'] for g in GB]
    gammas = [g['gamma'] for g in GB]
    types = [g['type'] for g in GB]
    colors = ['#4C72B0' if t == 'tilt_001' else '#DD8452' for t in types]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    bars = ax.bar(names, gammas, color=colors)
    for bar, v in zip(bars, gammas):
        ax.text(bar.get_x() + bar.get_width()/2, v + 8, f'{v:.1f}',
                ha='center', fontsize=9)
    ax.set_ylabel('grain boundary energy γ_GB (mJ/m²)')
    ax.set_title('β-Sn grain boundary energies (PFP/PBE+D3)')
    ax.set_ylim(0, max(gammas) * 1.18)
    # legend
    from matplotlib.patches import Patch
    leg = [Patch(color='#4C72B0', label='[001] symmetric tilt (rep. RBT min)'),
           Patch(color='#DD8452', label='compound twin (γ-surface phys. min)')]
    ax.legend(handles=leg, loc='upper right', fontsize=9)
    plt.xticks(rotation=20, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'gb_energy_bar.png'), dpi=150)
    plt.savefig(os.path.join(FIG, 'gb_energy_bar.tif'), dpi=300)
    print('wrote gb_energy_bar.{png,tif}')
    plt.close(fig)


# ---------- figure 3: γ_GB vs misorientation angle ---------- #
def fig_gb_misorientation():
    tilts = [g for g in GB if g['type'] == 'tilt_001']
    tilts_sorted = sorted(tilts, key=lambda g: g['theta'])
    theta = [g['theta'] for g in tilts_sorted]
    gamma = [g['gamma'] for g in tilts_sorted]
    names = [g['name'] for g in tilts_sorted]

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(theta, gamma, 'o-', color='#4C72B0', ms=10, lw=1.5)
    for t, g, n in zip(theta, gamma, names):
        ax.annotate(n, (t, g), xytext=(5, 8), textcoords='offset points', fontsize=8)
    ax.set_xlabel('misorientation angle θ (°)')
    ax.set_ylabel('γ_GB (mJ/m²)')
    ax.set_title('β-Sn symmetric tilt [001] grain boundaries vs θ (PFP/PBE+D3)')
    ax.set_xlim(15, 60)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'gb_vs_misorientation.png'), dpi=150)
    plt.savefig(os.path.join(FIG, 'gb_vs_misorientation.tif'), dpi=300)
    print('wrote gb_vs_misorientation.{png,tif}')
    plt.close(fig)


# ---------- markdown summary ---------- #
def write_markdown():
    md = []
    md.append('# β-Sn Surface & Grain Boundary Energy — Summary')
    md.append('')
    md.append('## 1. Surface energies (mJ/m²)')
    md.append('')
    md.append('| face | slab layers | γ_PFP (PBE+D3) | γ_DFT (OpenMX/PBE SP at PFP geom) | Δ (DFT–PFP) |')
    md.append('|------|-------------|---------------:|----------------------------------:|------------:|')
    for s in SURF:
        d = s['dft_sp'] - s['pfp']
        md.append(f"| {s['face']} | {s['L']} | {s['pfp']:.1f} | {s['dft_sp']:.1f} | {d:+.1f} |")
    md.append('')
    md.append(f"- Experimental β-Sn average: ~{LIT_SURF['avg_experimental']} mJ/m² "
              f"({LIT_SURF['src']}).")
    md.append('- DFT SP consistently under-predicts vs PFP+D3 (by ~20–30%).')
    md.append('- DFT fully relaxed values pending on supercomputer jobs.')
    md.append('')
    md.append('## 2. Grain boundary energies (PFP/PBE+D3, mJ/m²)')
    md.append('')
    md.append('| GB | type | θ (°) | γ_GB | method |')
    md.append('|----|------|------:|-----:|--------|')
    for g in GB:
        theta = f"{g['theta']:.2f}" if g.get('theta') is not None else '—'
        method = 'RBT 4-pt scan → min' if g['type'] == 'tilt_001' else 'γ-surface 10×10 → physical min'
        md.append(f"| {g['name']} | {g['type']} | {theta} | {g['gamma']:.1f} | {method} |")
    md.append('')
    md.append('Notes:')
    md.append('- [001] symmetric tilt convergence checked on Σ5(210) (N=3 gives γ within ±5 mJ/m² of N=8).')
    md.append('- Compound twins constructed via mirror; RBT space includes many "bulk-recovery" points (γ<20 mJ/m²) due to BCT body-center symmetry — only "physical" minima reported.')
    md.append('- (301) γ_GB = 21.6 mJ/m² is at the lower tail: consistent with β-Sn being the classic deformation-twin material ("tin cry").')
    md.append('- (012)[100] tilt skipped: a≠c tetragonality precludes orthogonal CSL; will require Type-II treatment.')
    md.append('')
    md.append('## 3. Method settings')
    md.append('')
    md.append('- **Engine**: PFP v8 (calc_mode `PBE_PLUS_D3`) on Matlantis')
    md.append('- **Bulk lattice** (PFP-optimized): a = 5.8457 Å, c = 3.1732 Å')
    md.append('- **Surface slabs**: 4 orientations, 15 Å vacuum, middle 30% atoms frozen during PFP relaxation')
    md.append('- **GB construction**: ase.build.make_supercell with CSL basis; symmetric tilt via mirror + RBT; dedup tol 0.5 Å (MIC)')
    md.append('- **Relaxation**: LBFGS, fmax 0.05 eV/Å, 200–300 step cap')
    md.append('- **γ formula**: γ = (E_cell − N·E_bulk/atom) / (2·A)')

    outp = os.path.join(BASE, 'summary_table.md')
    with open(outp, 'w') as f:
        f.write('\n'.join(md))
    print(f'wrote {outp}')


def write_json():
    outp = os.path.join(BASE, 'summary_all.json')
    with open(outp, 'w') as f:
        json.dump({
            'surface': SURF,
            'gb': GB,
            'references': {'exp_surface_avg_mJ_m2': LIT_SURF['avg_experimental'],
                           'source': LIT_SURF['src']},
        }, f, indent=2)
    print(f'wrote {outp}')


if __name__ == '__main__':
    fig_surface()
    fig_gb_bar()
    fig_gb_misorientation()
    write_markdown()
    write_json()
