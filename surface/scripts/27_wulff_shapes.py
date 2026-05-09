"""Build Wulff shapes per method using wulffpack, render grid figure."""
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from ase import Atoms
from wulffpack import SingleCrystal

BASE = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb')
FIG = os.path.join(BASE, 'figures')
os.makedirs(FIG, exist_ok=True)

# β-Sn primitive structure (I4₁/amd, tetragonal, a=5.97, c=3.218)
A, C = 5.970, 3.218
BETA_SN = Atoms(
    'Sn4',
    scaled_positions=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]],
    cell=[[A,0,0],[0,A,0],[0,0,C]],
    pbc=True,
)

FAMILIES = [(1,0,0), (0,0,1), (1,1,0), (1,0,1), (1,1,1), (2,1,1), (1,1,2)]

DATA = {
    'DFT/PBE (BFGS)':   {'100': 492.2, '001': 550.5, '110': 526.9, '101': 505.6, '111': 550.0, '211': 538.4, '112': 358.9},
    'PFP/PBE':          {'100': 385.5, '001': 460.1, '110': 487.0, '101': 412.7, '111': 417.1, '211': 424.0, '112': 442.5},
    'PFP/PBE+D3':       {'100': 642.2, '001': 734.4, '110': 737.5, '101': 832.5, '111': 685.1, '211': 704.1, '112': 668.2},
    'PFP/R²SCAN':       {'100': 784.6, '001': 527.5, '110': 599.3, '101': 476.6, '111': 538.8, '211': 529.2, '112': 483.7},
    'PFP/R²SCAN+D3':    {'100': 614.8, '001': 680.4, '110': 730.6, '101': 707.3, '111': 690.1, '211': 680.9, '112': 633.9},
    'MEAM/Ravelo 1997': {'100': 735.0, '001': 896.8, '110': 954.2, '101': 818.3, '111': 909.3, '211': 941.4, '112': 814.9},
    'MEAM/Etesami 2018':{'100': 1376.3,'001': 1299.7,'110': 1587.4,'101': 1354.7,'111': 1384.7,'211': 1224.8,'112': 1303.3},
}

COLORS = {
    (1,0,0): '#4C72B0',  # blue
    (0,0,1): '#55A868',  # green
    (1,1,0): '#C44E52',  # red
    (1,0,1): '#8172B3',  # purple
    (1,1,1): '#CCB974',  # yellow
    (2,1,1): '#64B5CD',  # cyan
    (1,1,2): '#FF7F0E',  # orange
}


def make_wulff(method_name, gam_dict):
    surface_energies = {}
    for hkl in FAMILIES:
        tag = ''.join(str(h) for h in hkl)
        g = gam_dict.get(tag)
        if g is not None:
            surface_energies[hkl] = g  # mJ/m² (relative values matter only for shape)
    if len(surface_energies) < 3:
        return None
    try:
        sc = SingleCrystal(surface_energies=surface_energies,
                           primitive_structure=BETA_SN,
                           natoms=2000)
        return sc
    except Exception as e:
        print(f'  {method_name}: SingleCrystal failed: {e}')
        return None


def plot_grid():
    methods = list(DATA.keys())
    n = len(methods)
    cols = 4
    rows = (n + cols - 1) // cols
    fig = plt.figure(figsize=(4.5 * cols, 4.5 * rows))

    for i, m in enumerate(methods, 1):
        ax = fig.add_subplot(rows, cols, i, projection='3d')
        sc = make_wulff(m, DATA[m])
        if sc is None:
            ax.text2D(0.5, 0.5, 'insufficient data', ha='center', transform=ax.transAxes)
            ax.set_title(m + ' (n/a)', fontsize=10)
            continue
        try:
            sc.make_plot(ax, alpha=0.85, linewidth=0.4, colors=COLORS)
        except Exception as e:
            print(f'  {m}: make_plot failed: {e}')
        # annotate faces present
        gam = DATA[m]
        txt = ', '.join(f"({tag}):{gam[tag]:.0f}" for tag in ['100','001','110','101','111','211','112']
                        if gam.get(tag) is not None)
        ax.set_title(f"{m}\n{txt}", fontsize=7)
        # consistent view & limits
        ax.view_init(elev=22, azim=35)
        ax.set_box_aspect([1, 1, 1])

    # shared legend from dummy handles
    from matplotlib.patches import Patch
    handles = [Patch(color=c, label='({},{},{})'.format(*k)) for k, c in COLORS.items()]
    fig.legend(handles=handles, loc='lower center', ncol=len(COLORS), fontsize=9,
               bbox_to_anchor=(0.5, -0.02))

    plt.suptitle('β-Sn Wulff shapes — 7 methods  (common primitive lattice a=5.97, c=3.218 Å)',
                 fontsize=12, y=1.00)
    plt.tight_layout(rect=[0, 0.04, 1, 0.98])
    outp_png = os.path.join(FIG, 'wulff_all_methods.png')
    outp_tif = os.path.join(FIG, 'wulff_all_methods.tif')
    plt.savefig(outp_png, dpi=150, bbox_inches='tight')
    plt.savefig(outp_tif, dpi=300, bbox_inches='tight')
    print(f'saved -> {outp_png}')


if __name__ == '__main__':
    plot_grid()
