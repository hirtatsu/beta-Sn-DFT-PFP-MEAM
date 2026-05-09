"""Plot γ-surface contour maps for (101) and (301) twins."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

TMP = os.path.expanduser('~/tmp')
FIG = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/figures')
os.makedirs(FIG, exist_ok=True)

TWINS = [
    ('twin_101', 'β-Sn (101) compound twin — mirror construction'),
    ('twin_301', 'β-Sn (301) compound twin — mirror construction'),
]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (name, title) in zip(axes, TWINS):
    data = np.load(os.path.join(TMP, f'{name}_gamma_surface.npz'))
    frac = data['frac']
    gamma = data['gamma']
    trivial = data['trivial']

    # mask: for visualization, replace trivial points (γ<20) with NaN or different color
    masked_gamma = gamma.copy()
    # keep original values but mark
    # Use periodic tiling for a nicer contour (2x2 copies)
    extended = np.tile(gamma, (2, 2))
    frac_ext = np.concatenate([frac, frac + 1.0])

    # contour
    X, Y = np.meshgrid(frac_ext, frac_ext, indexing='ij')
    vmin = max(1, np.nanmin(gamma))
    vmax = np.nanmax(gamma)
    cs = ax.contourf(X, Y, extended, levels=20, cmap='viridis',
                     vmin=vmin, vmax=vmax)

    # overlay bulk-recovery points
    for i, rb in enumerate(frac):
        for j, rc in enumerate(frac):
            if trivial[i, j]:
                for dx in (0, 1):
                    for dy in (0, 1):
                        ax.plot(rb + dx, rc + dy, 'rx', ms=8, mew=2)

    # mark best (min excluding trivial)
    mask_valid = ~trivial & ~np.isnan(gamma)
    if mask_valid.any():
        g_for_min = np.where(mask_valid, gamma, np.inf)
        idx = np.unravel_index(np.argmin(g_for_min), gamma.shape)
        gmin = gamma[idx]
        for dx in (0, 1):
            for dy in (0, 1):
                ax.plot(frac[idx[0]] + dx, frac[idx[1]] + dy,
                        'w*', ms=16, mec='black', mew=1.5)

    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)
    ax.set_xlabel('RBT along b′ (fractional, tiled 2×)')
    ax.set_ylabel('RBT along c′ (fractional, tiled 2×)')
    ax.set_aspect('equal')
    fig.colorbar(cs, ax=ax, label='γ (mJ/m²)', shrink=0.8)
    vmin_lbl = np.nanmin(gamma)
    vmax_lbl = np.nanmax(gamma)
    ax.set_title(f'{title}\nγ range: {vmin_lbl:.1f}–{vmax_lbl:.1f} mJ/m², '
                 f'γ_min(physical) = {gmin:.1f}\n'
                 f'red × : bulk-recovery (γ<20), white ★ : physical minimum')

plt.tight_layout()
outp = os.path.join(FIG, 'twin_gamma_surface.png')
plt.savefig(outp, dpi=150, bbox_inches='tight')
print(f'saved {outp}')
