"""
13_twin_gamma_surface.py (run on Matlantis)

γ-surface (10×10 RBT grid) for β-Sn (101) and (301) compound twins using the
mirror construction. Output: per-GB contour data and best γ_GB.

For each twin:
  - Build CSL supercell with mirror (Type-I / compound construction)
  - Scan rigid-body translation in (b', c') fractional grid 10×10
  - At each RBT, dedup + relax with PFP/PBE+D3
  - Flag "trivial bulk recovery" (γ < 20 mJ/m² suggests mirror + RBT reproduced bulk)
"""
import json
import os
import time
import numpy as np
from ase import Atoms
from ase.build import make_supercell
from ase.io import read, write
from ase.optimize import LBFGS
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

EV_A2_TO_MJ_M2 = 16021.766
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESDIR = os.path.join(ROOT, 'results_gb')
LOGDIR = os.path.join(ROOT, 'logs_gb')
os.makedirs(RESDIR, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

DEDUP_TOL = 0.50

TWINS = [
    {'name': 'twin_101', 'P': [[1,0,1],[0,1,0],[-1,0,1]], 'N': 3},
    {'name': 'twin_301', 'P': [[3,0,1],[0,1,0],[-1,0,3]], 'N': 2},
]

N_GRID = 10  # 10 x 10 RBT samples


def build_mirror_twin(bulk, P, N, rbt_b, rbt_c):
    csl = make_supercell(bulk, P)
    g1 = csl.repeat((N, 1, 1))
    combined_cell = g1.cell.array.copy()
    combined_cell[0] = 2.0 * combined_cell[0]

    scaled1 = g1.get_scaled_positions().copy()
    scaled1[:, 0] *= 0.5
    scaled2 = scaled1.copy()
    scaled2[:, 0] = 1.0 - scaled2[:, 0]
    scaled2[:, 1] = (scaled2[:, 1] + rbt_b) % 1.0
    scaled2[:, 2] = (scaled2[:, 2] + rbt_c) % 1.0

    symbols = list(g1.get_chemical_symbols()) * 2
    combined = Atoms(symbols=symbols,
                     scaled_positions=np.vstack([scaled1, scaled2]),
                     cell=combined_cell, pbc=True)
    combined.wrap()

    d = combined.get_all_distances(mic=True)
    keep = np.ones(len(combined), dtype=bool)
    for i in range(len(combined)):
        if not keep[i]:
            continue
        for j in range(i+1, len(combined)):
            if keep[j] and d[i, j] < DEDUP_TOL:
                keep[j] = False
    combined = combined[keep]

    bvec = combined_cell[1]
    cvec = combined_cell[2]
    area = float(np.linalg.norm(np.cross(bvec, cvec)))
    return combined, area


def main():
    bulk = read(os.path.join(ROOT, 'results', 'bulk_pfp_pbe_d3.traj'))
    with open(os.path.join(ROOT, 'results', 'bulk_pfp_pbe_d3.json')) as f:
        e_bulk_per_atom = json.load(f)['E_per_atom_eV']
    print(f'bulk E/atom = {e_bulk_per_atom:.6f} eV', flush=True)

    est = Estimator(calc_mode='PBE_PLUS_D3')
    calc = ASECalculator(est)

    all_out = []
    for twin in TWINS:
        name = twin['name']
        P = twin['P']
        N = twin['N']
        t0 = time.time()
        print(f'\n##### {name}  (N={N}) — γ-surface {N_GRID}×{N_GRID} #####', flush=True)

        grid = np.zeros((N_GRID, N_GRID))
        natoms = np.zeros((N_GRID, N_GRID), dtype=int)
        trivial = np.zeros((N_GRID, N_GRID), dtype=bool)

        frac_list = np.linspace(0.0, 1.0, N_GRID, endpoint=False)
        for i, rb in enumerate(frac_list):
            for j, rc in enumerate(frac_list):
                atoms, area = build_mirror_twin(bulk, P, N, rb, rc)
                atoms.calc = calc
                try:
                    opt = LBFGS(atoms, logfile=None)
                    opt.run(fmax=0.05, steps=200)
                except Exception as e:
                    grid[i, j] = np.nan
                    continue
                E = atoms.get_potential_energy()
                n = len(atoms)
                gamma = (E - n * e_bulk_per_atom) / (2 * area) * EV_A2_TO_MJ_M2
                grid[i, j] = gamma
                natoms[i, j] = n
                trivial[i, j] = gamma < 20.0  # "near-zero" bulk recovery
            print(f"  rb={rb:.2f} row done  (elapsed {time.time()-t0:.0f}s, γ range this row: "
                  f"{np.nanmin(grid[i,:j+1]):.0f}–{np.nanmax(grid[i,:j+1]):.0f})", flush=True)

        # summary
        valid = grid[~trivial & ~np.isnan(grid)]
        print(f"\n{name}: full grid γ range = {np.nanmin(grid):.1f} – {np.nanmax(grid):.1f} mJ/m²", flush=True)
        print(f"  number of trivial (bulk-recovery, γ<20) points: {trivial.sum()}", flush=True)
        if valid.size > 0:
            print(f"  γ_min (excluding trivial) = {valid.min():.1f} mJ/m²", flush=True)
            idx = np.unravel_index(np.argmin(np.where(trivial, np.inf, grid)), grid.shape)
            best_rbt = (frac_list[idx[0]], frac_list[idx[1]])
            print(f"  best RBT: ({best_rbt[0]:.2f}, {best_rbt[1]:.2f})", flush=True)
        else:
            best_rbt = (float('nan'), float('nan'))
            print(f"  WARNING: entire grid is trivial bulk recovery!", flush=True)

        np.savez(os.path.join(RESDIR, f'{name}_gamma_surface.npz'),
                 frac=frac_list, gamma=grid, natoms=natoms, trivial=trivial)
        all_out.append({
            'name': name, 'P': P, 'N': N,
            'grid_frac': frac_list.tolist(),
            'gamma_grid': grid.tolist(),
            'trivial_mask': trivial.tolist(),
            'gamma_min_physical': float(valid.min()) if valid.size > 0 else None,
            'best_rbt': list(best_rbt),
            'n_trivial': int(trivial.sum()),
        })
        with open(os.path.join(RESDIR, 'twin_gamma_surface.json'), 'w') as f:
            json.dump(all_out, f, indent=2)
        print(f'  wall: {time.time()-t0:.0f}s', flush=True)

    print('\n=== TWIN γ-SURFACE SUMMARY ===', flush=True)
    for e in all_out:
        print(f"{e['name']}: γ_min (physical) = {e['gamma_min_physical']:.1f} mJ/m² at RBT "
              f"({e['best_rbt'][0]:.2f}, {e['best_rbt'][1]:.2f}); "
              f"trivial pts: {e['n_trivial']}/{N_GRID**2}", flush=True)


if __name__ == '__main__':
    main()
