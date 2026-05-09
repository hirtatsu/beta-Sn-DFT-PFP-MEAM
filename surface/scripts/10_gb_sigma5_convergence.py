"""
10_gb_sigma5_convergence.py (run on Matlantis)

Σ5 (210) [001] symmetric tilt grain boundary for β-Sn — convergence test
of γ_GB vs. thickness perpendicular to GB plane.

Construction:
  - Cut bulk by CSL basis: a'=(2,1,0), b'=(-1,2,0), c'=(0,0,1)
  - Repeat N times along a' (grain 1)
  - Create combined cell of length 2N along a'
  - Grain 1 atoms at fractional [0, 0.5] along a
  - Grain 2 atoms = mirror of grain 1 through f_a = 0.5 (→ grain 2 at [0.5, 1])
  - Two GBs per periodic cell (f_a = 0 and 0.5)
  - Dedup with MIC-aware distance (handles non-orthorhombic cell correctly)
  - γ_GB = (E − N_atoms × E_bulk_per_atom) / (2 × A),  A = |b'| × c

Robust dedup: atoms within tol merged (keep lower index)
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

DEDUP_TOL = 0.50  # Å — merge atoms closer than this (PBC-aware)


def build_gb_cell(bulk, N):
    # CSL supercell (Σ=5) — use make_supercell (exact, preserves atom count)
    P = [[2, 1, 0], [-1, 2, 0], [0, 0, 1]]
    g1 = make_supercell(bulk, P)
    g1 = g1.repeat((N, 1, 1))

    combined_cell = g1.cell.array.copy()
    combined_cell[0] = 2.0 * combined_cell[0]  # double along a'

    scaled1 = g1.get_scaled_positions().copy()
    scaled1[:, 0] *= 0.5  # grain 1 → f_a ∈ [0, 0.5]
    scaled2 = scaled1.copy()
    scaled2[:, 0] = 1.0 - scaled2[:, 0]  # mirror through GB plane f_a = 0.5

    symbols = list(g1.get_chemical_symbols()) * 2
    combined = Atoms(symbols=symbols,
                     scaled_positions=np.vstack([scaled1, scaled2]),
                     cell=combined_cell, pbc=True)
    combined.wrap()

    # Proper MIC-aware deduplication
    dists = combined.get_all_distances(mic=True)
    keep = np.ones(len(combined), dtype=bool)
    for i in range(len(combined)):
        if not keep[i]:
            continue
        for j in range(i + 1, len(combined)):
            if keep[j] and dists[i, j] < DEDUP_TOL:
                keep[j] = False
    combined = combined[keep]

    # GB in-plane area: |b'| × c (perpendicular to GB normal a')
    bvec = combined_cell[1]
    cvec = combined_cell[2]
    area = float(np.linalg.norm(np.cross(bvec, cvec)))
    la = float(np.linalg.norm(combined_cell[0]))
    lb = float(np.linalg.norm(bvec))
    lc = float(np.linalg.norm(cvec))
    return combined, la, lb, lc, area


def main():
    bulk = read(os.path.join(ROOT, 'results', 'bulk_pfp_pbe_d3.traj'))
    with open(os.path.join(ROOT, 'results', 'bulk_pfp_pbe_d3.json')) as f:
        e_bulk_per_atom = json.load(f)['E_per_atom_eV']
    print(f'bulk: E_per_atom = {e_bulk_per_atom:.4f} eV', flush=True)

    est = Estimator(calc_mode='PBE_PLUS_D3')
    calc = ASECalculator(est)

    thickness_list = [2, 3, 4, 6, 8]
    results = []
    t0 = time.time()

    for N in thickness_list:
        label = f'sigma5_210_001_N{N}'
        print(f'\n=== {label}  (elapsed {time.time()-t0:.0f}s) ===', flush=True)
        gb_atoms, la, lb, lc, area = build_gb_cell(bulk, N)
        print(f'  built n={len(gb_atoms)}, cell(a,b,c)=({la:.2f},{lb:.2f},{lc:.2f}) Å, GB area={area:.2f}', flush=True)

        # report min pairwise distance (sanity)
        dmin = None
        if len(gb_atoms) > 1:
            d = gb_atoms.get_all_distances(mic=True)
            np.fill_diagonal(d, np.inf)
            dmin = float(d.min())
        print(f'  min pair dist: {dmin:.3f} Å', flush=True)

        gb_atoms.calc = calc
        E_init = gb_atoms.get_potential_energy()
        print(f'  E_init = {E_init:.4f} eV', flush=True)

        logfile = os.path.join(LOGDIR, f'{label}.log')
        opt = LBFGS(gb_atoms, logfile=logfile)
        try:
            opt.run(fmax=0.05, steps=300)
        except Exception as e:
            print(f'  opt failed: {e}', flush=True)
            results.append({'N': N, 'error': str(e), 'E_init': E_init})
            continue

        E = gb_atoms.get_potential_energy()
        n = len(gb_atoms)
        gamma = (E - n * e_bulk_per_atom) / (2.0 * area) * EV_A2_TO_MJ_M2
        print(f'  E_relax = {E:.4f} eV, γ_GB = {gamma:.1f} mJ/m²', flush=True)
        write(os.path.join(RESDIR, f'{label}.traj'), gb_atoms)
        results.append({
            'N': N, 'label': label, 'n_atoms': n,
            'cell_Ang': [la, lb, lc], 'area_A2': area,
            'E_init_eV': E_init, 'E_relax_eV': E, 'gamma_mJ_m2': gamma,
            'min_pair_dist_initial': dmin,
        })
        with open(os.path.join(RESDIR, 'sigma5_convergence.json'), 'w') as f:
            json.dump(results, f, indent=2)

    print('\n=== Σ5(210)[001] thickness convergence ===', flush=True)
    print(f"{'N':>3} {'n':>5} {'la(Å)':>8} {'γ_GB (mJ/m²)':>14}", flush=True)
    for r in results:
        if 'error' in r:
            print(f"{r['N']:>3} ERROR")
            continue
        print(f"{r['N']:>3} {r['n_atoms']:>5} {r['cell_Ang'][0]:>8.2f} {r['gamma_mJ_m2']:>14.1f}", flush=True)
    print(f'wall: {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
