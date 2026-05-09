"""
12_gb_main_calc.py (run on Matlantis)

Main GB energy calculation for 6 GBs:
  - Σ5(210)[001], Σ5(310)[001], Σ13(320)[001], Σ17(410)[001] symmetric tilts
  - (101), (301) "twins" — mirror construction (Type-I approximation for tetragonal a≠c)
  - (012)[100] tilt skipped (a≠c makes orthogonal CSL impossible; needs special handling)

For each GB, 4 representative RBT shifts are tried and the minimum-energy γ is reported.
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

DEDUP_TOL = 0.50  # Å

GB_LIST = [
    {'name': 'sigma5_210_001',  'P': [[2,1,0],[-1,2,0],[0,0,1]], 'N': 3, 'type': 'tilt_001'},
    {'name': 'sigma5_310_001',  'P': [[3,1,0],[-1,3,0],[0,0,1]], 'N': 2, 'type': 'tilt_001'},
    {'name': 'sigma13_320_001', 'P': [[3,2,0],[-2,3,0],[0,0,1]], 'N': 2, 'type': 'tilt_001'},
    {'name': 'sigma17_410_001', 'P': [[4,1,0],[-1,4,0],[0,0,1]], 'N': 2, 'type': 'tilt_001'},
    {'name': 'twin_101_mirror', 'P': [[1,0,1],[0,1,0],[-1,0,1]], 'N': 3, 'type': 'twin_mirror',
     'note': 'Type-I mirror approximation; true β-Sn (101) twin is Type-II'},
    {'name': 'twin_301_mirror', 'P': [[3,0,1],[0,1,0],[-1,0,3]], 'N': 2, 'type': 'twin_mirror',
     'note': 'Type-I mirror approximation; true β-Sn (301) twin is Type-II'},
]

RBT_LIST = [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)]


def build_gb_with_rbt(bulk, P, N, rbt_b, rbt_c):
    """Symmetric tilt / mirror-twin GB with rigid body translation (fractional in b', c')."""
    csl = make_supercell(bulk, P)
    g1 = csl.repeat((N, 1, 1))
    combined_cell = g1.cell.array.copy()
    combined_cell[0] = 2.0 * combined_cell[0]

    scaled1 = g1.get_scaled_positions().copy()
    scaled1[:, 0] *= 0.5  # grain 1 at [0, 0.5]
    scaled2 = scaled1.copy()
    scaled2[:, 0] = 1.0 - scaled2[:, 0]   # mirror through f_a = 0.5
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
    print(f'bulk E/atom = {e_bulk_per_atom:.4f} eV', flush=True)

    est = Estimator(calc_mode='PBE_PLUS_D3')
    calc = ASECalculator(est)

    all_results = []
    t0 = time.time()

    for gb in GB_LIST:
        name = gb['name']
        best = {'gamma_mJ_m2': 1e9}
        scan = []
        print(f'\n### {name} (type {gb["type"]}, N={gb["N"]}) ###', flush=True)

        for (rb, rc) in RBT_LIST:
            label = f"{name}_rbt{rb:.1f}_{rc:.1f}"
            atoms, area = build_gb_with_rbt(bulk, gb['P'], gb['N'], rb, rc)
            n0 = len(atoms)
            atoms.calc = calc
            logfile = os.path.join(LOGDIR, f'{label}.log')
            try:
                opt = LBFGS(atoms, logfile=logfile)
                opt.run(fmax=0.05, steps=300)
            except Exception as e:
                print(f'  [RBT {rb:.1f},{rc:.1f}] FAILED: {e}', flush=True)
                scan.append({'rbt': (rb, rc), 'error': str(e)})
                continue
            E = atoms.get_potential_energy()
            n = len(atoms)
            gamma = (E - n * e_bulk_per_atom) / (2 * area) * EV_A2_TO_MJ_M2
            row = {'rbt': [rb, rc], 'n_atoms': n, 'area_A2': area,
                   'E_eV': E, 'gamma_mJ_m2': gamma}
            scan.append(row)
            print(f'  RBT ({rb:.1f},{rc:.1f})  n={n}/{n0}  γ={gamma:7.1f} mJ/m^2 (elapsed {time.time()-t0:.0f}s)', flush=True)
            if gamma < best['gamma_mJ_m2']:
                best = dict(row)
                best['label'] = label
                # save best traj
                write(os.path.join(RESDIR, f'{name}_best.traj'), atoms)

        entry = {'name': name, 'type': gb['type'],
                 'P': gb['P'], 'N': gb['N'],
                 'note': gb.get('note', ''),
                 'scan': scan, 'best': best}
        all_results.append(entry)
        with open(os.path.join(RESDIR, 'main_results.json'), 'w') as f:
            json.dump(all_results, f, indent=2)

    # summary
    print('\n' + '=' * 66, flush=True)
    print(f"{'GB':>22} {'type':>14} {'γ_min':>8} {'best RBT':>12}", flush=True)
    print('-' * 66, flush=True)
    for e in all_results:
        b = e['best']
        if 'gamma_mJ_m2' not in b:
            print(f"{e['name']:>22}  FAILED all RBT")
            continue
        print(f"{e['name']:>22} {e['type']:>14} {b['gamma_mJ_m2']:>8.1f} "
              f"({b['rbt'][0]:.1f},{b['rbt'][1]:.1f})", flush=True)
    print(f'wall: {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
