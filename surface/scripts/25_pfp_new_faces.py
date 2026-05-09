"""PFP γ for additional faces (111), (211), (112) in 4 modes.
Uses cached mode-self-consistent bulk from previous run.
"""
import json
import os
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.io import read
from ase.optimize import LBFGS
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

EV_A2_TO_MJ_M2 = 16021.766
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUTDIR = os.path.join(ROOT, 'results_modes')

MODES = ['PBE', 'PBE_PLUS_D3', 'R2SCAN', 'R2SCAN_PLUS_D3']
NEW_FACES = [('111', 12), ('211', 12), ('112', 12)]
VACUUM = 15.0
FMAX = 0.015


def main():
    # load prior results
    inp = os.path.join(OUTDIR, 'modes_surface_energies.json')
    with open(inp) as f:
        all_results = json.load(f)

    for mode in MODES:
        print(f'\n##### {mode} #####', flush=True)
        bulk = read(os.path.join(OUTDIR, f'bulk_{mode}.traj'))
        e_per_atom = all_results[mode]['bulk']['E_per_atom_eV']
        est = Estimator(calc_mode=mode)
        calc = ASECalculator(est)

        for tag, L in NEW_FACES:
            miller = tuple(int(x) for x in tag)
            slab = ase_surface(bulk, miller, layers=L, vacuum=VACUUM, periodic=True)
            slab.center(vacuum=VACUUM, axis=2)
            slab.calc = calc
            n = len(slab)
            cell = slab.cell.array
            area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
            try:
                LBFGS(slab, logfile=None).run(fmax=FMAX, steps=300)
                E = float(slab.get_potential_energy())
                gamma = (E - n * e_per_atom) / (2 * area) * EV_A2_TO_MJ_M2
                print(f'  ({tag}) L={L}  n={n}  A={area:.2f}  γ={gamma:.1f} mJ/m²', flush=True)
                all_results[mode]['slabs'].append({
                    'face': tag, 'L': L, 'n_atoms': n, 'area_A2': area,
                    'E_slab_eV': E, 'gamma_mJ_m2': gamma})
            except Exception as e:
                print(f'  ({tag}) FAILED: {e}', flush=True)

    with open(inp, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f'\nupdated -> {inp}', flush=True)

    # summary
    print('\nγ (mJ/m²):')
    face_cols = ['001','110','100','101','111','211','112']
    print(f"{'mode':>15} " + ' '.join(f'{f:>7}' for f in face_cols))
    for mode, r in all_results.items():
        by_face = {s['face']: s['gamma_mJ_m2'] for s in r['slabs']}
        row = f"{mode:>15} " + ' '.join(f"{by_face.get(f, 0):>7.1f}" for f in face_cols)
        print(row)


if __name__ == '__main__':
    main()
