"""MEAM γ for additional faces (111), (211), (112) — Ravelo1997 + Etesami2018."""
import os
import json
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.calculators.lammpslib import LAMMPSlib
from ase.filters import ExpCellFilter
from ase.optimize import LBFGS

POT_DIR = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/meam_potentials')
EV_A2_TO_MJ_M2 = 16021.766
VACUUM = 15.0

POTENTIALS = {
    'Ravelo1997':  ('library.Sn.Ravelo1997.meam',  'Sn.Ravelo1997.meam'),
    'Etesami2018': ('library.Sn.Etesami2018.meam', 'Sn.Etesami2018.meam'),
}

EXP_LAT = {'a': 5.831, 'c': 3.182}
NEW_FACES = [('111', 12), ('211', 12), ('112', 12)]


def build_beta_sn(a, c):
    return Atoms('Sn4',
                 scaled_positions=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]],
                 cell=[[a,0,0],[0,a,0],[0,0,c]], pbc=True)


def meam_calc(name, lib, pot):
    cmds = [
        f'pair_style meam',
        f'pair_coeff * * {os.path.join(POT_DIR, lib)} Sn {os.path.join(POT_DIR, pot)} Sn',
    ]
    return LAMMPSlib(lmpcmds=cmds, atom_types={'Sn': 1},
                     log_file=f'/tmp/lammps_{name}.log', keep_alive=True)


def main():
    out = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/meam_surface_energies.json')
    with open(out) as f:
        all_results = json.load(f)

    for name, (lib, pot) in POTENTIALS.items():
        print(f'\n##### {name} #####', flush=True)
        # reuse relaxed bulk E_per_atom from previous run
        bulk_info = all_results[name]['bulk']
        a, c = bulk_info['a'], bulk_info['c']
        e_per_atom = bulk_info['E_per_atom_eV']
        print(f'  bulk: a={a:.4f} c={c:.4f} E/atom={e_per_atom:.4f} eV (cached)', flush=True)

        bulk = build_beta_sn(a, c)
        bulk.calc = meam_calc(name + '_b', lib, pot)
        _ = bulk.get_potential_energy()  # sanity

        for tag, L in NEW_FACES:
            miller = tuple(int(x) for x in tag)
            slab = ase_surface(bulk, miller, layers=L, vacuum=VACUUM, periodic=True)
            slab.center(vacuum=VACUUM, axis=2)
            slab.calc = meam_calc(name + '_s_' + tag, lib, pot)
            try:
                LBFGS(slab, logfile=None).run(fmax=0.01, steps=300)
                E = float(slab.get_potential_energy())
                n = len(slab)
                cell = slab.cell.array
                A = float(np.linalg.norm(np.cross(cell[0], cell[1])))
                gamma = (E - n * e_per_atom) / (2 * A) * EV_A2_TO_MJ_M2
                print(f'  ({tag}) L={L}  n={n}  A={A:.3f}  γ={gamma:.1f} mJ/m²', flush=True)
                all_results[name]['slabs'].append({
                    'face': tag, 'L': L, 'n': n, 'area_A2': A, 'E_eV': E, 'gamma_mJ_m2': gamma})
            except Exception as e:
                print(f'  ({tag}) FAILED: {e}', flush=True)

    with open(out, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f'\nupdated -> {out}', flush=True)

    # summary
    print('\nγ (mJ/m²):')
    print(f"{'method':>12} " + ' '.join(f'{f:>7}' for f in ['001','110','100','101','111','211','112']))
    for name, r in all_results.items():
        by_face = {s['face']: s['gamma_mJ_m2'] for s in r['slabs']}
        row = f"{name:>12} "
        for f in ['001','110','100','101','111','211','112']:
            row += f"{by_face.get(f, 0):>7.1f} "
        print(row)


if __name__ == '__main__':
    main()
