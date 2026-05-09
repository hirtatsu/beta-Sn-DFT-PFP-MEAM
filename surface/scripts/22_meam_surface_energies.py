"""MEAM surface energy benchmark for β-Sn — Ravelo1997 & Etesami2018."""
import os
import json
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.calculators.lammpslib import LAMMPSlib
from ase.filters import ExpCellFilter
from ase.io import write
from ase.optimize import LBFGS

POT_DIR = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/meam_potentials')
EV_A2_TO_MJ_M2 = 16021.766
VACUUM = 15.0

POTENTIALS = {
    'Ravelo1997':  ('library.Sn.Ravelo1997.meam',  'Sn.Ravelo1997.meam'),
    'Etesami2018': ('library.Sn.Etesami2018.meam', 'Sn.Etesami2018.meam'),
}

EXP_LAT = {'a': 5.831, 'c': 3.182}
FACES = [('001', 8), ('110', 8), ('100', 12), ('101', 12)]


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


def relax_bulk(name, lib, pot):
    atoms = build_beta_sn(EXP_LAT['a'], EXP_LAT['c'])
    atoms.calc = meam_calc(name, lib, pot)
    ecf = ExpCellFilter(atoms, hydrostatic_strain=False)
    LBFGS(ecf, logfile=None).run(fmax=0.001, steps=300)
    a, c = float(atoms.cell.lengths()[0]), float(atoms.cell.lengths()[2])
    e_per_atom = float(atoms.get_potential_energy()) / len(atoms)
    return atoms, a, c, e_per_atom


def relax_slab(name, lib, pot, bulk, miller, layers):
    slab = ase_surface(bulk, miller, layers=layers, vacuum=VACUUM, periodic=True)
    slab.center(vacuum=VACUUM, axis=2)
    slab.calc = meam_calc(name + '_slab_' + ''.join(str(x) for x in miller), lib, pot)
    LBFGS(slab, logfile=None).run(fmax=0.01, steps=300)
    E = float(slab.get_potential_energy())
    cell = slab.cell.array
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    return slab, E, area


def main():
    all_results = {}
    for name, (lib, pot) in POTENTIALS.items():
        print(f'\n##### {name} #####', flush=True)
        bulk, a, c, e_per_atom = relax_bulk(name, lib, pot)
        print(f'  bulk: a={a:.4f} c={c:.4f} c/a={c/a:.3f} E/atom={e_per_atom:.4f} eV', flush=True)
        faces_out = []
        for tag, L in FACES:
            miller = tuple(int(x) for x in tag)
            slab, E, area = relax_slab(name, lib, pot, bulk, miller, L)
            n = len(slab)
            gamma = (E - n * e_per_atom) / (2 * area) * EV_A2_TO_MJ_M2
            print(f'  ({tag}) L={L}  n={n}  A={area:.3f}  γ={gamma:.1f} mJ/m²', flush=True)
            faces_out.append({'face': tag, 'L': L, 'n': n, 'area_A2': area,
                              'E_eV': E, 'gamma_mJ_m2': gamma})
        all_results[name] = {'bulk': {'a': a, 'c': c, 'c_over_a': c/a, 'E_per_atom_eV': e_per_atom},
                             'slabs': faces_out}

    out = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/meam_surface_energies.json')
    with open(out, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f'\nsaved -> {out}', flush=True)

    # summary table
    print('\n' + '=' * 70)
    print('γ (mJ/m²) comparison')
    print(f"{'method':>12} {'(001)':>7} {'(110)':>7} {'(100)':>7} {'(101)':>7}")
    for name, r in all_results.items():
        by_face = {s['face']: s['gamma_mJ_m2'] for s in r['slabs']}
        print(f"{name:>12} {by_face['001']:>7.1f} {by_face['110']:>7.1f} "
              f"{by_face['100']:>7.1f} {by_face['101']:>7.1f}")


if __name__ == '__main__':
    main()
