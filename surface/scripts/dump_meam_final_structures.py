"""Re-run MEAM relaxation (fast, ~30s) and save final atoms as .xyz + .cif for the record."""
import os, json
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.calculators.lammpslib import LAMMPSlib
from ase.filters import ExpCellFilter
from ase.io import write
from ase.optimize import LBFGS

FD = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/final_deliverables/data/meam_lammps')
POT = os.path.join(FD, 'meam_potentials')
OUT = os.path.join(FD, 'structures')
os.makedirs(OUT, exist_ok=True)

POTENTIALS = {
    'Ravelo1997':  ('library.Sn.Ravelo1997.meam',  'Sn.Ravelo1997.meam'),
    'Etesami2018': ('library.Sn.Etesami2018.meam', 'Sn.Etesami2018.meam'),
}
EXP = {'a': 5.831, 'c': 3.182}
FACES = [('001', 8), ('110', 8), ('100', 12), ('101', 12), ('111', 12)]
VACUUM = 15.0


def build_bulk(a, c):
    return Atoms('Sn4',
                 scaled_positions=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]],
                 cell=[[a,0,0],[0,a,0],[0,0,c]], pbc=True)


def calc_for(lib, pot, tag):
    return LAMMPSlib(lmpcmds=[
        'pair_style meam',
        f'pair_coeff * * {os.path.join(POT, lib)} Sn {os.path.join(POT, pot)} Sn',
    ], atom_types={'Sn': 1},
       log_file=f'/tmp/lammps_dump_{tag}.log', keep_alive=True)


for name, (lib, pot) in POTENTIALS.items():
    print(f'\n=== {name} ===')
    # Bulk relax
    bulk = build_bulk(EXP['a'], EXP['c'])
    bulk.calc = calc_for(lib, pot, f'{name}_bulk')
    LBFGS(ExpCellFilter(bulk, hydrostatic_strain=False), logfile=None).run(fmax=0.001, steps=300)
    write(os.path.join(OUT, f'bulk_{name}.cif'), bulk)
    write(os.path.join(OUT, f'bulk_{name}.xyz'), bulk)
    print(f'  bulk: a={bulk.cell.lengths()[0]:.4f}, c={bulk.cell.lengths()[2]:.4f}, E/atom={bulk.get_potential_energy()/len(bulk):.4f} eV')

    # Slabs
    for tag, L in FACES:
        miller = tuple(int(x) for x in tag)
        slab = ase_surface(bulk, miller, layers=L, vacuum=VACUUM, periodic=True)
        slab.center(vacuum=VACUUM, axis=2)
        slab.calc = calc_for(lib, pot, f'{name}_{tag}')
        LBFGS(slab, logfile=None).run(fmax=0.01, steps=300)
        stem = f'slab_{tag}_L{L}_{name}'
        write(os.path.join(OUT, f'{stem}.cif'), slab)
        write(os.path.join(OUT, f'{stem}.xyz'), slab)
        print(f'  ({tag}) L={L}: n={len(slab)}, E={slab.get_potential_energy():.4f} eV')

print(f'\nsaved structures to {OUT}')
