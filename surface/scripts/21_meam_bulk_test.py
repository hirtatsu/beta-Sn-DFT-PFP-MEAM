"""Sanity-test MEAM Ravelo-Baskes (reconstructed) and Etesami on β-Sn bulk.
Uses LAMMPS through ASE to relax bulk β-Sn and report a, c, E/atom, c/a.
Compare with Vella 2017 Table III reference values for Ravelo-Baskes:
  β-Sn: a=5.86, c=3.09, Ec=3.08 (expected)
"""
import os
import numpy as np
from ase import Atoms
from ase.calculators.lammpslib import LAMMPSlib
from ase.optimize import LBFGS
from ase.filters import ExpCellFilter

POT_DIR = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/meam_potentials')

POTENTIALS = {
    'Ravelo1997': ('library.Sn.Ravelo1997.meam', 'Sn.Ravelo1997.meam'),
    'Etesami2018': ('library.Sn.Etesami2018.meam', 'Sn.Etesami2018.meam'),
}

EXP_LAT = {'a': 5.831, 'c': 3.182}

def build_beta_sn():
    a, c = EXP_LAT['a'], EXP_LAT['c']
    return Atoms('Sn4',
                 scaled_positions=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]],
                 cell=[[a,0,0],[0,a,0],[0,0,c]], pbc=True)


def run_relax(name, lib, pot):
    atoms = build_beta_sn()
    cmds = [
        f'pair_style meam',
        f'pair_coeff * * {os.path.join(POT_DIR, lib)} Sn {os.path.join(POT_DIR, pot)} Sn',
    ]
    calc = LAMMPSlib(lmpcmds=cmds, atom_types={'Sn': 1}, log_file=f'/tmp/lammps_{name}.log',
                     keep_alive=True)
    atoms.calc = calc
    E_init = atoms.get_potential_energy()
    # cell relax
    ecf = ExpCellFilter(atoms, hydrostatic_strain=False)
    opt = LBFGS(ecf, logfile=None)
    opt.run(fmax=0.001, steps=300)
    a_eq = float(atoms.cell.lengths()[0])
    c_eq = float(atoms.cell.lengths()[2])
    E = float(atoms.get_potential_energy())
    return a_eq, c_eq, E / len(atoms), E


print(f"{'potential':>15} {'a':>7} {'c':>7} {'c/a':>6} {'E/atom':>9}  {'Δa%':>6} {'Δc%':>6}")
print('-' * 68)
print(f"{'experiment':>15} {5.831:>7.4f} {3.182:>7.4f} {3.182/5.831:>6.3f} {'—':>9}")
for name, (lib, pot) in POTENTIALS.items():
    try:
        a, c, eper, etot = run_relax(name, lib, pot)
        print(f"{name:>15} {a:>7.4f} {c:>7.4f} {c/a:>6.3f} {eper:>9.4f} "
              f"{100*(a-EXP_LAT['a'])/EXP_LAT['a']:>6.2f} "
              f"{100*(c-EXP_LAT['c'])/EXP_LAT['c']:>6.2f}")
    except Exception as e:
        print(f"{name:>15} FAILED: {e}")
