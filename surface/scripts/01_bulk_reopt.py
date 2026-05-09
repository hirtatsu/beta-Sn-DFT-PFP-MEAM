"""
01_bulk_reopt.py
β-Sn bulk re-optimization with PFP/PBE+D3.
Writes optimized cell & atoms trajectory for downstream slab/GB construction.
"""
import json
import os
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import LBFGS
from ase.filters import ExpCellFilter
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

OUTDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
os.makedirs(OUTDIR, exist_ok=True)

EXP = {'a': 5.831, 'c': 3.182}


def build_beta_sn():
    a, c = EXP['a'], EXP['c']
    cell = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    scaled = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
              [0.0, 0.5, 0.25], [0.5, 0.0, 0.75]]
    return Atoms('Sn4', scaled_positions=scaled, cell=cell, pbc=True)


def main():
    est = Estimator(calc_mode='PBE_PLUS_D3')
    calc = ASECalculator(est)

    atoms = build_beta_sn()
    atoms.calc = calc

    logfile = os.path.join(OUTDIR, 'bulk_reopt.log')
    ecf = ExpCellFilter(atoms, hydrostatic_strain=False)
    LBFGS(ecf, logfile=logfile).run(fmax=0.001, steps=400)

    a_eq = float(atoms.cell.lengths()[0])
    b_eq = float(atoms.cell.lengths()[1])
    c_eq = float(atoms.cell.lengths()[2])
    vol = float(atoms.get_volume())
    e_tot = float(atoms.get_potential_energy())

    result = {
        'mode': 'PBE_PLUS_D3',
        'a': a_eq,
        'b': b_eq,
        'c': c_eq,
        'c_over_a': c_eq / a_eq,
        'volume_A3': vol,
        'E_total_eV': e_tot,
        'E_per_atom_eV': e_tot / len(atoms),
        'da_pct_vs_exp': 100 * (a_eq - EXP['a']) / EXP['a'],
        'dc_pct_vs_exp': 100 * (c_eq - EXP['c']) / EXP['c'],
        'n_atoms': len(atoms),
        'cell': atoms.cell.array.tolist(),
        'scaled_positions': atoms.get_scaled_positions().tolist(),
    }

    with open(os.path.join(OUTDIR, 'bulk_pfp_pbe_d3.json'), 'w') as f:
        json.dump(result, f, indent=2)

    write(os.path.join(OUTDIR, 'bulk_pfp_pbe_d3.traj'), atoms)
    write(os.path.join(OUTDIR, 'bulk_pfp_pbe_d3.cif'), atoms)

    print(f"a = {a_eq:.4f} A  (Δ = {result['da_pct_vs_exp']:+.2f}%)")
    print(f"c = {c_eq:.4f} A  (Δ = {result['dc_pct_vs_exp']:+.2f}%)")
    print(f"c/a = {result['c_over_a']:.4f}  (exp = {EXP['c']/EXP['a']:.4f})")
    print(f"E_per_atom = {result['E_per_atom_eV']:.4f} eV")
    print(f"saved -> {OUTDIR}")


if __name__ == '__main__':
    main()
