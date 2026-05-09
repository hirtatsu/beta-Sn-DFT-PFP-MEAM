"""
elastic_betasn.py

β-Sn (I4_1/amd, 4-atom cell) elastic constant tensor via PFP v8.
Tatsumi et al. (MSMSE 2026, in review) Table 2.

Method (per manuscript):
  - Experimental lattice (ICSD 40037): a=5.831 Å, c=3.182 Å
  - Full relaxation with ExpCellFilter, fmax=0.001 eV/Å
  - 6 Voigt strain modes at ±0.5%, atoms-only relaxation at fixed cell (fmax=0.005 eV/Å)
  - Central difference: C_ij = (σ_i(+ε_j) - σ_i(-ε_j)) / (2 Δε)
  - Symmetrize to tetragonal 4/mmm

Usage:
  python elastic_betasn.py PBE
  python elastic_betasn.py R2SCAN
"""

import sys
import json
import os
import numpy as np
from ase import Atoms
from ase.optimize import LBFGS
from ase.filters import ExpCellFilter

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

EV_PER_A3_TO_GPA = 160.2176634
VOIGT_LABELS = ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']
DELTA = 0.005  # ±0.5% strain
FMAX_CELL = 0.001
FMAX_ATOMS = 0.005

# Reference values (manuscript Table 2 and Table 1)
PAPER_LATTICE = {
    'PBE':    {'a': 5.929, 'c': 3.201},
    'R2SCAN': {'a': 5.882, 'c': 3.192},
}
PAPER_CIJ = {
    'PBE':    {'C11': 114.0, 'C12': 41.5, 'C13': 41.2, 'C33': 104.6, 'C44': 29.5, 'C66': 31.9, 'B': 64.5},
    'R2SCAN': {'C11': 118.0, 'C12': 15.2, 'C13': 40.7, 'C33': 129.7, 'C44': 29.7, 'C66': 24.8, 'B': 62.1},
}


def build_beta_sn():
    """β-Sn (I4_1/amd, Z=4) with experimental ICSD 40037 parameters."""
    a, c = 5.831, 3.182
    cell = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    # Wyckoff 4a of I4_1/amd (origin choice 2) in conventional cell
    scaled = [[0.0, 0.0, 0.0],
              [0.5, 0.5, 0.5],
              [0.0, 0.5, 0.25],
              [0.5, 0.0, 0.75]]
    return Atoms('Sn4', scaled_positions=scaled, cell=cell, pbc=True)


def apply_voigt_strain(atoms, eps):
    """Apply symmetric Voigt strain [e1,e2,e3,e4,e5,e6] (engineering shear)."""
    F = np.array([
        [1.0 + eps[0], 0.5 * eps[5], 0.5 * eps[4]],
        [0.5 * eps[5], 1.0 + eps[1], 0.5 * eps[3]],
        [0.5 * eps[4], 0.5 * eps[3], 1.0 + eps[2]],
    ])
    new = atoms.copy()
    new.set_cell(atoms.cell @ F.T, scale_atoms=True)
    return new


def stress_GPa(atoms):
    return atoms.get_stress(voigt=True) * EV_PER_A3_TO_GPA


def main(mode):
    print(f"=== β-Sn elastic tensor — independent recalc  [{mode}] ===")
    estimator = Estimator(calc_mode=mode)
    calc = ASECalculator(estimator)

    # 1. Equilibrium cell
    atoms = build_beta_sn()
    atoms.calc = calc
    a0 = atoms.cell.lengths()[0]
    c0 = atoms.cell.lengths()[2]
    print(f"Experimental start: a={a0:.4f} Å, c={c0:.4f} Å")

    print("Relaxing cell + atoms (ExpCellFilter, fmax=0.001)...")
    ecf = ExpCellFilter(atoms, hydrostatic_strain=False)
    LBFGS(ecf, logfile=f'opt_{mode.lower()}.log').run(fmax=FMAX_CELL)
    a_eq = atoms.cell.lengths()[0]
    c_eq = atoms.cell.lengths()[2]
    print(f"Equilibrium:        a={a_eq:.4f} Å, c={c_eq:.4f} Å, c/a={c_eq/a_eq:.4f}")
    print(f"Paper {mode}:        a={PAPER_LATTICE[mode]['a']:.4f} Å, c={PAPER_LATTICE[mode]['c']:.4f} Å")

    s0 = stress_GPa(atoms)
    print(f"Residual stress (GPa): {np.round(s0, 3)}")

    # 2. Finite-strain elastic tensor
    print(f"\nFinite-strain sweep, Δε = ±{DELTA*100:.1f}% (atoms-only relax per strain)")
    C = np.zeros((6, 6))
    for j in range(6):
        for sign, tag in [(+1, '+'), (-1, '-')]:
            eps = np.zeros(6)
            eps[j] = sign * DELTA
            a_strained = apply_voigt_strain(atoms, eps)
            a_strained.calc = calc
            LBFGS(a_strained, logfile=None).run(fmax=FMAX_ATOMS, steps=100)
            s = stress_GPa(a_strained)
            if sign == +1:
                s_plus = s
            else:
                s_minus = s
        C[:, j] = (s_plus - s_minus) / (2.0 * DELTA)
        print(f"  ε_{VOIGT_LABELS[j]}: dσ/dε → C_{{*,{j+1}}} = "
              + ", ".join(f"{C[i,j]:+7.2f}" for i in range(6)))

    # 3. Tetragonal (4/mmm) symmetrization
    C11 = 0.5 * (C[0, 0] + C[1, 1])
    C12 = 0.5 * (C[0, 1] + C[1, 0])
    C13 = 0.25 * (C[0, 2] + C[1, 2] + C[2, 0] + C[2, 1])
    C33 = C[2, 2]
    C44 = 0.5 * (C[3, 3] + C[4, 4])
    C66 = C[5, 5]
    B_voigt = (2 * C11 + C33 + 2 * C12 + 4 * C13) / 9.0
    mine = {'C11': C11, 'C12': C12, 'C13': C13, 'C33': C33, 'C44': C44, 'C66': C66, 'B': B_voigt}

    # 4. Comparison with manuscript Table 2
    paper = PAPER_CIJ[mode]
    print("\n" + "=" * 72)
    print(f"{'component':>10}  {'this run':>10}  {'paper':>10}  {'Δ (GPa)':>10}  {'Δ (%)':>10}")
    print("-" * 72)
    for k in ['C11', 'C12', 'C13', 'C33', 'C44', 'C66', 'B']:
        diff = mine[k] - paper[k]
        pct = 100 * diff / paper[k] if paper[k] != 0 else float('nan')
        print(f"{k:>10}  {mine[k]:>10.2f}  {paper[k]:>10.2f}  {diff:>+10.2f}  {pct:>+9.2f}%")
    print("=" * 72)

    os.makedirs('results', exist_ok=True)
    out = {
        'mode': mode,
        'lattice': {'a': a_eq, 'c': c_eq, 'c_over_a': c_eq/a_eq},
        'elastic_constants_GPa': {k: float(v) for k, v in mine.items()},
        'paper_Table2_GPa': paper,
        'full_Cij_GPa': C.tolist(),
        'strain_delta': DELTA,
        'fmax_cell': FMAX_CELL,
        'fmax_atoms': FMAX_ATOMS,
    }
    out_path = f'results/elastic_{mode.lower()}.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {out_path}")
    return mine


if __name__ == '__main__':
    mode = sys.argv[1].upper() if len(sys.argv) > 1 else 'PBE'
    if mode not in ('PBE', 'R2SCAN'):
        sys.exit(f"Usage: python {sys.argv[0]} [PBE|R2SCAN]")
    main(mode)
