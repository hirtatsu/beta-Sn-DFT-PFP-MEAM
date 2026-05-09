"""
elastic_all_modes.py
β-Sn elastic tensor for all crystal-relevant PFP v8 modes.
Ranks modes by MAPE vs experiment (Rayne & Chandrasekhar 1960).
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
DELTA = 0.005
MODES = ['CRYSTAL', 'CRYSTAL_PLUS_D3', 'CRYSTAL_U0', 'CRYSTAL_U0_PLUS_D3',
         'PBE', 'PBE_PLUS_D3', 'PBE_U', 'PBE_U_PLUS_D3',
         'R2SCAN', 'R2SCAN_PLUS_D3']

EXP = {'C11': 73.4, 'C12': 59.1, 'C13': 35.8, 'C33': 90.7, 'C44': 22.0, 'C66': 24.0}
EXP_LATTICE = {'a': 5.831, 'c': 3.182}


def build_beta_sn():
    a, c = EXP_LATTICE['a'], EXP_LATTICE['c']
    cell = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    scaled = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
              [0.0, 0.5, 0.25], [0.5, 0.0, 0.75]]
    return Atoms('Sn4', scaled_positions=scaled, cell=cell, pbc=True)


def apply_strain(atoms, eps):
    F = np.array([
        [1 + eps[0], eps[5] / 2, eps[4] / 2],
        [eps[5] / 2, 1 + eps[1], eps[3] / 2],
        [eps[4] / 2, eps[3] / 2, 1 + eps[2]],
    ])
    new = atoms.copy()
    new.set_cell(atoms.cell @ F.T, scale_atoms=True)
    return new


def stress_GPa(atoms):
    return atoms.get_stress(voigt=True) * EV_PER_A3_TO_GPA


def run_mode(mode):
    est = Estimator(calc_mode=mode)
    calc = ASECalculator(est)

    atoms = build_beta_sn()
    atoms.calc = calc
    try:
        ecf = ExpCellFilter(atoms, hydrostatic_strain=False)
        LBFGS(ecf, logfile=None).run(fmax=0.001, steps=300)
    except Exception as e:
        return {'error': f'opt failed: {e}'}

    a_eq = atoms.cell.lengths()[0]
    c_eq = atoms.cell.lengths()[2]

    C = np.zeros((6, 6))
    for j in range(6):
        s = {}
        for sign in (+1, -1):
            eps = np.zeros(6)
            eps[j] = sign * DELTA
            a2 = apply_strain(atoms, eps)
            a2.calc = calc
            LBFGS(a2, logfile=None).run(fmax=0.005, steps=100)
            s[sign] = stress_GPa(a2)
        C[:, j] = (s[+1] - s[-1]) / (2 * DELTA)

    cij = {
        'C11': 0.5 * (C[0, 0] + C[1, 1]),
        'C12': 0.5 * (C[0, 1] + C[1, 0]),
        'C13': 0.25 * (C[0, 2] + C[1, 2] + C[2, 0] + C[2, 1]),
        'C33': C[2, 2],
        'C44': 0.5 * (C[3, 3] + C[4, 4]),
        'C66': C[5, 5],
    }

    # Errors vs experiment
    diffs = {k: cij[k] - EXP[k] for k in EXP}
    pct = {k: 100 * diffs[k] / EXP[k] for k in EXP}
    mape = float(np.mean([abs(v) for v in pct.values()]))
    rmse = float(np.sqrt(np.mean([(cij[k] - EXP[k])**2 for k in EXP])))
    max_abs_pct = float(max(abs(v) for v in pct.values()))

    return {
        'mode': mode,
        'lattice': {'a': a_eq, 'c': c_eq,
                    'da_pct': 100 * (a_eq - EXP_LATTICE['a']) / EXP_LATTICE['a'],
                    'dc_pct': 100 * (c_eq - EXP_LATTICE['c']) / EXP_LATTICE['c']},
        'cij': cij,
        'diff_pct': pct,
        'MAPE': mape,
        'RMSE_GPa': rmse,
        'max_abs_pct': max_abs_pct,
        'full_C': C.tolist(),
    }


def main():
    os.makedirs('results', exist_ok=True)
    results = {}
    for mode in MODES:
        print(f"\n{'='*60}\nRunning {mode}\n{'='*60}", flush=True)
        try:
            r = run_mode(mode)
        except Exception as e:
            print(f"  !! failed: {e}", flush=True)
            r = {'error': str(e)}
        results[mode] = r
        if 'error' in r:
            print(f"  ERROR: {r['error']}")
            continue
        print(f"  a={r['lattice']['a']:.4f} c={r['lattice']['c']:.4f}  "
              f"(Δa={r['lattice']['da_pct']:+.2f}%, Δc={r['lattice']['dc_pct']:+.2f}%)")
        print(f"  Cij (GPa): " + "  ".join(f"{k}={v:.1f}" for k, v in r['cij'].items()))
        print(f"  MAPE={r['MAPE']:.1f}%  RMSE={r['RMSE_GPa']:.2f} GPa  max|%err|={r['max_abs_pct']:.1f}%")

        with open('results/all_modes.json', 'w') as f:
            json.dump(results, f, indent=2)

    # Rank
    ok = [(m, r) for m, r in results.items() if 'error' not in r]
    ok.sort(key=lambda x: x[1]['MAPE'])
    print("\n\n" + "#" * 74)
    print("  RANKING vs Rayne & Chandrasekhar 1960  (lower MAPE = better)")
    print("#" * 74)
    header = f"{'rank':>4} {'mode':>22} {'MAPE %':>8} {'RMSE':>7} {'max%':>6}   " \
             + "  ".join(f"{k:>5}" for k in ['C11', 'C12', 'C13', 'C33', 'C44', 'C66'])
    print(header)
    print("-" * 74)
    print(f"{'-':>4} {'EXPERIMENT':>22} {'-':>8} {'-':>7} {'-':>6}   "
          + "  ".join(f"{EXP[k]:>5.1f}" for k in ['C11', 'C12', 'C13', 'C33', 'C44', 'C66']))
    for i, (m, r) in enumerate(ok, 1):
        row = f"{i:>4} {m:>22} {r['MAPE']:>8.1f} {r['RMSE_GPa']:>7.2f} {r['max_abs_pct']:>6.1f}   " \
              + "  ".join(f"{r['cij'][k]:>5.1f}" for k in ['C11', 'C12', 'C13', 'C33', 'C44', 'C66'])
        print(row)


if __name__ == '__main__':
    main()
