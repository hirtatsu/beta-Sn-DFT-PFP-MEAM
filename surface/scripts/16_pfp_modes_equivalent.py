"""
16_pfp_modes_equivalent.py (run on Matlantis)

Compute β-Sn surface energies for 4 PFP modes under conditions matching DFT BFGS:
  - Bulk re-optimized in EACH mode (ExpCellFilter + LBFGS, fmax 0.015)
  - Slabs built from that mode's bulk (mode-self-consistent)
  - ALL slab atoms free (no middle freezing)
  - LBFGS, fmax = 0.015 eV/Å  (≈ DFT 3e-4 Ha/Bohr)
  - Layers: 8 for (001),(110); 12 for (100),(101); vacuum 15 Å
"""
import json
import os
import time
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.io import write
from ase.optimize import LBFGS
from ase.filters import ExpCellFilter
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

EV_A2_TO_MJ_M2 = 16021.766
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUTDIR = os.path.join(ROOT, 'results_modes')
LOGDIR = os.path.join(ROOT, 'logs_modes')
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

MODES = ['PBE', 'PBE_PLUS_D3', 'R2SCAN', 'R2SCAN_PLUS_D3']
FACES = [('001', 8), ('110', 8), ('100', 12), ('101', 12)]
VACUUM = 15.0
FMAX = 0.015  # eV/Å (DFT-equivalent)
EXP_LAT = {'a': 5.831, 'c': 3.182}  # experimental start


def build_bulk(a, c):
    cell = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    scaled = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
              [0.0, 0.5, 0.25], [0.5, 0.0, 0.75]]
    return Atoms('Sn4', scaled_positions=scaled, cell=cell, pbc=True)


def relax_bulk(calc, logfile):
    atoms = build_bulk(EXP_LAT['a'], EXP_LAT['c'])
    atoms.calc = calc
    ecf = ExpCellFilter(atoms, hydrostatic_strain=False)
    opt = LBFGS(ecf, logfile=logfile)
    opt.run(fmax=FMAX, steps=400)
    return atoms


def relax_slab(atoms, calc, logfile):
    atoms.calc = calc
    opt = LBFGS(atoms, logfile=logfile)
    opt.run(fmax=FMAX, steps=300)


def main():
    est0 = Estimator(calc_mode=MODES[0])  # warmup
    _ = ASECalculator(est0)

    all_results = {}
    t0 = time.time()

    for mode in MODES:
        print(f'\n########## MODE: {mode} ##########', flush=True)
        est = Estimator(calc_mode=mode)
        calc = ASECalculator(est)

        mode_res = {'bulk': {}, 'slabs': []}

        # Bulk
        bulk_log = os.path.join(LOGDIR, f'bulk_{mode}.log')
        bulk = relax_bulk(calc, bulk_log)
        a = float(bulk.cell.lengths()[0])
        c = float(bulk.cell.lengths()[2])
        e_bulk = float(bulk.get_potential_energy())
        e_bulk_per_atom = e_bulk / 4.0
        mode_res['bulk'] = {
            'a': a, 'c': c, 'c_over_a': c/a,
            'E_per_atom_eV': e_bulk_per_atom,
            'da_pct': 100*(a - EXP_LAT['a'])/EXP_LAT['a'],
            'dc_pct': 100*(c - EXP_LAT['c'])/EXP_LAT['c'],
        }
        print(f'  bulk: a={a:.4f}, c={c:.4f} (Δa={mode_res["bulk"]["da_pct"]:+.2f}%, '
              f'Δc={mode_res["bulk"]["dc_pct"]:+.2f}%), E/atom={e_bulk_per_atom:.4f} eV', flush=True)
        write(os.path.join(OUTDIR, f'bulk_{mode}.traj'), bulk)

        # Slabs
        for face, L in FACES:
            miller = tuple(int(x) for x in face)
            slab = ase_surface(bulk, miller, layers=L, vacuum=VACUUM, periodic=True)
            slab.center(vacuum=VACUUM, axis=2)
            n = len(slab)
            cell = slab.cell.array
            area = float(np.linalg.norm(np.cross(cell[0], cell[1])))

            slab_log = os.path.join(LOGDIR, f'slab_{face}_L{L}_{mode}.log')
            relax_slab(slab, calc, slab_log)
            E = float(slab.get_potential_energy())
            gamma = (E - n * e_bulk_per_atom) / (2 * area) * EV_A2_TO_MJ_M2

            mode_res['slabs'].append({
                'face': face, 'L': L, 'n_atoms': n,
                'area_A2': area, 'E_slab_eV': E, 'gamma_mJ_m2': gamma,
            })
            print(f'  ({face}) L={L} n={n} γ={gamma:.1f} mJ/m²   '
                  f'(elapsed {time.time()-t0:.0f}s)', flush=True)
            write(os.path.join(OUTDIR, f'slab_{face}_L{L}_{mode}.traj'), slab)

        all_results[mode] = mode_res
        with open(os.path.join(OUTDIR, 'modes_surface_energies.json'), 'w') as f:
            json.dump(all_results, f, indent=2)

    # summary
    print('\n' + '=' * 72, flush=True)
    print('Surface energy γ (mJ/m²) at mode-self-consistent conditions', flush=True)
    print(f"{'mode':>16} {'(001)':>8} {'(110)':>8} {'(100)':>8} {'(101)':>8} {'a':>7} {'c':>7}", flush=True)
    for mode, r in all_results.items():
        by_face = {s['face']: s['gamma_mJ_m2'] for s in r['slabs']}
        print(f"{mode:>16} "
              f"{by_face.get('001',0):>8.1f} {by_face.get('110',0):>8.1f} "
              f"{by_face.get('100',0):>8.1f} {by_face.get('101',0):>8.1f} "
              f"{r['bulk']['a']:>7.4f} {r['bulk']['c']:>7.4f}", flush=True)
    print(f'wall: {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
