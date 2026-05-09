"""
02_surface_convergence.py
Slab thickness & vacuum convergence for β-Sn surfaces on PFP/PBE+D3.

Phase A: 4 faces × thickness {8,12,16,20} @ vacuum 15 A
Phase B: (001) face × vacuum {10,15,20} @ thickness = converged from A

Surface energy: gamma = (E_slab - n * e_bulk_per_atom) / (2 * A)
"""
import json
import os
import time
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.io import read, write
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESDIR = os.path.join(ROOT, 'results')
LOGDIR = os.path.join(ROOT, 'logs')
os.makedirs(RESDIR, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

EV_A2_TO_MJ_M2 = 16021.766  # eV/A^2 -> mJ/m^2
FACES = [(1, 0, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1)]
THICKNESS_LIST = [8, 12, 16, 20]
VACUUM_FIXED = 15.0


def load_bulk():
    with open(os.path.join(RESDIR, 'bulk_pfp_pbe_d3.json')) as f:
        info = json.load(f)
    bulk = read(os.path.join(RESDIR, 'bulk_pfp_pbe_d3.traj'))
    return bulk, info['E_per_atom_eV']


def build_slab(bulk, miller, layers, vacuum):
    slab = ase_surface(bulk, miller, layers=layers, vacuum=vacuum, periodic=True)
    slab.center(vacuum=vacuum, axis=2)
    return slab


def relax_slab(slab, calc, logfile):
    slab.calc = calc
    z = slab.positions[:, 2]
    z_mid = 0.5 * (z.min() + z.max())
    thickness = z.max() - z.min()
    freeze_mask = np.abs(z - z_mid) < 0.15 * thickness
    slab.set_constraint(FixAtoms(mask=freeze_mask))
    opt = LBFGS(slab, logfile=logfile)
    opt.run(fmax=0.02, steps=200)
    slab.set_constraint()
    return slab


def surface_energy(slab, e_bulk_per_atom):
    e_slab = slab.get_potential_energy()
    n = len(slab)
    cell = slab.cell.array
    area = np.linalg.norm(np.cross(cell[0], cell[1]))
    gamma = (e_slab - n * e_bulk_per_atom) / (2.0 * area)
    return gamma * EV_A2_TO_MJ_M2, e_slab, n, area


def run_all():
    bulk, e_bulk = load_bulk()
    est = Estimator(calc_mode='PBE_PLUS_D3')
    calc = ASECalculator(est)

    all_results = []
    t0 = time.time()

    # Phase A
    for miller in FACES:
        tag = ''.join(str(i) for i in miller)
        for layers in THICKNESS_LIST:
            label = f"{tag}_L{layers}_V{int(VACUUM_FIXED)}"
            print(f"\n[A] {label}  (elapsed {time.time()-t0:.0f}s)", flush=True)
            slab = build_slab(bulk, miller, layers, VACUUM_FIXED)
            n0 = len(slab)
            log = os.path.join(LOGDIR, f'surf_{label}.log')
            try:
                slab = relax_slab(slab, calc, log)
                gamma, e_slab, n, area = surface_energy(slab, e_bulk)
                print(f"  n_atoms={n}  A={area:.3f} A^2  gamma={gamma:.1f} mJ/m^2", flush=True)
                all_results.append({
                    'phase': 'A',
                    'miller': miller,
                    'tag': tag,
                    'layers': layers,
                    'vacuum': VACUUM_FIXED,
                    'n_atoms_built': n0,
                    'n_atoms': n,
                    'area_A2': area,
                    'E_slab_eV': e_slab,
                    'gamma_mJ_m2': gamma,
                })
                write(os.path.join(RESDIR, f'slab_{label}.traj'), slab)
            except Exception as e:
                print(f"  FAILED: {e}", flush=True)
                all_results.append({
                    'phase': 'A', 'miller': miller, 'tag': tag,
                    'layers': layers, 'vacuum': VACUUM_FIXED,
                    'error': str(e),
                })
            # save incrementally
            with open(os.path.join(RESDIR, 'surface_convergence.json'), 'w') as f:
                json.dump(all_results, f, indent=2)

    # Determine converged thickness per face (first layer where |gamma-gamma_max_thickness| < 5 mJ/m^2)
    converged = {}
    for miller in FACES:
        tag = ''.join(str(i) for i in miller)
        rows = [r for r in all_results if r.get('tag') == tag and 'gamma_mJ_m2' in r]
        rows.sort(key=lambda r: r['layers'])
        if not rows:
            continue
        ref = rows[-1]['gamma_mJ_m2']
        for r in rows:
            if abs(r['gamma_mJ_m2'] - ref) < 5.0:
                converged[tag] = r['layers']
                break
        else:
            converged[tag] = rows[-1]['layers']
    print("\nConverged thickness per face:", converged, flush=True)

    # Phase B: vacuum convergence on (001) with its converged thickness
    n001 = converged.get('001', 16)
    for vac in [10.0, 15.0, 20.0]:
        label = f"001_L{n001}_V{int(vac)}"
        print(f"\n[B] {label}", flush=True)
        if abs(vac - VACUUM_FIXED) < 1e-6:
            # already computed in phase A
            prev = [r for r in all_results if r.get('tag') == '001'
                    and r.get('layers') == n001 and abs(r.get('vacuum', 0) - vac) < 1e-6
                    and 'gamma_mJ_m2' in r]
            if prev:
                entry = dict(prev[0])
                entry['phase'] = 'B'
                all_results.append(entry)
                print(f"  reused  gamma={entry['gamma_mJ_m2']:.1f}", flush=True)
                continue
        slab = build_slab(bulk, (0, 0, 1), n001, vac)
        log = os.path.join(LOGDIR, f'surf_{label}.log')
        try:
            slab = relax_slab(slab, calc, log)
            gamma, e_slab, n, area = surface_energy(slab, e_bulk)
            print(f"  gamma={gamma:.1f} mJ/m^2", flush=True)
            all_results.append({
                'phase': 'B',
                'miller': (0, 0, 1), 'tag': '001',
                'layers': n001, 'vacuum': vac,
                'n_atoms': n, 'area_A2': area,
                'E_slab_eV': e_slab, 'gamma_mJ_m2': gamma,
            })
        except Exception as e:
            print(f"  FAILED: {e}", flush=True)
            all_results.append({
                'phase': 'B', 'miller': (0, 0, 1), 'tag': '001',
                'layers': n001, 'vacuum': vac, 'error': str(e),
            })
        with open(os.path.join(RESDIR, 'surface_convergence.json'), 'w') as f:
            json.dump(all_results, f, indent=2)

    # Summary
    print("\n" + "=" * 70)
    print(f"{'face':>6} {'layers':>6} {'vac':>5} {'nat':>4} {'gamma (mJ/m^2)':>16}")
    print("-" * 70)
    for r in all_results:
        if 'gamma_mJ_m2' not in r:
            continue
        print(f"{r['tag']:>6} {r['layers']:>6} {r['vacuum']:>5.0f} "
              f"{r.get('n_atoms',0):>4} {r['gamma_mJ_m2']:>16.1f}")
    print(f"\nTotal wall time: {time.time()-t0:.0f}s")


if __name__ == '__main__':
    run_all()
