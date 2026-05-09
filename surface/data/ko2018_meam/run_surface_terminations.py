"""Explore alternative surface terminations for beta-Sn (001) and (111) with Ko 2018 MEAM.

Strategy: pre-shift the bulk along the surface-normal direction before slab cleavage.
beta-Sn has 4 atoms at z-fractions {0, 0.25, 0.5, 0.75}; 4 candidate cuts may collapse
to 2 inequivalent terminations by inversion symmetry.
"""
import os
import json
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import LBFGS

import os as _os
PROJECT_ROOT = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), "..", "..", ".."))
POT_DIR = _os.path.join(PROJECT_ROOT, "meam_potentials")
LAMMPS_BIN = _os.environ.get("LAMMPS_BIN", "lmp_kokkos")
EV_A2_TO_MJ_M2 = 16021.766
VACUUM = 18.0

A0, C0 = 5.8589, 3.2059
EPER_REF = -3.10195  # measured in our previous bulk run


def build_beta_sn(a, c, z_shift=0.0):
    """Conventional 4-atom tetragonal cell, optionally pre-shifted in z (fractional)."""
    sp = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.0, 0.5, 0.25],
        [0.5, 0.0, 0.75],
    ]
    sp = [[x, y, (z + z_shift) % 1.0] for x, y, z in sp]
    return Atoms(
        'Sn4',
        scaled_positions=sp,
        cell=[[a, 0, 0], [0, a, 0], [0, 0, c]],
        pbc=True,
    )


def make_calc(label):
    return LAMMPS(
        command=LAMMPS_BIN,
        pair_style='meam',
        pair_coeff=[
            f'* * {POT_DIR}/library.Sn.Ko2018.meam Sn {POT_DIR}/Sn.Ko2018.meam Sn'
        ],
        masses=['1 118.71'],
        files=[],
        keep_alive=False,
        tmp_dir=f'/tmp/lammps_{label}',
        keep_tmp_files=False,
        always_triclinic=True,
    )


def relax_slab(miller, layers, z_shift, label):
    bulk = build_beta_sn(A0, C0, z_shift=z_shift)
    slab = ase_surface(bulk, miller, layers=layers, vacuum=VACUUM, periodic=True)
    slab.center(vacuum=VACUUM, axis=2)
    slab.calc = make_calc(label)
    LBFGS(slab, logfile=None).run(fmax=0.005, steps=300)
    e = float(slab.get_potential_energy())
    cell = slab.cell.array
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    return slab, e, area


def main():
    print(f'POT_DIR = {POT_DIR}')
    print(f'LAMMPS  = {LAMMPS_BIN}')
    print(f'Bulk E/atom = {EPER_REF} eV', flush=True)

    cases = []

    # (001): try 4 z-shifts (0, 1/4, 1/2, 3/4)
    for shift in [0.0, 0.25, 0.50, 0.75]:
        cases.append({'face': '001', 'L': 12, 'z_shift': shift,
                      'label': f'001_L12_shift{int(shift*100):02d}'})
    # (001) parity: try L=11 and L=13 to expose different stacking ends
    cases.append({'face': '001', 'L': 11, 'z_shift': 0.0, 'label': '001_L11_shift00'})
    cases.append({'face': '001', 'L': 13, 'z_shift': 0.0, 'label': '001_L13_shift00'})

    # (111): try 4 z-shifts and parity
    for shift in [0.0, 0.25, 0.50, 0.75]:
        cases.append({'face': '111', 'L': 12, 'z_shift': shift,
                      'label': f'111_L12_shift{int(shift*100):02d}'})
    cases.append({'face': '111', 'L': 11, 'z_shift': 0.0, 'label': '111_L11_shift00'})
    cases.append({'face': '111', 'L': 13, 'z_shift': 0.0, 'label': '111_L13_shift00'})

    print()
    header = f'{"face":>5} {"L":>3} {"shift":>6} {"n":>4} {"area_A2":>10} {"E_eV":>14} {"gamma_mJ_m2":>14}'
    print(header, flush=True)
    print('-' * len(header), flush=True)

    results = {'potential': 'Ko 2018 MEAM (Sn)', 'E_per_atom_eV': EPER_REF, 'cases': []}

    for case in cases:
        miller = tuple(int(x) for x in case['face'])
        try:
            slab, e, area = relax_slab(miller, case['L'], case['z_shift'], case['label'])
            n = len(slab)
            gamma = (e - n * EPER_REF) / (2 * area) * EV_A2_TO_MJ_M2
            print(f"{case['face']:>5} {case['L']:>3} {case['z_shift']:>6.2f} "
                  f"{n:>4} {area:>10.3f} {e:>14.4f} {gamma:>14.2f}", flush=True)
            results['cases'].append({
                **case, 'n': n, 'area_A2': area, 'E_eV': e, 'gamma_mJ_m2': gamma
            })
        except Exception as exc:
            print(f"{case['face']:>5} {case['L']:>3} {case['z_shift']:>6.2f}  FAILED: {exc}", flush=True)
            results['cases'].append({**case, 'error': str(exc)})

    # Summary: best (lowest gamma) per face
    print()
    print('--- Best termination per face ---')
    for face in ['001', '111']:
        face_cases = [c for c in results['cases'] if c.get('face') == face and 'gamma_mJ_m2' in c]
        if face_cases:
            best = min(face_cases, key=lambda x: x['gamma_mJ_m2'])
            print(f"  ({face})  best γ = {best['gamma_mJ_m2']:.2f}  "
                  f"shift={best['z_shift']:.2f}  L={best['L']}")

    out_path = _os.path.join(PROJECT_ROOT, '/surface/data/ko2018_meam/result_surface_terminations.json'.lstrip('/'))
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nsaved -> {out_path}', flush=True)


if __name__ == '__main__':
    main()
