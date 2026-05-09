"""Surface energies of beta-Sn with Ko 2018 MEAM via lmp_kokkos.

Comparison target: NIST IPR predicted values for the same potential.
"""
import os
import json
import shutil
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
VACUUM = 18.0  # one-sided vacuum thickness (Ang); slab.center adds 2*VACUUM total

# Ko 2018 NIST IPR predicted bulk
A0, C0 = 5.8589, 3.2059
EPER_REF = -3.1019  # eV/atom

# Slab settings — use enough layers for convergence
FACES = [
    ('001', 12),
    ('100', 12),
    ('110', 8),
    ('101', 12),
    ('111', 12),
    ('112', 12),
    ('211', 12),
]


def build_beta_sn(a, c):
    """Conventional 4-atom tetragonal cell of beta-Sn (A5, I4_1/amd)."""
    return Atoms(
        'Sn4',
        scaled_positions=[
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.0, 0.5, 0.25],
            [0.5, 0.0, 0.75],
        ],
        cell=[[a, 0, 0], [0, a, 0], [0, 0, c]],
        pbc=True,
    )


def make_calc(label, files=None):
    """Configure ASE LAMMPSrun calculator using lmp_kokkos + Ko 2018 MEAM."""
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


def get_bulk_energy_per_atom():
    """Quick bulk single-point at NIST IPR equilibrium lattice."""
    bulk = build_beta_sn(A0, C0)
    bulk.calc = make_calc('bulk_sp')
    e = float(bulk.get_potential_energy())
    return e / len(bulk)


def relax_slab(miller, layers):
    bulk = build_beta_sn(A0, C0)
    slab = ase_surface(bulk, miller, layers=layers, vacuum=VACUUM, periodic=True)
    slab.center(vacuum=VACUUM, axis=2)
    tag = ''.join(str(x) for x in miller)
    slab.calc = make_calc(f'slab_{tag}_L{layers}')
    e_initial = float(slab.get_potential_energy())
    # Atoms-only relaxation, cell fixed (LAMMPS default for ase calc)
    LBFGS(slab, logfile=None).run(fmax=0.005, steps=300)
    e_final = float(slab.get_potential_energy())
    cell = slab.cell.array
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    return slab, e_initial, e_final, area


def main():
    os.makedirs('/tmp/lammps_output', exist_ok=True)
    print(f'POT_DIR = {POT_DIR}')
    print(f'LAMMPS  = {LAMMPS_BIN}')
    print(f'Bulk    = a={A0} c={C0}, ref E/atom = {EPER_REF} eV')
    print(flush=True)

    print('--- bulk single-point ---', flush=True)
    e_per_atom = get_bulk_energy_per_atom()
    print(f'  E/atom = {e_per_atom:.6f} eV (ref {EPER_REF})', flush=True)

    results = {
        'potential': 'Ko 2018 MEAM (Sn)',
        'lammps': '22 Jul 2025 Update 4 (lmp_kokkos)',
        'bulk': {'a_Ang': A0, 'c_Ang': C0, 'E_per_atom_eV': e_per_atom},
        'faces': []
    }
    print()
    print(f'{"face":>5} {"L":>3} {"n":>4} {"area_A2":>10} {"E_eV":>14} {"gamma_mJ_m2":>14}', flush=True)
    print('-' * 60, flush=True)

    for tag, L in FACES:
        miller = tuple(int(x) for x in tag)
        try:
            slab, e_init, e, area = relax_slab(miller, L)
            n = len(slab)
            gamma = (e - n * e_per_atom) / (2 * area) * EV_A2_TO_MJ_M2
            print(f'{tag:>5} {L:>3} {n:>4} {area:>10.3f} {e:>14.4f} {gamma:>14.2f}', flush=True)
            results['faces'].append({
                'face': tag, 'L': L, 'n_atoms': n, 'area_A2': area,
                'E_init_eV': e_init, 'E_relaxed_eV': e, 'gamma_mJ_m2': gamma
            })
        except Exception as exc:
            print(f'{tag:>5} {L:>3}  FAILED: {exc}', flush=True)
            results['faces'].append({'face': tag, 'L': L, 'error': str(exc)})

    out_path = _os.path.join(_os.path.dirname(__file__), 'result_surface_energies.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nsaved -> {out_path}', flush=True)


if __name__ == '__main__':
    main()
