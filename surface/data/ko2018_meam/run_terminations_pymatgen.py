"""Enumerate all unique surface terminations of beta-Sn (001) and (111) via pymatgen,
and compute gamma for each with Ko 2018 MEAM (lmp_kokkos)."""
import os
import json
import numpy as np
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import LBFGS
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

import os as _os
PROJECT_ROOT = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), "..", "..", ".."))
POT_DIR = _os.path.join(PROJECT_ROOT, "meam_potentials")
LAMMPS_BIN = _os.environ.get("LAMMPS_BIN", "lmp_kokkos")
EV_A2_TO_MJ_M2 = 16021.766
EPER_REF = -3.10195  # eV/atom from previous bulk relaxation
A0, C0 = 5.8589, 3.2059

# Faces to investigate alternative terminations
FACES = [(0, 0, 1), (1, 1, 1)]

# Slab parameters: enough thickness for convergence, big vacuum
MIN_SLAB_SIZE = 15.0  # Ang
MIN_VACUUM_SIZE = 18.0  # Ang


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


def build_bulk_pymatgen():
    lattice = Lattice.tetragonal(A0, C0)
    species = ['Sn'] * 4
    coords = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.0, 0.5, 0.25],
        [0.5, 0.0, 0.75],
    ]
    return Structure(lattice, species, coords)


def main():
    bulk = build_bulk_pymatgen()
    print(f'Bulk: a={A0}, c={C0}, formula={bulk.formula}, n_sites={len(bulk)}', flush=True)
    print(f'Bulk E/atom (ref) = {EPER_REF} eV', flush=True)
    print(flush=True)

    adaptor = AseAtomsAdaptor()
    results = {'potential': 'Ko 2018 MEAM (Sn)', 'E_per_atom_eV': EPER_REF, 'cases': []}

    for miller in FACES:
        face_tag = ''.join(str(x) for x in miller)
        print(f'\n##### face ({face_tag}) #####', flush=True)
        gen = SlabGenerator(
            bulk,
            miller,
            min_slab_size=MIN_SLAB_SIZE,
            min_vacuum_size=MIN_VACUUM_SIZE,
            center_slab=True,
            primitive=False,
            in_unit_planes=False,
        )
        # symmetrize=True → only return slabs with equivalent top & bottom terminations
        # symmetrize=False → list all distinct top-side terminations
        slabs = gen.get_slabs(symmetrize=False)
        print(f'  pymatgen returned {len(slabs)} candidate terminations', flush=True)

        for i, slab in enumerate(slabs):
            label = f'{face_tag}_term{i}'
            atoms = adaptor.get_atoms(slab)
            n = len(atoms)
            cell = atoms.cell.array
            area = float(np.linalg.norm(np.cross(cell[0], cell[1])))

            try:
                atoms.calc = make_calc(label)
                e_init = float(atoms.get_potential_energy())
                LBFGS(atoms, logfile=None).run(fmax=0.005, steps=300)
                e = float(atoms.get_potential_energy())
                gamma = (e - n * EPER_REF) / (2 * area) * EV_A2_TO_MJ_M2
                print(f'  [term {i}] n={n}  area={area:.3f}  E={e:.4f}  γ={gamma:.2f} mJ/m²',
                      flush=True)
                results['cases'].append({
                    'face': face_tag, 'term_index': i, 'n_atoms': n,
                    'area_A2': area, 'E_init_eV': e_init, 'E_relaxed_eV': e,
                    'gamma_mJ_m2': gamma, 'symmetric': False,
                })
            except Exception as exc:
                print(f'  [term {i}] FAILED: {exc}', flush=True)
                results['cases'].append({
                    'face': face_tag, 'term_index': i, 'error': str(exc)
                })

    # Best (lowest γ) per face
    print('\n--- Lowest γ per face ---')
    for miller in FACES:
        face_tag = ''.join(str(x) for x in miller)
        face_cases = [c for c in results['cases']
                      if c.get('face') == face_tag and 'gamma_mJ_m2' in c]
        if face_cases:
            best = min(face_cases, key=lambda x: x['gamma_mJ_m2'])
            print(f'  ({face_tag})  best γ = {best["gamma_mJ_m2"]:.2f} mJ/m²  '
                  f'(term_index={best["term_index"]}, n={best["n_atoms"]})')

    out_path = (_os.path.join(PROJECT_ROOT, 'surface', 'data')+'/'
                'ko2018_meam/result_terminations_pymatgen.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nsaved -> {out_path}', flush=True)


if __name__ == '__main__':
    main()
