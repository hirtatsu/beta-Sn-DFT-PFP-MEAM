"""More thorough termination scan: vary slab thickness, vacuum, and pymatgen tol."""
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
EPER_REF = -3.10195
A0, C0 = 5.8589, 3.2059


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


def build_bulk():
    return Structure(
        Lattice.tetragonal(A0, C0),
        ['Sn'] * 4,
        [[0,0,0], [0.5,0.5,0.5], [0,0.5,0.25], [0.5,0,0.75]]
    )


def relax_and_gamma(atoms, label, fmax=0.001, steps=500):
    n = len(atoms)
    cell = atoms.cell.array
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    atoms.calc = make_calc(label)
    LBFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
    e = float(atoms.get_potential_energy())
    return n, area, e, (e - n * EPER_REF) / (2 * area) * EV_A2_TO_MJ_M2


def main():
    bulk = build_bulk()
    adaptor = AseAtomsAdaptor()
    print(f'Bulk: a={A0}, c={C0}, n_sites={len(bulk)}', flush=True)
    print(f'NIST IPR (001) = 395.47, (111) = 458.84 mJ/m²', flush=True)
    print()

    results = {'cases': []}

    for face in [(0, 0, 1), (1, 1, 1)]:
        face_tag = ''.join(str(x) for x in face)
        print(f'\n##### face ({face_tag}) — varying slab parameters #####', flush=True)
        for min_slab in [12, 18, 24]:
            for min_vac in [18, 25]:
                for sym in [False, True]:
                    gen = SlabGenerator(
                        bulk, face,
                        min_slab_size=min_slab,
                        min_vacuum_size=min_vac,
                        center_slab=True,
                        primitive=False,
                        in_unit_planes=False,
                    )
                    try:
                        slabs = gen.get_slabs(symmetrize=sym, ftol=0.01)
                    except Exception as exc:
                        print(f'  slab=[{min_slab},{min_vac}] sym={sym}  FAILED to gen: {exc}', flush=True)
                        continue
                    for i, slab in enumerate(slabs):
                        atoms = adaptor.get_atoms(slab)
                        label = f'{face_tag}_s{min_slab}_v{min_vac}_sym{int(sym)}_t{i}'
                        try:
                            n, area, e, gamma = relax_and_gamma(atoms, label, fmax=0.001)
                            print(f'  slab={min_slab} vac={min_vac} sym={sym} term={i} '
                                  f'n={n:>3} γ={gamma:>7.2f} mJ/m²', flush=True)
                            results['cases'].append({
                                'face': face_tag, 'min_slab': min_slab, 'min_vac': min_vac,
                                'symmetric': sym, 'term_index': i,
                                'n_atoms': n, 'area_A2': area, 'E_eV': e,
                                'gamma_mJ_m2': gamma,
                            })
                        except Exception as exc:
                            print(f'  slab={min_slab} vac={min_vac} sym={sym} term={i}  FAILED: {exc}',
                                  flush=True)
                            results['cases'].append({
                                'face': face_tag, 'min_slab': min_slab, 'min_vac': min_vac,
                                'symmetric': sym, 'term_index': i, 'error': str(exc),
                            })

    # Summary: best per face
    print('\n--- Range of γ found ---')
    for face_tag in ['001', '111']:
        gammas = [c['gamma_mJ_m2'] for c in results['cases']
                  if c.get('face') == face_tag and 'gamma_mJ_m2' in c]
        if gammas:
            print(f'  ({face_tag})  γ range: {min(gammas):.2f} ~ {max(gammas):.2f} mJ/m²  '
                  f'(n_cases={len(gammas)})')

    out_path = (_os.path.join(PROJECT_ROOT, 'surface', 'data')+'/'
                'ko2018_meam/result_terminations_v2.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nsaved -> {out_path}', flush=True)


if __name__ == '__main__':
    main()
