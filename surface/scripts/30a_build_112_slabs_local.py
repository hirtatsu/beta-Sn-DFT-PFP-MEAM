"""Locally build (112) slab candidates via pymatgen SlabGenerator.
Emit ase .traj files for Matlantis to evaluate γ with PFP.
"""
import os
import json
from ase.io import write
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

OUT = os.path.expanduser('~/tmp/112_slabs')
os.makedirs(OUT, exist_ok=True)

# PFP/PBE-optimized β-Sn bulk (from previous run)
A, C = 5.9295, 3.2008


def build_bulk():
    latt = Lattice.from_parameters(A, A, C, 90, 90, 90)
    species = ['Sn']*4
    coords = [[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]]
    return Structure(latt, species, coords, coords_are_cartesian=False)


def slabs_for(miller, min_slab, min_vac, label_prefix):
    bulk = build_bulk()
    sg = SlabGenerator(bulk, miller_index=miller,
                       min_slab_size=min_slab, min_vacuum_size=min_vac,
                       lll_reduce=False, center_slab=True,
                       in_unit_planes=False, primitive=True,
                       reorient_lattice=True)
    slabs = sg.get_slabs(symmetrize=True, tol=0.01)
    print(f'{label_prefix}: got {len(slabs)} symmetric slabs (min_slab={min_slab}, vac={min_vac})')
    ada = AseAtomsAdaptor()
    records = []
    for i, s in enumerate(slabs):
        atoms = ada.get_atoms(s)
        fname = f'{label_prefix}_t{i}.traj'
        write(os.path.join(OUT, fname), atoms)
        print(f'  #{i}: n={len(atoms)}, cell={atoms.cell.lengths()[0]:.2f}x{atoms.cell.lengths()[1]:.2f}x{atoms.cell.lengths()[2]:.2f}')
        records.append({'label': label_prefix, 'term': i, 'file': fname,
                        'min_slab_A': min_slab, 'n_atoms': len(atoms),
                        'cell': atoms.cell.lengths().tolist()})
    return records


def main():
    all_rec = []
    # (a) terminations at 20 A slab
    all_rec += slabs_for((1,1,2), 20, 15, 'slab_112_term')
    # (b) thickness scan (using first termination, typically the low-γ one)
    for ms in [12, 18, 24, 30, 36]:
        all_rec += slabs_for((1,1,2), ms, 15, f'slab_112_thick{ms}')
    with open(os.path.join(OUT, 'manifest.json'), 'w') as f:
        json.dump(all_rec, f, indent=2)
    print(f'\nwrote {len(all_rec)} structures into {OUT}')


if __name__ == '__main__':
    main()
