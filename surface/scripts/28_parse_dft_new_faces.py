"""Parse the 3 new DFT Relax results (111, 211, 112) and compute γ."""
import os
import re
import numpy as np
from ase.io import read

HA_EV = 27.211386245988
EV_A2_TO_MJ_M2 = 16021.766

BASE = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/supercomputer')
E_bulk_Ha = -284.157741623577
E_bulk_per_atom_eV = (E_bulk_Ha / 4.0) * HA_EV


def parse_md_last(md_path):
    with open(md_path) as f:
        lines = f.readlines()
    last = -1
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].strip().isdigit():
            last = i
            break
    n = int(lines[last].strip())
    header = lines[last + 1]
    E = float(re.search(r'Energy=\s*(-?\d+\.\d+)', header).group(1))
    cv = [float(x) for x in re.search(r'Cell_Vectors=\s*([-\d.\s]+)', header).group(1).split()[:9]]
    cell = np.array(cv).reshape(3, 3)
    return E, cell, n


print(f"E_bulk/atom = {E_bulk_per_atom_eV:.4f} eV")
print(f"{'face':>5} {'n':>3} {'A(A²)':>8} {'E(Ha)':>14} {'dE(eV)':>8} {'γ (mJ/m²)':>11}")
print('-' * 62)

results = {}
for face in ['111', '211', '112']:
    md = f"{BASE}/slab_{face}_L12_Relax/slab_{face}_L12_Relax.md"
    E_Ha, cell, n = parse_md_last(md)
    A = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    E_eV = E_Ha * HA_EV
    dE = E_eV - n * E_bulk_per_atom_eV
    gamma = dE / (2 * A) * EV_A2_TO_MJ_M2
    print(f"({face}) {n:>3} {A:>8.3f} {E_Ha:>14.6f} {dE:>8.3f} {gamma:>11.1f}")
    results[face] = {'E_Ha': E_Ha, 'n': n, 'A_A2': A, 'gamma_mJ_m2': gamma}

import json
with open(os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/dft_new_faces.json'), 'w') as f:
    json.dump(results, f, indent=2)
