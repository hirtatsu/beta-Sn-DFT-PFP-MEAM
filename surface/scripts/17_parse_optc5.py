"""Parse OptC5 results and compute γ_surf compared to BFGS."""
import os
import re
import numpy as np

HA_EV = 27.211386245988
EV_A2_TO_MJ_M2 = 16021.766

BASE = os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/supercomputer')

E_bulk_Ha = -284.157741623577
E_bulk_per_atom_eV = (E_bulk_Ha / 4.0) * HA_EV

# BFGS results (from earlier parse)
bfgs = {
    '001': {'E_Ha': -2273.171923566892, 'A0_A2': 35.6409, 'n': 32, 'a0': 5.97, 'b0': 5.97},
    '110': {'E_Ha': -2273.196264822546, 'A0_A2': 27.1691, 'n': 32, 'a0': 8.442855, 'b0': 3.218},
    '100': {'E_Ha': -3409.849519041379, 'A0_A2': 19.2115, 'n': 48, 'a0': 5.97,    'b0': 3.218},
    '101': {'E_Ha': -3409.798996954964, 'A0_A2': 40.4889, 'n': 48, 'a0': 6.782066,'b0': 5.97},
}


def parse_md_last_frame(md_path):
    """Return (E_Ha, cell 3x3 Å, n_atoms)."""
    with open(md_path) as f:
        lines = f.readlines()
    # find the last natom line
    last_idx = -1
    for i in range(len(lines) - 1, -1, -1):
        s = lines[i].strip()
        if s.isdigit():
            last_idx = i
            break
    n = int(lines[last_idx].strip())
    header = lines[last_idx + 1]
    m_E = re.search(r'Energy=\s*(-?\d+\.\d+)', header)
    m_cv = re.search(r'Cell_Vectors=\s*([-\d.\s]+)', header)
    E = float(m_E.group(1))
    cv_vals = [float(x) for x in m_cv.group(1).split()[:9]]
    cell = np.array(cv_vals).reshape(3, 3)
    return E, cell, n


print(f"Bulk reference (DFT/PBE BFGS): E/atom = {E_bulk_per_atom_eV:.4f} eV\n")
print(f"{'face':>5}  {'method':>8}  {'a':>7} {'b':>7} {'A(A2)':>8}  "
      f"{'E_slab(eV)':>12}  {'dE(eV)':>7}  {'gamma':>8}")
print('-' * 78)

results = []
for face in ['001', '110', '100', '101']:
    bf = bfgs[face]
    # BFGS
    E_bfgs_eV = bf['E_Ha'] * HA_EV
    dE_bfgs = E_bfgs_eV - bf['n'] * E_bulk_per_atom_eV
    gamma_bfgs = dE_bfgs / (2 * bf['A0_A2']) * EV_A2_TO_MJ_M2
    print(f"({face})  BFGS    {bf['a0']:>7.3f} {bf['b0']:>7.3f} {bf['A0_A2']:>8.3f}  "
          f"{E_bfgs_eV:>12.2f}  {dE_bfgs:>7.3f}  {gamma_bfgs:>8.1f}")
    # OptC5
    name = f"slab_{face}_L{'8' if face in ['001','110'] else '12'}_Relax"
    md_path = f"{BASE}/{name}_CellRelax/{name}.md"
    E_Ha, cell, n = parse_md_last_frame(md_path)
    A = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    a_len = float(np.linalg.norm(cell[0]))
    b_len = float(np.linalg.norm(cell[1]))
    E_opt_eV = E_Ha * HA_EV
    dE_opt = E_opt_eV - n * E_bulk_per_atom_eV
    gamma_opt = dE_opt / (2 * A) * EV_A2_TO_MJ_M2
    dA_pct = 100 * (A - bf['A0_A2']) / bf['A0_A2']
    print(f"({face})  OptC5   {a_len:>7.3f} {b_len:>7.3f} {A:>8.3f}  "
          f"{E_opt_eV:>12.2f}  {dE_opt:>7.3f}  {gamma_opt:>8.1f}"
          f"  (ΔA {dA_pct:+.1f}%)")
    results.append({'face': face,
                    'BFGS': {'a': bf['a0'], 'b': bf['b0'], 'A': bf['A0_A2'], 'gamma': gamma_bfgs},
                    'OptC5': {'a': a_len, 'b': b_len, 'A': A, 'gamma': gamma_opt}})
    print()

print()
print(f"{'face':>5}  {'γ_BFGS':>8}  {'γ_OptC5':>8}  {'Δγ':>8}")
for r in results:
    dg = r['OptC5']['gamma'] - r['BFGS']['gamma']
    print(f"({r['face']})  {r['BFGS']['gamma']:>8.1f}  {r['OptC5']['gamma']:>8.1f}  {dg:>+8.1f}")

import json
with open(os.path.expanduser('~/projects/beta-sn-elasticity/surface_gb/optc5_results.json'), 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nsaved -> ~/projects/beta-sn-elasticity/surface_gb/optc5_results.json")
