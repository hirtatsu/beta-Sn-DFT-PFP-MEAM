"""Quick diagnostic: does a pristine CSL-cut cell have the same E_per_atom as the 4-atom bulk?"""
import json
import os
import numpy as np
from ase.build import cut
from ase.io import read
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
bulk = read(os.path.join(ROOT, 'results', 'bulk_pfp_pbe_d3.traj'))
with open(os.path.join(ROOT, 'results', 'bulk_pfp_pbe_d3.json')) as f:
    info = json.load(f)
e_ref = info['E_per_atom_eV']

est = Estimator(calc_mode='PBE_PLUS_D3')
calc = ASECalculator(est)

print(f'reference (4-atom bulk): E/atom = {e_ref:.6f} eV  (vol/atom = {bulk.get_volume()/4:.3f} Å³)')
print()

for basis in [(2, 1, 0), (1, 1, 0), (3, 1, 0)]:
    tag = ''.join(str(x) for x in basis)
    other = (-basis[1], basis[0], 0)
    csl = cut(bulk, a=basis, b=other, c=(0, 0, 1))
    csl.calc = calc
    n = len(csl)
    E = csl.get_potential_energy()
    dets = abs(np.linalg.det(np.array([basis, other, (0, 0, 1)])))
    print(f'cut({tag},{other},001): n={n}, expected={int(dets)*4}, '
          f'E/atom={E/n:.6f} eV, diff={E/n - e_ref:+.6f}, vol/atom={csl.get_volume()/n:.3f}')
