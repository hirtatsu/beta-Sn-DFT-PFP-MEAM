"""Evaluate γ for (112) slab candidates with PFP/PBE on Matlantis."""
import os, json, numpy as np
from ase.io import read
from ase.optimize import LBFGS
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, 'slabs_112_check')
EV_A2_TO_MJ_M2 = 16021.766

# PFP/PBE bulk e_per_atom (cached)
with open(os.path.join(ROOT, 'results_modes', 'modes_surface_energies.json')) as f:
    modes = json.load(f)
e_per_atom = modes['PBE']['bulk']['E_per_atom_eV']
print(f'PFP/PBE bulk: E/atom = {e_per_atom:.4f} eV')

est = Estimator(calc_mode='PBE')
calc = ASECalculator(est)

with open(os.path.join(SRC, 'manifest.json')) as f:
    manifest = json.load(f)

out = []
for rec in manifest:
    fpath = os.path.join(SRC, rec['file'])
    atoms = read(fpath)
    atoms.calc = calc
    try:
        LBFGS(atoms, logfile=None).run(fmax=0.015, steps=300)
    except Exception as e:
        print(f"  {rec['file']} FAILED: {e}")
        continue
    E = float(atoms.get_potential_energy())
    cell = atoms.cell.array
    A = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    n = len(atoms)
    gamma = (E - n * e_per_atom) / (2 * A) * EV_A2_TO_MJ_M2
    c_thick = float(cell[2][2]) - 15  # subtract vacuum approx
    row = dict(rec, E_eV=E, area_A2=A, slab_thick_A=c_thick, gamma_mJ_m2=gamma)
    out.append(row)
    print(f"  {rec['file']:30s}  n={n:>3}  A={A:>6.2f}  thickness~{c_thick:.1f}  γ={gamma:.1f} mJ/m²")

# Compare with ase.build.surface L=12 we used (48 atoms, A=76.28) → PFP/PBE γ=442.5 per earlier run
print()
print('reference (ase.build.surface, L=12): n=48, A=76.28, γ_PFP/PBE=442.5 mJ/m²')

with open(os.path.join(SRC, 'pfp_results.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(f'saved -> {SRC}/pfp_results.json')
