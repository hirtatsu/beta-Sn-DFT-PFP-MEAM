import os
from ase.io import read

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
atoms = read(os.path.join(PROJECT_ROOT, "beta-Sn.cif"))

print("=== β-Sn (white tin) structure ===")
print(f"Chemical formula  : {atoms.get_chemical_formula()}")
print(f"Number of atoms   : {len(atoms)}")
print(f"Cell parameters   : a={atoms.cell.lengths()[0]:.4f}, "
      f"b={atoms.cell.lengths()[1]:.4f}, "
      f"c={atoms.cell.lengths()[2]:.4f} Å")
print(f"Cell angles       : α={atoms.cell.angles()[0]:.2f}, "
      f"β={atoms.cell.angles()[1]:.2f}, "
      f"γ={atoms.cell.angles()[2]:.2f} deg")
print(f"Volume            : {atoms.get_volume():.4f} Å^3")
print(f"Density           : {sum(atoms.get_masses()) / atoms.get_volume() * 1.66054:.4f} g/cm^3")
print(f"PBC               : {atoms.pbc}")

print("\n--- Lattice vectors (Å) ---")
for i, v in enumerate(atoms.cell):
    print(f"  a{i+1} = [{v[0]:8.4f} {v[1]:8.4f} {v[2]:8.4f}]")

print("\n--- Atomic positions ---")
print(f"{'idx':>3} {'elem':>4} {'x':>10} {'y':>10} {'z':>10}   (fractional)")
for i, (sym, pos) in enumerate(zip(atoms.get_chemical_symbols(),
                                   atoms.get_scaled_positions())):
    print(f"{i:>3} {sym:>4} {pos[0]:10.5f} {pos[1]:10.5f} {pos[2]:10.5f}")
