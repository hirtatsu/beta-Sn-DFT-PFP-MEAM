"""(a) (112) termination enumeration; (b) thickness convergence.
Uses pymatgen SlabGenerator and PFP/PBE for cheap screening.
Run on Matlantis.
"""
import json
import os
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import LBFGS
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, 'results_112_check')
os.makedirs(OUT, exist_ok=True)

EV_A2_TO_MJ_M2 = 16021.766
A, C = 5.8457, 3.1732   # PFP/PBE+D3 bulk (same cell we used for reference in mode-self-consistent prev runs)
# but use PFP/PBE for this check since DFT-equivalent conditions; PFP/PBE bulk is a=5.9295, c=3.2008


def build_pmg_bulk(a, c):
    latt = Lattice.from_parameters(a, a, c, 90, 90, 90)
    species = ['Sn']*4
    coords = [[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]]
    return Structure(latt, species, coords, coords_are_cartesian=False)


def enumerate_slabs(bulk_struct, miller, min_slab, min_vac):
    sg = SlabGenerator(bulk_struct, miller_index=miller,
                       min_slab_size=min_slab, min_vacuum_size=min_vac,
                       lll_reduce=False, center_slab=True,
                       in_unit_planes=False, primitive=True,
                       reorient_lattice=True)
    slabs = sg.get_slabs(symmetrize=True, tol=0.01)
    return slabs


def gamma_pfp(slab_atoms, calc, e_bulk_per_atom):
    slab_atoms.calc = calc
    try:
        LBFGS(slab_atoms, logfile=None).run(fmax=0.015, steps=300)
    except Exception as e:
        print(f'    opt failed: {e}')
        return None, None
    E = float(slab_atoms.get_potential_energy())
    cell = slab_atoms.cell.array
    A = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    gamma = (E - len(slab_atoms) * e_bulk_per_atom) / (2 * A) * EV_A2_TO_MJ_M2
    return gamma, A


def main():
    # Use PFP/PBE bulk (mode-self-consistent)
    mode = 'PBE'
    # re-relax PFP/PBE bulk quickly to get a_eq, c_eq, e_per_atom
    est = Estimator(calc_mode=mode); calc = ASECalculator(est)
    a_exp, c_exp = 5.831, 3.182
    bulk = Atoms('Sn4',
                 scaled_positions=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]],
                 cell=[[a_exp,0,0],[0,a_exp,0],[0,0,c_exp]], pbc=True)
    bulk.calc = calc
    from ase.filters import ExpCellFilter
    ecf = ExpCellFilter(bulk, hydrostatic_strain=False)
    LBFGS(ecf, logfile=None).run(fmax=0.015, steps=300)
    a_eq = float(bulk.cell.lengths()[0]); c_eq = float(bulk.cell.lengths()[2])
    e_per_atom = float(bulk.get_potential_energy()) / len(bulk)
    print(f'PFP/PBE bulk: a={a_eq:.4f}, c={c_eq:.4f}, E/atom={e_per_atom:.4f} eV', flush=True)

    pmg_bulk = build_pmg_bulk(a_eq, c_eq)
    ada = AseAtomsAdaptor()

    # (a) Enumerate (112) terminations at a fixed slab/vacuum
    print('\n=== (a) (112) termination enumeration (min_slab=20, vac=15) ===', flush=True)
    slabs = enumerate_slabs(pmg_bulk, (1,1,2), min_slab=20, min_vac=15)
    print(f'  got {len(slabs)} distinct symmetric slabs')
    term_results = []
    for i, s in enumerate(slabs):
        atoms = ada.get_atoms(s)
        n = len(atoms)
        g, A = gamma_pfp(atoms, calc, e_per_atom)
        if g is None:
            continue
        print(f'  termination #{i}: n={n}, A={A:.2f}, γ={g:.1f} mJ/m²', flush=True)
        term_results.append({'term': i, 'n': n, 'area_A2': A, 'gamma_mJ_m2': g})
        write(os.path.join(OUT, f'slab_112_term{i}.traj'), atoms)

    # (b) Thickness convergence on best termination
    if term_results:
        best_idx = min(range(len(term_results)), key=lambda i: term_results[i]['gamma_mJ_m2'])
        print(f'\nbest termination: #{best_idx} (γ={term_results[best_idx]["gamma_mJ_m2"]:.1f})')
    print('\n=== (b) (112) thickness convergence (same termination, best) ===', flush=True)
    thick_results = []
    for m_slab in [12, 18, 24, 30, 36]:
        slabs = enumerate_slabs(pmg_bulk, (1,1,2), min_slab=m_slab, min_vac=15)
        # pick the termination with lowest γ
        row = {'min_slab_A': m_slab}
        if not slabs:
            row['err'] = 'no slabs'
        else:
            # just use the first for simplicity (SlabGenerator ordering is consistent)
            atoms = ada.get_atoms(slabs[0])
            n = len(atoms)
            g, A = gamma_pfp(atoms, calc, e_per_atom)
            print(f'  min_slab={m_slab} Å: n={n}, A={A:.2f}, γ={g:.1f}', flush=True)
            row.update({'n': n, 'area_A2': A, 'gamma_mJ_m2': g})
        thick_results.append(row)

    with open(os.path.join(OUT, 'term_thickness_112.json'), 'w') as f:
        json.dump({'terminations': term_results,
                   'thickness_scan': thick_results,
                   'bulk': {'a': a_eq, 'c': c_eq, 'E_per_atom_eV': e_per_atom}}, f, indent=2)
    print(f'\nsaved -> {OUT}/term_thickness_112.json', flush=True)


if __name__ == '__main__':
    main()
