"""
06_make_dft_slabs_and_inputs.py (run on Matlantis)

Build β-Sn slabs using the DFT/PBE-optimized bulk lattice (a=5.970, c=3.218 Å,
from OpenMX elastic paper) as the reference. Generate OpenMX inputs:
  - bulk_dft_SP.dat: single-point at DFT bulk lattice (reference for γ)
  - slab_XXX_LN_Relax.dat: full-atom relaxation, MD BFGS, maxIter 1000, fmax 3e-4

k-mesh target Δk = 0.15 rad/Å (Acta Materialia 2026 compliant)
SCF: rmm-diisk, Max=0.03, Kerker=10, History=15, StartPulay=60, criterion=1e-7
"""
import json
import math
import os
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.io import write

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESDIR = os.path.join(ROOT, 'results_dft')
OUTDIR = os.path.join(ROOT, 'openmx_inputs_dft')
os.makedirs(RESDIR, exist_ok=True)
os.makedirs(OUTDIR, exist_ok=True)

# DFT/PBE-optimized lattice from OpenMX elastic paper (intel_results2)
DFT_LATTICE = {'a': 5.970, 'c': 3.218}

SLABS = [('100', 12), ('001', 8), ('110', 8), ('101', 12)]
VACUUM = 15.0


def kgrid_for(L_A, is_vacuum=False):
    if is_vacuum:
        return 1
    return max(2, int(math.ceil(2 * math.pi / (0.15 * L_A))))


def build_bulk():
    a, c = DFT_LATTICE['a'], DFT_LATTICE['c']
    cell = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    scaled = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
              [0.0, 0.5, 0.25], [0.5, 0.0, 0.75]]
    return Atoms('Sn4', scaled_positions=scaled, cell=cell, pbc=True)


def build_slab(bulk, miller, layers, vacuum):
    slab = ase_surface(bulk, miller, layers=layers, vacuum=vacuum, periodic=True)
    slab.center(vacuum=vacuum, axis=2)
    return slab


# ---------- OpenMX template ----------
HEADER_COMMON = """#
# {title}
# β-Sn {kind}
# k-mesh target Δk ≈ 0.15 rad/Å (Acta Materialia 2026 compliant)
#
System.CurrrentDirectory         ./
System.Name                      {sysname}
DATA.PATH                        __DFT_DATA_PATH__
level.of.stdout                  1
level.of.fileout                 0

Species.Number                   1
<Definition.of.Atomic.Species
 Sn  Sn7.0-s2p2d3f1  Sn_PBE19
Definition.of.Atomic.Species>

Atoms.Number                     {natoms}
Atoms.SpeciesAndCoordinates.Unit Ang
<Atoms.SpeciesAndCoordinates
{atom_lines}Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit           Ang
<Atoms.UnitVectors
{cell_lines}Atoms.UnitVectors>

scf.XcType                       GGA-PBE
scf.SpinPolarization             off
scf.ElectronicTemperature        300.0
scf.energycutoff                 200.0
scf.maxIter                      600
scf.EigenvalueSolver             Band
scf.Kgrid                        {kx} {ky} {kz}
scf.Mixing.Type                  rmm-diisk
scf.Init.Mixing.Weight           0.005
scf.Min.Mixing.Weight            0.001
scf.Max.Mixing.Weight            0.03
scf.Mixing.History               15
scf.Mixing.StartPulay            60
scf.Kerker.factor                10.0
scf.criterion                    1.0e-7

"""

MD_SP = """MD.Type                          Nomd
MD.maxIter                       1
"""

MD_RELAX = """MD.Type                          BFGS
MD.maxIter                       1000
MD.Opt.criterion                 3.0e-4
"""


def fmt_atoms(atoms):
    lines = ''
    for i, (sym, xyz) in enumerate(zip(atoms.get_chemical_symbols(), atoms.get_positions()), 1):
        lines += f"  {i:3d}   {sym}   {xyz[0]: .8f}  {xyz[1]: .8f}  {xyz[2]: .8f}   7.0  7.0\n"
    return lines


def fmt_cell(atoms):
    lines = ''
    for row in atoms.cell.array:
        lines += f"  {row[0]: .6f}  {row[1]: .6f}  {row[2]: .6f}\n"
    return lines


def write_openmx(path, sysname, title, kind, atoms, kgrid, md_block):
    content = HEADER_COMMON.format(
        title=title, kind=kind, sysname=sysname, natoms=len(atoms),
        atom_lines=fmt_atoms(atoms), cell_lines=fmt_cell(atoms),
        kx=kgrid[0], ky=kgrid[1], kz=kgrid[2],
    ) + md_block
    with open(path, 'w') as f:
        f.write(content)


def main():
    bulk = build_bulk()
    write(os.path.join(RESDIR, 'bulk_dft.cif'), bulk)
    write(os.path.join(RESDIR, 'bulk_dft.traj'), bulk)

    manifest = {
        'bulk': {'a': DFT_LATTICE['a'], 'c': DFT_LATTICE['c'],
                 'n_atoms': 4, 'source': 'OpenMX/PBE elastic paper intel_results2'},
        'slabs': [],
    }

    # bulk SP
    kb = (kgrid_for(DFT_LATTICE['a']),
          kgrid_for(DFT_LATTICE['a']),
          kgrid_for(DFT_LATTICE['c']))
    write_openmx(
        os.path.join(OUTDIR, 'bulk_dft_SP.dat'),
        'bulk_dft_SP',
        'β-Sn bulk SP at DFT/PBE-optimized lattice',
        f'bulk 4 atoms, a={DFT_LATTICE["a"]:.3f}, c={DFT_LATTICE["c"]:.3f}',
        bulk, kb, MD_SP,
    )
    manifest['bulk']['kgrid'] = list(kb)
    print(f'bulk_dft_SP.dat  n=4  kgrid={kb}')

    # slab Relax inputs
    for tag, L in SLABS:
        miller = tuple(int(c) for c in tag)
        slab = build_slab(bulk, miller, L, VACUUM)
        Lx = float(np.linalg.norm(slab.cell.array[0]))
        Ly = float(np.linalg.norm(slab.cell.array[1]))
        Lz = float(np.linalg.norm(slab.cell.array[2]))
        A = float(np.linalg.norm(np.cross(slab.cell.array[0], slab.cell.array[1])))
        kg = (kgrid_for(Lx), kgrid_for(Ly), 1)

        sysname = f'slab_{tag}_L{L}_Relax_dft'
        write_openmx(
            os.path.join(OUTDIR, sysname + '.dat'),
            sysname,
            f'β-Sn ({tag}) slab L={L} full-atom relaxation at DFT lattice',
            f'slab ({tag}) L={L}, vacuum {VACUUM:.0f} Å',
            slab, kg, MD_RELAX,
        )
        write(os.path.join(RESDIR, f'{sysname}.traj'), slab)

        manifest['slabs'].append({
            'sysname': sysname, 'tag': tag, 'layers': L, 'n_atoms': len(slab),
            'cell_Ang': [Lx, Ly, Lz], 'area_A2': A, 'kgrid': list(kg),
        })
        print(f'{sysname}.dat  n={len(slab)}  cell=({Lx:.2f},{Ly:.2f},{Lz:.2f})  kgrid={kg}')

    with open(os.path.join(OUTDIR, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f'\nwrote manifest -> {OUTDIR}/manifest.json')


if __name__ == '__main__':
    main()
