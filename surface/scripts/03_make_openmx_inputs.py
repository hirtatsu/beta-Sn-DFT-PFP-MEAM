
import json
import os
import numpy as np
from ase.io import read
from ase import Atoms

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESDIR = os.path.join(ROOT, 'results')
OUTDIR = os.path.join(ROOT, 'openmx_inputs')
os.makedirs(OUTDIR, exist_ok=True)

# Target slab configs (thinner than production for DFT feasibility)
SLABS = [
    ('100', 12),
    ('001', 8),
    ('110', 8),
    ('101', 12),
]

# k-grid target: Δk = 0.15 rad/Å (within Acta Materialia 2026 paper range 0.10-0.54)
# k = ceil(2π / (Δk × L)) = ceil(41.888 / L)
def kgrid_for(L_A, is_vacuum=False):
    import math
    if is_vacuum:
        return 1
    return max(2, int(math.ceil(2 * math.pi / (0.15 * L_A))))


OPENMX_HEAD = """#
# {title}
# β-Sn {kind}: single-point DFT/PBE (MD.Type Nomd)
# Basis/pp identical to Tatsumi elastic manuscript.
#
System.CurrrentDirectory         ./
System.Name                      {sysname}
DATA.PATH                        /path/to/work_dir/openmx_betaSn/DFT_DATA
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

MD.Type                          Nomd
MD.maxIter                       1
"""


def write_openmx(path, sysname, title, kind, atoms, kgrid):
    n = len(atoms)
    pos = atoms.get_positions()
    cell = atoms.cell.array

    atom_lines = ''
    for i, (sym, xyz) in enumerate(zip(atoms.get_chemical_symbols(), pos), 1):
        atom_lines += f"  {i:3d}   {sym}   {xyz[0]: .8f}  {xyz[1]: .8f}  {xyz[2]: .8f}   7.0  7.0\n"

    cell_lines = ''
    for row in cell:
        cell_lines += f"  {row[0]: .6f}  {row[1]: .6f}  {row[2]: .6f}\n"

    content = OPENMX_HEAD.format(
        title=title,
        kind=kind,
        sysname=sysname,
        natoms=n,
        atom_lines=atom_lines,
        cell_lines=cell_lines,
        kx=kgrid[0], ky=kgrid[1], kz=kgrid[2],
    )
    with open(path, 'w') as f:
        f.write(content)


def main():
    # (1) Bulk SP at PFP lattice
    with open(os.path.join(RESDIR, 'bulk_pfp_pbe_d3.json')) as f:
        bulk_info = json.load(f)
    a_pfp = bulk_info['a']
    c_pfp = bulk_info['c']
    cell_b = [[a_pfp, 0, 0], [0, a_pfp, 0], [0, 0, c_pfp]]
    scaled = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
              [0.0, 0.5, 0.25], [0.5, 0.0, 0.75]]
    bulk_atoms = Atoms('Sn4', scaled_positions=scaled, cell=cell_b, pbc=True)
    kb = (kgrid_for(a_pfp), kgrid_for(a_pfp), kgrid_for(c_pfp))
    write_openmx(
        os.path.join(OUTDIR, 'bulk_pfp_SP.dat'),
        'bulk_pfp_SP',
        'β-Sn bulk SP at PFP/PBE+D3 lattice',
        f'bulk 4 atoms, a={a_pfp:.4f}, c={c_pfp:.4f}',
        bulk_atoms, kb,
    )
    print(f'bulk_pfp_SP.dat  n=4  kgrid={kb}')

    # (2) slabs
    manifest = []
    for tag, L in SLABS:
        trajp = os.path.join(RESDIR, f'slab_{tag}_L{L}_V15.traj')
        slab = read(trajp)
        Lx = float(np.linalg.norm(slab.cell.array[0]))
        Ly = float(np.linalg.norm(slab.cell.array[1]))
        Lz = float(np.linalg.norm(slab.cell.array[2]))
        kg = (kgrid_for(Lx), kgrid_for(Ly), 1)
        sysname = f'slab_{tag}_L{L}_SP'
        write_openmx(
            os.path.join(OUTDIR, sysname + '.dat'),
            sysname,
            f'β-Sn ({tag}) slab L={L} SP',
            f'slab ({tag}) L={L}, vacuum 15 Å',
            slab, kg,
        )
        print(f'{sysname}.dat  n={len(slab)}  cell=({Lx:.2f},{Ly:.2f},{Lz:.2f})  kgrid={kg}')
        manifest.append({
            'sysname': sysname, 'tag': tag, 'layers': L,
            'n_atoms': len(slab),
            'cell_Ang': [Lx, Ly, Lz],
            'area_A2': float(np.linalg.norm(np.cross(slab.cell.array[0], slab.cell.array[1]))),
            'kgrid': list(kg),
        })

    with open(os.path.join(OUTDIR, 'manifest.json'), 'w') as f:
        json.dump({'bulk': {'a': a_pfp, 'c': c_pfp, 'kgrid': list(kb)},
                   'slabs': manifest}, f, indent=2)
    print(f'\nwrote manifest -> {OUTDIR}/manifest.json')


if __name__ == '__main__':
    main()
