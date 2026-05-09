"""
04_make_openmx_relax_inputs.py (run on Matlantis)
Produces OpenMX atom-relaxation inputs for 4 slabs.
Cell is fixed to PFP-optimized lattice (consistent with bulk SP reference).
Middle 30% atoms (by z) are frozen to match PFP's slab relaxation setup.
"""
import json
import os
import numpy as np
from ase.io import read

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESDIR = os.path.join(ROOT, 'results')
OUTDIR = os.path.join(ROOT, 'openmx_inputs_relax')
os.makedirs(OUTDIR, exist_ok=True)

SLABS = [('100', 12), ('001', 8), ('110', 8), ('101', 12)]


def kgrid_for(L_A, is_vacuum=False):
    import math
    if is_vacuum:
        return 1
    return max(2, int(math.ceil(2 * math.pi / (0.15 * L_A))))


TEMPLATE = """#
# β-Sn ({tag}) slab L={L} atom relaxation (cell fixed to PFP lattice)
# middle 30% atoms frozen; BFGS until fmax < 2e-3 Ha/Bohr
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

MD.Type                          BFGS
MD.maxIter                       100
MD.Opt.criterion                 3.0e-4

<MD.Fixed.XYZ
{fix_lines}MD.Fixed.XYZ>
"""


def main():
    for tag, L in SLABS:
        trajp = os.path.join(RESDIR, f'slab_{tag}_L{L}_V15.traj')
        slab = read(trajp)
        pos = slab.get_positions()
        cell = slab.cell.array
        Lx = float(np.linalg.norm(cell[0]))
        Ly = float(np.linalg.norm(cell[1]))

        # middle 30% by z
        z = pos[:, 2]
        zmid = 0.5 * (z.min() + z.max())
        half_thick = 0.5 * (z.max() - z.min())
        freeze = np.abs(z - zmid) < 0.15 * 2 * half_thick  # |z-zmid|/thick < 0.15

        atom_lines = ''
        fix_lines = ''
        for i, (sym, xyz, f) in enumerate(zip(slab.get_chemical_symbols(), pos, freeze), 1):
            atom_lines += f"  {i:3d}   {sym}   {xyz[0]: .8f}  {xyz[1]: .8f}  {xyz[2]: .8f}   7.0  7.0\n"
            fx = 1 if f else 0
            fix_lines += f"  {i:3d}  {fx} {fx} {fx}\n"

        cell_lines = ''
        for row in cell:
            cell_lines += f"  {row[0]: .6f}  {row[1]: .6f}  {row[2]: .6f}\n"

        sysname = f'slab_{tag}_L{L}_Relax'
        content = TEMPLATE.format(
            tag=tag, L=L, sysname=sysname, natoms=len(slab),
            atom_lines=atom_lines, cell_lines=cell_lines, fix_lines=fix_lines,
            kx=kgrid_for(Lx), ky=kgrid_for(Ly), kz=1,
        )
        outp = os.path.join(OUTDIR, sysname + '.dat')
        with open(outp, 'w') as f:
            f.write(content)
        n_fix = int(freeze.sum())
        print(f'{sysname}.dat  n={len(slab)}  fixed={n_fix}  kgrid=({kgrid_for(Lx)},{kgrid_for(Ly)},1)')


if __name__ == '__main__':
    main()
