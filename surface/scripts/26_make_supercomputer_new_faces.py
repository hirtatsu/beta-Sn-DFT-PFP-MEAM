"""Prepare supercomputer/<sysname>/ dirs for new (111), (211), (112) slabs.
Same Acta-compliant settings as previous Relax jobs; all atoms free; BFGS; walltime 12h.
"""
import json
import math
import os
import re
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.io import write

HOME = os.path.expanduser('~')
BASE = os.path.join(HOME, 'projects/beta-sn-elasticity/surface_gb/supercomputer')
JOB_TEMPLATE = os.path.join(HOME, 'jobA.sh')

DATA_PATH_REL = '../../DFT_DATA19'
DFT_LATTICE = {'a': 5.970, 'c': 3.218}
NEW_FACES = [('111', 12), ('211', 12), ('112', 12)]
VACUUM = 15.0


def kgrid_for(L_A):
    return max(2, int(math.ceil(2 * math.pi / (0.15 * L_A))))


def build_bulk():
    a, c = DFT_LATTICE['a'], DFT_LATTICE['c']
    cell = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    scaled = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
              [0.0, 0.5, 0.25], [0.5, 0.0, 0.75]]
    return Atoms('Sn4', scaled_positions=scaled, cell=cell, pbc=True)


HEADER = """#
# β-Sn ({face}) slab L={L} full-atom relaxation at DFT lattice
# k-mesh target Δk ≈ 0.15 rad/Å (Acta Materialia 2026 compliant)
#
System.CurrrentDirectory         ./
System.Name                      {sysname}
DATA.PATH                        {datapath}
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
MD.maxIter                       1000
MD.Opt.criterion                 3.0e-4
"""


def fmt_atoms(atoms):
    s = ''
    for i, (sym, xyz) in enumerate(zip(atoms.get_chemical_symbols(),
                                        atoms.get_positions()), 1):
        s += f"  {i:3d}   {sym}   {xyz[0]: .8f}  {xyz[1]: .8f}  {xyz[2]: .8f}   7.0  7.0\n"
    return s


def fmt_cell(atoms):
    s = ''
    for row in atoms.cell.array:
        s += f"  {row[0]: .6f}  {row[1]: .6f}  {row[2]: .6f}\n"
    return s


def write_jobscript(dst_path, input_dat):
    with open(JOB_TEMPLATE) as f:
        text = f.read()
    # rewrite mpirun binary path and input
    text = re.sub(
        r'(mpirun\s+-np\s+\S+\s+)\S+(\s+)\S+\.dat(\s*.*)',
        lambda m: m.group(1) + '../openmx' + m.group(2) + input_dat + m.group(3),
        text,
    )
    # extend walltime to 12 h
    text = re.sub(r'#PBS -l walltime=\S+', '#PBS -l walltime=12:00:00', text)
    with open(dst_path, 'w') as f:
        f.write(text)
    os.chmod(dst_path, 0o755)


def main():
    bulk = build_bulk()
    for face, L in NEW_FACES:
        miller = tuple(int(x) for x in face)
        slab = ase_surface(bulk, miller, layers=L, vacuum=VACUUM, periodic=True)
        slab.center(vacuum=VACUUM, axis=2)
        Lx = float(np.linalg.norm(slab.cell.array[0]))
        Ly = float(np.linalg.norm(slab.cell.array[1]))
        A = float(np.linalg.norm(np.cross(slab.cell.array[0], slab.cell.array[1])))
        kg = (kgrid_for(Lx), kgrid_for(Ly), 1)
        sysname = f'slab_{face}_L{L}_Relax'
        dirpath = os.path.join(BASE, sysname)
        os.makedirs(dirpath, exist_ok=True)

        content = HEADER.format(
            face=face, L=L, sysname=sysname, datapath=DATA_PATH_REL,
            natoms=len(slab),
            atom_lines=fmt_atoms(slab), cell_lines=fmt_cell(slab),
            kx=kg[0], ky=kg[1], kz=kg[2])
        with open(os.path.join(dirpath, sysname + '.dat'), 'w') as f:
            f.write(content)
        write(os.path.join(dirpath, sysname + '.traj'), slab)
        write_jobscript(os.path.join(dirpath, 'jobA.sh'), sysname + '.dat')

        print(f'{sysname}/  n={len(slab)}  cell=({Lx:.2f},{Ly:.2f})  A={A:.2f}  kgrid={kg}')


if __name__ == '__main__':
    main()
