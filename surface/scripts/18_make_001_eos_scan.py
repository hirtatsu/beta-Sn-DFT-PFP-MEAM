"""
18_make_001_eos_scan.py (local WSL)

Prepare a 5-point in-plane EOS scan for β-Sn (001) slab, keeping a=b enforced.
Scan a values: [5.82, 5.90, 5.97, 6.05, 6.12] Å  (±3% around DFT bulk 5.97)

For each a:
  - in-plane cell = a × a × Lz (Lz fixed at original ≈ 54.94 Å)
  - rescale atomic x,y positions by (a / 5.97), keep z unchanged (same layer thicknesses)
  - atom-only BFGS, no freezing, fmax 3e-4 Ha/Bohr

Creates supercomputer/slab_001_L8_EOS_aX.XX/ dirs with:
  - .dat input
  - jobA.sh (copy of template)

DATA.PATH: ../../DFT_DATA19 (as before)
"""
import math
import os
import re
import shutil
import numpy as np
from ase.io import read
from ase import Atoms

HOME = os.path.expanduser('~')
BASE = os.path.join(HOME, 'projects/beta-sn-elasticity/surface_gb/supercomputer')
JOB_TEMPLATE = os.path.join(HOME, 'jobA.sh')
BFGS_DIR = os.path.join(BASE, 'slab_001_L8_Relax')

A_REF = 5.970
A_VALUES = [5.82, 5.90, 5.97, 6.05, 6.12]


def kgrid_for_xy(L_A):
    return max(2, int(math.ceil(2 * math.pi / (0.15 * L_A))))


DAT_HEADER = """#
# β-Sn (001) slab L=8 atom-only BFGS at a=b={a_val:.3f} Å (EOS scan)
# Built from BFGS-relaxed slab; positions rescaled isotropically in-plane.
#
System.CurrrentDirectory         ./
System.Name                      slab_001_L8_a{tag}
DATA.PATH                        ../../DFT_DATA19
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


def write_jobscript(dst_path, input_dat):
    with open(JOB_TEMPLATE) as f:
        text = f.read()
    new = re.sub(
        r'(mpirun\s+-np\s+\S+\s+)\S+(\s+)\S+\.dat(\s*.*)',
        lambda m: m.group(1) + '../openmx' + m.group(2) + input_dat + m.group(3),
        text,
    )
    with open(dst_path, 'w') as f:
        f.write(new)
    os.chmod(dst_path, 0o755)


def main():
    # source: BFGS-relaxed slab_001
    src_traj = os.path.join(BFGS_DIR, 'slab_001_L8_Relax.traj')
    base_slab = read(src_traj)
    pos0 = base_slab.get_positions().copy()
    cell0 = base_slab.cell.array.copy()
    n = len(base_slab)
    a0 = float(cell0[0, 0])
    b0 = float(cell0[1, 1])
    Lz = float(cell0[2, 2])
    if not np.isclose(a0, b0, atol=1e-3):
        print(f'WARN: BFGS traj has a≠b: {a0} vs {b0}, using a={a0}')
    print(f'Reference BFGS slab: a={a0:.3f}, b={b0:.3f}, Lz={Lz:.3f}, n={n}')

    for a_val in A_VALUES:
        tag = f"{a_val:.2f}".replace('.', 'p')
        dirname = f'slab_001_L8_EOS_a{tag}'
        dirpath = os.path.join(BASE, dirname)
        os.makedirs(dirpath, exist_ok=True)

        # rescale x,y by (a_val / a0), keep z
        scale = a_val / a0
        new_pos = pos0.copy()
        new_pos[:, 0] *= scale
        new_pos[:, 1] *= scale
        # new cell
        new_cell = np.array([
            [a_val, 0.0, 0.0],
            [0.0, a_val, 0.0],
            [0.0, 0.0, Lz],
        ])

        atom_lines = ''
        for i, (sym, xyz) in enumerate(zip(base_slab.get_chemical_symbols(), new_pos), 1):
            atom_lines += (f"  {i:3d}   {sym}   "
                           f"{xyz[0]: .8f}  {xyz[1]: .8f}  {xyz[2]: .8f}   7.0  7.0\n")
        cell_lines = ''
        for row in new_cell:
            cell_lines += f"  {row[0]: .6f}  {row[1]: .6f}  {row[2]: .6f}\n"

        kx = kgrid_for_xy(a_val)
        ky = kgrid_for_xy(a_val)
        sysname = f'slab_001_L8_a{tag}'
        content = DAT_HEADER.format(
            a_val=a_val, tag=tag, natoms=n,
            atom_lines=atom_lines, cell_lines=cell_lines,
            kx=kx, ky=ky, kz=1,
        )
        with open(os.path.join(dirpath, sysname + '.dat'), 'w') as f:
            f.write(content)
        write_jobscript(os.path.join(dirpath, 'jobA.sh'), sysname + '.dat')
        print(f'{dirname}/  a={a_val:.3f}  k={kx}x{ky}x1  jobA.sh')


if __name__ == '__main__':
    main()
