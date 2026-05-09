"""
07_build_supercomputer_tree.py (run locally on WSL)

Builds a per-job directory tree under supercomputer/:
  supercomputer/<sysname>/
      <sysname>.dat
      jobA.sh

 - OpenMX .dat: DATA.PATH = ../../DFT_DATA19
 - jobA.sh: copied from ~/jobA.sh, mpirun line rewritten to the correct input file
"""
import json
import math
import os
import re
import shutil
import numpy as np
from ase import Atoms
from ase.build import surface as ase_surface
from ase.io import write

HOME = os.path.expanduser('~')
BASE = os.path.join(HOME, 'projects/beta-sn-elasticity/surface_gb/supercomputer')
JOB_TEMPLATE = os.path.join(HOME, 'jobA.sh')

DATA_PATH_REL = '../../DFT_DATA19'
DFT_LATTICE = {'a': 5.970, 'c': 3.218}

JOBS = [
    {'sysname': 'bulk_dft_SP', 'kind': 'bulk',      'layers': None, 'face': None},
    {'sysname': 'slab_001_L8_Relax',  'kind': 'slab', 'layers': 8,  'face': '001'},
    {'sysname': 'slab_110_L8_Relax',  'kind': 'slab', 'layers': 8,  'face': '110'},
    {'sysname': 'slab_100_L12_Relax', 'kind': 'slab', 'layers': 12, 'face': '100'},
    {'sysname': 'slab_101_L12_Relax', 'kind': 'slab', 'layers': 12, 'face': '101'},
]
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
# {title}
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

"""

MD_SP = """MD.Type                          Nomd
MD.maxIter                       1
"""

MD_RELAX = """MD.Type                          BFGS
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


def write_openmx(path, sysname, title, atoms, kgrid, md_block):
    content = HEADER.format(
        title=title, sysname=sysname, datapath=DATA_PATH_REL,
        natoms=len(atoms),
        atom_lines=fmt_atoms(atoms), cell_lines=fmt_cell(atoms),
        kx=kgrid[0], ky=kgrid[1], kz=kgrid[2],
    ) + md_block
    with open(path, 'w') as f:
        f.write(content)


def write_jobscript(dst_path, input_dat):
    with open(JOB_TEMPLATE) as f:
        text = f.read()
    # rewrite mpirun line's input file (preserves flags / core count)
    new = re.sub(
        r'(mpirun\s+-np\s+\S+\s+\./openmx\s+)(\S+\.dat)(\s*.*)',
        lambda m: m.group(1) + input_dat + m.group(3),
        text,
    )
    with open(dst_path, 'w') as f:
        f.write(new)
    os.chmod(dst_path, 0o755)


def main():
    bulk = build_bulk()
    manifest = {
        'data_path_relative': DATA_PATH_REL,
        'bulk_lattice': DFT_LATTICE,
        'jobs': [],
    }

    for j in JOBS:
        sysname = j['sysname']
        dirpath = os.path.join(BASE, sysname)
        os.makedirs(dirpath, exist_ok=True)
        dat_path = os.path.join(dirpath, sysname + '.dat')

        if j['kind'] == 'bulk':
            kg = (kgrid_for(DFT_LATTICE['a']),
                  kgrid_for(DFT_LATTICE['a']),
                  kgrid_for(DFT_LATTICE['c']))
            write_openmx(
                dat_path, sysname,
                f'β-Sn bulk single-point at a={DFT_LATTICE["a"]:.3f}, c={DFT_LATTICE["c"]:.3f}',
                bulk, kg, MD_SP,
            )
            info = {'sysname': sysname, 'kind': 'bulk', 'n_atoms': 4,
                    'kgrid': list(kg)}
        else:
            miller = tuple(int(c) for c in j['face'])
            slab = ase_surface(bulk, miller, layers=j['layers'],
                               vacuum=VACUUM, periodic=True)
            slab.center(vacuum=VACUUM, axis=2)
            Lx = float(np.linalg.norm(slab.cell.array[0]))
            Ly = float(np.linalg.norm(slab.cell.array[1]))
            A = float(np.linalg.norm(np.cross(slab.cell.array[0],
                                              slab.cell.array[1])))
            kg = (kgrid_for(Lx), kgrid_for(Ly), 1)
            write_openmx(
                dat_path, sysname,
                f'β-Sn ({j["face"]}) slab L={j["layers"]} full-atom relaxation',
                slab, kg, MD_RELAX,
            )
            write(os.path.join(dirpath, sysname + '.traj'), slab)
            info = {'sysname': sysname, 'kind': 'slab', 'face': j['face'],
                    'layers': j['layers'], 'n_atoms': len(slab),
                    'cell_Ang': [Lx, Ly, float(np.linalg.norm(slab.cell.array[2]))],
                    'area_A2': A, 'kgrid': list(kg)}

        # job script
        job_path = os.path.join(dirpath, 'jobA.sh')
        write_jobscript(job_path, sysname + '.dat')
        manifest['jobs'].append(info)
        print(f'{sysname}/  ->  dat + jobA.sh (kgrid={kg})')

    with open(os.path.join(BASE, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f'\nwrote manifest -> {BASE}/manifest.json')


if __name__ == '__main__':
    main()
