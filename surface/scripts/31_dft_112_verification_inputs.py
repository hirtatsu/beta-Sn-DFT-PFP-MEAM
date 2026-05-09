"""Generate verification DFT inputs for (112):
  - use pymatgen SlabGenerator symmetric slab
  - build on DFT-optimized bulk a=5.970, c=3.218
  - try two thicknesses: 24 Å and 36 Å
  - dense k-grid (Δk ≈ 0.10 rad/Å)
"""
import math, os, re, numpy as np
from ase.io import write
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

HOME = os.path.expanduser('~')
BASE = os.path.join(HOME, 'projects/beta-sn-elasticity/surface_gb/supercomputer')
JOB_TEMPLATE = os.path.join(HOME, 'jobA.sh')

A, C = 5.970, 3.218
DATA_PATH_REL = '../../DFT_DATA19'
VACUUM = 15.0


def kgrid_dense(L_A):
    # Δk = 0.10 rad/Å (finer than prior 0.15)
    return max(2, int(math.ceil(2 * math.pi / (0.10 * L_A))))


def build_bulk_pmg():
    latt = Lattice.from_parameters(A, A, C, 90, 90, 90)
    return Structure(latt, ['Sn']*4,
                     [[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]])


HEADER = """#
# {title}
#   pymatgen symmetric slab, Δk ≈ 0.10 rad/Å (denser than previous)
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
    return ''.join(f"  {i+1:3d}   {sym}   {p[0]: .8f}  {p[1]: .8f}  {p[2]: .8f}   7.0  7.0\n"
                   for i, (sym, p) in enumerate(zip(atoms.get_chemical_symbols(), atoms.get_positions())))


def fmt_cell(atoms):
    return ''.join(f"  {r[0]: .6f}  {r[1]: .6f}  {r[2]: .6f}\n" for r in atoms.cell.array)


def write_jobA(dst, input_dat, walltime='18:00:00'):
    with open(JOB_TEMPLATE) as f:
        text = f.read()
    text = re.sub(r'(mpirun\s+-np\s+\S+\s+)\S+(\s+)\S+\.dat(\s*.*)',
                  lambda m: m.group(1) + '../openmx' + m.group(2) + input_dat + m.group(3), text)
    text = re.sub(r'#PBS -l walltime=\S+', f'#PBS -l walltime={walltime}', text)
    with open(dst, 'w') as f:
        f.write(text)
    os.chmod(dst, 0o755)


def main():
    bulk = build_bulk_pmg()
    ada = AseAtomsAdaptor()

    for min_slab in [24, 36]:
        sg = SlabGenerator(bulk, (1,1,2), min_slab_size=min_slab, min_vacuum_size=VACUUM,
                           lll_reduce=False, center_slab=True, in_unit_planes=False,
                           primitive=True, reorient_lattice=True)
        slabs = sg.get_slabs(symmetrize=True, tol=0.01)
        assert len(slabs) == 1
        atoms = ada.get_atoms(slabs[0])
        Lx = float(np.linalg.norm(atoms.cell.array[0]))
        Ly = float(np.linalg.norm(atoms.cell.array[1]))
        A_area = float(np.linalg.norm(np.cross(atoms.cell.array[0], atoms.cell.array[1])))
        kg = (kgrid_dense(Lx), kgrid_dense(Ly), 1)

        sysname = f'slab_112_pmg_t{min_slab}_denseK'
        dirpath = os.path.join(BASE, sysname)
        os.makedirs(dirpath, exist_ok=True)
        content = HEADER.format(
            title=f'β-Sn (112) slab (pymatgen, thick={min_slab} Å, dense k)',
            sysname=sysname, datapath=DATA_PATH_REL, natoms=len(atoms),
            atom_lines=fmt_atoms(atoms), cell_lines=fmt_cell(atoms),
            kx=kg[0], ky=kg[1], kz=kg[2])
        with open(os.path.join(dirpath, sysname + '.dat'), 'w') as f:
            f.write(content)
        write(os.path.join(dirpath, sysname + '.traj'), atoms)
        write_jobA(os.path.join(dirpath, 'jobA.sh'), sysname + '.dat')
        print(f'{sysname}/  n={len(atoms)}  |a|={Lx:.2f} |b|={Ly:.2f} A={A_area:.2f}  kgrid={kg}')


if __name__ == '__main__':
    main()
