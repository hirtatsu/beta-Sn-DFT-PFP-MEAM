"""Generate a dense-k bulk SP reference to match new (112) slabs.
k-grid: Δk ≈ 0.10 rad/Å  → 11×11×20 for a=5.970, c=3.218.
"""
import math, os, re

HOME = os.path.expanduser('~')
BASE = os.path.join(HOME, 'projects/beta-sn-elasticity/surface_gb/supercomputer')
JOB_TEMPLATE = os.path.join(HOME, 'jobA.sh')
A, C = 5.970, 3.218
DATA_PATH_REL = '../../DFT_DATA19'


def k_dense(L):
    return max(2, int(math.ceil(2 * math.pi / (0.10 * L))))


kx = ky = k_dense(A); kz = k_dense(C)

content = f"""#
# β-Sn bulk SP at DFT lattice, dense k (Δk≈0.10 rad/Å)
# Reference for (112) verification slabs with same k-density.
#
System.CurrrentDirectory         ./
System.Name                      bulk_dft_denseK_SP
DATA.PATH                        {DATA_PATH_REL}
level.of.stdout                  1
level.of.fileout                 0

Species.Number                   1
<Definition.of.Atomic.Species
 Sn  Sn7.0-s2p2d3f1  Sn_PBE19
Definition.of.Atomic.Species>

Atoms.Number                     4
Atoms.SpeciesAndCoordinates.Unit Ang
<Atoms.SpeciesAndCoordinates
    1   Sn    0.00000000  0.00000000  0.00000000  7.0  7.0
    2   Sn    2.98500000  2.98500000  1.60900000  7.0  7.0
    3   Sn    0.00000000  2.98500000  0.80450000  7.0  7.0
    4   Sn    2.98500000  0.00000000  2.41350000  7.0  7.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit           Ang
<Atoms.UnitVectors
  5.970000  0.000000  0.000000
  0.000000  5.970000  0.000000
  0.000000  0.000000  3.218000
Atoms.UnitVectors>

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

sysname = 'bulk_dft_denseK_SP'
dirpath = os.path.join(BASE, sysname)
os.makedirs(dirpath, exist_ok=True)
with open(os.path.join(dirpath, sysname + '.dat'), 'w') as f:
    f.write(content)

# job script
with open(JOB_TEMPLATE) as f:
    t = f.read()
t = re.sub(r'(mpirun\s+-np\s+\S+\s+)\S+(\s+)\S+\.dat(\s*.*)',
           lambda m: m.group(1) + '../openmx' + m.group(2) + sysname + '.dat' + m.group(3), t)
t = re.sub(r'#PBS -l walltime=\S+', '#PBS -l walltime=01:00:00', t)
with open(os.path.join(dirpath, 'jobA.sh'), 'w') as f:
    f.write(t)
os.chmod(os.path.join(dirpath, 'jobA.sh'), 0o755)
print(f'wrote {dirpath}/  (k={kx}×{ky}×{kz}, walltime 1h)')
