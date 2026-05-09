# β-Sn Surface Energy — Supercomputer Job Tree

## Directory layout

```
supercomputer/
├── README.md
├── manifest.json
├── bulk_dft_SP/
│   ├── bulk_dft_SP.dat
│   └── jobA.sh
├── slab_001_L8_Relax/
│   ├── slab_001_L8_Relax.dat        (32 atoms)
│   ├── slab_001_L8_Relax.traj       (initial geometry, ASE format)
│   └── jobA.sh
├── slab_110_L8_Relax/
│   ├── slab_110_L8_Relax.dat        (32 atoms)
│   ├── slab_110_L8_Relax.traj
│   └── jobA.sh
├── slab_100_L12_Relax/
│   ├── slab_100_L12_Relax.dat       (48 atoms)
│   ├── slab_100_L12_Relax.traj
│   └── jobA.sh
└── slab_101_L12_Relax/
    ├── slab_101_L12_Relax.dat       (48 atoms)
    ├── slab_101_L12_Relax.traj
    └── jobA.sh
```

Each subdirectory is self-contained: `cd` into it, `qsub jobA.sh`.

## Expected path layout on the supercomputer

Place the `openmx` executable at `supercomputer/openmx` (one level above each job directory). `jobA.sh` calls it via `../openmx`. The `DFT_DATA19` database sits one level above `supercomputer/`:

```
<workdir>/
├── DFT_DATA19/                ← OpenMX PAO + PP database
└── supercomputer/
    ├── openmx                 ← OpenMX executable
    ├── bulk_dft_SP/
    │   ├── bulk_dft_SP.dat
    │   └── jobA.sh            (runs ../openmx)
    ├── slab_001_L8_Relax/
    │   └── ...
    └── ...
```

From inside any `supercomputer/<sysname>/` directory:
- `../openmx` → the executable
- `../../DFT_DATA19` → the pseudopotential/basis database

## jobA.sh template

Copied from `~/jobA.sh`:
- PBS_P `NIFS26KISM025`
- Queue `A_S`
- 1 node × 256 cpus (128 MPI × 2 OpenMP threads)
- 752 GB memory
- **Walltime: 3 h** (may need extension for 48-atom slabs — see below)

Only the `mpirun` line has been rewritten to reference each directory's own `.dat` file.

## Walltime guidance

| job | atoms | suggested walltime |
|---|---:|---|
| bulk_dft_SP         | 4  | 1 h (default 3 h OK)    |
| slab_001_L8_Relax   | 32 | 6–12 h (extend `#PBS -l walltime=12:00:00`) |
| slab_110_L8_Relax   | 32 | 6–12 h |
| slab_100_L12_Relax  | 48 | 12–24 h |
| slab_101_L12_Relax  | 48 | 12–24 h |

The (100) surface historically needed the most SCF iterations per MD step (charge sloshing from short c-axis in plane). Consider queueing it first so it can start earlier if the queue is long.

## Computational setup (summary)

| item | value |
|---|---|
| Code | OpenMX v2019 |
| Basis | Sn7.0-s2p2d3f1 |
| Pseudopotential | Sn_PBE19 |
| XC functional | GGA-PBE |
| Spin polarization | off (β-Sn non-magnetic) |
| scf.energycutoff | 200 Ryd |
| ElectronicTemperature | 300 K |
| scf.criterion | 1.0 × 10⁻⁷ |
| Mixing | rmm-diisk, Init=0.005, Max=0.03, Kerker=10, History=15, StartPulay=60, maxIter=600 |
| k-mesh | Δk ≈ 0.15 rad/Å (uniform) |
| MD.Type (Relax) | BFGS |
| MD.maxIter (Relax) | **1000** |
| MD.Opt.criterion (Relax) | **3.0 × 10⁻⁴ Ha/Bohr** |
| Bulk lattice | a = 5.970 Å, c = 3.218 Å (from OpenMX elastic paper) |
| Vacuum (slabs) | 15 Å |
| Atom constraints | **none** (all atoms free) |

## Post-processing

After all jobs finish, from any working directory:
```bash
# collect final energies
for d in supercomputer/*/ ; do
  name=$(basename "$d")
  grep 'Utot' "$d/log.txt" | tail -1 | awk -v n="$name" '{printf "%-30s %s\n", n, $3}'
done
```

Surface energy (per face):
```
γ = (E_slab − N × E_bulk_per_atom) / (2 × A)
   where E_bulk_per_atom = Utot(bulk_dft_SP) / 4
         N = slab atom count
         A = in-plane area (see manifest.json for each slab)
   unit: 1 eV/Å² = 16021.77 mJ/m²
```

## Reference PFP/PBE+D3 values (for comparison)

Computed at PFP-optimized lattice (a = 5.8457, c = 3.1732 Å):

| face | γ_PFP (mJ/m²) |
|---|---:|
| (001) | 736.2 |
| (110) | 737.5 |
| (100) | 834.5 |
| (101) | 833.1 |
