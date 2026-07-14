# Comparison of Elastic Constants and Surface Energies of β-Sn

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21357962.svg)](https://doi.org/10.5281/zenodo.21357962)

Companion data and code repository for:

> H. Tatsumi, A. M. Ito, A. Takayama, H. Nishikawa.
> *Comparison of Elastic Constants and Surface Energies of β-Sn from Density Functional Theory, Universal Machine Learning Potential, and Empirical Potentials*.
> Submitted to *Modelling and Simulation in Materials Science and Engineering*
> (2026, in review).

This repository contains every input, output, and analysis script required to
reproduce the elastic constants, surface energies, and Wulff shapes presented
in the paper for **8 computational methods** applied to β-Sn (tetragonal,
I4₁/amd):

- **DFT/PBE** with OpenMX 3.9.9 (norm-conserving pseudopotentials, PBE19 basis)
- **PFP** (Preferred Potential v8) on Matlantis in 4 modes:
  - PFP/PBE
  - PFP/PBE+D3
  - PFP/r²SCAN
  - PFP/r²SCAN+D3
- **MEAM** (LAMMPS) with three classical potentials:
  - MEAM (Ravelo & Baskes 1997)
  - MEAM (Etesami et al. 2018)
  - MEAM (Ko et al. 2018)

## Citation

```
H. Tatsumi, A. M. Ito, A. Takayama, H. Nishikawa.
"Comparison of Elastic Constants and Surface Energies of β-Sn from Density Functional Theory, Universal Machine Learning Potential, and Empirical Potentials".
Modelling and Simulation in Materials Science and Engineering (2026, in review).
```

Archived at Zenodo — this version (v1.0.0, as reviewed): [doi:10.5281/zenodo.21357963](https://doi.org/10.5281/zenodo.21357963); all versions: [doi:10.5281/zenodo.21357962](https://doi.org/10.5281/zenodo.21357962).

## What's included

- **`elastic/`** — full elastic constant tensor calculations
  - `data/intel_results/` — DFT/PBE Cᵢⱼ run with all 6 Voigt strain trajectories
  - `data/ko2018_meam/` — MEAM (Ko 2018) reference + per-component error
  - `data/cij_*.csv` — consolidated Cᵢⱼ table and MAPE ranking vs experiment
  - `scripts/elastic_betasn.py` — reference implementation of the
    paper-matching elastic protocol (LBFGS, ε=±0.5%, internal relax,
    tetragonal 4/mmm symmetrization)
- **`surface/`** — surface energy calculations for 5 low-index faces
  - `data/surface_energies_5faces.{json,csv}` — master γ table (8 methods × 5
    faces; lowest-energy states, DFT (001)/(110) at L = 12)
  - `data/surface_energies_multistart.{json,csv}` — ideal-cleavage vs.
    lowest-energy surface states (SI Table S1)
  - `data/dft_openmx/`, `data/pfp_matlantis/`, `data/meam_lammps/`,
    `data/ko2018_meam/` — per-method raw data
  - `supercomputer/` — DFT/PBE OpenMX raw outputs (Plasma Simulator, NIFS):
    slab inputs, full SCF/MD logs, per-step `.ene`/`.md`/`.out`, including the
    L = 12/16 thickness runs and the seeded re-relaxations used in the revision
  - `make_figures.py` — regenerates the surface-energy figure previews
    (data of manuscript Figs. 4–5 and per-method Wulff shapes)
  - `scripts/27_wulff_shapes.py`, `30_surface_comparison_5faces_8methods.py`,
    `29_final_comparison_8methods.py`
- **`meam_potentials/`** — Sn.MEAM parameter files for the three classical
  potentials, exactly as used with LAMMPS in `data/meam_lammps/` and
  `data/ko2018_meam/`
- **`figures/`** — manuscript figures (150-dpi PNG previews and lossless
  300-dpi TIFs)
- **`beta-Sn.cif`** — reference crystal structure used as starting point,
  constructed from the published parameters corresponding to ICSD 40037
  (a = 5.831 Å, c = 3.182 Å); not a database export

## Headline results (from the paper)

### Bulk lattice constants (Å)
| | a | c | c/a |
|---|---:|---:|---:|
| Experiment (ICSD 40037) | 5.831 | 3.182 | 0.546 |
| DFT/PBE | 5.970 | 3.218 | 0.539 |
| PFP/PBE | 5.929 | 3.201 | 0.540 |
| PFP/PBE+D3 | 5.846 | 3.173 | 0.543 |

### Elastic constants — MAPE vs experiment (lower = better)
| Method | C₁₁ | C₃₃ | C₁₂ | C₁₃ | C₄₄ | C₆₆ | MAPE |
|---|---:|---:|---:|---:|---:|---:|---:|
| Experiment (Rayne–Chandrasekhar 1960) | 72.3 | 88.4 | 59.4 | 35.8 | 22.0 | 24.0 | — |
| **PFP/PBE+D3** | 98.5 | 121.1 | 36.2 | 36.5 | **22.9** | 16.2 | **25.1 %** ★ |
| DFT/PBE | 89.7 | 91.8 | 17.4 | 31.5 | 17.9 | 17.6 | 26.0 % |
| MEAM/Ko 2018 | 89.7 | 93.7 | 46.7 | 36.9 | 7.9 | 10.6 | 29.1 % |
| PFP/PBE | 114.0 | 104.6 | 41.5 | 41.2 | 29.5 | 31.9 | 31.4 % |

### Surface energies (mJ/m², lowest-energy states)
DFT/PBE predicts (100) to be the lowest-energy plane (492.2), and six of the
seven non-DFT methods reproduce this. Among them, PFP/r²SCAN agrees most
closely with DFT/PBE (MAPE 5.0%). All slab relaxations include
symmetry-broken starting configurations so that lower-lying surface states
are not missed; the ideal-cleavage vs. lowest-energy values for all 40
method–surface combinations are in
`surface/data/surface_energies_multistart.{json,csv}` (SI Table S1).

See `figures/fig5_surface_energies_8methods.png` and
`figures/supp_wulff_8methods.png` for the full picture.

## Reproducing the work

> **Note on paths.** Some legacy scripts in `elastic/scripts/`, `surface/scripts/`,
> and `surface/data/{ko2018_meam,…}/` contain hardcoded absolute paths from the
> original development environment. They are kept verbatim as a record of how
> the data were produced. To re-run on your own machine, either edit those paths
> to match your workspace or set the `LAMMPS_BIN`, `POT_DIR`, etc. variables at
> the top of each script. The numerical results they produced are already cached
> in JSON/CSV form under `*/data/`.

### DFT/PBE (OpenMX 3.9.9 on supercomputer)
- `surface/supercomputer/<run>/jobA.sh` — PBS job script (modify for your queue)
- `surface/supercomputer/<run>/<run>.dat` — OpenMX input file
- `meam_potentials/` — for LAMMPS comparison runs
- Reference structure: `beta-Sn.cif`

### PFP (Matlantis)
- `surface/data/pfp_matlantis/` — per-mode notebooks and energy tables
- Requires Matlantis (Preferred Networks) account

### MEAM (LAMMPS)
- `meam_potentials/Sn.{Ravelo1997, Etesami2018, Ko2018}.meam` (and matching `library.*.meam` files)
- LAMMPS scripts in `surface/data/{meam_lammps,ko2018_meam}/`

### Re-fitting the figures
```bash
python surface/make_figures.py
```

## License

Code: MIT.
Data (CIFs, JSON results, raw OpenMX/LAMMPS outputs): CC-BY-4.0.

**Exception:** the MEAM parameter files under `meam_potentials/` (and the
copies under `surface/data/meam_lammps/meam_potentials/`) reproduce
third-party published parameterizations and are *not* covered by the licenses
above — see `meam_potentials/NOTICE` for provenance and citation requirements.

## Contact

For questions: tatsumi.jwri@osaka-u.ac.jp
