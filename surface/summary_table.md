# β-Sn Surface & Grain Boundary Energy — Summary

## 1. Surface energies (mJ/m²)

| face | slab layers | γ_PFP (PBE+D3) | γ_DFT (OpenMX/PBE SP at PFP geom) | Δ (DFT–PFP) |
|------|-------------|---------------:|----------------------------------:|------------:|
| (001) | 8 | 736.2 | 579.6 | -156.6 |
| (110) | 8 | 737.5 | 532.5 | -205.0 |
| (100) | 12 | 834.5 | 619.7 | -214.8 |
| (101) | 12 | 833.1 | 624.1 | -209.0 |

- Experimental β-Sn average: ~670 mJ/m² (Tyson–Miller 1977, semi-empirical avg for β-Sn).
- DFT SP consistently under-predicts vs PFP+D3 (by ~20–30%).
- DFT fully relaxed values pending on supercomputer jobs.

## 2. Grain boundary energies (PFP/PBE+D3, mJ/m²)

| GB | type | θ (°) | γ_GB | method |
|----|------|------:|-----:|--------|
| Σ5(210)[001] | tilt_001 | 36.87 | 181.2 | RBT 4-pt scan → min |
| Σ17(410)[001] | tilt_001 | 28.07 | 300.4 | RBT 4-pt scan → min |
| Σ13(320)[001] | tilt_001 | 22.62 | 373.7 | RBT 4-pt scan → min |
| Σ5(310)[001] | tilt_001 | 53.13 | 369.7 | RBT 4-pt scan → min |
| (101) twin | twin | — | 126.4 | γ-surface 10×10 → physical min |
| (301) twin | twin | — | 21.6 | γ-surface 10×10 → physical min |

Notes:
- [001] symmetric tilt convergence checked on Σ5(210) (N=3 gives γ within ±5 mJ/m² of N=8).
- Compound twins constructed via mirror; RBT space includes many "bulk-recovery" points (γ<20 mJ/m²) due to BCT body-center symmetry — only "physical" minima reported.
- (301) γ_GB = 21.6 mJ/m² is at the lower tail: consistent with β-Sn being the classic deformation-twin material ("tin cry").
- (012)[100] tilt skipped: a≠c tetragonality precludes orthogonal CSL; will require Type-II treatment.

## 3. Method settings

- **Engine**: PFP v8 (calc_mode `PBE_PLUS_D3`) on Matlantis
- **Bulk lattice** (PFP-optimized): a = 5.8457 Å, c = 3.1732 Å
- **Surface slabs**: 4 orientations, 15 Å vacuum, middle 30% atoms frozen during PFP relaxation
- **GB construction**: ase.build.make_supercell with CSL basis; symmetric tilt via mirror + RBT; dedup tol 0.5 Å (MIC)
- **Relaxation**: LBFGS, fmax 0.05 eV/Å, 200–300 step cap
- **γ formula**: γ = (E_cell − N·E_bulk/atom) / (2·A)