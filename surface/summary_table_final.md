# β-Sn Surface & Wulff Shape — Final Summary (7 faces × 7 methods)

## 1. γ_surf table (mJ/m²)

| face | DFT/PBE BFGS | PFP/PBE | PFP/PBE+D3 | PFP/R²SCAN | PFP/R²SCAN+D3 | MEAM/Ravelo 1997 | MEAM/Etesami 2018 |
|------|-------------:|--------:|-----------:|-----------:|--------------:|-----------------:|------------------:|
| (100) | **492.2**  | 385.5  | 642.2  | 784.6 | 614.8  | 735.0  | 1376.3 |
| (101) | **505.6**  | 412.7  | 832.5  | 476.6 | 707.3  | 818.3  | 1354.7 |
| (110) | **526.9**  | 487.0  | 737.5  | 599.3 | 730.6  | 954.2  | 1587.4 |
| (001) | **550.5**  | 460.1  | 734.4  | 527.5 | 680.4  | 896.8  | 1299.7 |
| (111) | **pending** | 417.1 | 685.1  | 538.8 | 690.1  | 909.3  | 1384.7 |
| (211) | **pending** | 424.0 | 704.1  | 529.2 | 680.9  | 941.4  | 1224.8 |
| (112) | **pending** | 442.5 | 668.2  | 483.7 | 633.9  | 814.9  | 1303.3 |

Reference: experimental β-Sn γ_avg ≈ 670 mJ/m² [Tyson–Miller 1977]. DFT (001), (110), (100), (101) confirmed by OptC5/EOS checks; (111), (211), (112) jobs prepared in `supercomputer/slab_{111,211,112}_L12_Relax/`.

## 2. Face ordering (low → high γ)

| method | ordering |
|---|---|
| DFT/PBE  | (100) < (101) < (110) < (001) |
| PFP/PBE  | (100) < (101) < (001) < (110) ≈ (111) ≈ (211) ≈ (112) |
| PFP/PBE+D3 | (100) < (112) < (111) < (001) ≈ (110) < (211) < (101) |
| PFP/R²SCAN | (101) < (112) < (001) < (211) ≈ (111) < (110) < (100) |
| PFP/R²SCAN+D3 | (100) < (112) < (001) < (211) < (111) < (101) < (110) |
| MEAM/Ravelo 1997 | (100) < (112) < (101) < (001) < (111) < (211) < (110) |
| MEAM/Etesami 2018 | (211) < (001) < (112) < (101) < (100) < (111) < (110) |

## 3. MAE vs DFT/PBE (mJ/m², over 4 confirmed faces)

| rank | method | MAE |
|---:|---|---:|
| 1 | PFP/PBE | 82 |
| 2 | PFP/R²SCAN | 104 |
| 3 | PFP/R²SCAN+D3 | 164 |
| 4 | PFP/PBE+D3 | 218 |
| 5 | MEAM/Ravelo 1997 | 332 |
| 6 | MEAM/Etesami 2018 | 886 |

## 4. Wulff shape files

- `figures/wulff_all_methods.png` (+ `.tif`) — 7 手法の Wulff 形状グリッド
- DFT は (111), (211), (112) 未完でシンプルな多面体、他は 7 面から構築

## 5. Computational settings

- β-Sn bulk (DFT reference): a = 5.970 Å, c = 3.218 Å, I4₁/amd
- Slab: 15 Å vacuum, all atoms free
- DFT: OpenMX/PBE BFGS, SCF 1e-7, Δk=0.15 rad/Å
- PFP: Matlantis v8, fmax = 0.015 eV/Å
- MEAM: Ravelo-Baskes 1997 (Vella 2017 Table I 再構築) & Etesami 2018 (intel supplied)
- LAMMPS: `pair_style meam`, reference fcc, Sn 'fcc' 12 50 118.710
