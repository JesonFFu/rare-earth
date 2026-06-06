# v6 3D cylindrical-hole photonic-crystal random-particle RE summary

All files were generated inside `E:\xitugrating\sin451_qbic_3d_design_v6_cyl_phc_random_particles`; previous versions were not edited.

## Model choice

- Used a 3D square-lattice cylindrical air-hole SiN slab unit cell because cylindrical PhC BIC leakage is inherently 3D.
- SiO2 was truncated to `500 nm` plus PML to fit the 16 GB memory budget; this is an existence check, not the final full-substrate fabrication model.
- Rare-earth particles are deterministic random placements inside the periodic representative cell; they are not placed on a regular grid.
- Added a refinement set anchored by `bic0413.mph` scale (`period=356 nm`, `r/period=0.220`) as a range check only.

## Scan result

- Mode rows simulated: `196`.
- Strict `|lambda-451| <= 0.5 nm`, `Q >= 1e9`: `0`.
- Relaxed `|lambda-451| <= 1.0 nm`, `Q >= 1e8`: `0`.

## Best candidate

| item | value |
|---|---:|
| P (nm) | `250` |
| cylinder radius (nm) | `50` |
| axis gap (nm) | `150` |
| lambda451 (nm) | `449.683097` |
| Q451 | `6.98417e+06` |
| particle count | `11` |
| top / bottom / side count | `5 / 4 / 2` |

## Diagnostic notes

- Highest-Q solved mode: `fine_P250_G140`, lambda `443.271 nm`, Q `2.36209e+07`.
- Best mode within `451+/-0.5 nm`: `cyl_P260_G150`, lambda `450.701 nm`, Q `2.89232e+06`.
- Best `bic0413.mph`-anchored scale row: `bicref_P370_rho0p20`, lambda `449.354 nm`, Q `31505.8`.
- The random particulate RE breaks the clean in-plane symmetry enough that the cylindrical-hole PhC branch does not retain the 1D model's ultra-high-Q behavior in this memory-limited 3D scan.

Conclusion: this 3D cylindrical-hole random-particle scan did not meet the high-Q feasibility threshold; inspect the candidates and consider larger-radius/period refinements or a finite-array design.

## Band sanity check

- k=(0, 0) pi/P: lambda `449.683 nm`, Q `6.98417e+06`.
- k=(0.02, 0) pi/P: lambda `449.481 nm`, Q `2.81916e+06`.
- k=(0.05, 0) pi/P: lambda `448.438 nm`, Q `638338`.
- k=(0.1, 0) pi/P: lambda `444.945 nm`, Q `131468`.
- k=(0.2, 0) pi/P: lambda `446.842 nm`, Q `3470.8`.
- k=(0.5, 0) pi/P: lambda `452.642 nm`, Q `1258.83`.
- k=(0, 0.5) pi/P: lambda `452.739 nm`, Q `1281.57`.
- k=(0.5, 0.5) pi/P: lambda `448.824 nm`, Q `534.165`.

## Output files

- `sin451_qbic_v6_best_3d_cyl_hole_phc_random_particles.mph`
- `coarse_scan_451_3d_particle.csv`
- `best_candidates_451_3d_particle.csv`
- `band_451_3d_particle.csv`
- `random_particle_layout.csv`
- `scan_q_vs_wavelength.svg`
- `band_q_wavelength.svg`
