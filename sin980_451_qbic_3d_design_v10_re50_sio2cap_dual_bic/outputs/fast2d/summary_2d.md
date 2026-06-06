# v10 Fast 2D Dual-BIC Screening

Directory: `E:\xitugrating\sin980_451_qbic_3d_design_v10_re50_sio2cap_dual_bic`. This 2D screening keeps the v10 stack and only patterns the 300 nm SiN layer.

## Best 2D Candidate

- Design: `fine_primitive_P562_f0p72`
- Family: `primitive_1d_grating_fine`
- Pitches: `562.000 nm`
- Widths: `404.640 nm`
- Minimum feature/gap: `157.360 nm`
- 980 nm branch: lambda `980.115156 nm`, Q `3.54213e+09`
- 451 nm branch: lambda `453.291596 nm`, Q `1300.06`
- Meets requested Q thresholds: `False`

## Family Bests

| family | design | lambda980 | Q980 | lambda451 | Q451 | min feature |
|---|---|---:|---:|---:|---:|---:|
| primitive_1d_grating_fine | fine_primitive_P562_f0p72 | 980.1152 | 3.542e+09 | 453.2916 | 1300 | 157.4 |
| primitive_1d_grating | primitive_P294_f0p42 | 990.1350 | 7.055 | 445.0267 | 3.358e+09 | 123.5 |

## 980 nm Q(k)

- k=0 pi/L: lambda `980.115156 nm`, Q `3.54213e+09`
- k=0.002 pi/L: lambda `980.119995 nm`, Q `1.36112e+06`
- k=0.006 pi/L: lambda `980.158697 nm`, Q `151462`
- k=0.02 pi/L: lambda `980.597852 nm`, Q `13811.2`

## 451 nm Q(k)

- k=0 pi/L: lambda `453.291596 nm`, Q `1300.06`
- k=0.002 pi/L: lambda `453.290919 nm`, Q `1295.59`
- k=0.006 pi/L: lambda `453.285502 nm`, Q `1260.93`
- k=0.02 pi/L: lambda `453.223872 nm`, Q `966.922`

## Files

- `fast2d_best_980.mph` and `fast2d_best_451.mph`: same geometry with target-dependent material constants.
- `coarse_scan_2d.csv`, `fine_scan_2d.csv`, `dual_best_candidates_2d.csv`.
- `band_980_2d.csv`, `band_451_2d.csv`, `q_near_gamma_980_2d.csv`, `q_near_gamma_451_2d.csv`.

## Interpretation

- Treat this as a fast existence screen. The earlier 3D coarse scan is retained separately in `outputs/coarse_dual_scan.csv`.
