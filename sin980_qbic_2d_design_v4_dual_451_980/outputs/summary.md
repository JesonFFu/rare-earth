# v4 980/451 nm dual-resonance quasi-BIC design summary

All v4 files were generated inside `E:\xitugrating\sin980_qbic_2d_design_v4_dual_451_980`; v1/v2/v3 and the original inputs were not edited.

## Targets

- 980 nm pump quasi-BIC: `|lambda-980| <= 0.5 nm`, `Q >= 1e+08`.
- 451 nm emission out-coupling resonance: `|lambda-451| <= 5.0 nm`, `Q >= 100`, preferably `1e3-1e5`.
- Fabrication filter: minimum feature `>= 100 nm`; remaining air gap after 16 nm side RE layers `>= 100 nm`.

## Scan status

- Periodic/supercell candidates simulated: `55`.
- Wide top-groove backup candidates simulated: `9`.
- Folded-only supercell candidates excluded from recommendation: `4`.
- True geometry candidates satisfying both thresholds: `0`.

## Recommended candidate

| item | value |
|---|---:|
| design id | `periodic_N1_baseline` |
| family | `periodic_supercell` |
| supercell periods | `1` |
| pitches (nm) | `564.300` |
| ridge widths (nm) | `383.724` |
| top groove width/depth (nm) | `0.000 / 0.000` |
| min feature (nm) | `180.576` |
| min etched gap (nm) | `180.576` |
| min remaining air gap after RE (nm) | `148.576` |
| lambda980 (nm) | `980.009896` |
| Q980 | `1.53252e+11` |
| lambda451 (nm) | `451.625021` |
| Q451 | `79.367` |
| 451 RE overlap proxy | `0.0497597` |

Saved models:

- `outputs/sin980_451_dual_best.mph` uses the recommended geometry solved at 980 nm.
- `outputs/sin980_451_dual_best_451check.mph` uses the same geometry solved at 451 nm for field inspection.

Note: this is the best fabrication-filtered true-geometry proof-of-existence candidate found in the current scan, but it does not fully satisfy every target threshold.

## Folded-band evidence, not recommended as a standalone geometry

- Best folded-only case: `periodic_N2_dP40_dW0`.
- It gives 980 nm `980.01 nm`, Q `1.90148e+10` and 451 nm `450.618 nm`, Q `153.326`.
- Because its supercell does not introduce a real geometric modulation, it is treated as evidence that a 451 nm branch exists near the desired wavelength, not as the final dual-resonance structure.

## Wide top-groove comparison

- Best backup: `wide_N2_w140_d100_dW40`.
- 980 nm: `979.38 nm`, Q `21.1772`.
- 451 nm: `451.424 nm`, Q `115.988`.
- This backup is not preferred unless it clearly beats the periodic/supercell route because rare-earth coverage in a top groove is less reliable than in the large primary gaps.

## Output files

- `periodic_supercell_scan.csv`
- `wide_top_groove_scan.csv`
- `dual_best_candidates.csv`
- `band_980.csv`, `band_451.csv`
- `field_980_field_ez2.svg`, `field_451_field_ez2.svg`
- `scan_q980.svg`, `scan_q451.svg`, `scan_lambda980.svg`, `scan_lambda451.svg`