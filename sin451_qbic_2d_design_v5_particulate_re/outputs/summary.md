# v5 451 nm particulate rare-earth quasi-BIC summary

All files were generated inside `E:\xitugrating\sin451_qbic_2d_design_v5_particulate_re`; previous versions and original files were not edited.

## Scan result

- Coarse mode rows: `708`.
- Fine mode rows: `3648`.
- Candidates meeting `|lambda-451| <= 0.5 nm` and `Q >= 1e9`: `2`.
- Candidates meeting `|lambda-451| <= 1.0 nm` and `Q >= 1e8`: `2`.

## Best particulate-RE candidate

| item | value |
|---|---:|
| P (nm) | `294` |
| W (nm) | `159` |
| etched gap (nm) | `135` |
| lambda451 (nm) | `450.798604` |
| Q451 | `1.08366e+13` |
| mode index | `1` |
| particle count | `19` |
| top coverage | `0.7862` |
| bottom coverage | `0.4854` |
| sidewall coverage | `0.5` |
| particle RE field-overlap proxy | `0.189221` |

Conclusion: the particulate-RE 451 nm design meets the preferred `Q >= 1e9` target.

## RE style comparison

- Continuous RE at same P/W: `lambda=456.688 nm`, `Q=8.60008e+11`.
- No RE at same P/W: `lambda=444.38 nm`, `Q=1.7334e+10`.
- This comparison helps separate particle scattering/loading from the underlying SiN grating radiation loss.

## Output files

- `sin451_qbic_v5_best_particulate_re.mph`
- `coarse_scan_451_particle.csv`, `fine_scan_451_particle.csv`, `best_candidates_451_particle.csv`
- `band_451_particle.csv`, `band_451_particle_wavelength.svg`, `band_451_particle_q.svg`
- `field_451_particle_ez2.csv`, `field_451_particle_ez2.svg`
- `particle_geometry_check.csv`, `compare_re_styles.csv`