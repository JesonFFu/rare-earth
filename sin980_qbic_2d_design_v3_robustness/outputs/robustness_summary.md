# V3 robustness scan: large-gap SiN grating with 16 nm RE layer

## Scan setup

- Nominal center: `P=564.3 nm`, `W=383.724 nm`, top nominal gap `180.576 nm`.
- Period perturbations: `-5,-2,0,2,5` nm.
- Ridge-width perturbations: `-10,-5,0,5,10` nm.
- Sidewall taper angles: `0,2,-2` deg; positive means bottom is wider than top.
- Each result tracks the nearest 980 nm Gamma-point eigenmode.

## Summary

- Simulated cases: `75`.
- Cases with `Q >= 1e9`: `59/75`.
- Cases with minimum etched gap `>=140 nm`: `75/75`.
- Worst Q in scan: `533.76` at `P=566.3 nm`, `W=373.724 nm`, `taper=2 deg`.
- Largest wavelength detuning: `-25.1843 nm`.

Grouped by sidewall angle:

| taper (deg) | cases | Q >= 1e9 | min Q | median Q | max \|lambda-980\| (nm) |
|---:|---:|---:|---:|---:|---:|
| 0 | 25 | 25 | 8.17028e10 | 1.65610e11 | 12.3078 |
| +2 | 25 | 12 | 533.76 | 547.264 | 25.1843 |
| -2 | 25 | 22 | 541.425 | 3.73622e11 | 17.4036 |

## Key cases

| case | P (nm) | W_top (nm) | fill | taper (deg) | top gap (nm) | bottom gap (nm) | min gap (nm) | remaining air gap (nm) | lambda (nm) | Q |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| nominal | 564.3 | 383.724 | 0.68 | 0 | 180.576 | 180.576 | 180.576 | 148.576 | 980.010144 | 1.6561e+11 |
| worst Q | 566.3 | 373.724 | 0.65994 | 2 | 192.576 | 171.624 | 171.624 | 139.624 | 1004.49185 | 533.76 |
| largest detuning | 569.3 | 373.724 | 0.656462 | 2 | 195.576 | 174.624 | 174.624 | 142.624 | 954.815703 | 7.1698e+10 |

## Best candidates in robustness grid

| rank | P (nm) | W_top (nm) | fill | taper (deg) | top gap (nm) | bottom gap (nm) | min gap (nm) | remaining air gap (nm) | lambda (nm) | Q |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 564.3 | 383.724 | 0.68 | 0 | 180.576 | 180.576 | 180.576 | 148.576 | 980.010144 | 1.6561e+11 |
| 2 | 564.3 | 393.724 | 0.697721 | -2 | 170.576 | 191.528 | 170.576 | 138.576 | 979.416964 | 1.6563e+11 |
| 3 | 569.3 | 388.724 | 0.68281 | -2 | 180.576 | 201.528 | 180.576 | 148.576 | 980.773967 | 8.57952e+11 |
| 4 | 562.3 | 388.724 | 0.691311 | 0 | 173.576 | 173.576 | 173.576 | 141.576 | 981.38139 | 1.60941e+11 |
| 5 | 569.3 | 378.724 | 0.665245 | 0 | 190.576 | 190.576 | 190.576 | 158.576 | 981.412222 | 9.32942e+10 |
| 6 | 566.3 | 393.724 | 0.695257 | -2 | 172.576 | 193.528 | 172.576 | 140.576 | 981.452769 | 9.81827e+11 |
| 7 | 559.3 | 393.724 | 0.703959 | 0 | 165.576 | 165.576 | 165.576 | 133.576 | 981.45309 | 5.98089e+11 |
| 8 | 566.3 | 378.724 | 0.668769 | 0 | 187.576 | 187.576 | 187.576 | 155.576 | 978.419194 | 1.81618e+11 |
| 9 | 559.3 | 388.724 | 0.695019 | 0 | 170.576 | 170.576 | 170.576 | 138.576 | 978.226781 | 9.78185e+11 |
| 10 | 562.3 | 383.724 | 0.682419 | 0 | 178.576 | 178.576 | 178.576 | 146.576 | 977.955987 | 2.37939e+11 |
| 11 | 566.3 | 383.724 | 0.677598 | 0 | 182.576 | 182.576 | 182.576 | 150.576 | 982.051014 | 1.64866e+11 |
| 12 | 566.3 | 388.724 | 0.686428 | -2 | 177.576 | 198.528 | 177.576 | 145.576 | 977.791293 | 3.73622e+11 |

## Lowest-Q cases

| P (nm) | W_top (nm) | fill | taper (deg) | min gap (nm) | lambda (nm) | Q |
|---:|---:|---:|---:|---:|---:|---:|
| 566.3 | 373.724 | 0.65994 | 2 | 171.624 | 1004.49185 | 533.76 |
| 564.3 | 373.724 | 0.662279 | 2 | 169.624 | 1000.94772 | 537.573 |
| 564.3 | 378.724 | 0.671139 | 2 | 164.624 | 1000.94772 | 537.573 |
| 564.3 | 383.724 | 0.68 | 2 | 159.624 | 1000.94772 | 537.573 |
| 562.3 | 373.724 | 0.664635 | -2 | 188.576 | 997.403556 | 541.425 |
| 562.3 | 373.724 | 0.664635 | 2 | 167.624 | 997.403556 | 541.425 |
| 562.3 | 378.724 | 0.673527 | 2 | 162.624 | 997.403556 | 541.425 |
| 562.3 | 388.724 | 0.691311 | 2 | 152.624 | 997.403556 | 541.425 |

## Recommendation

The large-gap design is robust to the scanned period and top-width errors when the sidewalls remain close to vertical: all `0 deg` cases stay above `Q=1e9`. The dominant risk is a `+2 deg` bottom-wider trapezoid, which can move the high-Q branch away from 980 nm and leave only a low-Q 980-nearest mode. For fabrication, prioritize vertical or slightly bottom-narrow sidewalls; if the process is known to produce bottom-wider sidewalls, retune `P`/`W` with the `+2 deg` model rather than using the vertical-wall nominal directly.
