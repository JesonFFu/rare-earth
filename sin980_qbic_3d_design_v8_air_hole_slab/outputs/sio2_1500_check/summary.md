# v8 Air-Hole Slab SiO2 1500 nm Check

All files are inside `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sio2_1500_check`. Original v8 files were not overwritten.

## Original reference
- h_sio2_eff = 700 nm
- nearest 980: P = 540 nm, R = 75.6 nm, lambda = 980.077319 nm, Q = 199356
- high-Q near 980: P = 540 nm, R = 74 nm, lambda = 980.550636 nm, Q = 200672

## Same geometries with h_sio2_eff = 1500 nm
- nearest 980 model = `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sio2_1500_check\sin980_qbic_v8_air_hole_sio2_1500_nearest_980.mph`
- nearest 980: lambda = 980.089864 nm, Q = 2.12372e+09, shift = +0.012545 nm, Q ratio = 10652.9
- high-Q model = `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sio2_1500_check\sin980_qbic_v8_air_hole_sio2_1500_high_q.mph`
- high-Q: lambda = 980.562017 nm, Q = 2.73707e+09, shift = +0.011382 nm, Q ratio = 13639.5

## Local re-optimization
- model = `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sio2_1500_check\sin980_qbic_v8_air_hole_sio2_1500_optimized.mph`
- P = 540.000 nm
- R = 78.000 nm
- SiN bridge width = 384.000 nm
- lambda = 979.319900 nm
- Q = 3.58215e+09

## Field plots
- The saved .mph files include two 2D `|E|^2` slice plot groups: xy slab mid-plane and xz center-plane.
