# Gap >= 120 nm Recommendation

All files are inside `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sio2_1500_check`. The original v8 files were not overwritten.

## Recommended Exact-980 Model

- Model: `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sio2_1500_check\sin980_qbic_v8_air_hole_sio2_1500_exact980_gap120.mph`
- Geometry: 300 nm SiN air-hole slab on 1500 nm effective SiO2.
- Period `P = 542 nm`
- Air-hole radius `R = 74 nm`
- Gap / SiN bridge width `P - 2R = 394 nm`
- Resonance wavelength `lambda = 980.034515 nm`
- Q factor `Q = 2.52381e9`

This is the best recommendation when the priority is staying very close to 980 nm while keeping `gap >= 120 nm`.

## Higher-Q Nearby Candidate

- Model: `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sio2_1500_check\sin980_qbic_v8_air_hole_sio2_1500_optimized.mph`
- Period `P = 540 nm`
- Air-hole radius `R = 78 nm`
- Gap / SiN bridge width `P - 2R = 384 nm`
- Resonance wavelength `lambda = 979.319900 nm`
- Q factor `Q = 3.58215e9`

This is the higher-Q choice inside the local `|lambda-980| <= 1 nm` search window, but it is farther from exact 980 nm than the recommended model above.

## Field Plots

The saved `.mph` files include two 2D `|E|^2` slice plot groups:

- xy slice at the SiN slab mid-plane
- xz slice through the air-hole center
