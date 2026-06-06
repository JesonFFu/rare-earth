# 980 nm 3D Photonic-Crystal Comparison

All new files were written in the v7/v8 result folders only. Earlier v1-v6 models, source files, and `bic0413.mph` were not modified.

## v7 SiN Cylindrical Pillar Array

- Geometry: square-lattice 300 nm SiN pillars on SiO2/Si.
- Recommended practical-gap model: `E:\xitugrating\sin980_qbic_3d_design_v7_cyl_pillar_corrected_particles\outputs\sin980_qbic_v7_cyl_best_gapgt80_near_980.mph`
- P = 596.000 nm
- pillar radius R = 255.500 nm
- pillar-to-pillar gap = 85.000 nm
- wavelength = 979.934560 nm after saved-model re-solve
- Q = 18856.1

## v8 SiN Air-Hole Slab

- Geometry: square-lattice 300 nm continuous SiN slab with a through-etched circular air hole.
- Highest-Q near-980 saved model: `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sin980_qbic_v8_air_hole_focus_best.mph`
- P = 540.000 nm
- air-hole radius R = 74.000 nm
- SiN bridge width = 392.000 nm
- wavelength = 980.550636 nm
- Q = 200672

- Nearest-980 saved model: `E:\xitugrating\sin980_qbic_3d_design_v8_air_hole_slab\outputs\sin980_qbic_v8_air_hole_nearest_980.mph`
- P = 540.000 nm
- air-hole radius R = 75.600 nm
- SiN bridge width = 388.800 nm
- wavelength = 980.077319 nm
- Q = 199356

## Conclusion

The air-hole slab is the better 980 nm candidate from this pass. It improves Q by roughly one order of magnitude compared with the practical-gap cylindrical-pillar array while keeping a very large SiN bridge width, so it is also more fabrication-friendly.

Neither geometry reached `Q >= 1e9` in this 3D substrate/PML model. The result is therefore a BIC-like guided-resonance candidate rather than a confirmed ultra-high-Q BIC. The next check should be a small-k `Q(k)` sweep and mesh/PML convergence on the v8 air-hole model.
