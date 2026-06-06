# v7 input manifest

This v7 run is isolated from previous versions and writes only inside:

`E:\xitugrating\sin980_qbic_3d_design_v7_cyl_pillar_corrected_particles`

Read-only scale reference:

- `E:\xitugrating\bic0413.mph`
- Copied as `inputs\bic0413_readonly_scale_reference.mph`
- Used only as a broad parameter scale reference; the v7 target and geometry are different.

v7 target:

- 980 nm high-Q BIC-like eigenmode in a 3D square-lattice SiN cylindrical-pillar array.
- Clean model is optimized first.
- Corrected random rare-earth particles are added only as a perturbation after the clean BIC search.
- Main pillar gap constraint was relaxed for high-Q discovery; current focused scan allows `P-2R >= 10 nm`.
- This small-gap scan is exploratory and not yet a fabrication recommendation.

Corrected particle placement:

- Allowed: exposed SiO2 top surface, SiN pillar top surface, SiN pillar sidewall.
- Not allowed: suspended particles in air, particles embedded inside SiN, particles embedded inside SiO2.
- Particle index: n_re = 1.52.

Memory simplification:

- Effective SiO2 thickness is truncated to 700 nm.
- Effective Si buffer is 250 nm.
- Top and bottom PML thickness is 500 nm.
