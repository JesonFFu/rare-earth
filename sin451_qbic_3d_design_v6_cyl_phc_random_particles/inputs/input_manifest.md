# v6 input manifest

This v6 run is isolated from all previous versions.

Reference model copied read-only from v5:
- E:\xitugrating\sin451_qbic_2d_design_v5_particulate_re\outputs\sin451_qbic_v5_best_particulate_re.mph

Additional read-only scale reference:
- E:\xitugrating\bic0413.mph
- Extracted range anchor: period=356 nm, r/period=0.220, r≈78.32 nm
- Used only to add a parameter refinement range; the v6 geometry and objective remain different.

v6 target:
- 3D square-lattice SiN slab with cylindrical air-hole photonic-crystal array
- target wavelength 451 nm
- particulate rare-earth proxy n=1.52
- random-but-deterministic particles: size >=16 nm, mostly 24-44 nm diameter
- particles are denser on slab top and exposed hole bottom than on hole sidewall
- 3D unit cell selected because x/y periodicity and vertical radiation are essential for a cylindrical PhC BIC
- SiO2 thickness is truncated to 500 nm for the 16 GB memory budget
