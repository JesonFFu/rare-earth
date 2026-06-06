# v10 Continued Family Scan Report

## Safety And Memory Status

This continuation was run with a strict single-process policy. No parallel COMSOL model solves were launched.

Current machine behavior:

- Available memory before 3D solve attempts was about `7.2 GB`.
- Full 3D and lite 3D eigenfrequency solves both triggered the memory guard.
- The most conservative solve attempt used `profile=lite`, `cores=1`, `mesh=7`, `neigs=4`; available memory still dropped to about `5.10 GB`, and the guard stopped the process.

Conclusion: under the current RAM state, 3D eigenfrequency solving is not safe enough to continue unattended. I therefore completed geometry-only executable models for every requested structure family, and continued numerical screening with lower-memory 2D/1D models.

## Structure Families Covered

All requested families now have independent executable COMSOL `.mph` files under:

`E:\xitugrating\sin980_451_qbic_3d_design_v10_re50_sio2cap_dual_bic\outputs\family3d_geometry_only`

The files are geometry/study/material executable templates in the `balanced` profile:

- `single_air_hole_slab_reference`: square-lattice circular air-hole slab near v8.
- `rectangular_air_hole_slab`: rectangular-lattice circular air-hole slab.
- `dual_hole_radius_supercell`: 2x1 supercell with alternating hole radii.
- `dual_hole_shift_supercell`: 2x1 supercell with shifted hole positions.
- `folded_2x2_air_hole_radius`: 2x2 folded supercell with alternating radii.
- `folded_2x2_air_hole_shift`: 2x2 folded supercell with position offsets.
- `capsule_elliptic_air_hole`: manufacturable capsule/elliptical air-hole approximation.
- `dual_capsule_elliptic_supercell`: two different capsule/elliptical holes in a supercell.
- `main_hole_with_offset_satellite`: circular hole plus offset satellite hole.
- `wide_slot_cross_hole`: cross-like through-etched SiN air opening.
- `dual_period_1d_grating`: 1D two-period through-slot grating.
- `tetramer_hole_slab`: four-hole tetramer/quadrumer slab.
- `tetramer_pillar_array`: four-pillar tetramer/quadrumer array.

All generated files keep:

- `h_sio2_eff = 1500 nm`
- `h_re = 50 nm`
- `h_cap_sio2 = 12.5 nm`
- `h_sin = 300 nm`
- Patterning only in the SiN layer

## Numerical Screening Completed

The low-memory fast screening results are consolidated in:

- `combined_fast2d_candidates.csv`
- `combined_fast2d_q_distribution.svg`

The completed scans include:

- 75-point 2D 1D/coarse supercell grating scan.
- 70-point local scan around the 980 nm high-Q grating branch.
- 100-point wide-shallow-groove scan.

Best 980-prioritized 2D candidate:

- Design: `fine_primitive_P562_f0p72`
- Period: `562 nm`
- SiN width: `404.64 nm`
- Gap: `157.36 nm`
- 980 nm: `lambda = 980.115156 nm`, `Q = 3.54213e9`
- 451 nm: `lambda = 453.291596 nm`, `Q = 1300.06`

Best 451 high-Q control:

- Design: `primitive_P294_f0p42`
- Period: `294 nm`
- SiN width: `123.48 nm`
- Gap: `170.52 nm`
- 451 nm: `lambda = 445.026698 nm`, `Q = 3.35753e9`
- 980 nm: `lambda = 990.135049 nm`, `Q = 7.055`

Wide shallow groove result:

- The top-groove family preserved high 980 Q in some cases, but did not raise the 451 branch to `Q >= 1e5`.
- Best groove points still had 451 Q only in the `10^3` range or lower.

## Current Scientific Conclusion

The search has not yet found a manufacturable single-SiN-layer structure that simultaneously satisfies:

- `Q980 >= 1e8`
- `Q451 >= 1e5`

The evidence so far is consistent across three low-memory routes:

- A 980 nm BIC-like branch exists and can be tuned close to 980 nm with very high Q.
- A 451 nm high-Q branch exists in short-period 1D models.
- The two high-Q branches do not naturally overlap in the same simple single-layer SiN geometry with the continuous 50 nm rare-earth film and 12.5 nm cap.

## What Remains To Do

When memory permits, run guarded 3D eigenfrequency solves family by family, starting with:

1. `rectangular_air_hole_slab`
2. `dual_hole_radius_supercell`
3. `dual_hole_shift_supercell`
4. `folded_2x2_air_hole_radius`
5. `folded_2x2_air_hole_shift`
6. `capsule_elliptic_air_hole`
7. `dual_capsule_elliptic_supercell`
8. `main_hole_with_offset_satellite`
9. `wide_slot_cross_hole`
10. `tetramer_hole_slab`
11. `tetramer_pillar_array`

Recommended guarded command template:

```powershell
powershell -ExecutionPolicy Bypass -File .\sin980_451_qbic_3d_design_v10_re50_sio2cap_dual_bic\scripts\run_family3d_guarded.ps1 `
  -Family dual_hole_radius_supercell `
  -Profile lite `
  -MinFreeGB 6.0 `
  -Cores 1 `
  -MeshSize 7 `
  -Neigs 4 `
  -PollSeconds 10
```

If this still trips the guard, the next practical step is not to roughen the mesh further, but to reduce the number of 3D candidates using more 2D/2.5D screening, then reserve 3D for the top one or two structures.
