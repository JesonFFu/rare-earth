# v10 3D Family Scan Memory Guard Review

## Current Status

- No v10 COMSOL/Python scan process is currently running.
- Existing service process: `D:\Anaconda\envs\comsol_env\python.exe -m src.server`.
- A Vivado-related Java process is present and is not part of the COMSOL scan.
- Current checked memory after user cleanup: total `15.34 GB`, free about `8.0 GB`.

## Guarded Test Result

Attempted the lightest 3D pilot family:

- Family: `single_air_hole_slab_reference`
- Design: `mono_hole_P540_R74`
- Profile: original full stack
- Cores: `1`
- Mesh: `autoMeshSize=7`
- Eigenmodes per wavelength: `6`
- Guard threshold: `7.5 GB` free physical memory

Observed:

- Before COMSOL client: free memory about `8.08 GB`.
- Before solving first model: free memory about `7.73 GB`.
- During first 3D solve: free memory dropped to `5.74 GB`.
- Guard stopped the process automatically.

Conclusion: with the current available RAM, even the lightest full 3D v10 eigenfrequency model can exceed the safe memory envelope. Continuing full 3D family scans without simplification risks freezing the machine.

## Partial 3D Output Already Present

A prior interrupted run produced one single-wavelength 3D model:

- File: `outputs/family3d/capsule_hole_Px700_Py520_L320_W140_x_980.mph`
- Family: `capsule_elliptic_air_hole`
- Geometry: `Px=700 nm`, `Py=520 nm`, capsule/ellipse-approximation hole `L=320 nm`, `W=140 nm`
- 980 nm result: `lambda=990.996016 nm`, `Q=325.525`
- 451 nm counterpart was not completed.

This is a partial artifact only. It is useful as a geometry/debug proof that the capsule-hole model can be built and solved, but it is not a valid dual-wavelength result and not a high-Q candidate.

## Gradient Simplification Strategy

All profiles keep the core physical stack rule:

- SiO2 is kept at `1500 nm`.
- Rare-earth film is kept at `50 nm`.
- SiO2 cap is kept at `12.5 nm`.
- SiN layer remains `300 nm` and is the only patterned layer.

### Level 0: Full Profile

- Air: `900 nm`
- PML: `500 nm`
- Si buffer: `250 nm`
- Mesh: recommended `autoMeshSize=7` or finer only for final verification
- Eigenmodes: `6-10`
- Status: stopped by guard at first 3D model under current memory.

### Level 1: Balanced Profile

- Air: `650 nm`
- PML: `350 nm`
- Si buffer: `150 nm`
- Mesh: `autoMeshSize=7`
- Eigenmodes: `6`
- Purpose: first safe 3D existence test while preserving the 1500 nm SiO2 layer and nontrivial radiation region.
- Recommended next attempt if approved.

### Level 2: Lite Profile

- Air: `450 nm`
- PML: `250 nm`
- Si buffer: `100 nm`
- Mesh: `autoMeshSize=7`
- Eigenmodes: `4-6`
- Purpose: only for broad family triage when balanced still trips the guard.
- Any candidate found here must be rechecked with the balanced or full profile before being considered reliable.

### Level 3: 2D / 2.5D Proxy

- Continue expanding fast 2D and reduced-dimensional screening.
- Use this to rank structure families and avoid wasting 3D memory on weak candidates.
- 3D should then be used only for the best one or two candidates per family.

## Current Recommended Path

1. Try one `balanced` 3D pilot for `single_air_hole_slab_reference`.
2. If memory stays above threshold, continue one family at a time in this order:
   - `rectangular_air_hole_slab`
   - `dual_hole_radius_supercell`
   - `dual_hole_shift_supercell`
   - `folded_2x2_air_hole_radius`
   - `folded_2x2_air_hole_shift`
   - `capsule_elliptic_air_hole`
   - `dual_capsule_elliptic_supercell`
   - `main_hole_with_offset_satellite`
   - `wide_slot_cross_hole`
   - `dual_period_1d_grating`
   - `tetramer_hole_slab`
   - `tetramer_pillar_array`
3. If balanced still trips the guard, use `lite` only to identify promising families, then reserve balanced/full checks for final candidates.

## Guarded Command Template

```powershell
powershell -ExecutionPolicy Bypass -File .\sin980_451_qbic_3d_design_v10_re50_sio2cap_dual_bic\scripts\run_family3d_guarded.ps1 `
  -Family single_air_hole_slab_reference `
  -Profile balanced `
  -MinFreeGB 7.5 `
  -Cores 1 `
  -MeshSize 7 `
  -Neigs 6 `
  -PollSeconds 20
```
