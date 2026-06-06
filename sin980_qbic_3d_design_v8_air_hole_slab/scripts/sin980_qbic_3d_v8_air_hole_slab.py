from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
OUT = ROOT / "outputs"
V7_SCRIPT_DIR = ROOT.parent / "sin980_qbic_3d_design_v7_cyl_pillar_corrected_particles" / "scripts"
sys.path.insert(0, str(V7_SCRIPT_DIR))

import sin980_qbic_3d_v7_cyl_pillar_corrected_particles as base


TARGET_NM = base.TARGET_NM
N_SIN = base.N_SIN
N_SIO2 = base.N_SIO2
N_SI = base.N_SI
H_SIN_NM = base.H_SIN_NM
H_SIO2_NM = base.H_SIO2_NM
H_SI_BUF_NM = base.H_SI_BUF_NM
H_PML_NM = base.H_PML_NM
H_AIR_NM = base.H_AIR_NM
MIN_NECK_NM = 80.0
MIN_HOLE_DIAMETER_NM = 140.0


Design = base.Design
ModeResult = base.ModeResult


def build_model(client, design: Design, name: str):
    model = client.create(name)
    jmodel = model.java
    half = design.period_nm / 2.0
    z_si = H_PML_NM
    z_sio2 = H_PML_NM + H_SI_BUF_NM
    z_sin = z_sio2 + H_SIO2_NM
    z_top = z_sin + H_SIN_NM
    z_top_pml = z_top + H_AIR_NM
    total_height = 2 * H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_SIN_NM + H_AIR_NM

    base.set_param(model, "lambda_target", f"{TARGET_NM}[nm]", "Target wavelength")
    base.set_param(model, "P", f"{design.period_nm:.9g}[nm]", "Square lattice period")
    base.set_param(model, "R", f"{design.radius_nm:.9g}[nm]", "Air-hole radius")
    base.set_param(model, "axis_gap", f"{design.axis_gap_nm:.9g}[nm]", "Remaining SiN bridge width P-2R")
    base.set_param(model, "h_sin", "300[nm]", "SiN slab thickness and air-hole etch depth")
    base.set_param(model, "h_sio2_eff", f"{H_SIO2_NM:.9g}[nm]", "Memory-truncated SiO2 buffer")
    base.set_param(model, "z_sio2_top", f"{z_sin:.9g}[nm]", "Absolute SiO2 top / SiN base z")
    base.set_param(model, "z_slab_mid", f"{z_sin + 0.5 * H_SIN_NM:.9g}[nm]", "Absolute slab mid-plane z")
    base.set_param(model, "z_slab_top", f"{z_top:.9g}[nm]", "Absolute slab top z")
    base.set_param(model, "n_sin", f"{N_SIN:.9f}", "SiN index at 980 nm")
    base.set_param(model, "n_sio2", f"{N_SIO2:.9f}", "SiO2 index at 980 nm")
    base.set_param(model, "n_si", f"{N_SI:.9f}", "Effective Si index at 980 nm")
    base.set_param(model, "kx", "0[1/m]", "Floquet kx")
    base.set_param(model, "ky", "0[1/m]", "Floquet ky")

    jmodel.component().create("comp1", True)
    comp = jmodel.component("comp1")
    comp.label("3D square cell SiN air-hole slab photonic crystal")
    comp.geom().create("geom1", 3)
    geom = comp.geom("geom1")
    geom.lengthUnit("nm")

    base.create_block(geom, "si_pml", (-half, -half, 0.0), (design.period_nm, design.period_nm, H_PML_NM), "Si bottom PML")
    base.create_block(geom, "si_buf", (-half, -half, z_si), (design.period_nm, design.period_nm, H_SI_BUF_NM), "Si buffer")
    base.create_block(geom, "sio2", (-half, -half, z_sio2), (design.period_nm, design.period_nm, H_SIO2_NM), "SiO2")
    base.create_block(geom, "sin_slab_blk", (-half, -half, z_sin), (design.period_nm, design.period_nm, H_SIN_NM), "300 nm SiN slab before air-hole subtraction")
    base.create_cylinder(geom, "air_hole", 0.0, 0.0, z_sin, design.radius_nm, H_SIN_NM, "300 nm through-etched circular air hole")
    base.create_difference(geom, "sin_slab", ["sin_slab_blk"], ["air_hole"], "SiN slab after through air-hole etch", keep_tools=True)
    base.create_block(geom, "air_upper", (-half, -half, z_top), (design.period_nm, design.period_nm, H_AIR_NM), "air above slab")
    base.create_block(geom, "air_pml", (-half, -half, z_top_pml), (design.period_nm, design.period_nm, H_PML_NM), "air top PML")
    geom.run()

    base.create_union_selection(comp, "air_dom", "3", ["geom1_air_hole_dom", "geom1_air_upper_dom", "geom1_air_pml_dom"])
    base.create_union_selection(comp, "si_dom", "3", ["geom1_si_pml_dom", "geom1_si_buf_dom"])
    base.create_box_selection(comp, "x_left_bnd", "2", f"{-half-1}", f"{-half+1}", f"{-half-1}", f"{half+1}", "-1", f"{total_height+1}")
    base.create_box_selection(comp, "x_right_bnd", "2", f"{half-1}", f"{half+1}", f"{-half-1}", f"{half+1}", "-1", f"{total_height+1}")
    base.create_union_selection(comp, "x_periodic_bnd", "2", ["x_left_bnd", "x_right_bnd"])
    base.create_box_selection(comp, "y_front_bnd", "2", f"{-half-1}", f"{half+1}", f"{-half-1}", f"{-half+1}", "-1", f"{total_height+1}")
    base.create_box_selection(comp, "y_back_bnd", "2", f"{-half-1}", f"{half+1}", f"{half-1}", f"{half+1}", "-1", f"{total_height+1}")
    base.create_union_selection(comp, "y_periodic_bnd", "2", ["y_front_bnd", "y_back_bnd"])

    pml_bot = comp.coordSystem().create("pml_bot", "PML")
    pml_bot.label("Bottom Si PML")
    pml_bot.selection().named("geom1_si_pml_dom")
    pml_top = comp.coordSystem().create("pml_top", "PML")
    pml_top.label("Top air PML")
    pml_top.selection().named("geom1_air_pml_dom")

    base.create_refractive_material(comp, "mat_air", "Air in etched hole and top region", "air_dom", "1", "0")
    base.create_refractive_material(comp, "mat_sio2", "SiO2 980 nm", "geom1_sio2_dom", "n_sio2", "0")
    base.create_refractive_material(comp, "mat_si", "Si effective 980 nm", "si_dom", "n_si", "0")
    base.create_refractive_material(comp, "mat_sin", "SiN air-hole slab 980 nm", "geom1_sin_slab_dom", "n_sin", "0")

    ewfd = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
    ewfd.label("Electromagnetic Waves, Frequency Domain")
    pcx = ewfd.create("pcx", "PeriodicCondition", 2)
    pcx.label("Floquet periodic x")
    pcx.selection().named("x_periodic_bnd")
    pcx.set("PeriodicType", "Floquet")
    pcx.set("kFloquet", ["kx", "ky", "0"])
    pcy = ewfd.create("pcy", "PeriodicCondition", 2)
    pcy.label("Floquet periodic y")
    pcy.selection().named("y_periodic_bnd")
    pcy.set("PeriodicType", "Floquet")
    pcy.set("kFloquet", ["kx", "ky", "0"])

    mesh = comp.mesh().create("mesh1")
    mesh.label("Extra-coarse 3D mesh for memory-limited 980 nm air-hole scan")
    mesh.autoMeshSize(7)

    study = jmodel.study().create("std1")
    study.label("Eigenfrequency near 980 nm")
    eig = study.create("eig", "Eigenfrequency")
    eig.set("shift", "c_const/lambda_target")
    eig.set("neigs", "12")
    eig.set("neigsmanual", "12")
    eig.set("neigsactive", "on")
    eig.set("eigunit", "THz")
    return model


def solve_design(model, design: Design, kx_norm: float = 0.0, ky_norm: float = 0.0) -> list[ModeResult]:
    base.set_param(model, "kx", f"{kx_norm:.12g}*pi/P")
    base.set_param(model, "ky", f"{ky_norm:.12g}*pi/P")
    model.java.component("comp1").geom("geom1").run()
    model.java.component("comp1").mesh("mesh1").run()
    model.java.study("std1").run()
    modes = base.get_eigenfrequencies_thz(model)
    empty = base.ParticleLayout((), (), ())
    return [base.mode_result(design, empty, kx_norm, ky_norm, i + 1, freq) for i, freq in enumerate(modes)]


def best_mode(rows: list[ModeResult]) -> ModeResult:
    return sorted(rows, key=lambda r: (r.score, abs(r.wavelength_nm - TARGET_NM), -r.q_value))[0]


def scan_designs(client, designs: list[Design], path: Path, label: str) -> list[ModeResult]:
    rows = base.read_modes_csv(path)
    done = {r.design_id for r in rows}
    for i, design in enumerate(designs, 1):
        if design.design_id in done:
            continue
        print(f"[{label}] {i}/{len(designs)} P={design.period_nm:.3f} R={design.radius_nm:.3f} neck={design.axis_gap_nm:.3f}", flush=True)
        model = build_model(client, design, f"v8_{design.design_id}")
        try:
            current = solve_design(model, design, 0.0, 0.0)
            rows.extend(current)
            base.write_modes_csv(path, rows)
            best = best_mode(current)
            print(f"  best lambda={best.wavelength_nm:.6f} nm Q={best.q_value:.4g} mode={best.mode_index}", flush=True)
        except Exception as exc:
            print(f"  FAILED {design.design_id}: {exc}", flush=True)
        finally:
            client.remove(model.name())
    return rows


def coarse_designs() -> list[Design]:
    designs: list[Design] = []
    for period in [500.0, 540.0, 580.0, 620.0, 660.0, 700.0, 740.0]:
        for ratio in [0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42]:
            radius = period * ratio
            design = Design(
                f"hole_coarse_P{period:.0f}_rho{ratio:.2f}".replace(".", "p"),
                period,
                radius,
                base.stable_seed(f"hole_coarse_{period}_{ratio}"),
            )
            if design.axis_gap_nm >= MIN_NECK_NM and 2.0 * radius >= MIN_HOLE_DIAMETER_NM:
                designs.append(design)
    return designs


def fine_designs(seed: ModeResult) -> list[Design]:
    designs: list[Design] = []
    seen: set[tuple[float, float]] = set()
    for dp in [-20.0, -10.0, 0.0, 10.0, 20.0]:
        for dr in [-12.0, -6.0, 0.0, 6.0, 12.0]:
            period = round(seed.period_nm + dp, 3)
            radius = round(seed.radius_nm + dr, 3)
            key = (period, radius)
            if key in seen:
                continue
            seen.add(key)
            design = Design(
                f"hole_fine_P{period:.0f}_R{radius:.0f}",
                period,
                radius,
                base.stable_seed(f"hole_fine_{period}_{radius}"),
            )
            if design.axis_gap_nm >= MIN_NECK_NM and 2.0 * radius >= MIN_HOLE_DIAMETER_NM:
                designs.append(design)
    return designs


def choose_candidates(rows: list[ModeResult], limit: int = 8) -> list[ModeResult]:
    by_design = []
    for design_id in sorted({r.design_id for r in rows}):
        by_design.append(best_mode([r for r in rows if r.design_id == design_id]))
    return sorted(by_design, key=lambda r: (r.score, abs(r.wavelength_nm - TARGET_NM), -r.q_value))[:limit]


def save_model(client, row: ModeResult, path: Path, name: str) -> ModeResult:
    design = Design(row.design_id, row.period_nm, row.radius_nm, base.stable_seed(row.design_id))
    model = build_model(client, design, name)
    try:
        rows = solve_design(model, design, 0.0, 0.0)
        selected = min(rows, key=lambda r: abs(r.wavelength_nm - row.wavelength_nm) + 1e-6 * abs(math.log10(max(r.q_value, 1.0)) - math.log10(max(row.q_value, 1.0))))
        base.add_field_plot_groups(model, "air_hole")
        model.save(str(path))
        return selected
    finally:
        client.remove(model.name())


def run_band(client, row: ModeResult) -> list[ModeResult]:
    design = Design(row.design_id, row.period_nm, row.radius_nm, base.stable_seed(row.design_id))
    model = build_model(client, design, "sin980_qbic_v8_air_hole_band")
    kpath = [
        (0.0, 0.0), (0.02, 0.0), (0.05, 0.0), (0.1, 0.0), (0.2, 0.0), (0.35, 0.0), (0.5, 0.0),
        (0.5, 0.25), (0.5, 0.5), (0.35, 0.35), (0.2, 0.2), (0.1, 0.1), (0.05, 0.05), (0.0, 0.0),
    ]
    rows: list[ModeResult] = []
    try:
        for kx, ky in kpath:
            print(f"[band] k=({kx:.3g},{ky:.3g}) pi/P", flush=True)
            current = solve_design(model, design, kx, ky)
            rows.append(best_mode(current))
        base.write_modes_csv(OUT / "band_air_hole_980.csv", rows)
        base.write_modes_csv(OUT / "q_vs_k_air_hole_980.csv", [r for r in rows if r.ky_norm == 0.0 and r.kx_norm <= 0.2])
    finally:
        client.remove(model.name())
    return rows


def write_summary(clean_best: ModeResult, saved_best: ModeResult, band: list[ModeResult]) -> None:
    lines = [
        "# v8 980 nm SiN Air-Hole Slab BIC Search",
        "",
        f"All files were generated inside `{ROOT}`; prior versions were not modified.",
        "",
        "## Geometry",
        "- 3D square-lattice COMSOL unit cell.",
        "- 300 nm SiN slab on truncated SiO2/Si with a through-etched circular air hole.",
        "- `R` is air-hole radius; `axis_gap = P - 2R` is remaining SiN bridge width.",
        "",
        "## Best saved model",
        f"- P = {saved_best.period_nm:.3f} nm",
        f"- air-hole radius R = {saved_best.radius_nm:.3f} nm",
        f"- SiN bridge width = {saved_best.axis_gap_nm:.3f} nm",
        f"- wavelength = {saved_best.wavelength_nm:.6f} nm",
        f"- Q = {saved_best.q_value:.6g}",
        f"- model = `{OUT / 'sin980_qbic_v8_air_hole_best.mph'}`",
        "",
        "## Interpretation",
    ]
    if saved_best.q_value >= 1e9:
        lines.append("- A `Q >= 1e9` 980 nm BIC candidate was found in the air-hole slab geometry.")
    else:
        lines.append("- The scan did not find a `Q >= 1e9` mode; the saved model is the best near-980 candidate from this pass.")
    if band:
        lines.extend([
            "",
            "## Band / Q check",
            f"- Band points solved: {len(band)}",
            f"- Gamma point Q in saved branch: {saved_best.q_value:.6g}",
        ])
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_all(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = base.import_mph()
    client = mph.Client(cores=args.cores)
    try:
        coarse_rows = scan_designs(client, coarse_designs(), OUT / "air_hole_coarse_scan.csv", "air-hole-coarse")
        coarse_candidates = choose_candidates(coarse_rows, 8)
        base.write_modes_csv(OUT / "air_hole_coarse_best_candidates.csv", coarse_candidates)
        fine_rows: list[ModeResult] = []
        if not args.skip_fine and coarse_candidates:
            fine_rows = scan_designs(client, fine_designs(coarse_candidates[0]), OUT / "air_hole_fine_scan.csv", "air-hole-fine")
        all_rows = coarse_rows + fine_rows
        candidates = choose_candidates(all_rows, 10)
        base.write_modes_csv(OUT / "air_hole_best_candidates.csv", candidates)
        best = candidates[0]
        saved = save_model(client, best, OUT / "sin980_qbic_v8_air_hole_best.mph", "sin980_qbic_v8_air_hole_best")
        band: list[ModeResult] = []
        if not args.skip_band:
            band = run_band(client, saved)
        write_summary(best, saved, band)
        print(
            "BEST_AIR_HOLE "
            f"P={saved.period_nm:.3f} R={saved.radius_nm:.3f} neck={saved.axis_gap_nm:.3f} "
            f"lambda={saved.wavelength_nm:.6f} Q={saved.q_value:.6g}",
            flush=True,
        )
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="v8 980 nm SiN air-hole slab BIC search.")
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--skip-fine", action="store_true")
    parser.add_argument("--skip-band", action="store_true")
    args = parser.parse_args()
    run_all(args)


if __name__ == "__main__":
    main()
