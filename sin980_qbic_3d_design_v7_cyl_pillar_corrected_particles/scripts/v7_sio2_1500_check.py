from __future__ import annotations

from pathlib import Path

import sin980_qbic_3d_v7_cyl_pillar_corrected_particles as v7


OUT = Path(__file__).resolve().parent.parent / "outputs" / "sio2_1500_check"
SIO2_NM = 1500.0
BASE_LAMBDA_NM = 979.9332948457478
BASE_Q = 18856.037357115474


def set_thick_sio2() -> None:
    v7.H_SIO2_NM = SIO2_NM


def solve_and_save(client, design: v7.Design, filename: str, model_name: str) -> v7.ModeResult:
    set_thick_sio2()
    model, layout = v7.build_model(client, design, model_name, include_particles=False)
    try:
        rows = v7.solve_design(model, design, layout, 0.0, 0.0)
        selected = v7.best_mode(rows)
        add_2d_field_plots(model, "sio2_1500")
        model.save(str(OUT / filename))
        v7.write_modes_csv(OUT / f"{Path(filename).stem}_modes.csv", rows)
        v7.write_modes_csv(OUT / f"{Path(filename).stem}_selected.csv", [selected])
        return selected
    finally:
        client.remove(model.name())


def add_2d_field_plots(model, label: str) -> bool:
    try:
        res = model.java.result()
        pg_xy = res.create(f"pg_E2_xy_{label}", "PlotGroup3D")
        pg_xy.label(f"{label}: |E|^2 xy slice at pillar mid-plane")
        pg_xy.set("data", "dset1")
        slc_xy = pg_xy.feature().create("slc_xy", "Slice")
        slc_xy.label("|E|^2 xy slice")
        slc_xy.set("expr", "ewfd.normE^2")
        slc_xy.set("quickplane", "xy")
        slc_xy.set("quickz", "z_pillar_mid")
        slc_xy.set("resolution", "normal")

        pg_xz = res.create(f"pg_E2_xz_{label}", "PlotGroup3D")
        pg_xz.label(f"{label}: |E|^2 xz slice through pillar center")
        pg_xz.set("data", "dset1")
        slc_xz = pg_xz.feature().create("slc_xz", "Slice")
        slc_xz.label("|E|^2 xz slice")
        slc_xz.set("expr", "ewfd.normE^2")
        slc_xz.set("quickplane", "zx")
        slc_xz.set("quicky", "0")
        slc_xz.set("resolution", "normal")
        return True
    except Exception as exc:
        print(f"  2D field plot setup skipped: {exc}", flush=True)
        return False


def changed(row: v7.ModeResult) -> bool:
    lambda_shift = abs(row.wavelength_nm - BASE_LAMBDA_NM)
    q_rel = abs(row.q_value / BASE_Q - 1.0)
    return lambda_shift > 0.5 or q_rel > 0.20


def local_designs() -> list[v7.Design]:
    out: list[v7.Design] = []
    for period in [594.0, 596.0, 598.0]:
        for gap in [80.0, 85.0, 90.0]:
            radius = (period - gap) / 2.0
            out.append(v7.Design(
                f"v7_sio2_1500_P{period:.0f}_G{gap:.0f}",
                period,
                radius,
                v7.stable_seed(f"v7_sio2_1500_{period}_{gap}"),
            ))
    return out


def choose_best(rows: list[v7.ModeResult]) -> v7.ModeResult:
    near = [r for r in rows if abs(r.wavelength_nm - v7.TARGET_NM) <= 3.0]
    pool = near or rows
    return sorted(pool, key=lambda r: (-r.q_value, abs(r.wavelength_nm - v7.TARGET_NM)))[0]


def write_summary(base_1500: v7.ModeResult, optimized: v7.ModeResult | None) -> None:
    lambda_shift = base_1500.wavelength_nm - BASE_LAMBDA_NM
    q_ratio = base_1500.q_value / BASE_Q
    lines = [
        "# v7 Cylinder SiO2 1500 nm Check",
        "",
        f"All files are inside `{OUT}`. Original v7 files were not overwritten.",
        "",
        "## Original reference",
        "- h_sio2_eff = 700 nm",
        "- P = 596 nm, R = 255.5 nm, gap = 85 nm",
        f"- lambda = {BASE_LAMBDA_NM:.6f} nm",
        f"- Q = {BASE_Q:.6g}",
        "",
        "## Same geometry with h_sio2_eff = 1500 nm",
        f"- model = `{OUT / 'sin980_qbic_v7_cyl_sio2_1500_same_geometry.mph'}`",
        f"- lambda = {base_1500.wavelength_nm:.6f} nm",
        f"- Q = {base_1500.q_value:.6g}",
        f"- wavelength shift = {lambda_shift:+.6f} nm",
        f"- Q ratio vs 700 nm = {q_ratio:.6g}",
        "",
        "## Local re-optimization",
    ]
    if optimized is None:
        lines.append("- Not triggered: wavelength/Q change stayed below the set threshold.")
    else:
        lines.extend([
            f"- model = `{OUT / 'sin980_qbic_v7_cyl_sio2_1500_optimized.mph'}`",
            f"- P = {optimized.period_nm:.3f} nm",
            f"- R = {optimized.radius_nm:.3f} nm",
            f"- gap = {optimized.axis_gap_nm:.3f} nm",
            f"- lambda = {optimized.wavelength_nm:.6f} nm",
            f"- Q = {optimized.q_value:.6g}",
        ])
    lines.extend([
        "",
        "## Field plots",
        "- The saved .mph files include two 2D `|E|^2` slice plot groups: xy mid-plane and xz center-plane.",
    ])
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    set_thick_sio2()
    mph = v7.import_mph()
    client = mph.Client(cores=4)
    try:
        base_design = v7.Design("v7_sio2_1500_same_P596_R255p5", 596.0, 255.5, v7.stable_seed("v7_sio2_1500_same"))
        base_1500 = solve_and_save(
            client,
            base_design,
            "sin980_qbic_v7_cyl_sio2_1500_same_geometry.mph",
            "sin980_qbic_v7_cyl_sio2_1500_same_geometry",
        )
        optimized: v7.ModeResult | None = None
        if changed(base_1500):
            rows = v7.scan_designs(client, local_designs(), OUT / "v7_sio2_1500_local_scan.csv", "v7-sio2-1500", include_particles=False)
            best = choose_best(rows)
            optimized = solve_and_save(
                client,
                best if isinstance(best, v7.Design) else v7.Design(best.design_id, best.period_nm, best.radius_nm, v7.stable_seed(best.design_id)),
                "sin980_qbic_v7_cyl_sio2_1500_optimized.mph",
                "sin980_qbic_v7_cyl_sio2_1500_optimized",
            )
            v7.write_modes_csv(OUT / "v7_sio2_1500_optimized_selected.csv", [optimized])
        write_summary(base_1500, optimized)
        print(
            "V7_SAME_1500 "
            f"P={base_1500.period_nm:.3f} R={base_1500.radius_nm:.3f} gap={base_1500.axis_gap_nm:.3f} "
            f"lambda={base_1500.wavelength_nm:.6f} Q={base_1500.q_value:.6g}",
            flush=True,
        )
        if optimized:
            print(
                "V7_OPT_1500 "
                f"P={optimized.period_nm:.3f} R={optimized.radius_nm:.3f} gap={optimized.axis_gap_nm:.3f} "
                f"lambda={optimized.wavelength_nm:.6f} Q={optimized.q_value:.6g}",
                flush=True,
            )
    finally:
        client.clear()


if __name__ == "__main__":
    main()
