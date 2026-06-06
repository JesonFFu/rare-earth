from __future__ import annotations

from pathlib import Path

import sin980_qbic_3d_v8_air_hole_slab as v8


OUT = Path(__file__).resolve().parent.parent / "outputs" / "sio2_1500_check"
SIO2_NM = 1500.0

BASE_NEAREST = {
    "label": "nearest_980",
    "period_nm": 540.0,
    "radius_nm": 75.6,
    "lambda_nm": 980.0773190620838,
    "q": 199355.75282980665,
}

BASE_HIGH_Q = {
    "label": "high_q_near_980",
    "period_nm": 540.0,
    "radius_nm": 74.0,
    "lambda_nm": 980.5506358007459,
    "q": 200671.85712462344,
}


def set_thick_sio2() -> None:
    v8.H_SIO2_NM = SIO2_NM


def add_2d_field_plots(model, label: str) -> bool:
    try:
        res = model.java.result()
        pg_xy = res.create(f"pg_E2_xy_{label}", "PlotGroup3D")
        pg_xy.label(f"{label}: |E|^2 xy slice at slab mid-plane")
        pg_xy.set("data", "dset1")
        slc_xy = pg_xy.feature().create("slc_xy", "Slice")
        slc_xy.label("|E|^2 xy slice")
        slc_xy.set("expr", "ewfd.normE^2")
        slc_xy.set("quickplane", "xy")
        slc_xy.set("quickz", "z_slab_mid")
        slc_xy.set("resolution", "normal")

        pg_xz = res.create(f"pg_E2_xz_{label}", "PlotGroup3D")
        pg_xz.label(f"{label}: |E|^2 xz slice through air-hole center")
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


def solve_and_save(client, design: v8.Design, filename: str, model_name: str, label: str) -> v8.ModeResult:
    set_thick_sio2()
    model = v8.build_model(client, design, model_name)
    try:
        rows = v8.solve_design(model, design, 0.0, 0.0)
        selected = v8.best_mode(rows)
        add_2d_field_plots(model, label)
        model.save(str(OUT / filename))
        v8.base.write_modes_csv(OUT / f"{Path(filename).stem}_modes.csv", rows)
        v8.base.write_modes_csv(OUT / f"{Path(filename).stem}_selected.csv", [selected])
        return selected
    finally:
        client.remove(model.name())


def design_from_base(base: dict[str, float | str]) -> v8.Design:
    return v8.Design(
        f"v8_sio2_1500_{base['label']}",
        float(base["period_nm"]),
        float(base["radius_nm"]),
        v8.base.stable_seed(f"v8_sio2_1500_{base['label']}"),
    )


def changed(row: v8.ModeResult, base: dict[str, float | str]) -> bool:
    lambda_shift = abs(row.wavelength_nm - float(base["lambda_nm"]))
    q_rel = abs(row.q_value / float(base["q"]) - 1.0)
    return lambda_shift > 0.5 or q_rel > 0.20


def local_designs() -> list[v8.Design]:
    out: list[v8.Design] = []
    for period in [538.0, 540.0, 542.0]:
        for radius in [72.0, 74.0, 75.6, 78.0]:
            d = v8.Design(
                f"v8_sio2_1500_P{period:.0f}_R{str(radius).replace('.', 'p')}",
                period,
                radius,
                v8.base.stable_seed(f"v8_sio2_1500_{period}_{radius}"),
            )
            if d.axis_gap_nm >= 80.0 and 2.0 * radius >= 140.0:
                out.append(d)
    return out


def choose_best(rows: list[v8.ModeResult]) -> v8.ModeResult:
    near = [r for r in rows if abs(r.wavelength_nm - v8.TARGET_NM) <= 1.0]
    pool = near or rows
    return sorted(pool, key=lambda r: (-r.q_value, abs(r.wavelength_nm - v8.TARGET_NM)))[0]


def write_summary(nearest: v8.ModeResult, high_q: v8.ModeResult, optimized: v8.ModeResult | None) -> None:
    nearest_shift = nearest.wavelength_nm - float(BASE_NEAREST["lambda_nm"])
    high_q_shift = high_q.wavelength_nm - float(BASE_HIGH_Q["lambda_nm"])
    lines = [
        "# v8 Air-Hole Slab SiO2 1500 nm Check",
        "",
        f"All files are inside `{OUT}`. Original v8 files were not overwritten.",
        "",
        "## Original reference",
        "- h_sio2_eff = 700 nm",
        "- nearest 980: P = 540 nm, R = 75.6 nm, lambda = 980.077319 nm, Q = 199356",
        "- high-Q near 980: P = 540 nm, R = 74 nm, lambda = 980.550636 nm, Q = 200672",
        "",
        "## Same geometries with h_sio2_eff = 1500 nm",
        f"- nearest 980 model = `{OUT / 'sin980_qbic_v8_air_hole_sio2_1500_nearest_980.mph'}`",
        f"- nearest 980: lambda = {nearest.wavelength_nm:.6f} nm, Q = {nearest.q_value:.6g}, shift = {nearest_shift:+.6f} nm, Q ratio = {nearest.q_value / float(BASE_NEAREST['q']):.6g}",
        f"- high-Q model = `{OUT / 'sin980_qbic_v8_air_hole_sio2_1500_high_q.mph'}`",
        f"- high-Q: lambda = {high_q.wavelength_nm:.6f} nm, Q = {high_q.q_value:.6g}, shift = {high_q_shift:+.6f} nm, Q ratio = {high_q.q_value / float(BASE_HIGH_Q['q']):.6g}",
        "",
        "## Local re-optimization",
    ]
    if optimized is None:
        lines.append("- Not triggered: both checked models stayed below the wavelength/Q change threshold.")
    else:
        lines.extend([
            f"- model = `{OUT / 'sin980_qbic_v8_air_hole_sio2_1500_optimized.mph'}`",
            f"- P = {optimized.period_nm:.3f} nm",
            f"- R = {optimized.radius_nm:.3f} nm",
            f"- SiN bridge width = {optimized.axis_gap_nm:.3f} nm",
            f"- lambda = {optimized.wavelength_nm:.6f} nm",
            f"- Q = {optimized.q_value:.6g}",
        ])
    lines.extend([
        "",
        "## Field plots",
        "- The saved .mph files include two 2D `|E|^2` slice plot groups: xy slab mid-plane and xz center-plane.",
    ])
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    set_thick_sio2()
    mph = v8.base.import_mph()
    client = mph.Client(cores=4)
    try:
        nearest = solve_and_save(
            client,
            design_from_base(BASE_NEAREST),
            "sin980_qbic_v8_air_hole_sio2_1500_nearest_980.mph",
            "sin980_qbic_v8_air_hole_sio2_1500_nearest_980",
            "nearest_980",
        )
        high_q = solve_and_save(
            client,
            design_from_base(BASE_HIGH_Q),
            "sin980_qbic_v8_air_hole_sio2_1500_high_q.mph",
            "sin980_qbic_v8_air_hole_sio2_1500_high_q",
            "high_q",
        )
        optimized: v8.ModeResult | None = None
        if changed(nearest, BASE_NEAREST) or changed(high_q, BASE_HIGH_Q):
            rows = v8.scan_designs(client, local_designs(), OUT / "v8_sio2_1500_local_scan.csv", "v8-sio2-1500")
            best = choose_best(rows)
            opt_design = v8.Design(best.design_id, best.period_nm, best.radius_nm, v8.base.stable_seed(best.design_id))
            optimized = solve_and_save(
                client,
                opt_design,
                "sin980_qbic_v8_air_hole_sio2_1500_optimized.mph",
                "sin980_qbic_v8_air_hole_sio2_1500_optimized",
                "optimized",
            )
            v8.base.write_modes_csv(OUT / "v8_sio2_1500_optimized_selected.csv", [optimized])
        write_summary(nearest, high_q, optimized)
        print(
            "V8_NEAREST_1500 "
            f"P={nearest.period_nm:.3f} R={nearest.radius_nm:.3f} neck={nearest.axis_gap_nm:.3f} "
            f"lambda={nearest.wavelength_nm:.6f} Q={nearest.q_value:.6g}",
            flush=True,
        )
        print(
            "V8_HIGHQ_1500 "
            f"P={high_q.period_nm:.3f} R={high_q.radius_nm:.3f} neck={high_q.axis_gap_nm:.3f} "
            f"lambda={high_q.wavelength_nm:.6f} Q={high_q.q_value:.6g}",
            flush=True,
        )
        if optimized:
            print(
                "V8_OPT_1500 "
                f"P={optimized.period_nm:.3f} R={optimized.radius_nm:.3f} neck={optimized.axis_gap_nm:.3f} "
                f"lambda={optimized.wavelength_nm:.6f} Q={optimized.q_value:.6g}",
                flush=True,
            )
    finally:
        client.clear()


if __name__ == "__main__":
    main()
