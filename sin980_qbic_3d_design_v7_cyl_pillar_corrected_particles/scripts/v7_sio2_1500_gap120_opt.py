from __future__ import annotations

from pathlib import Path

import sin980_qbic_3d_v7_cyl_pillar_corrected_particles as v7
from v7_sio2_1500_check import add_2d_field_plots


OUT = Path(__file__).resolve().parent.parent / "outputs" / "sio2_1500_gap120_opt"
SIO2_NM = 1500.0
MIN_GAP_NM = 120.000001


def set_thick_sio2() -> None:
    v7.H_SIO2_NM = SIO2_NM


def designs() -> list[v7.Design]:
    out: list[v7.Design] = []
    for period in [610.0, 612.0, 614.0, 616.0, 618.0, 620.0, 622.0, 624.0]:
        for gap in [125.0, 130.0, 140.0, 150.0]:
            radius = (period - gap) / 2.0
            d = v7.Design(
                f"v7_sio2_1500_gap120_P{period:.0f}_G{gap:.0f}",
                period,
                radius,
                v7.stable_seed(f"v7_sio2_1500_gap120_{period}_{gap}"),
            )
            if d.axis_gap_nm >= MIN_GAP_NM and d.radius_nm >= 80.0:
                out.append(d)
    return out


def choose_best(rows: list[v7.ModeResult]) -> v7.ModeResult:
    feasible = [r for r in rows if r.axis_gap_nm >= MIN_GAP_NM]
    near = [r for r in feasible if abs(r.wavelength_nm - v7.TARGET_NM) <= 1.0]
    if near:
        return sorted(near, key=lambda r: (-r.q_value, abs(r.wavelength_nm - v7.TARGET_NM)))[0]
    near3 = [r for r in feasible if abs(r.wavelength_nm - v7.TARGET_NM) <= 3.0]
    if near3:
        return sorted(near3, key=lambda r: (-r.q_value, abs(r.wavelength_nm - v7.TARGET_NM)))[0]
    return sorted(feasible, key=lambda r: (abs(r.wavelength_nm - v7.TARGET_NM), -r.q_value))[0]


def save_model(client, row: v7.ModeResult, filename: str, model_name: str) -> v7.ModeResult:
    set_thick_sio2()
    design = v7.Design(row.design_id, row.period_nm, row.radius_nm, v7.stable_seed(row.design_id))
    model, layout = v7.build_model(client, design, model_name, include_particles=False)
    try:
        rows = v7.solve_design(model, design, layout, 0.0, 0.0)
        selected = min(rows, key=lambda r: abs(r.wavelength_nm - row.wavelength_nm))
        add_2d_field_plots(model, "gap120_opt")
        model.save(str(OUT / filename))
        v7.write_modes_csv(OUT / f"{Path(filename).stem}_modes.csv", rows)
        v7.write_modes_csv(OUT / f"{Path(filename).stem}_selected.csv", [selected])
        return selected
    finally:
        client.remove(model.name())


def write_summary(best: v7.ModeResult) -> None:
    lines = [
        "# v7 Cylinder SiO2 1500 nm Gap >= 120 nm Optimization",
        "",
        f"All files are inside `{OUT}`. Original v7 files were not overwritten.",
        "",
        "## Constraint",
        "- h_sio2_eff = 1500 nm",
        "- cylindrical SiN pillar array",
        "- gap = P - 2R > 120 nm",
        "",
        "## Recommended saved model",
        f"- model = `{OUT / 'sin980_qbic_v7_cyl_sio2_1500_gap120_best.mph'}`",
        f"- P = {best.period_nm:.3f} nm",
        f"- R = {best.radius_nm:.3f} nm",
        f"- gap = {best.axis_gap_nm:.3f} nm",
        f"- wavelength = {best.wavelength_nm:.6f} nm",
        f"- Q = {best.q_value:.6g}",
        "",
        "## Field plots",
        "- The saved .mph includes two 2D `|E|^2` slice plot groups: xy pillar mid-plane and xz center-plane.",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    set_thick_sio2()
    mph = v7.import_mph()
    client = mph.Client(cores=4)
    try:
        rows = v7.scan_designs(client, designs(), OUT / "v7_sio2_1500_gap120_scan.csv", "v7-gap120", include_particles=False)
        best = choose_best(rows)
        v7.write_modes_csv(OUT / "v7_sio2_1500_gap120_best.csv", [best])
        saved = save_model(
            client,
            best,
            "sin980_qbic_v7_cyl_sio2_1500_gap120_best.mph",
            "sin980_qbic_v7_cyl_sio2_1500_gap120_best",
        )
        v7.write_modes_csv(OUT / "v7_sio2_1500_gap120_saved_best.csv", [saved])
        write_summary(saved)
        print(
            "V7_GAP120_BEST "
            f"P={saved.period_nm:.3f} R={saved.radius_nm:.3f} gap={saved.axis_gap_nm:.3f} "
            f"lambda={saved.wavelength_nm:.6f} Q={saved.q_value:.6g}",
            flush=True,
        )
    finally:
        client.clear()


if __name__ == "__main__":
    main()
