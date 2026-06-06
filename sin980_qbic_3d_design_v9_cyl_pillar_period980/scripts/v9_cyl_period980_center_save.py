from __future__ import annotations

from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
OUT = ROOT / "outputs"
V7_SCRIPT_DIR = ROOT.parent / "sin980_qbic_3d_design_v7_cyl_pillar_corrected_particles" / "scripts"
sys.path.insert(0, str(V7_SCRIPT_DIR))

import sin980_qbic_3d_v7_cyl_pillar_corrected_particles as v7


MIN_GAP_NM = 80.0


def center_designs() -> list[v7.Design]:
    designs: list[v7.Design] = []
    for period in [960.0, 980.0, 1000.0]:
        for ratio in [0.20, 0.28, 0.36, 0.42, 0.45]:
            radius = period * ratio
            design = v7.Design(
                f"center_P{period:.0f}_rho{str(ratio).replace('.', 'p')}",
                period,
                radius,
                v7.stable_seed(f"center_period980_{period}_{ratio}"),
            )
            if design.axis_gap_nm >= MIN_GAP_NM and design.radius_nm >= 80.0:
                designs.append(design)
    return designs


def choose_near(rows: list[v7.ModeResult]) -> v7.ModeResult:
    practical = [r for r in rows if r.axis_gap_nm >= MIN_GAP_NM]
    near = [r for r in practical if abs(r.wavelength_nm - v7.TARGET_NM) <= 20.0]
    pool = near or practical
    return sorted(pool, key=lambda r: (abs(r.wavelength_nm - v7.TARGET_NM), -r.q_value))[0]


def choose_highq(rows: list[v7.ModeResult]) -> v7.ModeResult:
    practical = [r for r in rows if r.axis_gap_nm >= MIN_GAP_NM]
    return sorted(practical, key=lambda r: -r.q_value)[0]


def write_summary(saved_near: v7.ModeResult, saved_highq: v7.ModeResult) -> None:
    lines = [
        "# v9 Centered Period-980 Cylinder Check",
        "",
        f"All files were generated inside `{ROOT}`.",
        "",
        "## Search",
        "- Geometry: 300 nm SiN cylindrical pillars on SiO2/Si.",
        "- Periods tested: 960, 980, 1000 nm.",
        "- Radius ratios tested: 0.20, 0.28, 0.36, 0.42, 0.45.",
        "- Practical filter: `P - 2R >= 80 nm`.",
        "",
        "## Best near 980 nm",
        f"- model = `{OUT / 'sin980_qbic_v9_cyl_period980_center_best_near_980.mph'}`",
        f"- P = {saved_near.period_nm:.3f} nm",
        f"- R = {saved_near.radius_nm:.3f} nm",
        f"- gap = {saved_near.axis_gap_nm:.3f} nm",
        f"- wavelength = {saved_near.wavelength_nm:.6f} nm",
        f"- Q = {saved_near.q_value:.6g}",
        "",
        "## Highest-Q candidate",
        f"- model = `{OUT / 'sin980_qbic_v9_cyl_period980_center_highest_q.mph'}`",
        f"- P = {saved_highq.period_nm:.3f} nm",
        f"- R = {saved_highq.radius_nm:.3f} nm",
        f"- gap = {saved_highq.axis_gap_nm:.3f} nm",
        f"- wavelength = {saved_highq.wavelength_nm:.6f} nm",
        f"- Q = {saved_highq.q_value:.6g}",
        "",
        "## Interpretation",
    ]
    if max(saved_near.q_value, saved_highq.q_value) >= 1e9:
        lines.append("- A `Q >= 1e9` candidate was found.")
    else:
        lines.append("- No high-Q BIC emerged in this period-near-980 pillar check.")
    lines.append("- In this pass, the period-near-980 pillar geometry remains weaker than the v8 air-hole slab.")
    (OUT / "summary_center.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = v7.import_mph()
    client = mph.Client(cores=4)
    try:
        rows = v7.scan_designs(client, center_designs(), OUT / "period980_cyl_center_scan.csv", "period980-center", include_particles=False)
        best_near = choose_near(rows)
        best_q = choose_highq(rows)
        v7.write_modes_csv(OUT / "period980_cyl_center_best_near_980.csv", [best_near])
        v7.write_modes_csv(OUT / "period980_cyl_center_highest_q.csv", [best_q])
        saved_near = v7.save_model(
            client,
            best_near,
            OUT / "sin980_qbic_v9_cyl_period980_center_best_near_980.mph",
            "sin980_qbic_v9_cyl_period980_center_best_near_980",
            include_particles=False,
        )
        if best_q.design_id == best_near.design_id and best_q.mode_index == best_near.mode_index:
            saved_highq = saved_near
            v7.write_modes_csv(OUT / "period980_cyl_center_saved_highest_q.csv", [saved_highq])
        else:
            saved_highq = v7.save_model(
                client,
                best_q,
                OUT / "sin980_qbic_v9_cyl_period980_center_highest_q.mph",
                "sin980_qbic_v9_cyl_period980_center_highest_q",
                include_particles=False,
            )
        v7.write_modes_csv(OUT / "period980_cyl_center_saved_best_near_980.csv", [saved_near])
        v7.write_modes_csv(OUT / "period980_cyl_center_saved_highest_q.csv", [saved_highq])
        write_summary(saved_near, saved_highq)
        print(
            "BEST_NEAR_980 "
            f"P={saved_near.period_nm:.3f} R={saved_near.radius_nm:.3f} gap={saved_near.axis_gap_nm:.3f} "
            f"lambda={saved_near.wavelength_nm:.6f} Q={saved_near.q_value:.6g}",
            flush=True,
        )
        print(
            "HIGHEST_Q "
            f"P={saved_highq.period_nm:.3f} R={saved_highq.radius_nm:.3f} gap={saved_highq.axis_gap_nm:.3f} "
            f"lambda={saved_highq.wavelength_nm:.6f} Q={saved_highq.q_value:.6g}",
            flush=True,
        )
    finally:
        client.clear()


if __name__ == "__main__":
    main()
