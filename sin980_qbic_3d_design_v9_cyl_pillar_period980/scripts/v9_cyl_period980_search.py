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


def period980_designs() -> list[v7.Design]:
    designs: list[v7.Design] = []
    for period in [900.0, 940.0, 960.0, 980.0, 1000.0, 1020.0, 1060.0]:
        for ratio in [0.14, 0.20, 0.26, 0.32, 0.38, 0.42, 0.45]:
            radius = period * ratio
            design = v7.Design(
                f"period980_P{period:.0f}_rho{str(ratio).replace('.', 'p')}",
                period,
                radius,
                v7.stable_seed(f"period980_{period}_{ratio}"),
            )
            if design.axis_gap_nm >= MIN_GAP_NM and design.radius_nm >= 80.0:
                designs.append(design)
    return designs


def choose_near_target(rows: list[v7.ModeResult], window_nm: float = 10.0) -> v7.ModeResult:
    practical = [r for r in rows if r.axis_gap_nm >= MIN_GAP_NM]
    near = [r for r in practical if abs(r.wavelength_nm - v7.TARGET_NM) <= window_nm]
    if near:
        return sorted(near, key=lambda r: (-r.q_value, abs(r.wavelength_nm - v7.TARGET_NM)))[0]
    return sorted(practical, key=lambda r: (abs(r.wavelength_nm - v7.TARGET_NM), -r.q_value))[0]


def choose_highest_q(rows: list[v7.ModeResult]) -> v7.ModeResult:
    practical = [r for r in rows if r.axis_gap_nm >= MIN_GAP_NM]
    return sorted(practical, key=lambda r: -r.q_value)[0]


def focus_designs(seed: v7.ModeResult) -> list[v7.Design]:
    designs: list[v7.Design] = []
    seen: set[tuple[float, float]] = set()
    for dp in [-20.0, -10.0, 0.0, 10.0, 20.0]:
        for dr in [-18.0, -9.0, 0.0, 9.0, 18.0]:
            period = round(seed.period_nm + dp, 3)
            radius = round(seed.radius_nm + dr, 3)
            key = (period, radius)
            if key in seen:
                continue
            seen.add(key)
            design = v7.Design(
                f"period980_focus_P{period:.0f}_R{radius:.0f}",
                period,
                radius,
                v7.stable_seed(f"period980_focus_{period}_{radius}"),
            )
            if design.axis_gap_nm >= MIN_GAP_NM and design.radius_nm >= 80.0:
                designs.append(design)
    return designs


def save_named(client, row: v7.ModeResult, filename: str, model_name: str) -> v7.ModeResult:
    return v7.save_model(client, row, OUT / filename, model_name, include_particles=False)


def write_summary(best_near: v7.ModeResult, saved_near: v7.ModeResult, best_q: v7.ModeResult, saved_q: v7.ModeResult) -> None:
    same = (
        abs(saved_near.period_nm - saved_q.period_nm) < 1e-9
        and abs(saved_near.radius_nm - saved_q.radius_nm) < 1e-9
        and saved_near.mode_index == saved_q.mode_index
    )
    lines = [
        "# v9 Period-Scale SiN Cylindrical Pillar Search",
        "",
        f"All files were generated inside `{ROOT}`. Older versions were not modified.",
        "",
        "## Geometry",
        "- 3D square lattice of 300 nm high SiN cylindrical pillars on SiO2/Si.",
        "- Period is scanned near 980 nm; this is not an air-hole slab.",
        "- Practical filter: `P - 2R >= 80 nm`.",
        "",
        "## Best near 980 nm",
        f"- model = `{OUT / 'sin980_qbic_v9_cyl_period980_best_near_980.mph'}`",
        f"- P = {saved_near.period_nm:.3f} nm",
        f"- R = {saved_near.radius_nm:.3f} nm",
        f"- gap = {saved_near.axis_gap_nm:.3f} nm",
        f"- wavelength = {saved_near.wavelength_nm:.6f} nm",
        f"- Q = {saved_near.q_value:.6g}",
        "",
        "## Highest-Q saved candidate",
    ]
    if same:
        lines.append("- Same as the near-980 model.")
    else:
        lines.extend([
            f"- model = `{OUT / 'sin980_qbic_v9_cyl_period980_highest_q.mph'}`",
            f"- P = {saved_q.period_nm:.3f} nm",
            f"- R = {saved_q.radius_nm:.3f} nm",
            f"- gap = {saved_q.axis_gap_nm:.3f} nm",
            f"- wavelength = {saved_q.wavelength_nm:.6f} nm",
            f"- Q = {saved_q.q_value:.6g}",
        ])
    lines.extend([
        "",
        "## Interpretation",
    ])
    if saved_near.q_value >= 1e9 or saved_q.q_value >= 1e9:
        lines.append("- A `Q >= 1e9` long-period pillar candidate was found.")
    else:
        lines.append("- No `Q >= 1e9` mode was found in this period-near-980 pillar scan.")
    lines.append("- Compare this directly with the v8 air-hole slab, which reached about `2e5` Q near 980 nm.")
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = v7.import_mph()
    client = mph.Client(cores=4)
    try:
        coarse = v7.scan_designs(client, period980_designs(), OUT / "period980_cyl_coarse_scan.csv", "period980-cyl", include_particles=False)
        near_seed = choose_near_target(coarse)
        high_seed = choose_highest_q(coarse)
        v7.write_modes_csv(OUT / "period980_cyl_coarse_best_near_980.csv", [near_seed])
        v7.write_modes_csv(OUT / "period980_cyl_coarse_highest_q.csv", [high_seed])
        focus_seed = near_seed if abs(near_seed.wavelength_nm - v7.TARGET_NM) <= 15.0 else high_seed
        fine = v7.scan_designs(client, focus_designs(focus_seed), OUT / "period980_cyl_fine_scan.csv", "period980-cyl-fine", include_particles=False)
        all_rows = coarse + fine
        best_near = choose_near_target(all_rows)
        best_q = choose_highest_q(all_rows)
        v7.write_modes_csv(OUT / "period980_cyl_best_near_980.csv", [best_near])
        v7.write_modes_csv(OUT / "period980_cyl_highest_q.csv", [best_q])
        saved_near = save_named(client, best_near, "sin980_qbic_v9_cyl_period980_best_near_980.mph", "sin980_qbic_v9_cyl_period980_best_near_980")
        if best_q.design_id == best_near.design_id and best_q.mode_index == best_near.mode_index:
            saved_q = saved_near
        else:
            saved_q = save_named(client, best_q, "sin980_qbic_v9_cyl_period980_highest_q.mph", "sin980_qbic_v9_cyl_period980_highest_q")
        v7.write_modes_csv(OUT / "period980_cyl_saved_best_near_980.csv", [saved_near])
        v7.write_modes_csv(OUT / "period980_cyl_saved_highest_q.csv", [saved_q])
        write_summary(best_near, saved_near, best_q, saved_q)
        print(
            "BEST_NEAR_980 "
            f"P={saved_near.period_nm:.3f} R={saved_near.radius_nm:.3f} gap={saved_near.axis_gap_nm:.3f} "
            f"lambda={saved_near.wavelength_nm:.6f} Q={saved_near.q_value:.6g}",
            flush=True,
        )
        print(
            "HIGHEST_Q "
            f"P={saved_q.period_nm:.3f} R={saved_q.radius_nm:.3f} gap={saved_q.axis_gap_nm:.3f} "
            f"lambda={saved_q.wavelength_nm:.6f} Q={saved_q.q_value:.6g}",
            flush=True,
        )
    finally:
        client.clear()


if __name__ == "__main__":
    main()
