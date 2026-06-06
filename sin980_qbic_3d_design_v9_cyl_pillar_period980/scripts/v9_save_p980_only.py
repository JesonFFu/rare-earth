from __future__ import annotations

from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
OUT = ROOT / "outputs"
V7_SCRIPT_DIR = ROOT.parent / "sin980_qbic_3d_design_v7_cyl_pillar_corrected_particles" / "scripts"
sys.path.insert(0, str(V7_SCRIPT_DIR))

import sin980_qbic_3d_v7_cyl_pillar_corrected_particles as v7


def designs() -> list[v7.Design]:
    period = 980.0
    out: list[v7.Design] = []
    for ratio in [0.28, 0.36, 0.42, 0.45]:
        radius = period * ratio
        out.append(v7.Design(
            f"p980_only_rho{str(ratio).replace('.', 'p')}",
            period,
            radius,
            v7.stable_seed(f"p980_only_{ratio}"),
        ))
    return out


def choose_best_near(rows: list[v7.ModeResult]) -> v7.ModeResult:
    near = [r for r in rows if abs(r.wavelength_nm - v7.TARGET_NM) <= 30.0]
    pool = near or rows
    return sorted(pool, key=lambda r: (abs(r.wavelength_nm - v7.TARGET_NM), -r.q_value))[0]


def choose_highest_q(rows: list[v7.ModeResult]) -> v7.ModeResult:
    return sorted(rows, key=lambda r: -r.q_value)[0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = v7.import_mph()
    client = mph.Client(cores=4)
    try:
        rows = v7.scan_designs(client, designs(), OUT / "p980_only_cyl_scan.csv", "p980-only-cyl", include_particles=False)
        best_near = choose_best_near(rows)
        best_q = choose_highest_q(rows)
        v7.write_modes_csv(OUT / "p980_only_best_near_980.csv", [best_near])
        v7.write_modes_csv(OUT / "p980_only_highest_q.csv", [best_q])
        saved_near = v7.save_model(
            client,
            best_near,
            OUT / "sin980_qbic_v9_p980_cyl_best_near_980.mph",
            "sin980_qbic_v9_p980_cyl_best_near_980",
            include_particles=False,
        )
        if best_q.design_id == best_near.design_id and best_q.mode_index == best_near.mode_index:
            saved_q = saved_near
        else:
            saved_q = v7.save_model(
                client,
                best_q,
                OUT / "sin980_qbic_v9_p980_cyl_highest_q.mph",
                "sin980_qbic_v9_p980_cyl_highest_q",
                include_particles=False,
            )
        v7.write_modes_csv(OUT / "p980_only_saved_best_near_980.csv", [saved_near])
        v7.write_modes_csv(OUT / "p980_only_saved_highest_q.csv", [saved_q])
        lines = [
            "# v9 P=980 nm SiN Cylinder Check",
            "",
            f"All files were generated inside `{ROOT}`.",
            "",
            "## P=980 nm pillar scan",
            "- Radius ratios: 0.28, 0.36, 0.42, 0.45.",
            "",
            "## Best near 980 nm",
            f"- model = `{OUT / 'sin980_qbic_v9_p980_cyl_best_near_980.mph'}`",
            f"- P = {saved_near.period_nm:.3f} nm",
            f"- R = {saved_near.radius_nm:.3f} nm",
            f"- gap = {saved_near.axis_gap_nm:.3f} nm",
            f"- wavelength = {saved_near.wavelength_nm:.6f} nm",
            f"- Q = {saved_near.q_value:.6g}",
            "",
            "## Highest-Q candidate",
            f"- P = {saved_q.period_nm:.3f} nm",
            f"- R = {saved_q.radius_nm:.3f} nm",
            f"- gap = {saved_q.axis_gap_nm:.3f} nm",
            f"- wavelength = {saved_q.wavelength_nm:.6f} nm",
            f"- Q = {saved_q.q_value:.6g}",
            "",
            "## Interpretation",
        ]
        if max(saved_near.q_value, saved_q.q_value) >= 1e9:
            lines.append("- A `Q >= 1e9` P=980 nm pillar candidate was found.")
        else:
            lines.append("- No high-Q BIC candidate appeared for P=980 nm pillars in this scan.")
        (OUT / "summary_p980_only.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(
            "P980_BEST_NEAR "
            f"P={saved_near.period_nm:.3f} R={saved_near.radius_nm:.3f} gap={saved_near.axis_gap_nm:.3f} "
            f"lambda={saved_near.wavelength_nm:.6f} Q={saved_near.q_value:.6g}",
            flush=True,
        )
        print(
            "P980_HIGHEST_Q "
            f"P={saved_q.period_nm:.3f} R={saved_q.radius_nm:.3f} gap={saved_q.axis_gap_nm:.3f} "
            f"lambda={saved_q.wavelength_nm:.6f} Q={saved_q.q_value:.6g}",
            flush=True,
        )
    finally:
        client.clear()


if __name__ == "__main__":
    main()
