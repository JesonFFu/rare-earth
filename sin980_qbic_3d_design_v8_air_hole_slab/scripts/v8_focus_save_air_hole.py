from __future__ import annotations

from pathlib import Path

import sin980_qbic_3d_v8_air_hole_slab as v8


OUT = Path(__file__).resolve().parent.parent / "outputs"


def focus_designs() -> list[v8.Design]:
    designs: list[v8.Design] = []
    for period in [536.0, 538.0, 540.0, 542.0, 544.0]:
        for radius in [70.0, 72.0, 74.0, 75.6, 78.0, 82.0]:
            design = v8.Design(
                f"hole_focus_P{period:.0f}_R{str(radius).replace('.', 'p')}",
                period,
                radius,
                v8.base.stable_seed(f"hole_focus_{period}_{radius}"),
            )
            if design.axis_gap_nm >= v8.MIN_NECK_NM and 2.0 * radius >= v8.MIN_HOLE_DIAMETER_NM:
                designs.append(design)
    return designs


def choose_best(rows: list[v8.ModeResult]) -> v8.ModeResult:
    near = [r for r in rows if abs(r.wavelength_nm - v8.TARGET_NM) <= 1.0]
    if near:
        return sorted(near, key=lambda r: (-r.q_value, abs(r.wavelength_nm - v8.TARGET_NM)))[0]
    return sorted(rows, key=lambda r: (abs(r.wavelength_nm - v8.TARGET_NM), -r.q_value))[0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = v8.base.import_mph()
    client = mph.Client(cores=4)
    try:
        rows = v8.scan_designs(client, focus_designs(), OUT / "air_hole_focus_scan.csv", "air-hole-focus")
        best = choose_best(rows)
        v8.base.write_modes_csv(OUT / "air_hole_focus_best.csv", [best])
        saved = v8.save_model(
            client,
            best,
            OUT / "sin980_qbic_v8_air_hole_focus_best.mph",
            "sin980_qbic_v8_air_hole_focus_best",
        )
        v8.base.write_modes_csv(OUT / "air_hole_focus_saved_best.csv", [saved])
        v8.write_summary(best, saved, [])
        print(
            "BEST_FOCUS_AIR_HOLE "
            f"P={saved.period_nm:.3f} R={saved.radius_nm:.3f} neck={saved.axis_gap_nm:.3f} "
            f"lambda={saved.wavelength_nm:.6f} Q={saved.q_value:.6g}",
            flush=True,
        )
    finally:
        client.clear()


if __name__ == "__main__":
    main()
