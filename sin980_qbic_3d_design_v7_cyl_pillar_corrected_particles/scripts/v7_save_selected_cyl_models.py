from __future__ import annotations

from pathlib import Path

import sin980_qbic_3d_v7_cyl_pillar_corrected_particles as v7


OUT = Path(__file__).resolve().parent.parent / "outputs"


def key_designs() -> list[v7.Design]:
    designs: list[v7.Design] = []
    for period, gap in [
        (594.0, 80.0),
        (594.0, 82.0),
        (594.0, 85.0),
        (596.0, 80.0),
        (596.0, 82.0),
        (596.0, 85.0),
        (600.0, 90.0),
        (600.0, 100.0),
    ]:
        radius = (period - gap) / 2.0
        designs.append(v7.Design(
            f"selected_cyl_P{period:.0f}_G{gap:.0f}",
            period,
            radius,
            v7.stable_seed(f"selected_cyl_{period}_{gap}"),
        ))
    return designs


def choose(rows: list[v7.ModeResult], minimum_gap: float) -> v7.ModeResult:
    pool = [r for r in rows if r.axis_gap_nm >= minimum_gap and abs(r.wavelength_nm - v7.TARGET_NM) <= 5.0]
    if not pool:
        pool = [r for r in rows if r.axis_gap_nm >= minimum_gap]
    return sorted(pool, key=lambda r: (abs(r.wavelength_nm - v7.TARGET_NM), -r.q_value))[0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = v7.import_mph()
    client = mph.Client(cores=4)
    try:
        rows = v7.scan_designs(
            client,
            key_designs(),
            OUT / "clean_selected_cyl_980_scan.csv",
            "selected-cyl",
            include_particles=False,
        )
        best_gap80 = choose(rows, 80.0)
        best_gapgt80 = choose(rows, 80.000001)
        v7.write_modes_csv(OUT / "clean_selected_cyl_best_gap80.csv", [best_gap80])
        v7.write_modes_csv(OUT / "clean_selected_cyl_best_gapgt80.csv", [best_gapgt80])
        saved80 = v7.save_model(
            client,
            best_gap80,
            OUT / "sin980_qbic_v7_cyl_best_gap80_near_980.mph",
            "sin980_qbic_v7_cyl_best_gap80_near_980",
            include_particles=False,
        )
        savedgt80 = v7.save_model(
            client,
            best_gapgt80,
            OUT / "sin980_qbic_v7_cyl_best_gapgt80_near_980.mph",
            "sin980_qbic_v7_cyl_best_gapgt80_near_980",
            include_particles=False,
        )
        print(
            "BEST_GAP80 "
            f"P={saved80.period_nm:.3f} R={saved80.radius_nm:.3f} gap={saved80.axis_gap_nm:.3f} "
            f"lambda={saved80.wavelength_nm:.6f} Q={saved80.q_value:.6g}",
            flush=True,
        )
        print(
            "BEST_GAPGT80 "
            f"P={savedgt80.period_nm:.3f} R={savedgt80.radius_nm:.3f} gap={savedgt80.axis_gap_nm:.3f} "
            f"lambda={savedgt80.wavelength_nm:.6f} Q={savedgt80.q_value:.6g}",
            flush=True,
        )
    finally:
        client.clear()


if __name__ == "__main__":
    main()
