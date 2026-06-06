from __future__ import annotations

import argparse
import math
from pathlib import Path

import sin980_qbic_3d_v7_cyl_pillar_corrected_particles as v7


OUT = Path(__file__).resolve().parent.parent / "outputs"


def focused_designs() -> list[v7.Design]:
    designs: list[v7.Design] = []
    seen: set[tuple[float, float]] = set()
    for period in [592.0, 594.0, 596.0, 598.0, 600.0, 602.0, 604.0]:
        for gap in [80.0, 82.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0]:
            radius = (period - gap) / 2.0
            key = (period, radius)
            if key in seen:
                continue
            seen.add(key)
            design = v7.Design(
                f"gap80fine_P{period:.0f}_G{gap:.0f}".replace(".", "p"),
                period,
                radius,
                v7.stable_seed(f"gap80fine_{period}_{gap}"),
            )
            if design.axis_gap_nm >= 80.0 and design.radius_nm >= 80.0:
                designs.append(design)
    return designs


def best_near_980(rows: list[v7.ModeResult]) -> v7.ModeResult:
    practical = [r for r in rows if r.axis_gap_nm >= 80.0]
    near = [r for r in practical if abs(r.wavelength_nm - v7.TARGET_NM) <= 3.0]
    pool = near or practical
    return sorted(pool, key=lambda r: (abs(r.wavelength_nm - v7.TARGET_NM), -r.q_value))[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-model", action="store_true")
    parser.add_argument("--band", action="store_true")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    mph = v7.import_mph()
    client = mph.Client(cores=4)
    try:
        rows = v7.scan_designs(
            client,
            focused_designs(),
            OUT / "clean_gap80_fine_scan.csv",
            "gap80-fine",
            include_particles=False,
        )
        best_q = sorted([r for r in rows if r.axis_gap_nm >= 80.0], key=lambda r: -r.q_value)[0]
        best_980 = best_near_980(rows)
        v7.write_modes_csv(OUT / "clean_gap80_fine_best_by_q.csv", [best_q])
        v7.write_modes_csv(OUT / "clean_gap80_fine_best_near_980.csv", [best_980])
        print(
            "BEST_Q "
            f"P={best_q.period_nm:.3f} R={best_q.radius_nm:.3f} gap={best_q.axis_gap_nm:.3f} "
            f"lambda={best_q.wavelength_nm:.6f} Q={best_q.q_value:.6g}",
            flush=True,
        )
        print(
            "BEST_NEAR_980 "
            f"P={best_980.period_nm:.3f} R={best_980.radius_nm:.3f} gap={best_980.axis_gap_nm:.3f} "
            f"lambda={best_980.wavelength_nm:.6f} Q={best_980.q_value:.6g}",
            flush=True,
        )
        if args.save_model:
            saved = v7.save_model(
                client,
                best_980,
                OUT / "sin980_qbic_v7_cyl_gap80_best_near_980.mph",
                "sin980_qbic_v7_cyl_gap80_best_near_980",
                include_particles=False,
            )
            print(
                "SAVED_MODEL "
                f"lambda={saved.wavelength_nm:.6f} Q={saved.q_value:.6g}",
                flush=True,
            )
        if args.band:
            band = v7.run_band(client, best_980)
            print(f"BAND_POINTS {len(band)}", flush=True)
    finally:
        client.clear()


if __name__ == "__main__":
    main()
