from __future__ import annotations

from pathlib import Path

import sin980_qbic_3d_v8_air_hole_slab as v8


OUT = Path(__file__).resolve().parent.parent / "outputs"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    design = v8.Design(
        "hole_nearest980_P540_R75p6",
        540.0,
        75.6,
        v8.base.stable_seed("hole_nearest980_P540_R75p6"),
    )
    mph = v8.base.import_mph()
    client = mph.Client(cores=4)
    try:
        model = v8.build_model(client, design, "sin980_qbic_v8_air_hole_nearest_980")
        try:
            rows = v8.solve_design(model, design, 0.0, 0.0)
            best = v8.best_mode(rows)
            v8.base.add_field_plot_groups(model, "air_hole_nearest_980")
            model.save(str(OUT / "sin980_qbic_v8_air_hole_nearest_980.mph"))
            v8.base.write_modes_csv(OUT / "air_hole_nearest_980_saved.csv", [best])
            print(
                "SAVED_NEAREST_980 "
                f"P={best.period_nm:.3f} R={best.radius_nm:.3f} neck={best.axis_gap_nm:.3f} "
                f"lambda={best.wavelength_nm:.6f} Q={best.q_value:.6g}",
                flush=True,
            )
        finally:
            client.remove(model.name())
    finally:
        client.clear()


if __name__ == "__main__":
    main()
