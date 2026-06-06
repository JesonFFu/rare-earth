from __future__ import annotations

import csv
import math
from pathlib import Path

import sin980_qbic_3d_v7_cyl_pillar_corrected_particles as v7


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs" / "sio2_1500_gap120_opt" / "band_proof"
SOURCE_MPH = ROOT / "outputs" / "sio2_1500_gap120_opt" / "sin980_qbic_v7_cyl_sio2_1500_gap120_best.mph"
TARGET_NM = 980.0


DESIGN = v7.Design(
    "v7_cyl_sio2_1500_gap120_best_band",
    616.0,
    243.0,
    v7.stable_seed("v7_cyl_sio2_1500_gap120_best_band"),
)
LAYOUT = v7.ParticleLayout((), (), ())


def band_path() -> list[tuple[str, float, float]]:
    return [
        ("G", 0.0, 0.0),
        ("G-X", 0.02, 0.0),
        ("G-X", 0.05, 0.0),
        ("G-X", 0.10, 0.0),
        ("G-X", 0.20, 0.0),
        ("G-X", 0.35, 0.0),
        ("X", 0.50, 0.0),
        ("X-M", 0.50, 0.10),
        ("X-M", 0.50, 0.25),
        ("M", 0.50, 0.50),
        ("M-G", 0.35, 0.35),
        ("M-G", 0.20, 0.20),
        ("M-G", 0.10, 0.10),
        ("M-G", 0.05, 0.05),
        ("G", 0.0, 0.0),
    ]


def q_gamma_path() -> list[tuple[str, float, float]]:
    return [
        ("G", 0.0, 0.0),
        ("G-X", 0.002, 0.0),
        ("G-X", 0.005, 0.0),
        ("G-X", 0.010, 0.0),
        ("G-X", 0.020, 0.0),
        ("G-X", 0.050, 0.0),
        ("G-X", 0.080, 0.0),
        ("G-X", 0.100, 0.0),
    ]


def cumulative_s(points: list[tuple[str, float, float]]) -> list[float]:
    s = [0.0]
    for (_, x0, y0), (_, x1, y1) in zip(points, points[1:]):
        s.append(s[-1] + math.hypot(x1 - x0, y1 - y0))
    return s


def solve_cached(model, cache: dict[tuple[float, float], list[v7.ModeResult]], kx: float, ky: float) -> list[v7.ModeResult]:
    key = (round(kx, 9), round(ky, 9))
    if key not in cache:
        cache[key] = v7.solve_design(model, DESIGN, LAYOUT, kx, ky)
    return cache[key]


def select_branch(rows: list[v7.ModeResult], previous_lambda: float | None) -> v7.ModeResult:
    if previous_lambda is None:
        return min(rows, key=lambda r: abs(r.wavelength_nm - TARGET_NM))
    return min(rows, key=lambda r: abs(r.wavelength_nm - previous_lambda))


def write_all_modes(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "path_index", "segment", "s_norm", "kx_norm", "ky_norm", "mode_index",
        "wavelength_nm", "freq_thz_real", "freq_thz_imag", "q_value",
        "period_nm", "radius_nm", "axis_gap_nm",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_branch(path: Path, rows: list[dict[str, object]]) -> None:
    write_all_modes(path, rows)


def svg_plot_band(path: Path, all_rows: list[dict[str, object]], branch_rows: list[dict[str, object]], title: str) -> None:
    width, height = 980, 560
    margin_left, margin_right, margin_top, margin_bottom = 76, 30, 42, 74
    xs = [float(r["s_norm"]) for r in all_rows]
    ys = [float(r["wavelength_nm"]) for r in all_rows]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    pad = max(1.0, 0.05 * (ymax - ymin))
    ymin -= pad
    ymax += pad

    def sx(x: float) -> float:
        return margin_left + (x - xmin) / (xmax - xmin) * (width - margin_left - margin_right)

    def sy(y: float) -> float:
        return height - margin_bottom - (y - ymin) / (ymax - ymin) * (height - margin_top - margin_bottom)

    branch_points = " ".join(f"{sx(float(r['s_norm'])):.2f},{sy(float(r['wavelength_nm'])):.2f}" for r in branch_rows)
    circles = []
    for r in all_rows:
        circles.append(
            f"<circle cx='{sx(float(r['s_norm'])):.2f}' cy='{sy(float(r['wavelength_nm'])):.2f}' r='3.0' fill='#9aa6b2' opacity='0.55'/>"
        )
    branch_circles = []
    for r in branch_rows:
        branch_circles.append(
            f"<circle cx='{sx(float(r['s_norm'])):.2f}' cy='{sy(float(r['wavelength_nm'])):.2f}' r='4.5' fill='#c93232'/>"
        )
    tick_labels = [("Gamma", 0.0), ("X", 0.5), ("M", 1.0), ("Gamma", max(xs))]
    vlines = []
    for label, xpos in tick_labels:
        if xmin <= xpos <= xmax:
            x = sx(xpos)
            vlines.append(f"<line x1='{x:.2f}' y1='{margin_top}' x2='{x:.2f}' y2='{height-margin_bottom}' stroke='#d3d7de' stroke-dasharray='4 5'/>")
            vlines.append(f"<text x='{x:.2f}' y='{height-34}' text-anchor='middle' font-size='18'>{label}</text>")
    ygrid = []
    for i in range(6):
        y = ymin + i * (ymax - ymin) / 5
        yy = sy(y)
        ygrid.append(f"<line x1='{margin_left}' y1='{yy:.2f}' x2='{width-margin_right}' y2='{yy:.2f}' stroke='#edf0f4'/>")
        ygrid.append(f"<text x='{margin_left-12}' y='{yy+5:.2f}' text-anchor='end' font-size='14'>{y:.1f}</text>")

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2}" y="28" text-anchor="middle" font-size="22" font-family="Arial">{title}</text>
{''.join(ygrid)}
{''.join(vlines)}
<line x1="{margin_left}" y1="{height-margin_bottom}" x2="{width-margin_right}" y2="{height-margin_bottom}" stroke="#20242a" stroke-width="1.5"/>
<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height-margin_bottom}" stroke="#20242a" stroke-width="1.5"/>
<text x="{width/2}" y="{height-8}" text-anchor="middle" font-size="16" font-family="Arial">k path: Gamma-X-M-Gamma</text>
<text x="20" y="{height/2}" text-anchor="middle" font-size="16" font-family="Arial" transform="rotate(-90 20 {height/2})">wavelength (nm)</text>
{''.join(circles)}
<polyline points="{branch_points}" fill="none" stroke="#c93232" stroke-width="3"/>
{''.join(branch_circles)}
<rect x="{width-250}" y="54" width="214" height="58" fill="white" stroke="#d7dce2"/>
<circle cx="{width-230}" cy="76" r="4" fill="#9aa6b2" opacity="0.65"/><text x="{width-216}" y="81" font-size="14">all solved modes</text>
<line x1="{width-238}" y1="99" x2="{width-222}" y2="99" stroke="#c93232" stroke-width="3"/><text x="{width-216}" y="104" font-size="14">tracked 980 nm branch</text>
</svg>"""
    path.write_text(svg, encoding="utf-8")


def svg_plot_q(path: Path, rows: list[dict[str, object]], title: str) -> None:
    width, height = 900, 520
    ml, mr, mt, mb = 78, 28, 42, 72
    xs = [float(r["kx_norm"]) for r in rows]
    qs = [max(float(r["q_value"]), 1.0) for r in rows]
    ys = [math.log10(q) for q in qs]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    ymin = max(0.0, ymin - 0.4)
    ymax += 0.4

    def sx(x: float) -> float:
        return ml + (x - xmin) / (xmax - xmin) * (width - ml - mr)

    def sy(y: float) -> float:
        return height - mb - (y - ymin) / (ymax - ymin) * (height - mt - mb)

    points = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(xs, ys))
    circles = "".join(f"<circle cx='{sx(x):.2f}' cy='{sy(y):.2f}' r='4.5' fill='#1f6feb'/>" for x, y in zip(xs, ys))
    ygrid = []
    for exp in range(math.floor(ymin), math.ceil(ymax) + 1):
        yy = sy(exp)
        ygrid.append(f"<line x1='{ml}' y1='{yy:.2f}' x2='{width-mr}' y2='{yy:.2f}' stroke='#edf0f4'/>")
        ygrid.append(f"<text x='{ml-12}' y='{yy+5:.2f}' text-anchor='end' font-size='14'>1e{exp}</text>")
    xgrid = []
    for x in xs:
        xx = sx(x)
        xgrid.append(f"<text x='{xx:.2f}' y='{height-38}' text-anchor='middle' font-size='13'>{x:.3g}</text>")

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2}" y="28" text-anchor="middle" font-size="22" font-family="Arial">{title}</text>
{''.join(ygrid)}
<line x1="{ml}" y1="{height-mb}" x2="{width-mr}" y2="{height-mb}" stroke="#20242a" stroke-width="1.5"/>
<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{height-mb}" stroke="#20242a" stroke-width="1.5"/>
<polyline points="{points}" fill="none" stroke="#1f6feb" stroke-width="3"/>
{circles}
{''.join(xgrid)}
<text x="{width/2}" y="{height-8}" text-anchor="middle" font-size="16" font-family="Arial">kx normalized to pi/P along Gamma-X</text>
<text x="22" y="{height/2}" text-anchor="middle" font-size="16" font-family="Arial" transform="rotate(-90 22 {height/2})">Q factor (log scale)</text>
</svg>"""
    path.write_text(svg, encoding="utf-8")


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = v7.import_mph()
    client = mph.Client(cores=4)
    cache: dict[tuple[float, float], list[v7.ModeResult]] = {}
    model = client.load(str(SOURCE_MPH))
    try:
        all_rows: list[dict[str, object]] = []
        branch_rows: list[dict[str, object]] = []
        path = band_path()
        svals = cumulative_s(path)
        previous_lambda: float | None = None
        for idx, ((segment, kx, ky), sval) in enumerate(zip(path, svals)):
            print(f"[band] {idx+1}/{len(path)} {segment} k=({kx:.4g},{ky:.4g})", flush=True)
            modes = solve_cached(model, cache, kx, ky)
            for mode in modes:
                all_rows.append({
                    "path_index": idx,
                    "segment": segment,
                    "s_norm": sval,
                    "kx_norm": kx,
                    "ky_norm": ky,
                    "mode_index": mode.mode_index,
                    "wavelength_nm": mode.wavelength_nm,
                    "freq_thz_real": mode.freq_thz_real,
                    "freq_thz_imag": mode.freq_thz_imag,
                    "q_value": mode.q_value,
                    "period_nm": mode.period_nm,
                    "radius_nm": mode.radius_nm,
                    "axis_gap_nm": mode.axis_gap_nm,
                })
            selected = select_branch(modes, previous_lambda)
            previous_lambda = selected.wavelength_nm
            branch_rows.append({
                "path_index": idx,
                "segment": segment,
                "s_norm": sval,
                "kx_norm": kx,
                "ky_norm": ky,
                "mode_index": selected.mode_index,
                "wavelength_nm": selected.wavelength_nm,
                "freq_thz_real": selected.freq_thz_real,
                "freq_thz_imag": selected.freq_thz_imag,
                "q_value": selected.q_value,
                "period_nm": selected.period_nm,
                "radius_nm": selected.radius_nm,
                "axis_gap_nm": selected.axis_gap_nm,
            })
            print(f"  branch lambda={selected.wavelength_nm:.6f} nm Q={selected.q_value:.4g}", flush=True)

        q_rows: list[dict[str, object]] = []
        previous_lambda = None
        for idx, (segment, kx, ky) in enumerate(q_gamma_path()):
            print(f"[q-near-gamma] {idx+1}/{len(q_gamma_path())} k=({kx:.4g},{ky:.4g})", flush=True)
            modes = solve_cached(model, cache, kx, ky)
            selected = select_branch(modes, previous_lambda)
            previous_lambda = selected.wavelength_nm
            q_rows.append({
                "path_index": idx,
                "segment": segment,
                "s_norm": kx,
                "kx_norm": kx,
                "ky_norm": ky,
                "mode_index": selected.mode_index,
                "wavelength_nm": selected.wavelength_nm,
                "freq_thz_real": selected.freq_thz_real,
                "freq_thz_imag": selected.freq_thz_imag,
                "q_value": selected.q_value,
                "period_nm": selected.period_nm,
                "radius_nm": selected.radius_nm,
                "axis_gap_nm": selected.axis_gap_nm,
            })
            print(f"  q-branch lambda={selected.wavelength_nm:.6f} nm Q={selected.q_value:.4g}", flush=True)

        write_all_modes(OUT / "band_all_modes.csv", all_rows)
        write_branch(OUT / "band_tracked_branch.csv", branch_rows)
        write_branch(OUT / "q_near_gamma.csv", q_rows)
        svg_plot_band(OUT / "band_tracked_branch.svg", all_rows, branch_rows, "v7 SiN pillar gap>120 nm band near 980 nm")
        svg_plot_q(OUT / "q_near_gamma.svg", q_rows, "Q(k) near Gamma for tracked 980 nm branch")
        summary = [
            "# v7 gap>120 nm cylinder band / BIC proof",
            "",
            f"Source model: `{SOURCE_MPH}`",
            "",
            "## Geometry",
            "- P = 616 nm",
            "- R = 243 nm",
            "- gap = 130 nm",
            "- h_sio2_eff = 1500 nm",
            "",
            "## Outputs",
            f"- Band CSV: `{OUT / 'band_tracked_branch.csv'}`",
            f"- All solved modes CSV: `{OUT / 'band_all_modes.csv'}`",
            f"- Near-Gamma Q CSV: `{OUT / 'q_near_gamma.csv'}`",
            f"- Band plot: `{OUT / 'band_tracked_branch.svg'}`",
            f"- Q(k) plot: `{OUT / 'q_near_gamma.svg'}`",
            "",
            "## Interpretation rule",
            "- A BIC-like mode should show the largest Q at Gamma and a strong Q decrease after moving away from Gamma.",
        ]
        gamma = q_rows[0]
        last = q_rows[-1]
        summary.extend([
            "",
            "## Key result",
            f"- Gamma: lambda = {float(gamma['wavelength_nm']):.6f} nm, Q = {float(gamma['q_value']):.6g}",
            f"- kx = {float(last['kx_norm']):.3f} pi/P: lambda = {float(last['wavelength_nm']):.6f} nm, Q = {float(last['q_value']):.6g}",
        ])
        (OUT / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    finally:
        client.remove(model.name())
        client.clear()


if __name__ == "__main__":
    run()
