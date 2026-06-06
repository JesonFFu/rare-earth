from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parent
SIN_DATA = WORKSPACE / "氮化硅光波导 --ZJ-LPCVD+Annealing-0718-SiN on SiO2-(1)((1).txt"
OUT = ROOT / "outputs"

COMSOLROOT = r"D:\comsol63\Multiphysics"
JAVA_HOME = rf"{COMSOLROOT}\java\win64\jre"
CONDA_PYTHON = r"D:\Anaconda\envs\comsol_env\python.exe"

LAMBDA0_NM = 980.0
C0 = 299_792_458.0


@dataclass
class ModeResult:
    period_nm: float
    fill: float
    ridge_width_nm: float
    taper_deg: float
    gap_nm: float
    bottom_gap_nm: float
    min_gap_nm: float
    remaining_air_gap_nm: float
    kx_norm: float
    mode_index: int
    freq_thz_real: float
    freq_thz_imag: float
    wavelength_nm: float
    q_value: float
    score: float


def configure_process_environment() -> None:
    os.environ["COMSOLROOT"] = COMSOLROOT
    os.environ["JAVA_HOME"] = JAVA_HOME
    prefix = rf"{COMSOLROOT}\bin\win64;{JAVA_HOME}\bin"
    os.environ["PATH"] = prefix + ";" + os.environ.get("PATH", "")
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")


def load_sin_index(wavelength_nm: float) -> tuple[float, float]:
    rows: list[tuple[float, float, float]] = []
    with SIN_DATA.open("r", encoding="utf-8") as handle:
        next(handle)
        next(handle)
        for line in handle:
            parts = line.split()
            if len(parts) >= 3:
                rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
    if not rows:
        raise RuntimeError(f"No SiN data loaded from {SIN_DATA}")
    rows.sort()
    if wavelength_nm < rows[0][0] or wavelength_nm > rows[-1][0]:
        raise ValueError(f"{wavelength_nm} nm is outside SiN table range.")
    for (x0, n0, k0), (x1, n1, k1) in zip(rows, rows[1:]):
        if x0 <= wavelength_nm <= x1:
            t = (wavelength_nm - x0) / (x1 - x0)
            return n0 + t * (n1 - n0), k0 + t * (k1 - k0)
    return rows[-1][1], rows[-1][2]


def import_mph():
    configure_process_environment()
    import mph

    return mph


def set_param(model, name: str, value: str, description: str = "") -> None:
    if description:
        model.java.param().set(name, value, description)
    else:
        model.java.param().set(name, value)


def create_rectangle(geom, tag: str, pos: list[str], size: list[str], label: str) -> None:
    rect = geom.create(tag, "Rectangle")
    rect.label(label)
    rect.set("base", "corner")
    rect.set("pos", pos)
    rect.set("size", size)
    rect.set("selresult", "on")
    rect.set("selresultshow", "all")


def create_polygon(geom, tag: str, x: list[str], y: list[str], label: str) -> None:
    poly = geom.create(tag, "Polygon")
    poly.label(label)
    poly.set("x", x)
    poly.set("y", y)
    poly.set("type", "solid")
    poly.set("selresult", "on")
    poly.set("selresultshow", "all")


def create_box_selection(comp, tag: str, entitydim: str, xmin: str, xmax: str, ymin: str, ymax: str) -> None:
    sel = comp.selection().create(tag, "Box")
    sel.label(tag)
    sel.set("entitydim", entitydim)
    sel.set("xmin", xmin)
    sel.set("xmax", xmax)
    sel.set("ymin", ymin)
    sel.set("ymax", ymax)
    sel.set("condition", "inside")


def create_union_selection(comp, tag: str, entitydim: str, inputs: list[str]) -> None:
    sel = comp.selection().create(tag, "Union")
    sel.label(tag)
    sel.set("entitydim", entitydim)
    sel.set("input", inputs)


def create_refractive_material(comp, tag: str, label: str, selection: str, n: str, k: str = "0") -> None:
    mat = comp.material().create(tag, "Common")
    mat.label(label)
    mat.selection().named(selection)
    mat.propertyGroup().create("RefractiveIndex", "Refractive index")
    mat.propertyGroup("RefractiveIndex").set("n", [n])
    mat.propertyGroup("RefractiveIndex").set("ki", [k])


def build_model(client, name: str = "sin980_qbic_2d"):
    n_sin, k_sin = load_sin_index(LAMBDA0_NM)
    model = client.create(name)
    jmodel = model.java

    set_param(model, "lambda0", f"{LAMBDA0_NM}[nm]", "Target wavelength")
    set_param(model, "P", "564.3[nm]", "Grating period")
    set_param(model, "W", "383.724[nm]", "Top/mask SiN ridge width")
    set_param(model, "fill", "W/P", "SiN duty cycle")
    set_param(model, "taper_deg", "0", "Signed sidewall angle in degrees; positive means bottom wider")
    set_param(model, "taper_dx", "0[nm]", "Signed lateral sidewall offset per side")
    set_param(model, "h_sin", "300[nm]", "Etched-through SiN grating height")
    set_param(model, "h_re", "16[nm]", "Conformal rare-earth equivalent layer thickness")
    set_param(model, "h_sio2", "3000[nm]", "SiO2 layer thickness")
    set_param(model, "h_si_buf", "500[nm]", "Thin Si buffer before bottom PML")
    set_param(model, "h_pml", "800[nm]", "PML thickness")
    set_param(model, "h_air", "1200[nm]", "Air spacer above grating")
    set_param(model, "kx", "0[1/m]", "Floquet wave vector along grating period")
    set_param(model, "n_sin", f"{n_sin:.9f}", "SiN refractive index at 980 nm")
    set_param(model, "k_sin", f"{k_sin:.9g}", "SiN extinction coefficient at 980 nm")
    set_param(model, "n_re", "1.49", "Rare-earth equivalent conformal layer refractive index")
    set_param(model, "k_re", "0", "Rare-earth equivalent conformal layer extinction coefficient")
    set_param(model, "n_sio2", "1.45", "SiO2 refractive index")
    set_param(model, "n_si", "3.55", "Lossless effective Si refractive index at 980 nm")

    jmodel.component().create("comp1", True)
    comp = jmodel.component("comp1")
    comp.label("2D x-y unit cell")
    comp.geom().create("geom1", 2)
    geom = comp.geom("geom1")
    geom.lengthUnit("nm")

    y_bot_pml = "0"
    y_si = "h_pml"
    y_sio2 = "h_pml+h_si_buf"
    y_sin = "h_pml+h_si_buf+h_sio2"
    y_air = "h_pml+h_si_buf+h_sio2+h_sin"
    y_air_main = "h_pml+h_si_buf+h_sio2+h_sin+h_re"
    y_top_pml = "h_pml+h_si_buf+h_sio2+h_sin+h_air"
    total_height = "2*h_pml+h_si_buf+h_sio2+h_sin+h_air"

    create_rectangle(geom, "si_pml", ["-P/2", y_bot_pml], ["P", "h_pml"], "Si bottom PML")
    create_rectangle(geom, "si_buf", ["-P/2", y_si], ["P", "h_si_buf"], "Si buffer")
    create_rectangle(geom, "sio2", ["-P/2", y_sio2], ["P", "h_sio2"], "SiO2")
    create_polygon(
        geom,
        "sin_ridge",
        ["-W/2-taper_dx", "W/2+taper_dx", "W/2", "-W/2"],
        [y_sin, y_sin, y_air, y_air],
        "300 nm tapered SiN ridge",
    )
    create_rectangle(geom, "re_top", ["-W/2", y_air], ["W", "h_re"], "Rare-earth conformal top layer")
    create_polygon(
        geom,
        "re_side_l",
        ["-W/2-taper_dx-h_re", "-W/2-taper_dx", "-W/2", "-W/2-h_re"],
        [y_sin, y_sin, y_air, y_air],
        "Rare-earth left sidewall layer",
    )
    create_polygon(
        geom,
        "re_side_r",
        ["W/2+taper_dx", "W/2+taper_dx+h_re", "W/2+h_re", "W/2"],
        [y_sin, y_sin, y_air, y_air],
        "Rare-earth right sidewall layer",
    )
    create_rectangle(geom, "re_bottom_l", ["-P/2", y_sin], ["P/2-W/2-taper_dx-h_re", "h_re"], "Rare-earth left groove bottom layer")
    create_rectangle(geom, "re_bottom_r", ["W/2+taper_dx+h_re", y_sin], ["P/2-W/2-taper_dx-h_re", "h_re"], "Rare-earth right groove bottom layer")
    create_polygon(
        geom,
        "air_slot_l",
        ["-P/2", "-W/2-taper_dx-h_re+taper_dx*h_re/h_sin", "-W/2-h_re", "-P/2"],
        ["h_pml+h_si_buf+h_sio2+h_re", "h_pml+h_si_buf+h_sio2+h_re", y_air, y_air],
        "Left remaining air groove",
    )
    create_polygon(
        geom,
        "air_slot_r",
        ["W/2+taper_dx+h_re-taper_dx*h_re/h_sin", "P/2", "P/2", "W/2+h_re"],
        ["h_pml+h_si_buf+h_sio2+h_re", "h_pml+h_si_buf+h_sio2+h_re", y_air, y_air],
        "Right remaining air groove",
    )
    create_rectangle(geom, "air_cap_l", ["-P/2", y_air], ["P/2-W/2", "h_re"], "Left air cap over groove")
    create_rectangle(geom, "air_cap_r", ["W/2", y_air], ["P/2-W/2", "h_re"], "Right air cap over groove")
    create_rectangle(geom, "air", ["-P/2", y_air_main], ["P", "h_air-h_re"], "Air")
    create_rectangle(geom, "air_pml", ["-P/2", y_top_pml], ["P", "h_pml"], "Air top PML")
    geom.run()

    create_union_selection(comp, "air_dom", "2", ["geom1_air_slot_l_dom", "geom1_air_slot_r_dom", "geom1_air_cap_l_dom", "geom1_air_cap_r_dom", "geom1_air_dom", "geom1_air_pml_dom"])
    create_union_selection(comp, "rare_earth_dom", "2", ["geom1_re_top_dom", "geom1_re_side_l_dom", "geom1_re_side_r_dom", "geom1_re_bottom_l_dom", "geom1_re_bottom_r_dom"])
    create_union_selection(comp, "si_dom", "2", ["geom1_si_pml_dom", "geom1_si_buf_dom"])
    create_box_selection(comp, "left_bnd", "1", "-P/2-1[nm]", "-P/2+1[nm]", "-1[nm]", total_height + "+1[nm]")
    create_box_selection(comp, "right_bnd", "1", "P/2-1[nm]", "P/2+1[nm]", "-1[nm]", total_height + "+1[nm]")
    create_union_selection(comp, "periodic_bnd", "1", ["left_bnd", "right_bnd"])

    pml_bot = comp.coordSystem().create("pml_bot", "PML")
    pml_bot.label("Bottom Si PML")
    pml_bot.selection().named("geom1_si_pml_dom")
    pml_top = comp.coordSystem().create("pml_top", "PML")
    pml_top.label("Top air PML")
    pml_top.selection().named("geom1_air_pml_dom")

    create_refractive_material(comp, "mat_air", "Air", "air_dom", "1", "0")
    create_refractive_material(comp, "mat_sio2", "SiO2", "geom1_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_sin", "SiN measured", "geom1_sin_ridge_dom", "n_sin", "k_sin")
    create_refractive_material(comp, "mat_re", "Rare-earth equivalent n=1.49", "rare_earth_dom", "n_re", "k_re")
    create_refractive_material(comp, "mat_si", "Si lossless effective", "si_dom", "n_si", "0")

    ewfd = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
    ewfd.label("Electromagnetic Waves, Frequency Domain")
    pc = ewfd.create("pc1", "PeriodicCondition", 1)
    pc.label("Floquet periodic x")
    pc.selection().named("periodic_bnd")
    pc.set("PeriodicType", "Floquet")
    pc.set("kFloquet", ["kx", "0", "0"])

    mesh = comp.mesh().create("mesh1")
    mesh.label("Physics-controlled mesh")
    mesh.autoMeshSize(3)

    study = jmodel.study().create("std1")
    study.label("Eigenfrequency near 980 nm")
    eig = study.create("eig", "Eigenfrequency")
    eig.set("shift", "c_const/lambda0")
    eig.set("neigs", "8")
    eig.set("neigsmanual", "8")
    eig.set("neigsactive", "on")
    eig.set("eigunit", "THz")

    return model


def taper_dx_nm(taper_deg: float, h_sin_nm: float = 300.0) -> float:
    return h_sin_nm * math.tan(math.radians(taper_deg))


def set_design(model, period_nm: float, fill: float, kx_norm: float = 0.0, taper_deg: float = 0.0, ridge_width_nm: float | None = None) -> None:
    width_nm = ridge_width_nm if ridge_width_nm is not None else period_nm * fill
    set_param(model, "P", f"{period_nm:.9g}[nm]")
    set_param(model, "W", f"{width_nm:.9g}[nm]")
    set_param(model, "fill", "W/P")
    set_param(model, "taper_deg", f"{taper_deg:.9g}")
    set_param(model, "taper_dx", f"{taper_dx_nm(taper_deg):.9g}[nm]")
    kx_expr = f"{kx_norm:.12g}*pi/P"
    set_param(model, "kx", kx_expr)


def solve_current_design(model) -> None:
    model.java.component("comp1").geom("geom1").run()
    model.java.component("comp1").mesh("mesh1").run()
    model.java.study("std1").run()


def unique_complex(values: Iterable[complex | float], tol: float = 1e-6) -> list[complex]:
    result: list[complex] = []
    for raw in values:
        val = complex(raw)
        if not any(abs(val - old) < tol for old in result):
            result.append(val)
    return result


def get_eigenfrequencies_thz(model) -> list[complex]:
    values = model.evaluate("freq", unit="THz")
    flat = values.ravel() if hasattr(values, "ravel") else values
    modes = unique_complex(flat, tol=1e-7)
    modes = [m for m in modes if math.isfinite(m.real) and abs(m.real) > 1e-9]
    modes.sort(key=lambda z: abs(C0 / (z.real * 1e12) * 1e9 - LAMBDA0_NM) if z.real else float("inf"))
    return modes


def mode_result(period_nm: float, ridge_width_nm: float, taper_deg: float, kx_norm: float, idx: int, freq_thz: complex) -> ModeResult:
    freq_hz_real = freq_thz.real * 1e12
    wavelength_nm = C0 / freq_hz_real * 1e9 if freq_hz_real > 0 else float("inf")
    q = abs(freq_thz.real / (2.0 * freq_thz.imag)) if abs(freq_thz.imag) > 0 else float("inf")
    wl_error = abs(wavelength_nm - LAMBDA0_NM)
    q_penalty = 0.0 if q >= 1e9 else math.log10(1e9 / max(q, 1.0))
    fill = ridge_width_nm / period_nm
    top_gap_nm = period_nm - ridge_width_nm
    bottom_width_nm = ridge_width_nm + 2.0 * taper_dx_nm(taper_deg)
    bottom_gap_nm = period_nm - bottom_width_nm
    min_gap_nm = min(top_gap_nm, bottom_gap_nm)
    remaining_air_gap_nm = min_gap_nm - 32.0
    score = wl_error + 10.0 * q_penalty
    return ModeResult(period_nm, fill, ridge_width_nm, taper_deg, top_gap_nm, bottom_gap_nm, min_gap_nm, remaining_air_gap_nm, kx_norm, idx, freq_thz.real, freq_thz.imag, wavelength_nm, q, score)


def solve_and_collect(model, period_nm: float, fill: float, kx_norm: float = 0.0, taper_deg: float = 0.0, ridge_width_nm: float | None = None) -> list[ModeResult]:
    width_nm = ridge_width_nm if ridge_width_nm is not None else period_nm * fill
    set_design(model, period_nm, fill, kx_norm, taper_deg=taper_deg, ridge_width_nm=width_nm)
    solve_current_design(model)
    modes = get_eigenfrequencies_thz(model)
    return [mode_result(period_nm, width_nm, taper_deg, kx_norm, i + 1, freq) for i, freq in enumerate(modes)]


def write_csv(path: Path, rows: list[ModeResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "period_nm",
            "fill",
            "ridge_width_nm",
            "taper_deg",
            "top_gap_nm",
            "bottom_gap_nm",
            "min_gap_nm",
            "remaining_air_gap_nm",
            "kx_norm",
            "mode_index",
            "freq_thz_real",
            "freq_thz_imag",
            "wavelength_nm",
            "q_value",
            "score",
        ])
        for r in rows:
            writer.writerow([
                r.period_nm,
                r.fill,
                r.ridge_width_nm,
                r.taper_deg,
                r.gap_nm,
                r.bottom_gap_nm,
                r.min_gap_nm,
                r.remaining_air_gap_nm,
                r.kx_norm,
                r.mode_index,
                r.freq_thz_real,
                r.freq_thz_imag,
                r.wavelength_nm,
                r.q_value,
                r.score,
            ])


def read_csv(path: Path) -> list[ModeResult]:
    rows: list[ModeResult] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            period = float(row["period_nm"])
            fill = float(row["fill"])
            width = float(row.get("ridge_width_nm", period * fill))
            taper = float(row.get("taper_deg", 0.0))
            top_gap = float(row.get("top_gap_nm", row.get("gap_nm", period - width)))
            bottom_gap = float(row.get("bottom_gap_nm", period - (width + 2.0 * taper_dx_nm(taper))))
            min_gap = float(row.get("min_gap_nm", min(top_gap, bottom_gap)))
            rows.append(
                ModeResult(
                    period,
                    fill,
                    width,
                    taper,
                    top_gap,
                    bottom_gap,
                    min_gap,
                    float(row.get("remaining_air_gap_nm", min_gap - 32.0)),
                    float(row["kx_norm"]),
                    int(row["mode_index"]),
                    float(row["freq_thz_real"]),
                    float(row["freq_thz_imag"]),
                    float(row["wavelength_nm"]),
                    float(row["q_value"]),
                    float(row["score"]),
                )
            )
    return rows


def best_rows(rows: list[ModeResult], limit: int = 10) -> list[ModeResult]:
    return sorted(rows, key=lambda r: (r.score, abs(r.wavelength_nm - LAMBDA0_NM), -r.q_value))[:limit]


def fabrication_feasible_rows(rows: list[ModeResult], min_gap_nm: float = 140.0) -> list[ModeResult]:
    return [r for r in rows if r.min_gap_nm >= min_gap_nm]


def qualified_rows(rows: list[ModeResult], min_gap_nm: float = 140.0, min_q: float = 1e9) -> list[ModeResult]:
    return [r for r in rows if r.min_gap_nm >= min_gap_nm and r.q_value >= min_q]


def near_gap_rows(rows: list[ModeResult], soft_min_gap_nm: float = 138.0, hard_min_gap_nm: float = 140.0) -> list[ModeResult]:
    return [r for r in rows if soft_min_gap_nm <= r.min_gap_nm < hard_min_gap_nm]


def svg_line_plot(path: Path, title: str, xs: list[float], series: list[tuple[str, list[float]]], y_label: str, log_y: bool = False, x_label: str = "kx / (pi/P)") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 900, 540
    ml, mr, mt, mb = 90, 30, 55, 80
    plot_w, plot_h = width - ml - mr, height - mt - mb
    all_y = [y for _, ys in series for y in ys if math.isfinite(y) and y > 0]
    if not xs or not all_y:
        path.write_text("<svg xmlns='http://www.w3.org/2000/svg'></svg>", encoding="utf-8")
        return
    x_min, x_max = min(xs), max(xs)
    if log_y:
        y_vals = [math.log10(max(y, 1e-300)) for y in all_y]
        y_min, y_max = min(y_vals), max(y_vals)
    else:
        y_min, y_max = min(all_y), max(all_y)
    if abs(x_max - x_min) < 1e-12:
        x_max += 1
    if abs(y_max - y_min) < 1e-12:
        y_max += 1

    def sx(x: float) -> float:
        return ml + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        yy = math.log10(max(y, 1e-300)) if log_y else y
        return mt + (y_max - yy) / (y_max - y_min) * plot_h

    colors = ["#0b5fff", "#d22f27", "#14853b", "#7c3aed", "#cc7a00"]
    out: list[str] = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2}' y='30' text-anchor='middle' font-size='22' font-family='Arial'>{title}</text>",
        f"<line x1='{ml}' y1='{mt+plot_h}' x2='{ml+plot_w}' y2='{mt+plot_h}' stroke='black'/>",
        f"<line x1='{ml}' y1='{mt}' x2='{ml}' y2='{mt+plot_h}' stroke='black'/>",
        f"<text x='{width/2}' y='{height-25}' text-anchor='middle' font-size='16' font-family='Arial'>{x_label}</text>",
        f"<text x='24' y='{height/2}' text-anchor='middle' transform='rotate(-90 24 {height/2})' font-size='16' font-family='Arial'>{y_label}</text>",
    ]
    for i in range(6):
        tx = x_min + i * (x_max - x_min) / 5
        px = sx(tx)
        out.append(f"<line x1='{px:.2f}' y1='{mt+plot_h}' x2='{px:.2f}' y2='{mt+plot_h+5}' stroke='black'/>")
        out.append(f"<text x='{px:.2f}' y='{mt+plot_h+24}' text-anchor='middle' font-size='12' font-family='Arial'>{tx:.2g}</text>")
    for i in range(6):
        ty_raw = y_min + i * (y_max - y_min) / 5
        ty = 10**ty_raw if log_y else ty_raw
        py = mt + (5 - i) / 5 * plot_h
        label = f"1e{ty_raw:.0f}" if log_y else f"{ty:.4g}"
        out.append(f"<line x1='{ml-5}' y1='{py:.2f}' x2='{ml}' y2='{py:.2f}' stroke='black'/>")
        out.append(f"<text x='{ml-10}' y='{py+4:.2f}' text-anchor='end' font-size='12' font-family='Arial'>{label}</text>")
    for si, (label, ys) in enumerate(series):
        pts = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(xs, ys) if math.isfinite(y) and y > 0)
        color = colors[si % len(colors)]
        out.append(f"<polyline points='{pts}' fill='none' stroke='{color}' stroke-width='2.5'/>")
        out.append(f"<text x='{ml+15}' y='{mt+20+si*22}' font-size='14' font-family='Arial' fill='{color}'>{label}</text>")
    out.append("</svg>")
    path.write_text("\n".join(out), encoding="utf-8")


def plasma_color(t: float) -> str:
    t = max(0.0, min(1.0, t))
    stops = [
        (0.0, (13, 8, 135)),
        (0.25, (84, 3, 160)),
        (0.50, (182, 54, 121)),
        (0.75, (251, 136, 97)),
        (1.0, (240, 249, 33)),
    ]
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t0 <= t <= t1:
            u = (t - t0) / (t1 - t0)
            rgb = tuple(round(c0[i] + u * (c1[i] - c0[i])) for i in range(3))
            return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
    return "#f0f921"


def svg_field_plot(path: Path, title: str, x_nm: list[float], y_nm: list[float], value: list[float], period_nm: float, fill: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 920, 620
    ml, mr, mt, mb = 80, 35, 55, 75
    plot_w, plot_h = width - ml - mr, height - mt - mb

    y_focus_min = 800 + 500 + 3000 - 250
    y_focus_max = 800 + 500 + 3000 + 300 + 650
    points = [(x, y, v) for x, y, v in zip(x_nm, y_nm, value) if -period_nm / 2 <= x <= period_nm / 2 and y_focus_min <= y <= y_focus_max]
    if not points:
        points = list(zip(x_nm, y_nm, value))
    if len(points) > 5500:
        step = max(1, len(points) // 5500)
        points = points[::step]

    x_min, x_max = -period_nm / 2, period_nm / 2
    y_min, y_max = min(p[1] for p in points), max(p[1] for p in points)
    vmax = max(max(p[2] for p in points), 1e-300)

    def sx(x: float) -> float:
        return ml + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        return mt + (y_max - y) / (y_max - y_min) * plot_h

    def normalized(v: float) -> float:
        return max(0.0, min(1.0, (math.log10(v / vmax + 1e-12) + 12.0) / 12.0))

    out = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2}' y='30' text-anchor='middle' font-size='22' font-family='Arial'>{title}</text>",
        f"<rect x='{ml}' y='{mt}' width='{plot_w}' height='{plot_h}' fill='#05051f' stroke='black'/>",
    ]
    for x, y, v in points:
        out.append(f"<circle cx='{sx(x):.2f}' cy='{sy(y):.2f}' r='1.8' fill='{plasma_color(normalized(v))}' opacity='0.9'/>")

    y_sio2_top = 800 + 500 + 3000
    y_sin_top = y_sio2_top + 300
    for yy, label in [(y_sio2_top, "SiO2 / SiN"), (y_sin_top, "SiN / air")]:
        if y_min <= yy <= y_max:
            py = sy(yy)
            out.append(f"<line x1='{ml}' y1='{py:.2f}' x2='{ml+plot_w}' y2='{py:.2f}' stroke='white' stroke-width='1.5' stroke-dasharray='6 4'/>")
            out.append(f"<text x='{ml+8}' y='{py-6:.2f}' font-size='12' font-family='Arial' fill='white'>{label}</text>")

    ridge_w = fill * period_nm
    rx0, rx1 = -ridge_w / 2, ridge_w / 2
    if y_min <= y_sio2_top <= y_max or y_min <= y_sin_top <= y_max:
        out.append(
            f"<rect x='{sx(rx0):.2f}' y='{sy(y_sin_top):.2f}' width='{sx(rx1)-sx(rx0):.2f}' height='{sy(y_sio2_top)-sy(y_sin_top):.2f}' "
            "fill='none' stroke='white' stroke-width='2'/>"
        )
        h_re = 16.0
        re_rects = [
            (rx0, y_sin_top, rx1, y_sin_top + h_re),
            (rx0 - h_re, y_sio2_top, rx0, y_sin_top),
            (rx1, y_sio2_top, rx1 + h_re, y_sin_top),
            (x_min, y_sio2_top, rx0 - h_re, y_sio2_top + h_re),
            (rx1 + h_re, y_sio2_top, x_max, y_sio2_top + h_re),
        ]
        for x0, y0, x1, y1 in re_rects:
            out.append(
                f"<rect x='{sx(x0):.2f}' y='{sy(y1):.2f}' width='{sx(x1)-sx(x0):.2f}' height='{sy(y0)-sy(y1):.2f}' "
                "fill='none' stroke='#ff4fd8' stroke-width='1.6'/>"
            )
        out.append(f"<text x='{ml+8}' y='{mt+20}' font-size='12' font-family='Arial' fill='white'>white: SiN</text>")
        out.append(f"<text x='{ml+8}' y='{mt+38}' font-size='12' font-family='Arial' fill='#ff4fd8'>magenta: 16 nm RE layer</text>")

    out.extend(
        [
            f"<line x1='{ml}' y1='{mt+plot_h}' x2='{ml+plot_w}' y2='{mt+plot_h}' stroke='black'/>",
            f"<line x1='{ml}' y1='{mt}' x2='{ml}' y2='{mt+plot_h}' stroke='black'/>",
            f"<text x='{width/2}' y='{height-25}' text-anchor='middle' font-size='16' font-family='Arial'>x (nm)</text>",
            f"<text x='24' y='{height/2}' text-anchor='middle' transform='rotate(-90 24 {height/2})' font-size='16' font-family='Arial'>vertical coordinate y (nm)</text>",
            "</svg>",
        ]
    )
    path.write_text("\n".join(out), encoding="utf-8")


def find_solution_index_for_row(freqs: list[complex], target: ModeResult) -> int:
    def metric(item: tuple[int, complex]) -> float:
        _, freq = item
        wl = C0 / (freq.real * 1e12) * 1e9
        q = abs(freq.real / (2.0 * freq.imag)) if abs(freq.imag) > 0 else float("inf")
        return abs(wl - target.wavelength_nm) + 1e-6 * abs(math.log10(max(q, 1.0)) - math.log10(max(target.q_value, 1.0)))

    return min(enumerate(freqs, start=1), key=metric)[0]


def run_coarse(args) -> None:
    mph = import_mph()
    OUT.mkdir(parents=True, exist_ok=True)
    client = mph.Client(cores=args.cores)
    model = build_model(client, "sin980_qbic_2d")
    rows: list[ModeResult] = []
    periods = [float(p) for p in args.periods.split(",")]
    fills = [float(f) for f in args.fills.split(",")]
    try:
        for period in periods:
            for fill in fills:
                print(f"Solving Gamma design P={period:.3f} nm, fill={fill:.3f}", flush=True)
                try:
                    current = solve_and_collect(model, period, fill, 0.0)
                    rows.extend(current)
                    write_csv(OUT / "coarse_scan.csv", rows)
                    best = best_rows(current, 1)[0]
                    print(f"  best lambda={best.wavelength_nm:.6g} nm, Q={best.q_value:.4g}, mode={best.mode_index}", flush=True)
                except Exception as exc:
                    print(f"  FAILED: {exc}", flush=True)
        best_unconstrained = best_rows(rows, 1)[0]
        feasible = fabrication_feasible_rows(rows, args.min_gap)
        qualified = qualified_rows(rows, args.min_gap, 1e9)
        best_gap = best_rows(qualified or feasible or rows, 1)[0]

        set_design(model, best_unconstrained.period_nm, best_unconstrained.fill, 0.0)
        solve_current_design(model)
        model.save(str(OUT / "sin980_qbic_2d_v2_best_unconstrained.mph"))

        set_design(model, best_gap.period_nm, best_gap.fill, 0.0)
        solve_current_design(model)
        model.save(str(OUT / "sin980_qbic_2d_v2_best_gap_feasible.mph"))

        write_csv(OUT / "best_unconstrained.csv", [best_unconstrained])
        write_csv(OUT / "best_gap_feasible.csv", [best_gap])
        write_csv(OUT / "best_gamma.csv", [best_gap])
        write_csv(OUT / "gap_check.csv", best_rows(rows, min(30, len(rows))))
        print(
            f"Best unconstrained: P={best_unconstrained.period_nm:.6g} nm fill={best_unconstrained.fill:.6g} "
            f"gap={best_unconstrained.gap_nm:.6g} nm lambda={best_unconstrained.wavelength_nm:.9g} nm Q={best_unconstrained.q_value:.6g}",
            flush=True,
        )
        print(
            f"Best gap-feasible: P={best_gap.period_nm:.6g} nm fill={best_gap.fill:.6g} "
            f"gap={best_gap.gap_nm:.6g} nm lambda={best_gap.wavelength_nm:.9g} nm Q={best_gap.q_value:.6g}",
            flush=True,
        )
    finally:
        client.clear()


def run_band(args) -> None:
    mph = import_mph()
    OUT.mkdir(parents=True, exist_ok=True)
    client = mph.Client(cores=args.cores)
    model = build_model(client, "sin980_qbic_2d_v2_band")
    rows: list[ModeResult] = []
    k_values = [float(k) for k in args.kpoints.split(",")]
    try:
        for kx_norm in k_values:
            print(f"Solving band point kx={kx_norm:.5g} pi/P", flush=True)
            current = solve_and_collect(model, args.period, args.fill, kx_norm)
            rows.extend(current)
            write_csv(OUT / "band_scan.csv", rows)
        best_per_k: list[ModeResult] = []
        for kx_norm in k_values:
            subset = [r for r in rows if abs(r.kx_norm - kx_norm) < 1e-12]
            if subset:
                best_per_k.append(best_rows(subset, 1)[0])
        write_csv(OUT / "band_best_branch.csv", best_per_k)
        svg_line_plot(
            OUT / "band_wavelength.svg",
            "SiN 1D quasi-BIC band near 980 nm",
            [r.kx_norm for r in best_per_k],
            [("best mode", [r.wavelength_nm for r in best_per_k])],
            "wavelength (nm)",
        )
        svg_line_plot(
            OUT / "band_q.svg",
            "Q factor along Gamma-X",
            [r.kx_norm for r in best_per_k],
            [("best mode", [r.q_value for r in best_per_k])],
            "Q",
            log_y=True,
        )
        best = best_rows([r for r in rows if r.kx_norm == 0.0], 1)[0]
        set_design(model, best.period_nm, best.fill, 0.0)
        solve_current_design(model)
        model.save(str(OUT / "sin980_qbic_2d_v2_band.mph"))
    finally:
        client.clear()


def parse_float_list(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def run_robust(args) -> None:
    mph = import_mph()
    OUT.mkdir(parents=True, exist_ok=True)
    client = mph.Client(cores=args.cores)
    model = build_model(client, "sin980_qbic_2d_v3_robustness")
    rows: list[ModeResult] = []
    dps = parse_float_list(args.dperiods)
    dws = parse_float_list(args.dwidths)
    tapers = parse_float_list(args.tapers)
    try:
        for taper in tapers:
            for dp in dps:
                for dw in dws:
                    period = args.period + dp
                    width = args.width + dw
                    fill = width / period
                    print(
                        f"Solving robustness P={period:.3f} nm (dP={dp:+.3f}), W={width:.3f} nm (dW={dw:+.3f}), taper={taper:+.3f} deg",
                        flush=True,
                    )
                    try:
                        current = solve_and_collect(model, period, fill, 0.0, taper_deg=taper, ridge_width_nm=width)
                        best = best_rows(current, 1)[0]
                        rows.append(best)
                        write_csv(OUT / "robustness_scan.csv", rows)
                        print(
                            f"  lambda={best.wavelength_nm:.6g} nm, Q={best.q_value:.4g}, min_gap={best.min_gap_nm:.3f} nm",
                            flush=True,
                        )
                    except Exception as exc:
                        print(f"  FAILED: {exc}", flush=True)

        center = min(rows, key=lambda r: abs(r.period_nm - args.period) + abs(r.ridge_width_nm - args.width) + abs(r.taper_deg))
        set_design(model, center.period_nm, center.fill, 0.0, taper_deg=center.taper_deg, ridge_width_nm=center.ridge_width_nm)
        solve_current_design(model)
        model.save(str(OUT / "sin980_qbic_2d_v3_nominal_large_gap_re.mph"))

        worst_q = min(rows, key=lambda r: r.q_value)
        set_design(model, worst_q.period_nm, worst_q.fill, 0.0, taper_deg=worst_q.taper_deg, ridge_width_nm=worst_q.ridge_width_nm)
        solve_current_design(model)
        model.save(str(OUT / "sin980_qbic_2d_v3_worst_q_case.mph"))

        write_robustness_report(args, rows)
    finally:
        client.clear()


def write_robustness_report(args, rows: list[ModeResult]) -> None:
    if not rows:
        return
    center = min(rows, key=lambda r: abs(r.period_nm - args.period) + abs(r.ridge_width_nm - args.width) + abs(r.taper_deg))
    worst_q = min(rows, key=lambda r: r.q_value)
    worst_detune = max(rows, key=lambda r: abs(r.wavelength_nm - LAMBDA0_NM))
    passing_q = [r for r in rows if r.q_value >= 1e9]
    passing_gap = [r for r in rows if r.min_gap_nm >= 140.0]

    def fmt_row(label: str, r: ModeResult) -> str:
        return (
            f"| {label} | {r.period_nm:.6g} | {r.ridge_width_nm:.6g} | {r.fill:.6g} | {r.taper_deg:.6g} | "
            f"{r.gap_nm:.6g} | {r.bottom_gap_nm:.6g} | {r.min_gap_nm:.6g} | {r.remaining_air_gap_nm:.6g} | "
            f"{r.wavelength_nm:.9g} | {r.q_value:.6g} |"
        )

    top = best_rows(rows, min(12, len(rows)))
    low_q = sorted(rows, key=lambda r: r.q_value)[: min(8, len(rows))]
    lines = [
        "# V3 robustness scan: large-gap SiN grating with 16 nm RE layer",
        "",
        "## Scan setup",
        "",
        f"- Nominal center: `P={args.period} nm`, `W={args.width} nm`, top nominal gap `{args.period - args.width:.6g} nm`.",
        f"- Period perturbations: `{args.dperiods}` nm.",
        f"- Ridge-width perturbations: `{args.dwidths}` nm.",
        f"- Sidewall taper angles: `{args.tapers}` deg; positive means bottom is wider than top.",
        "- Each result tracks the nearest 980 nm Gamma-point eigenmode.",
        "",
        "## Summary",
        "",
        f"- Simulated cases: `{len(rows)}`.",
        f"- Cases with `Q >= 1e9`: `{len(passing_q)}/{len(rows)}`.",
        f"- Cases with minimum etched gap `>=140 nm`: `{len(passing_gap)}/{len(rows)}`.",
        f"- Worst Q in scan: `{worst_q.q_value:.6g}` at `P={worst_q.period_nm:.6g} nm`, `W={worst_q.ridge_width_nm:.6g} nm`, `taper={worst_q.taper_deg:.6g} deg`.",
        f"- Largest wavelength detuning: `{worst_detune.wavelength_nm - LAMBDA0_NM:.6g} nm`.",
        "",
        "## Key cases",
        "",
        "| case | P (nm) | W_top (nm) | fill | taper (deg) | top gap (nm) | bottom gap (nm) | min gap (nm) | remaining air gap (nm) | lambda (nm) | Q |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        fmt_row("nominal", center),
        fmt_row("worst Q", worst_q),
        fmt_row("largest detuning", worst_detune),
        "",
        "## Best candidates in robustness grid",
        "",
        "| rank | P (nm) | W_top (nm) | fill | taper (deg) | top gap (nm) | bottom gap (nm) | min gap (nm) | remaining air gap (nm) | lambda (nm) | Q |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for i, r in enumerate(top, 1):
        lines.append(
            f"| {i} | {r.period_nm:.6g} | {r.ridge_width_nm:.6g} | {r.fill:.6g} | {r.taper_deg:.6g} | "
            f"{r.gap_nm:.6g} | {r.bottom_gap_nm:.6g} | {r.min_gap_nm:.6g} | {r.remaining_air_gap_nm:.6g} | "
            f"{r.wavelength_nm:.9g} | {r.q_value:.6g} |"
        )
    lines.extend(
        [
            "",
            "## Lowest-Q cases",
            "",
            "| P (nm) | W_top (nm) | fill | taper (deg) | min gap (nm) | lambda (nm) | Q |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for r in low_q:
        lines.append(f"| {r.period_nm:.6g} | {r.ridge_width_nm:.6g} | {r.fill:.6g} | {r.taper_deg:.6g} | {r.min_gap_nm:.6g} | {r.wavelength_nm:.9g} | {r.q_value:.6g} |")

    if worst_q.q_value >= 1e9:
        lines.extend(
            [
                "",
                "## Recommendation",
                "",
                "Within this perturbation grid, the large-gap RE design remains robust: every simulated case keeps `Q >= 1e9`.",
                "The strongest effect is wavelength detuning, so experimentally the safer strategy is to fabricate a period/width matrix around the nominal design and locate the 980 nm resonance optically.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "## Recommendation",
                "",
                "At least one perturbation case falls below `Q = 1e9`; use the table above to bias the fabricated design away from the sensitive region.",
            ]
        )

    (OUT / "robustness_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(OUT / "robustness_summary.md")


def run_report(args) -> None:
    coarse = OUT / "coarse_scan.csv"
    band = OUT / "band_best_branch.csv"
    v1_best = WORKSPACE / "sin980_qbic_2d_design" / "outputs" / "best_gamma_final.csv"
    lines: list[str] = ["# SiN 980 nm 2D quasi-BIC v2 rare-earth perturbation summary", ""]
    if v1_best.exists() and (OUT / "best_gap_feasible.csv").exists():
        old = read_csv(v1_best)[0]
        new = read_csv(OUT / "best_gap_feasible.csv")[0]
        lines.append("## V1 to V2 comparison")
        lines.append("")
        lines.append("| model | P (nm) | fill | nominal gap (nm) | remaining air gap (nm) | lambda (nm) | Q |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        lines.append(f"| V1 no RE | {old.period_nm:.6g} | {old.fill:.6g} | {old.gap_nm:.6g} | n/a | {old.wavelength_nm:.9g} | {old.q_value:.6g} |")
        lines.append(f"| V2 16 nm RE | {new.period_nm:.6g} | {new.fill:.6g} | {new.gap_nm:.6g} | {new.remaining_air_gap_nm:.6g} | {new.wavelength_nm:.9g} | {new.q_value:.6g} |")
        lines.append("")
        lines.append(f"V2 wavelength shift from V1: {new.wavelength_nm - old.wavelength_nm:.6g} nm.")
        lines.append("")
    if coarse.exists():
        rows = read_csv(coarse)
        best = best_rows(rows, 10)
        grouped: dict[tuple[float, float], list[ModeResult]] = {}
        for row in rows:
            grouped.setdefault((row.period_nm, row.fill), []).append(row)
        scan_best = [best_rows(group, 1)[0] for group in grouped.values()]
        fills = sorted(set(r.fill for r in scan_best))
        if fills:
            xs = sorted(set(r.period_nm for r in scan_best))
            wl_series = []
            q_series = []
            for fill in fills:
                branch = [next((r for r in scan_best if abs(r.period_nm - x) < 1e-12 and abs(r.fill - fill) < 1e-12), None) for x in xs]
                wl_series.append((f"fill={fill:.3g}", [r.wavelength_nm if r else float("nan") for r in branch]))
                q_series.append((f"fill={fill:.3g}", [r.q_value if r else float("nan") for r in branch]))
            svg_line_plot(OUT / "gamma_scan_wavelength.svg", "Gamma-point resonance tuning", xs, wl_series, "wavelength (nm)", x_label="period (nm)")
            svg_line_plot(OUT / "gamma_scan_q.svg", "Gamma-point Q tuning", xs, q_series, "Q", log_y=True, x_label="period (nm)")
        lines.append("## Best Gamma candidates")
        lines.append("")
        lines.append("Generated Gamma plots: `gamma_scan_wavelength.svg`, `gamma_scan_q.svg`.")
        lines.append("")
        lines.append("| rank | P (nm) | fill | gap (nm) | remaining air gap (nm) | lambda (nm) | Q | mode |")
        lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
        for i, r in enumerate(best, 1):
            lines.append(f"| {i} | {r.period_nm:.6g} | {r.fill:.6g} | {r.gap_nm:.6g} | {r.remaining_air_gap_nm:.6g} | {r.wavelength_nm:.9g} | {r.q_value:.6g} | {r.mode_index} |")
        lines.append("")
        if (OUT / "best_unconstrained.csv").exists() and (OUT / "best_gap_feasible.csv").exists():
            unconstrained = read_csv(OUT / "best_unconstrained.csv")[0]
            feasible = read_csv(OUT / "best_gap_feasible.csv")[0]
            near = best_rows(near_gap_rows(rows), min(5, len(near_gap_rows(rows))))
            lines.append("## Gap post-check")
            lines.append("")
            lines.append("| candidate | P (nm) | fill | gap (nm) | remaining air gap (nm) | lambda (nm) | Q |")
            lines.append("|---|---:|---:|---:|---:|---:|---:|")
            lines.append(f"| unconstrained best | {unconstrained.period_nm:.6g} | {unconstrained.fill:.6g} | {unconstrained.gap_nm:.6g} | {unconstrained.remaining_air_gap_nm:.6g} | {unconstrained.wavelength_nm:.9g} | {unconstrained.q_value:.6g} |")
            lines.append(f"| recommended compact-gap | {feasible.period_nm:.6g} | {feasible.fill:.6g} | {feasible.gap_nm:.6g} | {feasible.remaining_air_gap_nm:.6g} | {feasible.wavelength_nm:.9g} | {feasible.q_value:.6g} |")
            lines.append("")
            if near:
                lines.append("Near-threshold candidates with `138 nm <= gap < 140 nm` are treated as acceptable backups if they clearly outperform hard-gap designs.")
                lines.append("")
                lines.append("| P (nm) | fill | gap (nm) | remaining air gap (nm) | lambda (nm) | Q |")
                lines.append("|---:|---:|---:|---:|---:|---:|")
                for r in near:
                    lines.append(f"| {r.period_nm:.6g} | {r.fill:.6g} | {r.gap_nm:.6g} | {r.remaining_air_gap_nm:.6g} | {r.wavelength_nm:.9g} | {r.q_value:.6g} |")
                lines.append("")
    if band.exists():
        rows = read_csv(band)
        lines.append("## Gamma-X branch")
        lines.append("")
        lines.append("| kx/(pi/P) | lambda (nm) | Q |")
        lines.append("|---:|---:|---:|")
        for r in rows:
            lines.append(f"| {r.kx_norm:.6g} | {r.wavelength_nm:.9g} | {r.q_value:.6g} |")
        lines.append("")
        lines.append("Generated plots: `band_wavelength.svg`, `band_q.svg`.")
        lines.append("")
    if (OUT / "field_ez2.svg").exists():
        lines.append("## Field plot")
        lines.append("")
        lines.append("Generated field output: `field_ez2.csv`, `field_ez2.svg`.")
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(OUT / "summary.md")


def run_field(args) -> None:
    mph = import_mph()
    import numpy as np

    model_path = Path(args.model)
    best_path = Path(args.best)
    rows = read_csv(best_path)
    target = best_rows(rows, 1)[0]
    client = mph.Client(cores=args.cores)
    try:
        model = client.load(str(model_path))
        freqs = [complex(v) for v in np.asarray(model.evaluate("freq", unit="THz")).ravel()]
        sol_idx = find_solution_index_for_row(freqs, target)
        ds = model.datasets()[0] if model.datasets() else None
        x = np.asarray(model.evaluate("x", unit="nm", dataset=ds, inner=[sol_idx])).ravel().astype(float)
        y = np.asarray(model.evaluate("y", unit="nm", dataset=ds, inner=[sol_idx])).ravel().astype(float)
        ez = np.asarray(model.evaluate("abs(ewfd.Ez)", dataset=ds, inner=[sol_idx])).ravel().astype(float)
        n = min(len(x), len(y), len(ez))
        x, y, ez = x[:n], y[:n], ez[:n]
        ez2 = ez * ez
        csv_path = OUT / "field_ez2.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["x_nm", "y_nm", "abs_ez", "abs_ez_squared"])
            for row in zip(x, y, ez, ez2):
                writer.writerow(row)
        svg_field_plot(
            OUT / "field_ez2.svg",
            f"|Ez|^2 field with 16 nm RE layer, P={target.period_nm:.3f} nm, fill={target.fill:.3f}",
            x.tolist(),
            y.tolist(),
            ez2.tolist(),
            target.period_nm,
            target.fill,
        )
        print(f"Matched solution index {sol_idx}; wrote {csv_path} and {OUT / 'field_ez2.svg'}")
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build and solve a 2D SiN 980 nm quasi-BIC grating model in COMSOL.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_coarse = sub.add_parser("coarse", help="Run Gamma-point period/fill scan.")
    p_coarse.add_argument("--cores", type=int, default=2)
    p_coarse.add_argument("--periods", default="540,546,552,558,564,570")
    p_coarse.add_argument("--fills", default="0.68,0.70,0.72,0.74,0.76,0.78")
    p_coarse.add_argument("--min-gap", type=float, default=140.0)
    p_coarse.set_defaults(func=run_coarse)

    p_band = sub.add_parser("band", help="Run Gamma-X band scan for one design.")
    p_band.add_argument("--cores", type=int, default=2)
    p_band.add_argument("--period", type=float, required=True)
    p_band.add_argument("--fill", type=float, required=True)
    p_band.add_argument("--kpoints", default="0,0.02,0.05,0.1,0.2,0.35,0.5,0.75,1")
    p_band.set_defaults(func=run_band)

    p_robust = sub.add_parser("robust", help="Run fabrication robustness scan around the large-gap RE design.")
    p_robust.add_argument("--cores", type=int, default=2)
    p_robust.add_argument("--period", type=float, default=564.3)
    p_robust.add_argument("--width", type=float, default=383.724)
    p_robust.add_argument("--dperiods", default="-5,-2,0,2,5")
    p_robust.add_argument("--dwidths", default="-10,-5,0,5,10")
    p_robust.add_argument("--tapers", default="0,2,-2")
    p_robust.set_defaults(func=run_robust)

    p_report = sub.add_parser("report", help="Generate Markdown summary from CSV outputs.")
    p_report.set_defaults(func=run_report)

    p_field = sub.add_parser("field", help="Export |Ez|^2 field CSV and SVG for the best Gamma model.")
    p_field.add_argument("--cores", type=int, default=2)
    p_field.add_argument("--model", default=str(OUT / "sin980_qbic_2d_v2_best_gap_feasible.mph"))
    p_field.add_argument("--best", default=str(OUT / "best_gap_feasible.csv"))
    p_field.set_defaults(func=run_field)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
