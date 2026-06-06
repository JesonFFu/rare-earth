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
    set_param(model, "P", "560[nm]", "Grating period")
    set_param(model, "fill", "0.50", "SiN duty cycle")
    set_param(model, "h_sin", "300[nm]", "Etched-through SiN grating height")
    set_param(model, "h_sio2", "3000[nm]", "SiO2 layer thickness")
    set_param(model, "h_si_buf", "500[nm]", "Thin Si buffer before bottom PML")
    set_param(model, "h_pml", "800[nm]", "PML thickness")
    set_param(model, "h_air", "1200[nm]", "Air spacer above grating")
    set_param(model, "kx", "0[1/m]", "Floquet wave vector along grating period")
    set_param(model, "n_sin", f"{n_sin:.9f}", "SiN refractive index at 980 nm")
    set_param(model, "k_sin", f"{k_sin:.9g}", "SiN extinction coefficient at 980 nm")
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
    y_top_pml = "h_pml+h_si_buf+h_sio2+h_sin+h_air"
    total_height = "2*h_pml+h_si_buf+h_sio2+h_sin+h_air"

    create_rectangle(geom, "si_pml", ["-P/2", y_bot_pml], ["P", "h_pml"], "Si bottom PML")
    create_rectangle(geom, "si_buf", ["-P/2", y_si], ["P", "h_si_buf"], "Si buffer")
    create_rectangle(geom, "sio2", ["-P/2", y_sio2], ["P", "h_sio2"], "SiO2")
    create_rectangle(geom, "sin_ridge", ["-fill*P/2", y_sin], ["fill*P", "h_sin"], "300 nm SiN ridge")
    create_rectangle(geom, "air_slot_l", ["-P/2", y_sin], ["(1-fill)*P/2", "h_sin"], "Left air groove")
    create_rectangle(geom, "air_slot_r", ["fill*P/2", y_sin], ["(1-fill)*P/2", "h_sin"], "Right air groove")
    create_rectangle(geom, "air", ["-P/2", y_air], ["P", "h_air"], "Air")
    create_rectangle(geom, "air_pml", ["-P/2", y_top_pml], ["P", "h_pml"], "Air top PML")
    geom.run()

    create_union_selection(comp, "air_dom", "2", ["geom1_air_slot_l_dom", "geom1_air_slot_r_dom", "geom1_air_dom", "geom1_air_pml_dom"])
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


def set_design(model, period_nm: float, fill: float, kx_norm: float = 0.0) -> None:
    set_param(model, "P", f"{period_nm:.9g}[nm]")
    set_param(model, "fill", f"{fill:.9g}")
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


def mode_result(period_nm: float, fill: float, kx_norm: float, idx: int, freq_thz: complex) -> ModeResult:
    freq_hz_real = freq_thz.real * 1e12
    wavelength_nm = C0 / freq_hz_real * 1e9 if freq_hz_real > 0 else float("inf")
    q = abs(freq_thz.real / (2.0 * freq_thz.imag)) if abs(freq_thz.imag) > 0 else float("inf")
    wl_error = abs(wavelength_nm - LAMBDA0_NM)
    q_penalty = 0.0 if q >= 1e9 else math.log10(1e9 / max(q, 1.0))
    score = wl_error + 10.0 * q_penalty
    return ModeResult(period_nm, fill, kx_norm, idx, freq_thz.real, freq_thz.imag, wavelength_nm, q, score)


def solve_and_collect(model, period_nm: float, fill: float, kx_norm: float = 0.0) -> list[ModeResult]:
    set_design(model, period_nm, fill, kx_norm)
    solve_current_design(model)
    modes = get_eigenfrequencies_thz(model)
    return [mode_result(period_nm, fill, kx_norm, i + 1, freq) for i, freq in enumerate(modes)]


def write_csv(path: Path, rows: list[ModeResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["period_nm", "fill", "kx_norm", "mode_index", "freq_thz_real", "freq_thz_imag", "wavelength_nm", "q_value", "score"])
        for r in rows:
            writer.writerow([r.period_nm, r.fill, r.kx_norm, r.mode_index, r.freq_thz_real, r.freq_thz_imag, r.wavelength_nm, r.q_value, r.score])


def read_csv(path: Path) -> list[ModeResult]:
    rows: list[ModeResult] = []
    with path.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                ModeResult(
                    float(row["period_nm"]),
                    float(row["fill"]),
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


def svg_line_plot(path: Path, title: str, xs: list[float], series: list[tuple[str, list[float]]], y_label: str, log_y: bool = False) -> None:
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
        f"<text x='{width/2}' y='{height-25}' text-anchor='middle' font-size='16' font-family='Arial'>kx / (pi/P)</text>",
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
        best = best_rows(rows, 1)[0]
        set_design(model, best.period_nm, best.fill, 0.0)
        solve_current_design(model)
        model.save(str(OUT / "sin980_qbic_2d_best_gamma.mph"))
        write_csv(OUT / "best_gamma.csv", [best])
        print(f"Best Gamma candidate: P={best.period_nm:.6g} nm fill={best.fill:.6g} lambda={best.wavelength_nm:.9g} nm Q={best.q_value:.6g}", flush=True)
    finally:
        client.clear()


def run_band(args) -> None:
    mph = import_mph()
    OUT.mkdir(parents=True, exist_ok=True)
    client = mph.Client(cores=args.cores)
    model = build_model(client, "sin980_qbic_2d_band")
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
        model.save(str(OUT / "sin980_qbic_2d_band.mph"))
    finally:
        client.clear()


def run_report(args) -> None:
    coarse = OUT / "coarse_scan.csv"
    band = OUT / "band_best_branch.csv"
    lines: list[str] = ["# SiN 980 nm 2D quasi-BIC simulation summary", ""]
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
            svg_line_plot(OUT / "gamma_scan_wavelength.svg", "Gamma-point resonance tuning", xs, wl_series, "wavelength (nm)")
            svg_line_plot(OUT / "gamma_scan_q.svg", "Gamma-point Q tuning", xs, q_series, "Q", log_y=True)
        lines.append("## Best Gamma candidates")
        lines.append("")
        lines.append("Generated Gamma plots: `gamma_scan_wavelength.svg`, `gamma_scan_q.svg`.")
        lines.append("")
        lines.append("| rank | P (nm) | fill | lambda (nm) | Q | mode |")
        lines.append("|---:|---:|---:|---:|---:|---:|")
        for i, r in enumerate(best, 1):
            lines.append(f"| {i} | {r.period_nm:.6g} | {r.fill:.6g} | {r.wavelength_nm:.9g} | {r.q_value:.6g} | {r.mode_index} |")
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
            f"|Ez|^2 field, P={target.period_nm:.3f} nm, fill={target.fill:.3f}",
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
    p_coarse.add_argument("--periods", default="480,520,560,600,640,680,700")
    p_coarse.add_argument("--fills", default="0.25,0.35,0.45,0.55,0.65,0.75")
    p_coarse.set_defaults(func=run_coarse)

    p_band = sub.add_parser("band", help="Run Gamma-X band scan for one design.")
    p_band.add_argument("--cores", type=int, default=2)
    p_band.add_argument("--period", type=float, required=True)
    p_band.add_argument("--fill", type=float, required=True)
    p_band.add_argument("--kpoints", default="0,0.02,0.05,0.1,0.2,0.35,0.5,0.75,1")
    p_band.set_defaults(func=run_band)

    p_report = sub.add_parser("report", help="Generate Markdown summary from CSV outputs.")
    p_report.set_defaults(func=run_report)

    p_field = sub.add_parser("field", help="Export |Ez|^2 field CSV and SVG for the best Gamma model.")
    p_field.add_argument("--cores", type=int, default=2)
    p_field.add_argument("--model", default=str(OUT / "sin980_qbic_2d_best_gamma.mph"))
    p_field.add_argument("--best", default=str(OUT / "best_gamma_final.csv"))
    p_field.set_defaults(func=run_field)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
