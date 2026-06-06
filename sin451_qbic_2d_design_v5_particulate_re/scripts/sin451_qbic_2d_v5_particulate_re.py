from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
OUT = ROOT / "outputs"

COMSOLROOT = r"D:\comsol63\Multiphysics"
JAVA_HOME = rf"{COMSOLROOT}\java\win64\jre"
C0 = 299_792_458.0

TARGET_NM = 451.0
N_SIN = 2.128247126
N_SIO2 = 1.466
N_SI = 4.7
N_RE = 1.52

H_SIN_NM = 300.0
H_RE_NM = 16.0
H_SIO2_NM = 3000.0
H_SI_BUF_NM = 500.0
H_PML_NM = 800.0
H_AIR_NM = 1200.0

PARTICLE_W_NM = 25.0
PARTICLE_GAP_NM = 8.0
SIDE_PARTICLE_H_NM = 25.0
SIDE_PARTICLE_GAP_NM = 25.0

MIN_GAP_NM = 130.0
PREFERRED_GAP_NM = 140.0
MIN_RIDGE_W_NM = 100.0
TARGET_Q = 1e9
ACCEPT_Q = 1e8


@dataclass(frozen=True)
class Design:
    design_id: str
    period_nm: float
    width_nm: float

    @property
    def gap_nm(self) -> float:
        return self.period_nm - self.width_nm


@dataclass
class ModeResult:
    design_id: str
    re_style: str
    period_nm: float
    width_nm: float
    gap_nm: float
    kx_norm: float
    mode_index: int
    freq_thz_real: float
    freq_thz_imag: float
    wavelength_nm: float
    q_value: float
    particle_count: int
    top_coverage: float
    bottom_coverage: float
    sidewall_coverage: float
    re_overlap: float
    score: float


def configure_process_environment() -> None:
    os.environ["COMSOLROOT"] = COMSOLROOT
    os.environ["JAVA_HOME"] = JAVA_HOME
    prefix = rf"{COMSOLROOT}\bin\win64;{JAVA_HOME}\bin"
    os.environ["PATH"] = prefix + ";" + os.environ.get("PATH", "")
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")


def import_mph():
    configure_process_environment()
    import mph

    return mph


def set_param(model, name: str, value: str, description: str = "") -> None:
    if description:
        model.java.param().set(name, value, description)
    else:
        model.java.param().set(name, value)


def fmt(value: float | str) -> str:
    if isinstance(value, str):
        return value
    if abs(value) < 1e-12:
        value = 0.0
    return f"{value:.12g}"


def create_rectangle(geom, tag: str, pos: tuple[float, float], size: tuple[float, float], label: str) -> None:
    if size[0] <= 1e-9 or size[1] <= 1e-9:
        return
    rect = geom.create(tag, "Rectangle")
    rect.label(f"{tag}: {label}")
    rect.set("base", "corner")
    rect.set("pos", [fmt(pos[0]), fmt(pos[1])])
    rect.set("size", [fmt(size[0]), fmt(size[1])])
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
    if not inputs:
        raise ValueError(f"Selection {tag} has no inputs.")
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


def tile_segments(start: float, end: float, particle_w: float = PARTICLE_W_NM, gap_w: float = PARTICLE_GAP_NM) -> list[tuple[float, float]]:
    length = end - start
    if length < particle_w - 1e-9:
        return []
    n = max(1, int(math.floor((length + gap_w) / (particle_w + gap_w))))
    occupied = n * particle_w + (n - 1) * gap_w
    offset = max(0.0, (length - occupied) / 2.0)
    segments = []
    x = start + offset
    for _ in range(n):
        segments.append((x, x + particle_w))
        x += particle_w + gap_w
    return segments


def complement_segments(start: float, end: float, occupied: list[tuple[float, float]]) -> list[tuple[float, float]]:
    pieces = []
    cursor = start
    for left, right in sorted(occupied):
        if left - cursor > 1e-9:
            pieces.append((cursor, left))
        cursor = max(cursor, right)
    if end - cursor > 1e-9:
        pieces.append((cursor, end))
    return pieces


def side_segments(y_start: float, y_end: float) -> list[tuple[float, float]]:
    return tile_segments(y_start, y_end, SIDE_PARTICLE_H_NM, SIDE_PARTICLE_GAP_NM)


def sidewall_coverage() -> float:
    segs = side_segments(0.0, H_SIN_NM)
    return sum(b - a for a, b in segs) / H_SIN_NM if segs else 0.0


def particle_stats(design: Design) -> tuple[int, float, float, float]:
    ridge_left, ridge_right = -design.width_nm / 2.0, design.width_nm / 2.0
    x_min, x_max = -design.period_nm / 2.0, design.period_nm / 2.0
    top = tile_segments(ridge_left, ridge_right)
    left_bottom = tile_segments(x_min, ridge_left - H_RE_NM)
    right_bottom = tile_segments(ridge_right + H_RE_NM, x_max)
    side = side_segments(0.0, H_SIN_NM)
    top_cov = sum(b - a for a, b in top) / design.width_nm if design.width_nm > 0 else 0.0
    bottom_available = max(0.0, (ridge_left - H_RE_NM) - x_min) + max(0.0, x_max - (ridge_right + H_RE_NM))
    bottom_cov = (sum(b - a for a, b in left_bottom + right_bottom) / bottom_available) if bottom_available > 0 else 0.0
    return len(top) + len(left_bottom) + len(right_bottom) + 2 * len(side), top_cov, bottom_cov, sidewall_coverage()


def is_fabricable(design: Design) -> bool:
    return design.gap_nm >= MIN_GAP_NM and design.width_nm >= MIN_RIDGE_W_NM


def unique_complex(values: Iterable[complex | float], tol: float = 1e-6) -> list[complex]:
    result: list[complex] = []
    for raw in values:
        val = complex(raw)
        if not any(abs(val - old) < tol for old in result):
            result.append(val)
    return result


def score_mode(wavelength_nm: float, q_value: float) -> float:
    detune = abs(wavelength_nm - TARGET_NM)
    if q_value >= TARGET_Q:
        q_penalty = 0.0
    elif q_value >= ACCEPT_Q:
        q_penalty = 0.5 * math.log10(TARGET_Q / q_value)
    else:
        q_penalty = 4.0 * math.log10(ACCEPT_Q / max(q_value, 1.0))
    return detune + q_penalty


def mode_result(design: Design, re_style: str, kx_norm: float, idx: int, freq_thz: complex, re_overlap: float = math.nan) -> ModeResult:
    freq_hz_real = freq_thz.real * 1e12
    wavelength_nm = C0 / freq_hz_real * 1e9 if freq_hz_real > 0 else float("inf")
    q = abs(freq_thz.real / (2.0 * freq_thz.imag)) if abs(freq_thz.imag) > 0 else float("inf")
    count, top_cov, bottom_cov, side_cov = particle_stats(design)
    return ModeResult(
        design.design_id,
        re_style,
        design.period_nm,
        design.width_nm,
        design.gap_nm,
        kx_norm,
        idx,
        freq_thz.real,
        freq_thz.imag,
        wavelength_nm,
        q,
        count if re_style == "particulate" else 0,
        top_cov if re_style == "particulate" else (1.0 if re_style == "continuous" else 0.0),
        bottom_cov if re_style == "particulate" else (1.0 if re_style == "continuous" else 0.0),
        side_cov if re_style == "particulate" else (1.0 if re_style == "continuous" else 0.0),
        re_overlap,
        score_mode(wavelength_nm, q),
    )


def get_eigenfrequencies_thz(model) -> list[complex]:
    values = model.evaluate("freq", unit="THz")
    flat = values.ravel() if hasattr(values, "ravel") else values
    modes = unique_complex(flat, tol=1e-7)
    modes = [m for m in modes if math.isfinite(m.real) and abs(m.real) > 1e-9]
    modes.sort(key=lambda z: abs(C0 / (z.real * 1e12) * 1e9 - TARGET_NM) if z.real else float("inf"))
    return modes


def build_model(client, design: Design, re_style: str, name: str):
    model = client.create(name)
    jmodel = model.java
    x_min, x_max = -design.period_nm / 2.0, design.period_nm / 2.0
    ridge_left, ridge_right = -design.width_nm / 2.0, design.width_nm / 2.0
    y_si = H_PML_NM
    y_sio2 = H_PML_NM + H_SI_BUF_NM
    y_sin = y_sio2 + H_SIO2_NM
    y_sin_top = y_sin + H_SIN_NM
    y_air_main = y_sin_top + H_RE_NM
    y_top_pml = H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_SIN_NM + H_AIR_NM
    total_height = 2 * H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_SIN_NM + H_AIR_NM

    set_param(model, "lambda_target", f"{TARGET_NM}[nm]", "Target wavelength")
    set_param(model, "P", f"{design.period_nm:.9g}[nm]", "Grating period")
    set_param(model, "W", f"{design.width_nm:.9g}[nm]", "SiN ridge width")
    set_param(model, "gap", f"{design.gap_nm:.9g}[nm]", "Etched gap")
    set_param(model, "h_sin", "300[nm]", "SiN height")
    set_param(model, "h_re", "16[nm]", "Particle thickness")
    set_param(model, "particle_w", "25[nm]", "Particle width")
    set_param(model, "particle_gap", "8[nm]", "Particle spacing")
    set_param(model, "h_sio2", "3000[nm]", "SiO2 thickness")
    set_param(model, "n_sin", f"{N_SIN:.9f}", "SiN index at 451 nm")
    set_param(model, "n_sio2", f"{N_SIO2:.9f}", "SiO2 index at 451 nm")
    set_param(model, "n_si", f"{N_SI:.9f}", "Effective Si index at 451 nm")
    set_param(model, "n_re", f"{N_RE:.9f}", "Particulate rare-earth proxy index")
    set_param(model, "kx", "0[1/m]", "Floquet wave vector along grating period")

    jmodel.component().create("comp1", True)
    comp = jmodel.component("comp1")
    comp.label("2D 451 nm particulate-RE unit cell")
    comp.geom().create("geom1", 2)
    geom = comp.geom("geom1")
    geom.lengthUnit("nm")

    create_rectangle(geom, "si_pml", (x_min, 0.0), (design.period_nm, H_PML_NM), "Si bottom PML")
    create_rectangle(geom, "si_buf", (x_min, y_si), (design.period_nm, H_SI_BUF_NM), "Si buffer")
    create_rectangle(geom, "sio2", (x_min, y_sio2), (design.period_nm, H_SIO2_NM), "SiO2")
    create_rectangle(geom, "sin_ridge", (ridge_left, y_sin), (design.width_nm, H_SIN_NM), "300 nm SiN ridge")

    re_tags: list[str] = []
    air_tags: list[str] = []

    def add_re(tag: str, left: float, bottom: float, width: float, height: float, label: str) -> None:
        create_rectangle(geom, tag, (left, bottom), (width, height), label)
        re_tags.append(f"geom1_{tag}_dom")

    def add_air(tag: str, left: float, bottom: float, width: float, height: float, label: str) -> None:
        create_rectangle(geom, tag, (left, bottom), (width, height), label)
        air_tags.append(f"geom1_{tag}_dom")

    if re_style == "continuous":
        add_re("re_top", ridge_left, y_sin_top, design.width_nm, H_RE_NM, "continuous top RE")
        add_re("re_side_l", ridge_left - H_RE_NM, y_sin, H_RE_NM, H_SIN_NM, "continuous left sidewall RE")
        add_re("re_side_r", ridge_right, y_sin, H_RE_NM, H_SIN_NM, "continuous right sidewall RE")
        add_re("re_bottom_l", x_min, y_sin, max(0.0, ridge_left - H_RE_NM - x_min), H_RE_NM, "continuous left floor RE")
        add_re("re_bottom_r", ridge_right + H_RE_NM, y_sin, max(0.0, x_max - (ridge_right + H_RE_NM)), H_RE_NM, "continuous right floor RE")
        add_air("air_top_l", x_min, y_sin_top, ridge_left - x_min, H_RE_NM, "air cap left")
        add_air("air_top_r", ridge_right, y_sin_top, x_max - ridge_right, H_RE_NM, "air cap right")
        add_air("air_slot_l", x_min, y_sin + H_RE_NM, ridge_left - H_RE_NM - x_min, H_SIN_NM - H_RE_NM, "left slot air")
        add_air("air_slot_r", ridge_right + H_RE_NM, y_sin + H_RE_NM, x_max - (ridge_right + H_RE_NM), H_SIN_NM - H_RE_NM, "right slot air")
    elif re_style == "none":
        add_air("air_top", x_min, y_sin_top, design.period_nm, H_RE_NM, "air over grating")
        add_air("air_slot_l", x_min, y_sin, ridge_left - x_min, H_SIN_NM, "left slot air")
        add_air("air_slot_r", ridge_right, y_sin, x_max - ridge_right, H_SIN_NM, "right slot air")
    else:
        top_particles = tile_segments(ridge_left, ridge_right)
        for i, (left, right) in enumerate(top_particles):
            add_re(f"re_top_{i}", left, y_sin_top, right - left, H_RE_NM, "top RE particle")
        for i, (left, right) in enumerate(complement_segments(ridge_left, ridge_right, top_particles)):
            add_air(f"air_top_between_{i}", left, y_sin_top, right - left, H_RE_NM, "top particle gap air")
        add_air("air_top_l", x_min, y_sin_top, ridge_left - x_min, H_RE_NM, "left top air")
        add_air("air_top_r", ridge_right, y_sin_top, x_max - ridge_right, H_RE_NM, "right top air")

        side_l_x0, side_l_x1 = ridge_left - H_RE_NM, ridge_left
        side_r_x0, side_r_x1 = ridge_right, ridge_right + H_RE_NM
        side_abs = [(y_sin + a, y_sin + b) for a, b in side_segments(0.0, H_SIN_NM)]
        for i, (bottom, top) in enumerate(side_abs):
            add_re(f"re_side_l_{i}", side_l_x0, bottom, H_RE_NM, top - bottom, "sparse left sidewall particle")
            add_re(f"re_side_r_{i}", side_r_x0, bottom, H_RE_NM, top - bottom, "sparse right sidewall particle")
        for i, (bottom, top) in enumerate(complement_segments(y_sin, y_sin_top, side_abs)):
            add_air(f"air_side_l_{i}", side_l_x0, bottom, H_RE_NM, top - bottom, "left sidewall air gap")
            add_air(f"air_side_r_{i}", side_r_x0, bottom, H_RE_NM, top - bottom, "right sidewall air gap")

        left_floor_start, left_floor_end = x_min, max(x_min, ridge_left - H_RE_NM)
        right_floor_start, right_floor_end = min(x_max, ridge_right + H_RE_NM), x_max
        left_bottom = tile_segments(left_floor_start, left_floor_end)
        right_bottom = tile_segments(right_floor_start, right_floor_end)
        for i, (left, right) in enumerate(left_bottom):
            add_re(f"re_bottom_l_{i}", left, y_sin, right - left, H_RE_NM, "left floor RE particle")
        for i, (left, right) in enumerate(right_bottom):
            add_re(f"re_bottom_r_{i}", left, y_sin, right - left, H_RE_NM, "right floor RE particle")
        for i, (left, right) in enumerate(complement_segments(left_floor_start, left_floor_end, left_bottom)):
            add_air(f"air_bottom_l_{i}", left, y_sin, right - left, H_RE_NM, "left floor particle gap air")
        for i, (left, right) in enumerate(complement_segments(right_floor_start, right_floor_end, right_bottom)):
            add_air(f"air_bottom_r_{i}", left, y_sin, right - left, H_RE_NM, "right floor particle gap air")

        add_air("air_slot_l_bulk", x_min, y_sin + H_RE_NM, max(0.0, ridge_left - H_RE_NM - x_min), H_SIN_NM - H_RE_NM, "left slot bulk air")
        add_air("air_slot_r_bulk", ridge_right + H_RE_NM, y_sin + H_RE_NM, max(0.0, x_max - (ridge_right + H_RE_NM)), H_SIN_NM - H_RE_NM, "right slot bulk air")

    add_air("air", x_min, y_air_main, design.period_nm, H_AIR_NM - H_RE_NM, "air")
    add_air("air_pml", x_min, y_top_pml, design.period_nm, H_PML_NM, "air top PML")
    geom.run()

    create_union_selection(comp, "air_dom", "2", air_tags)
    if re_tags:
        create_union_selection(comp, "rare_earth_dom", "2", re_tags)
    create_union_selection(comp, "si_dom", "2", ["geom1_si_pml_dom", "geom1_si_buf_dom"])
    create_box_selection(comp, "left_bnd", "1", f"{x_min - 1}", f"{x_min + 1}", "-1", f"{total_height + 1}")
    create_box_selection(comp, "right_bnd", "1", f"{x_max - 1}", f"{x_max + 1}", "-1", f"{total_height + 1}")
    create_union_selection(comp, "periodic_bnd", "1", ["left_bnd", "right_bnd"])

    pml_bot = comp.coordSystem().create("pml_bot", "PML")
    pml_bot.label("Bottom Si PML")
    pml_bot.selection().named("geom1_si_pml_dom")
    pml_top = comp.coordSystem().create("pml_top", "PML")
    pml_top.label("Top air PML")
    pml_top.selection().named("geom1_air_pml_dom")

    create_refractive_material(comp, "mat_air", "Air", "air_dom", "1", "0")
    create_refractive_material(comp, "mat_sio2", "SiO2 451 nm", "geom1_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_sin", "SiN 451 nm", "geom1_sin_ridge_dom", "n_sin", "0")
    if re_tags:
        create_refractive_material(comp, "mat_re", "Particulate rare-earth proxy n=1.52", "rare_earth_dom", "n_re", "0")
    create_refractive_material(comp, "mat_si", "Si effective 451 nm", "si_dom", "n_si", "0")

    ewfd = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
    ewfd.label("Electromagnetic Waves, Frequency Domain")
    pc = ewfd.create("pc1", "PeriodicCondition", 1)
    pc.label("Floquet periodic x")
    pc.selection().named("periodic_bnd")
    pc.set("PeriodicType", "Floquet")
    pc.set("kFloquet", ["kx", "0", "0"])

    mesh = comp.mesh().create("mesh1")
    mesh.label("Physics-controlled mesh")
    mesh.autoMeshSize(2)

    study = jmodel.study().create("std1")
    study.label("Eigenfrequency near 451 nm")
    eig = study.create("eig", "Eigenfrequency")
    eig.set("shift", "c_const/lambda_target")
    eig.set("neigs", "12")
    eig.set("neigsmanual", "12")
    eig.set("neigsactive", "on")
    eig.set("eigunit", "THz")
    return model


def solve_design(model, design: Design, re_style: str, kx_norm: float = 0.0) -> list[ModeResult]:
    set_param(model, "kx", f"{kx_norm:.12g}*pi/P")
    model.java.component("comp1").geom("geom1").run()
    model.java.component("comp1").mesh("mesh1").run()
    model.java.study("std1").run()
    modes = get_eigenfrequencies_thz(model)
    return [mode_result(design, re_style, kx_norm, i + 1, freq) for i, freq in enumerate(modes)]


def best_mode(rows: list[ModeResult]) -> ModeResult:
    return sorted(rows, key=lambda r: (r.score, abs(r.wavelength_nm - TARGET_NM), -r.q_value))[0]


def safe_model_name(text: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in text)[:55]


def coarse_designs() -> list[Design]:
    designs: dict[str, Design] = {}
    for period in [306, 340, 380, 430, 480, 540, 600, 680, 760, 840, 900]:
        for gap in [130, 150, 180, 220, 280, 340]:
            width = period - gap
            design = Design(f"coarse_P{period:.0f}_G{gap:.0f}", float(period), float(width))
            if is_fabricable(design):
                designs[design.design_id] = design
    return list(designs.values())


def fine_designs(seed: ModeResult) -> list[Design]:
    designs: dict[str, Design] = {}
    for dp in [-12, -8, -4, -2, 0, 2, 4, 8, 12]:
        for dg in [-20, -10, -5, 0, 5, 10, 20]:
            period = seed.period_nm + dp
            gap = max(MIN_GAP_NM, seed.gap_nm + dg)
            width = period - gap
            design = Design(f"fine_P{period:.3f}_G{gap:.3f}", period, width)
            if is_fabricable(design):
                designs[design.design_id] = design
    return list(designs.values())


def write_modes_csv(path: Path, rows: list[ModeResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "design_id",
                "re_style",
                "period_nm",
                "width_nm",
                "gap_nm",
                "kx_norm",
                "mode_index",
                "freq_thz_real",
                "freq_thz_imag",
                "wavelength_nm",
                "q_value",
                "particle_count",
                "top_coverage",
                "bottom_coverage",
                "sidewall_coverage",
                "re_overlap",
                "score",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r.design_id,
                    r.re_style,
                    r.period_nm,
                    r.width_nm,
                    r.gap_nm,
                    r.kx_norm,
                    r.mode_index,
                    r.freq_thz_real,
                    r.freq_thz_imag,
                    r.wavelength_nm,
                    r.q_value,
                    r.particle_count,
                    r.top_coverage,
                    r.bottom_coverage,
                    r.sidewall_coverage,
                    r.re_overlap,
                    r.score,
                ]
            )


def read_modes_csv(path: Path) -> list[ModeResult]:
    rows: list[ModeResult] = []
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                ModeResult(
                    raw["design_id"],
                    raw["re_style"],
                    float(raw["period_nm"]),
                    float(raw["width_nm"]),
                    float(raw["gap_nm"]),
                    float(raw["kx_norm"]),
                    int(raw["mode_index"]),
                    float(raw["freq_thz_real"]),
                    float(raw["freq_thz_imag"]),
                    float(raw["wavelength_nm"]),
                    float(raw["q_value"]),
                    int(float(raw["particle_count"])),
                    float(raw["top_coverage"]),
                    float(raw["bottom_coverage"]),
                    float(raw["sidewall_coverage"]),
                    float(raw["re_overlap"]) if raw["re_overlap"] not in ("", "nan") else math.nan,
                    float(raw["score"]),
                )
            )
    return rows


def design_from_result(row: ModeResult) -> Design:
    return Design(row.design_id, row.period_nm, row.width_nm)


def scan_designs(client, designs: list[Design], path: Path, label: str) -> list[ModeResult]:
    rows = read_modes_csv(path)
    done = {r.design_id for r in rows}
    for i, design in enumerate(designs, 1):
        if design.design_id in done:
            continue
        print(f"[{label}] {i}/{len(designs)} solving P={design.period_nm:.3f} nm W={design.width_nm:.3f} nm gap={design.gap_nm:.3f} nm", flush=True)
        model = build_model(client, design, "particulate", safe_model_name(f"v5_{design.design_id}"))
        try:
            current = solve_design(model, design, "particulate", 0.0)
            rows.extend(current)
            write_modes_csv(path, rows)
            best = best_mode(current)
            print(f"  best lambda={best.wavelength_nm:.6f} nm Q={best.q_value:.4g} mode={best.mode_index}", flush=True)
        except Exception as exc:
            print(f"  FAILED {design.design_id}: {exc}", flush=True)
        finally:
            client.remove(model.name())
    return rows


def choose_candidates(rows: list[ModeResult], limit: int = 12) -> list[ModeResult]:
    best_by_design: dict[str, ModeResult] = {}
    for design_id in sorted({r.design_id for r in rows}):
        subset = [r for r in rows if r.design_id == design_id]
        best_by_design[design_id] = best_mode(subset)
    ranked = sorted(best_by_design.values(), key=lambda r: (r.score, abs(r.wavelength_nm - TARGET_NM), -r.q_value))
    highest_q = sorted(best_by_design.values(), key=lambda r: -r.q_value)[:3]
    closest = sorted(best_by_design.values(), key=lambda r: abs(r.wavelength_nm - TARGET_NM))[:3]
    merged: dict[str, ModeResult] = {}
    for row in ranked[:limit] + highest_q + closest:
        merged[f"{row.design_id}:{row.mode_index}"] = row
    return sorted(merged.values(), key=lambda r: (r.score, -r.q_value))[:limit]


def find_solution_index_for_mode(freqs: list[complex], target: ModeResult) -> int:
    def metric(item: tuple[int, complex]) -> float:
        _, freq = item
        wl = C0 / (freq.real * 1e12) * 1e9
        q = abs(freq.real / (2.0 * freq.imag)) if abs(freq.imag) > 0 else float("inf")
        q_part = 0.0
        if math.isfinite(q) and math.isfinite(target.q_value):
            q_part = 1e-6 * abs(math.log10(max(q, 1.0)) - math.log10(max(target.q_value, 1.0)))
        return abs(wl - target.wavelength_nm) + q_part

    return min(enumerate(freqs, start=1), key=metric)[0]


def point_in_particle_re(design: Design, x: float, y: float) -> bool:
    ridge_left, ridge_right = -design.width_nm / 2.0, design.width_nm / 2.0
    x_min, x_max = -design.period_nm / 2.0, design.period_nm / 2.0
    y_sin = H_PML_NM + H_SI_BUF_NM + H_SIO2_NM
    y_top = y_sin + H_SIN_NM
    for left, right in tile_segments(ridge_left, ridge_right):
        if left <= x <= right and y_top <= y <= y_top + H_RE_NM:
            return True
    for left, right in tile_segments(x_min, ridge_left - H_RE_NM):
        if left <= x <= right and y_sin <= y <= y_sin + H_RE_NM:
            return True
    for left, right in tile_segments(ridge_right + H_RE_NM, x_max):
        if left <= x <= right and y_sin <= y <= y_sin + H_RE_NM:
            return True
    for bottom, top in [(y_sin + a, y_sin + b) for a, b in side_segments(0.0, H_SIN_NM)]:
        if ridge_left - H_RE_NM <= x <= ridge_left and bottom <= y <= top:
            return True
        if ridge_right <= x <= ridge_right + H_RE_NM and bottom <= y <= top:
            return True
    return False


def compute_re_overlap(model, design: Design, target: ModeResult) -> float:
    import numpy as np

    freqs = [complex(v) for v in np.asarray(model.evaluate("freq", unit="THz")).ravel()]
    sol_idx = find_solution_index_for_mode(freqs, target)
    ds = model.datasets()[0] if model.datasets() else None
    x = np.asarray(model.evaluate("x", unit="nm", dataset=ds, inner=[sol_idx])).ravel().astype(float)
    y = np.asarray(model.evaluate("y", unit="nm", dataset=ds, inner=[sol_idx])).ravel().astype(float)
    ez = np.asarray(model.evaluate("abs(ewfd.Ez)", dataset=ds, inner=[sol_idx])).ravel().astype(float)
    n = min(len(x), len(y), len(ez))
    x, y, e2 = x[:n], y[:n], ez[:n] * ez[:n]
    y_sin = H_PML_NM + H_SI_BUF_NM + H_SIO2_NM
    total = 0.0
    re_sum = 0.0
    for xx, yy, val in zip(x, y, e2):
        if y_sin - 200.0 <= yy <= y_sin + H_SIN_NM + 350.0:
            total += float(val)
            if point_in_particle_re(design, float(xx), float(yy)):
                re_sum += float(val)
    return re_sum / total if total > 0 else math.nan


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


def svg_line_plot(path: Path, title: str, xs: list[float], ys: list[float], y_label: str, x_label: str, log_y: bool = False) -> None:
    width, height = 900, 530
    ml, mr, mt, mb = 92, 35, 55, 75
    plot_w, plot_h = width - ml - mr, height - mt - mb
    valid = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y) and (not log_y or y > 0)]
    if not valid:
        path.write_text("<svg xmlns='http://www.w3.org/2000/svg'></svg>", encoding="utf-8")
        return
    x_min, x_max = min(x for x, _ in valid), max(x for x, _ in valid)
    y_vals = [math.log10(y) if log_y else y for _, y in valid]
    y_min, y_max = min(y_vals), max(y_vals)
    if abs(x_max - x_min) < 1e-12:
        x_max += 1.0
    if abs(y_max - y_min) < 1e-12:
        y_max += 1.0

    def sx(x: float) -> float:
        return ml + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        yy = math.log10(max(y, 1e-300)) if log_y else y
        return mt + (y_max - yy) / (y_max - y_min) * plot_h

    points = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in valid)
    out = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2}' y='30' text-anchor='middle' font-size='22' font-family='Arial'>{title}</text>",
        f"<line x1='{ml}' y1='{mt+plot_h}' x2='{ml+plot_w}' y2='{mt+plot_h}' stroke='black'/>",
        f"<line x1='{ml}' y1='{mt}' x2='{ml}' y2='{mt+plot_h}' stroke='black'/>",
        f"<text x='{width/2}' y='{height-25}' text-anchor='middle' font-size='16' font-family='Arial'>{x_label}</text>",
        f"<text x='24' y='{height/2}' text-anchor='middle' transform='rotate(-90 24 {height/2})' font-size='16' font-family='Arial'>{y_label}</text>",
        f"<polyline points='{points}' fill='none' stroke='#0b5fff' stroke-width='2.4'/>",
        "</svg>",
    ]
    path.write_text("\n".join(out), encoding="utf-8")


def svg_field_plot(path: Path, title: str, design: Design, x_nm: list[float], y_nm: list[float], values: list[float]) -> None:
    width, height = 980, 650
    ml, mr, mt, mb = 82, 38, 54, 74
    plot_w, plot_h = width - ml - mr, height - mt - mb
    y_sin = H_PML_NM + H_SI_BUF_NM + H_SIO2_NM
    y_min_focus = y_sin - 220.0
    y_max_focus = y_sin + H_SIN_NM + 420.0
    points = [(x, y, v) for x, y, v in zip(x_nm, y_nm, values) if -design.period_nm / 2 <= x <= design.period_nm / 2 and y_min_focus <= y <= y_max_focus]
    if len(points) > 6500:
        step = max(1, len(points) // 6500)
        points = points[::step]
    if not points:
        points = list(zip(x_nm, y_nm, values))
    x_min, x_max = -design.period_nm / 2.0, design.period_nm / 2.0
    y_min, y_max = min(y for _, y, _ in points), max(y for _, y, _ in points)
    vmax = max(max(v for _, _, v in points), 1e-300)

    def sx(x: float) -> float:
        return ml + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        return mt + (y_max - y) / (y_max - y_min) * plot_h

    def norm(v: float) -> float:
        return max(0.0, min(1.0, (math.log10(v / vmax + 1e-12) + 12.0) / 12.0))

    ridge_left, ridge_right = -design.width_nm / 2.0, design.width_nm / 2.0
    out = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2}' y='30' text-anchor='middle' font-size='21' font-family='Arial'>{title}</text>",
        f"<rect x='{ml}' y='{mt}' width='{plot_w}' height='{plot_h}' fill='#05051f' stroke='black'/>",
    ]
    for x, y, v in points:
        out.append(f"<circle cx='{sx(x):.2f}' cy='{sy(y):.2f}' r='1.5' fill='{plasma_color(norm(v))}' opacity='0.9'/>")
    out.append(f"<rect x='{sx(ridge_left):.2f}' y='{sy(y_sin+H_SIN_NM):.2f}' width='{sx(ridge_right)-sx(ridge_left):.2f}' height='{sy(y_sin)-sy(y_sin+H_SIN_NM):.2f}' fill='none' stroke='white' stroke-width='1.8'/>")
    for yy, label in [(y_sin, "SiO2 / SiN"), (y_sin + H_SIN_NM, "SiN / air")]:
        if y_min <= yy <= y_max:
            out.append(f"<line x1='{ml}' y1='{sy(yy):.2f}' x2='{ml+plot_w}' y2='{sy(yy):.2f}' stroke='white' stroke-width='1.2' stroke-dasharray='6 4'/>")
            out.append(f"<text x='{ml+8}' y='{sy(yy)-6:.2f}' font-size='12' font-family='Arial' fill='white'>{label}</text>")
    out.append(f"<text x='{ml+8}' y='{mt+18}' font-size='12' font-family='Arial' fill='white'>white: SiN outline; particle mask in CSV</text>")
    out.append("</svg>")
    path.write_text("\n".join(out), encoding="utf-8")


def export_field(model, design: Design, target: ModeResult) -> ModeResult:
    import numpy as np

    freqs = [complex(v) for v in np.asarray(model.evaluate("freq", unit="THz")).ravel()]
    sol_idx = find_solution_index_for_mode(freqs, target)
    ds = model.datasets()[0] if model.datasets() else None
    x = np.asarray(model.evaluate("x", unit="nm", dataset=ds, inner=[sol_idx])).ravel().astype(float)
    y = np.asarray(model.evaluate("y", unit="nm", dataset=ds, inner=[sol_idx])).ravel().astype(float)
    ez = np.asarray(model.evaluate("abs(ewfd.Ez)", dataset=ds, inner=[sol_idx])).ravel().astype(float)
    n = min(len(x), len(y), len(ez))
    x, y, e2 = x[:n], y[:n], ez[:n] * ez[:n]
    overlap = compute_re_overlap(model, design, target)
    target.re_overlap = overlap
    csv_path = OUT / "field_451_particle_ez2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x_nm", "y_nm", "abs_ez_squared", "in_particle_re"])
        for xx, yy, val in zip(x, y, e2):
            writer.writerow([xx, yy, val, point_in_particle_re(design, float(xx), float(yy))])
    svg_field_plot(OUT / "field_451_particle_ez2.svg", f"451 nm |Ez|^2, {design.design_id}", design, x.tolist(), y.tolist(), e2.tolist())
    return target


def save_best_model(client, best: ModeResult) -> ModeResult:
    design = design_from_result(best)
    model = build_model(client, design, "particulate", "sin451_qbic_v5_best_particulate_re")
    try:
        rows = solve_design(model, design, "particulate", 0.0)
        selected = min(rows, key=lambda r: abs(r.wavelength_nm - best.wavelength_nm) + 1e-6 * abs(math.log10(max(r.q_value, 1.0)) - math.log10(max(best.q_value, 1.0))))
        selected = export_field(model, design, selected)
        model.save(str(OUT / "sin451_qbic_v5_best_particulate_re.mph"))
        return selected
    finally:
        client.remove(model.name())


def run_compare(client, best: ModeResult) -> list[ModeResult]:
    design = design_from_result(best)
    rows: list[ModeResult] = []
    for style in ["continuous", "none"]:
        model = build_model(client, design, style, safe_model_name(f"compare_{style}_{design.design_id}"))
        try:
            current = solve_design(model, design, style, 0.0)
            selected = best_mode(current)
            rows.extend(current)
            if style == "continuous":
                model.save(str(OUT / "sin451_qbic_v5_compare_continuous_re.mph"))
            if style == "none":
                model.save(str(OUT / "sin451_qbic_v5_compare_no_re.mph"))
            print(f"[compare {style}] lambda={selected.wavelength_nm:.6f} nm Q={selected.q_value:.4g}", flush=True)
        finally:
            client.remove(model.name())
    write_modes_csv(OUT / "compare_re_styles.csv", rows)
    return rows


def run_band(client, best: ModeResult) -> list[ModeResult]:
    design = design_from_result(best)
    model = build_model(client, design, "particulate", "sin451_qbic_v5_band")
    rows: list[ModeResult] = []
    try:
        for kx in [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0]:
            print(f"[band] kx={kx:.3g} pi/P", flush=True)
            current = solve_design(model, design, "particulate", kx)
            rows.append(best_mode(current))
        write_modes_csv(OUT / "band_451_particle.csv", rows)
        svg_line_plot(OUT / "band_451_particle_wavelength.svg", "451 nm band", [r.kx_norm for r in rows], [r.wavelength_nm for r in rows], "wavelength (nm)", "kx / (pi/P)")
        svg_line_plot(OUT / "band_451_particle_q.svg", "451 nm Q along Gamma-X", [r.kx_norm for r in rows], [r.q_value for r in rows], "Q", "kx / (pi/P)", log_y=True)
    finally:
        client.remove(model.name())
    return rows


def write_particle_geometry_check(best: ModeResult) -> None:
    design = design_from_result(best)
    count, top_cov, bottom_cov, side_cov = particle_stats(design)
    with (OUT / "particle_geometry_check.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["period_nm", "width_nm", "gap_nm", "particle_w_nm", "particle_gap_nm", "particle_thickness_nm", "side_particle_h_nm", "side_particle_gap_nm", "particle_count", "top_coverage", "bottom_coverage", "sidewall_coverage"])
        writer.writerow([design.period_nm, design.width_nm, design.gap_nm, PARTICLE_W_NM, PARTICLE_GAP_NM, H_RE_NM, SIDE_PARTICLE_H_NM, SIDE_PARTICLE_GAP_NM, count, top_cov, bottom_cov, side_cov])


def plot_scan_summary(rows: list[ModeResult]) -> None:
    best_by_design = []
    for design_id in sorted({r.design_id for r in rows}):
        best_by_design.append(best_mode([r for r in rows if r.design_id == design_id]))
    ranked = sorted(best_by_design, key=lambda r: (r.score, -r.q_value))[:40]
    xs = list(range(1, len(ranked) + 1))
    svg_line_plot(OUT / "scan_451_wavelength.svg", "451 nm best candidates wavelength", xs, [r.wavelength_nm for r in ranked], "wavelength (nm)", "candidate rank")
    svg_line_plot(OUT / "scan_451_q.svg", "451 nm best candidates Q", xs, [r.q_value for r in ranked], "Q", "candidate rank", log_y=True)


def write_report(coarse: list[ModeResult], fine: list[ModeResult], best: ModeResult, compare: list[ModeResult]) -> None:
    all_rows = coarse + fine
    accepted_1e9 = [r for r in all_rows if abs(r.wavelength_nm - TARGET_NM) <= 0.5 and r.q_value >= TARGET_Q]
    accepted_1e8 = [r for r in all_rows if abs(r.wavelength_nm - TARGET_NM) <= 1.0 and r.q_value >= ACCEPT_Q]
    cont = best_mode([r for r in compare if r.re_style == "continuous"]) if any(r.re_style == "continuous" for r in compare) else None
    none = best_mode([r for r in compare if r.re_style == "none"]) if any(r.re_style == "none" for r in compare) else None
    lines = [
        "# v5 451 nm particulate rare-earth quasi-BIC summary",
        "",
        "All files were generated inside `E:\\xitugrating\\sin451_qbic_2d_design_v5_particulate_re`; previous versions and original files were not edited.",
        "",
        "## Scan result",
        "",
        f"- Coarse mode rows: `{len(coarse)}`.",
        f"- Fine mode rows: `{len(fine)}`.",
        f"- Candidates meeting `|lambda-451| <= 0.5 nm` and `Q >= 1e9`: `{len(accepted_1e9)}`.",
        f"- Candidates meeting `|lambda-451| <= 1.0 nm` and `Q >= 1e8`: `{len(accepted_1e8)}`.",
        "",
        "## Best particulate-RE candidate",
        "",
        "| item | value |",
        "|---|---:|",
        f"| P (nm) | `{best.period_nm:.9g}` |",
        f"| W (nm) | `{best.width_nm:.9g}` |",
        f"| etched gap (nm) | `{best.gap_nm:.9g}` |",
        f"| lambda451 (nm) | `{best.wavelength_nm:.9g}` |",
        f"| Q451 | `{best.q_value:.6g}` |",
        f"| mode index | `{best.mode_index}` |",
        f"| particle count | `{best.particle_count}` |",
        f"| top coverage | `{best.top_coverage:.4g}` |",
        f"| bottom coverage | `{best.bottom_coverage:.4g}` |",
        f"| sidewall coverage | `{best.sidewall_coverage:.4g}` |",
        f"| particle RE field-overlap proxy | `{best.re_overlap:.6g}` |",
        "",
    ]
    if best.q_value >= TARGET_Q and abs(best.wavelength_nm - TARGET_NM) <= 0.5:
        lines.append("Conclusion: the particulate-RE 451 nm design meets the preferred `Q >= 1e9` target.")
    elif best.q_value >= ACCEPT_Q and abs(best.wavelength_nm - TARGET_NM) <= 1.0:
        lines.append("Conclusion: the particulate-RE 451 nm design meets the relaxed `Q >= 1e8` feasibility target, but not the preferred `1e9` target.")
    else:
        lines.append("Conclusion: no scanned particulate-RE design met the relaxed `Q >= 1e8` feasibility target near 451 nm; the saved model is the best proof-of-search candidate.")
    lines.extend(["", "## RE style comparison", ""])
    if cont:
        lines.append(f"- Continuous RE at same P/W: `lambda={cont.wavelength_nm:.6g} nm`, `Q={cont.q_value:.6g}`.")
    if none:
        lines.append(f"- No RE at same P/W: `lambda={none.wavelength_nm:.6g} nm`, `Q={none.q_value:.6g}`.")
    if cont and none:
        lines.append("- This comparison helps separate particle scattering/loading from the underlying SiN grating radiation loss.")
    lines.extend(
        [
            "",
            "## Output files",
            "",
            "- `sin451_qbic_v5_best_particulate_re.mph`",
            "- `coarse_scan_451_particle.csv`, `fine_scan_451_particle.csv`, `best_candidates_451_particle.csv`",
            "- `band_451_particle.csv`, `band_451_particle_wavelength.svg`, `band_451_particle_q.svg`",
            "- `field_451_particle_ez2.csv`, `field_451_particle_ez2.svg`",
            "- `particle_geometry_check.csv`, `compare_re_styles.csv`",
        ]
    )
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def run_all(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = import_mph()
    client = mph.Client(cores=args.cores)
    try:
        coarse_candidates = coarse_designs()
        if args.limit_coarse:
            coarse_candidates = coarse_candidates[: args.limit_coarse]
        coarse_rows = scan_designs(client, coarse_candidates, OUT / "coarse_scan_451_particle.csv", "coarse")
        seeds = choose_candidates(coarse_rows, limit=args.fine_seeds)
        fine_map: dict[str, Design] = {}
        for seed in seeds:
            for design in fine_designs(seed):
                fine_map[design.design_id] = design
        fine_candidates = list(fine_map.values())
        if args.limit_fine:
            fine_candidates = fine_candidates[: args.limit_fine]
        fine_rows = scan_designs(client, fine_candidates, OUT / "fine_scan_451_particle.csv", "fine") if fine_candidates else []
        all_rows = coarse_rows + fine_rows
        candidates = choose_candidates(all_rows, limit=20)
        write_modes_csv(OUT / "best_candidates_451_particle.csv", candidates)
        plot_scan_summary(all_rows)
        best = candidates[0]
        best = save_best_model(client, best)
        write_particle_geometry_check(best)
        compare = run_compare(client, best)
        run_band(client, best)
        write_report(coarse_rows, fine_rows, best, compare)
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="451 nm particulate rare-earth SiN quasi-BIC scan.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-all", help="Run v5 scan and generate outputs.")
    p.add_argument("--cores", type=int, default=2)
    p.add_argument("--limit-coarse", type=int, default=0)
    p.add_argument("--limit-fine", type=int, default=0)
    p.add_argument("--fine-seeds", type=int, default=6)
    p.set_defaults(func=run_all)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
