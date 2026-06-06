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
WORKSPACE = ROOT.parent
INPUTS = ROOT / "inputs"
OUT = ROOT / "outputs"

COMSOLROOT = r"D:\comsol63\Multiphysics"
JAVA_HOME = rf"{COMSOLROOT}\java\win64\jre"
C0 = 299_792_458.0

TARGET_980_NM = 980.0
TARGET_451_NM = 451.0
N_SIN_980 = 2.069946018
N_SIN_451 = 2.128247126
N_RE = 1.49
N_SIO2_980 = 1.45
N_SIO2_451 = 1.466
N_SI_980 = 3.55
N_SI_451 = 4.7

P0_NM = 564.3
W0_NM = 383.724
H_SIN_NM = 300.0
H_RE_NM = 16.0
H_SIO2_NM = 3000.0
H_SI_BUF_NM = 500.0
H_PML_NM = 800.0
H_AIR_NM = 1200.0

MIN_FEATURE_NM = 100.0
MIN_REMAINING_GAP_NM = 100.0
MIN_Q_980 = 1e8
MIN_Q_451 = 100.0
MAX_DETUNE_980_NM = 0.5
MAX_DETUNE_451_NM = 5.0


@dataclass(frozen=True)
class Ridge:
    left_nm: float
    right_nm: float

    @property
    def width_nm(self) -> float:
        return self.right_nm - self.left_nm


@dataclass(frozen=True)
class Design:
    design_id: str
    family: str
    pitches_nm: tuple[float, ...]
    widths_nm: tuple[float, ...]
    groove_w_nm: float = 0.0
    groove_d_nm: float = 0.0

    @property
    def n_periods(self) -> int:
        return len(self.pitches_nm)

    @property
    def length_nm(self) -> float:
        return sum(self.pitches_nm)

    @property
    def has_top_groove(self) -> bool:
        return self.groove_w_nm > 0 and self.groove_d_nm > 0


@dataclass
class ModeResult:
    design_id: str
    family: str
    target_nm: float
    n_periods: int
    kx_norm: float
    mode_index: int
    freq_thz_real: float
    freq_thz_imag: float
    wavelength_nm: float
    q_value: float
    re_overlap: float
    score: float


@dataclass
class DualResult:
    design: Design
    mode980: ModeResult
    mode451: ModeResult
    min_feature_nm: float
    min_gap_nm: float
    min_remaining_gap_nm: float
    dual_score: float


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


def create_rectangle(geom, tag: str, pos: tuple[float, float] | tuple[str, str], size: tuple[float, float] | tuple[str, str], label: str) -> None:
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


def fmt(value: float | str) -> str:
    if isinstance(value, str):
        return value
    if abs(value) < 1e-12:
        value = 0.0
    return f"{value:.12g}"


def ridge_edges(design: Design) -> list[Ridge]:
    x = -design.length_nm / 2.0
    ridges: list[Ridge] = []
    for pitch, width in zip(design.pitches_nm, design.widths_nm):
        center = x + pitch / 2.0
        ridges.append(Ridge(center - width / 2.0, center + width / 2.0))
        x += pitch
    return ridges


def gap_widths(design: Design) -> list[float]:
    ridges = ridge_edges(design)
    gaps: list[float] = []
    for a, b in zip(ridges, ridges[1:]):
        gaps.append(b.left_nm - a.right_nm)
    wrap_gap = (design.length_nm / 2.0 - ridges[-1].right_nm) + (ridges[0].left_nm + design.length_nm / 2.0)
    gaps.append(wrap_gap)
    return gaps


def min_remaining_gap_nm(design: Design) -> float:
    return min(gap_widths(design)) - 2.0 * H_RE_NM


def min_feature_nm(design: Design) -> float:
    values = [min(design.widths_nm), min(gap_widths(design))]
    if design.has_top_groove:
        rib_values = [(w - design.groove_w_nm) / 2.0 for w in design.widths_nm]
        values.extend([design.groove_w_nm, design.groove_d_nm, H_SIN_NM - design.groove_d_nm, min(rib_values)])
    return min(values)


def is_fabricable(design: Design) -> bool:
    if min_feature_nm(design) < MIN_FEATURE_NM:
        return False
    if min_remaining_gap_nm(design) < MIN_REMAINING_GAP_NM:
        return False
    if design.has_top_groove and design.groove_w_nm - 2.0 * H_RE_NM < MIN_FEATURE_NM:
        return False
    return True


def is_folded_only(design: Design) -> bool:
    if design.n_periods <= 1:
        return False
    widths = design.widths_nm
    gaps = gap_widths(design)
    return max(widths) - min(widths) < 1e-6 and max(gaps) - min(gaps) < 1e-6 and not design.has_top_groove


def unique_complex(values: Iterable[complex | float], tol: float = 1e-6) -> list[complex]:
    result: list[complex] = []
    for raw in values:
        val = complex(raw)
        if not any(abs(val - old) < tol for old in result):
            result.append(val)
    return result


def get_eigenfrequencies_thz(model, target_nm: float) -> list[complex]:
    values = model.evaluate("freq", unit="THz")
    flat = values.ravel() if hasattr(values, "ravel") else values
    modes = unique_complex(flat, tol=1e-7)
    modes = [m for m in modes if math.isfinite(m.real) and abs(m.real) > 1e-9]
    modes.sort(key=lambda z: abs(C0 / (z.real * 1e12) * 1e9 - target_nm) if z.real else float("inf"))
    return modes


def mode_result(design: Design, target_nm: float, kx_norm: float, idx: int, freq_thz: complex, re_overlap: float = math.nan) -> ModeResult:
    freq_hz_real = freq_thz.real * 1e12
    wavelength_nm = C0 / freq_hz_real * 1e9 if freq_hz_real > 0 else float("inf")
    q = abs(freq_thz.real / (2.0 * freq_thz.imag)) if abs(freq_thz.imag) > 0 else float("inf")
    detune = abs(wavelength_nm - target_nm)
    if target_nm > 700:
        q_penalty = 0.0 if q >= MIN_Q_980 else math.log10(MIN_Q_980 / max(q, 1.0))
        score = detune + 5.0 * q_penalty
    else:
        if q < MIN_Q_451:
            q_penalty = math.log10(MIN_Q_451 / max(q, 1.0))
        elif q > 1e5:
            q_penalty = 0.2 * math.log10(q / 1e5)
        else:
            q_penalty = 0.0
        overlap_bonus = 0.0 if math.isnan(re_overlap) else min(re_overlap, 0.2)
        score = detune + 2.0 * q_penalty - overlap_bonus
    return ModeResult(
        design.design_id,
        design.family,
        target_nm,
        design.n_periods,
        kx_norm,
        idx,
        freq_thz.real,
        freq_thz.imag,
        wavelength_nm,
        q,
        re_overlap,
        score,
    )


def target_materials(target_nm: float) -> tuple[float, float, float]:
    if target_nm > 700:
        return N_SIN_980, N_SIO2_980, N_SI_980
    return N_SIN_451, N_SIO2_451, N_SI_451


def set_target(model, target_nm: float) -> None:
    n_sin, n_sio2, n_si = target_materials(target_nm)
    set_param(model, "lambda_target", f"{target_nm:.9g}[nm]")
    set_param(model, "n_sin", f"{n_sin:.9f}")
    set_param(model, "n_sio2", f"{n_sio2:.9f}")
    set_param(model, "n_si", f"{n_si:.9f}")


def build_model(client, design: Design, name: str):
    model = client.create(name)
    jmodel = model.java

    total_height = 2 * H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_SIN_NM + H_AIR_NM
    y_bot_pml = 0.0
    y_si = H_PML_NM
    y_sio2 = H_PML_NM + H_SI_BUF_NM
    y_sin = y_sio2 + H_SIO2_NM
    y_sin_top = y_sin + H_SIN_NM
    y_air_main = y_sin_top + H_RE_NM
    y_top_pml = H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_SIN_NM + H_AIR_NM
    x_min = -design.length_nm / 2.0
    x_max = design.length_nm / 2.0

    set_param(model, "lambda_980", "980[nm]", "Primary pump target wavelength")
    set_param(model, "lambda_451", "451[nm]", "Rare-earth emission target wavelength")
    set_param(model, "lambda_target", "980[nm]", "Current eigenfrequency search wavelength")
    set_param(model, "L", f"{design.length_nm:.9g}[nm]", "Supercell length")
    set_param(model, "kx", "0[1/m]", "Floquet wave vector along supercell")
    set_param(model, "h_sin", "300[nm]", "SiN height")
    set_param(model, "h_re", "16[nm]", "Rare-earth equivalent layer thickness")
    set_param(model, "h_sio2", "3000[nm]", "SiO2 thickness")
    set_param(model, "n_sin_980", f"{N_SIN_980:.9f}", "Measured SiN index at 980 nm")
    set_param(model, "n_sin_451", f"{N_SIN_451:.9f}", "Interpolated SiN index at 451 nm")
    set_param(model, "n_sin", f"{N_SIN_980:.9f}", "Active SiN refractive index")
    set_param(model, "n_re", f"{N_RE:.9f}", "Rare-earth equivalent index")
    set_param(model, "n_sio2", f"{N_SIO2_980:.9f}", "Active SiO2 refractive index")
    set_param(model, "n_si", f"{N_SI_980:.9f}", "Active effective Si refractive index")

    jmodel.component().create("comp1", True)
    comp = jmodel.component("comp1")
    comp.label("2D x-y supercell")
    comp.geom().create("geom1", 2)
    geom = comp.geom("geom1")
    geom.lengthUnit("nm")

    create_rectangle(geom, "si_pml", (x_min, y_bot_pml), (design.length_nm, H_PML_NM), "Si bottom PML")
    create_rectangle(geom, "si_buf", (x_min, y_si), (design.length_nm, H_SI_BUF_NM), "Si buffer")
    create_rectangle(geom, "sio2", (x_min, y_sio2), (design.length_nm, H_SIO2_NM), "SiO2")

    sin_tags: list[str] = []
    re_tags: list[str] = []
    air_tags: list[str] = []
    ridges = ridge_edges(design)

    for i, ridge in enumerate(ridges):
        width = ridge.width_nm
        center = 0.5 * (ridge.left_nm + ridge.right_nm)
        if design.has_top_groove:
            groove_left = center - design.groove_w_nm / 2.0
            groove_right = center + design.groove_w_nm / 2.0
            y_groove_bottom = y_sin_top - design.groove_d_nm
            bottom_h = H_SIN_NM - design.groove_d_nm
            left_rib_w = groove_left - ridge.left_nm
            right_rib_w = ridge.right_nm - groove_right
            pieces = [
                (f"sin{i}_bottom", ridge.left_nm, y_sin, width, bottom_h, "SiN bottom below wide top groove"),
                (f"sin{i}_rib_l", ridge.left_nm, y_groove_bottom, left_rib_w, design.groove_d_nm, "SiN left rib"),
                (f"sin{i}_rib_r", groove_right, y_groove_bottom, right_rib_w, design.groove_d_nm, "SiN right rib"),
            ]
            for tag, x0, y0, w, h, label in pieces:
                create_rectangle(geom, tag, (x0, y0), (w, h), label)
                sin_tags.append(f"geom1_{tag}_dom")
            for tag, x0, w in [
                (f"re{i}_top_l", ridge.left_nm, left_rib_w),
                (f"re{i}_top_r", groove_right, right_rib_w),
            ]:
                create_rectangle(geom, tag, (x0, y_sin_top), (w, H_RE_NM), "Rare-earth top on SiN rib")
                re_tags.append(f"geom1_{tag}_dom")
            create_rectangle(geom, f"re{i}_groove_floor", (groove_left + H_RE_NM, y_groove_bottom), (design.groove_w_nm - 2 * H_RE_NM, H_RE_NM), "Rare-earth on wide groove floor")
            create_rectangle(geom, f"re{i}_groove_wall_l", (groove_left, y_groove_bottom), (H_RE_NM, design.groove_d_nm), "Rare-earth groove left wall")
            create_rectangle(geom, f"re{i}_groove_wall_r", (groove_right - H_RE_NM, y_groove_bottom), (H_RE_NM, design.groove_d_nm), "Rare-earth groove right wall")
            re_tags.extend([f"geom1_re{i}_groove_floor_dom", f"geom1_re{i}_groove_wall_l_dom", f"geom1_re{i}_groove_wall_r_dom"])
            air_w = design.groove_w_nm - 2 * H_RE_NM
            air_h = design.groove_d_nm - H_RE_NM
            create_rectangle(geom, f"air{i}_groove", (groove_left + H_RE_NM, y_groove_bottom + H_RE_NM), (air_w, air_h), "Air in wide top groove")
            air_tags.append(f"geom1_air{i}_groove_dom")
        else:
            create_rectangle(geom, f"sin{i}", (ridge.left_nm, y_sin), (width, H_SIN_NM), "SiN ridge")
            sin_tags.append(f"geom1_sin{i}_dom")
            create_rectangle(geom, f"re{i}_top", (ridge.left_nm, y_sin_top), (width, H_RE_NM), "Rare-earth top layer")
            re_tags.append(f"geom1_re{i}_top_dom")

        create_rectangle(geom, f"re{i}_side_l", (ridge.left_nm - H_RE_NM, y_sin), (H_RE_NM, H_SIN_NM), "Rare-earth left sidewall")
        create_rectangle(geom, f"re{i}_side_r", (ridge.right_nm, y_sin), (H_RE_NM, H_SIN_NM), "Rare-earth right sidewall")
        re_tags.extend([f"geom1_re{i}_side_l_dom", f"geom1_re{i}_side_r_dom"])

    def add_gap(tag: str, left: float, right: float, shrink_left: bool, shrink_right: bool) -> None:
        gap = right - left
        x0 = left + (H_RE_NM if shrink_left else 0.0)
        x1 = right - (H_RE_NM if shrink_right else 0.0)
        if x1 <= x0:
            return
        create_rectangle(geom, f"re_gap_{tag}", (x0, y_sin), (x1 - x0, H_RE_NM), "Rare-earth on exposed SiO2 groove floor")
        create_rectangle(geom, f"air_slot_{tag}", (x0, y_sin + H_RE_NM), (x1 - x0, H_SIN_NM - H_RE_NM), "Remaining air gap")
        create_rectangle(geom, f"air_cap_{tag}", (left, y_sin_top), (gap, H_RE_NM), "Air cap above groove")
        re_tags.append(f"geom1_re_gap_{tag}_dom")
        air_tags.extend([f"geom1_air_slot_{tag}_dom", f"geom1_air_cap_{tag}_dom"])

    for i, (a, b) in enumerate(zip(ridges, ridges[1:])):
        add_gap(f"{i}", a.right_nm, b.left_nm, True, True)
    add_gap("left_edge", x_min, ridges[0].left_nm, False, True)
    add_gap("right_edge", ridges[-1].right_nm, x_max, True, False)

    create_rectangle(geom, "air", (x_min, y_air_main), (design.length_nm, H_AIR_NM - H_RE_NM), "Air")
    create_rectangle(geom, "air_pml", (x_min, y_top_pml), (design.length_nm, H_PML_NM), "Air top PML")
    air_tags.extend(["geom1_air_dom", "geom1_air_pml_dom"])
    geom.run()

    create_union_selection(comp, "sin_dom", "2", sin_tags)
    create_union_selection(comp, "rare_earth_dom", "2", re_tags)
    create_union_selection(comp, "air_dom", "2", air_tags)
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
    create_refractive_material(comp, "mat_sio2", "SiO2 target-dependent", "geom1_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_sin", "SiN target-dependent", "sin_dom", "n_sin", "0")
    create_refractive_material(comp, "mat_re", "Rare-earth equivalent n=1.49", "rare_earth_dom", "n_re", "0")
    create_refractive_material(comp, "mat_si", "Si effective target-dependent", "si_dom", "n_si", "0")

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
    study.label("Eigenfrequency near target wavelength")
    eig = study.create("eig", "Eigenfrequency")
    eig.set("shift", "c_const/lambda_target")
    eig.set("neigs", "10")
    eig.set("neigsmanual", "10")
    eig.set("neigsactive", "on")
    eig.set("eigunit", "THz")

    return model


def solve_target(model, design: Design, target_nm: float, kx_norm: float = 0.0) -> list[ModeResult]:
    set_target(model, target_nm)
    set_param(model, "kx", f"{kx_norm:.12g}*pi/L")
    model.java.component("comp1").geom("geom1").run()
    model.java.component("comp1").mesh("mesh1").run()
    model.java.study("std1").run()
    modes = get_eigenfrequencies_thz(model, target_nm)
    return [mode_result(design, target_nm, kx_norm, i + 1, freq) for i, freq in enumerate(modes)]


def best_mode_for_target(results: list[ModeResult], target_nm: float) -> ModeResult:
    return sorted([r for r in results if abs(r.target_nm - target_nm) < 1e-9], key=lambda r: (r.score, abs(r.wavelength_nm - target_nm), -r.q_value))[0]


def dual_score(mode980: ModeResult, mode451: ModeResult) -> float:
    detune980 = abs(mode980.wavelength_nm - TARGET_980_NM)
    detune451 = abs(mode451.wavelength_nm - TARGET_451_NM)
    q980_penalty = 0.0 if mode980.q_value >= MIN_Q_980 else 20.0 * math.log10(MIN_Q_980 / max(mode980.q_value, 1.0))
    q451_penalty = 0.0 if mode451.q_value >= MIN_Q_451 else 5.0 * math.log10(MIN_Q_451 / max(mode451.q_value, 1.0))
    high_q451_penalty = 0.0 if mode451.q_value <= 1e5 else 0.1 * math.log10(mode451.q_value / 1e5)
    return detune980 * 4.0 + detune451 + q980_penalty + q451_penalty + high_q451_penalty


def passes_dual_targets(row: DualResult) -> bool:
    return (
        not is_folded_only(row.design)
        and
        abs(row.mode980.wavelength_nm - TARGET_980_NM) <= MAX_DETUNE_980_NM
        and row.mode980.q_value >= MIN_Q_980
        and abs(row.mode451.wavelength_nm - TARGET_451_NM) <= MAX_DETUNE_451_NM
        and row.mode451.q_value >= MIN_Q_451
    )


def safe_model_name(raw: str) -> str:
    cleaned = "".join(ch if ch.isalnum() else "_" for ch in raw)
    return cleaned[:50]


def periodic_designs() -> list[Design]:
    designs: dict[str, Design] = {}

    def add(design: Design) -> None:
        if is_fabricable(design):
            designs[design.design_id] = design

    add(Design("periodic_N1_baseline", "periodic_supercell", (P0_NM,), (W0_NM,)))

    for period in [562.3, 563.3, 564.0, 564.2, 564.25, 564.35, 564.4, 565.3, 566.3]:
        for width_offset in [0.0, -5.0, 5.0]:
            add(
                Design(
                    f"primitive_P{period:.2f}_Woff{width_offset:+.0f}",
                    "primitive_tune",
                    (period,),
                    (0.68 * period + width_offset,),
                )
            )

    for d_p in [0.0, 40.0, 80.0]:
        for d_w in [0.0, 40.0, 60.0]:
            add(
                Design(
                    f"periodic_N2_dP{d_p:.0f}_dW{d_w:.0f}",
                    "periodic_supercell",
                    (P0_NM - d_p / 2.0, P0_NM + d_p / 2.0),
                    (W0_NM - d_w / 2.0, W0_NM + d_w / 2.0),
                )
            )
            add(
                Design(
                    f"periodic_N3_dP{d_p:.0f}_dW{d_w:.0f}",
                    "periodic_supercell",
                    (P0_NM - d_p / 2.0, P0_NM + d_p, P0_NM - d_p / 2.0),
                    (W0_NM - d_w / 2.0, W0_NM + d_w, W0_NM - d_w / 2.0),
                )
            )

    for d_p in [20.0, 40.0, 60.0]:
        for d_w in [20.0, 40.0, 60.0]:
            add(
                Design(
                    f"periodic_N5_dP{d_p:.0f}_dW{d_w:.0f}",
                    "periodic_supercell",
                    (P0_NM - d_p, P0_NM, P0_NM + 2.0 * d_p, P0_NM, P0_NM - d_p),
                    (W0_NM - d_w, W0_NM, W0_NM + 2.0 * d_w, W0_NM, W0_NM - d_w),
                )
            )

    return list(designs.values())


def wide_top_groove_designs() -> list[Design]:
    designs: dict[str, Design] = {}

    def add(design: Design) -> None:
        if is_fabricable(design):
            designs[design.design_id] = design

    for groove_w in [140.0, 180.0, 220.0]:
        for groove_d in [100.0, 140.0, 180.0]:
            add(
                Design(
                    f"wide_top_groove_w{groove_w:.0f}_d{groove_d:.0f}",
                    "wide_top_groove",
                    (P0_NM,),
                    (W0_NM,),
                    groove_w,
                    groove_d,
                )
            )
            add(
                Design(
                    f"wide_N2_w{groove_w:.0f}_d{groove_d:.0f}_dW40",
                    "wide_top_groove",
                    (P0_NM, P0_NM),
                    (W0_NM - 20.0, W0_NM + 20.0),
                    groove_w,
                    groove_d,
                )
            )
    return list(designs.values())


def write_dual_csv(path: Path, rows: list[DualResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "design_id",
                "family",
                "n_periods",
                "pitches_nm",
                "widths_nm",
                "groove_w_nm",
                "groove_d_nm",
                "min_feature_nm",
                "min_gap_nm",
                "min_remaining_gap_nm",
                "lambda980_nm",
                "q980",
                "mode980",
                "lambda451_nm",
                "q451",
                "mode451",
                "re_overlap451",
                "folded_only",
                "dual_score",
                "passes_dual_targets",
            ]
        )
        for row in rows:
            d = row.design
            writer.writerow(
                [
                    d.design_id,
                    d.family,
                    d.n_periods,
                    ";".join(f"{p:.9g}" for p in d.pitches_nm),
                    ";".join(f"{w:.9g}" for w in d.widths_nm),
                    d.groove_w_nm,
                    d.groove_d_nm,
                    row.min_feature_nm,
                    row.min_gap_nm,
                    row.min_remaining_gap_nm,
                    row.mode980.wavelength_nm,
                    row.mode980.q_value,
                    row.mode980.mode_index,
                    row.mode451.wavelength_nm,
                    row.mode451.q_value,
                    row.mode451.mode_index,
                    row.mode451.re_overlap,
                    is_folded_only(d),
                    row.dual_score,
                    passes_dual_targets(row),
                ]
            )


def write_mode_csv(path: Path, rows: list[ModeResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "design_id",
                "family",
                "target_nm",
                "n_periods",
                "kx_norm",
                "mode_index",
                "freq_thz_real",
                "freq_thz_imag",
                "wavelength_nm",
                "q_value",
                "re_overlap",
                "score",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r.design_id,
                    r.family,
                    r.target_nm,
                    r.n_periods,
                    r.kx_norm,
                    r.mode_index,
                    r.freq_thz_real,
                    r.freq_thz_imag,
                    r.wavelength_nm,
                    r.q_value,
                    r.re_overlap,
                    r.score,
                ]
            )


def read_dual_csv(path: Path) -> list[DualResult]:
    rows: list[DualResult] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            design = Design(
                raw["design_id"],
                raw["family"],
                tuple(float(x) for x in raw["pitches_nm"].split(";") if x),
                tuple(float(x) for x in raw["widths_nm"].split(";") if x),
                float(raw["groove_w_nm"]),
                float(raw["groove_d_nm"]),
            )
            mode980 = ModeResult(
                design.design_id,
                design.family,
                TARGET_980_NM,
                design.n_periods,
                0.0,
                int(raw["mode980"]),
                math.nan,
                math.nan,
                float(raw["lambda980_nm"]),
                float(raw["q980"]),
                math.nan,
                math.nan,
            )
            mode451 = ModeResult(
                design.design_id,
                design.family,
                TARGET_451_NM,
                design.n_periods,
                0.0,
                int(raw["mode451"]),
                math.nan,
                math.nan,
                float(raw["lambda451_nm"]),
                float(raw["q451"]),
                float(raw["re_overlap451"]) if raw["re_overlap451"] not in ("", "nan") else math.nan,
                math.nan,
            )
            rows.append(
                DualResult(
                    design,
                    mode980,
                    mode451,
                    float(raw["min_feature_nm"]),
                    float(raw["min_gap_nm"]),
                    float(raw["min_remaining_gap_nm"]),
                    float(raw["dual_score"]),
                )
            )
    return rows


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


def svg_line_plot(path: Path, title: str, xs: list[float], series: list[tuple[str, list[float]]], y_label: str, x_label: str, log_y: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 920, 540
    ml, mr, mt, mb = 95, 35, 58, 78
    plot_w, plot_h = width - ml - mr, height - mt - mb
    valid = [(x, y) for _, ys in series for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y) and (not log_y or y > 0)]
    if not valid:
        path.write_text("<svg xmlns='http://www.w3.org/2000/svg'></svg>", encoding="utf-8")
        return
    x_min, x_max = min(x for x, _ in valid), max(x for x, _ in valid)
    y_values = [math.log10(y) if log_y else y for _, y in valid]
    y_min, y_max = min(y_values), max(y_values)
    if abs(x_max - x_min) < 1e-12:
        x_max += 1.0
    if abs(y_max - y_min) < 1e-12:
        y_max += 1.0

    def sx(x: float) -> float:
        return ml + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        yy = math.log10(max(y, 1e-300)) if log_y else y
        return mt + (y_max - yy) / (y_max - y_min) * plot_h

    colors = ["#0b5fff", "#d22f27", "#14853b", "#7c3aed"]
    out = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2}' y='32' text-anchor='middle' font-size='22' font-family='Arial'>{title}</text>",
        f"<line x1='{ml}' y1='{mt+plot_h}' x2='{ml+plot_w}' y2='{mt+plot_h}' stroke='black'/>",
        f"<line x1='{ml}' y1='{mt}' x2='{ml}' y2='{mt+plot_h}' stroke='black'/>",
        f"<text x='{width/2}' y='{height-26}' text-anchor='middle' font-size='16' font-family='Arial'>{x_label}</text>",
        f"<text x='25' y='{height/2}' text-anchor='middle' transform='rotate(-90 25 {height/2})' font-size='16' font-family='Arial'>{y_label}</text>",
    ]
    for i in range(6):
        tx = x_min + i * (x_max - x_min) / 5.0
        out.append(f"<text x='{sx(tx):.2f}' y='{mt+plot_h+24}' text-anchor='middle' font-size='12' font-family='Arial'>{tx:.3g}</text>")
    for i in range(6):
        yy = y_min + i * (y_max - y_min) / 5.0
        label = f"1e{yy:.0f}" if log_y else f"{yy:.4g}"
        py = mt + (5 - i) / 5.0 * plot_h
        out.append(f"<text x='{ml-9}' y='{py+4:.2f}' text-anchor='end' font-size='12' font-family='Arial'>{label}</text>")
    for si, (label, ys) in enumerate(series):
        pts = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(xs, ys) if math.isfinite(y) and (not log_y or y > 0))
        color = colors[si % len(colors)]
        out.append(f"<polyline points='{pts}' fill='none' stroke='{color}' stroke-width='2.4'/>")
        out.append(f"<text x='{ml+12}' y='{mt+20+si*21}' font-size='14' font-family='Arial' fill='{color}'>{label}</text>")
    out.append("</svg>")
    path.write_text("\n".join(out), encoding="utf-8")


def svg_scatter(path: Path, title: str, rows: list[DualResult], y_getter, y_label: str, log_y: bool = False) -> None:
    xs = list(range(1, len(rows) + 1))
    ys = [float(y_getter(r)) for r in rows]
    svg_line_plot(path, title, xs, [("candidate", ys)], y_label, "candidate rank", log_y)


def point_in_rare_earth(design: Design, x: float, y: float) -> bool:
    y_sio2 = H_PML_NM + H_SI_BUF_NM
    y_sin = y_sio2 + H_SIO2_NM
    y_sin_top = y_sin + H_SIN_NM
    ridges = ridge_edges(design)
    for ridge in ridges:
        center = 0.5 * (ridge.left_nm + ridge.right_nm)
        if ridge.left_nm <= x <= ridge.right_nm and y_sin_top <= y <= y_sin_top + H_RE_NM:
            return True
        if ridge.left_nm - H_RE_NM <= x <= ridge.left_nm and y_sin <= y <= y_sin_top:
            return True
        if ridge.right_nm <= x <= ridge.right_nm + H_RE_NM and y_sin <= y <= y_sin_top:
            return True
        if design.has_top_groove:
            gl = center - design.groove_w_nm / 2.0
            gr = center + design.groove_w_nm / 2.0
            gb = y_sin_top - design.groove_d_nm
            if gl + H_RE_NM <= x <= gr - H_RE_NM and gb <= y <= gb + H_RE_NM:
                return True
            if gl <= x <= gl + H_RE_NM and gb <= y <= y_sin_top:
                return True
            if gr - H_RE_NM <= x <= gr and gb <= y <= y_sin_top:
                return True
    for a, b in zip(ridges, ridges[1:]):
        if a.right_nm + H_RE_NM <= x <= b.left_nm - H_RE_NM and y_sin <= y <= y_sin + H_RE_NM:
            return True
    if -design.length_nm / 2.0 <= x <= ridges[0].left_nm - H_RE_NM and y_sin <= y <= y_sin + H_RE_NM:
        return True
    if ridges[-1].right_nm + H_RE_NM <= x <= design.length_nm / 2.0 and y_sin <= y <= y_sin + H_RE_NM:
        return True
    return False


def compute_re_overlap_from_solution(model, design: Design, mode: ModeResult) -> float:
    import numpy as np

    freqs = [complex(v) for v in np.asarray(model.evaluate("freq", unit="THz")).ravel()]
    target_idx = find_solution_index_for_mode(freqs, mode)
    ds = model.datasets()[0] if model.datasets() else None
    x = np.asarray(model.evaluate("x", unit="nm", dataset=ds, inner=[target_idx])).ravel().astype(float)
    y = np.asarray(model.evaluate("y", unit="nm", dataset=ds, inner=[target_idx])).ravel().astype(float)
    ez = np.asarray(model.evaluate("abs(ewfd.Ez)", dataset=ds, inner=[target_idx])).ravel().astype(float)
    n = min(len(x), len(y), len(ez))
    x, y, e2 = x[:n], y[:n], ez[:n] * ez[:n]
    y_focus_min = H_PML_NM + H_SI_BUF_NM + H_SIO2_NM - 200.0
    y_focus_max = H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_SIN_NM + 300.0
    total = 0.0
    re_sum = 0.0
    for xx, yy, val in zip(x, y, e2):
        if y_focus_min <= yy <= y_focus_max:
            total += float(val)
            if point_in_rare_earth(design, float(xx), float(yy)):
                re_sum += float(val)
    return re_sum / total if total > 0 else math.nan


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


def export_field(model, design: Design, mode: ModeResult, prefix: str) -> None:
    import numpy as np

    freqs = [complex(v) for v in np.asarray(model.evaluate("freq", unit="THz")).ravel()]
    target_idx = find_solution_index_for_mode(freqs, mode)
    ds = model.datasets()[0] if model.datasets() else None
    x = np.asarray(model.evaluate("x", unit="nm", dataset=ds, inner=[target_idx])).ravel().astype(float)
    y = np.asarray(model.evaluate("y", unit="nm", dataset=ds, inner=[target_idx])).ravel().astype(float)
    ez = np.asarray(model.evaluate("abs(ewfd.Ez)", dataset=ds, inner=[target_idx])).ravel().astype(float)
    n = min(len(x), len(y), len(ez))
    x, y, e2 = x[:n], y[:n], ez[:n] * ez[:n]
    csv_path = OUT / f"{prefix}_field_ez2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x_nm", "y_nm", "abs_ez_squared", "in_rare_earth_mask"])
        for xx, yy, val in zip(x, y, e2):
            writer.writerow([xx, yy, val, point_in_rare_earth(design, float(xx), float(yy))])
    svg_field_plot(OUT / f"{prefix}_field_ez2.svg", f"{prefix} |Ez|^2, {design.design_id}", design, x.tolist(), y.tolist(), e2.tolist())


def svg_field_plot(path: Path, title: str, design: Design, x_nm: list[float], y_nm: list[float], value: list[float]) -> None:
    width, height = 980, 650
    ml, mr, mt, mb = 82, 38, 54, 74
    plot_w, plot_h = width - ml - mr, height - mt - mb
    y_sin = H_PML_NM + H_SI_BUF_NM + H_SIO2_NM
    y_min_focus = y_sin - 220.0
    y_max_focus = y_sin + H_SIN_NM + 420.0
    points = [(x, y, v) for x, y, v in zip(x_nm, y_nm, value) if -design.length_nm / 2 <= x <= design.length_nm / 2 and y_min_focus <= y <= y_max_focus]
    if not points:
        points = list(zip(x_nm, y_nm, value))
    if len(points) > 6500:
        step = max(1, len(points) // 6500)
        points = points[::step]
    x_min, x_max = -design.length_nm / 2.0, design.length_nm / 2.0
    y_min, y_max = min(y for _, y, _ in points), max(y for _, y, _ in points)
    vmax = max(max(v for _, _, v in points), 1e-300)

    def sx(x: float) -> float:
        return ml + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        return mt + (y_max - y) / (y_max - y_min) * plot_h

    def normalized(v: float) -> float:
        return max(0.0, min(1.0, (math.log10(v / vmax + 1e-12) + 12.0) / 12.0))

    out = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2}' y='30' text-anchor='middle' font-size='21' font-family='Arial'>{title}</text>",
        f"<rect x='{ml}' y='{mt}' width='{plot_w}' height='{plot_h}' fill='#05051f' stroke='black'/>",
    ]
    for x, y, v in points:
        out.append(f"<circle cx='{sx(x):.2f}' cy='{sy(y):.2f}' r='1.6' fill='{plasma_color(normalized(v))}' opacity='0.9'/>")
    for ridge in ridge_edges(design):
        out.append(f"<rect x='{sx(ridge.left_nm):.2f}' y='{sy(y_sin+H_SIN_NM):.2f}' width='{sx(ridge.right_nm)-sx(ridge.left_nm):.2f}' height='{sy(y_sin)-sy(y_sin+H_SIN_NM):.2f}' fill='none' stroke='white' stroke-width='1.8'/>")
    for yy, label in [(y_sin, "SiO2 / SiN"), (y_sin + H_SIN_NM, "SiN / air")]:
        if y_min <= yy <= y_max:
            out.append(f"<line x1='{ml}' y1='{sy(yy):.2f}' x2='{ml+plot_w}' y2='{sy(yy):.2f}' stroke='white' stroke-width='1.2' stroke-dasharray='6 4'/>")
            out.append(f"<text x='{ml+8}' y='{sy(yy)-6:.2f}' font-size='12' font-family='Arial' fill='white'>{label}</text>")
    out.append(f"<text x='{ml+10}' y='{mt+20}' font-size='12' font-family='Arial' fill='white'>white: SiN outline; RE mask is included in CSV</text>")
    out.extend(
        [
            f"<text x='{width/2}' y='{height-26}' text-anchor='middle' font-size='16' font-family='Arial'>x (nm)</text>",
            f"<text x='25' y='{height/2}' text-anchor='middle' transform='rotate(-90 25 {height/2})' font-size='16' font-family='Arial'>vertical coordinate y (nm)</text>",
            "</svg>",
        ]
    )
    path.write_text("\n".join(out), encoding="utf-8")


def solve_design(client, design: Design, compute_overlap: bool = False) -> DualResult:
    model_name = safe_model_name(f"v4_{design.design_id}")
    model = build_model(client, design, model_name)
    try:
        r980 = best_mode_for_target(solve_target(model, design, TARGET_980_NM, 0.0), TARGET_980_NM)
        r451 = best_mode_for_target(solve_target(model, design, TARGET_451_NM, 0.0), TARGET_451_NM)
        if compute_overlap:
            overlap = compute_re_overlap_from_solution(model, design, r451)
            r451.re_overlap = overlap
            r451.score = mode_result(design, TARGET_451_NM, 0.0, r451.mode_index, complex(r451.freq_thz_real, r451.freq_thz_imag), overlap).score
        return DualResult(
            design,
            r980,
            r451,
            min_feature_nm(design),
            min(gap_widths(design)),
            min_remaining_gap_nm(design),
            dual_score(r980, r451),
        )
    finally:
        client.remove(model.name())


def run_scan_family(client, designs: list[Design], csv_path: Path, label: str, compute_overlap_top: int = 5) -> list[DualResult]:
    rows: list[DualResult] = []
    for i, design in enumerate(designs, 1):
        print(f"[{label}] {i}/{len(designs)} solving {design.design_id}", flush=True)
        try:
            result = solve_design(client, design, compute_overlap=False)
            rows.append(result)
            write_dual_csv(csv_path, rows)
            print(
                f"  980: {result.mode980.wavelength_nm:.4f} nm Q={result.mode980.q_value:.3g}; "
                f"451: {result.mode451.wavelength_nm:.4f} nm Q={result.mode451.q_value:.3g}; score={result.dual_score:.4g}",
                flush=True,
            )
        except Exception as exc:
            print(f"  FAILED {design.design_id}: {exc}", flush=True)
    ranked = sorted(rows, key=lambda r: (0 if passes_dual_targets(r) else 1, 1 if is_folded_only(r.design) else 0, r.dual_score))
    enriched: list[DualResult] = []
    for result in ranked[: min(compute_overlap_top, len(ranked))]:
        print(f"[{label}] computing 451 RE overlap for {result.design.design_id}", flush=True)
        try:
            enriched_result = solve_design(client, result.design, compute_overlap=True)
            enriched.append(enriched_result)
            for idx, old in enumerate(rows):
                if old.design.design_id == result.design.design_id:
                    rows[idx] = enriched_result
                    break
            write_dual_csv(csv_path, rows)
        except Exception as exc:
            print(f"  overlap failed for {result.design.design_id}: {exc}", flush=True)
    return rows


def save_best_models(client, best: DualResult, wide_compare: DualResult | None = None) -> None:
    model = build_model(client, best.design, "sin980_451_dual_best")
    try:
        r980 = best_mode_for_target(solve_target(model, best.design, TARGET_980_NM, 0.0), TARGET_980_NM)
        model.save(str(OUT / "sin980_451_dual_best.mph"))
        export_field(model, best.design, r980, "field_980")
        r451 = best_mode_for_target(solve_target(model, best.design, TARGET_451_NM, 0.0), TARGET_451_NM)
        model.save(str(OUT / "sin980_451_dual_best_451check.mph"))
        export_field(model, best.design, r451, "field_451")
    finally:
        client.remove(model.name())

    if wide_compare is not None:
        model = build_model(client, wide_compare.design, "sin980_451_dual_wide_top_groove_compare")
        try:
            solve_target(model, wide_compare.design, TARGET_980_NM, 0.0)
            model.save(str(OUT / "sin980_451_dual_wide_top_groove_compare.mph"))
        finally:
            client.remove(model.name())


def run_band_for_best(client, best: DualResult) -> None:
    k_values = [0.0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0]
    model = build_model(client, best.design, "sin980_451_dual_band")
    rows980: list[ModeResult] = []
    rows451: list[ModeResult] = []
    try:
        for kx in k_values:
            print(f"[band] kx={kx:.3g} pi/L at 980 nm", flush=True)
            rows = solve_target(model, best.design, TARGET_980_NM, kx)
            rows980.append(best_mode_for_target(rows, TARGET_980_NM))
        for kx in k_values:
            print(f"[band] kx={kx:.3g} pi/L at 451 nm", flush=True)
            rows = solve_target(model, best.design, TARGET_451_NM, kx)
            rows451.append(best_mode_for_target(rows, TARGET_451_NM))
        write_mode_csv(OUT / "band_980.csv", rows980)
        write_mode_csv(OUT / "band_451.csv", rows451)
        svg_line_plot(OUT / "band_980_wavelength.svg", "980 nm folded band", [r.kx_norm for r in rows980], [("980 branch", [r.wavelength_nm for r in rows980])], "wavelength (nm)", "kx / (pi/L)")
        svg_line_plot(OUT / "band_980_q.svg", "980 nm Q along folded band", [r.kx_norm for r in rows980], [("980 branch", [r.q_value for r in rows980])], "Q", "kx / (pi/L)", log_y=True)
        svg_line_plot(OUT / "band_451_wavelength.svg", "451 nm folded band", [r.kx_norm for r in rows451], [("451 branch", [r.wavelength_nm for r in rows451])], "wavelength (nm)", "kx / (pi/L)")
        svg_line_plot(OUT / "band_451_q.svg", "451 nm Q along folded band", [r.kx_norm for r in rows451], [("451 branch", [r.q_value for r in rows451])], "Q", "kx / (pi/L)", log_y=True)
    finally:
        client.remove(model.name())


def plot_scan_outputs(periodic_rows: list[DualResult], wide_rows: list[DualResult]) -> None:
    ranked = sorted(periodic_rows + wide_rows, key=lambda r: (0 if passes_dual_targets(r) else 1, r.dual_score))
    if not ranked:
        return
    top = ranked[: min(30, len(ranked))]
    svg_scatter(OUT / "scan_lambda980.svg", "Top candidate 980 nm wavelength", top, lambda r: r.mode980.wavelength_nm, "lambda980 (nm)")
    svg_scatter(OUT / "scan_q980.svg", "Top candidate 980 nm Q", top, lambda r: r.mode980.q_value, "Q980", log_y=True)
    svg_scatter(OUT / "scan_lambda451.svg", "Top candidate 451 nm wavelength", top, lambda r: r.mode451.wavelength_nm, "lambda451 (nm)")
    svg_scatter(OUT / "scan_q451.svg", "Top candidate 451 nm Q", top, lambda r: r.mode451.q_value, "Q451", log_y=True)
    overlap_rows = [r for r in top if math.isfinite(r.mode451.re_overlap)]
    if overlap_rows:
        svg_scatter(OUT / "re_overlap451.svg", "451 nm rare-earth layer field-overlap proxy", overlap_rows, lambda r: r.mode451.re_overlap, "RE overlap proxy")


def write_report(periodic_rows: list[DualResult], wide_rows: list[DualResult], best: DualResult | None, wide_compare: DualResult | None, wide_was_run: bool) -> None:
    folded_rows = [r for r in periodic_rows + wide_rows if is_folded_only(r.design)]
    true_rows = [r for r in periodic_rows + wide_rows if not is_folded_only(r.design)]
    true_passing = [r for r in true_rows if passes_dual_targets(r)]
    lines: list[str] = [
        "# v4 980/451 nm dual-resonance quasi-BIC design summary",
        "",
        "All v4 files were generated inside `E:\\xitugrating\\sin980_qbic_2d_design_v4_dual_451_980`; v1/v2/v3 and the original inputs were not edited.",
        "",
        "## Targets",
        "",
        f"- 980 nm pump quasi-BIC: `|lambda-980| <= {MAX_DETUNE_980_NM} nm`, `Q >= {MIN_Q_980:.3g}`.",
        f"- 451 nm emission out-coupling resonance: `|lambda-451| <= {MAX_DETUNE_451_NM} nm`, `Q >= {MIN_Q_451:.3g}`, preferably `1e3-1e5`.",
        f"- Fabrication filter: minimum feature `>= {MIN_FEATURE_NM:.0f} nm`; remaining air gap after 16 nm side RE layers `>= {MIN_REMAINING_GAP_NM:.0f} nm`.",
        "",
        "## Scan status",
        "",
        f"- Periodic/supercell candidates simulated: `{len(periodic_rows)}`.",
        f"- Wide top-groove backup candidates simulated: `{len(wide_rows)}`." if wide_was_run else "- Wide top-groove backup scan was skipped because a periodic/supercell candidate was available for recommendation.",
        f"- Folded-only supercell candidates excluded from recommendation: `{len(folded_rows)}`.",
        f"- True geometry candidates satisfying both thresholds: `{len(true_passing)}`.",
        "",
    ]
    if best is None:
        lines.extend(
            [
                "## Result",
                "",
                "No candidate satisfied the dual-resonance thresholds under the fabrication filters.",
                "",
                "Recommended next relaxations are: allow a finite nonperiodic grating optimization, broaden the 451 nm detuning window, or accept a 451 nm mode with lower out-coupling Q as a proof of existence.",
            ]
        )
    else:
        d = best.design
        lines.extend(
            [
                "## Recommended candidate",
                "",
                "| item | value |",
                "|---|---:|",
                f"| design id | `{d.design_id}` |",
                f"| family | `{d.family}` |",
                f"| supercell periods | `{d.n_periods}` |",
                f"| pitches (nm) | `{'; '.join(f'{p:.3f}' for p in d.pitches_nm)}` |",
                f"| ridge widths (nm) | `{'; '.join(f'{w:.3f}' for w in d.widths_nm)}` |",
                f"| top groove width/depth (nm) | `{d.groove_w_nm:.3f} / {d.groove_d_nm:.3f}` |",
                f"| min feature (nm) | `{best.min_feature_nm:.3f}` |",
                f"| min etched gap (nm) | `{best.min_gap_nm:.3f}` |",
                f"| min remaining air gap after RE (nm) | `{best.min_remaining_gap_nm:.3f}` |",
                f"| lambda980 (nm) | `{best.mode980.wavelength_nm:.9g}` |",
                f"| Q980 | `{best.mode980.q_value:.6g}` |",
                f"| lambda451 (nm) | `{best.mode451.wavelength_nm:.9g}` |",
                f"| Q451 | `{best.mode451.q_value:.6g}` |",
                f"| 451 RE overlap proxy | `{best.mode451.re_overlap:.6g}` |",
                "",
                "Saved models:",
                "",
                "- `outputs/sin980_451_dual_best.mph` uses the recommended geometry solved at 980 nm.",
                "- `outputs/sin980_451_dual_best_451check.mph` uses the same geometry solved at 451 nm for field inspection.",
            ]
        )
        if not passes_dual_targets(best):
            lines.extend(
                [
                    "",
                    "Note: this is the best fabrication-filtered true-geometry proof-of-existence candidate found in the current scan, but it does not fully satisfy every target threshold.",
                ]
            )
    if folded_rows:
        folded_best = sorted(folded_rows, key=lambda r: r.dual_score)[0]
        lines.extend(
            [
                "",
                "## Folded-band evidence, not recommended as a standalone geometry",
                "",
                f"- Best folded-only case: `{folded_best.design.design_id}`.",
                f"- It gives 980 nm `{folded_best.mode980.wavelength_nm:.6g} nm`, Q `{folded_best.mode980.q_value:.6g}` and 451 nm `{folded_best.mode451.wavelength_nm:.6g} nm`, Q `{folded_best.mode451.q_value:.6g}`.",
                "- Because its supercell does not introduce a real geometric modulation, it is treated as evidence that a 451 nm branch exists near the desired wavelength, not as the final dual-resonance structure.",
            ]
        )
    if wide_compare is not None:
        lines.extend(
            [
                "",
                "## Wide top-groove comparison",
                "",
                f"- Best backup: `{wide_compare.design.design_id}`.",
                f"- 980 nm: `{wide_compare.mode980.wavelength_nm:.6g} nm`, Q `{wide_compare.mode980.q_value:.6g}`.",
                f"- 451 nm: `{wide_compare.mode451.wavelength_nm:.6g} nm`, Q `{wide_compare.mode451.q_value:.6g}`.",
                "- This backup is not preferred unless it clearly beats the periodic/supercell route because rare-earth coverage in a top groove is less reliable than in the large primary gaps.",
            ]
        )
    lines.extend(
        [
            "",
            "## Output files",
            "",
            "- `periodic_supercell_scan.csv`",
            "- `wide_top_groove_scan.csv`",
            "- `dual_best_candidates.csv`",
            "- `band_980.csv`, `band_451.csv`",
            "- `field_980_field_ez2.svg`, `field_451_field_ez2.svg`",
            "- `scan_q980.svg`, `scan_q451.svg`, `scan_lambda980.svg`, `scan_lambda451.svg`",
        ]
    )
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def run_all(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = import_mph()
    client = mph.Client(cores=args.cores)
    periodic_rows: list[DualResult] = []
    wide_rows: list[DualResult] = []
    wide_was_run = False
    try:
        periodic = periodic_designs()
        if args.limit_periodic:
            periodic = periodic[: args.limit_periodic]
        periodic_rows = run_scan_family(client, periodic, OUT / "periodic_supercell_scan.csv", "periodic", compute_overlap_top=args.overlap_top)
        passing_periodic = [r for r in periodic_rows if passes_dual_targets(r)]
        if args.force_wide or not passing_periodic:
            wide_was_run = True
            wide = wide_top_groove_designs()
            if args.limit_wide:
                wide = wide[: args.limit_wide]
            wide_rows = run_scan_family(client, wide, OUT / "wide_top_groove_scan.csv", "wide", compute_overlap_top=args.overlap_top)
        else:
            write_dual_csv(OUT / "wide_top_groove_scan.csv", [])

        all_rows = periodic_rows + wide_rows
        ranked = sorted(all_rows, key=lambda r: (0 if passes_dual_targets(r) else 1, 1 if is_folded_only(r.design) else 0, r.dual_score))
        best = ranked[0] if ranked else None
        write_dual_csv(OUT / "dual_best_candidates.csv", ranked[: min(20, len(ranked))])
        wide_candidates = sorted(wide_rows, key=lambda r: (0 if passes_dual_targets(r) else 1, 1 if is_folded_only(r.design) else 0, r.dual_score))
        wide_compare = wide_candidates[0] if wide_candidates else None
        plot_scan_outputs(periodic_rows, wide_rows)
        if best is not None:
            save_best_models(client, best, wide_compare)
            run_band_for_best(client, best)
        write_report(periodic_rows, wide_rows, best, wide_compare, wide_was_run)
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="v4 dual-resonance 980/451 nm SiN quasi-BIC design scan.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_all = sub.add_parser("run-all", help="Run the isolated v4 scan, save models, plots, and report.")
    p_all.add_argument("--cores", type=int, default=2)
    p_all.add_argument("--limit-periodic", type=int, default=0, help="Debug limit for periodic/supercell candidates; 0 means all.")
    p_all.add_argument("--limit-wide", type=int, default=0, help="Debug limit for wide-groove candidates; 0 means all.")
    p_all.add_argument("--force-wide", action="store_true", help="Always run the wide top-groove backup scan.")
    p_all.add_argument("--overlap-top", type=int, default=5, help="Compute 451 nm RE overlap for the top N candidates in each family.")
    p_all.set_defaults(func=run_all)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
