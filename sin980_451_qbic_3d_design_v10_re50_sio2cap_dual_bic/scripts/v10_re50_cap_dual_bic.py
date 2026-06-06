from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import shutil
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
Q_TARGET_980 = 1e8
Q_TARGET_451 = 1e5

N_SIN_980 = 2.069946018
N_SIN_451 = 2.128247126
N_SIO2_980 = 1.45
N_SIO2_451 = 1.466
N_RE_980 = 1.355
N_RE_451 = 1.365
N_SI_980 = 3.55
N_SI_451 = 4.7

H_SI_BUF_NM = 250.0
H_SIO2_NM = 1500.0
H_RE_NM = 50.0
H_CAP_NM = 12.5
H_SIN_NM = 300.0
H_AIR_NM = 900.0
H_PML_NM = 500.0

MIN_FEATURE_NM = 100.0
PREFERRED_GAP_NM = 120.0
MESH_AUTO_SIZE = 7
NEIGS = 10


@dataclass(frozen=True)
class CircleHole:
    x_nm: float
    y_nm: float
    radius_nm: float


@dataclass(frozen=True)
class RectHole:
    x_nm: float
    y_nm: float
    width_nm: float
    depth_nm: float


@dataclass(frozen=True)
class CapsuleHole:
    x_nm: float
    y_nm: float
    length_nm: float
    width_nm: float
    axis: str = "x"


@dataclass(frozen=True)
class Pillar:
    x_nm: float
    y_nm: float
    radius_nm: float


@dataclass(frozen=True)
class Design:
    design_id: str
    family: str
    px_nm: float
    py_nm: float
    circle_holes: tuple[CircleHole, ...] = ()
    rect_holes: tuple[RectHole, ...] = ()
    capsule_holes: tuple[CapsuleHole, ...] = ()
    pillars: tuple[Pillar, ...] = ()
    notes: str = ""

    @property
    def is_pillar(self) -> bool:
        return bool(self.pillars)

    @property
    def feature_min_nm(self) -> float:
        values: list[float] = [self.px_nm, self.py_nm]
        values.extend(2.0 * h.radius_nm for h in self.circle_holes)
        values.extend(min(h.width_nm, h.depth_nm) for h in self.rect_holes)
        values.extend(min(h.length_nm, h.width_nm) for h in self.capsule_holes)
        values.extend(2.0 * p.radius_nm for p in self.pillars)
        values.extend(self._edge_and_pair_clearances())
        return min(values) if values else min(self.px_nm, self.py_nm)

    @property
    def preferred_gap_ok(self) -> bool:
        return self.feature_min_nm >= PREFERRED_GAP_NM

    @property
    def fabricable(self) -> bool:
        return self.feature_min_nm >= MIN_FEATURE_NM

    def _edge_and_pair_clearances(self) -> list[float]:
        clearances: list[float] = []
        round_objects: list[tuple[float, float, float, str]] = []
        for h in self.circle_holes:
            round_objects.append((h.x_nm, h.y_nm, h.radius_nm, "hole"))
        for h in self.rect_holes:
            if h.width_nm < self.px_nm:
                clearances.extend([self.px_nm / 2.0 - abs(h.x_nm) - h.width_nm / 2.0])
            if h.depth_nm < self.py_nm:
                clearances.extend([self.py_nm / 2.0 - abs(h.y_nm) - h.depth_nm / 2.0])
        for h in self.capsule_holes:
            rx = h.length_nm / 2.0 if h.axis.lower() == "x" else h.width_nm / 2.0
            ry = h.width_nm / 2.0 if h.axis.lower() == "x" else h.length_nm / 2.0
            clearances.extend([
                self.px_nm / 2.0 - abs(h.x_nm) - rx,
                self.py_nm / 2.0 - abs(h.y_nm) - ry,
            ])
        for p in self.pillars:
            round_objects.append((p.x_nm, p.y_nm, p.radius_nm, "pillar"))

        # Edge clearance is a proxy for the remaining SiN bridge to periodic images.
        for x, y, r, _kind in round_objects:
            clearances.extend([
                self.px_nm / 2.0 - abs(x) - r,
                self.py_nm / 2.0 - abs(y) - r,
            ])

        for i, a in enumerate(round_objects):
            for b in round_objects[i + 1 :]:
                dx = abs(a[0] - b[0])
                dy = abs(a[1] - b[1])
                # Include nearest periodic image distance for a compact fabrication proxy.
                dx = min(dx, self.px_nm - dx)
                dy = min(dy, self.py_nm - dy)
                clearances.append(math.hypot(dx, dy) - a[2] - b[2])
        for i, a in enumerate(self.rect_holes):
            for b in self.rect_holes[i + 1 :]:
                dx = min(abs(a.x_nm - b.x_nm), self.px_nm - abs(a.x_nm - b.x_nm))
                dy = min(abs(a.y_nm - b.y_nm), self.py_nm - abs(a.y_nm - b.y_nm))
                gap_x = dx - (a.width_nm + b.width_nm) / 2.0
                gap_y = dy - (a.depth_nm + b.depth_nm) / 2.0
                if gap_x <= 0 and gap_y <= 0:
                    # Overlapping rectangles form one connected air opening, not a narrow bridge.
                    continue
                if gap_x <= 0:
                    clearances.append(gap_y)
                elif gap_y <= 0:
                    clearances.append(gap_x)
                else:
                    clearances.append(math.hypot(gap_x, gap_y))
        return [v for v in clearances if math.isfinite(v)]


@dataclass
class ModeResult:
    design_id: str
    family: str
    target_nm: float
    px_nm: float
    py_nm: float
    feature_min_nm: float
    fabricable: bool
    preferred_gap_ok: bool
    kx_norm: float
    ky_norm: float
    mode_index: int
    freq_thz_real: float
    freq_thz_imag: float
    wavelength_nm: float
    q_value: float
    score: float


@dataclass
class DualCandidate:
    design: Design
    mode980: ModeResult
    mode451: ModeResult
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


def stable_seed(text: str) -> int:
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16)


def fmt(value: float | str) -> str:
    if isinstance(value, str):
        return value
    if abs(value) < 1e-12:
        value = 0.0
    return f"{value:.12g}"


def target_indices(target_nm: float) -> tuple[float, float, float, float]:
    if abs(target_nm - TARGET_451_NM) < 10.0:
        return N_SIN_451, N_SIO2_451, N_RE_451, N_SI_451
    return N_SIN_980, N_SIO2_980, N_RE_980, N_SI_980


def set_param(model, name: str, value: str, description: str = "") -> None:
    if description:
        model.java.param().set(name, value, description)
    else:
        model.java.param().set(name, value)


def create_block(geom, tag: str, pos: tuple[float, float, float], size: tuple[float, float, float], label: str) -> None:
    blk = geom.create(tag, "Block")
    blk.label(f"{tag}: {label}")
    blk.set("base", "corner")
    blk.set("pos", [fmt(v) for v in pos])
    blk.set("size", [fmt(v) for v in size])
    blk.set("selresult", "on")
    blk.set("selresultshow", "all")


def create_cylinder(geom, tag: str, x: float, y: float, z: float, radius: float, height: float, label: str) -> None:
    cyl = geom.create(tag, "Cylinder")
    cyl.label(f"{tag}: {label}")
    cyl.set("r", fmt(radius))
    cyl.set("h", fmt(height))
    cyl.set("pos", [fmt(x), fmt(y), fmt(z)])
    cyl.set("selresult", "on")
    cyl.set("selresultshow", "all")


def create_difference(geom, tag: str, inputs: list[str], tools: list[str], label: str, keep_tools: bool = True) -> None:
    dif = geom.create(tag, "Difference")
    dif.label(f"{tag}: {label}")
    dif.selection("input").set(inputs)
    dif.selection("input2").set(tools)
    try:
        dif.set("keepadd", "off")
    except Exception:
        pass
    if keep_tools:
        try:
            dif.set("keepsubtract", "on")
        except Exception:
            pass
    dif.set("selresult", "on")
    dif.set("selresultshow", "all")


def create_union_selection(comp, tag: str, entitydim: str, inputs: list[str]) -> None:
    sel = comp.selection().create(tag, "Union")
    sel.label(tag)
    sel.set("entitydim", entitydim)
    sel.set("input", inputs)


def create_box_selection(comp, tag: str, entitydim: str, xmin: str, xmax: str, ymin: str, ymax: str, zmin: str, zmax: str) -> None:
    sel = comp.selection().create(tag, "Box")
    sel.label(tag)
    sel.set("entitydim", entitydim)
    sel.set("xmin", xmin)
    sel.set("xmax", xmax)
    sel.set("ymin", ymin)
    sel.set("ymax", ymax)
    sel.set("zmin", zmin)
    sel.set("zmax", zmax)
    sel.set("condition", "inside")


def create_refractive_material(comp, tag: str, label: str, selection: str, n: str, k: str = "0") -> None:
    mat = comp.material().create(tag, "Common")
    mat.label(label)
    mat.selection().named(selection)
    mat.propertyGroup().create("RefractiveIndex", "Refractive index")
    mat.propertyGroup("RefractiveIndex").set("n", [n])
    mat.propertyGroup("RefractiveIndex").set("ki", [k])


def unique_complex(values: Iterable[complex | float], tol: float = 1e-7) -> list[complex]:
    out: list[complex] = []
    for raw in values:
        val = complex(raw)
        if not any(abs(val - old) < tol for old in out):
            out.append(val)
    return out


def get_eigenfrequencies_thz(model) -> list[complex]:
    values = model.evaluate("freq", unit="THz")
    flat = values.ravel() if hasattr(values, "ravel") else values
    modes = unique_complex(flat)
    return [m for m in modes if math.isfinite(m.real) and abs(m.real) > 1e-9]


def mode_result(design: Design, target_nm: float, kx_norm: float, ky_norm: float, idx: int, freq_thz: complex) -> ModeResult:
    f_real = float(freq_thz.real)
    f_imag = float(freq_thz.imag)
    wavelength_nm = C0 / (abs(f_real) * 1e12) * 1e9 if abs(f_real) > 1e-15 else math.inf
    q_value = abs(f_real) / (2.0 * abs(f_imag)) if abs(f_imag) > 1e-18 else 1e18
    detune = abs(wavelength_nm - target_nm)
    target_q = Q_TARGET_451 if abs(target_nm - TARGET_451_NM) < 10 else Q_TARGET_980
    # The score favors near-target modes first, then high Q, with a penalty for tiny features.
    q_penalty = max(0.0, math.log10(target_q) - math.log10(max(q_value, 1.0)))
    fabrication_penalty = 5.0 if not design.fabricable else (0.5 if not design.preferred_gap_ok else 0.0)
    score = detune + 3.0 * q_penalty + fabrication_penalty
    return ModeResult(
        design.design_id,
        design.family,
        target_nm,
        design.px_nm,
        design.py_nm,
        design.feature_min_nm,
        design.fabricable,
        design.preferred_gap_ok,
        kx_norm,
        ky_norm,
        idx,
        f_real,
        f_imag,
        wavelength_nm,
        q_value,
        score,
    )


def best_mode(rows: list[ModeResult], target_nm: float) -> ModeResult:
    near = [r for r in rows if abs(r.wavelength_nm - target_nm) <= 6.0]
    pool = near or rows
    return sorted(pool, key=lambda r: (r.score, abs(r.wavelength_nm - target_nm), -r.q_value))[0]


def build_model(client, design: Design, target_nm: float, name: str):
    n_sin, n_sio2, n_re, n_si = target_indices(target_nm)
    model = client.create(name)
    jmodel = model.java

    half_x = design.px_nm / 2.0
    half_y = design.py_nm / 2.0
    z_si = H_PML_NM
    z_sio2 = H_PML_NM + H_SI_BUF_NM
    z_re = z_sio2 + H_SIO2_NM
    z_cap = z_re + H_RE_NM
    z_sin = z_cap + H_CAP_NM
    z_top = z_sin + H_SIN_NM
    z_top_pml = z_top + H_AIR_NM
    total_height = 2 * H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_RE_NM + H_CAP_NM + H_SIN_NM + H_AIR_NM

    set_param(model, "lambda_target", f"{target_nm:.9g}[nm]", "Eigenfrequency search wavelength")
    set_param(model, "Px", f"{design.px_nm:.9g}[nm]", "Unit-cell x period")
    set_param(model, "Py", f"{design.py_nm:.9g}[nm]", "Unit-cell y period")
    set_param(model, "h_sio2_eff", f"{H_SIO2_NM:.9g}[nm]", "Truncated SiO2 thickness")
    set_param(model, "h_re", f"{H_RE_NM:.9g}[nm]", "Continuous rare-earth film thickness")
    set_param(model, "h_cap_sio2", f"{H_CAP_NM:.9g}[nm]", "Continuous SiO2 cap thickness")
    set_param(model, "h_sin", f"{H_SIN_NM:.9g}[nm]", "Patterned SiN thickness")
    set_param(model, "z_re_mid", f"{z_re + H_RE_NM / 2.0:.9g}[nm]", "Rare-earth film mid-plane")
    set_param(model, "z_slab_mid", f"{z_sin + H_SIN_NM / 2.0:.9g}[nm]", "SiN mid-plane")
    set_param(model, "z_slab_top", f"{z_top:.9g}[nm]", "SiN top")
    set_param(model, "n_sin", f"{n_sin:.9f}", f"SiN index at {target_nm:.0f} nm")
    set_param(model, "n_sio2", f"{n_sio2:.9f}", f"SiO2 index at {target_nm:.0f} nm")
    set_param(model, "n_re", f"{n_re:.9f}", f"Rare-earth-film index at {target_nm:.0f} nm")
    set_param(model, "n_si", f"{n_si:.9f}", f"Effective Si index at {target_nm:.0f} nm")
    set_param(model, "kx", "0[1/m]", "Floquet kx")
    set_param(model, "ky", "0[1/m]", "Floquet ky")

    jmodel.component().create("comp1", True)
    comp = jmodel.component("comp1")
    comp.label(f"v10 {design.family} {target_nm:.0f} nm model")
    comp.geom().create("geom1", 3)
    geom = comp.geom("geom1")
    geom.lengthUnit("nm")

    create_block(geom, "si_pml", (-half_x, -half_y, 0.0), (design.px_nm, design.py_nm, H_PML_NM), "Si bottom PML")
    create_block(geom, "si_buf", (-half_x, -half_y, z_si), (design.px_nm, design.py_nm, H_SI_BUF_NM), "Si buffer")
    create_block(geom, "sio2", (-half_x, -half_y, z_sio2), (design.px_nm, design.py_nm, H_SIO2_NM), "1500 nm effective SiO2")
    create_block(geom, "re_film", (-half_x, -half_y, z_re), (design.px_nm, design.py_nm, H_RE_NM), "continuous rare-earth nanocrystal film")
    create_block(geom, "cap_sio2", (-half_x, -half_y, z_cap), (design.px_nm, design.py_nm, H_CAP_NM), "continuous SiO2 protection cap")

    air_domain_tags: list[str] = []
    sin_domain_tags: list[str] = []
    if design.is_pillar:
        tool_tags: list[str] = []
        for i, p in enumerate(design.pillars):
            tag = f"pillar{i}"
            create_cylinder(geom, tag, p.x_nm, p.y_nm, z_sin, p.radius_nm, H_SIN_NM, "SiN pillar")
            sin_domain_tags.append(f"geom1_{tag}_dom")
            tool_tags.append(tag)
        create_block(geom, "air_sin_blk", (-half_x, -half_y, z_sin), (design.px_nm, design.py_nm, H_SIN_NM), "air between SiN pillars")
        create_difference(geom, "air_sin", ["air_sin_blk"], tool_tags, "air surrounding SiN pillars", keep_tools=True)
        air_domain_tags.append("geom1_air_sin_dom")
    else:
        create_block(geom, "sin_slab_blk", (-half_x, -half_y, z_sin), (design.px_nm, design.py_nm, H_SIN_NM), "SiN slab before through-etch pattern")
        tool_tags = []
        for i, h in enumerate(design.circle_holes):
            tag = f"hole{i}"
            create_cylinder(geom, tag, h.x_nm, h.y_nm, z_sin, h.radius_nm, H_SIN_NM, "through-etched circular air hole in SiN only")
            tool_tags.append(tag)
            air_domain_tags.append(f"geom1_{tag}_dom")
        for i, h in enumerate(design.rect_holes):
            tag = f"slot{i}"
            create_block(
                geom,
                tag,
                (h.x_nm - h.width_nm / 2.0, h.y_nm - h.depth_nm / 2.0, z_sin),
                (h.width_nm, h.depth_nm, H_SIN_NM),
                "through-etched rectangular air slot in SiN only",
            )
            tool_tags.append(tag)
            air_domain_tags.append(f"geom1_{tag}_dom")
        for i, h in enumerate(design.capsule_holes):
            axis = h.axis.lower()
            cap_r = h.width_nm / 2.0
            mid_len = max(h.width_nm, h.length_nm - h.width_nm)
            if axis == "x":
                rect_w, rect_d = mid_len, h.width_nm
                c1 = (h.x_nm - mid_len / 2.0, h.y_nm)
                c2 = (h.x_nm + mid_len / 2.0, h.y_nm)
            else:
                rect_w, rect_d = h.width_nm, mid_len
                c1 = (h.x_nm, h.y_nm - mid_len / 2.0)
                c2 = (h.x_nm, h.y_nm + mid_len / 2.0)
            rect_tag = f"cap{i}_mid"
            cap1_tag = f"cap{i}_a"
            cap2_tag = f"cap{i}_b"
            create_block(
                geom,
                rect_tag,
                (h.x_nm - rect_w / 2.0, h.y_nm - rect_d / 2.0, z_sin),
                (rect_w, rect_d, H_SIN_NM),
                "through-etched capsule/ellipse-approximation mid slot in SiN only",
            )
            create_cylinder(geom, cap1_tag, c1[0], c1[1], z_sin, cap_r, H_SIN_NM, "capsule/ellipse-approximation round end")
            create_cylinder(geom, cap2_tag, c2[0], c2[1], z_sin, cap_r, H_SIN_NM, "capsule/ellipse-approximation round end")
            tool_tags.extend([rect_tag, cap1_tag, cap2_tag])
            air_domain_tags.extend([f"geom1_{rect_tag}_dom", f"geom1_{cap1_tag}_dom", f"geom1_{cap2_tag}_dom"])
        if tool_tags:
            create_difference(geom, "sin_slab", ["sin_slab_blk"], tool_tags, "patterned SiN slab", keep_tools=True)
            sin_domain_tags.append("geom1_sin_slab_dom")
        else:
            sin_domain_tags.append("geom1_sin_slab_blk_dom")

    create_block(geom, "air_upper", (-half_x, -half_y, z_top), (design.px_nm, design.py_nm, H_AIR_NM), "air above SiN")
    create_block(geom, "air_pml", (-half_x, -half_y, z_top_pml), (design.px_nm, design.py_nm, H_PML_NM), "top air PML")
    air_domain_tags.extend(["geom1_air_upper_dom", "geom1_air_pml_dom"])
    geom.run()

    create_union_selection(comp, "air_dom", "3", air_domain_tags)
    create_union_selection(comp, "sin_dom", "3", sin_domain_tags)
    create_union_selection(comp, "si_dom", "3", ["geom1_si_pml_dom", "geom1_si_buf_dom"])
    create_box_selection(comp, "x_left_bnd", "2", f"{-half_x-1}", f"{-half_x+1}", f"{-half_y-1}", f"{half_y+1}", "-1", f"{total_height+1}")
    create_box_selection(comp, "x_right_bnd", "2", f"{half_x-1}", f"{half_x+1}", f"{-half_y-1}", f"{half_y+1}", "-1", f"{total_height+1}")
    create_union_selection(comp, "x_periodic_bnd", "2", ["x_left_bnd", "x_right_bnd"])
    create_box_selection(comp, "y_front_bnd", "2", f"{-half_x-1}", f"{half_x+1}", f"{-half_y-1}", f"{-half_y+1}", "-1", f"{total_height+1}")
    create_box_selection(comp, "y_back_bnd", "2", f"{-half_x-1}", f"{half_x+1}", f"{half_y-1}", f"{half_y+1}", "-1", f"{total_height+1}")
    create_union_selection(comp, "y_periodic_bnd", "2", ["y_front_bnd", "y_back_bnd"])

    pml_bot = comp.coordSystem().create("pml_bot", "PML")
    pml_bot.label("Bottom Si PML")
    pml_bot.selection().named("geom1_si_pml_dom")
    pml_top = comp.coordSystem().create("pml_top", "PML")
    pml_top.label("Top air PML")
    pml_top.selection().named("geom1_air_pml_dom")

    create_refractive_material(comp, "mat_air", "Air in etched SiN pattern and top region", "air_dom", "1", "0")
    create_refractive_material(comp, "mat_si", "Effective Si", "si_dom", "n_si", "0")
    create_refractive_material(comp, "mat_sio2", "SiO2 substrate", "geom1_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_re", "50 nm rare-earth nanocrystal film", "geom1_re_film_dom", "n_re", "0")
    create_refractive_material(comp, "mat_cap", "12.5 nm SiO2 cap", "geom1_cap_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_sin", "Patterned SiN layer", "sin_dom", "n_sin", "0")

    ewfd = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
    ewfd.label("Electromagnetic Waves, Frequency Domain")
    pcx = ewfd.create("pcx", "PeriodicCondition", 2)
    pcx.label("Floquet periodic x")
    pcx.selection().named("x_periodic_bnd")
    pcx.set("PeriodicType", "Floquet")
    pcx.set("kFloquet", ["kx", "ky", "0"])
    pcy = ewfd.create("pcy", "PeriodicCondition", 2)
    pcy.label("Floquet periodic y")
    pcy.selection().named("y_periodic_bnd")
    pcy.set("PeriodicType", "Floquet")
    pcy.set("kFloquet", ["kx", "ky", "0"])

    mesh = comp.mesh().create("mesh1")
    mesh.label("Memory-controlled 3D mesh")
    mesh.autoMeshSize(MESH_AUTO_SIZE)

    study = jmodel.study().create("std1")
    study.label(f"Eigenfrequency near {target_nm:.0f} nm")
    eig = study.create("eig", "Eigenfrequency")
    eig.set("shift", "c_const/lambda_target")
    eig.set("neigs", str(NEIGS))
    eig.set("neigsmanual", str(NEIGS))
    eig.set("neigsactive", "on")
    eig.set("eigunit", "THz")
    return model


def solve_model(model, design: Design, target_nm: float, kx_norm: float = 0.0, ky_norm: float = 0.0) -> list[ModeResult]:
    set_param(model, "kx", f"{kx_norm:.12g}*pi/Px")
    set_param(model, "ky", f"{ky_norm:.12g}*pi/Py")
    model.java.component("comp1").geom("geom1").run()
    model.java.component("comp1").mesh("mesh1").run()
    model.java.study("std1").run()
    modes = get_eigenfrequencies_thz(model)
    return [mode_result(design, target_nm, kx_norm, ky_norm, i + 1, freq) for i, freq in enumerate(modes)]


def solve_design(client, design: Design, target_nm: float, kx_norm: float = 0.0, ky_norm: float = 0.0) -> list[ModeResult]:
    model = build_model(client, design, target_nm, f"v10_{design.design_id}_{int(target_nm)}")
    try:
        return solve_model(model, design, target_nm, kx_norm, ky_norm)
    finally:
        client.remove(model.name())


def add_field_plot_groups(model, label: str) -> None:
    res = model.java.result()
    for plane, tag, quick, fixed in [
        ("xy", f"field_{label}_xy", "xy", ("quickz", "z_re_mid")),
        ("xz", f"field_{label}_xz", "zx", ("quicky", "0")),
    ]:
        pg = res.create(tag, "PlotGroup3D")
        pg.label(f"{label}: |E|^2 {plane} slice")
        pg.set("data", "dset1")
        slc = pg.feature().create("slc", "Slice")
        slc.label(f"|E|^2 {plane}")
        slc.set("expr", "ewfd.normE^2")
        slc.set("quickplane", quick)
        slc.set(fixed[0], fixed[1])
        slc.set("resolution", "normal")


def save_solved_model(client, design: Design, target_nm: float, path: Path, model_name: str) -> ModeResult:
    model = build_model(client, design, target_nm, model_name)
    try:
        rows = solve_model(model, design, target_nm, 0.0, 0.0)
        selected = best_mode(rows, target_nm)
        add_field_plot_groups(model, f"{int(target_nm)}")
        model.save(str(path))
        write_modes_csv(path.with_suffix(".modes.csv"), rows)
        write_modes_csv(path.with_suffix(".selected.csv"), [selected])
        return selected
    finally:
        client.remove(model.name())


def dual_score(design: Design, m980: ModeResult, m451: ModeResult) -> float:
    det980 = abs(m980.wavelength_nm - TARGET_980_NM) / 2.0
    det451 = abs(m451.wavelength_nm - TARGET_451_NM) / 1.0
    q980_pen = max(0.0, math.log10(Q_TARGET_980) - math.log10(max(m980.q_value, 1.0)))
    q451_pen = max(0.0, math.log10(Q_TARGET_451) - math.log10(max(m451.q_value, 1.0)))
    fab_pen = 20.0 if not design.fabricable else (2.0 if not design.preferred_gap_ok else 0.0)
    return det980 + det451 + 6.0 * q980_pen + 3.0 * q451_pen + fab_pen


def scan_designs(client, designs: list[Design], path: Path, limit: int = 0) -> list[DualCandidate]:
    if limit:
        designs = designs[:limit]
    existing = read_dual_csv(path)
    done = {c.design.design_id for c in existing}
    out = existing[:]
    for i, design in enumerate(designs, 1):
        if design.design_id in done:
            continue
        print(
            f"[scan] {i}/{len(designs)} {design.family} {design.design_id} "
            f"Px={design.px_nm:.1f} Py={design.py_nm:.1f} feature_min={design.feature_min_nm:.1f}",
            flush=True,
        )
        try:
            rows980 = solve_design(client, design, TARGET_980_NM)
            m980 = best_mode(rows980, TARGET_980_NM)
            print(f"  980: lambda={m980.wavelength_nm:.6f} nm Q={m980.q_value:.4g} mode={m980.mode_index}", flush=True)
            rows451 = solve_design(client, design, TARGET_451_NM)
            m451 = best_mode(rows451, TARGET_451_NM)
            print(f"  451: lambda={m451.wavelength_nm:.6f} nm Q={m451.q_value:.4g} mode={m451.mode_index}", flush=True)
            out.append(DualCandidate(design, m980, m451, dual_score(design, m980, m451)))
            write_dual_csv(path, out)
        except Exception as exc:
            print(f"  FAILED {design.design_id}: {exc}", flush=True)
    return out


def write_modes_csv(path: Path, rows: list[ModeResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "design_id",
        "family",
        "target_nm",
        "px_nm",
        "py_nm",
        "feature_min_nm",
        "fabricable",
        "preferred_gap_ok",
        "kx_norm",
        "ky_norm",
        "mode_index",
        "freq_thz_real",
        "freq_thz_imag",
        "wavelength_nm",
        "q_value",
        "score",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({f: getattr(r, f) for f in fields})


def read_dual_csv(path: Path) -> list[DualCandidate]:
    if not path.exists():
        return []
    rows: list[DualCandidate] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            design = design_by_id(row["design_id"])
            if design is None:
                continue
            m980 = ModeResult(
                row["design_id"],
                row["family"],
                TARGET_980_NM,
                float(row["px_nm"]),
                float(row["py_nm"]),
                float(row["feature_min_nm"]),
                row["fabricable"].lower() == "true",
                row["preferred_gap_ok"].lower() == "true",
                0.0,
                0.0,
                int(row["mode980_index"]),
                float(row["freq980_thz_real"]),
                float(row["freq980_thz_imag"]),
                float(row["lambda980_nm"]),
                float(row["q980"]),
                0.0,
            )
            m451 = ModeResult(
                row["design_id"],
                row["family"],
                TARGET_451_NM,
                float(row["px_nm"]),
                float(row["py_nm"]),
                float(row["feature_min_nm"]),
                row["fabricable"].lower() == "true",
                row["preferred_gap_ok"].lower() == "true",
                0.0,
                0.0,
                int(row["mode451_index"]),
                float(row["freq451_thz_real"]),
                float(row["freq451_thz_imag"]),
                float(row["lambda451_nm"]),
                float(row["q451"]),
                0.0,
            )
            rows.append(DualCandidate(design, m980, m451, float(row["dual_score"])))
    return rows


def write_dual_csv(path: Path, rows: list[DualCandidate]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "design_id",
        "family",
        "notes",
        "px_nm",
        "py_nm",
        "feature_min_nm",
        "fabricable",
        "preferred_gap_ok",
        "lambda980_nm",
        "q980",
        "mode980_index",
        "freq980_thz_real",
        "freq980_thz_imag",
        "lambda451_nm",
        "q451",
        "mode451_index",
        "freq451_thz_real",
        "freq451_thz_imag",
        "dual_score",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for c in sorted(rows, key=lambda item: item.score):
            writer.writerow(
                {
                    "design_id": c.design.design_id,
                    "family": c.design.family,
                    "notes": c.design.notes,
                    "px_nm": c.design.px_nm,
                    "py_nm": c.design.py_nm,
                    "feature_min_nm": c.design.feature_min_nm,
                    "fabricable": c.design.fabricable,
                    "preferred_gap_ok": c.design.preferred_gap_ok,
                    "lambda980_nm": c.mode980.wavelength_nm,
                    "q980": c.mode980.q_value,
                    "mode980_index": c.mode980.mode_index,
                    "freq980_thz_real": c.mode980.freq_thz_real,
                    "freq980_thz_imag": c.mode980.freq_thz_imag,
                    "lambda451_nm": c.mode451.wavelength_nm,
                    "q451": c.mode451.q_value,
                    "mode451_index": c.mode451.mode_index,
                    "freq451_thz_real": c.mode451.freq_thz_real,
                    "freq451_thz_imag": c.mode451.freq_thz_imag,
                    "dual_score": c.score,
                }
            )


def baseline_and_literature_designs() -> list[Design]:
    designs: list[Design] = []

    # Air-hole slab near old 980 nm high-Q scale, used only as a reference branch.
    for p, r in [(520, 66), (540, 74), (560, 82), (620, 100), (700, 130)]:
        designs.append(
            Design(
                f"mono_hole_P{p}_R{r}",
                "single_air_hole_slab_reference",
                float(p),
                float(p),
                circle_holes=(CircleHole(0.0, 0.0, float(r)),),
                notes="single circular air-hole PhC slab reference",
            )
        )

    # Rectangular-lattice single-hole references around the high-Q v8 scale.
    for px, py, r in [
        (540.0, 460.0, 74.0),
        (560.0, 500.0, 82.0),
        (620.0, 520.0, 96.0),
        (700.0, 560.0, 120.0),
    ]:
        designs.append(
            Design(
                f"rect_lattice_hole_Px{px:.0f}_Py{py:.0f}_R{r:.0f}",
                "rectangular_air_hole_slab",
                px,
                py,
                circle_holes=(CircleHole(0.0, 0.0, r),),
                notes="rectangular-lattice single circular air-hole PhC slab",
            )
        )

    # Brillouin-zone folded two-hole supercells: alternating radii and/or shifted sites.
    for base_p in [450.0, 500.0, 540.0, 600.0, 700.0]:
        for r1, r2 in [(0.18, 0.28), (0.20, 0.32), (0.24, 0.34)]:
            rr1 = base_p * r1
            rr2 = base_p * r2
            if min(base_p - 2 * rr1, base_p - 2 * rr2) < MIN_FEATURE_NM:
                continue
            designs.append(
                Design(
                    f"dual_radius_P{base_p:.0f}_r{rr1:.0f}_{rr2:.0f}",
                    "dual_hole_radius_supercell",
                    2 * base_p,
                    base_p,
                    circle_holes=(CircleHole(-base_p / 2, 0.0, rr1), CircleHole(base_p / 2, 0.0, rr2)),
                    notes="2x1 supercell with alternating air-hole radii",
                )
            )
        for shift in [20.0, 35.0, 50.0]:
            r = base_p * 0.24
            designs.append(
                Design(
                    f"dual_shift_P{base_p:.0f}_R{r:.0f}_dx{shift:.0f}",
                    "dual_hole_shift_supercell",
                    2 * base_p,
                    base_p,
                    circle_holes=(CircleHole(-base_p / 2 - shift / 2, 0.0, r), CircleHole(base_p / 2 + shift / 2, 0.0, r)),
                    notes="2x1 supercell with position dimerization",
                )
            )

    # 2x2 folded supercells with alternating radii or deliberate position offsets.
    for base_p, r_a, r_b in [
        (420.0, 70.0, 95.0),
        (500.0, 88.0, 125.0),
        (540.0, 95.0, 135.0),
        (600.0, 110.0, 155.0),
    ]:
        holes = (
            CircleHole(-base_p / 2, -base_p / 2, r_a),
            CircleHole(base_p / 2, -base_p / 2, r_b),
            CircleHole(-base_p / 2, base_p / 2, r_b),
            CircleHole(base_p / 2, base_p / 2, r_a),
        )
        designs.append(
            Design(
                f"fold2x2_radius_P{base_p:.0f}_R{r_a:.0f}_{r_b:.0f}",
                "folded_2x2_air_hole_radius",
                2 * base_p,
                2 * base_p,
                circle_holes=holes,
                notes="2x2 BZ-folded air-hole slab with alternating radii",
            )
        )
    for base_p, r, shift in [
        (500.0, 105.0, 35.0),
        (540.0, 118.0, 45.0),
        (600.0, 135.0, 55.0),
    ]:
        holes = (
            CircleHole(-base_p / 2 - shift / 2, -base_p / 2, r),
            CircleHole(base_p / 2 + shift / 2, -base_p / 2, r),
            CircleHole(-base_p / 2 + shift / 2, base_p / 2, r),
            CircleHole(base_p / 2 - shift / 2, base_p / 2, r),
        )
        designs.append(
            Design(
                f"fold2x2_shift_P{base_p:.0f}_R{r:.0f}_d{shift:.0f}",
                "folded_2x2_air_hole_shift",
                2 * base_p,
                2 * base_p,
                circle_holes=holes,
                notes="2x2 BZ-folded air-hole slab with sublattice position offsets",
            )
        )

    # Elliptical-hole routes approximated by manufacturable capsule holes.
    for px, py, length, width, axis in [
        (620.0, 500.0, 270.0, 120.0, "x"),
        (700.0, 520.0, 320.0, 140.0, "x"),
        (700.0, 560.0, 300.0, 130.0, "y"),
        (840.0, 620.0, 380.0, 160.0, "x"),
    ]:
        designs.append(
            Design(
                f"capsule_hole_Px{px:.0f}_Py{py:.0f}_L{length:.0f}_W{width:.0f}_{axis}",
                "capsule_elliptic_air_hole",
                px,
                py,
                capsule_holes=(CapsuleHole(0.0, 0.0, length, width, axis),),
                notes="single elliptical/capsule air hole approximation through SiN",
            )
        )
    for base_p, length1, width1, length2, width2 in [
        (520.0, 230.0, 110.0, 280.0, 130.0),
        (600.0, 260.0, 120.0, 330.0, 150.0),
        (700.0, 310.0, 140.0, 390.0, 170.0),
    ]:
        designs.append(
            Design(
                f"dual_capsule_P{base_p:.0f}_L{length1:.0f}_{length2:.0f}",
                "dual_capsule_elliptic_supercell",
                2 * base_p,
                base_p,
                capsule_holes=(
                    CapsuleHole(-base_p / 2, 0.0, length1, width1, "x"),
                    CapsuleHole(base_p / 2, 0.0, length2, width2, "x"),
                ),
                notes="2x1 supercell with two different elliptical/capsule air holes",
            )
        )

    # Circular hole plus manufacturable offset satellite hole.
    for p, r_main, r_sat, dx, dy in [
        (620.0, 105.0, 55.0, 170.0, 0.0),
        (700.0, 125.0, 60.0, 200.0, 0.0),
        (760.0, 135.0, 65.0, 190.0, 120.0),
    ]:
        designs.append(
            Design(
                f"main_sat_hole_P{p:.0f}_R{r_main:.0f}_{r_sat:.0f}",
                "main_hole_with_offset_satellite",
                p,
                p,
                circle_holes=(CircleHole(0.0, 0.0, r_main), CircleHole(dx, dy, r_sat)),
                notes="circular air hole with an offset satellite hole, all features >=100 nm",
            )
        )

    # Slot and bow-tie-inspired wide rectangular cuts. These are deliberately not sub-100 nm grooves.
    for px, py, w, d in [
        (620.0, 420.0, 180.0, 260.0),
        (700.0, 500.0, 220.0, 300.0),
        (840.0, 540.0, 260.0, 340.0),
        (980.0, 600.0, 320.0, 380.0),
    ]:
        designs.append(
            Design(
                f"slot_cross_Px{px:.0f}_Py{py:.0f}",
                "wide_slot_cross_hole",
                px,
                py,
                rect_holes=(RectHole(0.0, 0.0, w, d), RectHole(0.0, 0.0, d, w)),
                notes="cross-like wide rectangular air opening through SiN",
            )
        )

    # 1D dual-period grating-waveguide approximated as a 3D unit cell with two through-slots.
    for p1, p2, slot in [(360.0, 430.0, 140.0), (420.0, 520.0, 160.0), (520.0, 650.0, 190.0), (620.0, 760.0, 220.0)]:
        px = p1 + p2
        py = 360.0
        x0 = -px / 2.0
        c1 = x0 + p1 / 2.0
        c2 = x0 + p1 + p2 / 2.0
        designs.append(
            Design(
                f"grating_dual_P{p1:.0f}_{p2:.0f}_S{slot:.0f}",
                "dual_period_1d_grating",
                px,
                py,
                rect_holes=(RectHole(c1, 0.0, slot, py + 20.0), RectHole(c2, 0.0, slot * 0.75, py + 20.0)),
                notes="1D two-pitch grating slab, implemented as through-slots",
            )
        )

    # Four-hole and four-pillar tetramer/quadrumer routes.
    for p, r_big, r_small, sep in [
        (900.0, 120.0, 78.0, 380.0),
        (1000.0, 135.0, 90.0, 430.0),
        (1200.0, 165.0, 110.0, 520.0),
    ]:
        holes = (
            CircleHole(-sep / 2, -sep / 2, r_big),
            CircleHole(sep / 2, -sep / 2, r_small),
            CircleHole(-sep / 2, sep / 2, r_small),
            CircleHole(sep / 2, sep / 2, r_big),
        )
        designs.append(
            Design(
                f"tetramer_hole_P{p:.0f}_R{r_big:.0f}_{r_small:.0f}",
                "tetramer_hole_slab",
                p,
                p,
                circle_holes=holes,
                notes="2x2-like alternating-radius four-hole slab",
            )
        )
    for p, r_big, r_small, sep in [
        (800.0, 110.0, 75.0, 330.0),
        (980.0, 145.0, 95.0, 420.0),
        (1180.0, 180.0, 120.0, 520.0),
    ]:
        pillars = (
            Pillar(-sep / 2, -sep / 2, r_big),
            Pillar(sep / 2, -sep / 2, r_small),
            Pillar(-sep / 2, sep / 2, r_small),
            Pillar(sep / 2, sep / 2, r_big),
        )
        designs.append(
            Design(
                f"tetramer_pillar_P{p:.0f}_R{r_big:.0f}_{r_small:.0f}",
                "tetramer_pillar_array",
                p,
                p,
                pillars=pillars,
                notes="four SiN pillars on intact cap/RE stack",
            )
        )

    filtered: list[Design] = []
    seen: set[str] = set()
    for d in designs:
        if d.design_id not in seen and d.fabricable:
            filtered.append(d)
            seen.add(d.design_id)
    return filtered


def design_by_id(design_id: str) -> Design | None:
    for d in baseline_and_literature_designs():
        if d.design_id == design_id:
            return d
    return None


def choose_best_dual(rows: list[DualCandidate]) -> DualCandidate:
    feasible = [r for r in rows if r.design.fabricable]
    pool = feasible or rows
    strict = [
        r
        for r in pool
        if abs(r.mode980.wavelength_nm - TARGET_980_NM) <= 2.0
        and r.mode980.q_value >= Q_TARGET_980
        and abs(r.mode451.wavelength_nm - TARGET_451_NM) <= 2.0
        and r.mode451.q_value >= Q_TARGET_451
    ]
    return sorted(strict or pool, key=lambda r: r.score)[0]


def choose_family_bests(rows: list[DualCandidate]) -> list[DualCandidate]:
    out: list[DualCandidate] = []
    for family in sorted({r.design.family for r in rows}):
        fam = [r for r in rows if r.design.family == family]
        out.append(sorted(fam, key=lambda r: r.score)[0])
    return sorted(out, key=lambda r: r.score)


def refine_designs(seed: DualCandidate) -> list[Design]:
    d = seed.design
    out: list[Design] = []
    if d.family.startswith("dual_hole"):
        base_p = d.py_nm
        holes = d.circle_holes
        if len(holes) >= 2:
            for dp in [-30.0, -15.0, 0.0, 15.0, 30.0]:
                for dr1 in [-12.0, 0.0, 12.0]:
                    for dr2 in [-12.0, 0.0, 12.0]:
                        p = max(360.0, base_p + dp)
                        r1 = max(50.0, holes[0].radius_nm + dr1)
                        r2 = max(50.0, holes[1].radius_nm + dr2)
                        nd = Design(
                            f"refine_{d.family}_P{p:.0f}_R{r1:.0f}_{r2:.0f}",
                            d.family,
                            2 * p,
                            p,
                            circle_holes=(CircleHole(-p / 2, 0.0, r1), CircleHole(p / 2, 0.0, r2)),
                            notes="local refinement around best dual-hole branch",
                        )
                        if nd.fabricable:
                            out.append(nd)
    elif d.family in {"wide_slot_cross_hole", "dual_period_1d_grating"}:
        for sx in [0.92, 1.0, 1.08]:
            for sw in [0.88, 1.0, 1.12]:
                rects = tuple(RectHole(h.x_nm * sx, h.y_nm, max(100.0, h.width_nm * sw), max(100.0, h.depth_nm * sw)) for h in d.rect_holes)
                nd = Design(
                    f"refine_{d.family}_Px{d.px_nm*sx:.0f}_S{sw:.2f}".replace(".", "p"),
                    d.family,
                    d.px_nm * sx,
                    d.py_nm,
                    rect_holes=rects,
                    notes="local refinement around best slot/grating branch",
                )
                if nd.fabricable:
                    out.append(nd)
    elif d.family.startswith("tetramer"):
        for sp in [0.92, 1.0, 1.08]:
            for sr in [0.88, 1.0, 1.12]:
                if d.circle_holes:
                    holes = tuple(CircleHole(h.x_nm * sp, h.y_nm * sp, max(50.0, h.radius_nm * sr)) for h in d.circle_holes)
                    nd = Design(
                        f"refine_{d.family}_P{d.px_nm*sp:.0f}_R{sr:.2f}".replace(".", "p"),
                        d.family,
                        d.px_nm * sp,
                        d.py_nm * sp,
                        circle_holes=holes,
                        notes="local refinement around tetramer holes",
                    )
                else:
                    pillars = tuple(Pillar(p.x_nm * sp, p.y_nm * sp, max(50.0, p.radius_nm * sr)) for p in d.pillars)
                    nd = Design(
                        f"refine_{d.family}_P{d.px_nm*sp:.0f}_R{sr:.2f}".replace(".", "p"),
                        d.family,
                        d.px_nm * sp,
                        d.py_nm * sp,
                        pillars=pillars,
                        notes="local refinement around tetramer pillars",
                    )
                if nd.fabricable:
                    out.append(nd)
    return out


def run_band(client, design: Design, target_nm: float, csv_path: Path, q_path: Path) -> list[ModeResult]:
    model = build_model(client, design, target_nm, f"v10_band_{design.design_id}_{int(target_nm)}")
    kpath = [
        (0.0, 0.0),
        (0.005, 0.0),
        (0.02, 0.0),
        (0.06, 0.0),
        (0.15, 0.0),
        (0.30, 0.0),
        (0.50, 0.0),
        (0.50, 0.25),
        (0.50, 0.50),
        (0.25, 0.25),
        (0.0, 0.0),
    ]
    rows: list[ModeResult] = []
    try:
        for kx, ky in kpath:
            print(f"[band {int(target_nm)}] k=({kx:g},{ky:g})", flush=True)
            current = solve_model(model, design, target_nm, kx, ky)
            rows.append(best_mode(current, target_nm))
        write_modes_csv(csv_path, rows)
        write_modes_csv(q_path, [r for r in rows if r.ky_norm == 0.0 and r.kx_norm <= 0.06])
        return rows
    finally:
        client.remove(model.name())


def svg_dual_candidates(rows: list[DualCandidate], path: Path) -> None:
    if not rows:
        return
    width, height = 920, 560
    left, right, top, bottom = 95, 40, 45, 85
    pw, ph = width - left - right, height - top - bottom
    l980 = [r.mode980.wavelength_nm for r in rows]
    q980 = [r.mode980.q_value for r in rows]
    q451 = [r.mode451.q_value for r in rows]
    xvals = list(range(len(rows)))
    qmin = max(1.0, min(min(q980), min(q451)) / 2)
    qmax = max(max(q980), max(q451)) * 2

    def sx(i: int) -> float:
        return left + i * pw / max(1, len(rows) - 1)

    def sy(q: float) -> float:
        return top + (math.log10(qmax) - math.log10(max(q, qmin))) / (math.log10(qmax) - math.log10(qmin)) * ph

    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    lines.append(f'<text x="{width/2}" y="27" text-anchor="middle" font-family="Arial" font-size="18">v10 dual-candidate Q comparison</text>')
    lines.append(f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#111"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#111"/>')
    for decade in range(int(math.floor(math.log10(qmin))), int(math.ceil(math.log10(qmax))) + 1):
        qv = 10**decade
        if qmin <= qv <= qmax:
            lines.append(f'<line x1="{left}" y1="{sy(qv):.2f}" x2="{left+pw}" y2="{sy(qv):.2f}" stroke="#eee"/>')
            lines.append(f'<text x="{left-8}" y="{sy(qv)+4:.2f}" text-anchor="end" font-family="Arial" font-size="12">1e{decade}</text>')
    for i, r in enumerate(rows):
        x = sx(i)
        lines.append(f'<circle cx="{x:.2f}" cy="{sy(r.mode980.q_value):.2f}" r="4" fill="#1f77b4"><title>{r.design.design_id} 980 lambda={r.mode980.wavelength_nm:.4g} Q={r.mode980.q_value:.4g}</title></circle>')
        lines.append(f'<circle cx="{x:.2f}" cy="{sy(r.mode451.q_value):.2f}" r="4" fill="#d62728"><title>{r.design.design_id} 451 lambda={r.mode451.wavelength_nm:.4g} Q={r.mode451.q_value:.4g}</title></circle>')
        lines.append(f'<text x="{x:.2f}" y="{top+ph+20}" text-anchor="end" transform="rotate(-35 {x:.2f} {top+ph+20})" font-family="Arial" font-size="10">{r.design.family[:18]}</text>')
    lines.append(f'<text x="{left+pw/2}" y="{height-20}" text-anchor="middle" font-family="Arial" font-size="14">candidate index, sorted by score</text>')
    lines.append(f'<text x="24" y="{top+ph/2}" transform="rotate(-90 24 {top+ph/2})" text-anchor="middle" font-family="Arial" font-size="14">Q (log)</text>')
    lines.append(f'<text x="{width-170}" y="55" font-family="Arial" font-size="13" fill="#1f77b4">blue: 980 nm branch</text>')
    lines.append(f'<text x="{width-170}" y="75" font-family="Arial" font-size="13" fill="#d62728">red: 451 nm branch</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def svg_band(rows: list[ModeResult], path: Path, title: str) -> None:
    if not rows:
        return
    width, height = 900, 520
    left, right, top, bottom = 80, 85, 40, 85
    pw, ph = width - left - right, height - top - bottom
    wl = [r.wavelength_nm for r in rows]
    qs = [r.q_value for r in rows]
    wmin, wmax = min(wl) - 2.0, max(wl) + 2.0
    qmin, qmax = max(1.0, min(qs) / 2.0), max(qs) * 2.0

    def sx(i: int) -> float:
        return left + i * pw / max(1, len(rows) - 1)

    def sy_w(v: float) -> float:
        return top + (wmax - v) / (wmax - wmin) * ph

    def sy_q(v: float) -> float:
        return top + (math.log10(qmax) - math.log10(max(v, qmin))) / (math.log10(qmax) - math.log10(qmin)) * ph

    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    lines.append(f'<text x="{width/2}" y="25" text-anchor="middle" font-family="Arial" font-size="18">{title}</text>')
    lines.append(f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#111"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#111"/>')
    lines.append(f'<line x1="{left+pw}" y1="{top}" x2="{left+pw}" y2="{top+ph}" stroke="#111"/>')
    wpts = " ".join(f"{sx(i):.2f},{sy_w(r.wavelength_nm):.2f}" for i, r in enumerate(rows))
    qpts = " ".join(f"{sx(i):.2f},{sy_q(r.q_value):.2f}" for i, r in enumerate(rows))
    lines.append(f'<polyline fill="none" stroke="#1f77b4" stroke-width="2" points="{wpts}"/>')
    lines.append(f'<polyline fill="none" stroke="#d62728" stroke-width="2" points="{qpts}"/>')
    for i, r in enumerate(rows):
        label = f"({r.kx_norm:g},{r.ky_norm:g})"
        lines.append(f'<text x="{sx(i):.2f}" y="{top+ph+20}" text-anchor="end" transform="rotate(-35 {sx(i):.2f} {top+ph+20})" font-family="Arial" font-size="10">{label}</text>')
    lines.append(f'<text x="22" y="{top+ph/2}" transform="rotate(-90 22 {top+ph/2})" text-anchor="middle" font-family="Arial" font-size="14" fill="#1f77b4">wavelength (nm)</text>')
    lines.append(f'<text x="{width-23}" y="{top+ph/2}" transform="rotate(90 {width-23} {top+ph/2})" text-anchor="middle" font-family="Arial" font-size="14" fill="#d62728">Q (log)</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def copy_reference_inputs() -> None:
    INPUTS.mkdir(parents=True, exist_ok=True)
    refs = [
        WORKSPACE / "reference" / "Symmetry-Breaking Papers Relevant to Dual-Wavelength High-Q Quasi-BIC Structures.pdf",
        WORKSPACE / "sin980_qbic_3d_design_v8_air_hole_slab" / "outputs" / "sio2_1500_check" / "sin980_qbic_v8_air_hole_sio2_1500_exact980_gap120.mph",
        WORKSPACE / "sin980_qbic_3d_design_v7_cyl_pillar_corrected_particles" / "outputs" / "sio2_1500_gap120_opt" / "sin980_qbic_v7_cyl_sio2_1500_gap120_best.mph",
    ]
    for src in refs:
        if src.exists():
            dst = INPUTS / src.name
            if not dst.exists():
                shutil.copy2(src, dst)
    (INPUTS / "v10_material_stack_notes.txt").write_text(
        "\n".join(
            [
                "v10 material stack and design constraints",
                "Effective Si / 1500 nm SiO2 / 50 nm rare-earth film / 12.5 nm SiO2 cap / 300 nm patterned SiN / air.",
                "Only the SiN layer is patterned. The rare-earth film and SiO2 cap remain continuous.",
                "n_re(980)=1.355, n_re(451)=1.365.",
                "n_sin(980)=2.069946018, n_sin(451)=2.128247126.",
                "n_sio2(980)=1.45, n_sio2(451)=1.466.",
                "Targets: Q980>=1e8, Q451>=1e5.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def write_summary(best: DualCandidate, family_bests: list[DualCandidate], band980: list[ModeResult], band451: list[ModeResult]) -> None:
    passes = best.mode980.q_value >= Q_TARGET_980 and best.mode451.q_value >= Q_TARGET_451
    lines = [
        "# v10 50 nm RE Film + 12.5 nm SiO2 Cap Dual-BIC Search",
        "",
        f"All v10 files were generated inside `{ROOT}`. Earlier v1-v9 folders and original inputs were not overwritten.",
        "",
        "## Fixed Stack",
        "",
        "- Effective Si / 1500 nm SiO2 / 50 nm rare-earth nanocrystal film / 12.5 nm SiO2 cap / 300 nm patterned SiN / air.",
        "- Only the SiN layer is patterned. The rare-earth film and cap are continuous and unetched.",
        f"- RE index: n(980) = {N_RE_980}, n(451) = {N_RE_451}; lossless.",
        f"- SiN index: n(980) = {N_SIN_980}, n(451) = {N_SIN_451}; lossless.",
        f"- SiO2 index: n(980) = {N_SIO2_980}, n(451) = {N_SIO2_451}; lossless.",
        "",
        "## Best Candidate",
        "",
        f"- Design: `{best.design.design_id}`",
        f"- Family: `{best.design.family}`",
        f"- Px = {best.design.px_nm:.3f} nm, Py = {best.design.py_nm:.3f} nm",
        f"- Minimum SiN feature / bridge proxy = {best.design.feature_min_nm:.3f} nm",
        f"- Fabricable hard check >=100 nm: `{best.design.fabricable}`",
        f"- Preferred gap/bridge >=120 nm: `{best.design.preferred_gap_ok}`",
        f"- 980 branch: lambda = {best.mode980.wavelength_nm:.6f} nm, Q = {best.mode980.q_value:.6g}",
        f"- 451 branch: lambda = {best.mode451.wavelength_nm:.6f} nm, Q = {best.mode451.q_value:.6g}",
        f"- Meets requested Q thresholds: `{passes}`",
        "",
        "## Family Bests",
        "",
        "| family | design | lambda980 (nm) | Q980 | lambda451 (nm) | Q451 | min feature (nm) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for c in family_bests:
        lines.append(
            f"| {c.design.family} | {c.design.design_id} | {c.mode980.wavelength_nm:.4f} | {c.mode980.q_value:.4g} | "
            f"{c.mode451.wavelength_nm:.4f} | {c.mode451.q_value:.4g} | {c.design.feature_min_nm:.1f} |"
        )
    lines.extend(
        [
            "",
            "## BIC Evidence",
            "",
            "- The generated band and Q(k) CSVs use a short Gamma-X-M-Gamma path plus small-k points near Gamma.",
            "- A BIC-like branch should show highest Q at Gamma and a rapid Q decrease for small nonzero k.",
        ]
    )
    if band980:
        lines.append("- 980 Q(k) near Gamma:")
        for r in band980[:4]:
            lines.append(f"  - k=({r.kx_norm:g},{r.ky_norm:g}) pi/P: lambda {r.wavelength_nm:.6f} nm, Q {r.q_value:.6g}")
    if band451:
        lines.append("- 451 Q(k) near Gamma:")
        for r in band451[:4]:
            lines.append(f"  - k=({r.kx_norm:g},{r.ky_norm:g}) pi/P: lambda {r.wavelength_nm:.6f} nm, Q {r.q_value:.6g}")
    lines.extend(
        [
            "",
            "## Output Files",
            "",
            "- `overall_best_dual_bic.mph`: 980 nm material setting for the final geometry, with xy/xz field plot groups.",
            "- `overall_best_dual_bic_451check.mph`: same final geometry at 451 nm material setting.",
            "- Per-family best `.mph` files are saved as `family_best_<family>_980.mph` and `family_best_<family>_451.mph`.",
            "- `dual_best_candidates.csv`, `coarse_dual_scan.csv`, `fine_dual_scan.csv`, `band_980.csv`, `band_451.csv`, `q_near_gamma_980.csv`, `q_near_gamma_451.csv`.",
        ]
    )
    if not passes:
        lines.extend(
            [
                "",
                "## Current Limitation",
                "",
                "- This run did not yet prove a simultaneous high-Q pair at the requested thresholds.",
                "- The likely reason is that a single 300 nm SiN layer with an intact low-index 50 nm film and cap can separately support high-Q branches, but aligning a 980 nm folded branch and a 451 nm visible branch in the same manufacturable in-plane geometry is strongly constrained.",
                "- Next useful relaxations would be allowing a second SiN etch depth, a second patterned SiN level, or a larger supercell with more in-plane degrees of freedom.",
            ]
        )
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_all(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    INPUTS.mkdir(parents=True, exist_ok=True)
    copy_reference_inputs()
    mph = import_mph()
    client = mph.Client(cores=args.cores)
    try:
        designs = baseline_and_literature_designs()
        if args.family:
            designs = [d for d in designs if d.family == args.family]
        coarse_rows = scan_designs(client, designs, OUT / "coarse_dual_scan.csv", args.limit)
        if not coarse_rows:
            (OUT / "summary.md").write_text("# v10 dual BIC search\n\nNo candidates solved successfully.\n", encoding="utf-8")
            return
        family_bests = choose_family_bests(coarse_rows)
        write_dual_csv(OUT / "family_best_candidates.csv", family_bests)
        refine_seed = choose_best_dual(coarse_rows)
        fine_designs = refine_designs(refine_seed)
        fine_rows = scan_designs(client, fine_designs, OUT / "fine_dual_scan.csv", args.fine_limit) if fine_designs and not args.skip_fine else []
        all_rows = coarse_rows + fine_rows
        best = choose_best_dual(all_rows)
        family_bests = choose_family_bests(all_rows)
        write_dual_csv(OUT / "dual_best_candidates.csv", sorted(all_rows, key=lambda r: r.score)[:20])
        write_dual_csv(OUT / "family_best_candidates.csv", family_bests)
        svg_dual_candidates(sorted(all_rows, key=lambda r: r.score)[:30], OUT / "dual_q_candidates.svg")

        print(f"[save] overall best {best.design.design_id}", flush=True)
        save_solved_model(client, best.design, TARGET_980_NM, OUT / "overall_best_dual_bic.mph", "overall_best_dual_bic_980")
        save_solved_model(client, best.design, TARGET_451_NM, OUT / "overall_best_dual_bic_451check.mph", "overall_best_dual_bic_451")

        for fam_best in family_bests:
            safe_family = fam_best.design.family.replace("/", "_")
            save_solved_model(client, fam_best.design, TARGET_980_NM, OUT / f"family_best_{safe_family}_980.mph", f"family_best_{safe_family}_980")
            save_solved_model(client, fam_best.design, TARGET_451_NM, OUT / f"family_best_{safe_family}_451.mph", f"family_best_{safe_family}_451")

        band980: list[ModeResult] = []
        band451: list[ModeResult] = []
        if not args.skip_band:
            band980 = run_band(client, best.design, TARGET_980_NM, OUT / "band_980.csv", OUT / "q_near_gamma_980.csv")
            band451 = run_band(client, best.design, TARGET_451_NM, OUT / "band_451.csv", OUT / "q_near_gamma_451.csv")
            svg_band(band980, OUT / "band_980.svg", "v10 980 nm band and Q(k)")
            svg_band(band451, OUT / "band_451.svg", "v10 451 nm band and Q(k)")
        write_summary(best, family_bests, band980, band451)
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="v10 RE50 + SiO2 cap dual-wavelength SiN BIC search.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-all")
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--fine-limit", type=int, default=0)
    p.add_argument("--family", default="")
    p.add_argument("--skip-fine", action="store_true")
    p.add_argument("--skip-band", action="store_true")
    p.set_defaults(func=run_all)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
