from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
OUT = ROOT / "outputs"
INPUTS = ROOT / "inputs"

COMSOLROOT = r"D:\comsol63\Multiphysics"
JAVA_HOME = rf"{COMSOLROOT}\java\win64\jre"
C0 = 299_792_458.0

TARGET_NM = 980.0
N_SIN = 2.069946018
N_SIO2 = 1.45
N_SI = 3.55
N_RE = 1.52

H_SIN_NM = 300.0
H_RE_NM = 16.0
H_SIO2_NM = 700.0
H_SI_BUF_NM = 250.0
H_PML_NM = 500.0
H_AIR_NM = 900.0

PARTICLE_D_MIN_NM = 24.0
PARTICLE_D_MAX_NM = 55.0
MIN_AXIS_GAP_NM = 10.0
TARGET_Q = 1e9
ACCEPT_Q = 1e8


@dataclass(frozen=True)
class Disk:
    x_nm: float
    y_nm: float
    diameter_nm: float

    @property
    def radius_nm(self) -> float:
        return self.diameter_nm / 2.0


@dataclass(frozen=True)
class SideParticle:
    theta: float
    z_local_nm: float
    diameter_nm: float

    @property
    def radius_nm(self) -> float:
        return self.diameter_nm / 2.0


@dataclass(frozen=True)
class ParticleLayout:
    pillar_top: tuple[Disk, ...]
    sio2_top: tuple[Disk, ...]
    pillar_side: tuple[SideParticle, ...]

    @property
    def count(self) -> int:
        return len(self.pillar_top) + len(self.sio2_top) + len(self.pillar_side)


@dataclass(frozen=True)
class Design:
    design_id: str
    period_nm: float
    radius_nm: float
    seed: int = 0

    @property
    def axis_gap_nm(self) -> float:
        return self.period_nm - 2.0 * self.radius_nm

    @property
    def fill_area(self) -> float:
        return math.pi * self.radius_nm**2 / self.period_nm**2


@dataclass
class ModeResult:
    design_id: str
    period_nm: float
    radius_nm: float
    axis_gap_nm: float
    kx_norm: float
    ky_norm: float
    mode_index: int
    freq_thz_real: float
    freq_thz_imag: float
    wavelength_nm: float
    q_value: float
    particle_count: int
    top_count: int
    sio2_count: int
    side_count: int
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


def create_sphere(geom, tag: str, x: float, y: float, z: float, radius: float, label: str) -> None:
    sph = geom.create(tag, "Sphere")
    sph.label(f"{tag}: {label}")
    sph.set("r", fmt(radius))
    sph.set("pos", [fmt(x), fmt(y), fmt(z)])
    sph.set("selresult", "on")
    sph.set("selresultshow", "all")


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


def stable_seed(text: str) -> int:
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def nonoverlap_disks_in_circle(rng: random.Random, radius: float, count: int) -> list[Disk]:
    disks: list[Disk] = []
    attempts = 0
    while len(disks) < count and attempts < 5000:
        attempts += 1
        d = rng.uniform(PARTICLE_D_MIN_NM, PARTICLE_D_MAX_NM)
        rr = d / 2.0
        rho = math.sqrt(rng.random()) * max(1.0, radius - rr - 2.0)
        theta = rng.uniform(0.0, 2.0 * math.pi)
        x = rho * math.cos(theta)
        y = rho * math.sin(theta)
        if all(math.hypot(x - old.x_nm, y - old.y_nm) >= rr + old.radius_nm + 5.0 for old in disks):
            disks.append(Disk(x, y, d))
    return disks


def nonoverlap_disks_in_square_outside_pillar(rng: random.Random, period: float, pillar_radius: float, count: int) -> list[Disk]:
    disks: list[Disk] = []
    half = period / 2.0
    attempts = 0
    while len(disks) < count and attempts < 9000:
        attempts += 1
        d = rng.uniform(PARTICLE_D_MIN_NM, PARTICLE_D_MAX_NM)
        rr = d / 2.0
        x = rng.uniform(-half + rr + 3.0, half - rr - 3.0)
        y = rng.uniform(-half + rr + 3.0, half - rr - 3.0)
        if math.hypot(x, y) < pillar_radius + rr + 3.0:
            continue
        if all(math.hypot(x - old.x_nm, y - old.y_nm) >= rr + old.radius_nm + 5.0 for old in disks):
            disks.append(Disk(x, y, d))
    return disks


def side_particles(rng: random.Random, count: int) -> list[SideParticle]:
    parts: list[SideParticle] = []
    used: list[tuple[float, float, float]] = []
    attempts = 0
    while len(parts) < count and attempts < 5000:
        attempts += 1
        d = rng.uniform(PARTICLE_D_MIN_NM, PARTICLE_D_MAX_NM)
        rr = d / 2.0
        z_local = rng.uniform(rr + 8.0, H_SIN_NM - rr - 8.0)
        theta = rng.uniform(0.0, 2.0 * math.pi)
        ok = True
        for old_theta, old_z, old_d in used:
            arc = abs(math.atan2(math.sin(theta - old_theta), math.cos(theta - old_theta)))
            if arc < 0.40 and abs(z_local - old_z) < 0.5 * (d + old_d) + 12.0:
                ok = False
                break
        if ok:
            parts.append(SideParticle(theta, z_local, d))
            used.append((theta, z_local, d))
    return parts


def particle_layout(design: Design, enabled: bool) -> ParticleLayout:
    if not enabled:
        return ParticleLayout((), (), ())
    rng = random.Random(design.seed or stable_seed(design.design_id))
    pillar_area = math.pi * design.radius_nm**2
    exposed_area = design.period_nm**2 - pillar_area
    top_count = max(4, min(8, int(pillar_area / 9000.0)))
    sio2_count = max(4, min(10, int(exposed_area / 26000.0)))
    side_count = 4
    return ParticleLayout(
        tuple(nonoverlap_disks_in_circle(rng, design.radius_nm, top_count)),
        tuple(nonoverlap_disks_in_square_outside_pillar(rng, design.period_nm, design.radius_nm, sio2_count)),
        tuple(side_particles(rng, side_count)),
    )


def is_fabricable(design: Design) -> bool:
    return design.axis_gap_nm >= MIN_AXIS_GAP_NM and design.radius_nm >= 80.0 and design.period_nm <= 850.0


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
        penalty = 0.0
    elif q_value >= ACCEPT_Q:
        penalty = 0.5 * math.log10(TARGET_Q / max(q_value, 1.0))
    else:
        penalty = 4.0 * math.log10(ACCEPT_Q / max(q_value, 1.0))
    return detune + penalty


def mode_result(design: Design, layout: ParticleLayout, kx_norm: float, ky_norm: float, idx: int, freq_thz: complex) -> ModeResult:
    freq_hz = freq_thz.real * 1e12
    wavelength = C0 / freq_hz * 1e9 if freq_hz > 0 else float("inf")
    q = abs(freq_thz.real / (2.0 * freq_thz.imag)) if abs(freq_thz.imag) > 0 else float("inf")
    return ModeResult(
        design.design_id,
        design.period_nm,
        design.radius_nm,
        design.axis_gap_nm,
        kx_norm,
        ky_norm,
        idx,
        freq_thz.real,
        freq_thz.imag,
        wavelength,
        q,
        layout.count,
        len(layout.pillar_top),
        len(layout.sio2_top),
        len(layout.pillar_side),
        score_mode(wavelength, q),
    )


def get_eigenfrequencies_thz(model) -> list[complex]:
    values = model.evaluate("freq", unit="THz")
    flat = values.ravel() if hasattr(values, "ravel") else values
    modes = unique_complex(flat, tol=1e-7)
    modes = [m for m in modes if math.isfinite(m.real) and abs(m.real) > 1e-9]
    modes.sort(key=lambda z: abs(C0 / (z.real * 1e12) * 1e9 - TARGET_NM) if z.real else float("inf"))
    return modes


def add_field_plot_groups(model, label: str) -> bool:
    try:
        res = model.java.result()
        pg = res.create(f"pg_E2_{label}", "PlotGroup3D")
        pg.label(f"{label}: |E|^2 slices near 980 nm")
        pg.set("data", "dset1")
        slc = pg.feature().create("slc1", "Slice")
        slc.label("|E|^2 horizontal and vertical slices")
        slc.set("expr", "ewfd.normE^2")
        slc.set("quickplane", "xy")
        slc.set("quickz", "z_pillar_mid")
        slc.set("resolution", "normal")
        iso = pg.feature().create("iso1", "Isosurface")
        iso.label("|E|^2 isosurface")
        iso.set("expr", "ewfd.normE^2")
        iso.set("number", "5")
        iso.set("resolution", "normal")
        return True
    except Exception as exc:
        print(f"  plot group setup skipped: {exc}", flush=True)
        return False


def build_model(client, design: Design, name: str, include_particles: bool = False):
    layout = particle_layout(design, include_particles)
    model = client.create(name)
    jmodel = model.java
    half = design.period_nm / 2.0
    z_si = H_PML_NM
    z_sio2 = H_PML_NM + H_SI_BUF_NM
    z_sin = z_sio2 + H_SIO2_NM
    z_top = z_sin + H_SIN_NM
    z_top_pml = z_top + H_AIR_NM
    total_height = 2 * H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_SIN_NM + H_AIR_NM

    set_param(model, "lambda_target", f"{TARGET_NM}[nm]", "Target wavelength")
    set_param(model, "P", f"{design.period_nm:.9g}[nm]", "Square lattice period")
    set_param(model, "R", f"{design.radius_nm:.9g}[nm]", "SiN pillar radius")
    set_param(model, "axis_gap", f"{design.axis_gap_nm:.9g}[nm]", "Axis-aligned gap between SiN pillars")
    set_param(model, "h_sin", "300[nm]", "SiN pillar height")
    set_param(model, "h_re", "16[nm]", "Rare-earth particle reference height")
    set_param(model, "h_sio2_eff", f"{H_SIO2_NM:.9g}[nm]", "Memory-truncated SiO2 buffer")
    set_param(model, "z_sio2_top", f"{z_sin:.9g}[nm]", "Absolute SiO2 top / pillar base z")
    set_param(model, "z_pillar_mid", f"{z_sin + 0.5 * H_SIN_NM:.9g}[nm]", "Absolute pillar mid-plane z")
    set_param(model, "z_pillar_top", f"{z_top:.9g}[nm]", "Absolute pillar top z")
    set_param(model, "n_sin", f"{N_SIN:.9f}", "SiN index at 980 nm")
    set_param(model, "n_sio2", f"{N_SIO2:.9f}", "SiO2 index at 980 nm")
    set_param(model, "n_si", f"{N_SI:.9f}", "Effective Si index at 980 nm")
    set_param(model, "n_re", f"{N_RE:.9f}", "Random particulate RE index")
    set_param(model, "kx", "0[1/m]", "Floquet kx")
    set_param(model, "ky", "0[1/m]", "Floquet ky")

    jmodel.component().create("comp1", True)
    comp = jmodel.component("comp1")
    comp.label("3D square cell SiN cylindrical pillar photonic crystal")
    comp.geom().create("geom1", 3)
    geom = comp.geom("geom1")
    geom.lengthUnit("nm")

    create_block(geom, "si_pml", (-half, -half, 0.0), (design.period_nm, design.period_nm, H_PML_NM), "Si bottom PML")
    create_block(geom, "si_buf", (-half, -half, z_si), (design.period_nm, design.period_nm, H_SI_BUF_NM), "Si buffer")
    create_block(geom, "sio2", (-half, -half, z_sio2), (design.period_nm, design.period_nm, H_SIO2_NM), "SiO2")
    create_block(geom, "air_lower_blk", (-half, -half, z_sin), (design.period_nm, design.period_nm, H_SIN_NM), "air block before pillar subtraction")
    create_cylinder(geom, "sin_pillar", 0.0, 0.0, z_sin, design.radius_nm, H_SIN_NM, "300 nm SiN cylindrical pillar")
    create_difference(geom, "air_lower", ["air_lower_blk"], ["sin_pillar"], "air surrounding SiN pillar after Boolean subtraction")
    create_block(geom, "air_upper", (-half, -half, z_top), (design.period_nm, design.period_nm, H_AIR_NM), "air above pillar")
    create_block(geom, "air_pml", (-half, -half, z_top_pml), (design.period_nm, design.period_nm, H_PML_NM), "air top PML")

    re_tags: list[str] = []
    for i, p in enumerate(layout.pillar_top):
        tag = f"re_top_{i}"
        create_sphere(geom, tag, p.x_nm, p.y_nm, z_top + p.radius_nm, p.radius_nm, "RE particle touching SiN pillar top")
        re_tags.append(f"geom1_{tag}_dom")
    for i, p in enumerate(layout.sio2_top):
        tag = f"re_sio2_{i}"
        create_sphere(geom, tag, p.x_nm, p.y_nm, z_sin + p.radius_nm, p.radius_nm, "RE particle touching exposed SiO2 top")
        re_tags.append(f"geom1_{tag}_dom")
    for i, p in enumerate(layout.pillar_side):
        rr = p.radius_nm
        radial = design.radius_nm + rr
        x = radial * math.cos(p.theta)
        y = radial * math.sin(p.theta)
        tag = f"re_side_{i}"
        create_sphere(geom, tag, x, y, z_sin + p.z_local_nm, rr, "RE particle tangent to SiN pillar sidewall")
        re_tags.append(f"geom1_{tag}_dom")
    geom.run()

    create_union_selection(comp, "air_dom", "3", ["geom1_air_lower_dom", "geom1_air_upper_dom", "geom1_air_pml_dom"])
    create_union_selection(comp, "si_dom", "3", ["geom1_si_pml_dom", "geom1_si_buf_dom"])
    if re_tags:
        create_union_selection(comp, "rare_earth_dom", "3", re_tags)
    create_box_selection(comp, "x_left_bnd", "2", f"{-half-1}", f"{-half+1}", f"{-half-1}", f"{half+1}", "-1", f"{total_height+1}")
    create_box_selection(comp, "x_right_bnd", "2", f"{half-1}", f"{half+1}", f"{-half-1}", f"{half+1}", "-1", f"{total_height+1}")
    create_union_selection(comp, "x_periodic_bnd", "2", ["x_left_bnd", "x_right_bnd"])
    create_box_selection(comp, "y_front_bnd", "2", f"{-half-1}", f"{half+1}", f"{-half-1}", f"{-half+1}", "-1", f"{total_height+1}")
    create_box_selection(comp, "y_back_bnd", "2", f"{-half-1}", f"{half+1}", f"{half-1}", f"{half+1}", "-1", f"{total_height+1}")
    create_union_selection(comp, "y_periodic_bnd", "2", ["y_front_bnd", "y_back_bnd"])

    pml_bot = comp.coordSystem().create("pml_bot", "PML")
    pml_bot.label("Bottom Si PML")
    pml_bot.selection().named("geom1_si_pml_dom")
    pml_top = comp.coordSystem().create("pml_top", "PML")
    pml_top.label("Top air PML")
    pml_top.selection().named("geom1_air_pml_dom")

    create_refractive_material(comp, "mat_air", "Air surrounding pillar", "air_dom", "1", "0")
    create_refractive_material(comp, "mat_sio2", "SiO2 980 nm", "geom1_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_si", "Si effective 980 nm", "si_dom", "n_si", "0")
    create_refractive_material(comp, "mat_sin", "SiN cylindrical pillar 980 nm", "geom1_sin_pillar_dom", "n_sin", "0")
    if re_tags:
        create_refractive_material(comp, "mat_re", "Corrected random particulate RE n=1.52", "rare_earth_dom", "n_re", "0")

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
    mesh.label("Extra-coarse 3D mesh for memory-limited 980 nm BIC scan")
    mesh.autoMeshSize(7)

    study = jmodel.study().create("std1")
    study.label("Eigenfrequency near 980 nm")
    eig = study.create("eig", "Eigenfrequency")
    eig.set("shift", "c_const/lambda_target")
    eig.set("neigs", "12")
    eig.set("neigsmanual", "12")
    eig.set("neigsactive", "on")
    eig.set("eigunit", "THz")
    return model, layout


def solve_design(model, design: Design, layout: ParticleLayout, kx_norm: float = 0.0, ky_norm: float = 0.0) -> list[ModeResult]:
    set_param(model, "kx", f"{kx_norm:.12g}*pi/P")
    set_param(model, "ky", f"{ky_norm:.12g}*pi/P")
    model.java.component("comp1").geom("geom1").run()
    model.java.component("comp1").mesh("mesh1").run()
    model.java.study("std1").run()
    modes = get_eigenfrequencies_thz(model)
    return [mode_result(design, layout, kx_norm, ky_norm, i + 1, freq) for i, freq in enumerate(modes)]


def best_mode(rows: list[ModeResult]) -> ModeResult:
    return sorted(rows, key=lambda r: (r.score, abs(r.wavelength_nm - TARGET_NM), -r.q_value))[0]


def write_modes_csv(path: Path, rows: list[ModeResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "design_id", "period_nm", "radius_nm", "axis_gap_nm", "kx_norm", "ky_norm", "mode_index",
            "freq_thz_real", "freq_thz_imag", "wavelength_nm", "q_value", "particle_count",
            "top_count", "sio2_count", "side_count", "score",
        ])
        for r in rows:
            writer.writerow([
                r.design_id, r.period_nm, r.radius_nm, r.axis_gap_nm, r.kx_norm, r.ky_norm, r.mode_index,
                r.freq_thz_real, r.freq_thz_imag, r.wavelength_nm, r.q_value, r.particle_count,
                r.top_count, r.sio2_count, r.side_count, r.score,
            ])


def read_modes_csv(path: Path) -> list[ModeResult]:
    if not path.exists():
        return []
    rows: list[ModeResult] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            rows.append(ModeResult(
                raw["design_id"], float(raw["period_nm"]), float(raw["radius_nm"]), float(raw["axis_gap_nm"]),
                float(raw["kx_norm"]), float(raw["ky_norm"]), int(raw["mode_index"]), float(raw["freq_thz_real"]),
                float(raw["freq_thz_imag"]), float(raw["wavelength_nm"]), float(raw["q_value"]),
                int(float(raw["particle_count"])), int(float(raw["top_count"])), int(float(raw["sio2_count"])),
                int(float(raw["side_count"])), float(raw["score"]),
            ))
    return rows


def design_from_result(row: ModeResult) -> Design:
    return Design(row.design_id, row.period_nm, row.radius_nm, stable_seed(row.design_id))


def coarse_designs() -> list[Design]:
    out: list[Design] = []
    # High-Q priority pass: previous v7 scans showed Q increasing as the
    # inter-pillar gap shrank. This final search moves the gap=10-30 nm
    # branch toward 980 nm by reducing the period.
    for period in [572.0, 574.0, 576.0, 578.0, 580.0, 582.0, 584.0]:
        for gap in [10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]:
            radius = (period - gap) / 2.0
            d = Design(f"targetgap_P{period:.0f}_G{gap:.0f}", period, radius, stable_seed(f"targetgap_{period}_{gap}"))
            if is_fabricable(d):
                out.append(d)
    return out


def fine_designs(best: ModeResult) -> list[Design]:
    out: list[Design] = []
    seen: set[tuple[float, float]] = set()
    for dp in [-20.0, -10.0, 0.0, 10.0, 20.0]:
        for dr in [-10.0, -5.0, 0.0, 5.0, 10.0]:
            period = round(best.period_nm + dp, 3)
            radius = round(best.radius_nm + dr, 3)
            key = (period, radius)
            if key in seen:
                continue
            seen.add(key)
            d = Design(f"fine_P{period:.0f}_R{radius:.0f}", period, radius, stable_seed(f"fine_{period}_{radius}"))
            if is_fabricable(d):
                out.append(d)
    return out


def particle_designs(best: ModeResult) -> list[Design]:
    out: list[Design] = []
    for dp, dr in [(0, 0), (-10, 0), (10, 0), (0, -5), (0, 5), (-10, -5), (10, 5)]:
        period = round(best.period_nm + dp, 3)
        radius = round(best.radius_nm + dr, 3)
        d = Design(f"particle_P{period:.0f}_R{radius:.0f}", period, radius, stable_seed(f"particle_{period}_{radius}"))
        if is_fabricable(d):
            out.append(d)
    return out


def scan_designs(client, designs: list[Design], path: Path, label: str, include_particles: bool = False) -> list[ModeResult]:
    rows = read_modes_csv(path)
    done = {r.design_id for r in rows}
    for i, design in enumerate(designs, 1):
        if design.design_id in done:
            continue
        print(f"[{label}] {i}/{len(designs)} P={design.period_nm:.3f} R={design.radius_nm:.3f} gap={design.axis_gap_nm:.3f} particles={include_particles}", flush=True)
        model, layout = build_model(client, design, f"v7_{design.design_id}", include_particles=include_particles)
        try:
            current = solve_design(model, design, layout, 0.0, 0.0)
            rows.extend(current)
            write_modes_csv(path, rows)
            best = best_mode(current)
            print(f"  best lambda={best.wavelength_nm:.6f} nm Q={best.q_value:.4g} mode={best.mode_index} particles={best.particle_count}", flush=True)
        except Exception as exc:
            print(f"  FAILED {design.design_id}: {exc}", flush=True)
        finally:
            client.remove(model.name())
    return rows


def choose_candidates(rows: list[ModeResult], limit: int = 10) -> list[ModeResult]:
    best_by_design = []
    for design_id in sorted({r.design_id for r in rows}):
        best_by_design.append(best_mode([r for r in rows if r.design_id == design_id]))
    ranked = sorted(best_by_design, key=lambda r: (r.score, -r.q_value))
    high_q = sorted(best_by_design, key=lambda r: -r.q_value)[:3]
    merged: dict[str, ModeResult] = {}
    for r in ranked[:limit] + high_q:
        merged[f"{r.design_id}:{r.mode_index}"] = r
    return sorted(merged.values(), key=lambda r: (r.score, -r.q_value))[:limit]


def save_model(client, best: ModeResult, path: Path, name: str, include_particles: bool) -> ModeResult:
    design = design_from_result(best)
    model, layout = build_model(client, design, name, include_particles=include_particles)
    try:
        rows = solve_design(model, design, layout, 0.0, 0.0)
        selected = min(rows, key=lambda r: abs(r.wavelength_nm - best.wavelength_nm) + 1e-6 * abs(math.log10(max(r.q_value, 1.0)) - math.log10(max(best.q_value, 1.0))))
        add_field_plot_groups(model, "particle" if include_particles else "clean")
        model.save(str(path))
        if include_particles:
            write_particle_geometry_csv(design, layout)
        return selected
    finally:
        client.remove(model.name())


def run_band(client, best: ModeResult) -> list[ModeResult]:
    design = design_from_result(best)
    model, layout = build_model(client, design, "sin980_qbic_v7_clean_band", include_particles=False)
    rows: list[ModeResult] = []
    path = [
        (0.0, 0.0), (0.02, 0.0), (0.05, 0.0), (0.1, 0.0), (0.2, 0.0), (0.35, 0.0), (0.5, 0.0),
        (0.5, 0.25), (0.5, 0.5), (0.35, 0.35), (0.2, 0.2), (0.1, 0.1), (0.05, 0.05), (0.0, 0.0),
    ]
    try:
        for kx, ky in path:
            print(f"[band] k=({kx:.3g},{ky:.3g}) pi/P", flush=True)
            current = solve_design(model, design, layout, kx, ky)
            rows.append(best_mode(current))
        write_modes_csv(OUT / "band_clean_980.csv", rows)
        q_near = [r for r in rows if r.kx_norm <= 0.2 and abs(r.ky_norm) < 1e-12]
        write_modes_csv(OUT / "q_vs_k_clean_980.csv", q_near)
    finally:
        client.remove(model.name())
    return rows


def nearest_gap(x: float, y: float, disks: list[tuple[float, float, float]]) -> float:
    gaps = []
    for ox, oy, orad in disks:
        gaps.append(math.hypot(x - ox, y - oy) - orad)
    return min(gaps) if gaps else float("nan")


def write_particle_geometry_csv(design: Design, layout: ParticleLayout) -> None:
    rows = []
    placed: list[tuple[float, float, float]] = []
    z_si = H_PML_NM
    z_sio2 = H_PML_NM + H_SI_BUF_NM
    z_sin = z_sio2 + H_SIO2_NM
    z_top = z_sin + H_SIN_NM

    def add_row(region: str, x: float, y: float, z: float, diameter: float, valid: bool, contact: str) -> None:
        rr = diameter / 2.0
        rows.append([region, x, y, z, diameter, contact, nearest_gap(x, y, placed), valid])
        placed.append((x, y, rr))

    for p in layout.pillar_top:
        add_row("pillar_top", p.x_nm, p.y_nm, z_top + p.radius_nm, p.diameter_nm, math.hypot(p.x_nm, p.y_nm) + p.radius_nm <= design.radius_nm + 1e-6, "SiN top")
    for p in layout.sio2_top:
        add_row("sio2_top", p.x_nm, p.y_nm, z_sin + p.radius_nm, p.diameter_nm, math.hypot(p.x_nm, p.y_nm) >= design.radius_nm + p.radius_nm - 1e-6, "SiO2 top")
    for p in layout.pillar_side:
        rr = p.radius_nm
        x = (design.radius_nm + rr) * math.cos(p.theta)
        y = (design.radius_nm + rr) * math.sin(p.theta)
        z = z_sin + p.z_local_nm
        valid = z_sin + rr <= z <= z_top - rr and abs(math.hypot(x, y) - (design.radius_nm + rr)) < 1e-6
        add_row("pillar_side", x, y, z, p.diameter_nm, valid, "SiN sidewall")

    with (OUT / "particle_geometry_check.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["region", "x_nm", "y_nm", "z_nm", "diameter_nm", "contact_surface", "nearest_prior_center_gap_nm", "valid"])
        writer.writerows(rows)


def svg_scatter(rows: list[ModeResult], path: Path, title: str) -> None:
    pts = [(r.wavelength_nm, r.q_value, r.design_id) for r in rows if math.isfinite(r.wavelength_nm) and math.isfinite(r.q_value)]
    if not pts:
        return
    width, height = 900, 540
    left, right, top, bottom = 80, 30, 40, 70
    pw, ph = width - left - right, height - top - bottom
    xmin, xmax = min(x for x, _, _ in pts) - 5.0, max(x for x, _, _ in pts) + 5.0
    qmin, qmax = max(1.0, min(q for _, q, _ in pts) / 1.5), max(q for _, q, _ in pts) * 1.5

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * pw

    def sy(q: float) -> float:
        return top + (math.log10(qmax) - math.log10(max(q, qmin))) / (math.log10(qmax) - math.log10(qmin)) * ph

    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    lines.append(f'<text x="{width/2}" y="24" text-anchor="middle" font-family="Arial" font-size="18">{title}</text>')
    lines.append(f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#111"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#111"/>')
    for xv in [940, 960, 980, 1000, 1020]:
        if xmin <= xv <= xmax:
            color = "#111" if xv == TARGET_NM else "#ddd"
            dash = ' stroke-dasharray="5,4"' if xv == TARGET_NM else ""
            lines.append(f'<line x1="{sx(xv):.2f}" y1="{top}" x2="{sx(xv):.2f}" y2="{top+ph}" stroke="{color}"{dash}/>')
            lines.append(f'<text x="{sx(xv):.2f}" y="{top+ph+20}" text-anchor="middle" font-family="Arial" font-size="12">{xv}</text>')
    for decade in range(int(math.floor(math.log10(qmin))), int(math.ceil(math.log10(qmax))) + 1):
        qv = 10**decade
        if qmin <= qv <= qmax:
            lines.append(f'<line x1="{left}" y1="{sy(qv):.2f}" x2="{left+pw}" y2="{sy(qv):.2f}" stroke="#eee"/>')
            lines.append(f'<text x="{left-8}" y="{sy(qv)+4:.2f}" text-anchor="end" font-family="Arial" font-size="12">1e{decade}</text>')
    for x, q, did in pts:
        color = "#1f77b4" if not did.startswith("particle") else "#d62728"
        lines.append(f'<circle cx="{sx(x):.2f}" cy="{sy(q):.2f}" r="3.5" fill="{color}" opacity="0.75"><title>{did} lambda={x:.6g} Q={q:.6g}</title></circle>')
    lines.append(f'<text x="{left+pw/2}" y="{height-22}" text-anchor="middle" font-family="Arial" font-size="14">Wavelength (nm)</text>')
    lines.append(f'<text x="22" y="{top+ph/2}" transform="rotate(-90 22 {top+ph/2})" text-anchor="middle" font-family="Arial" font-size="14">Q (log)</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def svg_band(rows: list[ModeResult], path: Path) -> None:
    if not rows:
        return
    width, height = 920, 540
    left, right, top, bottom = 80, 90, 40, 90
    pw, ph = width - left - right, height - top - bottom
    labels = [f"({r.kx_norm:g},{r.ky_norm:g})" for r in rows]
    xcoords = [left + i * pw / max(1, len(rows) - 1) for i in range(len(rows))]
    wl = [r.wavelength_nm for r in rows]
    q = [r.q_value for r in rows]
    wmin, wmax = min(wl) - 5.0, max(wl) + 5.0
    qmin, qmax = max(1.0, min(q) / 1.3), max(q) * 1.3

    def sy_w(y: float) -> float:
        return top + (wmax - y) / (wmax - wmin) * ph

    def sy_q(y: float) -> float:
        return top + (math.log10(qmax) - math.log10(max(y, qmin))) / (math.log10(qmax) - math.log10(qmin)) * ph

    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    lines.append(f'<text x="{width/2}" y="24" text-anchor="middle" font-family="Arial" font-size="18">v7 clean band and Q(k)</text>')
    lines.append(f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#111"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#111"/>')
    lines.append(f'<line x1="{left+pw}" y1="{top}" x2="{left+pw}" y2="{top+ph}" stroke="#111"/>')
    wpts = " ".join(f"{xcoords[i]:.2f},{sy_w(wl[i]):.2f}" for i in range(len(rows)))
    qpts = " ".join(f"{xcoords[i]:.2f},{sy_q(q[i]):.2f}" for i in range(len(rows)))
    lines.append(f'<polyline fill="none" stroke="#1f77b4" stroke-width="2" points="{wpts}"/>')
    lines.append(f'<polyline fill="none" stroke="#d62728" stroke-width="2" points="{qpts}"/>')
    for i, label in enumerate(labels):
        lines.append(f'<text x="{xcoords[i]:.2f}" y="{top+ph+22}" text-anchor="end" transform="rotate(-35 {xcoords[i]:.2f} {top+ph+22})" font-family="Arial" font-size="10">{label}</text>')
    lines.append(f'<text x="22" y="{top+ph/2}" transform="rotate(-90 22 {top+ph/2})" text-anchor="middle" font-family="Arial" font-size="14" fill="#1f77b4">Wavelength (nm)</text>')
    lines.append(f'<text x="{width-24}" y="{top+ph/2}" transform="rotate(90 {width-24} {top+ph/2})" text-anchor="middle" font-family="Arial" font-size="14" fill="#d62728">Q (log)</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_report(clean_rows: list[ModeResult], fine_rows: list[ModeResult], particle_rows: list[ModeResult], clean_best: ModeResult, particle_best: ModeResult | None, band: list[ModeResult]) -> None:
    all_clean = clean_rows + fine_rows
    strict = [r for r in all_clean if abs(r.wavelength_nm - TARGET_NM) <= 0.5 and r.q_value >= TARGET_Q]
    relaxed = [r for r in all_clean if abs(r.wavelength_nm - TARGET_NM) <= 1.0 and r.q_value >= ACCEPT_Q]
    highest_clean = max(all_clean, key=lambda r: r.q_value)
    near_clean = [r for r in all_clean if abs(r.wavelength_nm - TARGET_NM) <= 1.0]
    best_near = max(near_clean, key=lambda r: r.q_value) if near_clean else None
    particle_text = "not run"
    if particle_best:
        particle_text = f"lambda `{particle_best.wavelength_nm:.6g} nm`, Q `{particle_best.q_value:.6g}`"
    lines = [
        "# v7 980 nm SiN cylindrical-pillar BIC with corrected random particles",
        "",
        f"All files were generated inside `{ROOT}`; previous versions and source .mph files were not edited.",
        "",
        "## Model",
        "",
        "- Geometry: 3D square-lattice SiN cylindrical pillars on SiO2, not cylindrical holes.",
        f"- Materials: n_sin `{N_SIN:.9f}`, n_sio2 `{N_SIO2:.3f}`, n_si `{N_SI:.3f}`, n_re `{N_RE:.3f}`.",
        f"- Stack is memory-reduced: SiO2 `{H_SIO2_NM:.0f} nm`, Si buffer `{H_SI_BUF_NM:.0f} nm`, PML `{H_PML_NM:.0f} nm`.",
        "- Corrected particles are only on exposed SiO2 top, SiN pillar top, and SiN pillar sidewall.",
        "",
        "## Scan result",
        "",
        f"- Clean rows simulated: `{len(all_clean)}`.",
        f"- Particle rows simulated: `{len(particle_rows)}`.",
        f"- Strict clean target `|lambda-980| <= 0.5 nm`, `Q >= 1e9`: `{len(strict)}`.",
        f"- Relaxed clean target `|lambda-980| <= 1.0 nm`, `Q >= 1e8`: `{len(relaxed)}`.",
        "",
        "## Best clean candidate",
        "",
        "| item | value |",
        "|---|---:|",
        f"| P (nm) | `{clean_best.period_nm:.9g}` |",
        f"| R (nm) | `{clean_best.radius_nm:.9g}` |",
        f"| axis gap (nm) | `{clean_best.axis_gap_nm:.9g}` |",
        f"| lambda980 (nm) | `{clean_best.wavelength_nm:.9g}` |",
        f"| Q980 | `{clean_best.q_value:.6g}` |",
        "",
        "## Diagnostics",
        "",
        f"- Highest-Q clean mode: `{highest_clean.design_id}`, lambda `{highest_clean.wavelength_nm:.6g} nm`, Q `{highest_clean.q_value:.6g}`.",
        f"- Best clean mode within `980+/-1 nm`: `{best_near.design_id if best_near else 'none'}`" + (f", lambda `{best_near.wavelength_nm:.6g} nm`, Q `{best_near.q_value:.6g}`." if best_near else "."),
        f"- Best particle-perturbed result: {particle_text}.",
    ]
    if clean_best.q_value > 1e8:
        lines.append("- Compared with the 451 nm 3D random-particle run, the 980 nm cylindrical-pillar branch is substantially more favorable.")
    else:
        lines.append("- In this first pillar-array parameter family, 980 nm did not automatically recover ultra-high Q; geometry family and particle perturbation still matter.")
    if particle_best and clean_best.q_value > 0:
        lines.append(f"- Particle perturbation Q ratio: `{particle_best.q_value / clean_best.q_value:.6g}`.")
    if band:
        lines.extend(["", "## Band / Q(k)", ""])
        for r in band:
            lines.append(f"- k=({r.kx_norm:.3g}, {r.ky_norm:.3g}) pi/P: lambda `{r.wavelength_nm:.6g} nm`, Q `{r.q_value:.6g}`.")
    lines.extend([
        "",
        "## Output files",
        "",
        "- `sin980_qbic_v7_clean_best.mph`",
        "- `sin980_qbic_v7_particle_perturbed_best.mph`",
        "- `clean_coarse_scan.csv`",
        "- `clean_fine_scan.csv`",
        "- `particle_perturbation_scan.csv`",
        "- `band_clean_980.csv`",
        "- `q_vs_k_clean_980.csv`",
        "- `particle_geometry_check.csv`",
        "- `scan_q_vs_wavelength.svg`",
        "- `band_q_wavelength.svg`",
    ])
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def run_all(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = import_mph()
    client = mph.Client(cores=args.cores)
    try:
        coarse = coarse_designs()
        if args.limit:
            coarse = coarse[: args.limit]
        coarse_rows = scan_designs(client, coarse, OUT / "clean_coarse_scan.csv", "clean-coarse", include_particles=False)
        coarse_candidates = choose_candidates(coarse_rows, 8)
        write_modes_csv(OUT / "clean_coarse_best_candidates.csv", coarse_candidates)
        if not coarse_candidates:
            (OUT / "summary.md").write_text("# v7 980 nm pillar scan\n\nNo clean coarse candidates solved successfully.", encoding="utf-8")
            return
        fine = fine_designs(coarse_candidates[0])
        fine_rows = scan_designs(client, fine, OUT / "clean_fine_scan.csv", "clean-fine", include_particles=False)
        all_clean = coarse_rows + fine_rows
        clean_candidates = choose_candidates(all_clean, 10)
        write_modes_csv(OUT / "clean_best_candidates.csv", clean_candidates)
        clean_best = clean_candidates[0]
        clean_best = save_model(client, clean_best, OUT / "sin980_qbic_v7_clean_best.mph", "sin980_qbic_v7_clean_best", include_particles=False)
        particle_rows: list[ModeResult] = []
        particle_best: ModeResult | None = None
        if not args.skip_particles:
            particle_rows = scan_designs(client, particle_designs(clean_best), OUT / "particle_perturbation_scan.csv", "particle-perturb", include_particles=True)
            if particle_rows:
                particle_candidates = choose_candidates(particle_rows, 8)
                write_modes_csv(OUT / "particle_best_candidates.csv", particle_candidates)
                particle_best = save_model(client, particle_candidates[0], OUT / "sin980_qbic_v7_particle_perturbed_best.mph", "sin980_qbic_v7_particle_best", include_particles=True)
        band: list[ModeResult] = []
        if not args.skip_band:
            band = run_band(client, clean_best)
        svg_scatter(all_clean + particle_rows, OUT / "scan_q_vs_wavelength.svg", "v7 980 nm pillar scan")
        svg_band(band, OUT / "band_q_wavelength.svg")
        write_report(coarse_rows, fine_rows, particle_rows, clean_best, particle_best, band)
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="v7 980 nm SiN cylindrical-pillar BIC with corrected random particles.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-all")
    p.add_argument("--cores", type=int, default=2)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--skip-band", action="store_true")
    p.add_argument("--skip-particles", action="store_true")
    p.set_defaults(func=run_all)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
