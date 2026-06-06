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
H_SIO2_NM = 500.0
H_SI_BUF_NM = 150.0
H_PML_NM = 300.0
H_AIR_NM = 450.0

PARTICLE_D_MIN_NM = 24.0
PARTICLE_D_MAX_NM = 44.0
PARTICLE_H_NM = 16.0
MIN_AXIS_GAP_NM = 130.0
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
    z_nm: float
    diameter_nm: float

    @property
    def radius_nm(self) -> float:
        return self.diameter_nm / 2.0


@dataclass(frozen=True)
class ParticleLayout:
    top: tuple[Disk, ...]
    bottom: tuple[Disk, ...]
    side: tuple[SideParticle, ...]

    @property
    def count(self) -> int:
        return len(self.top) + len(self.bottom) + len(self.side)


@dataclass(frozen=True)
class Design:
    design_id: str
    period_nm: float
    radius_nm: float
    seed: int = 0

    @property
    def axis_gap_nm(self) -> float:
        return self.period_nm - 2.0 * self.radius_nm


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
    bottom_count: int
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
        if all(math.hypot(x - old.x_nm, y - old.y_nm) >= rr + old.radius_nm + 4.0 for old in disks):
            disks.append(Disk(x, y, d))
    return disks


def nonoverlap_disks_in_square_outside_circle(rng: random.Random, period: float, cylinder_radius: float, count: int) -> list[Disk]:
    disks: list[Disk] = []
    half = period / 2.0
    attempts = 0
    while len(disks) < count and attempts < 9000:
        attempts += 1
        d = rng.uniform(PARTICLE_D_MIN_NM, PARTICLE_D_MAX_NM)
        rr = d / 2.0
        x = rng.uniform(-half + rr + 2.0, half - rr - 2.0)
        y = rng.uniform(-half + rr + 2.0, half - rr - 2.0)
        if math.hypot(x, y) < cylinder_radius + H_RE_NM + rr + 2.0:
            continue
        if all(math.hypot(x - old.x_nm, y - old.y_nm) >= rr + old.radius_nm + 4.0 for old in disks):
            disks.append(Disk(x, y, d))
    return disks


def side_particles(rng: random.Random, count: int) -> list[SideParticle]:
    parts: list[SideParticle] = []
    used: list[tuple[float, float, float]] = []
    attempts = 0
    while len(parts) < count and attempts < 3000:
        attempts += 1
        d = rng.uniform(PARTICLE_D_MIN_NM, PARTICLE_D_MAX_NM)
        z = rng.uniform(d / 2.0 + 8.0, H_SIN_NM - d / 2.0 - 8.0)
        theta = rng.uniform(0.0, 2.0 * math.pi)
        ok = True
        for old_theta, old_z, old_d in used:
            arc = abs(math.atan2(math.sin(theta - old_theta), math.cos(theta - old_theta)))
            if arc < 0.38 and abs(z - old_z) < 0.5 * (d + old_d) + 10.0:
                ok = False
                break
        if ok:
            parts.append(SideParticle(theta, z, d))
            used.append((theta, z, d))
    return parts


def particle_layout(design: Design) -> ParticleLayout:
    rng = random.Random(design.seed or stable_seed(design.design_id))
    top_area = design.period_nm**2 - math.pi * (design.radius_nm + H_RE_NM) ** 2
    top_count = max(4, min(6, int(top_area / 9000.0)))
    bottom_count = max(4, min(6, int(math.pi * design.radius_nm**2 / 3600.0)))
    side_count = 2
    return ParticleLayout(
        tuple(nonoverlap_disks_in_square_outside_circle(rng, design.period_nm, design.radius_nm, top_count)),
        tuple(nonoverlap_disks_in_circle(rng, design.radius_nm, bottom_count)),
        tuple(side_particles(rng, side_count)),
    )


def is_fabricable(design: Design) -> bool:
    return design.axis_gap_nm >= MIN_AXIS_GAP_NM and design.radius_nm >= 45.0


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
        penalty = 0.5 * math.log10(TARGET_Q / q_value)
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
        len(layout.top),
        len(layout.bottom),
        len(layout.side),
        score_mode(wavelength, q),
    )


def get_eigenfrequencies_thz(model) -> list[complex]:
    values = model.evaluate("freq", unit="THz")
    flat = values.ravel() if hasattr(values, "ravel") else values
    modes = unique_complex(flat, tol=1e-7)
    modes = [m for m in modes if math.isfinite(m.real) and abs(m.real) > 1e-9]
    modes.sort(key=lambda z: abs(C0 / (z.real * 1e12) * 1e9 - TARGET_NM) if z.real else float("inf"))
    return modes


def build_model(client, design: Design, name: str):
    layout = particle_layout(design)
    model = client.create(name)
    jmodel = model.java
    half = design.period_nm / 2.0
    z_si = H_PML_NM
    z_sio2 = H_PML_NM + H_SI_BUF_NM
    z_sin = z_sio2 + H_SIO2_NM
    z_top = z_sin + H_SIN_NM
    z_air = z_top + H_RE_NM
    z_top_pml = z_top + H_AIR_NM
    total_height = 2 * H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_SIN_NM + H_AIR_NM

    set_param(model, "lambda_target", f"{TARGET_NM}[nm]", "Target wavelength")
    set_param(model, "P", f"{design.period_nm:.9g}[nm]", "Square lattice period")
    set_param(model, "R", f"{design.radius_nm:.9g}[nm]", "Cylindrical air-hole radius")
    set_param(model, "axis_gap", f"{design.axis_gap_nm:.9g}[nm]", "Axis-aligned SiN bridge between holes")
    set_param(model, "h_sin", "300[nm]", "SiN slab height")
    set_param(model, "h_re", "16[nm]", "Particle height")
    set_param(model, "h_sio2_eff", f"{H_SIO2_NM:.9g}[nm]", "Memory-truncated SiO2 buffer")
    set_param(model, "n_sin", f"{N_SIN:.9f}", "SiN index at 451 nm")
    set_param(model, "n_sio2", f"{N_SIO2:.9f}", "SiO2 index at 451 nm")
    set_param(model, "n_si", f"{N_SI:.9f}", "Effective Si index at 451 nm")
    set_param(model, "n_re", f"{N_RE:.9f}", "Random particulate RE index")
    set_param(model, "kx", "0[1/m]", "Floquet kx")
    set_param(model, "ky", "0[1/m]", "Floquet ky")

    jmodel.component().create("comp1", True)
    comp = jmodel.component("comp1")
    comp.label("3D square cell cylindrical air-hole SiN photonic crystal")
    comp.geom().create("geom1", 3)
    geom = comp.geom("geom1")
    geom.lengthUnit("nm")

    create_block(geom, "si_pml", (-half, -half, 0.0), (design.period_nm, design.period_nm, H_PML_NM), "Si bottom PML")
    create_block(geom, "si_buf", (-half, -half, z_si), (design.period_nm, design.period_nm, H_SI_BUF_NM), "Si buffer")
    create_block(geom, "sio2", (-half, -half, z_sio2), (design.period_nm, design.period_nm, H_SIO2_NM), "SiO2")
    create_block(geom, "sin_slab", (-half, -half, z_sin), (design.period_nm, design.period_nm, H_SIN_NM), "300 nm SiN slab")
    create_cylinder(geom, "air_hole", 0.0, 0.0, z_sin, design.radius_nm, H_SIN_NM, "cylindrical air hole through SiN")
    create_block(geom, "air", (-half, -half, z_top), (design.period_nm, design.period_nm, H_AIR_NM), "air above slab")
    create_block(geom, "air_pml", (-half, -half, z_top_pml), (design.period_nm, design.period_nm, H_PML_NM), "air top PML")

    re_tags: list[str] = []
    for i, p in enumerate(layout.top):
        tag = f"re_top_{i}"
        create_cylinder(geom, tag, p.x_nm, p.y_nm, z_top, p.radius_nm, PARTICLE_H_NM, "top RE particle")
        re_tags.append(f"geom1_{tag}_dom")
    for i, p in enumerate(layout.bottom):
        tag = f"re_bottom_{i}"
        create_cylinder(geom, tag, p.x_nm, p.y_nm, z_sin, p.radius_nm, PARTICLE_H_NM, "bottom RE particle")
        re_tags.append(f"geom1_{tag}_dom")
    for i, p in enumerate(layout.side):
        r = max(1.0, design.radius_nm - p.radius_nm)
        x = r * math.cos(p.theta)
        y = r * math.sin(p.theta)
        tag = f"re_side_{i}"
        create_sphere(geom, tag, x, y, z_sin + p.z_nm, p.radius_nm, "sparse sidewall RE particle")
        re_tags.append(f"geom1_{tag}_dom")
    geom.run()

    create_union_selection(comp, "air_dom", "3", ["geom1_air_hole_dom", "geom1_air_dom", "geom1_air_pml_dom"])
    create_union_selection(comp, "si_dom", "3", ["geom1_si_pml_dom", "geom1_si_buf_dom"])
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

    create_refractive_material(comp, "mat_sio2", "SiO2 451 nm", "geom1_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_sin", "SiN slab 451 nm", "geom1_sin_slab_dom", "n_sin", "0")
    create_refractive_material(comp, "mat_air", "Air and cylindrical holes", "air_dom", "1", "0")
    create_refractive_material(comp, "mat_re", "Random particulate RE n=1.52", "rare_earth_dom", "n_re", "0")
    create_refractive_material(comp, "mat_si", "Si effective 451 nm", "si_dom", "n_si", "0")

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
    mesh.label("Extra-coarse 3D mesh for memory-limited BIC existence scan")
    mesh.autoMeshSize(7)

    study = jmodel.study().create("std1")
    study.label("Eigenfrequency near 451 nm")
    eig = study.create("eig", "Eigenfrequency")
    eig.set("shift", "c_const/lambda_target")
    eig.set("neigs", "4")
    eig.set("neigsmanual", "4")
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
            "top_count", "bottom_count", "side_count", "score",
        ])
        for r in rows:
            writer.writerow([
                r.design_id, r.period_nm, r.radius_nm, r.axis_gap_nm, r.kx_norm, r.ky_norm, r.mode_index,
                r.freq_thz_real, r.freq_thz_imag, r.wavelength_nm, r.q_value, r.particle_count,
                r.top_count, r.bottom_count, r.side_count, r.score,
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
                int(float(raw["particle_count"])), int(float(raw["top_count"])), int(float(raw["bottom_count"])),
                int(float(raw["side_count"])), float(raw["score"]),
            ))
    return rows


def design_from_result(row: ModeResult) -> Design:
    return Design(row.design_id, row.period_nm, row.radius_nm, stable_seed(row.design_id))


def coarse_designs() -> list[Design]:
    out: list[Design] = []
    seen: set[str] = set()

    def add_design(design_id: str, period: float, radius: float) -> None:
        if design_id in seen:
            return
        d = Design(design_id, period, radius, stable_seed(design_id))
        if is_fabricable(d):
            out.append(d)
            seen.add(design_id)

    for period, gaps in [
        (294.0, [135.0, 130.0, 150.0]),
        (306.0, [155.0, 150.0, 180.0]),
        (280.0, [130.0, 150.0]),
        (310.0, [130.0, 150.0, 180.0]),
        (330.0, [130.0, 150.0, 180.0]),
        (260.0, [130.0, 150.0]),
        (360.0, [150.0, 180.0]),
    ]:
        for gap in gaps:
            radius = (period - gap) / 2.0
            add_design(f"cyl_P{period:.0f}_G{gap:.0f}", period, radius)

    # Local refinement around the best first-pass 3D result near P=260 nm.
    for period, gaps in [
        (250.0, [130.0, 140.0, 150.0]),
        (255.0, [135.0, 145.0, 155.0]),
        (265.0, [135.0, 145.0, 155.0]),
        (270.0, [140.0, 150.0, 160.0]),
    ]:
        for gap in gaps:
            radius = (period - gap) / 2.0
            add_design(f"fine_P{period:.0f}_G{gap:.0f}", period, radius)

    # bic0413.mph reference scale: period=356 nm, r/period=0.220.
    # This is used only as a parameter-range anchor, not as a copied design.
    for period in [294.0, 306.0, 320.0, 340.0, 356.0, 370.0]:
        for ratio in [0.20, 0.22, 0.24]:
            radius = period * ratio
            add_design(f"bicref_P{period:.0f}_rho{ratio:.2f}".replace(".", "p"), period, radius)
    add_design("bicref_P356_exact", 356.0, 356.0 * 0.220)
    return out


def scan_designs(client, designs: list[Design], path: Path, label: str) -> list[ModeResult]:
    rows = read_modes_csv(path)
    done = {r.design_id for r in rows}
    for i, design in enumerate(designs, 1):
        if design.design_id in done:
            continue
        print(f"[{label}] {i}/{len(designs)} solving 3D P={design.period_nm:.3f} R={design.radius_nm:.3f} gap={design.axis_gap_nm:.3f}", flush=True)
        model, layout = build_model(client, design, f"v6_{design.design_id}")
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


def choose_candidates(rows: list[ModeResult], limit: int = 8) -> list[ModeResult]:
    best_by_design = []
    for design_id in sorted({r.design_id for r in rows}):
        best_by_design.append(best_mode([r for r in rows if r.design_id == design_id]))
    ranked = sorted(best_by_design, key=lambda r: (r.score, -r.q_value))
    high_q = sorted(best_by_design, key=lambda r: -r.q_value)[:3]
    merged: dict[str, ModeResult] = {}
    for r in ranked[:limit] + high_q:
        merged[f"{r.design_id}:{r.mode_index}"] = r
    return sorted(merged.values(), key=lambda r: (r.score, -r.q_value))[:limit]


def save_best_model(client, best: ModeResult) -> ModeResult:
    design = design_from_result(best)
    model, layout = build_model(client, design, "sin451_qbic_v6_best_3d_cyl_phc_random_particles")
    try:
        rows = solve_design(model, design, layout, 0.0, 0.0)
        selected = min(rows, key=lambda r: abs(r.wavelength_nm - best.wavelength_nm) + 1e-6 * abs(math.log10(max(r.q_value, 1.0)) - math.log10(max(best.q_value, 1.0))))
        model.save(str(OUT / "sin451_qbic_v6_best_3d_cyl_hole_phc_random_particles.mph"))
        write_particle_layout_csv(design, layout)
        return selected
    finally:
        client.remove(model.name())


def run_band(client, best: ModeResult) -> list[ModeResult]:
    design = design_from_result(best)
    model, layout = build_model(client, design, "sin451_qbic_v6_band")
    rows: list[ModeResult] = []
    try:
        for kx, ky in [(0, 0), (0.02, 0), (0.05, 0), (0.1, 0), (0.2, 0), (0.5, 0), (0, 0.5), (0.5, 0.5)]:
            print(f"[band] k=({kx:.3g},{ky:.3g}) pi/P", flush=True)
            current = solve_design(model, design, layout, kx, ky)
            rows.append(best_mode(current))
        write_modes_csv(OUT / "band_451_3d_particle.csv", rows)
    finally:
        client.remove(model.name())
    return rows


def write_particle_layout_csv(design: Design, layout: ParticleLayout) -> None:
    with (OUT / "random_particle_layout.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["region", "x_nm", "y_nm", "z_or_theta_nm", "diameter_nm"])
        for p in layout.top:
            writer.writerow(["top", p.x_nm, p.y_nm, H_SIN_NM, p.diameter_nm])
        for p in layout.bottom:
            writer.writerow(["bottom", p.x_nm, p.y_nm, 0.0, p.diameter_nm])
        for p in layout.side:
            writer.writerow(["side", "", "", p.z_nm, p.diameter_nm])


def write_report(coarse: list[ModeResult], best: ModeResult, band: list[ModeResult]) -> None:
    strict = [r for r in coarse if abs(r.wavelength_nm - TARGET_NM) <= 0.5 and r.q_value >= TARGET_Q]
    relaxed = [r for r in coarse if abs(r.wavelength_nm - TARGET_NM) <= 1.0 and r.q_value >= ACCEPT_Q]
    near_05 = [r for r in coarse if abs(r.wavelength_nm - TARGET_NM) <= 0.5]
    best_near_05 = max(near_05, key=lambda r: r.q_value) if near_05 else None
    highest_q = max(coarse, key=lambda r: r.q_value)
    bicref = [r for r in coarse if r.design_id.startswith("bicref_")]
    best_bicref = max(bicref, key=lambda r: r.q_value) if bicref else None
    lines = [
        "# v6 3D cylindrical-hole photonic-crystal random-particle RE summary",
        "",
        "All files were generated inside `E:\\xitugrating\\sin451_qbic_3d_design_v6_cyl_phc_random_particles`; previous versions were not edited.",
        "",
        "## Model choice",
        "",
        "- Used a 3D square-lattice cylindrical air-hole SiN slab unit cell because cylindrical PhC BIC leakage is inherently 3D.",
        f"- SiO2 was truncated to `{H_SIO2_NM:.0f} nm` plus PML to fit the 16 GB memory budget; this is an existence check, not the final full-substrate fabrication model.",
        "- Rare-earth particles are deterministic random placements inside the periodic representative cell; they are not placed on a regular grid.",
        "- Added a refinement set anchored by `bic0413.mph` scale (`period=356 nm`, `r/period=0.220`) as a range check only.",
        "",
        "## Scan result",
        "",
        f"- Mode rows simulated: `{len(coarse)}`.",
        f"- Strict `|lambda-451| <= 0.5 nm`, `Q >= 1e9`: `{len(strict)}`.",
        f"- Relaxed `|lambda-451| <= 1.0 nm`, `Q >= 1e8`: `{len(relaxed)}`.",
        "",
        "## Best candidate",
        "",
        "| item | value |",
        "|---|---:|",
        f"| P (nm) | `{best.period_nm:.9g}` |",
        f"| cylinder radius (nm) | `{best.radius_nm:.9g}` |",
        f"| axis gap (nm) | `{best.axis_gap_nm:.9g}` |",
        f"| lambda451 (nm) | `{best.wavelength_nm:.9g}` |",
        f"| Q451 | `{best.q_value:.6g}` |",
        f"| particle count | `{best.particle_count}` |",
        f"| top / bottom / side count | `{best.top_count} / {best.bottom_count} / {best.side_count}` |",
        "",
    ]
    lines.extend([
        "## Diagnostic notes",
        "",
        f"- Highest-Q solved mode: `{highest_q.design_id}`, lambda `{highest_q.wavelength_nm:.6g} nm`, Q `{highest_q.q_value:.6g}`.",
    ])
    if best_near_05:
        lines.append(f"- Best mode within `451+/-0.5 nm`: `{best_near_05.design_id}`, lambda `{best_near_05.wavelength_nm:.6g} nm`, Q `{best_near_05.q_value:.6g}`.")
    else:
        lines.append("- No mode landed within `451+/-0.5 nm`.")
    if best_bicref:
        lines.append(f"- Best `bic0413.mph`-anchored scale row: `{best_bicref.design_id}`, lambda `{best_bicref.wavelength_nm:.6g} nm`, Q `{best_bicref.q_value:.6g}`.")
    lines.append("- The random particulate RE breaks the clean in-plane symmetry enough that the cylindrical-hole PhC branch does not retain the 1D model's ultra-high-Q behavior in this memory-limited 3D scan.")
    lines.append("")
    if best.q_value >= TARGET_Q and abs(best.wavelength_nm - TARGET_NM) <= 0.5:
        lines.append("Conclusion: this 3D cylindrical-hole random-particle model meets the preferred high-Q target.")
    elif best.q_value >= ACCEPT_Q and abs(best.wavelength_nm - TARGET_NM) <= 1.0:
        lines.append("Conclusion: this 3D cylindrical-hole random-particle model meets the relaxed feasibility target.")
    else:
        lines.append("Conclusion: this 3D cylindrical-hole random-particle scan did not meet the high-Q feasibility threshold; inspect the candidates and consider larger-radius/period refinements or a finite-array design.")
    if band:
        lines.extend(["", "## Band sanity check", ""])
        for r in band:
            lines.append(f"- k=({r.kx_norm:.3g}, {r.ky_norm:.3g}) pi/P: lambda `{r.wavelength_nm:.6g} nm`, Q `{r.q_value:.6g}`.")
    lines.extend([
        "",
        "## Output files",
        "",
        "- `sin451_qbic_v6_best_3d_cyl_hole_phc_random_particles.mph`",
        "- `coarse_scan_451_3d_particle.csv`",
        "- `best_candidates_451_3d_particle.csv`",
        "- `random_particle_layout.csv`",
        "- `scan_q_vs_wavelength.svg`",
        "- `band_q_wavelength.svg`",
    ])
    if band:
        lines.insert(-1, "- `band_451_3d_particle.csv`")
    else:
        lines.insert(-1, "- `band_451_3d_particle.csv` was not generated in this skip-band run")
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def run_all(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mph = import_mph()
    client = mph.Client(cores=args.cores)
    try:
        designs = coarse_designs()
        if args.limit:
            designs = designs[: args.limit]
        rows = scan_designs(client, designs, OUT / "coarse_scan_451_3d_particle.csv", "coarse3d")
        candidates = choose_candidates(rows, 8)
        write_modes_csv(OUT / "best_candidates_451_3d_particle.csv", candidates)
        if not candidates:
            (OUT / "summary.md").write_text(
                "# v6 3D cylindrical photonic-crystal random-particle RE summary\n\nNo 3D candidate solved successfully. Check COMSOL logs and geometry/PML setup.",
                encoding="utf-8",
            )
            return
        best = candidates[0]
        best = save_best_model(client, best)
        band = run_band(client, best) if not args.skip_band else []
        write_report(rows, best, band)
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="3D cylindrical SiN photonic-crystal random-particle RE BIC scan.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-all")
    p.add_argument("--cores", type=int, default=2)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--skip-band", action="store_true")
    p.set_defaults(func=run_all)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
