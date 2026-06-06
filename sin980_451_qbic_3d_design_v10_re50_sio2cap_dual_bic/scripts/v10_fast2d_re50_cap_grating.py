from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
OUT = ROOT / "outputs" / "fast2d"
INPUTS = ROOT / "inputs"

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

H_PML_NM = 600.0
H_SI_BUF_NM = 250.0
H_SIO2_NM = 1500.0
H_RE_NM = 50.0
H_CAP_NM = 12.5
H_SIN_NM = 300.0
H_AIR_NM = 1000.0

MIN_FEATURE_NM = 100.0
PREFERRED_GAP_NM = 120.0


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
    notes: str = ""
    groove_widths_nm: tuple[float, ...] = ()
    groove_depths_nm: tuple[float, ...] = ()

    @property
    def length_nm(self) -> float:
        return sum(self.pitches_nm)

    @property
    def n_periods(self) -> int:
        return len(self.pitches_nm)

    @property
    def ridges(self) -> list[Ridge]:
        x = -self.length_nm / 2.0
        out: list[Ridge] = []
        for pitch, width in zip(self.pitches_nm, self.widths_nm):
            center = x + pitch / 2.0
            out.append(Ridge(center - width / 2.0, center + width / 2.0))
            x += pitch
        return out

    @property
    def gaps_nm(self) -> list[float]:
        ridges = self.ridges
        out: list[float] = []
        for a, b in zip(ridges, ridges[1:]):
            out.append(b.left_nm - a.right_nm)
        wrap_gap = (self.length_nm / 2.0 - ridges[-1].right_nm) + (ridges[0].left_nm + self.length_nm / 2.0)
        out.append(wrap_gap)
        return out

    @property
    def min_feature_nm(self) -> float:
        values = [min(self.widths_nm), min(self.gaps_nm)]
        for i, width in enumerate(self.widths_nm):
            gw = self.groove_widths_nm[i] if i < len(self.groove_widths_nm) else 0.0
            gd = self.groove_depths_nm[i] if i < len(self.groove_depths_nm) else 0.0
            if gw > 0.0 and gd > 0.0:
                values.extend([gw, gd, H_SIN_NM - gd, (width - gw) / 2.0])
        return min(values)

    @property
    def fabricable(self) -> bool:
        return self.min_feature_nm >= MIN_FEATURE_NM

    @property
    def preferred_gap_ok(self) -> bool:
        return min(self.gaps_nm) >= PREFERRED_GAP_NM

    def groove(self, index: int) -> tuple[float, float]:
        gw = self.groove_widths_nm[index] if index < len(self.groove_widths_nm) else 0.0
        gd = self.groove_depths_nm[index] if index < len(self.groove_depths_nm) else 0.0
        return gw, gd


@dataclass
class ModeResult:
    design_id: str
    family: str
    target_nm: float
    length_nm: float
    n_periods: int
    min_feature_nm: float
    min_gap_nm: float
    kx_norm: float
    mode_index: int
    freq_thz_real: float
    freq_thz_imag: float
    wavelength_nm: float
    q_value: float
    score: float


@dataclass
class DualResult:
    design: Design
    mode980: ModeResult
    mode451: ModeResult
    score: float


def configure_process_environment() -> None:
    os.environ["COMSOLROOT"] = COMSOLROOT
    os.environ["JAVA_HOME"] = JAVA_HOME
    os.environ["PATH"] = rf"{COMSOLROOT}\bin\win64;{JAVA_HOME}\bin;" + os.environ.get("PATH", "")
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")


def import_mph():
    configure_process_environment()
    import mph

    return mph


def fmt(value: float | str) -> str:
    if isinstance(value, str):
        return value
    if abs(value) < 1e-12:
        value = 0.0
    return f"{value:.12g}"


def set_param(model, name: str, value: str, description: str = "") -> None:
    if description:
        model.java.param().set(name, value, description)
    else:
        model.java.param().set(name, value)


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


def target_indices(target_nm: float) -> tuple[float, float, float, float]:
    if abs(target_nm - TARGET_451_NM) < 10:
        return N_SIN_451, N_SIO2_451, N_RE_451, N_SI_451
    return N_SIN_980, N_SIO2_980, N_RE_980, N_SI_980


def unique_complex(values: Iterable[complex | float], tol: float = 1e-7) -> list[complex]:
    out: list[complex] = []
    for raw in values:
        val = complex(raw)
        if not any(abs(val - old) < tol for old in out):
            out.append(val)
    return out


def get_eigenfrequencies_thz(model, target_nm: float) -> list[complex]:
    values = model.evaluate("freq", unit="THz")
    flat = values.ravel() if hasattr(values, "ravel") else values
    modes = [m for m in unique_complex(flat) if math.isfinite(m.real) and abs(m.real) > 1e-9]
    modes.sort(key=lambda z: abs(C0 / (z.real * 1e12) * 1e9 - target_nm) if z.real else float("inf"))
    return modes


def mode_result(design: Design, target_nm: float, kx_norm: float, idx: int, freq_thz: complex) -> ModeResult:
    f_real = float(freq_thz.real)
    f_imag = float(freq_thz.imag)
    wavelength_nm = C0 / (f_real * 1e12) * 1e9 if f_real > 0 else math.inf
    q_value = abs(f_real) / (2.0 * abs(f_imag)) if abs(f_imag) > 0 else 1e18
    target_q = Q_TARGET_451 if target_nm < 700 else Q_TARGET_980
    detune = abs(wavelength_nm - target_nm)
    q_penalty = max(0.0, math.log10(target_q) - math.log10(max(q_value, 1.0)))
    score = detune + 4.0 * q_penalty
    return ModeResult(
        design.design_id,
        design.family,
        target_nm,
        design.length_nm,
        design.n_periods,
        design.min_feature_nm,
        min(design.gaps_nm),
        kx_norm,
        idx,
        f_real,
        f_imag,
        wavelength_nm,
        q_value,
        score,
    )


def best_mode(rows: list[ModeResult], target_nm: float) -> ModeResult:
    near = [r for r in rows if abs(r.wavelength_nm - target_nm) <= 8.0]
    return sorted(near or rows, key=lambda r: (r.score, abs(r.wavelength_nm - target_nm), -r.q_value))[0]


def dual_score(design: Design, m980: ModeResult, m451: ModeResult) -> float:
    d980 = abs(m980.wavelength_nm - TARGET_980_NM)
    d451 = abs(m451.wavelength_nm - TARGET_451_NM)
    q980_pen = max(0.0, math.log10(Q_TARGET_980) - math.log10(max(m980.q_value, 1.0)))
    q451_pen = max(0.0, math.log10(Q_TARGET_451) - math.log10(max(m451.q_value, 1.0)))
    gap_pen = 0.0 if design.preferred_gap_ok else 1.5
    return d980 * 2.0 + d451 + 7.0 * q980_pen + 4.0 * q451_pen + gap_pen


def build_model(client, design: Design, target_nm: float, name: str):
    n_sin, n_sio2, n_re, n_si = target_indices(target_nm)
    model = client.create(name)
    jmodel = model.java

    total_h = 2 * H_PML_NM + H_SI_BUF_NM + H_SIO2_NM + H_RE_NM + H_CAP_NM + H_SIN_NM + H_AIR_NM
    z_si = H_PML_NM
    z_sio2 = H_PML_NM + H_SI_BUF_NM
    z_re = z_sio2 + H_SIO2_NM
    z_cap = z_re + H_RE_NM
    z_sin = z_cap + H_CAP_NM
    z_top = z_sin + H_SIN_NM
    z_top_pml = z_top + H_AIR_NM
    x_min = -design.length_nm / 2.0
    x_max = design.length_nm / 2.0

    set_param(model, "lambda_target", f"{target_nm:.9g}[nm]", "target wavelength")
    set_param(model, "L", f"{design.length_nm:.9g}[nm]", "supercell length")
    set_param(model, "kx", "0[1/m]", "Floquet wave vector")
    set_param(model, "h_sio2_eff", f"{H_SIO2_NM:.9g}[nm]", "effective SiO2 thickness")
    set_param(model, "h_re", f"{H_RE_NM:.9g}[nm]", "continuous rare-earth film")
    set_param(model, "h_cap_sio2", f"{H_CAP_NM:.9g}[nm]", "continuous SiO2 cap")
    set_param(model, "h_sin", f"{H_SIN_NM:.9g}[nm]", "patterned SiN")
    set_param(model, "z_re_mid", f"{z_re + H_RE_NM / 2:.9g}[nm]", "rare-earth film mid-plane")
    set_param(model, "z_slab_mid", f"{z_sin + H_SIN_NM / 2:.9g}[nm]", "SiN mid-plane")
    set_param(model, "n_sin", f"{n_sin:.9f}", "target-dependent SiN index")
    set_param(model, "n_sio2", f"{n_sio2:.9f}", "target-dependent SiO2 index")
    set_param(model, "n_re", f"{n_re:.9f}", "target-dependent rare-earth-film index")
    set_param(model, "n_si", f"{n_si:.9f}", "target-dependent Si index")

    jmodel.component().create("comp1", True)
    comp = jmodel.component("comp1")
    comp.label("v10 fast 2D x-z supercell")
    comp.geom().create("geom1", 2)
    geom = comp.geom("geom1")
    geom.lengthUnit("nm")

    create_rectangle(geom, "si_pml", (x_min, 0.0), (design.length_nm, H_PML_NM), "Si bottom PML")
    create_rectangle(geom, "si_buf", (x_min, z_si), (design.length_nm, H_SI_BUF_NM), "Si buffer")
    create_rectangle(geom, "sio2", (x_min, z_sio2), (design.length_nm, H_SIO2_NM), "SiO2")
    create_rectangle(geom, "re_film", (x_min, z_re), (design.length_nm, H_RE_NM), "continuous rare-earth film")
    create_rectangle(geom, "cap_sio2", (x_min, z_cap), (design.length_nm, H_CAP_NM), "continuous SiO2 cap")

    sin_tags: list[str] = []
    air_tags: list[str] = []
    ridges = design.ridges
    for i, ridge in enumerate(ridges):
        gw, gd = design.groove(i)
        if gw > 0.0 and gd > 0.0:
            groove_left = (ridge.left_nm + ridge.right_nm - gw) / 2.0
            groove_right = groove_left + gw
            bottom_h = H_SIN_NM - gd
            create_rectangle(geom, f"sin{i}_bottom", (ridge.left_nm, z_sin), (ridge.width_nm, bottom_h), "SiN below wide shallow groove")
            create_rectangle(geom, f"sin{i}_rib_l", (ridge.left_nm, z_sin + bottom_h), (groove_left - ridge.left_nm, gd), "left SiN rib beside groove")
            create_rectangle(geom, f"sin{i}_rib_r", (groove_right, z_sin + bottom_h), (ridge.right_nm - groove_right, gd), "right SiN rib beside groove")
            create_rectangle(geom, f"air_groove{i}", (groove_left, z_sin + bottom_h), (gw, gd), "wide shallow air groove in SiN only")
            sin_tags.extend([f"geom1_sin{i}_bottom_dom", f"geom1_sin{i}_rib_l_dom", f"geom1_sin{i}_rib_r_dom"])
            air_tags.append(f"geom1_air_groove{i}_dom")
        else:
            create_rectangle(geom, f"sin{i}", (ridge.left_nm, z_sin), (ridge.width_nm, H_SIN_NM), "patterned SiN ridge")
            sin_tags.append(f"geom1_sin{i}_dom")

    def add_air_gap(tag: str, left: float, right: float) -> None:
        if right - left <= 1e-9:
            return
        create_rectangle(geom, f"air_gap_{tag}", (left, z_sin), (right - left, H_SIN_NM), "air in through-etched SiN gap")
        air_tags.append(f"geom1_air_gap_{tag}_dom")

    for i, (a, b) in enumerate(zip(ridges, ridges[1:])):
        add_air_gap(str(i), a.right_nm, b.left_nm)
    add_air_gap("left", x_min, ridges[0].left_nm)
    add_air_gap("right", ridges[-1].right_nm, x_max)
    create_rectangle(geom, "air", (x_min, z_top), (design.length_nm, H_AIR_NM), "air above SiN")
    create_rectangle(geom, "air_pml", (x_min, z_top_pml), (design.length_nm, H_PML_NM), "air top PML")
    air_tags.extend(["geom1_air_dom", "geom1_air_pml_dom"])
    geom.run()

    create_union_selection(comp, "sin_dom", "2", sin_tags)
    create_union_selection(comp, "air_dom", "2", air_tags)
    create_union_selection(comp, "si_dom", "2", ["geom1_si_pml_dom", "geom1_si_buf_dom"])
    create_box_selection(comp, "left_bnd", "1", f"{x_min-1}", f"{x_min+1}", "-1", f"{total_h+1}")
    create_box_selection(comp, "right_bnd", "1", f"{x_max-1}", f"{x_max+1}", "-1", f"{total_h+1}")
    create_union_selection(comp, "periodic_bnd", "1", ["left_bnd", "right_bnd"])

    pml_bot = comp.coordSystem().create("pml_bot", "PML")
    pml_bot.label("Bottom Si PML")
    pml_bot.selection().named("geom1_si_pml_dom")
    pml_top = comp.coordSystem().create("pml_top", "PML")
    pml_top.label("Top air PML")
    pml_top.selection().named("geom1_air_pml_dom")

    create_refractive_material(comp, "mat_air", "Air", "air_dom", "1", "0")
    create_refractive_material(comp, "mat_si", "Effective Si", "si_dom", "n_si", "0")
    create_refractive_material(comp, "mat_sio2", "SiO2 substrate", "geom1_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_re", "50 nm rare-earth film", "geom1_re_film_dom", "n_re", "0")
    create_refractive_material(comp, "mat_cap", "12.5 nm SiO2 cap", "geom1_cap_sio2_dom", "n_sio2", "0")
    create_refractive_material(comp, "mat_sin", "patterned SiN", "sin_dom", "n_sin", "0")

    ewfd = comp.physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
    ewfd.label("Electromagnetic Waves, Frequency Domain")
    pc = ewfd.create("pc1", "PeriodicCondition", 1)
    pc.label("Floquet periodic x")
    pc.selection().named("periodic_bnd")
    pc.set("PeriodicType", "Floquet")
    pc.set("kFloquet", ["kx", "0", "0"])

    mesh = comp.mesh().create("mesh1")
    mesh.label("fast 2D physics mesh")
    mesh.autoMeshSize(4)

    study = jmodel.study().create("std1")
    study.label(f"Eigenfrequency near {target_nm:.0f} nm")
    eig = study.create("eig", "Eigenfrequency")
    eig.set("shift", "c_const/lambda_target")
    eig.set("neigs", "12")
    eig.set("neigsmanual", "12")
    eig.set("neigsactive", "on")
    eig.set("eigunit", "THz")
    return model


def solve_model(model, design: Design, target_nm: float, kx_norm: float = 0.0) -> list[ModeResult]:
    set_param(model, "kx", f"{kx_norm:.12g}*pi/L")
    model.java.component("comp1").geom("geom1").run()
    model.java.component("comp1").mesh("mesh1").run()
    model.java.study("std1").run()
    modes = get_eigenfrequencies_thz(model, target_nm)
    return [mode_result(design, target_nm, kx_norm, i + 1, f) for i, f in enumerate(modes)]


def solve_design(client, design: Design, target_nm: float) -> list[ModeResult]:
    model = build_model(client, design, target_nm, f"v10_fast2d_{design.design_id}_{int(target_nm)}")
    try:
        return solve_model(model, design, target_nm, 0.0)
    finally:
        client.remove(model.name())


def add_field_plot_groups(model, label: str) -> None:
    res = model.java.result()
    for tag, expr in [(f"field_{label}_Ez2_xz", "ewfd.Ez^2"), (f"field_{label}_E2_xz", "ewfd.normE^2")]:
        pg = res.create(tag, "PlotGroup2D")
        pg.label(f"{label}: {expr} xz map")
        pg.set("data", "dset1")
        surf = pg.feature().create("surf", "Surface")
        surf.set("expr", expr)
        surf.set("resolution", "normal")


def save_solved_model(client, design: Design, target_nm: float, path: Path, name: str) -> ModeResult:
    model = build_model(client, design, target_nm, name)
    try:
        rows = solve_model(model, design, target_nm, 0.0)
        selected = best_mode(rows, target_nm)
        add_field_plot_groups(model, str(int(target_nm)))
        model.save(str(path))
        write_modes_csv(path.with_suffix(".modes.csv"), rows)
        write_modes_csv(path.with_suffix(".selected.csv"), [selected])
        return selected
    finally:
        client.remove(model.name())


def candidates() -> list[Design]:
    out: dict[str, Design] = {}

    def add(d: Design) -> None:
        if d.fabricable:
            out[d.design_id] = d

    for p in [294.0, 320.0, 420.0, 500.0, 540.0, 560.0, 600.0, 700.0, 840.0, 980.0]:
        for fill in [0.42, 0.55, 0.68, 0.75]:
            add(Design(f"primitive_P{p:.0f}_f{fill:.2f}".replace(".", "p"), "primitive_1d_grating", (p,), (p * fill,), "single-period 1D reference"))

    # Combine known 451-scale and 980-scale branches in a single supercell.
    for p1, w1, p2, w2 in [
        (294.0, 159.0, 554.0, 414.0),
        (300.0, 165.0, 560.0, 400.0),
        (320.0, 180.0, 600.0, 420.0),
        (360.0, 200.0, 650.0, 455.0),
        (420.0, 240.0, 760.0, 520.0),
    ]:
        add(Design(f"hybrid_451_980_P{p1:.0f}_{p2:.0f}", "hybrid_451_980_supercell", (p1, p2), (w1, w2), "one visible-scale ridge plus one 980-scale ridge"))

    # Brillouin-zone folded two-ridge cells around literature dimerization logic.
    for p in [420.0, 500.0, 540.0, 600.0, 700.0]:
        for dw in [40.0, 70.0, 100.0]:
            w0 = 0.62 * p
            add(Design(f"dimer_P{p:.0f}_dW{dw:.0f}", "dimerized_same_pitch_grating", (p, p), (w0 - dw / 2, w0 + dw / 2), "alternating SiN ridge width"))
        for dp in [40.0, 80.0, 120.0]:
            add(Design(f"dimer_P{p:.0f}_dP{dp:.0f}", "dimerized_pitch_grating", (p - dp / 2, p + dp / 2), (0.62 * (p - dp / 2), 0.62 * (p + dp / 2)), "alternating pitch"))

    # Three-ridge supercells with one visible-tuned segment and two pump-tuned segments.
    for scale in [1.0, 1.08, 0.92]:
        add(
            Design(
                f"trimer_visible_pump_s{scale:.2f}".replace(".", "p"),
                "trimer_visible_pump_supercell",
                (294.0 * scale, 540.0 * scale, 540.0 * scale),
                (159.0 * scale, 385.0 * scale, 410.0 * scale),
                "three-ridge supercell for folded dual branches",
            )
        )
    return list(out.values())


def write_modes_csv(path: Path, rows: list[ModeResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "design_id",
        "family",
        "target_nm",
        "length_nm",
        "n_periods",
        "min_feature_nm",
        "min_gap_nm",
        "kx_norm",
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


def write_dual_csv(path: Path, rows: list[DualResult]) -> None:
    fields = [
        "design_id",
        "family",
        "notes",
        "length_nm",
        "n_periods",
        "pitches_nm",
        "widths_nm",
        "groove_widths_nm",
        "groove_depths_nm",
        "min_feature_nm",
        "min_gap_nm",
        "lambda980_nm",
        "q980",
        "mode980_index",
        "lambda451_nm",
        "q451",
        "mode451_index",
        "dual_score",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for r in sorted(rows, key=lambda item: item.score):
            writer.writerow(
                {
                    "design_id": r.design.design_id,
                    "family": r.design.family,
                    "notes": r.design.notes,
                    "length_nm": r.design.length_nm,
                    "n_periods": r.design.n_periods,
                    "pitches_nm": ";".join(f"{v:.6g}" for v in r.design.pitches_nm),
                    "widths_nm": ";".join(f"{v:.6g}" for v in r.design.widths_nm),
                    "groove_widths_nm": ";".join(f"{v:.6g}" for v in r.design.groove_widths_nm),
                    "groove_depths_nm": ";".join(f"{v:.6g}" for v in r.design.groove_depths_nm),
                    "min_feature_nm": r.design.min_feature_nm,
                    "min_gap_nm": min(r.design.gaps_nm),
                    "lambda980_nm": r.mode980.wavelength_nm,
                    "q980": r.mode980.q_value,
                    "mode980_index": r.mode980.mode_index,
                    "lambda451_nm": r.mode451.wavelength_nm,
                    "q451": r.mode451.q_value,
                    "mode451_index": r.mode451.mode_index,
                    "dual_score": r.score,
                }
            )


def scan(client, designs: list[Design], path: Path, limit: int = 0) -> list[DualResult]:
    if limit:
        designs = designs[:limit]
    rows: list[DualResult] = []
    for i, d in enumerate(designs, 1):
        print(f"[2D] {i}/{len(designs)} {d.family} {d.design_id} L={d.length_nm:.1f} min={d.min_feature_nm:.1f}", flush=True)
        try:
            m980 = best_mode(solve_design(client, d, TARGET_980_NM), TARGET_980_NM)
            print(f"  980 lambda={m980.wavelength_nm:.6f} Q={m980.q_value:.4g}", flush=True)
            m451 = best_mode(solve_design(client, d, TARGET_451_NM), TARGET_451_NM)
            print(f"  451 lambda={m451.wavelength_nm:.6f} Q={m451.q_value:.4g}", flush=True)
            rows.append(DualResult(d, m980, m451, dual_score(d, m980, m451)))
            write_dual_csv(path, rows)
        except Exception as exc:
            print(f"  FAILED {d.design_id}: {exc}", flush=True)
    return rows


def choose_best(rows: list[DualResult]) -> DualResult:
    strict = [
        r
        for r in rows
        if abs(r.mode980.wavelength_nm - TARGET_980_NM) <= 1.0
        and r.mode980.q_value >= Q_TARGET_980
        and abs(r.mode451.wavelength_nm - TARGET_451_NM) <= 2.0
        and r.mode451.q_value >= Q_TARGET_451
    ]
    return sorted(strict or rows, key=lambda r: r.score)[0]


def family_bests(rows: list[DualResult]) -> list[DualResult]:
    out: list[DualResult] = []
    for fam in sorted({r.design.family for r in rows}):
        out.append(sorted([r for r in rows if r.design.family == fam], key=lambda r: r.score)[0])
    return sorted(out, key=lambda r: r.score)


def refine_designs(seed: DualResult) -> list[Design]:
    d = seed.design
    out: list[Design] = []
    for sp in [0.94, 0.98, 1.0, 1.02, 1.06]:
        for sw in [0.92, 0.98, 1.0, 1.02, 1.08]:
            nd = Design(
                f"refine_{d.design_id}_sp{sp:.2f}_sw{sw:.2f}".replace(".", "p"),
                d.family,
                tuple(p * sp for p in d.pitches_nm),
                tuple(w * sp * sw for w in d.widths_nm),
                "local refinement from best 2D candidate",
            )
            if nd.fabricable:
                out.append(nd)
    return out


def run_band(client, design: Design, target_nm: float, csv_path: Path, q_path: Path) -> list[ModeResult]:
    model = build_model(client, design, target_nm, f"v10_fast2d_band_{int(target_nm)}")
    kpath = [0.0, 0.002, 0.006, 0.02, 0.06, 0.15, 0.30, 0.50]
    rows: list[ModeResult] = []
    try:
        for k in kpath:
            current = solve_model(model, design, target_nm, k)
            rows.append(best_mode(current, target_nm))
        write_modes_csv(csv_path, rows)
        write_modes_csv(q_path, rows[:4])
        return rows
    finally:
        client.remove(model.name())


def svg_band(rows: list[ModeResult], path: Path, title: str) -> None:
    if not rows:
        return
    width, height = 820, 500
    left, right, top, bottom = 75, 80, 40, 70
    pw, ph = width - left - right, height - top - bottom
    wmin, wmax = min(r.wavelength_nm for r in rows) - 2, max(r.wavelength_nm for r in rows) + 2
    qmin, qmax = max(1.0, min(r.q_value for r in rows) / 2), max(r.q_value for r in rows) * 2

    def sx(i: int) -> float:
        return left + i * pw / max(1, len(rows) - 1)

    def sy_w(v: float) -> float:
        return top + (wmax - v) / (wmax - wmin) * ph

    def sy_q(v: float) -> float:
        return top + (math.log10(qmax) - math.log10(max(v, qmin))) / (math.log10(qmax) - math.log10(qmin)) * ph

    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    lines.append(f'<text x="{width/2}" y="25" text-anchor="middle" font-family="Arial" font-size="17">{title}</text>')
    lines.append(f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#111"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#111"/>')
    lines.append(f'<line x1="{left+pw}" y1="{top}" x2="{left+pw}" y2="{top+ph}" stroke="#111"/>')
    wpts = " ".join(f"{sx(i):.2f},{sy_w(r.wavelength_nm):.2f}" for i, r in enumerate(rows))
    qpts = " ".join(f"{sx(i):.2f},{sy_q(r.q_value):.2f}" for i, r in enumerate(rows))
    lines.append(f'<polyline fill="none" stroke="#1f77b4" stroke-width="2" points="{wpts}"/>')
    lines.append(f'<polyline fill="none" stroke="#d62728" stroke-width="2" points="{qpts}"/>')
    for i, r in enumerate(rows):
        lines.append(f'<text x="{sx(i):.2f}" y="{top+ph+20}" text-anchor="middle" font-family="Arial" font-size="10">{r.kx_norm:g}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_summary(best: DualResult, rows: list[DualResult], band980: list[ModeResult], band451: list[ModeResult]) -> None:
    passes = best.mode980.q_value >= Q_TARGET_980 and best.mode451.q_value >= Q_TARGET_451
    lines = [
        "# v10 Fast 2D Dual-BIC Screening",
        "",
        f"Directory: `{ROOT}`. This 2D screening keeps the v10 stack and only patterns the 300 nm SiN layer.",
        "",
        "## Best 2D Candidate",
        "",
        f"- Design: `{best.design.design_id}`",
        f"- Family: `{best.design.family}`",
        f"- Pitches: `{'; '.join(f'{v:.3f}' for v in best.design.pitches_nm)} nm`",
        f"- Widths: `{'; '.join(f'{v:.3f}' for v in best.design.widths_nm)} nm`",
        f"- Minimum feature/gap: `{best.design.min_feature_nm:.3f} nm`",
        f"- 980 nm branch: lambda `{best.mode980.wavelength_nm:.6f} nm`, Q `{best.mode980.q_value:.6g}`",
        f"- 451 nm branch: lambda `{best.mode451.wavelength_nm:.6f} nm`, Q `{best.mode451.q_value:.6g}`",
        f"- Meets requested Q thresholds: `{passes}`",
        "",
        "## Family Bests",
        "",
        "| family | design | lambda980 | Q980 | lambda451 | Q451 | min feature |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for r in family_bests(rows):
        lines.append(
            f"| {r.design.family} | {r.design.design_id} | {r.mode980.wavelength_nm:.4f} | {r.mode980.q_value:.4g} | "
            f"{r.mode451.wavelength_nm:.4f} | {r.mode451.q_value:.4g} | {r.design.min_feature_nm:.1f} |"
        )
    if band980:
        lines.extend(["", "## 980 nm Q(k)", ""])
        for r in band980[:4]:
            lines.append(f"- k={r.kx_norm:g} pi/L: lambda `{r.wavelength_nm:.6f} nm`, Q `{r.q_value:.6g}`")
    if band451:
        lines.extend(["", "## 451 nm Q(k)", ""])
        for r in band451[:4]:
            lines.append(f"- k={r.kx_norm:g} pi/L: lambda `{r.wavelength_nm:.6f} nm`, Q `{r.q_value:.6g}`")
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `fast2d_best_980.mph` and `fast2d_best_451.mph`: same geometry with target-dependent material constants.",
            "- `coarse_scan_2d.csv`, `fine_scan_2d.csv`, `dual_best_candidates_2d.csv`.",
            "- `band_980_2d.csv`, `band_451_2d.csv`, `q_near_gamma_980_2d.csv`, `q_near_gamma_451_2d.csv`.",
            "",
            "## Interpretation",
            "",
            "- Treat this as a fast existence screen. The earlier 3D coarse scan is retained separately in `outputs/coarse_dual_scan.csv`.",
        ]
    )
    (OUT / "summary_2d.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_all(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    INPUTS.mkdir(parents=True, exist_ok=True)
    mph = import_mph()
    client = mph.Client(cores=args.cores)
    try:
        coarse = scan(client, candidates(), OUT / "coarse_scan_2d.csv", args.limit)
        if not coarse:
            return
        seed = choose_best(coarse)
        fine = scan(client, refine_designs(seed), OUT / "fine_scan_2d.csv", args.fine_limit) if not args.skip_fine else []
        rows = coarse + fine
        best = choose_best(rows)
        write_dual_csv(OUT / "dual_best_candidates_2d.csv", sorted(rows, key=lambda r: r.score)[:30])
        save_solved_model(client, best.design, TARGET_980_NM, OUT / "fast2d_best_980.mph", "v10_fast2d_best_980")
        save_solved_model(client, best.design, TARGET_451_NM, OUT / "fast2d_best_451.mph", "v10_fast2d_best_451")
        # Provide the plan-level filename as a 2D fallback model while 3D refinement remains expensive.
        if not (ROOT / "outputs" / "overall_best_dual_bic.mph").exists():
            shutil.copy2(OUT / "fast2d_best_980.mph", ROOT / "outputs" / "overall_best_dual_bic.mph")
        band980: list[ModeResult] = []
        band451: list[ModeResult] = []
        if not args.skip_band:
            band980 = run_band(client, best.design, TARGET_980_NM, OUT / "band_980_2d.csv", OUT / "q_near_gamma_980_2d.csv")
            band451 = run_band(client, best.design, TARGET_451_NM, OUT / "band_451_2d.csv", OUT / "q_near_gamma_451_2d.csv")
            svg_band(band980, OUT / "band_980_2d.svg", "v10 fast 2D 980 nm band and Q(k)")
            svg_band(band451, OUT / "band_451_2d.svg", "v10 fast 2D 451 nm band and Q(k)")
        write_summary(best, rows, band980, band451)
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="v10 fast 2D screening for dual 980/451 BIC branches.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-all")
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--fine-limit", type=int, default=0)
    p.add_argument("--skip-fine", action="store_true")
    p.add_argument("--skip-band", action="store_true")
    p.set_defaults(func=run_all)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
