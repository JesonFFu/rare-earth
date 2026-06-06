from __future__ import annotations

import argparse
import csv
import ctypes
from pathlib import Path

import v10_re50_cap_dual_bic as v


OUT = v.OUT / "family3d"


PREFERRED_IDS = {
    "single_air_hole_slab_reference": "mono_hole_P540_R74",
    "rectangular_air_hole_slab": "rect_lattice_hole_Px540_Py460_R74",
    "dual_hole_radius_supercell": "dual_radius_P540_r97_151",
    "dual_hole_shift_supercell": "dual_shift_P540_R130_dx35",
    "folded_2x2_air_hole_radius": "fold2x2_radius_P540_R95_135",
    "folded_2x2_air_hole_shift": "fold2x2_shift_P540_R118_d45",
    "capsule_elliptic_air_hole": "capsule_hole_Px700_Py520_L320_W140_x",
    "dual_capsule_elliptic_supercell": "dual_capsule_P600_L260_330",
    "main_hole_with_offset_satellite": "main_sat_hole_P700_R125_60",
    "wide_slot_cross_hole": "slot_cross_Px700_Py500",
    "dual_period_1d_grating": "grating_dual_P420_520_S160",
    "tetramer_hole_slab": "tetramer_hole_P1000_R135_90",
    "tetramer_pillar_array": "tetramer_pillar_P980_R145_95",
}


class MemoryStatusEx(ctypes.Structure):
    _fields_ = [
        ("dwLength", ctypes.c_ulong),
        ("dwMemoryLoad", ctypes.c_ulong),
        ("ullTotalPhys", ctypes.c_ulonglong),
        ("ullAvailPhys", ctypes.c_ulonglong),
        ("ullTotalPageFile", ctypes.c_ulonglong),
        ("ullAvailPageFile", ctypes.c_ulonglong),
        ("ullTotalVirtual", ctypes.c_ulonglong),
        ("ullAvailVirtual", ctypes.c_ulonglong),
        ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
    ]


def memory_status() -> tuple[float, float, float]:
    stat = MemoryStatusEx()
    stat.dwLength = ctypes.sizeof(MemoryStatusEx)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
    total_gb = stat.ullTotalPhys / (1024**3)
    avail_gb = stat.ullAvailPhys / (1024**3)
    used_percent = 100.0 * (1.0 - avail_gb / total_gb)
    return total_gb, avail_gb, used_percent


def require_memory(min_free_gb: float, label: str) -> None:
    total_gb, avail_gb, used_percent = memory_status()
    print(f"[memory] {label}: total={total_gb:.2f} GB free={avail_gb:.2f} GB used={used_percent:.1f}%", flush=True)
    if avail_gb < min_free_gb:
        raise RuntimeError(
            f"Available memory {avail_gb:.2f} GB is below guard threshold {min_free_gb:.2f} GB. "
            "Stop for review instead of launching another COMSOL model."
        )


def chosen_designs() -> list[v.Design]:
    designs = v.baseline_and_literature_designs()
    by_id = {d.design_id: d for d in designs}
    by_family: dict[str, list[v.Design]] = {}
    for d in designs:
        by_family.setdefault(d.family, []).append(d)
    out: list[v.Design] = []
    for family in sorted(by_family):
        preferred = PREFERRED_IDS.get(family)
        if preferred and preferred in by_id:
            out.append(by_id[preferred])
        else:
            out.append(by_family[family][0])
    return out


def result_fields() -> list[str]:
    return [
        "design_id",
        "family",
        "notes",
        "px_nm",
        "py_nm",
        "feature_min_nm",
        "lambda980_nm",
        "q980",
        "mode980_index",
        "lambda451_nm",
        "q451",
        "mode451_index",
        "dual_score",
        "mph_980",
        "mph_451",
    ]


def read_done(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open("r", newline="", encoding="utf-8") as handle:
        return {row["design_id"] for row in csv.DictReader(handle)}


def append_result(path: Path, design: v.Design, m980: v.ModeResult, m451: v.ModeResult, mph980: Path, mph451: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    score = v.dual_score(design, m980, m451)
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=result_fields())
        if not exists:
            writer.writeheader()
        writer.writerow(
            {
                "design_id": design.design_id,
                "family": design.family,
                "notes": design.notes,
                "px_nm": design.px_nm,
                "py_nm": design.py_nm,
                "feature_min_nm": design.feature_min_nm,
                "lambda980_nm": m980.wavelength_nm,
                "q980": m980.q_value,
                "mode980_index": m980.mode_index,
                "lambda451_nm": m451.wavelength_nm,
                "q451": m451.q_value,
                "mode451_index": m451.mode_index,
                "dual_score": score,
                "mph_980": str(mph980),
                "mph_451": str(mph451),
            }
        )


def run(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "family3d_pilot_results.csv"
    done = read_done(csv_path)
    designs = chosen_designs()
    if args.family:
        designs = [d for d in designs if d.family == args.family]
    if args.limit:
        designs = designs[: args.limit]
    if args.profile == "full":
        v.H_AIR_NM = 900.0
        v.H_PML_NM = 500.0
        v.H_SI_BUF_NM = 250.0
    elif args.profile == "balanced":
        v.H_AIR_NM = 650.0
        v.H_PML_NM = 350.0
        v.H_SI_BUF_NM = 150.0
    elif args.profile == "lite":
        v.H_AIR_NM = 450.0
        v.H_PML_NM = 250.0
        v.H_SI_BUF_NM = 100.0
    else:
        raise ValueError(f"Unknown profile {args.profile}")
    v.MESH_AUTO_SIZE = args.mesh_size
    v.NEIGS = args.neigs
    print(
        f"[profile] {args.profile}: h_air={v.H_AIR_NM} nm h_pml={v.H_PML_NM} nm "
        f"h_si_buf={v.H_SI_BUF_NM} nm h_sio2={v.H_SIO2_NM} nm mesh={v.MESH_AUTO_SIZE} neigs={v.NEIGS}",
        flush=True,
    )
    require_memory(args.min_free_gb, "before COMSOL client")
    if args.dry_run:
        print("[dry-run] COMSOL client will not be started.", flush=True)
        for d in designs:
            state = "done" if d.design_id in done else "pending"
            print(f"[dry-run] {state} {d.family} {d.design_id} min={d.feature_min_nm:.1f}", flush=True)
        return
    mph = v.import_mph()
    client = mph.Client(cores=args.cores)
    try:
        for index, design in enumerate(designs, 1):
            if design.design_id in done and not args.force:
                print(f"[skip] {design.family} {design.design_id}", flush=True)
                continue
            require_memory(args.min_free_gb, f"before {design.design_id}")
            print(
                f"[pilot] {index}/{len(designs)} {design.family} {design.design_id} "
                f"Px={design.px_nm:.1f} Py={design.py_nm:.1f} min={design.feature_min_nm:.1f}",
                flush=True,
            )
            safe = design.design_id.replace("/", "_")
            mph980 = OUT / f"{safe}_980.mph"
            mph451 = OUT / f"{safe}_451.mph"
            m980 = v.save_solved_model(client, design, v.TARGET_980_NM, mph980, f"pilot_{safe}_980")
            print(f"  980 lambda={m980.wavelength_nm:.6f} Q={m980.q_value:.4g}", flush=True)
            require_memory(args.min_free_gb, f"after {design.design_id} 980")
            m451 = v.save_solved_model(client, design, v.TARGET_451_NM, mph451, f"pilot_{safe}_451")
            print(f"  451 lambda={m451.wavelength_nm:.6f} Q={m451.q_value:.4g}", flush=True)
            require_memory(args.min_free_gb, f"after {design.design_id} 451")
            append_result(csv_path, design, m980, m451, mph980, mph451)
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="Save one 3D COMSOL pilot model per v10 structure family.")
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--family", default="")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--mesh-size", type=int, default=7, help="COMSOL autoMeshSize, larger is coarser/lower memory.")
    parser.add_argument("--neigs", type=int, default=6, help="Number of eigenmodes to request per target wavelength.")
    parser.add_argument("--min-free-gb", type=float, default=7.5, help="Stop before a model if free physical memory is below this.")
    parser.add_argument("--dry-run", action="store_true", help="Check memory and list selected designs without starting COMSOL.")
    parser.add_argument(
        "--profile",
        choices=["full", "balanced", "lite"],
        default="full",
        help="Layer simplification profile. SiO2 remains 1500 nm in all profiles.",
    )
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
