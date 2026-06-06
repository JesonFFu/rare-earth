from __future__ import annotations

import argparse
import csv
import ctypes
from pathlib import Path

import v10_re50_cap_dual_bic as v
from v10_family3d_pilot import PREFERRED_IDS, chosen_designs


OUT = v.OUT / "family3d_geometry_only"


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
        raise RuntimeError(f"Free memory {avail_gb:.2f} GB is below guard threshold {min_free_gb:.2f} GB.")


def apply_profile(profile: str) -> None:
    if profile == "full":
        v.H_AIR_NM = 900.0
        v.H_PML_NM = 500.0
        v.H_SI_BUF_NM = 250.0
    elif profile == "balanced":
        v.H_AIR_NM = 650.0
        v.H_PML_NM = 350.0
        v.H_SI_BUF_NM = 150.0
    elif profile == "lite":
        v.H_AIR_NM = 450.0
        v.H_PML_NM = 250.0
        v.H_SI_BUF_NM = 100.0
    else:
        raise ValueError(profile)


def write_manifest(path: Path, rows: list[dict[str, str | float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "design_id",
        "family",
        "px_nm",
        "py_nm",
        "feature_min_nm",
        "profile",
        "mph_980",
        "mph_451",
        "status",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(args) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    apply_profile(args.profile)
    v.MESH_AUTO_SIZE = args.mesh_size
    v.NEIGS = args.neigs
    designs = chosen_designs()
    if args.family:
        designs = [d for d in designs if d.family == args.family]
    if args.limit:
        designs = designs[: args.limit]
    print(
        f"[profile] {args.profile}: h_air={v.H_AIR_NM} h_pml={v.H_PML_NM} h_si_buf={v.H_SI_BUF_NM} "
        f"h_sio2={v.H_SIO2_NM} mesh={v.MESH_AUTO_SIZE} neigs={v.NEIGS}",
        flush=True,
    )
    require_memory(args.min_free_gb, "before COMSOL client")
    mph = v.import_mph()
    client = mph.Client(cores=args.cores)
    rows: list[dict[str, str | float]] = []
    try:
        for i, design in enumerate(designs, 1):
            require_memory(args.min_free_gb, f"before geometry {design.design_id}")
            safe = design.design_id.replace("/", "_")
            path980 = OUT / f"{safe}_{args.profile}_geometry_980.mph"
            path451 = OUT / f"{safe}_{args.profile}_geometry_451.mph"
            status = "ok"
            try:
                for target, path in [(v.TARGET_980_NM, path980), (v.TARGET_451_NM, path451)]:
                    if path.exists() and not args.force:
                        print(f"[skip-file] {path.name}", flush=True)
                        continue
                    print(f"[geometry] {i}/{len(designs)} {design.family} {design.design_id} target={target:.0f}", flush=True)
                    model = v.build_model(client, design, target, f"geomonly_{safe}_{int(target)}")
                    try:
                        model.save(str(path))
                    finally:
                        client.remove(model.name())
                    require_memory(args.min_free_gb, f"after geometry {design.design_id} {target:.0f}")
            except Exception as exc:
                status = f"failed: {exc}"
                print(f"[failed] {design.design_id}: {exc}", flush=True)
            rows.append(
                {
                    "design_id": design.design_id,
                    "family": design.family,
                    "px_nm": design.px_nm,
                    "py_nm": design.py_nm,
                    "feature_min_nm": design.feature_min_nm,
                    "profile": args.profile,
                    "mph_980": str(path980),
                    "mph_451": str(path451),
                    "status": status,
                }
            )
            write_manifest(OUT / "geometry_manifest.csv", rows)
    finally:
        client.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build unsolved but executable 3D geometry-only MPH files for v10 families.")
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--family", default="")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--profile", choices=["full", "balanced", "lite"], default="balanced")
    parser.add_argument("--mesh-size", type=int, default=7)
    parser.add_argument("--neigs", type=int, default=6)
    parser.add_argument("--min-free-gb", type=float, default=6.5)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
