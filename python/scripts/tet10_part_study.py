#!/usr/bin/env python3
"""Linear against quadratic tetrahedra on a part meshed from CAD.

usage:
    python3 python/scripts/tet10_part_study.py [--output results/tet10_part_study]
        [--bin build/bin] [--sizes 8 6 4 3] [--tet4-sizes 2 1.5]

The engine mount of ``make_meshes.py`` (a slab with four bolt holes and an
upright with a pin hole, drawn in millimetres) is meshed by Gmsh at each
element size of ``--sizes`` twice: with linear tetrahedra, and with quadratic
ones whose edge nodes Gmsh places on the curved hole surfaces. Each mesh is
solved by ``sparlab_solve`` as

* ``Tet4``            - the linear mesh as it is,
* ``Tet10 (straight)`` - the linear mesh elevated in SparLab (``mesh.order:
  2``): an edge node at the midpoint of every straight edge, the same
  faceted geometry as the Tet4 mesh,
* ``Tet10 (curved)``   - the quadratic Gmsh mesh,

plus linear meshes at the extra sizes of ``--tet4-sizes``, to see how far a
linear mesh has to be refined to catch up. Loads: the bolt-hole surfaces are
clamped and the pin-hole surface carries a uniform traction - vertical, and
lateral in a second load case - scaled so that its resultant on the exact
cylinder is 12 kN (5 kN); each element integrates it consistently, so the
problem is the same continuous one on every mesh.

The mesh sequence refines uniformly: Gmsh's curvature rule (elements per
full circle, which sets the size on the holes) scales with the element size,
12 at 4 mm - the committed benchmark mesh - so the size on the bolt holes
shrinks with h instead of staying at 2 pi r / 12 once h drops below it.

The measured quantity is the compliance f.u of each load case, the energy
norm the displacement solution converges in. A displacement-based element is
too stiff: its compliance approaches the exact value from below, and the
shortfall is the "reads too stiff" error of a coarse part mesh. The reference
is the finest curved Tet10 compliance; with three curved Tet10 sizes in a
geometric sequence the script also records a Richardson estimate of the
limit, and reports every value against both.

Output: ``tet10_part_study.csv`` (one row per mesh), ``summary.json``, and
every run's own result directory under ``runs/``. Nothing is invented: each
number is read from the summary a run wrote.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import _bootstrap  # noqa: F401

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

import make_meshes  # noqa: E402  (after the path setup)

#: Pin-hole geometry of make_meshes.engine_mount_3d [mm].
PIN_RADIUS_MM = 10.0
PIN_LENGTH_MM = 60.0
#: Exact volume of the part [m^3]: slab + upright - bolt holes - pin hole.
EXACT_VOLUME_M3 = (160 * 60 * 12 + 80 * 60 * 78 - 4 * math.pi * 5.0 ** 2 * 12
                   - math.pi * PIN_RADIUS_MM ** 2 * PIN_LENGTH_MM) * 1e-9
#: Resultants of the two load cases [N].
LOADS = {"vertical": (0.0, 0.0, -12000.0), "lateral": (5000.0, 0.0, 0.0)}


def curvature_for(size_mm: float) -> float:
    """Elements per full circle on the holes: 12 at 4 mm, scaled with 1/h."""
    return 12.0 * 4.0 / size_mm


def pin_area_m2() -> float:
    return 2.0 * math.pi * PIN_RADIUS_MM * PIN_LENGTH_MM * 1e-6


def deck(name: str, mesh_file: str, order: Optional[int]) -> Dict:
    area = pin_area_m2()
    mesh = {"type": "file", "path": mesh_file, "scale": 0.001}
    if order is not None:
        mesh["order"] = order
    return {
        "name": name,
        "description": "engine mount, static analysis for the Tet4 / Tet10 study",
        "units": "SI (m, N, Pa, kg)",
        "mesh": mesh,
        "material": {"name": "Al 7075-T6", "youngs_modulus": 71.7e9,
                     "poisson_ratio": 0.33, "density": 2810.0},
        "boundary_conditions": [
            {"name": "bolts", "fix": ["x", "y", "z"], "region": {"group": "bolt_holes"}}],
        "load_cases": [
            {"name": case,
             "tractions": [{"name": "pin",
                            "traction": [f / area for f in force],
                            "region": {"group": "pin_hole"}}]}
            for case, force in LOADS.items()],
        "solver": {"linear": {"type": "auto", "iterative_tolerance": 1e-11,
                              "residual_tolerance": 1e-8}},
        "modal": {"enabled": False},
        "output": {"csv": False, "vtk": False},
    }


def run_case(binary: Path, directory: Path, name: str, mesh_file: Path,
             order: Optional[int]) -> Dict:
    directory.mkdir(parents=True, exist_ok=True)
    deck_path = directory / "deck.json"
    deck_path.write_text(json.dumps(deck(name, os.path.relpath(mesh_file, directory), order),
                                    indent=2) + "\n")
    out = directory / "out"
    start = time.perf_counter()
    proc = subprocess.run([str(binary), "--config", str(deck_path), "--output", str(out),
                           "--verbosity", "warn"],
                          capture_output=True, text=True, check=False)
    wall = time.perf_counter() - start
    (directory / "log.txt").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"sparlab_solve failed on {deck_path}:\n{proc.stdout[-3000:]}"
                           f"\n{proc.stderr[-3000:]}")
    summary = json.loads((out / "summary.json").read_text())
    return {"summary": summary, "wall_seconds": wall}


def richardson(values: List[float], sizes: List[float]) -> Optional[Dict]:
    """Limit and order from three values on sizes h1 > h2 > h3 (any ratio)."""
    if len(values) < 3:
        return None
    (c1, c2, c3), (h1, h2, h3) = values[-3:], sizes[-3:]
    d12, d23 = c2 - c1, c3 - c2
    if d12 == 0.0 or d23 == 0.0 or (d12 > 0) != (d23 > 0):
        return None
    # Solve (c2 - c1) / (c3 - c2) = (h1^p - h2^p) / (h2^p - h3^p) for p.
    target = d12 / d23

    def ratio(p: float) -> float:
        return (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p)

    lo, hi = 0.05, 8.0
    if (ratio(lo) - target) * (ratio(hi) - target) > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (ratio(lo) - target) * (ratio(mid) - target) <= 0:
            hi = mid
        else:
            lo = mid
    p = 0.5 * (lo + hi)
    limit = c3 + d23 * h3 ** p / (h2 ** p - h3 ** p)
    return {"limit": limit, "order": p, "from_sizes_mm": [h1, h2, h3]}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", default="results/tet10_part_study")
    parser.add_argument("--bin", default="build/bin")
    parser.add_argument("--sizes", type=float, nargs="+", default=[8.0, 6.0, 4.0, 3.0],
                        help="element sizes of the Tet4 and Tet10 meshes [mm]")
    parser.add_argument("--tet4-sizes", type=float, nargs="*", default=[2.0, 1.5],
                        help="extra, finer linear meshes [mm]")
    args = parser.parse_args(argv)
    out = Path(args.output)
    meshes = out / "meshes"
    meshes.mkdir(parents=True, exist_ok=True)
    binary = Path(args.bin) / "sparlab_solve"
    if not binary.is_file():
        print(f"{binary} not found; build the project first", file=sys.stderr)
        return 2

    sizes = sorted(set(args.sizes), reverse=True)
    extra = sorted(set(args.tet4_sizes) - set(sizes), reverse=True)
    runs = []  # (variant, size, mesh file, order key)
    make_meshes.gmsh.initialize()
    make_meshes.gmsh.option.setNumber("General.Terminal", 0)
    mesh_records = {}
    try:
        for h in sizes + extra:
            for order in ((1, 2) if h in sizes else (1,)):
                path = meshes / f"engine_mount_h{h:g}_o{order}.inp"
                make_meshes.gmsh.clear()
                rec = make_meshes.engine_mount_3d(path, h, order=order,
                                                  curvature=curvature_for(h))
                mesh_records[path.name] = rec
                print(f"  meshed {path.name}: {rec['cells']} cells, {rec['nodes']} nodes")
                if order == 1:
                    runs.append(("Tet4", h, path, None))
                    if h in sizes:
                        runs.append(("Tet10 (straight)", h, path, 2))
                else:
                    runs.append(("Tet10 (curved)", h, path, None))
    finally:
        make_meshes.gmsh.finalize()

    rows = []
    for variant, h, path, order in runs:
        tag = {"Tet4": "tet4", "Tet10 (straight)": "tet10_straight",
               "Tet10 (curved)": "tet10_curved"}[variant]
        name = f"engine_mount_{tag}_h{h:g}"
        print(f"  solving {name} ...", flush=True)
        result = run_case(binary, out / "runs" / name, name, path, order)
        s = result["summary"]
        row = {
            "variant": variant,
            "size_mm": h,
            "element_type": s["mesh"]["element_type"],
            "num_elements": s["mesh"]["num_elements"],
            "num_nodes": s["mesh"]["num_nodes"],
            "num_dofs": s["mesh"]["num_dofs"],
            "domain_volume_m3": s["mesh"]["domain_volume_m3"],
            "volume_error": s["mesh"]["domain_volume_m3"] / EXACT_VOLUME_M3 - 1.0,
            "wall_seconds": result["wall_seconds"],
        }
        for lc in s["load_cases"]:
            row[f"compliance_{lc['name']}_J"] = lc["compliance_J"]
            row[f"max_u_{lc['name']}_m"] = lc["max_displacement_magnitude_m"]
            row[f"resultant_{lc['name']}_N"] = math.sqrt(
                sum(v * v for v in lc["equilibrium"]["applied_force_N"]))
            row[f"solver_iterations_{lc['name']}"] = lc.get("solver_iterations")
        row["linear_solver"] = s["load_cases"][0]["linear_solver"]
        rows.append(row)
        print(f"    {row['num_dofs']} DOFs, compliance vertical "
              f"{row['compliance_vertical_J']:.6g} J, lateral "
              f"{row['compliance_lateral_J']:.6g} J ({row['wall_seconds']:.1f} s)")

    # References: the finest curved Tet10 run, and a Richardson estimate.
    curved = sorted((r for r in rows if r["variant"] == "Tet10 (curved)"),
                    key=lambda r: -r["size_mm"])
    references = {}
    for case in LOADS:
        key = f"compliance_{case}_J"
        finest = curved[-1][key]
        rich = richardson([r[key] for r in curved], [r["size_mm"] for r in curved])
        references[case] = {"finest_tet10_curved_J": finest,
                            "finest_size_mm": curved[-1]["size_mm"],
                            "richardson": rich}
        for r in rows:
            r[f"rel_diff_{case}_vs_finest_tet10"] = r[key] / finest - 1.0
            if rich is not None:
                r[f"rel_diff_{case}_vs_richardson"] = r[key] / rich["limit"] - 1.0

    columns = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in columns:
                columns.append(k)
    with open(out / "tet10_part_study.csv", "w", encoding="utf-8") as handle:
        handle.write(",".join(columns) + "\n")
        for r in rows:
            handle.write(",".join("" if r.get(c) is None else
                                  (f"{r[c]:.12g}" if isinstance(r[c], float) else str(r[c]))
                                  for c in columns) + "\n")
    summary = {
        "about": ("Compliance of the engine mount on linear and quadratic tetrahedral "
                  "meshes of the same CAD geometry (Gmsh), solved by sparlab_solve"),
        "exact_volume_m3": EXACT_VOLUME_M3,
        "pin_traction_resultants_N": LOADS,
        "gmsh_version": make_meshes.gmsh.__version__,
        "meshes": mesh_records,
        "references": references,
        "rows": rows,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"  wrote {out / 'tet10_part_study.csv'} and {out / 'summary.json'}")
    for case, ref in references.items():
        rich = ref["richardson"]
        extra_text = (f", Richardson limit {rich['limit']:.6g} J (order {rich['order']:.2f})"
                      if rich else "")
        print(f"  {case}: finest curved Tet10 {ref['finest_tet10_curved_J']:.6g} J{extra_text}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
