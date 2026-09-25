#!/usr/bin/env python3
"""Generate the unstructured meshes of the real-geometry benchmarks with Gmsh.

Two parts, both drawn in millimetres (the decks read them with
``mesh.scale = 0.001``), both meshed with linear simplices and carrying the
physical groups the decks address by name:

* ``lug_bracket_2d.msh`` (MSH 4.1, triangles): a 200 x 100 mm plate with
  rounded corners, two bolt holes on the left and a lug hole on the right.
  Groups: ``bolt_holes`` and ``load_hole`` (the hole edges), ``hole_rims``
  (an 8 mm ring of material kept around every hole) and ``design`` (the rest).
* ``engine_mount_3d.inp`` (Abaqus / CalculiX input, tetrahedra): a base slab
  with four vertical bolt holes and an upright block with a horizontal pin
  hole. Groups: ``bolt_holes`` and ``pin_hole`` (the cylindrical hole faces),
  ``rings`` (material kept around every hole) and ``design``.

Meshing runs single-threaded with fixed options, so the same Gmsh version
writes the same files; the version and the options are printed and recorded
in the header comment of each file's companion ``.json``. Run from the
repository root::

    python3 python/scripts/make_meshes.py [--size-2d 1.5] [--size-3d 4.0]
        [--output configs/meshes] [--fine]

``--fine`` additionally writes ``engine_mount_3d_fine.inp`` at half the 3-D
element size (several hundred thousand tetrahedra; not committed).
``--tet10`` additionally writes ``engine_mount_3d_tet10.inp``: the same part
meshed with quadratic tetrahedra (C3D10) at ``--size-tet10`` (6 mm), their
edge nodes placed on the curved hole surfaces by Gmsh.

``engine_mount_3d(path, size, order)`` is also what
``tet10_part_study.py`` calls to mesh the part at several sizes and orders.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

try:
    import gmsh
except ImportError:  # pragma: no cover - reported to the user
    sys.exit("make_meshes.py needs the gmsh Python package: pip install gmsh")


def _options(size: float, dim: int, order: int = 1, curvature: float = 12) -> dict:
    opts = {
        "General.NumThreads": 1,
        "Mesh.ElementOrder": order,
        "Mesh.MeshSizeMin": 0.3 * size,
        "Mesh.MeshSizeMax": size,
        # Elements per full circle along curved edges: the size on the holes.
        "Mesh.MeshSizeFromCurvature": curvature,
        "Mesh.Algorithm": 6,  # Frontal-Delaunay (2-D)
        "Mesh.Algorithm3D": 1,  # Delaunay (3-D)
        "Mesh.Optimize": 1,
        "Mesh.RecombineAll": 0,
        "Mesh.SaveAll": 0,
        "Mesh.Binary": 0,
    }
    if dim == 2:
        opts["Mesh.MshFileVersion"] = 4.1
    if order == 2:
        # Edge nodes on the CAD geometry (curved cells along the holes); the
        # high-order optimiser untangles any cell the curving inverts.
        opts["Mesh.SecondOrderLinear"] = 0
        opts["Mesh.HighOrderOptimize"] = 2
    return opts


def _apply(opts: dict) -> None:
    for key, value in opts.items():
        gmsh.option.setNumber(key, value)


def _near(point, centre, tol=1e-6) -> bool:
    return all(abs(a - b) <= tol for a, b in zip(point, centre))


def lug_bracket_2d(path: Path, size: float) -> dict:
    gmsh.model.add("lug_bracket_2d")
    occ = gmsh.model.occ
    plate = occ.addRectangle(0, 0, 0, 200, 100, roundedRadius=10)
    holes = {"bolt": [(25.0, 20.0, 8.0), (25.0, 80.0, 8.0)], "load": [(170.0, 50.0, 10.0)]}
    all_holes = holes["bolt"] + holes["load"]
    disks = [occ.addDisk(x, y, 0, r, r) for x, y, r in all_holes]
    body, _ = occ.cut([(2, plate)], [(2, d) for d in disks])
    rims = []
    for x, y, r in all_holes:
        outer = occ.addDisk(x, y, 0, r + 8.0, r + 8.0)
        inner = occ.addDisk(x, y, 0, r, r)
        ring, _ = occ.cut([(2, outer)], [(2, inner)])
        rims.extend(ring)
    pieces, _ = occ.fragment(body, rims)
    occ.synchronize()

    rim_surfaces, design_surfaces = [], []
    for dim, tag in gmsh.model.getEntities(2):
        cx, cy, _ = occ.getCenterOfMass(dim, tag)
        in_rim = any(math.hypot(cx - x, cy - y) < r + 8.0 + 1e-6 for x, y, r in all_holes)
        (rim_surfaces if in_rim else design_surfaces).append(tag)
    bolt_curves, load_curves = [], []
    for dim, tag in gmsh.model.getEntities(1):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
        centre = (0.5 * (xmin + xmax), 0.5 * (ymin + ymax))
        radius = 0.5 * (xmax - xmin)
        for x, y, r in holes["bolt"]:
            if _near(centre, (x, y), 1e-4) and abs(radius - r) < 1e-4:
                bolt_curves.append(tag)
        for x, y, r in holes["load"]:
            if _near(centre, (x, y), 1e-4) and abs(radius - r) < 1e-4:
                load_curves.append(tag)
    assert len(bolt_curves) == 2 and len(load_curves) == 1, (bolt_curves, load_curves)
    groups = {
        "bolt_holes": (1, bolt_curves),
        "load_hole": (1, load_curves),
        "hole_rims": (2, rim_surfaces),
        "design": (2, design_surfaces),
    }
    for name, (dim, tags) in groups.items():
        gmsh.model.setPhysicalName(dim, gmsh.model.addPhysicalGroup(dim, tags), name)

    opts = _options(size, 2)
    _apply(opts)
    gmsh.model.mesh.generate(2)
    gmsh.write(str(path))
    return _summary(path, opts, groups, 2)


def engine_mount_3d(path: Path, size: float, order: int = 1, curvature: float = 12) -> dict:
    gmsh.model.add("engine_mount_3d")
    occ = gmsh.model.occ
    base = occ.addBox(0, 0, 0, 160, 60, 12)
    tower = occ.addBox(40, 0, 12, 80, 60, 78)
    domain, _ = occ.fuse([(3, base)], [(3, tower)])
    bolts = [(15.0, 15.0), (15.0, 45.0), (145.0, 15.0), (145.0, 45.0)]
    bolt_r, bolt_ring = 5.0, 11.0
    pin = (80.0, 70.0)  # (x, z), axis along y
    pin_r, pin_ring = 10.0, 18.0

    cutters = [occ.addCylinder(x, y, 0, 0, 0, 12, bolt_r) for x, y in bolts]
    cutters.append(occ.addCylinder(pin[0], 0, pin[1], 0, 60, 0, pin_r))
    body, _ = occ.cut(domain, [(3, c) for c in cutters])
    rings = []
    for x, y in bolts:
        outer = occ.addCylinder(x, y, 0, 0, 0, 12, bolt_ring)
        inner = occ.addCylinder(x, y, 0, 0, 0, 12, bolt_r)
        ring, _ = occ.cut([(3, outer)], [(3, inner)])
        rings.extend(ring)
    outer = occ.addCylinder(pin[0], 0, pin[1], 0, 60, 0, pin_ring)
    inner = occ.addCylinder(pin[0], 0, pin[1], 0, 60, 0, pin_r)
    ring, _ = occ.cut([(3, outer)], [(3, inner)])
    rings.extend(ring)
    occ.fragment(body, rings)
    occ.synchronize()

    ring_volumes, design_volumes = [], []
    for dim, tag in gmsh.model.getEntities(3):
        cx, cy, cz = occ.getCenterOfMass(dim, tag)
        near_bolt = any(math.hypot(cx - x, cy - y) < bolt_ring and cz < 12.0 for x, y in bolts)
        near_pin = math.hypot(cx - pin[0], cz - pin[1]) < pin_ring
        (ring_volumes if near_bolt or near_pin else design_volumes).append(tag)
    bolt_faces, pin_faces = [], []
    for dim, tag in gmsh.model.getEntities(2):
        if gmsh.model.getType(dim, tag) != "Cylinder":
            continue
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
        for x, y in bolts:
            if abs(0.5 * (xmin + xmax) - x) < 1e-3 and abs(0.5 * (ymin + ymax) - y) < 1e-3 \
                    and abs(dx - 2 * bolt_r) < 1e-3:
                bolt_faces.append(tag)
        if abs(0.5 * (xmin + xmax) - pin[0]) < 1e-3 and abs(0.5 * (zmin + zmax) - pin[1]) < 1e-3 \
                and abs(dx - 2 * pin_r) < 1e-3 and dy > 50:
            pin_faces.append(tag)
    assert len(bolt_faces) == 4 and len(pin_faces) >= 1, (bolt_faces, pin_faces)
    groups = {
        "bolt_holes": (2, bolt_faces),
        "pin_hole": (2, pin_faces),
        "rings": (3, ring_volumes),
        "design": (3, design_volumes),
    }
    for name, (dim, tags) in groups.items():
        gmsh.model.setPhysicalName(dim, gmsh.model.addPhysicalGroup(dim, tags), name)

    opts = _options(size, 3, order, curvature)
    _apply(opts)
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
    gmsh.model.mesh.generate(3)
    gmsh.write(str(path))
    # Gmsh writes the output path into the *Heading; keep only the file name
    # so the file does not depend on where it was generated.
    lines = path.read_text().splitlines(keepends=True)
    if len(lines) > 1 and lines[0].strip().lower() == "*heading":
        lines[1] = f" {path.name} (Gmsh {gmsh.__version__}, units: mm)\n"
        path.write_text("".join(lines))
    return _summary(path, opts, groups, 3)


def _summary(path: Path, opts: dict, groups: dict, dim: int) -> dict:
    types, tags, _ = gmsh.model.mesh.getElements(dim)
    cells = sum(len(t) for t in tags)
    nodes = len(gmsh.model.mesh.getNodes()[0])
    return {
        "file": path.name,
        "gmsh_version": gmsh.__version__,
        "options": opts,
        "units": "millimetres (read with mesh.scale = 0.001)",
        "dimension": dim,
        "cells": cells,
        "nodes": nodes,
        "groups": {name: {"dim": d, "entities": t} for name, (d, t) in groups.items()},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", default="configs/meshes")
    parser.add_argument("--size-2d", type=float, default=1.5, help="2-D element size [mm]")
    parser.add_argument("--size-3d", type=float, default=4.0, help="3-D element size [mm]")
    parser.add_argument("--fine", action="store_true",
                        help="also write engine_mount_3d_fine.inp at half the 3-D size")
    parser.add_argument("--tet10", action="store_true",
                        help="also write engine_mount_3d_tet10.inp (quadratic tetrahedra)")
    parser.add_argument("--size-tet10", type=float, default=6.0,
                        help="element size of the quadratic mesh [mm]")
    parser.add_argument("--only-tet10", action="store_true",
                        help="write only engine_mount_3d_tet10.inp")
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    records = []
    try:
        if not args.only_tet10:
            records.append(lug_bracket_2d(out / "lug_bracket_2d.msh", args.size_2d))
            gmsh.clear()
            records.append(engine_mount_3d(out / "engine_mount_3d.inp", args.size_3d))
        if args.fine and not args.only_tet10:
            gmsh.clear()
            records.append(engine_mount_3d(out / "engine_mount_3d_fine.inp", 0.5 * args.size_3d))
        if args.tet10 or args.only_tet10:
            gmsh.clear()
            records.append(engine_mount_3d(out / "engine_mount_3d_tet10.inp", args.size_tet10,
                                           order=2))
    finally:
        gmsh.finalize()
    for rec in records:
        (out / (Path(rec["file"]).stem + ".json")).write_text(json.dumps(rec, indent=2) + "\n")
        print(f"{rec['file']}: {rec['cells']} cells, {rec['nodes']} nodes "
              f"(gmsh {rec['gmsh_version']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
