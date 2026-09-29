#!/usr/bin/env python3
"""Write the shell mesh of the shell cross-validation decks.

Structured and deterministic (numpy only), an Abaqus / CalculiX ``.inp`` of
S4R cells - which SparLab reads as its MITC4 shell - with named sets for the
supports, loads and sections the deck addresses by name:

* ``shell_box_beam.inp``: a cantilever of square box section, mid-surface
  side 0.1 m and 1 m long along x, 8 cells per wall and 40 along. Its four
  walls meet at 90 degree folds (each cell keeps its own normal there). The
  cells run round the section anticlockwise seen from +x, so every normal
  points out of the box. Sets: ``root`` (the nodes at x = 0), ``tip`` (the
  nodes at x = 1 m), ``tip_corner`` (the tip node at y = z = 0.05 m),
  ``flanges`` (the cells of the walls z = +-0.05 m; the deck gives them
  their own thickness) and ``webs`` (the walls y = +-0.05 m).

Run from the repository root::

    python3 python/scripts/make_shell_meshes.py [--output configs/meshes]
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np


def box_beam(side: float, length: float, per_wall: int, along: int):
    """Nodes (n x 3) and S4 cells (m x 4, 0-based) of the tube, the wall
    index of each cell, and the node positions round the section."""
    h = 0.5 * side
    corners = np.array([[-h, -h], [h, -h], [h, h], [-h, h]])
    around = 4 * per_wall
    ring = np.array([corners[i // per_wall] + (i % per_wall) / per_wall
                     * (corners[(i // per_wall + 1) % 4] - corners[i // per_wall])
                     for i in range(around)])
    nodes = np.array([[length * j / along, p[0], p[1]]
                      for j in range(along + 1) for p in ring])
    cells: List[List[int]] = []
    walls: List[int] = []
    for j in range(along):
        for i in range(around):
            i1 = (i + 1) % around
            cells.append([j * around + i, j * around + i1, (j + 1) * around + i1,
                          (j + 1) * around + i])
            walls.append(i // per_wall)
    return nodes, np.array(cells), np.array(walls)


def write_set(out, keyword: str, name: str, ids: List[int]) -> None:
    out.write(f"*{keyword}, {keyword}={name}\n")
    for k in range(0, len(ids), 12):
        out.write(", ".join(str(i + 1) for i in ids[k:k + 12]) + "\n")


def write_box_beam(path: Path) -> Dict[str, int]:
    side, length = 0.1, 1.0
    nodes, cells, walls = box_beam(side, length, 8, 40)
    with open(path, "w", encoding="utf-8") as out:
        out.write("** SparLab shell cross-validation mesh (python/scripts/make_shell_meshes.py)\n")
        out.write("** square box cantilever, S4R cells (read as SparLab's MITC4)\n")
        out.write("*NODE\n")
        for i, x in enumerate(nodes):
            out.write(f"{i + 1}, {x[0]:.12g}, {x[1]:.12g}, {x[2]:.12g}\n")
        out.write("*ELEMENT, TYPE=S4R, ELSET=ALL\n")
        for e, c in enumerate(cells):
            out.write(f"{e + 1}, " + ", ".join(str(n + 1) for n in c) + "\n")
        tol = 1e-12
        write_set(out, "NSET", "root", [i for i, x in enumerate(nodes) if abs(x[0]) < tol])
        write_set(out, "NSET", "tip", [i for i, x in enumerate(nodes) if abs(x[0] - length) < tol])
        write_set(out, "NSET", "tip_corner",
                  [i for i, x in enumerate(nodes)
                   if abs(x[0] - length) < tol and abs(x[1] - 0.5 * side) < tol
                   and abs(x[2] - 0.5 * side) < tol])
        # Walls 0 and 2 are z = -h and z = +h; 1 and 3 are y = +h and y = -h.
        write_set(out, "ELSET", "flanges", [e for e, w in enumerate(walls) if w in (0, 2)])
        write_set(out, "ELSET", "webs", [e for e, w in enumerate(walls) if w in (1, 3)])
    return {"nodes": len(nodes), "cells": len(cells)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output", default="configs/meshes")
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    counts = write_box_beam(out / "shell_box_beam.inp")
    print(f"wrote {out / 'shell_box_beam.inp'}: {counts['nodes']} nodes, {counts['cells']} S4R cells")


if __name__ == "__main__":
    main()
