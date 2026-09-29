#!/usr/bin/env python3
"""Write the two-body meshes of the contact cross-validation decks.

Structured and deterministic (numpy only), each an Abaqus / CalculiX ``.inp``
with a named node set (``*NSET``) for every support, load and contact
surface the decks address by name:

* ``contact_blocks_hex.inp``: two Hex8 blocks, 0.4 x 0.2 x 0.3 m below and
  0.4 x 0.15 x 0.3 m above, touching at y = 0.2 m with non-matching meshes
  (8 x 4 x 6 and 10 x 3 x 5 cells). Sets: ``lower_bottom``, ``lower_top``,
  ``upper_bottom``, ``upper_top``, ``load`` (the upper top with x <= 0.2 m),
  ``sym_x`` and ``sym_z`` (the planes x = 0 and z = 0, both blocks).
* ``contact_blocks_gap_hex.inp``: the same blocks, the upper one 20 um above
  the lower one. Sets as above, and ``lower_sym_x`` (the lower block's plane
  x = 0) and ``punch`` (the upper top with x <= 0.2 m).
* ``contact_cap_q4.inp``: a plane half model (x >= 0) of a cylinder of
  radius 0.5 m, its cap 20 mm wide and 20 mm high meshed with 26 x 10 Q4
  cells whose bottom row follows the circle, touching the top of a block
  20 x 20 mm meshed with 20 x 10 cells at x = 0 - a master surface meshed
  finer than its slave, and not matching it. Sets: ``block_top``,
  ``block_bottom``, ``cap_bottom``, ``cap_top``, ``sym_x``.

Run from the repository root::

    python3 python/scripts/make_contact_meshes.py [--output configs/meshes]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

TOL = 1e-12


def hex_block(nx: int, ny: int, nz: int, box: Sequence[float],
              first: int) -> Tuple[np.ndarray, List[List[int]]]:
    """Nodes and C3D8 cells of a structured block; `box` = (x0, x1, y0, y1,
    z0, z1). Node ids from `first` (0-based)."""
    x0, x1, y0, y1, z0, z1 = box
    xs, ys, zs = np.linspace(x0, x1, nx + 1), np.linspace(y0, y1, ny + 1), np.linspace(z0, z1, nz + 1)
    nodes = np.array([(x, y, z) for z in zs for y in ys for x in xs])

    def idx(i: int, j: int, k: int) -> int:
        return first + i + (nx + 1) * (j + (ny + 1) * k)

    cells = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # C3D8: the face 1-2-3-4 at z_k, counter-clockwise seen from +z
                # (a positive Jacobian), then 5-8 above it.
                cells.append([idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k),
                              idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1),
                              idx(i, j + 1, k + 1)])
    return nodes, cells


def quad_grid(xs: Sequence[float], y_of: callable, ny: int, first: int) -> Tuple[np.ndarray, List[List[int]]]:
    """Nodes and CPE4 cells of a plane grid: columns at `xs`, each spanning
    y_of(x) = (bottom, top) in `ny` equal cells."""
    nx = len(xs) - 1
    nodes = []
    for j in range(ny + 1):
        for x in xs:
            bottom, top = y_of(x)
            nodes.append((x, bottom + (top - bottom) * j / ny, 0.0))

    def idx(i: int, j: int) -> int:
        return first + i + (nx + 1) * j

    cells = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)]
             for j in range(ny) for i in range(nx)]
    return np.array(nodes), cells


def write_inp(path: Path, title: str, nodes: np.ndarray, blocks: List[Tuple[str, str, List[List[int]]]],
              sets: Dict[str, List[int]]) -> None:
    lines = [f"** {title}", "** written by python/scripts/make_contact_meshes.py", "*NODE"]
    lines += [f"{i + 1}, {p[0]:.17g}, {p[1]:.17g}, {p[2]:.17g}" for i, p in enumerate(nodes)]
    number = 1
    for kind, elset, cells in blocks:
        lines.append(f"*ELEMENT, TYPE={kind}, ELSET={elset}")
        for cell in cells:
            lines.append(f"{number}, " + ", ".join(str(n + 1) for n in cell))
            number += 1
    for name, ids in sets.items():
        if not ids:
            raise SystemExit(f"node set {name} of {path.name} is empty")
        lines.append(f"*NSET, NSET={name}")
        for i in range(0, len(ids), 16):
            lines.append(", ".join(str(n + 1) for n in ids[i:i + 16]))
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {path} ({len(nodes)} nodes, {number - 1} cells)")


def select(nodes: np.ndarray, ids: range, test) -> List[int]:
    return [n for n in ids if test(nodes[n])]


def blocks(gap: float) -> Tuple[np.ndarray, List, Dict[str, List[int]]]:
    lx, lz, h1, h2 = 0.4, 0.3, 0.2, 0.15
    low_nodes, low_cells = hex_block(8, 4, 6, (0.0, lx, 0.0, h1, 0.0, lz), 0)
    up_nodes, up_cells = hex_block(10, 3, 5, (0.0, lx, h1 + gap, h1 + gap + h2, 0.0, lz), len(low_nodes))
    nodes = np.vstack([low_nodes, up_nodes])
    low = range(len(low_nodes))
    up = range(len(low_nodes), len(nodes))
    everything = range(len(nodes))
    top = h1 + gap + h2
    sets = {
        "lower_bottom": select(nodes, low, lambda p: abs(p[1]) < TOL),
        "lower_top": select(nodes, low, lambda p: abs(p[1] - h1) < TOL),
        "upper_bottom": select(nodes, up, lambda p: abs(p[1] - h1 - gap) < TOL),
        "upper_top": select(nodes, up, lambda p: abs(p[1] - top) < TOL),
        "sym_x": select(nodes, everything, lambda p: abs(p[0]) < TOL),
        "sym_z": select(nodes, everything, lambda p: abs(p[2]) < TOL),
    }
    half = select(nodes, up, lambda p: abs(p[1] - top) < TOL and p[0] <= 0.2 + TOL)
    if gap > 0.0:
        sets["lower_sym_x"] = select(nodes, low, lambda p: abs(p[0]) < TOL)
        sets["punch"] = half
    else:
        sets["load"] = half
    return nodes, [("C3D8", "LOWER", low_cells), ("C3D8", "UPPER", up_cells)], sets


def cap_on_block() -> Tuple[np.ndarray, List, Dict[str, List[int]]]:
    radius, width, height = 0.5, 0.02, 0.02
    block_nodes, block_cells = quad_grid(np.linspace(0.0, width, 21), lambda x: (-height, 0.0), 10, 0)

    def cap_span(x: float) -> Tuple[float, float]:
        return radius - math.sqrt(radius * radius - x * x), height

    cap_nodes, cap_cells = quad_grid(np.linspace(0.0, width, 27), cap_span, 10, len(block_nodes))
    nodes = np.vstack([block_nodes, cap_nodes])
    block = range(len(block_nodes))
    cap = range(len(block_nodes), len(nodes))
    sets = {
        "block_top": select(nodes, block, lambda p: abs(p[1]) < TOL),
        "block_bottom": select(nodes, block, lambda p: abs(p[1] + height) < TOL),
        "cap_bottom": select(nodes, cap,
                             lambda p: abs(p[1] - cap_span(p[0])[0]) < TOL * radius),
        "cap_top": select(nodes, cap, lambda p: abs(p[1] - height) < TOL),
        "sym_x": select(nodes, range(len(nodes)), lambda p: abs(p[0]) < TOL),
    }
    return nodes, [("CPE4", "BLOCK", block_cells), ("CPE4", "CAP", cap_cells)], sets


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", default="configs/meshes")
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    nodes, cells, sets = blocks(0.0)
    write_inp(out / "contact_blocks_hex.inp", "two touching Hex8 blocks, non-matching meshes",
              nodes, cells, sets)
    nodes, cells, sets = blocks(2.0e-5)
    write_inp(out / "contact_blocks_gap_hex.inp",
              "two Hex8 blocks 20 um apart, non-matching meshes", nodes, cells, sets)
    nodes, cells, sets = cap_on_block()
    write_inp(out / "contact_cap_q4.inp", "a cylinder cap (Q4) touching a block, half model",
              nodes, cells, sets)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
