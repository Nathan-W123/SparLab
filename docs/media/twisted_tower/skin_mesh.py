"""MITC4 shell mesh of the tower's skin and its SparLab topology deck.

usage: python3 skin_mesh.py <out dir> [n] [nz] [volume] [filter elements]

Four walls of n x nz S4 cells round a square of side L (corners shared) and
an n x n cap on top, 6 mm thick. Cells are written wall by wall (faces
y = 0, x = L, y = L, x = 0, each column-major in (along, up)) and then the cap;
skin_cells.json records each cell's face and (along, up) or (i, j) index so a
density field can be drawn on the same squares. The base edge is clamped; the
cap and the top two rows are kept solid (the twist goes in through them), as
is the bottom row.
"""
import json
import sys

import numpy as np

out = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 32
nz = int(sys.argv[3]) if len(sys.argv) > 3 else 80
vf = float(sys.argv[4]) if len(sys.argv) > 4 else 0.35
rf = float(sys.argv[5]) if len(sys.argv) > 5 else 2.5
L, H, T = 0.1, 0.25, 0.006
h, hz = L / n, H / nz

nodes, index = [], {}


def node(x, y, z):
    key = (round(x / h), round(y / h), round(z / hz))
    if key not in index:
        index[key] = len(nodes) + 1
        nodes.append((x, y, z))
    return index[key]


# Perimeter positions, counter-clockwise seen from above, starting at (0, 0).
def perim(s):                      # s in [0, 4n]
    side, t = divmod(s, n)
    if side == 0: return (t * h, 0.0)
    if side == 1: return (L, t * h)
    if side == 2: return (L - t * h, L)
    return (0.0, L - t * h)


cells, meta = [], []
for face in range(4):
    for a in range(n):
        for k in range(nz):
            s0, s1 = face * n + a, face * n + a + 1
            p0, p1 = perim(s0), perim(s1 % (4 * n))
            # Outward normal: counter-clockwise along the perimeter, then up.
            q = [node(p0[0], p0[1], k * hz), node(p1[0], p1[1], k * hz),
                 node(p1[0], p1[1], (k + 1) * hz), node(p0[0], p0[1], (k + 1) * hz)]
            cells.append(q)
            meta.append(["wall", face, a, k])
for i in range(n):
    for j in range(n):
        q = [node(i * h, j * h, H), node((i + 1) * h, j * h, H), node((i + 1) * h, (j + 1) * h, H), node(i * h, (j + 1) * h, H)]
        cells.append(q)
        meta.append(["cap", -1, i, j])

lines = ["*HEADING", "Tower skin, MITC4", "*NODE"]
lines += [f"{k + 1}, {x:.10g}, {y:.10g}, {z:.10g}" for k, (x, y, z) in enumerate(nodes)]
lines.append("*ELEMENT, TYPE=S4, ELSET=SKIN")
lines += [f"{e + 1}, " + ", ".join(map(str, c)) for e, c in enumerate(cells)]
open(f"{out}/skin.inp", "w").write("\n".join(lines) + "\n")
json.dump({"n": n, "nz": nz, "L": L, "H": H, "T": T, "cells": meta}, open(f"{out}/skin_cells.json", "w"))

top = lambda **k: {"box": dict(zmin=H, **k)}
def twist(sense):
    F = 500.0 * sense
    return [{"force": [0.0, F, 0.0, ], "distribution": "total", "region": top(xmin=L)},
            {"force": [0.0, -F, 0.0], "distribution": "total", "region": top(xmax=0.0)},
            {"force": [-F, 0.0, 0.0], "distribution": "total", "region": top(ymin=L)},
            {"force": [F, 0.0, 0.0], "distribution": "total", "region": top(ymax=0.0)},
            {"force": [0.0, 0.0, -1000.0], "distribution": "total", "region": top()}]
cfg = {
    "name": "tower_skin", "units": "SI (m, N, Pa, kg)",
    "description": "Square tower skin (MITC4) clamped at the base, pressed and twisted both ways at the top",
    "mesh": {"type": "file", "path": "skin.inp"},
    "material": {"name": "Ti-6Al-4V", "youngs_modulus": 113.8e9, "poisson_ratio": 0.342, "density": 4430.0},
    "model": {"thickness": T},
    "boundary_conditions": [{"name": "footing", "fix": ["x", "y", "z", "rx", "ry", "rz"], "region": {"box": {"zmax": 0.0}}}],
    "load_cases": [{"name": "twist_ccw", "weight": 1.0, "point_loads": twist(1.0)},
                   {"name": "twist_cw", "weight": 1.0, "point_loads": twist(-1.0)}],
    "solver": {"linear": {"type": "simplicial_ldlt", "residual_tolerance": 1e-8}},
    "modal": {"enabled": False},
    "topology": {
        "enabled": True, "volume_fraction": vf,
        "simp": {"penalty": 3.0, "emin_ratio": 1e-9},
        "filter": {"type": "density", "radius_elements": rf},
        "passive_regions": [
            {"name": "cap_and_ring", "type": "solid", "region": {"box": {"zmin": H - 2 * hz}}},
            {"name": "footing_ring", "type": "solid", "region": {"box": {"zmax": hz}}}],
        "projection": {"enabled": True, "eta": 0.5, "beta_start": 1.0, "beta_max": 16.0,
                       "beta_factor": 2.0, "beta_interval": 25},
        "optimizer": {"max_iterations": 180, "change_tolerance": 1e-2, "objective_tolerance": 1e-3,
                      "objective_window": 10, "move_limit": 0.1, "history_stride": 2}},
    "output": {"csv": True, "vtk": False, "density_history": True},
}
json.dump(cfg, open(f"{out}/config.json", "w"), indent=1)
print(f"{len(cells)} shells ({4 * n * nz} wall, {n * n} cap), {len(nodes)} nodes")
