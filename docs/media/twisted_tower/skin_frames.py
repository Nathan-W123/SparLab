"""Contour meshes for the tower-skin animation, from a sparlab_topopt run on
the MITC4 skin of skin_mesh.py.

usage: python3 skin_frames.py <run dir> <skin dir> <out dir> [meshes] [display scale]

Each shell cell's density is painted onto the voxels under it, on a box
grid with the shell grid's spacing; the surface drawn is that density field
cut to a fixed 6.25 mm skin (the four walls meeting at the corners, the cap
on top), so the walls keep their thickness and only the density decides
which parts of them are there. Mesh 0
is the whole skin (a closed tower); the last is the engine's final physical
density at the interpretation threshold 0.5. In between, mesh k shows the
optimiser's skin density at iteration 1 + (n - 1) s^1.6 (snapshots blended),
contoured at the level that keeps a smoothly falling share of the skin: what
it rates lowest goes first, cells it has already emptied ranked by their
distance to the nearest kept cell.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(here))
sys.path.insert(0, "/tmp/claude-0/-home-user-SparLab/d56879da-faaf-5259-b306-3c582599460d/scratchpad/knot")
from box_surface import BoxGrid                   # noqa: E402
from meshio_bin import write_bin, write_index     # noqa: E402

LOW = 0.01
run, skin, out = sys.argv[1], sys.argv[2], sys.argv[3]
count = int(sys.argv[4]) if len(sys.argv) > 4 else 211
scale = float(sys.argv[5]) if len(sys.argv) > 5 else 10.0
W = 6                                   # depth each shell's density is painted to
T_VIS = 0.00625                         # drawn wall thickness [m] (the shell is 6 mm)
os.makedirs(out, exist_ok=True)
info = json.load(open(f"{skin}/skin_cells.json"))
n, nz, L, H = info["n"], info["nz"], info["L"], info["H"]
cells = info["cells"]
g = BoxGrid(n, n, nz, L, L, H)

# Voxel targets of every shell cell (flat indices into the (nz, n, n) grid).
targets = []
for kind, face, a, b in cells:
    vox = []
    for w in range(W):
        if kind == "cap":
            vox.append(((nz - 1 - w), b, a))                    # (k, j, i)
        elif face == 0:
            vox.append((b, w, a))
        elif face == 1:
            vox.append((b, a, n - 1 - w))
        elif face == 2:
            vox.append((b, n - 1 - w, n - 1 - a))
        else:
            vox.append((b, n - 1 - a, w))
    targets.append([k * n * n + j * n + i for k, j, i in vox])
tgt = np.array(targets)                                        # (cells, W)
skin_vox = np.zeros(n * n * nz, bool)
skin_vox[tgt.ravel()] = True


def paint(shell_field):
    v = np.zeros(n * n * nz)
    np.maximum.at(v, tgt.ravel(), np.repeat(np.asarray(shell_field, float), W))
    return v


final = pd.read_csv(f"{run}/density_final.csv")
assert len(final) == len(cells)
cols = [c for c in final.columns if c.split("[")[0] in ("cx", "cy", "cz")]
cz = final[cols[2]].to_numpy()
mine = np.array([(H if k == "cap" else (b + 0.5) * H / nz) for k, f, a, b in cells])
assert np.abs(cz - mine).max() < 1e-9, "cell order differs from skin_cells.json"
rho_final = final["physical_density[-]"].to_numpy()
hist = pd.read_csv(f"{run}/density_history.csv")
iterations = hist["iteration"].to_numpy()
snaps = hist.drop(columns=["iteration"]).to_numpy()
n_iter = int(json.load(open(f"{run}/summary.json"))["optimization_result"]["iterations"])


def ranked(vfield):
    keep = (vfield >= LOW) & skin_vox
    if not keep.any():
        return vfield
    dist = distance_transform_edt(~keep.reshape(nz, n, n), sampling=g.h[::-1]).ravel()
    return np.where(skin_vox & ~keep, LOW - dist, np.where(skin_vox, vfield, -1.0))


# The drawn skin is a fixed slab T_VIS deep under the four walls and the cap
# (signed, scaled to 1 at the outer faces); the density only decides where
# the slab is kept, so the walls never thicken as the level changes.
zz, yy, xx = np.meshgrid(np.arange(nz + 1) * H / nz, np.arange(n + 1) * L / n, np.arange(n + 1) * L / n, indexing="ij")
d_wall = np.minimum(np.minimum(xx, L - xx), np.minimum(yy, L - yy))
slab = np.maximum(T_VIS - d_wall, T_VIS - (H - zz)) / T_VIS
# Nodes on the outer faces: the drawn pattern's area is measured there, on
# the same smoothed field that is drawn, so every frame's level matches it.
eps = 1e-12
face_nodes = ((xx < eps) | (xx > L - eps) | (yy < eps) | (yy > L - eps) | (zz > H - eps)).ravel()
f_final = float(np.mean(g.nodal(paint(rho_final), 1).ravel()[face_nodes] >= 0.5))
D_final = None
meta = {"iterations": n_iter, "shells": len(cells), "frames": []}
t0 = time.time()
for k in range(count):
    s = k / (count - 1)
    if k == 0:
        vfield, level, it, share = np.where(skin_vox, 1.0, 0.0), 0.5, 0.0, 1.0
    elif k == count - 1:
        vfield, level, it = paint(rho_final), 0.5, float(n_iter)
        share = f_final
    else:
        it = 1.0 + (n_iter - 1.0) * s ** 1.6
        i = min(max(int(np.searchsorted(iterations, it, side="right")) - 1, 0), len(iterations) - 2)
        w = min(max((it - iterations[i]) / (iterations[i + 1] - iterations[i]), 0.0), 1.0)
        vfield = ranked(paint((1 - w) * snaps[i] + w * snaps[i + 1]))
        share = 1.0 - (1.0 - f_final) * (s * s * (3 - 2 * s))
        vals = np.sort(g.nodal(vfield, 1).ravel()[face_nodes])[::-1]
        level = vals[min(int(share * len(vals)), len(vals) - 1)] - 1e-9
    D = g.nodal(vfield, 1)
    X, F = g.contour_nodal(np.minimum(D - level, slab), 0.0, up=2, taubin=40,
                           min_area=4 * float(np.min(g.h)) ** 2)
    write_bin(f"{out}/{k:04d}.bin", X, F)
    meta["frames"].append({"iteration": it, "share": share, "triangles": int(len(F))})
    if k % 30 == 0 or k == count - 1:
        print(f"mesh {k:3d}: iteration {it:6.1f}, share {share:.3f}, level {level:+.3f}, {len(F)} triangles ({time.time() - t0:.0f} s)", flush=True)
write_index(out, count, scale, center=[L / 2, L / 2, H / 2], radius=float(np.linalg.norm([L, L, H]) / 2), **meta)
