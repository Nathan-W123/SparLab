"""Rendering solid (Hex8 and Tet4) meshes: boundary surfaces, feature edges and views.

A solid mesh is drawn through its boundary faces only - the interior is not
visible - and a field is shown on that surface: a nodal field averaged over the
corners of each face (quadrilaterals for Hex8, triangles for Tet4), an element
field taken from the element that owns the face. Faces are flat-shaded from one fixed light direction so the geometry
can be told apart at all; the colour itself is the data, and the shading factor
is the same on every figure so two panels remain comparable.

matplotlib's 3-D axes sort whole polygons by depth and do not clip one artist
against another. That is adequate for the closed, mostly convex surfaces here;
markers and arrows are drawn on top of the surface, which each figure says.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_rgba
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from . import style as st

#: Default camera: x to the right and slightly towards the viewer, z up.
VIEW = {"elev": 24.0, "azim": -58.0}

#: Direction *towards* the light used for flat shading, and the ambient floor,
#: so a face turned away from the light keeps 58 % of its colour.
LIGHT = np.array([-0.35, -0.55, 0.76])
AMBIENT = 0.58


def face_normals(polys: np.ndarray) -> np.ndarray:
    """Unit normals of (n, 4, 3) quads (from their diagonals) or (n, 3, 3)
    triangles (from two edges)."""
    polys = np.asarray(polys, dtype=float)
    if polys.shape[1] == 3:
        normals = np.cross(polys[:, 1] - polys[:, 0], polys[:, 2] - polys[:, 0])
    else:
        normals = np.cross(polys[:, 2] - polys[:, 0], polys[:, 3] - polys[:, 1])
    length = np.linalg.norm(normals, axis=1)
    length[length == 0.0] = 1.0
    return normals / length[:, None]


def shade(colors: np.ndarray, normals: np.ndarray, ambient: float = AMBIENT) -> np.ndarray:
    """Darken RGBA face colours by their angle to the light."""
    light = LIGHT / np.linalg.norm(LIGHT)
    lambert = np.clip(normals @ light, 0.0, 1.0)
    factor = ambient + (1.0 - ambient) * lambert
    out = np.array(colors, dtype=float, copy=True)
    out[:, :3] *= factor[:, None]
    return out


def surface(ax, polys, values: Optional[Sequence[float]] = None,
            cmap: str = st.FIELD_CMAP, vmin: Optional[float] = None,
            vmax: Optional[float] = None, norm=None, color: str = "#d9d7cf",
            alpha: float = 1.0, edge_color=(0.0, 0.0, 0.0, 0.12),
            linewidth: float = 0.15, shading: bool = True):
    """Draw quads or triangles as a shaded surface.

    Returns `(collection, mappable)`; the mappable carries the colour scale
    for a colour bar and is None for a flat-coloured surface.
    """
    polys = np.asarray(polys, dtype=float)
    if polys.size == 0:
        return None, None
    if polys.shape[1] == 3 and edge_color == (0.0, 0.0, 0.0, 0.12):
        # Tetrahedral surfaces have several times as many (smaller) faces as a
        # hexahedral one; the default edge weight would grey the field out.
        edge_color = (0.0, 0.0, 0.0, 0.05)
    normals = face_normals(polys)
    mappable = None
    if values is not None:
        values = np.asarray(values, dtype=float)
        if norm is None:
            norm = Normalize(
                vmin=float(values.min()) if vmin is None else vmin,
                vmax=float(values.max()) if vmax is None else vmax,
            )
        mappable = ScalarMappable(norm=norm, cmap=cmap)
        mappable.set_array(values)
        colors = mappable.to_rgba(values)
    else:
        colors = np.tile(to_rgba(color), (len(polys), 1))
    if shading:
        colors = shade(colors, normals)
    colors[:, 3] = alpha
    collection = Poly3DCollection(
        polys, facecolors=colors,
        edgecolors=edge_color if edge_color is not None else "none",
        linewidths=linewidth, zsort="average",
    )
    ax.add_collection3d(collection)
    return collection, mappable


def feature_edges(faces, coords, angle_deg: float = 30.0) -> np.ndarray:
    """(m, 2, 3) segments along the sharp edges of a quad or triangle surface.

    An edge is kept when its two faces meet at more than `angle_deg`, or when
    it belongs to one face only (an open boundary). On a box this is its twelve
    edges; on a thresholded topology it outlines every step of the surface.
    """
    faces = np.asarray(faces, dtype=int)
    coords = np.asarray(coords, dtype=float)
    if faces.size == 0:
        return np.empty((0, 2, 3))
    normals = face_normals(coords[faces])
    corners = np.stack([faces, np.roll(faces, -1, axis=1)], axis=2).reshape(-1, 2)
    face_of = np.repeat(np.arange(len(faces)), faces.shape[1])
    keys = np.sort(corners, axis=1)
    order = np.lexsort((keys[:, 1], keys[:, 0]))
    keys = keys[order]
    face_of = face_of[order]
    unique, start, counts = np.unique(keys, axis=0, return_index=True,
                                      return_counts=True)
    sharp = counts != 2
    pairs = np.flatnonzero(counts == 2)
    if pairs.size:
        first = face_of[start[pairs]]
        second = face_of[start[pairs] + 1]
        cosines = np.einsum("ij,ij->i", normals[first], normals[second])
        sharp[pairs] = cosines < np.cos(np.radians(angle_deg))
    return coords[unique[sharp]]


def edges(ax, segments, color: str = st.INK_MUTED, linewidth: float = 0.8,
          linestyle: str = "-", label: Optional[str] = None, alpha: float = 1.0):
    """Draw line segments (from `feature_edges`) on a 3-D axes."""
    segments = np.asarray(segments, dtype=float)
    if segments.size == 0:
        return None
    collection = Line3DCollection(
        segments, colors=color, linewidths=linewidth, linestyles=linestyle,
        label=label, alpha=alpha,
    )
    ax.add_collection3d(collection)
    return collection


def set_solid_axes(ax, *coord_sets, margin: float = 0.06, view=None,
                   labels: bool = True, zoom: float = 1.3) -> None:
    """Fit a 3-D axes to the given point sets with a true (equal) aspect.

    `zoom` enlarges the drawing inside the axes box: matplotlib's 3-D axes
    leave a wide margin around the data at zoom 1.
    """
    points = [np.asarray(c, dtype=float).reshape(-1, 3)
              for c in coord_sets if c is not None and len(c)]
    if not points:
        raise ValueError("no coordinates to fit the axes to")
    stacked = np.vstack(points)
    lo = stacked.min(axis=0)
    hi = stacked.max(axis=0)
    span = np.where(hi - lo > 0.0, hi - lo, 1.0)
    lo = lo - margin * span
    hi = hi + margin * span
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(tuple(hi - lo), zoom=zoom)
    ax.view_init(**(VIEW if view is None else view))
    ax.set_facecolor(st.SURFACE)
    longest = float(span.max())
    for axis, extent in zip((ax.xaxis, ax.yaxis, ax.zaxis), span):
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        try:  # private, but the only handle on the pane grid style
            axis._axinfo["grid"]["color"] = st.GRID
            axis._axinfo["grid"]["linewidth"] = 0.5
        except (AttributeError, KeyError):  # pragma: no cover
            pass
        # A short axis of a slender body gets two ticks: four collide.
        axis.set_major_locator(MaxNLocator(4 if extent >= 0.3 * longest else 2))
    ax.tick_params(labelsize=6.8, colors=st.INK_SECONDARY, pad=1)
    if labels:
        ax.set_xlabel("x [m]", fontsize=8, labelpad=5, color=st.INK_SECONDARY)
        ax.set_ylabel("y [m]", fontsize=8, labelpad=5, color=st.INK_SECONDARY)
        ax.set_zlabel("z [m]", fontsize=8, labelpad=5, color=st.INK_SECONDARY)


def bare_solid_axes(ax, *coord_sets, margin: float = 0.04, view=None,
                    zoom: float = 1.15) -> None:
    """A geometry panel without axes, for small multiples."""
    set_solid_axes(ax, *coord_sets, margin=margin, view=view, labels=False, zoom=zoom)
    ax.set_axis_off()


def face_values(faces, nodal) -> np.ndarray:
    """A nodal field averaged over the corners of each face."""
    return np.asarray(nodal, dtype=float)[np.asarray(faces, dtype=int)].mean(axis=1)


def diagonal(coords) -> float:
    """Length of the bounding-box diagonal of a point set [m]."""
    coords = np.asarray(coords, dtype=float)
    return float(np.linalg.norm(coords.max(axis=0) - coords.min(axis=0)))


def auto_scale(coords, magnitude, target_fraction: float = 0.08) -> float:
    """Displacement exaggeration that makes the peak a fraction of the diagonal."""
    peak = float(np.max(np.abs(np.asarray(magnitude, dtype=float))))
    if peak <= 0.0:
        return 1.0
    return target_fraction * diagonal(coords) / peak
