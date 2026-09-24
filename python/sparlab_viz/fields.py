"""Turning SparLab element and nodal fields into matplotlib artists.

Element fields are drawn as a `PolyCollection` over the actual element polygons
rather than as an image, so distorted meshes and sub-meshes render correctly and
no interpolation is invented between cells. Nodal fields are drawn with
`tripcolor` using Gouraud shading, which is the honest representation of a field
that really is defined at nodes.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.tri import Triangulation

from .loaders import Mesh
from .style import (
    GRID,
    INK_MUTED,
    INK_SECONDARY,
    SURFACE,
    symmetric_norm,
)


def element_collection(
    ax,
    mesh: Mesh,
    values: Optional[Sequence[float]] = None,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    norm=None,
    edge_color: Optional[str] = None,
    edge_width: float = 0.0,
    alpha: float = 1.0,
    mask: Optional[Sequence[bool]] = None,
) -> PolyCollection:
    """Draw one value per element as filled polygons.

    `mask` selects which elements are drawn (used to show only the material
    above a density threshold).
    """
    polygons = mesh.polygons
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        polygons = polygons[mask]
        if values is not None:
            values = np.asarray(values, dtype=float)[mask]

    # Triangles are drawn with an edge in their own colour: without it the
    # anti-aliased seams between thousands of small cells show as a texture.
    seamless = edge_color is None and mesh.is_simplex
    collection = PolyCollection(
        [poly for poly in polygons],
        cmap=cmap,
        alpha=alpha,
        edgecolors="face" if seamless else (edge_color if edge_color else "none"),
        linewidths=0.25 if seamless else edge_width,
        antialiased=True,
    )
    if values is not None:
        values = np.asarray(values, dtype=float)
        collection.set_array(values)
        if norm is not None:
            collection.set_norm(norm)
        else:
            collection.set_clim(
                float(values.min()) if vmin is None else vmin,
                float(values.max()) if vmax is None else vmax,
            )
    else:
        collection.set_facecolor("none")
    ax.add_collection(collection)
    return collection


def signed_element_collection(ax, mesh: Mesh, values, **kwargs) -> PolyCollection:
    """Element field with a diverging ramp centred exactly on zero."""
    kwargs.setdefault("cmap", "coolwarm")
    return element_collection(ax, mesh, values, norm=symmetric_norm(values), **kwargs)


def _quad_triangulation(mesh: Mesh, coords: np.ndarray) -> Triangulation:
    """Triangles for nodal shading: a Tri3 mesh's own cells, or each
    quadrilateral split into two."""
    cells = mesh.elements
    if cells.shape[1] == 3:
        triangles = cells
    elif cells.shape[1] == 4:
        triangles = np.vstack([cells[:, [0, 1, 2]], cells[:, [0, 2, 3]]])
    else:
        raise ValueError("nodal shading needs a plane (Quad4 or Tri3) mesh")
    return Triangulation(coords[:, 0], coords[:, 1], triangles)


def nodal_tripcolor(
    ax,
    mesh: Mesh,
    values: Sequence[float],
    cmap: str = "viridis",
    coords: Optional[np.ndarray] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
):
    """Draw a nodal field with Gouraud shading."""
    coords = mesh.nodes if coords is None else np.asarray(coords, dtype=float)
    triangulation = _quad_triangulation(mesh, coords)
    values = np.asarray(values, dtype=float)
    return ax.tripcolor(
        triangulation, values, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax
    )


def deformed_coordinates(mesh: Mesh, ux, uy, scale: float) -> np.ndarray:
    """Nodal coordinates displaced by `scale` times the displacement field."""
    coords = mesh.nodes.copy()
    coords[:, 0] += scale * np.asarray(ux, dtype=float)
    coords[:, 1] += scale * np.asarray(uy, dtype=float)
    return coords


def auto_deformation_scale(mesh: Mesh, magnitude, target_fraction: float = 0.08) -> float:
    """Scale that makes the peak displacement `target_fraction` of the diagonal.

    Deformation plots are always exaggerated; the factor is returned so it can
    be printed on the figure instead of left implicit.
    """
    xmin, xmax, ymin, ymax = mesh.extent
    diagonal = float(np.hypot(xmax - xmin, ymax - ymin))
    peak = float(np.max(np.abs(np.asarray(magnitude, dtype=float))))
    if peak <= 0.0:
        return 1.0
    return target_fraction * diagonal / peak


def deformed_collection(
    ax,
    mesh: Mesh,
    ux,
    uy,
    values: Optional[Sequence[float]] = None,
    scale: float = 1.0,
    **kwargs,
) -> PolyCollection:
    """Element polygons drawn on the displaced configuration."""
    coords = deformed_coordinates(mesh, ux, uy, scale)
    displaced = Mesh(
        nodes=coords,
        elements=mesh.elements,
        element_type=mesh.element_type,
        prescribed=mesh.prescribed,
        load_cases=mesh.load_cases,
    )
    return element_collection(ax, displaced, values, **kwargs)


def domain_fill(ax, mesh: Mesh, color: str = "#eceae2", alpha: float = 1.0,
                zorder: float = 0.5) -> PolyCollection:
    """Fill the whole domain with a flat tint.

    Used on boundary-condition figures, where the point is to show *where* the
    design domain is, not a field over it: a colormap at a constant value gives
    a fill that depends on the map's endpoint and can come out white.
    """
    collection = PolyCollection(
        [poly for poly in mesh.polygons], facecolors=color, edgecolors="none",
        alpha=alpha, zorder=zorder, antialiased=False,
    )
    ax.add_collection(collection)
    return collection


def mesh_outline(ax, mesh: Mesh, coords: Optional[np.ndarray] = None,
                 color: str = INK_MUTED, linewidth: float = 0.8,
                 linestyle: str = "-", label: Optional[str] = None,
                 zorder: float = 2.0):
    """Draw only the boundary of the mesh (edges owned by one element) - the
    outline and any holes - for quadrilateral or triangular cells."""
    coords = mesh.nodes if coords is None else np.asarray(coords, dtype=float)
    cells = np.asarray(mesh.elements, dtype=int)
    k = cells.shape[1]
    starts = cells.reshape(-1)
    ends = np.roll(cells, -1, axis=1).reshape(-1)
    keys = np.stack([np.minimum(starts, ends), np.maximum(starts, ends)], axis=1)
    unique, counts = np.unique(keys, axis=0, return_counts=True)
    boundary = unique[counts == 1]
    del k
    segments = [[coords[a], coords[b]] for a, b in boundary]
    from matplotlib.collections import LineCollection

    collection = LineCollection(
        segments, colors=color, linewidths=linewidth, linestyles=linestyle,
        zorder=zorder, label=label,
    )
    ax.add_collection(collection)
    return collection


def density_outline(ax, mesh: Mesh, density, threshold: float = 0.5, **kwargs):
    """Boundary of the material above `threshold`.

    This is the *interpretation* of a density field as geometry, so any figure
    using it must say which threshold was applied.
    """
    density = np.asarray(density, dtype=float)
    keep = density >= threshold
    if not np.any(keep):
        return None
    sub = Mesh(
        nodes=mesh.nodes,
        elements=mesh.elements[keep],
        element_type=mesh.element_type,
    )
    return mesh_outline(ax, sub, **kwargs)


def add_colorbar(fig, mappable, ax, label: str, pad: float = 0.02,
                 fraction: float = 0.035, orientation: str = "vertical"):
    """Attach a recessive colour bar with a units label."""
    bar = fig.colorbar(
        mappable, ax=ax, pad=pad, fraction=fraction, orientation=orientation
    )
    bar.set_label(label, fontsize=8.5, color=INK_SECONDARY)
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, colors=INK_SECONDARY, length=2)
    return bar




def set_domain_limits(ax, mesh: Mesh, coords: Optional[np.ndarray] = None,
                      margin: float = 0.05, equal: bool = True) -> None:
    """Fit the axes to the domain with a small margin and an equal aspect."""
    coords = mesh.nodes if coords is None else np.asarray(coords, dtype=float)
    xmin, xmax = float(coords[:, 0].min()), float(coords[:, 0].max())
    ymin, ymax = float(coords[:, 1].min()), float(coords[:, 1].max())
    dx = (xmax - xmin) or 1.0
    dy = (ymax - ymin) or 1.0
    ax.set_xlim(xmin - margin * dx, xmax + margin * dx)
    ax.set_ylim(ymin - margin * dy, ymax + margin * dy)
    if equal:
        ax.set_aspect("equal", adjustable="box")


def bare_axes(ax) -> None:
    """Strip the frame from a geometry panel, where axes carry no information."""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def geometry_axes(ax, mesh: Mesh, coords: Optional[np.ndarray] = None,
                  units: bool = True, margin: float = 0.05) -> None:
    """Configure a geometry panel: equal aspect, light grid, metre labels."""
    set_domain_limits(ax, mesh, coords, margin=margin)
    ax.grid(True, color=GRID, linewidth=0.6)
    if units:
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
    ax.set_facecolor(SURFACE)


def structured_grid_shape(mesh: Mesh) -> Optional[Tuple[int, int]]:
    """(nx, ny) when the element centroids form a regular grid, else None.

    Used only for image-style density montages, where a regular grid allows a
    much cheaper `imshow` than one polygon collection per frame.
    """
    centroids = mesh.element_centroids
    xs = np.unique(np.round(centroids[:, 0], 12))
    ys = np.unique(np.round(centroids[:, 1], 12))
    if xs.size * ys.size != centroids.shape[0]:
        return None
    return int(xs.size), int(ys.size)


def density_image(density, shape: Tuple[int, int]) -> np.ndarray:
    """Reshape an element density vector into (ny, nx) image order.

    SparLab numbers structured elements with x fastest, so a C-order reshape to
    (ny, nx) puts row 0 at the bottom of the domain; `origin="lower"` then draws
    it the right way up.
    """
    nx, ny = shape
    return np.asarray(density, dtype=float).reshape(ny, nx)
