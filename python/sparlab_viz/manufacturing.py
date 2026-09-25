"""Figures for the manufacturing-oriented features: the robust formulation
with its length-scale check, and the additive-manufacturing overhang filter
with its overhang check.

The numbers on every figure come from the runs' `summary.json`. The one thing
computed here is which elements the overhang check counted as unsupported, so
they can be coloured; the count is then compared with the run's own and a
mismatch is an error, not a figure.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from . import fields as fld
from . import solid as sd
from . import style as st
from .loaders import CaseResults

#: Colour of an unsupported solid element: outside the density ramp, and
#: carried by a legend entry so it is never colour alone.
UNSUPPORTED = "#e34948"
SUPPORTED = "#cccac2"


# ---------------------------------------------------------------------------
# Structured-grid helpers
# ---------------------------------------------------------------------------
def grid_indices(table) -> Tuple[np.ndarray, Tuple[int, ...]]:
    """Integer cell indices (i, j[, k]) of each element of a structured grid,
    from the element centroids in density_final.csv."""
    columns = [c for c in ("cx[m]", "cy[m]", "cz[m]") if c in table.columns]
    idx, shape = [], []
    for c in columns:
        values = table[c].to_numpy()
        levels = np.unique(np.round(values, 12))
        idx.append(np.searchsorted(levels, np.round(values, 12)))
        shape.append(levels.size)
    return np.stack(idx, axis=1), tuple(shape)


def unsupported_mask(table, density: np.ndarray, direction: str,
                     threshold: float) -> np.ndarray:
    """Solid elements (density >= threshold) off the build plate with no
    solid element in their support stencil: the element below it and, in the
    same layer below, its face neighbours (3 in 2-D, a cross of 5 in 3-D)."""
    index, shape = grid_indices(table)
    dim = len(shape)
    sense = -1 if direction.startswith("-") else 1
    axis = "xyz".index(direction.lstrip("+-"))
    solid = np.zeros(shape, dtype=bool)
    solid[tuple(index.T)] = density >= threshold
    out = np.zeros(len(density), dtype=bool)
    lateral = [k for k in range(dim) if k != axis]
    for e, cell in enumerate(index):
        if not solid[tuple(cell)]:
            continue
        layer = cell[axis] if sense > 0 else shape[axis] - 1 - cell[axis]
        if layer == 0:
            continue
        below = cell.copy()
        below[axis] -= sense
        supports = [below]
        for k in lateral:
            for step in (-1, 1):
                side = below.copy()
                side[k] += step
                if 0 <= side[k] < shape[k]:
                    supports.append(side)
        out[e] = not any(solid[tuple(s)] for s in supports)
    return out


def _overhang_block(case: CaseResults) -> Dict:
    block = case.summary.get("manufacturing_checks", {}).get("overhang")
    if not block:
        raise ValueError(f"{case.directory} has no overhang check")
    return block


def _checked_unsupported(case: CaseResults, table, density) -> np.ndarray:
    block = _overhang_block(case)
    mask = unsupported_mask(table, density, block["build_direction"],
                            float(block.get("threshold", 0.5)))
    expected = int(block["unsupported_elements"])
    if int(mask.sum()) != expected:
        raise ValueError(
            f"{case.directory}: the figure's overhang check finds {int(mask.sum())} "
            f"unsupported elements, the run's {expected}; the two stencils disagree")
    return mask


def _plate_marker(ax, mesh, direction: str) -> None:
    """A thick line along the build plate and an arrow in the build direction."""
    xmin, xmax, ymin, ymax = mesh.extent
    axis = "xyz".index(direction.lstrip("+-"))
    sense = -1 if direction.startswith("-") else 1
    kw = dict(color=st.INK_PRIMARY, linewidth=3.0, solid_capstyle="butt")
    if axis == 1:
        y = ymin if sense > 0 else ymax
        ax.plot([xmin, xmax], [y, y], **kw)
        span = 0.14 * (ymax - ymin)
        ax.annotate("", xy=(xmax + 0.03 * (xmax - xmin), y + sense * span),
                    xytext=(xmax + 0.03 * (xmax - xmin), y),
                    arrowprops=dict(arrowstyle="-|>", color=st.INK_PRIMARY, lw=1.2),
                    annotation_clip=False)
    else:
        x = xmin if sense > 0 else xmax
        ax.plot([x, x], [ymin, ymax], **kw)


# ---------------------------------------------------------------------------
# Robust formulation
# ---------------------------------------------------------------------------
def _three_designs(case: CaseResults):
    result = case.summary.get("optimization_result", {})
    block = result.get("robust") or result.get("erosion_check")
    if not block:
        raise ValueError(f"{case.directory} records neither a robust run nor an "
                         "erosion check")
    table = case.density()
    if table is None or "eroded_density[-]" not in table.columns:
        raise FileNotFoundError(f"{case.directory} has no eroded/dilated densities")
    designs = {d["design"].split(" ")[0]: d for d in block["designs"]}
    fields = {"eroded": table["eroded_density[-]"].to_numpy(),
              "intermediate": table["physical_density[-]"].to_numpy(),
              "dilated": table["dilated_density[-]"].to_numpy()}
    return designs, fields


def _length_scale_text(case: CaseResults) -> str:
    block = case.summary.get("manufacturing_checks", {}).get("length_scale")
    if not block:
        return ""
    cell = float(case.summary.get("mesh", {}).get("mean_element_size_m", 0.0)) or 1.0
    solid = block["solid_min_size_m"] / cell
    void = block["void_min_size_m"] / cell
    more = " (scan cap)" if block.get("solid_bound_reached_scan_cap") else ""
    return (f"smallest member ~{solid:g} cells{more}, narrowest gap ~{void:g} cells")


def plot_robust_comparison(plain: CaseResults, robust: CaseResults, path: str) -> str:
    """Eroded, blueprint and dilated designs of a plain (erosion-checked) and a
    robust run on one grid, with their compliances and length scales."""
    if plain.dim != 2 or robust.dim != 2:
        raise ValueError("the robust comparison figure is drawn for plane cases")
    runs = [("without the robust formulation", plain),
            ("robust formulation", robust)]
    data = [_three_designs(case) for _, case in runs]
    mesh = plain.mesh
    xmin, xmax, ymin, ymax = mesh.extent
    aspect = (ymax - ymin) / max(xmax - xmin, 1.0e-12)
    field = float(np.clip(3.55 * aspect + 0.45, 0.9, 3.0))
    fig, axes = st.figure(7.4, 3 * field + 1.9, nrows=3, ncols=2)
    rows = [("eroded", "eroded (part comes out thinner)"),
            ("intermediate", "blueprint (as drawn)"),
            ("dilated", "dilated (part comes out thicker)")]
    mappable = None
    for col, ((label, case), (designs, fields)) in enumerate(zip(runs, data)):
        blueprint_c = float(designs["intermediate"]["compliance_J"])
        for row, (key, text) in enumerate(rows):
            ax = axes[row, col]
            mappable = fld.element_collection(ax, case.mesh, fields[key],
                                              cmap=st.DENSITY_CMAP_SURFACE, vmin=0.0, vmax=1.0)
            fld.mesh_outline(ax, case.mesh, color=st.INK_MUTED, linewidth=0.6)
            fld.set_domain_limits(ax, case.mesh, margin=0.02)
            fld.bare_axes(ax)
            d = designs[key]
            c = float(d["compliance_J"])
            change = "" if key == "intermediate" else f" ({100.0 * (c / blueprint_c - 1.0):+.1f} %)"
            head = f"{label}\n" if row == 0 else ""
            ax.set_title(f"{head}{text}, eta {float(d['eta']):g}:\ncompliance "
                         f"{st.format_si(c)} J{change}, volume {float(d['volume_fraction']):.3f}",
                         loc="left", fontsize=8.2)
        axes[2, col].annotate(_length_scale_text(case), xy=(0.0, -0.08),
                              xycoords="axes fraction", va="top", fontsize=8,
                              color=st.INK_SECONDARY)
    fld.add_colorbar(fig, mappable, list(axes[:, 1]), "physical density [-]")
    st.figure_title(
        fig, f"{robust.name}: a design that survives erosion",
        "one filtered field projected at three thresholds (Wang, Lazarov and Sigmund "
        "2011); the robust run minimised the eroded design's compliance and exported "
        "the blueprint, the plain run was evaluated at the same thresholds once at "
        "the end (projection.erosion_check)",
    )
    st.annotate_note(
        fig,
        "Erosion and dilation stand for a uniform manufacturing error of the whole "
        "part; they are not a model of a particular process. The length scales are "
        "measured on the thresholded blueprint by morphological opening and closing "
        "(to half a cell), a measurement of the design, not a guarantee.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Overhang filter
# ---------------------------------------------------------------------------
def plot_overhang_comparison(designs: Sequence[Tuple[str, CaseResults]], path: str,
                             name: Optional[str] = None) -> str:
    """Plane designs with their unsupported elements marked and the build plate
    drawn, one row per design. `name` titles the figure."""
    count = len(designs)
    mesh = designs[0][1].mesh
    xmin, xmax, ymin, ymax = mesh.extent
    aspect = (ymax - ymin) / max(xmax - xmin, 1.0e-12)
    field = float(np.clip(6.6 * aspect + 0.55, 1.0, 3.2))
    fig, axes = st.figure(7.4, count * field + 1.6, nrows=count, ncols=1)
    axes = np.atleast_1d(axes)
    for ax, (label, case) in zip(axes, designs):
        table = case.density()
        if table is None:
            raise FileNotFoundError(f"{case.directory} has no density_final.csv")
        density = table["physical_density[-]"].to_numpy()
        block = _overhang_block(case)
        mask = _checked_unsupported(case, table, density)
        fld.element_collection(ax, case.mesh, density, cmap=st.DENSITY_CMAP_SURFACE,
                               vmin=0.0, vmax=1.0)
        if mask.any():
            ax.add_collection(PolyCollection(list(case.mesh.polygons[mask]),
                                             facecolors=UNSUPPORTED, edgecolors=UNSUPPORTED,
                                             linewidths=0.3))
        fld.mesh_outline(ax, case.mesh, color=st.INK_MUTED, linewidth=0.6)
        _plate_marker(ax, case.mesh, block["build_direction"])
        fld.set_domain_limits(ax, case.mesh, margin=0.04)
        fld.bare_axes(ax)
        result = case.summary.get("optimization_result", {})
        ax.set_title(
            f"{label}: compliance {st.format_si(result.get('compliance_J', 0.0))} J; "
            f"built {block['build_direction']}: {int(block['unsupported_elements'])} of "
            f"{int(block['solid_elements'])} solid elements unsupported "
            f"({100.0 * float(block['unsupported_fraction']):.2f} % of the solid)",
            loc="left", fontsize=8.4)
    handles = [Patch(facecolor=UNSUPPORTED, edgecolor=UNSUPPORTED,
                     label="unsupported solid element (density >= 0.5, nothing solid "
                           "below or diagonally below)"),
               Line2D([0], [0], color=st.INK_PRIMARY, linewidth=3.0,
                      label="build plate (arrow: build direction)")]
    axes[-1].legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.03),
                    ncol=1, fontsize=8, frameon=False)
    first = designs[0][1]
    st.figure_title(
        fig, f"{name or first.name}: printability under a 45-degree overhang limit",
        "the overhang filter (Langelaar 2016) between the density filter and the "
        "projection lets an element hold no more material than the smooth maximum of "
        "its supports in the layer below; the check counts the thresholded design's "
        "solid elements with no solid support",
    )
    st.annotate_note(
        fig,
        "Only the overhang rule is modelled - square cells, supports directly or "
        "diagonally below, the plate as the first layer. Support removal, thermal "
        "distortion, minimum wall thickness and anisotropic material are not, so "
        "'0 unsupported' means printable under this rule, nothing more.",
    )
    return st.save_figure(fig, path)


def plot_overhang_comparison_3d(designs: Sequence[Tuple[str, CaseResults]], path: str,
                                name: Optional[str] = None) -> str:
    """Solid designs, their unsupported elements' faces in red, drawn with the
    build direction up. `name` titles the figure."""
    count = len(designs)
    fig, axes = st.figure(8.0, 4.6, nrows=1, ncols=count, subplot_kw={"projection": "3d"})
    axes = np.atleast_1d(axes)
    for ax, (label, case) in zip(axes, designs):
        table = case.density()
        if table is None:
            raise FileNotFoundError(f"{case.directory} has no density_final.csv")
        density = table["physical_density[-]"].to_numpy()
        block = _overhang_block(case)
        threshold = float(block.get("threshold", 0.5))
        mask = _checked_unsupported(case, table, density)
        mesh = case.mesh
        solid = density >= threshold
        faces, owners = mesh.boundary_faces(solid, return_owners=True)
        sd.surface(ax, mesh.nodes[faces], mask[owners].astype(float),
                   cmap=ListedColormap([SUPPORTED, UNSUPPORTED]), vmin=0.0, vmax=1.0,
                   edge_color=(0.0, 0.0, 0.0, 0.10), linewidth=0.15)
        direction = block["build_direction"]
        view = dict(sd.VIEW, vertical_axis=direction.lstrip("+-"))
        if direction.startswith("-"):
            view["elev"] = -view["elev"]
        sd.set_solid_axes(ax, mesh.nodes, margin=0.06, view=view)
        result = case.summary.get("optimization_result", {})
        ax.set_title(
            f"{label}\ncompliance {st.format_si(result.get('compliance_J', 0.0))} J\n"
            f"{int(block['unsupported_elements'])} of {int(block['solid_elements'])} solid "
            f"elements unsupported",
            loc="left", fontsize=8.4)
    handles = [Patch(facecolor=UNSUPPORTED, label="faces of unsupported elements"),
               Patch(facecolor=SUPPORTED, label="supported material")]
    fig.legend(handles=handles, loc="outside lower left", fontsize=8, frameon=False)
    first = designs[0][1]
    direction = _overhang_block(first)["build_direction"]
    st.figure_title(
        fig, f"{name or first.name}: printability in 3-D",
        f"built {direction}, drawn with the build direction up (the plate at the "
        "bottom). Supports of an element: the element below it and that element's four "
        "face neighbours in the layer, a 45-degree limit in both lateral directions on "
        "cubic cells; the surface is the thresholded design (density >= 0.5)",
    )
    return st.save_figure(fig, path)
