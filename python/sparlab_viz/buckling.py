"""Figures for the linear buckling check and the buckling constraint.

Every load factor on these figures is read from a run's `summary.json` or its
buckling CSV and VTK files; nothing is recomputed. Two kinds of load factor
appear and are always labelled apart:

* the **SIMP model's** - what the constraint acts on during the optimisation,
  grey material included (for a robust run, the eroded design's);
* the **exported part's** - the density field thresholded at the
  interpretation threshold, largest face-connected group, full material,
  re-analysed after the run (the `buckling_check` block). That is the number
  that says whether the part meets the requirement.

Where the run also checked the part with the non-linear analysis along a
load path (the `nonlinear_check` block), the bracket of its critical load
factor is drawn beside the linear one.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

from . import fields as fld
from . import style as st
from .loaders import CaseResults, VtkGrid


def part_load_factor(case: CaseResults, load_case: Optional[str] = None) -> Optional[float]:
    """Lowest positive load factor of the exported part, or None."""
    block = case.summary.get("buckling_check", {}).get("interpreted_structure")
    return _first_factor(block, load_case)


def solid_load_factor(case: CaseResults, load_case: Optional[str] = None) -> Optional[float]:
    """Lowest positive load factor of the full solid design domain, or None."""
    block = case.summary.get("buckling_check", {}).get("full_solid")
    return _first_factor(block, load_case)


def simp_load_factor(case: CaseResults) -> Optional[float]:
    """Lowest constrained load factor of the SIMP model at the final design."""
    block = case.summary.get("optimization_result", {}).get("buckling")
    if not block:
        return None
    value = block.get("min_load_factor")
    return float(value) if value is not None else None


def nonlinear_bracket(case: CaseResults,
                      load_case: Optional[str] = None) -> Optional[Tuple[float, float]]:
    """The load factors between which the non-linear check of the exported
    part met its critical point, or None (no check, or none met)."""
    block = case.summary.get("nonlinear_check", {})
    for entry in block.get("load_cases", []):
        if load_case is None or entry.get("load_case") == load_case:
            bracket = entry.get("critical_load_factor_bracket")
            if bracket and len(bracket) == 2:
                return float(bracket[0]), float(bracket[1])
            return None
    return None


def plane_mesh(case: CaseResults):
    """The mesh a plane figure draws: a plane mesh itself, a flat shell panel
    in its plane (its buckling modes then leave the plane and show as colour
    alone); anything else is refused."""
    if case.dim == 2:
        return case.mesh
    flat = case.mesh.flat_shell_as_plane()
    if flat is None:
        raise ValueError("the buckling comparison figure is drawn for plane cases and flat "
                         "shell panels")
    return flat


def required_load_factor(case: CaseResults) -> Optional[float]:
    block = case.summary.get("optimization_setup", {}).get("buckling_constraint")
    return float(block["min_load_factor"]) if block else None


def _first_factor(block, load_case: Optional[str]) -> Optional[float]:
    if not block:
        return None
    for entry in block.get("load_cases", []):
        if load_case is None or entry.get("load_case") == load_case:
            factors = entry.get("load_factors") or []
            return float(factors[0]) if factors else None
    return None


def _boundary_segments(grid: VtkGrid, coords: np.ndarray) -> np.ndarray:
    """Edges of a quadrilateral grid that belong to one cell only."""
    count = {}
    for cell in grid.cells:
        n = len(cell)
        for a in range(n):
            edge = tuple(sorted((int(cell[a]), int(cell[(a + 1) % n]))))
            count[edge] = count.get(edge, 0) + 1
    edges = np.array([e for e, c in count.items() if c == 1], dtype=int)
    return coords[edges][:, :, :2]


def _mode_panel(ax, grid: VtkGrid, amplitude: float) -> PolyCollection:
    """The exported part displaced along a buckling mode (max |phi| = 1 in the
    file), coloured by the local mode magnitude, over its undeformed outline."""
    phi = grid.point_data["buckling_mode"]
    magnitude = grid.point_data.get("buckling_mode_magnitude",
                                    np.linalg.norm(phi, axis=1))
    coords = grid.points + amplitude * phi
    polygons = [coords[cell][:, :2] for cell in grid.cells]
    values = np.array([magnitude[cell].mean() for cell in grid.cells])
    collection = PolyCollection(polygons, cmap=st.FIELD_CMAP, edgecolors="face",
                                linewidths=0.2)
    collection.set_array(values)
    collection.set_clim(0.0, 1.0)
    ax.add_collection(collection)
    outline = LineCollection(_boundary_segments(grid, grid.points), colors=st.INK_MUTED,
                             linewidths=0.6, linestyles="--")
    ax.add_collection(outline)
    return collection


def plot_buckling_comparison(designs: Sequence[Tuple[str, CaseResults]], path: str,
                             load_case: Optional[str] = None,
                             name: Optional[str] = None) -> str:
    """Designs side by side, each with the first buckling mode of its exported
    part, and their load factors against the requirement. `name` titles the
    figure (default: the last design's case)."""
    if not designs:
        raise ValueError("no designs to compare")
    meshes = [plane_mesh(case) for _, case in designs]
    shell = designs[0][1].mesh.element_type == "Shell4"
    count = len(designs)
    mesh = meshes[0]
    xmin, xmax, ymin, ymax = mesh.extent
    height = ymax - ymin
    aspect = height / max(xmax - xmin, 1.0e-12)
    panel_w = 7.4 / count
    field = float(np.clip(panel_w * aspect * 0.78, 1.6, 3.6))
    fig = st.figure(7.4, 2 * field + 3.9, nrows=3, ncols=count,
                    height_ratios=[field, field, 2.2])[0]
    axes = np.array(fig.axes).reshape(3, count)
    for ax in axes[2, 1:]:
        ax.remove()
    chart = axes[2, 0]
    chart.remove()
    chart = fig.add_subplot(fig.axes[0].get_gridspec()[2, :])

    mode_mappable = None
    simp, part, labels, brackets = [], [], [], []
    required = None
    # An out-of-plane mode of a flat shell has no in-plane motion to draw.
    amplitude = 0.0 if shell else 0.06 * height
    for k, (label, case) in enumerate(designs):
        table = case.density()
        if table is None:
            raise FileNotFoundError(f"{case.directory} has no density_final.csv")
        density = table["physical_density[-]"].to_numpy()
        ax = axes[0, k]
        fld.element_collection(ax, meshes[k], density, cmap=st.DENSITY_CMAP_SURFACE,
                               vmin=0.0, vmax=1.0)
        fld.mesh_outline(ax, meshes[k], color=st.INK_MUTED, linewidth=0.6)
        fld.set_domain_limits(ax, meshes[k], margin=0.03)
        fld.bare_axes(ax)
        result = case.summary.get("optimization_result", {})
        ax.set_title(f"{textwrap_label(label, 26)}\n"
                     f"compliance {st.format_si(result.get('compliance_J', 0.0))} J",
                     loc="left", fontsize=8.6)

        lam_part = part_load_factor(case, load_case)
        lam_simp = simp_load_factor(case)
        required = required_load_factor(case) or required
        simp.append(lam_simp)
        part.append(lam_part)
        labels.append(label)
        brackets.append(nonlinear_bracket(case, load_case))

        ax = axes[1, k]
        cases = case.summary.get("buckling_check", {}).get("interpreted_structure", {})
        case_name = load_case or (cases.get("load_cases") or [{}])[0].get("load_case", "")
        grid = case.buckling_mode("topology", case_name, 1) if case_name else None
        if grid is None:
            ax.text(0.5, 0.5, "no mode-shape file\n(output.vtk and\noutput.mode_shapes)",
                    ha="center", va="center", fontsize=8, color=st.INK_MUTED,
                    transform=ax.transAxes)
            ax.axis("off")
            continue
        mode_mappable = _mode_panel(ax, grid, amplitude)
        pts = np.vstack([grid.points[:, :2],
                         (grid.points + amplitude * grid.point_data["buckling_mode"])[:, :2]])
        pad = 0.04 * height
        ax.set_xlim(min(xmin, pts[:, 0].min()) - pad, max(xmax, pts[:, 0].max()) + pad)
        ax.set_ylim(min(ymin, pts[:, 1].min()) - pad, max(ymax, pts[:, 1].max()) + pad)
        ax.set_aspect("equal")
        fld.bare_axes(ax)
        ax.set_title("exported part, mode 1" +
                     (f"\nlambda = {lam_part:.3f}" if lam_part is not None else ""),
                     loc="left", fontsize=8.6)
    if mode_mappable is not None:
        fld.add_colorbar(fig, mode_mappable, list(axes[1, :]), "mode magnitude [-]")

    # Load factors against the requirement: the SIMP model the optimiser saw
    # and the part it exported.
    x = np.arange(count)
    chart.plot(x - 0.08, [np.nan if v is None else v for v in simp], "o", markersize=8,
               color=st.series_color(0), label="SIMP model (the constraint's view)")
    chart.plot(x + 0.08, [np.nan if v is None else v for v in part], "s", markersize=8,
               color=st.series_color(1), label="exported part (re-analysed)")
    if any(b is not None for b in brackets):
        # The bracket of the non-linear check's critical load factor, at the
        # part's position: the two states it lies between.
        for i, b in enumerate(brackets):
            if b is None:
                continue
            chart.plot([x[i] + 0.24, x[i] + 0.24], [b[0], b[1]], color=st.series_color(2),
                       linewidth=2.2, solid_capstyle="butt")
            chart.plot([x[i] + 0.24], [0.5 * (b[0] + b[1])], "^", markersize=7,
                       color=st.series_color(2),
                       label="non-linear check: critical bracket" if i == min(
                           j for j, v in enumerate(brackets) if v is not None) else None)
            chart.annotate(f"{0.5 * (b[0] + b[1]):.3f}", (x[i] + 0.24, b[1]),
                           xytext=(6, 4), textcoords="offset points", ha="left",
                           va="bottom", fontsize=7.6, color=st.INK_SECONDARY)
    for i in range(count):
        if simp[i] is not None:
            chart.annotate(f"{simp[i]:.2f}", (x[i] - 0.08, simp[i]), xytext=(-8, 7),
                           textcoords="offset points", ha="right", va="bottom",
                           fontsize=8, color=st.INK_SECONDARY)
        else:
            chart.annotate("no constraint", (x[i] - 0.08, 0.0), xytext=(-8, 6),
                           textcoords="offset points", ha="right", va="bottom",
                           fontsize=7.6, color=st.INK_MUTED)
        if part[i] is not None:
            chart.annotate(f"{part[i]:.2f}", (x[i] + 0.08, part[i]), xytext=(8, -7),
                           textcoords="offset points", ha="left", va="top",
                           fontsize=8, color=st.INK_SECONDARY)
    if required is not None:
        # The label sits just right of the first design, clear of the markers
        # and their labels (the first design is the unconstrained one).
        chart.axhline(required, color=st.INK_SECONDARY, linewidth=1.0, linestyle="--")
        chart.annotate(f"required lambda >= {required:g}",
                       (0.25 if count > 1 else -0.45, required), xytext=(0, 4),
                       textcoords="offset points", ha="center" if count > 1 else "left",
                       va="bottom", fontsize=8, color=st.INK_SECONDARY)
    top = max([v for v in simp + part if v is not None] + [required or 0.0] +
              [b[1] for b in brackets if b is not None])
    chart.set_ylim(0.0, 1.45 * top)
    chart.set_xlim(-0.5, count - 0.5)
    chart.set_xticks(x)
    chart.set_xticklabels([textwrap_label(lbl) for lbl in labels], fontsize=8)
    chart.set_ylabel("lowest load factor [-]")
    st.title(chart, "Lowest buckling load factor",
             "multiples of the design load at which the linear (bifurcation) "
             "analysis predicts buckling")
    st.legend(chart, loc="lower right")

    first = designs[-1][1]
    setup = first.summary.get("optimization_setup", {})
    drawn = ("coloured by the mode's out-of-plane magnitude on the sheet (the mode "
             "leaves the plane)" if shell else
             "drawn at 6 % of the height on its dashed undeformed outline")
    st.figure_title(
        fig, f"{name or first.name}: minimum compliance with and without a buckling "
        "constraint",
        f"same mesh, load, filter and volume fraction "
        f"({setup.get('volume_fraction_target', float('nan')):g}); top row: blueprint "
        "density; middle row: the first buckling mode of the thresholded part "
        f"(density >= 0.5, largest face-connected group, full material), {drawn}",
    )
    note = ("Linear buckling is bifurcation of the ideal geometry under the linear "
            "static stress state: an upper bound for a real, imperfect part, with no "
            "post-buckling. The SIMP model's load factor includes the stress stiffness "
            "of intermediate densities; the exported part has none, which is why the "
            "two differ and why the part's value is the one to judge.")
    if any(b is not None for b in brackets):
        note += (" The non-linear check follows the part along a load path with large "
                 "displacement; its bracket is where the tangent stops being positive "
                 "definite. A member can buckle before that into a stable post-buckled "
                 "state, which the check reports as softening instead.")
    st.annotate_note(fig, note)
    return st.save_figure(fig, path)


def textwrap_label(label: str, width: int = 24) -> str:
    import textwrap

    return "\n".join(textwrap.wrap(label, width))


def plot_buckling_history(case: CaseResults, path: str) -> str:
    """Compliance and lowest SIMP load factor per iteration of a constrained run."""
    history = case.history()
    if history is None or "min_load_factor[-]" not in history.columns:
        raise FileNotFoundError(f"{case.directory} has no buckling history")
    iterations = history["iteration"].to_numpy()
    required = required_load_factor(case)
    robust = "eroded_volume_fraction[-]" in history.columns
    fig, axes = st.stacked_panels(3, width=7.4, panel_height=1.9)

    ax = axes[0]
    ax.plot(iterations, history["compliance[J]"], color=st.series_color(0), linewidth=1.6)
    ax.set_yscale("log")
    ax.set_ylabel("compliance [J]")
    st.title(ax, "Compliance" + (" of the eroded design (the robust objective)"
                                 if robust else ""))

    ax = axes[1]
    ax.plot(iterations, history["min_load_factor[-]"], color=st.series_color(1),
            linewidth=1.6)
    if required is not None:
        ax.axhline(required, color=st.INK_SECONDARY, linewidth=1.0, linestyle="--")
        ax.annotate(f"required {required:g}", (iterations[-1], required), xytext=(0, 4),
                    textcoords="offset points", ha="right", va="bottom", fontsize=8,
                    color=st.INK_SECONDARY)
    ax.set_ylabel("lowest load factor [-]")
    finite = history["min_load_factor[-]"].replace([np.inf, -np.inf], np.nan).dropna()
    if not finite.empty:
        ax.set_ylim(0.0, max(float(finite.max()), required or 0.0) * 1.1)
    st.title(ax, "Lowest buckling load factor of the SIMP model"
             + (" (eroded design)" if robust else ""))

    ax = axes[2]
    if "beta[-]" in history.columns:
        ax.step(iterations, history["beta[-]"], where="post", color=st.series_color(2),
                linewidth=1.6)
        ax.set_yscale("log", base=2)
        ax.set_ylabel("projection beta [-]")
        st.title(ax, "Projection sharpness (continuation)")
    ax.set_xlabel("iteration")
    result = case.summary.get("optimization_result", {})
    st.figure_title(
        fig, f"{case.name}: convergence with the buckling constraint",
        f"{result.get('iterations', 0)} iterations, stop: "
        f"{result.get('stop_reason', 'n/a')}; each iteration solves the buckling "
        "eigenproblem warm-started from the previous subspace",
    )
    return st.save_figure(fig, path)
