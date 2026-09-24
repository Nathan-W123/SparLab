"""The figures the SparLab documentation uses.

Every function takes loaded results and a target path, writes one figure, and
returns that path. None of them recompute physics: values come from the CSV
tables and summaries a solver run produced.

Where a figure shows something that is an *interpretation* rather than a
computed field - a deformation scale factor, a density threshold read as
geometry - the figure says so on its face.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import fields as fld
from . import plots3d as p3
from . import style as st
from .loaders import CaseResults


# ---------------------------------------------------------------------------
# Mesh, boundary conditions and loads
# ---------------------------------------------------------------------------
def plot_mesh_and_bcs(case: CaseResults, path: str,
                      max_mesh_lines: int = 120) -> str:
    """Mesh, prescribed DOFs, applied loads and any passive design regions.

    Passive regions are part of the problem definition, not of the answer, so
    they belong on this figure rather than only on the density plot.
    """
    if case.dim == 3:
        return p3.plot_mesh_and_bcs(case, path)
    mesh = case.mesh
    # The panel keeps the domain's aspect and the legend sits outside it, so the
    # height that leaves no slack follows from the domain: a tall figure around
    # a wide domain is a page of white space.
    xmin, xmax, ymin, ymax = mesh.extent
    aspect = (xmax - xmin) / max(ymax - ymin, 1.0e-12)
    fig, ax = st.figure(7.2, min(6.6, max(3.1, 3.4 / aspect + 1.55)))

    # Drawing every element edge is unreadable above a few thousand cells, so
    # the interior grid is shown only when it is legible and the boundary is
    # always drawn.
    nx_ny = fld.structured_grid_shape(mesh)
    show_grid = nx_ny is not None and max(nx_ny) <= max_mesh_lines
    fld.domain_fill(ax, mesh)
    if show_grid:
        fld.element_collection(ax, mesh, None, edge_color=st.GRID, edge_width=0.4)
    fld.mesh_outline(ax, mesh, color=st.INK_SECONDARY, linewidth=1.1)

    # Passive design regions, read from the density table's tag column.
    passive_note = ""
    table = case.density() if case.is_topology_run else None
    if table is not None and "passive_tag[-]" in table:
        tags = table["passive_tag[-]"].to_numpy()
        for value, color, label in (
            (1.0, "#4a3aa7", "passive solid (density pinned at 1)"),
            (2.0, "#e34948", "passive void (density pinned at 0)"),
        ):
            selected = tags == value
            if not selected.any():
                continue
            fld.element_collection(ax, mesh, None, mask=selected)
            # A flat tint rather than a colormap value, so the meaning is
            # categorical rather than a magnitude.
            from matplotlib.collections import PolyCollection

            ax.add_collection(PolyCollection(
                [poly for poly in mesh.polygons[selected]], facecolors=color,
                edgecolors="none", alpha=0.55, zorder=1.5, label=label,
            ))
        counts = [int((tags == 1.0).sum()), int((tags == 2.0).sum())]
        if any(counts):
            passive_note = (
                f" {counts[0]} passive solid and {counts[1]} passive void "
                "elements are pinned and not optimised."
            )

    # Prescribed DOFs, split by component so a roller reads differently from a
    # clamp.
    both = np.intersect1d(mesh.constrained_nodes("x"), mesh.constrained_nodes("y"))
    only_x = np.setdiff1d(mesh.constrained_nodes("x"), both)
    only_y = np.setdiff1d(mesh.constrained_nodes("y"), both)
    for nodes, marker, label, slot in (
        (both, "s", "fixed in x and y", 0),
        (only_x, ">", "fixed in x only", 2),
        (only_y, "^", "fixed in y only", 3),
    ):
        if nodes.size == 0:
            continue
        ax.plot(mesh.nodes[nodes, 0], mesh.nodes[nodes, 1], marker,
                color=st.series_color(slot), markersize=3.4, linestyle="none",
                label=f"{label} ({nodes.size} nodes)", zorder=5)

    # Applied loads. Drawing one arrow per loaded node turns into an unreadable
    # blob on a slender domain, so each load case shows where it is applied
    # (markers) and its resultant (one arrow from the loaded centroid). The
    # arrow length is referenced to the SHORTER domain dimension so the scaling
    # works at any aspect ratio.
    extent = mesh.extent
    reference = min(extent[1] - extent[0], extent[3] - extent[2])
    load_slots = [1, 4, 5, 7]
    for index, name in enumerate(mesh.load_case_names):
        forces = mesh.nodal_forces(name)
        if forces.size == 0:
            continue
        nodes = forces[:, 0].astype(int)
        fx, fy = forces[:, 1], forces[:, 2]
        resultant = np.array([fx.sum(), fy.sum()])
        magnitude = float(np.hypot(*resultant))
        color = st.series_color(load_slots[index % len(load_slots)])
        ax.plot(mesh.nodes[nodes, 0], mesh.nodes[nodes, 1], "o", color=color,
                markersize=3.0, linestyle="none", zorder=6,
                label=f"load case '{name}': {st.format_si(magnitude)} N over "
                      f"{nodes.size} node(s)")
        if magnitude <= 0.0:
            continue
        origin = mesh.nodes[nodes].mean(axis=0)
        direction = resultant / magnitude * 0.22 * reference
        ax.annotate(
            "", xy=(origin[0] + direction[0], origin[1] + direction[1]),
            xytext=(origin[0], origin[1]),
            arrowprops=dict(arrowstyle="-|>", color=color, linewidth=1.8,
                            shrinkA=0.0, shrinkB=0.0),
            zorder=7, annotation_clip=False,
        )

    fld.geometry_axes(ax, mesh, margin=0.14)
    st.title(
        ax, f"{case.name}: mesh and boundary conditions",
        f"{mesh.num_elements} {mesh.element_type} elements, {mesh.num_nodes} nodes, "
        f"{len(mesh.prescribed)} prescribed DOFs",
    )
    st.legend(ax, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    st.annotate_note(
        fig,
        ("Interior element edges are omitted above 120 cells per direction. "
         if not show_grid else "")
        + "The tinted region is the design domain, the dark line is its "
          "boundary, and each load case is shown as markers at the loaded nodes "
          "plus one arrow for its resultant. Arrow length is a fixed fraction of "
          "the domain, not proportional to the load." + passive_note,
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Static fields
# ---------------------------------------------------------------------------
def plot_deformed_shape(case: CaseResults, load_case: str, path: str,
                        density: Optional[np.ndarray] = None) -> str:
    """Undeformed outline against the exaggerated deformed shape."""
    if case.dim == 3:
        return p3.plot_deformed_shape(case, load_case, path, density=density)
    mesh = case.mesh
    disp = case.displacement(load_case)
    ux = disp["ux[m]"].to_numpy()
    uy = disp["uy[m]"].to_numpy()
    magnitude = disp["umag[m]"].to_numpy()
    scale = fld.auto_deformation_scale(mesh, magnitude)

    fig, ax = st.figure(7.2, 4.2)
    fld.mesh_outline(ax, mesh, color=st.INK_MUTED, linewidth=0.9, linestyle="--",
                     label="undeformed")
    coords = fld.deformed_coordinates(mesh, ux, uy, scale)
    nodal = np.asarray(magnitude, dtype=float)
    mappable = fld.nodal_tripcolor(ax, mesh, nodal, cmap=st.FIELD_CMAP, coords=coords)
    if density is not None:
        fld.density_outline(ax, mesh, density, 0.5, color=st.INK_SECONDARY,
                            linewidth=0.5)
    fld.add_colorbar(fig, mappable, ax, "|u| [m]")

    fld.geometry_axes(ax, mesh, coords)
    st.title(
        ax, f"{case.name}: deformed shape, load case '{load_case}'",
        f"displacements exaggerated {scale:.4g}x; peak |u| = "
        f"{st.format_si(float(nodal.max()))} m",
    )
    st.legend(ax, loc="lower right")
    st.annotate_note(
        fig,
        f"The deformed outline is the undeformed mesh plus {scale:.4g} times the "
        "computed displacement. The exaggeration is for visibility only; the "
        "analysis itself is small-strain linear elasticity.",
    )
    return st.save_figure(fig, path)


def plot_displacement_magnitude(case: CaseResults, load_case: str, path: str,
                                density: Optional[np.ndarray] = None) -> str:
    """Nodal displacement magnitude on the undeformed configuration.

    `density` is used by the solid renderer only, which draws the retained
    structure rather than the whole domain; the plane figure shows the field
    over the full domain with the design outline drawn on top elsewhere.
    """
    if case.dim == 3:
        return p3.plot_displacement_magnitude(case, load_case, path, density=density)
    mesh = case.mesh
    disp = case.displacement(load_case)
    magnitude = disp["umag[m]"].to_numpy()

    fig, axes = st.figure(7.4, 5.6, nrows=2, ncols=1)
    apply = [("umag[m]", "|u| [m]", st.FIELD_CMAP, None),
             ("ux[m]", "u_x [m]", st.SIGNED_CMAP, "signed")]
    for ax, (column, label, cmap, kind) in zip(axes, apply):
        values = disp[column].to_numpy()
        if kind == "signed":
            mappable = fld.nodal_tripcolor(ax, mesh, values, cmap=cmap)
            mappable.set_norm(st.symmetric_norm(values))
        else:
            mappable = fld.nodal_tripcolor(ax, mesh, values, cmap=cmap)
        fld.mesh_outline(ax, mesh, color=st.INK_SECONDARY, linewidth=0.7)
        fld.add_colorbar(fig, mappable, ax, label)
        fld.geometry_axes(ax, mesh)
    st.figure_title(
        fig, f"{case.name}: displacement, load case '{load_case}'",
        f"peak |u| = {st.format_si(float(magnitude.max()))} m. "
        "Upper panel: magnitude |u|, sequential scale",
        ax=axes[0],
    )
    # A caption, not a heading: styled like the subtitle so the figure has one
    # title, and short enough not to run past the panel.
    axes[1].set_title("lower panel: horizontal component u_x, diverging scale "
                      "centred on zero", loc="left", fontsize=8.5,
                      fontweight="normal", color=st.INK_SECONDARY)
    return st.save_figure(fig, path)


def plot_stress_fields(case: CaseResults, load_case: str, path: str,
                       density: Optional[np.ndarray] = None,
                       threshold: float = 0.5) -> str:
    """Normal, shear and von Mises stress on one page."""
    if case.dim == 3:
        return p3.plot_stress_fields(case, load_case, path, density=density,
                                     threshold=threshold)
    mesh = case.mesh
    stress = case.stress(load_case)
    mask = None
    note_extra = ""
    if density is not None:
        mask = np.asarray(density, dtype=float) >= threshold
        note_extra = (
            f" Only elements with density >= {threshold} are shown: below that "
            "the SIMP stress is a macroscopic average over a near-void cell and "
            "is not a material stress."
        )

    panels = [
        ("sxx[Pa]", r"$\sigma_{xx}$ [Pa]", st.SIGNED_CMAP, True),
        ("syy[Pa]", r"$\sigma_{yy}$ [Pa]", st.SIGNED_CMAP, True),
        ("sxy[Pa]", r"$\sigma_{xy}$ [Pa]", st.SIGNED_CMAP, True),
        ("von_mises[Pa]", r"von Mises $\sigma_{vm}$ [Pa]", st.FIELD_CMAP, False),
    ]
    # The grid follows the domain. A wide domain (the MBB beam is 3:1) gets one
    # column of short panels; a compact one gets a 2 x 2 block, because stacking
    # four 3:2 panels in a single column leaves most of the page white.
    xmin, xmax, ymin, ymax = mesh.extent
    wide = (xmax - xmin) / max(ymax - ymin, 1.0e-12) >= 2.2
    fig, grid = (st.figure(7.6, 8.8, nrows=4, ncols=1) if wide
                 else st.figure(7.6, 5.6, nrows=2, ncols=2))
    axes = list(np.atleast_1d(grid).ravel())
    for ax, (column, label, cmap, signed) in zip(axes, panels):
        values = stress[column].to_numpy()
        shown = values[mask] if mask is not None else values
        if signed:
            collection = fld.element_collection(
                ax, mesh, values, cmap=cmap, norm=st.symmetric_norm(shown), mask=mask
            )
        else:
            collection = fld.element_collection(
                ax, mesh, values, cmap=cmap, vmin=0.0,
                vmax=float(np.max(shown)) if shown.size else 1.0, mask=mask
            )
        fld.mesh_outline(ax, mesh, color=st.INK_SECONDARY, linewidth=0.7)
        fld.add_colorbar(fig, collection, ax, label)
        fld.geometry_axes(ax, mesh)
    # No panel titles: each colour bar is labelled with its own component, so a
    # title above the panel would repeat it. The subtitle gives the reading
    # order instead, and axes[0]'s title slot carries that subtitle.

    peak = float(stress["von_mises[Pa]"].to_numpy()[mask].max()) if mask is not None \
        else float(stress["von_mises[Pa]"].max())
    st.figure_title(
        fig, f"{case.name}: stress components, load case '{load_case}'",
        f"peak von Mises {st.format_si(peak)} Pa; element values are the average "
        "over the 2x2 stiffness quadrature points. Panels, in order: normal "
        "sigma_xx, normal sigma_yy, shear sigma_xy, von Mises - each labelled "
        "on its own colour bar",
        ax=axes[0],
    )
    st.annotate_note(
        fig,
        "Normal and shear components use a diverging scale with zero exactly at "
        "the neutral midpoint and symmetric limits; von Mises is non-negative "
        "and uses a sequential scale." + note_extra,
    )
    return st.save_figure(fig, path)


def plot_reactions(case: CaseResults, load_case: str, path: str) -> str:
    """Support reactions against the applied load, with the balance printed."""
    if case.dim == 3:
        return p3.plot_reactions(case, load_case, path)
    mesh = case.mesh
    reactions = case.reactions(load_case)
    summary = case.summary
    entry = next(
        (lc for lc in summary.get("load_cases", []) if lc["name"] == load_case), None
    )

    fig, ax = st.figure(7.2, 4.4)
    fld.domain_fill(ax, mesh)
    fld.mesh_outline(ax, mesh, color=st.INK_SECONDARY, linewidth=1.0)

    extent = mesh.extent
    reference = min(extent[1] - extent[0], extent[3] - extent[2])
    rx = reactions["rx[N]"].to_numpy()
    ry = reactions["ry[N]"].to_numpy()
    peak = float(np.max(np.hypot(rx, ry))) or 1.0
    arrow = 0.22 * reference / peak
    ax.quiver(reactions["x[m]"], reactions["y[m]"], rx * arrow, ry * arrow,
              color=st.series_color(0), angles="xy", scale_units="xy", scale=1.0,
              width=0.004, label="support reactions (per node)", zorder=5)

    forces = mesh.nodal_forces(load_case)
    if forces.size:
        nodes = forces[:, 0].astype(int)
        fx, fy = forces[:, 1], forces[:, 2]
        resultant = np.array([fx.sum(), fy.sum()])
        magnitude = float(np.hypot(*resultant)) or 1.0
        origin = mesh.nodes[nodes].mean(axis=0)
        direction = resultant / magnitude * 0.22 * reference
        ax.annotate(
            "", xy=(origin[0] + direction[0], origin[1] + direction[1]),
            xytext=(origin[0], origin[1]),
            arrowprops=dict(arrowstyle="-|>", color=st.series_color(1),
                            linewidth=2.0, shrinkA=0.0, shrinkB=0.0),
            zorder=7, annotation_clip=False,
        )
        ax.plot([], [], "-", color=st.series_color(1), linewidth=2.0,
                label=f"applied resultant {st.format_si(magnitude)} N")

    subtitle = f"{len(reactions)} constrained nodes"
    if entry is not None:
        eq = entry["equilibrium"]
        subtitle = (
            f"applied ({st.format_si(eq['applied_force_N'][0])}, "
            f"{st.format_si(eq['applied_force_N'][1])}) N, "
            f"reactions ({st.format_si(eq['reaction_force_N'][0])}, "
            f"{st.format_si(eq['reaction_force_N'][1])}) N, "
            f"relative force-balance error {eq['relative_force_error']:.2e}"
        )
    fld.geometry_axes(ax, mesh, margin=0.14)
    st.title(ax, f"{case.name}: reactions, load case '{load_case}'", subtitle)
    st.legend(ax, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    st.annotate_note(
        fig,
        "Reaction arrows are per node and scaled to the largest nodal reaction; "
        "the applied load is drawn as a single resultant from the centroid of "
        "the loaded nodes. Arrow length is therefore not comparable between the "
        "two sets - the force balance printed above is the quantitative check.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Modal
# ---------------------------------------------------------------------------
def plot_mode_shapes(case: CaseResults, path: str, tag: str = "",
                     num_modes: int = 6) -> str:
    """Small multiples of the lowest mode shapes."""
    if case.dim == 3:
        return p3.plot_mode_shapes(case, path, tag=tag, num_modes=num_modes)
    modes = case.modes(tag)
    shapes = case.mode_shapes(tag)
    if modes is None or shapes is None:
        raise FileNotFoundError(
            f"{case.directory} has no modal output"
            + (f" for tag '{tag}'" if tag else "")
        )

    # The mode-shape table may belong to a sub-mesh (the interpreted topology),
    # so the mesh is rebuilt from its own coordinates when the counts differ.
    mesh = case.mesh
    if len(shapes) != mesh.num_nodes:
        raise ValueError(
            f"mode_shapes has {len(shapes)} rows but mesh.json has "
            f"{mesh.num_nodes} nodes; the two belong to different meshes"
        )

    count = min(num_modes, len(modes))
    rows = (count + 1) // 2
    # Size each row from the domain aspect ratio: a 10:1 panel does not need the
    # same height as a square one.
    xmin, xmax, ymin, ymax = mesh.extent
    aspect = (ymax - ymin) / max(xmax - xmin, 1e-30)
    row_height = float(np.clip(3.6 * aspect + 0.55, 0.85, 2.6))
    fig, axes = st.figure(7.6, row_height * rows + 0.5, nrows=rows, ncols=2)
    axes = np.atleast_1d(axes).ravel()

    for index in range(count):
        ax = axes[index]
        ux = shapes[f"ux_mode{index}[m]"].to_numpy()
        uy = shapes[f"uy_mode{index}[m]"].to_numpy()
        magnitude = np.hypot(ux, uy)
        scale = fld.auto_deformation_scale(mesh, magnitude, 0.12)
        coords = fld.deformed_coordinates(mesh, ux, uy, scale)
        fld.mesh_outline(ax, mesh, color=st.INK_MUTED, linewidth=0.6, linestyle="--")
        fld.nodal_tripcolor(ax, mesh, magnitude, cmap=st.FIELD_CMAP, coords=coords)
        fld.set_domain_limits(ax, mesh, coords, margin=0.12)
        fld.bare_axes(ax)
        row = modes.iloc[index]
        ax.set_title(
            f"mode {index + 1}: {row['frequency[Hz]']:.1f} Hz", loc="left",
            fontsize=9.5,
        )
    for index in range(count, len(axes)):
        axes[index].axis("off")

    label = f" ({tag})" if tag else ""
    fig.suptitle(
        f"{case.name}: mode shapes{label}", x=0.01, ha="left", fontsize=11,
        fontweight="bold", color=st.INK_PRIMARY,
    )
    st.annotate_note(
        fig,
        "Mode shapes are M-orthonormal (phi^T M phi = 1) and are drawn on an "
        "exaggerated displaced configuration; colour is the local mode-shape "
        "magnitude. Dashed outline is the undeformed structure.",
    )
    return st.save_figure(fig, path)


def plot_modal_comparison(summary: Dict, path: str, case_name: str,
                          dim: int = 2) -> str:
    """Frequencies of the solid domain and the optimised topology, and mass."""
    solid = summary.get("modal_initial_solid")
    topology = summary.get("modal_optimised_topology")
    if solid is None or topology is None:
        raise KeyError("summary.json has no modal comparison block")

    f_solid = np.asarray(solid["frequencies_hz"], dtype=float)
    f_topo = np.asarray(topology["frequencies_hz"], dtype=float)
    count = min(f_solid.size, f_topo.size)
    index = np.arange(1, count + 1)
    width = 0.38

    fig, axes = st.stacked_panels(2, width=7.0, panel_height=2.5)
    ax = axes[0]
    ax.bar(index - width / 2, f_solid[:count], width, color=st.series_color(0),
           label=f"full solid domain ({solid['total_mass_kg']:.3g} kg)")
    ax.bar(index + width / 2, f_topo[:count], width, color=st.series_color(1),
           label=f"optimised topology ({topology['total_mass_kg']:.3g} kg)")
    for x, value in zip(index - width / 2, f_solid[:count]):
        ax.text(x, value, f"{value:.0f}", ha="center", va="bottom", fontsize=7.4,
                color=st.INK_SECONDARY)
    for x, value in zip(index + width / 2, f_topo[:count]):
        ax.text(x, value, f"{value:.0f}", ha="center", va="bottom", fontsize=7.4,
                color=st.INK_SECONDARY)
    ax.set_ylabel("frequency [Hz]")
    st.title(
        ax, f"{case_name}: natural frequencies before and after optimisation",
        "the optimised structure is the density field thresholded at "
        f"{summary.get('solid_interpretation', {}).get('threshold', 0.5)} and "
        "reduced to its largest connected group",
    )
    st.legend(ax, loc="upper left")

    ax = axes[1]
    ratio = f_topo[:count] / f_solid[:count]
    ax.axhline(1.0, color=st.INK_MUTED, linewidth=0.9, linestyle="--")
    ax.bar(index, ratio, 0.5, color=st.series_color(2),
           label="optimised / solid frequency ratio")
    mass_ratio = topology["total_mass_kg"] / solid["total_mass_kg"]
    ax.set_ylabel("frequency ratio [-]")
    ax.set_xlabel("mode number")
    ax.set_xticks(index)
    st.title(
        ax, "frequency ratio at a mass ratio of " f"{mass_ratio:.3f}",
        "above 1 means the optimised design is stiffer per unit mass in that mode",
    )
    st.legend(ax, loc="upper right")
    if dim == 2:
        note = (
            "In this 2-D idealisation K and M both scale linearly with thickness, "
            "so a uniformly thinned plate of the same mass has exactly the "
            "full-solid frequencies. The full-solid bars are therefore also the "
            "equal-mass uniform baseline."
        )
    else:
        note = (
            "A solid has no thickness to scale, so there is no uniform equal-mass "
            "baseline with the full-solid frequencies: the full-solid bars are the "
            "reference structure, and the mass ratio above says how much lighter "
            "the optimised design is."
        )
    st.annotate_note(fig, note)
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Topology optimisation
# ---------------------------------------------------------------------------
def plot_final_topology(case: CaseResults, path: str) -> str:
    """Final density field plus the thresholded interpretation."""
    if case.dim == 3:
        return p3.plot_final_topology(case, path)
    density_table = case.density()
    if density_table is None:
        raise FileNotFoundError(f"{case.directory} has no density_final.csv")
    density = density_table["physical_density[-]"].to_numpy()
    threshold = (
        case.summary.get("solid_interpretation", {}).get("threshold", 0.5)
    )
    result = case.summary.get("optimization_result", {})

    fig, axes = st.figure(7.4, 5.6, nrows=2, ncols=1)
    collection = fld.element_collection(axes[0], case.mesh, density,
                                        cmap=st.DENSITY_CMAP_SURFACE,
                                        vmin=0.0, vmax=1.0)
    fld.mesh_outline(axes[0], case.mesh, color=st.INK_MUTED, linewidth=0.8)
    # One bar for both panels: it is the same density scale, and a bar attached
    # to the upper panel alone would make the two panels different widths.
    fld.add_colorbar(fig, collection, list(axes), "density [-]")
    fld.geometry_axes(axes[0], case.mesh)
    st.figure_title(
        fig, f"{case.name}: optimised density",
        "SIMP design variable field. Compliance "
        f"{st.format_si(result.get('compliance_J', float('nan')))} J, volume "
        f"fraction {result.get('volume_fraction', float('nan')):.4f}, grey level "
        f"{result.get('grey_level', float('nan')):.3f}, "
        f"{result.get('iterations', 0)} iterations",
        ax=axes[0],
    )

    # A solid fill over the retained material only, plus the original domain
    # boundary as a dashed reference.
    mask = density >= threshold
    fld.element_collection(axes[1], case.mesh,
                           np.full(case.mesh.num_elements, 1.0),
                           cmap=st.DENSITY_CMAP_SURFACE, vmin=0.0, vmax=1.0,
                           mask=mask)
    fld.mesh_outline(axes[1], case.mesh, color=st.INK_MUTED, linewidth=0.8,
                     linestyle="--")
    fld.geometry_axes(axes[1], case.mesh)
    interpretation = case.summary.get("solid_interpretation", {})
    st.title(
        axes[1], f"interpretation as solid material at density >= {threshold}",
        f"{interpretation.get('elements_retained', int(mask.sum()))} of "
        f"{case.mesh.num_elements} elements retained in "
        f"{interpretation.get('connected_groups_above_threshold', 1)} connected "
        f"group(s); {st.format_si(interpretation.get('volume_discarded_as_islands_m3', 0.0))}"
        " m^3 discarded as disconnected islands",
    )
    st.annotate_note(
        fig,
        "The upper panel is the SIMP design variable field: a relative material "
        "distribution, not a solid body. The lower panel is one explicit "
        "interpretation of it, and the threshold and connectivity rule are "
        "stated because the answer depends on both.",
    )
    return st.save_figure(fig, path)


def plot_convergence_history(case: CaseResults, path: str) -> str:
    """Compliance, volume fraction, design change and grey level vs iteration.

    Four quantities of different scale, so four stacked panels sharing the x
    axis. A single chart with two y-scales would make the reader guess which
    curve belongs to which axis.
    """
    history = case.history()
    if history is None:
        raise FileNotFoundError(f"{case.directory} has no history.csv")
    result = case.summary.get("optimization_result", {})
    setup = case.summary.get("optimization_setup", {})
    iterations = history["iteration"].to_numpy()
    method = str(setup.get("method", "oc")).lower()
    stress_constrained = "max_stress_ratio[-]" in history and bool(result.get("stress"))

    fig, axes = st.stacked_panels(5 if stress_constrained else 4, width=7.0,
                                  panel_height=1.75)

    ax = axes[0]
    # A linear scale reads better here: the objective falls by a factor of a few
    # over the run, and a log axis only crowds the tick labels.
    ax.plot(iterations, history["compliance[J]"], color=st.series_color(0))
    ax.set_ylabel("compliance [J]")
    ax.set_ylim(bottom=0.0)
    st.limit_ticks(ax, x=8, y=5)
    st.title(
        ax, f"{case.name}: optimisation history",
        f"stopped by '{result.get('stop_reason', 'unknown')}' after "
        f"{result.get('iterations', len(iterations))} iterations; final compliance "
        f"{st.format_si(result.get('compliance_J', float('nan')))} J",
    )
    # One series, so the title names it and no legend box is needed.

    ax = axes[1]
    target = setup.get("volume_fraction_target")
    violation = result.get("volume_constraint_relative_violation", 0.0)
    if target and method == "mma":
        # MMA treats the volume as an explicit constraint that its subproblem
        # satisfies only approximately, so the interesting reading is how far
        # each iterate sits from the target and which iterates overshoot it.
        relative = (history["volume_fraction[-]"].to_numpy() - target) / target
        magnitude = np.maximum(np.abs(relative), 1.0e-16)
        ax.plot(iterations, magnitude, color=st.series_color(2),
                label="|volume - target| / target")
        over = relative > 0.0
        if over.any():
            ax.plot(iterations[over], magnitude[over], "o", color=st.series_color(7),
                    markersize=3.2, label=f"over the target ({int(over.sum())} iterates)")
        feasibility = setup.get("constraint_tolerance")
        if feasibility:
            ax.axhline(feasibility, color=st.INK_MUTED, linewidth=0.9, linestyle="--",
                       label=f"feasibility tolerance {feasibility:g}")
        ax.set_yscale("log")
        ax.set_ylabel("relative distance [-]")
        st.title(
            ax, f"volume constraint, target fraction {target:g}",
            f"met to {abs(violation):.1e} relative at the returned iterate; an "
            "explicit MMA constraint, so iterates may sit on either side of the target",
        )
        st.legend(ax, loc="upper right")
    elif target:
        relative = (history["volume_fraction[-]"].to_numpy() - target) / target
        ax.plot(iterations, relative, color=st.series_color(2))
        ax.axhline(0.0, color=st.INK_MUTED, linewidth=0.9, linestyle="--")
        ax.set_yscale("symlog", linthresh=1.0e-12)
        ax.set_ylabel("relative violation [-]")
        st.limit_ticks(ax, x=8, y=5)
        st.title(
            ax, f"volume constraint, target fraction {target:g}",
            f"met to {abs(violation):.1e} relative at the final iterate, which is "
            "the multiplier bisection tolerance, not a drift",
        )
    else:
        ax.plot(iterations, history["volume_fraction[-]"], color=st.series_color(2))
        ax.set_ylabel("volume fraction [-]")
        st.title(ax, "volume constraint")
    # One series against a zero reference line; the title carries the reading.

    ax = axes[2]
    ax.plot(iterations, history["max_design_change[-]"], color=st.series_color(1),
            label=r"max $|\Delta x|$")
    tol = setup.get("change_tolerance")
    if tol:
        ax.axhline(tol, color=st.INK_MUTED, linewidth=0.9, linestyle="--",
                   label=f"tolerance {tol:g}")
    ax.set_yscale("log")
    ax.set_ylabel("design change [-]")
    st.title(
        ax, "design change per iteration",
        "on a fine mesh this measure stalls above its tolerance because thin "
        "members migrate one cell at a time long after the objective has settled",
    )
    st.legend(ax, loc="best")

    grey_panel = 3
    if stress_constrained:
        ax = axes[3]
        ax.plot(iterations, history["max_stress_ratio[-]"], color=st.series_color(3),
                label="max relaxed von Mises / limit")
        ax.axhline(1.0, color=st.INK_MUTED, linewidth=0.9, linestyle="--",
                   label="stress limit")
        ax.set_ylabel("stress ratio [-]")
        st.limit_ticks(ax, x=8, y=5)
        stress_block = result.get("stress", {})
        st.title(
            ax, "stress constraint",
            "final max relaxed ratio "
            f"{stress_block.get('max_relaxed_stress_ratio', float('nan')):.4f}. MMA "
            "constrains the p-norm aggregate, whose scale is re-fitted to the true "
            "maximum every iteration, so the ratio hovers about 1 while the design "
            "moves and settles below it",
        )
        st.legend(ax, loc="best")
        grey_panel = 4

    ax = axes[grey_panel]
    ax.plot(iterations, history["grey_level[-]"], color=st.series_color(6))
    if "penalty[-]" in history:
        penalty = history["penalty[-]"].to_numpy()
        if penalty.max() > penalty.min():
            changes = np.flatnonzero(np.diff(penalty) > 0) + 1
            for change in changes:
                ax.axvline(iterations[change], color=st.INK_MUTED, linewidth=0.8,
                           linestyle=":")
            ax.text(
                iterations[changes[0]], ax.get_ylim()[1],
                f" SIMP penalty raised (p: {penalty.min():g} -> {penalty.max():g})",
                fontsize=7.4, color=st.INK_MUTED, va="top",
            )
    ax.set_ylabel("grey level [-]")
    ax.set_xlabel("iteration")
    ax.set_ylim(0.0, 1.0)
    st.title(
        ax, "how binary the design is",
        "grey level 4/n sum rho (1 - rho): 0 = pure 0/1, 1 = every element at 0.5",
    )

    if method == "mma":
        why = ("MMA minimises a convex approximation with moving asymptotes and a "
               "move limit, and a stress constraint's p-norm scale changes between "
               "iterations")
    else:
        why = ("optimality criteria is a fixed-point update with a move limit, and "
               "each continuation step raises the SIMP penalty, which raises "
               "compliance at fixed density")
    st.annotate_note(
        fig,
        f"Compliance is not guaranteed to fall monotonically: {why}. Every value is "
        "plotted rather than smoothed.",
    )
    return st.save_figure(fig, path)


def plot_density_evolution(case: CaseResults, path: str, frames: int = 8) -> str:
    """Small multiples of the density field as the optimisation proceeds."""
    if case.dim == 3:
        return p3.plot_density_evolution(case, path, frames=min(frames, 6))
    table = case.density_history()
    if table is None:
        raise FileNotFoundError(f"{case.directory} has no density_history.csv")
    columns = [c for c in table.columns if c.startswith("rho_")]
    iterations = table["iteration"].to_numpy()
    picks = np.unique(np.linspace(0, len(table) - 1, frames).round().astype(int))

    rows = (len(picks) + 1) // 2
    fig, axes = st.figure(7.6, 1.55 * rows, nrows=rows, ncols=2)
    axes = np.atleast_1d(axes).ravel()
    shape = fld.structured_grid_shape(case.mesh)
    extent = case.mesh.extent

    for slot, pick in enumerate(picks):
        ax = axes[slot]
        density = table.iloc[pick][columns].to_numpy(dtype=float)
        if shape is not None:
            ax.imshow(
                fld.density_image(density, shape), cmap=st.DENSITY_CMAP_SURFACE,
                vmin=0.0, vmax=1.0, origin="lower",
                extent=(extent[0], extent[1], extent[2], extent[3]),
                interpolation="nearest", aspect="equal",
            )
        else:
            fld.element_collection(ax, case.mesh, density,
                                   cmap=st.DENSITY_CMAP_SURFACE, vmin=0.0, vmax=1.0)
            fld.set_domain_limits(ax, case.mesh)
        fld.bare_axes(ax)
        ax.set_title(f"iteration {int(iterations[pick])}", loc="left", fontsize=9)
    for slot in range(len(picks), len(axes)):
        axes[slot].axis("off")

    fig.suptitle(
        f"{case.name}: density evolution", x=0.01, ha="left", fontsize=11,
        fontweight="bold", color=st.INK_PRIMARY,
    )
    st.annotate_note(
        fig,
        "Black is density 1 (full material), the page colour is density 0. "
        f"Snapshots taken from {len(table)} recorded frames.",
    )
    return st.save_figure(fig, path)


def animate_density_evolution(case: CaseResults, path: str, fps: int = 8,
                              hold_frames: int = 10, width: float = 5.6,
                              dpi: int = 110) -> str:
    """Animated GIF of the density field from the initial domain to the final
    topology.

    `hold_frames` repeats of the last frame keep the finished design on screen
    when the GIF loops.
    """
    if case.dim == 3:
        return p3.animate_density_evolution(case, path, fps=fps, hold_frames=hold_frames,
                                            width=width, dpi=dpi)
    import matplotlib.animation as animation

    table = case.density_history()
    if table is None:
        raise FileNotFoundError(f"{case.directory} has no density_history.csv")
    columns = [c for c in table.columns if c.startswith("rho_")]
    iterations = table["iteration"].to_numpy()
    shape = fld.structured_grid_shape(case.mesh)
    extent = case.mesh.extent
    setup = case.summary.get("optimization_setup", {})

    st.apply_style()
    import matplotlib.pyplot as plt

    # A GIF is committed to the repository, so it is rendered smaller and at a
    # lower DPI than the still figures: legibility on a README page, not print
    # resolution, is what it needs.
    aspect = (extent[3] - extent[2]) / max(extent[1] - extent[0], 1e-30)
    fig, ax = plt.subplots(
        figsize=(width, max(2.0, width * aspect * 0.82 + 1.05)), dpi=dpi
    )
    fig.patch.set_facecolor(st.SURFACE)

    first = table.iloc[0][columns].to_numpy(dtype=float)
    if shape is not None:
        image = ax.imshow(
            fld.density_image(first, shape), cmap=st.DENSITY_CMAP_SURFACE,
            vmin=0.0, vmax=1.0, origin="lower",
            extent=(extent[0], extent[1], extent[2], extent[3]),
            interpolation="nearest", aspect="equal",
        )
        collection = None
    else:
        image = None
        collection = fld.element_collection(ax, case.mesh, first,
                                            cmap=st.DENSITY_CMAP_SURFACE,
                                            vmin=0.0, vmax=1.0)
        fld.set_domain_limits(ax, case.mesh)
    fld.bare_axes(ax)

    header = ax.set_title("", loc="left", fontsize=9.5, color=st.INK_PRIMARY)
    caption = (
        f"{case.name}   vf {setup.get('volume_fraction_target', float('nan')):g}"
        f"   p {setup.get('simp_penalty', float('nan')):g}"
        f"   filter r {setup.get('filter_radius_m', float('nan')):.3g} m"
        "\nblack = material, page colour = void"
    )
    fig.text(0.012, 0.012, caption, fontsize=6.4, color=st.INK_MUTED,
             ha="left", va="bottom", linespacing=1.5)

    order = list(range(len(table))) + [len(table) - 1] * hold_frames
    history = case.history()

    def draw(frame_index: int):
        row = order[frame_index]
        density = table.iloc[row][columns].to_numpy(dtype=float)
        if image is not None:
            image.set_data(fld.density_image(density, shape))
        else:
            collection.set_array(density)
        iteration = int(iterations[row])
        compliance = float("nan")
        if history is not None:
            match = history.loc[history["iteration"] == iteration, "compliance[J]"]
            if not match.empty:
                compliance = float(match.iloc[0])
        header.set_text(
            f"iteration {iteration}   compliance "
            f"{st.format_si(compliance, digits=5)} J"
        )
        return [a for a in (image, collection, header) if a is not None]

    anim = animation.FuncAnimation(
        fig, draw, frames=len(order), interval=1000 / max(fps, 1), blit=False
    )
    directory = os.path.dirname(os.path.abspath(path))
    if directory:
        os.makedirs(directory, exist_ok=True)
    anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    st.shrink_gif_palette(path, fps)
    return path


# ---------------------------------------------------------------------------
# Stress-constrained vs unconstrained design
# ---------------------------------------------------------------------------
def plot_stress_constraint_comparison(constrained: CaseResults,
                                      unconstrained: CaseResults, path: str,
                                      load_case: Optional[str] = None) -> str:
    """Von Mises stress of the two designs on one colour scale with the limit.

    Both fields are the solid-material von Mises stress recovered from the SIMP
    model's displacement field, shown on the elements at or above the
    interpretation threshold; the re-solve peak of the thresholded structure,
    which is the number the summary judges the limit by, is quoted alongside.
    """
    if constrained.dim != 2 or unconstrained.dim != 2:
        raise ValueError("the stress comparison figure is drawn for plane cases")
    name = load_case or constrained.load_case_names[0]
    result = constrained.summary.get("optimization_result", {})
    interpreted = constrained.summary.get("interpreted_solid_analysis", {})
    stress_block = result.get("stress")
    if not stress_block:
        raise KeyError("the constrained run has no stress-constraint block in summary.json")
    limit = interpreted.get("stress_limit_Pa")
    if not limit:
        limit = stress_block["max_relaxed_stress_Pa"] / stress_block["max_relaxed_stress_ratio"]
    threshold = float(constrained.summary.get("solid_interpretation", {})
                      .get("threshold", 0.5))
    column = "solid_von_mises[Pa]"

    def field(case: CaseResults):
        table = case.density()
        if table is None:
            raise FileNotFoundError(f"{case.directory} has no density_final.csv")
        density = table["physical_density[-]"].to_numpy()
        stress = case.stress(name)
        values = stress[column if column in stress else "von_mises[Pa]"].to_numpy()
        return density >= threshold, values

    cases = [("unconstrained", unconstrained), ("stress-constrained", constrained)]
    fields = [field(c) for _, c in cases]
    vmax = max(float(values[mask].max()) for mask, values in fields)
    vmax = max(vmax, float(limit))

    # Two stacked panels: each carries a subtitle with its own numbers, which
    # side by side would not fit a half-width panel.
    xmin, xmax, ymin, ymax = constrained.mesh.extent
    aspect = (ymax - ymin) / max(xmax - xmin, 1.0e-12)
    fig, grid = st.figure(7.4, float(np.clip(2.0 * (5.2 * aspect + 1.3), 5.0, 9.0)),
                          nrows=2, ncols=1)
    axes = list(np.atleast_1d(grid).ravel())
    collection = None
    for ax, (label, case), (mask, values) in zip(axes, cases, fields):
        collection = fld.element_collection(ax, case.mesh, values, cmap=st.FIELD_CMAP,
                                            vmin=0.0, vmax=vmax, mask=mask)
        fld.mesh_outline(ax, case.mesh, color=st.INK_MUTED, linewidth=0.7, linestyle="--")
        fld.geometry_axes(ax, case.mesh)
        summary = case.summary
        peak_field = float(values[mask].max())
        resolve = summary.get("interpreted_solid_analysis", {}).get("max_von_mises_Pa")
        compliance = summary.get("optimization_result", {}).get("compliance_J", float("nan"))
        st.title(
            ax, label,
            f"compliance {st.format_si(compliance)} J; field peak "
            f"{st.format_si(peak_field)} Pa; re-solve peak "
            f"{st.format_si(resolve) if resolve else 'n/a'} Pa"
            + (f" = {resolve / limit:.3f} x limit" if resolve else ""),
        )
    bar = fld.add_colorbar(fig, collection, axes,
                           r"solid-material von Mises $\sigma_{vm}$ [Pa]")
    bar.ax.axhline(limit, color=st.series_color(7), linewidth=1.6)

    st.figure_title(
        fig, f"{constrained.name}: effect of the stress constraint, load case '{name}'",
        f"limit {st.format_si(limit)} Pa on the relaxed p-norm aggregate (the red "
        "line on the colour bar); constrained design: max relaxed ratio "
        f"{stress_block.get('max_relaxed_stress_ratio', float('nan')):.4f}, "
        f"re-solved structure at {interpreted.get('max_von_mises_over_limit', float('nan')):.3f} "
        "of the limit",
        ax=axes[0],
    )
    st.annotate_note(
        fig,
        f"Both panels show elements with density >= {threshold} only, coloured by the "
        "von Mises stress the material would carry at full density, recovered from "
        "the SIMP model's displacement field and averaged over each element's "
        "quadrature points. The re-solve peak in each title comes from analysing "
        "the thresholded structure as solid material, which is the check the "
        "summary reports against the limit. Neither is a manufacturability claim.",
    )
    return st.save_figure(fig, path)
