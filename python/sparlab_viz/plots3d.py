"""Figures for solid (Hex8) results.

`plots` dispatches here when a case is three-dimensional. Each figure mirrors
its 2-D counterpart: the same title, the same numbers from the same summary
fields, and the same rule that an interpretation - a deformation scale, a
density threshold read as geometry - is stated on the figure itself.

A solid is shown through its boundary surface. For a topology run the surface
drawn is that of the elements at or above the interpretation threshold, which
is exactly the `structure_after.stl` the solver exported; the design domain is
kept as a dashed outline so the reader sees what was removed.
"""

from __future__ import annotations

import os
from functools import reduce
from typing import Optional

import numpy as np

from . import fields as fld
from . import solid as sd
from . import style as st
from .loaders import CaseResults, ResultError, load_json

DOMAIN_TINT = "#dcdad2"
PASSIVE_SOLID = "#4a3aa7"
PASSIVE_VOID = "#e34948"
DEFAULT_THRESHOLD = 0.5


def _axes3d(width: float, height: float, nrows: int = 1, ncols: int = 1):
    return st.figure(width, height, nrows=nrows, ncols=ncols,
                     subplot_kw={"projection": "3d"})


def _threshold(case: CaseResults) -> float:
    return float(case.summary.get("solid_interpretation", {})
                 .get("threshold", DEFAULT_THRESHOLD))


def _structure(mesh, density, threshold):
    """Boundary faces (and owning elements) of the material above `threshold`.

    With `density` None the whole mesh is the structure.
    """
    mask = None if density is None else np.asarray(density, dtype=float) >= threshold
    faces, owners = mesh.boundary_faces(mask, return_owners=True)
    return faces, owners, mask


def _domain_outline(ax, mesh, coords=None, max_segments: Optional[int] = None, **kwargs):
    """Dashed feature edges of the mesh as the undeformed reference.

    A box has twelve; an interpreted topology has hundreds, which drawn
    through a coloured surface is clutter, so above `max_segments` the
    bounding box is drawn instead (the caller's caption says which).
    """
    faces = mesh.boundary_faces()
    points = mesh.nodes if coords is None else coords
    segments = sd.feature_edges(faces, points)
    if max_segments is not None and len(segments) > max_segments:
        lo, hi = points.min(axis=0), points.max(axis=0)
        corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1])
                            for z in (lo[2], hi[2])])
        pairs = [(0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (1, 3), (4, 6), (5, 7),
                 (0, 4), (1, 5), (2, 6), (3, 7)]
        segments = corners[np.asarray(pairs)]
    return sd.edges(ax, segments, **kwargs)


def _view(mesh) -> dict:
    """Camera for a case: one fixed view, z up, for every figure of a case.

    A camera chosen per case (say, the loaded axis up) was tried and read
    worse: for the cantilevered bracket it put the tip face nearest the viewer
    and hid the tapered top and the cavities that make the topology legible.
    A single fixed view also keeps the before/after, evolution and field
    figures of one case directly comparable. `mesh` is accepted so the choice
    can be revisited in one place.
    """
    del mesh
    return dict(sd.VIEW)


def _horizontal_colorbar(fig, mappable, ax, label: str, shrink: float = 0.72):
    """A colour bar under a 3-D panel: a tall vertical bar beside a squat
    solid dwarfs it and collides with the z-axis labels."""
    from matplotlib.ticker import FuncFormatter

    bar = fig.colorbar(mappable, ax=ax, orientation="horizontal", fraction=0.045,
                       pad=0.02, shrink=shrink)
    # Stresses in MPa on the ticks: a "1e7" offset label under a short bar
    # collides with the unit label.
    vmin, vmax = mappable.get_clim()
    if "[Pa]" in label and max(abs(vmin), abs(vmax)) >= 1.0e5:
        bar.formatter = FuncFormatter(lambda value, _pos: f"{value / 1.0e6:g}")
        bar.update_ticks()
        label = label.replace("[Pa]", "[MPa]")
    bar.set_label(label, fontsize=8.5, color=st.INK_SECONDARY)
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, colors=st.INK_SECONDARY, length=2)
    return bar


def _density_note(mask, threshold: float) -> str:
    if mask is None:
        return ""
    return (f" Only the {int(mask.sum())} elements with density >= {threshold} are "
            "drawn: the field on a near-void SIMP cell is a macroscopic average, "
            "not a material response, and the surface shown is the exported "
            "structure_after geometry.")


# ---------------------------------------------------------------------------
# Mesh, boundary conditions and loads
# ---------------------------------------------------------------------------
def _element_value_note(mesh) -> str:
    """How an element's stress value was obtained, for the figure text."""
    if mesh.element_type == "Tet4":
        return "element values are constant over each linear tetrahedron"
    return "element values are the average over the 2x2x2 stiffness quadrature points"


def plot_mesh_and_bcs(case: CaseResults, path: str) -> str:
    mesh = case.mesh
    # A tall part (the engine mount's upright) draws a taller box, whose axis
    # labels the layout engine does not see; the figure grows with the
    # height-to-length ratio so they stay clear of the footnote.
    extent = mesh.extent
    tallness = (extent[5] - extent[4]) / max(extent[1] - extent[0], 1.0e-30)
    fig, ax = _axes3d(7.8, 4.8 + 2.0 * float(np.clip(tallness - 0.25, 0.0, 0.75)))
    # Draw in the order added rather than by mplot3d's per-artist depth sort,
    # so markers inside holes and passive regions stay visible on top.
    ax.computed_zorder = False
    faces = mesh.boundary_faces()
    sd.surface(ax, mesh.nodes[faces], color=DOMAIN_TINT,
               edge_color=(0.0, 0.0, 0.0, 0.07), linewidth=0.12)
    _domain_outline(ax, mesh, color=st.INK_SECONDARY, linewidth=1.0)

    passive_note = ""
    table = case.density() if case.is_topology_run else None
    if table is not None and "passive_tag[-]" in table:
        tags = table["passive_tag[-]"].to_numpy()
        for value, color, label in (
            (1.0, PASSIVE_SOLID, "passive solid (density pinned at 1)"),
            (2.0, PASSIVE_VOID, "passive void (density pinned at 0)"),
        ):
            selected = tags == value
            if not selected.any():
                continue
            sd.surface(ax, mesh.nodes[mesh.boundary_faces(selected)], color=color,
                       edge_color="none")
            ax.plot([], [], [], "s", color=color, label=label)
        counts = [int((tags == 1.0).sum()), int((tags == 2.0).sum())]
        if any(counts):
            passive_note = (
                f" {counts[0]} passive solid and {counts[1]} passive void "
                "elements are pinned and not optimised."
            )

    fixed = [mesh.constrained_nodes(c) for c in "xyz"]
    clamped = reduce(np.intersect1d, fixed)
    partial = np.setdiff1d(np.unique(np.concatenate(fixed)) if any(f.size for f in fixed)
                           else np.empty(0, dtype=int), clamped)
    for nodes, marker, label, slot in (
        (clamped, "s", "fixed in x, y and z", 0),
        (partial, "^", "fixed in one or two components", 2),
    ):
        if nodes.size == 0:
            continue
        ax.scatter(mesh.nodes[nodes, 0], mesh.nodes[nodes, 1], mesh.nodes[nodes, 2],
                   marker=marker, s=9, color=st.series_color(slot), depthshade=False,
                   label=f"{label} ({nodes.size} nodes)", zorder=5)

    reference = max(extent[1] - extent[0], extent[3] - extent[2], extent[5] - extent[4])
    load_slots = [1, 4, 5, 7]
    for index, name in enumerate(mesh.load_case_names):
        forces = mesh.nodal_forces(name)
        if forces.size == 0:
            continue
        nodes = forces[:, 0].astype(int)
        resultant = forces[:, 1:4].sum(axis=0)
        magnitude = float(np.linalg.norm(resultant))
        color = st.series_color(load_slots[index % len(load_slots)])
        ax.scatter(mesh.nodes[nodes, 0], mesh.nodes[nodes, 1], mesh.nodes[nodes, 2],
                   marker="o", s=4 if nodes.size > 100 else 8, color=color,
                   depthshade=False, zorder=6,
                   label=f"load case '{name}': {st.format_si(magnitude)} N over "
                         f"{nodes.size} node(s)")
        if magnitude <= 0.0:
            continue
        origin = mesh.nodes[nodes].mean(axis=0)
        direction = resultant / magnitude * 0.3 * reference
        # The arrow starts on the loaded nodes and points along the resultant;
        # it is drawn last so it sits on top of the markers.
        ax.quiver(origin[0], origin[1], origin[2], direction[0], direction[1],
                  direction[2], color=color, linewidth=2.2, arrow_length_ratio=0.25,
                  zorder=8)

    sd.set_solid_axes(ax, mesh.nodes, margin=0.10, view=_view(mesh))
    st.title(
        ax, f"{case.name}: mesh and boundary conditions",
        f"{mesh.num_elements} {mesh.element_type} elements, {mesh.num_nodes} nodes, "
        f"{len(mesh.prescribed)} prescribed DOFs",
    )
    # The legend sits over the top of the box, so it gets an opaque backing
    # in the page colour rather than letting the edges run through its text.
    st.legend(ax, loc="upper right", fontsize=8, frameon=True, facecolor=st.SURFACE,
              edgecolor="none", framealpha=0.92)
    if ax.get_legend() is not None:
        ax.get_legend().set_zorder(20)  # above the outline (drawn in order added)
    st.annotate_note(
        fig,
        "The tinted body is the design domain, drawn through its boundary faces "
        "with flat shading from a fixed light; the dark lines are its edges. "
        "Each load case is shown as markers at the loaded nodes plus one arrow "
        "for its resultant; arrow length is a fixed fraction of the domain, not "
        "proportional to the load. Markers and passive regions are drawn on top "
        "of the surface, so nodes on hidden faces show through." + passive_note,
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Static fields on the surface
# ---------------------------------------------------------------------------
def plot_deformed_shape(case: CaseResults, load_case: str, path: str,
                        density: Optional[np.ndarray] = None,
                        threshold: Optional[float] = None) -> str:
    mesh = case.mesh
    threshold = _threshold(case) if threshold is None else threshold
    disp = case.displacement(load_case)
    u = disp[["ux[m]", "uy[m]", "uz[m]"]].to_numpy(dtype=float)
    magnitude = disp["umag[m]"].to_numpy(dtype=float)
    scale = sd.auto_scale(mesh.nodes, magnitude)
    coords = mesh.nodes + scale * u
    faces, _owners, mask = _structure(mesh, density, threshold)

    fig, ax = _axes3d(7.8, 5.0)
    _domain_outline(ax, mesh, color=st.INK_MUTED, linewidth=0.9, linestyle="--")
    _collection, mappable = sd.surface(
        ax, coords[faces], sd.face_values(faces, magnitude), cmap=st.FIELD_CMAP,
        vmin=0.0, vmax=float(magnitude.max()) or 1.0,
    )
    if mappable is not None:
        _horizontal_colorbar(fig, mappable, ax, "|u| [m]")
    sd.set_solid_axes(ax, mesh.nodes, coords[np.unique(faces)] if faces.size else None,
                      margin=0.08, view=_view(mesh))
    st.title(
        ax, f"{case.name}: deformed shape, load case '{load_case}'",
        f"displacements exaggerated {scale:.4g}x; peak |u| = "
        f"{st.format_si(float(magnitude.max()))} m; dashed lines are the "
        "undeformed domain edges",
    )
    st.annotate_note(
        fig,
        f"The surface is the boundary of the structure displaced by {scale:.4g} "
        "times the computed displacement and coloured by |u| averaged over each "
        "face's corners. The exaggeration is for visibility only; the analysis "
        "itself is small-strain linear elasticity." + _density_note(mask, threshold),
    )
    return st.save_figure(fig, path)


def plot_displacement_magnitude(case: CaseResults, load_case: str, path: str,
                                density: Optional[np.ndarray] = None,
                                threshold: Optional[float] = None) -> str:
    mesh = case.mesh
    threshold = _threshold(case) if threshold is None else threshold
    disp = case.displacement(load_case)
    magnitude = disp["umag[m]"].to_numpy(dtype=float)
    faces, _owners, mask = _structure(mesh, density, threshold)

    fig, ax = _axes3d(7.8, 4.8)
    _domain_outline(ax, mesh, color=st.INK_MUTED, linewidth=0.8, linestyle="--")
    _collection, mappable = sd.surface(
        ax, mesh.nodes[faces], sd.face_values(faces, magnitude), cmap=st.FIELD_CMAP,
        vmin=0.0, vmax=float(magnitude.max()) or 1.0,
    )
    if mappable is not None:
        _horizontal_colorbar(fig, mappable, ax, "|u| [m]")
    sd.set_solid_axes(ax, mesh.nodes, margin=0.08, view=_view(mesh))
    st.title(
        ax, f"{case.name}: displacement, load case '{load_case}'",
        f"peak |u| = {st.format_si(float(magnitude.max()))} m on the undeformed "
        "configuration; sequential scale",
    )
    st.annotate_note(
        fig,
        "Nodal displacement magnitude averaged over the corners of each boundary "
        "face; interior nodes are not visible." + _density_note(mask, threshold),
    )
    return st.save_figure(fig, path)


def plot_stress_fields(case: CaseResults, load_case: str, path: str,
                       density: Optional[np.ndarray] = None,
                       threshold: Optional[float] = None) -> str:
    mesh = case.mesh
    threshold = _threshold(case) if threshold is None else threshold
    stress = case.stress(load_case)
    faces, owners, mask = _structure(mesh, density, threshold)

    panels = [
        ("von_mises[Pa]", r"von Mises $\sigma_{vm}$ [Pa]", st.FIELD_CMAP, False),
        ("principal_max[Pa]", r"max principal $\sigma_1$ [Pa]", st.SIGNED_CMAP, True),
        ("principal_min[Pa]", r"min principal $\sigma_3$ [Pa]", st.SIGNED_CMAP, True),
        ("sxx[Pa]", r"$\sigma_{xx}$ [Pa]", st.SIGNED_CMAP, True),
    ]
    fig, grid = _axes3d(8.0, 6.0, nrows=2, ncols=2)
    axes = list(np.atleast_1d(grid).ravel())
    for ax, (column, label, cmap, signed) in zip(axes, panels):
        values = stress[column].to_numpy(dtype=float)
        shown = values[mask] if mask is not None else values
        on_faces = values[owners]
        _domain_outline(ax, mesh, color=st.INK_MUTED, linewidth=0.7, linestyle="--")
        if signed:
            _c, mappable = sd.surface(ax, mesh.nodes[faces], on_faces, cmap=cmap,
                                      norm=st.symmetric_norm(shown))
        else:
            _c, mappable = sd.surface(ax, mesh.nodes[faces], on_faces, cmap=cmap,
                                      vmin=0.0,
                                      vmax=float(np.max(shown)) if shown.size else 1.0)
        if mappable is not None:
            _horizontal_colorbar(fig, mappable, ax, label, shrink=0.8)
        sd.set_solid_axes(ax, mesh.nodes, margin=0.06, view=_view(mesh))

    peak = (float(stress["von_mises[Pa]"].to_numpy()[mask].max()) if mask is not None
            else float(stress["von_mises[Pa]"].max()))
    st.figure_title(
        fig, f"{case.name}: stress on the surface, load case '{load_case}'",
        f"peak von Mises {st.format_si(peak)} Pa; {_element_value_note(mesh)}, "
        "shown through the element that owns each boundary face. Panels: von Mises, "
        "max principal, min principal, sigma_xx - each labelled on its own colour bar",
        ax=axes[0],
    )
    st.annotate_note(
        fig,
        "Principal and normal components use a diverging scale with zero exactly "
        "at the neutral midpoint and symmetric limits; von Mises is non-negative "
        "and uses a sequential scale. Interior elements are not visible."
        + _density_note(mask, threshold),
    )
    return st.save_figure(fig, path)


def plot_reactions(case: CaseResults, load_case: str, path: str) -> str:
    """Applied load against the support reactions, component by component.

    A 3-D arrow plot of hundreds of nodal reactions is unreadable, so the
    quantitative check - the force and moment balance the solver reported - is
    the figure, with the nodal distribution summarised in the subtitle.
    """
    reactions = case.reactions(load_case)
    entry = next((lc for lc in case.summary.get("load_cases", [])
                  if lc["name"] == load_case), None)
    if entry is None:
        raise ResultError(f"summary.json has no load case '{load_case}'")
    eq = entry["equilibrium"]
    applied = np.asarray(eq["applied_force_N"], dtype=float)
    reaction = np.asarray(eq["reaction_force_N"], dtype=float)
    components = ["x", "y", "z"][: applied.size]
    index = np.arange(len(components))
    width = 0.38

    fig, ax = st.figure(7.0, 3.6)
    ax.axhline(0.0, color=st.INK_MUTED, linewidth=0.9)
    ax.bar(index - width / 2, applied, width, color=st.series_color(0),
           label="applied load")
    ax.bar(index + width / 2, reaction, width, color=st.series_color(1),
           label="sum of support reactions")
    for x, value in zip(index - width / 2, applied):
        ax.text(x, value, f" {st.format_si(value)}", ha="center",
                va="bottom" if value >= 0 else "top", fontsize=7.4,
                color=st.INK_SECONDARY)
    for x, value in zip(index + width / 2, reaction):
        ax.text(x, value, f" {st.format_si(value)}", ha="center",
                va="bottom" if value >= 0 else "top", fontsize=7.4,
                color=st.INK_SECONDARY)
    ax.set_xticks(index)
    ax.set_xticklabels([f"F{c}" for c in components])
    ax.set_ylabel("force [N]")
    peak = float(reactions["rmag[N]"].max()) if "rmag[N]" in reactions else float("nan")
    moment_text = ""
    if "relative_moment_error" in eq:
        moment_text = f", moment balance error {eq['relative_moment_error']:.2e}"
    st.title(
        ax, f"{case.name}: equilibrium, load case '{load_case}'",
        f"{len(reactions)} constrained nodes, largest nodal reaction "
        f"{st.format_si(peak)} N; relative force-balance error "
        f"{eq['relative_force_error']:.2e}{moment_text}",
    )
    st.legend(ax, loc="best")
    st.annotate_note(
        fig,
        "Each pair of bars should be equal and opposite. The relative errors are "
        "|F_applied + F_reaction| / |F_applied| (and the same for moments about "
        "the origin), computed by the solver after the solve.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Modal
# ---------------------------------------------------------------------------
def plot_mode_shapes(case: CaseResults, path: str, tag: str = "",
                     num_modes: int = 6) -> str:
    modes = case.modes(tag)
    shapes = case.mode_shapes(tag)
    if modes is None or shapes is None:
        raise FileNotFoundError(
            f"{case.directory} has no modal output" + (f" for tag '{tag}'" if tag else "")
        )
    mesh = case.mesh
    if len(shapes) != mesh.num_nodes:
        raise ValueError(
            f"mode_shapes has {len(shapes)} rows but the mesh has {mesh.num_nodes} "
            "nodes; the two belong to different meshes"
        )
    faces = mesh.boundary_faces()
    count = min(num_modes, len(modes))
    rows = (count + 1) // 2
    # Row height from the domain's shape: a slender bar needs less than a cube.
    extent = mesh.extent
    spans = np.array([extent[1] - extent[0], extent[3] - extent[2], extent[5] - extent[4]])
    slenderness = float(spans.max() / max(np.sort(spans)[1], 1e-30))
    row_height = float(np.clip(2.6 / max(slenderness, 1.0) ** 0.5, 1.4, 2.6))
    fig, grid = _axes3d(7.8, row_height * rows + 0.5, nrows=rows, ncols=2)
    axes = np.atleast_1d(grid).ravel()

    for index in range(count):
        ax = axes[index]
        phi = np.column_stack([shapes[f"u{c}_mode{index}[m]"].to_numpy(dtype=float)
                               for c in "xyz"])
        magnitude = np.linalg.norm(phi, axis=1)
        scale = sd.auto_scale(mesh.nodes, magnitude, 0.10)
        coords = mesh.nodes + scale * phi
        _domain_outline(ax, mesh, max_segments=60, color=st.INK_MUTED, linewidth=0.6,
                        linestyle="--")
        sd.surface(ax, coords[faces], sd.face_values(faces, magnitude),
                   cmap=st.FIELD_CMAP, vmin=0.0, edge_color="none")
        sd.bare_solid_axes(ax, mesh.nodes, coords, margin=0.06, view=_view(mesh))
        row = modes.iloc[index]
        ax.set_title(f"mode {index + 1}: {row['frequency[Hz]']:.1f} Hz", loc="left",
                     fontsize=9.5)
    for index in range(count, len(axes)):
        axes[index].axis("off")

    label = f" ({tag})" if tag else ""
    fig.suptitle(f"{case.name}: mode shapes{label}", x=0.01, ha="left", fontsize=11,
                 fontweight="bold", color=st.INK_PRIMARY)
    st.annotate_note(
        fig,
        "Mode shapes are M-orthonormal (phi^T M phi = 1) and are drawn as the "
        "boundary surface on an exaggerated displaced configuration; colour is "
        "the local mode-shape magnitude. Dashed lines are the undeformed edges, "
        "or the undeformed bounding box when the structure has too many edges to "
        "show through the surface.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Topology optimisation
# ---------------------------------------------------------------------------
def initial_density(case: CaseResults) -> float:
    """The uniform value the free design variables started from.

    Newer summaries record it; for older ones the echoed deck says whether
    `topology.initial_density` overrode the default start at the volume
    fraction.
    """
    setup = case.summary.get("optimization_setup", {})
    if setup.get("initial_density") is not None:
        return float(setup["initial_density"])
    deck = os.path.join(case.directory, "config.json")
    if os.path.isfile(deck):
        try:
            value = (load_json(deck).get("topology") or {}).get("initial_density")
            if value is not None and float(value) >= 0.0:
                return float(value)
        except (ResultError, ValueError, TypeError):
            pass
    return float(setup.get("volume_fraction_target", float("nan")))


def density_field_label(case: CaseResults) -> str:
    """What the plotted density is: filtered, and projected when the run used
    the Heaviside projection (with its final sharpness)."""
    projection = case.summary.get("optimization_result", {}).get("projection")
    if projection:
        return (f"Physical density: filtered, then Heaviside-projected (beta = "
                f"{projection.get('final_beta', 0):g}, eta = {projection.get('eta', 0.5):g})")
    return "Physical (filtered) density"


def plot_final_topology(case: CaseResults, path: str) -> str:
    """The design domain before and the interpreted structure after.

    The 'after' surface is the boundary of the elements at or above the
    threshold - the same geometry as the exported structure_after.stl - with
    each face coloured by its element's density, so the reader can see where
    the interpretation cut through intermediate material.
    """
    table = case.density()
    if table is None:
        raise FileNotFoundError(f"{case.directory} has no density_final.csv")
    density = table["physical_density[-]"].to_numpy(dtype=float)
    threshold = _threshold(case)
    result = case.summary.get("optimization_result", {})
    setup = case.summary.get("optimization_setup", {})
    interpretation = case.summary.get("solid_interpretation", {})
    mesh = case.mesh

    fig, (before, after) = _axes3d(8.0, 4.2, nrows=1, ncols=2)
    domain_faces = mesh.boundary_faces()
    sd.surface(before, mesh.nodes[domain_faces], color=DOMAIN_TINT,
               edge_color=(0.0, 0.0, 0.0, 0.07), linewidth=0.12)
    _domain_outline(before, mesh, color=st.INK_SECONDARY, linewidth=0.9)
    sd.set_solid_axes(before, mesh.nodes, margin=0.06, view=_view(mesh))
    start = initial_density(case)
    target = float(setup.get("volume_fraction_target", float("nan")))
    before.set_title(
        f"before: design domain, uniform density {start:g}"
        if abs(start - target) < 1.0e-12
        else f"before: uniform start {start:g}, target {target:g}",
        loc="left", fontsize=9,
    )

    mask = density >= threshold
    faces, owners = mesh.boundary_faces(mask, return_owners=True)
    _domain_outline(after, mesh, color=st.INK_MUTED, linewidth=0.8, linestyle="--")
    _collection, mappable = sd.surface(after, mesh.nodes[faces], density[owners],
                                       cmap=st.DENSITY_CMAP_SURFACE, vmin=0.0, vmax=1.0,
                                       edge_color=(0.0, 0.0, 0.0, 0.10), linewidth=0.15)
    if mappable is not None:
        _horizontal_colorbar(fig, mappable, after, "density [-]", shrink=0.8)
    sd.set_solid_axes(after, mesh.nodes, margin=0.06, view=_view(mesh))
    after.set_title(
        f"after: {int(mask.sum())} of {mesh.num_elements} elements at density >= "
        f"{threshold}", loc="left", fontsize=9,
    )

    st.figure_title(
        fig, f"{case.name}: structure before and after optimisation",
        f"{density_field_label(case)} at the last iterate. Compliance "
        f"{st.format_si(result.get('compliance_J', float('nan')))} J, volume "
        f"fraction {result.get('volume_fraction', float('nan')):.4f}, grey level "
        f"{result.get('grey_level', float('nan')):.3f}, "
        f"{result.get('iterations', 0)} iterations; "
        f"{interpretation.get('elements_retained', int(mask.sum()))} elements "
        f"retained in {interpretation.get('connected_groups_above_threshold', 1)} "
        "connected group(s), "
        f"{st.format_si(interpretation.get('volume_discarded_as_islands_m3', 0.0))} "
        "m^3 discarded as islands",
        ax=before,
    )
    st.annotate_note(
        fig,
        "The density field is a relative material distribution, not a solid "
        "body. The right-hand panel is one explicit interpretation of it - the "
        "boundary of the elements at or above the threshold, which is what "
        "structure_after.stl contains - and its faces are coloured by the "
        "density of the element they belong to, so a light face is material the "
        "threshold only just kept. Threshold and connectivity rule are stated "
        "because the answer depends on both.",
    )
    return st.save_figure(fig, path)


def plot_density_evolution(case: CaseResults, path: str, frames: int = 6) -> str:
    table = case.density_history()
    if table is None:
        raise FileNotFoundError(f"{case.directory} has no density_history.csv")
    columns = [c for c in table.columns if c.startswith("rho_")]
    iterations = table["iteration"].to_numpy()
    picks = np.unique(np.linspace(0, len(table) - 1, frames).round().astype(int))
    threshold = _threshold(case)
    mesh = case.mesh

    cols = 3 if len(picks) > 4 else 2
    rows = (len(picks) + cols - 1) // cols
    fig, grid = _axes3d(7.8, 2.3 * rows + 0.4, nrows=rows, ncols=cols)
    axes = np.atleast_1d(grid).ravel()
    empty = 0
    for slot, pick in enumerate(picks):
        ax = axes[slot]
        density = table.iloc[pick][columns].to_numpy(dtype=float)
        mask = density >= threshold
        _domain_outline(ax, mesh, color=st.INK_MUTED, linewidth=0.6, linestyle="--")
        if mask.any():
            faces, owners = mesh.boundary_faces(mask, return_owners=True)
            sd.surface(ax, mesh.nodes[faces], density[owners],
                       cmap=st.DENSITY_CMAP_SURFACE, vmin=0.0, vmax=1.0,
                       edge_color="none")
        else:
            empty += 1
        sd.bare_solid_axes(ax, mesh.nodes, margin=0.04, view=_view(mesh))
        ax.set_title(f"iteration {int(iterations[pick])}", loc="left", fontsize=9)
    for slot in range(len(picks), len(axes)):
        axes[slot].axis("off")

    fig.suptitle(f"{case.name}: density evolution", x=0.01, ha="left", fontsize=11,
                 fontweight="bold", color=st.INK_PRIMARY)
    st.annotate_note(
        fig,
        f"Each frame is the boundary of the elements at density >= {threshold}, "
        "coloured by density (black is 1). The dashed box is the design domain."
        + (f" {empty} frame(s) show no material because every density was still "
           "below the threshold - the design starts uniform." if empty else "")
        + f" Snapshots taken from {len(table)} recorded frames.",
    )
    return st.save_figure(fig, path)


def animate_density_evolution(case: CaseResults, path: str, fps: int = 8,
                              hold_frames: int = 10, width: float = 5.6,
                              dpi: int = 110) -> str:
    import matplotlib.animation as animation

    table = case.density_history()
    if table is None:
        raise FileNotFoundError(f"{case.directory} has no density_history.csv")
    columns = [c for c in table.columns if c.startswith("rho_")]
    iterations = table["iteration"].to_numpy()
    threshold = _threshold(case)
    mesh = case.mesh
    setup = case.summary.get("optimization_setup", {})

    st.apply_style()
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(width, width * 0.78), dpi=dpi)
    fig.patch.set_facecolor(st.SURFACE)
    ax = fig.add_subplot(projection="3d")
    _domain_outline(ax, mesh, color=st.INK_MUTED, linewidth=0.7, linestyle="--")
    sd.bare_solid_axes(ax, mesh.nodes, margin=0.03, view=_view(mesh))
    header = fig.text(0.02, 0.95, "", fontsize=9.5, color=st.INK_PRIMARY, ha="left",
                      va="top")
    caption = (
        f"{case.name}   vf {setup.get('volume_fraction_target', float('nan')):g}"
        f"   p {setup.get('simp_penalty', float('nan')):g}"
        f"   filter r {setup.get('filter_radius_m', float('nan')):.3g} m"
        f"\nsurface of the elements at density >= {threshold}, black = density 1"
    )
    fig.text(0.02, 0.02, caption, fontsize=6.4, color=st.INK_MUTED, ha="left",
             va="bottom", linespacing=1.5)

    order = list(range(len(table))) + [len(table) - 1] * hold_frames
    history = case.history()
    state = {"collection": None}

    def draw(frame_index: int):
        row = order[frame_index]
        density = table.iloc[row][columns].to_numpy(dtype=float)
        if state["collection"] is not None:
            state["collection"].remove()
            state["collection"] = None
        mask = density >= threshold
        if mask.any():
            faces, owners = mesh.boundary_faces(mask, return_owners=True)
            state["collection"], _ = sd.surface(
                ax, mesh.nodes[faces], density[owners], cmap=st.DENSITY_CMAP_SURFACE,
                vmin=0.0, vmax=1.0, edge_color="none",
            )
        iteration = int(iterations[row])
        compliance = float("nan")
        if history is not None:
            match = history.loc[history["iteration"] == iteration, "compliance[J]"]
            if not match.empty:
                compliance = float(match.iloc[0])
        header.set_text(f"iteration {iteration}   compliance "
                        f"{st.format_si(compliance, digits=5)} J"
                        + ("" if mask.any() else "   (no element above the threshold yet)"))
        return [a for a in (state["collection"], header) if a is not None]

    anim = animation.FuncAnimation(fig, draw, frames=len(order),
                                   interval=1000 / max(fps, 1), blit=False)
    directory = os.path.dirname(os.path.abspath(path))
    if directory:
        os.makedirs(directory, exist_ok=True)
    anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    st.shrink_gif_palette(path, fps)
    return path
