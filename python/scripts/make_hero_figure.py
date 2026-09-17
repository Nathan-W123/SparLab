#!/usr/bin/env python3
"""One standalone figure of an optimised part carrying its own stress field.

usage:
    python3 python/scripts/make_hero_figure.py \
        --case results/aerospace_bracket --load-case down_limit \
        --output docs/figures/hero_aerospace_bracket.png

Single panel, no comparison: the retained material of the optimised design,
coloured by the von Mises stress *in the material*, with the supports and the
applied load drawn on top. Everything is read from the run directory - the
geometry from `density_final.csv`, the stress from `stress_<case>.csv`, the
support nodes and load vectors from `mesh.json` - so the picture cannot drift
from the solver output.

Two deliberate choices worth knowing about:

* the field plotted is `solid_von_mises`, the von Mises of the *solid material*
  stress. The other column, `von_mises`, is the macroscopic SIMP stress, which
  in a partly dense element is an average over a cell that is part void and is
  not a material stress. Over the retained material the two differ by the SIMP
  factor;
* only elements at or above the interpretation threshold are drawn, and the
  original design domain is outlined so what was removed is still legible.
  This is an interpretation of a density field as geometry, so the threshold is
  printed on the figure.
"""

from __future__ import annotations

import argparse
import sys

import _bootstrap  # noqa: F401

import numpy as np
from matplotlib.patches import Circle, Rectangle

from sparlab_viz import fields as fld
from sparlab_viz import style as st
from sparlab_viz.loaders import ResultError
from sparlab_viz.loaders import load_case as load_results


def _clusters(coords, gap: float = 0.02):
    """Group support nodes by proximity and return (centre, radius) per group.

    The supports of these decks are circular bolt holes, so a single-link
    grouping on the node coordinates recovers one group per hole without the
    figure needing to know the deck's geometry.
    """
    remaining = [tuple(p) for p in coords]
    groups = []
    while remaining:
        seed = remaining.pop()
        group = [seed]
        changed = True
        while changed:
            changed = False
            for point in list(remaining):
                if any(abs(point[0] - q[0]) <= gap and abs(point[1] - q[1]) <= gap
                       for q in group):
                    group.append(point)
                    remaining.remove(point)
                    changed = True
        points = np.asarray(group)
        centre = points.mean(axis=0)
        radius = float(np.max(np.hypot(points[:, 0] - centre[0],
                                       points[:, 1] - centre[1])))
        groups.append((centre, radius))
    return groups


def _resultant(load_case: dict) -> tuple:
    """Total applied force and its application point, from the nodal forces."""
    forces = load_case.get("nodal_forces") or []
    if not forces:
        return None
    fx = float(sum(f["fx_N"] for f in forces))
    fy = float(sum(f["fy_N"] for f in forces))
    return fx, fy, [f["node"] for f in forces]


def build(case_dir: str, load_case: str, path: str, threshold: float = None,
          width: float = 10.0) -> str:
    case = load_results(case_dir)
    mesh = case.mesh
    summary = case.summary

    density_table = case.density()
    if density_table is None:
        raise ResultError(f"{case_dir} has no density_final.csv; this figure is "
                          "for a topology-optimisation run")
    density = density_table["physical_density[-]"].to_numpy()
    stress = case.stress(load_case)
    if "solid_von_mises[Pa]" not in stress.columns:
        raise ResultError("stress table has no solid_von_mises[Pa] column")

    if threshold is None:
        threshold = float(
            summary.get("solid_interpretation", {}).get("threshold", 0.5)
        )
    keep = density >= threshold
    values = stress["solid_von_mises[Pa]"].to_numpy() / 1.0e6  # MPa

    # Colour limits from the retained material only, clipped at the 99th
    # percentile: a handful of re-entrant corner cells sit far above the bulk of
    # the structure and would otherwise flatten the whole field to one colour.
    shown = values[keep]
    vmax = float(np.percentile(shown, 99.0))
    peak = float(shown.max())

    xmin, xmax, ymin, ymax = mesh.extent
    aspect = (xmax - xmin) / max(ymax - ymin, 1.0e-12)
    fig, ax = st.figure(width, width / aspect * 0.86 + 1.9)

    # The original design domain, so the removed material stays legible. Drawn
    # as one rectangle rather than through mesh_outline, whose per-element edge
    # segments are too short for a dash pattern to read.
    ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                           facecolor="none", edgecolor=st.INK_MUTED,
                           linewidth=1.1, linestyle=(0, (6, 4)), zorder=1.5,
                           label="original design domain"))
    collection = fld.element_collection(
        ax, mesh, values, cmap=st.FIELD_CMAP, vmin=0.0, vmax=vmax, mask=keep
    )
    collection.set_clim(0.0, vmax)

    # --- supports ---------------------------------------------------------
    # The constrained nodes fill each bolt hole, so a marker per node paints a
    # solid disc that competes with the stress ramp. Each cluster is drawn as
    # its own outline instead, which is how a support is normally annotated.
    fixed = mesh.constrained_nodes()
    if fixed.size:
        coords = mesh.nodes[fixed]
        for centre, radius in _clusters(coords):
            ax.add_patch(Circle(centre, radius, facecolor="none",
                                edgecolor=st.series_color(4), linewidth=2.2,
                                zorder=4))
            ax.plot([centre[0]], [centre[1]], "x", color=st.series_color(4),
                    markersize=7, markeredgewidth=1.8, zorder=4)
        ax.plot([], [], "o", markerfacecolor="none",
                markeredgecolor=st.series_color(4), markeredgewidth=2.2,
                markersize=9, linestyle="none",
                label=f"bolted: {fixed.size} nodes fixed in x and y")

    # --- applied load -----------------------------------------------------
    entry = next((lc for lc in (mesh.load_cases or [])
                  if lc.get("name") == load_case), None)
    if entry is not None:
        result = _resultant(entry)
        if result is not None:
            fx, fy, nodes = result
            magnitude = float(np.hypot(fx, fy))
            loaded = mesh.nodes[np.asarray(nodes, dtype=int)]
            tail = loaded.mean(axis=0)
            span = min(xmax - xmin, ymax - ymin)
            length = 0.30 * span
            scale = length / max(magnitude, 1.0e-30)
            ax.annotate(
                "", xy=(tail[0] + scale * fx, tail[1] + scale * fy),
                xytext=(tail[0], tail[1]), annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>", linewidth=2.6,
                                color=st.series_color(1), shrinkA=0, shrinkB=0,
                                mutation_scale=22),
            )
            ax.plot([], [], "-", color=st.series_color(1), linewidth=2.6,
                    label=f"applied load {st.format_si(magnitude)} N "
                          f"over {len(nodes)} nodes")

    fld.add_colorbar(fig, collection, ax,
                     "von Mises stress in the material [MPa]",
                     fraction=0.030, pad=0.015)
    fld.geometry_axes(ax, mesh)
    st.legend(ax, loc="lower left")

    result_block = summary.get("optimization_result", {})
    interp = summary.get("solid_interpretation", {})
    mesh_block = summary.get("mesh", {})
    names = summary.get("optimization_setup", {}).get("load_case_names") or []

    # Headline comparison, straight from the summary: the equal-mass uniform
    # plate is the only fair baseline in this 2-D idealisation, and the run
    # verifies the thickness scaling behind it rather than assuming it.
    comparison = summary.get("mass_stiffness_comparison", {})
    gain = comparison.get("stiffness_gain_over_equal_mass_plate")
    solid_modal = summary.get("modal_initial_solid") or {}
    topo_modal = summary.get("modal_optimised_topology") or {}
    headline = "stiffness and frequency compared in the summary"
    if gain:
        headline = f"{gain - 1.0:+.1%} stiffer than an equal-mass uniform plate"
        f1s = (solid_modal.get("frequencies_hz") or [None])[0]
        f1t = (topo_modal.get("frequencies_hz") or [None])[0]
        if f1s and f1t:
            headline += f" and {f1t / f1s - 1.0:+.0%} in first natural frequency"

    st.figure_title(
        fig,
        f"{case.name}: optimised topology carrying its own stress field",
        f"{mesh_block.get('num_elements', mesh.num_elements)} Q4 elements, "
        f"{mesh_block.get('num_dofs', 0)} DOFs, {len(names)} load cases, "
        f"{result_block.get('volume_fraction', float('nan')):.0%} volume "
        f"constraint - {headline}. Material shown is the "
        f"{interp.get('elements_retained', int(keep.sum()))} elements at density "
        f">= {threshold:g}, one connected group; the stress is the "
        f"'{load_case}' case",
    )
    st.annotate_note(
        fig,
        f"Colour is the von Mises stress of the solid material, clipped at the "
        f"99th percentile ({vmax:.0f} MPa) so the bulk of the field stays "
        f"readable; the peak is {peak:.0f} MPa at a re-entrant corner and is "
        f"mesh-sensitive. No stress constraint is modelled, so this is a "
        f"post-hoc field, not a substantiation. The dashed rectangle is the "
        f"design domain the optimiser started from.",
    )
    fig.savefig(path, dpi=200)
    import matplotlib.pyplot as plt
    plt.close(fig)
    return path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--case", default="results/aerospace_bracket")
    parser.add_argument("--load-case", default="down_limit")
    parser.add_argument("--output", default="docs/figures/hero_aerospace_bracket.png")
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--width", type=float, default=10.0)
    args = parser.parse_args(argv)

    try:
        print(f"  wrote {build(args.case, args.load_case, args.output, args.threshold, args.width)}")
    except (ResultError, FileNotFoundError, KeyError) as error:
        print(f"could not build the figure: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
