#!/usr/bin/env python3
"""Generate every figure for one result directory.

usage:
    python3 python/scripts/plot_case.py --case results/cantilever_beam \
                                        --figures docs/figures

Produces, as applicable to what the run wrote:
    <case>_mesh_bcs.png             mesh, prescribed DOFs and applied loads
    <case>_deformed_<lc>.png        undeformed vs exaggerated deformed shape
    <case>_displacement_<lc>.png    displacement magnitude and u_x
    <case>_stress_<lc>.png          sigma_xx, sigma_yy, sigma_xy, von Mises
    <case>_reactions_<lc>.png       support reactions vs applied load
    <case>_modes.png                lowest mode shapes
    <case>_modes_topology.png       mode shapes of the interpreted topology
    <case>_modal_comparison.png     frequencies before and after optimisation
    <case>_topology.png             final density field and its interpretation
    <case>_convergence.png          compliance, volume, change and grey level
    <case>_density_evolution.png    density small multiples
    <case>_evolution.gif            density animation

With --minimal only <case>_topology.png and <case>_convergence.png are drawn.
"""

from __future__ import annotations

import argparse
import os
import sys

import _bootstrap  # noqa: F401  (sys.path side effect)

import numpy as np

from sparlab_viz import load_case, load_mesh, load_csv
from sparlab_viz import plots
from sparlab_viz.loaders import CaseResults, ResultError


def _figure_name(figures_dir: str, case_name: str, suffix: str,
                 extension: str = "png") -> str:
    return os.path.join(figures_dir, f"{case_name}_{suffix}.{extension}")


def _topology_case(case: CaseResults) -> CaseResults:
    """A CaseResults view whose mesh is the interpreted topology sub-mesh.

    Modal output tagged `topology` was computed on the extracted sub-mesh, so it
    must be plotted against that mesh rather than the design domain. The
    sub-mesh is not written separately; instead the mode-shape table carries its
    own nodal coordinates, which is what is used here.
    """
    table = case.mode_shapes("topology")
    if table is None:
        raise ResultError("no topology-tagged mode shapes in this run")
    from sparlab_viz.loaders import Mesh

    columns = [c for c in ("x[m]", "y[m]", "z[m]") if c in table]
    nodes = np.column_stack([table[c].to_numpy() for c in columns])
    # Rebuild connectivity by matching the retained elements of the parent mesh
    # against the sub-mesh node coordinates.
    parent = case.mesh
    if nodes.shape[1] != parent.dim:
        raise ResultError("the topology mode-shape table and the mesh differ in dimension")
    lookup = {tuple(np.round(row, 12)): i for i, row in enumerate(nodes)}
    elements = []
    for element in parent.elements:
        mapped = []
        for node in element:
            key = tuple(np.round(parent.nodes[node], 12))
            index = lookup.get(key)
            if index is None:
                mapped = []
                break
            mapped.append(index)
        if len(mapped) == parent.elements.shape[1]:
            elements.append(mapped)
    if not elements:
        raise ResultError(
            "could not reconstruct the interpreted topology sub-mesh from the "
            "mode-shape coordinates"
        )
    sub_mesh = Mesh(nodes=nodes, elements=np.asarray(elements, dtype=int),
                    element_type=parent.element_type)
    return CaseResults(directory=case.directory, mesh=sub_mesh,
                       summary=case.summary)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", required=True, help="result directory of one run")
    parser.add_argument("--figures", default="docs/figures",
                        help="directory to write figures into")
    parser.add_argument("--name", default=None,
                        help="figure-name prefix (default: the case name)")
    parser.add_argument("--no-animation", action="store_true",
                        help="skip the GIF, which is the slowest artefact")
    parser.add_argument("--load-cases", default=None,
                        help="comma-separated subset of load cases to plot")
    parser.add_argument("--minimal", action="store_true",
                        help="only the final topology and the convergence history "
                             "(for a variant of a deck whose other figures already "
                             "exist, such as the same deck with the projection on)")
    args = parser.parse_args(argv)

    case = load_case(args.case)
    prefix = args.name or os.path.basename(os.path.normpath(args.case))
    os.makedirs(args.figures, exist_ok=True)
    written = []

    if args.minimal:
        if not case.is_topology_run:
            print(f"  --minimal draws topology figures; {args.case} is not a topology run",
                  file=sys.stderr)
            return 1
        written.append(plots.plot_final_topology(
            case, _figure_name(args.figures, prefix, "topology")))
        written.append(plots.plot_convergence_history(
            case, _figure_name(args.figures, prefix, "convergence")))
        for path in written:
            print(f"  wrote {path}")
        return 0

    density = None
    if case.is_topology_run:
        table = case.density()
        if table is not None:
            density = table["physical_density[-]"].to_numpy()

    written.append(plots.plot_mesh_and_bcs(
        case, _figure_name(args.figures, prefix, "mesh_bcs")))

    names = case.load_case_names
    if args.load_cases:
        wanted = {n.strip() for n in args.load_cases.split(",")}
        names = [n for n in names if n in wanted]

    for name in names:
        if not case.has(f"displacement_{name}.csv"):
            continue
        written.append(plots.plot_deformed_shape(
            case, name, _figure_name(args.figures, prefix, f"deformed_{name}"),
            density=density))
        written.append(plots.plot_displacement_magnitude(
            case, name, _figure_name(args.figures, prefix, f"displacement_{name}"),
            density=density))
        if case.has(f"stress_{name}.csv"):
            written.append(plots.plot_stress_fields(
                case, name, _figure_name(args.figures, prefix, f"stress_{name}"),
                density=density))
        if case.has(f"reactions_{name}.csv"):
            written.append(plots.plot_reactions(
                case, name, _figure_name(args.figures, prefix, f"reactions_{name}")))

    # --- modal ---------------------------------------------------------------
    for tag, suffix in (("", "modes"), ("solid", "modes")):
        if case.modes(tag) is not None and case.mode_shapes(tag) is not None:
            try:
                written.append(plots.plot_mode_shapes(
                    case, _figure_name(args.figures, prefix, suffix), tag=tag))
            except (ValueError, ResultError) as error:
                print(f"  skipped {suffix}: {error}", file=sys.stderr)
            break

    if case.modes("topology") is not None:
        try:
            written.append(plots.plot_mode_shapes(
                _topology_case(case),
                _figure_name(args.figures, prefix, "modes_topology"),
                tag="topology"))
        except (ValueError, ResultError) as error:
            print(f"  skipped topology mode shapes: {error}", file=sys.stderr)

    if "modal_optimised_topology" in case.summary:
        written.append(plots.plot_modal_comparison(
            case.summary, _figure_name(args.figures, prefix, "modal_comparison"),
            case.name, dim=case.dim))

    # --- topology ------------------------------------------------------------
    if case.is_topology_run:
        written.append(plots.plot_final_topology(
            case, _figure_name(args.figures, prefix, "topology")))
        written.append(plots.plot_convergence_history(
            case, _figure_name(args.figures, prefix, "convergence")))
        if case.density_history() is not None:
            written.append(plots.plot_density_evolution(
                case, _figure_name(args.figures, prefix, "density_evolution")))
            if not args.no_animation:
                written.append(plots.animate_density_evolution(
                    case, _figure_name(args.figures, prefix, "evolution", "gif")))

    for path in written:
        print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
