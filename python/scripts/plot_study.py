#!/usr/bin/env python3
"""Aggregate and plot the aerospace parametric design study.

usage:
    python3 python/scripts/plot_study.py --study results/study \
                                         --figures docs/figures \
                                         --table docs/results/aerospace_study.csv

Writes one CSV with every study point and the figures the design-study
documentation uses.
"""

from __future__ import annotations

import argparse
import os
import sys

import _bootstrap  # noqa: F401

from sparlab_viz import studies
from sparlab_viz.loaders import ResultError


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--study", default="results/study")
    parser.add_argument("--figures", default="docs/figures")
    parser.add_argument("--table", default="docs/results/aerospace_study.csv")
    args = parser.parse_args(argv)

    table = studies.collect_study(args.study)
    os.makedirs(args.figures, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.table)), exist_ok=True)
    table.sort_values(["arm", "tag"]).to_csv(args.table, index=False)
    print(f"  wrote {args.table} ({len(table)} study points)")

    written = []
    failures = []
    tasks = [
        ("mass-stiffness Pareto",
         lambda: studies.plot_mass_stiffness_pareto(
             table, os.path.join(args.figures, "study_pareto.png"))),
        ("numerical settings",
         lambda: studies.plot_numerical_settings(
             table, os.path.join(args.figures, "study_settings.png"))),
        ("mesh dependence",
         lambda: studies.plot_mesh_dependence(
             table, os.path.join(args.figures, "study_mesh_dependence.png"))),
        ("load weighting",
         lambda: studies.plot_load_weighting(
             table, os.path.join(args.figures, "study_load_weighting.png"))),
        ("volume-fraction topologies",
         lambda: studies.plot_study_topologies(
             args.study, table, "volume_fraction",
             os.path.join(args.figures, "study_topologies_volume.png"),
             "Aerospace bracket: optimised topology vs volume fraction",
             "Black is density 1, the page colour is density 0. Every panel uses "
             "the same domain, loads, supports and passive regions; only the "
             "volume-fraction target differs.",
             order_by="volume_fraction_target")),
        ("penalty topologies",
         lambda: studies.plot_study_topologies(
             args.study, table, "penalty",
             os.path.join(args.figures, "study_topologies_penalty.png"),
             "Aerospace bracket: optimised topology vs SIMP penalty p",
             "A low penalty leaves large grey regions because intermediate "
             "density is not penalised enough to be driven out; a high penalty "
             "sharpens the design but makes the problem more non-convex and more "
             "dependent on the starting point.",
             order_by="penalty")),
        ("filter-radius topologies",
         lambda: studies.plot_study_topologies(
             args.study, table, "filter_radius",
             os.path.join(args.figures, "study_topologies_filter.png"),
             "Aerospace bracket: optimised topology vs filter radius",
             "The filter radius sets the minimum member size. A small radius "
             "gives many thin members and a lower compliance; a large radius "
             "gives fewer, thicker members that a manufacturing process is more "
             "likely to be able to produce, at a stiffness cost.",
             order_by="filter_radius_m")),
        ("load-weighting topologies",
         lambda: studies.plot_study_topologies(
             args.study, table, "load_weighting",
             os.path.join(args.figures, "study_topologies_loads.png"),
             "Aerospace bracket: optimised topology vs lateral load weight",
             "Tags are (down, reversal, lateral) weights. As the lateral case "
             "gains weight the design adds chordwise material to carry it, at "
             "the cost of vertical stiffness.",
             order_by="weight_lateral")),
    ]

    for label, task in tasks:
        try:
            written.append(task())
        except (ResultError, FileNotFoundError, KeyError, ValueError) as error:
            failures.append(f"{label}: {error}")

    for path in written:
        print(f"  wrote {path}")
    for failure in failures:
        print(f"  skipped {failure}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
