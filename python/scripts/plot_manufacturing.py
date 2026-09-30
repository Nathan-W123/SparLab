#!/usr/bin/env python3
"""Comparison figures for the buckling-constrained column and shell panel,
the robust MBB beam, the overhang-filtered designs, the bridge under its own
weight and the heated clamped beam.

usage:
    python3 python/scripts/plot_manufacturing.py --results results \
                                                 --figures docs/figures

Each figure sets a benchmark beside the comparison runs that
scripts/run_all_benchmarks.sh makes of the same deck with one feature
switched off or a load scaled. A missing run is reported and its figure skipped; the exit
code is non-zero if any figure could not be drawn.
"""

from __future__ import annotations

import argparse
import os
import sys

import _bootstrap  # noqa: F401

from sparlab_viz import buckling, manufacturing, physics
from sparlab_viz.loaders import ResultError, load_case


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results")
    parser.add_argument("--figures", default="docs/figures")
    args = parser.parse_args(argv)
    os.makedirs(args.figures, exist_ok=True)

    def case(name: str):
        return load_case(os.path.join(args.results, name))

    def figure(name: str) -> str:
        return os.path.join(args.figures, name)

    tasks = [
        ("buckling comparison", lambda: buckling.plot_buckling_comparison(
            [("compliance only", case("column_buckling_unconstrained")),
             ("lambda >= 6, plain projection", case("column_buckling_nonrobust")),
             ("lambda >= 6, robust projection", case("column_buckling"))],
            figure("column_buckling_comparison.png"), name="column_buckling")),
        ("buckling-constrained convergence", lambda: buckling.plot_buckling_history(
            case("column_buckling"), figure("column_buckling_history.png"))),
        ("shell buckling comparison", lambda: buckling.plot_buckling_comparison(
            [("compliance only", case("shell_panel_buckling_unconstrained")),
             ("lambda >= 10, plain projection", case("shell_panel_buckling_nonrobust")),
             ("lambda >= 10, robust projection", case("shell_panel_buckling"))],
            figure("shell_panel_buckling_comparison.png"), name="shell_panel_buckling")),
        ("self-weight comparison", lambda: physics.plot_self_weight_comparison(
            [(0.0, case("bridge_self_weight_g0")), (1.0, case("bridge_self_weight")),
             (5.0, case("bridge_self_weight_g5"))],
            figure("bridge_self_weight_comparison.png"), name="bridge_self_weight")),
        ("thermal sweep", lambda: physics.plot_thermal_sweep(
            [(float(dt), case(f"clamped_beam_thermal_dT{dt}")) for dt in (0, 1, 2, 3, 5, 20, 40)]
            + [(10.0, case("clamped_beam_thermal"))],
            figure("clamped_beam_thermal_sweep.png"), name="clamped_beam_thermal")),
        ("robust formulation", lambda: manufacturing.plot_robust_comparison(
            case("mbb_beam_robust_off"), case("mbb_beam_robust"),
            figure("mbb_beam_robust_comparison.png"))),
        ("overhang filter (2-D)", lambda: manufacturing.plot_overhang_comparison(
            [("no filter", case("mbb_beam_overhang_off")),
             ("filter, built +y", case("mbb_beam_overhang")),
             ("filter, built -y", case("mbb_beam_overhang_down"))],
            figure("mbb_beam_overhang_comparison.png"), name="mbb_beam_overhang")),
        ("overhang filter (3-D)", lambda: manufacturing.plot_overhang_comparison_3d(
            [("no filter", case("bracket_3d_overhang_off")),
             ("overhang filter", case("bracket_3d_overhang"))],
            figure("bracket_3d_overhang_comparison.png"), name="bracket_3d_overhang")),
    ]
    failures = []
    for label, task in tasks:
        try:
            print(f"  wrote {task()}")
        except (ResultError, FileNotFoundError, KeyError, ValueError) as error:
            failures.append(f"{label}: {error}")
    for failure in failures:
        print(f"  skipped {failure}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
