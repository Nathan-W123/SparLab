#!/usr/bin/env python3
"""Regenerate every figure and the animation from the results on disk.

usage:
    python3 python/scripts/make_all_figures.py --results results \
                                               --figures docs/figures

Each step is independent: a missing run produces a skip message and a non-zero
exit code at the end, never a half-written figure. This is the script behind
`make figures`.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

import _bootstrap  # noqa: F401

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
CASES = [
    "cantilever_analysis",
    "block_3d_analysis",
    "cantilever_beam",
    "mbb_beam",
    "mbb_beam_projected",
    "aerospace_bracket",
    "wing_rib",
    "l_bracket_stress",
    "bracket_3d",
    "bracket_3d_projected",
    "lug_bracket_2d",
    "engine_mount_3d",
    "bracket_3d_large",
    "column_buckling",
    "mbb_beam_robust",
    "mbb_beam_overhang",
    "bracket_3d_overhang",
]

#: Variants of another benchmark (the same problem with the projection, a
#: buckling constraint, the robust formulation or the overhang filter on):
#: only their topology and convergence figures are drawn, the rest would repeat
#: the parent; plot_manufacturing.py draws their comparisons.
MINIMAL_CASES = {"mbb_beam_projected", "bracket_3d_projected", "column_buckling",
                 "mbb_beam_robust", "mbb_beam_overhang", "bracket_3d_overhang"}

#: Cases whose density animation is skipped: 110k hexahedra per frame make the
#: GIF slow to render and large, and the small-multiples figure shows the same.
NO_ANIMATION_CASES = {"bracket_3d_large"}

#: The stress-constrained case and the unconstrained run of the same deck
#: (sparlab_topopt --no-stress) that the comparison figure sets beside it.
STRESS_COMPARISON = ("l_bracket_stress", "l_bracket_unconstrained")


def run(label: str, argv: list) -> bool:
    print(f"\n--- {label}")
    completed = subprocess.run([sys.executable] + argv, cwd=os.getcwd())
    return completed.returncode == 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results")
    parser.add_argument("--figures", default="docs/figures")
    parser.add_argument("--no-animation", action="store_true")
    parser.add_argument("--only", default=None,
                        help="comma-separated subset of cases")
    args = parser.parse_args(argv)

    cases = CASES
    if args.only:
        wanted = {c.strip() for c in args.only.split(",")}
        cases = [c for c in CASES if c in wanted]

    failures = []
    for case in cases:
        directory = os.path.join(args.results, case)
        if not os.path.isdir(directory):
            failures.append(f"{case}: {directory} does not exist")
            continue
        argv_case = [os.path.join(SCRIPTS, "plot_case.py"),
                     "--case", directory, "--figures", args.figures]
        if case in MINIMAL_CASES:
            argv_case.append("--minimal")
        elif args.no_animation or case in NO_ANIMATION_CASES:
            argv_case.append("--no-animation")
        if not run(f"figures for {case}", argv_case):
            failures.append(f"{case}: plot_case.py failed")

    if not run("verification, benchmark and cross-validation figures",
               [os.path.join(SCRIPTS, "plot_verification.py"),
                "--verification", os.path.join(args.results, "verification"),
                "--benchmark", os.path.join(args.results, "benchmark"),
                "--cross-validation", os.path.join(args.results, "cross_validation"),
                "--tet10-study", os.path.join(args.results, "tet10_part_study"),
                "--figures", args.figures]):
        failures.append("verification/benchmark figures failed")

    constrained, unconstrained = (os.path.join(args.results, c) for c in STRESS_COMPARISON)
    if not args.only or STRESS_COMPARISON[0] in cases:
        if os.path.isdir(constrained) and os.path.isdir(unconstrained):
            if not run("stress-constraint comparison figure",
                       [os.path.join(SCRIPTS, "plot_comparison.py"),
                        "--constrained", constrained, "--unconstrained", unconstrained,
                        "--figures", args.figures]):
                failures.append("stress-constraint comparison figure failed")
        else:
            failures.append(f"stress comparison: {constrained} and {unconstrained} "
                            "must both exist (scripts/run_all_benchmarks.sh writes them)")

    if not args.only or "mbb_beam_projected" in cases:
        if not run("Heaviside projection figures",
                   [os.path.join(SCRIPTS, "plot_projection.py"),
                    "--results", args.results, "--figures", args.figures]):
            failures.append("projection figures failed")

    if not args.only or any(c in cases for c in ("column_buckling", "mbb_beam_robust",
                                                 "mbb_beam_overhang",
                                                 "bracket_3d_overhang")):
        if not run("buckling, robust and overhang comparison figures",
                   [os.path.join(SCRIPTS, "plot_manufacturing.py"),
                    "--results", args.results, "--figures", args.figures]):
            failures.append("buckling / robust / overhang figures failed")

    study_dir = os.path.join(args.results, "study")
    if os.path.isdir(study_dir):
        if not run("design-study figures",
                   [os.path.join(SCRIPTS, "plot_study.py"),
                    "--study", study_dir, "--figures", args.figures]):
            failures.append("design-study figures failed")
    else:
        failures.append(f"study: {study_dir} does not exist")

    print("\n--- summary")
    if failures:
        for failure in failures:
            print(f"  incomplete: {failure}", file=sys.stderr)
        print(f"  {len(failures)} step(s) did not produce output", file=sys.stderr)
        return 1
    print("  every figure regenerated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
