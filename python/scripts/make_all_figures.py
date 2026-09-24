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
    "aerospace_bracket",
    "wing_rib",
    "l_bracket_stress",
    "bracket_3d",
]

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
        if args.no_animation:
            argv_case.append("--no-animation")
        if not run(f"figures for {case}", argv_case):
            failures.append(f"{case}: plot_case.py failed")

    if not run("verification, benchmark and cross-validation figures",
               [os.path.join(SCRIPTS, "plot_verification.py"),
                "--verification", os.path.join(args.results, "verification"),
                "--benchmark", os.path.join(args.results, "benchmark"),
                "--cross-validation", os.path.join(args.results, "cross_validation"),
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
