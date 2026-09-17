#!/usr/bin/env python3
"""Figures for the verification / validation studies and the runtime benchmark.

usage:
    python3 python/scripts/plot_verification.py \
        --verification results/verification \
        --benchmark results/benchmark \
        --figures docs/figures
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
    parser.add_argument("--verification", default="results/verification")
    parser.add_argument("--benchmark", default="results/benchmark")
    parser.add_argument("--figures", default="docs/figures")
    args = parser.parse_args(argv)

    os.makedirs(args.figures, exist_ok=True)
    written = []
    failures = []

    tasks = [
        ("sensitivity finite-difference comparison",
         lambda: studies.plot_sensitivity_check(
             args.verification, os.path.join(args.figures, "verify_sensitivity.png"))),
        ("mesh convergence",
         lambda: studies.plot_mesh_convergence(
             args.verification,
             os.path.join(args.figures, "verify_mesh_convergence.png"))),
        ("modal convergence",
         lambda: studies.plot_modal_convergence(
             args.verification,
             os.path.join(args.figures, "verify_modal_convergence.png"))),
        ("runtime scaling",
         lambda: studies.plot_runtime_scaling(
             args.benchmark, os.path.join(args.figures, "runtime_scaling.png"))),
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
    return 1 if failures and not written else 0


if __name__ == "__main__":
    raise SystemExit(main())
