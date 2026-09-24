#!/usr/bin/env python3
"""Figures for the verification / validation studies, the runtime benchmarks
and the cross-validation against independent codes.

usage:
    python3 python/scripts/plot_verification.py \
        --verification results/verification \
        --benchmark results/benchmark \
        --cross-validation results/cross_validation \
        --figures docs/figures

Each figure is independent: a missing input is reported and skipped, and the
exit code is non-zero only when nothing at all could be drawn.
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
    parser.add_argument("--cross-validation", default="results/cross_validation")
    parser.add_argument("--figures", default="docs/figures")
    args = parser.parse_args(argv)

    os.makedirs(args.figures, exist_ok=True)
    written = []
    failures = []

    def figure(name: str) -> str:
        return os.path.join(args.figures, name)

    tasks = [
        ("sensitivity finite-difference comparison",
         lambda: studies.plot_sensitivity_check(
             args.verification, figure("verify_sensitivity.png"))),
        ("mesh convergence",
         lambda: studies.plot_mesh_convergence(
             args.verification, figure("verify_mesh_convergence.png"))),
        ("modal convergence",
         lambda: studies.plot_modal_convergence(
             args.verification, figure("verify_modal_convergence.png"))),
        ("sensitivity finite-difference comparison (3-D)",
         lambda: studies.plot_sensitivity_check_3d(
             args.verification, figure("verify_sensitivity_3d.png"))),
        ("mesh convergence (3-D)",
         lambda: studies.plot_mesh_convergence_3d(
             args.verification, figure("verify_mesh_convergence_3d.png"))),
        ("modal convergence (3-D)",
         lambda: studies.plot_modal_convergence_3d(
             args.verification, figure("verify_modal_convergence_3d.png"))),
        ("runtime scaling",
         lambda: studies.plot_runtime_scaling(
             args.benchmark, figure("runtime_scaling.png"))),
        ("runtime scaling (3-D)",
         lambda: studies.plot_runtime_scaling(
             args.benchmark, figure("runtime_scaling_3d.png"),
             stem="runtime_scaling_3d", dim=3)),
        ("cross-validation",
         lambda: studies.plot_cross_validation(
             args.cross_validation, figure("cross_validation.png"))),
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
