#!/usr/bin/env python3
"""Figures for the Heaviside projection.

usage:
    python3 python/scripts/plot_projection.py --results results \
                                              --figures docs/figures

Writes, from the runs that exist under --results:
    mbb_beam_projection.png   the MBB deck without and with the projection
                              (results/mbb_beam and results/mbb_beam_projected,
                              which scripts/run_all_benchmarks.sh writes)
    projection_summary.png    grey level and thresholded-vs-SIMP compliance of
                              every projected run beside its unprojected twin

Each figure is independent: a missing run is reported and skipped, and the
exit code is non-zero only when nothing could be drawn.
"""

from __future__ import annotations

import argparse
import os
import sys

import _bootstrap  # noqa: F401

from sparlab_viz import load_case, plots, studies
from sparlab_viz.loaders import ResultError


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results")
    parser.add_argument("--figures", default="docs/figures")
    parser.add_argument("--plain", default="mbb_beam",
                        help="run without the projection (a directory under --results)")
    parser.add_argument("--projected", default="mbb_beam_projected",
                        help="the same deck with the projection")
    args = parser.parse_args(argv)
    os.makedirs(args.figures, exist_ok=True)

    written, failures = [], []

    def comparison():
        plain = load_case(os.path.join(args.results, args.plain))
        projected = load_case(os.path.join(args.results, args.projected))
        return plots.plot_projection_comparison(
            plain, projected, os.path.join(args.figures, f"{args.plain}_projection.png"))

    tasks = [
        ("projection comparison", comparison),
        ("projection summary", lambda: studies.plot_projection_summary(
            args.results, os.path.join(args.figures, "projection_summary.png"))),
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
