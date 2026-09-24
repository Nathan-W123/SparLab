#!/usr/bin/env python3
"""The stress-constrained design beside the unconstrained run of the same deck.

usage:
    python3 python/scripts/plot_comparison.py \
        --constrained results/l_bracket_stress \
        --unconstrained results/l_bracket_unconstrained \
        --figures docs/figures

Writes <name>_stress_comparison.png, where <name> defaults to the constrained
case's directory name. The unconstrained run is the same configuration run
with `sparlab_topopt --no-stress`, which scripts/run_all_benchmarks.sh does.
"""

from __future__ import annotations

import argparse
import os

import _bootstrap  # noqa: F401

from sparlab_viz import load_case
from sparlab_viz import plots


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--constrained", required=True,
                        help="result directory of the stress-constrained run")
    parser.add_argument("--unconstrained", required=True,
                        help="result directory of the same deck run without the constraint")
    parser.add_argument("--figures", default="docs/figures")
    parser.add_argument("--name", default=None, help="figure-name prefix")
    parser.add_argument("--load-case", default=None,
                        help="load case to show (default: the first)")
    args = parser.parse_args(argv)

    constrained = load_case(args.constrained)
    unconstrained = load_case(args.unconstrained)
    prefix = args.name or os.path.basename(os.path.normpath(args.constrained))
    os.makedirs(args.figures, exist_ok=True)
    path = plots.plot_stress_constraint_comparison(
        constrained, unconstrained,
        os.path.join(args.figures, f"{prefix}_stress_comparison.png"),
        load_case=args.load_case,
    )
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
