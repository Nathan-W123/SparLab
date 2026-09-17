#!/usr/bin/env python3
"""Turn the run summaries into the small Markdown/CSV tables the docs quote.

usage:
    python3 python/scripts/write_result_tables.py --results results \
                                                  --output docs/results

Every number written here is read straight out of a `summary.json`, so the
documentation and the solver output cannot drift apart. Re-run this after any
benchmark run and the tables update.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional

import _bootstrap  # noqa: F401

import pandas as pd

from sparlab_viz.loaders import ResultError, load_json


BENCHMARK_CASES = ["cantilever_beam", "mbb_beam", "aerospace_bracket", "wing_rib"]


def _fmt(value, digits: int = 4) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int,)):
        return str(value)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number != number:  # NaN
        return "n/a"
    if number == 0:
        return "0"
    if 1e-3 <= abs(number) < 1e5:
        return f"{number:.{digits}g}"
    return f"{number:.{digits}g}"


def _markdown_table(headers: List[str], rows: List[List[str]]) -> str:
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))
    lines = [
        "| " + " | ".join(h.ljust(widths[i]) for i, h in enumerate(headers)) + " |",
        "|" + "|".join("-" * (w + 2) for w in widths) + "|",
    ]
    for row in rows:
        lines.append(
            "| " + " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)) + " |"
        )
    return "\n".join(lines)


def benchmark_table(results_dir: str) -> Optional[str]:
    rows = []
    frames = []
    for case in BENCHMARK_CASES:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        mesh = doc.get("mesh", {})
        setup = doc.get("optimization_setup", {})
        result = doc.get("optimization_result", {})
        comparison = doc.get("mass_stiffness_comparison", {})
        interp = doc.get("solid_interpretation", {})
        record = {
            "case": case,
            "elements": mesh.get("num_elements"),
            "dofs": mesh.get("num_dofs"),
            "volume_fraction_target": setup.get("volume_fraction_target"),
            "volume_fraction_achieved": result.get("volume_fraction"),
            "volume_violation": result.get("volume_constraint_relative_violation"),
            "iterations": result.get("iterations"),
            "stop_reason": result.get("stop_reason"),
            "compliance_J": result.get("compliance_J"),
            "initial_compliance_J": result.get("initial_compliance_J"),
            "full_solid_compliance_J": comparison.get("full_solid_compliance_J"),
            "equal_mass_plate_compliance_J":
                comparison.get("equal_mass_uniform_plate_compliance_J"),
            "stiffness_gain_over_equal_mass_plate":
                comparison.get("stiffness_gain_over_equal_mass_plate"),
            "grey_level": result.get("grey_level"),
            "retained_elements": interp.get("elements_retained"),
            "connected_groups": interp.get("connected_groups_above_threshold"),
            "islands_discarded_m3": interp.get("volume_discarded_as_islands_m3"),
            "runtime_s": result.get("total_seconds"),
            "s_per_iteration": result.get("seconds_per_iteration"),
        }
        solid = doc.get("modal_initial_solid")
        topology = doc.get("modal_optimised_topology")
        if solid and solid.get("frequencies_hz"):
            record["f1_solid_Hz"] = solid["frequencies_hz"][0]
            record["mass_solid_kg"] = solid.get("total_mass_kg")
        if topology and topology.get("frequencies_hz"):
            record["f1_topology_Hz"] = topology["frequencies_hz"][0]
            record["mass_topology_kg"] = topology.get("total_mass_kg")
        frames.append(record)
        rows.append([
            case,
            _fmt(record["elements"]),
            _fmt(record["volume_fraction_target"]),
            _fmt(record["volume_fraction_achieved"], 6),
            _fmt(record["volume_violation"], 2),
            _fmt(record["iterations"]),
            str(record["stop_reason"]),
            _fmt(record["compliance_J"], 6),
            _fmt(record["equal_mass_plate_compliance_J"], 6),
            _fmt(record["stiffness_gain_over_equal_mass_plate"], 4),
            _fmt(record["grey_level"], 3),
            _fmt(record["runtime_s"], 4),
        ])

    if not rows:
        return None
    pd.DataFrame(frames).to_csv(
        os.path.join(_OUTPUT, "benchmarks.csv"), index=False
    )
    headers = [
        "case", "elements", "vf target", "vf achieved", "vf violation",
        "iterations", "stop reason", "compliance [J]",
        "equal-mass plate [J]", "stiffness gain", "grey", "runtime [s]",
    ]
    return _markdown_table(headers, rows)


#: Arm name -> (what it varies, the summary field that records it).
STUDY_ARMS = {
    "volume_fraction": ("volume-fraction target", "volume_fraction_target"),
    "load_weighting": ("lateral load-case weight", "weight_lateral"),
    "mesh_fixed_r": ("elements, filter radius fixed in metres", "num_elements"),
    "mesh_fixed_cells": ("elements, filter radius fixed in cells", "num_elements"),
    "penalty": ("SIMP penalty p", "penalty"),
    "filter_radius": ("filter radius [m]", "filter_radius_m"),
    "youngs_modulus": ("Young's modulus [Pa]", "youngs_modulus_Pa"),
}


def study_table(study_dir: str) -> Optional[str]:
    """One row per design-study arm, straight from the study summaries."""
    from sparlab_viz import studies

    try:
        table = studies.collect_study(study_dir)
    except (ResultError, FileNotFoundError):
        return None

    rows = []
    for arm, (what, column) in STUDY_ARMS.items():
        sub = table[table["arm"] == arm]
        if sub.empty:
            continue
        span = "n/a"
        if column in sub.columns and sub[column].notna().any():
            span = f"{_fmt(sub[column].min())} - {_fmt(sub[column].max())}"
        rows.append([
            arm,
            what,
            str(len(sub)),
            span,
            f"{int(sub['converged'].sum())}/{len(sub)}",
            f"{_fmt(sub['compliance_J'].min(), 6)} - {_fmt(sub['compliance_J'].max(), 6)}",
            f"{_fmt(sub['grey_level'].min(), 3)} - {_fmt(sub['grey_level'].max(), 3)}",
            str(int(sub["connected_groups"].max())),
            _fmt(sub["seconds"].sum(), 4),
        ])
    if not rows:
        return None
    headers = ["arm", "varies", "points", "range", "converged", "compliance [J]",
               "grey", "worst group count", "total runtime [s]"]
    return _markdown_table(headers, rows)


def verification_table(results_dir: str) -> Optional[str]:
    path = os.path.join(results_dir, "verification", "summary.json")
    if not os.path.isfile(path):
        return None
    doc = load_json(path)
    rows = []
    frames = []
    for outcome in doc.get("outcomes", []):
        rows.append([
            outcome["study"],
            outcome["kind"],
            outcome["metric"],
            _fmt(outcome["value"], 4),
            _fmt(outcome["tolerance"], 3),
            "PASS" if outcome["passed"] else "FAIL",
        ])
        frames.append(outcome)
    if not rows:
        return None
    pd.DataFrame(frames).to_csv(
        os.path.join(_OUTPUT, "verification.csv"), index=False
    )
    return _markdown_table(
        ["study", "kind", "metric", "value", "tolerance", "result"], rows
    )


def modal_table(results_dir: str) -> Optional[str]:
    rows = []
    for case in BENCHMARK_CASES + ["cantilever_analysis"]:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        for key, label in (
            ("modal", "static model"),
            ("modal_initial_solid", "full solid domain"),
            ("modal_mass_matched_baseline", "equal-mass uniform plate"),
            ("modal_optimised_topology", "optimised topology"),
        ):
            block = doc.get(key)
            if not block or not block.get("frequencies_hz"):
                continue
            freqs = block["frequencies_hz"][:4]
            rows.append([
                case, label,
                _fmt(block.get("total_mass_kg"), 5),
                *[_fmt(f, 5) for f in freqs],
                *["n/a"] * (4 - len(freqs)),
                _fmt(max(block.get("eigenpair_residuals", [0.0])), 2),
                "yes" if block.get("converged") else "no",
            ])
    if not rows:
        return None
    return _markdown_table(
        ["case", "structure", "mass [kg]", "f1 [Hz]", "f2 [Hz]", "f3 [Hz]",
         "f4 [Hz]", "max residual", "converged"],
        rows,
    )


def scaling_table(results_dir: str) -> Optional[str]:
    path = os.path.join(results_dir, "benchmark", "runtime_scaling.json")
    if not os.path.isfile(path):
        return None
    doc = load_json(path)
    rows = []
    for record in doc.get("records", []):
        rows.append([
            f"{record['nx']} x {record['ny']}",
            _fmt(record["num_elements"]),
            _fmt(record["num_dofs"]),
            _fmt(record["stiffness_nonzeros"]),
            _fmt(record["assemble_s"], 3),
            _fmt(record["factorize_s"], 3),
            _fmt(record["solve_s"], 3),
            _fmt(record["objective_gradient_s"], 3),
        ])
    if not rows:
        return None
    exponents = doc.get("scaling_exponent_vs_dofs", {})
    table = _markdown_table(
        ["mesh", "elements", "DOFs", "nonzeros", "assemble [s]",
         "factorise [s]", "solve [s]", "obj+grad [s]"],
        rows,
    )
    slope_line = ", ".join(
        f"{k} {v:.2f}" for k, v in exponents.items()
    )
    return table + f"\n\nFitted slopes of log(time) vs log(DOFs): {slope_line}.\n"


_OUTPUT = "docs/results"


def main(argv=None) -> int:
    global _OUTPUT
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results")
    parser.add_argument("--output", default="docs/results")
    args = parser.parse_args(argv)
    _OUTPUT = args.output
    os.makedirs(args.output, exist_ok=True)

    sections = [
        ("Verification and validation", verification_table(args.results),
         "Generated by `python/scripts/write_result_tables.py` from "
         "`results/verification/summary.json`. Verification compares against "
         "exact answers for the discrete problem; validation compares against an "
         "independent theory, where a finite gap is expected."),
        ("Benchmark results", benchmark_table(args.results),
         "Generated from each `results/<case>/summary.json`. `stiffness gain` is "
         "the compliance of an equal-mass uniform plate divided by the optimised "
         "compliance, so values above 1 mean the optimised design is stiffer at "
         "the same mass."),
        ("Modal results", modal_table(args.results),
         "`max residual` is the largest relative eigenpair residual "
         "||K phi - lambda M phi|| / ||lambda M phi|| over the reported modes."),
        ("Runtime scaling", scaling_table(args.results),
         "Minimum over repeats at each size, which is the standard estimator for "
         "wall-clock benchmarks."),
        ("Design study", study_table(os.path.join(args.results, "study")),
         "One row per arm of the aerospace parametric study; the per-point table "
         "is `aerospace_study.csv` beside this file and the interpretation is in "
         "`docs/aerospace_study.md`. A worst group count above 1 means that arm "
         "contains a point whose thresholded density field is not a single "
         "connected structure."),
    ]

    lines = [
        "# SparLab result tables",
        "",
        "**Generated file - do not edit by hand.** Refresh with `make results` "
        "after a benchmark run.",
        "",
    ]
    missing = []
    for title, table, note in sections:
        lines.append(f"## {title}")
        lines.append("")
        if table is None:
            missing.append(title)
            lines.append("_Not available: the corresponding run has not been made._")
        else:
            lines.append(note)
            lines.append("")
            lines.append(table)
        lines.append("")

    path = os.path.join(args.output, "README.md")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    print(f"  wrote {path}")
    for name in missing:
        print(f"  missing section: {name}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
