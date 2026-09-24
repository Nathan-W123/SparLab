#!/usr/bin/env python3
"""Cross-validate SparLab's static solution against two independent FE codes.

usage:
    python3 python/scripts/cross_validate.py \
        --case results/cantilever_analysis --case results/block_3d_analysis \
        --output results/cross_validation

For every result directory given (a `sparlab_solve` run made with
`--export-calculix`), the same discrete problem - identical mesh, element type,
material, consistent nodal loads and prescribed DOFs, all read from the run's
own `mesh.json` and CalculiX decks - is solved by

  * CalculiX (`ccx`), from the exported `calculix_<load case>.inp`, and
  * scikit-fem, assembled in this script from `mesh.json` with the linear
    isotropic elasticity form on ElementQuad1 / ElementHex1,

and the nodal displacements are compared with SparLab's `displacement_<load
case>.csv` node by node. Reported per code and load case:

  * `max_abs_diff_m`, `max_rel_diff` = max |u_ref - u| / max |u|, and the RMS
    of the same ratio,
  * `ref_max_abs_m` and `sparlab_max_abs_m`,
  * pass/fail against the tolerance stated in the output.

Two things are deliberately measured rather than assumed. scikit-fem uses the
same element formulation, so it agrees to the level the two linear solvers
agree (1e-10 and better); the tolerance is 1e-7. CalculiX writes its nodal
results to the .frd file with six significant digits, so even for C3D8 - the
same element as SparLab's Hex8 - the comparison is bounded below by that
rounding: half a unit in the sixth digit of the largest displacement, 5e-6
relative. The CalculiX tolerance is therefore 1e-5 for both element families,
and the summary records the rounding floor next to the measured difference.
CalculiX's CPS4 plane-stress element is *expanded* into a solid through the
thickness with the plane-stress condition imposed on the expanded element,
which is a different discretisation of the plane problem; how close it lands
is measured, not assumed.

Nothing here recomputes SparLab's numbers: they are read from the run.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional

import _bootstrap  # noqa: F401

import numpy as np

from sparlab_viz.loaders import Mesh, ResultError, load_case


# ---------------------------------------------------------------------------
# CalculiX
# ---------------------------------------------------------------------------
def parse_frd_displacements(path: str, num_nodes: int) -> np.ndarray:
    """Nodal displacements (num_nodes, 3) from a CalculiX .frd file (last step)."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    blocks: List[np.ndarray] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(" -4") and "DISP" in line:
            values = np.full((num_nodes, 3), np.nan)
            i += 1
            # component definition lines start with " -5"
            while i < len(lines) and lines[i].startswith(" -5"):
                i += 1
            while i < len(lines) and lines[i].startswith(" -1"):
                row = lines[i]
                node = int(row[3:13])
                comps = [float(row[13 + 12 * k: 25 + 12 * k]) for k in range(3)]
                values[node - 1] = comps
                i += 1
            blocks.append(values)
            continue
        i += 1
    if not blocks:
        raise ResultError(f"{path} holds no displacement block")
    return blocks[-1]


def run_calculix(deck: str, workdir: str) -> str:
    """Run ccx on `deck` inside `workdir` and return the .frd path."""
    job = os.path.splitext(os.path.basename(deck))[0]
    shutil.copy(deck, os.path.join(workdir, job + ".inp"))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(["ccx", "-i", job], cwd=workdir, capture_output=True,
                          text=True, env=env, check=False)
    frd = os.path.join(workdir, job + ".frd")
    if proc.returncode != 0 or not os.path.isfile(frd):
        raise ResultError(f"ccx failed on {deck}:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    return frd


def calculix_version() -> str:
    """The ccx version number ("2.21"), or its raw banner if it cannot be parsed."""
    proc = subprocess.run(["ccx", "-v"], capture_output=True, text=True, check=False)
    text = (proc.stdout + proc.stderr).strip()
    match = re.search(r"[Vv]ersion\s+([0-9][0-9.]*)", text)
    if match:
        return match.group(1)
    lines = text.splitlines()
    return lines[-1].strip() if lines else "unknown"


# ---------------------------------------------------------------------------
# scikit-fem
# ---------------------------------------------------------------------------
def solve_with_skfem(mesh: Mesh, youngs: float, poisson: float, thickness: float,
                     stress_state: str, forces: np.ndarray, prescribed: List[Dict],
                     ) -> np.ndarray:
    """Solve the same problem with scikit-fem; returns (num_nodes, dim)."""
    import skfem
    from skfem import Basis, ElementHex1, ElementQuad1, ElementVector, MeshHex, MeshQuad
    from skfem import asm, condense, solve
    from skfem.models.elasticity import lame_parameters, linear_elasticity

    dim = mesh.dim
    # scikit-fem wants (dim, nodes) and (corners, elements) arrays, C-ordered.
    p = np.ascontiguousarray(mesh.nodes.T.astype(float))
    if dim == 2:
        # scikit-fem's quad corners are ordered (0,0), (0,1), (1,1), (1,0) on
        # its reference square; SparLab stores them counter-clockwise, so the
        # second and fourth corners swap.
        t = np.ascontiguousarray(mesh.elements[:, [0, 3, 2, 1]].T)
        m = MeshQuad(p, t)
        element = ElementVector(ElementQuad1())
    else:
        # scikit-fem hex corners: (0,0,0), (0,1,0), (1,0,0), (0,0,1), (1,1,0),
        # (0,1,1), (1,0,1), (1,1,1) in reference coordinates; SparLab uses the
        # VTK order (bottom face counter-clockwise, then the top face).
        t = np.ascontiguousarray(mesh.elements[:, [0, 3, 1, 4, 2, 7, 5, 6]].T)
        m = MeshHex(p, t)
        element = ElementVector(ElementHex1())
    basis = Basis(m, element)

    lam, mu = lame_parameters(youngs, poisson)
    if stress_state == "plane_stress":
        lam = 2.0 * lam * mu / (lam + 2.0 * mu)
    elif stress_state not in ("plane_strain", "three_dimensional"):
        raise ResultError(f"unsupported stress state {stress_state}")
    k = asm(linear_elasticity(lam, mu), basis)
    if dim == 2:
        k = k * thickness

    # Nodal loads and prescribed DOFs mapped through skfem's DOF numbering.
    nodal_dofs = basis.nodal_dofs  # (dim, num_nodes)
    f = np.zeros(basis.N)
    for row in forces:
        node = int(row[0])
        for c in range(dim):
            f[nodal_dofs[c, node]] += row[1 + c]
    x = np.zeros(basis.N)
    fixed = []
    for entry in prescribed:
        c = "xyz".index(entry["component"])
        dof = nodal_dofs[c, int(entry["node"])]
        fixed.append(dof)
        x[dof] = float(entry["value_m"])
    fixed = np.unique(np.asarray(fixed, dtype=int))
    u = solve(*condense(k, f, x=x, D=fixed))
    out = np.zeros((mesh.num_nodes, dim))
    for c in range(dim):
        out[:, c] = u[nodal_dofs[c]]
    return out


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
def compare(reference: np.ndarray, ours: np.ndarray) -> Dict[str, float]:
    diff = reference - ours
    scale = float(np.abs(ours).max())
    per_node = np.linalg.norm(diff, axis=1) / max(scale, 1e-300)
    return {
        "max_abs_diff_m": float(np.abs(diff).max()),
        "max_rel_diff": float(np.abs(diff).max() / max(scale, 1e-300)),
        "rms_rel_diff": float(np.sqrt(np.mean(per_node ** 2))),
        "ref_max_abs_m": float(np.abs(reference).max()),
        "sparlab_max_abs_m": scale,
    }


def sparlab_displacement(case, load_case: str) -> np.ndarray:
    table = case.displacement(load_case)
    cols = ["ux[m]", "uy[m]", "uz[m]"][: case.dim]
    return table[cols].to_numpy()


def cross_validate_case(case_dir: str, tolerances: Dict[str, float],
                        skip_calculix: bool) -> Dict:
    case = load_case(case_dir)
    mesh = case.mesh
    summary = case.summary
    material = summary["material"]
    thickness = float(summary["mesh"].get("thickness_m", 1.0))
    stress_state = material["stress_state"]
    element_type = summary["mesh"]["element_type"]
    report = {
        "case": case.name,
        "directory": case_dir,
        "dim": mesh.dim,
        "element_type": element_type,
        "num_nodes": mesh.num_nodes,
        "num_elements": mesh.num_elements,
        "stress_state": stress_state,
        "load_cases": [],
    }
    ccx_type = {"Quad4": "CPE4" if stress_state == "plane_strain" else "CPS4",
                "Hex8": "C3D8"}[element_type]

    for name in mesh.load_case_names:
        ours = sparlab_displacement(case, name)
        entry = {"load_case": name, "codes": {}}

        # --- scikit-fem ---
        forces = mesh.nodal_forces(name)
        sk = solve_with_skfem(mesh, float(material["youngs_modulus_Pa"]),
                              float(material["poisson_ratio"]), thickness, stress_state,
                              forces, mesh.prescribed)
        stats = compare(sk, ours)
        stats["tolerance"] = tolerances["skfem"]
        stats["passed"] = stats["max_rel_diff"] <= tolerances["skfem"]
        stats["element"] = {"Quad4": "ElementQuad1", "Hex8": "ElementHex1"}[element_type]
        entry["codes"]["scikit-fem"] = stats

        # --- CalculiX ---
        if not skip_calculix:
            deck = case.path(f"calculix_{_safe(name)}.inp")
            if not os.path.isfile(deck):
                raise ResultError(f"{deck} is missing; rerun sparlab_solve with "
                                  "--export-calculix")
            with tempfile.TemporaryDirectory(prefix="sparlab_ccx_") as work:
                frd = run_calculix(deck, work)
                ref = parse_frd_displacements(frd, mesh.num_nodes)[:, : mesh.dim]
            if np.isnan(ref).any():
                raise ResultError(f"CalculiX returned no displacement for some nodes of {name}")
            key = "calculix_solid" if ccx_type == "C3D8" else "calculix_plane"
            stats = compare(ref, ours)
            stats["tolerance"] = tolerances[key]
            stats["passed"] = stats["max_rel_diff"] <= tolerances[key]
            stats["element"] = ccx_type
            # Six significant digits in the .frd file: half a unit in the last
            # digit of the largest value, relative to that value.
            stats["frd_rounding_floor_rel"] = 5.0e-6
            stats["within_frd_rounding"] = stats["max_rel_diff"] <= 5.0e-6
            entry["codes"]["calculix"] = stats
        report["load_cases"].append(entry)
    return report


def _safe(name: str) -> str:
    return "".join(c if (c.isalnum() or c in "-_") else "_" for c in name) or "unnamed"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", action="append", required=True,
                        help="result directory of a sparlab_solve run (repeatable)")
    parser.add_argument("--output", default="results/cross_validation")
    parser.add_argument("--skip-calculix", action="store_true")
    parser.add_argument("--tol-skfem", type=float, default=1e-7)
    parser.add_argument("--tol-calculix-solid", type=float, default=1e-5)
    parser.add_argument("--tol-calculix-plane", type=float, default=1e-5)
    args = parser.parse_args(argv)

    if not args.skip_calculix and shutil.which("ccx") is None:
        print("ccx (CalculiX) is not on PATH; pass --skip-calculix to compare with "
              "scikit-fem only", file=sys.stderr)
        return 2

    tolerances = {"skfem": args.tol_skfem, "calculix_solid": args.tol_calculix_solid,
                  "calculix_plane": args.tol_calculix_plane}
    import skfem
    summary = {
        "about": ("Node-by-node comparison of SparLab nodal displacements with "
                  "independent solutions of the same discrete problem"),
        "codes": {
            "scikit-fem": {"version": skfem.__version__,
                           "note": "same element formulation (bilinear/trilinear "
                                   "isoparametric, full integration); differences are "
                                   "linear-solver round-off"},
            "calculix": {"version": None if args.skip_calculix else calculix_version(),
                         "note": "C3D8 is the same element as SparLab's Hex8; CPS4/CPE4 "
                                 "are plane elements CalculiX expands through the "
                                 "thickness, a different discretisation of the plane "
                                 "problem. Nodal results are read from the .frd file, "
                                 "which carries six significant digits, so differences "
                                 "below 5e-6 relative are its rounding, not a "
                                 "disagreement"},
        },
        "tolerances": tolerances,
        "cases": [],
    }
    all_passed = True
    for case_dir in args.case:
        try:
            report = cross_validate_case(case_dir, tolerances, args.skip_calculix)
        except (ResultError, FileNotFoundError, KeyError) as error:
            print(f"cross-validation failed for {case_dir}: {error}", file=sys.stderr)
            return 1
        summary["cases"].append(report)
        for lc in report["load_cases"]:
            for code, stats in lc["codes"].items():
                all_passed = all_passed and stats["passed"]
                print(f"  {report['case']:<24} {lc['load_case']:<14} {code:<11} "
                      f"{stats['element']:<12} max rel diff {stats['max_rel_diff']:.3e} "
                      f"(tol {stats['tolerance']:.0e}) {'PASS' if stats['passed'] else 'FAIL'}")
    summary["all_passed"] = all_passed
    os.makedirs(args.output, exist_ok=True)
    path = os.path.join(args.output, "summary.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(f"  wrote {path}")
    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
