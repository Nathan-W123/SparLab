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
    isotropic elasticity form on ElementQuad1 / ElementHex1 / ElementTriP1 /
    ElementTetP1, and on ElementTetP2 over an isoparametric MeshTet2 for the
    quadratic tetrahedra (curved cells included),

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
relative. The same holds for C3D4, the linear tetrahedron, and for C3D10, the
quadratic one (four-point stiffness integration, as in SparLab). The CalculiX
tolerance is therefore 1e-5 for every element family, and the summary records
the rounding floor next to the measured difference. CalculiX's plane elements
(CPS4, CPS3) are *expanded* into solids through the thickness with the
plane-stress condition imposed on the expanded element, which is a different
discretisation of the plane problem. It coincides with plane stress only
for a Poisson ratio of zero: measured on the Gmsh lug bracket, CPS3 agrees to
the .frd rounding floor at nu = 0 (2.8e-6) and differs by about 1e-3 at
nu = 0.33. A plane comparison with nu != 0 is therefore reported as an
informational comparison between two idealisations, and scikit-fem - the same
plane element - is the verification for those cases.

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


def parse_dat_buckling_factors(path: str) -> np.ndarray:
    """Buckling factors from the "B U C K L I N G   F A C T O R" block of a
    CalculiX .dat file."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    factors: List[float] = []
    for i, line in enumerate(lines):
        if "B U C K L I N G" in line:
            for row in lines[i + 1:]:
                parts = row.split()
                if len(parts) == 2 and parts[0].isdigit():
                    factors.append(float(parts[1]))
                elif factors and not row.strip():
                    break
            break
    if not factors:
        raise ResultError(f"{path} holds no buckling factors")
    return np.asarray(factors)


def run_calculix_buckling(static_deck: str, workdir: str, num_modes: int) -> np.ndarray:
    """The static deck turned into a *BUCKLE step (same mesh, supports and
    nodal loads, which define the load pattern) and solved by ccx."""
    with open(static_deck, "r", encoding="utf-8") as handle:
        text = handle.read()
    if "*STATIC\n" not in text:
        raise ResultError(f"{static_deck} has no *STATIC step to turn into *BUCKLE")
    text = text.replace("*STATIC\n", f"*BUCKLE\n{num_modes}, 1.E-12, "
                                      f"{max(4 * num_modes, 40)}, 10000\n", 1)
    deck = os.path.join(workdir, "buckle.inp")
    with open(deck, "w", encoding="utf-8") as handle:
        handle.write(text)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(["ccx", "-i", "buckle"], cwd=workdir, capture_output=True,
                          text=True, env=env, check=False)
    dat = os.path.join(workdir, "buckle.dat")
    if proc.returncode != 0 or not os.path.isfile(dat):
        raise ResultError(f"ccx *BUCKLE failed on {static_deck}:\n{proc.stdout[-2000:]}")
    return parse_dat_buckling_factors(dat)


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
class SkfemProblem:
    """The discrete problem of a SparLab run rebuilt in scikit-fem: basis,
    stiffness, nodal loads and prescribed DOFs in scikit-fem's numbering."""

    def __init__(self, mesh: Mesh, youngs: float, poisson: float, thickness: float,
                 stress_state: str, forces: np.ndarray, prescribed: List[Dict]):
        from skfem import (Basis, ElementHex1, ElementQuad1, ElementTetP1, ElementTetP2,
                           ElementTriP1, ElementVector, MeshHex, MeshQuad, MeshTet, MeshTet2,
                           MeshTri)
        from skfem import asm
        from skfem.models.elasticity import lame_parameters, linear_elasticity

        dim = mesh.dim
        self.mesh = mesh
        self.dim = dim
        self.thickness = thickness if dim == 2 else 1.0
        # scikit-fem wants (dim, nodes) and (corners, elements) arrays, C-ordered.
        p = np.ascontiguousarray(mesh.nodes.T.astype(float))
        # SparLab node -> scikit-fem DOF per component; the identity map of the
        # nodal DOFs except for quadratic meshes, which scikit-fem renumbers.
        dof_of_node = None
        if mesh.element_type == "Tet10":
            # Quadratic tetrahedra, curved when the file's edge nodes are: an
            # isoparametric MeshTet2 built from all ten nodes of every cell.
            # With sort_t=False scikit-fem keeps the connectivity rows, and its
            # edge order - (0,1), (1,2), (2,0), (0,3), (1,3), (2,3) - is
            # SparLab's (VTK's), which the element_dofs check below confirms.
            # intorder=2 is the 4-point rule SparLab (and CalculiX's C3D10)
            # integrate the stiffness with, so curved cells are the same
            # discrete problem too.
            m = MeshTet2(p, np.ascontiguousarray(mesh.elements.T), sort_t=False)
            element = ElementVector(ElementTetP2())
            basis = Basis(m, element, intorder=2)
            scalar = Basis(m, ElementTetP2(), intorder=2)
            located = m.doflocs[:, scalar.element_dofs].transpose(2, 1, 0)
            expected = mesh.nodes[mesh.elements]
            if not np.allclose(located, expected, rtol=0.0,
                               atol=1e-12 * float(np.abs(mesh.nodes).max())):
                raise ResultError("scikit-fem's Tet10 node order does not match SparLab's")
            corners = np.unique(mesh.elements[:, :4])
            dof_of_node = np.full((mesh.num_nodes, dim), -1, dtype=int)
            dof_of_node[corners] = basis.nodal_dofs.T
            for j in range(6):
                dof_of_node[mesh.elements[:, 4 + j]] = basis.edge_dofs[:, m.t2e[j]].T
            if (dof_of_node < 0).any():
                raise ResultError("some Tet10 nodes have no scikit-fem DOF")
        elif mesh.element_type == "Tri3":
            # Linear simplices: any vertex order describes the same element.
            m = MeshTri(p, np.ascontiguousarray(mesh.elements.T))
            element = ElementVector(ElementTriP1())
        elif mesh.element_type == "Tet4":
            m = MeshTet(p, np.ascontiguousarray(mesh.elements.T))
            element = ElementVector(ElementTetP1())
        elif dim == 2:
            # scikit-fem's quad corners are ordered (0,0), (0,1), (1,1), (1,0) on
            # its reference square; SparLab stores them counter-clockwise, so
            # the second and fourth corners swap.
            t = np.ascontiguousarray(mesh.elements[:, [0, 3, 2, 1]].T)
            m = MeshQuad(p, t)
            element = ElementVector(ElementQuad1())
        else:
            # scikit-fem hex corners: (0,0,0), (0,1,0), (1,0,0), (0,0,1),
            # (1,1,0), (0,1,1), (1,0,1), (1,1,1) in reference coordinates;
            # SparLab uses the VTK order (bottom face counter-clockwise, then
            # the top face).
            t = np.ascontiguousarray(mesh.elements[:, [0, 3, 1, 4, 2, 7, 5, 6]].T)
            m = MeshHex(p, t)
            element = ElementVector(ElementHex1())
        if dof_of_node is None:
            basis = Basis(m, element)
            dof_of_node = basis.nodal_dofs.T
        self.basis = basis
        self.dof_of_node = dof_of_node

        lam, mu = lame_parameters(youngs, poisson)
        if stress_state == "plane_stress":
            lam = 2.0 * lam * mu / (lam + 2.0 * mu)
        elif stress_state not in ("plane_strain", "three_dimensional"):
            raise ResultError(f"unsupported stress state {stress_state}")
        self.lam, self.mu = lam, mu
        self.k = asm(linear_elasticity(lam, mu), basis) * self.thickness

        # Nodal loads and prescribed DOFs mapped through skfem's DOF numbering.
        self.f = np.zeros(basis.N)
        for row in forces:
            node = int(row[0])
            for c in range(dim):
                self.f[dof_of_node[node, c]] += row[1 + c]
        self.x = np.zeros(basis.N)
        fixed = []
        for entry in prescribed:
            c = "xyz".index(entry["component"])
            dof = dof_of_node[int(entry["node"]), c]
            fixed.append(dof)
            self.x[dof] = float(entry["value_m"])
        self.fixed = np.unique(np.asarray(fixed, dtype=int))
        self.free = np.setdiff1d(np.arange(basis.N), self.fixed)

    def static(self) -> np.ndarray:
        """The static solution in scikit-fem's numbering."""
        from skfem import condense, solve

        return solve(*condense(self.k, self.f, x=self.x, D=self.fixed))

    def nodal(self, u: np.ndarray) -> np.ndarray:
        out = np.zeros((self.mesh.num_nodes, self.dim))
        for c in range(self.dim):
            out[:, c] = u[self.dof_of_node[:, c]]
        return out

    def buckling(self, u: np.ndarray, num_modes: int) -> np.ndarray:
        """Lowest positive load factors of (K + lambda K_G(u)) phi = 0 on the
        free DOFs, K_G assembled here from the stress of `u` at the same
        quadrature points, by a dense generalised eigensolve."""
        import scipy.linalg
        from skfem import BilinearForm, asm
        from skfem.helpers import eye, grad, sym_grad, trace

        dim = self.dim
        eps = sym_grad(self.basis.interpolate(u))
        sigma = self.lam * eye(trace(eps), dim) + 2.0 * self.mu * eps

        @BilinearForm
        def geometric(du, dv, w):
            # sum_ijk sigma_ij d(du_k)/dx_i d(dv_k)/dx_j
            return np.einsum("ij...,ki...,kj...->...", w["s"], grad(du), grad(dv))

        kg = asm(geometric, self.basis, s=sigma) * self.thickness
        free = self.free
        kff = self.k[free][:, free].toarray()
        gff = kg[free][:, free].toarray()
        # 1/lambda from (-K_G) phi = (1/lambda) K phi, K positive definite.
        inverse = scipy.linalg.eigh(-gff, kff, eigvals_only=True)
        positive = inverse[inverse > 0.0]
        return np.sort(1.0 / positive)[:num_modes]


def solve_with_skfem(mesh: Mesh, youngs: float, poisson: float, thickness: float,
                     stress_state: str, forces: np.ndarray, prescribed: List[Dict],
                     ) -> np.ndarray:
    """Solve the same problem with scikit-fem; returns (num_nodes, dim)."""
    problem = SkfemProblem(mesh, youngs, poisson, thickness, stress_state, forces,
                           prescribed)
    return problem.nodal(problem.static())


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
    plane_strain = stress_state == "plane_strain"
    ccx_type = {"Quad4": "CPE4" if plane_strain else "CPS4",
                "Tri3": "CPE3" if plane_strain else "CPS3",
                "Hex8": "C3D8", "Tet4": "C3D4", "Tet10": "C3D10"}[element_type]

    buckling = {entry["load_case"]: entry
                for entry in (summary.get("buckling") or {}).get("load_cases", [])}
    for name in mesh.load_case_names:
        ours = sparlab_displacement(case, name)
        entry = {"load_case": name, "codes": {}}

        # --- scikit-fem ---
        forces = mesh.nodal_forces(name)
        problem = SkfemProblem(mesh, float(material["youngs_modulus_Pa"]),
                               float(material["poisson_ratio"]), thickness, stress_state,
                               forces, mesh.prescribed)
        u_sk = problem.static()
        sk = problem.nodal(u_sk)
        stats = compare(sk, ours)
        stats["tolerance"] = tolerances["skfem"]
        stats["passed"] = stats["max_rel_diff"] <= tolerances["skfem"]
        stats["element"] = {"Quad4": "ElementQuad1", "Hex8": "ElementHex1",
                            "Tri3": "ElementTriP1", "Tet4": "ElementTetP1",
                            "Tet10": "ElementTetP2"}[element_type]
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
            key = ("calculix_solid" if ccx_type in ("C3D8", "C3D4", "C3D10")
                   else "calculix_plane")
            stats = compare(ref, ours)
            stats["tolerance"] = tolerances[key]
            stats["element"] = ccx_type
            poisson = float(material["poisson_ratio"])
            if key == "calculix_plane" and poisson != 0.0:
                # Not the same discrete problem (see the module docstring):
                # recorded, not judged.
                stats["passed"] = None
                stats["comparison"] = (
                    f"different idealisation: CalculiX expands {ccx_type} into a 3-D "
                    f"layer, which matches plane stress only for nu = 0 (here nu = "
                    f"{poisson:g}); informational")
            else:
                stats["passed"] = stats["max_rel_diff"] <= tolerances[key]
                stats["comparison"] = "same discrete problem"
            # Six significant digits in the .frd file: half a unit in the last
            # digit of the largest value, relative to that value.
            stats["frd_rounding_floor_rel"] = 5.0e-6
            stats["within_frd_rounding"] = stats["max_rel_diff"] <= 5.0e-6
            entry["codes"]["calculix"] = stats

        # --- linear buckling, when the run computed it ---
        ours_lf = np.asarray(buckling.get(name, {}).get("load_factors") or [], dtype=float)
        if ours_lf.size:
            count = int(ours_lf.size)
            if problem.free.size <= tolerances["skfem_buckling_max_dofs"]:
                ref = problem.buckling(u_sk, count)
                stats = compare_factors(
                    ref, ours_lf, tolerances["skfem_buckling"],
                    "K_G assembled by scikit-fem from its own static solution at the "
                    "same quadrature points; dense generalised eigensolve")
                stats["element"] = entry["codes"]["scikit-fem"]["element"] + " K_G"
                entry["codes"]["scikit-fem buckling"] = stats
            else:
                entry["skipped"] = (f"scikit-fem buckling skipped: {problem.free.size} free "
                                    "DOFs exceed the dense-eigensolve limit")
            if not skip_calculix and mesh.dim == 3:
                deck = case.path(f"calculix_{_safe(name)}.inp")
                with tempfile.TemporaryDirectory(prefix="sparlab_ccx_buckle_") as work:
                    ref = run_calculix_buckling(deck, work, count)
                stats = compare_factors(
                    ref, ours_lf, tolerances["calculix_buckling"],
                    f"CalculiX *BUCKLE on the exported {ccx_type} deck (same mesh, "
                    "supports and nodal loads); its own stress-stiffness evaluation, so "
                    "an independent implementation rather than the same discrete problem")
                stats["element"] = f"{ccx_type} *BUCKLE"
                entry["codes"]["calculix buckling"] = stats
        report["load_cases"].append(entry)
    return report


def compare_factors(reference: np.ndarray, ours: np.ndarray, tolerance: float,
                    comparison: str) -> Dict:
    """Mode-by-mode relative difference of two sets of load factors."""
    count = min(reference.size, ours.size)
    if count == 0:
        raise ResultError("no load factors to compare")
    rel = np.abs(reference[:count] - ours[:count]) / np.abs(ours[:count])
    return {
        "element": "buckling",
        "reference_load_factors": reference[:count].tolist(),
        "sparlab_load_factors": ours[:count].tolist(),
        "max_rel_diff": float(rel.max()),
        "per_mode_rel_diff": rel.tolist(),
        "tolerance": tolerance,
        "passed": bool(rel.max() <= tolerance),
        "comparison": comparison,
    }


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
    parser.add_argument("--tol-skfem-buckling", type=float, default=1e-7)
    parser.add_argument("--tol-calculix-buckling", type=float, default=1e-4)
    parser.add_argument("--skfem-buckling-max-dofs", type=int, default=6000,
                        help="largest free-DOF count for the dense buckling eigensolve")
    args = parser.parse_args(argv)

    if not args.skip_calculix and shutil.which("ccx") is None:
        print("ccx (CalculiX) is not on PATH; pass --skip-calculix to compare with "
              "scikit-fem only", file=sys.stderr)
        return 2

    tolerances = {"skfem": args.tol_skfem, "calculix_solid": args.tol_calculix_solid,
                  "calculix_plane": args.tol_calculix_plane,
                  "skfem_buckling": args.tol_skfem_buckling,
                  "calculix_buckling": args.tol_calculix_buckling,
                  "skfem_buckling_max_dofs": args.skfem_buckling_max_dofs}
    import skfem
    summary = {
        "about": ("Node-by-node comparison of SparLab nodal displacements with "
                  "independent solutions of the same discrete problem"),
        "codes": {
            "scikit-fem": {"version": skfem.__version__,
                           "note": "same element formulations (bilinear/trilinear "
                                   "isoparametric with full integration, linear "
                                   "simplices, isoparametric quadratic tetrahedra with "
                                   "the 4-point rule); differences are linear-solver "
                                   "round-off"},
            "calculix": {"version": None if args.skip_calculix else calculix_version(),
                         "note": "C3D8, C3D4 and C3D10 are the same elements as "
                                 "SparLab's Hex8, Tet4 and Tet10; CPS4/CPS3 (CPE4/CPE3) are plane elements "
                                 "CalculiX expands through the thickness, a different "
                                 "discretisation of the plane problem. Nodal results are "
                                 "read from the .frd file, which carries six significant "
                                 "digits, so differences below 5e-6 relative are its "
                                 "rounding, not a disagreement"},
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
                if stats["passed"] is None:
                    verdict = "INFO (different idealisation)"
                else:
                    all_passed = all_passed and stats["passed"]
                    verdict = "PASS" if stats["passed"] else "FAIL"
                print(f"  {report['case']:<28} {lc['load_case']:<14} {code:<11} "
                      f"{stats['element']:<12} max rel diff {stats['max_rel_diff']:.3e} "
                      f"(tol {stats['tolerance']:.0e}) {verdict}")
    summary["all_passed"] = all_passed
    os.makedirs(args.output, exist_ok=True)
    path = os.path.join(args.output, "summary.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(f"  wrote {path}")
    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
