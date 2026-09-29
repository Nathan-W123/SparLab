"""Cross-validation of SparLab's MITC4 shell, used by cross_validate.py for a
run of Shell4 elements.

The same discrete problem is solved again by an independent implementation
of the element in NumPy, written from its equations (Dvorkin and Bathe's
MITC4 as Shell4.hpp states it), not translated from the C++:

  * the geometry X = sum N_k (x_k + zeta t/2 V_k) and the displacement
    u = sum N_k (u_k + zeta t/2 theta_k x V_k), with the directors V_k and
    thicknesses SparLab used (mesh.json `shell`), as 3 x 24 interpolation
    matrices and their derivatives in r, s and zeta;
  * the covariant strains 1/2 (g_i . u_,j + g_j . u_,i), the transverse
    shears replaced by the MITC interpolation of their values at the edge
    midpoints A, C (r-zeta) and B, D (s-zeta), turned into the local frame
    (e3 along g_zeta, e1 the projection of global x) by T = E^T G^-T;
  * plane stress with k G (k = 5/6) in transverse shear, 2 x 2 x 2 Gauss
    points; the drilling penalty 1/2 alpha G t int (n . theta - omega)^2 dA
    at the 2 x 2 points, omega = 1/2 n . (a^a x u_,a);
  * the consistent mass (mass_points x mass_points x 3 points) and the
    geometric stiffness of the in-plane stresses on the gradients of the
    whole displacement field;
  * a pressure's and self-weight's nodal forces integrated from the deck's
    definitions (whole-model or box regions), against SparLab's.

It then solves the static problem with the loads and prescribed DOFs of
mesh.json, the modal problem (the lowest frequencies, scipy's shift-invert
Lanczos) and the buckling problem of each static solution, and compares
with SparLab's results: displacements and rotations node by node,
frequencies and load factors mode by mode. The agreement is that of two
solves of one discrete problem, their round-off: judge_fields weighs it
against the system's round-off scale kappa_1 eps and the backward error of
SparLab's solution in the NumPy system.

CalculiX's S4 is a different discretisation - it expands each shell into a
layer of incompatible-mode solids (C3D8I) over normals it averages itself -
so its displacements are compared for information: the difference measures
two discretisations of the same shell, not the implementation.
"""

from __future__ import annotations

import json
import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

SHEAR_FACTOR = 5.0 / 6.0
CORNERS = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
_G2 = 1.0 / math.sqrt(3.0)


def _start(n: int) -> np.ndarray:
    """A fixed random start vector for ARPACK, whose default one is drawn
    afresh on every call: the results repeat bit for bit, and no symmetry of
    the model hides a mode from the Krylov space."""
    return np.random.default_rng(20240917).standard_normal(n)


def seeded_onenormest(operator) -> float:
    """scipy's onenormest (Hager and Higham's block estimate of the 1-norm)
    with its random starting block drawn from a fixed seed: it draws from
    numpy's global random state, which is restored afterwards, so the
    estimate repeats bit for bit from run to run."""
    state = np.random.get_state()
    np.random.seed(20240917)
    try:
        return float(spla.onenormest(operator))
    finally:
        np.random.set_state(state)


def gauss(n: int) -> Tuple[np.ndarray, np.ndarray]:
    return np.polynomial.legendre.leggauss(n)


def bilinear(r: float, s: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """N_k and their derivatives in r and s at (r, s)."""
    n = 0.25 * (1.0 + CORNERS[:, 0] * r) * (1.0 + CORNERS[:, 1] * s)
    dr = 0.25 * CORNERS[:, 0] * (1.0 + CORNERS[:, 1] * s)
    ds = 0.25 * CORNERS[:, 1] * (1.0 + CORNERS[:, 0] * r)
    return n, dr, ds


def skew(v: np.ndarray) -> np.ndarray:
    """[v]x with [v]x w = v x w."""
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


class Mitc4:
    """One MITC4 element: its corner coordinates x (4 x 3), unit directors
    v (4 x 3), thickness, plane-stress matrix and drilling factor."""

    def __init__(self, x: np.ndarray, v: np.ndarray, thickness: float, youngs: float,
                 poisson: float, drilling: float):
        self.x = np.asarray(x, dtype=float)
        self.v = np.asarray(v, dtype=float) / np.linalg.norm(v, axis=1)[:, None]
        self.h = 0.5 * thickness
        self.t = thickness
        g = youngs / (2.0 * (1.0 + poisson))
        self.shear = g
        c = youngs / (1.0 - poisson * poisson)
        self.d = np.array([[c, c * poisson, 0.0], [c * poisson, c, 0.0], [0.0, 0.0, g]])
        # (11, 22, 12, 23, 13): plane stress, then k G on both shears.
        self.c = np.zeros((5, 5))
        self.c[:3, :3] = self.d
        self.c[3, 3] = self.c[4, 4] = SHEAR_FACTOR * g
        self.kd = drilling * g * thickness
        # theta_k x V_k = -[V_k]x theta_k
        self.psi = [-skew(self.v[k]) for k in range(4)]

    # -- interpolation --------------------------------------------------------
    def interpolation(self, r: float, s: float, z: float):
        """The 3 x 24 matrices of u, u_,r, u_,s, u_,zeta at (r, s, z), and the
        base vectors g_r, g_s, g_zeta there."""
        n, dr, ds = bilinear(r, s)
        u = np.zeros((3, 24))
        ur = np.zeros((3, 24))
        us = np.zeros((3, 24))
        uz = np.zeros((3, 24))
        eye = np.eye(3)
        for k in range(4):
            cols = slice(6 * k, 6 * k + 3)
            rots = slice(6 * k + 3, 6 * k + 6)
            u[:, cols] = n[k] * eye
            u[:, rots] = z * self.h * n[k] * self.psi[k]
            ur[:, cols] = dr[k] * eye
            ur[:, rots] = z * self.h * dr[k] * self.psi[k]
            us[:, cols] = ds[k] * eye
            us[:, rots] = z * self.h * ds[k] * self.psi[k]
            uz[:, rots] = self.h * n[k] * self.psi[k]
        pos = self.x + z * self.h * self.v
        gr = dr @ pos
        gs = ds @ pos
        gz = self.h * (n @ self.v)
        return u, (ur, us, uz), (gr, gs, gz)

    @staticmethod
    def frame(gz: np.ndarray) -> np.ndarray:
        """Columns e1, e2, e3: e3 along g_zeta, e1 the projection of global x
        (global z when x is within 0.1 degree of e3)."""
        e3 = gz / np.linalg.norm(gz)
        e1 = np.array([1.0, 0.0, 0.0]) - e3[0] * e3
        if np.linalg.norm(e1) < math.sin(math.radians(0.1)):
            e1 = np.array([0.0, 0.0, 1.0]) - e3[2] * e3
        e1 /= np.linalg.norm(e1)
        return np.column_stack([e1, np.cross(e3, e1), e3])

    def covariant(self, r: float, s: float, z: float) -> Dict[str, np.ndarray]:
        _, (ur, us, uz), (gr, gs, gz) = self.interpolation(r, s, z)
        return {
            "rr": gr @ ur,
            "ss": gs @ us,
            "rs": 0.5 * (gr @ us + gs @ ur),
            "rz": 0.5 * (gr @ uz + gz @ ur),
            "sz": 0.5 * (gs @ uz + gz @ us),
        }

    def point(self, r: float, s: float, z: float, tying: Dict[float, Dict]):
        """The local strain operator (11, 22, 2*12, 2*23, 2*13; 5 x 24), the
        in-plane gradient operators D_1, D_2 (3 x 24), and det J."""
        u, grads, (gr, gs, gz) = self.interpolation(r, s, z)
        base = np.column_stack([gr, gs, gz])
        det = float(np.linalg.det(base))
        if det <= 0.0:
            raise ValueError("non-positive volume Jacobian")
        cov = self.covariant(r, s, z)
        a, c_, d_, b_ = tying[z]
        cov["rz"] = 0.5 * (1.0 + s) * a["rz"] + 0.5 * (1.0 - s) * c_["rz"]
        cov["sz"] = 0.5 * (1.0 + r) * d_["sz"] + 0.5 * (1.0 - r) * b_["sz"]
        tmat = self.frame(gz).T @ np.linalg.inv(base).T  # T(a, i) = e_a . g^i
        eps = [[cov["rr"], cov["rs"], cov["rz"]],
               [cov["rs"], cov["ss"], cov["sz"]],
               [cov["rz"], cov["sz"], None]]

        def local(a_: int, b: int) -> np.ndarray:
            out = np.zeros(24)
            for i in range(3):
                for j in range(3):
                    if eps[i][j] is not None:
                        out += tmat[a_, i] * tmat[b, j] * eps[i][j]
            return out

        bmat = np.vstack([local(0, 0), local(1, 1), 2.0 * local(0, 1), 2.0 * local(1, 2),
                          2.0 * local(0, 2)])
        d1 = sum(tmat[0, i] * grads[i] for i in range(3))
        d2 = sum(tmat[1, i] * grads[i] for i in range(3))
        return bmat, (d1, d2), det, u

    def tying_rows(self) -> Dict[float, Tuple]:
        rows = {}
        for z in (-_G2, _G2):
            rows[z] = (self.covariant(0.0, 1.0, z), self.covariant(0.0, -1.0, z),
                       self.covariant(1.0, 0.0, z), self.covariant(-1.0, 0.0, z))
        return rows

    # -- matrices -------------------------------------------------------------
    def drilling_row(self, r: float, s: float) -> Tuple[np.ndarray, float]:
        n, dr, ds = bilinear(r, s)
        gr = dr @ self.x
        gs = ds @ self.x
        cross = np.cross(gr, gs)
        area = float(np.linalg.norm(cross))
        normal = cross / area
        metric = np.array([[gr @ gr, gr @ gs], [gs @ gr, gs @ gs]])
        inv = np.linalg.inv(metric)
        ar = inv[0, 0] * gr + inv[0, 1] * gs
        as_ = inv[1, 0] * gr + inv[1, 1] * gs
        row = np.zeros(24)
        for k in range(4):
            row[6 * k + 3:6 * k + 6] = n[k] * normal
            # omega = 1/2 n . (a^r x u_,r + a^s x u_,s) = 1/2 u_,a . (n x a^a)
            row[6 * k:6 * k + 3] = -0.5 * (dr[k] * np.cross(normal, ar)
                                           + ds[k] * np.cross(normal, as_))
        return row, area

    def stiffness(self) -> np.ndarray:
        tying = self.tying_rows()
        ke = np.zeros((24, 24))
        for r in (-_G2, _G2):
            for s in (-_G2, _G2):
                for z in (-_G2, _G2):
                    bmat, _, det, _ = self.point(r, s, z, tying)
                    ke += det * bmat.T @ self.c @ bmat
        if self.kd > 0.0:
            for r in (-_G2, _G2):
                for s in (-_G2, _G2):
                    row, area = self.drilling_row(r, s)
                    ke += self.kd * area * np.outer(row, row)
        return 0.5 * (ke + ke.T)

    def mass(self, density: float, mass_points: int) -> np.ndarray:
        me = np.zeros((24, 24))
        xs, ws = gauss(mass_points)
        zs, wz = gauss(3)
        for r, wr in zip(xs, ws):
            for s, wsv in zip(xs, ws):
                for z, wzv in zip(zs, wz):
                    u, _, (gr, gs, gz) = self.interpolation(r, s, z)
                    det = float(np.linalg.det(np.column_stack([gr, gs, gz])))
                    me += density * wr * wsv * wzv * det * u.T @ u
        return 0.5 * (me + me.T)

    def geometric(self, ue: np.ndarray) -> np.ndarray:
        tying = self.tying_rows()
        kg = np.zeros((24, 24))
        for r in (-_G2, _G2):
            for s in (-_G2, _G2):
                for z in (-_G2, _G2):
                    bmat, (d1, d2), det, _ = self.point(r, s, z, tying)
                    sigma = self.c @ (bmat @ ue)
                    kg += det * (sigma[0] * d1.T @ d1 + sigma[1] * d2.T @ d2
                                 + sigma[2] * (d1.T @ d2 + d2.T @ d1))
        return 0.5 * (kg + kg.T)

    def pressure(self, p: float) -> np.ndarray:
        """-p int N_k (g_r x g_s) dr ds on the translations."""
        f = np.zeros(24)
        for r in (-_G2, _G2):
            for s in (-_G2, _G2):
                n, dr, ds = bilinear(r, s)
                area_vector = np.cross(dr @ self.x, ds @ self.x)
                for k in range(4):
                    f[6 * k:6 * k + 3] -= p * n[k] * area_vector
        return f

    def body(self, b: np.ndarray, mass_points: int) -> np.ndarray:
        """int U^T b dV for a constant body force density b."""
        f = np.zeros(24)
        xs, ws = gauss(mass_points)
        zs, wz = gauss(3)
        for r, wr in zip(xs, ws):
            for s, wsv in zip(xs, ws):
                for z, wzv in zip(zs, wz):
                    u, _, (gr, gs, gz) = self.interpolation(r, s, z)
                    det = float(np.linalg.det(np.column_stack([gr, gs, gz])))
                    f += wr * wsv * wzv * det * (u.T @ b)
        return f

    def centroid(self) -> np.ndarray:
        return self.x.mean(axis=0)


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------
class ShellProblem:
    def __init__(self, mesh_json: Dict, summary: Dict, deck: Dict):
        self.nodes = np.asarray(mesh_json["nodes_m"], dtype=float)
        self.elements = np.asarray(mesh_json["elements"], dtype=int)
        shell = mesh_json["shell"]
        self.thickness = np.asarray(shell["thickness_m"], dtype=float)
        self.directors = np.asarray(shell["directors"], dtype=float)  # ne x 4 x 3
        drilling = float(shell["drilling_factor"])
        materials = summary.get("materials") or [summary["material"]]
        per_element = mesh_json.get("element_materials") or [0] * len(self.elements)
        self.density = [float(materials[m].get("density_kg_per_m3", 0.0)) for m in per_element]
        integration = (deck.get("model") or {}).get("integration") or {}
        self.mass_points = int(integration.get("mass_points", 3))
        self.elems = [Mitc4(self.nodes[conn], self.directors[e], self.thickness[e],
                            float(materials[per_element[e]]["youngs_modulus_Pa"]),
                            float(materials[per_element[e]]["poisson_ratio"]), drilling)
                      for e, conn in enumerate(self.elements)]
        self.ndof = 6 * len(self.nodes)
        self.prescribed = {}
        component = {"x": 0, "y": 1, "z": 2, "rx": 3, "ry": 4, "rz": 5}
        for p in mesh_json["prescribed_dofs"]:
            value = p.get("value_m", p.get("value_rad", 0.0))
            self.prescribed[6 * int(p["node"]) + component[p["component"]]] = float(value)
        fixed = np.array(sorted(self.prescribed), dtype=int)
        self.fixed = fixed
        self.free = np.setdiff1d(np.arange(self.ndof), fixed)
        self.k = self._assemble(lambda e, el: el.stiffness())

    def _dofs(self, conn: np.ndarray) -> np.ndarray:
        return np.concatenate([np.arange(6 * n, 6 * n + 6) for n in conn])

    def _assemble(self, element_matrix) -> sp.csr_matrix:
        """Sum of element_matrix(e, element) over the elements."""
        rows, cols, vals = [], [], []
        for e, (el, conn) in enumerate(zip(self.elems, self.elements)):
            me = element_matrix(e, el)
            dofs = self._dofs(conn)
            rows.append(np.repeat(dofs, 24))
            cols.append(np.tile(dofs, 24))
            vals.append(me.ravel())
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.ndof, self.ndof))

    def condition_estimate(self) -> float:
        """The 1-norm condition number of K_ff, estimated as for the solid
        decks (Hager and Higham, scipy's onenormest, on a factorisation of
        K_ff): two backward-stable solves can differ by about kappa * eps of
        the largest value."""
        kff = self.k[self.free][:, self.free].tocsc()
        solve = spla.factorized(kff)
        n = kff.shape[0]
        inverse = spla.LinearOperator((n, n), matvec=solve, rmatvec=solve, dtype=float)
        return float(abs(kff).sum(axis=0).max() * seeded_onenormest(inverse))

    def backward_error(self, u: np.ndarray, f: np.ndarray) -> float:
        """||K u - f|| / (||K|| ||u|| + ||f||) on the free rows, ||K|| the
        largest absolute row sum: how far u is from solving this system."""
        kf = self.k[self.free]
        r = kf @ u - f[self.free]
        norm = float(np.abs(kf).sum(axis=1).max())
        return float(np.linalg.norm(r) / (norm * np.linalg.norm(u) + np.linalg.norm(f[self.free])))

    def static(self, f: np.ndarray) -> np.ndarray:
        u = np.zeros(self.ndof)
        for d, v in self.prescribed.items():
            u[d] = v
        kff = self.k[self.free][:, self.free].tocsc()
        rhs = f[self.free] - self.k[self.free][:, self.fixed] @ u[self.fixed]
        u[self.free] = spla.spsolve(kff, rhs)
        return u

    def modal(self, count: int) -> np.ndarray:
        """The lowest `count` angular frequencies of the consistent mass."""
        m = self._assemble(lambda e, el: el.mass(self.density[e], self.mass_points))
        kff = self.k[self.free][:, self.free].tocsc()
        mff = m[self.free][:, self.free].tocsc()
        vals = spla.eigsh(kff, k=count, M=mff, sigma=0.0, which="LM",
                          return_eigenvectors=False, tol=1e-13, v0=_start(kff.shape[0]))
        return np.sqrt(np.sort(vals))

    def buckling(self, u: np.ndarray, count: int) -> np.ndarray:
        """The lowest positive load factors of K phi = lambda (-K_G) phi, from
        the reversed pencil (-K_G) phi = mu K phi, K positive definite."""
        kg = self._assemble(lambda e, el: el.geometric(u[self._dofs(self.elements[e])]))
        kff = self.k[self.free][:, self.free].tocsc()
        gff = -kg[self.free][:, self.free].tocsc()
        vals = spla.eigsh(gff, k=count + 4, M=kff, which="LA", return_eigenvectors=False,
                          tol=1e-13, v0=_start(kff.shape[0]))
        mu = np.sort(vals[vals > 0.0])[::-1]
        return 1.0 / mu[:count]

    # -- loads from the deck --------------------------------------------------
    def region_elements(self, region: Optional[Dict]) -> Optional[List[int]]:
        """Elements whose centroid lies in a whole-model or box region; None
        for any other selector (not reproduced here)."""
        if region is None or region == {} or region.get("all") is True:
            return list(range(len(self.elems)))
        box = region.get("box")
        if box is None or set(region) - {"box", "name"}:
            return None
        diag = float(np.linalg.norm(self.nodes.max(axis=0) - self.nodes.min(axis=0)))
        tol = 1e-9 * diag
        lo = np.array([box.get(k, -np.inf) for k in ("xmin", "ymin", "zmin")]) - tol
        hi = np.array([box.get(k, np.inf) for k in ("xmax", "ymax", "zmax")]) + tol
        return [e for e, el in enumerate(self.elems)
                if np.all(el.centroid() >= lo) and np.all(el.centroid() <= hi)]

    def deck_loads(self, load_case: Dict) -> Tuple[Optional[np.ndarray], List[str]]:
        """Pressures and self-weight integrated here; None when the case has
        another kind of load or a region this script does not select."""
        others = {"point_loads", "tractions", "body_forces", "centrifugal", "temperature"}
        if any(load_case.get(k) for k in others):
            return None, []
        f = np.zeros(self.ndof)
        kinds = []
        for p in load_case.get("pressures") or []:
            elements = self.region_elements(p.get("region"))
            if elements is None:
                return None, []
            for e in elements:
                f[self._dofs(self.elements[e])] += self.elems[e].pressure(float(p["pressure"]))
            kinds.append("pressure")
        gravity = load_case.get("gravity")
        if gravity is not None and np.linalg.norm(gravity) > 0.0:
            g = np.asarray(gravity, dtype=float)
            for e, el in enumerate(self.elems):
                f[self._dofs(self.elements[e])] += el.body(self.density[e] * g, self.mass_points)
            kinds.append("self-weight")
        return (f, kinds) if kinds else (None, [])


def _relative(reference: np.ndarray, ours: np.ndarray) -> Dict[str, float]:
    diff = np.abs(reference - ours)
    scale = max(float(np.abs(reference).max()), 1e-300)
    return {"max_abs_diff": float(diff.max()), "max_rel_diff": float(diff.max() / scale),
            "ref_max_abs": scale}


def compare_fields(reference: np.ndarray, ours: np.ndarray, size: float) -> Dict:
    """Translations over their largest value; rotations over the larger of
    their largest value and the largest translation over the model's size
    (so that a field without rotations - an in-plane load - is not judged on
    its round-off)."""
    ref = reference.reshape(-1, 6)
    our = ours.reshape(-1, 6)
    t = _relative(ref[:, :3], our[:, :3])
    rot_scale = max(float(np.abs(ref[:, 3:]).max()), t["ref_max_abs"] / size, 1e-300)
    rot_diff = float(np.abs(ref[:, 3:] - our[:, 3:]).max())
    return {"max_abs_diff_m": t["max_abs_diff"],
            "max_rel_diff": max(t["max_rel_diff"], rot_diff / rot_scale),
            "translation_rel_diff": t["max_rel_diff"], "rotation_rel_diff": rot_diff / rot_scale,
            "max_abs_diff_rad": rot_diff, "ref_max_abs_m": t["ref_max_abs"],
            "rotation_scale_rad": rot_scale}


def judge_fields(stats: Dict, tolerance: float, round_off: float, backward: float) -> None:
    """Translations within the tolerance; rotations within the tolerance or
    ten times the round-off scale kappa eps of the system; and SparLab's
    solution solving the NumPy system to round-off. The drilling rotations
    are held by a penalty 1e-3 G t against the E t of the stiffest DOFs, so
    they are where the system's conditioning shows: each backward-stable
    solution lies within about kappa eps of the exact one there, two of them
    within twice that, and kappa is an estimate - hence the factor ten. The
    backward error is the proof that both solve one system."""
    stats["round_off_scale"] = round_off
    stats["sparlab_backward_error"] = backward
    stats["passed"] = (stats["translation_rel_diff"] <= tolerance
                       and stats["rotation_rel_diff"] <= max(tolerance, 10.0 * round_off)
                       and backward <= 1.0e-13)


def compare_values(reference: np.ndarray, ours: np.ndarray) -> Dict:
    n = min(reference.size, ours.size)
    rel = np.abs(reference[:n] - ours[:n]) / np.abs(reference[:n])
    return {"modes": int(n), "max_rel_diff": float(rel.max()),
            "reference": reference[:n].tolist(), "sparlab": ours[:n].tolist()}


def shell_case_report(case, case_dir: str, tolerances: Dict[str, float], skip_calculix: bool,
                      run_calculix, parse_frd_displacements, safe,
                      run_calculix_buckling=None) -> Dict:
    """The comparisons of a shell run (see the module docstring)."""
    mesh_json = json.load(open(os.path.join(case_dir, "mesh.json"), encoding="utf-8"))
    summary = case.summary
    deck = json.load(open(os.path.join(case_dir, "config.json"), encoding="utf-8"))
    problem = ShellProblem(mesh_json, summary, deck)
    size = float(np.linalg.norm(problem.nodes.max(axis=0) - problem.nodes.min(axis=0)))
    kappa = problem.condition_estimate()
    round_off = kappa * float(np.finfo(float).eps)
    report = {
        "case": case.name, "directory": case_dir, "dim": 3, "element_type": "Shell4",
        "num_nodes": int(problem.nodes.shape[0]), "num_elements": len(problem.elems),
        "stress_state": "shell", "load_cases": [],
    }
    load_decks = {lc.get("name"): lc for lc in deck.get("load_cases", [])}
    buckling = {entry["load_case"]: entry
                for entry in (summary.get("buckling") or {}).get("load_cases", [])}
    for lc in mesh_json["load_cases"]:
        name = lc["name"]
        f = np.zeros(problem.ndof)
        for nf in lc["nodal_forces"]:
            n = int(nf["node"])
            for k, key in enumerate(("fx_N", "fy_N", "fz_N", "mx_Nm", "my_Nm", "mz_Nm")):
                f[6 * n + k] += float(nf.get(key, 0.0))
        table = case.table(f"displacement_{safe(name)}.csv")
        ours = np.column_stack([table[c].to_numpy() for c in
                                ("ux[m]", "uy[m]", "uz[m]", "rx[rad]", "ry[rad]", "rz[rad]")])
        u_np = problem.static(f)
        entry = {"load_case": name, "codes": {}}
        stats = compare_fields(u_np, ours.ravel(), size)
        stats["tolerance"] = tolerances["shell_numpy"]
        stats["condition_estimate"] = kappa
        judge_fields(stats, tolerances["shell_numpy"], round_off,
                     problem.backward_error(ours.ravel(), f))
        stats["element"] = "NumPy MITC4"
        stats["comparison"] = ("an independent MITC4 in NumPy, SparLab's assembled load vector "
                               "(mesh.json); translations and rotations each over their largest "
                               "value")
        entry["codes"]["numpy mitc4"] = stats

        f_deck, kinds = problem.deck_loads(load_decks.get(name, {}))
        if f_deck is not None:
            stats = compare_fields(problem.static(f_deck), ours.ravel(), size)
            load_diff = float(np.abs(f_deck - f).max() / max(np.abs(f).max(), 1e-300))
            stats["load_vector_rel_diff"] = load_diff
            stats["tolerance"] = tolerances["shell_numpy"]
            judge_fields(stats, tolerances["shell_numpy"], round_off,
                         problem.backward_error(ours.ravel(), f_deck))
            stats["passed"] = stats["passed"] and load_diff <= tolerances["shell_numpy"]
            stats["element"] = "NumPy MITC4"
            stats["loads"] = kinds
            stats["comparison"] = ("the deck's " + " and ".join(kinds) + " integrated in NumPy "
                                   "(the mid-surface area vector; the body force through the "
                                   "shell's volume)")
            entry["codes"]["numpy mitc4 loads"] = stats

        ours_lf = np.asarray(buckling.get(name, {}).get("load_factors") or [], dtype=float)
        if ours_lf.size:
            ref = problem.buckling(u_np, int(ours_lf.size))
            stats = compare_values(ref, ours_lf)
            stats["tolerance"] = tolerances["shell_numpy_eigen"]
            stats["passed"] = stats["max_rel_diff"] <= tolerances["shell_numpy_eigen"]
            stats["element"] = "NumPy MITC4 K_G"
            stats["comparison"] = ("K_G of the NumPy static solution; scipy Lanczos on the "
                                   "reversed pencil")
            entry["codes"]["numpy mitc4 buckling"] = stats

        if not skip_calculix:
            ccx_deck = case.path(f"calculix_{safe(name)}.inp")
            if os.path.isfile(ccx_deck):
                import tempfile
                with tempfile.TemporaryDirectory(prefix="sparlab_ccx_shell_") as work:
                    frd = run_calculix(ccx_deck, work)
                    ref = parse_frd_displacements(frd, problem.nodes.shape[0])[:, :3]
                if np.isnan(ref).any():
                    raise RuntimeError(f"CalculiX returned no displacement for some nodes of {name}")
                t = _relative(ref, ours[:, :3])
                stats = {"max_abs_diff_m": t["max_abs_diff"], "max_rel_diff": t["max_rel_diff"],
                         "ref_max_abs_m": t["ref_max_abs"], "tolerance": None, "passed": None,
                         "element": "S4 (C3D8I expansion)",
                         "info": "different discretisation",
                         "comparison": ("a different discretisation: CalculiX expands each S4 "
                                        "into a layer of incompatible-mode solids over its own "
                                        "averaged normals, and applies nodal moments and held "
                                        "rotations through rigid knots; informational")}
                entry["codes"]["calculix"] = stats
                if ours_lf.size and run_calculix_buckling is not None:
                    with tempfile.TemporaryDirectory(prefix="sparlab_ccx_shell_buckle_") as work:
                        ref = run_calculix_buckling(ccx_deck, work, int(ours_lf.size))
                    stats = compare_values(ref, ours_lf)
                    stats["tolerance"] = None
                    stats["passed"] = None
                    stats["element"] = "S4 *BUCKLE"
                    stats["info"] = "different discretisation"
                    stats["comparison"] = ("CalculiX *BUCKLE on its S4 expansion: a different "
                                           "discretisation; informational")
                    entry["codes"]["calculix buckling"] = stats
        report["load_cases"].append(entry)

    frequencies = summary.get("modal")
    if frequencies:
        ours_f = np.asarray(frequencies.get("frequencies_hz") or [], dtype=float)
        if ours_f.size:
            ref = problem.modal(int(ours_f.size)) / (2.0 * math.pi)
            stats = compare_values(ref, ours_f)
            stats["tolerance"] = tolerances["shell_numpy_eigen"]
            stats["passed"] = stats["max_rel_diff"] <= tolerances["shell_numpy_eigen"]
            stats["element"] = "NumPy MITC4 mass"
            stats["comparison"] = ("the frequencies of the model (listed with its first load "
                                   "case): consistent mass, scipy shift-invert Lanczos at "
                                   "sigma = 0")
            report["load_cases"][0]["codes"]["numpy mitc4 modal"] = stats
    return report
