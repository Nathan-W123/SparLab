"""Cross-validation of SparLab's Timoshenko beam, used by cross_validate.py for
a run of Beam2 elements.

The same discrete frame is solved again by an independent implementation of
the element in NumPy, written from the beam's equations as Beam2.hpp states
them, not translated from the C++:

  * the local axes from the node coordinates and each section's orientation
    vector (y' its component normal to the axis; without one z' is the
    projection of global Z, or of X within 0.1 degree of vertical, and
    y' = z' x x');
  * the interdependent interpolation of Reddy (1997), written out here, and
    from it by 6-point Gauss integration the stiffness, the consistent mass
    (rotary inertia rho I included) and the geometric stiffness of the axial
    force N int [u'^2 + v'^2 + w'^2 + (I_p/A) theta_x'^2 + (I_y/A)
    theta_y'^2 + (I_z/A) theta_z'^2], N = E A (u_1 - u_0) / L;
  * the stiffness and the consistent mass checked against Przemieniecki's
    closed forms of the Timoshenko beam (Phi = 12 E I / (k G A L^2), 0
    without shear deformation);
  * the lumped mass by its definition: rho A L / 2 on each node's
    translations and half the sections' inertia tensor, rho L / 2
    R^T diag(I_p, I_y, I_z) R, on its rotations;
  * each element's distributed load, int H^T q, from mesh.json's `beam`
    block, with the nodal loads listed apart from it.

It then solves the static problem with the loads and prescribed DOFs of
mesh.json, the modal problem with the run's mass, the buckling problem of
each static solution and the steady harmonic response of a run with
`frequency_response`, and compares with SparLab's results: displacements
and rotations node by node, the end resultants element by element,
frequencies and load factors mode by mode, the complex monitors frequency by
frequency. The agreement is that of two solves of one discrete problem.

CalculiX: its linear beam B31 is expanded into bricks over the rectangle - a
3-D solid model of each member, not the Timoshenko beam - so its
displacements, *FREQUENCY and *BUCKLE are compared for information. Its U1
user element is a Timoshenko beam for static analysis, but in version 2.21
its shear term stiffens the beam instead of softening it (the deflection
falls below Euler-Bernoulli's as the shear coefficient falls); with a
coefficient of 1e12 it is the Euler-Bernoulli beam, with the torsion
constant I_y + I_z. A run whose sections have no shear deformation and that
torsion constant (a tube) is solved by a U1 deck written here and compared
to the 7 significant digits of CalculiX's printed output.
"""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import shell_xval

_PARALLEL = math.cos(0.1 * math.pi / 180.0)
_GAUSS = np.polynomial.legendre.leggauss(6)


def local_axes(x0: np.ndarray, x1: np.ndarray, orientation: np.ndarray) -> Tuple[float, np.ndarray]:
    """The length and the rotation R (rows x', y', z'; local = R global)."""
    d = x1 - x0
    length = float(np.linalg.norm(d))
    ex = d / length
    if np.linalg.norm(orientation) > 0.0:
        p = orientation - orientation.dot(ex) * ex
        ey = p / np.linalg.norm(p)
        ez = np.cross(ex, ey)
    else:
        up = np.array([1.0, 0.0, 0.0]) if abs(ex[2]) > _PARALLEL else np.array([0.0, 0.0, 1.0])
        ez = up - up.dot(ex) * ex
        ez /= np.linalg.norm(ez)
        ey = np.cross(ez, ex)
    return length, np.vstack([ex, ey, ez])


def bending_functions(xi: float, length: float, phi: float):
    """Reddy's interdependent interpolation of one bending plane in the DOFs
    (v_0, theta_0, v_1, theta_1), theta = v' without shear: the deflection
    and rotation functions and their derivatives in x."""
    mu = 1.0 / (1.0 + phi)
    l = length
    x2, x3 = xi * xi, xi ** 3
    nv = mu * np.array([1.0 - 3.0 * x2 + 2.0 * x3 + phi * (1.0 - xi),
                        l * (xi - 2.0 * x2 + x3 + 0.5 * phi * (xi - x2)),
                        3.0 * x2 - 2.0 * x3 + phi * xi,
                        l * (-x2 + x3 - 0.5 * phi * (xi - x2))])
    dnv = mu * np.array([-6.0 * xi + 6.0 * x2 - phi,
                         l * (1.0 - 4.0 * xi + 3.0 * x2 + 0.5 * phi * (1.0 - 2.0 * xi)),
                         6.0 * xi - 6.0 * x2 + phi,
                         l * (-2.0 * xi + 3.0 * x2 - 0.5 * phi * (1.0 - 2.0 * xi))]) / l
    nt = mu * np.array([6.0 * (x2 - xi) / l,
                        1.0 - 4.0 * xi + 3.0 * x2 + phi * (1.0 - xi),
                        6.0 * (xi - x2) / l,
                        3.0 * x2 - 2.0 * xi + phi * xi])
    dnt = mu * np.array([6.0 * (2.0 * xi - 1.0) / l,
                         -4.0 + 6.0 * xi - phi,
                         6.0 * (1.0 - 2.0 * xi) / l,
                         6.0 * xi - 2.0 + phi]) / l
    return nv, dnv, nt, dnt


class Beam:
    """One Timoshenko beam element from mesh.json's per-element data."""

    def __init__(self, x0: np.ndarray, x1: np.ndarray, section: Dict):
        self.length, self.r = local_axes(np.asarray(x0, float), np.asarray(x1, float),
                                         np.asarray(section["orientation"], float))
        self.a = float(section["area_m2"])
        self.iy = float(section["iy_m4"])
        self.iz = float(section["iz_m4"])
        self.j = float(section["torsion_m4"])
        self.ky = float(section["shear_y"])
        self.kz = float(section["shear_z"])
        self.e = float(section["youngs_modulus_Pa"])
        self.g = float(section["shear_modulus_Pa"])
        self.rho = float(section["density_kg_m3"])
        l2 = self.length ** 2
        self.phi_y = 12.0 * self.e * self.iz / (self.ky * self.g * self.a * l2) if self.ky > 0 else 0.0
        self.phi_z = 12.0 * self.e * self.iy / (self.kz * self.g * self.a * l2) if self.kz > 0 else 0.0
        self.t = np.zeros((12, 12))
        for b in range(4):
            self.t[3 * b:3 * b + 3, 3 * b:3 * b + 3] = self.r

    # -- the interpolation ----------------------------------------------------
    def fields(self, xi: float) -> Tuple[np.ndarray, np.ndarray]:
        """H (6 x 12): u, v, w, theta_x, theta_y, theta_z at xi in [0, 1] from
        the local DOFs; and its derivative in x."""
        h = np.zeros((6, 12))
        dh = np.zeros((6, 12))
        l = self.length
        for row, (i0, i1) in ((0, (0, 6)), (3, (3, 9))):
            h[row, i0], h[row, i1] = 1.0 - xi, xi
            dh[row, i0], dh[row, i1] = -1.0 / l, 1.0 / l
        # x'-y' plane: (v, theta_z) at DOFs 1, 5, 7, 11.
        nv, dnv, nt, dnt = bending_functions(xi, l, self.phi_y)
        cols = [1, 5, 7, 11]
        h[1, cols], dh[1, cols] = nv, dnv
        h[5, cols], dh[5, cols] = nt, dnt
        # x'-z' plane: (w, -theta_y) at DOFs 2, 4, 8, 10.
        nv, dnv, nt, dnt = bending_functions(xi, l, self.phi_z)
        sign = np.array([1.0, -1.0, 1.0, -1.0])
        cols = [2, 4, 8, 10]
        h[2, cols], dh[2, cols] = sign * nv, sign * dnv
        h[4, cols], dh[4, cols] = -sign * nt, -sign * dnt
        return h, dh

    def _integrate(self, integrand) -> np.ndarray:
        out = np.zeros((12, 12))
        for p, w in zip(*_GAUSS):
            out += 0.5 * w * self.length * integrand(0.5 * (p + 1.0))
        return 0.5 * (out + out.T)

    # -- local matrices -------------------------------------------------------
    def stiffness_local(self) -> np.ndarray:
        d = np.diag([self.e * self.a, self.ky * self.g * self.a, self.kz * self.g * self.a,
                     self.g * self.j, self.e * self.iy, self.e * self.iz])

        def strains(xi):
            h, dh = self.fields(xi)
            b = dh.copy()
            b[1] = dh[1] - h[5]  # gamma_y = v' - theta_z
            b[2] = dh[2] + h[4]  # gamma_z = w' + theta_y
            return b.T @ d @ b

        return self._integrate(strains)

    def mass_local(self) -> np.ndarray:
        m = self.rho * np.diag([self.a, self.a, self.a, self.iy + self.iz, self.iy, self.iz])
        return self._integrate(lambda xi: self.fields(xi)[0].T @ m @ self.fields(xi)[0])

    def geometric_local(self, n: float) -> np.ndarray:
        w = np.diag([1.0, 1.0, 1.0, (self.iy + self.iz) / self.a, self.iy / self.a,
                     self.iz / self.a])
        return n * self._integrate(lambda xi: self.fields(xi)[1].T @ w @ self.fields(xi)[1])

    def load_local(self, q_local: np.ndarray) -> np.ndarray:
        out = np.zeros(12)
        for p, w in zip(*_GAUSS):
            h, _ = self.fields(0.5 * (p + 1.0))
            out += 0.5 * w * self.length * (h[:3].T @ q_local)
        return out

    # -- closed forms (Przemieniecki) ------------------------------------------
    def closed_form_stiffness(self) -> np.ndarray:
        l = self.length
        k = np.zeros((12, 12))
        ax, tq = self.e * self.a / l, self.g * self.j / l
        k[np.ix_([0, 6], [0, 6])] = ax * np.array([[1, -1], [-1, 1]])
        k[np.ix_([3, 9], [3, 9])] = tq * np.array([[1, -1], [-1, 1]])
        for cols, inertia, phi, s in (([1, 5, 7, 11], self.iz, self.phi_y, 1.0),
                                      ([2, 4, 8, 10], self.iy, self.phi_z, -1.0)):
            c = self.e * inertia / ((1.0 + phi) * l ** 3)
            block = np.array([[12, 6 * l, -12, 6 * l],
                              [6 * l, (4 + phi) * l * l, -6 * l, (2 - phi) * l * l],
                              [-12, -6 * l, 12, -6 * l],
                              [6 * l, (2 - phi) * l * l, -6 * l, (4 + phi) * l * l]], float)
            flip = np.array([1.0, s, 1.0, s])
            k[np.ix_(cols, cols)] = c * block * np.outer(flip, flip)
        return k

    def closed_form_mass(self) -> np.ndarray:
        """Przemieniecki's consistent mass of the Timoshenko beam with the
        rotary inertia, in terms of Phi."""
        l = self.length
        m = np.zeros((12, 12))
        rho = self.rho
        m[np.ix_([0, 6], [0, 6])] = rho * self.a * l / 6.0 * np.array([[2, 1], [1, 2]])
        m[np.ix_([3, 9], [3, 9])] = rho * (self.iy + self.iz) * l / 6.0 * np.array([[2, 1], [1, 2]])
        for cols, inertia, phi, s in (([1, 5, 7, 11], self.iz, self.phi_y, 1.0),
                                      ([2, 4, 8, 10], self.iy, self.phi_z, -1.0)):
            f = 1.0 / (1.0 + phi) ** 2
            p2 = phi * phi
            t11 = 13 / 35 + 7 * phi / 10 + p2 / 3
            t12 = (11 / 210 + 11 * phi / 120 + p2 / 24) * l
            t13 = 9 / 70 + 3 * phi / 10 + p2 / 6
            t14 = -(13 / 420 + 3 * phi / 40 + p2 / 24) * l
            t22 = (1 / 105 + phi / 60 + p2 / 120) * l * l
            t24 = -(1 / 140 + phi / 60 + p2 / 120) * l * l
            trans = rho * self.a * l * f * np.array([[t11, t12, t13, t14],
                                                     [t12, t22, -t14, t24],
                                                     [t13, -t14, t11, -t12],
                                                     [t14, t24, -t12, t22]])
            r11 = 6 / 5
            r12 = (1 / 10 - phi / 2) * l
            r22 = (2 / 15 + phi / 6 + p2 / 3) * l * l
            r24 = (-1 / 30 - phi / 6 + p2 / 6) * l * l
            rot = rho * inertia / l * f * np.array([[r11, r12, -r11, r12],
                                                    [r12, r22, -r12, r24],
                                                    [-r11, -r12, r11, -r12],
                                                    [r12, r24, -r12, r22]])
            flip = np.array([1.0, s, 1.0, s])
            m[np.ix_(cols, cols)] = (trans + rot) * np.outer(flip, flip)
        return m

    # -- global matrices ------------------------------------------------------
    def stiffness(self) -> np.ndarray:
        return self.t.T @ self.stiffness_local() @ self.t

    def mass(self, lumped: bool) -> np.ndarray:
        if not lumped:
            return self.t.T @ self.mass_local() @ self.t
        half = 0.5 * self.rho * self.length
        out = np.zeros((12, 12))
        inertia = half * self.r.T @ np.diag([self.iy + self.iz, self.iy, self.iz]) @ self.r
        for a in range(2):
            out[6 * a:6 * a + 3, 6 * a:6 * a + 3] = half * self.a * np.eye(3)
            out[6 * a + 3:6 * a + 6, 6 * a + 3:6 * a + 6] = inertia
        return out

    def axial_force(self, ue: np.ndarray) -> float:
        d = self.t @ ue
        return self.e * self.a * (d[6] - d[0]) / self.length

    def geometric(self, ue: np.ndarray) -> np.ndarray:
        return self.t.T @ self.geometric_local(self.axial_force(ue)) @ self.t

    def load(self, q: np.ndarray) -> np.ndarray:
        return self.t.T @ self.load_local(self.r @ np.asarray(q, float))

    def end_forces(self, ue: np.ndarray, q: np.ndarray) -> np.ndarray:
        """(start, end) resultants N, Q_y, Q_z, T, M_y, M_z in local axes."""
        f = self.stiffness_local() @ (self.t @ ue) - self.load_local(self.r @ np.asarray(q, float))
        return np.concatenate([-f[:6], f[6:]])


class BeamProblem:
    def __init__(self, mesh_json: Dict):
        self.nodes = np.asarray(mesh_json["nodes_m"], dtype=float)
        self.elements = np.asarray(mesh_json["elements"], dtype=int)
        beam = mesh_json["beam"]
        self.sections = beam["elements"]
        self.elems = [Beam(self.nodes[c[0]], self.nodes[c[1]], s)
                      for c, s in zip(self.elements, self.sections)]
        self.cases = {lc["name"]: lc for lc in beam["load_cases"]}
        self.ndof = 6 * len(self.nodes)
        self.prescribed = {}
        component = {"x": 0, "y": 1, "z": 2, "rx": 3, "ry": 4, "rz": 5}
        for p in mesh_json["prescribed_dofs"]:
            value = p.get("value_m", p.get("value_rad", 0.0))
            self.prescribed[6 * int(p["node"]) + component[p["component"]]] = float(value)
        self.fixed = np.array(sorted(self.prescribed), dtype=int)
        self.free = np.setdiff1d(np.arange(self.ndof), self.fixed)
        self.k = self._assemble(lambda e, el: el.stiffness())
        # The closed forms, against the integrated interpolation.
        self.closed_form_diff = {"stiffness": 0.0, "mass": 0.0}
        for el in self.elems:
            for key, closed, integrated in (
                    ("stiffness", el.closed_form_stiffness(), el.stiffness_local()),
                    ("mass", el.closed_form_mass(), el.mass_local())):
                scale = np.abs(integrated).max()
                self.closed_form_diff[key] = max(self.closed_form_diff[key],
                                                 float(np.abs(closed - integrated).max() / scale))

    def dofs(self, e: int) -> np.ndarray:
        return np.concatenate([np.arange(6 * n, 6 * n + 6) for n in self.elements[e]])

    def _assemble(self, element_matrix) -> sp.csr_matrix:
        rows, cols, vals = [], [], []
        for e, el in enumerate(self.elems):
            me = element_matrix(e, el)
            d = self.dofs(e)
            rows.append(np.repeat(d, 12))
            cols.append(np.tile(d, 12))
            vals.append(me.ravel())
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.ndof, self.ndof))

    def mass(self, lumped: bool) -> sp.csr_matrix:
        return self._assemble(lambda e, el: el.mass(lumped))

    def load(self, name: str) -> Tuple[np.ndarray, np.ndarray]:
        """The load vector of a case: the nodal loads plus int H^T q; and the
        distributed part alone."""
        lc = self.cases[name]
        nodal = np.zeros(self.ndof)
        for row in lc["nodal_loads_node_fx_fy_fz_mx_my_mz"]:
            n = int(row[0])
            nodal[6 * n:6 * n + 6] += np.asarray(row[1:], dtype=float)
        distributed = np.zeros(self.ndof)
        for e, q in enumerate(lc["distributed_N_per_m"]):
            if np.any(np.asarray(q) != 0.0):
                distributed[self.dofs(e)] += self.elems[e].load(q)
        return nodal + distributed, distributed

    def condition_estimate(self) -> float:
        kff = self.k[self.free][:, self.free].tocsc()
        solve = spla.factorized(kff)
        n = kff.shape[0]
        inverse = spla.LinearOperator((n, n), matvec=solve, rmatvec=solve, dtype=float)
        return float(abs(kff).sum(axis=0).max() * shell_xval.seeded_onenormest(inverse))

    def backward_error(self, u: np.ndarray, f: np.ndarray) -> float:
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

    def modal(self, count: int, lumped: bool) -> np.ndarray:
        m = self.mass(lumped)
        kff = self.k[self.free][:, self.free].tocsc()
        mff = m[self.free][:, self.free].tocsc()
        vals = spla.eigsh(kff, k=count, M=mff, sigma=0.0, which="LM", return_eigenvectors=False,
                          tol=1e-13, v0=shell_xval._start(kff.shape[0]))
        return np.sqrt(np.sort(vals))

    def buckling(self, u: np.ndarray, count: int) -> np.ndarray:
        kg = self._assemble(lambda e, el: el.geometric(u[self.dofs(e)]))
        kff = self.k[self.free][:, self.free].tocsc()
        gff = -kg[self.free][:, self.free].tocsc()
        vals = spla.eigsh(gff, k=count + 4, M=kff, which="LA", return_eigenvectors=False,
                          tol=1e-13, v0=shell_xval._start(kff.shape[0]))
        mu = np.sort(vals[vals > 0.0])[::-1]
        return 1.0 / mu[:count]

    def end_forces(self, u: np.ndarray, name: str) -> np.ndarray:
        q = self.cases[name]["distributed_N_per_m"]
        return np.array([el.end_forces(u[self.dofs(e)], q[e]) for e, el in enumerate(self.elems)])

    def harmonic(self, f: np.ndarray, frequencies: List[float], lumped: bool, eta: float,
                 a: float, b: float) -> Tuple[np.ndarray, np.ndarray, sp.csr_matrix]:
        """U (frequencies x N) of [K (1 + i eta) - omega^2 M + i omega (a M + b K)] U = f,
        the prescribed DOFs held at their values; and the reactions."""
        m = self.mass(lumped).tocsr()
        k = self.k.tocsr()
        out = np.zeros((len(frequencies), self.ndof), dtype=complex)
        reactions = np.zeros_like(out)
        for j, freq in enumerate(frequencies):
            omega = 2.0 * math.pi * freq
            dyn = ((1.0 + 1j * (eta + omega * b)) * k + (-omega * omega + 1j * omega * a) * m).tocsr()
            u = np.zeros(self.ndof, dtype=complex)
            for d, v in self.prescribed.items():
                u[d] = v
            rhs = f[self.free] - dyn[self.free][:, self.fixed] @ u[self.fixed]
            u[self.free] = spla.spsolve(dyn[self.free][:, self.free].tocsc(), rhs)
            out[j] = u
            full = dyn @ u - f
            reactions[j, self.fixed] = full[self.fixed]
        return out, reactions, m


# ---------------------------------------------------------------------------
# CalculiX
# ---------------------------------------------------------------------------
def run_calculix_frequency(static_deck: str, workdir: str, num_modes: int) -> np.ndarray:
    """The static deck turned into a *FREQUENCY step: CalculiX's lowest
    frequencies [Hz] of its model (for B31, the expansion into bricks)."""
    text = open(static_deck, encoding="utf-8").read()
    if "*STATIC\n" not in text:
        raise RuntimeError(f"{static_deck} has no *STATIC step to turn into *FREQUENCY")
    text = text.replace("*STATIC\n", f"*FREQUENCY\n{num_modes}\n", 1)
    text = re.sub(r"\*CLOAD\n(?:[^*].*\n)*", "", text)
    text = re.sub(r"\*DLOAD\n(?:[^*].*\n)*", "", text)
    deck = os.path.join(workdir, "frequency.inp")
    with open(deck, "w", encoding="utf-8") as out:
        out.write(text)
    proc = subprocess.run(["ccx", "-i", "frequency"], cwd=workdir, capture_output=True, text=True,
                          timeout=1800)
    dat = os.path.join(workdir, "frequency.dat")
    if proc.returncode != 0 or not os.path.isfile(dat):
        raise RuntimeError(f"ccx *FREQUENCY failed on {static_deck}:\n{proc.stdout[-2000:]}")
    values = []
    lines = open(dat, encoding="utf-8").read().splitlines()
    for i, line in enumerate(lines):
        if "E I G E N V A L U E   O U T P U T" in line:
            for row in lines[i + 1:]:
                parts = row.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    values.append(float(parts[3]))  # cycles per time
                elif values and not row.strip():
                    break
            break
    if len(values) < num_modes:
        raise RuntimeError(f"ccx *FREQUENCY returned {len(values)} of {num_modes} frequencies")
    return np.asarray(values[:num_modes])


def write_u1_deck(path: str, problem: BeamProblem, f: np.ndarray, poisson: float) -> None:
    """A CalculiX deck of U1 elements (A, I11 = I_y, I12 = 0, I22 = I_z,
    k = 1e12: the Euler-Bernoulli limit), one section per element with its
    y' as the 1-direction, the nodal loads f and the prescribed DOFs."""
    def fmt(v: float) -> str:
        return repr(float(v))

    lines = ["*HEADING", "SparLab cross-validation: U1 Euler-Bernoulli frame", "*NODE, NSET=NALL"]
    for n, x in enumerate(problem.nodes):
        lines.append(f"{n + 1}, {fmt(x[0])}, {fmt(x[1])}, {fmt(x[2])}")
    lines.append("*USER ELEMENT, TYPE=U1, NODES=2, INTEGRATION POINTS=2, MAXDOF=6")
    lines.append("*ELEMENT, TYPE=U1, ELSET=EALL")
    for e, c in enumerate(problem.elements):
        lines.append(f"{e + 1}, {c[0] + 1}, {c[1] + 1}")
    materials = {}
    for e, (el, s) in enumerate(zip(problem.elems, problem.sections)):
        key = (s["youngs_modulus_Pa"], s["shear_modulus_Pa"])
        materials.setdefault(key, f"MAT{len(materials) + 1}")
        lines.append(f"*ELSET, ELSET=E{e + 1}\n{e + 1}")
    for (young, _), name in materials.items():
        lines += [f"*MATERIAL, NAME={name}", "*ELASTIC", f"{fmt(young)}, {fmt(poisson)}"]
    for e, (el, s) in enumerate(zip(problem.elems, problem.sections)):
        name = materials[(s["youngs_modulus_Pa"], s["shear_modulus_Pa"])]
        y = el.r[1]
        lines += [f"*BEAM SECTION, ELSET=E{e + 1}, MATERIAL={name}, SECTION=GENERAL",
                  f"{fmt(el.a)}, {fmt(el.iy)}, 0., {fmt(el.iz)}, 1.e12",
                  f"{fmt(y[0])}, {fmt(y[1])}, {fmt(y[2])}"]
    lines += ["*STEP", "*STATIC", "*BOUNDARY"]
    for d, v in sorted(problem.prescribed.items()):
        lines.append(f"{d // 6 + 1}, {d % 6 + 1}, {d % 6 + 1}, {fmt(v)}")
    lines.append("*CLOAD")
    for d in np.nonzero(f)[0]:
        lines.append(f"{d // 6 + 1}, {d % 6 + 1}, {fmt(f[d])}")
    lines += ["*NODE PRINT, NSET=NALL", "U", "*END STEP"]
    with open(path, "w", encoding="utf-8") as out:
        out.write("\n".join(lines) + "\n")


def run_u1(problem: BeamProblem, f: np.ndarray, poisson: float, workdir: str) -> np.ndarray:
    """CalculiX's U1 displacements and rotations (N x 6) from its .dat file
    (7 significant digits)."""
    write_u1_deck(os.path.join(workdir, "u1.inp"), problem, f, poisson)
    proc = subprocess.run(["ccx", "-i", "u1"], cwd=workdir, capture_output=True, text=True,
                          timeout=1800)
    dat = os.path.join(workdir, "u1.dat")
    if proc.returncode != 0 or not os.path.isfile(dat):
        raise RuntimeError(f"ccx U1 failed:\n{proc.stdout[-2000:]}")
    out = np.full((len(problem.nodes), 6), np.nan)
    started = False
    for line in open(dat, encoding="utf-8"):
        if "displacements" in line:
            started = True
            continue
        parts = line.split()
        if started and len(parts) == 7 and parts[0].isdigit():
            out[int(parts[0]) - 1] = [float(v) for v in parts[1:]]
    if np.isnan(out).any():
        raise RuntimeError("CalculiX U1 printed no displacement for some nodes")
    return out


def u1_obstacle(problem: BeamProblem) -> str:
    """Why the run's beams are not U1's Euler-Bernoulli beam, or ''."""
    for s, el in zip(problem.sections, problem.elems):
        if float(s["shear_y"]) != 0.0 or float(s["shear_z"]) != 0.0:
            return "a section deforms in shear, which U1 2.21 gets wrong"
        if abs(el.j - (el.iy + el.iz)) > 1e-12 * el.j:
            return "a torsion constant differs from I_y + I_z, U1's"
    return ""


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------
def beam_case_report(case, case_dir: str, tolerances: Dict[str, float], skip_calculix: bool,
                     run_calculix, parse_frd_displacements, safe,
                     run_calculix_buckling=None) -> Dict:
    """The comparisons of a beam run (see the module docstring)."""
    mesh_json = json.load(open(os.path.join(case_dir, "mesh.json"), encoding="utf-8"))
    summary = case.summary
    problem = BeamProblem(mesh_json)
    size = float(np.linalg.norm(problem.nodes.max(axis=0) - problem.nodes.min(axis=0)))
    kappa = problem.condition_estimate()
    round_off = kappa * float(np.finfo(float).eps)
    poisson = float(summary["material"]["poisson_ratio"])
    report = {
        "case": case.name, "directory": case_dir, "dim": 3, "element_type": "Beam2",
        "num_nodes": int(problem.nodes.shape[0]), "num_elements": len(problem.elems),
        "stress_state": "beam", "load_cases": [],
        "closed_form_rel_diff": problem.closed_form_diff,
    }
    buckling = {entry["load_case"]: entry
                for entry in (summary.get("buckling") or {}).get("load_cases", [])}
    harmonic = {entry["load_case"]: entry
                for entry in (summary.get("frequency_response") or {}).get("load_cases", [])}
    tol = tolerances["beam_numpy"]
    tol_eigen = tolerances["beam_numpy_eigen"]
    for lc in mesh_json["load_cases"]:
        name = lc["name"]
        f_sparlab = np.zeros(problem.ndof)
        for nf in lc["nodal_forces"]:
            n = int(nf["node"])
            for k, key in enumerate(("fx_N", "fy_N", "fz_N", "mx_Nm", "my_Nm", "mz_Nm")):
                f_sparlab[6 * n + k] += float(nf.get(key, 0.0))
        f, _ = problem.load(name)
        load_diff = float(np.abs(f - f_sparlab).max() / max(np.abs(f_sparlab).max(), 1e-300))
        table = case.table(f"displacement_{safe(name)}.csv")
        ours = np.column_stack([table[c].to_numpy() for c in
                                ("ux[m]", "uy[m]", "uz[m]", "rx[rad]", "ry[rad]", "rz[rad]")])
        u_np = problem.static(f)
        entry = {"load_case": name, "codes": {}}
        stats = shell_xval.compare_fields(u_np, ours.ravel(), size)
        stats["tolerance"] = tol
        stats["condition_estimate"] = kappa
        stats["load_vector_rel_diff"] = load_diff
        shell_xval.judge_fields(stats, tol, round_off, problem.backward_error(ours.ravel(), f))
        stats["passed"] = stats["passed"] and load_diff <= tol
        stats["element"] = "NumPy Timoshenko"
        stats["comparison"] = ("an independent Timoshenko frame in NumPy, its own load vector "
                               "(the nodal loads and int H^T q of mesh.json's distributed "
                               "loads); translations and rotations each over their largest "
                               "value")
        entry["codes"]["numpy timoshenko"] = stats

        # The end resultants, from the NumPy solution.
        forces = case.table(f"beam_{safe(name)}.csv")
        cols = [f"{q}{end}[{u}]" for end in "01" for q, u in
                (("N", "N"), ("Qy", "N"), ("Qz", "N"), ("T", "Nm"), ("My", "Nm"), ("Mz", "Nm"))]
        ours_f = np.column_stack([forces[c].to_numpy() for c in cols])
        ref_f = problem.end_forces(u_np, name)
        force_scale = max(float(np.abs(ref_f[:, [0, 1, 2, 6, 7, 8]]).max()), 1e-300)
        moment_scale = max(float(np.abs(ref_f[:, [3, 4, 5, 9, 10, 11]]).max()),
                           force_scale * size * 1e-3)
        diff = np.abs(ref_f - ours_f)
        rel = max(float(diff[:, [0, 1, 2, 6, 7, 8]].max()) / force_scale,
                  float(diff[:, [3, 4, 5, 9, 10, 11]].max()) / moment_scale)
        entry["codes"]["numpy timoshenko end forces"] = {
            "max_rel_diff": rel, "force_scale_N": force_scale, "moment_scale_Nm": moment_scale,
            "tolerance": tol, "passed": rel <= max(tol, 10.0 * round_off),
            "element": "NumPy Timoshenko",
            "comparison": ("K_e u_e - f_q of the NumPy solution in each element's axes, the "
                           "start's negated, against beam_<case>.csv; forces over the largest "
                           "force, moments over the largest moment"),
        }

        ours_lf = np.asarray(buckling.get(name, {}).get("load_factors") or [], dtype=float)
        if ours_lf.size:
            ref = problem.buckling(u_np, int(ours_lf.size))
            stats = shell_xval.compare_values(ref, ours_lf)
            stats["tolerance"] = tol_eigen
            stats["passed"] = stats["max_rel_diff"] <= tol_eigen
            stats["element"] = "NumPy Timoshenko K_G"
            stats["comparison"] = ("K_G of the NumPy static solution's axial forces; scipy "
                                   "Lanczos on the reversed pencil")
            entry["codes"]["numpy timoshenko buckling"] = stats

        fr = harmonic.get(name)
        if fr is not None:
            block = summary["frequency_response"]
            frequencies = [float(x) for x in block["frequencies_Hz"]]
            damping = block["damping"]
            u, reactions, _ = problem.harmonic(
                f, frequencies, block["mass"] == "lumped", float(damping["structural_loss_factor"]),
                float(damping["mass_1_per_s"]), float(damping["stiffness_s"]))
            omega = 2.0 * math.pi * np.asarray(frequencies)[:, None]
            table_fr = case.table(f"frequency_response_{safe(name)}.csv")
            monitor_diff = 0.0
            for monitor in fr.get("monitors", []):
                dofs = 6 * np.asarray(monitor["nodes"], dtype=int) + int(monitor["component"])
                source = {"displacement": u, "velocity": 1j * omega * u,
                          "acceleration": -omega * omega * u,
                          "reaction": reactions}[monitor["quantity"]]
                ref = source[:, dofs].sum(axis=1)
                if monitor["quantity"] != "reaction":
                    ref = ref / len(dofs)
                label, unit = monitor["name"], monitor["unit"]
                ours_m = (table_fr[f"{label}_re[{unit}]"].to_numpy()
                          + 1j * table_fr[f"{label}_im[{unit}]"].to_numpy())
                monitor_diff = max(monitor_diff, float(np.abs(ref - ours_m).max()
                                                       / max(np.abs(ref).max(), 1e-300)))
            entry["codes"]["numpy timoshenko harmonic"] = {
                "max_rel_diff": monitor_diff, "frequencies": len(frequencies),
                "tolerance": tol_eigen, "passed": monitor_diff <= tol_eigen,
                "element": "NumPy Timoshenko" + (" lumped" if block["mass"] == "lumped" else ""),
                "comparison": ("a direct complex solve per frequency of the NumPy K and M with "
                               "structural and Rayleigh damping; the complex monitors at every "
                               "frequency"),
            }

        if not skip_calculix:
            ccx_deck = case.path(f"calculix_{safe(name)}.inp")
            if os.path.isfile(ccx_deck):
                with tempfile.TemporaryDirectory(prefix="sparlab_ccx_beam_") as work:
                    frd = run_calculix(ccx_deck, work)
                    ref = parse_frd_displacements(frd, problem.nodes.shape[0])[:, :3]
                if np.isnan(ref).any():
                    raise RuntimeError(f"CalculiX returned no displacement for some nodes of {name}")
                diff = np.abs(ref - ours[:, :3])
                scale = max(float(np.abs(ref).max()), 1e-300)
                entry["codes"]["calculix"] = {
                    "max_abs_diff_m": float(diff.max()), "max_rel_diff": float(diff.max() / scale),
                    "ref_max_abs_m": scale, "tolerance": None, "passed": None,
                    "element": "B31 (brick expansion)", "info": "different model",
                    "comparison": ("a different model: CalculiX expands each B31 into bricks "
                                   "over its rectangle - 3-D elasticity of the section, no "
                                   "shear coefficient, Saint-Venant torsion of its own mesh; "
                                   "informational"),
                }
                if ours_lf.size and run_calculix_buckling is not None:
                    with tempfile.TemporaryDirectory(prefix="sparlab_ccx_beam_buckle_") as work:
                        ref = run_calculix_buckling(ccx_deck, work, int(ours_lf.size))
                    stats = shell_xval.compare_values(ref, ours_lf)
                    stats.update({"tolerance": None, "passed": None, "element": "B31 *BUCKLE",
                                  "info": "different model",
                                  "comparison": ("CalculiX *BUCKLE on its brick expansion of "
                                                 "B31: a different model; informational")})
                    entry["codes"]["calculix buckling"] = stats
            obstacle = u1_obstacle(problem)
            if not obstacle:
                with tempfile.TemporaryDirectory(prefix="sparlab_ccx_u1_") as work:
                    ref = run_u1(problem, f_sparlab, poisson, work)
                t = np.abs(ref[:, :3] - ours[:, :3]).max() / np.abs(ref[:, :3]).max()
                r = np.abs(ref[:, 3:] - ours[:, 3:]).max() / np.abs(ref[:, 3:]).max()
                entry["codes"]["calculix u1"] = {
                    "max_rel_diff": float(max(t, r)), "translation_rel_diff": float(t),
                    "rotation_rel_diff": float(r), "tolerance": tolerances["beam_calculix_u1"],
                    "passed": bool(max(t, r) <= tolerances["beam_calculix_u1"]),
                    "element": "U1 (k = 1e12)",
                    "comparison": ("CalculiX's U1 beam with a shear coefficient of 1e12 - the "
                                   "Euler-Bernoulli beam, torsion constant I_y + I_z - the same "
                                   "nodal loads; judged at the 7 digits of its printed output"),
                }
        report["load_cases"].append(entry)

    modal = summary.get("modal")
    if modal:
        ours_hz = np.asarray(modal.get("frequencies_hz") or [], dtype=float)
        lumped = modal.get("mass_type") == "lumped"
        if ours_hz.size:
            ref = problem.modal(int(ours_hz.size), lumped) / (2.0 * math.pi)
            stats = shell_xval.compare_values(ref, ours_hz)
            stats["tolerance"] = tol_eigen
            stats["passed"] = stats["max_rel_diff"] <= tol_eigen
            stats["element"] = "NumPy Timoshenko " + ("lumped" if lumped else "consistent") + " mass"
            stats["comparison"] = ("the frequencies of the model (listed with its first load "
                                   "case): scipy shift-invert Lanczos at sigma = 0")
            report["load_cases"][0]["codes"]["numpy timoshenko modal"] = stats
            ccx_deck = case.path(f"calculix_{safe(report['load_cases'][0]['load_case'])}.inp")
            if not skip_calculix and os.path.isfile(ccx_deck) and not lumped:
                with tempfile.TemporaryDirectory(prefix="sparlab_ccx_beam_freq_") as work:
                    ref = run_calculix_frequency(ccx_deck, work, int(ours_hz.size))
                stats = shell_xval.compare_values(ref, ours_hz)
                stats.update({"tolerance": None, "passed": None, "element": "B31 *FREQUENCY",
                              "info": "different model",
                              "comparison": ("CalculiX *FREQUENCY on its brick expansion of "
                                             "B31: a different model; informational")})
                report["load_cases"][0]["codes"]["calculix modal"] = stats
    return report
