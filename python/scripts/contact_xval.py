"""Cross-validation of SparLab's contact, used by cross_validate.py for a run
with a `contact` block.

The same discrete problem is solved again, independently, by scikit-fem:

  * the stiffness from scikit-fem (SkfemProblem), the loads and prescribed
    displacements of mesh.json, scaled by SparLab's load factors;
  * the contact discretisation computed here from the pairs' faces
    (mesh.json `contact_surfaces`) and definitions (summary.json): the weights
    D_j = int N_j dA of the slave nodes, the dual basis of every slave face
    (A = D M^-1 per face), the mortar integrals M_jl and the weighted initial
    gaps - by exact quadrature on the overlaps of the slave faces and the
    master faces projected onto them (segments in 2-D, clipped polygons in
    3-D) - and, against a rigid obstacle, its distance and normal at the
    slave nodes. The nodes SparLab leaves out are left out by the same
    rules (coverage below 0.99 or above 1.01, every component prescribed,
    prescribed along the normal);
  * the solution by a semismooth Newton method on the uncondensed problem -
    the displacements and, per slave node, the contact pressure and the
    tangential traction as unknowns - with the complementarity functions of
    Alart and Curnier (the tangential traction the projection of
    tau - c u_slip / D onto the disc of radius mu p), where SparLab condenses
    the pressure out and switches rows. Coulomb friction acts on the slip of
    each step, measured from its start, so the steps are SparLab's.

The mortar integrals are exact for a flat slave surface (the projection
along its constant normal is affine); a curved slave surface - whose normal
field SparLab interpolates - is outside this comparison and is reported as
such. Compared: the displacement field and, node by node, the contact
status, the pressure and the tangential traction, at the end of the run.
"""

from __future__ import annotations

import json
import math
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# SparLab's rules (Contact.cpp): coverage bounds of a mortar row, the
# touching tolerance of the first iteration.
MIN_COVERAGE = 0.99
MAX_COVERAGE = 1.01
TOUCHING = 1.0e-9


# ---------------------------------------------------------------------------
# Quadrature
# ---------------------------------------------------------------------------
def gauss_segment(n: int) -> Tuple[np.ndarray, np.ndarray]:
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def gauss_triangle(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Collapsed Gauss rule on the unit triangle (points (xi, eta), weights
    summing to 1/2): exact for polynomials of degree 2n - 2."""
    s, ws = gauss_segment(n)
    pts, wts = [], []
    for a, wa in zip(s, ws):
        for b, wb in zip(s, ws):
            pts.append((a * (1.0 - b), b))
            wts.append(wa * wb * (1.0 - b))
    return np.array(pts), np.array(wts)


def shape(corners: int, xi: float, eta: float = 0.0) -> np.ndarray:
    """Linear shape functions of a segment, triangle or quadrilateral face in
    the order of its corners (the quadrilateral on [-1, 1]^2, counter-
    clockwise from (-1, -1))."""
    if corners == 2:
        return np.array([1.0 - xi, xi])
    if corners == 3:
        return np.array([1.0 - xi - eta, xi, eta])
    return 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                            (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])


def shape_derivatives(corners: int, xi: float, eta: float) -> np.ndarray:
    if corners == 3:
        return np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
    return 0.25 * np.array([[-(1 - eta), -(1 - xi)], [(1 - eta), -(1 + xi)],
                            [(1 + eta), (1 + xi)], [-(1 + eta), (1 - xi)]])


def face_rule(corners: int) -> Tuple[np.ndarray, np.ndarray]:
    """Points and weights on the reference face: exact for the products of
    two shape functions and the area element of a flat face."""
    if corners == 2:
        s, w = gauss_segment(4)
        return s[:, None], w
    if corners == 3:
        return gauss_triangle(5)
    g, w = np.polynomial.legendre.leggauss(4)
    return (np.array([(a, b) for a in g for b in g]),
            np.array([wa * wb for wa in w for wb in w]))


# ---------------------------------------------------------------------------
# Slave faces: weights, dual basis, normals
# ---------------------------------------------------------------------------
class Face:
    """A boundary face (edge in 2-D) with its element's outward normal."""

    def __init__(self, nodes: Sequence[int], coords: np.ndarray, centroid: np.ndarray,
                 thickness: float):
        self.nodes = list(nodes)
        self.x = np.asarray(coords, dtype=float)       # (corners, 3)
        self.corners = len(self.nodes)
        self.thickness = thickness
        if self.corners == 2:
            t = self.x[1] - self.x[0]
            n = np.array([t[1], -t[0], 0.0])
        elif self.corners == 3:
            n = np.cross(self.x[1] - self.x[0], self.x[2] - self.x[0])
        else:
            n = np.cross(self.x[2] - self.x[0], self.x[3] - self.x[1])
        n /= np.linalg.norm(n)
        if np.dot(n, self.x.mean(axis=0) - centroid) < 0.0:
            n = -n
        self.normal = n                                  # outward, unit

    def point(self, xi: float, eta: float = 0.0) -> np.ndarray:
        return shape(self.corners, xi, eta) @ self.x

    def area_element(self, xi: float, eta: float = 0.0) -> float:
        if self.corners == 2:
            return float(np.linalg.norm(self.x[1] - self.x[0])) * self.thickness
        d = shape_derivatives(self.corners, xi, eta)
        a = d[:, 0] @ self.x
        b = d[:, 1] @ self.x
        return float(np.linalg.norm(np.cross(a, b)))

    def weights_and_dual(self) -> Tuple[np.ndarray, np.ndarray]:
        """D_a = int N_a dA and the dual-basis coefficients A (psi = A N)."""
        pts, wts = face_rule(self.corners)
        d = np.zeros(self.corners)
        m = np.zeros((self.corners, self.corners))
        for p, w in zip(pts, wts):
            xi, eta = (p[0], 0.0) if self.corners == 2 else (p[0], p[1])
            n = shape(self.corners, xi, eta)
            da = w * self.area_element(xi, eta)
            d += n * da
            m += np.outer(n, n) * da
        return d, np.diag(d) @ np.linalg.inv(m)


def plane_frame(normal: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    k = int(np.argmin(np.abs(normal)))
    e1 = np.zeros(3)
    e1[k] = 1.0
    e1 -= e1.dot(normal) * normal
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(normal, e1)


def clip(subject: List[np.ndarray], clipper: List[np.ndarray]) -> List[np.ndarray]:
    """Sutherland-Hodgman: `subject` clipped by the convex, counter-clockwise
    polygon `clipper` (2-D points)."""
    out = subject
    for i in range(len(clipper)):
        a, b = clipper[i], clipper[(i + 1) % len(clipper)]
        edge = b - a
        inside = lambda p: edge[0] * (p[1] - a[1]) - edge[1] * (p[0] - a[0]) >= -1e-14 * (
            np.dot(edge, edge))
        points, out = out, []
        if not points:
            break
        for j in range(len(points)):
            p, q = points[j], points[(j + 1) % len(points)]
            pin, qin = inside(p), inside(q)
            if pin:
                out.append(p)
            if pin != qin:
                d = q - p
                den = edge[0] * d[1] - edge[1] * d[0]
                if den != 0.0:
                    t = (edge[0] * (a[1] - p[1]) - edge[1] * (a[0] - p[0])) / den
                    out.append(p + min(max(t, 0.0), 1.0) * d)
    return out


def ccw(points: List[np.ndarray]) -> List[np.ndarray]:
    area = 0.0
    for i in range(len(points)):
        p, q = points[i], points[(i + 1) % len(points)]
        area += p[0] * q[1] - q[0] * p[1]
    return points if area > 0.0 else points[::-1]


def invert_map(corners: int, x2: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Local coordinates of `target` on a face given by its 2-D corner
    coordinates `x2`: exact for a triangle, Newton for a quadrilateral."""
    if corners == 3:
        m = np.column_stack([x2[1] - x2[0], x2[2] - x2[0]])
        return np.linalg.solve(m, target - x2[0])
    xi = np.zeros(2)
    for _ in range(50):
        r = shape(4, xi[0], xi[1]) @ x2 - target
        j = (shape_derivatives(4, xi[0], xi[1]).T @ x2).T
        step = np.linalg.solve(j, r)
        xi -= step
        if np.abs(step).max() < 1e-15:
            break
    return xi


# ---------------------------------------------------------------------------
# The contact discretisation
# ---------------------------------------------------------------------------
class ContactNode:
    def __init__(self):
        self.node = 0
        self.pair = 0
        self.weight = 0.0
        self.normal = np.zeros(3)
        self.initial_gap = 0.0                 # weighted
        self.masters: List[Tuple[int, float]] = []
        self.free: List[int] = []
        self.tangents: List[np.ndarray] = []
        self.c = 0.0
        self.size = 0.0
        self.friction = 0.0
        self.motion = np.zeros(3)


def obstacle_gap(ob: Dict, x: np.ndarray, dim: int) -> Tuple[float, np.ndarray]:
    kind = ob["type"]
    p = np.asarray(ob["point_m"], dtype=float)
    if kind == "plane":
        n = np.asarray(ob["normal"], dtype=float)
        n /= np.linalg.norm(n)
        return float(n.dot(x - p)), n
    sign = -1.0 if ob.get("inside", False) else 1.0
    r = x - p
    if kind == "cylinder":
        axis = np.array([0.0, 0.0, 1.0]) if dim == 2 else np.asarray(ob["axis"], dtype=float)
        axis = axis / np.linalg.norm(axis)
        r = r - r.dot(axis) * axis
    dist = float(np.linalg.norm(r))
    return sign * (dist - float(ob["radius_m"])), sign * r / dist


def build_contact(mesh, summary: Dict, surfaces: List[Dict], prescribed: set,
                  element_youngs: np.ndarray, thickness: float) -> Tuple[List[ContactNode], List[str]]:
    """The slave nodes that take part, and the reasons for any part of the
    comparison that cannot be made. `prescribed` holds (node, component)."""
    dim = mesh.dim
    nodes3 = np.zeros((mesh.num_nodes, 3))
    nodes3[:, :dim] = mesh.nodes
    centroids = nodes3[mesh.elements].mean(axis=1)
    block = summary["nonlinear"]["contact"]
    complementarity = float(block.get("complementarity", 1.0))
    search_factor = float(block.get("search_factor", 2.0))
    by_name = {s["name"]: s for s in surfaces}
    out: List[ContactNode] = []
    problems: List[str] = []

    def faces_of(entries: List[Dict]) -> List[Face]:
        faces = []
        for f in entries:
            ids = [int(n) for n in f["nodes"]]
            faces.append(Face(ids, nodes3[ids], centroids[int(f["element"])],
                              thickness if dim == 2 else 1.0))
        return faces

    for k, pair in enumerate(block["pairs"]):
        rigid = pair["kind"] == "rigid obstacle"
        slave = faces_of(by_name[pair["name"]]["slave_faces"])
        master = [] if rigid else faces_of(by_name[pair["name"]]["master_faces"])
        weight: Dict[int, float] = {}
        normal_sum: Dict[int, np.ndarray] = {}
        modulus: Dict[int, float] = {}
        rows: Dict[int, Dict[int, float]] = {}
        gaps: Dict[int, float] = {}
        flat = True
        for f in slave:
            flat = flat and abs(abs(f.normal.dot(slave[0].normal)) - 1.0) < 1e-12
        if not rigid and not flat:
            problems.append(f"contact pair '{pair['name']}': the slave surface is curved; the "
                            "mortar integrals here assume a flat one")
            continue
        element_of_face = [int(f["element"]) for f in by_name[pair["name"]]["slave_faces"]]
        for f, e in zip(slave, element_of_face):
            d, a = f.weights_and_dual()
            for i, n in enumerate(f.nodes):
                weight[n] = weight.get(n, 0.0) + d[i]
                normal_sum[n] = normal_sum.get(n, np.zeros(3)) + f.normal
                modulus[n] = max(modulus.get(n, 0.0), float(element_youngs[e]))
            if rigid:
                continue
            # Mortar integrals over the overlaps with the master faces within
            # the search distance (SparLab's: the slave face's bounding box
            # grown by search_factor times its diagonal).
            ns = f.normal
            lo, hi = f.x.min(axis=0), f.x.max(axis=0)
            reach = search_factor * float(np.linalg.norm(hi - lo))
            near = [m for m in master
                    if not (np.any(m.x.min(axis=0) > hi + reach)
                            or np.any(m.x.max(axis=0) < lo - reach))]
            if dim == 2:
                e_t = f.x[1] - f.x[0]
                length = float(np.linalg.norm(e_t))
                e_t /= length
                s, w = gauss_segment(5)
                for m in near:
                    sc = (m.x[0] - f.x[0]).dot(e_t) / length
                    sd = (m.x[1] - f.x[0]).dot(e_t) / length
                    lo, hi = max(0.0, min(sc, sd)), min(1.0, max(sc, sd))
                    if hi - lo <= 1e-14:
                        continue
                    for si, wi in zip(lo + (hi - lo) * s, (hi - lo) * w):
                        xs = f.point(si)
                        r = (si - sc) / (sd - sc)
                        xm = (1.0 - r) * m.x[0] + r * m.x[1]
                        psi = a @ shape(2, si)
                        nm = shape(2, r)
                        g0 = float((xm - xs).dot(ns))
                        da = wi * length * f.thickness
                        for i, j in enumerate(f.nodes):
                            row = rows.setdefault(j, {})
                            for q, l in enumerate(m.nodes):
                                row[l] = row.get(l, 0.0) + psi[i] * nm[q] * da
                            gaps[j] = gaps.get(j, 0.0) + psi[i] * g0 * da
            else:
                e1, e2 = plane_frame(ns)
                origin = f.x[0]
                to2 = lambda x: np.array([(x - origin).dot(e1), (x - origin).dot(e2)])
                s2 = ccw([to2(x) for x in f.x])
                s2_order = [to2(x) for x in f.x]
                tri_pts, tri_wts = gauss_triangle(6)
                for m in near:
                    m2 = [to2(x) for x in m.x]
                    poly = clip(ccw(list(m2)), s2)
                    if len(poly) < 3:
                        continue
                    for t in range(1, len(poly) - 1):
                        p0, p1, p2 = poly[0], poly[t], poly[t + 1]
                        area2 = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])
                        if abs(area2) <= 1e-30:
                            continue
                        for (u, v), wq in zip(tri_pts, tri_wts):
                            q2 = p0 + u * (p1 - p0) + v * (p2 - p0)
                            da = wq * abs(area2)
                            xi_s = invert_map(f.corners, np.array(s2_order), q2)
                            xi_m = invert_map(m.corners, np.array(m2), q2)
                            psi = a @ shape(f.corners, xi_s[0], xi_s[1])
                            nm = shape(m.corners, xi_m[0], xi_m[1])
                            xs = f.point(xi_s[0], xi_s[1])
                            xm = m.point(xi_m[0], xi_m[1])
                            g0 = float((xm - xs).dot(ns))
                            for i, j in enumerate(f.nodes):
                                row = rows.setdefault(j, {})
                                for q, l in enumerate(m.nodes):
                                    row[l] = row.get(l, 0.0) + psi[i] * nm[q] * da
                                gaps[j] = gaps.get(j, 0.0) + psi[i] * g0 * da
        for j in sorted(weight):
            c = ContactNode()
            c.node, c.pair, c.weight = j, k, weight[j]
            c.friction = float(pair.get("friction", 0.0))
            x = nodes3[j]
            if rigid:
                g, nu = obstacle_gap(pair["obstacle"], x, dim)
                c.normal = nu
                c.initial_gap = c.weight * g
                c.motion = np.asarray(pair["obstacle"]["motion_m"], dtype=float)
                if dim == 2:
                    c.motion[2] = 0.0
            else:
                row = rows.get(j)
                covered = sum(row.values()) if row else 0.0
                ratio = covered / c.weight
                if row is None or ratio < MIN_COVERAGE or ratio > MAX_COVERAGE:
                    continue
                factor = c.weight / covered
                c.masters = [(l, v * factor) for l, v in row.items() if abs(v) > 1e-14 * c.weight]
                c.initial_gap = gaps.get(j, 0.0) * factor
                ns = normal_sum[j] / np.linalg.norm(normal_sum[j])
                c.normal = -ns
            c.free = [q for q in range(dim) if (j, q) not in prescribed]
            if not c.free:
                continue
            fn = np.zeros(3)
            for q in c.free:
                fn[q] = c.normal[q]
            if np.linalg.norm(fn) < 1e-6:
                continue
            fnu = fn / np.linalg.norm(fn)
            axes = sorted((abs(fnu[q]), q) for q in c.free)[:-1]
            for _, q in axes:
                t = np.zeros(3)
                t[q] = 1.0
                t -= t.dot(fnu) * fnu
                for prev in c.tangents:
                    t -= t.dot(prev) * prev
                if np.linalg.norm(t) > 1e-10:
                    c.tangents.append(t / np.linalg.norm(t))
            c.size = c.weight / thickness if dim == 2 else math.sqrt(c.weight)
            c.c = complementarity * modulus[j] / c.size
            out.append(c)
    return out, problems


# ---------------------------------------------------------------------------
# The semismooth Newton method (uncondensed, Alart-Curnier)
# ---------------------------------------------------------------------------
def sparlab_rule_stiffness(problem) -> sp.csr_matrix:
    """The stiffness integrated with SparLab's rule - 2 x 2 (x 2) points for
    the Q4 and Hex8, as the non-linear comparisons do - where scikit-fem's
    default takes 3 x 3 (4 x 4 x 4): the two agree on parallelogram cells
    only, and a curved cap's cells are not."""
    from skfem import Basis, BilinearForm, asm
    from skfem.helpers import ddot, eye, sym_grad, trace

    intorder = {"Quad4": 3, "Hex8": 3, "Tri3": 1, "Tet4": 1, "Tet10": 2}[problem.mesh.element_type]
    basis = Basis(problem.skfem_mesh, problem.vector_element, intorder=intorder)
    if not np.array_equal(basis.nodal_dofs, problem.basis.nodal_dofs):
        raise RuntimeError("the scikit-fem DOF numbering changed with the quadrature")
    nqp = basis.X.shape[-1]
    lam = np.repeat(np.asarray(problem.lam_e)[:, None], nqp, axis=1)
    mu = np.repeat(np.asarray(problem.mu_e)[:, None], nqp, axis=1)
    dim = problem.dim

    @BilinearForm
    def elasticity(u, v, w):
        eps = sym_grad(u)
        return ddot(w["lam"] * eye(trace(eps), dim) + 2.0 * w["mu"] * eps, sym_grad(v))

    return (asm(elasticity, basis, lam=lam, mu=mu) * problem.thickness).tocsr()


class ContactSolver:
    def __init__(self, problem, contact: List[ContactNode], dim: int):
        self.p = problem
        self.nodes = contact
        self.dim = dim
        free = problem.free
        self.position = -np.ones(problem.basis.N, dtype=int)
        self.position[free] = np.arange(free.size)
        k = sparlab_rule_stiffness(problem)
        self.kff = k[free][:, free].tocsc()
        self.kfp = k[free][:, problem.fixed]
        # The kinematic rows of every node: normal, then tangents, over the
        # free DOFs (C), and over all DOFs (for the values).
        rows_c, cols_c, vals_c = [], [], []
        self.offsets = []
        count = 0
        self.dof = problem.dof_of_node
        for c in contact:
            directions = [c.normal] + (c.tangents if c.friction > 0.0 else [])
            self.offsets.append(count)
            for d_i, w in enumerate(directions):
                for q in range(dim):
                    if w[q] == 0.0:
                        continue
                    col = self.position[self.dof[c.node, q]]
                    if col >= 0:
                        rows_c.append(count + d_i)
                        cols_c.append(col)
                        vals_c.append(c.weight * w[q])
                    for l, m in c.masters:
                        col = self.position[self.dof[l, q]]
                        if col >= 0:
                            rows_c.append(count + d_i)
                            cols_c.append(col)
                            vals_c.append(-m * w[q])
            count += len(directions)
        self.nz = count
        self.c_matrix = sp.csr_matrix((vals_c, (rows_c, cols_c)), shape=(count, free.size))

    def kinematic(self, c: ContactNode, u: np.ndarray) -> np.ndarray:
        """w = D u_j - sum_l M_jl u_l (full-length u in scikit-fem numbering)."""
        w = np.zeros(3)
        for q in range(self.dim):
            w[q] = c.weight * u[self.dof[c.node, q]]
            for l, m in c.masters:
                w[q] -= m * u[self.dof[l, q]]
        return w

    def solve(self, factors: Sequence[float], tolerance: float = 1e-12,
              max_iterations: int = 100) -> Dict:
        pr = self.p
        u = pr.x * 0.0
        z = np.zeros(self.nz)
        lam0 = 0.0
        u0 = u.copy()
        history = []
        for step, lam in enumerate(factors):
            u = u.copy()
            u[pr.fixed] = lam * pr.x[pr.fixed]
            status_prev = None
            converged = False
            for it in range(1, max_iterations + 1):
                first = step == 0 and it == 1
                f_vec, jac, status, scale = self.system(u, z, lam, u0, lam0, first)
                with np.errstate(over="ignore", divide="ignore"):
                    # (Before any contact force the pressures' scale is 0 and
                    # an unconverged iterate's measure infinite.)
                    norm = np.linalg.norm(f_vec[: pr.free.size]) / scale[0]
                    norm_c = (np.linalg.norm(f_vec[pr.free.size:]) / scale[1]) if self.nz else 0.0
                if status_prev is not None and status == status_prev and max(norm, norm_c) <= tolerance:
                    converged = True
                    break
                status_prev = status
                dx = spla.spsolve(jac.tocsc(), -f_vec)
                if not np.all(np.isfinite(dx)):
                    raise RuntimeError("the scikit-fem contact Newton step is not finite")
                u[pr.free] += dx[: pr.free.size]
                z += dx[pr.free.size:]
            if not converged:
                raise RuntimeError(f"the scikit-fem contact Newton method did not converge at "
                                   f"lambda = {lam} in {max_iterations} iterations")
            history.append({"load_factor": lam, "iterations": it})
            start, lam_start = u0, lam0
            u0, lam0 = u.copy(), lam
        return {"u": u, "z": z, "status": status, "history": history, "u_start": start,
                "load_factor": lam0, "load_factor_start": lam_start}

    def system(self, u, z, lam, u0, lam0, first):
        pr = self.p
        nf = pr.free.size
        r = self.kff @ u[pr.free] + self.kfp @ u[pr.fixed] - lam * pr.f[pr.free]
        contact_force = self.c_matrix.T @ z if self.nz else np.zeros(nf)
        f_u = r - contact_force
        rows, cols, vals = [], [], []
        f_c = np.zeros(self.nz)
        status = []
        pmax = 0.0
        for idx, c in enumerate(self.nodes):
            o = self.offsets[idx]
            m = len(c.tangents) if c.friction > 0.0 else 0
            d = c.weight
            w = self.kinematic(c, u)
            gap = c.initial_gap + c.normal.dot(w) - d * lam * c.normal.dot(c.motion)
            p = z[o]
            pmax = max(pmax, abs(p))
            trial = p - c.c * gap / d
            touching = first and p >= 0.0 and abs(gap / d) <= TOUCHING * c.size
            active = trial > 0.0 or touching
            cn = self.c_matrix[o]
            if active:
                f_c[o] = c.c * gap / d
                rows.extend([o] * cn.nnz)
                cols.extend(cn.indices.tolist())
                vals.extend((c.c / d * cn.data).tolist())
            else:
                f_c[o] = p
                rows.append(o)
                cols.append(nf + o)
                vals.append(1.0)
            if m == 0:
                # Frictionless, or with no free tangential direction (a
                # frictional node that cannot slip sticks, as SparLab reports).
                status.append(("stick" if c.friction > 0.0 else "slip") if active else "open")
                continue
            # The weighted slip of the step and the Coulomb disc.
            w0 = self.kinematic(c, u0)
            rel = w - w0 - d * (lam - lam0) * c.motion
            slip = np.array([t.dot(rel) for t in c.tangents])
            tau = z[o + 1: o + 1 + m]
            q = tau - c.c * slip / d
            p_hat = max(trial, 0.0) if active else 0.0
            radius = c.friction * p_hat
            qn = float(np.linalg.norm(q))
            ct = self.c_matrix[o + 1: o + 1 + m]
            if qn <= radius:
                f_c[o + 1: o + 1 + m] = tau - q
                coo = ct.tocoo()
                rows.extend((o + 1 + coo.row).tolist())
                cols.extend(coo.col.tolist())
                vals.extend((c.c / d * coo.data).tolist())
                status.append("stick" if active else "open")
            else:
                qh = q / qn
                proj = np.eye(m) - np.outer(qh, qh)
                f_c[o + 1: o + 1 + m] = tau - radius * qh
                # d/d tau
                dt = np.eye(m) - radius * proj / qn
                for a_ in range(m):
                    for b_ in range(m):
                        if dt[a_, b_] != 0.0:
                            rows.append(o + 1 + a_)
                            cols.append(nf + o + 1 + b_)
                            vals.append(dt[a_, b_])
                # d/du: -radius proj/|q| dq/du - qh dradius/du
                du = (radius / qn) * proj @ (c.c / d * ct.toarray())
                if p_hat > 0.0:
                    du += np.outer(qh, c.friction * c.c / d * cn.toarray().ravel())
                    for a_ in range(m):
                        rows.append(o + 1 + a_)
                        cols.append(nf + o)
                        vals.append(-qh[a_] * c.friction)
                nzr, nzc = np.nonzero(du)
                rows.extend((o + 1 + nzr).tolist())
                cols.extend(nzc.tolist())
                vals.extend(du[nzr, nzc].tolist())
                status.append("slip" if active else "open")
        n = nf + self.nz
        top = sp.hstack([self.kff, -self.c_matrix.T]) if self.nz else self.kff
        bottom = sp.csr_matrix((vals, (rows, cols)), shape=(self.nz, n))
        jac = sp.vstack([top, bottom]).tocsr()
        scale_u = max(np.linalg.norm(lam * pr.f[pr.free]), np.linalg.norm(self.kfp @ u[pr.fixed]),
                      np.linalg.norm(contact_force), 1e-300)
        return np.concatenate([f_u, f_c]), jac, status, (scale_u, max(pmax, 1e-300))


# ---------------------------------------------------------------------------
# Comparison with SparLab's run
# ---------------------------------------------------------------------------
def contact_comparisons(case, summary: Dict, name: str, nl_case: Dict, problem, tolerance: float,
                        load_factors: Sequence[float], nonlinear_displacement: np.ndarray) -> Dict:
    """The scikit-fem solution of load case `name` compared with SparLab's."""
    mesh = case.mesh
    dim = mesh.dim
    with open(case.path("mesh.json"), "r", encoding="utf-8") as handle:
        surfaces = json.load(handle).get("contact_surfaces")
    if not surfaces:
        raise ValueError("mesh.json holds no contact_surfaces; rerun sparlab_solve")
    prescribed = set()
    for entry in mesh.prescribed:
        prescribed.add((int(entry["node"]), "xyz".index(entry["component"])))
    materials = summary.get("materials") or [summary["material"]]
    youngs = np.array([float(m["youngs_modulus_Pa"]) for m in materials])
    element_youngs = youngs[mesh.element_materials] if mesh.element_materials is not None \
        else np.full(mesh.num_elements, youngs[0])
    thickness = float(summary["mesh"].get("thickness_m", 1.0))
    contact, problems = build_contact(mesh, summary, surfaces, prescribed, element_youngs,
                                      thickness)
    if problems:
        return {"codes": {}, "notes": problems}
    solver = ContactSolver(problem, contact, dim)
    result = solver.solve(load_factors)
    ours = nonlinear_displacement
    theirs = problem.nodal(result["u"])
    diff = np.abs(theirs - ours)
    ref = float(np.abs(ours).max())
    stats = {
        "max_abs_diff_m": float(diff.max()),
        "max_rel_diff": float(diff.max() / max(ref, 1e-300)),
        "rms_rel_diff": float(np.sqrt((diff ** 2).mean()) / max(ref, 1e-300)),
        "ref_max_abs_m": float(np.abs(theirs).max()),
        "sparlab_max_abs_m": ref,
        "load_factors": len(load_factors),
        "newton_iterations": [h["iterations"] for h in result["history"]],
    }
    # Node by node: status, pressure, tangential traction.
    table = case.table(f"contact_{_safe(name)}.csv")
    by_node = {int(r["node"]): r for _, r in table.iterrows()}
    pressures, sp_pressures, mismatches, boundary = [], [], [], []
    trac_diff, trac_ref = 0.0, 0.0
    p_scale = max(float(np.abs(result["z"]).max()) if result["z"].size else 0.0, 1e-300)
    u_scale = max(ref, 1e-300)
    u_last = result["u"]
    u_start = result["u_start"]
    for idx, c in enumerate(contact):
        o = solver.offsets[idx]
        row = by_node.get(c.node)
        if row is None:
            mismatches.append(f"node {c.node}: taking part here, not in SparLab's run")
            continue
        p = float(result["z"][o])
        status = result["status"][idx]
        if status != row["status"]:
            # A node on the edge of its status - touching with no pressure, or
            # at the Coulomb limit without slip - is classed by round-off
            # alone: both classes describe its state.
            w = solver.kinematic(c, u_last)
            gap = (c.initial_gap + c.normal.dot(w)
                   - c.weight * result["load_factor"] * c.normal.dot(c.motion)) / c.weight
            on_edge = False
            if "open" in (status, row["status"]):
                on_edge = abs(p) <= 1e-9 * p_scale and abs(gap) <= 1e-9 * u_scale
            elif c.friction > 0.0:
                m = len(c.tangents)
                tau = np.linalg.norm(result["z"][o + 1: o + 1 + m])
                rel = w - solver.kinematic(c, u_start) - c.weight * (
                    result["load_factor"] - result["load_factor_start"]) * c.motion
                slip = np.linalg.norm([t.dot(rel) for t in c.tangents]) / c.weight
                on_edge = (abs(tau - c.friction * p) <= 1e-9 * p_scale
                           and slip <= 1e-9 * u_scale)
            text = f"node {c.node}: {status} here, {row['status']} in SparLab"
            (boundary if on_edge else mismatches).append(text)
        pressures.append(p)
        sp_pressures.append(float(row["pressure[Pa]"]))
        if c.friction > 0.0:
            m = len(c.tangents)
            tau = sum(t * v for t, v in zip(c.tangents, result["z"][o + 1: o + 1 + m]))
            theirs_t = np.array([float(row[f"t{a}[Pa]"]) for a in "xyz"[:dim]] + [0.0] * (3 - dim))
            trac_diff = max(trac_diff, float(np.abs(tau - theirs_t).max()))
            trac_ref = max(trac_ref, float(np.abs(theirs_t).max()))
    if len(by_node) != len(contact):
        mismatches.append(f"{len(contact)} slave nodes take part here, {len(by_node)} in SparLab")
    pmax = max(max(np.abs(sp_pressures)) if sp_pressures else 0.0, 1e-300)
    stats["pressure_max_rel_diff"] = float(np.abs(np.array(pressures) - np.array(sp_pressures)).max()
                                           / pmax) if pressures else 0.0
    if trac_ref > 0.0:
        stats["traction_max_rel_diff"] = trac_diff / max(pmax, 1e-300)
    stats["status_mismatches"] = mismatches
    stats["status_on_edge"] = boundary
    counts = {}
    for s in result["status"]:
        counts[s] = counts.get(s, 0) + 1
    stats["nodes"] = counts
    stats["tolerance"] = tolerance
    stats["passed"] = (stats["max_rel_diff"] <= tolerance
                       and stats["pressure_max_rel_diff"] <= tolerance
                       and stats.get("traction_max_rel_diff", 0.0) <= tolerance
                       and not mismatches)
    stats["comparison"] = (
        "the same discrete problem solved independently: scikit-fem's stiffness, the dual "
        "mortar weights, integrals and gaps computed here from the pairs' faces, and a "
        "semismooth Newton method on the uncondensed Alart-Curnier complementarity functions, "
        "through SparLab's load factors")
    return {"codes": {"scikit-fem contact": stats}, "notes": []}


def _safe(name: str) -> str:
    return "".join(ch if (ch.isalnum() or ch in "-_") else "_" for ch in name) or "unnamed"
