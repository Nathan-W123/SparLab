"""Isosurfaces of a density field on a SparLab structured_hex box mesh.

contour(element_field, level, up=2, smooth=2, taubin=40, min_area=0)
    -> (vertices [m], triangles)

The element field (engine order, x fastest) is averaged to the nodes,
smoothed by a few Laplacian passes, padded with a void layer on every face so
the surface closes there, contoured by marching cubes (upsampled `up` times),
and mapped to metres; the void layer maps onto the box faces themselves, so
material never shows beyond the domain. Taubin smoothing then evens out the
cell-scale steps, with vertices on the box faces kept on them.
"""
import numpy as np
from scipy.ndimage import map_coordinates, uniform_filter
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from skimage.measure import marching_cubes


class BoxGrid:
    def __init__(self, nx, ny, nz, lx, ly, lz):
        self.n = np.array([nx, ny, nz])
        self.L = np.array([lx, ly, lz], dtype=float)
        self.h = self.L / self.n

    def centroids(self):
        nx, ny, nz = self.n
        k, j, i = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
        return np.column_stack([(i.ravel() + 0.5) * self.h[0], (j.ravel() + 0.5) * self.h[1],
                                (k.ravel() + 0.5) * self.h[2]])

    def nodal(self, element_field, smooth=2):
        nx, ny, nz = self.n
        E = np.asarray(element_field, dtype=float).reshape(nz, ny, nx)
        P = np.pad(E, 1, mode="constant")
        C = np.pad(np.ones_like(E), 1, mode="constant")
        s = sum(P[a:a + nz + 1, b:b + ny + 1, c:c + nx + 1] for a in (0, 1) for b in (0, 1) for c in (0, 1))
        w = sum(C[a:a + nz + 1, b:b + ny + 1, c:c + nx + 1] for a in (0, 1) for b in (0, 1) for c in (0, 1))
        V = s / w
        for _ in range(smooth):
            V = 0.5 * V + 0.5 * uniform_filter(V, size=3, mode="nearest")
        return V                                   # (nz+1, ny+1, nx+1)

    def contour(self, element_field, level, up=2, smooth=2, taubin=40, min_area=0.0):
        return self.contour_nodal(self.nodal(element_field, smooth), level, up, taubin, min_area)

    def contour_nodal(self, nodal_field, level, up=2, taubin=40, min_area=0.0):
        """As contour(), for a field already on the nodes, (nz+1, ny+1, nx+1)."""
        V = np.pad(nodal_field, 1, mode="constant", constant_values=-1.0e3)
        if up > 1:
            shape = [(d - 1) * up + 1 for d in V.shape]
            grid = np.meshgrid(*[np.linspace(0, d - 1, s) for d, s in zip(V.shape, shape)], indexing="ij")
            V = map_coordinates(V, grid, order=1)
        v, F, _, _ = marching_cubes(V, level)
        v = v / up - 1.0                            # node index coordinates (z, y, x)
        nzyx = self.n[::-1]
        on_face = (v <= 0.0) | (v >= nzyx)
        v = np.clip(v, 0.0, nzyx)
        X = (v * self.h[::-1])[:, ::-1].astype(float)          # metres, (x, y, z)
        face = on_face[:, ::-1]
        # marching_cubes already winds this volume outward
        if min_area > 0:
            e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]]])
            adj = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(len(X), len(X)))
            _, label = connected_components(adj, directed=False)
            area = 0.5 * np.linalg.norm(np.cross(X[F[:, 1]] - X[F[:, 0]], X[F[:, 2]] - X[F[:, 0]]), axis=1)
            piece = np.bincount(label[F[:, 0]], weights=area, minlength=label.max() + 1)
            F = F[piece[label[F[:, 0]]] >= min_area]
            used = np.unique(F)
            remap = -np.ones(len(X), dtype=np.int64)
            remap[used] = np.arange(len(used))
            X, F, face = X[used], remap[F], face[used]
        if taubin > 0 and len(F):
            e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
            A = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(len(X), len(X))).tocsr()
            A = ((A + A.T) > 0).astype(float)
            W = diags(1.0 / np.maximum(A.sum(axis=1).A1, 1.0)) @ A
            lo, hi = X.min(axis=0), X.max(axis=0)
            target = np.where(X <= 1e-12, 0.0, self.L)            # the face each pinned coordinate is on
            for _ in range(taubin):
                for fct in (0.5, -0.53):
                    X = X + fct * (W @ X - X)
                    X[face] = target[face]          # keep face vertices on their faces
        return X, F
