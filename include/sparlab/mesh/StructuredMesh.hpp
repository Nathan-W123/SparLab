/// \file StructuredMesh.hpp
/// \brief Generator for structured rectangular Q4 grids.
///
/// Node numbering (nx = elements along x, ny = elements along y):
/// \code
///   node(i, j) = j * (nx + 1) + i,      i in [0, nx], j in [0, ny]
///   elem(i, j) = j * nx + i,            i in [0, nx), j in [0, ny)
/// \endcode
/// The element's nodes are stored counter-clockwise starting from the lower
/// left corner:
/// \code
///   3 --------- 2          eta
///   |           |           ^
///   |    (e)    |           |
///   |           |           +---> xi
///   0 --------- 1
/// \endcode
#pragma once

#include "sparlab/mesh/Mesh.hpp"

namespace sparlab {

/// Parameters of a structured rectangular domain.
struct StructuredMeshSpec {
  Index nx = 1;         ///< number of elements along x (>= 1)
  Index ny = 1;         ///< number of elements along y (>= 1)
  Scalar lx = 1.0;      ///< domain length along x [m] (> 0)
  Scalar ly = 1.0;      ///< domain length along y [m] (> 0)
  Scalar x0 = 0.0;      ///< x coordinate of the lower-left corner [m]
  Scalar y0 = 0.0;      ///< y coordinate of the lower-left corner [m]
};

/// Build a uniform structured Q4 mesh. The resulting mesh is validated before
/// being returned and carries StructuredGridInfo so downstream components can
/// exploit the identical-cell property.
/// \throws ConfigError for non-positive counts or lengths.
Mesh make_structured_quad_mesh(const StructuredMeshSpec& spec);

/// Build a structured Q4 mesh whose *interior* nodes are randomly displaced by
/// up to `perturbation` times the local cell size. Boundary nodes are left in
/// place so the domain shape is preserved.
///
/// Distorted meshes are what make a patch test meaningful: on a uniform grid a
/// Q4 element reproduces linear fields for trivial reasons, whereas on a
/// distorted grid only a correctly formed isoparametric mapping does. The
/// generator is deterministic for a given seed.
/// \param perturbation fraction of the cell size, in [0, 0.45).
/// \throws ConfigError for an out-of-range perturbation.
Mesh make_perturbed_quad_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                              unsigned int seed = 12345u);

/// Convenience accessors for structured grids (used by tests and selectors).
/// \{
Index structured_node_index(const StructuredGridInfo& info, Index i, Index j);
Index structured_element_index(const StructuredGridInfo& info, Index i, Index j);
/// \}

}  // namespace sparlab
