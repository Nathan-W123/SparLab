/// \file StructuredMesh.hpp
/// \brief Generators for structured Q4 grids and structured Hex8 grids.
///
/// 2-D node numbering (nx = elements along x, ny = elements along y):
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
/// 3-D numbering adds a slowest-varying z index:
/// \code
///   node(i, j, k) = k * (nx + 1) * (ny + 1) + j * (nx + 1) + i
///   elem(i, j, k) = k * nx * ny + j * nx + i
/// \endcode
/// and each hexahedron lists its bottom face (z = z_k) counter-clockwise seen
/// from +z followed by the top face in the same order, which is the VTK
/// hexahedron convention.
#pragma once

#include "sparlab/mesh/Mesh.hpp"

namespace sparlab {

/// Parameters of a structured box domain. `nz`/`lz` are ignored by the 2-D
/// generator and required by the 3-D one.
struct StructuredMeshSpec {
  Index nx = 1;         ///< number of elements along x (>= 1)
  Index ny = 1;         ///< number of elements along y (>= 1)
  Index nz = 1;         ///< number of elements along z (>= 1, 3-D only)
  Scalar lx = 1.0;      ///< domain length along x [m] (> 0)
  Scalar ly = 1.0;      ///< domain length along y [m] (> 0)
  Scalar lz = 1.0;      ///< domain length along z [m] (> 0, 3-D only)
  Scalar x0 = 0.0;      ///< x coordinate of the lower corner [m]
  Scalar y0 = 0.0;      ///< y coordinate of the lower corner [m]
  Scalar z0 = 0.0;      ///< z coordinate of the lower corner [m] (3-D only)
};

/// Build a uniform structured Q4 mesh. The resulting mesh is validated before
/// being returned and carries StructuredGridInfo so downstream components can
/// exploit the identical-cell property.
/// \throws ConfigError for non-positive counts or lengths.
Mesh make_structured_quad_mesh(const StructuredMeshSpec& spec);

/// Build a uniform structured Hex8 mesh; otherwise as the Q4 generator.
Mesh make_structured_hex_mesh(const StructuredMeshSpec& spec);

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

/// 3-D counterpart of `make_perturbed_quad_mesh`: interior nodes of a
/// structured Hex8 grid are displaced by up to `perturbation` of the cell
/// size in every direction, boundary nodes stay in place.
/// \param perturbation fraction of the cell size, in [0, 0.35).
Mesh make_perturbed_hex_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                             unsigned int seed = 12345u);

/// Convenience accessors for structured grids (used by tests and selectors).
/// The two-index forms address a 2-D grid; the three-index forms a 3-D grid.
/// \{
Index structured_node_index(const StructuredGridInfo& info, Index i, Index j);
Index structured_element_index(const StructuredGridInfo& info, Index i, Index j);
Index structured_node_index(const StructuredGridInfo& info, Index i, Index j, Index k);
Index structured_element_index(const StructuredGridInfo& info, Index i, Index j, Index k);
/// \}

}  // namespace sparlab
