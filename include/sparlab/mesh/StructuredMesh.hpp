/// \file StructuredMesh.hpp
/// \brief Generators for structured Q4 and Hex8 grids, the triangle and
///        tetrahedron meshes obtained by splitting their cells, and structured
///        shell surfaces (plates, cylinders, spheres).
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
///
/// The simplex generators split every cell of the corresponding grid: a
/// quadrilateral into two triangles along alternating diagonals (so the mesh
/// has no preferred direction), a hexahedron into the six Kuhn tetrahedra
/// that share its main diagonal from corner 0 to corner 6. Both splits are
/// conforming across cells, and the elements of cell `c` are numbered
/// `2 c, 2 c + 1` and `6 c ... 6 c + 5` respectively. They exist for
/// verification (the same box as the Q4 / Hex8 studies, deterministic) and
/// for decks that want simplices on a box; real geometry comes from a mesh
/// file (MeshReader.hpp).
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

/// Structured Tri3 mesh: every Q4 cell of `make_structured_quad_mesh(spec)`
/// split into two counter-clockwise triangles, along the diagonal 0-2 in cells
/// with i + j even and 1-3 otherwise.
Mesh make_structured_tri_mesh(const StructuredMeshSpec& spec);

/// Structured Tet4 mesh: every Hex8 cell of `make_structured_hex_mesh(spec)`
/// split into six positively oriented Kuhn tetrahedra.
Mesh make_structured_tet_mesh(const StructuredMeshSpec& spec);

/// Structured Tet10 mesh: `make_structured_tet_mesh(spec)` with a node at the
/// midpoint of every edge (`elevate_to_tet10`).
Mesh make_structured_tet10_mesh(const StructuredMeshSpec& spec);

/// The perturbed quad / hex meshes above, split the same way (identical node
/// positions for the same seed), for patch tests on distorted simplices.
/// \{
Mesh make_perturbed_tri_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                             unsigned int seed = 12345u);
Mesh make_perturbed_tet_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                             unsigned int seed = 12345u);
/// Straight-sided Tet10 cells on the perturbed Tet4 mesh (edge nodes at the
/// midpoints of the distorted edges), so every cell map stays affine and the
/// element space contains every quadratic field.
Mesh make_perturbed_tet10_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                               unsigned int seed = 12345u);
/// \}

/// Surfaces of the structured shell generator.
enum class ShellShape {
  Plate,     ///< a rectangle in the plane z = origin.z
  Cylinder,  ///< a cylinder, or a panel of one, about a global axis
  Sphere     ///< a zone of a sphere between two polar angles
};

std::string to_string(ShellShape shape);
ShellShape parse_shell_shape(const std::string& text);

/// A structured MITC4 shell surface. Direction 1 runs along x (plate),
/// around the axis (cylinder) or in longitude (sphere); direction 2 along y,
/// along the axis, or in latitude, from the larger polar angle to the
/// smaller. Nodes are numbered `j * (n1 + 1) + i` (`j * n1 + i` when
/// direction 1 closes on itself, a full 360 degrees), and every element runs
/// (i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1), so its normal
/// \f$g_1\times g_2\f$ points along +z on a plate and away from the axis or
/// the centre on a cylinder or a sphere - the side a positive pressure
/// presses on. The mesh carries the surface's exact normals at its nodes.
struct ShellMeshSpec {
  ShellShape shape = ShellShape::Plate;
  Index n1 = 1;  ///< elements along direction 1 (>= 1; >= 3 around a closed circle)
  Index n2 = 1;  ///< elements along direction 2 (>= 1)
  /// Plate: its lower corner. Cylinder: the point of its axis where it
  /// starts. Sphere: its centre [m].
  Vector3 origin = Vector3::Zero();
  Scalar lx = 1.0;      ///< plate extent along x [m]
  Scalar ly = 1.0;      ///< plate extent along y [m]
  /// Plate: interior nodes moved in the plane by up to this fraction of the
  /// cell size along x and along y, in [0, 0.25), which keeps every cell
  /// convex (distorted meshes for patch tests).
  Scalar perturbation = 0.0;
  unsigned int seed = 12345u;
  Scalar radius = 1.0;  ///< cylinder and sphere [m]
  Scalar length = 1.0;  ///< cylinder, along its axis [m]
  /// Cylinder axis: 0 = x, 1 = y, 2 = z. The angle around it starts at the
  /// next axis and turns towards the one after (about z: from +x towards +y;
  /// about x: from +y towards +z; about y: from +z towards +x).
  int axis = 2;
  /// Cylinder: the angles around the axis; sphere: the longitudes about +z
  /// from +x [degrees]. A span of 360 closes the surface.
  Scalar angle_start = 0.0;
  Scalar angle_end = 360.0;
  /// Sphere: the polar angles from +z [degrees], 0 < start < end < 180 (a
  /// pole would collapse the quadrilaterals around it).
  Scalar polar_start = 30.0;
  Scalar polar_end = 150.0;
};

/// Build a structured shell mesh (Shell4 cells in 3-D) with exact node normals.
/// \throws ConfigError for out-of-range counts, extents or angles.
Mesh make_structured_shell_mesh(const ShellMeshSpec& spec);

/// Convenience accessors for structured grids (used by tests and selectors).
/// The two-index forms address a 2-D grid; the three-index forms a 3-D grid.
/// \{
Index structured_node_index(const StructuredGridInfo& info, Index i, Index j);
Index structured_element_index(const StructuredGridInfo& info, Index i, Index j);
Index structured_node_index(const StructuredGridInfo& info, Index i, Index j, Index k);
Index structured_element_index(const StructuredGridInfo& info, Index i, Index j, Index k);
/// \}

}  // namespace sparlab
