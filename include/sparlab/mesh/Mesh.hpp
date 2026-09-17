/// \file Mesh.hpp
/// \brief Storage-only mesh container plus validation.
///
/// The mesh owns nodal coordinates and element connectivity and nothing else:
/// materials, loads and degrees of freedom live in separate components. A mesh
/// carries a single element type today (Q4) but the connectivity is stored as a
/// flat array with an explicit stride, so mixed-topology meshes only require
/// replacing the stride with a per-element offset table.
#pragma once

#include "sparlab/core/Types.hpp"

#include <optional>
#include <string>
#include <vector>

namespace sparlab {

/// Supported element topologies.
enum class ElementType {
  Quad4  ///< Four-node bilinear isoparametric quadrilateral.
};

/// Number of nodes carried by an element topology.
int nodes_per_element(ElementType type);

/// Human-readable name, also used in VTK output.
std::string to_string(ElementType type);

/// Metadata retained when a mesh came from the structured generator. Presence
/// of this record allows the assembler to reuse a single element stiffness
/// matrix for the whole grid (all cells are geometrically identical).
struct StructuredGridInfo {
  Index nx = 0;       ///< elements along x
  Index ny = 0;       ///< elements along y
  Scalar lx = 0.0;    ///< domain length along x [m]
  Scalar ly = 0.0;    ///< domain length along y [m]
  bool uniform = false;  ///< true when every cell is the same rectangle
};

/// Node coordinates + element connectivity for a 2-D mesh.
class Mesh {
 public:
  Mesh() = default;

  /// \param coords 2 x num_nodes matrix of nodal coordinates [m].
  /// \param connectivity flat, row-major connectivity with
  ///        `nodes_per_element(type)` entries per element, counter-clockwise.
  Mesh(Eigen::Matrix2Xd coords, std::vector<Index> connectivity, ElementType type);

  Index num_nodes() const { return static_cast<Index>(coords_.cols()); }
  Index num_elements() const {
    return static_cast<Index>(connectivity_.size()) / nodes_per_element_;
  }
  int nodes_per_elem() const { return nodes_per_element_; }
  int dim() const { return kDim; }
  ElementType element_type() const { return type_; }

  /// Coordinates of node `n` [m].
  Vector2 node(Index n) const { return coords_.col(n); }

  const Eigen::Matrix2Xd& coordinates() const { return coords_; }

  /// Global node indices of element `e` (pointer to `nodes_per_elem()` values).
  const Index* element_nodes(Index e) const {
    return connectivity_.data() + static_cast<std::size_t>(e) * nodes_per_element_;
  }

  const std::vector<Index>& connectivity() const { return connectivity_; }

  /// 2 x nodes_per_elem matrix of the element's nodal coordinates [m].
  Eigen::Matrix<Scalar, 2, Eigen::Dynamic> element_coordinates(Index e) const;

  /// Geometric centroid of element `e` (average of its nodes) [m].
  Vector2 element_centroid(Index e) const;

  /// Signed area of element `e` [m^2] via the shoelace formula. Positive for
  /// counter-clockwise node ordering.
  Scalar element_area(Index e) const;

  /// Axis-aligned bounding box of the mesh: (xmin, ymin, xmax, ymax) [m].
  Eigen::Vector4d bounding_box() const;

  void set_structured_info(const StructuredGridInfo& info) { structured_ = info; }
  const std::optional<StructuredGridInfo>& structured_info() const { return structured_; }

  /// Validate connectivity ranges, element areas and node degeneracy.
  /// \throws MeshError with an actionable message on the first problem found.
  void validate() const;

  /// Boundary edge of the mesh: the two node indices plus the owning element.
  struct BoundaryEdge {
    Index node_a = 0;
    Index node_b = 0;
    Index element = 0;
    int local_edge = 0;  ///< local edge id within the element (0..3 for Q4)
  };

  /// Edges referenced by exactly one element, i.e. the mesh boundary.
  /// Computed on demand; O(num_elements * nodes_per_elem).
  std::vector<BoundaryEdge> boundary_edges() const;

 private:
  Eigen::Matrix2Xd coords_;
  std::vector<Index> connectivity_;
  ElementType type_ = ElementType::Quad4;
  int nodes_per_element_ = 4;
  std::optional<StructuredGridInfo> structured_;
};

}  // namespace sparlab
