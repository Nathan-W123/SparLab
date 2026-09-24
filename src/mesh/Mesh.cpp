#include "sparlab/mesh/Mesh.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/Hex8.hpp"
#include "sparlab/elements/Quad4.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

namespace sparlab {

int nodes_per_element(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return 4;
    case ElementType::Hex8: return 8;
  }
  throw MeshError("unhandled element type");
}

int element_dimension(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return 2;
    case ElementType::Hex8: return 3;
  }
  throw MeshError("unhandled element type");
}

std::string to_string(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return "Quad4";
    case ElementType::Hex8: return "Hex8";
  }
  return "Unknown";
}

const std::vector<std::vector<int>>& element_local_faces(ElementType type) {
  // Q4 edges run counter-clockwise; Hex8 faces follow the VTK hexahedron
  // convention (bottom, top, front, right, back, left) with each face's nodes
  // listed so that the face parametrisation is well defined.
  static const std::vector<std::vector<int>> quad_edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  static const std::vector<std::vector<int>> hex_faces = {
      {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  switch (type) {
    case ElementType::Quad4: return quad_edges;
    case ElementType::Hex8: return hex_faces;
  }
  throw MeshError("unhandled element type");
}

Mesh::Mesh(Matrix coords, std::vector<Index> connectivity, ElementType type)
    : coords_(std::move(coords)),
      connectivity_(std::move(connectivity)),
      type_(type),
      nodes_per_element_(nodes_per_element(type)),
      dim_(element_dimension(type)) {
  if (coords_.rows() != dim_) {
    std::ostringstream os;
    os << "coordinate array has " << coords_.rows() << " rows but " << to_string(type)
       << " elements live in " << dim_ << " dimensions";
    throw MeshError(os.str());
  }
  if (connectivity_.size() % static_cast<std::size_t>(nodes_per_element_) != 0) {
    std::ostringstream os;
    os << "connectivity length " << connectivity_.size() << " is not a multiple of "
       << nodes_per_element_ << " nodes per " << to_string(type) << " element";
    throw MeshError(os.str());
  }
}

Matrix Mesh::element_coordinates(Index e) const {
  Matrix xe(dim_, nodes_per_element_);
  const Index* nodes = element_nodes(e);
  for (int a = 0; a < nodes_per_element_; ++a) xe.col(a) = coords_.col(nodes[a]);
  return xe;
}

Vector3 Mesh::element_centroid(Index e) const {
  Vector3 c = Vector3::Zero();
  const Index* nodes = element_nodes(e);
  for (int a = 0; a < nodes_per_element_; ++a) c.head(dim_) += coords_.col(nodes[a]);
  return c / static_cast<Scalar>(nodes_per_element_);
}

Scalar Mesh::element_measure(Index e) const {
  const Index* nodes = element_nodes(e);
  if (dim_ == 2) {
    Scalar twice_area = 0.0;
    for (int a = 0; a < nodes_per_element_; ++a) {
      const Vector2 p = coords_.col(nodes[a]);
      const Vector2 q = coords_.col(nodes[(a + 1) % nodes_per_element_]);
      twice_area += p.x() * q.y() - q.x() * p.y();
    }
    return 0.5 * twice_area;
  }
  Scalar min_det = 0.0;
  return hex8_volume(element_coordinates(e), &min_det);
}

BoundingBox Mesh::bounding_box() const {
  if (num_nodes() == 0) throw MeshError("bounding box requested for an empty mesh");
  BoundingBox bb;
  for (int i = 0; i < dim_; ++i) {
    bb.lower(i) = coords_.row(i).minCoeff();
    bb.upper(i) = coords_.row(i).maxCoeff();
  }
  return bb;
}

void Mesh::validate() const {
  if (num_nodes() == 0) throw MeshError("mesh contains no nodes");
  if (num_elements() == 0) throw MeshError("mesh contains no elements");
  if (!coords_.allFinite()) throw MeshError("nodal coordinates contain NaN or Inf");

  const Index nn = num_nodes();
  for (std::size_t k = 0; k < connectivity_.size(); ++k) {
    const Index node = connectivity_[k];
    if (node < 0 || node >= nn) {
      std::ostringstream os;
      os << "element " << (k / nodes_per_element_) << " references node " << node
         << " which is outside the valid range [0, " << nn - 1 << "]";
      throw MeshError(os.str());
    }
  }

  const Index ne = num_elements();
  const char* unit = dim_ == 2 ? " m^2" : " m^3";
  Scalar min_measure = std::numeric_limits<Scalar>::max();
  Scalar max_measure = 0.0;
  for (Index e = 0; e < ne; ++e) {
    // Duplicated node within one element => degenerate topology.
    const Index* nodes = element_nodes(e);
    for (int a = 0; a < nodes_per_element_; ++a) {
      for (int b = a + 1; b < nodes_per_element_; ++b) {
        if (nodes[a] == nodes[b]) {
          std::ostringstream os;
          os << "element " << e << " uses node " << nodes[a] << " at local positions "
             << a << " and " << b << " (collapsed/degenerate element)";
          throw MeshError(os.str());
        }
      }
    }
    if (dim_ == 2) {
      const Scalar area = element_measure(e);
      if (!(area > 0.0)) {
        std::ostringstream os;
        os << "element " << e << " has signed area " << area
           << " m^2; nodes must be listed counter-clockwise and the element must not be "
              "inverted or collapsed";
        throw MeshError(os.str());
      }
      min_measure = std::min(min_measure, area);
      max_measure = std::max(max_measure, area);
    } else {
      // A hexahedron can have positive volume and still be folded at a
      // corner, so the Jacobian is checked at every integration point rather
      // than only through the volume.
      Scalar min_det = 0.0;
      const Scalar volume = hex8_volume(element_coordinates(e), &min_det);
      if (!(min_det > 0.0) || !(volume > 0.0)) {
        std::ostringstream os;
        os << "element " << e << " has volume " << volume
           << " m^3 and a minimum Jacobian determinant of " << min_det
           << " m^3 over its 2 x 2 x 2 Gauss points; nodes must follow the VTK "
              "hexahedron ordering (bottom face counter-clockwise seen from +z, then "
              "the top face) and the element must not be inverted or folded";
        throw MeshError(os.str());
      }
      min_measure = std::min(min_measure, volume);
      max_measure = std::max(max_measure, volume);
    }
  }

  if (max_measure / min_measure > 1.0e8) {
    log::warn("extreme element ", (dim_ == 2 ? "area" : "volume"), " ratio ",
              max_measure / min_measure, " (min ", min_measure, unit, ", max ",
              max_measure, unit, "); conditioning of the stiffness matrix may suffer");
  }

  // Nodes referenced by no element contribute empty rows/columns to K and are a
  // frequent cause of singular systems.
  std::vector<char> used(static_cast<std::size_t>(nn), 0);
  for (Index node : connectivity_) used[static_cast<std::size_t>(node)] = 1;
  const auto orphans = std::count(used.begin(), used.end(), 0);
  if (orphans > 0) {
    std::ostringstream os;
    os << orphans << " of " << nn
       << " nodes are not referenced by any element; their degrees of freedom would "
          "make the stiffness matrix singular";
    throw MeshError(os.str());
  }
}

std::vector<Mesh::BoundaryFace> Mesh::boundary_faces() const {
  const std::vector<std::vector<int>>& local = element_local_faces(type_);

  // Faces are keyed by their sorted node list, so orientation does not matter.
  std::map<std::vector<Index>, std::vector<BoundaryFace>> face_map;
  const Index ne = num_elements();
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = element_nodes(e);
    for (std::size_t lf = 0; lf < local.size(); ++lf) {
      BoundaryFace face;
      face.element = e;
      face.local_face = static_cast<int>(lf);
      for (int a : local[lf]) face.nodes.push_back(nodes[a]);
      std::vector<Index> key = face.nodes;
      std::sort(key.begin(), key.end());
      face_map[key].push_back(std::move(face));
    }
  }

  std::vector<BoundaryFace> boundary;
  for (const auto& entry : face_map) {
    if (entry.second.size() == 1) boundary.push_back(entry.second.front());
  }
  return boundary;
}

Scalar Mesh::face_measure(const BoundaryFace& face) const {
  if (face.nodes.size() == 2) {
    return (coords_.col(face.nodes[1]) - coords_.col(face.nodes[0])).norm();
  }
  if (face.nodes.size() == 4 && dim_ == 3) {
    Matrix xf(3, 4);
    for (int a = 0; a < 4; ++a) xf.col(a) = coords_.col(face.nodes[static_cast<std::size_t>(a)]);
    return hex8_face_area(xf);
  }
  throw MeshError("face_measure received a face with an unsupported node count");
}

}  // namespace sparlab
