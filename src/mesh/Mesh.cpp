#include "sparlab/mesh/Mesh.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

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
  }
  throw MeshError("unhandled element type");
}

std::string to_string(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return "Quad4";
  }
  return "Unknown";
}

Mesh::Mesh(Eigen::Matrix2Xd coords, std::vector<Index> connectivity, ElementType type)
    : coords_(std::move(coords)),
      connectivity_(std::move(connectivity)),
      type_(type),
      nodes_per_element_(nodes_per_element(type)) {
  if (connectivity_.size() % static_cast<std::size_t>(nodes_per_element_) != 0) {
    std::ostringstream os;
    os << "connectivity length " << connectivity_.size() << " is not a multiple of "
       << nodes_per_element_ << " nodes per " << to_string(type) << " element";
    throw MeshError(os.str());
  }
}

Eigen::Matrix<Scalar, 2, Eigen::Dynamic> Mesh::element_coordinates(Index e) const {
  Eigen::Matrix<Scalar, 2, Eigen::Dynamic> xe(2, nodes_per_element_);
  const Index* nodes = element_nodes(e);
  for (int a = 0; a < nodes_per_element_; ++a) xe.col(a) = coords_.col(nodes[a]);
  return xe;
}

Vector2 Mesh::element_centroid(Index e) const {
  Vector2 c = Vector2::Zero();
  const Index* nodes = element_nodes(e);
  for (int a = 0; a < nodes_per_element_; ++a) c += coords_.col(nodes[a]);
  return c / static_cast<Scalar>(nodes_per_element_);
}

Scalar Mesh::element_area(Index e) const {
  const Index* nodes = element_nodes(e);
  Scalar twice_area = 0.0;
  for (int a = 0; a < nodes_per_element_; ++a) {
    const Vector2 p = coords_.col(nodes[a]);
    const Vector2 q = coords_.col(nodes[(a + 1) % nodes_per_element_]);
    twice_area += p.x() * q.y() - q.x() * p.y();
  }
  return 0.5 * twice_area;
}

Eigen::Vector4d Mesh::bounding_box() const {
  if (num_nodes() == 0) throw MeshError("bounding box requested for an empty mesh");
  Eigen::Vector4d bb;
  bb << coords_.row(0).minCoeff(), coords_.row(1).minCoeff(),
      coords_.row(0).maxCoeff(), coords_.row(1).maxCoeff();
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
  Scalar min_area = std::numeric_limits<Scalar>::max();
  Scalar max_area = 0.0;
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
    const Scalar area = element_area(e);
    if (!(area > 0.0)) {
      std::ostringstream os;
      os << "element " << e << " has signed area " << area
         << " m^2; nodes must be listed counter-clockwise and the element must not be "
            "inverted or collapsed";
      throw MeshError(os.str());
    }
    min_area = std::min(min_area, area);
    max_area = std::max(max_area, area);
  }

  if (max_area / min_area > 1.0e8) {
    log::warn("extreme element area ratio ", max_area / min_area,
              " (min ", min_area, " m^2, max ", max_area,
              " m^2); conditioning of the stiffness matrix may suffer");
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

std::vector<Mesh::BoundaryEdge> Mesh::boundary_edges() const {
  // Local edge node pairs for the supported topologies.
  static const int quad_edges[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  if (type_ != ElementType::Quad4) throw MeshError("boundary_edges supports Quad4 only");

  std::map<std::pair<Index, Index>, std::vector<BoundaryEdge>> edge_map;
  const Index ne = num_elements();
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = element_nodes(e);
    for (int le = 0; le < 4; ++le) {
      BoundaryEdge edge;
      edge.node_a = nodes[quad_edges[le][0]];
      edge.node_b = nodes[quad_edges[le][1]];
      edge.element = e;
      edge.local_edge = le;
      const auto key = std::minmax(edge.node_a, edge.node_b);
      edge_map[{key.first, key.second}].push_back(edge);
    }
  }

  std::vector<BoundaryEdge> boundary;
  for (const auto& entry : edge_map) {
    if (entry.second.size() == 1) boundary.push_back(entry.second.front());
  }
  return boundary;
}

}  // namespace sparlab
