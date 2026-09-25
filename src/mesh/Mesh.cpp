#include "sparlab/mesh/Mesh.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/Hex8.hpp"
#include "sparlab/elements/Quad4.hpp"
#include "sparlab/elements/Tet10.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <utility>

namespace sparlab {

int nodes_per_element(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return 4;
    case ElementType::Hex8: return 8;
    case ElementType::Tri3: return 3;
    case ElementType::Tet4: return 4;
    case ElementType::Tet10: return 10;
  }
  throw MeshError("unhandled element type");
}

int element_dimension(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return 2;
    case ElementType::Hex8: return 3;
    case ElementType::Tri3: return 2;
    case ElementType::Tet4: return 3;
    case ElementType::Tet10: return 3;
  }
  throw MeshError("unhandled element type");
}

bool is_simplex(ElementType type) {
  return type == ElementType::Tri3 || type == ElementType::Tet4 ||
         type == ElementType::Tet10;
}

int corner_nodes(ElementType type) {
  return type == ElementType::Tet10 ? 4 : nodes_per_element(type);
}

int face_corner_nodes(ElementType type) {
  switch (type) {
    case ElementType::Quad4:
    case ElementType::Tri3: return 2;
    case ElementType::Hex8: return 4;
    case ElementType::Tet4:
    case ElementType::Tet10: return 3;
  }
  throw MeshError("unhandled element type");
}

std::string to_string(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return "Quad4";
    case ElementType::Hex8: return "Hex8";
    case ElementType::Tri3: return "Tri3";
    case ElementType::Tet4: return "Tet4";
    case ElementType::Tet10: return "Tet10";
  }
  return "Unknown";
}

const std::vector<std::vector<int>>& element_local_faces(ElementType type) {
  // Q4 and Tri3 edges run counter-clockwise; Hex8 faces follow the VTK
  // hexahedron convention (bottom, top, front, right, back, left); Tet4 faces
  // are listed opposite nodes 3, 2, 0 and 1, and a Tet10 face adds the edge
  // nodes of its corner pairs (0,1), (1,2), (2,0) in that order. Every face
  // is wound so that its right-hand normal points out of a positively
  // oriented element, which the STL export relies on.
  static const std::vector<std::vector<int>> quad_edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  static const std::vector<std::vector<int>> hex_faces = {
      {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  static const std::vector<std::vector<int>> tri_edges = {{0, 1}, {1, 2}, {2, 0}};
  static const std::vector<std::vector<int>> tet_faces = {
      {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}};
  static const std::vector<std::vector<int>> tet10_faces = {
      {0, 2, 1, 6, 5, 4}, {0, 1, 3, 4, 8, 7}, {1, 2, 3, 5, 9, 8}, {0, 3, 2, 7, 9, 6}};
  switch (type) {
    case ElementType::Quad4: return quad_edges;
    case ElementType::Hex8: return hex_faces;
    case ElementType::Tri3: return tri_edges;
    case ElementType::Tet4: return tet_faces;
    case ElementType::Tet10: return tet10_faces;
  }
  throw MeshError("unhandled element type");
}

Scalar tet4_volume(const Matrix& x) {
  if (x.rows() != 3 || x.cols() != 4) {
    std::ostringstream os;
    os << "tet4_volume expects a 3 x 4 coordinate matrix, received " << x.rows() << " x "
       << x.cols();
    throw MeshError(os.str());
  }
  const Vector3 a = x.col(1) - x.col(0);
  const Vector3 b = x.col(2) - x.col(0);
  const Vector3 c = x.col(3) - x.col(0);
  return a.dot(b.cross(c)) / 6.0;
}

Mesh elevate_to_tet10(const Mesh& linear) {
  if (linear.element_type() != ElementType::Tet4) {
    throw MeshError("only a Tet4 mesh can be elevated to Tet10; this mesh holds " +
                    to_string(linear.element_type()) + " cells");
  }
  static const int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
  const Index nn = linear.num_nodes();
  const Index ne = linear.num_elements();

  // Unique edges keyed by their sorted end nodes, in order of first
  // appearance; `edge_of[6 e + k]` is the edge of local edge k of cell e.
  struct EdgeRecord {
    Index a;
    Index b;
    std::size_t slot;  // 6 e + k
  };
  std::vector<EdgeRecord> records(static_cast<std::size_t>(ne) * 6);
  for (Index e = 0; e < ne; ++e) {
    const Index* n = linear.element_nodes(e);
    for (int k = 0; k < 6; ++k) {
      const Index p = n[edges[k][0]];
      const Index q = n[edges[k][1]];
      const std::size_t slot = static_cast<std::size_t>(e) * 6 + static_cast<std::size_t>(k);
      records[slot] = {std::min(p, q), std::max(p, q), slot};
    }
  }
  std::vector<EdgeRecord> sorted = records;
  std::sort(sorted.begin(), sorted.end(), [](const EdgeRecord& x, const EdgeRecord& y) {
    if (x.a != y.a) return x.a < y.a;
    if (x.b != y.b) return x.b < y.b;
    return x.slot < y.slot;
  });
  // Group equal edges; each group's first slot (the smallest) fixes the
  // numbering by first appearance.
  std::vector<Index> group_of(records.size(), -1);
  std::vector<std::size_t> group_first;
  for (std::size_t i = 0; i < sorted.size();) {
    std::size_t j = i;
    while (j < sorted.size() && sorted[j].a == sorted[i].a && sorted[j].b == sorted[i].b) {
      group_of[sorted[j].slot] = static_cast<Index>(group_first.size());
      ++j;
    }
    group_first.push_back(sorted[i].slot);
    i = j;
  }
  std::vector<Index> order(group_first.size());
  std::iota(order.begin(), order.end(), Index{0});
  std::sort(order.begin(), order.end(), [&](Index x, Index y) {
    return group_first[static_cast<std::size_t>(x)] < group_first[static_cast<std::size_t>(y)];
  });
  std::vector<Index> edge_node(group_first.size());
  for (std::size_t r = 0; r < order.size(); ++r) {
    edge_node[static_cast<std::size_t>(order[r])] = nn + static_cast<Index>(r);
  }
  const Index num_edges = static_cast<Index>(group_first.size());

  Matrix coords(3, nn + num_edges);
  coords.leftCols(nn) = linear.coordinates();
  std::vector<std::array<Index, 2>> ends(static_cast<std::size_t>(num_edges));
  for (std::size_t g = 0; g < group_first.size(); ++g) {
    const EdgeRecord& r = records[group_first[g]];
    const Index node = edge_node[g];
    coords.col(node) = 0.5 * (linear.coordinates().col(r.a) + linear.coordinates().col(r.b));
    ends[static_cast<std::size_t>(node - nn)] = {r.a, r.b};
  }
  std::vector<Index> connectivity;
  connectivity.reserve(static_cast<std::size_t>(ne) * 10);
  for (Index e = 0; e < ne; ++e) {
    const Index* n = linear.element_nodes(e);
    connectivity.insert(connectivity.end(), n, n + 4);
    for (int k = 0; k < 6; ++k) {
      const std::size_t slot = static_cast<std::size_t>(e) * 6 + static_cast<std::size_t>(k);
      connectivity.push_back(edge_node[static_cast<std::size_t>(group_of[slot])]);
    }
  }
  Mesh quadratic(std::move(coords), std::move(connectivity), ElementType::Tet10);

  for (const auto& entry : linear.element_sets()) {
    quadratic.set_element_set(entry.first, entry.second);
  }
  if (!linear.node_sets().empty()) {
    // Edge node of the edge (a, b), from the quadratic cells.
    const auto edge_node_of = [&](Index cell, int k) {
      return quadratic.element_nodes(cell)[4 + k];
    };
    const std::vector<Mesh::BoundaryFace> boundary = linear.boundary_faces();
    // Local edge index of the face corner pair (i, j) of a Tet4 face.
    const auto local_edge = [](int i, int j) {
      for (int k = 0; k < 6; ++k) {
        if ((edges[k][0] == i && edges[k][1] == j) || (edges[k][0] == j && edges[k][1] == i)) {
          return k;
        }
      }
      return -1;
    };
    for (const auto& entry : linear.node_sets()) {
      std::vector<char> in_set(static_cast<std::size_t>(nn), 0);
      for (Index n : entry.second) in_set[static_cast<std::size_t>(n)] = 1;
      std::vector<Index> nodes = entry.second;
      // (a) the edges of the cells of an element set with the same name.
      const auto cells = linear.element_sets().find(entry.first);
      if (cells != linear.element_sets().end()) {
        for (Index e : cells->second) {
          for (int k = 0; k < 6; ++k) nodes.push_back(edge_node_of(e, k));
        }
      }
      // (b) the edges of boundary faces lying wholly in the set.
      bool covers_face = false;
      std::vector<Index> boundary_edge_nodes;
      for (const Mesh::BoundaryFace& face : boundary) {
        const std::vector<int>& local = element_local_faces(ElementType::Tet4)
            [static_cast<std::size_t>(face.local_face)];
        const bool whole = in_set[static_cast<std::size_t>(face.nodes[0])] &&
                           in_set[static_cast<std::size_t>(face.nodes[1])] &&
                           in_set[static_cast<std::size_t>(face.nodes[2])];
        for (int c = 0; c < 3; ++c) {
          const int i = local[static_cast<std::size_t>(c)];
          const int j = local[static_cast<std::size_t>((c + 1) % 3)];
          const Index node = edge_node_of(face.element, local_edge(i, j));
          if (whole) nodes.push_back(node);
          const std::array<Index, 2>& ab = ends[static_cast<std::size_t>(node - nn)];
          if (in_set[static_cast<std::size_t>(ab[0])] && in_set[static_cast<std::size_t>(ab[1])]) {
            boundary_edge_nodes.push_back(node);
          }
        }
        covers_face = covers_face || whole;
      }
      // (c) a curve or point group: boundary edges with both ends in it.
      if (!covers_face) nodes.insert(nodes.end(), boundary_edge_nodes.begin(), boundary_edge_nodes.end());
      quadratic.set_node_set(entry.first, std::move(nodes));
    }
  }
  quadratic.validate();
  return quadratic;
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
  if (type_ == ElementType::Tet4) return tet4_volume(element_coordinates(e));
  if (type_ == ElementType::Tet10) return tet10_volume(element_coordinates(e));
  Scalar min_det = 0.0;
  return hex8_volume(element_coordinates(e), &min_det);
}

Scalar Mesh::mean_element_size() const {
  const Index ne = num_elements();
  if (ne == 0) throw MeshError("mean element size requested for a mesh with no elements");
  if (!is_simplex(type_)) {
    // Side of the square (2-D) or cube (3-D) with the mean cell measure.
    Scalar measure_sum = 0.0;
    for (Index e = 0; e < ne; ++e) measure_sum += element_measure(e);
    const Scalar mean_measure = measure_sum / static_cast<Scalar>(ne);
    return dim_ == 2 ? std::sqrt(mean_measure) : std::cbrt(mean_measure);
  }
  // Mean edge length, the quantity a mesh generator's size field controls;
  // the corners span a Tet10's edges.
  const int npe = corner_nodes(type_);
  Scalar total = 0.0;
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = element_nodes(e);
    Scalar sum = 0.0;
    int count = 0;
    for (int a = 0; a < npe; ++a) {
      for (int b = a + 1; b < npe; ++b) {
        sum += (coords_.col(nodes[a]) - coords_.col(nodes[b])).norm();
        ++count;
      }
    }
    total += sum / static_cast<Scalar>(count);
  }
  return total / static_cast<Scalar>(ne);
}

MeshQuality Mesh::quality() const {
  MeshQuality q;
  const Index ne = num_elements();
  if (ne == 0) return q;
  static const int hex_corner[8][3] = {{1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7},
                                       {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};
  switch (type_) {
    case ElementType::Quad4:
    case ElementType::Hex8: q.metric = "minimum scaled Jacobian at the corners"; break;
    case ElementType::Tri3: q.metric = "4 sqrt(3) A / sum of squared edge lengths"; break;
    case ElementType::Tet4: q.metric = "6 sqrt(2) V / rms edge length cubed"; break;
    case ElementType::Tet10:
      q.metric = "6 sqrt(2) V / rms edge length cubed of the corners, times the nodal "
                 "Jacobian ratio";
      break;
  }
  q.min = std::numeric_limits<Scalar>::max();
  Scalar sum = 0.0;
  for (Index e = 0; e < ne; ++e) {
    const Index* n = element_nodes(e);
    Scalar value = 0.0;
    if (type_ == ElementType::Quad4) {
      value = std::numeric_limits<Scalar>::max();
      for (int a = 0; a < 4; ++a) {
        const Vector2 p = coords_.col(n[a]);
        const Vector2 e1 = Vector2(coords_.col(n[(a + 1) % 4])) - p;
        const Vector2 e2 = Vector2(coords_.col(n[(a + 3) % 4])) - p;
        const Scalar denom = e1.norm() * e2.norm();
        const Scalar sj =
            denom > 0.0 ? (e1.x() * e2.y() - e1.y() * e2.x()) / denom : -1.0;
        value = std::min(value, sj);
      }
    } else if (type_ == ElementType::Hex8) {
      value = std::numeric_limits<Scalar>::max();
      for (int a = 0; a < 8; ++a) {
        const Vector3 p = coords_.col(n[a]);
        const Vector3 e1 = Vector3(coords_.col(n[hex_corner[a][0]])) - p;
        const Vector3 e2 = Vector3(coords_.col(n[hex_corner[a][1]])) - p;
        const Vector3 e3 = Vector3(coords_.col(n[hex_corner[a][2]])) - p;
        const Scalar denom = e1.norm() * e2.norm() * e3.norm();
        const Scalar sj = denom > 0.0 ? e1.dot(e2.cross(e3)) / denom : -1.0;
        value = std::min(value, sj);
      }
    } else if (type_ == ElementType::Tri3) {
      Scalar l2 = 0.0;
      for (int a = 0; a < 3; ++a) {
        l2 += (coords_.col(n[a]) - coords_.col(n[(a + 1) % 3])).squaredNorm();
      }
      value = l2 > 0.0 ? 4.0 * std::sqrt(3.0) * element_measure(e) / l2 : -1.0;
    } else {
      Scalar l2 = 0.0;
      for (int a = 0; a < 4; ++a) {
        for (int b = a + 1; b < 4; ++b) {
          l2 += (coords_.col(n[a]) - coords_.col(n[b])).squaredNorm();
        }
      }
      const Scalar lrms = std::sqrt(l2 / 6.0);
      const Matrix xe = element_coordinates(e);
      const Scalar corner_volume = tet4_volume(xe.leftCols(4));
      value = lrms > 0.0 ? 6.0 * std::sqrt(2.0) * corner_volume / (lrms * lrms * lrms)
                         : -1.0;
      if (type_ == ElementType::Tet10) value *= std::max(tet10_jacobian_ratio(xe), -1.0);
    }
    sum += value;
    if (value < q.min) {
      q.min = value;
      q.worst_element = e;
    }
    if (value < q.poor_threshold) ++q.poor_elements;
  }
  q.mean = sum / static_cast<Scalar>(ne);
  return q;
}

void Mesh::set_node_set(const std::string& name, std::vector<Index> nodes) {
  if (name.empty()) throw MeshError("a node set needs a non-empty name");
  std::sort(nodes.begin(), nodes.end());
  nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
  if (!nodes.empty() && (nodes.front() < 0 || nodes.back() >= num_nodes())) {
    std::ostringstream os;
    os << "node set '" << name << "' references node "
       << (nodes.front() < 0 ? nodes.front() : nodes.back()) << " outside [0, "
       << num_nodes() - 1 << "]";
    throw MeshError(os.str());
  }
  node_sets_[name] = std::move(nodes);
}

void Mesh::set_element_set(const std::string& name, std::vector<Index> elements) {
  if (name.empty()) throw MeshError("an element set needs a non-empty name");
  std::sort(elements.begin(), elements.end());
  elements.erase(std::unique(elements.begin(), elements.end()), elements.end());
  if (!elements.empty() && (elements.front() < 0 || elements.back() >= num_elements())) {
    std::ostringstream os;
    os << "element set '" << name << "' references element "
       << (elements.front() < 0 ? elements.front() : elements.back()) << " outside [0, "
       << num_elements() - 1 << "]";
    throw MeshError(os.str());
  }
  element_sets_[name] = std::move(elements);
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
    } else if (type_ == ElementType::Tet4) {
      const Scalar volume = tet4_volume(element_coordinates(e));
      if (!(volume > 0.0)) {
        std::ostringstream os;
        os << "element " << e << " has signed volume " << volume
           << " m^3; Tet4 nodes must be ordered so that (x1 - x0) x (x2 - x0) . (x3 - x0) "
              "> 0 (swapping two nodes of an inverted cell fixes it) and the element "
              "must not be flat";
        throw MeshError(os.str());
      }
      min_measure = std::min(min_measure, volume);
      max_measure = std::max(max_measure, volume);
    } else if (type_ == ElementType::Tet10) {
      // The corners must be positively oriented, and the quadratic map must
      // not fold at an integration point: an edge node pushed far off its
      // edge can invert part of a cell whose corner volume is positive.
      const Matrix xe = element_coordinates(e);
      const Scalar corner_volume = tet4_volume(xe.leftCols(4));
      if (!(corner_volume > 0.0)) {
        std::ostringstream os;
        os << "element " << e << " has a corner-tetrahedron volume of " << corner_volume
           << " m^3; Tet10 corners must be ordered so that (x1 - x0) x (x2 - x0) . "
              "(x3 - x0) > 0 and the element must not be flat";
        throw MeshError(os.str());
      }
      Scalar min_det = 0.0;
      const Scalar volume = tet10_volume(xe, &min_det);
      if (!(min_det > 0.0) || !(volume > 0.0)) {
        std::ostringstream os;
        os << "element " << e << " has volume " << volume
           << " m^3 and a minimum Jacobian determinant of " << min_det
           << " m^3 over its four integration points; an edge node lies so far from "
              "the midpoint of its edge that the curved element folds. Move the edge "
              "nodes towards their edge midpoints (Gmsh: Mesh.HighOrderOptimize = 2) or "
              "refine where the boundary curves sharply";
        throw MeshError(os.str());
      }
      min_measure = std::min(min_measure, volume);
      max_measure = std::max(max_measure, volume);
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
  const std::size_t corners = static_cast<std::size_t>(face_corner_nodes(type_));

  // Faces are keyed by their sorted corner nodes, so orientation does not
  // matter (the corners alone identify a Tet10 face).
  // Sorting one flat array of keys is far cheaper than a map of vectors on a
  // large mesh, and the boundary comes out in ascending key order - the order
  // the earlier map-based search produced, which fixes the summation order of
  // traction loads and the triangle order of an STL file.
  struct Record {
    std::array<Index, 4> key;
    Index element;
    int local_face;
  };
  const Index ne = num_elements();
  // Sized up front and filled by index: besides saving the reallocations,
  // this keeps GCC from assuming a one-record buffer inside std::sort and
  // raising a spurious -Warray-bounds.
  std::vector<Record> records(static_cast<std::size_t>(ne) * local.size());
  std::size_t next = 0;
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = element_nodes(e);
    for (std::size_t lf = 0; lf < local.size(); ++lf) {
      Record& r = records[next++];
      r.key.fill(-1);
      for (std::size_t a = 0; a < corners; ++a) r.key[a] = nodes[local[lf][a]];
      std::sort(r.key.begin(), r.key.begin() + static_cast<long>(corners));
      r.element = e;
      r.local_face = static_cast<int>(lf);
    }
  }
  std::sort(records.begin(), records.end(), [](const Record& a, const Record& b) {
    if (a.key != b.key) return a.key < b.key;
    if (a.element != b.element) return a.element < b.element;
    return a.local_face < b.local_face;
  });

  std::vector<BoundaryFace> boundary;
  std::size_t i = 0;
  while (i < records.size()) {
    std::size_t j = i + 1;
    while (j < records.size() && records[j].key == records[i].key) ++j;
    if (j - i == 1) {
      BoundaryFace face;
      face.element = records[i].element;
      face.local_face = records[i].local_face;
      const Index* nodes = element_nodes(face.element);
      for (int a : local[static_cast<std::size_t>(face.local_face)]) {
        face.nodes.push_back(nodes[a]);
      }
      boundary.push_back(std::move(face));
    }
    i = j;
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
  if (face.nodes.size() == 3 && dim_ == 3) {
    const Vector3 a = coords_.col(face.nodes[0]);
    const Vector3 b = coords_.col(face.nodes[1]);
    const Vector3 c = coords_.col(face.nodes[2]);
    return 0.5 * (b - a).cross(c - a).norm();
  }
  if (face.nodes.size() == 6 && dim_ == 3) {
    Matrix xf(3, 6);
    for (int a = 0; a < 6; ++a) xf.col(a) = coords_.col(face.nodes[static_cast<std::size_t>(a)]);
    return tri6_face_area(xf);
  }
  throw MeshError("face_measure received a face with an unsupported node count");
}

}  // namespace sparlab
