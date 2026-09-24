#include "sparlab/io/StlWriter.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

namespace sparlab {
namespace {

void emit_quad(TriangleSurface& out, Index p0, Index p1, Index p2, Index p3) {
  // Two triangles with the same winding as the quad (p0 p1 p2 p3).
  out.faces.push_back({p0, p1, p2});
  out.faces.push_back({p0, p2, p3});
}

void write_float(std::ostream& out, float value) {
  out.write(reinterpret_cast<const char*>(&value), sizeof(float));
}

void write_vec(std::ostream& out, const Vector3& v) {
  write_float(out, static_cast<float>(v.x()));
  write_float(out, static_cast<float>(v.y()));
  write_float(out, static_cast<float>(v.z()));
}

float read_float(std::istream& in) {
  float value = 0.0f;
  in.read(reinterpret_cast<char*>(&value), sizeof(float));
  return value;
}

}  // namespace

Vector3 Triangle::normal() const {
  const Vector3 n = (b - a).cross(c - a);
  const Scalar length = n.norm();
  return length > 0.0 ? Vector3(n / length) : Vector3::Zero();
}

Scalar Triangle::area() const { return 0.5 * (b - a).cross(c - a).norm(); }

Triangle TriangleSurface::triangle(Index i) const {
  const std::array<Index, 3>& f = faces[static_cast<std::size_t>(i)];
  Triangle t;
  t.a = points[static_cast<std::size_t>(f[0])];
  t.b = points[static_cast<std::size_t>(f[1])];
  t.c = points[static_cast<std::size_t>(f[2])];
  return t;
}

std::vector<Triangle> TriangleSurface::triangles() const {
  std::vector<Triangle> out;
  out.reserve(faces.size());
  for (Index i = 0; i < num_triangles(); ++i) out.push_back(triangle(i));
  return out;
}

TriangleSurface boundary_surface(const Mesh& mesh, Scalar thickness) {
  TriangleSurface surface;
  const std::vector<Mesh::BoundaryFace> faces = mesh.boundary_faces();

  if (mesh.dim() == 3) {
    // Every mesh node becomes a surface point (unused interior points are
    // harmless to the STL, which stores triangles only) and the outward-wound
    // face tables can be used as they are.
    surface.points.reserve(static_cast<std::size_t>(mesh.num_nodes()));
    for (Index n = 0; n < mesh.num_nodes(); ++n) surface.points.push_back(mesh.node(n));
    surface.faces.reserve(2 * faces.size());
    for (const Mesh::BoundaryFace& f : faces) {
      if (f.nodes.size() != 4) throw IoError("expected quadrilateral boundary faces");
      emit_quad(surface, f.nodes[0], f.nodes[1], f.nodes[2], f.nodes[3]);
    }
    return surface;
  }

  if (!(thickness > 0.0)) {
    std::ostringstream os;
    os << "a 2-D mesh is extruded by the model thickness, which must be positive (got "
       << thickness << " m)";
    throw ConfigError(os.str());
  }
  const Index nn = mesh.num_nodes();
  const Index ne = mesh.num_elements();
  const int npe = mesh.nodes_per_elem();
  // Points 0..nn-1 are the bottom cap (z = 0), nn..2nn-1 the top (z = t).
  surface.points.reserve(static_cast<std::size_t>(2 * nn));
  for (Index n = 0; n < nn; ++n) surface.points.push_back(mesh.node(n));
  for (Index n = 0; n < nn; ++n) {
    surface.points.push_back(mesh.node(n) + Vector3(0.0, 0.0, thickness));
  }
  surface.faces.reserve(static_cast<std::size_t>(ne) * 4 + faces.size() * 2);

  // Caps: the counter-clockwise element polygon gives an outward (+z) normal
  // on the top; the bottom cap is wound the other way.
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 1; a + 1 < npe; ++a) {
      surface.faces.push_back({nodes[0] + nn, nodes[a] + nn, nodes[a + 1] + nn});
      surface.faces.push_back({nodes[0], nodes[a + 1], nodes[a]});
    }
  }
  // Side walls: for a counter-clockwise element the outward normal of edge
  // a -> b is (b - a) rotated clockwise, which the quad (a, b, b + t, a + t)
  // reproduces by the right-hand rule.
  for (const Mesh::BoundaryFace& f : faces) {
    if (f.nodes.size() != 2) throw IoError("expected two-node boundary edges");
    emit_quad(surface, f.nodes[0], f.nodes[1], f.nodes[1] + nn, f.nodes[0] + nn);
  }
  return surface;
}

SurfaceStats surface_stats(const TriangleSurface& surface) {
  SurfaceStats stats;
  stats.num_triangles = surface.num_triangles();
  stats.num_points = static_cast<Index>(surface.points.size());
  if (surface.faces.empty()) return stats;
  stats.bounds.lower = Vector3::Constant(std::numeric_limits<Scalar>::max());
  stats.bounds.upper = Vector3::Constant(-std::numeric_limits<Scalar>::max());

  // Directed edges: a closed, consistently oriented surface uses every edge
  // as often in one direction as in the other. A 2-manifold uses it exactly
  // once each way; cells touching only along an edge use it twice each way,
  // which is counted separately because the surface still closes.
  std::map<std::pair<Index, Index>, Index> edges;
  for (const std::array<Index, 3>& f : surface.faces) {
    for (int k = 0; k < 3; ++k) ++edges[{f[k], f[(k + 1) % 3]}];
  }
  Index unmatched = 0;
  Index non_manifold = 0;
  for (const auto& entry : edges) {
    const auto reverse = edges.find({entry.first.second, entry.first.first});
    const Index reversed = reverse == edges.end() ? 0 : reverse->second;
    if (entry.second != reversed) {
      unmatched += std::abs(entry.second - reversed);
    } else if (entry.second > 1 && entry.first.first < entry.first.second) {
      ++non_manifold;  // count each undirected edge once
    }
  }
  stats.unmatched_edges = unmatched;
  stats.non_manifold_edges = non_manifold;
  stats.closed = unmatched == 0;

  for (Index i = 0; i < surface.num_triangles(); ++i) {
    const Triangle t = surface.triangle(i);
    stats.area += t.area();
    // Divergence theorem: V = 1/6 sum a . (b x c) over outward triangles.
    stats.enclosed_volume += t.a.dot(t.b.cross(t.c)) / 6.0;
    for (const Vector3* p : {&t.a, &t.b, &t.c}) {
      stats.bounds.lower = stats.bounds.lower.cwiseMin(*p);
      stats.bounds.upper = stats.bounds.upper.cwiseMax(*p);
    }
  }
  return stats;
}

void write_stl(const std::string& path, const TriangleSurface& surface,
               const std::string& header) {
  if (surface.faces.size() >
      static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw IoError("too many triangles for the STL facet counter");
  }
  std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
  if (!out) throw IoError("cannot open '" + path + "' for writing");

  // A binary STL must not begin with "solid" or some readers take it for
  // ASCII; the header is padded to exactly 80 bytes.
  char head[80];
  std::memset(head, 0, sizeof(head));
  const std::string text = "SparLab binary STL: " + header;
  std::memcpy(head, text.data(), std::min<std::size_t>(text.size(), 79));
  out.write(head, sizeof(head));
  const std::uint32_t count = static_cast<std::uint32_t>(surface.faces.size());
  out.write(reinterpret_cast<const char*>(&count), sizeof(count));
  const std::uint16_t attribute = 0;
  for (Index i = 0; i < surface.num_triangles(); ++i) {
    const Triangle t = surface.triangle(i);
    write_vec(out, t.normal());
    write_vec(out, t.a);
    write_vec(out, t.b);
    write_vec(out, t.c);
    out.write(reinterpret_cast<const char*>(&attribute), sizeof(attribute));
  }
  out.flush();
  if (!out) throw IoError("failed while writing '" + path + "'");
}

std::vector<Triangle> read_stl(const std::string& path) {
  std::ifstream in(path, std::ios::in | std::ios::binary);
  if (!in) throw IoError("cannot open '" + path + "' for reading");
  char head[80];
  in.read(head, sizeof(head));
  std::uint32_t count = 0;
  in.read(reinterpret_cast<char*>(&count), sizeof(count));
  if (!in) throw IoError("'" + path + "' is too short to be a binary STL");
  std::vector<Triangle> triangles;
  triangles.reserve(count);
  for (std::uint32_t i = 0; i < count; ++i) {
    Triangle t;
    for (int k = 0; k < 3; ++k) (void)read_float(in);  // stored normal
    for (Vector3* p : {&t.a, &t.b, &t.c}) {
      const float x = read_float(in);
      const float y = read_float(in);
      const float z = read_float(in);
      *p = Vector3(x, y, z);
    }
    std::uint16_t attribute = 0;
    in.read(reinterpret_cast<char*>(&attribute), sizeof(attribute));
    if (!in) throw IoError("'" + path + "' ended before its declared facet count");
    triangles.push_back(t);
  }
  return triangles;
}

}  // namespace sparlab
