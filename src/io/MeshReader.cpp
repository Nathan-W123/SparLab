#include "sparlab/io/MeshReader.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/Hex8.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_map>
#include <utility>

namespace sparlab {
namespace {

// ---------------------------------------------------------------------------
// Element kinds shared by both formats
// ---------------------------------------------------------------------------

enum class Kind { Point, Line, Tri3, Quad4, Tet4, Hex8, Unsupported };

struct KindInfo {
  Kind kind = Kind::Unsupported;
  int dim = -1;
  int nodes = -1;         ///< expected node count; -1 when unknown
  std::string name;       ///< human-readable, used in messages
  std::string advice;     ///< how to avoid it, for unsupported cells
};

const char* kSecondOrderAdvice =
    "SparLab's elements are first order; re-export the mesh with linear elements "
    "(Gmsh: Mesh.ElementOrder = 1; Abaqus/CalculiX: C3D4/C3D8, CPS3/CPS4)";
const char* kPrismAdvice =
    "SparLab supports tetrahedra or hexahedra only; re-mesh the volume with "
    "tetrahedra (Gmsh: Mesh.RecombineAll = 0, no Recombine in extrusions)";

KindInfo make_kind(Kind kind, int dim, int nodes, std::string name, std::string advice = "") {
  KindInfo info;
  info.kind = kind;
  info.dim = dim;
  info.nodes = nodes;
  info.name = std::move(name);
  info.advice = std::move(advice);
  return info;
}

/// Gmsh element type codes (see the Gmsh reference manual, "MSH file format").
KindInfo gmsh_kind(int code) {
  switch (code) {
    case 1: return make_kind(Kind::Line, 1, 2, "2-node line");
    case 2: return make_kind(Kind::Tri3, 2, 3, "3-node triangle");
    case 3: return make_kind(Kind::Quad4, 2, 4, "4-node quadrilateral");
    case 4: return make_kind(Kind::Tet4, 3, 4, "4-node tetrahedron");
    case 5: return make_kind(Kind::Hex8, 3, 8, "8-node hexahedron");
    case 6: return make_kind(Kind::Unsupported, 3, 6, "6-node prism", kPrismAdvice);
    case 7: return make_kind(Kind::Unsupported, 3, 5, "5-node pyramid", kPrismAdvice);
    case 8: return make_kind(Kind::Line, 1, 3, "3-node line");
    case 9: return make_kind(Kind::Unsupported, 2, 6, "6-node triangle", kSecondOrderAdvice);
    case 10: return make_kind(Kind::Unsupported, 2, 9, "9-node quadrilateral", kSecondOrderAdvice);
    case 11: return make_kind(Kind::Unsupported, 3, 10, "10-node tetrahedron", kSecondOrderAdvice);
    case 12: return make_kind(Kind::Unsupported, 3, 27, "27-node hexahedron", kSecondOrderAdvice);
    case 13: return make_kind(Kind::Unsupported, 3, 18, "18-node prism", kPrismAdvice);
    case 14: return make_kind(Kind::Unsupported, 3, 14, "14-node pyramid", kPrismAdvice);
    case 15: return make_kind(Kind::Point, 0, 1, "point");
    case 16: return make_kind(Kind::Unsupported, 2, 8, "8-node quadrilateral", kSecondOrderAdvice);
    case 17: return make_kind(Kind::Unsupported, 3, 20, "20-node hexahedron", kSecondOrderAdvice);
    case 18: return make_kind(Kind::Unsupported, 3, 15, "15-node prism", kPrismAdvice);
    case 19: return make_kind(Kind::Unsupported, 3, 13, "13-node pyramid", kPrismAdvice);
    case 20: return make_kind(Kind::Unsupported, 2, 9, "9-node triangle", kSecondOrderAdvice);
    case 21: return make_kind(Kind::Unsupported, 2, 10, "10-node triangle", kSecondOrderAdvice);
    case 22: return make_kind(Kind::Unsupported, 2, 12, "12-node triangle", kSecondOrderAdvice);
    case 23: return make_kind(Kind::Unsupported, 2, 15, "15-node triangle", kSecondOrderAdvice);
    case 24: return make_kind(Kind::Unsupported, 2, 15, "15-node triangle", kSecondOrderAdvice);
    case 25: return make_kind(Kind::Unsupported, 2, 21, "21-node triangle", kSecondOrderAdvice);
    case 26: return make_kind(Kind::Line, 1, 4, "4-node line");
    case 27: return make_kind(Kind::Line, 1, 5, "5-node line");
    case 28: return make_kind(Kind::Line, 1, 6, "6-node line");
    case 29: return make_kind(Kind::Unsupported, 3, 20, "20-node tetrahedron", kSecondOrderAdvice);
    case 30: return make_kind(Kind::Unsupported, 3, 35, "35-node tetrahedron", kSecondOrderAdvice);
    case 31: return make_kind(Kind::Unsupported, 3, 56, "56-node tetrahedron", kSecondOrderAdvice);
    case 92: return make_kind(Kind::Unsupported, 3, 64, "64-node hexahedron", kSecondOrderAdvice);
    case 93: return make_kind(Kind::Unsupported, 3, 125, "125-node hexahedron", kSecondOrderAdvice);
    default: break;
  }
  std::ostringstream os;
  os << "Gmsh element type " << code << " is not recognised";
  throw IoError(os.str());
}

std::string upper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return s;
}

/// Abaqus / CalculiX element type names: a family prefix followed by the node
/// count and optional formulation letters (C3D8R, CPS4I, S4R, ...).
KindInfo abaqus_kind(const std::string& type_in) {
  const std::string type = upper(type_in);
  struct Family {
    const char* prefix;
    int dim;
    bool solid;  // true: C3D-like (dim 3), false: plane / surface / line
  };
  // Longest prefixes first so "SFM3D" is not read as "S".
  static const Family families[] = {
      {"CONN3D", 1, false}, {"CONN2D", 1, false}, {"SPRING", 1, false},
      {"DASHPOT", 1, false}, {"SFM3D", 2, false}, {"DC3D", 3, true},
      {"DC2D", 2, false},   {"RB3D", 1, false},   {"RB2D", 1, false},
      {"C3D", 3, true},     {"CPS", 2, false},    {"CPE", 2, false},
      {"CAX", 2, false},    {"CGAX", 2, false},   {"M3D", 2, false},
      {"R3D", 2, false},    {"R2D", 1, false},    {"T2D", 1, false},
      {"T3D", 1, false},    {"STRI", 2, false},   {"SC", 2, false},
      {"B2", 1, false},     {"B3", 1, false},     {"S", 2, false}};
  for (const Family& f : families) {
    const std::string prefix = f.prefix;
    if (type.compare(0, prefix.size(), prefix) != 0) continue;
    std::size_t k = prefix.size();
    std::string digits;
    while (k < type.size() && std::isdigit(static_cast<unsigned char>(type[k]))) {
      digits.push_back(type[k++]);
    }
    if (digits.empty()) continue;
    const int n = std::atoi(digits.c_str());
    if (f.dim == 1) return make_kind(Kind::Line, 1, -1, type + " line element");
    if (f.dim == 3) {
      if (n == 4) return make_kind(Kind::Tet4, 3, 4, type + " (4-node tetrahedron)");
      if (n == 8) return make_kind(Kind::Hex8, 3, 8, type + " (8-node hexahedron)");
      if (n == 6 || n == 15 || n == 5 || n == 13) {
        return make_kind(Kind::Unsupported, 3, n, type + " (prism / pyramid)", kPrismAdvice);
      }
      return make_kind(Kind::Unsupported, 3, n, type + " (second-order solid)",
                       kSecondOrderAdvice);
    }
    // Plane, axisymmetric, shell, membrane and rigid surface elements.
    if (n == 3) return make_kind(Kind::Tri3, 2, 3, type + " (3-node triangle)");
    if (n == 4) return make_kind(Kind::Quad4, 2, 4, type + " (4-node quadrilateral)");
    return make_kind(Kind::Unsupported, 2, n, type + " (second-order surface)",
                     kSecondOrderAdvice);
  }
  throw IoError("Abaqus element type '" + type_in + "' is not recognised");
}

std::string dimension_prefix(int dim) {
  switch (dim) {
    case 0: return "point:";
    case 1: return "curve:";
    case 2: return "surface:";
    default: return "volume:";
  }
}

// ---------------------------------------------------------------------------
// Format-independent intermediate mesh
// ---------------------------------------------------------------------------

struct RawElement {
  long long tag = 0;
  KindInfo info;
  std::vector<long long> nodes;
  std::vector<std::string> groups;
};

struct RawMesh {
  std::string source;
  std::vector<long long> node_tags;
  std::vector<Vector3> node_coords;
  std::vector<RawElement> elements;
  std::map<std::string, std::vector<long long>> node_set_tags;
  std::map<std::string, std::vector<long long>> element_set_tags;
  std::map<std::string, Index> ignored;
  std::vector<std::string> warnings;
};

void trim(std::string& s) {
  const auto not_space = [](unsigned char c) { return !std::isspace(c); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
  s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
}

std::vector<std::string> split_whitespace(const std::string& line) {
  std::vector<std::string> out;
  std::istringstream is(line);
  std::string token;
  while (is >> token) out.push_back(token);
  return out;
}

long long to_integer(const std::string& token, const std::string& where) {
  char* end = nullptr;
  const long long value = std::strtoll(token.c_str(), &end, 10);
  if (end == token.c_str() || *end != '\0') {
    throw IoError(where + ": expected an integer, found '" + token + "'");
  }
  return value;
}

Scalar to_real(const std::string& token, const std::string& where) {
  char* end = nullptr;
  const Scalar value = std::strtod(token.c_str(), &end);
  if (end == token.c_str() || *end != '\0' || !std::isfinite(value)) {
    throw IoError(where + ": expected a finite number, found '" + token + "'");
  }
  return value;
}

std::string context(const std::string& source, long long line) {
  std::ostringstream os;
  os << source << ":" << line;
  return os.str();
}

/// Turn the raw file content into a validated mesh with its named sets.
Mesh build_mesh(RawMesh& raw, const MeshReadOptions& options, MeshReadReport& report) {
  const std::string& src = raw.source;
  report.scale = options.scale;
  report.nodes_in_file = static_cast<Index>(raw.node_tags.size());
  report.elements_in_file = static_cast<Index>(raw.elements.size());
  report.ignored.insert(raw.ignored.begin(), raw.ignored.end());
  report.warnings.insert(report.warnings.end(), raw.warnings.begin(), raw.warnings.end());
  if (!(options.scale > 0.0) || !std::isfinite(options.scale)) {
    std::ostringstream os;
    os << "mesh scale must be a positive factor (got " << options.scale << ")";
    throw ConfigError(os.str());
  }

  if (raw.elements.empty()) {
    throw MeshError(src + " contains no elements");
  }
  int dim = -1;
  for (const RawElement& e : raw.elements) dim = std::max(dim, e.info.dim);
  if (dim < 2) {
    throw MeshError(src + " contains only points and lines; SparLab needs cells: "
                          "triangles or quadrilaterals for a plane model, tetrahedra or "
                          "hexahedra for a solid");
  }

  // The cells are the elements of the highest dimension, all of one type.
  std::map<std::string, Index> cell_types;
  for (const RawElement& e : raw.elements) {
    if (e.info.dim == dim) ++cell_types[e.info.name];
  }
  if (cell_types.size() > 1) {
    std::ostringstream os;
    os << src << " mixes " << cell_types.size() << " kinds of " << dim << "-D cell (";
    bool first = true;
    for (const auto& entry : cell_types) {
      os << (first ? "" : ", ") << entry.second << " x " << entry.first;
      first = false;
    }
    os << "); a SparLab mesh carries one cell type. ";
    if (dim == 2) {
      os << "Mesh it all-triangle (Gmsh: Mesh.RecombineAll = 0) or all-quadrilateral "
            "(Gmsh: Mesh.RecombineAll = 1 with Mesh.SubdivisionAlgorithm = 1)";
    } else {
      os << "Mesh the volume with tetrahedra only, or with hexahedra only (extruded or "
            "transfinite volumes)";
    }
    throw MeshError(os.str());
  }
  const RawElement* sample = nullptr;
  for (const RawElement& e : raw.elements) {
    if (e.info.dim == dim) {
      sample = &e;
      break;
    }
  }
  const KindInfo cell = sample->info;
  if (cell.kind == Kind::Unsupported) {
    throw MeshError(src + " is made of " + cell.name + " cells, which SparLab does not "
                    "implement. " + cell.advice);
  }
  ElementType type = ElementType::Tri3;
  switch (cell.kind) {
    case Kind::Tri3: type = ElementType::Tri3; break;
    case Kind::Quad4: type = ElementType::Quad4; break;
    case Kind::Tet4: type = ElementType::Tet4; break;
    case Kind::Hex8: type = ElementType::Hex8; break;
    default: throw MeshError(src + ": unexpected cell kind");
  }
  const int npe = nodes_per_element(type);
  report.dimension = dim;
  report.element_type = type;

  // Node tags -> file positions.
  std::unordered_map<long long, Index> tag_to_file;
  tag_to_file.reserve(raw.node_tags.size() * 2);
  for (std::size_t i = 0; i < raw.node_tags.size(); ++i) {
    if (!tag_to_file.emplace(raw.node_tags[i], static_cast<Index>(i)).second) {
      std::ostringstream os;
      os << src << " defines node tag " << raw.node_tags[i] << " twice";
      throw IoError(os.str());
    }
  }
  const auto file_index = [&](long long tag, long long element_tag) {
    const auto it = tag_to_file.find(tag);
    if (it == tag_to_file.end()) {
      std::ostringstream os;
      os << src << ": element " << element_tag << " references node " << tag
         << ", which the file does not define";
      throw IoError(os.str());
    }
    return it->second;
  };

  // Cells, in file order.
  std::vector<std::array<Index, 8>> cells;
  std::vector<long long> cell_tags;
  std::vector<char> used(raw.node_tags.size(), 0);
  for (const RawElement& e : raw.elements) {
    if (e.info.dim != dim) continue;
    if (static_cast<int>(e.nodes.size()) != npe) {
      std::ostringstream os;
      os << src << ": element " << e.tag << " (" << e.info.name << ") lists "
         << e.nodes.size() << " nodes, expected " << npe;
      throw IoError(os.str());
    }
    std::array<Index, 8> c{};
    for (int a = 0; a < npe; ++a) {
      c[static_cast<std::size_t>(a)] = file_index(e.nodes[static_cast<std::size_t>(a)], e.tag);
      used[static_cast<std::size_t>(c[static_cast<std::size_t>(a)])] = 1;
    }
    cells.push_back(c);
    cell_tags.push_back(e.tag);
  }
  report.cells = static_cast<Index>(cells.size());
  report.boundary_elements = report.elements_in_file - report.cells;

  // Compact node numbering in file order, keeping only referenced nodes.
  std::vector<Index> file_to_mesh(raw.node_tags.size(), -1);
  Index next = 0;
  for (std::size_t i = 0; i < used.size(); ++i) {
    if (used[i]) file_to_mesh[i] = next++;
  }
  report.unreferenced_nodes = report.nodes_in_file - next;

  Matrix coords(dim, next);
  for (std::size_t i = 0; i < used.size(); ++i) {
    if (file_to_mesh[i] < 0) continue;
    const Vector3 x = options.scale * raw.node_coords[i];
    coords.col(file_to_mesh[i]) = x.head(dim);
  }
  if (!coords.allFinite()) throw MeshError(src + ": node coordinates are not finite");

  Vector3 lo = Vector3::Constant(std::numeric_limits<Scalar>::max());
  Vector3 hi = Vector3::Constant(-std::numeric_limits<Scalar>::max());
  for (std::size_t i = 0; i < used.size(); ++i) {
    if (!used[i]) continue;
    const Vector3 x = options.scale * raw.node_coords[i];
    lo = lo.cwiseMin(x);
    hi = hi.cwiseMax(x);
  }
  const Scalar diagonal = (hi - lo).norm();
  if (dim == 2) {
    const Scalar zspan = std::max(std::abs(lo.z()), std::abs(hi.z()));
    if (zspan > 1.0e-9 * std::max(diagonal, 1.0e-300)) {
      std::ostringstream os;
      os << src << " is a 2-D mesh but its nodes reach z = " << (std::abs(lo.z()) > std::abs(hi.z()) ? lo.z() : hi.z())
         << " m; a plane model must lie in the z = 0 plane. Draw the geometry in the x-y "
            "plane, or mesh the part as a solid";
      throw MeshError(os.str());
    }
  }
  if (diagonal > 20.0) {
    std::ostringstream os;
    os << "the mesh spans " << diagonal << " m corner to corner; if the geometry was drawn "
          "in millimetres, set mesh.scale = 0.001";
    report.warnings.push_back(os.str());
  } else if (diagonal < 1.0e-4) {
    std::ostringstream os;
    os << "the mesh spans only " << diagonal << " m corner to corner; check mesh.scale";
    report.warnings.push_back(os.str());
  }

  // Coincident nodes.
  const Scalar tol =
      options.duplicate_tolerance > 0.0 ? options.duplicate_tolerance : 1.0e-9 * diagonal;
  std::vector<Index> representative(static_cast<std::size_t>(next));
  std::iota(representative.begin(), representative.end(), Index{0});
  if (tol > 0.0 && next > 1) {
    std::unordered_map<std::size_t, std::vector<Index>> buckets;
    buckets.reserve(static_cast<std::size_t>(next) * 2);
    const auto cell_of = [&](const Vector3& x, int k) {
      return static_cast<long long>(std::floor(x(k) / tol));
    };
    const auto key = [](long long i, long long j, long long k) {
      std::size_t h = static_cast<std::size_t>(i) * 73856093u;
      h ^= static_cast<std::size_t>(j) * 19349663u;
      h ^= static_cast<std::size_t>(k) * 83492791u;
      return h;
    };
    for (Index n = 0; n < next; ++n) {
      Vector3 x = Vector3::Zero();
      x.head(dim) = coords.col(n);
      const long long ci = cell_of(x, 0);
      const long long cj = cell_of(x, 1);
      const long long ck = cell_of(x, 2);
      Index found = -1;
      for (long long di = -1; di <= 1 && found < 0; ++di) {
        for (long long dj = -1; dj <= 1 && found < 0; ++dj) {
          for (long long dk = -1; dk <= 1 && found < 0; ++dk) {
            const auto it = buckets.find(key(ci + di, cj + dj, ck + dk));
            if (it == buckets.end()) continue;
            for (Index m : it->second) {
              Vector3 y = Vector3::Zero();
              y.head(dim) = coords.col(m);
              if ((x - y).norm() <= tol) {
                found = m;
                break;
              }
            }
          }
        }
      }
      if (found >= 0) {
        representative[static_cast<std::size_t>(n)] = representative[static_cast<std::size_t>(found)];
        ++report.duplicate_nodes;
      }
      buckets[key(ci, cj, ck)].push_back(n);
    }
  }
  if (report.duplicate_nodes > 0) {
    std::ostringstream os;
    os << report.duplicate_nodes << " node(s) coincide with another node within " << tol
       << " m; cells on either side of such a seam are not connected";
    if (options.merge_duplicate_nodes) {
      os << ". They were merged (mesh.merge_duplicate_nodes = true)";
      report.duplicates_merged = true;
    } else {
      os << ". Merge them in the mesher (Gmsh: Coherence Mesh) or set "
            "mesh.merge_duplicate_nodes = true";
    }
    report.warnings.push_back(os.str());
  }

  // Connectivity with orientation repair.
  std::vector<Index> connectivity;
  connectivity.reserve(cells.size() * static_cast<std::size_t>(npe));
  Matrix xe(dim, npe);
  for (std::size_t c = 0; c < cells.size(); ++c) {
    std::array<Index, 8> n{};
    for (int a = 0; a < npe; ++a) {
      Index m = file_to_mesh[static_cast<std::size_t>(cells[c][static_cast<std::size_t>(a)])];
      if (options.merge_duplicate_nodes) m = representative[static_cast<std::size_t>(m)];
      n[static_cast<std::size_t>(a)] = m;
    }
    for (int a = 0; a < npe; ++a) xe.col(a) = coords.col(n[static_cast<std::size_t>(a)]);
    bool flipped = false;
    if (type == ElementType::Tri3) {
      const Scalar twice = (xe(0, 1) - xe(0, 0)) * (xe(1, 2) - xe(1, 0)) -
                           (xe(0, 2) - xe(0, 0)) * (xe(1, 1) - xe(1, 0));
      if (twice < 0.0) {
        std::swap(n[1], n[2]);
        flipped = true;
      }
    } else if (type == ElementType::Quad4) {
      Scalar twice = 0.0;
      for (int a = 0; a < 4; ++a) {
        const int b = (a + 1) % 4;
        twice += xe(0, a) * xe(1, b) - xe(0, b) * xe(1, a);
      }
      if (twice < 0.0) {
        std::swap(n[1], n[3]);
        flipped = true;
      }
    } else if (type == ElementType::Tet4) {
      const Vector3 a = xe.col(1) - xe.col(0);
      const Vector3 b = xe.col(2) - xe.col(0);
      const Vector3 d = xe.col(3) - xe.col(0);
      if (a.dot(b.cross(d)) < 0.0) {
        std::swap(n[1], n[2]);
        flipped = true;
      }
    } else {
      Scalar min_det = 0.0;
      const Scalar volume = hex8_volume(xe, &min_det);
      if (volume < 0.0) {
        // A mirrored hexahedron lists both faces clockwise; reversing each
        // face (0 3 2 1 / 4 7 6 5) restores the VTK ordering.
        std::swap(n[1], n[3]);
        std::swap(n[5], n[7]);
        flipped = true;
      }
    }
    if (flipped) ++report.reoriented;
    for (int a = 0; a < npe; ++a) connectivity.push_back(n[static_cast<std::size_t>(a)]);
  }
  if (report.reoriented > 0) {
    std::ostringstream os;
    os << report.reoriented << " of " << cells.size()
       << " cell(s) were listed with the opposite orientation and have been re-ordered";
    log::info(src, ": ", os.str());
  }

  // Merging can leave nodes unreferenced and cells collapsed.
  if (options.merge_duplicate_nodes && report.duplicate_nodes > 0) {
    Index collapsed = 0;
    for (std::size_t c = 0; c < cells.size(); ++c) {
      std::set<Index> distinct(connectivity.begin() + static_cast<long>(c) * npe,
                               connectivity.begin() + static_cast<long>(c + 1) * npe);
      if (static_cast<int>(distinct.size()) < npe) ++collapsed;
    }
    if (collapsed > 0) {
      std::ostringstream os;
      os << src << ": merging coincident nodes within " << tol << " m collapsed "
         << collapsed << " cell(s); the tolerance is larger than the smallest cell. "
            "Lower mesh.duplicate_tolerance";
      throw MeshError(os.str());
    }
    std::vector<Index> remap(static_cast<std::size_t>(next), -1);
    Index kept = 0;
    for (Index node : connectivity) {
      if (remap[static_cast<std::size_t>(node)] < 0) remap[static_cast<std::size_t>(node)] = 0;
    }
    for (Index n = 0; n < next; ++n) {
      if (remap[static_cast<std::size_t>(n)] == 0) remap[static_cast<std::size_t>(n)] = kept++;
    }
    Matrix merged(dim, kept);
    for (Index n = 0; n < next; ++n) {
      if (remap[static_cast<std::size_t>(n)] >= 0) merged.col(remap[static_cast<std::size_t>(n)]) = coords.col(n);
    }
    for (Index& node : connectivity) node = remap[static_cast<std::size_t>(node)];
    for (Index& m : file_to_mesh) {
      if (m >= 0) m = remap[static_cast<std::size_t>(representative[static_cast<std::size_t>(m)])];
    }
    coords = std::move(merged);
    next = kept;
  }
  report.nodes_used = next;

  Mesh mesh(std::move(coords), std::move(connectivity), type);

  // Named sets.
  std::map<std::string, std::set<Index>> node_sets;
  std::map<std::string, std::set<Index>> element_sets;
  std::unordered_map<long long, Index> cell_of_tag;
  for (std::size_t c = 0; c < cell_tags.size(); ++c) {
    cell_of_tag.emplace(cell_tags[c], static_cast<Index>(c));
  }
  Index dropped_set_nodes = 0;
  const auto add_nodes = [&](const std::string& name, const std::vector<long long>& tags,
                             long long element_tag) {
    std::set<Index>& target = node_sets[name];
    for (long long tag : tags) {
      const Index f = file_index(tag, element_tag);
      const Index m = file_to_mesh[static_cast<std::size_t>(f)];
      if (m >= 0) {
        target.insert(m);
      } else {
        ++dropped_set_nodes;
      }
    }
  };
  {
    Index c = 0;
    for (const RawElement& e : raw.elements) {
      const bool is_cell = e.info.dim == dim;
      for (const std::string& g : e.groups) {
        add_nodes(g, e.nodes, e.tag);
        if (is_cell) element_sets[g].insert(c);
      }
      if (is_cell) ++c;
    }
  }
  std::unordered_map<long long, const RawElement*> element_of_tag;
  for (const RawElement& e : raw.elements) element_of_tag.emplace(e.tag, &e);
  for (const auto& entry : raw.node_set_tags) {
    std::set<Index>& target = node_sets[entry.first];
    for (long long tag : entry.second) {
      const auto it = tag_to_file.find(tag);
      if (it == tag_to_file.end()) {
        std::ostringstream os;
        os << src << ": node set '" << entry.first << "' lists node " << tag
           << ", which the file does not define";
        throw IoError(os.str());
      }
      const Index m = file_to_mesh[static_cast<std::size_t>(it->second)];
      if (m >= 0) {
        target.insert(m);
      } else {
        ++dropped_set_nodes;
      }
    }
  }
  for (const auto& entry : raw.element_set_tags) {
    for (long long tag : entry.second) {
      const auto cell_it = cell_of_tag.find(tag);
      if (cell_it != cell_of_tag.end()) {
        element_sets[entry.first].insert(cell_it->second);
        add_nodes(entry.first, element_of_tag.at(tag)->nodes, tag);
        continue;
      }
      const auto el = element_of_tag.find(tag);
      if (el == element_of_tag.end()) {
        std::ostringstream os;
        os << src << ": element set '" << entry.first << "' lists element " << tag
           << ", which the file does not define";
        throw IoError(os.str());
      }
      add_nodes(entry.first, el->second->nodes, tag);
    }
  }
  for (const auto& entry : node_sets) {
    mesh.set_node_set(entry.first, std::vector<Index>(entry.second.begin(), entry.second.end()));
    report.node_sets[entry.first] = static_cast<Index>(entry.second.size());
    if (entry.second.empty()) {
      report.warnings.push_back("group '" + entry.first +
                                "' contains no node of the mesh cells and selects nothing");
    }
  }
  for (const auto& entry : element_sets) {
    mesh.set_element_set(entry.first,
                         std::vector<Index>(entry.second.begin(), entry.second.end()));
    report.element_sets[entry.first] = static_cast<Index>(entry.second.size());
  }
  if (dropped_set_nodes > 0) {
    std::ostringstream os;
    os << dropped_set_nodes
       << " group membership(s) referred to nodes no cell uses (for example a geometry "
          "point) and were dropped";
    report.warnings.push_back(os.str());
  }

  mesh.validate();
  report.quality = mesh.quality();
  if (report.quality.min < 0.05) {
    std::ostringstream os;
    os << "the worst cell (element " << report.quality.worst_element << ") has quality "
       << report.quality.min << " (" << report.quality.metric << ", 1 = ideal); "
       << report.quality.poor_elements << " cell(s) fall below "
       << report.quality.poor_threshold
       << ". Results in and near them are unreliable: refine or smooth the mesh there";
    report.warnings.push_back(os.str());
  }
  for (const std::string& w : report.warnings) log::warn(src, ": ", w);
  log::info("read ", src, " (", report.format, (report.version.empty() ? "" : " "),
            report.version, "): ", mesh.num_nodes(), " nodes, ", mesh.num_elements(), " ",
            to_string(type), " cells, ", report.boundary_elements,
            " boundary element(s), ", report.node_sets.size(), " node set(s), ",
            report.element_sets.size(), " element set(s); quality min ",
            report.quality.min, ", mean ", report.quality.mean);
  return mesh;
}

// ---------------------------------------------------------------------------
// Gmsh
// ---------------------------------------------------------------------------

class LineReader {
 public:
  LineReader(std::istream& in, std::string source) : in_(in), source_(std::move(source)) {}

  /// Next non-empty line (trimmed); false at end of input.
  bool next(std::string& line) {
    while (std::getline(in_, line)) {
      ++number_;
      trim(line);
      if (!line.empty()) return true;
    }
    return false;
  }

  std::string require(const std::string& what) {
    std::string line;
    if (!next(line)) throw IoError(source_ + ": unexpected end of file while reading " + what);
    return line;
  }

  std::string where() const { return context(source_, number_); }
  const std::string& source() const { return source_; }

 private:
  std::istream& in_;
  std::string source_;
  long long number_ = 0;
};

void expect_end(LineReader& reader, const std::string& section) {
  const std::string line = reader.require("$End" + section);
  if (line != "$End" + section) {
    throw IoError(reader.where() + ": expected $End" + section + ", found '" + line + "'");
  }
}

void skip_section(LineReader& reader, const std::string& section) {
  std::string line;
  while (reader.next(line)) {
    if (line == "$End" + section) return;
  }
  throw IoError(reader.source() + ": section $" + section + " is not terminated");
}

}  // namespace

Mesh read_gmsh(std::istream& in, const std::string& source, const MeshReadOptions& options,
               MeshReadReport* report_out) {
  MeshReadReport local;
  MeshReadReport& report = report_out ? *report_out : local;
  report = MeshReadReport();
  report.path = source;
  report.format = "gmsh";

  RawMesh raw;
  raw.source = source;
  LineReader reader(in, source);

  std::string version;
  std::map<std::pair<int, long long>, std::string> physical_names;
  std::map<std::pair<int, long long>, std::vector<long long>> entity_physicals;
  bool have_nodes = false;
  bool have_elements = false;

  const auto group_name = [&](int dim, long long tag) {
    const auto it = physical_names.find({dim, tag});
    if (it != physical_names.end()) return it->second;
    return dimension_prefix(dim) + std::to_string(tag);
  };

  std::string line;
  while (reader.next(line)) {
    if (line == "$MeshFormat") {
      const std::vector<std::string> t = split_whitespace(reader.require("the format line"));
      if (t.size() < 3) throw IoError(reader.where() + ": malformed $MeshFormat line");
      version = t[0];
      if (t[1] != "0") {
        throw IoError(source + " is a binary Gmsh file; export it as ASCII "
                               "(Gmsh: Mesh.Binary = 0)");
      }
      if (version != "2.2" && version.rfind("4.1", 0) != 0 && version != "2.1" &&
          version != "2.0") {
        throw IoError(source + " uses MSH format version " + version +
                      "; SparLab reads versions 2.2 and 4.1 (Gmsh: Mesh.MshFileVersion = "
                      "4.1 or 2.2)");
      }
      expect_end(reader, "MeshFormat");
    } else if (line == "$PhysicalNames") {
      const long long n = to_integer(reader.require("the physical-name count"), reader.where());
      for (long long i = 0; i < n; ++i) {
        const std::string entry = reader.require("a physical name");
        std::istringstream is(entry);
        int dim = 0;
        long long tag = 0;
        if (!(is >> dim >> tag)) throw IoError(reader.where() + ": malformed physical name");
        std::string rest;
        std::getline(is, rest);
        trim(rest);
        if (rest.size() >= 2 && rest.front() == '"' && rest.back() == '"') {
          rest = rest.substr(1, rest.size() - 2);
        }
        if (rest.empty()) rest = dimension_prefix(dim) + std::to_string(tag);
        physical_names[{dim, tag}] = rest;
      }
      expect_end(reader, "PhysicalNames");
    } else if (line == "$Entities") {
      if (version.rfind("4", 0) != 0) {
        skip_section(reader, "Entities");
        continue;
      }
      const std::vector<std::string> counts = split_whitespace(reader.require("entity counts"));
      if (counts.size() < 4) throw IoError(reader.where() + ": malformed $Entities header");
      for (int dim = 0; dim <= 3; ++dim) {
        const long long n = to_integer(counts[static_cast<std::size_t>(dim)], reader.where());
        for (long long i = 0; i < n; ++i) {
          const std::vector<std::string> t = split_whitespace(reader.require("an entity"));
          // points: tag x y z nphys phys...; others: tag 6 bbox values nphys phys...
          const std::size_t phys_at = dim == 0 ? 4 : 7;
          if (t.size() <= phys_at) throw IoError(reader.where() + ": malformed entity line");
          const long long tag = to_integer(t[0], reader.where());
          const long long nphys = to_integer(t[phys_at], reader.where());
          if (t.size() < phys_at + 1 + static_cast<std::size_t>(nphys)) {
            throw IoError(reader.where() + ": entity lists fewer physical tags than declared");
          }
          std::vector<long long> phys;
          for (long long k = 0; k < nphys; ++k) {
            phys.push_back(std::llabs(to_integer(t[phys_at + 1 + static_cast<std::size_t>(k)], reader.where())));
          }
          entity_physicals[{dim, tag}] = phys;
        }
      }
      expect_end(reader, "Entities");
    } else if (line == "$Nodes") {
      if (version.empty()) throw IoError(source + ": $Nodes appears before $MeshFormat");
      have_nodes = true;
      if (version.rfind("4", 0) == 0) {
        const std::vector<std::string> h = split_whitespace(reader.require("the node header"));
        if (h.size() < 4) throw IoError(reader.where() + ": malformed $Nodes header");
        const long long blocks = to_integer(h[0], reader.where());
        const long long total = to_integer(h[1], reader.where());
        raw.node_tags.reserve(static_cast<std::size_t>(total));
        raw.node_coords.reserve(static_cast<std::size_t>(total));
        for (long long b = 0; b < blocks; ++b) {
          const std::vector<std::string> bh = split_whitespace(reader.require("a node block"));
          if (bh.size() < 4) throw IoError(reader.where() + ": malformed node block header");
          const int edim = static_cast<int>(to_integer(bh[0], reader.where()));
          const bool parametric = to_integer(bh[2], reader.where()) != 0;
          const long long count = to_integer(bh[3], reader.where());
          const std::size_t first = raw.node_tags.size();
          for (long long k = 0; k < count; ++k) {
            raw.node_tags.push_back(to_integer(reader.require("a node tag"), reader.where()));
          }
          for (long long k = 0; k < count; ++k) {
            const std::vector<std::string> c = split_whitespace(reader.require("node coordinates"));
            const std::size_t need = 3 + (parametric ? static_cast<std::size_t>(std::max(edim, 0)) : 0);
            if (c.size() < need) throw IoError(reader.where() + ": malformed node coordinates");
            raw.node_coords.emplace_back(to_real(c[0], reader.where()),
                                         to_real(c[1], reader.where()),
                                         to_real(c[2], reader.where()));
          }
          (void)first;
        }
        if (static_cast<long long>(raw.node_tags.size()) != total) {
          throw IoError(source + ": $Nodes declares " + std::to_string(total) +
                        " nodes but its blocks hold " + std::to_string(raw.node_tags.size()));
        }
      } else {
        const long long total = to_integer(reader.require("the node count"), reader.where());
        raw.node_tags.reserve(static_cast<std::size_t>(total));
        raw.node_coords.reserve(static_cast<std::size_t>(total));
        for (long long k = 0; k < total; ++k) {
          const std::vector<std::string> c = split_whitespace(reader.require("a node"));
          if (c.size() < 4) throw IoError(reader.where() + ": malformed node line");
          raw.node_tags.push_back(to_integer(c[0], reader.where()));
          raw.node_coords.emplace_back(to_real(c[1], reader.where()),
                                       to_real(c[2], reader.where()),
                                       to_real(c[3], reader.where()));
        }
      }
      expect_end(reader, "Nodes");
    } else if (line == "$Elements") {
      if (version.empty()) throw IoError(source + ": $Elements appears before $MeshFormat");
      have_elements = true;
      if (version.rfind("4", 0) == 0) {
        const std::vector<std::string> h = split_whitespace(reader.require("the element header"));
        if (h.size() < 4) throw IoError(reader.where() + ": malformed $Elements header");
        const long long blocks = to_integer(h[0], reader.where());
        const long long total = to_integer(h[1], reader.where());
        raw.elements.reserve(static_cast<std::size_t>(total));
        for (long long b = 0; b < blocks; ++b) {
          const std::vector<std::string> bh = split_whitespace(reader.require("an element block"));
          if (bh.size() < 4) throw IoError(reader.where() + ": malformed element block header");
          const int edim = static_cast<int>(to_integer(bh[0], reader.where()));
          const long long etag = to_integer(bh[1], reader.where());
          const KindInfo info = gmsh_kind(static_cast<int>(to_integer(bh[2], reader.where())));
          const long long count = to_integer(bh[3], reader.where());
          std::vector<std::string> groups;
          const auto phys = entity_physicals.find({edim, etag});
          if (phys != entity_physicals.end()) {
            for (long long p : phys->second) groups.push_back(group_name(edim, p));
          }
          for (long long k = 0; k < count; ++k) {
            const std::vector<std::string> t = split_whitespace(reader.require("an element"));
            if (t.size() != static_cast<std::size_t>(1 + info.nodes)) {
              throw IoError(reader.where() + ": " + info.name + " element lists " +
                            std::to_string(t.size() - 1) + " nodes");
            }
            RawElement e;
            e.tag = to_integer(t[0], reader.where());
            e.info = info;
            e.groups = groups;
            e.nodes.reserve(t.size() - 1);
            for (std::size_t a = 1; a < t.size(); ++a) {
              e.nodes.push_back(to_integer(t[a], reader.where()));
            }
            raw.elements.push_back(std::move(e));
          }
        }
      } else {
        const long long total = to_integer(reader.require("the element count"), reader.where());
        raw.elements.reserve(static_cast<std::size_t>(total));
        for (long long k = 0; k < total; ++k) {
          const std::vector<std::string> t = split_whitespace(reader.require("an element"));
          if (t.size() < 3) throw IoError(reader.where() + ": malformed element line");
          RawElement e;
          e.tag = to_integer(t[0], reader.where());
          e.info = gmsh_kind(static_cast<int>(to_integer(t[1], reader.where())));
          const std::size_t ntags = static_cast<std::size_t>(to_integer(t[2], reader.where()));
          if (t.size() != 3 + ntags + static_cast<std::size_t>(e.info.nodes)) {
            throw IoError(reader.where() + ": " + e.info.name + " element has " +
                          std::to_string(t.size()) + " fields, expected " +
                          std::to_string(3 + ntags + static_cast<std::size_t>(e.info.nodes)));
          }
          if (ntags >= 1) {
            const long long phys = to_integer(t[3], reader.where());
            if (phys != 0) e.groups.push_back(group_name(e.info.dim, phys));
          }
          for (std::size_t a = 3 + ntags; a < t.size(); ++a) {
            e.nodes.push_back(to_integer(t[a], reader.where()));
          }
          raw.elements.push_back(std::move(e));
        }
      }
      expect_end(reader, "Elements");
    } else if (!line.empty() && line.front() == '$') {
      const std::string section = line.substr(1);
      ++raw.ignored["$" + section + " section"];
      skip_section(reader, section);
    } else {
      throw IoError(reader.where() + ": unexpected content '" + line.substr(0, 40) +
                    "' outside a $Section; is this a Gmsh .msh file?");
    }
  }

  if (version.empty()) {
    throw IoError(source + " has no $MeshFormat section; it is not a Gmsh mesh file");
  }
  if (!have_nodes || !have_elements) {
    throw IoError(source + " lacks a $Nodes or $Elements section");
  }
  report.version = version;
  // Gmsh saves only the elements of physical groups once any group is defined
  // (unless Mesh.SaveAll = 1), which silently loses the volume when only the
  // boundary groups were declared.
  {
    int top = -1;
    for (const RawElement& e : raw.elements) top = std::max(top, e.info.dim);
    bool any_group = false;
    for (const RawElement& e : raw.elements) any_group = any_group || !e.groups.empty();
    if (top >= 0 && top < 2 && any_group) {
      throw MeshError(source + " holds only boundary elements: Gmsh writes only the "
                               "elements of physical groups once any group is defined. "
                               "Add a physical surface (2-D) or volume (3-D), or set "
                               "Mesh.SaveAll = 1");
    }
  }
  return build_mesh(raw, options, report);
}

// ---------------------------------------------------------------------------
// Abaqus / CalculiX
// ---------------------------------------------------------------------------

namespace {

struct InpState {
  enum class Block { None, Node, Element, Nset, Elset, Ignored } block = Block::None;
  std::string keyword;  // current keyword, upper case
  std::string set_name;
  bool generate = false;
  KindInfo element_kind;
  std::string pending;  // element definition continued on the next line
  long long pending_line = 0;
};

std::vector<std::string> split_commas(const std::string& line) {
  std::vector<std::string> out;
  std::string item;
  std::istringstream is(line);
  while (std::getline(is, item, ',')) {
    trim(item);
    out.push_back(item);
  }
  if (!line.empty() && line.back() == ',') out.push_back("");
  return out;
}

std::map<std::string, std::string> keyword_parameters(const std::vector<std::string>& parts) {
  std::map<std::string, std::string> params;
  for (std::size_t i = 1; i < parts.size(); ++i) {
    const std::string& p = parts[i];
    if (p.empty()) continue;
    const std::size_t eq = p.find('=');
    std::string key = upper(p.substr(0, eq));
    trim(key);
    std::string value = eq == std::string::npos ? "" : p.substr(eq + 1);
    trim(value);
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
      value = value.substr(1, value.size() - 2);
    }
    params[key] = value;
  }
  return params;
}

bool is_number(const std::string& token) {
  if (token.empty()) return false;
  char* end = nullptr;
  std::strtod(token.c_str(), &end);
  return end != token.c_str() && *end == '\0';
}

void parse_inp(std::istream& in, const std::string& source, const std::string& base_dir,
               RawMesh& raw, int depth);

void finish_element(InpState& state, RawMesh& raw, const std::string& source) {
  if (state.pending.empty()) return;
  const std::vector<std::string> t = split_commas(state.pending);
  const std::string where = context(source, state.pending_line);
  state.pending.clear();
  std::vector<std::string> fields;
  for (const std::string& f : t) {
    if (!f.empty()) fields.push_back(f);
  }
  if (fields.size() < 2) throw IoError(where + ": element line has no nodes");
  RawElement e;
  e.tag = to_integer(fields[0], where);
  e.info = state.element_kind;
  for (std::size_t a = 1; a < fields.size(); ++a) e.nodes.push_back(to_integer(fields[a], where));
  if (e.info.nodes > 0 && static_cast<int>(e.nodes.size()) != e.info.nodes) {
    throw IoError(where + ": " + e.info.name + " element " + std::to_string(e.tag) + " lists " +
                  std::to_string(e.nodes.size()) + " nodes, expected " +
                  std::to_string(e.info.nodes));
  }
  if (!state.set_name.empty()) raw.element_set_tags[state.set_name].push_back(e.tag);
  raw.elements.push_back(std::move(e));
}

void parse_inp(std::istream& in, const std::string& source, const std::string& base_dir,
               RawMesh& raw, int depth) {
  if (depth > 8) throw IoError(source + ": *INCLUDE nesting deeper than 8 levels");
  InpState state;
  std::string line;
  long long number = 0;
  while (std::getline(in, line)) {
    ++number;
    trim(line);
    if (line.empty() || line.rfind("**", 0) == 0) continue;
    const std::string where = context(source, number);
    if (line.front() == '*') {
      finish_element(state, raw, source);
      const std::vector<std::string> parts = split_commas(line);
      std::string keyword = upper(parts[0]);
      trim(keyword);
      state.keyword = keyword;
      const std::map<std::string, std::string> params = keyword_parameters(parts);
      state.generate = params.count("GENERATE") > 0;
      state.set_name.clear();
      if (keyword == "*NODE") {
        state.block = InpState::Block::Node;
        if (params.count("NSET")) state.set_name = params.at("NSET");
      } else if (keyword == "*ELEMENT") {
        if (!params.count("TYPE")) throw IoError(where + ": *ELEMENT without TYPE=");
        state.block = InpState::Block::Element;
        state.element_kind = abaqus_kind(params.at("TYPE"));
        if (params.count("ELSET")) state.set_name = params.at("ELSET");
      } else if (keyword == "*NSET") {
        if (!params.count("NSET")) throw IoError(where + ": *NSET without NSET=");
        state.block = InpState::Block::Nset;
        state.set_name = params.at("NSET");
        raw.node_set_tags[state.set_name];
      } else if (keyword == "*ELSET") {
        if (!params.count("ELSET")) throw IoError(where + ": *ELSET without ELSET=");
        state.block = InpState::Block::Elset;
        state.set_name = params.at("ELSET");
        raw.element_set_tags[state.set_name];
      } else if (keyword == "*INCLUDE") {
        if (!params.count("INPUT")) throw IoError(where + ": *INCLUDE without INPUT=");
        std::string path = params.at("INPUT");
        if (!path.empty() && path.front() != '/') path = base_dir + "/" + path;
        std::ifstream inc(path);
        if (!inc) throw IoError(where + ": cannot open included file '" + path + "'");
        const std::size_t slash = path.find_last_of('/');
        parse_inp(inc, path, slash == std::string::npos ? "." : path.substr(0, slash), raw,
                  depth + 1);
        state.block = InpState::Block::None;
      } else {
        state.block = InpState::Block::Ignored;
        ++raw.ignored[keyword + " cards"];
      }
      continue;
    }

    switch (state.block) {
      case InpState::Block::Node: {
        const std::vector<std::string> t = split_commas(line);
        std::vector<std::string> f;
        for (const std::string& s : t) {
          if (!s.empty()) f.push_back(s);
        }
        if (f.size() < 3) throw IoError(where + ": a node needs a tag and at least x, y");
        const long long tag = to_integer(f[0], where);
        raw.node_tags.push_back(tag);
        raw.node_coords.emplace_back(to_real(f[1], where), to_real(f[2], where),
                                     f.size() > 3 ? to_real(f[3], where) : 0.0);
        if (!state.set_name.empty()) raw.node_set_tags[state.set_name].push_back(tag);
        break;
      }
      case InpState::Block::Element: {
        // A definition continues on the next line when it ends with a comma
        // and still lacks nodes; some writers end every line with a comma.
        if (state.pending.empty()) state.pending_line = number;
        state.pending += line;
        std::size_t fields = 0;
        for (const std::string& f : split_commas(state.pending)) fields += f.empty() ? 0 : 1;
        const int expected = state.element_kind.nodes;
        if (line.back() != ',' ||
            (expected > 0 && fields >= 1 + static_cast<std::size_t>(expected))) {
          finish_element(state, raw, source);
        }
        break;
      }
      case InpState::Block::Nset:
      case InpState::Block::Elset: {
        auto& sets = state.block == InpState::Block::Nset ? raw.node_set_tags
                                                          : raw.element_set_tags;
        std::vector<long long>& target = sets[state.set_name];
        const std::vector<std::string> t = split_commas(line);
        std::vector<std::string> f;
        for (const std::string& s : t) {
          if (!s.empty()) f.push_back(s);
        }
        if (state.generate) {
          if (f.size() < 2) throw IoError(where + ": GENERATE needs start, end[, step]");
          const long long a = to_integer(f[0], where);
          const long long b = to_integer(f[1], where);
          const long long step = f.size() > 2 ? to_integer(f[2], where) : 1;
          if (step <= 0 || b < a) throw IoError(where + ": invalid GENERATE range");
          for (long long v = a; v <= b; v += step) target.push_back(v);
        } else {
          for (const std::string& token : f) {
            if (is_number(token)) {
              target.push_back(to_integer(token, where));
            } else {
              // A set may list other sets by name.
              const auto it = sets.find(token);
              if (it == sets.end()) {
                throw IoError(where + ": set '" + state.set_name + "' refers to set '" +
                              token + "', which is not defined before it");
              }
              target.insert(target.end(), it->second.begin(), it->second.end());
            }
          }
        }
        break;
      }
      case InpState::Block::Ignored:
        if (state.keyword == "*INSTANCE") {
          throw IoError(where + ": *INSTANCE carries a translation or rotation, which "
                        "SparLab does not apply; export the model as a flat (non-assembly) "
                        "input file");
        }
        break;
      case InpState::Block::None:
        throw IoError(where + ": data line outside any keyword block");
    }
  }
  finish_element(state, raw, source);
}

}  // namespace

Mesh read_abaqus_inp(std::istream& in, const std::string& source,
                     const std::string& base_directory, const MeshReadOptions& options,
                     MeshReadReport* report_out) {
  MeshReadReport local;
  MeshReadReport& report = report_out ? *report_out : local;
  report = MeshReadReport();
  report.path = source;
  report.format = "abaqus";
  RawMesh raw;
  raw.source = source;
  parse_inp(in, source, base_directory, raw, 0);
  if (raw.node_tags.empty()) throw IoError(source + " defines no *NODE block");
  const auto parts = raw.ignored.find("*PART cards");
  if (parts != raw.ignored.end() && parts->second > 1) {
    std::ostringstream os;
    os << source << " defines " << parts->second
       << " *PART blocks; SparLab reads one flat mesh with a single node numbering. "
          "Export the model as a flat input file (Abaqus/CAE: uncheck \"Do not use parts "
          "and assemblies in input files\" is the default to change) or merge the parts "
          "in the mesher";
    throw IoError(os.str());
  }
  if (raw.ignored.count("*SURFACE cards")) {
    raw.warnings.push_back("*SURFACE definitions are not imported; name a load or support "
                           "region with an *NSET (or *ELSET of boundary elements) instead");
  }
  for (const char* card : {"*BOUNDARY cards", "*CLOAD cards", "*DLOAD cards",
                           "*MATERIAL cards", "*STEP cards"}) {
    if (raw.ignored.count(card)) {
      raw.warnings.push_back(std::string(card) +
                             " in the file are not imported; the SparLab deck defines "
                             "supports, loads and materials");
      break;
    }
  }
  return build_mesh(raw, options, report);
}

Mesh read_mesh_file(const std::string& path, const MeshReadOptions& options,
                    MeshReadReport* report) {
  std::string format = options.format;
  if (format.empty() || format == "auto") {
    const std::size_t dot = path.find_last_of('.');
    const std::string ext = dot == std::string::npos ? "" : upper(path.substr(dot + 1));
    if (ext == "MSH") {
      format = "gmsh";
    } else if (ext == "INP") {
      format = "abaqus";
    } else {
      throw ConfigError("cannot infer the mesh format of '" + path +
                        "' from its extension; set mesh.format to \"gmsh\" or \"abaqus\"");
    }
  }
  std::ifstream in(path);
  if (!in) throw IoError("cannot open mesh file '" + path + "'");
  if (format == "gmsh") return read_gmsh(in, path, options, report);
  if (format == "abaqus" || format == "calculix" || format == "inp") {
    const std::size_t slash = path.find_last_of('/');
    return read_abaqus_inp(in, path, slash == std::string::npos ? "." : path.substr(0, slash),
                           options, report);
  }
  throw ConfigError("mesh.format must be \"auto\", \"gmsh\" or \"abaqus\", got \"" + format +
                    "\"");
}

}  // namespace sparlab
