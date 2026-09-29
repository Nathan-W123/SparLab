#include "sparlab/mesh/StructuredMesh.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <Eigen/Geometry>

#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <utility>

namespace sparlab {
namespace {

void check_counts_2d(const StructuredMeshSpec& spec) {
  if (spec.nx < 1 || spec.ny < 1) {
    std::ostringstream os;
    os << "structured mesh needs at least one element per direction (got nx=" << spec.nx
       << ", ny=" << spec.ny << ")";
    throw ConfigError(os.str());
  }
  if (!(spec.lx > 0.0) || !(spec.ly > 0.0)) {
    std::ostringstream os;
    os << "structured mesh needs positive extents (got lx=" << spec.lx
       << " m, ly=" << spec.ly << " m)";
    throw ConfigError(os.str());
  }
}

void check_counts_3d(const StructuredMeshSpec& spec) {
  check_counts_2d(spec);
  if (spec.nz < 1) {
    std::ostringstream os;
    os << "structured hex mesh needs at least one element along z (got nz=" << spec.nz
       << ")";
    throw ConfigError(os.str());
  }
  if (!(spec.lz > 0.0)) {
    std::ostringstream os;
    os << "structured hex mesh needs a positive extent along z (got lz=" << spec.lz
       << " m)";
    throw ConfigError(os.str());
  }
}

/// Split every Q4 of a structured grid into two triangles. Alternating the
/// diagonal between neighbouring cells removes the directional bias a single
/// diagonal direction would give the stiffness.
Mesh split_into_triangles(const Mesh& quads, const StructuredMeshSpec& spec) {
  std::vector<Index> connectivity;
  connectivity.reserve(static_cast<std::size_t>(quads.num_elements()) * 6);
  for (Index j = 0; j < spec.ny; ++j) {
    for (Index i = 0; i < spec.nx; ++i) {
      const Index* n = quads.element_nodes(j * spec.nx + i);
      if ((i + j) % 2 == 0) {
        connectivity.insert(connectivity.end(), {n[0], n[1], n[2], n[0], n[2], n[3]});
      } else {
        connectivity.insert(connectivity.end(), {n[0], n[1], n[3], n[1], n[2], n[3]});
      }
    }
  }
  Mesh mesh(quads.coordinates(), std::move(connectivity), ElementType::Tri3);
  mesh.validate();
  return mesh;
}

/// Split every Hex8 into the six Kuhn tetrahedra around the diagonal from
/// local corner 0 to corner 6. Each tetrahedron follows one monotone path
/// 0 -> 6 along the cell edges; every face of the cube is then cut along the
/// diagonal joining its lowest and highest corners, which is the same
/// diagonal the neighbouring cell uses, so the split is conforming. Paths
/// that are odd permutations of (x, y, z) come out negatively oriented and
/// get two nodes swapped.
Mesh split_into_tetrahedra(const Mesh& hexes) {
  static const int paths[6][4] = {{0, 1, 2, 6}, {0, 1, 5, 6}, {0, 3, 2, 6},
                                  {0, 3, 7, 6}, {0, 4, 5, 6}, {0, 4, 7, 6}};
  const Matrix& x = hexes.coordinates();
  std::vector<Index> connectivity;
  connectivity.reserve(static_cast<std::size_t>(hexes.num_elements()) * 24);
  for (Index e = 0; e < hexes.num_elements(); ++e) {
    const Index* n = hexes.element_nodes(e);
    for (const auto& path : paths) {
      Index t[4] = {n[path[0]], n[path[1]], n[path[2]], n[path[3]]};
      const Vector3 a = x.col(t[1]) - x.col(t[0]);
      const Vector3 b = x.col(t[2]) - x.col(t[0]);
      const Vector3 c = x.col(t[3]) - x.col(t[0]);
      if (a.dot(b.cross(c)) < 0.0) std::swap(t[1], t[2]);
      connectivity.insert(connectivity.end(), {t[0], t[1], t[2], t[3]});
    }
  }
  Mesh mesh(hexes.coordinates(), std::move(connectivity), ElementType::Tet4);
  mesh.validate();
  return mesh;
}

StructuredGridInfo grid_info(const StructuredMeshSpec& spec, bool three_d, bool uniform) {
  StructuredGridInfo info;
  info.nx = spec.nx;
  info.ny = spec.ny;
  info.nz = three_d ? spec.nz : 0;
  info.lx = spec.lx;
  info.ly = spec.ly;
  info.lz = three_d ? spec.lz : 0.0;
  info.uniform = uniform;
  return info;
}

}  // namespace

Index structured_node_index(const StructuredGridInfo& info, Index i, Index j) {
  if (i < 0 || i > info.nx || j < 0 || j > info.ny) {
    std::ostringstream os;
    os << "structured node (" << i << ", " << j << ") is outside the grid [0.."
       << info.nx << "] x [0.." << info.ny << "]";
    throw MeshError(os.str());
  }
  return j * (info.nx + 1) + i;
}

Index structured_element_index(const StructuredGridInfo& info, Index i, Index j) {
  if (i < 0 || i >= info.nx || j < 0 || j >= info.ny) {
    std::ostringstream os;
    os << "structured element (" << i << ", " << j << ") is outside the grid [0.."
       << info.nx - 1 << "] x [0.." << info.ny - 1 << "]";
    throw MeshError(os.str());
  }
  return j * info.nx + i;
}

Index structured_node_index(const StructuredGridInfo& info, Index i, Index j, Index k) {
  if (i < 0 || i > info.nx || j < 0 || j > info.ny || k < 0 || k > info.nz) {
    std::ostringstream os;
    os << "structured node (" << i << ", " << j << ", " << k
       << ") is outside the grid [0.." << info.nx << "] x [0.." << info.ny << "] x [0.."
       << info.nz << "]";
    throw MeshError(os.str());
  }
  return k * (info.nx + 1) * (info.ny + 1) + j * (info.nx + 1) + i;
}

Index structured_element_index(const StructuredGridInfo& info, Index i, Index j,
                               Index k) {
  if (i < 0 || i >= info.nx || j < 0 || j >= info.ny || k < 0 || k >= info.nz) {
    std::ostringstream os;
    os << "structured element (" << i << ", " << j << ", " << k
       << ") is outside the grid [0.." << info.nx - 1 << "] x [0.." << info.ny - 1
       << "] x [0.." << info.nz - 1 << "]";
    throw MeshError(os.str());
  }
  return k * info.nx * info.ny + j * info.nx + i;
}

Mesh make_structured_quad_mesh(const StructuredMeshSpec& spec) {
  check_counts_2d(spec);

  const Index nnx = spec.nx + 1;
  const Index nny = spec.ny + 1;
  const Scalar dx = spec.lx / static_cast<Scalar>(spec.nx);
  const Scalar dy = spec.ly / static_cast<Scalar>(spec.ny);

  Matrix coords(2, nnx * nny);
  for (Index j = 0; j < nny; ++j) {
    for (Index i = 0; i < nnx; ++i) {
      const Index n = j * nnx + i;
      coords(0, n) = spec.x0 + static_cast<Scalar>(i) * dx;
      coords(1, n) = spec.y0 + static_cast<Scalar>(j) * dy;
    }
  }

  std::vector<Index> connectivity;
  connectivity.reserve(static_cast<std::size_t>(spec.nx) * spec.ny * 4);
  for (Index j = 0; j < spec.ny; ++j) {
    for (Index i = 0; i < spec.nx; ++i) {
      const Index n0 = j * nnx + i;
      connectivity.push_back(n0);
      connectivity.push_back(n0 + 1);
      connectivity.push_back(n0 + nnx + 1);
      connectivity.push_back(n0 + nnx);
    }
  }

  Mesh mesh(std::move(coords), std::move(connectivity), ElementType::Quad4);
  mesh.set_structured_info(grid_info(spec, false, true));
  mesh.validate();

  log::debug("generated structured Q4 mesh: ", spec.nx, " x ", spec.ny, " = ",
             mesh.num_elements(), " elements, ", mesh.num_nodes(), " nodes, cell ", dx,
             " x ", dy, " m");
  return mesh;
}

Mesh make_structured_hex_mesh(const StructuredMeshSpec& spec) {
  check_counts_3d(spec);

  const Index nnx = spec.nx + 1;
  const Index nny = spec.ny + 1;
  const Index nnz = spec.nz + 1;
  const Scalar dx = spec.lx / static_cast<Scalar>(spec.nx);
  const Scalar dy = spec.ly / static_cast<Scalar>(spec.ny);
  const Scalar dz = spec.lz / static_cast<Scalar>(spec.nz);

  Matrix coords(3, nnx * nny * nnz);
  for (Index k = 0; k < nnz; ++k) {
    for (Index j = 0; j < nny; ++j) {
      for (Index i = 0; i < nnx; ++i) {
        const Index n = k * nnx * nny + j * nnx + i;
        coords(0, n) = spec.x0 + static_cast<Scalar>(i) * dx;
        coords(1, n) = spec.y0 + static_cast<Scalar>(j) * dy;
        coords(2, n) = spec.z0 + static_cast<Scalar>(k) * dz;
      }
    }
  }

  std::vector<Index> connectivity;
  connectivity.reserve(static_cast<std::size_t>(spec.nx) * spec.ny * spec.nz * 8);
  const Index layer = nnx * nny;
  for (Index k = 0; k < spec.nz; ++k) {
    for (Index j = 0; j < spec.ny; ++j) {
      for (Index i = 0; i < spec.nx; ++i) {
        const Index n0 = k * layer + j * nnx + i;
        // bottom face, counter-clockwise seen from +z
        connectivity.push_back(n0);
        connectivity.push_back(n0 + 1);
        connectivity.push_back(n0 + nnx + 1);
        connectivity.push_back(n0 + nnx);
        // top face, same order
        connectivity.push_back(n0 + layer);
        connectivity.push_back(n0 + layer + 1);
        connectivity.push_back(n0 + layer + nnx + 1);
        connectivity.push_back(n0 + layer + nnx);
      }
    }
  }

  Mesh mesh(std::move(coords), std::move(connectivity), ElementType::Hex8);
  mesh.set_structured_info(grid_info(spec, true, true));
  mesh.validate();

  log::debug("generated structured Hex8 mesh: ", spec.nx, " x ", spec.ny, " x ", spec.nz,
             " = ", mesh.num_elements(), " elements, ", mesh.num_nodes(),
             " nodes, cell ", dx, " x ", dy, " x ", dz, " m");
  return mesh;
}

Mesh make_perturbed_quad_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                              unsigned int seed) {
  if (!(perturbation >= 0.0 && perturbation < 0.45)) {
    std::ostringstream os;
    os << "mesh perturbation must lie in [0, 0.45) to keep elements convex (got "
       << perturbation << ")";
    throw ConfigError(os.str());
  }
  Mesh mesh = make_structured_quad_mesh(spec);
  if (perturbation == 0.0) return mesh;

  const Index nnx = spec.nx + 1;
  const Index nny = spec.ny + 1;
  const Scalar dx = spec.lx / static_cast<Scalar>(spec.nx);
  const Scalar dy = spec.ly / static_cast<Scalar>(spec.ny);

  Matrix coords = mesh.coordinates();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
  for (Index j = 1; j < nny - 1; ++j) {
    for (Index i = 1; i < nnx - 1; ++i) {
      const Index n = j * nnx + i;
      coords(0, n) += perturbation * dx * dist(rng);
      coords(1, n) += perturbation * dy * dist(rng);
    }
  }

  Mesh perturbed(std::move(coords), mesh.connectivity(), ElementType::Quad4);
  // The grid is no longer uniform, so the element-matrix cache must not be
  // used; recording nx/ny keeps the index helpers usable.
  perturbed.set_structured_info(grid_info(spec, false, false));
  perturbed.validate();
  log::debug("perturbed structured mesh: ", spec.nx, " x ", spec.ny,
             " with interior nodes displaced by up to ", perturbation,
             " of the cell size");
  return perturbed;
}

Mesh make_perturbed_hex_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                             unsigned int seed) {
  if (!(perturbation >= 0.0 && perturbation < 0.35)) {
    std::ostringstream os;
    os << "hex mesh perturbation must lie in [0, 0.35) to keep cells valid (got "
       << perturbation << ")";
    throw ConfigError(os.str());
  }
  Mesh mesh = make_structured_hex_mesh(spec);
  if (perturbation == 0.0) return mesh;

  const Index nnx = spec.nx + 1;
  const Index nny = spec.ny + 1;
  const Index nnz = spec.nz + 1;
  const Scalar dx = spec.lx / static_cast<Scalar>(spec.nx);
  const Scalar dy = spec.ly / static_cast<Scalar>(spec.ny);
  const Scalar dz = spec.lz / static_cast<Scalar>(spec.nz);

  Matrix coords = mesh.coordinates();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
  for (Index k = 1; k < nnz - 1; ++k) {
    for (Index j = 1; j < nny - 1; ++j) {
      for (Index i = 1; i < nnx - 1; ++i) {
        const Index n = k * nnx * nny + j * nnx + i;
        coords(0, n) += perturbation * dx * dist(rng);
        coords(1, n) += perturbation * dy * dist(rng);
        coords(2, n) += perturbation * dz * dist(rng);
      }
    }
  }

  Mesh perturbed(std::move(coords), mesh.connectivity(), ElementType::Hex8);
  perturbed.set_structured_info(grid_info(spec, true, false));
  perturbed.validate();
  log::debug("perturbed structured hex mesh: ", spec.nx, " x ", spec.ny, " x ", spec.nz,
             " with interior nodes displaced by up to ", perturbation,
             " of the cell size");
  return perturbed;
}

Mesh make_structured_tri_mesh(const StructuredMeshSpec& spec) {
  return split_into_triangles(make_structured_quad_mesh(spec), spec);
}

Mesh make_structured_tet_mesh(const StructuredMeshSpec& spec) {
  return split_into_tetrahedra(make_structured_hex_mesh(spec));
}

Mesh make_structured_tet10_mesh(const StructuredMeshSpec& spec) {
  return elevate_to_tet10(make_structured_tet_mesh(spec));
}

Mesh make_perturbed_tri_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                             unsigned int seed) {
  return split_into_triangles(make_perturbed_quad_mesh(spec, perturbation, seed), spec);
}

Mesh make_perturbed_tet_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                             unsigned int seed) {
  return split_into_tetrahedra(make_perturbed_hex_mesh(spec, perturbation, seed));
}

std::string to_string(ShellShape shape) {
  switch (shape) {
    case ShellShape::Plate: return "plate";
    case ShellShape::Cylinder: return "cylinder";
    case ShellShape::Sphere: return "sphere";
  }
  return "unknown";
}

ShellShape parse_shell_shape(const std::string& text) {
  if (text == "plate") return ShellShape::Plate;
  if (text == "cylinder") return ShellShape::Cylinder;
  if (text == "sphere") return ShellShape::Sphere;
  throw ConfigError("shell shape must be \"plate\", \"cylinder\" or \"sphere\", got '" +
                    text + "'");
}

Mesh make_structured_shell_mesh(const ShellMeshSpec& spec) {
  constexpr Scalar kPi = 3.14159265358979323846;
  constexpr Scalar kDegree = kPi / 180.0;
  const auto fail = [](const std::string& text) { throw ConfigError("shell mesh: " + text); };
  if (spec.n1 < 1 || spec.n2 < 1) {
    std::ostringstream os;
    os << "at least one element per direction is needed (got n1 = " << spec.n1
       << ", n2 = " << spec.n2 << ")";
    fail(os.str());
  }
  const Scalar span = spec.angle_end - spec.angle_start;
  bool closed = false;
  if (spec.shape == ShellShape::Plate) {
    if (!(spec.lx > 0.0) || !(spec.ly > 0.0)) {
      std::ostringstream os;
      os << "a plate needs positive extents (got lx = " << spec.lx << " m, ly = " << spec.ly
         << " m)";
      fail(os.str());
    }
    // Moving the corners of a rectangle by up to p of its sides along each
    // axis keeps it convex while p < 1/4 (at 1/4 three corners can line up).
    if (!(spec.perturbation >= 0.0 && spec.perturbation < 0.25)) {
      std::ostringstream os;
      os << "the perturbation must lie in [0, 0.25) to keep every cell convex (got "
         << spec.perturbation << ")";
      fail(os.str());
    }
  } else {
    if (!(spec.radius > 0.0)) {
      std::ostringstream os;
      os << "a " << to_string(spec.shape) << " needs a positive radius (got " << spec.radius
         << " m)";
      fail(os.str());
    }
    if (!(span > 0.0) || span > 360.0 + 1.0e-12) {
      std::ostringstream os;
      os << "the angles must run upwards over at most 360 degrees (got " << spec.angle_start
         << " to " << spec.angle_end << ")";
      fail(os.str());
    }
    closed = std::abs(span - 360.0) <= 1.0e-12;
    if (closed && spec.n1 < 3) {
      std::ostringstream os;
      os << "a closed circle needs at least 3 elements around it (got n1 = " << spec.n1 << ")";
      fail(os.str());
    }
    if (spec.shape == ShellShape::Cylinder) {
      if (!(spec.length > 0.0)) {
        std::ostringstream os;
        os << "a cylinder needs a positive length (got " << spec.length << " m)";
        fail(os.str());
      }
      if (spec.axis < 0 || spec.axis > 2) {
        std::ostringstream os;
        os << "the cylinder axis must be 0 (x), 1 (y) or 2 (z) (got " << spec.axis << ")";
        fail(os.str());
      }
    } else if (!(spec.polar_start > 0.0 && spec.polar_start < spec.polar_end &&
                 spec.polar_end < 180.0)) {
      std::ostringstream os;
      os << "the polar angles must satisfy 0 < start < end < 180 degrees (got "
         << spec.polar_start << " to " << spec.polar_end
         << "); at a pole the quadrilaterals would collapse into triangles, so leave an "
            "opening there (as the pinched hemisphere does)";
      fail(os.str());
    }
  }

  const Index m1 = closed ? spec.n1 : spec.n1 + 1;  // node columns along direction 1
  const Index m2 = spec.n2 + 1;
  Matrix coords(3, m1 * m2);
  Matrix normals(3, m1 * m2);
  const Scalar d1 = 1.0 / static_cast<Scalar>(spec.n1);
  const Scalar d2 = 1.0 / static_cast<Scalar>(spec.n2);
  for (Index j = 0; j < m2; ++j) {
    for (Index i = 0; i < m1; ++i) {
      const Scalar a = static_cast<Scalar>(i) * d1;  // in [0, 1] along direction 1
      const Scalar b = static_cast<Scalar>(j) * d2;  // along direction 2
      const Index n = j * m1 + i;
      switch (spec.shape) {
        case ShellShape::Plate:
          coords.col(n) = spec.origin + Vector3(a * spec.lx, b * spec.ly, 0.0);
          normals.col(n) = Vector3::UnitZ();
          break;
        case ShellShape::Cylinder: {
          const Vector3 e = Vector3::Unit(spec.axis);
          const Vector3 p = Vector3::Unit((spec.axis + 1) % 3);
          const Vector3 q = Vector3::Unit((spec.axis + 2) % 3);
          const Scalar theta = (spec.angle_start + a * span) * kDegree;
          const Vector3 radial = std::cos(theta) * p + std::sin(theta) * q;
          coords.col(n) = spec.origin + spec.radius * radial + (b * spec.length) * e;
          normals.col(n) = radial;
          break;
        }
        case ShellShape::Sphere: {
          const Scalar lambda = (spec.angle_start + a * span) * kDegree;
          const Scalar phi =
              (spec.polar_end - b * (spec.polar_end - spec.polar_start)) * kDegree;
          const Vector3 radial(std::sin(phi) * std::cos(lambda), std::sin(phi) * std::sin(lambda),
                               std::cos(phi));
          coords.col(n) = spec.origin + spec.radius * radial;
          normals.col(n) = radial;
          break;
        }
      }
    }
  }
  if (spec.shape == ShellShape::Plate && spec.perturbation > 0.0) {
    std::mt19937 rng(spec.seed);
    std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
    const Scalar dx = spec.lx * d1;
    const Scalar dy = spec.ly * d2;
    for (Index j = 1; j + 1 < m2; ++j) {
      for (Index i = 1; i + 1 < m1; ++i) {
        const Index n = j * m1 + i;
        coords(0, n) += spec.perturbation * dx * dist(rng);
        coords(1, n) += spec.perturbation * dy * dist(rng);
      }
    }
  }
  std::vector<Index> connectivity;
  connectivity.reserve(static_cast<std::size_t>(spec.n1 * spec.n2) * 4);
  for (Index j = 0; j < spec.n2; ++j) {
    for (Index i = 0; i < spec.n1; ++i) {
      const Index i1 = closed ? (i + 1) % spec.n1 : i + 1;
      connectivity.insert(connectivity.end(),
                          {j * m1 + i, j * m1 + i1, (j + 1) * m1 + i1, (j + 1) * m1 + i});
    }
  }
  Mesh mesh(std::move(coords), std::move(connectivity), ElementType::Shell4);
  mesh.set_node_normals(std::move(normals));
  mesh.validate();
  log::debug("structured shell mesh: ", to_string(spec.shape), ", ", spec.n1, " x ", spec.n2,
             " elements", closed ? " (closed around)" : "");
  return mesh;
}

Mesh make_perturbed_tet10_mesh(const StructuredMeshSpec& spec, Scalar perturbation,
                               unsigned int seed) {
  return elevate_to_tet10(make_perturbed_tet_mesh(spec, perturbation, seed));
}

}  // namespace sparlab
