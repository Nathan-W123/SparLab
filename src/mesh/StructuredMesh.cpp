#include "sparlab/mesh/StructuredMesh.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <random>
#include <sstream>

namespace sparlab {

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

Mesh make_structured_quad_mesh(const StructuredMeshSpec& spec) {
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

  const Index nnx = spec.nx + 1;
  const Index nny = spec.ny + 1;
  const Scalar dx = spec.lx / static_cast<Scalar>(spec.nx);
  const Scalar dy = spec.ly / static_cast<Scalar>(spec.ny);

  Eigen::Matrix2Xd coords(2, nnx * nny);
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
  StructuredGridInfo info;
  info.nx = spec.nx;
  info.ny = spec.ny;
  info.lx = spec.lx;
  info.ly = spec.ly;
  info.uniform = true;
  mesh.set_structured_info(info);
  mesh.validate();

  log::debug("generated structured Q4 mesh: ", spec.nx, " x ", spec.ny, " = ",
             mesh.num_elements(), " elements, ", mesh.num_nodes(), " nodes, cell ", dx,
             " x ", dy, " m");
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

  Eigen::Matrix2Xd coords = mesh.coordinates();
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
  StructuredGridInfo info;
  info.nx = spec.nx;
  info.ny = spec.ny;
  info.lx = spec.lx;
  info.ly = spec.ly;
  info.uniform = false;
  perturbed.set_structured_info(info);
  perturbed.validate();
  log::debug("perturbed structured mesh: ", spec.nx, " x ", spec.ny,
             " with interior nodes displaced by up to ", perturbation,
             " of the cell size");
  return perturbed;
}

}  // namespace sparlab
