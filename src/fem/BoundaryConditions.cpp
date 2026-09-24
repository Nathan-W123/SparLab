#include "sparlab/fem/BoundaryConditions.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
#include <set>
#include <sstream>

namespace sparlab {
namespace {

void require_in_plane(const Vector3& v, int dim, const std::string& what) {
  if (dim == 2 && v.z() != 0.0) {
    std::ostringstream os;
    os << what << " has a z component of " << v.z()
       << " but the model is two-dimensional; a plane model carries no out-of-plane "
          "component (use a 3-D mesh or drop it)";
    throw ConfigError(os.str());
  }
}

}  // namespace

Index apply_constraints(const Mesh& mesh,
                        const std::vector<DisplacementConstraint>& constraints,
                        DofManager& dofs) {
  const int dim = mesh.dim();
  std::set<Index> touched;
  for (const DisplacementConstraint& bc : constraints) {
    if (!bc.fix_x && !bc.fix_y && !bc.fix_z) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' constrains no component; remove it or list the "
                        "components to fix");
    }
    if (dim == 2 && bc.fix_z) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' fixes z, but the model is two-dimensional and has no z "
                        "degree of freedom");
    }
    const std::vector<Index> nodes = bc.region.select_nodes(mesh);
    if (nodes.empty()) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' selected no nodes; check its coordinates against the mesh "
                        "extents (an unconstrained model yields a singular stiffness "
                        "matrix)");
    }
    for (Index n : nodes) {
      for (int k = 0; k < dim; ++k) {
        if (!bc.fixes(k)) continue;
        dofs.prescribe(n, k, bc.value(k));
        touched.insert(dofs.dof(n, k));
      }
    }
    log::debug("boundary condition '", bc.region.name, "': ", nodes.size(),
               " nodes, fix_x=", bc.fix_x, " fix_y=", bc.fix_y,
               (dim == 3 ? " fix_z=" : ""), (dim == 3 ? (bc.fix_z ? "1" : "0") : ""));
  }
  return static_cast<Index>(touched.size());
}

Vector assemble_load_vector(const Mesh& mesh, const Element& element,
                            const LoadCaseSpec& load_case, Scalar thickness,
                            const IntegrationOptions& integration) {
  const int dim = mesh.dim();
  Vector f = Vector::Zero(mesh.num_nodes() * dim);

  for (const PointLoadSpec& load : load_case.point_loads) {
    require_in_plane(load.force, dim,
                     "point load '" + load.region.name + "' in load case '" +
                         load_case.name + "'");
    const std::vector<Index> nodes = load.region.select_nodes(mesh);
    if (nodes.empty()) {
      throw ConfigError("point load region '" + load.region.name +
                        "' in load case '" + load_case.name +
                        "' selected no nodes; check its coordinates");
    }
    const Scalar scale =
        load.distribute_total ? 1.0 / static_cast<Scalar>(nodes.size()) : 1.0;
    for (Index n : nodes) {
      for (int k = 0; k < dim; ++k) f(n * dim + k) += scale * load.force(k);
    }
    log::debug("point load '", load.region.name, "' in case '", load_case.name, "': ",
               nodes.size(), " nodes, resultant (", load.force.x(), ", ",
               load.force.y(), (dim == 3 ? ", " : ""), (dim == 3 ? load.force.z() : 0.0),
               ") N", load.distribute_total ? " distributed" : " per node");
  }

  if (!load_case.tractions.empty()) {
    const std::vector<Mesh::BoundaryFace> faces = mesh.boundary_faces();
    for (const TractionLoadSpec& load : load_case.tractions) {
      require_in_plane(load.traction, dim,
                       "traction '" + load.region.name + "' in load case '" +
                           load_case.name + "'");
      const std::vector<Index> nodes = load.region.select_nodes(mesh);
      std::vector<char> in_region(static_cast<std::size_t>(mesh.num_nodes()), 0);
      for (Index n : nodes) in_region[static_cast<std::size_t>(n)] = 1;

      Index matched = 0;
      Scalar total_measure = 0.0;
      for (const Mesh::BoundaryFace& face : faces) {
        const bool inside =
            std::all_of(face.nodes.begin(), face.nodes.end(), [&](Index n) {
              return in_region[static_cast<std::size_t>(n)] != 0;
            });
        if (!inside) continue;
        const Matrix coords = mesh.element_coordinates(face.element);
        const Vector fe = element.boundary_traction(coords, face.local_face,
                                                    load.traction, thickness, integration);
        const Index* enodes = mesh.element_nodes(face.element);
        for (int a = 0; a < mesh.nodes_per_elem(); ++a) {
          for (int k = 0; k < dim; ++k) f(enodes[a] * dim + k) += fe(dim * a + k);
        }
        total_measure += mesh.face_measure(face);
        ++matched;
      }
      if (matched == 0) {
        throw ConfigError("traction region '" + load.region.name + "' in load case '" +
                          load_case.name + "' matched no boundary " +
                          (dim == 2 ? "edge" : "face") +
                          "; the region must contain every node of at least one " +
                          (dim == 2 ? "edge" : "face") + " on the mesh boundary");
      }
      log::debug("traction '", load.region.name, "' in case '", load_case.name, "': ",
                 matched, (dim == 2 ? " boundary edges, length " : " boundary faces, area "),
                 total_measure, (dim == 2 ? " m, traction (" : " m^2, traction ("),
                 load.traction.x(), ", ", load.traction.y(),
                 (dim == 3 ? ", " : ""), (dim == 3 ? load.traction.z() : 0.0), ") Pa");
    }
  }

  return f;
}

}  // namespace sparlab
