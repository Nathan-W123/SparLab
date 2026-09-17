#include "sparlab/fem/BoundaryConditions.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
#include <set>
#include <sstream>

namespace sparlab {

Index apply_constraints(const Mesh& mesh,
                        const std::vector<DisplacementConstraint>& constraints,
                        DofManager& dofs) {
  std::set<Index> touched;
  for (const DisplacementConstraint& bc : constraints) {
    if (!bc.fix_x && !bc.fix_y) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' constrains neither x nor y; remove it or list the "
                        "components to fix");
    }
    const std::vector<Index> nodes = bc.region.select_nodes(mesh);
    if (nodes.empty()) {
      throw ConfigError("boundary condition '" + bc.region.name +
                        "' selected no nodes; check its coordinates against the mesh "
                        "extents (an unconstrained model yields a singular stiffness "
                        "matrix)");
    }
    for (Index n : nodes) {
      if (bc.fix_x) {
        dofs.prescribe(n, 0, bc.value_x);
        touched.insert(dofs.dof(n, 0));
      }
      if (bc.fix_y) {
        dofs.prescribe(n, 1, bc.value_y);
        touched.insert(dofs.dof(n, 1));
      }
    }
    log::debug("boundary condition '", bc.region.name, "': ", nodes.size(),
               " nodes, fix_x=", bc.fix_x, " fix_y=", bc.fix_y);
  }
  return static_cast<Index>(touched.size());
}

Vector assemble_load_vector(const Mesh& mesh, const Element& element,
                            const LoadCaseSpec& load_case, Scalar thickness,
                            const IntegrationOptions& integration) {
  Vector f = Vector::Zero(mesh.num_nodes() * kDofsPerNode);

  for (const PointLoadSpec& load : load_case.point_loads) {
    const std::vector<Index> nodes = load.region.select_nodes(mesh);
    if (nodes.empty()) {
      throw ConfigError("point load region '" + load.region.name +
                        "' in load case '" + load_case.name +
                        "' selected no nodes; check its coordinates");
    }
    const Scalar scale =
        load.distribute_total ? 1.0 / static_cast<Scalar>(nodes.size()) : 1.0;
    for (Index n : nodes) {
      f(n * kDofsPerNode + 0) += scale * load.force.x();
      f(n * kDofsPerNode + 1) += scale * load.force.y();
    }
    log::debug("point load '", load.region.name, "' in case '", load_case.name, "': ",
               nodes.size(), " nodes, resultant (", load.force.x(), ", ",
               load.force.y(), ") N",
               load.distribute_total ? " distributed" : " per node");
  }

  if (!load_case.tractions.empty()) {
    const std::vector<Mesh::BoundaryEdge> edges = mesh.boundary_edges();
    for (const TractionLoadSpec& load : load_case.tractions) {
      const std::vector<Index> nodes = load.region.select_nodes(mesh);
      std::vector<char> in_region(static_cast<std::size_t>(mesh.num_nodes()), 0);
      for (Index n : nodes) in_region[static_cast<std::size_t>(n)] = 1;

      Index matched = 0;
      Scalar total_length = 0.0;
      for (const Mesh::BoundaryEdge& edge : edges) {
        if (!in_region[static_cast<std::size_t>(edge.node_a)] ||
            !in_region[static_cast<std::size_t>(edge.node_b)]) {
          continue;
        }
        const auto coords = mesh.element_coordinates(edge.element);
        const Vector fe = element.edge_traction(coords, edge.local_edge, load.traction,
                                                thickness, integration);
        const Index* enodes = mesh.element_nodes(edge.element);
        for (int a = 0; a < mesh.nodes_per_elem(); ++a) {
          f(enodes[a] * kDofsPerNode + 0) += fe(kDofsPerNode * a + 0);
          f(enodes[a] * kDofsPerNode + 1) += fe(kDofsPerNode * a + 1);
        }
        total_length += (mesh.node(edge.node_b) - mesh.node(edge.node_a)).norm();
        ++matched;
      }
      if (matched == 0) {
        throw ConfigError("traction region '" + load.region.name + "' in load case '" +
                          load_case.name +
                          "' matched no boundary edge; the region must contain both end "
                          "nodes of at least one edge on the mesh boundary");
      }
      log::debug("traction '", load.region.name, "' in case '", load_case.name, "': ",
                 matched, " boundary edges, length ", total_length, " m, traction (",
                 load.traction.x(), ", ", load.traction.y(), ") Pa");
    }
  }

  return f;
}

}  // namespace sparlab
