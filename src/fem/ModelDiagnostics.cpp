#include "sparlab/fem/ModelDiagnostics.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <numeric>
#include <sstream>
#include <vector>

namespace sparlab {
namespace {

/// Union-find over element indices.
class DisjointSet {
 public:
  explicit DisjointSet(std::size_t n) : parent_(n) {
    std::iota(parent_.begin(), parent_.end(), static_cast<Index>(0));
  }
  Index find(Index a) {
    while (parent_[static_cast<std::size_t>(a)] != a) {
      parent_[static_cast<std::size_t>(a)] =
          parent_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(a)])];
      a = parent_[static_cast<std::size_t>(a)];
    }
    return a;
  }
  void unite(Index a, Index b) {
    a = find(a);
    b = find(b);
    if (a != b) parent_[static_cast<std::size_t>(b)] = a;
  }

 private:
  std::vector<Index> parent_;
};

}  // namespace

std::vector<std::vector<Index>> element_connected_components(const Mesh& mesh) {
  const Index ne = mesh.num_elements();
  const int npe = mesh.nodes_per_elem();

  // First element seen at each node, used to link elements sharing that node.
  std::vector<Index> first_owner(static_cast<std::size_t>(mesh.num_nodes()), -1);
  DisjointSet ds(static_cast<std::size_t>(ne));
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 0; a < npe; ++a) {
      const Index n = nodes[a];
      Index& owner = first_owner[static_cast<std::size_t>(n)];
      if (owner < 0) {
        owner = e;
      } else {
        ds.unite(owner, e);
      }
    }
  }

  std::vector<Index> root_to_component(static_cast<std::size_t>(ne), -1);
  std::vector<std::vector<Index>> components;
  for (Index e = 0; e < ne; ++e) {
    const Index root = ds.find(e);
    Index& cid = root_to_component[static_cast<std::size_t>(root)];
    if (cid < 0) {
      cid = static_cast<Index>(components.size());
      components.emplace_back();
    }
    components[static_cast<std::size_t>(cid)].push_back(e);
  }
  return components;
}

ModelDiagnostics diagnose_model(const FemModel& model) {
  const Mesh& mesh = model.mesh();
  const DofManager& dofs = model.dofs();

  ModelDiagnostics diag;
  diag.num_dofs = dofs.num_dofs();
  diag.num_free_dofs = dofs.num_free();
  diag.num_prescribed_dofs = dofs.num_constrained();

  if (diag.num_prescribed_dofs == 0) {
    diag.problems.push_back(
        "no displacement boundary condition is prescribed: the model can translate and "
        "rotate freely, so the stiffness matrix has a three-dimensional null space. "
        "Add at least three independent constraints.");
  }

  const auto components = element_connected_components(mesh);
  const int npe = mesh.nodes_per_elem();

  for (std::size_t c = 0; c < components.size(); ++c) {
    MeshComponent comp;
    comp.elements = components[c];

    // Collect the component's unique nodes.
    std::vector<Index> nodes;
    nodes.reserve(comp.elements.size() * static_cast<std::size_t>(npe));
    for (Index e : comp.elements) {
      const Index* en = mesh.element_nodes(e);
      for (int a = 0; a < npe; ++a) nodes.push_back(en[a]);
    }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    comp.nodes = nodes;

    // Reference point for the rotation mode: the component's node centroid.
    Vector2 centroid = Vector2::Zero();
    for (Index n : nodes) centroid += mesh.node(n);
    centroid /= static_cast<Scalar>(nodes.size());

    // Rows: prescribed DOFs of this component. Columns: the three rigid modes.
    std::vector<Eigen::Vector3d> rows;
    for (Index n : nodes) {
      const Vector2 x = mesh.node(n) - centroid;
      for (int k = 0; k < kDofsPerNode; ++k) {
        const Index d = dofs.dof(n, k);
        if (!dofs.is_constrained(d)) continue;
        Eigen::Vector3d row;
        if (k == 0) {
          row << 1.0, 0.0, -x.y();  // u_x of (tx, ty, theta)
        } else {
          row << 0.0, 1.0, x.x();   // u_y of (tx, ty, theta)
        }
        rows.push_back(row);
      }
    }
    comp.prescribed_dofs = static_cast<Index>(rows.size());

    if (rows.empty()) {
      comp.rigid_null_dimension = 3;
    } else {
      Eigen::MatrixXd r(static_cast<Eigen::Index>(rows.size()), 3);
      for (std::size_t i = 0; i < rows.size(); ++i) r.row(static_cast<Eigen::Index>(i)) = rows[i];
      // Scale the rotation column so the rank test is dimensionally sensible:
      // translations are O(1) while the rotation column is O(length).
      const Scalar rot_scale = r.col(2).cwiseAbs().maxCoeff();
      if (rot_scale > 0.0) r.col(2) /= rot_scale;
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(r);
      const Eigen::VectorXd sv = svd.singularValues();
      const Scalar tol = 1.0e-10 * std::max(sv(0), 1.0);
      int rank = 0;
      for (Eigen::Index i = 0; i < sv.size(); ++i) {
        if (sv(i) > tol) ++rank;
      }
      comp.rigid_null_dimension = 3 - rank;
    }

    if (comp.rigid_null_dimension > 0) {
      std::ostringstream os;
      os << "element group " << c << " (" << comp.elements.size() << " elements, "
         << comp.nodes.size() << " nodes, centroid at (" << centroid.x() << ", "
         << centroid.y() << ") m) retains " << comp.rigid_null_dimension
         << " rigid-body degree(s) of freedom";
      if (comp.prescribed_dofs == 0) {
        os << " because it has no prescribed DOF at all (a floating region)";
      } else {
        os << " despite " << comp.prescribed_dofs
           << " prescribed DOF(s); the prescribed DOFs do not suppress translation in "
              "both directions plus rotation";
      }
      os << ". Add displacement constraints to that region.";
      diag.problems.push_back(os.str());
    }

    diag.components.push_back(std::move(comp));
  }

  if (components.size() > 1) {
    std::ostringstream os;
    os << "mesh consists of " << components.size()
       << " disconnected element groups; each one needs its own constraints";
    log::warn(os.str());
  }

  return diag;
}

void require_well_posed(const FemModel& model) {
  const ModelDiagnostics diag = diagnose_model(model);
  if (diag.well_posed()) {
    log::debug("model diagnostics passed: ", diag.num_free_dofs, " free DOFs, ",
               diag.num_prescribed_dofs, " prescribed, ", diag.components.size(),
               " connected element group(s)");
    return;
  }
  std::ostringstream os;
  os << "the model is not well posed for a static solve (" << diag.problems.size()
     << " problem(s)):";
  for (const std::string& p : diag.problems) os << "\n  - " << p;
  throw ModelError(os.str());
}

}  // namespace sparlab
