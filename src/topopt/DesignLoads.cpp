#include "sparlab/topopt/DesignLoads.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Loads.hpp"

#include <sstream>

namespace sparlab {
namespace {

/// Element vectors as the columns of an edofs x ne matrix (zero columns for
/// elements without one); a 0 x 0 matrix when no element has one.
Matrix as_columns(const std::vector<Vector>& vectors, int edofs) {
  bool any = false;
  for (const Vector& v : vectors) any = any || v.size() > 0;
  if (!any) return Matrix();
  Matrix out = Matrix::Zero(edofs, static_cast<Index>(vectors.size()));
  for (std::size_t e = 0; e < vectors.size(); ++e) {
    if (vectors[e].size() > 0) out.col(static_cast<Index>(e)) = vectors[e];
  }
  return out;
}

}  // namespace

DesignLoads::DesignLoads(const FemModel& model) : model_(model) {
  if (!model.finalized()) {
    throw ModelError("design loads requested before FemModel::finalize() was called");
  }
  const int edofs = model.mesh().nodes_per_elem() * model.dofs_per_node();
  const std::vector<LoadCaseSpec>& specs = model.load_case_specs();
  for (std::size_t l = 0; l < specs.size(); ++l) {
    const LoadCaseData& data = model.load_case_data(l);
    fixed_.push_back(data.mechanical);
    body_.push_back(specs[l].has_body_loads()
                        ? as_columns(element_body_loads(model, specs[l]), edofs)
                        : Matrix());
    if (data.temperature.size() > 0) {
      const ElementThermalLoads thermal = element_thermal_loads(model, data.temperature);
      thermal_.push_back(as_columns(thermal.force, edofs));
      self_energy_.push_back(thermal_.back().size() > 0 ? thermal.self_energy : Vector());
    } else {
      thermal_.push_back(Matrix());
      self_energy_.emplace_back();
    }
    dependent_ = dependent_ || has_body(l) || has_thermal(l);
  }
}

const Vector& DesignLoads::temperature(std::size_t l) const {
  return model_.load_case_data(l).temperature;
}

Vector DesignLoads::scatter(const Matrix& columns, const Vector& factor) const {
  const Mesh& mesh = model_.mesh();
  const int npe = mesh.nodes_per_elem();
  Vector out = Vector::Zero(model_.dofs().num_dofs());
  Vector fe(columns.rows());
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    if (factor(e) == 0.0) continue;
    fe = columns.col(e);
    model_.dofs().scatter_add(mesh.element_nodes(e), npe, fe, out, factor(e));
  }
  return out;
}

Vector DesignLoads::load(std::size_t l, const Vector& rho, const SimpOptions& simp) const {
  if (l >= num_cases()) throw ModelError("design load requested for a load case out of range");
  if (!has_body(l) && !has_thermal(l)) return model_.load_vectors()[l];
  const Index ne = model_.mesh().num_elements();
  if (rho.size() != ne) {
    std::ostringstream os;
    os << "design loads: density vector has length " << rho.size() << " but the mesh has "
       << ne << " elements";
    throw ModelError(os.str());
  }
  // The same order as FemModel::finalize: mechanical, then body, then thermal,
  // so a full-density design reproduces the model's load vector.
  Vector f = fixed_[l];
  if (has_body(l)) {
    Vector gamma(ne);
    for (Index e = 0; e < ne; ++e) gamma(e) = body_load_factor(rho(e), simp);
    f += scatter(body_[l], gamma);
  }
  if (has_thermal(l)) f += scatter(thermal_[l], simp_stiffness_factors(rho, simp));
  return f;
}

Vector DesignLoads::contract(std::size_t l, const Vector& rho, const SimpOptions& simp,
                             const Vector& a) const {
  const Mesh& mesh = model_.mesh();
  const Index ne = mesh.num_elements();
  Vector out = Vector::Zero(ne);
  if (!has_body(l) && !has_thermal(l)) return out;
  const int npe = mesh.nodes_per_elem();
  Vector ae(npe * model_.dofs_per_node());
  for (Index e = 0; e < ne; ++e) {
    model_.dofs().gather(mesh.element_nodes(e), npe, a, ae);
    if (has_body(l)) out(e) += body_load_derivative(rho(e), simp) * ae.dot(body_[l].col(e));
    if (has_thermal(l)) {
      out(e) += simp_stiffness_derivative(rho(e), simp) * ae.dot(thermal_[l].col(e));
    }
  }
  return out;
}

Vector DesignLoads::thermal_work(std::size_t l, const Vector& u) const {
  const Mesh& mesh = model_.mesh();
  const Index ne = mesh.num_elements();
  Vector out = Vector::Zero(ne);
  if (!has_thermal(l)) return out;
  const int npe = mesh.nodes_per_elem();
  Vector ue(npe * model_.dofs_per_node());
  for (Index e = 0; e < ne; ++e) {
    model_.dofs().gather(mesh.element_nodes(e), npe, u, ue);
    out(e) = ue.dot(thermal_[l].col(e));
  }
  return out;
}

}  // namespace sparlab
