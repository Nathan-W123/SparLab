#include "sparlab/fem/FemModel.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <numeric>
#include <sstream>
#include <utility>

namespace sparlab {

std::string to_string(MassType type) {
  switch (type) {
    case MassType::Consistent: return "consistent";
    case MassType::Lumped: return "lumped";
  }
  return "unknown";
}

FemModel::FemModel(Mesh mesh, IsotropicMaterial material, Scalar thickness,
                   StressState stress_state, IntegrationOptions integration)
    : mesh_(std::move(mesh)),
      material_(std::move(material)),
      thickness_(thickness),
      stress_state_(stress_state),
      integration_(integration),
      element_(make_element(mesh_.element_type())),
      d_(material_.constitutive(stress_state)),
      dofs_(mesh_.num_nodes()) {
  if (!(thickness_ > 0.0)) {
    std::ostringstream os;
    os << "model thickness must be positive (got " << thickness_ << " m)";
    throw ConfigError(os.str());
  }
  mesh_.validate();
}

void FemModel::set_material(const IsotropicMaterial& material) {
  material_ = material;
  d_ = material_.constitutive(stress_state_);
}

void FemModel::finalize(bool require_load_cases) {
  if (finalized_) return;
  if (load_case_specs_.empty() && require_load_cases) {
    throw ConfigError("model has no load cases; define at least one");
  }
  const Index constrained = apply_constraints(mesh_, constraints_, dofs_);
  log::info("applied ", constraints_.size(), " boundary condition group(s): ",
            constrained, " of ", dofs_.num_dofs(), " DOFs prescribed, ",
            dofs_.num_free(), " free");

  load_vectors_.clear();
  load_vectors_.reserve(load_case_specs_.size());
  for (const LoadCaseSpec& spec : load_case_specs_) {
    if (!(spec.weight >= 0.0)) {
      std::ostringstream os;
      os << "load case '" << spec.name << "' has a negative weight (" << spec.weight
         << "); weights must be non-negative";
      throw ConfigError(os.str());
    }
    load_vectors_.push_back(
        assemble_load_vector(mesh_, *element_, spec, thickness_, integration_));
    // A zero load vector is a mistake *unless* the case is driven by prescribed
    // displacements, which is how the patch test and any enforced-deflection
    // study work.
    if (load_vectors_.back().norm() == 0.0 && !dofs_.has_nonzero_prescribed() &&
        !spec.prescribed_displacement_only) {
      log::warn("load case '", spec.name,
                "' has a zero resultant force vector and no non-zero prescribed "
                "displacement, so its solution is identically zero");
    }
  }

  if (!load_case_specs_.empty()) {
    Scalar weight_sum = 0.0;
    for (const LoadCaseSpec& spec : load_case_specs_) weight_sum += spec.weight;
    if (!(weight_sum > 0.0)) {
      throw ConfigError(
          "the sum of load-case weights is zero; the compliance objective would be "
          "identically zero");
    }
  }
  finalized_ = true;
}

const std::vector<Vector>& FemModel::load_vectors() const {
  if (!finalized_) {
    throw ModelError("load vectors requested before FemModel::finalize() was called");
  }
  return load_vectors_;
}

std::vector<Scalar> FemModel::normalised_weights() const {
  Scalar sum = 0.0;
  for (const LoadCaseSpec& spec : load_case_specs_) sum += spec.weight;
  if (!(sum > 0.0)) return {};
  std::vector<Scalar> w;
  w.reserve(load_case_specs_.size());
  for (const LoadCaseSpec& spec : load_case_specs_) w.push_back(spec.weight / sum);
  return w;
}

Vector FemModel::element_volumes() const {
  const Index ne = mesh_.num_elements();
  Vector v(ne);
  for (Index e = 0; e < ne; ++e) v(e) = mesh_.element_area(e) * thickness_;
  return v;
}

Scalar FemModel::domain_volume() const { return element_volumes().sum(); }

}  // namespace sparlab
