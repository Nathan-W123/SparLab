#include "sparlab/fem/FemModel.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/Shell4.hpp"
#include "sparlab/fem/HeatConduction.hpp"
#include "sparlab/fem/LinearSolver.hpp"
#include "sparlab/fem/Loads.hpp"

#include <array>
#include <cmath>
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
      materials_{std::move(material)},
      thickness_(thickness),
      stress_state_(stress_state),
      integration_(integration),
      element_(make_element(mesh_.element_type())),
      d_{materials_.front().constitutive(stress_state)},
      dofs_(mesh_.num_nodes(), element_->dofs_per_node()) {
  const bool shell_mesh = sparlab::is_shell(mesh_.element_type());
  if (shell_mesh != (stress_state_ == StressState::Shell)) {
    std::ostringstream os;
    os << (shell_mesh ? "a shell mesh (Shell4 elements) needs the stress state 'shell'"
                      : "the stress state 'shell' needs a shell mesh (Shell4 elements)")
       << ", got '" << to_string(stress_state_) << "' with " << to_string(mesh_.element_type())
       << " elements";
    throw ConfigError(os.str());
  }
  if (stress_state_dimension(stress_state_) != mesh_.dim()) {
    std::ostringstream os;
    os << "stress state '" << to_string(stress_state_) << "' belongs to a "
       << stress_state_dimension(stress_state_) << "-D model but the mesh is "
       << mesh_.dim() << "-D (" << to_string(mesh_.element_type()) << " elements); use "
       << (mesh_.dim() == 3 ? "'three_dimensional'" : "'plane_stress' or 'plane_strain'");
    throw ConfigError(os.str());
  }
  if (!(thickness_ > 0.0)) {
    std::ostringstream os;
    os << "model thickness must be positive (got " << thickness_ << " m)";
    throw ConfigError(os.str());
  }
  if (mesh_.dim() == 3 && !shell_mesh && thickness_ != 1.0) {
    std::ostringstream os;
    os << "a 3-D model has no thickness (got " << thickness_
       << " m); leave model.thickness at its default of 1 for a solid mesh";
    throw ConfigError(os.str());
  }
  mesh_.validate();
  if (shell_mesh) {
    element_ = std::make_unique<Shell4Element>(shell_options_.drilling_factor);
    compute_directors();
  }
}

void FemModel::assign_thickness(Scalar t, const std::vector<Index>& elements) {
  if (!is_shell()) throw ModelError("only a shell model takes element thicknesses");
  if (finalized_) throw ModelError("thicknesses must be assigned before FemModel::finalize()");
  if (!(t > 0.0) || !std::isfinite(t)) {
    std::ostringstream os;
    os << "a shell thickness must be positive (got " << t << " m)";
    throw ConfigError(os.str());
  }
  const Index ne = mesh_.num_elements();
  if (element_thickness_.empty()) element_thickness_.assign(static_cast<std::size_t>(ne), thickness_);
  for (Index e : elements) {
    if (e < 0 || e >= ne) {
      std::ostringstream os;
      os << "thickness assigned to element " << e << ", outside [0, " << ne - 1 << "]";
      throw ModelError(os.str());
    }
    element_thickness_[static_cast<std::size_t>(e)] = t;
  }
}

void FemModel::set_shell_options(const ShellOptions& options) {
  if (!is_shell()) throw ModelError("shell options given to a model that is not a shell");
  if (finalized_) throw ModelError("shell options must be set before FemModel::finalize()");
  if (!(options.fold_angle_deg > 0.0) || !(options.fold_angle_deg < 90.0)) {
    std::ostringstream os;
    os << "the shell fold angle must lie in (0, 90) degrees (got " << options.fold_angle_deg
       << ")";
    throw ConfigError(os.str());
  }
  shell_options_ = options;
  element_ = std::make_unique<Shell4Element>(options.drilling_factor);
  compute_directors();
}

void FemModel::compute_directors() {
  // Per node, the element normals there are grouped: a normal joins the first
  // group whose mean lies within the fold angle of it (either way round, as
  // two elements of one surface may be numbered with opposite senses), and
  // each element takes its group's mean, turned to its own side.
  static const std::array<std::array<Scalar, 2>, 4> corners{
      {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}}};
  const Index ne = mesh_.num_elements();
  const Scalar cos_fold = std::cos(shell_options_.fold_angle_deg * 3.14159265358979323846 / 180.0);
  std::vector<Eigen::Matrix<Scalar, 3, 4>> own(static_cast<std::size_t>(ne));
  for (Index e = 0; e < ne; ++e) {
    const Matrix x = mesh_.element_coordinates(e);
    for (int k = 0; k < 4; ++k) {
      own[static_cast<std::size_t>(e)].col(k) = Shell4Element::normal(
          x, corners[static_cast<std::size_t>(k)][0], corners[static_cast<std::size_t>(k)][1]);
    }
  }
  struct Use {
    Index element;
    int corner;
  };
  std::vector<std::vector<Use>> at_node(static_cast<std::size_t>(mesh_.num_nodes()));
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = mesh_.element_nodes(e);
    for (int k = 0; k < 4; ++k) at_node[static_cast<std::size_t>(nodes[k])].push_back({e, k});
  }
  directors_ = own;
  if (mesh_.has_node_normals()) {
    // The surface's own normals, turned to each element's side. A facet of
    // a coarse mesh of a curved surface departs from the surface's normal by
    // half its angle (22.5 degrees for 8 cells round a circle), so no fold
    // test applies; a normal within 15 degrees of the element's plane,
    // though, tilts the element's fibres flat (its volume Jacobian goes to
    // zero at 90 degrees): the mesh does not follow those normals.
    const Scalar cos_limit = std::cos(75.0 * 3.14159265358979323846 / 180.0);
    const Matrix& normals = mesh_.node_normals();
    for (Index e = 0; e < ne; ++e) {
      const Index* nodes = mesh_.element_nodes(e);
      for (int k = 0; k < 4; ++k) {
        const Vector3 v = normals.col(nodes[k]);
        const Vector3 n = own[static_cast<std::size_t>(e)].col(k);
        if (std::abs(v.dot(n)) < cos_limit) {
          std::ostringstream os;
          os << "the normal the mesh gives at node " << nodes[k] << " lies "
             << std::acos(std::min(1.0, std::abs(v.dot(n)))) * 180.0 / 3.14159265358979323846
             << " degrees from the normal of element " << e
             << " there (at most 75 allowed); the mesh does not follow the surface its "
                "normals describe";
          throw MeshError(os.str());
        }
        directors_[static_cast<std::size_t>(e)].col(k) = v.dot(n) >= 0.0 ? v : Vector3(-v);
      }
    }
    return;
  }
  for (const std::vector<Use>& uses : at_node) {
    std::vector<Vector3> sums;
    std::vector<std::size_t> group(uses.size());
    for (std::size_t i = 0; i < uses.size(); ++i) {
      const Vector3 n =
          own[static_cast<std::size_t>(uses[i].element)].col(uses[i].corner);
      std::size_t g = 0;
      for (; g < sums.size(); ++g) {
        const Vector3 mean = sums[g].normalized();
        if (std::abs(mean.dot(n)) >= cos_fold) {
          sums[g] += mean.dot(n) >= 0.0 ? n : Vector3(-n);
          break;
        }
      }
      if (g == sums.size()) sums.push_back(n);
      group[i] = g;
    }
    for (std::size_t i = 0; i < uses.size(); ++i) {
      const Vector3 n =
          own[static_cast<std::size_t>(uses[i].element)].col(uses[i].corner);
      const Scalar len = sums[group[i]].norm();
      if (!(len > 1.0e-12)) continue;  // opposite normals cancel: keep the element's own
      const Vector3 mean = sums[group[i]] / len;
      directors_[static_cast<std::size_t>(uses[i].element)].col(uses[i].corner) =
          mean.dot(n) >= 0.0 ? mean : Vector3(-mean);
    }
  }
}

Matrix FemModel::element_geometry(Index e) const {
  if (directors_.empty()) return mesh_.element_coordinates(e);
  Matrix g(6, 4);
  g.topRows(3) = mesh_.element_coordinates(e);
  g.bottomRows(3) = directors_[static_cast<std::size_t>(e)];
  return g;
}

void FemModel::set_material(const IsotropicMaterial& material) {
  materials_.front() = material;
  d_.front() = material.constitutive(stress_state_);
}

void FemModel::assign_material(const IsotropicMaterial& material,
                               const std::vector<Index>& elements) {
  if (finalized_) {
    throw ModelError("materials must be assigned before FemModel::finalize()");
  }
  const Index ne = mesh_.num_elements();
  for (Index e : elements) {
    if (e < 0 || e >= ne) {
      std::ostringstream os;
      os << "material '" << material.name() << "' assigned to element " << e
         << ", outside [0, " << ne - 1 << "]";
      throw ModelError(os.str());
    }
  }
  materials_.push_back(material);
  d_.push_back(material.constitutive(stress_state_));
  if (element_material_.empty()) element_material_.assign(static_cast<std::size_t>(ne), 0);
  const int id = static_cast<int>(materials_.size()) - 1;
  for (Index e : elements) element_material_[static_cast<std::size_t>(e)] = id;
}

void FemModel::set_conduction_solver(const LinearSolverOptions& options) {
  conduction_solver_ = std::make_shared<LinearSolverOptions>(options);
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
  load_data_.clear();
  load_vectors_.reserve(load_case_specs_.size());
  load_data_.reserve(load_case_specs_.size());
  for (const LoadCaseSpec& spec : load_case_specs_) {
    if (!(spec.weight >= 0.0)) {
      std::ostringstream os;
      os << "load case '" << spec.name << "' has a negative weight (" << spec.weight
         << "); weights must be non-negative";
      throw ConfigError(os.str());
    }
    LoadCaseData data;
    if (is_shell()) {
      if (spec.has_temperature()) {
        throw ConfigError("load case '" + spec.name + "': a shell model takes no temperature "
                          "field - the shell element has no thermal strain");
      }
      if (spec.centrifugal.enabled) {
        throw ConfigError("load case '" + spec.name + "': a shell model takes no steady "
                          "rotation (its centrifugal load through the thickness is not "
                          "integrated)");
      }
      data.mechanical = assemble_shell_load_vector(*this, spec);
    } else {
      data.mechanical = assemble_load_vector(mesh_, *element_, spec, thickness_, integration_);
    }
    Vector total = data.mechanical;
    if (spec.has_body_loads()) {
      data.body = assemble_body_load_vector(*this, spec);
      for (Index n = 0; n < mesh_.num_nodes(); ++n) {
        for (int k = 0; k < mesh_.dim(); ++k) {
          data.body_resultant(k) += data.body(n * dofs_.dofs_per_node() + k);
        }
      }
      total += data.body;
    }
    if (spec.has_temperature()) {
      if (spec.temperature.source == TemperatureSpec::Source::Conduction) {
        LinearSolverOptions linear;
        if (conduction_solver_) linear = *conduction_solver_;
        const ConductionResult solved = solve_conduction(*this, spec.temperature.conduction, linear);
        data.temperature = solved.temperature;
        data.conduction = solved.summary;
        data.conduction_solved = true;
        log::info("load case '", spec.name, "': steady conduction, temperature ",
                  solved.summary.min_temperature, " to ", solved.summary.max_temperature,
                  " K, heat in ", solved.summary.applied_heat, " W, out through prescribed "
                  "temperatures ", solved.summary.prescribed_heat, " W (relative balance "
                  "error ", solved.summary.relative_balance_error, ")");
      } else {
        data.temperature = resolve_region_temperatures(mesh_, spec.temperature);
      }
      const ThermalLoad thermal = assemble_thermal_load(*this, data.temperature);
      data.thermal = thermal.force;
      data.thermal_self_energy = thermal.self_energy;
      total += data.thermal;
    }
    // A zero load vector is a mistake *unless* the case is driven by prescribed
    // displacements, which is how the patch test and any enforced-deflection
    // study work.
    if (total.norm() == 0.0 && !dofs_.has_nonzero_prescribed() &&
        !spec.prescribed_displacement_only) {
      log::warn("load case '", spec.name,
                "' has a zero resultant force vector and no non-zero prescribed "
                "displacement, so its solution is identically zero");
    }
    load_vectors_.push_back(std::move(total));
    load_data_.push_back(std::move(data));
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

const LoadCaseData& FemModel::load_case_data(std::size_t l) const {
  if (!finalized_) {
    throw ModelError("load-case data requested before FemModel::finalize() was called");
  }
  if (l >= load_data_.size()) {
    std::ostringstream os;
    os << "load case " << l << " requested, the model has " << load_data_.size();
    throw ModelError(os.str());
  }
  return load_data_[l];
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
  if (is_shell()) {
    for (Index e = 0; e < ne; ++e) v(e) = mesh_.element_measure(e) * thickness_of(e);
  } else if (mesh_.dim() == 2) {
    for (Index e = 0; e < ne; ++e) v(e) = mesh_.element_measure(e) * thickness_;
  } else {
    for (Index e = 0; e < ne; ++e) v(e) = mesh_.element_measure(e);
  }
  return v;
}

Scalar FemModel::domain_volume() const { return element_volumes().sum(); }

}  // namespace sparlab
