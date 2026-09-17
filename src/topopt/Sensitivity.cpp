#include "sparlab/topopt/Sensitivity.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/core/Timer.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {

ComplianceObjective::ComplianceObjective(const FemModel& model,
                                         const Assembler& assembler,
                                         const DensityFilter& filter,
                                         const DesignDomain& domain, SimpOptions simp,
                                         StaticAnalysisOptions analysis_options)
    : model_(model),
      assembler_(assembler),
      filter_(filter),
      domain_(domain),
      simp_(simp),
      analysis_options_(analysis_options),
      weights_(model.normalised_weights()) {
  simp_.validate();
  if (filter_.num_elements() != model_.mesh().num_elements()) {
    throw ConfigError("the filter was built for a different mesh than the model");
  }
  // The SIMP stiffness floor keeps K_ff positive definite everywhere, so the
  // rigid-body diagnostics are run once by the caller on the solid model
  // instead of on every SIMP iterate.
  analysis_options_.check_model = false;
}

ObjectiveEvaluation ComplianceObjective::evaluate(const Vector& x, bool need_gradients) {
  const Index ne = model_.mesh().num_elements();
  if (x.size() != ne) {
    std::ostringstream os;
    os << "design vector has length " << x.size() << " but the mesh has " << ne
       << " elements";
    throw ConfigError(os.str());
  }
  if (!x.allFinite()) throw ConfigError("design vector contains NaN or Inf");

  ObjectiveEvaluation out;
  out.physical_density = filter_.to_physical(x);
  // Filtering is a convex combination, so the result stays in [0,1] up to
  // round-off; clamp to keep pow() well defined.
  out.physical_density = out.physical_density.cwiseMax(0.0).cwiseMin(1.0);
  out.stiffness_factors = simp_stiffness_factors(out.physical_density, simp_);

  Timer timer;
  StaticAnalysis analysis(model_, assembler_, analysis_options_);
  analysis.prepare(&out.stiffness_factors);
  out.assemble_seconds = timer.elapsed_seconds();

  timer.reset();
  const std::vector<Vector>& loads = model_.load_vectors();
  out.displacements.reserve(loads.size());
  out.load_case_compliance.reserve(loads.size());

  const int npe = model_.mesh().nodes_per_elem();
  const int edofs = npe * kDofsPerNode;
  out.element_strain_energy.setZero(ne);
  Vector dc_dphys = Vector::Zero(ne);
  const Vector dfactor = need_gradients
                             ? simp_stiffness_derivatives(out.physical_density, simp_)
                             : Vector();

  Vector ue(edofs);
  for (std::size_t l = 0; l < loads.size(); ++l) {
    const Vector u = analysis.solve_load_vector(loads[l]);
    ++num_solves_;
    const Scalar c = loads[l].dot(u);
    out.load_case_compliance.push_back(c);
    out.compliance += weights_[l] * c;
    out.displacements.push_back(u);

    for (Index e = 0; e < ne; ++e) {
      const Index* nodes = model_.mesh().element_nodes(e);
      for (int a = 0; a < npe; ++a) {
        ue(kDofsPerNode * a + 0) = u(nodes[a] * kDofsPerNode + 0);
        ue(kDofsPerNode * a + 1) = u(nodes[a] * kDofsPerNode + 1);
      }
      // u_e^T K_e^0 u_e: the unit-density element strain energy times two.
      const Scalar quad = ue.dot(assembler_.element_stiffness(e) * ue);
      out.element_strain_energy(e) +=
          weights_[l] * 0.5 * out.stiffness_factors(e) * quad;
      if (need_gradients) dc_dphys(e) -= weights_[l] * dfactor(e) * quad;
    }
  }
  out.solve_seconds = timer.elapsed_seconds();

  out.volume = domain_.volume_of(out.physical_density);
  out.volume_fraction = out.volume / domain_.domain_volume();

  if (need_gradients) {
    out.dc_dphysical = dc_dphys;
    out.dc_dx = filter_.transform_gradient(x, dc_dphys);
    out.dv_dx = filter_.pull_back(domain_.element_volumes());

    if (out.dc_dx.maxCoeff() > 0.0) {
      // For minimum compliance the gradient must be non-positive everywhere.
      // A positive entry indicates a bug or a pathological filter; report it
      // rather than clipping silently.
      log::warn("compliance gradient has a positive entry (max = ",
                out.dc_dx.maxCoeff(),
                "), which is not expected for minimum-compliance SIMP; the optimality "
                "criteria update will clamp it");
    }
  }

  if (!std::isfinite(out.compliance)) {
    throw SolverError("the compliance objective evaluated to a non-finite value");
  }
  return out;
}

Scalar ComplianceObjective::compliance_at(const Vector& x) {
  return evaluate(x, /*need_gradients=*/false).compliance;
}

SensitivityCheckResult verify_sensitivities(ComplianceObjective& objective,
                                            const DesignDomain& domain, const Vector& x,
                                            const std::vector<Index>& elements,
                                            Scalar step, Scalar tolerance) {
  if (!(step > 0.0)) {
    throw ConfigError("finite-difference step must be positive");
  }
  SensitivityCheckResult result;
  result.step = step;
  result.tolerance = tolerance;

  const ObjectiveEvaluation base = objective.evaluate(x, /*need_gradients=*/true);
  result.gradient_infinity_norm = base.dc_dx.cwiseAbs().maxCoeff();

  std::vector<Index> targets = elements;
  if (targets.empty()) {
    targets.resize(static_cast<std::size_t>(domain.num_elements()));
    for (Index e = 0; e < domain.num_elements(); ++e) {
      targets[static_cast<std::size_t>(e)] = e;
    }
  }

  Scalar sum_sq_rel = 0.0;
  Vector xp = x;
  Vector xm = x;
  for (Index e : targets) {
    SensitivityCheckEntry entry;
    entry.element = e;
    entry.analytical = base.dc_dx(e);

    const Scalar lo = domain.lower_bounds()(e);
    const Scalar hi = domain.upper_bounds()(e);
    if (!domain.is_free(e)) {
      entry.excluded = true;
      entry.exclusion_reason =
          "passive variable (lower bound equals upper bound), so no admissible "
          "perturbation exists";
    } else if (x(e) - step < lo || x(e) + step > hi) {
      std::ostringstream os;
      os << "design variable " << x(e) << " is within the step " << step
         << " of its bounds [" << lo << ", " << hi
         << "], so a central difference would leave the feasible box";
      entry.excluded = true;
      entry.exclusion_reason = os.str();
    }

    if (entry.excluded) {
      ++result.num_excluded;
      result.entries.push_back(std::move(entry));
      continue;
    }

    xp(e) = x(e) + step;
    xm(e) = x(e) - step;
    const Scalar cp = objective.compliance_at(xp);
    const Scalar cm = objective.compliance_at(xm);
    xp(e) = x(e);
    xm(e) = x(e);

    entry.finite_difference = (cp - cm) / (2.0 * step);
    entry.absolute_error = std::abs(entry.analytical - entry.finite_difference);
    const Scalar scale =
        std::max({std::abs(entry.analytical), std::abs(entry.finite_difference),
                  1.0e-30});
    entry.relative_error = entry.absolute_error / scale;

    result.max_absolute_error = std::max(result.max_absolute_error, entry.absolute_error);
    result.max_relative_error = std::max(result.max_relative_error, entry.relative_error);
    sum_sq_rel += entry.relative_error * entry.relative_error;
    ++result.num_tested;
    result.entries.push_back(std::move(entry));
  }

  if (result.num_tested == 0) {
    throw ConfigError(
        "every candidate element was excluded from the finite-difference check; move "
        "the design point away from its bounds or reduce the step size");
  }
  result.rms_relative_error =
      std::sqrt(sum_sq_rel / static_cast<Scalar>(result.num_tested));

  // Directional derivative along the gradient restricted to tested elements.
  {
    Vector direction = Vector::Zero(x.size());
    for (const SensitivityCheckEntry& entry : result.entries) {
      if (!entry.excluded) direction(entry.element) = entry.analytical;
    }
    const Scalar dnorm = direction.norm();
    if (dnorm > 0.0) {
      direction /= dnorm;
      // Shrink the step so no component leaves the box.
      Scalar h = step;
      for (Index e = 0; e < domain.num_elements(); ++e) {
        if (direction(e) == 0.0) continue;
        const Scalar room = direction(e) > 0.0
                                ? (domain.upper_bounds()(e) - x(e)) / direction(e)
                                : (domain.lower_bounds()(e) - x(e)) / direction(e);
        h = std::min(h, std::abs(room));
      }
      if (h > 0.0) {
        const Scalar cp = objective.compliance_at((x + h * direction).eval());
        const Scalar cm = objective.compliance_at((x - h * direction).eval());
        const Scalar fd = (cp - cm) / (2.0 * h);
        const Scalar an = base.dc_dx.dot(direction);
        result.directional_relative_error =
            std::abs(fd - an) / std::max({std::abs(fd), std::abs(an), 1.0e-30});
      }
    }
  }

  result.passed = result.max_relative_error <= tolerance;
  log::info("sensitivity check: step ", step, ", ", result.num_tested, " elements (",
            result.num_excluded, " excluded), max relative error ",
            result.max_relative_error, ", RMS ", result.rms_relative_error,
            ", directional ", result.directional_relative_error, " -> ",
            (result.passed ? "PASS" : "FAIL"), " at tolerance ", tolerance);
  return result;
}

}  // namespace sparlab
