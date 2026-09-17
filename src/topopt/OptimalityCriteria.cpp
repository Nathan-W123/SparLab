#include "sparlab/topopt/OptimalityCriteria.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

/// Candidate design for a trial multiplier.
Vector oc_candidate(const DesignDomain& domain, const Vector& x, const Vector& dc_dx,
                    const Vector& dv_dx, Scalar lambda,
                    const OptimalityCriteriaOptions& options) {
  const Index ne = domain.num_elements();
  Vector xnew(ne);
  for (Index e = 0; e < ne; ++e) {
    const Scalar lo = std::max(domain.lower_bounds()(e), x(e) - options.move_limit);
    const Scalar hi = std::min(domain.upper_bounds()(e), x(e) + options.move_limit);
    if (!domain.is_free(e)) {
      xnew(e) = domain.lower_bounds()(e);
      continue;
    }
    // B_e = -dc/dx / (lambda dv/dx). Clamp at zero: a non-negative objective
    // gradient carries no information about where to add material, and a
    // negative base would make pow() undefined.
    const Scalar denom = lambda * dv_dx(e);
    Scalar b = 0.0;
    if (denom > 0.0) b = std::max(0.0, -dc_dx(e) / denom);
    const Scalar target = x(e) * std::pow(b, options.damping);
    xnew(e) = std::clamp(target, lo, hi);
  }
  return xnew;
}

}  // namespace

OptimalityCriteriaStep optimality_criteria_update(
    const DesignDomain& domain, const Vector& x, const Vector& dc_dx,
    const Vector& dv_dx, const std::function<Scalar(const Vector&)>& physical_volume_of,
    const OptimalityCriteriaOptions& options) {
  const Index ne = domain.num_elements();
  if (x.size() != ne || dc_dx.size() != ne || dv_dx.size() != ne) {
    throw ConfigError("optimality criteria received vectors of inconsistent length");
  }
  if (!(options.move_limit > 0.0 && options.move_limit <= 1.0)) {
    std::ostringstream os;
    os << "optimizer.move_limit must lie in (0, 1] (got " << options.move_limit << ")";
    throw ConfigError(os.str());
  }
  if (!(options.damping > 0.0 && options.damping <= 1.0)) {
    std::ostringstream os;
    os << "optimizer.damping must lie in (0, 1] (got " << options.damping << ")";
    throw ConfigError(os.str());
  }
  if (!dc_dx.allFinite() || !dv_dx.allFinite()) {
    throw SolverError("optimality criteria received a non-finite gradient");
  }

  const Scalar target = domain.volume_target();

  // Volume is monotonically non-increasing in lambda, so widen the bracket
  // until it straddles the target.
  Scalar lo = options.lambda_lower;
  Scalar hi = options.lambda_upper;
  Scalar vol_lo = physical_volume_of(oc_candidate(domain, x, dc_dx, dv_dx, lo, options));
  Scalar vol_hi = physical_volume_of(oc_candidate(domain, x, dc_dx, dv_dx, hi, options));
  int widen = 0;
  while (vol_lo < target && widen < 60) {
    lo *= 1.0e-3;
    vol_lo = physical_volume_of(oc_candidate(domain, x, dc_dx, dv_dx, lo, options));
    ++widen;
  }
  while (vol_hi > target && widen < 120) {
    hi *= 1.0e3;
    vol_hi = physical_volume_of(oc_candidate(domain, x, dc_dx, dv_dx, hi, options));
    ++widen;
  }
  if (vol_lo < target || vol_hi > target) {
    std::ostringstream os;
    os << "the optimality-criteria bisection cannot bracket the volume target "
       << target << " m^3: the move-limited box yields volumes in [" << vol_hi << ", "
       << vol_lo << "] m^3. This happens when the move limit is too small to reach "
          "the target from the current design, or when passive regions make the target "
          "unreachable";
    throw ConvergenceError(os.str());
  }

  OptimalityCriteriaStep step;
  Vector candidate = oc_candidate(domain, x, dc_dx, dv_dx, lo, options);
  Scalar volume = vol_lo;
  for (int it = 0; it < options.max_bisections; ++it) {
    const Scalar mid = std::sqrt(lo * hi);  // geometric bisection: lambda spans decades
    candidate = oc_candidate(domain, x, dc_dx, dv_dx, mid, options);
    volume = physical_volume_of(candidate);
    step.bisections = it + 1;
    step.lambda = mid;
    if (std::abs(volume - target) <= options.volume_tolerance * target) {
      step.volume_converged = true;
      break;
    }
    if (volume > target) {
      lo = mid;
    } else {
      hi = mid;
    }
    if (hi / lo - 1.0 < 1.0e-15) {
      // The bracket collapsed: the volume is a step function of lambda here
      // (many variables pinned at bounds). Accept the current candidate and
      // report the achieved volume so the caller can see the gap.
      step.volume_converged =
          std::abs(volume - target) <= 1.0e-6 * std::max(target, 1.0e-30);
      break;
    }
  }

  step.achieved_volume = volume;
  step.max_change = (candidate - x).cwiseAbs().maxCoeff();
  step.x = std::move(candidate);

  if (!step.volume_converged) {
    log::warn("optimality-criteria bisection stopped after ", step.bisections,
              " steps with volume ", step.achieved_volume, " m^3 versus target ",
              target, " m^3 (relative gap ",
              std::abs(step.achieved_volume - target) / std::max(target, 1.0e-30), ")");
  }
  return step;
}

}  // namespace sparlab
