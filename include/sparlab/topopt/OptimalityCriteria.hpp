/// \file OptimalityCriteria.hpp
/// \brief Optimality-criteria update with bisection on the volume multiplier.
///
/// Minimum-compliance SIMP with a single linear volume constraint has a
/// separable KKT system whose stationarity condition is
/// \f[
///   B_e \equiv \frac{-\partial c/\partial x_e}
///                   {\lambda\, \partial g/\partial x_e} = 1
///   \quad\text{for every variable strictly inside its bounds.}
/// \f]
/// The classical fixed-point update (Bendsoe & Sigmund) is therefore
/// \f[
///   x_e^{k+1} = \mathrm{clip}\!\left(x_e^k\, B_e^{\eta},\;
///     \max(\underline{x}_e,\, x_e^k - m),\;
///     \min(\overline{x}_e,\, x_e^k + m)\right),
/// \f]
/// with damping \f$\eta = 1/2\f$ and move limit \f$m\f$. \f$\lambda\f$ is found
/// by bisection so the volume constraint is met to a tight tolerance; because
/// the mapped volume is monotonically decreasing in \f$\lambda\f$, bisection is
/// globally convergent.
///
/// Why OC rather than MMA: the problem has one constraint, the objective is
/// separable in the SIMP sense, and OC converges in tens of iterations with no
/// tuning. Its limitation - only one inequality constraint - is stated in
/// docs/limitations.md.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/topopt/DesignDomain.hpp"

#include <functional>

namespace sparlab {

struct OptimalityCriteriaOptions {
  Scalar move_limit = 0.2;       ///< m, maximum per-iteration density change [-]
  Scalar damping = 0.5;          ///< eta [-]
  Scalar lambda_lower = 1.0e-12; ///< initial bisection bracket lower bound
  Scalar lambda_upper = 1.0e12;  ///< initial bisection bracket upper bound
  int max_bisections = 200;      ///< iteration cap for the multiplier search
  /// Relative volume tolerance of the bisection: |V(x)-V*| <= tol * V*.
  Scalar volume_tolerance = 1.0e-10;
};

struct OptimalityCriteriaStep {
  Vector x;                      ///< updated design variables
  Scalar lambda = 0.0;           ///< converged Lagrange multiplier
  Scalar achieved_volume = 0.0;  ///< \f$\tilde\rho^T v\f$ after the update [m^3]
  Scalar max_change = 0.0;       ///< \f$\|x^{k+1}-x^k\|_\infty\f$
  int bisections = 0;
  bool volume_converged = false;
};

/// One optimality-criteria update.
/// \param domain design bounds, passive tags and element volumes.
/// \param x current design variables.
/// \param dc_dx objective gradient (expected non-positive).
/// \param dv_dx volume-constraint gradient (positive).
/// \param physical_volume_of maps design variables to physical volume; this is
///        `DesignDomain::volume_of(filter.to_physical(x))` and is passed as a
///        callable so the update stays independent of the filter type.
/// \throws ConvergenceError when the bisection cannot bracket the target
///         volume, which means the bounds make the constraint unreachable.
OptimalityCriteriaStep optimality_criteria_update(
    const DesignDomain& domain, const Vector& x, const Vector& dc_dx,
    const Vector& dv_dx, const std::function<Scalar(const Vector&)>& physical_volume_of,
    const OptimalityCriteriaOptions& options);

}  // namespace sparlab
