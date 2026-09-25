/// \file BucklingConstraint.hpp
/// \brief Aggregated lower bound on the linear buckling load factors of a
///        SIMP design, with adjoint sensitivities (MMA only).
///
/// For each constrained load case the lowest `num_modes` positive load
/// factors \f$\lambda_i\f$ of
/// \f$(K(\rho) + \lambda K_G(\rho, u))\phi_i = 0\f$ (Buckling.hpp) must stay
/// above \f$\lambda_{req}\f$. With \f$r_i = \lambda_{req}/\lambda_i\f$ the
/// constraint is the Kreisselmeier-Steinhauser aggregate
/// \f[
///   g = \mathrm{KS}_P(r) - 1 = r_{max} + \frac{1}{P}\ln\sum_i e^{P(r_i - r_{max})}
///       - 1 \le 0,
/// \f]
/// a smooth upper bound on \f$\max_i r_i - 1\f$ that exceeds it by at most
/// \f$\ln(m)/P\f$: the constraint is conservative, and it stays
/// differentiable when two load factors cross or coincide, where a single
/// \f$\lambda_1\f$ is not.
///
/// **Interpolation.** The stiffness is the SIMP law
/// \f$E_K(\rho) = E_{min} + \rho^p(E_0 - E_{min})\f$; the stress entering
/// \f$K_G\f$ uses \f$E_G(\rho) = \rho^p E_0\f$ without the floor (Gao and Ma
/// 2015). With the same law in both, a void element's geometric stiffness
/// keeps pace with its elastic stiffness, and near-void regions produce
/// spurious localised "pseudo" buckling modes at very low load factors. The
/// fraction of each mode's strain energy in elements with
/// \f$\rho \ge 0.5\f$ is recorded so any mode that still localises in void
/// is visible.
///
/// **Sensitivity.** With \f$\phi^T K \phi = 1\f$, so
/// \f$\phi^T K_G \phi = -1/\lambda\f$,
/// \f[
///   \frac{d\lambda}{d\rho_e} =
///     \lambda \frac{E_K'}{E_0}\phi_e^T K_e^0 \phi_e
///   + \lambda^2 \frac{E_G'}{E_0}\phi_e^T K_{G,e}^0(u_e)\phi_e
///   - \lambda^2 \frac{E_K'}{E_0} a_e^T K_e^0 u_e,
///   \qquad K a = \sum_e \frac{E_G}{E_0} A_e^T g_e(\phi_e),
/// \f]
/// where \f$g_e\f$ is `Element::geometric_stiffness_derivative` (the
/// derivative of \f$\phi_e^T K_{G,e}^0\phi_e\f$ with respect to \f$u_e\f$)
/// and the adjoint solve accounts for the dependence of the stress on the
/// design through \f$u\f$. Each mode costs one adjoint solve with the
/// factorisation already in hand. The formula holds for a simple eigenvalue;
/// at a repeated one the KS aggregate still has a gradient, which the
/// finite-difference check in the tests confirms away from crossings.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/topopt/Sensitivity.hpp"

#include <string>
#include <vector>

namespace sparlab {

struct BucklingConstraintOptions {
  bool enabled = false;
  Scalar min_load_factor = 1.0;    ///< lambda_req [-]
  int num_modes = 6;               ///< load factors aggregated per load case
  Scalar ks_parameter = 40.0;      ///< KS aggregation parameter P
  /// Load cases to constrain, by name; empty constrains every case.
  std::vector<std::string> load_cases;
  /// Density at and above which an element counts as solid in the
  /// pseudo-mode check.
  Scalar solid_threshold = 0.5;
  /// Subspace-iteration settings of the eigensolve (num_modes is set from
  /// the option above).
  BucklingOptions eigen;

  void validate() const;
};

struct BucklingEvaluation {
  Vector load_factors;             ///< smallest positive lambda_i, ascending [-]
  Vector solid_energy_fraction;    ///< per mode
  Scalar ks = 0.0;                 ///< KS of lambda_req / lambda_i
  Scalar constraint = 0.0;         ///< ks - 1
  bool no_positive_load_factor = false;
  int iterations = 0;
  bool transformed = false;        ///< the eigensolve used the spectral transformation
  Scalar sigma = 0.0;              ///< its shift
  Vector dg_dphysical;             ///< d(constraint)/d rho (when requested)
  Vector dg_dx;                    ///< d(constraint)/d x through filter and projection
  /// Per mode, d lambda_i / d rho (when requested); used by the tests.
  std::vector<Vector> dlambda_dphysical;
};

/// Stress interpolation factor \f$E_G(\rho)/E_0 = \rho^p\f$ and its derivative.
/// \{
Vector buckling_stress_factors(const Vector& rho, Scalar penalty);
Vector buckling_stress_derivatives(const Vector& rho, Scalar penalty);
/// \}

class BucklingConstraint {
 public:
  BucklingConstraint(const FemModel& model, const Assembler& assembler,
                     BucklingConstraintOptions options);

  const BucklingConstraintOptions& options() const { return options_; }

  /// Indices of the constrained load cases.
  const std::vector<std::size_t>& load_cases() const { return cases_; }

  /// Evaluate the constraint of load case `load_case` at the design behind
  /// `eval`, which must be the latest `objective.evaluate` (its factorisation
  /// is reused). The eigensolve starts from the subspace of the previous call
  /// for the same case, which keeps an optimisation iteration to a few
  /// subspace iterations.
  BucklingEvaluation evaluate(ComplianceObjective& objective, const ObjectiveEvaluation& eval,
                              std::size_t load_case, bool need_gradients);

 private:
  const FemModel& model_;
  const Assembler& assembler_;
  BucklingConstraintOptions options_;
  std::vector<std::size_t> cases_;
  std::vector<Matrix> subspace_;   ///< warm start per load case
};

}  // namespace sparlab
