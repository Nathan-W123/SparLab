/// \file Sensitivity.hpp
/// \brief Weighted-compliance objective, its analytical sensitivities, and the
///        finite-difference verification of those sensitivities.
///
/// **Objective.** For load cases \f$l\f$ with normalised weights \f$w_l\f$,
/// \f[
///   c(x) = \sum_l w_l\, f_l^T u_l, \qquad K(\tilde\rho)\, u_l = f_l .
/// \f]
///
/// **Sensitivity.** Compliance is self-adjoint, so differentiating
/// \f$K u = f\f$ with \f$f\f$ independent of the design gives
/// \f$ \partial u / \partial \tilde\rho_e = -K^{-1} (\partial K/\partial
/// \tilde\rho_e) u \f$ and therefore
/// \f[
///   \frac{\partial c}{\partial \tilde\rho_e}
///     = -\sum_l w_l\, u_{l,e}^T
///       \frac{\partial K_e}{\partial \tilde\rho_e} u_{l,e}
///     = -\sum_l w_l\,
///       \frac{\mathrm{d}}{\mathrm{d}\tilde\rho_e}
///       \!\left(\frac{E(\tilde\rho_e)}{E_0}\right)
///       u_{l,e}^T K_e^0 u_{l,e},
/// \f]
/// which is non-positive: adding material never increases compliance. The
/// design-variable gradient follows from the filter chain rule,
/// \f$ \nabla_x c = \hat{H}^T \nabla_{\tilde\rho} c \f$.
///
/// **Volume constraint.** \f$ g(x) = \tilde\rho^T v - \nu V \f$ with
/// \f$ \nabla_x g = \hat{H}^T v \f$.
///
/// No adjoint solve is required, so one linear solve per load case yields both
/// the objective and the exact gradient.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/Projection.hpp"
#include "sparlab/topopt/SimpInterpolation.hpp"

#include <memory>
#include <vector>

namespace sparlab {

/// Objective, constraint and gradients at one design point.
struct ObjectiveEvaluation {
  /// Physical density [-]: the filtered density, projected when the Heaviside
  /// projection is on (\f$\bar\rho\f$), otherwise \f$\tilde\rho\f$ itself.
  Vector physical_density;
  Vector filtered_density;        ///< \f$\tilde\rho\f$ [-]
  /// \f$d\bar\rho/d\tilde\rho\f$; empty when the projection is off.
  Vector projection_derivative;
  Vector stiffness_factors;       ///< \f$E(\tilde\rho)/E_0\f$ [-]
  Scalar compliance = 0.0;        ///< weighted compliance [J]
  std::vector<Scalar> load_case_compliance;  ///< per load case [J]
  Scalar volume = 0.0;            ///< \f$\tilde\rho^T v\f$ [m^3]
  Scalar volume_fraction = 0.0;   ///< volume / domain volume [-]
  Vector dc_dphysical;            ///< \f$\partial c/\partial\tilde\rho\f$ [J]
  Vector dc_dx;                   ///< \f$\partial c/\partial x\f$ [J]
  Vector dv_dx;                   ///< \f$\partial g/\partial x\f$ [m^3]
  Vector element_strain_energy;   ///< per element, summed over weighted cases [J]
  std::vector<Vector> displacements;  ///< one per load case [m]
  Scalar max_scaled_residual = 0.0;
  Scalar solve_seconds = 0.0;
  Scalar assemble_seconds = 0.0;
  /// Iterations of the iterative linear solver summed over the load cases
  /// (0 for a direct solver).
  int linear_iterations = 0;
};

/// Evaluates the weighted compliance objective and its gradients.
class ComplianceObjective {
 public:
  ComplianceObjective(const FemModel& model, const Assembler& assembler,
                      const DensityFilter& filter, const DesignDomain& domain,
                      SimpOptions simp, StaticAnalysisOptions analysis_options);

  /// Objective + exact gradients at design point `x`.
  ObjectiveEvaluation evaluate(const Vector& x, bool need_gradients = true);

  /// Objective only (no gradients). Used by the finite-difference check.
  Scalar compliance_at(const Vector& x);

  /// Current SIMP options (the optimizer mutates the penalty during
  /// continuation).
  SimpOptions& simp() { return simp_; }
  const SimpOptions& simp() const { return simp_; }

  /// Switch the Heaviside projection on with sharpness `beta` and threshold
  /// `eta` (the optimiser raises beta during continuation), or off.
  /// \throws ConfigError with the sensitivity filter, which has no chain rule.
  void set_projection(Scalar beta, Scalar eta);
  void disable_projection() { projection_ = false; }
  bool projection() const { return projection_; }
  Scalar projection_beta() const { return beta_; }
  Scalar projection_eta() const { return eta_; }

  /// Physical density of a design: filtered, clamped to [0, 1], projected.
  Vector physical_density(const Vector& x) const;

  /// Chain a derivative with respect to the physical density back to the
  /// design variables through the projection and the filter, using the
  /// state of evaluation `eval`.
  Vector chain_to_design(const ObjectiveEvaluation& eval, const Vector& d_dphysical) const;

  Index num_solves() const { return num_solves_; }

  /// Solve \f$K(\tilde\rho)\,\lambda = r\f$ with the factorisation of the last
  /// `evaluate` call (homogeneous conditions at prescribed DOFs). This is the
  /// adjoint solve behind every non-self-adjoint sensitivity, e.g. the stress
  /// constraint; it costs one back-substitution.
  /// \throws ModelError before the first evaluation.
  /// `slot` >= 0 names a recurring adjoint problem (one per load case, say):
  /// its previous solution seeds an iterative solver.
  Vector solve_adjoint(const Vector& rhs, int slot = -1);

  /// Iterations of the iterative solver in the adjoint solves since the last
  /// `evaluate` call.
  int adjoint_iterations() const { return adjoint_iterations_; }

  /// Iterative-solver iterations of every solve so far (0 for a direct solver).
  Index total_linear_iterations() const { return total_linear_iterations_; }

  /// The linear solver kept across design iterations (its multigrid
  /// hierarchy and statistics), or nullptr before the first evaluation.
  const LinearSolver* solver() const { return solver_.get(); }

  /// The global stiffness matrix \f$K(\tilde\rho)\f$ of the last `evaluate`
  /// call (full size), e.g. for a buckling eigenproblem on the same design.
  /// \throws ModelError before the first evaluation.
  const SparseMatrix& stiffness() const;

  const FemModel& model() const { return model_; }
  const Assembler& assembler() const { return assembler_; }
  const DensityFilter& filter() const { return filter_; }

 private:
  const FemModel& model_;
  const Assembler& assembler_;
  const DensityFilter& filter_;
  const DesignDomain& domain_;
  SimpOptions simp_;
  StaticAnalysisOptions analysis_options_;
  std::vector<Scalar> weights_;
  bool projection_ = false;
  Scalar beta_ = 1.0;
  Scalar eta_ = 0.5;
  Index num_solves_ = 0;
  /// One solver for the whole optimisation (a multigrid hierarchy is reused;
  /// a direct factorisation is simply recomputed).
  std::unique_ptr<LinearSolver> solver_;
  /// Factorisation of the last evaluated design, kept for adjoint solves.
  std::unique_ptr<StaticAnalysis> analysis_;
  /// Last solutions, the initial guesses of the next iterative solves.
  std::vector<Vector> previous_displacements_;
  std::vector<Vector> previous_adjoints_;
  int adjoint_iterations_ = 0;
  Index total_linear_iterations_ = 0;
};

/// Outcome of comparing an analytical gradient entry with central differences.
struct SensitivityCheckEntry {
  Index element = -1;
  Scalar analytical = 0.0;
  Scalar finite_difference = 0.0;
  Scalar absolute_error = 0.0;
  Scalar relative_error = 0.0;
  bool excluded = false;
  std::string exclusion_reason;
};

struct SensitivityCheckResult {
  Scalar step = 0.0;                  ///< perturbation h used [-]
  Index num_tested = 0;
  Index num_excluded = 0;
  Scalar max_absolute_error = 0.0;
  Scalar max_relative_error = 0.0;
  /// Largest error with each entry judged against max(|analytical|, |FD|,
  /// 1e-3 ||gradient||_inf): entries far below the gradient's own scale -
  /// numerous under a sharp projection, whose derivative is ~0 away from the
  /// threshold - are compared with the absolute round-off floor of a central
  /// difference rather than with themselves.
  Scalar max_scaled_error = 0.0;
  Scalar rms_relative_error = 0.0;
  Scalar gradient_infinity_norm = 0.0;
  /// Relative error of the directional derivative along the full gradient,
  /// which is the most sensitive single scalar test.
  Scalar directional_relative_error = 0.0;
  std::vector<SensitivityCheckEntry> entries;
  bool passed = false;
  Scalar tolerance = 0.0;
};

/// Compare \f$\partial c/\partial x\f$ against central differences
/// \f[
///   \frac{\partial c}{\partial x_e} \approx
///   \frac{c(x + h e_e) - c(x - h e_e)}{2h}.
/// \f]
/// Elements are excluded when a perturbation would cross a bound (passive
/// variables, or free variables within `h` of 0 or 1), because the one-sided
/// value is then not a valid central difference. Every exclusion is reported.
///
/// \param elements elements to test; empty means "every element".
/// \param step perturbation \f$h\f$ on the design variable [-].
/// \param tolerance pass threshold on the maximum relative error.
SensitivityCheckResult verify_sensitivities(ComplianceObjective& objective,
                                            const DesignDomain& domain, const Vector& x,
                                            const std::vector<Index>& elements,
                                            Scalar step, Scalar tolerance);

}  // namespace sparlab
