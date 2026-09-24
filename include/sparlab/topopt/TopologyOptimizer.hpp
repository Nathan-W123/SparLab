/// \file TopologyOptimizer.hpp
/// \brief SIMP minimum-compliance optimisation loop, with optimality criteria
///        or MMA as the update, and optional aggregated stress constraints.
///
/// Each iteration performs exactly one objective/gradient evaluation (one
/// factorisation of \f$K_{ff}\f$ shared by all load cases) followed by one
/// design update:
///   * **optimality criteria** (default) for the single volume constraint,
///     bit-for-bit the classic fixed-point update with bisection on the
///     multiplier;
///   * **MMA** for any number of constraints. The volume constraint is
///     \f$V/V^* - 1 \le 0\f$ and, when enabled, one aggregated stress
///     constraint per load case \f$c_l\, g_{PN,l} - 1 \le 0\f$ (see
///     StressConstraint.hpp), each needing one extra adjoint solve with the
///     factorisation already in hand. The objective is scaled by the initial
///     compliance so MMA sees O(1) numbers.
///
/// **Convergence.** Two independent indicators are tracked and *either* can
/// declare convergence:
///   * the design change \f$\|x^{k+1}-x^k\|_\infty\f$ falls below
///     `change_tolerance`, or
///   * `objective_tolerance` is enabled and the relative compliance change over
///     the last `objective_window` iterations falls below it.
/// With MMA the design must also be feasible: no constraint above
/// `constraint_tolerance`. The second criterion matters on fine meshes, where
/// a handful of elements can keep oscillating between bounds long after the
/// objective has settled; a design-change-only rule would then never
/// terminate. The result records which criterion fired and the final value of
/// both, so a reader can see whether the design or only the objective had
/// settled. Reaching `max_iterations` without meeting either criterion is
/// reported as non-convergence: the history and the final design are still
/// returned, `converged` is false and a warning names the stalled quantity.
/// Nothing is hidden.
///
/// **Continuation.** The SIMP penalty can be ramped from `penalty_start` to
/// `simp.penalty` in `continuation_steps` stages, each lasting
/// `continuation_iterations` iterations or until the stage converges. Starting
/// from a low penalty makes the early iterations nearly convex and reduces the
/// dependence of the final topology on the starting design; the penalty history
/// is recorded per iteration.
///
/// **Monotonicity.** Compliance is *not* guaranteed to decrease monotonically:
/// neither update is a line-search method, and each continuation step raises
/// the penalty, which raises compliance at fixed density. The history records
/// every value so the behaviour can be inspected rather than asserted.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/Mma.hpp"
#include "sparlab/topopt/OptimalityCriteria.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/SimpInterpolation.hpp"
#include "sparlab/topopt/StressConstraint.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Design update used by the loop.
enum class OptimizerMethod {
  OptimalityCriteria,  ///< single volume constraint, bisection on the multiplier
  MMA                  ///< method of moving asymptotes, any number of constraints
};

std::string to_string(OptimizerMethod method);
OptimizerMethod parse_optimizer_method(const std::string& text);

struct TopologyOptimizerOptions {
  SimpOptions simp;
  OptimalityCriteriaOptions oc;
  MmaOptions mma;
  StressConstraintOptions stress;
  StaticAnalysisOptions analysis;
  OptimizerMethod method = OptimizerMethod::OptimalityCriteria;

  int max_iterations = 200;
  /// Stop when max |dx| <= this value.
  Scalar change_tolerance = 1.0e-3;
  /// Alternative stop criterion on the relative compliance change over
  /// `objective_window` iterations; <= 0 disables it.
  Scalar objective_tolerance = 0.0;
  /// Window (in iterations) over which the relative compliance change is
  /// measured.
  int objective_window = 5;
  /// MMA only: largest constraint value (V/V* - 1, c g_PN - 1) accepted as
  /// feasible when declaring convergence.
  Scalar constraint_tolerance = 1.0e-4;

  /// Continuation on the SIMP penalty. `continuation_steps <= 1` disables it.
  int continuation_steps = 1;
  Scalar penalty_start = 1.0;
  int continuation_iterations = 25;

  /// Record the full density field every `history_stride` iterations for the
  /// evolution animation (1 = every iteration, 0 = never).
  int history_stride = 1;

  /// Threshold used when interpreting the final density field as solid
  /// geometry (reported, never silently applied to the optimisation itself).
  Scalar interpretation_threshold = 0.5;
};

/// One iteration of the optimisation history.
struct TopologyIteration {
  int iteration = 0;
  Scalar penalty = 0.0;
  Scalar compliance = 0.0;
  Scalar volume = 0.0;
  Scalar volume_fraction = 0.0;
  Scalar max_change = 0.0;
  Scalar lambda = 0.0;      ///< volume-constraint multiplier (OC bisection or MMA)
  Scalar gray_level = 0.0;  ///< 4/n * sum rho(1-rho): 0 = black/white, 1 = all grey
  Scalar seconds = 0.0;
  int bisections = 0;       ///< OC bisections or MMA subproblem Newton iterations
  bool volume_converged = false;
  // MMA-only records (zero for OC runs).
  Scalar max_stress_ratio = 0.0;      ///< max over cases of max_e sigma_rel / sigma_lim
  Scalar stress_constraint = 0.0;     ///< max over cases of c g_PN - 1
  Scalar constraint_violation = 0.0;  ///< max over all constraints, > 0 = violated
};

/// Stress state of one load case at the final design.
struct StressConstraintRecord {
  std::string load_case;
  Scalar max_relaxed_stress_Pa = 0.0;
  Scalar max_relaxed_ratio = 0.0;     ///< max_e rho^q sigma_vm / sigma_lim
  Index max_element = -1;
  Scalar p_norm_ratio = 0.0;
  Scalar scale = 1.0;
  Scalar constraint = 0.0;            ///< c g_PN - 1 at the final design
  /// Largest unrelaxed solid-material stress over the elements at or above the
  /// interpretation threshold, divided by the limit: what the constraint does
  /// not bound but a designer wants to know.
  Scalar max_solid_ratio_retained = 0.0;
};

struct TopologyOptimizationResult {
  Vector design;                 ///< final design variables x
  Vector physical_density;       ///< final filtered density
  Vector stiffness_factors;      ///< final E(rho)/E_0
  Scalar compliance = 0.0;       ///< final weighted compliance [J]
  std::vector<Scalar> load_case_compliance;
  Scalar volume = 0.0;           ///< final physical volume [m^3]
  Scalar volume_fraction = 0.0;
  Scalar volume_constraint_violation = 0.0;  ///< (V - V*) / V*
  Scalar gray_level = 0.0;
  bool converged = false;
  /// Which criterion ended the loop: "design_change", "objective_stall" or
  /// "iteration_cap".
  std::string stop_reason = "iteration_cap";
  /// Final values of both convergence indicators.
  Scalar final_design_change = 0.0;
  Scalar final_objective_change = 0.0;
  int iterations = 0;
  Scalar total_seconds = 0.0;
  Index linear_solves = 0;
  std::vector<TopologyIteration> history;
  /// Density snapshots for the evolution animation, paired with their iteration.
  std::vector<int> snapshot_iterations;
  std::vector<Vector> snapshots;
  std::vector<std::string> warnings;
  /// Displacements of the final design, one per load case.
  std::vector<Vector> displacements;
  Vector element_strain_energy;

  OptimizerMethod method = OptimizerMethod::OptimalityCriteria;
  bool stress_constrained = false;
  std::vector<StressConstraintRecord> stress;  ///< one per load case when constrained
  Scalar max_stress_ratio = 0.0;               ///< max over cases of max_relaxed_ratio
  Scalar constraint_violation = 0.0;           ///< max over all constraints at the end
  bool feasible = true;                        ///< constraint_violation <= tolerance
};

class TopologyOptimizer {
 public:
  TopologyOptimizer(const FemModel& model, const Assembler& assembler,
                    const DensityFilter& filter, const DesignDomain& domain,
                    TopologyOptimizerOptions options);

  /// Run the loop from the domain's initial design.
  TopologyOptimizationResult run();

  const TopologyOptimizerOptions& options() const { return options_; }

 private:
  Scalar penalty_for_iteration(int iteration) const;
  TopologyOptimizationResult run_oc();
  TopologyOptimizationResult run_mma();
  void finish(TopologyOptimizationResult& result, ComplianceObjective& objective,
              const Vector& x, int performed, const std::vector<Scalar>& stress_scales,
              const StressConstraint* stress) const;

  const FemModel& model_;
  const Assembler& assembler_;
  const DensityFilter& filter_;
  const DesignDomain& domain_;
  TopologyOptimizerOptions options_;
};

/// Measure of how binary a density field is:
/// \f$ M_{nd} = \frac{4}{n}\sum_e \rho_e (1-\rho_e) \f$, which is 0 for a pure
/// 0/1 design and 1 when every element sits at 0.5.
Scalar gray_level(const Vector& density);

}  // namespace sparlab
