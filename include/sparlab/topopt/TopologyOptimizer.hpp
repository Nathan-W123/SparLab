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
///     factorisation already in hand, and one aggregated buckling constraint
///     per constrained load case, \f$\mathrm{KS}(\lambda_{req}/\lambda_i) - 1
///     \le 0\f$ (BucklingConstraint.hpp), which adds a buckling eigensolve and
///     one adjoint solve per aggregated mode. The objective is scaled by the
///     initial compliance so MMA sees O(1) numbers.
///
/// **Convergence.** Two independent indicators are tracked and *either* can
/// declare convergence:
///   * the design change \f$\|x^{k+1}-x^k\|_\infty\f$ falls below
///     `change_tolerance`, or
///   * `objective_tolerance` is enabled and the relative spread
///     (max - min) / |c_k| of the last `objective_window` + 1 compliance
///     values falls below it. For a monotone history that is the change over
///     the window; unlike a two-point difference it is not fooled by an
///     oscillation whose period divides the window.
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
#include "sparlab/topopt/BucklingConstraint.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/Mma.hpp"
#include "sparlab/topopt/Projection.hpp"
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
  BucklingConstraintOptions buckling;
  StaticAnalysisOptions analysis;
  OptimizerMethod method = OptimizerMethod::OptimalityCriteria;

  int max_iterations = 200;
  /// Stop when max |dx| <= this value.
  Scalar change_tolerance = 1.0e-3;
  /// Alternative stop criterion on the relative spread (max - min) / |c_k| of
  /// the last `objective_window` + 1 compliance values; <= 0 disables it.
  Scalar objective_tolerance = 0.0;
  /// Window (in iterations) over which the compliance spread is measured.
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

  /// Heaviside projection of the filtered density with beta continuation
  /// (and, optionally, the robust formulation).
  ProjectionOptions projection;

  /// Additive-manufacturing overhang filter (applied when `overhang.filter`).
  OverhangOptions overhang;
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
  Scalar beta = 0.0;        ///< projection sharpness (0 when the projection is off)
  int linear_iterations = 0;  ///< iterative-solver iterations, all solves of the iteration
  // MMA-only records (zero for OC runs).
  Scalar max_stress_ratio = 0.0;      ///< max over cases of max_e sigma_rel / sigma_lim
  Scalar stress_constraint = 0.0;     ///< max over cases of c g_PN - 1
  Scalar min_load_factor = 0.0;       ///< min over constrained cases of lambda_1
  Scalar buckling_constraint = 0.0;   ///< max over cases of KS - 1
  int buckling_iterations = 0;        ///< subspace iterations, all cases
  // Robust formulation (zero otherwise): `compliance` is then the eroded
  // design's and `volume_fraction` the blueprint's.
  Scalar eroded_volume_fraction = 0.0;
  Scalar dilated_volume_fraction = 0.0;
  Scalar dilated_target_fraction = 0.0;  ///< rescaled target of the dilated design
  Scalar constraint_violation = 0.0;  ///< max over all constraints, > 0 = violated
};

/// The three designs of a robust run at its final iterate.
struct RobustRecord {
  Scalar eta_eroded = 0.0;
  Scalar eta_intermediate = 0.0;
  Scalar eta_dilated = 0.0;
  Scalar compliance_eroded = 0.0;        ///< the objective [J]
  Scalar compliance_intermediate = 0.0;  ///< the blueprint [J]
  Scalar compliance_dilated = 0.0;       ///< [J]
  Scalar volume_fraction_eroded = 0.0;
  Scalar volume_fraction_intermediate = 0.0;
  Scalar volume_fraction_dilated = 0.0;
  Scalar dilated_target_fraction = 0.0;
  Vector eroded_density;
  Vector dilated_density;
};

/// Buckling state of one constrained load case at the final design.
struct BucklingConstraintRecord {
  std::string load_case;
  Vector load_factors;              ///< aggregated lambda_i, ascending [-]
  Vector solid_energy_fraction;     ///< per mode
  Scalar ks = 0.0;                  ///< KS of lambda_req / lambda_i
  Scalar constraint = 0.0;          ///< ks - 1
  bool no_positive_load_factor = false;
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

  /// Heaviside projection at the end of the run (filtered_density is the
  /// density before projection; physical_density after it).
  bool projected = false;
  Scalar final_beta = 0.0;
  Scalar projection_eta = 0.0;
  Vector filtered_density;

  /// Linear solver actually used and its work over the whole run.
  std::string linear_solver;
  Index linear_iterations = 0;       ///< iterative-solver iterations, all solves
  bool has_amg_stats = false;
  AmgStats amg_stats;                ///< hierarchy of the last factorisation

  OptimizerMethod method = OptimizerMethod::OptimalityCriteria;
  bool stress_constrained = false;
  std::vector<StressConstraintRecord> stress;  ///< one per load case when constrained
  Scalar max_stress_ratio = 0.0;               ///< max over cases of max_relaxed_ratio
  bool buckling_constrained = false;
  std::vector<BucklingConstraintRecord> buckling;  ///< one per constrained load case

  /// Robust formulation: the reported density, compliance, volume and
  /// displacements above are the blueprint's (the intermediate design);
  /// the eroded and dilated ones are here.
  bool robust = false;
  RobustRecord robust_record;

  /// The overhang filter was part of the density chain; `printable_density`
  /// is its output at the final design (before projection).
  bool overhang_filtered = false;
  Vector printable_density;
  Scalar min_load_factor = 0.0;                ///< min over constrained cases of lambda_1
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
              const StressConstraint* stress, BucklingConstraint* buckling) const;

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

/// The objective-stall measure: the relative spread (max - min) / |c_k| of
/// the last `window` + 1 values of `history` (infinity while it is shorter).
/// For a monotone history it equals |c_k - c_{k-w}| / |c_k|; an oscillation
/// shows up as its full amplitude, whatever its period.
Scalar objective_spread(const std::vector<Scalar>& history, int window);

}  // namespace sparlab
