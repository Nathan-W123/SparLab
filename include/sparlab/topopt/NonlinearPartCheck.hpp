/// \file NonlinearPartCheck.hpp
/// \brief The non-linear check of a part designed by linear analysis: the part
///        analysed again with large displacement and, where its material has a
///        yield stress, J2 plasticity (NonlinearStatic.hpp), and set beside its
///        linear analysis.
///
/// The optimiser designs for linear elasticity, and the part a topology run
/// exports - the density field thresholded, its largest connected group, full
/// material - is what would be built. Whether the linear analysis describes
/// that part at its design load is a question the linear analysis cannot
/// answer. This check answers it per load case:
///
/// * **Does the part carry the load?** The non-linear solver follows the
///   equilibrium path from the unloaded part. The part is shown to carry the
///   load up to the *stable load factor*: the largest load factor of the
///   converged states before the first whose tangent is indefinite (a negative
///   pivot of its \f$LDL^T\f$ factor, which by Sylvester's law of inertia means
///   a limit or bifurcation point has been passed). It carries the design load
///   - the load case as the optimiser saw it, \f$\lambda = 1\f$ - when that
///   factor is at least 1. A run that stops short brackets the critical load
///   factor: a limit or bifurcation point (load control meeting an unstable
///   tangent or a jump to another branch), or a limit load such as a plastic
///   collapse (no step converging beyond a load factor), or, along an
///   arc-length path, the change of inertia between two converged states. A
///   tangent without inertia (a follower pressure's, factorised by LU) leaves
///   stability unassessed, and the verdict undetermined.
/// * **How far from linear is it?** At the final state, load factor
///   \f$\lambda\f$, the non-linear response over the linear one at the same
///   \f$\lambda\f$ (which is \f$\lambda\f$ times the linear response at 1):
///   the largest nodal displacement, the end compliance \f$f^T u\f$ with
///   \f$f\f$ the load vector of the case at \f$\lambda = 1\f$ (the measure
///   the optimiser minimised, thermal and body loads included), and the
///   largest element von Mises stress (the Cauchy stress, averaged over the
///   stiffness points, against the linear stress averaged over the stress
///   points: the same average for elements whose strain is linear in space).
///   Each ratio tends to 1 as the load tends to 0.
/// * **Does it yield?** For a material with a yield stress, the linear
///   estimate of the load factor of first yield - the yield stress over the
///   largest linear element von Mises stress - beside the non-linear run's
///   plastic points and largest plastic strain.
/// * **Does the path soften?** A member can buckle into a stable
///   post-buckled state, or yielding can spread, while every tangent stays
///   positive definite: the part keeps carrying load, on a far softer path.
///   The check measures the path's incremental stiffness - the load factor
///   gained per unit of the largest nodal displacement, over the loading
///   increments - against its first increment's, and brackets where it first
///   falls below half. On an imperfect part the member buckling of a linear
///   buckling analysis shows here, near its load factor, and not as an
///   unstable tangent.
///
/// **The model's range.** A verdict rests on the laws the run used, so it is
/// withheld ("undetermined") when the final state lies outside their range: an
/// integration point compressed to a volume ratio \f$J < 1/\sqrt{3}\f$ (an
/// inverted element among them) under the Saint Venant-Kirchhoff form, which
/// an elastoplastic material also takes - there its compressive force falls
/// again, so neither an equilibrium nor an instability found there is the
/// material's. Strains beyond the small strains that form and the
/// elastoplastic law assume (the solver's 5 %) are reported with the verdict,
/// which then holds only qualitatively where they occur. Small-strain
/// kinematics models no geometric non-linearity: its verdict covers yielding
/// and plastic collapse, not buckling.
///
/// The SIMP design itself is not analysed: its void elements, with a small
/// fraction of the stiffness, distort without bound under a large-displacement
/// analysis, and it is not what would be built.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"

#include <limits>
#include <string>
#include <vector>

namespace sparlab {

/// A ratio of the non-linear to the linear response deviating from 1 by more
/// than this is warned about.
constexpr Scalar kLinearDeviationWarning = 0.05;
/// The strain beyond which the solver warns that the Saint Venant-Kirchhoff
/// form, the elastoplastic law or small-strain kinematics is outside its range.
constexpr Scalar kSmallStrainRange = 0.05;
/// A path whose incremental stiffness falls below this fraction of its first
/// increment's has softened.
constexpr Scalar kSofteningRatio = 0.5;

/// The check of one load case.
struct NonlinearPartCase {
  std::string load_case;
  /// Why the non-linear analysis could not be run at all (the fields below
  /// are then unset); empty when it ran, whether or not it completed.
  std::string failure;
  NonlinearResult nonlinear;

  /// The linear analysis of the part, at lambda = 1.
  Scalar linear_compliance = 0.0;        ///< f^T u [J]
  Scalar linear_max_displacement = 0.0;  ///< largest nodal displacement [m]
  Scalar linear_max_von_mises = 0.0;     ///< largest element von Mises stress [Pa]

  /// The non-linear over the linear response at the final load factor; NaN
  /// when that factor is 0 or the linear response is round-off (a freely
  /// expanding part's stresses, a fully held part's displacements).
  Scalar displacement_ratio = std::numeric_limits<Scalar>::quiet_NaN();
  Scalar compliance_ratio = std::numeric_limits<Scalar>::quiet_NaN();
  Scalar von_mises_ratio = std::numeric_limits<Scalar>::quiet_NaN();

  /// Every converged state reported the inertia of its tangent.
  bool stability_assessed = false;
  /// The final state lies within the range of the laws (see the file
  /// comment); false withholds the verdict.
  bool within_law_range = true;
  /// The largest load factor of the converged states before the first
  /// unstable one (0 when none converged).
  Scalar stable_load_factor = 0.0;
  /// Shown to carry the design load: stability assessed, and the stable
  /// load factor at least 1.
  bool carries_design_load = false;
  /// Where the run brackets a critical point: the load factors of the states
  /// it lies between (along the path); NaN when there is none.
  Scalar critical_lower = std::numeric_limits<Scalar>::quiet_NaN();
  Scalar critical_upper = std::numeric_limits<Scalar>::quiet_NaN();
  /// The path's incremental stiffness over its first increment's: the
  /// smallest along the loading increments, and the load factors of the
  /// increment in which it first fell below kSofteningRatio. NaN where the
  /// path has no second loading increment, or never softened so far.
  Scalar min_stiffness_ratio = std::numeric_limits<Scalar>::quiet_NaN();
  Scalar softening_lower = std::numeric_limits<Scalar>::quiet_NaN();
  Scalar softening_upper = std::numeric_limits<Scalar>::quiet_NaN();
  /// "carries" (shown to carry the design load), "fails" (an unstable
  /// state met below it) or "undetermined" (neither shown: the run stopped
  /// without meeting an unstable state, ended below the design load, or
  /// could not assess stability), and the reason with its numbers.
  std::string verdict;
  std::string assessment;

  /// The initial yield stress of the part's material [Pa] and the load
  /// factor at which the linear analysis first reaches it; NaN for an
  /// elastic material.
  Scalar yield_stress = std::numeric_limits<Scalar>::quiet_NaN();
  Scalar linear_first_yield_load_factor = std::numeric_limits<Scalar>::quiet_NaN();

  std::vector<std::string> warnings;
};

struct NonlinearPartCheck {
  std::vector<NonlinearPartCase> cases;
  /// The options the runs used: the given ones less the monitors that select
  /// no node of the part (named in `dropped_monitors`).
  NonlinearOptions options;
  std::vector<std::string> dropped_monitors;
};

/// Check load cases `load_cases` of the finalised `part` (a continuum model:
/// plane or solid cells). `linear` holds the part's linear solutions of all
/// its load cases, indexed like them. A load case whose non-linear analysis
/// cannot be run is recorded with its `failure`, not thrown.
/// \throws ModelError for a shell or beam part or mismatched solutions.
NonlinearPartCheck check_part_nonlinear(const FemModel& part, const Assembler& assembler,
                                        const NonlinearOptions& options,
                                        const std::vector<std::size_t>& load_cases,
                                        const std::vector<StaticSolution>& linear);

/// Classify a finished run (used by check_part_nonlinear, exposed for the
/// tests): the stability fields, the critical bracket and the verdict.
void assess_nonlinear_path(NonlinearPartCase& c);

}  // namespace sparlab
