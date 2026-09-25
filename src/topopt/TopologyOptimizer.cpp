#include "sparlab/topopt/TopologyOptimizer.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/core/Timer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>

namespace sparlab {
namespace {

/// Beta continuation of the Heaviside projection: stage k uses
/// beta_start * factor^k (capped), and a stage ends after `beta_interval`
/// iterations or, optionally, as soon as its design change has settled.
class BetaSchedule {
 public:
  explicit BetaSchedule(const ProjectionOptions& options)
      : options_(options), stages_(options.enabled ? options.num_stages() : 1) {}

  bool enabled() const { return options_.enabled; }
  Scalar beta() const { return options_.enabled ? options_.beta_of_stage(stage_) : 0.0; }
  bool at_final() const { return !options_.enabled || stage_ + 1 >= stages_; }

  /// Called after iteration `iteration`; true when beta was raised.
  bool advance(int iteration, Scalar design_change, Scalar tolerance) {
    if (at_final()) return false;
    const int spent = iteration - start_ + 1;
    if (spent >= options_.beta_interval ||
        (options_.advance_on_convergence && design_change <= tolerance)) {
      ++stage_;
      start_ = iteration + 1;
      return true;
    }
    return false;
  }

 private:
  ProjectionOptions options_;
  int stages_ = 1;
  int stage_ = 0;
  int start_ = 1;
};

/// The overhang filter of a run, built when the options ask for it; passive
/// elements keep their density.
std::unique_ptr<OverhangFilter> make_overhang(const OverhangOptions& options, const Mesh& mesh,
                                              const DesignDomain& domain) {
  if (!options.filter) return nullptr;
  std::vector<char> passive(static_cast<std::size_t>(domain.num_elements()), 0);
  for (Index e = 0; e < domain.num_elements(); ++e) {
    passive[static_cast<std::size_t>(e)] = domain.is_free(e) ? 0 : 1;
  }
  return std::make_unique<OverhangFilter>(mesh, options, passive);
}

/// Volumes of the robust formulation at one design: the dilated design's
/// volume and gradient (the constrained quantity) and the blueprint's.
struct RobustVolumes {
  Scalar dilated = 0.0;
  Scalar intermediate = 0.0;
  Scalar eroded = 0.0;
  Vector dilated_gradient;
};

RobustVolumes robust_volumes(const ComplianceObjective& objective, const DesignDomain& domain,
                             const ObjectiveEvaluation& eval, const Vector& x,
                             const ProjectionOptions& projection) {
  RobustVolumes out;
  out.dilated = domain.volume_of(objective.physical_density_at(x, projection.dilated_eta()));
  out.intermediate = domain.volume_of(objective.physical_density_at(x, projection.eta));
  out.eroded = eval.volume;
  out.dilated_gradient =
      objective.chain_at(eval, projection.dilated_eta(), domain.element_volumes());
  return out;
}

/// ", beta = ..., CG iterations = ..." for the iteration log, when present.
std::string iteration_extras(const TopologyIteration& record) {
  std::ostringstream os;
  if (record.beta > 0.0) os << ", beta = " << record.beta;
  if (record.linear_iterations > 0) os << ", CG iterations = " << record.linear_iterations;
  return os.str();
}

}  // namespace

std::string to_string(OptimizerMethod method) {
  switch (method) {
    case OptimizerMethod::OptimalityCriteria: return "oc";
    case OptimizerMethod::MMA: return "mma";
  }
  return "unknown";
}

OptimizerMethod parse_optimizer_method(const std::string& text) {
  if (text == "oc" || text == "optimality_criteria") return OptimizerMethod::OptimalityCriteria;
  if (text == "mma") return OptimizerMethod::MMA;
  throw ConfigError("unknown optimizer method '" + text + "' (expected oc|mma)");
}

Scalar gray_level(const Vector& density) {
  if (density.size() == 0) return 0.0;
  Scalar sum = 0.0;
  for (Eigen::Index e = 0; e < density.size(); ++e) {
    sum += density(e) * (1.0 - density(e));
  }
  return 4.0 * sum / static_cast<Scalar>(density.size());
}

/// The objective-stall measure: the relative spread (max - min) / |latest| of
/// the last `window` + 1 compliance values, or infinity while the history is
/// shorter than that.
///
/// For a monotone history this is exactly |c_k - c_{k-w}| / |c_k|. Unlike that
/// two-point difference it cannot be fooled by an oscillation whose period
/// divides the window: a design cycling at the move limit (as the adaptive
/// stress-constraint scaling can make it) returns to the same compliance
/// every period, and a two-point test then declares a stall mid-cycle.
Scalar objective_spread(const std::vector<Scalar>& history, int window) {
  const std::size_t w = static_cast<std::size_t>(std::max(window, 1));
  if (history.size() <= w) return std::numeric_limits<Scalar>::infinity();
  const auto first = history.end() - static_cast<std::ptrdiff_t>(w + 1);
  const auto [lowest, highest] = std::minmax_element(first, history.end());
  return (*highest - *lowest) / std::max(std::abs(history.back()), 1.0e-30);
}

TopologyOptimizer::TopologyOptimizer(const FemModel& model, const Assembler& assembler,
                                     const DensityFilter& filter,
                                     const DesignDomain& domain,
                                     TopologyOptimizerOptions options)
    : model_(model),
      assembler_(assembler),
      filter_(filter),
      domain_(domain),
      options_(options) {
  options_.simp.validate();
  if (options_.max_iterations < 1) {
    throw ConfigError("optimizer.max_iterations must be at least 1");
  }
  if (!(options_.change_tolerance > 0.0)) {
    throw ConfigError("optimizer.change_tolerance must be positive");
  }
  if (options_.continuation_steps < 1) {
    throw ConfigError("optimizer.continuation_steps must be at least 1");
  }
  if (options_.continuation_steps > 1) {
    if (!(options_.penalty_start >= 1.0)) {
      throw ConfigError("optimizer.penalty_start must be at least 1");
    }
    if (options_.penalty_start > options_.simp.penalty) {
      throw ConfigError(
          "optimizer.penalty_start must not exceed the final simp.penalty; "
          "continuation ramps the penalty upward");
    }
    if (options_.continuation_iterations < 1) {
      throw ConfigError("optimizer.continuation_iterations must be at least 1");
    }
  }
  if (!(options_.interpretation_threshold > 0.0 &&
        options_.interpretation_threshold < 1.0)) {
    throw ConfigError("optimizer.interpretation_threshold must lie in (0, 1)");
  }
  if (options_.stress.enabled && options_.method != OptimizerMethod::MMA) {
    throw ConfigError(
        "stress constraints need optimizer.method = \"mma\"; the optimality-criteria "
        "update handles the volume constraint only");
  }
  if (options_.buckling.enabled && options_.method != OptimizerMethod::MMA) {
    throw ConfigError(
        "a buckling constraint needs optimizer.method = \"mma\"; the optimality-criteria "
        "update handles the volume constraint only");
  }
  if (options_.buckling.enabled && filter_.type() == FilterType::Sensitivity) {
    throw ConfigError(
        "a buckling constraint needs an exact design-to-density map; use filter.type = "
        "density (or none), not the heuristic sensitivity filter");
  }
  options_.buckling.validate();
  if (options_.method == OptimizerMethod::MMA) {
    options_.mma.validate();
    if (!(options_.constraint_tolerance > 0.0)) {
      throw ConfigError("optimizer.constraint_tolerance must be positive");
    }
  }
  options_.stress.validate(options_.simp);
  options_.projection.validate();
  if (options_.overhang.filter) {
    options_.overhang.validate();
    if (filter_.type() == FilterType::Sensitivity) {
      throw ConfigError(
          "the overhang filter needs the density filter (or none): the sensitivity filter "
          "has no chain rule to extend through it");
    }
  }
  if (options_.projection.enabled && filter_.type() == FilterType::Sensitivity) {
    throw ConfigError(
        "topology.projection needs the density filter (or none): the sensitivity filter "
        "changes the gradient heuristically, so the projection's chain rule cannot be "
        "applied to it");
  }
}

Scalar TopologyOptimizer::penalty_for_iteration(int iteration) const {
  const Scalar p_final = options_.simp.penalty;
  if (options_.continuation_steps <= 1) return p_final;
  const int stage = std::min(iteration / std::max(options_.continuation_iterations, 1),
                             options_.continuation_steps - 1);
  const Scalar t = static_cast<Scalar>(stage) /
                   static_cast<Scalar>(options_.continuation_steps - 1);
  return options_.penalty_start + t * (p_final - options_.penalty_start);
}

TopologyOptimizationResult TopologyOptimizer::run() {
  return options_.method == OptimizerMethod::MMA ? run_mma() : run_oc();
}

// ----------------------------------------------------------------------------
// Shared final evaluation
// ----------------------------------------------------------------------------

void TopologyOptimizer::finish(TopologyOptimizationResult& result,
                               ComplianceObjective& objective, const Vector& x,
                               int performed, const std::vector<Scalar>& stress_scales,
                               const StressConstraint* stress,
                               BucklingConstraint* buckling) const {
  result.iterations = performed;

  // Final evaluation at the converged design so the reported fields match x.
  const ObjectiveEvaluation eval = objective.evaluate(x, /*need_gradients=*/true);
  result.design = x;
  result.physical_density = eval.physical_density;
  result.stiffness_factors = eval.stiffness_factors;
  result.compliance = eval.compliance;
  result.load_case_compliance = eval.load_case_compliance;
  result.volume = eval.volume;
  result.volume_fraction = eval.volume_fraction;
  result.gray_level = gray_level(eval.physical_density);
  result.displacements = eval.displacements;
  result.element_strain_energy = eval.element_strain_energy;

  const Scalar target = domain_.volume_target();
  result.volume_constraint_violation = (result.volume - target) / target;
  result.constraint_violation = result.volume_constraint_violation;

  if (stress != nullptr) {
    result.stress_constrained = true;
    const std::vector<LoadCaseSpec>& specs = model_.load_case_specs();
    for (std::size_t l = 0; l < eval.displacements.size(); ++l) {
      const StressEvaluation se =
          stress->evaluate(objective, eval, l, stress_scales[l], /*need_gradients=*/false);
      StressConstraintRecord rec;
      rec.load_case = specs[l].name;
      rec.max_relaxed_stress_Pa = se.max_relaxed_ratio * options_.stress.limit;
      rec.max_relaxed_ratio = se.max_relaxed_ratio;
      rec.max_element = se.max_element;
      rec.p_norm_ratio = se.p_norm_ratio;
      rec.scale = se.scale;
      rec.constraint = se.constraint;
      for (Index e = 0; e < model_.mesh().num_elements(); ++e) {
        if (eval.physical_density(e) >= options_.interpretation_threshold) {
          rec.max_solid_ratio_retained = std::max(
              rec.max_solid_ratio_retained, se.solid_von_mises(e) / options_.stress.limit);
        }
      }
      result.max_stress_ratio = std::max(result.max_stress_ratio, rec.max_relaxed_ratio);
      result.constraint_violation = std::max(result.constraint_violation, rec.constraint);
      result.stress.push_back(rec);
    }
  }
  if (buckling != nullptr) {
    result.buckling_constrained = true;
    result.min_load_factor = std::numeric_limits<Scalar>::infinity();
    const std::vector<LoadCaseSpec>& specs = model_.load_case_specs();
    for (std::size_t l : buckling->load_cases()) {
      const BucklingEvaluation be =
          buckling->evaluate(objective, eval, l, /*need_gradients=*/false);
      BucklingConstraintRecord rec;
      rec.load_case = specs[l].name;
      rec.load_factors = be.load_factors;
      rec.solid_energy_fraction = be.solid_energy_fraction;
      rec.ks = be.ks;
      rec.constraint = be.constraint;
      rec.no_positive_load_factor = be.no_positive_load_factor;
      if (be.load_factors.size() > 0) {
        result.min_load_factor = std::min(result.min_load_factor, be.load_factors(0));
        if (be.solid_energy_fraction(0) < 0.5) {
          std::ostringstream os;
          os << "the lowest buckling mode of load case '" << rec.load_case << "' keeps only "
             << be.solid_energy_fraction(0)
             << " of its strain energy in solid elements (density >= "
             << options_.buckling.solid_threshold
             << "): it may be a pseudo mode of near-void material rather than a real "
                "instability. Check it on the interpreted structure";
          result.warnings.push_back(os.str());
          log::warn(os.str());
        }
      }
      result.constraint_violation = std::max(result.constraint_violation, rec.constraint);
      result.buckling.push_back(rec);
    }
  }
  // Robust formulation: the objective and the constraints above act on the
  // eroded design; what is reported and exported is the blueprint, with the
  // dilated design recorded next to it.
  if (options_.projection.enabled && options_.projection.robust) {
    result.robust = true;
    RobustRecord& rr = result.robust_record;
    const Scalar beta = objective.projection_beta();
    rr.eta_eroded = options_.projection.eroded_eta();
    rr.eta_intermediate = options_.projection.eta;
    rr.eta_dilated = options_.projection.dilated_eta();
    rr.compliance_eroded = eval.compliance;
    rr.volume_fraction_eroded = eval.volume_fraction;
    rr.eroded_density = eval.physical_density;
    objective.set_projection(beta, rr.eta_dilated);
    const ObjectiveEvaluation dilated = objective.evaluate(x, /*need_gradients=*/false);
    rr.compliance_dilated = dilated.compliance;
    rr.volume_fraction_dilated = dilated.volume_fraction;
    rr.dilated_density = dilated.physical_density;
    objective.set_projection(beta, rr.eta_intermediate);
    const ObjectiveEvaluation blueprint = objective.evaluate(x, /*need_gradients=*/false);
    rr.compliance_intermediate = blueprint.compliance;
    rr.volume_fraction_intermediate = blueprint.volume_fraction;
    rr.dilated_target_fraction =
        result.history.empty() ? 0.0 : result.history.back().dilated_target_fraction;
    objective.set_projection(beta, rr.eta_eroded);
    result.physical_density = blueprint.physical_density;
    result.stiffness_factors = blueprint.stiffness_factors;
    result.compliance = blueprint.compliance;
    result.load_case_compliance = blueprint.load_case_compliance;
    result.volume = blueprint.volume;
    result.volume_fraction = blueprint.volume_fraction;
    result.gray_level = gray_level(blueprint.physical_density);
    result.displacements = blueprint.displacements;
    result.element_strain_energy = blueprint.element_strain_energy;
    result.volume_constraint_violation = (result.volume - target) / target;
    // The volume constraint of the loop is the dilated design's against its
    // rescaled target; the blueprint follows V* only through the ratio of the
    // two volumes, as of the last rescaling, so it misses V* by however much
    // that ratio has moved since - which the summary records.
    Scalar violation = rr.dilated_target_fraction > 0.0
                           ? rr.volume_fraction_dilated / rr.dilated_target_fraction - 1.0
                           : 0.0;
    for (const StressConstraintRecord& rec : result.stress) {
      violation = std::max(violation, rec.constraint);
    }
    for (const BucklingConstraintRecord& rec : result.buckling) {
      violation = std::max(violation, rec.constraint);
    }
    result.constraint_violation = violation;
  }
  // Erosion check of a non-robust design: the same filtered field projected
  // at eta +- delta, as if the part came out uniformly thinner or thicker.
  if (options_.projection.enabled && !options_.projection.robust &&
      options_.projection.erosion_check) {
    result.erosion_checked = true;
    RobustRecord& rr = result.robust_record;
    const Scalar beta = objective.projection_beta();
    const Scalar eta = objective.projection_eta();
    rr.eta_eroded = eta + options_.projection.robust_delta;
    rr.eta_intermediate = eta;
    rr.eta_dilated = eta - options_.projection.robust_delta;
    rr.compliance_intermediate = eval.compliance;
    rr.volume_fraction_intermediate = eval.volume_fraction;
    objective.set_projection(beta, rr.eta_eroded);
    const ObjectiveEvaluation eroded = objective.evaluate(x, /*need_gradients=*/false);
    rr.compliance_eroded = eroded.compliance;
    rr.volume_fraction_eroded = eroded.volume_fraction;
    rr.eroded_density = eroded.physical_density;
    objective.set_projection(beta, rr.eta_dilated);
    const ObjectiveEvaluation dilated = objective.evaluate(x, /*need_gradients=*/false);
    rr.compliance_dilated = dilated.compliance;
    rr.volume_fraction_dilated = dilated.volume_fraction;
    rr.dilated_density = dilated.physical_density;
    objective.set_projection(beta, eta);
  }
  result.overhang_filtered = objective.overhang() != nullptr;
  if (result.overhang_filtered) result.printable_density = eval.printable_density;
  result.feasible = result.constraint_violation <= options_.constraint_tolerance;
  result.linear_solves = objective.num_solves();
  result.projected = objective.projection();
  result.final_beta = objective.projection() ? objective.projection_beta() : 0.0;
  result.projection_eta = objective.projection() ? objective.projection_eta() : 0.0;
  result.filtered_density = eval.filtered_density;
  result.linear_iterations = objective.total_linear_iterations();
  if (objective.solver() != nullptr) {
    result.linear_solver = objective.solver()->name();
    if (const AmgStats* stats = objective.solver()->amg_stats()) {
      result.has_amg_stats = true;
      result.amg_stats = *stats;
    }
  }
  if (options_.projection.enabled &&
      result.final_beta < options_.projection.beta_max) {
    std::ostringstream os;
    os << "the projection continuation reached beta = " << result.final_beta
       << " of beta_max = " << options_.projection.beta_max
       << " before the iteration cap; the design is less crisp than configured. Raise "
          "optimizer.max_iterations or lower projection.beta_interval";
    result.warnings.push_back(os.str());
    log::warn(os.str());
  }

  if (options_.history_stride > 0 &&
      (result.snapshot_iterations.empty() ||
       result.snapshot_iterations.back() != performed)) {
    result.snapshot_iterations.push_back(performed);
    result.snapshots.push_back(result.physical_density);
  }

  if (!result.converged) {
    std::ostringstream os;
    os << "the optimiser reached the iteration cap (" << options_.max_iterations
       << ") without meeting either convergence criterion: the last design change was "
       << result.final_design_change << " against a tolerance of "
       << options_.change_tolerance << ", and ";
    if (std::isfinite(result.final_objective_change)) {
      os << "the relative compliance spread over the last "
         << options_.objective_window + 1 << " iterations was "
         << result.final_objective_change << " against a tolerance of "
         << options_.objective_tolerance;
    } else {
      // Fewer compliances than the window holds, counted from the start or
      // from the last projection step (the window restarts there).
      os << "the relative compliance spread could not be measured: it needs "
         << options_.objective_window + 1
         << " iterations since the start or the last projection step, and the run "
            "ended before it had them";
    }
    if (result.method == OptimizerMethod::MMA) {
      os << " (largest constraint value " << result.constraint_violation
         << " against the feasibility tolerance " << options_.constraint_tolerance << ")";
    }
    os << ". The returned design is the last iterate, not a converged optimum";
    result.warnings.push_back(os.str());
    log::warn(os.str());
  }

  const Scalar volume_tolerance =
      result.method == OptimizerMethod::MMA ? options_.constraint_tolerance : 1.0e-6;
  if (result.robust) {
    const RobustRecord& rr = result.robust_record;
    const Scalar dilated_violation =
        rr.dilated_target_fraction > 0.0
            ? rr.volume_fraction_dilated / rr.dilated_target_fraction - 1.0
            : 0.0;
    if (dilated_violation > std::max(volume_tolerance, 1.0e-6)) {
      std::ostringstream os;
      os << "the robust volume constraint is violated: the dilated design's volume "
            "fraction "
         << rr.volume_fraction_dilated << " exceeds its rescaled target "
         << rr.dilated_target_fraction << " by a relative " << dilated_violation;
      result.warnings.push_back(os.str());
      log::warn(os.str());
    }
  } else if (result.volume_constraint_violation > volume_tolerance) {
    std::ostringstream os;
    os << "the volume constraint is violated: final volume " << result.volume
       << " m^3 exceeds the target " << target << " m^3 by a relative "
       << result.volume_constraint_violation;
    result.warnings.push_back(os.str());
    log::warn(os.str());
  }
  if (result.buckling_constrained && !result.feasible) {
    std::ostringstream os;
    os << "a constraint is violated at the final design (largest value "
       << result.constraint_violation << "); the smallest aggregated buckling load factor "
       << "is " << result.min_load_factor << " against the required "
       << options_.buckling.min_load_factor;
    result.warnings.push_back(os.str());
    log::warn(os.str());
  }
  if (result.stress_constrained && !result.feasible) {
    std::ostringstream os;
    os << "the stress constraint is violated at the final design: the largest aggregated "
          "constraint value is "
       << result.constraint_violation << " (max relaxed stress ratio "
       << result.max_stress_ratio << "); the design is infeasible, not optimal";
    result.warnings.push_back(os.str());
    log::warn(os.str());
  }
}

// ----------------------------------------------------------------------------
// Optimality criteria
// ----------------------------------------------------------------------------

TopologyOptimizationResult TopologyOptimizer::run_oc() {
  Timer total_timer;
  TopologyOptimizationResult result;
  result.method = OptimizerMethod::OptimalityCriteria;

  ComplianceObjective objective(model_, assembler_, filter_, domain_, options_.simp,
                                options_.analysis);
  const std::unique_ptr<OverhangFilter> overhang =
      make_overhang(options_.overhang, model_.mesh(), domain_);
  objective.set_overhang(overhang.get());
  BetaSchedule beta_schedule(options_.projection);
  const bool robust = options_.projection.enabled && options_.projection.robust;
  // The objective acts on the eroded design in a robust run.
  const Scalar objective_eta =
      robust ? options_.projection.eroded_eta() : options_.projection.eta;
  if (beta_schedule.enabled()) objective.set_projection(beta_schedule.beta(), objective_eta);

  // The bisection's volume is that of the physical (filtered and, when on,
  // projected) density the update would produce - of the dilated design in a
  // robust run, whose target is rescaled so the blueprint meets V*.
  const auto volume_of_design = [this, &objective, robust](const Vector& xi) {
    return domain_.volume_of(robust ? objective.physical_density_at(
                                          xi, options_.projection.dilated_eta())
                                    : objective.physical_density(xi));
  };
  Scalar dilated_target = domain_.volume_target();

  Vector x = domain_.initial_design();
  domain_.clamp(x);

  log::info("topology optimisation: ", domain_.describe());
  log::info("  filter ", to_string(filter_.type()), " radius ", filter_.radius(),
            " m (support ", filter_.average_support(), " elements), penalty ",
            options_.simp.penalty, ", emin_ratio ", options_.simp.emin_ratio,
            ", move limit ", options_.oc.move_limit);
  if (options_.continuation_steps > 1) {
    log::info("  continuation: penalty ", options_.penalty_start, " -> ",
              options_.simp.penalty, " in ", options_.continuation_steps, " stages of ",
              options_.continuation_iterations, " iterations");
  }
  if (beta_schedule.enabled()) {
    log::info("  Heaviside projection: eta ", options_.projection.eta, ", beta ",
              options_.projection.beta_start, " -> ", options_.projection.beta_max,
              " (x", options_.projection.beta_factor, " every ",
              options_.projection.beta_interval, " iterations or on convergence)");
  }

  std::vector<Scalar> compliance_history;
  ObjectiveEvaluation eval;
  int iteration = 0;
  Scalar previous_penalty = -1.0;

  for (iteration = 1; iteration <= options_.max_iterations; ++iteration) {
    Timer iter_timer;

    const Scalar penalty = penalty_for_iteration(iteration - 1);
    if (penalty != previous_penalty) {
      objective.simp().penalty = penalty;
      if (previous_penalty > 0.0) {
        log::info("continuation: SIMP penalty raised to ", penalty, " at iteration ",
                  iteration);
      }
      previous_penalty = penalty;
    }

    eval = objective.evaluate(x, /*need_gradients=*/true);

    Vector dv_dx = eval.dv_dx;
    RobustVolumes rv;
    if (robust) {
      rv = robust_volumes(objective, domain_, eval, x, options_.projection);
      if ((iteration - 1) % options_.projection.robust_volume_interval == 0 &&
          rv.intermediate > 0.0) {
        dilated_target = domain_.volume_target() * rv.dilated / rv.intermediate;
      }
      dv_dx = rv.dilated_gradient;
    }
    const OptimalityCriteriaStep step =
        optimality_criteria_update(domain_, x, eval.dc_dx, dv_dx, volume_of_design,
                                   options_.oc, robust ? dilated_target : 0.0,
                                   /*step_toward_unreachable_target=*/robust);

    TopologyIteration record;
    record.iteration = iteration;
    record.penalty = penalty;
    record.compliance = eval.compliance;
    record.volume = eval.volume;
    record.volume_fraction = eval.volume_fraction;
    record.max_change = step.max_change;
    record.lambda = step.lambda;
    record.gray_level = gray_level(eval.physical_density);
    record.bisections = step.bisections;
    record.volume_converged = step.volume_converged;
    record.beta = beta_schedule.beta();
    record.linear_iterations = eval.linear_iterations;
    if (robust) {
      record.volume = rv.intermediate;
      record.volume_fraction = rv.intermediate / domain_.domain_volume();
      record.eroded_volume_fraction = eval.volume_fraction;
      record.dilated_volume_fraction = rv.dilated / domain_.domain_volume();
      record.dilated_target_fraction = dilated_target / domain_.domain_volume();
    }
    record.seconds = iter_timer.elapsed_seconds();
    result.history.push_back(record);
    compliance_history.push_back(eval.compliance);

    if (options_.history_stride > 0 &&
        (iteration == 1 || iteration % options_.history_stride == 0)) {
      result.snapshot_iterations.push_back(iteration);
      result.snapshots.push_back(eval.physical_density);
    }

    log::info("iter ", iteration, ": c = ", eval.compliance, " J, vf = ",
              eval.volume_fraction, ", dx = ", step.max_change, ", p = ", penalty,
              ", grey = ", record.gray_level, iteration_extras(record));

    x = step.x;
    domain_.clamp(x);

    // Convergence: design change, plus an optional objective-stall criterion.
    // During continuation, only the final penalty stage may terminate the loop,
    // otherwise the run would stop before the penalty has been ramped.
    const bool at_final_penalty =
        options_.continuation_steps <= 1 ||
        iteration >= options_.continuation_iterations *
                         (options_.continuation_steps - 1);
    const bool change_ok = step.max_change <= options_.change_tolerance;
    const Scalar objective_change =
        objective_spread(compliance_history, options_.objective_window);
    const bool objective_ok = options_.objective_tolerance > 0.0 &&
                              objective_change <= options_.objective_tolerance;
    result.final_design_change = step.max_change;
    result.final_objective_change = objective_change;

    if (at_final_penalty && beta_schedule.at_final() && (change_ok || objective_ok)) {
      result.converged = true;
      result.stop_reason = change_ok ? "design_change" : "objective_stall";
      if (change_ok) {
        log::info("converged after ", iteration, " iterations: max |dx| = ",
                  step.max_change, " <= ", options_.change_tolerance);
      } else {
        log::info("converged after ", iteration, " iterations: relative compliance "
                  "spread over the last ", options_.objective_window + 1,
                  " iterations = ", objective_change, " <= ",
                  options_.objective_tolerance, " (design change ", step.max_change,
                  " is still above ", options_.change_tolerance,
                  ", so a few variables are still oscillating between bounds)");
      }
      break;
    }
    if (beta_schedule.advance(iteration, step.max_change, options_.change_tolerance)) {
      objective.set_projection(beta_schedule.beta(), objective_eta);
      compliance_history.clear();  // the objective changes with beta
      log::info("continuation: projection beta raised to ", beta_schedule.beta(),
                " after iteration ", iteration);
    }
  }

  const int performed = std::min(iteration, options_.max_iterations);
  finish(result, objective, x, performed, {}, nullptr, nullptr);

  result.total_seconds = total_timer.elapsed_seconds();
  log::info("topology optimisation finished: ", result.iterations, " iterations, ",
            result.linear_solves, " linear solves, compliance ", result.compliance,
            " J, volume fraction ", result.volume_fraction, " (target ",
            domain_.volume_fraction(), "), grey level ", result.gray_level, ", ",
            result.total_seconds, " s");
  return result;
}

// ----------------------------------------------------------------------------
// MMA
// ----------------------------------------------------------------------------

TopologyOptimizationResult TopologyOptimizer::run_mma() {
  Timer total_timer;
  TopologyOptimizationResult result;
  result.method = OptimizerMethod::MMA;

  ComplianceObjective objective(model_, assembler_, filter_, domain_, options_.simp,
                                options_.analysis);
  const std::unique_ptr<OverhangFilter> overhang =
      make_overhang(options_.overhang, model_.mesh(), domain_);
  objective.set_overhang(overhang.get());
  BetaSchedule beta_schedule(options_.projection);
  const bool robust = options_.projection.enabled && options_.projection.robust;
  const Scalar objective_eta =
      robust ? options_.projection.eroded_eta() : options_.projection.eta;
  if (beta_schedule.enabled()) objective.set_projection(beta_schedule.beta(), objective_eta);
  Scalar dilated_target = domain_.volume_target();
  std::unique_ptr<StressConstraint> stress;
  if (options_.stress.enabled) {
    stress = std::make_unique<StressConstraint>(model_, assembler_, filter_,
                                                options_.stress);
  }
  std::unique_ptr<BucklingConstraint> buckling;
  if (options_.buckling.enabled) {
    buckling = std::make_unique<BucklingConstraint>(model_, assembler_, options_.buckling);
  }

  // MMA works on the free variables only; passive ones keep their value.
  std::vector<Index> free_ids;
  for (Index e = 0; e < domain_.num_elements(); ++e) {
    if (domain_.is_free(e)) free_ids.push_back(e);
  }
  const Index nf = static_cast<Index>(free_ids.size());
  Vector lower(nf);
  Vector upper(nf);
  for (Index i = 0; i < nf; ++i) {
    lower(i) = domain_.lower_bounds()(free_ids[static_cast<std::size_t>(i)]);
    upper(i) = domain_.upper_bounds()(free_ids[static_cast<std::size_t>(i)]);
  }
  const auto pack = [&](const Vector& full) {
    Vector out(nf);
    for (Index i = 0; i < nf; ++i) out(i) = full(free_ids[static_cast<std::size_t>(i)]);
    return out;
  };

  const std::size_t num_cases = model_.load_case_specs().size();
  const Index stress_rows = stress ? static_cast<Index>(num_cases) : 0;
  const Index buckling_rows = buckling ? static_cast<Index>(buckling->load_cases().size()) : 0;
  const Index m = 1 + stress_rows + buckling_rows;
  MmaOptimizer mma(nf, m, lower, upper, options_.mma);
  std::vector<Scalar> scales(num_cases, 1.0);

  Vector x = domain_.initial_design();
  domain_.clamp(x);
  const Scalar target = domain_.volume_target();

  log::info("topology optimisation (MMA): ", domain_.describe());
  log::info("  filter ", to_string(filter_.type()), " radius ", filter_.radius(),
            " m (support ", filter_.average_support(), " elements), penalty ",
            options_.simp.penalty, ", emin_ratio ", options_.simp.emin_ratio,
            ", move limit ", options_.mma.move_limit, ", ", m, " constraint(s)");
  if (stress) {
    log::info("  stress constraint: limit ", options_.stress.limit, " Pa, P = ",
              options_.stress.p_norm, ", q = ", options_.stress.relaxation,
              ", one aggregated constraint per load case");
  }
  if (buckling) {
    log::info("  buckling constraint: lambda >= ", options_.buckling.min_load_factor,
              " for the lowest ", options_.buckling.num_modes, " modes of ",
              buckling->load_cases().size(), " load case(s), KS parameter ",
              options_.buckling.ks_parameter, ", stress interpolation rho^p without E_min");
  }
  if (options_.continuation_steps > 1) {
    log::info("  continuation: penalty ", options_.penalty_start, " -> ",
              options_.simp.penalty, " in ", options_.continuation_steps, " stages of ",
              options_.continuation_iterations, " iterations");
  }
  if (beta_schedule.enabled()) {
    log::info("  Heaviside projection: eta ", options_.projection.eta, ", beta ",
              options_.projection.beta_start, " -> ", options_.projection.beta_max,
              " (x", options_.projection.beta_factor, " every ",
              options_.projection.beta_interval, " iterations or on convergence)");
  }

  std::vector<Scalar> compliance_history;
  int iteration = 0;
  Scalar previous_penalty = -1.0;
  Scalar c_ref = 0.0;

  for (iteration = 1; iteration <= options_.max_iterations; ++iteration) {
    Timer iter_timer;

    const Scalar penalty = penalty_for_iteration(iteration - 1);
    if (penalty != previous_penalty) {
      objective.simp().penalty = penalty;
      if (previous_penalty > 0.0) {
        log::info("continuation: SIMP penalty raised to ", penalty, " at iteration ",
                  iteration);
      }
      previous_penalty = penalty;
    }

    const ObjectiveEvaluation eval = objective.evaluate(x, /*need_gradients=*/true);
    if (iteration == 1) {
      c_ref = eval.compliance;
      if (!(c_ref > 0.0)) {
        throw SolverError("the initial compliance is not positive; MMA cannot scale it");
      }
    }

    // Objective and constraints in MMA form (all O(1)).
    const Scalar f0 = eval.compliance / c_ref;
    const Vector df0 = pack(eval.dc_dx) / c_ref;
    Vector fval(m);
    Matrix dfdx(m, nf);
    RobustVolumes rv;
    if (robust) {
      // The volume constraint acts on the dilated design against a target
      // rescaled so that the blueprint meets the volume fraction.
      rv = robust_volumes(objective, domain_, eval, x, options_.projection);
      if ((iteration - 1) % options_.projection.robust_volume_interval == 0 &&
          rv.intermediate > 0.0) {
        dilated_target = target * rv.dilated / rv.intermediate;
      }
      fval(0) = rv.dilated / dilated_target - 1.0;
      dfdx.row(0) = (pack(rv.dilated_gradient) / dilated_target).transpose();
    } else {
      fval(0) = eval.volume / target - 1.0;
      dfdx.row(0) = (pack(eval.dv_dx) / target).transpose();
    }

    std::vector<StressEvaluation> stress_evals;
    Scalar max_stress_ratio = 0.0;
    Scalar max_stress_constraint = -std::numeric_limits<Scalar>::infinity();
    if (stress) {
      for (std::size_t l = 0; l < num_cases; ++l) {
        stress_evals.push_back(
            stress->evaluate(objective, eval, l, scales[l], /*need_gradients=*/true));
        const StressEvaluation& se = stress_evals.back();
        fval(1 + static_cast<Index>(l)) = se.constraint;
        dfdx.row(1 + static_cast<Index>(l)) = pack(se.dg_dx).transpose();
        max_stress_ratio = std::max(max_stress_ratio, se.max_relaxed_ratio);
        max_stress_constraint = std::max(max_stress_constraint, se.constraint);
      }
    }
    Scalar min_load_factor = std::numeric_limits<Scalar>::infinity();
    Scalar max_buckling_constraint = -std::numeric_limits<Scalar>::infinity();
    int buckling_iterations = 0;
    if (buckling) {
      for (std::size_t b = 0; b < buckling->load_cases().size(); ++b) {
        const BucklingEvaluation be = buckling->evaluate(
            objective, eval, buckling->load_cases()[b], /*need_gradients=*/true);
        const Index row = 1 + stress_rows + static_cast<Index>(b);
        fval(row) = be.constraint;
        dfdx.row(row) = pack(be.dg_dx).transpose();
        if (be.load_factors.size() > 0) {
          min_load_factor = std::min(min_load_factor, be.load_factors(0));
        }
        max_buckling_constraint = std::max(max_buckling_constraint, be.constraint);
        buckling_iterations += be.iterations;
      }
    }

    const MmaStep step = mma.update(pack(x), f0, df0, fval, dfdx);
    Vector x_new = x;
    for (Index i = 0; i < nf; ++i) x_new(free_ids[static_cast<std::size_t>(i)]) = step.x(i);
    domain_.clamp(x_new);

    TopologyIteration record;
    record.iteration = iteration;
    record.penalty = penalty;
    record.compliance = eval.compliance;
    record.volume = eval.volume;
    record.volume_fraction = eval.volume_fraction;
    record.max_change = step.max_change;
    record.lambda = step.lambda(0);
    record.gray_level = gray_level(eval.physical_density);
    record.bisections = step.subproblem_iterations;
    record.volume_converged = fval(0) <= options_.constraint_tolerance;
    record.max_stress_ratio = max_stress_ratio;
    record.stress_constraint = stress ? max_stress_constraint : 0.0;
    record.min_load_factor =
        buckling && std::isfinite(min_load_factor) ? min_load_factor : 0.0;
    record.buckling_constraint = buckling ? max_buckling_constraint : 0.0;
    record.buckling_iterations = buckling_iterations;
    record.constraint_violation = fval.maxCoeff();
    record.beta = beta_schedule.beta();
    record.linear_iterations = eval.linear_iterations + objective.adjoint_iterations();
    if (robust) {
      record.volume = rv.intermediate;
      record.volume_fraction = rv.intermediate / domain_.domain_volume();
      record.eroded_volume_fraction = eval.volume_fraction;
      record.dilated_volume_fraction = rv.dilated / domain_.domain_volume();
      record.dilated_target_fraction = dilated_target / domain_.domain_volume();
    }
    record.seconds = iter_timer.elapsed_seconds();
    result.history.push_back(record);
    compliance_history.push_back(eval.compliance);

    if (options_.history_stride > 0 &&
        (iteration == 1 || iteration % options_.history_stride == 0)) {
      result.snapshot_iterations.push_back(iteration);
      result.snapshots.push_back(eval.physical_density);
    }

    log::info("iter ", iteration, ": c = ", eval.compliance, " J, vf = ",
              eval.volume_fraction, ", dx = ", step.max_change, ", p = ", penalty,
              ", grey = ", record.gray_level, ", g_max = ", record.constraint_violation,
              (stress ? ", stress ratio = " : ""), (stress ? max_stress_ratio : 0.0),
              (buckling ? ", lambda_1 = " : ""), (buckling ? record.min_load_factor : 0.0),
              ", mma iters = ", step.subproblem_iterations, iteration_extras(record));

    const bool at_final_penalty =
        options_.continuation_steps <= 1 ||
        iteration >= options_.continuation_iterations *
                         (options_.continuation_steps - 1);
    const bool feasible = fval.maxCoeff() <= options_.constraint_tolerance;
    const bool change_ok = step.max_change <= options_.change_tolerance;
    const Scalar objective_change =
        objective_spread(compliance_history, options_.objective_window);
    const bool objective_ok = options_.objective_tolerance > 0.0 &&
                              objective_change <= options_.objective_tolerance;
    result.final_design_change = step.max_change;
    result.final_objective_change = objective_change;

    // Convergence is judged on the iterate that was just evaluated, so that
    // iterate is the design returned: the feasibility and stationarity it
    // reports would not be guaranteed for the (unevaluated) update, which under
    // the objective-stall criterion can still differ by up to the move limit.
    if (at_final_penalty && beta_schedule.at_final() && feasible &&
        (change_ok || objective_ok)) {
      result.converged = true;
      result.stop_reason = change_ok ? "design_change" : "objective_stall";
      log::info("converged after ", iteration, " iterations (", result.stop_reason,
                "): max |dx| = ", step.max_change, ", relative compliance spread = ",
                objective_change, ", largest constraint value = ", fval.maxCoeff(),
                "; the returned design is this evaluated iterate");
      break;
    }

    // Not converged: advance the design and the adaptive stress scales.
    for (std::size_t l = 0; l < stress_evals.size(); ++l) {
      scales[l] = stress->next_scale(scales[l], stress_evals[l]);
    }
    x = x_new;
    if (beta_schedule.advance(iteration, step.max_change, options_.change_tolerance)) {
      objective.set_projection(beta_schedule.beta(), objective_eta);
      compliance_history.clear();  // the objective changes with beta
      log::info("continuation: projection beta raised to ", beta_schedule.beta(),
                " after iteration ", iteration);
    }
  }

  // Converged runs return the last evaluated iterate; a run that hit the
  // iteration cap returns its last update, which finish() evaluates.
  const int performed = std::min(iteration, options_.max_iterations);
  finish(result, objective, x, performed, scales, stress.get(), buckling.get());

  result.total_seconds = total_timer.elapsed_seconds();
  log::info("topology optimisation finished (MMA): ", result.iterations, " iterations, ",
            result.linear_solves, " linear solves, compliance ", result.compliance,
            " J, volume fraction ", result.volume_fraction, " (target ",
            domain_.volume_fraction(), "), grey level ", result.gray_level,
            ", largest constraint value ", result.constraint_violation, ", ",
            result.total_seconds, " s");
  return result;
}

}  // namespace sparlab
