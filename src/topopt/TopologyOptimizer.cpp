#include "sparlab/topopt/TopologyOptimizer.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/core/Timer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace sparlab {

Scalar gray_level(const Vector& density) {
  if (density.size() == 0) return 0.0;
  Scalar sum = 0.0;
  for (Eigen::Index e = 0; e < density.size(); ++e) {
    sum += density(e) * (1.0 - density(e));
  }
  return 4.0 * sum / static_cast<Scalar>(density.size());
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
  Timer total_timer;
  TopologyOptimizationResult result;

  ComplianceObjective objective(model_, assembler_, filter_, domain_, options_.simp,
                                options_.analysis);

  const auto volume_of_design = [this](const Vector& xi) {
    Vector rho = filter_.to_physical(xi);
    rho = rho.cwiseMax(0.0).cwiseMin(1.0);
    return domain_.volume_of(rho);
  };

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

    const OptimalityCriteriaStep step = optimality_criteria_update(
        domain_, x, eval.dc_dx, eval.dv_dx, volume_of_design, options_.oc);

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
              ", grey = ", record.gray_level);

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
    bool objective_ok = false;
    Scalar objective_change = std::numeric_limits<Scalar>::infinity();
    {
      const int w = std::max(options_.objective_window, 1);
      if (static_cast<int>(compliance_history.size()) > w) {
        const Scalar recent = compliance_history.back();
        const Scalar older = compliance_history[compliance_history.size() - 1 - w];
        objective_change =
            std::abs(recent - older) / std::max(std::abs(recent), 1.0e-30);
        objective_ok = options_.objective_tolerance > 0.0 &&
                       objective_change <= options_.objective_tolerance;
      }
    }
    result.final_design_change = step.max_change;
    result.final_objective_change = objective_change;

    if (at_final_penalty && (change_ok || objective_ok)) {
      result.converged = true;
      result.stop_reason = change_ok ? "design_change" : "objective_stall";
      if (change_ok) {
        log::info("converged after ", iteration, " iterations: max |dx| = ",
                  step.max_change, " <= ", options_.change_tolerance);
      } else {
        log::info("converged after ", iteration, " iterations: relative compliance "
                  "change over ", options_.objective_window, " iterations = ",
                  objective_change, " <= ", options_.objective_tolerance,
                  " (design change ", step.max_change, " is still above ",
                  options_.change_tolerance,
                  ", so a few variables are still oscillating between bounds)");
      }
      break;
    }
  }

  const int performed = std::min(iteration, options_.max_iterations);
  result.iterations = performed;

  // Final evaluation at the converged design so the reported fields match x.
  eval = objective.evaluate(x, /*need_gradients=*/true);
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
  result.linear_solves = objective.num_solves();

  const Scalar target = domain_.volume_target();
  result.volume_constraint_violation = (result.volume - target) / target;

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
       << options_.change_tolerance << ", and the relative compliance change over "
       << options_.objective_window << " iterations was "
       << result.final_objective_change << " against a tolerance of "
       << options_.objective_tolerance
       << ". The returned design is the last iterate, not a converged optimum";
    result.warnings.push_back(os.str());
    log::warn(os.str());
  }

  if (result.volume_constraint_violation > 1.0e-6) {
    std::ostringstream os;
    os << "the volume constraint is violated: final volume " << result.volume
       << " m^3 exceeds the target " << target << " m^3 by a relative "
       << result.volume_constraint_violation;
    result.warnings.push_back(os.str());
    log::warn(os.str());
  }

  result.total_seconds = total_timer.elapsed_seconds();
  log::info("topology optimisation finished: ", result.iterations, " iterations, ",
            result.linear_solves, " linear solves, compliance ", result.compliance,
            " J, volume fraction ", result.volume_fraction, " (target ",
            domain_.volume_fraction(), "), grey level ", result.gray_level, ", ",
            result.total_seconds, " s");
  return result;
}

}  // namespace sparlab
