#include "sparlab/topopt/NonlinearPartCheck.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/StressRecovery.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

/// The largest nodal displacement magnitude of a full-length vector.
Scalar max_nodal_displacement(const Vector& u, Index nodes, int dim) {
  Scalar m = 0.0;
  for (Index n = 0; n < nodes; ++n) m = std::max(m, u.segment(n * dim, dim).norm());
  return m;
}

std::string number(Scalar v) {
  std::ostringstream os;
  os.precision(6);
  os << v;
  return os.str();
}

/// A linear response below this fraction of its scale is round-off: the
/// displacement against the largest strain times the part's size, f^T u
/// against |f| |u|.
constexpr Scalar kRoundoffResponse = 1.0e-9;
/// A linear von Mises stress below this fraction of Young's modulus times the
/// largest strain is a stress-free state's round-off.
constexpr Scalar kStressResponse = 1.0e-6;

/// numerator / (lambda * denominator), NaN where that is undefined.
Scalar ratio(Scalar numerator, Scalar lambda, Scalar denominator) {
  const Scalar d = lambda * denominator;
  if (!(std::abs(d) > 0.0) || !std::isfinite(d)) return std::numeric_limits<Scalar>::quiet_NaN();
  return numerator / d;
}

}  // namespace

void assess_nonlinear_path(NonlinearPartCase& c) {
  const NonlinearResult& r = c.nonlinear;
  // The load factor lambda = 1 is reached exactly (load control lands on the
  // ends of its path; the arc-length method's last step switches to load
  // control), so only round-off is allowed for.
  constexpr Scalar kDesignLoad = 1.0 - 1.0e-9;

  c.stability_assessed = std::all_of(r.steps.begin(), r.steps.end(),
                                     [](const NonlinearStep& s) { return s.negative_pivots >= 0; });
  Scalar stable = 0.0;    // the unloaded state is stable
  Scalar previous = 0.0;  // the load factor of the state before the current one
  bool unstable = false;
  for (const NonlinearStep& s : r.steps) {
    if (s.negative_pivots > 0) {
      c.critical_lower = previous;
      c.critical_upper = s.load_factor;
      unstable = true;
      break;
    }
    stable = std::max(stable, s.load_factor);
    previous = s.load_factor;
  }
  c.stable_load_factor = stable;

  // The incremental stiffness of the loading increments that move the part,
  // against the first one's.
  {
    Scalar initial = std::numeric_limits<Scalar>::quiet_NaN();
    Scalar lambda0 = 0.0;
    Scalar d0 = 0.0;
    for (const NonlinearStep& s : r.steps) {
      const Scalar dl = s.load_factor - lambda0;
      const Scalar dd = s.max_displacement - d0;
      const Scalar from = lambda0;
      lambda0 = s.load_factor;
      d0 = s.max_displacement;
      if (!(dl > 0.0) || !(dd > 0.0)) continue;
      const Scalar stiffness = dl / dd;
      if (std::isnan(initial)) {
        initial = stiffness;
        continue;
      }
      const Scalar ratio = stiffness / initial;
      c.min_stiffness_ratio =
          std::isnan(c.min_stiffness_ratio) ? ratio : std::min(c.min_stiffness_ratio, ratio);
      if (std::isnan(c.softening_lower) && ratio < kSofteningRatio) {
        c.softening_lower = from;
        c.softening_upper = s.load_factor;
      }
    }
  }
  const bool bound = !r.completed && !unstable && !std::isnan(r.critical_bound);
  const bool unreached = !r.completed && !unstable && !bound &&
                         !std::isnan(r.unreached_load_factor);
  if (bound) {
    c.critical_lower = r.load_factor;
    c.critical_upper = r.critical_bound;
  } else if (unreached) {
    c.critical_lower = r.load_factor;
    c.critical_upper = r.unreached_load_factor;
  }

  // The range of the laws: the Saint Venant-Kirchhoff form (which an
  // elastoplastic material takes with finite kinematics) softens in strong
  // compression, and it, the elastoplastic law and small-strain kinematics
  // assume small strains.
  const bool small = r.kinematics == to_string(Kinematics::SmallStrain);
  const bool svk_form =
      !small && (r.plastic || r.law == to_string(HyperelasticModel::SaintVenantKirchhoff));
  c.within_law_range = !(svk_form && !r.steps.empty() && r.min_jacobian < 1.0 / std::sqrt(3.0));

  std::ostringstream os;
  if (!c.within_law_range) {
    c.verdict = "undetermined";
    os << "the verdict is withheld: an integration point is compressed to a volume ratio J = "
       << number(r.min_jacobian) << (r.min_jacobian <= 0.0 ? " (an inverted element)" : "")
       << " < 1/sqrt(3), where the Saint Venant-Kirchhoff form's compressive force falls "
          "again, so the state is not the material's";
    if (r.plastic) {
      os << " (where that lies under a concentrated load, the yielding cells under it are "
            "crushing: spread the load over a larger face)";
    }
    os << "; the run ";
    if (unstable || bound || unreached) {
      os << "brackets a critical point between lambda = " << number(c.critical_lower)
         << " and " << number(c.critical_upper);
    } else {
      os << "reaches lambda = " << number(stable)
         << (c.stability_assessed ? " through states with positive definite tangents"
                                  : " (stability not assessed)");
    }
  } else if (stable >= kDesignLoad && c.stability_assessed) {
    c.carries_design_load = true;
    c.verdict = "carries";
    os << "stable equilibrium states from the unloaded part to lambda = " << number(stable)
       << ": the part carries the design load (lambda = 1)";
  } else if (stable >= kDesignLoad) {
    c.verdict = "undetermined";
    os << "equilibrium states reach lambda = " << number(stable)
       << ", but a non-symmetric tangent (a follower pressure, factorised by LU) reveals no "
          "inertia, so their stability is not assessed; with \"follower_pressure\": false the "
          "pressure keeps to the undeformed faces and the tangent is symmetric";
  } else if (unstable) {
    c.verdict = "fails";
    os << "the path passes a limit or bifurcation point below the design load: the tangent "
          "turns indefinite between the states at lambda = "
       << number(c.critical_lower) << " and " << number(c.critical_upper);
  } else if (bound) {
    c.verdict = "fails";
    os << "load control meets a critical point below the design load - a limit or "
          "bifurcation point"
       << (r.plastic ? ", or a plastic collapse," : "") << " between lambda = "
       << number(c.critical_lower) << " and " << number(c.critical_upper) << " ("
       << r.termination << ")";
  } else if (unreached) {
    c.verdict = "undetermined";
    os << "no step converged beyond lambda = " << number(c.critical_lower) << " (up to "
       << number(c.critical_upper)
       << "): a limit load below the design load - a limit point or, for a yielding "
          "material, a plastic collapse - lies there, or Newton's method failed there ("
       << r.termination << ")";
  } else if (!r.completed) {
    c.verdict = "undetermined";
    os << "the run stopped at lambda = " << number(r.load_factor) << " below the design load: "
       << r.termination;
  } else {
    c.verdict = "undetermined";
    os << "the run's load path ends at lambda = " << number(stable)
       << " at most, below the design load (lambda = 1)";
  }
  if (c.within_law_range && !std::isnan(c.softening_lower)) {
    os << "; the path softens: its incremental stiffness (load factor per largest "
          "displacement) fell below half its initial value between lambda = "
       << number(c.softening_lower) << " and " << number(c.softening_upper) << ", to "
       << number(c.min_stiffness_ratio)
       << " of it at least - a member buckling into a stable post-buckled state, or "
          "yielding spreading";
  }
  if (c.within_law_range && (small || svk_form) && r.max_green_strain > kSmallStrainRange) {
    os << "; strains reach " << number(r.max_green_strain)
       << ", beyond the small strains the law assumes, so the result holds only "
          "qualitatively where they occur";
  }
  if (small) {
    os << " (small-strain kinematics: no geometric non-linearity, so buckling is not "
          "checked)";
  }
  c.assessment = os.str();
}

NonlinearPartCheck check_part_nonlinear(const FemModel& part, const Assembler& assembler,
                                        const NonlinearOptions& options,
                                        const std::vector<std::size_t>& load_cases,
                                        const std::vector<StaticSolution>& linear) {
  const Mesh& mesh = part.mesh();
  const int dim = mesh.dim();
  if (part.dofs_per_node() != dim) {
    throw ModelError("the non-linear check of a part is formulated for plane and solid cells; "
                     "a shell or beam part has only its linear analyses (static, modal, "
                     "buckling)");
  }
  if (options.contact.enabled) {
    throw ModelError("the non-linear check of a part does not model contact: the part is "
                     "held by the supports the optimiser designed it for");
  }
  const std::size_t num_cases = part.load_case_specs().size();
  if (linear.size() != num_cases) {
    std::ostringstream os;
    os << "the non-linear check of a part needs the linear solutions of all its " << num_cases
       << " load case(s), not " << linear.size();
    throw ModelError(os.str());
  }

  NonlinearPartCheck check;
  check.options = options;
  check.options.monitors.clear();
  for (const NonlinearMonitor& m : options.monitors) {
    // A monitor region can lie wholly in the material the design removed.
    if (m.region.select_nodes(mesh).empty()) {
      check.dropped_monitors.push_back(m.name);
    } else {
      check.options.monitors.push_back(m);
    }
  }

  const Index nodes = mesh.num_nodes();
  Scalar min_yield = std::numeric_limits<Scalar>::infinity();
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const PlasticityParameters& p = part.material_of(e).plasticity();
    if (p.enabled()) min_yield = std::min(min_yield, p.yield_stress);
  }
  const bool yielding = std::isfinite(min_yield);

  NonlinearStaticAnalysis analysis(part, assembler, check.options);
  for (std::size_t l : load_cases) {
    if (l >= num_cases) throw ModelError("the non-linear check names a load case the part lacks");
    const StaticSolution& lin = linear[l];
    NonlinearPartCase c;
    c.load_case = part.load_case_specs()[l].name;
    if (lin.load_case_name != c.load_case) {
      throw ModelError("the linear solutions given to the non-linear check are not in the "
                       "order of the part's load cases");
    }
    c.linear_compliance = lin.compliance;
    c.linear_max_displacement = max_nodal_displacement(lin.displacement, nodes, dim);
    const Vector& temperature = part.load_case_data(l).temperature;
    const StressField stress = recover_stresses(part, assembler, lin.displacement, nullptr,
                                                temperature.size() > 0 ? &temperature : nullptr);
    c.linear_max_von_mises = stress.element_von_mises.maxCoeff();
    // What counts as a response rather than round-off: the largest strain of
    // the linear state, total or thermal, sets the scales. A freely expanding
    // part is stress-free, and its stresses are round-off; a fully held one
    // does not move.
    Scalar strain_scale = 0.0;
    Scalar stress_scale = 0.0;
    for (Index e = 0; e < mesh.num_elements(); ++e) {
      const IsotropicMaterial& m = part.material_of(e);
      Scalar strain = stress.element_strain.col(e).cwiseAbs().maxCoeff();
      if (temperature.size() > 0) {
        strain = std::max(strain, std::abs(m.thermal_expansion() *
                                           (stress.element_temperature(e) -
                                            m.reference_temperature())));
      }
      strain_scale = std::max(strain_scale, strain);
      stress_scale = std::max(stress_scale, m.youngs_modulus() * strain);
    }
    const Scalar length = mesh.bounding_box().extent().norm();
    const bool moves = c.linear_max_displacement > kRoundoffResponse * strain_scale * length;
    const bool stressed = c.linear_max_von_mises > kStressResponse * stress_scale;
    if (yielding) {
      c.yield_stress = min_yield;
      Scalar first = std::numeric_limits<Scalar>::infinity();
      for (Index e = 0; e < mesh.num_elements(); ++e) {
        const PlasticityParameters& p = part.material_of(e).plasticity();
        const Scalar vm = stress.element_von_mises(e);
        if (p.enabled() && vm > 0.0) first = std::min(first, p.yield_stress / vm);
      }
      if (std::isfinite(first)) c.linear_first_yield_load_factor = first;
    }

    try {
      c.nonlinear = analysis.solve(l);
    } catch (const std::exception& error) {
      c.failure = error.what();
      c.verdict = "undetermined";
      c.assessment = std::string("the non-linear analysis could not be run: ") + error.what();
      c.warnings.push_back(c.assessment);
      check.cases.push_back(std::move(c));
      continue;
    }
    const NonlinearResult& r = c.nonlinear;
    const Scalar lambda = r.load_factor;
    const Vector& f = part.load_vectors()[l];
    if (moves) {
      c.displacement_ratio = ratio(max_nodal_displacement(r.displacement, nodes, dim),
                                   std::abs(lambda), c.linear_max_displacement);
      if (std::abs(c.linear_compliance) >
          kRoundoffResponse * f.norm() * lin.displacement.norm()) {
        c.compliance_ratio = ratio(f.dot(r.displacement), lambda, c.linear_compliance);
      }
    }
    if (stressed) {
      c.von_mises_ratio =
          ratio(r.element_von_mises.size() > 0 ? r.element_von_mises.maxCoeff() : 0.0,
                std::abs(lambda), c.linear_max_von_mises);
    }
    assess_nonlinear_path(c);

    if (c.verdict != "carries") c.warnings.push_back(c.assessment);
    const auto deviation = [&](Scalar value, const std::string& what) {
      if (std::isnan(value) || !(std::abs(value - 1.0) > kLinearDeviationWarning)) return;
      std::ostringstream os;
      os << "at lambda = " << number(lambda) << " the " << what << " is " << number(value)
         << " times the linear prediction: the linear analysis the design rests on is off by "
         << number(100.0 * std::abs(value - 1.0)) << " %";
      c.warnings.push_back(os.str());
    };
    deviation(c.displacement_ratio, "largest displacement");
    deviation(c.compliance_ratio, "end compliance f^T u");
    deviation(c.von_mises_ratio, "largest von Mises stress");
    if (c.within_law_range && !std::isnan(c.softening_upper) &&
        c.softening_upper <= 1.0 + 1.0e-9) {
      std::ostringstream os;
      os << "the path softens below the design load: its incremental stiffness fell below "
            "half its initial value between lambda = "
         << number(c.softening_lower) << " and " << number(c.softening_upper)
         << ", so at the design load a member has buckled or material has yielded";
      c.warnings.push_back(os.str());
    }
    if (r.plastic && r.plastic_points > 0) {
      std::ostringstream os;
      os << r.plastic_points << " of " << r.total_points
         << " integration points have yielded (largest equivalent plastic strain "
         << number(r.max_plastic_strain)
         << "): the part, designed for linear elasticity, is not elastic at lambda = "
         << number(lambda);
      c.warnings.push_back(os.str());
    }
    check.cases.push_back(std::move(c));
  }
  return check;
}

}  // namespace sparlab
