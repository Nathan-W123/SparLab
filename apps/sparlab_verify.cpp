/// \file sparlab_verify.cpp
/// \brief Verification and validation studies.
///
/// **Verification** asks "does the code solve the equations it claims to
/// solve?" and is judged against exact answers for the discrete problem:
/// the patch test, solver agreement, and the finite-difference check of the
/// analytical sensitivities. These must pass to tight tolerances.
///
/// **Validation** asks "are the modelling assumptions appropriate?" and is
/// judged against an independent theory whose assumptions differ from the
/// model's: Euler-Bernoulli and Timoshenko beam deflection, and
/// Euler-Bernoulli natural frequencies. A finite gap there is expected and is
/// reported as a modelling difference, not an error. The studies quantify it
/// rather than asserting agreement.
///
/// Studies implemented:
///   * `sensitivity`      analytical vs central-difference topology gradients
///                        over a range of perturbation sizes;
///   * `mesh-convergence` cantilever tip deflection and compliance vs mesh
///                        size, with observed convergence order, against
///                        Euler-Bernoulli and Timoshenko theory at two
///                        Poisson ratios;
///   * `modal`            cantilever bending frequencies vs mesh size against
///                        Euler-Bernoulli theory;
///   * `solver-agreement` all five linear solvers on one small model;
///   * `patch-test`       constant-strain patch test on a distorted mesh.

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>

using namespace sparlab;

namespace {

// ---------------------------------------------------------------------------
// Shared cantilever definition used by the validation studies.
//
// L/h = 10 keeps the beam slender enough for beam theory to be a reasonable
// reference while remaining a genuinely two-dimensional model. The tip load is
// specified as a resultant spread over the tip edge nodes so the total force is
// mesh independent.
// ---------------------------------------------------------------------------
struct CantileverSpec {
  Scalar length = 1.0;        ///< L [m]
  Scalar height = 0.1;        ///< h [m]
  Scalar thickness = 0.01;    ///< t [m]
  Scalar youngs = 70.0e9;     ///< E [Pa]
  Scalar poisson = 0.3;       ///< nu [-]
  Scalar density = 2700.0;    ///< rho [kg/m^3]
  Scalar tip_load = -1000.0;  ///< P [N], negative = downward
};

struct CantileverResult {
  Index nx = 0;
  Index ny = 0;
  Index num_elements = 0;
  Index num_dofs = 0;
  Scalar element_size = 0.0;
  Scalar tip_deflection_mean = 0.0;   ///< mean u_y over the tip edge [m]
  Scalar tip_deflection_mid = 0.0;    ///< u_y at the tip mid-height node [m]
  Scalar compliance = 0.0;
  Scalar strain_energy = 0.0;
  Scalar max_von_mises = 0.0;
  Scalar reaction_error = 0.0;
  Scalar solve_seconds = 0.0;
};

/// Euler-Bernoulli tip deflection of a cantilever with an end load:
/// delta = P L^3 / (3 E I), I = t h^3 / 12.
Scalar euler_bernoulli_tip(const CantileverSpec& s) {
  const Scalar i = s.thickness * s.height * s.height * s.height / 12.0;
  return s.tip_load * s.length * s.length * s.length / (3.0 * s.youngs * i);
}

/// Timoshenko tip deflection: adds the shear term P L / (k G A) with the
/// rectangular-section shear coefficient k = 5/6.
Scalar timoshenko_tip(const CantileverSpec& s) {
  const Scalar i = s.thickness * s.height * s.height * s.height / 12.0;
  const Scalar a = s.thickness * s.height;
  const Scalar g = s.youngs / (2.0 * (1.0 + s.poisson));
  const Scalar k = 5.0 / 6.0;
  return s.tip_load * s.length * s.length * s.length / (3.0 * s.youngs * i) +
         s.tip_load * s.length / (k * g * a);
}

FemModel build_cantilever(const CantileverSpec& spec, Index nx, Index ny) {
  StructuredMeshSpec mesh_spec;
  mesh_spec.nx = nx;
  mesh_spec.ny = ny;
  mesh_spec.lx = spec.length;
  mesh_spec.ly = spec.height;

  IsotropicMaterial material(spec.youngs, spec.poisson, spec.density, "verification");
  FemModel model(make_structured_quad_mesh(mesh_spec), material, spec.thickness,
                 StressState::PlaneStress, IntegrationOptions());

  DisplacementConstraint root;
  root.region.name = "root";
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = true;
  root.fix_y = true;
  model.constraints().push_back(root);

  LoadCaseSpec load;
  load.name = "tip_load";
  PointLoadSpec tip;
  tip.region.name = "tip_edge";
  Selector tip_box;
  tip_box.kind = SelectorKind::Box;
  tip_box.xmin = spec.length;
  tip.region.members.push_back(tip_box);
  tip.force = Vector2(0.0, spec.tip_load);
  tip.distribute_total = true;
  load.point_loads.push_back(tip);
  model.load_case_specs().push_back(load);

  model.finalize();
  return model;
}

CantileverResult solve_cantilever(const CantileverSpec& spec, Index nx, Index ny) {
  CantileverResult out;
  out.nx = nx;
  out.ny = ny;
  FemModel model = build_cantilever(spec, nx, ny);
  Assembler assembler(model);

  Timer timer;
  StaticAnalysisOptions options;
  StaticAnalysis analysis(model, assembler, options);
  const std::vector<StaticSolution> solutions = analysis.solve_all();
  out.solve_seconds = timer.elapsed_seconds();

  const StaticSolution& sol = solutions.front();
  const StressField field = recover_stresses(model, assembler, sol.displacement);

  out.num_elements = model.mesh().num_elements();
  out.num_dofs = model.dofs().num_dofs();
  out.element_size = spec.length / static_cast<Scalar>(nx);
  out.compliance = sol.compliance;
  out.strain_energy = sol.strain_energy;
  out.max_von_mises = field.element_von_mises.maxCoeff();
  out.reaction_error = sol.equilibrium.relative_force_error;

  // Average u_y over the tip edge, plus the mid-height value.
  const StructuredGridInfo& info = *model.mesh().structured_info();
  Scalar sum = 0.0;
  for (Index j = 0; j <= ny; ++j) {
    const Index node = structured_node_index(info, nx, j);
    sum += sol.displacement(node * kDofsPerNode + 1);
  }
  out.tip_deflection_mean = sum / static_cast<Scalar>(ny + 1);
  const Index mid = structured_node_index(info, nx, ny / 2);
  out.tip_deflection_mid = sol.displacement(mid * kDofsPerNode + 1);
  return out;
}

/// Observed convergence order from three successive refinements of a quantity
/// converging to `reference`: p = log(e1/e2) / log(h1/h2).
Scalar observed_order(Scalar h1, Scalar e1, Scalar h2, Scalar e2) {
  if (!(e1 > 0.0) || !(e2 > 0.0) || h1 == h2) return 0.0;
  return std::log(e1 / e2) / std::log(h1 / h2);
}

// ---------------------------------------------------------------------------
// Studies
// ---------------------------------------------------------------------------

struct StudyOutcome {
  std::string name;
  bool passed = true;
  std::string metric;
  Scalar value = 0.0;
  Scalar tolerance = 0.0;
  std::string kind;  ///< "verification" or "validation"
  std::string note;
};

/// Analytical topology sensitivities vs central differences.
StudyOutcome study_sensitivity(const std::string& out_dir, json::Value& summary,
                               Scalar tolerance) {
  // A deliberately small mesh: the check needs two extra linear solves per
  // element per step size.
  CantileverSpec spec;
  spec.length = 0.6;
  spec.height = 0.3;
  spec.tip_load = -500.0;
  FemModel model = build_cantilever(spec, 12, 6);
  Assembler assembler(model);

  // Radius of 1.5 cells: large enough for the filter to couple neighbours (so
  // the chain rule is genuinely exercised), small enough to stay cheap.
  const Scalar cell = spec.length / 12.0;
  DensityFilter filter(model.mesh(), FilterType::Density, 1.5 * cell);

  std::vector<PassiveRegionSpec> passive;
  {
    // One passive solid and one passive void patch, so the exclusion logic is
    // exercised and reported.
    PassiveRegionSpec solid;
    solid.region.name = "passive_solid_patch";
    Selector s;
    s.kind = SelectorKind::Box;
    s.xmax = 0.05;
    s.ymin = 0.10;
    s.ymax = 0.20;
    solid.region.members.push_back(s);
    solid.solid = true;
    passive.push_back(solid);

    PassiveRegionSpec hole;
    hole.region.name = "passive_void_patch";
    Selector c;
    c.kind = SelectorKind::Circle;
    c.center = Vector2(0.30, 0.15);
    c.radius = 0.04;
    hole.region.members.push_back(c);
    hole.solid = false;
    passive.push_back(hole);
  }

  DesignDomain domain(model, /*volume_fraction=*/0.5, /*initial_density=*/0.5, passive);

  SimpOptions simp;
  simp.penalty = 3.0;
  simp.emin_ratio = 1.0e-9;
  StaticAnalysisOptions analysis_options;
  // The FD check differences two compliance values that agree to ~10 digits, so
  // the linear solves must be as accurate as the direct solver allows.
  analysis_options.linear.type = LinearSolverType::SimplicialLdlt;
  analysis_options.linear.residual_tolerance = 1.0e-9;

  ComplianceObjective objective(model, assembler, filter, domain, simp,
                                analysis_options);

  // A non-uniform starting design: a uniform field makes every sensitivity
  // nearly identical and hides indexing mistakes.
  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    if (!domain.is_free(e)) continue;
    const Vector2 c = model.mesh().element_centroid(e);
    x(e) = 0.35 + 0.30 * std::sin(7.0 * c.x()) * std::cos(5.0 * c.y());
  }
  domain.clamp(x);

  // Every element is tested: the mesh is small enough that two extra linear
  // solves per element per step size is still fast, and full coverage is what
  // makes an indexing mistake impossible to miss.
  std::vector<Index> targets;
  for (Index e = 0; e < domain.num_elements(); ++e) targets.push_back(e);

  const std::vector<Scalar> steps = {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4,
                                     1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8};

  CsvWriter step_csv(path_join(out_dir, "sensitivity_steps.csv"),
                     {"step[-]", "num_tested", "num_excluded", "max_absolute_error[J]",
                      "max_relative_error[-]", "rms_relative_error[-]",
                      "directional_relative_error[-]", "gradient_inf_norm[J]",
                      "passed[-]"});
  CsvWriter element_csv(
      path_join(out_dir, "sensitivity_elements.csv"),
      {"step[-]", "element", "cx[m]", "cy[m]", "analytical[J]", "finite_difference[J]",
       "absolute_error[J]", "relative_error[-]", "excluded[-]"});

  json::Value step_records = json::Value::make_array();
  Scalar best_error = std::numeric_limits<Scalar>::max();
  Scalar best_step = 0.0;
  Index excluded_count = 0;
  std::vector<std::string> exclusion_reasons;

  for (Scalar step : steps) {
    const SensitivityCheckResult check =
        verify_sensitivities(objective, domain, x, targets, step, tolerance);
    step_csv.row({step, static_cast<Scalar>(check.num_tested),
                  static_cast<Scalar>(check.num_excluded), check.max_absolute_error,
                  check.max_relative_error, check.rms_relative_error,
                  check.directional_relative_error, check.gradient_infinity_norm,
                  check.passed ? 1.0 : 0.0});
    for (const SensitivityCheckEntry& entry : check.entries) {
      const Vector2 c = model.mesh().element_centroid(entry.element);
      element_csv.row({step, static_cast<Scalar>(entry.element), c.x(), c.y(),
                       entry.analytical, entry.finite_difference, entry.absolute_error,
                       entry.relative_error, entry.excluded ? 1.0 : 0.0});
    }
    if (check.max_relative_error < best_error) {
      best_error = check.max_relative_error;
      best_step = step;
    }
    excluded_count = check.num_excluded;
    if (exclusion_reasons.empty()) {
      for (const SensitivityCheckEntry& entry : check.entries) {
        if (!entry.excluded) continue;
        if (std::find(exclusion_reasons.begin(), exclusion_reasons.end(),
                      entry.exclusion_reason) == exclusion_reasons.end()) {
          exclusion_reasons.push_back(entry.exclusion_reason);
        }
      }
    }

    json::Value rec = json::Value::make_object();
    rec.set("step", json::Value::make_number(step));
    rec.set("num_tested", json::Value::make_number(check.num_tested));
    rec.set("num_excluded", json::Value::make_number(check.num_excluded));
    rec.set("max_absolute_error_J",
            json::Value::make_number(check.max_absolute_error));
    rec.set("max_relative_error", json::Value::make_number(check.max_relative_error));
    rec.set("rms_relative_error", json::Value::make_number(check.rms_relative_error));
    rec.set("directional_relative_error",
            json::Value::make_number(check.directional_relative_error));
    rec.set("passed", json::Value::make_bool(check.passed));
    step_records.push_back(rec);
  }
  step_csv.close();
  element_csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("mesh", json::Value::make_string("12 x 6 Q4, 72 elements"));
  block.set("filter", json::Value::make_string("density, radius = 1.5 cells"));
  block.set("simp_penalty", json::Value::make_number(simp.penalty));
  block.set("design_point",
            json::Value::make_string(
                "non-uniform: x = 0.35 + 0.30 sin(7x) cos(5y) on free elements"));
  block.set("tested_elements", json::Value::make_number(targets.size()));
  block.set("excluded_elements", json::Value::make_number(excluded_count));
  block.set("exclusion_reasons", json::array_of(exclusion_reasons));
  block.set("tolerance", json::Value::make_number(tolerance));
  block.set("best_step", json::Value::make_number(best_step));
  block.set("best_max_relative_error", json::Value::make_number(best_error));
  block.set("steps", step_records);
  summary.set("sensitivity", block);

  StudyOutcome outcome;
  outcome.name = "sensitivity (analytical vs central differences)";
  outcome.kind = "verification";
  outcome.metric = "min over steps of max relative error";
  outcome.value = best_error;
  outcome.tolerance = tolerance;
  outcome.passed = best_error <= tolerance;
  std::ostringstream note;
  note << "best at h = " << best_step << "; " << excluded_count << " of "
       << targets.size() << " candidates excluded at bounds";
  outcome.note = note.str();
  return outcome;
}

/// Cantilever tip deflection and compliance vs mesh size.
StudyOutcome study_mesh_convergence(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "mesh_convergence.csv"),
                {"poisson[-]", "nx", "ny", "num_elements", "num_dofs", "h[m]",
                 "tip_mean[m]", "tip_mid[m]", "euler_bernoulli[m]", "timoshenko[m]",
                 "rel_error_eb[-]", "rel_error_timoshenko[-]", "compliance[J]",
                 "strain_energy[J]", "compliance_over_2U[-]", "max_von_mises[Pa]",
                 "reaction_rel_error[-]", "solve_seconds[s]"});

  json::Value records = json::Value::make_array();
  const std::vector<Scalar> poissons = {0.0, 0.3};
  const std::vector<std::pair<Index, Index>> meshes = {
      {10, 2}, {20, 4}, {40, 8}, {80, 16}, {160, 32}, {320, 64}};

  Scalar worst_richardson = 0.0;
  std::map<Scalar, Scalar> order_by_poisson;
  std::map<Scalar, Scalar> finest_error_by_poisson;

  for (Scalar nu : poissons) {
    CantileverSpec spec;
    spec.poisson = nu;
    const Scalar eb = euler_bernoulli_tip(spec);
    const Scalar ti = timoshenko_tip(spec);

    std::vector<CantileverResult> results;
    for (const auto& m : meshes) {
      const CantileverResult r = solve_cantilever(spec, m.first, m.second);
      results.push_back(r);
      csv.row({nu, static_cast<Scalar>(r.nx), static_cast<Scalar>(r.ny),
               static_cast<Scalar>(r.num_elements), static_cast<Scalar>(r.num_dofs),
               r.element_size, r.tip_deflection_mean, r.tip_deflection_mid, eb, ti,
               std::abs(r.tip_deflection_mean - eb) / std::abs(eb),
               std::abs(r.tip_deflection_mean - ti) / std::abs(ti), r.compliance,
               r.strain_energy,
               r.strain_energy != 0.0 ? r.compliance / (2.0 * r.strain_energy) : 0.0,
               r.max_von_mises, r.reaction_error, r.solve_seconds});

      json::Value rec = json::Value::make_object();
      rec.set("poisson", json::Value::make_number(nu));
      rec.set("nx", json::Value::make_number(r.nx));
      rec.set("ny", json::Value::make_number(r.ny));
      rec.set("num_dofs", json::Value::make_number(r.num_dofs));
      rec.set("h_m", json::Value::make_number(r.element_size));
      rec.set("tip_mean_m", json::Value::make_number(r.tip_deflection_mean));
      rec.set("compliance_J", json::Value::make_number(r.compliance));
      rec.set("relative_error_vs_timoshenko",
              json::Value::make_number(std::abs(r.tip_deflection_mean - ti) /
                                       std::abs(ti)));
      records.push_back(rec);
    }

    // Self-convergence: use the finest mesh as the reference for the observed
    // order, which measures discretisation error alone (verification) rather
    // than the beam-theory gap (validation).
    const CantileverResult& finest = results.back();
    std::vector<Scalar> errors;
    std::vector<Scalar> sizes;
    for (std::size_t i = 0; i + 1 < results.size(); ++i) {
      errors.push_back(std::abs(results[i].tip_deflection_mean -
                                finest.tip_deflection_mean) /
                       std::abs(finest.tip_deflection_mean));
      sizes.push_back(results[i].element_size);
    }
    Scalar order = 0.0;
    if (errors.size() >= 2) {
      order = observed_order(sizes[errors.size() - 2], errors[errors.size() - 2],
                             sizes[errors.size() - 1], errors[errors.size() - 1]);
    }
    order_by_poisson[nu] = order;
    finest_error_by_poisson[nu] =
        std::abs(finest.tip_deflection_mean - ti) / std::abs(ti);
    worst_richardson = std::max(worst_richardson, errors.empty() ? 0.0 : errors.back());
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (self-convergence) + "
                                             "validation (beam theory)"));
  block.set("geometry", json::Value::make_string(
                            "L = 1 m, h = 0.1 m (L/h = 10), t = 0.01 m, "
                            "E = 70 GPa, tip resultant -1000 N"));
  block.set("records", records);
  json::Value orders = json::Value::make_object();
  for (const auto& kv : order_by_poisson) {
    orders.set(app::format(kv.first, 3), json::Value::make_number(kv.second));
  }
  block.set("observed_convergence_order_tip_deflection", orders);
  json::Value finest = json::Value::make_object();
  for (const auto& kv : finest_error_by_poisson) {
    finest.set(app::format(kv.first, 3), json::Value::make_number(kv.second));
  }
  block.set("finest_mesh_relative_error_vs_timoshenko", finest);
  block.set("note",
            json::Value::make_string(
                "The plane-stress model and beam theory are different models. The "
                "Timoshenko value (bending + shear) is the closer reference; the "
                "residual gap at nu = 0.3 is the anticlastic/Poisson effect that "
                "1-D beam theory omits, and is a modelling difference rather than a "
                "code error. Self-convergence against the finest mesh isolates the "
                "discretisation error."));
  summary.set("mesh_convergence", block);

  StudyOutcome outcome;
  outcome.name = "mesh convergence (cantilever tip deflection)";
  outcome.kind = "verification + validation";
  outcome.metric = "relative error vs Timoshenko on the finest mesh (nu = 0)";
  outcome.value = finest_error_by_poisson[0.0];
  outcome.tolerance = 0.02;
  outcome.passed = outcome.value <= outcome.tolerance;
  std::ostringstream note;
  note << "observed order " << app::format(order_by_poisson[0.0], 3) << " (nu = 0), "
       << app::format(order_by_poisson[0.3], 3) << " (nu = 0.3)";
  outcome.note = note.str();
  return outcome;
}

/// Cantilever bending frequencies vs Euler-Bernoulli theory.
StudyOutcome study_modal(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "modal_convergence.csv"),
                {"nx", "ny", "num_dofs", "h[m]", "total_mass[kg]", "expected_mass[kg]",
                 "mass_error[-]", "f1_fem[Hz]", "f1_theory[Hz]", "f1_rel_error[-]",
                 "f2_fem[Hz]", "f2_theory[Hz]", "f2_rel_error[-]",
                 "f_axial_fem[Hz]", "f_axial_theory[Hz]", "f_axial_rel_error[-]",
                 "max_eigenpair_residual[-]", "seconds[s]"});

  CantileverSpec spec;
  spec.poisson = 0.0;  // removes the Poisson coupling beam theory omits
  const std::vector<std::pair<Index, Index>> meshes = {
      {20, 4}, {40, 8}, {80, 16}, {160, 32}};

  json::Value records = json::Value::make_array();
  Scalar finest_error = 0.0;
  Scalar finest_axial_error = 0.0;
  Scalar worst_mass_error = 0.0;
  for (const auto& m : meshes) {
    FemModel model = build_cantilever(spec, m.first, m.second);
    Assembler assembler(model);
    ModalAnalysisOptions options;
    options.num_modes = 8;  // enough to contain the first axial mode
    options.mass_type = MassType::Consistent;

    Timer timer;
    const ModalResult modal = solve_modal(model, assembler, options);
    const Scalar seconds = timer.elapsed_seconds();

    const Scalar f1_theory = cantilever_bending_frequency(
        spec.youngs, spec.density, spec.length, spec.height, spec.thickness, 1);
    const Scalar f2_theory = cantilever_bending_frequency(
        spec.youngs, spec.density, spec.length, spec.height, spec.thickness, 2);
    const Scalar expected_mass =
        spec.density * spec.length * spec.height * spec.thickness;
    const Scalar mass_error =
        std::abs(modal.total_mass - expected_mass) / expected_mass;
    worst_mass_error = std::max(worst_mass_error, mass_error);

    const Scalar f1 = modal.frequencies_hz(0);
    // The second *bending* mode: in a 2-D model the axial mode can appear
    // between the bending modes, so the comparison uses the next frequency
    // whose value is closest to the theoretical second bending frequency.
    Scalar f2 = 0.0;
    Scalar best = std::numeric_limits<Scalar>::max();
    for (Eigen::Index i = 1; i < modal.frequencies_hz.size(); ++i) {
      const Scalar d = std::abs(modal.frequencies_hz(i) - f2_theory);
      if (d < best) {
        best = d;
        f2 = modal.frequencies_hz(i);
      }
    }
    // The first axial mode: identified as the mode whose kinetic energy is
    // overwhelmingly in the x direction, then compared with fixed-free rod
    // theory. This is an independent analytical reference from the bending
    // series and is exact at nu = 0.
    const Scalar f_axial_theory =
        rod_axial_frequency(spec.youngs, spec.density, spec.length, 1);
    Scalar f_axial = 0.0;
    Scalar best_axial_ratio = 0.0;
    for (Eigen::Index mode = 0; mode < modal.mode_shapes.cols(); ++mode) {
      Scalar ex = 0.0;
      Scalar ey = 0.0;
      for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
        const Scalar ux = modal.mode_shapes(n * kDofsPerNode + 0, mode);
        const Scalar uy = modal.mode_shapes(n * kDofsPerNode + 1, mode);
        ex += ux * ux;
        ey += uy * uy;
      }
      const Scalar ratio = ex / std::max(ex + ey, 1.0e-300);
      if (ratio > best_axial_ratio) {
        best_axial_ratio = ratio;
        f_axial = modal.frequencies_hz(mode);
      }
    }
    const Scalar axial_error =
        std::abs(f_axial - f_axial_theory) / f_axial_theory;

    csv.row({static_cast<Scalar>(m.first), static_cast<Scalar>(m.second),
             static_cast<Scalar>(model.dofs().num_dofs()),
             spec.length / static_cast<Scalar>(m.first), modal.total_mass,
             expected_mass, mass_error, f1, f1_theory,
             std::abs(f1 - f1_theory) / f1_theory, f2, f2_theory,
             std::abs(f2 - f2_theory) / f2_theory, f_axial, f_axial_theory,
             axial_error, modal.modal_residuals.maxCoeff(), seconds});

    json::Value rec = json::Value::make_object();
    rec.set("nx", json::Value::make_number(m.first));
    rec.set("ny", json::Value::make_number(m.second));
    rec.set("f1_hz", json::Value::make_number(f1));
    rec.set("f1_theory_hz", json::Value::make_number(f1_theory));
    rec.set("f1_relative_error",
            json::Value::make_number(std::abs(f1 - f1_theory) / f1_theory));
    rec.set("f2_hz", json::Value::make_number(f2));
    rec.set("f2_theory_hz", json::Value::make_number(f2_theory));
    rec.set("f_axial_hz", json::Value::make_number(f_axial));
    rec.set("f_axial_theory_hz", json::Value::make_number(f_axial_theory));
    rec.set("f_axial_relative_error", json::Value::make_number(axial_error));
    rec.set("mass_relative_error", json::Value::make_number(mass_error));
    rec.set("max_eigenpair_residual",
            json::Value::make_number(modal.modal_residuals.maxCoeff()));
    records.push_back(rec);
    finest_error = std::abs(f1 - f1_theory) / f1_theory;
    finest_axial_error = axial_error;
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("validation (Euler-Bernoulli frequencies) "
                                             "+ verification (mass conservation)"));
  block.set("geometry",
            json::Value::make_string("L = 1 m, h = 0.1 m, t = 0.01 m, E = 70 GPa, "
                                     "rho = 2700 kg/m^3, nu = 0"));
  block.set("bending_reference",
            json::Value::make_string(
                "f_n = (beta_n L)^2 / (2 pi L^2) sqrt(E h^2 / (12 rho)), "
                "beta_1 L = 1.87510407, beta_2 L = 4.69409113"));
  block.set("axial_reference",
            json::Value::make_string(
                "fixed-free rod, f_n = (2n-1)/(4L) sqrt(E/rho); exact for a bar "
                "in uniaxial stress, which plane stress reproduces at nu = 0"));
  block.set("records", records);
  block.set("max_mass_relative_error", json::Value::make_number(worst_mass_error));
  block.set("note",
            json::Value::make_string(
                "The 2-D model includes shear deformation and rotary inertia that "
                "Euler-Bernoulli theory omits, so the computed frequencies are "
                "expected to sit slightly below the theoretical values, and the gap "
                "grows with mode number. Conservation of total mass "
                "(sum(M)/2 = rho V) is a verification check and holds to round-off."));
  summary.set("modal_validation", block);

  StudyOutcome outcome;
  outcome.name = "modal frequencies (cantilever vs Euler-Bernoulli)";
  outcome.kind = "validation";
  outcome.metric = "f1 relative error on the finest mesh";
  outcome.value = finest_error;
  outcome.tolerance = 0.02;
  outcome.passed = finest_error <= outcome.tolerance;
  std::ostringstream note;
  note << "mass conservation error <= " << app::format(worst_mass_error, 3)
       << ", axial-mode error " << app::format(finest_axial_error, 3)
       << " vs fixed-free rod theory";
  outcome.note = note.str();
  return outcome;
}

/// All linear solvers on one small model.
StudyOutcome study_solver_agreement(const std::string& out_dir, json::Value& summary) {
  CantileverSpec spec;
  FemModel model = build_cantilever(spec, 12, 4);
  Assembler assembler(model);

  const std::vector<LinearSolverType> types = {
      LinearSolverType::SimplicialLdlt, LinearSolverType::SimplicialLlt,
      LinearSolverType::SparseLu, LinearSolverType::ConjugateGradient,
      LinearSolverType::DenseLu};

  CsvWriter csv(path_join(out_dir, "solver_agreement.csv"),
                {"solver", "num_dofs", "max_displacement[m]", "compliance[J]",
                 "scaled_residual[-]", "relative_difference_vs_dense[-]",
                 "iterations", "seconds[s]"});

  Vector reference;
  std::map<std::string, Scalar> differences;
  std::vector<std::vector<std::string>> rows;
  json::Value records = json::Value::make_array();

  // Solve with the dense solver first so it can serve as the reference.
  std::vector<LinearSolverType> ordered = {LinearSolverType::DenseLu};
  for (LinearSolverType t : types) {
    if (t != LinearSolverType::DenseLu) ordered.push_back(t);
  }

  for (LinearSolverType type : ordered) {
    StaticAnalysisOptions options;
    options.linear.type = type;
    options.linear.residual_tolerance = 1.0e-7;
    options.linear.iterative_tolerance = 1.0e-14;
    Timer timer;
    StaticAnalysis analysis(model, assembler, options);
    const std::vector<StaticSolution> solutions = analysis.solve_all();
    const Scalar seconds = timer.elapsed_seconds();
    const StaticSolution& sol = solutions.front();

    if (reference.size() == 0) reference = sol.displacement;
    const Scalar diff = (sol.displacement - reference).cwiseAbs().maxCoeff() /
                        std::max(reference.cwiseAbs().maxCoeff(), 1.0e-300);
    differences[to_string(type)] = diff;

    rows.push_back({to_string(type),
                    app::format(static_cast<Scalar>(model.dofs().num_dofs()), 9),
                    app::format(sol.max_displacement_magnitude, 17),
                    app::format(sol.compliance, 17),
                    app::format(sol.scaled_residual, 6), app::format(diff, 6),
                    app::format(static_cast<Scalar>(sol.solver_iterations), 6),
                    app::format(seconds, 4)});

    json::Value rec = json::Value::make_object();
    rec.set("solver", json::Value::make_string(to_string(type)));
    rec.set("max_displacement_m",
            json::Value::make_number(sol.max_displacement_magnitude));
    rec.set("compliance_J", json::Value::make_number(sol.compliance));
    rec.set("scaled_residual", json::Value::make_number(sol.scaled_residual));
    rec.set("relative_difference_vs_dense", json::Value::make_number(diff));
    rec.set("iterations", json::Value::make_number(sol.solver_iterations));
    records.push_back(rec);
  }
  for (const auto& row : rows) csv.raw_row(row);
  csv.close();

  Scalar worst = 0.0;
  for (const auto& kv : differences) worst = std::max(worst, kv.second);

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("mesh", json::Value::make_string("12 x 4 Q4, 130 DOFs"));
  block.set("reference_solver", json::Value::make_string("dense_lu"));
  block.set("records", records);
  block.set("max_relative_difference", json::Value::make_number(worst));
  summary.set("solver_agreement", block);

  StudyOutcome outcome;
  outcome.name = "solver agreement (sparse vs dense)";
  outcome.kind = "verification";
  outcome.metric = "max relative displacement difference vs dense LU";
  outcome.value = worst;
  outcome.tolerance = 1.0e-8;
  outcome.passed = worst <= outcome.tolerance;
  return outcome;
}

/// Constant-strain patch test on a distorted mesh.
StudyOutcome study_patch_test(const std::string& out_dir, json::Value& summary) {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 4;
  spec.lx = 2.0;
  spec.ly = 1.0;

  // The exact linear displacement field: u = a + G x with a constant strain.
  const Vector2 offset(1.0e-4, -2.0e-4);
  Matrix2 gradient;
  gradient << 3.0e-4, 1.0e-4, 1.0e-4, -2.0e-4;
  const Vector3 exact_strain(gradient(0, 0), gradient(1, 1),
                             gradient(0, 1) + gradient(1, 0));

  CsvWriter csv(path_join(out_dir, "patch_test.csv"),
                {"perturbation[-]", "num_elements", "max_displacement_error[m]",
                 "relative_displacement_error[-]", "max_strain_error[-]",
                 "relative_strain_error[-]", "max_stress_error[Pa]",
                 "relative_stress_error[-]"});

  json::Value records = json::Value::make_array();
  Scalar worst_relative = 0.0;
  for (Scalar perturbation : {0.0, 0.15, 0.30}) {
    Mesh mesh = make_perturbed_quad_mesh(spec, perturbation, 20240917u);
    IsotropicMaterial material(200.0e9, 0.3, 7850.0, "patch_steel");
    FemModel model(std::move(mesh), material, 0.02, StressState::PlaneStress,
                   IntegrationOptions());

    // Prescribe the exact field on every boundary node.
    const std::vector<Mesh::BoundaryEdge> edges = model.mesh().boundary_edges();
    std::vector<char> on_boundary(
        static_cast<std::size_t>(model.mesh().num_nodes()), 0);
    for (const Mesh::BoundaryEdge& e : edges) {
      on_boundary[static_cast<std::size_t>(e.node_a)] = 1;
      on_boundary[static_cast<std::size_t>(e.node_b)] = 1;
    }
    std::vector<Index> boundary_nodes;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      if (on_boundary[static_cast<std::size_t>(n)]) boundary_nodes.push_back(n);
    }

    DisplacementConstraint bc;
    bc.region.name = "patch_boundary";
    Selector sel;
    sel.kind = SelectorKind::NodeIds;
    sel.ids = boundary_nodes;
    bc.region.members.push_back(sel);
    bc.fix_x = true;
    bc.fix_y = true;
    model.constraints().push_back(bc);

    // A zero load case: the patch test is driven purely by the prescribed
    // boundary displacements.
    // The patch test is driven entirely by the prescribed boundary field.
    LoadCaseSpec load;
    load.name = "patch";
    load.prescribed_displacement_only = true;
    model.load_case_specs().push_back(load);
    model.finalize();

    // Overwrite the prescribed values with the exact field.
    for (Index n : boundary_nodes) {
      const Vector2 x = model.mesh().node(n);
      const Vector2 u = offset + gradient * x;
      model.dofs().prescribe(n, 0, u.x());
      model.dofs().prescribe(n, 1, u.y());
    }

    Assembler assembler(model);
    StaticAnalysisOptions options;
    options.linear.residual_tolerance = 1.0e-9;
    StaticAnalysis analysis(model, assembler, options);
    const std::vector<StaticSolution> solutions = analysis.solve_all();
    const Vector& u = solutions.front().displacement;

    Scalar max_u_error = 0.0;
    Scalar u_scale = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Vector2 x = model.mesh().node(n);
      const Vector2 expected = offset + gradient * x;
      max_u_error = std::max(
          max_u_error, std::max(std::abs(u(n * kDofsPerNode + 0) - expected.x()),
                                std::abs(u(n * kDofsPerNode + 1) - expected.y())));
      u_scale = std::max(u_scale, expected.cwiseAbs().maxCoeff());
    }

    const StressField field = recover_stresses(model, assembler, u);
    const Vector3 exact_stress = material.plane_stress_matrix() * exact_strain;
    Scalar max_strain_error = 0.0;
    Scalar max_stress_error = 0.0;
    for (Index e = 0; e < model.mesh().num_elements(); ++e) {
      max_strain_error = std::max(
          max_strain_error,
          (field.element_strain.col(e) - exact_strain).cwiseAbs().maxCoeff());
      max_stress_error = std::max(
          max_stress_error,
          (field.element_stress.col(e) - exact_stress).cwiseAbs().maxCoeff());
    }

    const Scalar rel_u = max_u_error / std::max(u_scale, 1.0e-300);
    const Scalar rel_e =
        max_strain_error / std::max(exact_strain.cwiseAbs().maxCoeff(), 1.0e-300);
    const Scalar rel_s =
        max_stress_error / std::max(exact_stress.cwiseAbs().maxCoeff(), 1.0e-300);
    worst_relative = std::max({worst_relative, rel_u, rel_e, rel_s});

    csv.row({perturbation, static_cast<Scalar>(model.mesh().num_elements()),
             max_u_error, rel_u, max_strain_error, rel_e, max_stress_error, rel_s});

    json::Value rec = json::Value::make_object();
    rec.set("perturbation", json::Value::make_number(perturbation));
    rec.set("relative_displacement_error", json::Value::make_number(rel_u));
    rec.set("relative_strain_error", json::Value::make_number(rel_e));
    rec.set("relative_stress_error", json::Value::make_number(rel_s));
    records.push_back(rec);
  }
  csv.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("field", json::Value::make_string(
                         "u = (1e-4, -2e-4) + [[3e-4, 1e-4], [1e-4, -2e-4]] x, "
                         "constant strain (3e-4, -2e-4, 2e-4)"));
  block.set("records", records);
  block.set("max_relative_error", json::Value::make_number(worst_relative));
  block.set("note", json::Value::make_string(
                        "The interior displacement, the recovered strain and the "
                        "recovered stress must all reproduce the exact linear field "
                        "to round-off, including on distorted meshes."));
  summary.set("patch_test", block);

  StudyOutcome outcome;
  outcome.name = "patch test (constant strain, distorted mesh)";
  outcome.kind = "verification";
  outcome.metric = "max relative error in u, strain and stress";
  outcome.value = worst_relative;
  outcome.tolerance = 1.0e-10;
  outcome.passed = worst_relative <= outcome.tolerance;
  return outcome;
}

}  // namespace

int main(int argc, char** argv) {
  return app::run_guarded([&]() -> int {
    const std::vector<std::string> known = {"study", "output", "verbosity",
                                            "sensitivity-tolerance", "help"};
    app::CommandLine cli(argc, argv, known);
    if (cli.has("help")) {
      return app::print_usage(
          "sparlab_verify", "[--study <name>] [--output <dir>]",
          {{"--study <name>",
            "all (default) | sensitivity | mesh-convergence | modal | "
            "solver-agreement | patch-test"},
           {"--output <dir>", "output directory (default results/verification)"},
           {"--sensitivity-tolerance <t>",
            "pass threshold on the max relative gradient error (default 1e-5)"},
           {"--verbosity <lvl>", "trace|debug|info|warn|error|silent"},
           {"--help", "show this message"}});
    }
    app::apply_verbosity(cli);

    const std::string study = cli.value("study", "all");
    const std::string out_dir = cli.value("output", "results/verification");
    const Scalar sensitivity_tolerance = cli.number("sensitivity-tolerance", 1.0e-5);
    ensure_directory(out_dir);

    const bool all = study == "all";
    json::Value summary = json::Value::make_object();
    {
      json::Value header = json::Value::make_object();
      header.set("study", json::Value::make_string(study));
      header.set("units", json::Value::make_string("SI: m, N, Pa, kg, Hz, J"));
      header.set("terminology",
                 json::Value::make_string(
                     "verification = the code solves the stated discrete equations "
                     "correctly (judged against exact answers); validation = the "
                     "modelling assumptions are appropriate (judged against an "
                     "independent theory, where a finite gap is expected)"));
      summary.set("about", header);
    }

    Timer wall;
    std::vector<StudyOutcome> outcomes;
    if (all || study == "patch-test") {
      outcomes.push_back(study_patch_test(out_dir, summary));
    }
    if (all || study == "solver-agreement") {
      outcomes.push_back(study_solver_agreement(out_dir, summary));
    }
    if (all || study == "sensitivity") {
      outcomes.push_back(study_sensitivity(out_dir, summary, sensitivity_tolerance));
    }
    if (all || study == "mesh-convergence") {
      outcomes.push_back(study_mesh_convergence(out_dir, summary));
    }
    if (all || study == "modal") {
      outcomes.push_back(study_modal(out_dir, summary));
    }
    if (outcomes.empty()) {
      throw ConfigError("unknown study '" + study +
                        "'; run with --help to list the available studies");
    }

    json::Value results = json::Value::make_array();
    bool all_passed = true;
    for (const StudyOutcome& o : outcomes) {
      json::Value rec = json::Value::make_object();
      rec.set("study", json::Value::make_string(o.name));
      rec.set("kind", json::Value::make_string(o.kind));
      rec.set("metric", json::Value::make_string(o.metric));
      rec.set("value", json::Value::make_number(o.value));
      rec.set("tolerance", json::Value::make_number(o.tolerance));
      rec.set("passed", json::Value::make_bool(o.passed));
      if (!o.note.empty()) rec.set("note", json::Value::make_string(o.note));
      results.push_back(rec);
      all_passed = all_passed && o.passed;
    }
    summary.set("outcomes", results);
    summary.set("all_passed", json::Value::make_bool(all_passed));
    summary.set("total_seconds", json::Value::make_number(wall.elapsed_seconds()));

    std::ofstream out(path_join(out_dir, "summary.json"));
    if (!out) throw IoError("cannot write the verification summary");
    out << json::dump(summary, 2) << '\n';
    out.close();

    std::cout << "\nverification / validation results\n";
    std::cout << std::string(96, '-') << "\n";
    for (const StudyOutcome& o : outcomes) {
      std::cout << (o.passed ? "  PASS  " : "  FAIL  ") << o.name << "\n";
      std::cout << "          " << o.kind << " | " << o.metric << " = "
                << app::format(o.value, 4) << " (tolerance "
                << app::format(o.tolerance, 3) << ")\n";
      if (!o.note.empty()) std::cout << "          " << o.note << "\n";
    }
    std::cout << std::string(96, '-') << "\n";
    std::cout << (all_passed ? "all studies passed" : "AT LEAST ONE STUDY FAILED")
              << " in " << app::format(wall.elapsed_seconds(), 4) << " s\n";
    std::cout << "results: " << out_dir << "\n";
    return all_passed ? 0 : 1;
  });
}
