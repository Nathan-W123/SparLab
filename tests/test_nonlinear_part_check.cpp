/// \file test_nonlinear_part_check.cpp
/// \brief The non-linear check of an exported part (NonlinearPartCheck.hpp):
///        the classification of a run's path, the small-load limit of the
///        ratios, the exact first yield and plastic collapse of a bar, the
///        bifurcation of a column against its linear buckling factor, and free
///        thermal expansion, which both analyses reproduce exactly.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/topopt/NonlinearPartCheck.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <limits>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {

constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

Selector box_selector(Scalar xmin, Scalar xmax, Scalar ymin = -kInf, Scalar ymax = kInf) {
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  s.ymin = ymin;
  s.ymax = ymax;
  return s;
}

StructuredMeshSpec strip(Index nx, Index ny, Scalar lx, Scalar ly) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.lx = lx;
  spec.ly = ly;
  return spec;
}

/// A run with converged steps at `factors` whose tangents had `pivots`
/// negative pivots, in the range of its law.
NonlinearPartCase path(const std::vector<Scalar>& factors, const std::vector<int>& pivots,
                       bool completed) {
  NonlinearPartCase c;
  NonlinearResult& r = c.nonlinear;
  r.kinematics = to_string(Kinematics::Finite);
  r.law = to_string(HyperelasticModel::SaintVenantKirchhoff);
  r.completed = completed;
  for (std::size_t i = 0; i < factors.size(); ++i) {
    NonlinearStep s;
    s.index = static_cast<int>(i) + 1;
    s.load_factor = factors[i];
    s.negative_pivots = pivots[i];
    r.steps.push_back(s);
  }
  r.load_factor = factors.empty() ? 0.0 : factors.back();
  r.min_jacobian = 0.999;
  r.max_green_strain = 1.0e-3;
  return c;
}

std::vector<StaticSolution> linear_solutions(const FemModel& model, const Assembler& assembler) {
  StaticAnalysis analysis(model, assembler);
  return analysis.solve_all();
}

/// A plane-stress bar along x held by rollers at x = 0 and one pinned
/// corner, so that a load along x or a uniform temperature leaves it in a
/// uniform state.
FemModel make_bar(Index nx, Index ny, Scalar length, Scalar height, const IsotropicMaterial& m) {
  FemModel model(make_structured_quad_mesh(strip(nx, ny, length, height)), m, 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint rollers;
  rollers.region.members.push_back(box_selector(-1.0, 0.0));
  rollers.fix_x = true;
  model.constraints().push_back(rollers);
  DisplacementConstraint pin;
  pin.region.members.push_back(box_selector(-1.0, 0.0, -1.0, 0.0));
  pin.fix_y = true;
  model.constraints().push_back(pin);
  return model;
}

}  // namespace

TEST_CASE("a run's path is classified from its converged states", "[part_check]") {
  SECTION("stable states to the design load carry it") {
    NonlinearPartCase c = path({0.25, 0.5, 0.75, 1.0}, {0, 0, 0, 0}, true);
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "carries");
    REQUIRE(c.carries_design_load);
    REQUIRE(c.stability_assessed);
    REQUIRE(c.within_law_range);
    REQUIRE(c.stable_load_factor == 1.0);
    REQUIRE(std::isnan(c.critical_lower));
  }
  SECTION("an unloading path that reached the design load carries it") {
    NonlinearPartCase c = path({0.5, 1.0, 0.5, 0.0}, {0, 0, 0, 0}, true);
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "carries");
    REQUIRE(c.stable_load_factor == 1.0);
  }
  SECTION("an arc-length path past a limit point fails, bracketed by the change of inertia") {
    // The load factor peaks between the third and fourth states, then the
    // path snaps through to lambda = 1 on an unstable branch.
    NonlinearPartCase c = path({0.3, 0.6, 0.8, 0.78, 0.9, 1.0}, {0, 0, 0, 1, 1, 0}, true);
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "fails");
    REQUIRE_FALSE(c.carries_design_load);
    REQUIRE(c.stable_load_factor == 0.8);
    REQUIRE(c.critical_lower == 0.8);
    REQUIRE(c.critical_upper == 0.78);
  }
  SECTION("load control meeting an unstable tangent fails") {
    NonlinearPartCase c = path({0.2, 0.4, 0.6, 0.65}, {0, 0, 0, 0}, false);
    c.nonlinear.critical_bound = 0.6625;
    c.nonlinear.termination = "the tangent stiffness loses positive definiteness";
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "fails");
    REQUIRE(c.critical_lower == 0.65);
    REQUIRE(c.critical_upper == 0.6625);
    REQUIRE_THAT(c.assessment, ContainsSubstring("limit or bifurcation point"));
    REQUIRE_THAT(c.assessment, !ContainsSubstring("plastic collapse"));
    c.nonlinear.plastic = true;
    assess_nonlinear_path(c);
    REQUIRE_THAT(c.assessment, ContainsSubstring("plastic collapse"));
  }
  SECTION("a step that no halving converges leaves the verdict undetermined") {
    NonlinearPartCase c = path({0.2, 0.4, 0.6, 0.65}, {0, 0, 0, 0}, false);
    c.nonlinear.unreached_load_factor = 0.65001;
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "undetermined");
    REQUIRE_FALSE(c.carries_design_load);
    REQUIRE(c.critical_lower == 0.65);
    REQUIRE(c.critical_upper == 0.65001);
    REQUIRE_THAT(c.assessment, ContainsSubstring("plastic collapse"));
  }
  SECTION("a tangent without inertia leaves stability, and the verdict, undetermined") {
    NonlinearPartCase c = path({0.5, 1.0}, {-1, -1}, true);
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "undetermined");
    REQUIRE_FALSE(c.stability_assessed);
    REQUIRE(c.stable_load_factor == 1.0);
    REQUIRE_THAT(c.assessment, ContainsSubstring("follower_pressure"));
  }
  SECTION("a path that ends below the design load does not show it") {
    NonlinearPartCase c = path({0.25, 0.5}, {0, 0}, true);
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "undetermined");
    REQUIRE(c.stable_load_factor == 0.5);
  }
  SECTION("a state compressed beyond the Saint Venant-Kirchhoff range withholds the verdict") {
    for (const Scalar j : {0.5, -0.2}) {
      NonlinearPartCase c = path({0.5, 1.0}, {0, 0}, true);
      c.nonlinear.min_jacobian = j;
      assess_nonlinear_path(c);
      REQUIRE(c.verdict == "undetermined");
      REQUIRE_FALSE(c.within_law_range);
      REQUIRE_FALSE(c.carries_design_load);
      REQUIRE_THAT(c.assessment, ContainsSubstring("withheld"));
    }
    // The neo-Hookean law has no such range (and never inverts).
    NonlinearPartCase c = path({0.5, 1.0}, {0, 0}, true);
    c.nonlinear.law = to_string(HyperelasticModel::NeoHookean);
    c.nonlinear.min_jacobian = 0.5;
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "carries");
  }
  SECTION("strains beyond the small-strain range and small-strain kinematics are qualified") {
    NonlinearPartCase c = path({0.5, 1.0}, {0, 0}, true);
    c.nonlinear.max_green_strain = 0.08;
    assess_nonlinear_path(c);
    REQUIRE(c.verdict == "carries");
    REQUIRE_THAT(c.assessment, ContainsSubstring("qualitatively"));
    NonlinearPartCase s = path({0.5, 1.0}, {0, 0}, true);
    s.nonlinear.kinematics = to_string(Kinematics::SmallStrain);
    assess_nonlinear_path(s);
    REQUIRE(s.verdict == "carries");
    REQUIRE_THAT(s.assessment, ContainsSubstring("buckling is not checked"));
  }
}

TEST_CASE("the ratios of a part tend to 1 as its load tends to 0", "[part_check]") {
  // A cantilever strip with a transverse tip resultant. Reversing the load
  // mirrors the deflection, so the end compliance P v deviates from linear at
  // second order in the load, and so does the largest displacement - or
  // faster, where the shortening that adds to it cancels part of the
  // deflection it loses. The peak stress, of the fibre whose Cauchy stress
  // the rotation raises, deviates at first order.
  std::vector<Scalar> displacement;
  std::vector<Scalar> compliance;
  std::vector<Scalar> stress;
  for (const Scalar load : {-10.0e3, -5.0e3, -2.5e3}) {
    CantileverCase cc;
    cc.tip_load = load;
    cc.poisson = 0.3;
    FemModel model = make_cantilever(cc, 40, 4);
    Assembler assembler(model);
    const std::vector<StaticSolution> linear = linear_solutions(model, assembler);
    NonlinearOptions options;
    const NonlinearPartCheck check = check_part_nonlinear(model, assembler, options, {0}, linear);
    REQUIRE(check.cases.size() == 1);
    const NonlinearPartCase& c = check.cases[0];
    INFO("tip load " << load << " N: " << c.assessment);
    REQUIRE(c.failure.empty());
    REQUIRE(c.verdict == "carries");
    REQUIRE(c.nonlinear.load_factor == 1.0);
    REQUIRE(c.linear_compliance == Approx(linear[0].compliance));
    REQUIRE_FALSE(std::isnan(c.displacement_ratio));
    REQUIRE_FALSE(std::isnan(c.compliance_ratio));
    REQUIRE_FALSE(std::isnan(c.von_mises_ratio));
    displacement.push_back(std::abs(c.displacement_ratio - 1.0));
    compliance.push_back(std::abs(c.compliance_ratio - 1.0));
    stress.push_back(std::abs(c.von_mises_ratio - 1.0));
  }
  INFO("displacement " << displacement[0] << " " << displacement[1] << " " << displacement[2]);
  INFO("compliance " << compliance[0] << " " << compliance[1] << " " << compliance[2]);
  INFO("stress " << stress[0] << " " << stress[1] << " " << stress[2]);
  REQUIRE(displacement[0] > 1.0e-4);  // a measurable non-linearity at 10 kN
  for (std::size_t k = 0; k + 1 < 3; ++k) {
    const Scalar order_u = std::log2(displacement[k] / displacement[k + 1]);
    const Scalar order_c = std::log2(compliance[k] / compliance[k + 1]);
    const Scalar order_s = std::log2(stress[k] / stress[k + 1]);
    REQUIRE(order_u >= 1.85);
    REQUIRE(order_c == Approx(2.0).margin(0.1));
    REQUIRE(order_s >= 0.75);
    REQUIRE(order_s <= 1.25);
  }
}

TEST_CASE("a bar's first yield and plastic collapse are bracketed exactly", "[part_check]") {
  // A uniform bar in tension, elastic-perfectly plastic, loaded to 1.5 times
  // its yield load: the stress is uniform, so first yield and collapse
  // coincide at lambda = 1/1.5, which the linear estimate gives exactly and
  // the small-strain run must bracket. At the collapse the tangent turns
  // singular, which load control meets as an unstable tangent.
  const Scalar yield = 250.0e6;
  IsotropicMaterial steel(200.0e9, 0.3, 7850.0, "steel");
  PlasticityParameters p;
  p.yield_stress = yield;
  steel.set_plasticity(p);
  FemModel model = make_bar(10, 2, 1.0, 0.1, steel);
  LoadCaseSpec pull;
  pull.name = "pull";
  TractionLoadSpec end;
  end.region.members.push_back(box_selector(1.0, 2.0));
  end.traction = Vector3(1.5 * yield, 0.0, 0.0);
  pull.tractions.push_back(end);
  model.load_case_specs().push_back(pull);
  model.finalize();
  Assembler assembler(model);
  const std::vector<StaticSolution> linear = linear_solutions(model, assembler);

  NonlinearOptions options;
  options.kinematics = Kinematics::SmallStrain;
  const NonlinearPartCheck check = check_part_nonlinear(model, assembler, options, {0}, linear);
  const NonlinearPartCase& c = check.cases.at(0);
  INFO(c.assessment);
  REQUIRE(c.failure.empty());
  REQUIRE(c.yield_stress == yield);
  REQUIRE(c.linear_first_yield_load_factor == Approx(1.0 / 1.5).epsilon(1.0e-10));
  REQUIRE(c.verdict == "fails");
  REQUIRE_FALSE(c.carries_design_load);
  REQUIRE(c.nonlinear.plastic);
  REQUIRE_THAT(c.assessment, ContainsSubstring("plastic collapse"));
  REQUIRE(c.critical_lower <= 1.0 / 1.5);
  REQUIRE(c.critical_upper > 1.0 / 1.5);
  REQUIRE(c.critical_upper - c.critical_lower < 1.0e-3);
  REQUIRE_THAT(c.assessment, ContainsSubstring("small-strain kinematics"));
}

TEST_CASE("a column's bifurcation agrees with its linear buckling factor", "[part_check]") {
  // A slender cantilever column under an axial dead load beyond its critical
  // load. The straight path turns unstable at the bifurcation; load control
  // brackets it, and it differs from the linear buckling factor of the same
  // mesh only by the pre-buckling strain, here P / (E A) ~ 1e-3.
  const Scalar length = 1.0;
  const Scalar h = 0.04;
  const Scalar t = 0.01;
  FemModel model(make_structured_quad_mesh(strip(40, 2, length, h)),
                 IsotropicMaterial(200.0e9, 0.3, 7850.0, "steel"), t, StressState::PlaneStress,
                 IntegrationOptions());
  DisplacementConstraint root;
  root.region.members.push_back(box_selector(-1.0, 0.0));
  root.fix_x = root.fix_y = true;
  model.constraints().push_back(root);
  LoadCaseSpec axial;
  axial.name = "axial";
  TractionLoadSpec top;
  top.region.members.push_back(box_selector(length, 2.0));
  top.traction = Vector3(-80.0e3 / (h * t), 0.0, 0.0);
  axial.tractions.push_back(top);
  model.load_case_specs().push_back(axial);
  model.finalize();
  Assembler assembler(model);

  BucklingOptions buckling;
  buckling.num_modes = 1;
  buckling.tolerance = 1.0e-12;
  const Scalar lambda_linear = analyse_buckling(model, assembler, 0, buckling).load_factors(0);
  REQUIRE(lambda_linear < 1.0);  // the load exceeds the critical one

  const std::vector<StaticSolution> linear = linear_solutions(model, assembler);
  const NonlinearPartCheck check =
      check_part_nonlinear(model, assembler, NonlinearOptions(), {0}, linear);
  const NonlinearPartCase& c = check.cases.at(0);
  INFO(c.assessment << "; linear buckling factor " << lambda_linear);
  REQUIRE(c.verdict == "fails");
  REQUIRE(c.stable_load_factor < 1.0);
  const Scalar strain = 80.0e3 / (200.0e9 * h * t);
  REQUIRE(c.critical_lower == Approx(lambda_linear).epsilon(5.0 * strain));
  REQUIRE(c.critical_upper == Approx(lambda_linear).epsilon(5.0 * strain));
  REQUIRE(c.critical_upper - c.critical_lower < 1.0e-3 * lambda_linear);
}

TEST_CASE("free thermal expansion is exact in the linear and the non-linear analysis",
          "[part_check]") {
  // The finite-strain law splits the free thermal stretch off
  // multiplicatively, so a freely heated bar takes u = alpha dT x exactly -
  // the linear solution - and stays stress-free: its stresses are round-off,
  // and no stress ratio is formed from them.
  IsotropicMaterial steel(200.0e9, 0.3, 7850.0, "steel");
  steel.set_thermal(1.2e-5, 20.0, 0.0);
  FemModel model = make_bar(10, 2, 1.0, 0.1, steel);
  LoadCaseSpec heat;
  heat.name = "heat";
  heat.temperature.source = TemperatureSpec::Source::Uniform;
  heat.temperature.uniform = 120.0;
  model.load_case_specs().push_back(heat);
  model.finalize();
  Assembler assembler(model);
  const std::vector<StaticSolution> linear = linear_solutions(model, assembler);
  const NonlinearPartCheck check =
      check_part_nonlinear(model, assembler, NonlinearOptions(), {0}, linear);
  const NonlinearPartCase& c = check.cases.at(0);
  INFO(c.assessment);
  REQUIRE(c.verdict == "carries");
  REQUIRE(c.displacement_ratio == Approx(1.0).margin(1.0e-8));
  REQUIRE(c.compliance_ratio == Approx(1.0).margin(1.0e-8));
  REQUIRE(std::isnan(c.von_mises_ratio));
  REQUIRE(c.linear_max_displacement ==
          Approx(1.2e-5 * 100.0 * std::hypot(1.0, 0.1)).epsilon(1.0e-10));
}

TEST_CASE("the part check drops monitors off the part and records failed runs",
          "[part_check]") {
  CantileverCase cc;
  FemModel model = make_cantilever(cc, 20, 2);
  Assembler assembler(model);
  const std::vector<StaticSolution> linear = linear_solutions(model, assembler);

  NonlinearOptions options;
  NonlinearMonitor tip;
  tip.name = "tip";
  tip.region.members.push_back(box_selector(cc.length, 2.0));
  tip.component = 1;
  options.monitors.push_back(tip);
  NonlinearMonitor gone = tip;
  gone.name = "gone";
  gone.region.members.clear();
  gone.region.members.push_back(box_selector(5.0, 6.0));
  options.monitors.push_back(gone);
  const NonlinearPartCheck check = check_part_nonlinear(model, assembler, options, {0}, linear);
  REQUIRE(check.dropped_monitors == std::vector<std::string>{"gone"});
  REQUIRE(check.options.monitors.size() == 1);
  REQUIRE(check.cases.at(0).nonlinear.monitor_names == std::vector<std::string>{"tip"});

  // A run the solver refuses is recorded with its reason, not thrown.
  NonlinearOptions bad;
  NonlinearMonitor z = tip;
  z.component = 2;
  bad.monitors.push_back(z);
  const NonlinearPartCheck failed = check_part_nonlinear(model, assembler, bad, {0}, linear);
  REQUIRE_THAT(failed.cases.at(0).failure, ContainsSubstring("component"));
  REQUIRE(failed.cases.at(0).verdict == "undetermined");

  // What the check cannot do is refused.
  NonlinearOptions contact;
  contact.contact.enabled = true;
  REQUIRE_THROWS_AS(check_part_nonlinear(model, assembler, contact, {0}, linear), ModelError);
  REQUIRE_THROWS_AS(check_part_nonlinear(model, assembler, options, {0}, {}), ModelError);
}
