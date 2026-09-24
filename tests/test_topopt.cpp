/// \file test_topopt.cpp
/// \brief SIMP interpolation, filters, sensitivity verification (the
///        finite-difference regression test), the optimality-criteria update
///        and end-to-end optimiser behaviour.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/OptimalityCriteria.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/SimpInterpolation.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

/// The small design problem used by the sensitivity and optimiser tests.
struct TopoptFixture {
  FemModel model;
  Scalar cell = 0.0;

  explicit TopoptFixture(Index nx = 12, Index ny = 6, Scalar lx = 0.6, Scalar ly = 0.3)
      : model(make_design_model(nx, ny, lx, ly)), cell(lx / static_cast<Scalar>(nx)) {}

  static FemModel make_design_model(Index nx, Index ny, Scalar lx, Scalar ly) {
    StructuredMeshSpec spec;
    spec.nx = nx;
    spec.ny = ny;
    spec.lx = lx;
    spec.ly = ly;

    FemModel model(make_structured_quad_mesh(spec), default_material(), 0.01,
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
    load.name = "tip";
    PointLoadSpec tip;
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector3(lx, 0.0, 0.0);
    tip.region.members.push_back(nearest);
    tip.force = Vector3(0.0, -500.0, 0.0);
    load.point_loads.push_back(tip);
    model.load_case_specs().push_back(load);
    model.finalize();
    return model;
  }
};

/// A spatially varying design well inside the bounds, so central differences
/// are admissible almost everywhere.
Vector wavy_design(const FemModel& model, const DesignDomain& domain) {
  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    if (!domain.is_free(e)) continue;
    const Vector3 c = model.mesh().element_centroid(e);
    x(e) = 0.40 + 0.25 * std::sin(9.0 * c.x()) * std::cos(7.0 * c.y());
  }
  domain.clamp(x);
  return x;
}

}  // namespace

TEST_CASE("modified SIMP interpolation and its derivative are consistent",
          "[topopt][simp]") {
  SimpOptions options;
  options.penalty = 3.0;
  options.emin_ratio = 1.0e-9;

  REQUIRE(simp_stiffness_factor(0.0, options) == Approx(options.emin_ratio));
  REQUIRE(simp_stiffness_factor(1.0, options) == Approx(1.0));
  REQUIRE(simp_stiffness_factor(0.5, options) ==
          Approx(options.emin_ratio + (1.0 - options.emin_ratio) * 0.125));

  // Monotone increasing.
  for (Scalar r = 0.0; r < 1.0; r += 0.05) {
    REQUIRE(simp_stiffness_factor(r + 0.05, options) >
            simp_stiffness_factor(r, options));
  }

  // The derivative matches central differences.
  const Scalar h = 1.0e-7;
  for (Scalar r : {0.1, 0.3, 0.5, 0.75, 0.95}) {
    const Scalar fd = (simp_stiffness_factor(r + h, options) -
                       simp_stiffness_factor(r - h, options)) /
                      (2.0 * h);
    REQUIRE(simp_stiffness_derivative(r, options) == Approx(fd).epsilon(1.0e-6));
  }

  // Values outside [0, 1] are clamped rather than producing nan from pow().
  REQUIRE(simp_stiffness_factor(-0.5, options) == Approx(options.emin_ratio));
  REQUIRE(simp_stiffness_factor(1.5, options) == Approx(1.0));
}

TEST_CASE("mass interpolation laws behave as documented", "[topopt][simp]") {
  SimpOptions linear;
  linear.mass_law = MassInterpolation::Linear;
  REQUIRE(simp_mass_factor(0.0, linear) == Approx(0.0));
  REQUIRE(simp_mass_factor(0.4, linear) == Approx(0.4));
  REQUIRE(simp_mass_factor(1.0, linear) == Approx(1.0));

  SimpOptions matched;
  matched.mass_law = MassInterpolation::PenaltyMatched;
  matched.penalty = 3.0;
  matched.mass_floor = 1.0e-9;
  REQUIRE(simp_mass_factor(0.0, matched) == Approx(matched.mass_floor));
  REQUIRE(simp_mass_factor(1.0, matched) == Approx(1.0));

  // The stiffness/mass ratio stays bounded at low density under the matched
  // law, which is what suppresses spurious localised modes.
  const Scalar rho = 1.0e-3;
  const Scalar ratio_matched =
      simp_stiffness_factor(rho, matched) / simp_mass_factor(rho, matched);
  const Scalar ratio_linear =
      simp_stiffness_factor(rho, linear) / simp_mass_factor(rho, linear);
  REQUIRE(ratio_matched > ratio_linear);
  REQUIRE(ratio_matched == Approx(1.0).epsilon(0.05));

  REQUIRE(parse_mass_interpolation("linear") == MassInterpolation::Linear);
  REQUIRE(parse_mass_interpolation("penalty_matched") ==
          MassInterpolation::PenaltyMatched);
  REQUIRE_THROWS_AS(parse_mass_interpolation("nonsense"), ConfigError);
}

TEST_CASE("SIMP options are validated", "[topopt][simp][diagnostics]") {
  SimpOptions options;
  options.penalty = 0.5;
  REQUIRE_THROWS_AS(options.validate(), ConfigError);
  options.penalty = 3.0;
  options.emin_ratio = 0.0;
  REQUIRE_THROWS_AS(options.validate(), ConfigError);
  options.emin_ratio = 1.5;
  REQUIRE_THROWS_AS(options.validate(), ConfigError);
  options.emin_ratio = 1.0e-9;
  options.mass_floor = 0.0;
  REQUIRE_THROWS_AS(options.validate(), ConfigError);
  options.mass_floor = 1.0e-9;
  REQUIRE_NOTHROW(options.validate());
}

TEST_CASE("the density filter is a normalised averaging operator",
          "[topopt][filter]") {
  StructuredMeshSpec spec;
  spec.nx = 10;
  spec.ny = 8;
  spec.lx = 1.0;
  spec.ly = 0.8;
  const Mesh mesh = make_structured_quad_mesh(spec);
  const Scalar cell = 0.1;

  const DensityFilter filter(mesh, FilterType::Density, 2.0 * cell);
  REQUIRE(filter.num_elements() == mesh.num_elements());
  REQUIRE(filter.average_support() > 4.0);

  // Rows sum to one, so a constant field is preserved exactly.
  const SparseMatrix& h = filter.weights();
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    Scalar row_sum = 0.0;
    for (SparseMatrix::InnerIterator it(h, 0); it; ++it) (void)it;
    for (Eigen::Index col = 0; col < h.outerSize(); ++col) {
      for (SparseMatrix::InnerIterator it(h, col); it; ++it) {
        if (it.row() == e) row_sum += it.value();
      }
    }
    REQUIRE(row_sum == Approx(1.0).epsilon(1.0e-12));
  }
  const Vector constant = Vector::Constant(mesh.num_elements(), 0.37);
  REQUIRE((filter.to_physical(constant) - constant).cwiseAbs().maxCoeff() ==
          Approx(0.0).margin(1.0e-14));

  // Filtering never leaves [min, max] of the input (it is a convex average).
  Vector x(mesh.num_elements());
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    x(e) = (e % 2 == 0) ? 1.0 : 0.0;
  }
  const Vector rho = filter.to_physical(x);
  REQUIRE(rho.minCoeff() >= 0.0);
  REQUIRE(rho.maxCoeff() <= 1.0);
  // A checkerboard is smoothed: the filtered field has a smaller range.
  REQUIRE(rho.maxCoeff() - rho.minCoeff() < 1.0);

  // pull_back is the exact adjoint: <H x, v> == <x, H^T v>.
  Vector v(mesh.num_elements());
  for (Index e = 0; e < mesh.num_elements(); ++e) v(e) = std::sin(0.7 * e);
  REQUIRE(filter.to_physical(x).dot(v) ==
          Approx(x.dot(filter.pull_back(v))).epsilon(1.0e-12));
}

TEST_CASE("filter variants and radius handling", "[topopt][filter]") {
  StructuredMeshSpec spec;
  spec.nx = 6;
  spec.ny = 6;
  const Mesh mesh = make_structured_quad_mesh(spec);

  SECTION("no filtering is the identity") {
    const DensityFilter none(mesh, FilterType::None, 0.0);
    const Vector x = Vector::LinSpaced(mesh.num_elements(), 0.0, 1.0);
    REQUIRE(none.to_physical(x).isApprox(x));
    REQUIRE(none.pull_back(x).isApprox(x));
    REQUIRE(none.transform_gradient(x, x).isApprox(x));
  }

  SECTION("the sensitivity filter leaves the density untouched") {
    const DensityFilter sens(mesh, FilterType::Sensitivity, 0.3);
    const Vector x = Vector::Constant(mesh.num_elements(), 0.5);
    REQUIRE(sens.to_physical(x).isApprox(x));
    // It does change the gradient, and preserves its sign.
    Vector dc = Vector::Constant(mesh.num_elements(), -1.0);
    dc(0) = -5.0;
    const Vector filtered = sens.transform_gradient(x, dc);
    REQUIRE(filtered.maxCoeff() < 0.0);
    REQUIRE_FALSE(filtered.isApprox(dc));
  }

  SECTION("an invalid radius is rejected") {
    REQUIRE_THROWS_AS(DensityFilter(mesh, FilterType::Density, 0.0), ConfigError);
    REQUIRE_THROWS_AS(DensityFilter(mesh, FilterType::Density, -1.0), ConfigError);
  }

  SECTION("length mismatches are rejected") {
    const DensityFilter filter(mesh, FilterType::Density, 0.3);
    REQUIRE_THROWS_AS(filter.to_physical(Vector::Ones(3)), ConfigError);
    REQUIRE_THROWS_AS(filter.pull_back(Vector::Ones(3)), ConfigError);
  }

  REQUIRE(parse_filter_type("density") == FilterType::Density);
  REQUIRE(parse_filter_type("none") == FilterType::None);
  REQUIRE(parse_filter_type("sensitivity") == FilterType::Sensitivity);
  REQUIRE_THROWS_AS(parse_filter_type("blur"), ConfigError);
}

TEST_CASE("design domain bounds, passive regions and feasibility checks",
          "[topopt][domain]") {
  TopoptFixture fixture;

  SECTION("a plain domain has free variables everywhere") {
    const DesignDomain domain(fixture.model, 0.4, -1.0, {});
    REQUIRE(domain.num_elements() == fixture.model.mesh().num_elements());
    REQUIRE(domain.num_free_variables() == domain.num_elements());
    REQUIRE(domain.lower_bounds().maxCoeff() == Approx(0.0));
    REQUIRE(domain.upper_bounds().minCoeff() == Approx(1.0));
    // Uniform start at the volume fraction.
    REQUIRE(domain.initial_design().minCoeff() == Approx(0.4));
    REQUIRE(domain.volume_target() ==
            Approx(0.4 * fixture.model.domain_volume()));
    REQUIRE(domain.fraction_of(domain.initial_design()) == Approx(0.4));
  }

  SECTION("passive regions pin their variables") {
    PassiveRegionSpec solid;
    solid.region.name = "solid_strip";
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.05;
    solid.region.members.push_back(box);
    solid.solid = true;

    PassiveRegionSpec hole;
    hole.region.name = "hole";
    Selector circle;
    circle.kind = SelectorKind::Circle;
    circle.center = Vector3(0.3, 0.15, 0.0);
    circle.radius = 0.05;
    hole.region.members.push_back(circle);
    hole.solid = false;

    const DesignDomain domain(fixture.model, 0.5, -1.0, {solid, hole});
    REQUIRE(domain.num_passive_solid() > 0);
    REQUIRE(domain.num_passive_void() > 0);
    REQUIRE(domain.num_free_variables() ==
            domain.num_elements() - domain.num_passive_solid() -
                domain.num_passive_void());
    REQUIRE(domain.passive_solid_volume() > 0.0);

    for (Index e = 0; e < domain.num_elements(); ++e) {
      if (domain.tags()[static_cast<std::size_t>(e)] == PassiveTag::Solid) {
        REQUIRE(domain.lower_bounds()(e) == Approx(1.0));
        REQUIRE(domain.upper_bounds()(e) == Approx(1.0));
        REQUIRE(domain.initial_design()(e) == Approx(1.0));
        REQUIRE_FALSE(domain.is_free(e));
      } else if (domain.tags()[static_cast<std::size_t>(e)] == PassiveTag::Void) {
        REQUIRE(domain.upper_bounds()(e) == Approx(0.0));
        REQUIRE_FALSE(domain.is_free(e));
      }
    }

    // Clamping respects the pinned values.
    Vector x = Vector::Constant(domain.num_elements(), 0.5);
    domain.clamp(x);
    for (Index e = 0; e < domain.num_elements(); ++e) {
      REQUIRE(x(e) >= domain.lower_bounds()(e) - 1.0e-15);
      REQUIRE(x(e) <= domain.upper_bounds()(e) + 1.0e-15);
    }
  }

  SECTION("infeasible specifications are rejected with an explanation") {
    // Volume fraction outside (0, 1].
    REQUIRE_THROWS_AS(DesignDomain(fixture.model, 0.0, -1.0, {}), ConfigError);
    REQUIRE_THROWS_AS(DesignDomain(fixture.model, 1.5, -1.0, {}), ConfigError);

    // A passive solid region larger than the volume budget.
    PassiveRegionSpec solid;
    solid.region.name = "too_big";
    Selector all;
    all.kind = SelectorKind::Box;
    all.xmax = 0.4;
    solid.region.members.push_back(all);
    solid.solid = true;
    REQUIRE_THROWS_AS(DesignDomain(fixture.model, 0.1, -1.0, {solid}), ConfigError);

    // A passive void region so large the target cannot be reached.
    PassiveRegionSpec hole;
    hole.region.name = "huge_hole";
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmin = 0.1;
    hole.region.members.push_back(box);
    hole.solid = false;
    REQUIRE_THROWS_AS(DesignDomain(fixture.model, 0.9, -1.0, {hole}), ConfigError);

    // An empty region is an authoring mistake.
    PassiveRegionSpec empty;
    empty.region.name = "nowhere";
    Selector far;
    far.kind = SelectorKind::Circle;
    far.center = Vector3(100.0, 100.0, 0.0);
    far.radius = 0.001;
    empty.region.members.push_back(far);
    REQUIRE_THROWS_AS(DesignDomain(fixture.model, 0.5, -1.0, {empty}), ConfigError);

    // Overlapping solid and void claims on the same element.
    PassiveRegionSpec a;
    a.region.name = "a";
    Selector ba;
    ba.kind = SelectorKind::Box;
    ba.xmax = 0.1;
    a.region.members.push_back(ba);
    a.solid = true;
    PassiveRegionSpec b = a;
    b.region.name = "b";
    b.solid = false;
    REQUIRE_THROWS_AS(DesignDomain(fixture.model, 0.5, -1.0, {a, b}), ConfigError);
  }
}

TEST_CASE("analytical compliance sensitivities match central differences",
          "[topopt][sensitivity][verification][regression]") {
  // This is the sensitivity regression test required of the project: the
  // analytical gradient of the *filtered* SIMP compliance is compared with a
  // central difference of the same objective on a small mesh, including passive
  // regions so the exclusion logic is exercised.
  TopoptFixture fixture;
  const DensityFilter filter(fixture.model.mesh(), FilterType::Density, 1.5 * fixture.cell);

  PassiveRegionSpec solid;
  solid.region.name = "root_pad";
  Selector sbox;
  sbox.kind = SelectorKind::Box;
  sbox.xmax = 0.05;
  sbox.ymin = 0.10;
  sbox.ymax = 0.20;
  solid.region.members.push_back(sbox);
  solid.solid = true;

  PassiveRegionSpec hole;
  hole.region.name = "bolt_hole";
  Selector circle;
  circle.kind = SelectorKind::Circle;
  circle.center = Vector3(0.30, 0.15, 0.0);
  circle.radius = 0.04;
  hole.region.members.push_back(circle);
  hole.solid = false;

  DesignDomain domain(fixture.model, 0.5, 0.5, {solid, hole});
  Assembler assembler(fixture.model);

  SimpOptions simp;
  simp.penalty = 3.0;
  StaticAnalysisOptions options;
  options.linear.residual_tolerance = 1.0e-9;
  ComplianceObjective objective(fixture.model, assembler, filter, domain, simp,
                                options);

  const Vector x = wavy_design(fixture.model, domain);
  std::vector<Index> targets;
  for (Index e = 0; e < domain.num_elements(); e += 5) targets.push_back(e);

  const Scalar tolerance = 1.0e-5;
  const SensitivityCheckResult check =
      verify_sensitivities(objective, domain, x, targets, 1.0e-4, tolerance);

  INFO("max relative error " << check.max_relative_error << ", RMS "
                             << check.rms_relative_error << ", directional "
                             << check.directional_relative_error << ", tested "
                             << check.num_tested << ", excluded "
                             << check.num_excluded);
  REQUIRE(check.num_tested > 5);
  REQUIRE(check.num_excluded > 0);  // the passive elements must be excluded
  REQUIRE(check.passed);
  REQUIRE(check.max_relative_error < tolerance);
  REQUIRE(check.rms_relative_error < 0.1 * tolerance);
  REQUIRE(check.directional_relative_error < tolerance);

  // Every excluded entry carries a reason.
  for (const SensitivityCheckEntry& entry : check.entries) {
    if (entry.excluded) REQUIRE_FALSE(entry.exclusion_reason.empty());
  }
}

TEST_CASE("sensitivities are exact without a filter too",
          "[topopt][sensitivity][verification]") {
  TopoptFixture fixture(8, 4, 0.4, 0.2);
  const DensityFilter filter(fixture.model.mesh(), FilterType::None, 0.0);
  DesignDomain domain(fixture.model, 0.5, 0.5, {});
  Assembler assembler(fixture.model);

  SimpOptions simp;
  simp.penalty = 2.5;
  StaticAnalysisOptions options;
  options.linear.residual_tolerance = 1.0e-9;
  ComplianceObjective objective(fixture.model, assembler, filter, domain, simp,
                                options);

  const Vector x = wavy_design(fixture.model, domain);
  const SensitivityCheckResult check =
      verify_sensitivities(objective, domain, x, {}, 1.0e-4, 1.0e-5);
  REQUIRE(check.passed);
  // Without a filter the physical and design gradients coincide.
  const ObjectiveEvaluation eval = objective.evaluate(x, true);
  REQUIRE(eval.dc_dx.isApprox(eval.dc_dphysical, 1.0e-14));
}

TEST_CASE("compliance gradient is non-positive and volume gradient positive",
          "[topopt][sensitivity]") {
  TopoptFixture fixture(10, 5, 0.5, 0.25);
  const DensityFilter filter(fixture.model.mesh(), FilterType::Density, 1.5 * fixture.cell);
  DesignDomain domain(fixture.model, 0.5, 0.5, {});
  Assembler assembler(fixture.model);
  ComplianceObjective objective(fixture.model, assembler, filter, domain,
                                SimpOptions(), StaticAnalysisOptions());

  const ObjectiveEvaluation eval = objective.evaluate(domain.initial_design(), true);
  REQUIRE(eval.dc_dphysical.maxCoeff() <= 0.0);
  REQUIRE(eval.dc_dx.maxCoeff() <= 1.0e-30);
  REQUIRE(eval.dv_dx.minCoeff() > 0.0);
  REQUIRE(eval.compliance > 0.0);
  REQUIRE(eval.volume_fraction == Approx(0.5).epsilon(1.0e-12));

  // Adding material lowers compliance.
  const Vector denser = (domain.initial_design().array() + 0.1).matrix();
  REQUIRE(objective.compliance_at(denser) < eval.compliance);

  // Total strain energy equals half the compliance for this single load case.
  REQUIRE(eval.element_strain_energy.sum() ==
          Approx(0.5 * eval.compliance).epsilon(1.0e-10));

  REQUIRE_THROWS_AS(objective.evaluate(Vector::Ones(3), true), ConfigError);
}

TEST_CASE("multi-load-case compliance is the weighted sum",
          "[topopt][sensitivity]") {
  TopoptFixture fixture(10, 5, 0.5, 0.25);
  FemModel& model = fixture.model;
  REQUIRE(model.load_case_specs().size() == 1);

  // Rebuild with two load cases of unequal weight.
  StructuredMeshSpec spec;
  spec.nx = 10;
  spec.ny = 5;
  spec.lx = 0.5;
  spec.ly = 0.25;
  FemModel two(make_structured_quad_mesh(spec), default_material(), 0.01,
               StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint root;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = true;
  root.fix_y = true;
  two.constraints().push_back(root);

  for (int k = 0; k < 2; ++k) {
    LoadCaseSpec load;
    load.name = k == 0 ? "down" : "sideways";
    load.weight = k == 0 ? 3.0 : 1.0;
    PointLoadSpec p;
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector3(spec.lx, k == 0 ? 0.0 : spec.ly, 0.0);
    p.region.members.push_back(nearest);
    p.force = k == 0 ? Vector3(0.0, -500.0, 0.0) : Vector3(300.0, 0.0, 0.0);
    load.point_loads.push_back(p);
    two.load_case_specs().push_back(load);
  }
  two.finalize();

  const std::vector<Scalar> weights = two.normalised_weights();
  REQUIRE(weights.size() == 2);
  REQUIRE(weights[0] == Approx(0.75));
  REQUIRE(weights[1] == Approx(0.25));

  const DensityFilter filter(two.mesh(), FilterType::Density, 1.5 * fixture.cell);
  DesignDomain domain(two, 0.5, 0.5, {});
  Assembler assembler(two);
  ComplianceObjective objective(two, assembler, filter, domain, SimpOptions(),
                                StaticAnalysisOptions());
  const ObjectiveEvaluation eval = objective.evaluate(domain.initial_design(), true);

  REQUIRE(eval.load_case_compliance.size() == 2);
  REQUIRE(eval.compliance ==
          Approx(0.75 * eval.load_case_compliance[0] +
                 0.25 * eval.load_case_compliance[1])
              .epsilon(1.0e-12));

  // And the gradient of the weighted objective is verified by FD as well.
  const SensitivityCheckResult check =
      verify_sensitivities(objective, domain, domain.initial_design(), {}, 1.0e-4,
                           1.0e-5);
  REQUIRE(check.passed);
}

TEST_CASE("optimality-criteria update satisfies the volume constraint",
          "[topopt][optimizer][verification]") {
  TopoptFixture fixture(16, 8, 0.8, 0.4);
  const DensityFilter filter(fixture.model.mesh(), FilterType::Density, 1.5 * fixture.cell);
  DesignDomain domain(fixture.model, 0.35, 0.35, {});
  Assembler assembler(fixture.model);
  ComplianceObjective objective(fixture.model, assembler, filter, domain,
                                SimpOptions(), StaticAnalysisOptions());

  const auto volume_of = [&](const Vector& xi) {
    Vector rho = filter.to_physical(xi);
    rho = rho.cwiseMax(0.0).cwiseMin(1.0);
    return domain.volume_of(rho);
  };

  Vector x = domain.initial_design();
  OptimalityCriteriaOptions oc;
  for (int iteration = 0; iteration < 5; ++iteration) {
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    const OptimalityCriteriaStep step =
        optimality_criteria_update(domain, x, eval.dc_dx, eval.dv_dx, volume_of, oc);

    REQUIRE(step.volume_converged);
    REQUIRE(step.achieved_volume ==
            Approx(domain.volume_target()).epsilon(1.0e-8));
    REQUIRE(step.lambda > 0.0);
    REQUIRE(step.max_change <= oc.move_limit + 1.0e-12);

    // Bounds are respected.
    for (Index e = 0; e < domain.num_elements(); ++e) {
      REQUIRE(step.x(e) >= domain.lower_bounds()(e) - 1.0e-12);
      REQUIRE(step.x(e) <= domain.upper_bounds()(e) + 1.0e-12);
      REQUIRE(std::abs(step.x(e) - x(e)) <= oc.move_limit + 1.0e-12);
    }
    x = step.x;
  }
}

TEST_CASE("optimality criteria validates its inputs",
          "[topopt][optimizer][diagnostics]") {
  TopoptFixture fixture(6, 3, 0.3, 0.15);
  const DensityFilter filter(fixture.model.mesh(), FilterType::None, 0.0);
  DesignDomain domain(fixture.model, 0.5, 0.5, {});
  const auto volume_of = [&](const Vector& xi) { return domain.volume_of(xi); };

  const Vector x = domain.initial_design();
  const Vector dc = Vector::Constant(domain.num_elements(), -1.0);
  const Vector dv = domain.element_volumes();

  OptimalityCriteriaOptions bad_move;
  bad_move.move_limit = 0.0;
  REQUIRE_THROWS_AS(
      optimality_criteria_update(domain, x, dc, dv, volume_of, bad_move), ConfigError);

  OptimalityCriteriaOptions bad_damping;
  bad_damping.damping = 0.0;
  REQUIRE_THROWS_AS(
      optimality_criteria_update(domain, x, dc, dv, volume_of, bad_damping),
      ConfigError);

  OptimalityCriteriaOptions ok;
  REQUIRE_THROWS_AS(
      optimality_criteria_update(domain, x, Vector::Ones(3), dv, volume_of, ok),
      ConfigError);

  Vector nan_gradient = dc;
  nan_gradient(0) = std::numeric_limits<Scalar>::quiet_NaN();
  REQUIRE_THROWS_AS(
      optimality_criteria_update(domain, x, nan_gradient, dv, volume_of, ok),
      SolverError);
}

TEST_CASE("the objective-stall measure is the spread over its window", "[topopt]") {
  // Too short a history: no verdict.
  REQUIRE(std::isinf(objective_spread({1.0, 0.9}, 5)));

  // Monotone: exactly the two-point change over the window.
  std::vector<Scalar> monotone;
  for (int k = 0; k < 30; ++k) monotone.push_back(1.0 + 1.0 / (1.0 + k));
  const Scalar two_point = std::abs(monotone.back() - monotone[monotone.size() - 1 - 20]) /
                           std::abs(monotone.back());
  REQUIRE(objective_spread(monotone, 20) == two_point);

  // A cycle of period 5 returns to the same compliance every 20 iterations: a
  // two-point test over a window of 20 sees no change at all, the spread sees
  // the full amplitude of the cycle.
  const std::array<Scalar, 5> cycle = {0.0900, 0.0915, 0.0940, 0.0910, 0.0897};
  std::vector<Scalar> oscillating;
  for (int k = 0; k < 60; ++k) oscillating.push_back(cycle[static_cast<std::size_t>(k % 5)]);
  const Scalar last = oscillating.back();
  REQUIRE(std::abs(last - oscillating[oscillating.size() - 1 - 20]) == 0.0);
  REQUIRE(objective_spread(oscillating, 20) == Approx((0.0940 - 0.0897) / last));
  REQUIRE(objective_spread(oscillating, 20) > 1.0e-2);
}

TEST_CASE("grey level measures how binary a density field is", "[topopt]") {
  REQUIRE(gray_level(Vector::Zero(10)) == Approx(0.0));
  REQUIRE(gray_level(Vector::Ones(10)) == Approx(0.0));
  REQUIRE(gray_level(Vector::Constant(10, 0.5)) == Approx(1.0));
  REQUIRE(gray_level(Vector::Constant(10, 0.25)) == Approx(0.75));
  REQUIRE(gray_level(Vector()) == Approx(0.0));
}

TEST_CASE("end-to-end optimisation converges and respects the volume constraint",
          "[topopt][optimizer][verification]") {
  TopoptFixture fixture(30, 15, 0.6, 0.3);
  const DensityFilter filter(fixture.model.mesh(), FilterType::Density, 1.5 * fixture.cell);
  DesignDomain domain(fixture.model, 0.4, -1.0, {});
  Assembler assembler(fixture.model);

  TopologyOptimizerOptions options;
  options.max_iterations = 250;
  options.change_tolerance = 5.0e-3;
  options.history_stride = 10;
  TopologyOptimizer optimizer(fixture.model, assembler, filter, domain, options);
  const TopologyOptimizationResult result = optimizer.run();

  INFO("iterations " << result.iterations << ", compliance " << result.compliance
                     << ", volume fraction " << result.volume_fraction << ", grey "
                     << result.gray_level);
  REQUIRE(result.converged);
  REQUIRE(result.warnings.empty());

  // The volume constraint holds to the bisection tolerance.
  REQUIRE(result.volume_fraction == Approx(0.4).epsilon(1.0e-6));
  REQUIRE(std::abs(result.volume_constraint_violation) < 1.0e-6);

  // Compliance improved substantially over the uniform start.
  REQUIRE(result.history.front().compliance > result.compliance);
  REQUIRE(result.history.front().compliance / result.compliance > 1.5);

  // The design is mostly black and white, and inside its bounds. A density
  // filter of radius 1.5 cells leaves a genuinely grey boundary layer, and on a
  // 30 x 15 mesh that layer is a large fraction of the domain, so the measure
  // settles near 0.3 rather than near 0. docs/topology_optimization.md explains
  // how the grey level should be read.
  REQUIRE(result.gray_level < 0.40);
  REQUIRE(result.physical_density.minCoeff() >= 0.0);
  REQUIRE(result.physical_density.maxCoeff() <= 1.0 + 1.0e-12);

  // Every volume record along the way also satisfies the constraint.
  for (const TopologyIteration& it : result.history) {
    REQUIRE(it.volume_fraction == Approx(0.4).epsilon(1.0e-5));
    REQUIRE(it.volume_converged);
  }

  // The compliance history is eventually monotone decreasing: OC is a
  // fixed-point update with a move limit, so early iterations may rise. Check
  // that the last third is non-increasing to within a small tolerance.
  const std::size_t start = 2 * result.history.size() / 3;
  for (std::size_t i = start + 1; i < result.history.size(); ++i) {
    REQUIRE(result.history[i].compliance <=
            result.history[i - 1].compliance * (1.0 + 1.0e-3));
  }

  // Snapshots were recorded for the animation.
  REQUIRE_FALSE(result.snapshots.empty());
  REQUIRE(result.snapshots.size() == result.snapshot_iterations.size());
  REQUIRE(result.snapshot_iterations.back() == result.iterations);
  REQUIRE(result.linear_solves >= static_cast<Index>(result.iterations));
}

TEST_CASE("passive regions are honoured by the optimiser",
          "[topopt][optimizer][verification]") {
  TopoptFixture fixture(24, 12, 0.6, 0.3);
  const DensityFilter filter(fixture.model.mesh(), FilterType::Density, 1.5 * fixture.cell);

  PassiveRegionSpec solid;
  solid.region.name = "attachment_pad";
  Selector sbox;
  sbox.kind = SelectorKind::Box;
  sbox.xmin = 0.55;
  sbox.ymin = 0.10;
  sbox.ymax = 0.20;
  solid.region.members.push_back(sbox);
  solid.solid = true;

  PassiveRegionSpec hole;
  hole.region.name = "bolt_hole";
  Selector circle;
  circle.kind = SelectorKind::Circle;
  circle.center = Vector3(0.15, 0.15, 0.0);
  circle.radius = 0.05;
  hole.region.members.push_back(circle);
  hole.solid = false;

  DesignDomain domain(fixture.model, 0.45, -1.0, {solid, hole});
  Assembler assembler(fixture.model);

  TopologyOptimizerOptions options;
  options.max_iterations = 60;
  options.change_tolerance = 1.0e-2;
  TopologyOptimizer optimizer(fixture.model, assembler, filter, domain, options);
  const TopologyOptimizationResult result = optimizer.run();

  for (Index e = 0; e < domain.num_elements(); ++e) {
    const PassiveTag tag = domain.tags()[static_cast<std::size_t>(e)];
    if (tag == PassiveTag::Solid) {
      REQUIRE(result.design(e) == Approx(1.0));
    } else if (tag == PassiveTag::Void) {
      REQUIRE(result.design(e) == Approx(0.0));
    }
  }
  REQUIRE(result.volume_fraction == Approx(0.45).epsilon(1.0e-5));
}

TEST_CASE("SIMP continuation ramps the penalty and is recorded",
          "[topopt][optimizer]") {
  TopoptFixture fixture(20, 10, 0.6, 0.3);
  const DensityFilter filter(fixture.model.mesh(), FilterType::Density, 1.5 * fixture.cell);
  DesignDomain domain(fixture.model, 0.4, -1.0, {});
  Assembler assembler(fixture.model);

  TopologyOptimizerOptions options;
  options.max_iterations = 45;
  options.continuation_steps = 3;
  options.penalty_start = 1.0;
  options.continuation_iterations = 10;
  options.simp.penalty = 3.0;
  options.change_tolerance = 1.0e-3;

  TopologyOptimizer optimizer(fixture.model, assembler, filter, domain, options);
  const TopologyOptimizationResult result = optimizer.run();

  REQUIRE(result.history.front().penalty == Approx(1.0));
  REQUIRE(result.history.back().penalty == Approx(3.0));
  // The penalty is non-decreasing along the history.
  for (std::size_t i = 1; i < result.history.size(); ++i) {
    REQUIRE(result.history[i].penalty >= result.history[i - 1].penalty);
  }
  REQUIRE(result.volume_fraction == Approx(0.4).epsilon(1.0e-5));
}

TEST_CASE("non-convergence is reported instead of hidden",
          "[topopt][optimizer][diagnostics]") {
  TopoptFixture fixture(20, 10, 0.6, 0.3);
  const DensityFilter filter(fixture.model.mesh(), FilterType::Density, 1.5 * fixture.cell);
  DesignDomain domain(fixture.model, 0.4, -1.0, {});
  Assembler assembler(fixture.model);

  TopologyOptimizerOptions options;
  options.max_iterations = 3;          // far too few
  options.change_tolerance = 1.0e-12;  // unreachable
  TopologyOptimizer optimizer(fixture.model, assembler, filter, domain, options);
  const TopologyOptimizationResult result = optimizer.run();

  REQUIRE_FALSE(result.converged);
  REQUIRE_FALSE(result.warnings.empty());
  REQUIRE(result.warnings.front().find("iteration cap") != std::string::npos);
  REQUIRE(result.iterations == 3);
  // Three compliances cannot fill the default six-value window, and the
  // warning says so instead of quoting an infinite spread.
  REQUIRE_FALSE(std::isfinite(result.final_objective_change));
  REQUIRE(result.warnings.front().find("could not be measured") != std::string::npos);
  REQUIRE(result.warnings.front().find("inf") == std::string::npos);
  // The volume constraint still holds on the returned iterate.
  REQUIRE(result.volume_fraction == Approx(0.4).epsilon(1.0e-5));
}

TEST_CASE("optimizer options are validated", "[topopt][optimizer][diagnostics]") {
  TopoptFixture fixture(6, 3, 0.3, 0.15);
  const DensityFilter filter(fixture.model.mesh(), FilterType::None, 0.0);
  DesignDomain domain(fixture.model, 0.5, -1.0, {});
  Assembler assembler(fixture.model);

  const auto make = [&](TopologyOptimizerOptions options) {
    return TopologyOptimizer(fixture.model, assembler, filter, domain, options);
  };

  TopologyOptimizerOptions bad;
  bad.max_iterations = 0;
  REQUIRE_THROWS_AS(make(bad), ConfigError);

  bad = TopologyOptimizerOptions();
  bad.change_tolerance = 0.0;
  REQUIRE_THROWS_AS(make(bad), ConfigError);

  bad = TopologyOptimizerOptions();
  bad.continuation_steps = 0;
  REQUIRE_THROWS_AS(make(bad), ConfigError);

  bad = TopologyOptimizerOptions();
  bad.continuation_steps = 3;
  bad.penalty_start = 5.0;
  bad.simp.penalty = 3.0;
  REQUIRE_THROWS_AS(make(bad), ConfigError);

  bad = TopologyOptimizerOptions();
  bad.interpretation_threshold = 1.0;
  REQUIRE_THROWS_AS(make(bad), ConfigError);
}

TEST_CASE("filter radius controls the minimum feature size",
          "[topopt][filter][optimizer]") {
  // A larger filter radius must produce a smoother (less fine-grained) design.
  // "Smoothness" is measured by the total variation of the density field
  // between neighbouring cells in the structured grid.
  const auto run = [](Scalar radius_cells) {
    TopoptFixture fixture(30, 15, 0.6, 0.3);
    const DensityFilter filter(fixture.model.mesh(), FilterType::Density,
                               radius_cells * fixture.cell);
    DesignDomain domain(fixture.model, 0.4, -1.0, {});
    Assembler assembler(fixture.model);
    TopologyOptimizerOptions options;
    options.max_iterations = 80;
    options.change_tolerance = 5.0e-3;
    options.history_stride = 0;
    TopologyOptimizer optimizer(fixture.model, assembler, filter, domain, options);
    const TopologyOptimizationResult result = optimizer.run();

    const StructuredGridInfo& info = *fixture.model.mesh().structured_info();
    Scalar variation = 0.0;
    for (Index j = 0; j < info.ny; ++j) {
      for (Index i = 0; i + 1 < info.nx; ++i) {
        variation += std::abs(result.physical_density(
                                  structured_element_index(info, i + 1, j)) -
                              result.physical_density(
                                  structured_element_index(info, i, j)));
      }
    }
    return variation;
  };

  const Scalar sharp = run(1.2);
  const Scalar smooth = run(3.0);
  INFO("total variation: r = 1.2 cells -> " << sharp << ", r = 3.0 cells -> "
                                            << smooth);
  REQUIRE(smooth < sharp);
}
