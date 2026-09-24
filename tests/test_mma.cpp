/// \file test_mma.cpp
/// \brief MMA on problems with known solutions, the von Mises matrix, the
///        p-norm aggregation, the adjoint stress sensitivities against
///        central differences (2-D and 3-D), and stress-constrained runs.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/Mma.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/StressConstraint.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <random>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

/// Central-difference check of a constraint gradient on every free element.
struct GradientCheck {
  Scalar max_relative_error = 0.0;
  Scalar directional_relative_error = 0.0;
  Index tested = 0;
};

GradientCheck check_stress_gradient(ComplianceObjective& objective,
                                    const StressConstraint& stress,
                                    const DesignDomain& domain, const Vector& x,
                                    std::size_t load_case, Scalar scale, Scalar step) {
  const ObjectiveEvaluation base = objective.evaluate(x, true);
  const StressEvaluation se = stress.evaluate(objective, base, load_case, scale, true);
  const auto value_at = [&](const Vector& xx) {
    const ObjectiveEvaluation e = objective.evaluate(xx, false);
    return stress.evaluate(objective, e, load_case, scale, false).constraint;
  };

  GradientCheck out;
  Vector xp = x;
  Vector xm = x;
  Vector direction = Vector::Zero(x.size());
  // Entries far below the gradient's own scale are judged against that scale:
  // a central difference of two O(1) constraint values carries an absolute
  // round-off floor, so a per-entry relative error would mean nothing there.
  const Scalar floor = 1.0e-3 * se.dg_dx.cwiseAbs().maxCoeff();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    if (!domain.is_free(e)) continue;
    if (x(e) - step < domain.lower_bounds()(e) || x(e) + step > domain.upper_bounds()(e)) {
      continue;
    }
    xp(e) = x(e) + step;
    xm(e) = x(e) - step;
    const Scalar fd = (value_at(xp) - value_at(xm)) / (2.0 * step);
    xp(e) = x(e);
    xm(e) = x(e);
    const Scalar an = se.dg_dx(e);
    const Scalar rel = std::abs(fd - an) / std::max({std::abs(fd), std::abs(an), floor});
    out.max_relative_error = std::max(out.max_relative_error, rel);
    direction(e) = an;
    ++out.tested;
  }
  const Scalar dnorm = direction.norm();
  if (dnorm > 0.0) {
    direction /= dnorm;
    const Scalar h = step;
    const Scalar fd = (value_at((x + h * direction).eval()) -
                       value_at((x - h * direction).eval())) / (2.0 * h);
    const Scalar an = se.dg_dx.dot(direction);
    out.directional_relative_error =
        std::abs(fd - an) / std::max({std::abs(fd), std::abs(an), 1.0e-30});
  }
  return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// MMA
// ---------------------------------------------------------------------------

TEST_CASE("MMA solves a separable problem with a known optimum",
          "[mma][optimizer][verification]") {
  // min sum x_i^2 subject to 1 - sum x_i <= 0 on [0, 1]^n: x_i = 1/n.
  const Index n = 12;
  MmaOptimizer mma(n, 1, Vector::Zero(n), Vector::Ones(n));
  Vector x = Vector::Constant(n, 0.9);
  MmaStep step;
  REQUIRE(mma.options().constraint_scale_cap == Approx(10.0));
  for (int it = 0; it < 60; ++it) {
    const Scalar f0 = x.squaredNorm();
    const Vector df0 = 2.0 * x;
    Vector fval(1);
    fval(0) = 1.0 - x.sum();
    Matrix dfdx = Matrix::Constant(1, n, -1.0);
    step = mma.update(x, f0, df0, fval, dfdx);
    REQUIRE(step.subproblem_converged);
    REQUIRE(step.lambda(0) >= 0.0);
    x = step.x;
    if (step.max_change < 1.0e-9) break;
  }
  REQUIRE((x.array() - 1.0 / n).abs().maxCoeff() < 1.0e-6);
  REQUIRE(x.sum() == Approx(1.0).epsilon(1.0e-6));
  REQUIRE(step.constraint_scale(0) == 1.0);  // well scaled: untouched
  // KKT: the multiplier equals 2/n (gradient balance 2 x_i = lambda).
  REQUIRE(step.lambda(0) == Approx(2.0 / n).epsilon(1.0e-3));
  REQUIRE(mma.iteration() > 2);
}

TEST_CASE("MMA handles two active constraints and non-uniform bounds",
          "[mma][optimizer][verification]") {
  // min (x1-1)^2 + (x2-1)^2 s.t. x1 + x2 <= 1 and x1 <= 0.3 on [0, 2]^2.
  // Solution (0.3, 0.7), both constraints active, f = 0.58.
  MmaOptimizer mma(2, 2, Vector::Zero(2), Vector::Constant(2, 2.0));
  Vector x(2);
  x << 1.5, 0.2;
  MmaStep step;
  for (int it = 0; it < 80; ++it) {
    const Scalar f0 = std::pow(x(0) - 1.0, 2) + std::pow(x(1) - 1.0, 2);
    Vector df0(2);
    df0 << 2.0 * (x(0) - 1.0), 2.0 * (x(1) - 1.0);
    Vector fval(2);
    fval << x(0) + x(1) - 1.0, x(0) - 0.3;
    Matrix dfdx(2, 2);
    dfdx << 1.0, 1.0, 1.0, 0.0;
    step = mma.update(x, f0, df0, fval, dfdx);
    x = step.x;
    if (step.max_change < 1.0e-10) break;
  }
  REQUIRE(x(0) == Approx(0.3).margin(1.0e-5));
  REQUIRE(x(1) == Approx(0.7).margin(1.0e-5));
  REQUIRE(std::pow(x(0) - 1.0, 2) + std::pow(x(1) - 1.0, 2) == Approx(0.58).epsilon(1.0e-4));
  // Both multipliers positive at a vertex solution; KKT: 2(x1-1) + l1 + l2 = 0,
  // 2(x2-1) + l1 = 0 -> l1 = 0.6, l2 = 0.8.
  REQUIRE(step.lambda(0) == Approx(0.6).epsilon(2.0e-2));
  REQUIRE(step.lambda(1) == Approx(0.8).epsilon(2.0e-2));
  REQUIRE(step.y.maxCoeff() < 1.0e-6);
}

TEST_CASE("MMA scales a badly violated constraint down for its subproblem",
          "[mma][optimizer]") {
  // min x1 + x2 s.t. 1e5 (1 - x1 - x2) <= 0 on [0, 1]^2: the constraint is the
  // volume-like one scaled by 1e5, which would defeat an absolute residual
  // target. The cap keeps the subproblem solvable and the answer unchanged:
  // x1 + x2 = 1 with the multiplier 1e-5 (gradient balance 1 = 1e5 lambda).
  MmaOptimizer mma(2, 1, Vector::Zero(2), Vector::Ones(2));
  Vector x = Vector::Constant(2, 0.9);
  MmaStep step;
  for (int it = 0; it < 80; ++it) {
    Vector fval(1);
    fval(0) = 1.0e5 * (1.0 - x.sum());
    step = mma.update(x, x.sum(), Vector::Ones(2), fval, Matrix::Constant(1, 2, -1.0e5));
    x = step.x;
    if (step.max_change < 1.0e-10) break;
  }
  REQUIRE(step.constraint_scale(0) < 1.0);
  REQUIRE(x.sum() == Approx(1.0).epsilon(1.0e-6));
  REQUIRE(step.lambda(0) == Approx(1.0e-5).epsilon(0.05));
}

TEST_CASE("MMA rejects malformed input", "[mma][diagnostics]") {
  REQUIRE_THROWS_AS(MmaOptimizer(3, 1, Vector::Zero(3), Vector::Zero(3)), ConfigError);
  REQUIRE_THROWS_AS(MmaOptimizer(0, 1, Vector(), Vector()), ConfigError);
  MmaOptions bad;
  bad.move_limit = 0.0;
  REQUIRE_THROWS_AS(MmaOptimizer(2, 1, Vector::Zero(2), Vector::Ones(2), bad), ConfigError);
  MmaOptimizer mma(2, 1, Vector::Zero(2), Vector::Ones(2));
  REQUIRE_THROWS_AS(mma.update(Vector::Zero(3), 0.0, Vector::Zero(2), Vector::Zero(1),
                               Matrix::Zero(1, 2)),
                    ConfigError);
  Vector nan_grad = Vector::Zero(2);
  nan_grad(0) = std::numeric_limits<Scalar>::quiet_NaN();
  REQUIRE_THROWS_AS(mma.update(Vector::Constant(2, 0.5), 0.0, nan_grad, Vector::Zero(1),
                               Matrix::Zero(1, 2)),
                    SolverError);
  REQUIRE_THROWS_AS(parse_optimizer_method("simplex"), ConfigError);
  REQUIRE(parse_optimizer_method("mma") == OptimizerMethod::MMA);
  REQUIRE(to_string(OptimizerMethod::OptimalityCriteria) == "oc");
}

// ---------------------------------------------------------------------------
// Stress constraint ingredients
// ---------------------------------------------------------------------------

TEST_CASE("the von Mises matrix reproduces von_mises() in every stress state",
          "[stress][mma][verification]") {
  std::mt19937 rng(3u);
  std::uniform_real_distribution<Scalar> dist(-1.0e8, 1.0e8);
  const Scalar nu = 0.3;
  for (StressState state : {StressState::PlaneStress, StressState::PlaneStrain,
                            StressState::ThreeDimensional}) {
    const Matrix v = von_mises_matrix(state, nu);
    const int nv = voigt_components(stress_state_dimension(state));
    REQUIRE(v.rows() == nv);
    REQUIRE((v - v.transpose()).cwiseAbs().maxCoeff() == Approx(0.0));
    for (int trial = 0; trial < 20; ++trial) {
      Vector s(nv);
      for (int i = 0; i < nv; ++i) s(i) = dist(rng);
      REQUIRE(std::sqrt(std::max(s.dot(v * s), 0.0)) ==
              Approx(von_mises(s, state, nu)).epsilon(1.0e-12));
    }
  }
}

TEST_CASE("stress constraint options are validated against the SIMP law",
          "[stress][diagnostics]") {
  SimpOptions simp;
  simp.penalty = 3.0;
  StressConstraintOptions o;
  o.enabled = true;
  o.limit = 0.0;
  REQUIRE_THROWS_AS(o.validate(simp), ConfigError);
  o.limit = 1.0e8;
  o.relaxation = 3.0;  // must be below p
  REQUIRE_THROWS_AS(o.validate(simp), ConfigError);
  o.relaxation = 0.5;
  o.p_norm = 1.0;
  REQUIRE_THROWS_AS(o.validate(simp), ConfigError);
  o.p_norm = 8.0;
  REQUIRE_NOTHROW(o.validate(simp));
  o.enabled = false;
  o.limit = -1.0;
  REQUIRE_NOTHROW(o.validate(simp));  // disabled options are not checked
}

TEST_CASE("the p-norm aggregation bounds the maximum relaxed stress from above",
          "[stress][mma][verification]") {
  FemModel model = make_small_plate(8, 4);
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density, 0.15);
  DesignDomain domain(model, 0.6, 0.6, {});
  StressConstraintOptions so;
  so.enabled = true;
  so.limit = 5.0e7;
  so.p_norm = 8.0;
  const StressConstraint stress(model, assembler, filter, so);
  ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                StaticAnalysisOptions());

  const ObjectiveEvaluation eval = objective.evaluate(domain.initial_design(), true);
  const StressEvaluation se = stress.evaluate(objective, eval, 0, 1.0, false);
  const Index ne = model.mesh().num_elements();
  REQUIRE(se.p_norm_ratio >= se.max_relaxed_ratio);
  REQUIRE(se.p_norm_ratio <= std::pow(static_cast<Scalar>(ne), 1.0 / 8.0) * se.max_relaxed_ratio);
  REQUIRE(se.constraint == Approx(se.p_norm_ratio - 1.0));
  // Relaxed stress is rho^q times the solid stress.
  for (Index e = 0; e < ne; ++e) {
    REQUIRE(se.relaxed_von_mises(e) ==
            Approx(std::pow(eval.physical_density(e), 0.5) * se.solid_von_mises(e)));
  }
  // The solid stress at the centre agrees with the recovery's element average
  // on a uniform (parallelogram) mesh, where B is linear in xi and eta.
  const StressField field = recover_stresses(model, assembler, eval.displacements[0]);
  for (Index e = 0; e < ne; ++e) {
    REQUIRE(se.solid_von_mises(e) == Approx(field.element_solid_von_mises(e)).epsilon(1.0e-10));
  }
  // The adaptive scale moves towards max / p-norm.
  const Scalar next = stress.next_scale(1.0, se);
  const Scalar target = se.max_relaxed_ratio / se.p_norm_ratio;
  REQUIRE(next == Approx(0.5 * target + 0.5));
  REQUIRE(next < 1.0);
}

TEST_CASE("stress constraint gradients match central differences in 2-D",
          "[stress][mma][sensitivity][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 10;
  spec.ny = 5;
  spec.lx = 0.5;
  spec.ly = 0.25;
  FemModel model(make_structured_quad_mesh(spec), default_material(), 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint root;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = root.fix_y = true;
  model.constraints().push_back(root);
  for (int k = 0; k < 2; ++k) {
    LoadCaseSpec load;
    load.name = k == 0 ? "down" : "side";
    PointLoadSpec p;
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector3(spec.lx, k == 0 ? 0.0 : spec.ly, 0.0);
    p.region.members.push_back(nearest);
    p.force = k == 0 ? Vector3(0.0, -800.0, 0.0) : Vector3(400.0, 0.0, 0.0);
    load.point_loads.push_back(p);
    model.load_case_specs().push_back(load);
  }
  model.finalize();
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density, 1.5 * 0.05);

  PassiveRegionSpec solid;
  solid.region.name = "pad";
  Selector pad;
  pad.kind = SelectorKind::Box;
  pad.xmax = 0.05;
  pad.ymax = 0.05;
  solid.region.members.push_back(pad);
  solid.solid = true;
  DesignDomain domain(model, 0.5, 0.5, {solid});

  SimpOptions simp;
  simp.penalty = 3.0;
  StaticAnalysisOptions options;
  options.linear.residual_tolerance = 1.0e-9;
  ComplianceObjective objective(model, assembler, filter, domain, simp, options);
  StressConstraintOptions so;
  so.enabled = true;
  so.limit = 4.0e7;
  so.p_norm = 6.0;
  so.relaxation = 0.5;
  const StressConstraint stress(model, assembler, filter, so);

  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    if (!domain.is_free(e)) continue;
    const Vector3 c = model.mesh().element_centroid(e);
    x(e) = 0.45 + 0.30 * std::sin(11.0 * c.x()) * std::cos(7.0 * c.y());
  }
  domain.clamp(x);

  for (std::size_t l = 0; l < 2; ++l) {
    const GradientCheck check = check_stress_gradient(objective, stress, domain, x, l, 0.7, 1.0e-5);
    INFO("load case " << l << ": max relative error " << check.max_relative_error
                      << ", directional " << check.directional_relative_error
                      << ", tested " << check.tested);
    REQUIRE(check.tested > 30);
    REQUIRE(check.max_relative_error < 1.0e-5);
    REQUIRE(check.directional_relative_error < 1.0e-6);
  }

  // Without a filter the physical and design gradients coincide.
  const DensityFilter none(model.mesh(), FilterType::None, 0.0);
  ComplianceObjective plain(model, assembler, none, domain, simp, options);
  const StressConstraint stress_plain(model, assembler, none, so);
  const ObjectiveEvaluation e0 = plain.evaluate(x, true);
  const StressEvaluation s0 = stress_plain.evaluate(plain, e0, 0, 1.0, true);
  REQUIRE(s0.dg_dx.isApprox(s0.dg_dphysical, 1.0e-14));
  const GradientCheck plain_check =
      check_stress_gradient(plain, stress_plain, domain, x, 0, 1.0, 1.0e-5);
  REQUIRE(plain_check.max_relative_error < 1.0e-5);

  // The sensitivity filter is refused: it has no exact chain rule.
  const DensityFilter heuristic(model.mesh(), FilterType::Sensitivity, 0.075);
  REQUIRE_THROWS_AS(StressConstraint(model, assembler, heuristic, so), ConfigError);
}

TEST_CASE("stress constraint gradients match central differences in 3-D",
          "[stress][mma][sensitivity][solid][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 6;
  spec.ny = 3;
  spec.nz = 2;
  spec.lx = 0.6;
  spec.ly = 0.3;
  spec.lz = 0.2;
  FemModel model(make_structured_hex_mesh(spec), default_material(), 1.0,
                 StressState::ThreeDimensional, IntegrationOptions());
  DisplacementConstraint root;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = root.fix_y = root.fix_z = true;
  model.constraints().push_back(root);
  LoadCaseSpec load;
  PointLoadSpec p;
  Selector tip;
  tip.kind = SelectorKind::Box;
  tip.xmin = spec.lx;
  tip.ymax = 0.0;
  p.region.members.push_back(tip);
  p.force = Vector3(0.0, -500.0, 150.0);
  load.point_loads.push_back(p);
  model.load_case_specs().push_back(load);
  model.finalize();
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density, 0.15);
  DesignDomain domain(model, 0.5, 0.5, {});
  StaticAnalysisOptions options;
  options.linear.residual_tolerance = 1.0e-9;
  ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(), options);
  StressConstraintOptions so;
  so.enabled = true;
  so.limit = 2.0e7;
  so.p_norm = 8.0;
  const StressConstraint stress(model, assembler, filter, so);

  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    const Vector3 c = model.mesh().element_centroid(e);
    x(e) = 0.5 + 0.25 * std::sin(9.0 * c.x()) * std::cos(7.0 * c.y()) * std::cos(5.0 * c.z());
  }
  domain.clamp(x);
  const GradientCheck check = check_stress_gradient(objective, stress, domain, x, 0, 0.6, 1.0e-5);
  INFO("max relative error " << check.max_relative_error << ", directional "
                             << check.directional_relative_error);
  REQUIRE(check.tested == domain.num_elements());
  REQUIRE(check.max_relative_error < 1.0e-5);
  REQUIRE(check.directional_relative_error < 1.0e-6);
}

// ---------------------------------------------------------------------------
// Optimisation runs
// ---------------------------------------------------------------------------

TEST_CASE("MMA reproduces the optimality-criteria optimum for compliance alone",
          "[mma][optimizer][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 24;
  spec.ny = 12;
  spec.lx = 0.6;
  spec.ly = 0.3;
  const auto build = [&]() {
    FemModel model(make_structured_quad_mesh(spec), default_material(), 0.01,
                   StressState::PlaneStress, IntegrationOptions());
    DisplacementConstraint root;
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.0;
    root.region.members.push_back(box);
    root.fix_x = root.fix_y = true;
    model.constraints().push_back(root);
    LoadCaseSpec load;
    PointLoadSpec p;
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector3(spec.lx, 0.0, 0.0);
    p.region.members.push_back(nearest);
    p.force = Vector3(0.0, -500.0, 0.0);
    load.point_loads.push_back(p);
    model.load_case_specs().push_back(load);
    model.finalize();
    return model;
  };
  FemModel model = build();
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density, 1.5 * 0.025);
  DesignDomain domain(model, 0.4, -1.0, {});

  TopologyOptimizerOptions oc;
  oc.max_iterations = 200;
  oc.change_tolerance = 5.0e-3;
  oc.history_stride = 0;
  const TopologyOptimizationResult ref =
      TopologyOptimizer(model, assembler, filter, domain, oc).run();
  REQUIRE(ref.converged);

  TopologyOptimizerOptions mma = oc;
  mma.method = OptimizerMethod::MMA;
  mma.objective_tolerance = 1.0e-4;
  mma.objective_window = 10;
  const TopologyOptimizationResult got =
      TopologyOptimizer(model, assembler, filter, domain, mma).run();
  INFO("OC compliance " << ref.compliance << " J in " << ref.iterations
                        << " iterations; MMA " << got.compliance << " J in "
                        << got.iterations << " (" << got.stop_reason << ")");
  REQUIRE(got.method == OptimizerMethod::MMA);
  REQUIRE(got.converged);
  REQUIRE(got.feasible);
  REQUIRE(std::abs(got.volume_constraint_violation) <= mma.constraint_tolerance);
  // Same problem, same local-optimum neighbourhood: within a few percent.
  REQUIRE(std::abs(got.compliance - ref.compliance) / ref.compliance < 0.05);
  REQUIRE(got.stress.empty());
  REQUIRE_FALSE(got.stress_constrained);
  for (const TopologyIteration& it : got.history) {
    REQUIRE(it.bisections > 0);  // MMA Newton iterations
    REQUIRE(it.max_stress_ratio == 0.0);
  }
}

namespace {

/// The L-bracket: a square domain whose upper-right quadrant is passive void,
/// clamped along the top of the left arm and loaded at the tip of the lower
/// arm through a passive solid pad. The re-entrant corner is the classic
/// design-dependent stress concentration; a point load through a solid pad
/// keeps the load-application stress design independent.
struct LBracket {
  FemModel model;
  DensityFilter filter;
  DesignDomain domain;

  static FemModel make_model() {
    StructuredMeshSpec spec;
    spec.nx = 32;
    spec.ny = 32;
    spec.lx = 0.4;
    spec.ly = 0.4;
    FemModel m(make_structured_quad_mesh(spec), default_material(), 0.01,
               StressState::PlaneStress, IntegrationOptions());
    DisplacementConstraint top;
    Selector box;
    box.kind = SelectorKind::Box;
    box.ymin = 0.4;
    box.xmax = 0.16;
    top.region.members.push_back(box);
    top.fix_x = top.fix_y = true;
    m.constraints().push_back(top);
    LoadCaseSpec load;
    load.name = "tip";
    PointLoadSpec p;
    Selector tip;
    tip.kind = SelectorKind::Box;
    tip.xmin = 0.4;
    tip.ymin = 0.1;
    tip.ymax = 0.1375;
    p.region.members.push_back(tip);
    p.force = Vector3(0.0, -600.0, 0.0);
    load.point_loads.push_back(p);
    m.load_case_specs().push_back(load);
    m.finalize();
    return m;
  }

  static std::vector<PassiveRegionSpec> passive() {
    PassiveRegionSpec cutout;
    cutout.region.name = "cutout";
    Selector q;
    q.kind = SelectorKind::Box;
    q.xmin = 0.16;
    q.ymin = 0.16;
    cutout.region.members.push_back(q);
    cutout.solid = false;
    PassiveRegionSpec pad;
    pad.region.name = "load_pad";
    Selector b;
    b.kind = SelectorKind::Box;
    b.xmin = 0.3625;
    b.ymin = 0.0875;
    b.ymax = 0.15;
    pad.region.members.push_back(b);
    pad.solid = true;
    return {cutout, pad};
  }

  LBracket()
      : model(make_model()),
        filter(model.mesh(), FilterType::Density, 1.5 * 0.0125),
        domain(model, 0.35, -1.0, passive()) {}
};

}  // namespace

TEST_CASE("a stress-constrained run rounds the re-entrant corner at the price of compliance",
          "[stress][mma][optimizer][verification]") {
  LBracket lb;
  Assembler assembler(lb.model);

  TopologyOptimizerOptions base;
  base.method = OptimizerMethod::MMA;
  base.max_iterations = 150;
  base.change_tolerance = 5.0e-3;
  base.objective_tolerance = 1.0e-4;
  base.objective_window = 10;
  base.history_stride = 0;
  const TopologyOptimizationResult free_run =
      TopologyOptimizer(lb.model, assembler, lb.filter, lb.domain, base).run();
  REQUIRE(free_run.converged);

  // Where and how high the unconstrained design's relaxed stress peaks.
  ComplianceObjective objective(lb.model, assembler, lb.filter, lb.domain, base.simp,
                                base.analysis);
  const ObjectiveEvaluation eval = objective.evaluate(free_run.design, true);
  StressConstraintOptions probe;
  probe.enabled = true;
  probe.limit = 1.0;
  probe.p_norm = 8.0;
  const StressConstraint probe_constraint(lb.model, assembler, lb.filter, probe);
  const StressEvaluation free_stress = probe_constraint.evaluate(objective, eval, 0, 1.0, false);
  const Scalar free_peak = free_stress.max_relaxed_ratio;  // Pa, since limit = 1
  const Vector3 hotspot = lb.model.mesh().element_centroid(free_stress.max_element);
  INFO("unconstrained peak " << free_peak << " Pa at (" << hotspot.x() << ", "
                             << hotspot.y() << ")");
  // The hotspot is the re-entrant corner region, not the load pad.
  REQUIRE(hotspot.x() < 0.3);

  TopologyOptimizerOptions constrained = base;
  constrained.stress.enabled = true;
  constrained.stress.limit = 0.7 * free_peak;
  constrained.stress.p_norm = 8.0;
  constrained.stress.relaxation = 0.5;
  constrained.max_iterations = 250;
  const TopologyOptimizationResult got =
      TopologyOptimizer(lb.model, assembler, lb.filter, lb.domain, constrained).run();

  INFO("constrained: max relaxed ratio " << got.max_stress_ratio << ", compliance "
                                         << got.compliance << " J vs " << free_run.compliance
                                         << " J unconstrained, " << got.iterations
                                         << " iterations, " << got.stop_reason
                                         << ", largest constraint " << got.constraint_violation);
  REQUIRE(got.stress_constrained);
  REQUIRE(got.stress.size() == 1);
  REQUIRE(got.stress.front().load_case == "tip");
  REQUIRE(got.feasible);
  REQUIRE(got.constraint_violation <= constrained.constraint_tolerance);
  REQUIRE(got.volume_fraction <= 0.35 * (1.0 + constrained.constraint_tolerance));
  // The aggregated constraint holds and the true relaxed peak sits within the
  // p-norm's slack of the limit; the peak came down by at least a fifth and the
  // compliance went up.
  REQUIRE(got.max_stress_ratio < 1.0 + 0.05);
  REQUIRE(got.max_stress_ratio * constrained.stress.limit < 0.8 * free_peak);
  REQUIRE(got.compliance > free_run.compliance);
  REQUIRE(got.stress.front().max_solid_ratio_retained > 0.0);
  REQUIRE(got.history.back().max_stress_ratio > 0.0);
}

TEST_CASE("an unattainable stress limit is reported as infeasible, not hidden",
          "[stress][mma][optimizer][diagnostics]") {
  LBracket lb;
  Assembler assembler(lb.model);
  TopologyOptimizerOptions options;
  options.method = OptimizerMethod::MMA;
  options.max_iterations = 25;
  options.change_tolerance = 1.0e-2;
  options.history_stride = 0;
  options.stress.enabled = true;
  options.stress.limit = 1.0e3;  // Pa: three orders below any achievable stress
  const TopologyOptimizationResult got =
      TopologyOptimizer(lb.model, assembler, lb.filter, lb.domain, options).run();
  REQUIRE_FALSE(got.feasible);
  REQUIRE_FALSE(got.converged);
  REQUIRE(got.constraint_violation > 1.0);
  REQUIRE(got.max_stress_ratio > 1.0);
  bool mentioned = false;
  for (const std::string& w : got.warnings) {
    if (w.find("stress constraint is violated") != std::string::npos) mentioned = true;
  }
  REQUIRE(mentioned);
}

TEST_CASE("stress and method keys parse and are validated", "[stress][mma][io][config]") {
  const std::string head = R"({
    "name": "stress_case",
    "mesh": { "nx": 8, "ny": 4, "lx": 0.4, "ly": 0.2 },
    "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3, "density": 2700 },
    "model": { "thickness": 0.01 },
    "boundary_conditions": [ { "fix": ["x","y"], "region": { "box": { "xmax": 0.0 } } } ],
    "load_cases": [ { "point_loads": [ { "force": [0,-100], "region": { "box": { "xmin": 0.4 } } } ] } ],
    "topology": { "enabled": true, "volume_fraction": 0.5,)";
  const auto parse = [&](const std::string& tail) {
    return parse_configuration(json::parse(head + tail, "s"), "s", true);
  };
  const Configuration ok = parse(R"(
      "optimizer": { "method": "mma", "move_limit": 0.15,
                     "mma": { "asymptote_init": 0.4, "constraint_penalty": 500 },
                     "constraint_tolerance": 2e-4 },
      "stress": { "enabled": true, "limit": 1.2e8, "p_norm": 10, "relaxation": 0.6 } } })");
  REQUIRE(ok.topology.optimizer.method == OptimizerMethod::MMA);
  REQUIRE(ok.topology.optimizer.mma.move_limit == Approx(0.15));
  REQUIRE(ok.topology.optimizer.mma.asymptote_init == Approx(0.4));
  REQUIRE(ok.topology.optimizer.mma.c == Approx(500.0));
  REQUIRE(ok.topology.optimizer.constraint_tolerance == Approx(2.0e-4));
  REQUIRE(ok.topology.optimizer.stress.enabled);
  REQUIRE(ok.topology.optimizer.stress.limit == Approx(1.2e8));
  REQUIRE(ok.topology.optimizer.stress.p_norm == Approx(10.0));
  REQUIRE(ok.topology.optimizer.stress.relaxation == Approx(0.6));
  REQUIRE(ok.topology.optimizer.mma.max_inner_iterations == 500);  // the default

  // The Newton budget of the subproblem is a deck key, validated at load.
  const Configuration budget = parse(R"(
      "optimizer": { "method": "mma", "mma": { "max_newton_iterations": 750 } } } })");
  REQUIRE(budget.topology.optimizer.mma.max_inner_iterations == 750);
  REQUIRE_THROWS_AS(parse(R"(
      "optimizer": { "method": "mma", "mma": { "max_newton_iterations": 5 } } } })"),
                    ConfigError);

  // Stress without MMA, with the sensitivity filter, or with q >= p is refused.
  REQUIRE_THROWS_AS(parse(R"("stress": { "enabled": true, "limit": 1e8 } } })"), ConfigError);
  REQUIRE_THROWS_AS(parse(R"("optimizer": { "method": "mma" },
      "filter": { "type": "sensitivity" },
      "stress": { "enabled": true, "limit": 1e8 } } })"), ConfigError);
  REQUIRE_THROWS_AS(parse(R"("optimizer": { "method": "mma" },
      "stress": { "enabled": true, "limit": 1e8, "relaxation": 3.5 } } })"), ConfigError);
  REQUIRE_THROWS_AS(parse(R"("optimizer": { "method": "newton" } } })"), ConfigError);
  // Defaults: OC, no stress.
  const Configuration plain = parse(R"("optimizer": { "max_iterations": 3 } } })");
  REQUIRE(plain.topology.optimizer.method == OptimizerMethod::OptimalityCriteria);
  REQUIRE_FALSE(plain.topology.optimizer.stress.enabled);
}
