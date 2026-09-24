/// \file test_projection.cpp
/// \brief Heaviside projection: the map itself, its continuation schedule, the
///        chain rule through it (compliance, volume and stress constraint
///        against central differences, 2-D and 3-D), and optimisation runs
///        with it.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/Projection.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/StressConstraint.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {

FemModel make_block(const Mesh& mesh) {
  const int dim = mesh.dim();
  const Scalar length = mesh.bounding_box().upper.x();
  FemModel model(Mesh(mesh), default_material(), dim == 2 ? 0.01 : 1.0,
                 dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                 IntegrationOptions());
  DisplacementConstraint root;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = root.fix_y = true;
  root.fix_z = dim == 3;
  model.constraints().push_back(root);
  LoadCaseSpec load;
  load.name = "tip";
  PointLoadSpec p;
  Selector corner;
  corner.kind = SelectorKind::NearestNode;
  corner.point = Vector3(length, 0.0, 0.0);
  p.region.members.push_back(corner);
  p.force = Vector3(50.0, -800.0, dim == 3 ? 100.0 : 0.0);
  load.point_loads.push_back(p);
  model.load_case_specs().push_back(load);
  model.finalize();
  return model;
}

Vector wavy_design(const FemModel& model, const DesignDomain& domain) {
  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    if (!domain.is_free(e)) continue;
    const Vector3 c = model.mesh().element_centroid(e);
    x(e) = 0.5 + 0.3 * std::sin(11.0 * c.x()) * std::cos(7.0 * c.y() + 5.0 * c.z());
  }
  domain.clamp(x);
  return x;
}

/// Max over free elements of the relative error between an analytical
/// gradient and central differences of `value`, with a floor at 1e-3 of the
/// gradient's scale (a central difference has an absolute round-off floor).
template <typename Value>
Scalar gradient_error(const Vector& analytical, const DesignDomain& domain, const Vector& x,
                      Scalar step, const Value& value) {
  const Scalar floor = 1.0e-3 * analytical.cwiseAbs().maxCoeff();
  Scalar worst = 0.0;
  Vector xp = x;
  Vector xm = x;
  for (Index e = 0; e < domain.num_elements(); ++e) {
    if (!domain.is_free(e)) continue;
    xp(e) = x(e) + step;
    xm(e) = x(e) - step;
    const Scalar fd = (value(xp) - value(xm)) / (2.0 * step);
    xp(e) = x(e);
    xm(e) = x(e);
    worst = std::max(worst, std::abs(fd - analytical(e)) /
                                std::max({std::abs(fd), std::abs(analytical(e)), floor}));
  }
  return worst;
}

StructuredMeshSpec box(Index nx, Index ny, Index nz, Scalar lx, Scalar ly, Scalar lz) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.lx = lx;
  spec.ly = ly;
  spec.lz = lz;
  return spec;
}

}  // namespace

TEST_CASE("the tanh projection keeps its end points, is monotone and has the right slope",
          "[projection][topopt]") {
  for (Scalar eta : {0.3, 0.5, 0.7}) {
    for (Scalar beta : {0.5, 1.0, 8.0, 64.0}) {
      INFO("eta " << eta << ", beta " << beta);
      REQUIRE(heaviside_project(0.0, beta, eta) == Approx(0.0).margin(1.0e-15));
      REQUIRE(heaviside_project(1.0, beta, eta) == Approx(1.0).epsilon(1.0e-15));
      Scalar previous = -1.0;
      for (int i = 0; i <= 100; ++i) {
        const Scalar r = 0.01 * i;
        const Scalar value = heaviside_project(r, beta, eta);
        REQUIRE(value >= previous);
        previous = value;
        const Scalar h = 1.0e-6;
        const Scalar fd =
            (heaviside_project(r + h, beta, eta) - heaviside_project(r - h, beta, eta)) /
            (2.0 * h);
        REQUIRE(heaviside_derivative(r, beta, eta) == Approx(fd).epsilon(1.0e-6).margin(1e-9));
      }
    }
  }
  // Small beta is nearly the identity; large beta nearly a step at eta.
  for (int i = 0; i <= 10; ++i) {
    const Scalar r = 0.1 * i;
    REQUIRE(std::abs(heaviside_project(r, 1.0e-3, 0.5) - r) < 1.0e-6);
  }
  REQUIRE(heaviside_project(0.45, 256.0, 0.5) < 1.0e-5);
  REQUIRE(heaviside_project(0.55, 256.0, 0.5) > 1.0 - 1.0e-5);
  const Vector v = heaviside_project(Vector::LinSpaced(5, 0.0, 1.0), 4.0, 0.5);
  REQUIRE(v(2) == Approx(0.5));
}

TEST_CASE("projection options validate and schedule beta", "[projection][topopt][config]") {
  ProjectionOptions o;
  o.enabled = true;
  REQUIRE_NOTHROW(o.validate());
  REQUIRE(o.num_stages() == 6);  // 1, 2, 4, 8, 16, 32
  REQUIRE(o.beta_of_stage(0) == 1.0);
  REQUIRE(o.beta_of_stage(3) == 8.0);
  REQUIRE(o.beta_of_stage(9) == 32.0);
  o.beta_max = 20.0;
  REQUIRE(o.num_stages() == 6);
  REQUIRE(o.beta_of_stage(5) == 20.0);
  ProjectionOptions bad = o;
  bad.eta = 1.0;
  REQUIRE_THROWS_WITH(bad.validate(), ContainsSubstring("eta"));
  bad = o;
  bad.beta_max = 0.5;
  REQUIRE_THROWS_WITH(bad.validate(), ContainsSubstring("beta_max"));
  bad = o;
  bad.beta_factor = 1.0;
  REQUIRE_THROWS_WITH(bad.validate(), ContainsSubstring("beta_factor"));
  bad = o;
  bad.beta_max = 5000.0;
  REQUIRE_THROWS_AS(bad.validate(), ConfigError);
  bad.enabled = false;
  REQUIRE_NOTHROW(bad.validate());  // ignored when off
}

TEST_CASE("compliance and volume gradients stay exact through the projection",
          "[projection][topopt][sensitivity][verification]") {
  const std::vector<Mesh> meshes = {
      make_perturbed_quad_mesh(box(10, 5, 1, 0.5, 0.25, 1.0), 0.2, 2u),
      make_structured_hex_mesh(box(6, 3, 2, 0.6, 0.3, 0.2)),
      make_structured_tet_mesh(box(4, 2, 2, 0.4, 0.2, 0.2))};
  for (const Mesh& mesh : meshes) {
    INFO(to_string(mesh.element_type()));
    FemModel model = make_block(mesh);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.6 * model.mesh().mean_element_size());
    PassiveRegionSpec pad;
    Selector corner;
    corner.kind = SelectorKind::Box;
    corner.xmax = 0.1;
    corner.ymax = 0.1;
    pad.region.members.push_back(corner);
    pad.solid = true;
    DesignDomain domain(model, 0.5, 0.5, {pad});
    StaticAnalysisOptions options;
    options.linear.type = LinearSolverType::SimplicialLdlt;
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(), options);
    for (Scalar beta : {1.0, 8.0}) {
      objective.set_projection(beta, 0.45);
      const Vector x = wavy_design(model, domain);
      const SensitivityCheckResult check =
          verify_sensitivities(objective, domain, x, {}, 1.0e-5, 1.0e-5);
      INFO("beta " << beta << ": compliance max relative error " << check.max_relative_error
                   << ", directional " << check.directional_relative_error);
      REQUIRE(check.passed);
      REQUIRE(check.directional_relative_error < 1.0e-6);

      const ObjectiveEvaluation eval = objective.evaluate(x, true);
      REQUIRE(eval.projection_derivative.size() == domain.num_elements());
      const Scalar volume_error =
          gradient_error(eval.dv_dx, domain, x, 1.0e-6, [&](const Vector& xx) {
            return domain.volume_of(objective.physical_density(xx));
          });
      INFO("volume gradient max relative error " << volume_error);
      REQUIRE(volume_error < 1.0e-6);
      // The physical density is the projection of the filtered one.
      REQUIRE((eval.physical_density -
               heaviside_project(eval.filtered_density, beta, 0.45))
                  .cwiseAbs()
                  .maxCoeff() == 0.0);
    }
  }
}

TEST_CASE("stress constraint gradients stay exact through the projection",
          "[projection][stress][sensitivity][verification]") {
  for (int dim : {2, 3}) {
    INFO("dim " << dim);
    const Mesh mesh = dim == 2 ? make_structured_quad_mesh(box(10, 5, 1, 0.5, 0.25, 1.0))
                               : make_structured_hex_mesh(box(6, 3, 2, 0.6, 0.3, 0.2));
    FemModel model = make_block(mesh);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    StaticAnalysisOptions options;
    options.linear.type = LinearSolverType::SimplicialLdlt;
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(), options);
    objective.set_projection(6.0, 0.5);
    StressConstraintOptions so;
    so.enabled = true;
    so.limit = 3.0e7;
    so.p_norm = 6.0;
    const StressConstraint stress(model, assembler, filter, so);
    const Vector x = wavy_design(model, domain);
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    const StressEvaluation se = stress.evaluate(objective, eval, 0, 0.8, true);
    const Scalar error = gradient_error(se.dg_dx, domain, x, 1.0e-5, [&](const Vector& xx) {
      const ObjectiveEvaluation e = objective.evaluate(xx, false);
      return stress.evaluate(objective, e, 0, 0.8, false).constraint;
    });
    INFO("stress gradient max relative error " << error);
    REQUIRE(error < 1.0e-5);
  }
}

TEST_CASE("the projection refuses the sensitivity filter", "[projection][config]") {
  FemModel model = make_block(make_structured_quad_mesh(box(6, 3, 1, 0.6, 0.3, 1.0)));
  Assembler assembler(model);
  const DensityFilter heuristic(model.mesh(), FilterType::Sensitivity, 0.15);
  DesignDomain domain(model, 0.5, 0.5, {});
  ComplianceObjective objective(model, assembler, heuristic, domain, SimpOptions(),
                                StaticAnalysisOptions());
  REQUIRE_THROWS_WITH(objective.set_projection(4.0, 0.5), ContainsSubstring("density filter"));
  TopologyOptimizerOptions options;
  options.projection.enabled = true;
  REQUIRE_THROWS_AS(TopologyOptimizer(model, assembler, heuristic, domain, options),
                    ConfigError);
}

TEST_CASE("projected optimisation runs end crisper at the final beta",
          "[projection][topopt]") {
  // A small MBB half-beam: support at the left edge (symmetry) and a roller
  // at the bottom right, load at the top left.
  StructuredMeshSpec spec = box(60, 20, 1, 3.0, 1.0, 1.0);
  FemModel model(make_structured_quad_mesh(spec), default_material(), 1.0,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint symmetry;
  Selector left;
  left.kind = SelectorKind::Box;
  left.xmax = 0.0;
  symmetry.region.members.push_back(left);
  symmetry.fix_x = true;
  model.constraints().push_back(symmetry);
  DisplacementConstraint roller;
  Selector corner;
  corner.kind = SelectorKind::NearestNode;
  corner.point = Vector3(3.0, 0.0, 0.0);
  roller.region.members.push_back(corner);
  roller.fix_y = true;
  model.constraints().push_back(roller);
  LoadCaseSpec load;
  load.name = "top";
  PointLoadSpec p;
  Selector top;
  top.kind = SelectorKind::NearestNode;
  top.point = Vector3(0.0, 1.0, 0.0);
  p.region.members.push_back(top);
  p.force = Vector3(0.0, -1.0, 0.0);
  load.point_loads.push_back(p);
  model.load_case_specs().push_back(load);
  model.finalize();
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density, 0.1 * 1.5 * 1.0);
  DesignDomain domain(model, 0.5, 0.5, {});

  for (OptimizerMethod method : {OptimizerMethod::OptimalityCriteria, OptimizerMethod::MMA}) {
    INFO("method " << to_string(method));
    TopologyOptimizerOptions options;
    options.method = method;
    options.max_iterations = 160;
    options.change_tolerance = 0.01;
    options.history_stride = 0;
    TopologyOptimizer plain(model, assembler, filter, domain, options);
    const TopologyOptimizationResult base = plain.run();

    options.projection.enabled = true;
    options.projection.beta_start = 1.0;
    options.projection.beta_max = 8.0;
    options.projection.beta_interval = 25;
    TopologyOptimizer projected(model, assembler, filter, domain, options);
    const TopologyOptimizationResult sharp = projected.run();
    INFO("grey level " << base.gray_level << " -> " << sharp.gray_level << ", iterations "
                       << base.iterations << " / " << sharp.iterations);
    REQUIRE(sharp.projected);
    REQUIRE(sharp.final_beta == 8.0);
    REQUIRE(sharp.filtered_density.size() == model.mesh().num_elements());
    REQUIRE(sharp.gray_level < 0.5 * base.gray_level);
    REQUIRE(sharp.volume_fraction == Approx(0.5).epsilon(1.0e-3));
    // Beta is non-decreasing through the history and ends at beta_max.
    Scalar last = 0.0;
    for (const TopologyIteration& it : sharp.history) {
      REQUIRE(it.beta >= last);
      last = it.beta;
    }
    REQUIRE(last == 8.0);
    if (sharp.converged) REQUIRE(sharp.history.back().beta == 8.0);
  }
}

TEST_CASE("a deck switches the projection on with its schedule",
          "[projection][io][config]") {
  const json::Value doc = json::parse(R"({
    "mesh": { "type": "structured_quad", "nx": 8, "ny": 4, "lx": 1, "ly": 0.5 },
    "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3 },
    "boundary_conditions": [ { "fix": ["x", "y"], "region": { "box": { "xmax": 0 } } } ],
    "load_cases": [ { "name": "c", "point_loads": [
        { "region": { "nearest_node": [1, 0] }, "force": [0, -1] } ] } ],
    "topology": { "enabled": true, "volume_fraction": 0.4,
                  "projection": { "enabled": true, "eta": 0.4, "beta_max": 16,
                                  "beta_interval": 30 } }
  })");
  const Configuration config = parse_configuration(doc, "inline", true);
  const ProjectionOptions& p = config.topology.optimizer.projection;
  REQUIRE(p.enabled);
  REQUIRE(p.eta == 0.4);
  REQUIRE(p.beta_max == 16.0);
  REQUIRE(p.beta_interval == 30);
  REQUIRE(p.beta_start == 1.0);
  const json::Value sensitivity = json::parse(R"({
    "mesh": { "type": "structured_quad", "nx": 8, "ny": 4, "lx": 1, "ly": 0.5 },
    "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3 },
    "boundary_conditions": [ { "fix": ["x", "y"], "region": { "box": { "xmax": 0 } } } ],
    "load_cases": [ { "name": "c", "point_loads": [
        { "region": { "nearest_node": [1, 0] }, "force": [0, -1] } ] } ],
    "topology": { "enabled": true, "filter": { "type": "sensitivity" },
                  "projection": { "enabled": true } }
  })");
  REQUIRE_THROWS_WITH(parse_configuration(sensitivity, "inline", true),
                      ContainsSubstring("density filter"));
}
