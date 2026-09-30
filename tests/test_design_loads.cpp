/// \file test_design_loads.cpp
/// \brief Design-dependent loads in topology optimisation: the body-load
///        interpolation, the load vectors of a design, the compliance, stress
///        and buckling gradients under self-weight, rotation, body forces and
///        temperature against central differences, the bounded parasitic
///        load of near-void material, and what the optimiser refuses.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Loads.hpp"
#include "sparlab/topopt/BucklingConstraint.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/DesignLoads.hpp"
#include "sparlab/topopt/StressConstraint.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <functional>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {

StructuredMeshSpec box_spec(Index nx, Index ny, Index nz, Scalar lx, Scalar ly, Scalar lz) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.lx = lx;
  spec.ly = ly;
  spec.lz = lz;
  return spec;
}

Selector box_selector(Scalar xmin, Scalar xmax) {
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  return s;
}

IsotropicMaterial heated_steel() {
  IsotropicMaterial m(200.0e9, 0.3, 7850.0, "steel");
  m.set_thermal(1.2e-5, 20.0, 50.0);
  return m;
}

/// The load cases of the plate: its weight with a point load, a rotation with
/// a body force density on the right half, a uniform temperature rise with the
/// point load, and two temperature regions.
std::vector<LoadCaseSpec> plate_cases(Scalar length, Scalar point_force, int dim) {
  PointLoadSpec point;
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector3(0.5 * length, 0.0, 0.0);
  point.region.members.push_back(nearest);
  point.force = Vector3(0.0, -point_force, 0.0);

  std::vector<LoadCaseSpec> cases(4);
  cases[0].name = "weight";
  cases[0].gravity = Vector3(0.0, -9.81, 0.0);
  cases[0].point_loads.push_back(point);
  cases[1].name = "spin";
  cases[1].weight = 0.5;
  cases[1].centrifugal.enabled = true;
  cases[1].centrifugal.angular_velocity = 20.0;
  cases[1].centrifugal.axis = Vector3::UnitZ();
  cases[1].centrifugal.point = Vector3(0.0, 0.0, 0.0);
  BodyForceSpec body;
  body.whole_model = false;
  body.region.members.push_back(box_selector(0.5 * length, length));
  body.force_density = Vector3(0.0, -1.0e5, 0.0);
  cases[1].body_forces.push_back(body);
  cases[2].name = "heat";
  cases[2].temperature.source = TemperatureSpec::Source::Uniform;
  cases[2].temperature.uniform = 70.0;
  cases[2].point_loads.push_back(point);
  cases[3].name = "patch";
  cases[3].weight = 2.0;
  cases[3].temperature.source = TemperatureSpec::Source::Regions;
  cases[3].temperature.uniform = 20.0;
  RegionValue hot;
  hot.region.members.push_back(box_selector(0.5 * length, length));
  hot.value = 120.0;
  cases[3].temperature.regions.push_back(hot);
  (void)dim;
  return cases;
}

/// A plate (Q4) or block (Hex8) of length `length` along x, clamped at x = 0
/// and held along x at x = length, so that heating it compresses it; with
/// `hold_end` false, a plate on statically determinate supports (a corner
/// pinned, the corner above it on a roller), free to expand. Not finalised.
FemModel plate_model(Mesh mesh, Scalar length, Scalar thickness, bool hold_end = true) {
  const int dim = mesh.dim();
  const Scalar height = mesh.bounding_box().upper.y();
  FemModel model(std::move(mesh), heated_steel(), thickness,
                 dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                 IntegrationOptions());
  if (!hold_end) {
    DisplacementConstraint pin;
    Selector corner;
    corner.kind = SelectorKind::NearestNode;
    corner.point = Vector3::Zero();
    pin.region.members.push_back(corner);
    pin.fix_x = pin.fix_y = true;
    model.constraints().push_back(pin);
    DisplacementConstraint roller;
    Selector above;
    above.kind = SelectorKind::NearestNode;
    above.point = Vector3(0.0, height, 0.0);
    roller.region.members.push_back(above);
    roller.fix_x = true;
    model.constraints().push_back(roller);
    return model;
  }
  DisplacementConstraint root;
  root.region.members.push_back(box_selector(-1.0, 0.0));
  root.fix_x = root.fix_y = true;
  root.fix_z = dim == 3;
  model.constraints().push_back(root);
  DisplacementConstraint end;
  end.region.members.push_back(box_selector(length, 2.0 * length));
  end.fix_x = true;
  model.constraints().push_back(end);
  return model;
}

/// The plate with all four load cases, finalised.
FemModel make_plate(Mesh mesh, Scalar length, Scalar thickness, Scalar point_force) {
  const int dim = mesh.dim();
  FemModel model = plate_model(std::move(mesh), length, thickness);
  model.load_case_specs() = plate_cases(length, point_force, dim);
  model.finalize();
  return model;
}

/// The plate with the uniform temperature rise alone (no point load), held
/// at its end or free to expand.
FemModel make_heated_plate(Mesh mesh, Scalar length, Scalar thickness, bool hold_end) {
  const int dim = mesh.dim();
  FemModel model = plate_model(std::move(mesh), length, thickness, hold_end);
  LoadCaseSpec heat = plate_cases(length, 0.0, dim)[2];
  heat.point_loads.clear();
  model.load_case_specs() = {heat};
  model.finalize();
  return model;
}

/// A smooth design in [0.02, 0.98] with a low-density region, so that some
/// filtered densities fall below the body-load threshold.
Vector graded_design(const FemModel& model, const DesignDomain& domain) {
  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    const Vector3 c = model.mesh().element_centroid(e);
    const Scalar q = 0.5 + 0.5 * std::sin(11.0 * c.x()) * std::cos(7.0 * c.y() + 5.0 * c.z());
    x(e) = 0.02 + 0.96 * q * q;
  }
  domain.clamp(x);
  return x;
}

/// Max relative error of `analytical` against central differences of `value`
/// - of second order, and of fourth, (-v(x+2h) + 8 v(x+h) - 8 v(x-h) +
/// v(x-2h)) / 12h, whose larger steps an eigenvalue's round-off needs - at
/// the best of four steps, over every `stride`-th variable, with entries
/// below 1e-3 of the gradient's scale judged against that scale. Fourth order
/// goes first, largest step first: an eigenvalue's round-off favours the
/// larger reach. The sweep stops at the first order and step whose error is
/// at or below `good_enough`: one step showing the agreement settles a check
/// against a tolerance above it.
Scalar best_gradient_error(const Vector& analytical, const Vector& x,
                           const std::function<Scalar(const Vector&)>& value,
                           Eigen::Index stride = 1, Scalar good_enough = 0.0) {
  const Scalar floor = 1.0e-3 * analytical.cwiseAbs().maxCoeff();
  Scalar best = std::numeric_limits<Scalar>::infinity();
  Vector xs = x;
  const auto at = [&](Eigen::Index e, Scalar offset) {
    xs(e) = x(e) + offset;
    const Scalar v = value(xs);
    xs(e) = x(e);
    return v;
  };
  for (const int order : {4, 2}) {
    for (const Scalar step : {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6}) {
      const Scalar reach = order == 4 ? 2.0 * step : step;
      Scalar worst = 0.0;
      for (Eigen::Index e = 0; e < x.size(); e += stride) {
        if (x(e) - reach < 0.0 || x(e) + reach > 1.0) continue;
        const Scalar fd =
            order == 4 ? (-at(e, 2.0 * step) + 8.0 * at(e, step) - 8.0 * at(e, -step) +
                          at(e, -2.0 * step)) / (12.0 * step)
                       : (at(e, step) - at(e, -step)) / (2.0 * step);
        worst = std::max(worst, std::abs(fd - analytical(e)) /
                                    std::max({std::abs(fd), std::abs(analytical(e)), floor}));
      }
      best = std::min(best, worst);
      if (best <= good_enough) return best;
    }
  }
  return best;
}

StaticAnalysisOptions direct_solver() {
  StaticAnalysisOptions options;
  options.linear.type = LinearSolverType::SimplicialLdlt;
  return options;
}

}  // namespace

TEST_CASE("the body-load interpolation is the density above the threshold and bounded "
          "against the stiffness below it",
          "[topopt][loads][simp]") {
  for (const Scalar p : {1.0, 2.0, 3.0, 4.5}) {
    SimpOptions simp;
    simp.penalty = p;
    simp.body_load_threshold = 0.1;
    INFO("p = " << p);
    // The physical load at and above the threshold, zero at zero.
    for (const Scalar rho : {0.1, 0.25, 0.5, 1.0}) {
      REQUIRE(body_load_factor(rho, simp) == rho);
      REQUIRE(body_load_derivative(rho, simp) == 1.0);
    }
    REQUIRE(body_load_factor(0.0, simp) == 0.0);
    // Continuous with a continuous slope at the threshold.
    const Scalar below = 0.1 * (1.0 - 1.0e-12);
    REQUIRE(body_load_factor(below, simp) == Approx(0.1).epsilon(1.0e-11));
    REQUIRE(body_load_derivative(below, simp) == Approx(1.0).epsilon(1.0e-10));
    // Increasing, its derivative the central difference of the factor, and
    // the load never more than p rho_t^(1-p) times the stiffness's share.
    Scalar previous = -1.0;
    bool increasing = true;
    Scalar derivative_error = 0.0;
    Scalar ratio = 0.0;
    for (int k = 1; k < 1000; ++k) {
      const Scalar rho = 0.001 * k;
      const Scalar gamma = body_load_factor(rho, simp);
      increasing = increasing && gamma > previous;
      previous = gamma;
      const Scalar h = 1.0e-7;
      if (std::abs(rho - 0.1) > 2.0 * h) {
        const Scalar fd =
            (body_load_factor(rho + h, simp) - body_load_factor(rho - h, simp)) / (2.0 * h);
        derivative_error = std::max(derivative_error,
                                    std::abs(body_load_derivative(rho, simp) - fd) /
                                        std::max(std::abs(fd), 1.0e-3));
      }
      ratio = std::max(ratio, gamma / simp_stiffness_factor(rho, simp));
    }
    REQUIRE(increasing);
    REQUIRE(derivative_error < 1.0e-6);
    REQUIRE(ratio <= p * std::pow(0.1, 1.0 - p) / (1.0 - simp.emin_ratio) * (1.0 + 1.0e-12));
    // At p = 1 the load is the density everywhere.
    if (p == 1.0) {
      REQUIRE(body_load_factor(0.03, simp) == Approx(0.03).epsilon(1.0e-14));
    }
  }
  // A zero threshold keeps gamma = rho; outside [0, 1) it is refused.
  SimpOptions plain;
  plain.body_load_threshold = 0.0;
  REQUIRE(body_load_factor(1.0e-3, plain) == 1.0e-3);
  plain.body_load_threshold = 1.0;
  REQUIRE_THROWS_WITH(plain.validate(), ContainsSubstring("body_load_threshold"));
  plain.body_load_threshold = -0.1;
  REQUIRE_THROWS_AS(plain.validate(), ConfigError);
}

TEST_CASE("the design's load vectors: the model's own at full density, scaled parts "
          "at a uniform one",
          "[topopt][loads]") {
  FemModel model = make_plate(make_structured_quad_mesh(box_spec(12, 6, 1, 0.6, 0.3, 1.0)),
                              0.6, 0.01, 100.0);
  const DesignLoads loads(model);
  REQUIRE(loads.design_dependent());
  REQUIRE(loads.num_cases() == 4);
  REQUIRE(loads.has_body(0));
  REQUIRE_FALSE(loads.has_thermal(0));
  REQUIRE(loads.has_body(1));
  REQUIRE(loads.has_thermal(2));
  REQUIRE(loads.has_thermal(3));
  REQUIRE(loads.temperature(3).size() == model.mesh().num_nodes());

  const Index ne = model.mesh().num_elements();
  SimpOptions simp;
  for (std::size_t l = 0; l < 4; ++l) {
    INFO("case " << l);
    const Vector full = loads.load(l, Vector::Ones(ne), simp);
    REQUIRE((full - model.load_vectors()[l]).cwiseAbs().maxCoeff() <=
            1.0e-15 * model.load_vectors()[l].cwiseAbs().maxCoeff());
    // At a uniform density the body part scales with gamma and the thermal
    // part with E(rho) / E0; the mechanical part stays.
    for (const Scalar rho : {0.5, 0.05}) {
      const LoadCaseData& data = model.load_case_data(l);
      Vector expected = data.mechanical;
      if (data.body.size() > 0) expected += body_load_factor(rho, simp) * data.body;
      if (data.thermal.size() > 0) expected += simp_stiffness_factor(rho, simp) * data.thermal;
      const Vector f = loads.load(l, Vector::Constant(ne, rho), simp);
      REQUIRE((f - expected).cwiseAbs().maxCoeff() <=
              1.0e-14 * model.load_vectors()[l].cwiseAbs().maxCoeff());
    }
  }
  // The element loads sum to the assembled ones.
  const std::vector<Vector> body = element_body_loads(model, model.load_case_specs()[1]);
  Vector sum = Vector::Zero(model.dofs().num_dofs());
  for (Index e = 0; e < ne; ++e) {
    if (body[static_cast<std::size_t>(e)].size() > 0) {
      model.dofs().scatter_add(model.mesh().element_nodes(e), 4, body[static_cast<std::size_t>(e)],
                               sum);
    }
  }
  REQUIRE((sum - model.load_case_data(1).body).cwiseAbs().maxCoeff() == 0.0);
  const ElementThermalLoads thermal =
      element_thermal_loads(model, model.load_case_data(3).temperature);
  REQUIRE(thermal.expanding == ne);
  REQUIRE(thermal.total_self_energy == model.load_case_data(3).thermal_self_energy);
  REQUIRE(thermal.self_energy.sum() ==
          Approx(model.load_case_data(3).thermal_self_energy).epsilon(1.0e-13));
}

TEST_CASE("compliance gradients under design-dependent loads match central differences",
          "[topopt][loads][sensitivity][verification]") {
  struct Case {
    std::string name;
    Mesh mesh;
    Scalar thickness;
    Scalar force;
    bool projection;
  };
  std::vector<Case> cases;
  cases.push_back({"Q4", make_structured_quad_mesh(box_spec(12, 6, 1, 0.6, 0.3, 1.0)), 0.01,
                   100.0, false});
  cases.push_back({"Q4 projected", make_structured_quad_mesh(box_spec(12, 6, 1, 0.6, 0.3, 1.0)),
                   0.01, 100.0, true});
  cases.push_back({"Hex8", make_structured_hex_mesh(box_spec(6, 3, 2, 0.6, 0.3, 0.1)), 1.0,
                   1000.0, false});
  for (Case& c : cases) {
    INFO(c.name);
    FemModel model = make_plate(std::move(c.mesh), 0.6, c.thickness, c.force);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                  direct_solver());
    if (c.projection) objective.set_projection(4.0, 0.5);
    const Vector x = graded_design(model, domain);
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    // Some densities lie below the body-load threshold, and the load term
    // makes some entries of the gradient positive.
    REQUIRE(eval.physical_density.minCoeff() < 0.1);
    REQUIRE(eval.dc_dx.maxCoeff() > 0.0);
    // The compliance is the work of this design's loads.
    Scalar c_sum = 0.0;
    const std::vector<Scalar> w = model.normalised_weights();
    for (std::size_t l = 0; l < eval.loads.size(); ++l) {
      REQUIRE(eval.load_case_compliance[l] ==
              Approx(eval.loads[l].dot(eval.displacements[l])).epsilon(1.0e-14));
      c_sum += w[l] * eval.load_case_compliance[l];
    }
    REQUIRE(eval.compliance == Approx(c_sum).epsilon(1.0e-14));
    const Scalar error = best_gradient_error(eval.dc_dx, x, [&](const Vector& xx) {
      return objective.compliance_at(xx);
    });
    INFO("compliance gradient error " << error);
    REQUIRE(error < 1.0e-6);
  }
}

TEST_CASE("the thermoelastic element energy is the elastic strain energy of the state",
          "[topopt][loads]") {
  // One load case: a uniform temperature rise of a plate held at both ends.
  FemModel model =
      make_heated_plate(make_structured_quad_mesh(box_spec(8, 4, 1, 0.6, 0.3, 1.0)), 0.6, 0.01,
                        /*hold_end=*/true);
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::None, 0.0);
  DesignDomain domain(model, 0.5, 0.5, {});
  ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                direct_solver());
  // At full density: 1/2 u^T K u - u^T f_th + 1/2 int eps0 D eps0.
  const ObjectiveEvaluation eval =
      objective.evaluate(Vector::Ones(model.mesh().num_elements()), false);
  const Vector& u = eval.displacements[0];
  const SparseMatrix k = assembler.assemble_stiffness();
  const LoadCaseData& data = model.load_case_data(0);
  const Scalar exact = 0.5 * u.dot(k * u) - u.dot(data.thermal) + data.thermal_self_energy;
  REQUIRE(eval.element_strain_energy.sum() == Approx(exact).epsilon(1.0e-11));
  REQUIRE(exact > 0.0);
  // Free expansion - the end released - stores no energy at any density.
  FemModel free_model =
      make_heated_plate(make_structured_quad_mesh(box_spec(8, 4, 1, 0.6, 0.3, 1.0)), 0.6, 0.01,
                        /*hold_end=*/false);
  Assembler free_assembler(free_model);
  DesignDomain free_domain(free_model, 0.5, 0.5, {});
  ComplianceObjective free_objective(free_model, free_assembler, filter, free_domain,
                                     SimpOptions(), direct_solver());
  const ObjectiveEvaluation free_eval =
      free_objective.evaluate(Vector::Constant(free_model.mesh().num_elements(), 0.7), false);
  REQUIRE(std::abs(free_eval.element_strain_energy.sum()) <=
          1.0e-9 * data.thermal_self_energy);
}

TEST_CASE("stress-constraint gradients under design-dependent loads match central "
          "differences",
          "[topopt][loads][stress][sensitivity][verification]") {
  FemModel model = make_plate(make_structured_quad_mesh(box_spec(10, 5, 1, 0.6, 0.3, 1.0)), 0.6,
                              0.01, 100.0);
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density,
                             1.5 * model.mesh().mean_element_size());
  DesignDomain domain(model, 0.5, 0.5, {});
  ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(), direct_solver());
  StressConstraintOptions so;
  so.enabled = true;
  so.limit = 5.0e7;
  so.p_norm = 6.0;
  so.relaxation = 0.5;
  const StressConstraint stress(model, assembler, filter, so);
  const Vector x = graded_design(model, domain);
  for (std::size_t l = 0; l < 4; ++l) {
    INFO("load case " << model.load_case_specs()[l].name);
    // The adjoint solve uses the factorisation of the last evaluation, which
    // the central differences of the previous case replaced.
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    const StressEvaluation se = stress.evaluate(objective, eval, l, 0.8, true);
    const Scalar error = best_gradient_error(se.dg_dx, x, [&](const Vector& xx) {
      const ObjectiveEvaluation e = objective.evaluate(xx, false);
      return stress.evaluate(objective, e, l, 0.8, false).constraint;
    });
    INFO("stress gradient error " << error);
    REQUIRE(error < 1.0e-6);
  }
  // With a temperature field the stress is that of the strain less the free
  // thermal strain: a uniformly heated plate free to expand is stress-free.
  FemModel free_model =
      make_heated_plate(make_structured_quad_mesh(box_spec(6, 3, 1, 0.6, 0.3, 1.0)), 0.6, 0.01,
                        /*hold_end=*/false);
  Assembler free_assembler(free_model);
  const DensityFilter none(free_model.mesh(), FilterType::None, 0.0);
  DesignDomain free_domain(free_model, 0.5, 0.5, {});
  ComplianceObjective free_objective(free_model, free_assembler, none, free_domain,
                                     SimpOptions(), direct_solver());
  const StressConstraint free_stress(free_model, free_assembler, none, so);
  const ObjectiveEvaluation free_eval =
      free_objective.evaluate(Vector::Ones(free_model.mesh().num_elements()), false);
  const StressEvaluation fs = free_stress.evaluate(free_objective, free_eval, 0, 1.0, false);
  REQUIRE(fs.solid_von_mises.maxCoeff() <= 1.0e-6 * 200.0e9 * 1.2e-5 * 50.0);
}

TEST_CASE("buckling-constraint gradients under self-weight and heating match central "
          "differences",
          "[topopt][loads][buckling][sensitivity][verification]") {
  // A column along x clamped at x = 0: (a) compressed by an end traction and
  // its own weight along -x; (b) held along x at the far end and heated.
  for (const bool heated : {false, true}) {
    INFO((heated ? "heated" : "weight"));
    FemModel model(make_structured_quad_mesh(box_spec(20, 4, 1, 0.8, 0.08, 1.0)),
                   heated_steel(), 0.01, StressState::PlaneStress, IntegrationOptions());
    DisplacementConstraint root;
    root.region.members.push_back(box_selector(-1.0, 0.0));
    root.fix_x = root.fix_y = true;
    model.constraints().push_back(root);
    LoadCaseSpec load;
    load.name = heated ? "heat" : "weight";
    if (heated) {
      DisplacementConstraint end;
      end.region.members.push_back(box_selector(0.8, 1.6));
      end.fix_x = true;
      model.constraints().push_back(end);
      load.temperature.source = TemperatureSpec::Source::Uniform;
      load.temperature.uniform = 30.0;
    } else {
      load.gravity = Vector3(-9.81e3, 0.0, 0.0);
      TractionLoadSpec top;
      top.region.members.push_back(box_selector(0.8, 1.6));
      top.traction = Vector3(-2.0e6, 0.0, 0.0);
      load.tractions.push_back(top);
    }
    model.load_case_specs().push_back(load);
    model.finalize();
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                  direct_solver());
    BucklingConstraintOptions bo;
    bo.enabled = true;
    bo.min_load_factor = 5.0;
    bo.num_modes = 3;
    bo.ks_parameter = 20.0;
    bo.eigen.tolerance = 1.0e-14;
    bo.eigen.residual_tolerance = 1.0e-11;
    bo.eigen.max_iterations = 2000;
    BucklingConstraint buckling(model, assembler, bo);
    Vector x = domain.initial_design();
    for (Index e = 0; e < domain.num_elements(); ++e) {
      const Vector3 c = model.mesh().element_centroid(e);
      x(e) = 0.55 + 0.3 * std::sin(9.0 * c.x()) * std::cos(40.0 * c.y());
    }
    domain.clamp(x);
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    const BucklingEvaluation be = buckling.evaluate(objective, eval, 0, true);
    REQUIRE(be.load_factors.size() == 3);
    REQUIRE(be.load_factors(0) > 0.0);
    REQUIRE(be.load_factors(1) > 1.01 * be.load_factors(0));
    INFO("lambda = " << be.load_factors.transpose());
    // Every eigenvalue costs a solve; the sweep stops once a step agrees
    // within half the tolerance.
    const Eigen::Index stride = std::max<Eigen::Index>(1, x.size() / 20);
    const Vector dlambda_dx = objective.chain_to_design(eval, be.dlambda_dphysical[0]);
    const Scalar error = best_gradient_error(
        dlambda_dx, x,
        [&](const Vector& xx) {
          const ObjectiveEvaluation e = objective.evaluate(xx, false);
          return buckling.evaluate(objective, e, 0, false).load_factors(0);
        },
        stride, 5.0e-6);
    INFO("lambda_1 gradient error " << error);
    REQUIRE(error < 1.0e-5);
  }
}

TEST_CASE("near-void material under self-weight: the threshold bounds its parasitic load",
          "[topopt][loads][simp]") {
  // A cantilever whose right half is near void under its own weight and a
  // tip load at the end of the solid half.
  StructuredMeshSpec spec = box_spec(16, 4, 1, 1.6, 0.2, 1.0);
  FemModel model(make_structured_quad_mesh(spec), heated_steel(), 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint root;
  root.region.members.push_back(box_selector(-1.0, 0.0));
  root.fix_x = root.fix_y = true;
  model.constraints().push_back(root);
  LoadCaseSpec weight;
  weight.name = "weight";
  weight.gravity = Vector3(0.0, -9.81, 0.0);
  PointLoadSpec tip;
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector3(0.8, 0.0, 0.0);
  tip.region.members.push_back(nearest);
  tip.force = Vector3(0.0, -200.0, 0.0);
  weight.point_loads.push_back(tip);
  model.load_case_specs().push_back(weight);
  model.finalize();
  Assembler assembler(model);
  const DensityFilter none(model.mesh(), FilterType::None, 0.0);
  DesignDomain domain(model, 0.5, 0.5, {});

  const auto compliance = [&](Scalar threshold, Scalar void_density) {
    SimpOptions simp;
    simp.body_load_threshold = threshold;
    ComplianceObjective objective(model, assembler, none, domain, simp, direct_solver());
    Vector x = Vector::Ones(model.mesh().num_elements());
    for (Index e = 0; e < x.size(); ++e) {
      if (model.mesh().element_centroid(e).x() > 0.8) x(e) = void_density;
    }
    return objective.compliance_at(x);
  };
  // The solid half alone: the same model with the right half weightless and
  // (nearly) stiffness-free.
  const Scalar solid_half = compliance(0.1, 0.0);
  Scalar previous = solid_half;
  for (const Scalar rho_v : {1.0e-2, 3.0e-3, 1.0e-3}) {
    const Scalar with_threshold = compliance(0.1, rho_v);
    const Scalar without = compliance(0.0, rho_v);
    INFO("void density " << rho_v << ": " << with_threshold << " J with the threshold, "
                         << without << " J without, " << solid_half << " J the solid half");
    // Bounded: within a few per cent of the solid half.
    REQUIRE(std::abs(with_threshold / solid_half - 1.0) < 0.02);
    // Unbounded without it: the near-void half sags under its own weight,
    // more the emptier it is.
    REQUIRE(without > 2.0 * solid_half);
    REQUIRE(without > previous);
    previous = without;
  }
  REQUIRE(previous > 10.0 * solid_half);
}

TEST_CASE("what the optimiser refuses with design-dependent loads",
          "[topopt][loads][diagnostics]") {
  FemModel model = make_plate(make_structured_quad_mesh(box_spec(8, 4, 1, 0.6, 0.3, 1.0)), 0.6,
                              0.01, 100.0);
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density, 0.1);
  DesignDomain domain(model, 0.5, 0.5, {});
  TopologyOptimizerOptions options;
  options.method = OptimizerMethod::OptimalityCriteria;
  // The OC update cannot follow a positive gradient.
  REQUIRE_THROWS_WITH(TopologyOptimizer(model, assembler, filter, domain, options),
                      ContainsSubstring("optimality-criteria") &&
                          ContainsSubstring("mma"));
  options.method = OptimizerMethod::MMA;
  REQUIRE_NOTHROW(TopologyOptimizer(model, assembler, filter, domain, options));
  // The heuristic sensitivity filter assumes a non-positive gradient.
  const DensityFilter heuristic(model.mesh(), FilterType::Sensitivity, 0.1);
  REQUIRE_THROWS_WITH(TopologyOptimizer(model, assembler, heuristic, domain, options),
                      ContainsSubstring("sensitivity filter"));

  // A conducted temperature depends on the design.
  FemModel conducted =
      plate_model(make_structured_quad_mesh(box_spec(8, 4, 1, 0.6, 0.3, 1.0)), 0.6, 0.01);
  LoadCaseSpec heat;
  heat.name = "conducted";
  heat.temperature.source = TemperatureSpec::Source::Conduction;
  RegionValue cold;
  cold.region.members.push_back(box_selector(-1.0, 0.0));
  cold.value = 20.0;
  RegionValue hot;
  hot.region.members.push_back(box_selector(0.6, 1.2));
  hot.value = 80.0;
  heat.temperature.conduction.prescribed = {cold, hot};
  conducted.load_case_specs() = {heat};
  conducted.finalize();
  Assembler conducted_assembler(conducted);
  DesignDomain conducted_domain(conducted, 0.5, 0.5, {});
  REQUIRE_THROWS_WITH(
      TopologyOptimizer(conducted, conducted_assembler, filter, conducted_domain, options),
      ContainsSubstring("conducted temperature field"));

  // Non-zero prescribed displacements break the self-adjoint compliance.
  FemModel pushed =
      plate_model(make_structured_quad_mesh(box_spec(8, 4, 1, 0.6, 0.3, 1.0)), 0.6, 0.01);
  pushed.constraints().back().value_x = 1.0e-4;
  pushed.load_case_specs() = plate_cases(0.6, 100.0, 2);
  pushed.finalize();
  Assembler pushed_assembler(pushed);
  DesignDomain pushed_domain(pushed, 0.5, 0.5, {});
  REQUIRE_THROWS_WITH(ComplianceObjective(pushed, pushed_assembler, filter, pushed_domain,
                                          SimpOptions(), direct_solver()),
                      ContainsSubstring("homogeneous supports"));
}

TEST_CASE("MMA designs a plate under its own weight and a point load", "[topopt][loads][mma]") {
  StructuredMeshSpec spec = box_spec(24, 12, 1, 1.2, 0.6, 1.0);
  FemModel model(make_structured_quad_mesh(spec), heated_steel(), 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint left;
  Selector corner_left;
  corner_left.kind = SelectorKind::NearestNode;
  corner_left.point = Vector3(0.0, 0.0, 0.0);
  left.region.members.push_back(corner_left);
  left.fix_x = left.fix_y = true;
  model.constraints().push_back(left);
  DisplacementConstraint right;
  Selector corner_right;
  corner_right.kind = SelectorKind::NearestNode;
  corner_right.point = Vector3(1.2, 0.0, 0.0);
  right.region.members.push_back(corner_right);
  right.fix_y = true;
  model.constraints().push_back(right);
  LoadCaseSpec load;
  load.name = "weight_and_load";
  load.gravity = Vector3(0.0, -9.81, 0.0);
  PointLoadSpec mid;
  Selector centre;
  centre.kind = SelectorKind::NearestNode;
  centre.point = Vector3(0.6, 0.6, 0.0);
  mid.region.members.push_back(centre);
  mid.force = Vector3(0.0, -300.0, 0.0);
  load.point_loads.push_back(mid);
  model.load_case_specs().push_back(load);
  model.finalize();
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density, 1.5 * 0.05);
  DesignDomain domain(model, 0.4, 0.4, {});
  TopologyOptimizerOptions options;
  options.method = OptimizerMethod::MMA;
  options.max_iterations = 40;
  options.analysis = direct_solver();
  options.history_stride = 0;
  TopologyOptimizer optimizer(model, assembler, filter, domain, options);
  const TopologyOptimizationResult result = optimizer.run();
  REQUIRE(result.history.size() >= 2);
  // The compliance falls and the final one is the work of the final design's
  // own loads: its weight is that of its densities.
  REQUIRE(result.compliance < 0.8 * result.history.front().compliance);
  REQUIRE(result.volume_fraction <= 0.4 * (1.0 + options.constraint_tolerance) + 1.0e-12);
  const DesignLoads loads(model);
  const Vector f = loads.load(0, result.physical_density, options.simp);
  REQUIRE(result.compliance == Approx(f.dot(result.displacements[0])).epsilon(1.0e-12));
  // The design's weight is the weight of its material above the threshold.
  Scalar weight = 0.0;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    weight += f(2 * n + 1) - model.load_case_data(0).mechanical(2 * n + 1);
  }
  Scalar expected = 0.0;
  const Vector volumes = domain.element_volumes();
  for (Index e = 0; e < volumes.size(); ++e) {
    expected += body_load_factor(result.physical_density(e), options.simp) * volumes(e);
  }
  expected *= -9.81 * 7850.0;
  REQUIRE(weight == Approx(expected).epsilon(1.0e-12));
}

