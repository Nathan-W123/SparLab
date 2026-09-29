/// \file verify_topopt.cpp
/// \brief Verification of the topology optimiser's design-dependent loads.
///
/// Study `design-loads`:
///   * the load vector of a design against its definition - at full density
///     the model's own, at a uniform one its mechanical part plus gamma(rho)
///     times the body part and E(rho)/E0 times the thermal part;
///   * the compliance gradient under self-weight, rotation with a body force,
///     a uniform temperature and two temperature regions (Q4, with and without
///     the projection; Hex8), the aggregated stress constraint's gradient
///     under the same four cases (Q4), and the lowest buckling load factor's
///     gradient under self-weight and under heating (Q4) - every one against
///     central differences at a design whose filtered densities reach below
///     the body-load threshold;
///   * the parasitic load of near-void material: a cantilever whose outer
///     half is near void, under its own weight, with the body-load threshold
///     and without it.
#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/fem/Assembler.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"
#include "sparlab/topopt/BucklingConstraint.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/DesignLoads.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/StressConstraint.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

std::string fmt(Scalar v, int digits = 4) { return app::format(v, digits); }

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

Selector x_range(Scalar xmin, Scalar xmax) {
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  return s;
}

Selector nearest(const Vector3& point) {
  Selector s;
  s.kind = SelectorKind::NearestNode;
  s.point = point;
  return s;
}

/// Steel with thermal expansion 1.2e-5 1/K from 20 degrees.
IsotropicMaterial steel() {
  IsotropicMaterial m(200.0e9, 0.3, 7850.0, "steel");
  m.set_thermal(1.2e-5, 20.0, 50.0);
  return m;
}

StaticAnalysisOptions direct_solver() {
  StaticAnalysisOptions options;
  options.linear.type = LinearSolverType::SimplicialLdlt;
  return options;
}

/// A plate (Q4) or block (Hex8) 0.6 m long, clamped at x = 0 and held along
/// x at x = 0.6 so that heating compresses it, with four load cases: its
/// weight and a point load; a rotation about z with a body force density on
/// the outer half; a uniform rise of 50 K and the point load; and the outer
/// half at 120 degrees against 20 elsewhere.
FemModel loaded_plate(Mesh mesh, Scalar thickness, Scalar point_force) {
  const int dim = mesh.dim();
  const Scalar length = 0.6;
  FemModel model(std::move(mesh), steel(), thickness,
                 dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                 IntegrationOptions());
  DisplacementConstraint root;
  root.region.name = "root";
  root.region.members.push_back(x_range(-1.0, 0.0));
  root.fix_x = root.fix_y = true;
  root.fix_z = dim == 3;
  model.constraints().push_back(root);
  DisplacementConstraint end;
  end.region.name = "end";
  end.region.members.push_back(x_range(length, 2.0 * length));
  end.fix_x = true;
  model.constraints().push_back(end);

  PointLoadSpec point;
  point.region.members.push_back(nearest(Vector3(0.5 * length, 0.0, 0.0)));
  point.force = Vector3(0.0, -point_force, 0.0);
  std::vector<LoadCaseSpec> cases(4);
  cases[0].name = "weight";
  cases[0].gravity = Vector3(0.0, -9.81, 0.0);
  cases[0].point_loads.push_back(point);
  cases[1].name = "spin";
  cases[1].centrifugal.enabled = true;
  cases[1].centrifugal.angular_velocity = 20.0;
  BodyForceSpec body;
  body.whole_model = false;
  body.region.members.push_back(x_range(0.5 * length, length));
  body.force_density = Vector3(0.0, -1.0e5, 0.0);
  cases[1].body_forces.push_back(body);
  cases[2].name = "heat";
  cases[2].temperature.source = TemperatureSpec::Source::Uniform;
  cases[2].temperature.uniform = 70.0;
  cases[2].point_loads.push_back(point);
  cases[3].name = "hot_half";
  cases[3].temperature.source = TemperatureSpec::Source::Regions;
  cases[3].temperature.uniform = 20.0;
  RegionValue hot;
  hot.region.members.push_back(x_range(0.5 * length, length));
  hot.value = 120.0;
  cases[3].temperature.regions.push_back(hot);
  model.load_case_specs() = cases;
  model.finalize();
  return model;
}

/// A smooth design in [0.02, 0.98] whose filtered densities reach below the
/// body-load threshold.
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

/// Max scaled error of `analytical` against central differences of `value`
/// at one step, over every `stride`-th variable: each entry judged against
/// max(|analytical|, |FD|, 1e-3 ||analytical||_inf). Of second order,
/// (v(x+h) - v(x-h)) / 2h, or fourth, (-v(x+2h) + 8 v(x+h) - 8 v(x-h) +
/// v(x-2h)) / 12h: an eigenvalue's round-off, which a small step amplifies,
/// needs the larger steps the fourth-order difference affords.
Scalar scaled_error(const Vector& analytical, const Vector& x, Scalar step, int order,
                    const std::function<Scalar(const Vector&)>& value, Eigen::Index stride,
                    Index& tested) {
  const Scalar floor = 1.0e-3 * analytical.cwiseAbs().maxCoeff();
  const Scalar reach = order == 4 ? 2.0 * step : step;
  Scalar worst = 0.0;
  tested = 0;
  Vector xs = x;
  const auto at = [&](Eigen::Index e, Scalar offset) {
    xs(e) = x(e) + offset;
    const Scalar v = value(xs);
    xs(e) = x(e);
    return v;
  };
  for (Eigen::Index e = 0; e < x.size(); e += stride) {
    if (x(e) - reach < 0.0 || x(e) + reach > 1.0) continue;
    const Scalar fd =
        order == 4 ? (-at(e, 2.0 * step) + 8.0 * at(e, step) - 8.0 * at(e, -step) +
                      at(e, -2.0 * step)) / (12.0 * step)
                   : (at(e, step) - at(e, -step)) / (2.0 * step);
    worst = std::max(worst, std::abs(fd - analytical(e)) /
                                std::max({std::abs(fd), std::abs(analytical(e)), floor}));
    ++tested;
  }
  return worst;
}

}  // namespace

StudyOutcome study_design_loads(const std::string& out_dir, json::Value& summary) {
  constexpr Scalar kTolerance = 1.0e-5;
  CsvWriter csv(path_join(out_dir, "design_loads.csv"),
                {"quantity", "mesh", "load_case", "beta[-]", "difference_order", "step[-]",
                 "num_tested", "max_scaled_error[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst_gradient = 0.0;
  bool passed = true;
  const auto record = [&](const std::string& quantity, const std::string& mesh,
                          const std::string& load_case, Scalar beta, const Vector& analytical,
                          const Vector& x, const std::function<Scalar(const Vector&)>& value,
                          Eigen::Index stride) {
    Scalar best = std::numeric_limits<Scalar>::infinity();
    for (const int order : {2, 4}) {
      for (const Scalar step : {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6}) {
        Index tested = 0;
        const Scalar error = scaled_error(analytical, x, step, order, value, stride, tested);
        csv.raw_row({quantity, mesh, load_case, fmt(beta, 3), std::to_string(order),
                     fmt(step, 3), std::to_string(tested), fmt(error)});
        best = std::min(best, error);
      }
    }
    worst_gradient = std::max(worst_gradient, best);
    json::Value rec = json::Value::make_object();
    rec.set("quantity", json::Value::make_string(quantity));
    rec.set("mesh", json::Value::make_string(mesh));
    rec.set("load_case", json::Value::make_string(load_case));
    rec.set("beta", json::Value::make_number(beta));
    rec.set("best_max_scaled_error", json::Value::make_number(best));
    records.push_back(rec);
  };

  // --- the load vectors of a design ------------------------------------------
  Scalar load_error = 0.0;
  {
    FemModel model = loaded_plate(make_structured_quad_mesh(box_spec(12, 6, 1, 0.6, 0.3, 1.0)),
                                  0.01, 100.0);
    const DesignLoads loads(model);
    const Index ne = model.mesh().num_elements();
    const SimpOptions simp;
    for (std::size_t l = 0; l < loads.num_cases(); ++l) {
      const Scalar scale = model.load_vectors()[l].cwiseAbs().maxCoeff();
      load_error = std::max(load_error, (loads.load(l, Vector::Ones(ne), simp) -
                                         model.load_vectors()[l]).cwiseAbs().maxCoeff() / scale);
      for (const Scalar rho : {0.5, 0.05}) {
        const LoadCaseData& data = model.load_case_data(l);
        Vector expected = data.mechanical;
        if (data.body.size() > 0) expected += body_load_factor(rho, simp) * data.body;
        if (data.thermal.size() > 0) expected += simp_stiffness_factor(rho, simp) * data.thermal;
        load_error = std::max(load_error,
                              (loads.load(l, Vector::Constant(ne, rho), simp) - expected)
                                      .cwiseAbs()
                                      .maxCoeff() /
                                  scale);
      }
    }
    passed = passed && load_error <= 1.0e-14;
  }

  // --- compliance gradients ----------------------------------------------------
  Scalar min_density = 1.0;
  {
    struct Case {
      std::string name;
      Mesh mesh;
      Scalar thickness;
      Scalar force;
      Scalar beta;
    };
    std::vector<Case> cases;
    cases.push_back({"Q4 12 x 6", make_structured_quad_mesh(box_spec(12, 6, 1, 0.6, 0.3, 1.0)),
                     0.01, 100.0, 0.0});
    cases.push_back({"Q4 12 x 6", make_structured_quad_mesh(box_spec(12, 6, 1, 0.6, 0.3, 1.0)),
                     0.01, 100.0, 4.0});
    cases.push_back({"Hex8 6 x 3 x 2",
                     make_structured_hex_mesh(box_spec(6, 3, 2, 0.6, 0.3, 0.1)), 1.0, 1000.0,
                     0.0});
    for (Case& c : cases) {
      FemModel model = loaded_plate(std::move(c.mesh), c.thickness, c.force);
      Assembler assembler(model);
      const DensityFilter filter(model.mesh(), FilterType::Density,
                                 1.5 * model.mesh().mean_element_size());
      DesignDomain domain(model, 0.5, 0.5, {});
      ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                    direct_solver());
      if (c.beta > 0.0) objective.set_projection(c.beta, 0.5);
      const Vector x = graded_design(model, domain);
      const ObjectiveEvaluation eval = objective.evaluate(x, true);
      min_density = std::min(min_density, eval.physical_density.minCoeff());
      record("compliance (four cases, weighted)", c.name, "all", c.beta, eval.dc_dx, x,
             [&](const Vector& xx) { return objective.compliance_at(xx); }, 1);
    }
  }

  // --- stress-constraint gradients ---------------------------------------------
  {
    FemModel model = loaded_plate(make_structured_quad_mesh(box_spec(10, 5, 1, 0.6, 0.3, 1.0)),
                                  0.01, 100.0);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                  direct_solver());
    StressConstraintOptions so;
    so.enabled = true;
    so.limit = 5.0e7;
    so.p_norm = 6.0;
    so.relaxation = 0.5;
    const StressConstraint stress(model, assembler, filter, so);
    const Vector x = graded_design(model, domain);
    for (std::size_t l = 0; l < model.load_case_specs().size(); ++l) {
      // A fresh evaluation: the adjoint uses the last factorisation.
      const ObjectiveEvaluation eval = objective.evaluate(x, true);
      const StressEvaluation se = stress.evaluate(objective, eval, l, 0.8, true);
      record("stress constraint (P = 6, q = 0.5)", "Q4 10 x 5",
             model.load_case_specs()[l].name, 0.0, se.dg_dx, x,
             [&](const Vector& xx) {
               const ObjectiveEvaluation e = objective.evaluate(xx, false);
               return stress.evaluate(objective, e, l, 0.8, false).constraint;
             },
             1);
    }
  }

  // --- buckling gradients --------------------------------------------------------
  for (const bool heated : {false, true}) {
    FemModel model(make_structured_quad_mesh(box_spec(20, 4, 1, 0.8, 0.08, 1.0)), steel(),
                   0.01, StressState::PlaneStress, IntegrationOptions());
    DisplacementConstraint root;
    root.region.members.push_back(x_range(-1.0, 0.0));
    root.fix_x = root.fix_y = true;
    model.constraints().push_back(root);
    LoadCaseSpec load;
    load.name = heated ? "heated, both ends held" : "weight and end compression";
    if (heated) {
      DisplacementConstraint end;
      end.region.members.push_back(x_range(0.8, 1.6));
      end.fix_x = true;
      model.constraints().push_back(end);
      load.temperature.source = TemperatureSpec::Source::Uniform;
      load.temperature.uniform = 30.0;
    } else {
      load.gravity = Vector3(-9.81e3, 0.0, 0.0);
      TractionLoadSpec top;
      top.region.members.push_back(x_range(0.8, 1.6));
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
    const Vector dlambda = objective.chain_to_design(eval, be.dlambda_dphysical.front());
    record("lowest buckling load factor", "Q4 20 x 4", load.name, 0.0, dlambda, x,
           [&](const Vector& xx) {
             const ObjectiveEvaluation e = objective.evaluate(xx, false);
             return buckling.evaluate(objective, e, 0, false).load_factors(0);
           },
           std::max<Eigen::Index>(1, x.size() / 20));
  }
  csv.close();

  // --- the parasitic load of near-void material ----------------------------------
  CsvWriter parasitic(path_join(out_dir, "design_loads_parasitic.csv"),
                      {"void_density[-]", "compliance_threshold_0.1[J]",
                       "compliance_no_threshold[J]", "solid_half[J]"});
  Scalar worst_bounded = 0.0;
  Scalar unbounded_ratio = 0.0;
  {
    FemModel model(make_structured_quad_mesh(box_spec(16, 4, 1, 1.6, 0.2, 1.0)), steel(), 0.01,
                   StressState::PlaneStress, IntegrationOptions());
    DisplacementConstraint root;
    root.region.members.push_back(x_range(-1.0, 0.0));
    root.fix_x = root.fix_y = true;
    model.constraints().push_back(root);
    LoadCaseSpec weight;
    weight.name = "weight";
    weight.gravity = Vector3(0.0, -9.81, 0.0);
    PointLoadSpec tip;
    tip.region.members.push_back(nearest(Vector3(0.8, 0.0, 0.0)));
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
    const Scalar solid_half = compliance(0.1, 0.0);
    for (const Scalar rho_v : {1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3}) {
      const Scalar bounded = compliance(0.1, rho_v);
      const Scalar unbounded = compliance(0.0, rho_v);
      parasitic.raw_row({fmt(rho_v, 3), fmt(bounded, 8), fmt(unbounded, 8), fmt(solid_half, 8)});
      if (rho_v <= 1.0e-2) worst_bounded = std::max(worst_bounded, std::abs(bounded / solid_half - 1.0));
      if (rho_v == 1.0e-3) unbounded_ratio = unbounded / solid_half;
    }
    passed = passed && worst_bounded <= 0.02 && unbounded_ratio >= 10.0;
  }
  parasitic.close();

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("records", records);
  block.set("load_vector_max_relative_error", json::Value::make_number(load_error));
  block.set("min_filtered_density", json::Value::make_number(min_density));
  block.set("near_void_bounded_max_relative_excess", json::Value::make_number(worst_bounded));
  block.set("near_void_unbounded_ratio_at_1e-3", json::Value::make_number(unbounded_ratio));
  block.set("note",
            json::Value::make_string(
                "Design-dependent loads (topopt/DesignLoads.hpp): the body loads scale with "
                "gamma(rho) (rho above the threshold 0.1, t [p x^p - (p - 1) x^(p+1)] below "
                "it), the thermal loads with E(rho)/E0. The load vector of a design against "
                "its definition; the compliance, stress-constraint and buckling-load-factor "
                "gradients against central differences of second and fourth order (best of "
                "the steps 1e-3 ... 1e-6; entries judged against max(|analytical|, |FD|, 1e-3 ||gradient||_inf)); and "
                "a cantilever with a near-void outer half under its own weight, whose "
                "compliance with the threshold stays within 2 % of the solid half alone and "
                "without it grows without bound as the half empties."));
  summary.set("design_loads", block);

  StudyOutcome outcome;
  outcome.name = "design-dependent loads: compliance, stress and buckling gradients";
  outcome.kind = "verification";
  outcome.metric =
      "worst case of the best-step max scaled gradient error (the load vectors to 1e-14, the "
      "near-void compliance within 2 % of the solid half and x10 without the threshold also "
      "required)";
  outcome.value = worst_gradient;
  outcome.tolerance = kTolerance;
  outcome.passed = passed && worst_gradient <= kTolerance;
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
