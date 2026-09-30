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
///
/// Study `shell-topology`: a density design on MITC4 shells - the compliance
/// gradient of a plate and a cylinder panel under a pressure, an edge load
/// and their weight, and the gradients of the lowest out-of-plane buckling
/// load factor and its KS aggregate on a compressed plate, against central
/// differences; a buckling-constrained run on that plate; and the solid an
/// exported shell part is thickened into.
///
/// Study `part-check`: the non-linear check of an exported part - its ratios
/// of the non-linear to the linear response against Euler's elastica for a
/// slender strip (with the tip section's rigid rotation for the largest
/// displacement) over a mesh ladder, and their order in the small-load limit;
/// the linear first-yield estimate and the collapse bracket of a uniform bar
/// against the exact collapse load; the bifurcation bracket of a column
/// against the linear buckling factor of the same mesh; and free thermal
/// expansion, which both analyses reproduce exactly.
#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/StlWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"
#include "sparlab/topopt/BucklingConstraint.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/DesignLoads.hpp"
#include "sparlab/topopt/NonlinearPartCheck.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/StressConstraint.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

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

constexpr Scalar kPiTopo = 3.14159265358979323846;

Selector plane_box(Scalar xmin, Scalar xmax, Scalar ymin, Scalar ymax) {
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  s.ymin = ymin;
  s.ymax = ymax;
  return s;
}

/// A plate 0.6 x 0.3 m (8 x 4 cells) or a cylinder panel of radius 1 m, 40
/// degrees round and 0.3 m long, of 5 mm aluminium, clamped along its first
/// straight edge, under a pressure of 2 kPa, an edge traction and its weight
/// (at 100 g) with a point force.
FemModel loaded_shell(bool curved) {
  ShellMeshSpec spec;
  spec.n1 = 8;
  spec.n2 = 4;
  if (curved) {
    spec.shape = ShellShape::Cylinder;
    spec.radius = 1.0;
    spec.length = 0.3;
    spec.axis = 1;
    spec.angle_start = 0.0;
    spec.angle_end = 40.0;
  } else {
    spec.lx = 0.6;
    spec.ly = 0.3;
  }
  FemModel model(make_structured_shell_mesh(spec),
                 IsotropicMaterial(70.0e9, 0.3, 2700.0, "aluminium"), 0.005,
                 StressState::Shell, IntegrationOptions());
  const Scalar angle = 40.0 * kPiTopo / 180.0;
  DisplacementConstraint root;
  root.region.name = "root";
  root.region.members.push_back(plane_box(-1.0, curved ? 1.0e-9 : 0.0, -1.0, 1.0));
  if (curved) root.region.members.back().zmin = 1.0 - 1.0e-9;
  for (int k = 0; k < 6; ++k) root.set(k, true);
  model.constraints().push_back(root);
  const Scalar far_x = curved ? std::sin(angle) : 0.6;
  LoadCaseSpec pressure;
  pressure.name = "pressure";
  PressureLoadSpec p;
  p.region.members.push_back(Selector());
  p.pressure = 2.0e3;
  pressure.pressures.push_back(p);
  LoadCaseSpec edge;
  edge.name = "edge";
  edge.weight = 0.5;
  TractionLoadSpec t;
  t.region.members.push_back(plane_box(far_x - 1.0e-6, 2.0, -1.0, 1.0));
  if (curved) t.region.members.back().zmax = std::cos(angle) + 1.0e-6;
  t.traction = curved ? Vector3(0.0, 2.0e5, 0.0) : Vector3(1.0e6, 2.0e5, 0.0);
  edge.tractions.push_back(t);
  LoadCaseSpec weight;
  weight.name = "weight";
  weight.gravity = Vector3(0.0, 0.0, -9.81e2);
  PointLoadSpec point;
  point.region.members.push_back(nearest(curved ? Vector3(far_x, 0.3, std::cos(angle))
                                                : Vector3(0.6, 0.3, 0.0)));
  point.force = Vector3(0.0, 0.0, -50.0);
  weight.point_loads.push_back(point);
  model.load_case_specs() = {pressure, edge, weight};
  model.finalize();
  return model;
}

/// A plate 0.6 x 0.3 m, 3 mm thick, simply supported and compressed along x
/// by 1 MPa: it buckles out of its plane.
FemModel compressed_plate(Index n1, Index n2) {
  ShellMeshSpec spec;
  spec.lx = 0.6;
  spec.ly = 0.3;
  spec.n1 = n1;
  spec.n2 = n2;
  FemModel model(make_structured_shell_mesh(spec),
                 IsotropicMaterial(70.0e9, 0.3, 2700.0, "aluminium"), 0.003,
                 StressState::Shell, IntegrationOptions());
  DisplacementConstraint edges;
  edges.region.members = {plane_box(-1.0, 0.0, -1.0, 1.0), plane_box(0.6, 1.0, -1.0, 1.0),
                          plane_box(-1.0, 1.0, -1.0, 0.0), plane_box(-1.0, 1.0, 0.3, 1.0)};
  edges.fix_z = true;
  DisplacementConstraint held_x;
  held_x.region.members = {plane_box(-1.0, 0.0, -1.0, 1.0)};
  held_x.fix_x = true;
  DisplacementConstraint held_y;
  held_y.region.members = {plane_box(-1.0, 1.0, -1.0, 0.0)};
  held_y.fix_y = true;
  DisplacementConstraint drill;
  drill.region.members = {nearest(Vector3::Zero())};
  drill.fix_rz = true;
  model.constraints() = {edges, held_x, held_y, drill};
  LoadCaseSpec compression;
  compression.name = "compression";
  TractionLoadSpec t;
  t.region.members.push_back(plane_box(0.6, 1.0, -1.0, 1.0));
  t.traction = Vector3(-1.0e6, 0.0, 0.0);
  compression.tractions.push_back(t);
  model.load_case_specs() = {compression};
  model.finalize();
  return model;
}

/// A smooth design in [0.05, 0.95].
Vector shell_design(const FemModel& model, const DesignDomain& domain) {
  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    const Vector3 c = model.mesh().element_centroid(e);
    x(e) = 0.5 + 0.45 * std::sin(9.0 * c.x() + 1.0) * std::cos(13.0 * c.y() + 3.0 * c.z());
  }
  domain.clamp(x);
  return x;
}

/// A plane-stress strip `length` x `height` (nx x ny Q4 cells) of a material
/// with nu = 0, clamped at x = 0, with a dead shear traction on its free end
/// of resultant -force along y: the continuum counterpart of Euler's elastica.
FemModel elastica_strip(Index nx, Index ny, Scalar length, Scalar height, Scalar thickness,
                        Scalar youngs, Scalar force) {
  FemModel model(make_structured_quad_mesh(box_spec(nx, ny, 1, length, height, 1.0)),
                 IsotropicMaterial(youngs, 0.0, 7850.0, "strip"), thickness,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint root;
  root.region.name = "root";
  root.region.members.push_back(x_range(-1.0, 0.0));
  root.fix_x = root.fix_y = true;
  model.constraints().push_back(root);
  LoadCaseSpec tip;
  tip.name = "tip";
  TractionLoadSpec shear;
  shear.region.members.push_back(x_range(length, 2.0 * length));
  shear.traction = Vector3(0.0, -force / (height * thickness), 0.0);
  tip.tractions.push_back(shear);
  model.load_case_specs().push_back(tip);
  model.finalize();
  return model;
}

/// A plane-stress bar 1 x 0.1 m (10 x 2 cells) held by rollers at x = 0 and
/// a pinned corner, so that an end load along x or a uniform temperature
/// leaves it in a uniform state.
FemModel uniform_bar(const IsotropicMaterial& material, const LoadCaseSpec& load) {
  FemModel model(make_structured_quad_mesh(box_spec(10, 2, 1, 1.0, 0.1, 1.0)), material, 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint rollers;
  rollers.region.members.push_back(x_range(-1.0, 0.0));
  rollers.fix_x = true;
  model.constraints().push_back(rollers);
  DisplacementConstraint pin;
  pin.region.members.push_back(plane_box(-1.0, 0.0, -1.0, 0.0));
  pin.fix_y = true;
  model.constraints().push_back(pin);
  model.load_case_specs().push_back(load);
  model.finalize();
  return model;
}

/// The part check of load case 0 of a finalised model.
NonlinearPartCase check_case(const FemModel& model, const NonlinearOptions& options) {
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler, direct_solver());
  const std::vector<StaticSolution> linear = analysis.solve_all();
  return check_part_nonlinear(model, assembler, options, {0}, linear).cases.at(0);
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
                "the steps 1e-3 ... 1e-6; entries judged against max(|analytical|, |FD|, "
                "1e-3 ||gradient||_inf)); and a cantilever with a near-void outer half under its own weight, whose "
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

StudyOutcome study_shell_topology(const std::string& out_dir, json::Value& summary) {
  constexpr Scalar kTolerance = 1.0e-5;
  CsvWriter csv(path_join(out_dir, "shell_topology.csv"),
                {"quantity", "model", "difference_order", "step[-]", "num_tested",
                 "max_scaled_error[-]"});
  json::Value records = json::Value::make_array();
  Scalar worst_gradient = 0.0;
  bool passed = true;
  const auto record = [&](const std::string& quantity, const std::string& model_name,
                          const Vector& analytical, const Vector& x,
                          const std::function<Scalar(const Vector&)>& value,
                          Eigen::Index stride) {
    Scalar best = std::numeric_limits<Scalar>::infinity();
    for (const int order : {2, 4}) {
      for (const Scalar step : {3.0e-3, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6}) {
        Index tested = 0;
        const Scalar error = scaled_error(analytical, x, step, order, value, stride, tested);
        csv.raw_row({quantity, model_name, std::to_string(order), fmt(step, 3),
                     std::to_string(tested), fmt(error)});
        best = std::min(best, error);
      }
    }
    worst_gradient = std::max(worst_gradient, best);
    json::Value rec = json::Value::make_object();
    rec.set("quantity", json::Value::make_string(quantity));
    rec.set("model", json::Value::make_string(model_name));
    rec.set("best_max_scaled_error", json::Value::make_number(best));
    records.push_back(rec);
  };

  // --- compliance gradients: a plate and a cylinder panel ---------------------
  for (const bool curved : {false, true}) {
    FemModel model = loaded_shell(curved);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                  direct_solver());
    const Vector x = shell_design(model, domain);
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    record("compliance (pressure, edge load, weight)",
           curved ? "cylinder panel, 8 x 4 MITC4" : "plate, 8 x 4 MITC4", eval.dc_dx, x,
           [&](const Vector& xx) { return objective.compliance_at(xx); }, 1);
  }

  // --- out-of-plane buckling of a compressed plate ----------------------------
  Scalar solid_ratio = 0.0;
  Scalar out_of_plane_share = 0.0;
  {
    FemModel model = compressed_plate(10, 5);
    Assembler assembler(model);
    BucklingOptions plain;
    plain.num_modes = 2;
    const BucklingResult solid = analyse_buckling(model, assembler, 0, plain);
    const Scalar d = 70.0e9 * std::pow(0.003, 3) / (12.0 * (1.0 - 0.09));
    solid_ratio = solid.load_factors(0) / (4.0 * kPiTopo * kPiTopo * d / (0.09 * 0.003) / 1.0e6);
    Scalar in_plane = 0.0;
    Scalar normal = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      in_plane += solid.mode_shapes.col(0).segment<2>(6 * n).squaredNorm();
      normal += std::pow(solid.mode_shapes(6 * n + 2, 0), 2);
    }
    out_of_plane_share = normal / (normal + in_plane);
    passed = passed && out_of_plane_share > 1.0 - 1.0e-6;

    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                  direct_solver());
    BucklingConstraintOptions bo;
    bo.enabled = true;
    bo.min_load_factor = 1.0;
    bo.num_modes = 3;
    bo.ks_parameter = 20.0;
    bo.eigen.tolerance = 1.0e-14;
    bo.eigen.residual_tolerance = 1.0e-11;
    bo.eigen.max_iterations = 2000;
    BucklingConstraint buckling(model, assembler, bo);
    const Vector x = shell_design(model, domain);
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    const BucklingEvaluation be = buckling.evaluate(objective, eval, 0, true);
    const Vector dlambda = objective.chain_to_design(eval, be.dlambda_dphysical.front());
    const Eigen::Index stride = std::max<Eigen::Index>(1, x.size() / 24);
    record("lowest out-of-plane buckling load factor", "compressed plate, 10 x 5 MITC4", dlambda,
           x,
           [&](const Vector& xx) {
             const ObjectiveEvaluation e = objective.evaluate(xx, false);
             return buckling.evaluate(objective, e, 0, false).load_factors(0);
           },
           stride);
    record("KS aggregate of the three lowest", "compressed plate, 10 x 5 MITC4", be.dg_dx, x,
           [&](const Vector& xx) {
             const ObjectiveEvaluation e = objective.evaluate(xx, false);
             return buckling.evaluate(objective, e, 0, false).constraint;
           },
           stride);
  }
  csv.close();

  // --- a buckling-constrained run -------------------------------------------------
  Scalar lambda_start = 0.0;
  Scalar lambda_required = 0.0;
  Scalar lambda_final = 0.0;
  Scalar volume_fraction = 0.0;
  int iterations = 0;
  {
    FemModel model = compressed_plate(12, 6);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.6, 0.6, {});
    TopologyOptimizerOptions options;
    options.method = OptimizerMethod::MMA;
    options.max_iterations = 30;
    options.analysis = direct_solver();
    options.history_stride = 0;
    options.mma.move_limit = 0.1;
    options.buckling.enabled = true;
    options.buckling.num_modes = 4;
    options.buckling.ks_parameter = 40.0;
    ComplianceObjective probe(model, assembler, filter, domain, options.simp, direct_solver());
    BucklingConstraintOptions bo = options.buckling;
    bo.min_load_factor = 1.0;
    BucklingConstraint probe_buckling(model, assembler, bo);
    const ObjectiveEvaluation start = probe.evaluate(domain.initial_design(), false);
    lambda_start = probe_buckling.evaluate(probe, start, 0, false).load_factors(0);
    lambda_required = 1.3 * lambda_start;
    options.buckling.min_load_factor = lambda_required;
    TopologyOptimizer optimizer(model, assembler, filter, domain, options);
    const TopologyOptimizationResult result = optimizer.run();
    lambda_final = result.min_load_factor;
    volume_fraction = result.volume_fraction;
    iterations = result.iterations;
    passed = passed && lambda_final >= 0.99 * lambda_required &&
             volume_fraction <= 0.6 * (1.0 + options.constraint_tolerance) + 1.0e-9;
  }

  // --- the solid an exported shell part stands for ------------------------------
  Scalar plate_volume_error = 0.0;
  Scalar panel_volume_error = 0.0;
  bool closed = true;
  for (const bool curved : {false, true}) {
    const FemModel model = loaded_shell(curved);
    const Vector thickness = Vector::Constant(model.mesh().num_elements(), 0.005);
    Scalar volume = 0.0;
    for (Index e = 0; e < model.mesh().num_elements(); ++e) {
      volume += model.mesh().element_measure(e) * thickness(e);
    }
    const SurfaceStats stats = surface_stats(shell_surface(model.mesh(), thickness));
    closed = closed && stats.closed && stats.non_manifold_edges == 0;
    (curved ? panel_volume_error : plate_volume_error) =
        std::abs(stats.enclosed_volume / volume - 1.0);
  }
  passed = passed && closed && plate_volume_error <= 1.0e-12 && panel_volume_error <= 1.0e-2;

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("records", records);
  block.set("solid_plate_lambda_over_thin_plate_k4", json::Value::make_number(solid_ratio));
  block.set("lowest_mode_out_of_plane_share", json::Value::make_number(out_of_plane_share));
  block.set("constrained_run_lambda_start", json::Value::make_number(lambda_start));
  block.set("constrained_run_lambda_required", json::Value::make_number(lambda_required));
  block.set("constrained_run_lambda_final", json::Value::make_number(lambda_final));
  block.set("constrained_run_volume_fraction", json::Value::make_number(volume_fraction));
  block.set("constrained_run_iterations", json::Value::make_number(iterations));
  block.set("thickened_plate_volume_relative_error", json::Value::make_number(plate_volume_error));
  block.set("thickened_panel_volume_relative_error", json::Value::make_number(panel_volume_error));
  block.set("thickened_surfaces_closed", json::Value::make_bool(closed));
  block.set("note",
            json::Value::make_string(
                "Topology optimisation of MITC4 shells: each cell's density scales its "
                "membrane, bending, shear and drilling stiffness. The compliance gradient "
                "of a plate and a cylinder panel under a pressure, an edge load and their "
                "weight, and the gradients of the lowest out-of-plane buckling load factor "
                "and of the KS aggregate of the three lowest on a simply supported plate "
                "compressed in its plane, against central differences of second and fourth "
                "order (best of the steps 3e-3 ... 1e-6); a buckling-constrained MMA run "
                "asked for 1.3 times the uniform design's load factor at 0.6 of the "
                "volume; and the solid an exported part is thickened into, against A t."));
  summary.set("shell_topology", block);

  StudyOutcome outcome;
  outcome.name = "shell topology: compliance and out-of-plane buckling gradients";
  outcome.kind = "verification";
  outcome.metric =
      "worst case of the best-step max scaled gradient error (the constrained run meeting "
      "its load factor within the volume, the lowest mode out of plane and the thickened "
      "part closed with the plate's volume also required)";
  outcome.value = worst_gradient;
  outcome.tolerance = kTolerance;
  outcome.passed = passed && worst_gradient <= kTolerance;
  return outcome;
}

StudyOutcome study_part_check(const std::string& out_dir, json::Value& summary) {
  // The continuum strip and the elastica differ by its shear and the second
  // order of its bending strain: (h/L)^2 + (k h / 2L)^2 <= 8e-4 here.
  constexpr Scalar kTolerance = 1.0e-3;
  bool passed = true;

  // --- the ratios against Euler's elastica ------------------------------------
  // A strip 1 m long and 20 mm deep: its shear and the difference between the
  // continuum and the beam are of order (h/L)^2 = 4e-4, and the largest
  // bending strain, k h / (2 L), stays below 2 %, where the Saint Venant-
  // Kirchhoff moment differs from E I kappa at second order in it. The
  // fully integrated Q4 cells lock in bending, so the ratios converge to the
  // elastica with the cell size: the three finest meshes give the observed
  // order and the Richardson extrapolation that is compared.
  const Scalar length = 1.0;
  const Scalar height = 0.02;
  const Scalar thickness = 0.01;
  const Scalar youngs = 100.0e9;
  const Scalar inertia = thickness * height * height * height / 12.0;
  const std::vector<std::pair<Index, Index>> ladder = {{50, 2}, {100, 4}, {200, 8}, {400, 16}};
  const std::vector<Scalar> ks = {0.025, 0.05, 0.1, 0.5, 1.0, 2.0};
  const int rk_steps = 20000;
  CsvWriter csv(path_join(out_dir, "part_check_elastica.csv"),
                {"k[-]", "mesh", "compliance_ratio[-]", "compliance_ratio_elastica[-]",
                 "compliance_ratio_error[-]", "displacement_ratio[-]",
                 "displacement_ratio_elastica[-]", "displacement_ratio_error[-]",
                 "stable_load_factor[-]", "verdict"});
  json::Value records = json::Value::make_array();
  json::Value extrapolated = json::Value::make_array();
  Scalar worst = 0.0;        // the extrapolated ratios
  Scalar worst_finest = 0.0; // the finest mesh's
  Scalar lowest_order = std::numeric_limits<Scalar>::infinity();
  Scalar reference_error = 0.0;
  std::vector<Scalar> small_load_deviation;  // finest mesh, k = 0.025, 0.05, 0.1
  for (const Scalar k : ks) {
    const ElasticaEnd end = euler_elastica(k, rk_steps);
    const ElasticaEnd check = euler_elastica(k, 2 * rk_steps);
    reference_error = std::max({reference_error, std::abs(check.deflection - end.deflection),
                                std::abs(check.shortening - end.shortening)});
    // The end compliance P v: the elastica's deflection over the linear
    // k L / 3. The largest displacement is a corner's of the end section,
    // which turns rigidly by theta about the centre line.
    const Scalar compliance_reference = end.deflection / (k / 3.0);
    const Scalar a = 0.5 * height / length;
    const Scalar top = std::hypot(end.shortening - a * std::sin(end.rotation),
                                  end.deflection + a * (1.0 - std::cos(end.rotation)));
    const Scalar bottom = std::hypot(end.shortening + a * std::sin(end.rotation),
                                     end.deflection - a * (1.0 - std::cos(end.rotation)));
    const Scalar displacement_reference =
        std::max(top, bottom) / std::hypot(k / 3.0, a * k / 2.0);
    std::vector<Scalar> compliance_ratios;
    std::vector<Scalar> displacement_ratios;
    for (std::size_t m = 0; m < ladder.size(); ++m) {
      const auto [nx, ny] = ladder[m];
      const FemModel model = elastica_strip(nx, ny, length, height, thickness, youngs,
                                            k * youngs * inertia / (length * length));
      const NonlinearPartCase c = check_case(model, NonlinearOptions());
      const Scalar ec = c.compliance_ratio - compliance_reference;
      const Scalar ed = c.displacement_ratio - displacement_reference;
      const std::string mesh = std::to_string(nx) + " x " + std::to_string(ny);
      csv.raw_row({fmt(k, 6), mesh, fmt(c.compliance_ratio, 10), fmt(compliance_reference, 10),
                   fmt(ec, 4), fmt(c.displacement_ratio, 10), fmt(displacement_reference, 10),
                   fmt(ed, 4), fmt(c.stable_load_factor, 6), c.verdict});
      json::Value rec = json::Value::make_object();
      rec.set("k", json::Value::make_number(k));
      rec.set("mesh", json::Value::make_string(mesh));
      rec.set("compliance_ratio", json::Value::make_number(c.compliance_ratio));
      rec.set("compliance_ratio_elastica", json::Value::make_number(compliance_reference));
      rec.set("displacement_ratio", json::Value::make_number(c.displacement_ratio));
      rec.set("displacement_ratio_elastica", json::Value::make_number(displacement_reference));
      rec.set("verdict", json::Value::make_string(c.verdict));
      records.push_back(rec);
      passed = passed && c.verdict == "carries";
      compliance_ratios.push_back(c.compliance_ratio);
      displacement_ratios.push_back(c.displacement_ratio);
      if (m + 1 == ladder.size()) {
        worst_finest = std::max({worst_finest, std::abs(ec), std::abs(ed)});
        if (k <= 0.1) small_load_deviation.push_back(std::abs(c.compliance_ratio - 1.0));
      }
    }
    // The observed order of the three finest meshes and the Richardson
    // extrapolation with it (second order is the Q4 displacement's; on these
    // meshes it is still approached from below). Where the load is large
    // enough for the mesh error to dominate, it must fall on every
    // refinement.
    const std::size_t n = ladder.size();
    json::Value x = json::Value::make_object();
    x.set("k", json::Value::make_number(k));
    for (const bool compliance : {true, false}) {
      const std::vector<Scalar>& r = compliance ? compliance_ratios : displacement_ratios;
      const Scalar reference = compliance ? compliance_reference : displacement_reference;
      const Scalar coarse_change = std::abs(r[n - 2] - r[n - 3]);
      const Scalar fine_change = std::abs(r[n - 1] - r[n - 2]);
      const Scalar order =
          fine_change > 0.0 && coarse_change > 0.0 ? std::log2(coarse_change / fine_change) : 0.0;
      const bool consistent = order >= 1.5 && order <= 2.5;
      const Scalar limit =
          r[n - 1] + (r[n - 1] - r[n - 2]) / (std::pow(2.0, consistent ? order : 2.0) - 1.0);
      worst = std::max(worst, std::abs(limit - reference));
      const std::string name = compliance ? "compliance" : "displacement";
      x.set(name + "_ratio_extrapolated", json::Value::make_number(limit));
      x.set(name + "_ratio_extrapolated_error", json::Value::make_number(limit - reference));
      x.set(name + "_ratio_observed_order", json::Value::make_number(order));
      if (k >= 0.5) {
        lowest_order = std::min(lowest_order, order);
        passed = passed && consistent;
        for (std::size_t m = 0; m + 1 < n; ++m) {
          passed = passed && std::abs(r[m + 1] - reference) < std::abs(r[m] - reference);
        }
      }
    }
    extrapolated.push_back(x);
  }
  csv.close();
  // The small-load limit: the end compliance deviates from linear at second
  // order in the load (reversing the load mirrors the deflection).
  std::vector<Scalar> small_load_orders;
  for (std::size_t i = 0; i + 1 < small_load_deviation.size(); ++i) {
    const Scalar order = std::log2(small_load_deviation[i + 1] / small_load_deviation[i]);
    small_load_orders.push_back(order);
    passed = passed && std::abs(order - 2.0) <= 0.1;
  }

  // --- first yield and plastic collapse of a uniform bar ------------------------
  const Scalar yield = 250.0e6;
  IsotropicMaterial plastic(200.0e9, 0.3, 7850.0, "steel");
  {
    PlasticityParameters p;
    p.yield_stress = yield;
    plastic.set_plasticity(p);
  }
  LoadCaseSpec pull;
  pull.name = "pull";
  TractionLoadSpec end_traction;
  end_traction.region.members.push_back(x_range(1.0, 2.0));
  end_traction.traction = Vector3(1.5 * yield, 0.0, 0.0);
  pull.tractions.push_back(end_traction);
  NonlinearOptions small_strain;
  small_strain.kinematics = Kinematics::SmallStrain;
  const NonlinearPartCase bar = check_case(uniform_bar(plastic, pull), small_strain);
  const Scalar collapse = 1.0 / 1.5;
  const bool bar_ok = bar.verdict == "fails" && bar.critical_lower <= collapse &&
                      bar.critical_upper > collapse &&
                      bar.critical_upper - bar.critical_lower < 1.0e-3 &&
                      std::abs(bar.linear_first_yield_load_factor - collapse) <= 1.0e-10;
  passed = passed && bar_ok;

  // --- the bifurcation of a column against its linear buckling factor -----------
  // A cantilever column 1 m x 40 mm (40 x 2 cells) under an axial dead load
  // of 80 kN, beyond its critical load. The bifurcation of the straight path
  // differs from the linear buckling factor only by the pre-buckling
  // deformation, of the order of the axial strain P / (E A).
  Scalar lambda_linear = 0.0;
  NonlinearPartCase column;
  const Scalar column_force = 80.0e3;
  const Scalar column_h = 0.04;
  const Scalar column_t = 0.01;
  {
    FemModel model(make_structured_quad_mesh(box_spec(40, 2, 1, 1.0, column_h, 1.0)),
                   IsotropicMaterial(200.0e9, 0.3, 7850.0, "steel"), column_t,
                   StressState::PlaneStress, IntegrationOptions());
    DisplacementConstraint root;
    root.region.members.push_back(x_range(-1.0, 0.0));
    root.fix_x = root.fix_y = true;
    model.constraints().push_back(root);
    LoadCaseSpec axial;
    axial.name = "axial";
    TractionLoadSpec top;
    top.region.members.push_back(x_range(1.0, 2.0));
    top.traction = Vector3(-column_force / (column_h * column_t), 0.0, 0.0);
    axial.tractions.push_back(top);
    model.load_case_specs().push_back(axial);
    model.finalize();
    Assembler assembler(model);
    BucklingOptions options;
    options.num_modes = 1;
    options.tolerance = 1.0e-12;
    options.linear.type = LinearSolverType::SimplicialLdlt;
    lambda_linear = analyse_buckling(model, assembler, 0, options).load_factors(0);
    column = check_case(model, NonlinearOptions());
  }
  const Scalar prebuckling_strain = lambda_linear * column_force / (200.0e9 * column_h * column_t);
  const Scalar column_gap =
      std::max(std::abs(column.critical_lower / lambda_linear - 1.0),
               std::abs(column.critical_upper / lambda_linear - 1.0));
  const bool column_ok = column.verdict == "fails" && column_gap <= 5.0 * prebuckling_strain;
  passed = passed && column_ok;

  // --- free thermal expansion -----------------------------------------------------
  // The finite-strain law splits the free thermal stretch off
  // multiplicatively, so both analyses give u = alpha dT x exactly.
  LoadCaseSpec heat;
  heat.name = "heat";
  heat.temperature.source = TemperatureSpec::Source::Uniform;
  heat.temperature.uniform = 120.0;
  const NonlinearPartCase thermal = check_case(uniform_bar(steel(), heat), NonlinearOptions());
  const Scalar thermal_error = std::max(std::abs(thermal.displacement_ratio - 1.0),
                                        std::abs(thermal.compliance_ratio - 1.0));
  const bool thermal_ok = thermal.verdict == "carries" && thermal_error <= 1.0e-8 &&
                          std::isnan(thermal.von_mises_ratio);
  passed = passed && thermal_ok;

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification"));
  block.set("elastica_records", records);
  block.set("elastica_extrapolated", extrapolated);
  block.set("elastica_reference_error", json::Value::make_number(reference_error));
  block.set("finest_mesh_worst_ratio_error", json::Value::make_number(worst_finest));
  block.set("extrapolated_worst_ratio_error", json::Value::make_number(worst));
  block.set("lowest_observed_order_k_ge_0.5", json::Value::make_number(lowest_order));
  block.set("small_load_compliance_deviation_orders", json::array_of(small_load_orders));
  json::Value b = json::Value::make_object();
  b.set("verdict", json::Value::make_string(bar.verdict));
  b.set("exact_collapse_load_factor", json::Value::make_number(collapse));
  b.set("critical_lower", json::Value::make_number(bar.critical_lower));
  b.set("critical_upper", json::Value::make_number(bar.critical_upper));
  b.set("linear_first_yield_load_factor", json::Value::make_number(bar.linear_first_yield_load_factor));
  b.set("passed", json::Value::make_bool(bar_ok));
  block.set("bar_collapse", b);
  json::Value col = json::Value::make_object();
  col.set("verdict", json::Value::make_string(column.verdict));
  col.set("linear_buckling_load_factor", json::Value::make_number(lambda_linear));
  col.set("critical_lower", json::Value::make_number(column.critical_lower));
  col.set("critical_upper", json::Value::make_number(column.critical_upper));
  col.set("relative_gap", json::Value::make_number(column_gap));
  col.set("prebuckling_axial_strain", json::Value::make_number(prebuckling_strain));
  col.set("passed", json::Value::make_bool(column_ok));
  block.set("column_bifurcation", col);
  json::Value th = json::Value::make_object();
  th.set("verdict", json::Value::make_string(thermal.verdict));
  th.set("ratio_error", json::Value::make_number(thermal_error));
  th.set("stress_ratio_formed", json::Value::make_bool(!std::isnan(thermal.von_mises_ratio)));
  th.set("passed", json::Value::make_bool(thermal_ok));
  block.set("free_thermal_expansion", th);
  block.set("note",
            json::Value::make_string(
                "The non-linear check of an exported part (topopt/NonlinearPartCheck.hpp): "
                "its end-compliance and largest-displacement ratios of the non-linear to the "
                "linear response for a strip 1 m x 20 mm (Q4, nu = 0) under a dead end "
                "force, k = P L^2 / (E I) from 0.025 to 2, against Euler's elastica solved by "
                "shooting (the largest displacement with the end section turning rigidly), "
                "on four meshes; their second order in the small-load limit; the linear "
                "first-yield estimate and the collapse bracket of a uniform "
                "elastic-perfectly plastic bar at 1.5 times its yield load (small strain) "
                "against 1/1.5; the bifurcation bracket of a cantilever column against the "
                "linear buckling factor of the same mesh; and free thermal expansion, "
                "exact in both analyses."));
  summary.set("part_check", block);

  StudyOutcome outcome;
  outcome.name = "non-linear check of the exported part vs the elastica and exact limits";
  outcome.kind = "verification";
  outcome.metric =
      "largest error of the compliance and displacement ratios, Richardson-extrapolated with "
      "the observed order of the three finest meshes, against the elastica (for k >= 0.5 "
      "that order within 1.5 ... 2.5 and the error falling on every refinement; the "
      "small-load order, the bar's collapse bracket, the column's bifurcation and free "
      "thermal expansion also required)";
  outcome.value = worst;
  outcome.tolerance = kTolerance;
  outcome.passed = passed && worst <= kTolerance;
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
