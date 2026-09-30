/// \file test_shell_topology.cpp
/// \brief Topology optimisation of MITC4 shell surfaces: the compliance
///        gradient on a plate and a curved panel, the out-of-plane buckling
///        load factor's gradient against central differences, the resultants
///        of a density design, the exported part's normals and thickened
///        surface, and a buckling-constrained run.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/StlWriter.hpp"
#include "sparlab/mesh/SubMesh.hpp"
#include "sparlab/topopt/BucklingConstraint.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <vector>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {

Selector box(Scalar xmin, Scalar xmax, Scalar ymin, Scalar ymax) {
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  s.ymin = ymin;
  s.ymax = ymax;
  return s;
}

IsotropicMaterial aluminium() { return IsotropicMaterial(70.0e9, 0.3, 2700.0, "aluminium"); }

StaticAnalysisOptions direct_solver() {
  StaticAnalysisOptions options;
  options.linear.type = LinearSolverType::SimplicialLdlt;
  return options;
}

/// A plate 0.6 x 0.3 m in z = 0 (or a cylinder panel of radius 1 m, 40 degrees
/// round and 0.3 m long) of 5 mm aluminium, clamped along x = 0 (along its
/// first straight edge), with three load cases: a pressure, an in-plane (or
/// along-the-axis) edge traction, and its weight with a point force.
FemModel loaded_shell(bool curved) {
  ShellMeshSpec spec;
  if (curved) {
    spec.shape = ShellShape::Cylinder;
    spec.radius = 1.0;
    spec.length = 0.3;
    spec.axis = 1;  // about y: the angle turns from +z towards +x
    spec.angle_start = 0.0;
    spec.angle_end = 40.0;
    spec.n1 = 8;
    spec.n2 = 4;
  } else {
    spec.lx = 0.6;
    spec.ly = 0.3;
    spec.n1 = 8;
    spec.n2 = 4;
  }
  FemModel model(make_structured_shell_mesh(spec), aluminium(), 0.005, StressState::Shell,
                 IntegrationOptions());
  DisplacementConstraint root;
  root.region.name = "root";
  root.region.members.push_back(curved ? box(-1.0, 1.0e-9, -1.0, 1.0)
                                       : box(-1.0, 0.0, -1.0, 1.0));
  if (curved) {
    // The panel's first straight edge is the line x = 0, z = 1.
    root.region.members.back().zmin = 1.0 - 1.0e-9;
  }
  for (int k = 0; k < 6; ++k) root.set(k, true);
  model.constraints().push_back(root);

  const Scalar far_x = curved ? std::sin(40.0 * 3.14159265358979323846 / 180.0) : 0.6;
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
  t.region.members.push_back(box(far_x - 1.0e-6, 2.0, -1.0, 1.0));
  if (curved) t.region.members.back().zmax = 1.0 - (1.0 - std::cos(40.0 * 3.14159265358979323846 / 180.0)) + 1.0e-6;
  t.traction = curved ? Vector3(0.0, 2.0e5, 0.0) : Vector3(1.0e6, 2.0e5, 0.0);
  edge.tractions.push_back(t);
  LoadCaseSpec weight;
  weight.name = "weight";
  weight.gravity = Vector3(0.0, 0.0, -9.81e2);
  PointLoadSpec point;
  Selector corner;
  corner.kind = SelectorKind::NearestNode;
  corner.point = curved ? Vector3(far_x, 0.3, std::cos(40.0 * 3.14159265358979323846 / 180.0))
                        : Vector3(0.6, 0.3, 0.0);
  point.region.members.push_back(corner);
  point.force = Vector3(0.0, 0.0, -50.0);
  weight.point_loads.push_back(point);
  model.load_case_specs() = {pressure, edge, weight};
  model.finalize();
  return model;
}

/// A smooth design in [0.05, 0.95].
Vector graded_design(const FemModel& model, const DesignDomain& domain) {
  Vector x = domain.initial_design();
  for (Index e = 0; e < domain.num_elements(); ++e) {
    const Vector3 c = model.mesh().element_centroid(e);
    x(e) = 0.5 + 0.45 * std::sin(9.0 * c.x() + 1.0) * std::cos(13.0 * c.y() + 3.0 * c.z());
  }
  domain.clamp(x);
  return x;
}

/// Best (over second- and fourth-order central differences and five steps)
/// max relative error of each gradient in `analytical`, entries below 1e-3 of
/// its scale judged against that scale; `values` returns the functionals at a
/// design, one per gradient, so the checks share their evaluations. A thin
/// shell's bending and membrane stiffness differ by orders of magnitude, so
/// its solves carry more round-off than a continuum's and the best step is a
/// larger one, and fourth order, largest step first, goes first. The sweep
/// stops at the first order and step at which every error is at or below
/// `good_enough`: one step showing the agreement settles a check against a
/// tolerance above it.
std::vector<Scalar> best_gradient_errors(const std::vector<Vector>& analytical, const Vector& x,
                                         const std::function<Vector(const Vector&)>& values,
                                         Eigen::Index stride = 1, Scalar good_enough = 0.0) {
  const std::size_t count = analytical.size();
  std::vector<Scalar> floor(count);
  for (std::size_t k = 0; k < count; ++k) {
    floor[k] = 1.0e-3 * analytical[k].cwiseAbs().maxCoeff();
  }
  std::vector<Scalar> best(count, std::numeric_limits<Scalar>::infinity());
  Vector xs = x;
  const auto at = [&](Eigen::Index e, Scalar offset) {
    xs(e) = x(e) + offset;
    const Vector v = values(xs);
    xs(e) = x(e);
    return v;
  };
  for (const int order : {4, 2}) {
    for (const Scalar step : {3.0e-3, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6}) {
      const Scalar reach = order == 4 ? 2.0 * step : step;
      std::vector<Scalar> worst(count, 0.0);
      for (Eigen::Index e = 0; e < x.size(); e += stride) {
        if (x(e) - reach < 0.0 || x(e) + reach > 1.0) continue;
        const Vector fd =
            order == 4 ? Vector((-at(e, 2.0 * step) + 8.0 * at(e, step) - 8.0 * at(e, -step) +
                                 at(e, -2.0 * step)) / (12.0 * step))
                       : Vector((at(e, step) - at(e, -step)) / (2.0 * step));
        for (std::size_t k = 0; k < count; ++k) {
          const Eigen::Index i = static_cast<Eigen::Index>(k);
          const Scalar a = analytical[k](e);
          worst[k] = std::max(worst[k], std::abs(fd(i) - a) /
                                            std::max({std::abs(fd(i)), std::abs(a), floor[k]}));
        }
      }
      bool settled = true;
      for (std::size_t k = 0; k < count; ++k) {
        best[k] = std::min(best[k], worst[k]);
        settled = settled && best[k] <= good_enough;
      }
      if (settled) return best;
    }
  }
  return best;
}

/// One gradient's error, as above.
Scalar best_gradient_error(const Vector& analytical, const Vector& x,
                           const std::function<Scalar(const Vector&)>& value,
                           Eigen::Index stride = 1, Scalar good_enough = 0.0) {
  return best_gradient_errors(
      {analytical}, x, [&](const Vector& xx) { return Vector::Constant(1, value(xx)); }, stride,
      good_enough)[0];
}

/// A plate 0.6 x 0.3 m, 3 mm thick, simply supported (w held on its edges,
/// u at x = 0, v at y = 0) and compressed along x by a traction on x = 0.6:
/// it buckles out of its plane.
FemModel compressed_plate(Index n1, Index n2) {
  ShellMeshSpec spec;
  spec.lx = 0.6;
  spec.ly = 0.3;
  spec.n1 = n1;
  spec.n2 = n2;
  FemModel model(make_structured_shell_mesh(spec), aluminium(), 0.003, StressState::Shell,
                 IntegrationOptions());
  DisplacementConstraint edges;
  edges.region.name = "edges";
  edges.region.members = {box(-1.0, 0.0, -1.0, 1.0), box(0.6, 1.0, -1.0, 1.0),
                          box(-1.0, 1.0, -1.0, 0.0), box(-1.0, 1.0, 0.3, 1.0)};
  edges.fix_z = true;
  DisplacementConstraint held_x;
  held_x.region.members = {box(-1.0, 0.0, -1.0, 1.0)};
  held_x.fix_x = true;
  DisplacementConstraint held_y;
  held_y.region.members = {box(-1.0, 1.0, -1.0, 0.0)};
  held_y.fix_y = true;
  // The drilling rotation of one corner, so the in-plane problem has no
  // rotation left free on a plate held only by translations.
  DisplacementConstraint drill;
  Selector corner;
  corner.kind = SelectorKind::NearestNode;
  corner.point = Vector3::Zero();
  drill.region.members = {corner};
  drill.fix_rz = true;
  model.constraints() = {edges, held_x, held_y, drill};
  LoadCaseSpec compression;
  compression.name = "compression";
  TractionLoadSpec t;
  t.region.members.push_back(box(0.6, 1.0, -1.0, 1.0));
  t.traction = Vector3(-1.0e6, 0.0, 0.0);
  compression.tractions.push_back(t);
  model.load_case_specs() = {compression};
  model.finalize();
  return model;
}

}  // namespace

TEST_CASE("shell topology: compliance gradients on a plate and a curved panel match central "
          "differences",
          "[topopt][shell][sensitivity][verification]") {
  for (const bool curved : {false, true}) {
    INFO((curved ? "cylinder panel" : "plate"));
    FemModel model = loaded_shell(curved);
    Assembler assembler(model);
    const DensityFilter filter(model.mesh(), FilterType::Density,
                               1.5 * model.mesh().mean_element_size());
    DesignDomain domain(model, 0.5, 0.5, {});
    ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                  direct_solver());
    const Vector x = graded_design(model, domain);
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    REQUIRE(eval.compliance > 0.0);
    // Every case loads the shell: out of its plane, in it, and by its weight.
    for (const Scalar c : eval.load_case_compliance) REQUIRE(c > 0.0);
    const Scalar error = best_gradient_error(eval.dc_dx, x, [&](const Vector& xx) {
      return objective.compliance_at(xx);
    });
    INFO("compliance gradient error " << error);
    REQUIRE(error < 1.0e-5);
  }
}

TEST_CASE("shell topology: the out-of-plane buckling load factor's gradient matches central "
          "differences",
          "[topopt][shell][buckling][sensitivity][verification]") {
  // 10 x 5 cells: under 400 free unknowns, the dense eigensolve.
  FemModel model = compressed_plate(10, 5);
  Assembler assembler(model);
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
  const Vector x = graded_design(model, domain);
  const ObjectiveEvaluation eval = objective.evaluate(x, true);
  const BucklingEvaluation be = buckling.evaluate(objective, eval, 0, true);
  REQUIRE(be.load_factors.size() == 3);
  REQUIRE(be.load_factors(0) > 0.0);
  REQUIRE(be.load_factors(1) > 1.01 * be.load_factors(0));
  INFO("lambda = " << be.load_factors.transpose());
  // The solid plate buckles out of its plane near the thin plate's load,
  // k pi^2 D / (b^2 t) with k = 4 (two half-waves along a = 2b): 5.7 % above
  // it on this coarse mesh of five cells a half-wave (the shell-plate-buckling
  // study measures the convergence). Its lowest mode is all w.
  {
    BucklingOptions options;
    options.num_modes = 2;
    const BucklingResult solid = analyse_buckling(model, assembler, 0, options);
    const Scalar d = 70.0e9 * std::pow(0.003, 3) / (12.0 * (1.0 - 0.09));
    const Scalar exact = 4.0 * 9.869604401089358 * d / (0.09 * 0.003) / 1.0e6;
    INFO("solid plate lambda_1 " << solid.load_factors(0) << " against " << exact);
    REQUIRE(solid.load_factors(0) > exact);
    REQUIRE(solid.load_factors(0) < 1.1 * exact);
    Scalar in_plane = 0.0;
    Scalar out_of_plane = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      in_plane += solid.mode_shapes.col(0).segment<2>(6 * n).squaredNorm();
      out_of_plane += std::pow(solid.mode_shapes(6 * n + 2, 0), 2);
    }
    REQUIRE(in_plane <= 1.0e-6 * out_of_plane);
  }
  // lambda_1 and its KS aggregate from the same perturbed solves; the sweep
  // stops once a step agrees within half the tolerance.
  const Vector dlambda = objective.chain_to_design(eval, be.dlambda_dphysical[0]);
  const std::vector<Scalar> errors = best_gradient_errors(
      {dlambda, be.dg_dx}, x,
      [&](const Vector& xx) {
        const ObjectiveEvaluation e = objective.evaluate(xx, false);
        const BucklingEvaluation b = buckling.evaluate(objective, e, 0, false);
        return Vector((Vector(2) << b.load_factors(0), b.constraint).finished());
      },
      std::max<Eigen::Index>(1, x.size() / 24), 5.0e-6);
  INFO("lambda_1 gradient error " << errors[0] << ", KS gradient error " << errors[1]);
  REQUIRE(errors[0] < 1.0e-5);
  REQUIRE(errors[1] < 1.0e-5);
}

TEST_CASE("a shell design's resultants, its exported part's normals and thickened surface",
          "[topopt][shell][geometry]") {
  // Resultants of a density design are those of each element's scaled
  // material.
  FemModel model = loaded_shell(true);
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler, direct_solver());
  const std::vector<StaticSolution> sols = analysis.solve_all();
  Vector scale(model.mesh().num_elements());
  for (Index e = 0; e < scale.size(); ++e) scale(e) = 0.1 + 0.8 * (e % 7) / 6.0;
  const ShellField plain = recover_shell_resultants(model, assembler, sols[0].displacement);
  const ShellField scaled =
      recover_shell_resultants(model, assembler, sols[0].displacement, &scale);
  for (Index e = 0; e < scale.size(); ++e) {
    const ShellResultants& a = plain.element[static_cast<std::size_t>(e)];
    const ShellResultants& b = scaled.element[static_cast<std::size_t>(e)];
    REQUIRE((b.membrane - scale(e) * a.membrane).norm() <= 1.0e-12 * a.membrane.norm() + 1e-300);
    REQUIRE((b.moment - scale(e) * a.moment).norm() <= 1.0e-12 * a.moment.norm() + 1e-300);
    REQUIRE(scaled.element_von_mises(e) ==
            Approx(scale(e) * plain.element_von_mises(e)).epsilon(1.0e-12));
    REQUIRE(scaled.element_strain_energy(e) ==
            Approx(scale(e) * plain.element_strain_energy(e)).epsilon(1.0e-12));
  }
  const Vector wrong = Vector::Ones(3);
  REQUIRE_THROWS_AS(recover_shell_resultants(model, assembler, sols[0].displacement, &wrong),
                    ModelError);

  // The exported part keeps the surface's exact normals at its nodes.
  std::vector<Index> kept;
  for (Index e = 0; e < model.mesh().num_elements(); e += 2) kept.push_back(e);
  const SubMeshResult sub = extract_element_subset(model.mesh(), kept);
  REQUIRE(sub.mesh.has_node_normals());
  for (std::size_t n = 0; n < sub.node_map.size(); ++n) {
    REQUIRE((sub.mesh.node_normals().col(static_cast<Eigen::Index>(n)) -
             model.mesh().node_normals().col(sub.node_map[n]))
                .norm() == 0.0);
  }

  // The thickened surface of a plate is closed and holds exactly A t; of a
  // cylinder panel, closed, and A t to the flat facets' error.
  const FemModel plate = loaded_shell(false);
  for (const bool curved : {false, true}) {
    INFO((curved ? "cylinder panel" : "plate"));
    const FemModel& m = curved ? model : plate;
    Vector thickness = Vector::Constant(m.mesh().num_elements(), 0.005);
    Scalar volume = 0.0;
    for (Index e = 0; e < m.mesh().num_elements(); ++e) {
      volume += m.mesh().element_measure(e) * thickness(e);
    }
    const TriangleSurface surface = shell_surface(m.mesh(), thickness);
    const SurfaceStats stats = surface_stats(surface);
    REQUIRE(stats.closed);
    REQUIRE(stats.non_manifold_edges == 0);
    REQUIRE(stats.enclosed_volume > 0.0);
    if (curved) {
      REQUIRE(std::abs(stats.enclosed_volume / volume - 1.0) < 1.0e-2);
    } else {
      REQUIRE(stats.enclosed_volume == Approx(volume).epsilon(1.0e-12));
    }
  }
  REQUIRE_THROWS_AS(shell_surface(model.mesh(), Vector::Ones(2)), IoError);
}

TEST_CASE("the optimiser designs a compressed plate against out-of-plane buckling",
          "[topopt][shell][buckling][mma]") {
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
  // The uniform starting design's lowest load factor, and a requirement above
  // it that the optimiser has to meet by placing material.
  ComplianceObjective probe(model, assembler, filter, domain, options.simp, direct_solver());
  BucklingConstraintOptions bo = options.buckling;
  bo.min_load_factor = 1.0;
  BucklingConstraint probe_buckling(model, assembler, bo);
  const ObjectiveEvaluation start = probe.evaluate(domain.initial_design(), false);
  const Scalar lambda_start = probe_buckling.evaluate(probe, start, 0, false).load_factors(0);
  options.buckling.min_load_factor = 1.3 * lambda_start;
  TopologyOptimizer optimizer(model, assembler, filter, domain, options);
  const TopologyOptimizationResult result = optimizer.run();
  INFO("lambda start " << lambda_start << ", required " << options.buckling.min_load_factor
                       << ", final " << result.min_load_factor);
  REQUIRE(result.buckling_constrained);
  REQUIRE(result.history.size() >= 2);
  // The design meets the requirement (to the KS aggregate's conservatism)
  // within the volume allowance.
  REQUIRE(result.min_load_factor >= 0.99 * options.buckling.min_load_factor);
  REQUIRE(result.volume_fraction <= 0.6 * (1.0 + options.constraint_tolerance) + 1.0e-9);

  // A shell refuses the overhang filter and the stress constraint.
  TopologyOptimizerOptions overhang = options;
  overhang.overhang.filter = true;
  REQUIRE_THROWS_WITH(TopologyOptimizer(model, assembler, filter, domain, overhang),
                      ContainsSubstring("overhang"));
  StressConstraintOptions so;
  so.enabled = true;
  so.limit = 1.0e8;
  REQUIRE_THROWS_WITH(StressConstraint(model, assembler, filter, so),
                      ContainsSubstring("shell"));
}

