/// \file test_manufacturing.cpp
/// \brief The additive-manufacturing overhang filter and check, the robust
///        (eroded / intermediate / dilated) projection, and the minimum
///        length-scale check.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/LengthScale.hpp"
#include "sparlab/topopt/OverhangFilter.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <random>

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

OverhangOptions build(const std::string& direction, int dim) {
  OverhangOptions o;
  o.filter = true;
  o.direction = parse_build_direction(direction, dim);
  return o;
}

Vector random_density(Index n, unsigned int seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<Scalar> dist(0.05, 0.95);
  Vector v(n);
  for (Index i = 0; i < n; ++i) v(i) = dist(rng);
  return v;
}

/// A clamped-left cantilever block with a tip load, for the optimisation tests.
FemModel make_block(const Mesh& mesh) {
  const int dim = mesh.dim();
  FemModel model(Mesh(mesh), IsotropicMaterial(70.0e9, 0.3, 2700.0, "al"),
                 dim == 2 ? 0.01 : 1.0,
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
  PointLoadSpec tip;
  Selector corner;
  corner.kind = SelectorKind::NearestNode;
  const BoundingBox bb = mesh.bounding_box();
  corner.point = Vector3(bb.upper.x(), bb.lower.y(), bb.lower.z());
  tip.region.members.push_back(corner);
  tip.force = Vector3(0.0, -1.0e3, 0.0);
  load.point_loads.push_back(tip);
  model.load_case_specs().push_back(load);
  model.finalize();
  return model;
}

}  // namespace

TEST_CASE("build directions parse and refuse what the model cannot have", "[manufacturing]") {
  REQUIRE(parse_build_direction("+y", 2).label() == "+y");
  REQUIRE(parse_build_direction("y", 2).label() == "+y");
  REQUIRE(parse_build_direction("-x", 2).label() == "-x");
  REQUIRE(parse_build_direction("-z", 3).axis == 2);
  REQUIRE(parse_build_direction("-z", 3).sense == -1);
  REQUIRE_THROWS_WITH(parse_build_direction("+z", 2), ContainsSubstring("plane"));
  REQUIRE_THROWS_WITH(parse_build_direction("up", 3), ContainsSubstring("+z, -z"));
  // Only structured Q4 / Hex8 grids have layers.
  REQUIRE_THROWS_WITH(OverhangFilter(make_structured_tri_mesh(box_spec(4, 4, 1, 1, 1, 1)),
                                     build("+y", 2)),
                      ContainsSubstring("structured Q4 or Hex8"));
}

TEST_CASE("the overhang filter keeps supported material and removes overhangs",
          "[manufacturing][overhang]") {
  // A 6 x 5 grid built along +y: layer j is row j.
  const Mesh mesh = make_structured_quad_mesh(box_spec(6, 5, 1, 0.6, 0.5, 1.0));
  const OverhangFilter am(mesh, build("+y", 2));
  REQUIRE(am.overhang_angle_degrees() == Approx(45.0));
  const auto elem = [](Index i, Index j) { return j * 6 + i; };
  REQUIRE(am.layer(elem(3, 0)) == 0);
  REQUIRE(am.supports(elem(3, 0)).empty());
  REQUIRE(am.supports(elem(3, 2)) == std::vector<Index>({elem(3, 1), elem(2, 1), elem(4, 1)}));
  REQUIRE(am.supports(elem(0, 2)).size() == 2);  // clipped at the domain edge

  // Solid everywhere: everything is printable (up to the smooth minimum's
  // offset of sqrt(epsilon)/2).
  const Vector full = Vector::Ones(mesh.num_elements());
  const Vector xi_full = am.apply(full);
  REQUIRE((xi_full - full).cwiseAbs().maxCoeff() < 0.01);

  // A floating bar in row 3 over void rows 1-2 (row 0 solid): it is removed,
  // and so is everything the bar would have supported.
  Vector floating = Vector::Zero(mesh.num_elements());
  for (Index i = 0; i < 6; ++i) {
    floating(elem(i, 0)) = 1.0;
    floating(elem(i, 3)) = 1.0;
    floating(elem(i, 4)) = 1.0;
  }
  const Vector xi_float = am.apply(floating);
  for (Index i = 0; i < 6; ++i) {
    REQUIRE(xi_float(elem(i, 0)) == Approx(1.0));
    REQUIRE(xi_float(elem(i, 3)) < 0.02);
    REQUIRE(xi_float(elem(i, 4)) < 0.02);
  }

  // A 45-degree staircase is self-supporting: every step rests diagonally on
  // the one below.
  Vector stairs = Vector::Zero(mesh.num_elements());
  for (Index j = 0; j < 5; ++j) stairs(elem(j, j)) = 1.0;
  const Vector xi_stairs = am.apply(stairs);
  for (Index j = 0; j < 5; ++j) REQUIRE(xi_stairs(elem(j, j)) > 0.95);
  // Built the other way (-y), the same staircase stands on its top step,
  // which rests on the plate, and so is still supported; built along +x it
  // stands on its first step.
  const OverhangFilter down(mesh, build("-y", 2));
  REQUIRE(down.layer(elem(0, 4)) == 0);
  const Vector xi_down = down.apply(stairs);
  for (Index j = 0; j < 5; ++j) REQUIRE(xi_down(elem(j, j)) > 0.95);

  // Passive elements keep their density even when unsupported.
  std::vector<char> passive;
  for (Index e = 0; e < mesh.num_elements(); ++e) passive.push_back(e / 6 == 3 ? 1 : 0);
  const OverhangFilter fixture(mesh, build("+y", 2), passive);
  const Vector xi_fixture = fixture.apply(floating);
  for (Index i = 0; i < 6; ++i) {
    REQUIRE(xi_fixture(elem(i, 3)) == 1.0);
    REQUIRE(xi_fixture(elem(i, 4)) > 0.95);  // supported by the fixture
  }

  // The check counts the unsupported solid elements.
  const Vector volumes = Vector::Constant(mesh.num_elements(), 0.01);
  const OverhangReport report = check_overhang(am, floating, volumes, 0.5);
  REQUIRE(report.solid_elements == 18);
  REQUIRE(report.unsupported_elements == 6);  // row 3; row 4 rests on row 3
  REQUIRE(report.unsupported_fraction == Approx(6.0 / 18.0));
  REQUIRE(report.lowest_unsupported_layer == 3);
  REQUIRE(check_overhang(am, stairs, volumes, 0.5).unsupported_elements == 0);
}

TEST_CASE("the overhang filter's adjoint recursion matches central differences",
          "[manufacturing][overhang][sensitivity][verification]") {
  struct Sample {
    Mesh mesh;
    std::string direction;
  };
  std::vector<Sample> samples;
  samples.push_back({make_structured_quad_mesh(box_spec(7, 6, 1, 0.7, 0.6, 1.0)), "+y"});
  samples.push_back({make_structured_quad_mesh(box_spec(7, 6, 1, 0.7, 0.6, 1.0)), "-x"});
  samples.push_back({make_structured_hex_mesh(box_spec(4, 4, 4, 0.4, 0.4, 0.4)), "+z"});
  samples.push_back({make_structured_hex_mesh(box_spec(4, 3, 5, 0.4, 0.3, 0.5)), "-y"});
  for (const Sample& s : samples) {
    INFO(to_string(s.mesh.element_type()) << " built " << s.direction);
    const OverhangFilter am(s.mesh, build(s.direction, s.mesh.dim()));
    const Index n = s.mesh.num_elements();
    const Vector x = random_density(n, 17u);
    const Vector w = random_density(n, 29u);  // f(xi) = w . xi
    const Vector analytical = am.pull_back(x, w);
    const Scalar h = 1.0e-6;
    Scalar worst = 0.0;
    for (Index e = 0; e < n; ++e) {
      Vector xp = x;
      Vector xm = x;
      xp(e) += h;
      xm(e) -= h;
      const Scalar fd = (w.dot(am.apply(xp)) - w.dot(am.apply(xm))) / (2.0 * h);
      worst = std::max(worst, std::abs(fd - analytical(e)) /
                                  std::max({std::abs(fd), std::abs(analytical(e)), 1.0e-3}));
    }
    // Measured: at most 1.04e-6, the truncation floor of the steep smooth
    // maximum (P = 40) at this step.
    INFO("worst relative error " << worst);
    REQUIRE(worst < 1.0e-5);
  }
}

TEST_CASE("the overhang filter stays finite on near-void densities",
          "[manufacturing][overhang]") {
  // Densities of 1e-8 raised to P = 40 underflow a direct sum of powers, and
  // the derivative's sum^(1/Q - 1) then overflows: an MBB run met exactly
  // this in its void region. Mixed with larger values, the scaled evaluation
  // must also match the textbook formula wherever that one is representable.
  for (const char* direction : {"+y", "-y"}) {
    const Mesh mesh = make_structured_quad_mesh(box_spec(8, 6, 1, 0.8, 0.6, 1.0));
    const OverhangFilter am(mesh, build(direction, 2));
    const Index n = mesh.num_elements();
    // Near void from the plate up (a solid first layer would lift the next
    // layer to about sqrt(eps)/2 through the smooth minimum), so the supports
    // of most elements are all 1e-8: a direct sum of three 1e-320 terms. One
    // column of random density keeps the filter's other branches busy, and
    // zeros and 1e-300 are scattered through the void.
    const Vector r = random_density(n, 41u);
    Vector x = Vector::Constant(n, 1.0e-8);
    for (Index e = 5; e < n; e += 8) x(e) = r(e);
    for (Index e = 1; e < n; e += 5) x(e) = 0.0;
    for (Index e = 3; e < n; e += 7) x(e) = 1.0e-300;
    const Vector xi = am.apply(x);
    const Vector g = am.pull_back(x, Vector::Ones(n));
    INFO("built " << direction);
    REQUIRE(xi.allFinite());
    REQUIRE(g.allFinite());
    REQUIRE(xi.minCoeff() >= 0.0);
    REQUIRE(xi.maxCoeff() <= 1.0 + 1.0e-12);
  }
  // Where nothing underflows, the scaled smooth maximum is the plain one:
  // two supports at 0.5 and one at 0.3 (Q from n = 3) above an element at 1.
  const Mesh mesh = make_structured_quad_mesh(box_spec(3, 2, 1, 0.3, 0.2, 1.0));
  const OverhangFilter am(mesh, build("+y", 2));
  Vector x(6);
  x << 0.5, 0.5, 0.3, 1.0, 1.0, 1.0;
  const Scalar p = 40.0;
  const Scalar q = p + std::log(3.0) / std::log(0.5);
  const Scalar smax = std::pow(2.0 * std::pow(0.5, p) + std::pow(0.3, p), 1.0 / q);
  const Scalar eps = 1.0e-4;
  const Scalar expected =
      0.5 * (1.0 + smax - std::sqrt((1.0 - smax) * (1.0 - smax) + eps) + std::sqrt(eps));
  REQUIRE(am.apply(x)(4) == Approx(expected).epsilon(1e-13));
}

TEST_CASE("robust projection: eroded <= blueprint <= dilated, and exact gradients",
          "[manufacturing][robust][sensitivity]") {
  const Mesh mesh = make_structured_quad_mesh(box_spec(12, 6, 1, 0.6, 0.3, 1.0));
  FemModel model = make_block(mesh);
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density,
                             1.5 * model.mesh().mean_element_size());
  DesignDomain domain(model, 0.5, 0.5, {});
  ComplianceObjective objective(model, assembler, filter, domain, SimpOptions(),
                                StaticAnalysisOptions());
  objective.set_projection(6.0, 0.6);
  const OverhangFilter am(mesh, build("+y", 2));
  for (const bool with_overhang : {false, true}) {
    INFO("overhang filter " << with_overhang);
    objective.set_overhang(with_overhang ? &am : nullptr);
    Vector x = random_density(mesh.num_elements(), 5u);
    const Vector eroded = objective.physical_density_at(x, 0.6);
    const Vector blueprint = objective.physical_density_at(x, 0.5);
    const Vector dilated = objective.physical_density_at(x, 0.4);
    REQUIRE((blueprint - eroded).minCoeff() >= 0.0);
    REQUIRE((dilated - blueprint).minCoeff() >= 0.0);
    REQUIRE((objective.physical_density(x) - eroded).cwiseAbs().maxCoeff() == 0.0);
    // The dilated volume's gradient through projection, overhang filter and
    // density filter.
    const ObjectiveEvaluation eval = objective.evaluate(x, true);
    const Vector grad = objective.chain_at(eval, 0.4, domain.element_volumes());
    Scalar worst = 0.0;
    const Scalar h = 1.0e-6;
    for (Index e = 0; e < mesh.num_elements(); e += 3) {
      Vector xp = x;
      Vector xm = x;
      xp(e) += h;
      xm(e) -= h;
      const Scalar fd = (domain.volume_of(objective.physical_density_at(xp, 0.4)) -
                         domain.volume_of(objective.physical_density_at(xm, 0.4))) /
                        (2.0 * h);
      worst = std::max(worst, std::abs(fd - grad(e)) /
                                  std::max({std::abs(fd), std::abs(grad(e)),
                                            1.0e-3 * grad.cwiseAbs().maxCoeff()}));
    }
    INFO("dilated volume gradient error " << worst);
    REQUIRE(worst < 1.0e-6);
    // And the compliance gradient through the overhang filter.
    Scalar c_worst = 0.0;
    for (Index e = 0; e < mesh.num_elements(); e += 5) {
      Vector xp = x;
      Vector xm = x;
      // A step of 1e-4: below it the solver's round-off in the compliance
      // dominates (measured: 8e-5 at 1e-5, 2e-3 at 1e-6), above it the
      // curvature of the beta = 6 projection.
      xp(e) += 1.0e-4;
      xm(e) -= 1.0e-4;
      const Scalar fd = (objective.compliance_at(xp) - objective.compliance_at(xm)) / 2.0e-4;
      c_worst = std::max(c_worst, std::abs(fd - eval.dc_dx(e)) /
                                      std::max({std::abs(fd), std::abs(eval.dc_dx(e)),
                                                1.0e-3 * eval.dc_dx.cwiseAbs().maxCoeff()}));
    }
    INFO("compliance gradient error " << c_worst);
    REQUIRE(c_worst < 3.0e-5);
  }
}

TEST_CASE("robust and overhang-filtered optimisations report their designs",
          "[manufacturing][robust][overhang][topopt]") {
  const Mesh mesh = make_structured_quad_mesh(box_spec(30, 15, 1, 0.6, 0.3, 1.0));
  FemModel model = make_block(mesh);
  Assembler assembler(model);
  const DensityFilter filter(model.mesh(), FilterType::Density,
                             2.0 * model.mesh().mean_element_size());
  DesignDomain domain(model, 0.4, -1.0, {});
  for (const OptimizerMethod method :
       {OptimizerMethod::OptimalityCriteria, OptimizerMethod::MMA}) {
    INFO("method " << to_string(method));
    TopologyOptimizerOptions options;
    options.method = method;
    options.max_iterations = 60;
    options.change_tolerance = 1.0e-3;
    options.projection.enabled = true;
    options.projection.beta_start = 2.0;
    options.projection.beta_max = 8.0;
    options.projection.beta_interval = 15;
    options.projection.robust = true;
    options.projection.robust_delta = 0.1;
    options.overhang = build("+y", 2);
    options.history_stride = 0;
    TopologyOptimizer optimizer(model, assembler, filter, domain, options);
    const TopologyOptimizationResult result = optimizer.run();
    REQUIRE(result.robust);
    REQUIRE(result.overhang_filtered);
    const RobustRecord& rr = result.robust_record;
    // Thinner parts are softer: the eroded design is the worst case.
    REQUIRE(rr.compliance_eroded > rr.compliance_intermediate);
    REQUIRE(rr.compliance_intermediate > rr.compliance_dilated);
    REQUIRE(rr.volume_fraction_eroded < rr.volume_fraction_intermediate);
    REQUIRE(rr.volume_fraction_intermediate < rr.volume_fraction_dilated);
    // The blueprint is what is reported, and it meets the volume fraction up
    // to the lag of the target rescaling.
    REQUIRE(result.volume_fraction == Approx(rr.volume_fraction_intermediate));
    REQUIRE(result.compliance == Approx(rr.compliance_intermediate));
    REQUIRE(std::abs(result.volume_fraction - 0.4) < 0.02);
    REQUIRE(result.printable_density.size() == mesh.num_elements());
    REQUIRE(result.history.back().dilated_volume_fraction > result.history.back().volume_fraction);
  }
  // The erosion check of a non-robust projected run: the same three
  // thresholds, evaluated once at the end, with the design itself untouched.
  {
    TopologyOptimizerOptions options;
    options.max_iterations = 40;
    options.change_tolerance = 1.0e-3;
    options.projection.enabled = true;
    options.projection.beta_start = 2.0;
    options.projection.beta_max = 8.0;
    options.projection.beta_interval = 10;
    options.projection.erosion_check = true;
    options.projection.robust_delta = 0.15;
    options.history_stride = 0;
    TopologyOptimizer optimizer(model, assembler, filter, domain, options);
    const TopologyOptimizationResult result = optimizer.run();
    REQUIRE_FALSE(result.robust);
    REQUIRE(result.erosion_checked);
    const RobustRecord& rr = result.robust_record;
    REQUIRE(rr.eta_eroded == Approx(0.65));
    REQUIRE(rr.eta_dilated == Approx(0.35));
    REQUIRE(rr.compliance_intermediate == Approx(result.compliance));
    REQUIRE(rr.compliance_eroded > rr.compliance_intermediate);
    REQUIRE(rr.compliance_intermediate > rr.compliance_dilated);
    REQUIRE(rr.volume_fraction_eroded < result.volume_fraction);
    REQUIRE(rr.eroded_density.size() == mesh.num_elements());
    // Evaluating the variants leaves the projection where the run ended.
    REQUIRE(result.projection_eta == Approx(0.5));
  }
  // The robust formulation and the erosion check need the projection; the
  // overhang filter needs a structured grid and the density filter.
  TopologyOptimizerOptions bad;
  bad.projection.robust = true;
  REQUIRE_THROWS_WITH(TopologyOptimizer(model, assembler, filter, domain, bad),
                      ContainsSubstring("projection"));
  TopologyOptimizerOptions bad_check;
  bad_check.projection.erosion_check = true;
  REQUIRE_THROWS_WITH(TopologyOptimizer(model, assembler, filter, domain, bad_check),
                      ContainsSubstring("erosion_check"));
  const DensityFilter heuristic(model.mesh(), FilterType::Sensitivity, 0.03);
  TopologyOptimizerOptions heuristic_options;
  heuristic_options.overhang = build("+y", 2);
  REQUIRE_THROWS_WITH(TopologyOptimizer(model, assembler, heuristic, domain, heuristic_options),
                      ContainsSubstring("overhang"));
}

TEST_CASE("the length-scale check flags thin members and narrow gaps",
          "[manufacturing][length_scale]") {
  // 40 x 20 cells of 1 mm.
  const Mesh mesh = make_structured_quad_mesh(box_spec(40, 20, 1, 0.04, 0.02, 1.0));
  const Vector volumes = Vector::Constant(mesh.num_elements(), 1.0e-6);
  const auto elem = [](Index i, Index j) { return j * 40 + i; };
  // A vertical bar 1 cell wide and one 6 cells wide, on void.
  Vector density = Vector::Zero(mesh.num_elements());
  for (Index j = 2; j < 18; ++j) {
    density(elem(8, j)) = 1.0;
    for (Index i = 20; i < 26; ++i) density(elem(i, j)) = 1.0;
  }
  const LengthScaleReport thin = check_length_scale(mesh, density, volumes, 0.5, 1.5e-3);
  REQUIRE(thin.solid_elements == 16 * 7);
  // The opening keeps the wide bar and removes the thin one.
  REQUIRE(thin.solid_violations == 16);
  REQUIRE(thin.solid_violation_fraction == Approx(16.0 / 112.0));
  // A 1-cell gap between two blocks is closed; the open field is not.
  Vector blocks = Vector::Zero(mesh.num_elements());
  for (Index j = 5; j < 15; ++j) {
    for (Index i = 5; i < 15; ++i) blocks(elem(i, j)) = 1.0;
    for (Index i = 16; i < 26; ++i) blocks(elem(i, j)) = 1.0;
  }
  const LengthScaleReport gap = check_length_scale(mesh, blocks, volumes, 0.5, 1.5e-3);
  REQUIRE(gap.solid_violations == 0);
  REQUIRE(gap.void_violations == 10);
  // A probe of 4.5 cells keeps the 10-cell blocks but rounds their convex
  // corners: three cells at each of the eight corners.
  REQUIRE(check_length_scale(mesh, blocks, volumes, 0.5, 4.5e-3).solid_violations == 24);
  REQUIRE_THROWS_AS(check_length_scale(mesh, blocks, volumes, 0.5, 0.0), ConfigError);

  // The scan in half-cell steps: the one-cell bar survives the half-cell
  // probe and not the one-cell probe, so the members measure 1 mm; the
  // one-cell gap measures 1 mm the same way.
  const LengthScaleScan bars =
      scan_length_scale(mesh, density, volumes, 0.5, 0.5e-3, 6.0e-3, 0.02);
  REQUIRE(bars.solid_min_size == Approx(1.0e-3));
  REQUIRE_FALSE(bars.solid_bound_reached_cap);
  const LengthScaleScan gaps =
      scan_length_scale(mesh, blocks, volumes, 0.5, 0.5e-3, 6.0e-3, 0.01);
  REQUIRE(gaps.void_min_size == Approx(1.0e-3));  // the gap is 1.7 % of the void
  // Opening also rounds convex corners: the one-cell probe (a plus-shaped
  // neighbourhood) takes the four corner cells of each block, 4 % of the
  // solid, so at 1 % the square blocks measure only one cell. At 5 % the
  // corners pass until the 2-cell probe takes three cells per corner (12 %).
  REQUIRE(gaps.solid_min_size == Approx(1.0e-3));
  REQUIRE(scan_length_scale(mesh, blocks, volumes, 0.5, 0.5e-3, 6.0e-3, 0.05).solid_min_size ==
          Approx(3.0e-3));
}
