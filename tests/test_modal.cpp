/// \file test_modal.cpp
/// \brief Mass-matrix properties, eigen-solver verification and frequency
///        validation against beam theory.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

TEST_CASE("global mass matrix is symmetric, positive definite and conserves mass",
          "[modal][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 6;
  spec.ny = 4;
  spec.lx = 1.2;
  spec.ly = 0.5;
  const Scalar thickness = 0.01;
  const Scalar density = 2700.0;

  FemModel model(make_perturbed_quad_mesh(spec, 0.2, 11u),
                 IsotropicMaterial(70.0e9, 0.3, density), thickness,
                 StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint bc;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  bc.region.members.push_back(box);
  bc.fix_x = true;
  bc.fix_y = true;
  model.constraints().push_back(bc);
  LoadCaseSpec load;
  load.name = "unused";
  PointLoadSpec p;
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector3(spec.lx, 0.0, 0.0);
  p.region.members.push_back(nearest);
  p.force = Vector3(0.0, -1.0, 0.0);
  load.point_loads.push_back(p);
  model.load_case_specs().push_back(load);
  model.finalize();

  Assembler assembler(model);
  const Scalar expected_mass = density * spec.lx * spec.ly * thickness;

  for (MassType type : {MassType::Consistent, MassType::Lumped}) {
    const SparseMatrix m = assembler.assemble_mass(type);
    const Matrix dense(m);
    REQUIRE((dense - dense.transpose()).cwiseAbs().maxCoeff() /
                dense.cwiseAbs().maxCoeff() == Approx(0.0).margin(1.0e-15));

    // sum(M) = 2 * total mass, because each translation direction contributes.
    REQUIRE(m.sum() == Approx(2.0 * expected_mass).epsilon(1.0e-11));

    // Positive definite on the free DOFs.
    const SparseMatrix mff = assembler.reduce_free_free(m);
    const Matrix dense_mff(mff);
    Eigen::SelfAdjointEigenSolver<Matrix> es(dense_mff);
    REQUIRE(es.eigenvalues().minCoeff() > 0.0);

    // Rigid translation kinetic energy equals the total mass.
    Vector tx = Vector::Zero(model.dofs().num_dofs());
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) tx(n * 2) = 1.0;
    REQUIRE(tx.dot(m * tx) == Approx(expected_mass).epsilon(1.0e-11));
  }

  // total_mass() reproduces the same value analytically.
  REQUIRE(assembler.total_mass() == Approx(expected_mass).epsilon(1.0e-12));
}

TEST_CASE("mass assembly requires a positive density", "[modal][diagnostics]") {
  StructuredMeshSpec spec;
  spec.nx = 2;
  spec.ny = 2;
  FemModel model(make_structured_quad_mesh(spec), IsotropicMaterial(1.0e9, 0.3, 0.0),
                 0.01, StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint bc;
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  bc.region.members.push_back(box);
  bc.fix_x = true;
  bc.fix_y = true;
  model.constraints().push_back(bc);
  model.finalize(/*require_load_cases=*/false);

  Assembler assembler(model);
  REQUIRE_THROWS_AS(assembler.assemble_mass(MassType::Consistent), ModelError);
}

TEST_CASE("subspace iteration agrees with a dense generalised eigensolve",
          "[modal][verification]") {
  // A mesh large enough to take the iterative path (> 400 free DOFs) but small
  // enough for a dense reference.
  CantileverCase c;
  c.poisson = 0.3;
  FemModel model = make_cantilever(c, 40, 8);
  Assembler assembler(model);

  ModalAnalysisOptions options;
  options.num_modes = 5;
  options.tolerance = 1.0e-12;
  const ModalResult modal = solve_modal(model, assembler, options);
  REQUIRE(modal.converged);
  REQUIRE(model.dofs().num_free() > 400);

  const SparseMatrix k =
      assembler.reduce_free_free(assembler.assemble_stiffness());
  const SparseMatrix m =
      assembler.reduce_free_free(assembler.assemble_mass(MassType::Consistent));
  const Matrix dense_k(k);
  const Matrix dense_m(m);
  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> ges(dense_k, dense_m);
  REQUIRE(ges.info() == Eigen::Success);

  for (int i = 0; i < options.num_modes; ++i) {
    const Scalar reference = ges.eigenvalues()(i);
    INFO("mode " << i << ": subspace " << modal.eigenvalues(i) << ", dense "
                 << reference);
    REQUIRE(modal.eigenvalues(i) == Approx(reference).epsilon(1.0e-8));
  }

  // Eigenvalues ascending, residuals small, modes M-orthonormal.
  for (Eigen::Index i = 1; i < modal.eigenvalues.size(); ++i) {
    REQUIRE(modal.eigenvalues(i) >= modal.eigenvalues(i - 1));
  }
  REQUIRE(modal.modal_residuals.maxCoeff() < 1.0e-8);
  for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
    const Vector phi = model.dofs().restrict_to_free(modal.mode_shapes.col(i));
    REQUIRE(phi.dot(m * phi) == Approx(1.0).epsilon(1.0e-8));
  }
  // Mode shapes vanish at prescribed DOFs.
  for (Index d : model.dofs().constrained_dofs()) {
    for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
      REQUIRE(modal.mode_shapes(d, i) == Approx(0.0));
    }
  }
}

TEST_CASE("cantilever bending frequencies match Euler-Bernoulli theory",
          "[modal][validation]") {
  CantileverCase c;
  c.poisson = 0.0;
  FemModel model = make_cantilever(c, 120, 12);
  Assembler assembler(model);

  ModalAnalysisOptions options;
  options.num_modes = 4;
  const ModalResult modal = solve_modal(model, assembler, options);
  REQUIRE(modal.converged);

  const Scalar f1_theory = cantilever_bending_frequency(
      c.youngs, c.density, c.length, c.height, c.thickness, 1);
  const Scalar f2_theory = cantilever_bending_frequency(
      c.youngs, c.density, c.length, c.height, c.thickness, 2);

  const Scalar f1 = modal.frequencies_hz(0);
  INFO("f1 FEM " << f1 << " Hz, theory " << f1_theory << " Hz");
  REQUIRE(std::abs(f1 - f1_theory) / f1_theory < 0.01);

  // The 2-D model includes shear deformation and rotary inertia, which
  // Euler-Bernoulli omits, so the computed frequency must be *below* theory.
  REQUIRE(f1 < f1_theory);

  // Find the second bending frequency among the computed set.
  Scalar f2 = 0.0;
  Scalar best = std::numeric_limits<Scalar>::max();
  for (Eigen::Index i = 1; i < modal.frequencies_hz.size(); ++i) {
    const Scalar d = std::abs(modal.frequencies_hz(i) - f2_theory);
    if (d < best) {
      best = d;
      f2 = modal.frequencies_hz(i);
    }
  }
  INFO("f2 FEM " << f2 << " Hz, theory " << f2_theory << " Hz");
  REQUIRE(std::abs(f2 - f2_theory) / f2_theory < 0.06);
  REQUIRE(f2 < f2_theory);

  // All frequencies are strictly positive and ascending.
  REQUIRE(modal.frequencies_hz(0) > 0.0);
  for (Eigen::Index i = 1; i < modal.frequencies_hz.size(); ++i) {
    REQUIRE(modal.frequencies_hz(i) >= modal.frequencies_hz(i - 1));
  }
  REQUIRE(modal.warnings.empty());
}

TEST_CASE("the first axial mode matches fixed-free rod theory",
          "[modal][validation]") {
  // A plane-stress cantilever carries axial modes interleaved with its bending
  // modes. With nu = 0 the plane-stress model is exactly a uniaxial bar, so the
  // axial wave speed sqrt(E/rho) applies and f_1 = c / (4 L) is an *exact*
  // reference - a different analytical check from the bending series.
  CantileverCase c;
  c.poisson = 0.0;
  FemModel model = make_cantilever(c, 160, 16);
  Assembler assembler(model);

  ModalAnalysisOptions options;
  options.num_modes = 8;
  const ModalResult modal = solve_modal(model, assembler, options);
  REQUIRE(modal.converged);

  const Scalar expected =
      rod_axial_frequency(c.youngs, c.density, c.length, 1);

  // Identify the axial mode as the one whose motion is predominantly along x.
  Eigen::Index axial = -1;
  Scalar best_ratio = 0.0;
  for (Eigen::Index m = 0; m < modal.mode_shapes.cols(); ++m) {
    Scalar ex = 0.0;
    Scalar ey = 0.0;
    for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
      const Scalar ux = modal.mode_shapes(n * 2 + 0, m);
      const Scalar uy = modal.mode_shapes(n * 2 + 1, m);
      ex += ux * ux;
      ey += uy * uy;
    }
    const Scalar ratio = ex / (ex + ey);
    if (ratio > best_ratio) {
      best_ratio = ratio;
      axial = m;
    }
  }
  REQUIRE(axial >= 0);
  REQUIRE(best_ratio > 0.9);  // overwhelmingly axial motion

  const Scalar computed = modal.frequencies_hz(axial);
  INFO("axial mode " << axial << ": FEM " << computed << " Hz, rod theory "
                     << expected << " Hz");
  REQUIRE(std::abs(computed - expected) / expected < 0.01);
}

TEST_CASE("the axial-frequency helper is self-consistent", "[modal]") {
  const Scalar base = rod_axial_frequency(70.0e9, 2700.0, 1.0, 1);
  // f = c / (4L) with c = sqrt(E/rho).
  REQUIRE(base == Approx(std::sqrt(70.0e9 / 2700.0) / 4.0));
  REQUIRE(rod_axial_frequency(70.0e9, 2700.0, 2.0, 1) == Approx(0.5 * base));
  REQUIRE(rod_axial_frequency(4.0 * 70.0e9, 2700.0, 1.0, 1) == Approx(2.0 * base));
  // Odd harmonics: f_n / f_1 = 2n - 1.
  REQUIRE(rod_axial_frequency(70.0e9, 2700.0, 1.0, 2) == Approx(3.0 * base));
  REQUIRE(rod_axial_frequency(70.0e9, 2700.0, 1.0, 3) == Approx(5.0 * base));
  REQUIRE_THROWS_AS(rod_axial_frequency(70.0e9, 2700.0, 1.0, 0), ConfigError);
  REQUIRE_THROWS_AS(rod_axial_frequency(-1.0, 2700.0, 1.0, 1), ConfigError);
}

TEST_CASE("frequencies are invariant under uniform thickness scaling",
          "[modal][verification]") {
  // K and M both scale linearly with thickness, so omega must not change.
  // This property is what makes the full solid domain a valid equal-mass
  // uniform baseline in the topology study.
  CantileverCase c;
  ModalAnalysisOptions options;
  options.num_modes = 3;

  CantileverCase thin = c;
  thin.thickness = c.thickness * 0.37;

  FemModel a = make_cantilever(c, 30, 6);
  FemModel b = make_cantilever(thin, 30, 6);
  Assembler aa(a);
  Assembler ab(b);
  const ModalResult ma = solve_modal(a, aa, options);
  const ModalResult mb = solve_modal(b, ab, options);

  for (int i = 0; i < options.num_modes; ++i) {
    REQUIRE(mb.frequencies_hz(i) == Approx(ma.frequencies_hz(i)).epsilon(1.0e-9));
  }
  // The mass, however, scales.
  REQUIRE(mb.total_mass == Approx(0.37 * ma.total_mass).epsilon(1.0e-10));
}

TEST_CASE("lumped mass gives slightly lower frequencies than consistent mass",
          "[modal]") {
  CantileverCase c;
  c.poisson = 0.0;
  FemModel model = make_cantilever(c, 60, 6);
  Assembler assembler(model);

  ModalAnalysisOptions consistent;
  consistent.num_modes = 3;
  consistent.mass_type = MassType::Consistent;
  ModalAnalysisOptions lumped = consistent;
  lumped.mass_type = MassType::Lumped;

  const ModalResult mc = solve_modal(model, assembler, consistent);
  const ModalResult ml = solve_modal(model, assembler, lumped);

  // Row-sum lumping moves mass towards the diagonal, lowering the frequencies;
  // the two must still agree to a few percent on a converged mesh.
  for (int i = 0; i < 3; ++i) {
    INFO("mode " << i << ": consistent " << mc.frequencies_hz(i) << " Hz, lumped "
                 << ml.frequencies_hz(i) << " Hz");
    REQUIRE(ml.frequencies_hz(i) < mc.frequencies_hz(i));
    REQUIRE(std::abs(ml.frequencies_hz(i) - mc.frequencies_hz(i)) /
                mc.frequencies_hz(i) < 0.05);
  }
}

TEST_CASE("an unconstrained model produces rigid-body eigenvalues that are flagged",
          "[modal][diagnostics]") {
  StructuredMeshSpec spec;
  spec.nx = 6;
  spec.ny = 6;
  spec.lx = 0.5;
  spec.ly = 0.5;
  FemModel model(make_structured_quad_mesh(spec), default_material(), 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  model.finalize(/*require_load_cases=*/false);

  Assembler assembler(model);
  ModalAnalysisOptions options;
  options.num_modes = 4;

  // K_ff is singular here. The eigen solver takes the dense path for this size
  // and must either flag the near-zero eigenvalues or report the failure; it
  // must never return them silently as physical frequencies.
  try {
    const ModalResult modal = solve_modal(model, assembler, options);
    REQUIRE_FALSE(modal.warnings.empty());
    // The first three eigenvalues are the rigid-body modes.
    REQUIRE(modal.eigenvalues(0) < 1.0e-6 * modal.eigenvalues(3));
  } catch (const SolverError&) {
    SUCCEED("the singular pair was reported as a solver error");
  }
}

TEST_CASE("modal analysis validates its own inputs", "[modal][diagnostics]") {
  FemModel model = make_small_plate();
  Assembler assembler(model);
  ModalAnalysisOptions options;
  options.num_modes = 0;
  REQUIRE_THROWS_AS(solve_modal(model, assembler, options), ConfigError);

  // Requesting more modes than there are free DOFs warns and truncates.
  options.num_modes = model.dofs().num_free() + 5;
  const ModalResult modal = solve_modal(model, assembler, options);
  REQUIRE(modal.eigenvalues.size() == model.dofs().num_free());
  REQUIRE_FALSE(modal.warnings.empty());
}

TEST_CASE("the analytical cantilever frequency helper is self-consistent",
          "[modal]") {
  // f scales as sqrt(E), as 1/L^2, as h and as 1/sqrt(rho).
  const Scalar base = cantilever_bending_frequency(70.0e9, 2700.0, 1.0, 0.1, 0.01, 1);
  REQUIRE(cantilever_bending_frequency(4.0 * 70.0e9, 2700.0, 1.0, 0.1, 0.01, 1) ==
          Approx(2.0 * base));
  REQUIRE(cantilever_bending_frequency(70.0e9, 4.0 * 2700.0, 1.0, 0.1, 0.01, 1) ==
          Approx(0.5 * base));
  REQUIRE(cantilever_bending_frequency(70.0e9, 2700.0, 2.0, 0.1, 0.01, 1) ==
          Approx(0.25 * base));
  REQUIRE(cantilever_bending_frequency(70.0e9, 2700.0, 1.0, 0.2, 0.01, 1) ==
          Approx(2.0 * base));
  // Independent of thickness, as expected for in-plane bending.
  REQUIRE(cantilever_bending_frequency(70.0e9, 2700.0, 1.0, 0.1, 0.05, 1) ==
          Approx(base));
  // Mode ratios follow (beta_n L)^2.
  const Scalar f2 = cantilever_bending_frequency(70.0e9, 2700.0, 1.0, 0.1, 0.01, 2);
  REQUIRE(f2 / base == Approx((4.69409113 * 4.69409113) / (1.87510407 * 1.87510407))
                           .epsilon(1.0e-9));

  REQUIRE_THROWS_AS(cantilever_bending_frequency(70.0e9, 2700.0, 1.0, 0.1, 0.01, 0),
                    ConfigError);
  REQUIRE_THROWS_AS(cantilever_bending_frequency(70.0e9, 2700.0, 1.0, 0.1, 0.01, 6),
                    ConfigError);
  REQUIRE_THROWS_AS(cantilever_bending_frequency(-1.0, 2700.0, 1.0, 0.1, 0.01, 1),
                    ConfigError);
}
