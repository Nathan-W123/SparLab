/// \file test_assembly_solver.cpp
/// \brief Global assembly, Dirichlet partitioning, solver agreement, reaction
///        equilibrium and ill-posed-model detection.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/LinearSolver.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

TEST_CASE("global stiffness is symmetric and annihilates rigid-body modes",
          "[assembly][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 5;
  spec.ny = 4;
  spec.lx = 1.2;
  spec.ly = 0.8;

  FemModel model(make_structured_quad_mesh(spec), default_material(), 0.01,
                 StressState::PlaneStress, IntegrationOptions());
  LoadCaseSpec load;
  load.name = "dummy";
  PointLoadSpec p;
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector3(spec.lx, spec.ly, 0.0);
  p.region.members.push_back(nearest);
  p.force = Vector3(1.0, 0.0, 0.0);
  load.point_loads.push_back(p);
  model.load_case_specs().push_back(load);
  // Deliberately unconstrained so the free-free block is the full matrix.
  model.finalize();

  Assembler assembler(model);
  const SparseMatrix k = assembler.assemble_stiffness();
  REQUIRE(k.rows() == model.dofs().num_dofs());

  const Matrix dense(k);
  const Scalar scale = dense.cwiseAbs().maxCoeff();
  REQUIRE((dense - dense.transpose()).cwiseAbs().maxCoeff() / scale ==
          Approx(0.0).margin(1.0e-14));

  const Matrix modes = rigid_body_modes(model.mesh());
  for (int m = 0; m < 3; ++m) {
    const Vector residual = k * modes.col(m);
    REQUIRE(residual.cwiseAbs().maxCoeff() / (scale * modes.col(m).norm()) ==
            Approx(0.0).margin(1.0e-13));
  }

  // An unconstrained model is not positive definite.
  REQUIRE_FALSE(is_positive_definite(k));
}

TEST_CASE("the element-matrix cache agrees with per-element integration",
          "[assembly][verification]") {
  // A uniform structured mesh triggers the cache; a perturbed one does not.
  // Assembling the *same* uniform mesh both ways must give identical matrices,
  // which is checked by comparing against a mesh with the uniform flag cleared.
  StructuredMeshSpec spec;
  spec.nx = 6;
  spec.ny = 5;
  spec.lx = 1.0;
  spec.ly = 0.6;

  Mesh uniform = make_structured_quad_mesh(spec);
  Mesh generic = uniform;
  StructuredGridInfo info = *uniform.structured_info();
  info.uniform = false;  // force the generic path
  generic.set_structured_info(info);

  const auto build = [&](Mesh mesh) {
    FemModel model(std::move(mesh), default_material(), 0.004, StressState::PlaneStress,
                   IntegrationOptions());
    DisplacementConstraint bc;
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.0;
    bc.region.members.push_back(box);
    bc.fix_x = true;
    bc.fix_y = true;
    model.constraints().push_back(bc);
    LoadCaseSpec load;
    load.name = "tip";
    PointLoadSpec p;
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector3(spec.lx, 0.0, 0.0);
    p.region.members.push_back(nearest);
    p.force = Vector3(0.0, -100.0, 0.0);
    load.point_loads.push_back(p);
    model.load_case_specs().push_back(load);
    model.finalize();
    return model;
  };

  FemModel cached_model = build(uniform);
  FemModel generic_model = build(generic);
  Assembler cached(cached_model);
  Assembler plain(generic_model);
  REQUIRE(cached.uses_element_cache());
  REQUIRE_FALSE(plain.uses_element_cache());

  const Matrix a(cached.assemble_stiffness());
  const Matrix b(plain.assemble_stiffness());
  REQUIRE((a - b).cwiseAbs().maxCoeff() / a.cwiseAbs().maxCoeff() ==
          Approx(0.0).margin(1.0e-15));

  const Matrix ma(cached.assemble_mass(MassType::Consistent));
  const Matrix mb(plain.assemble_mass(MassType::Consistent));
  REQUIRE((ma - mb).cwiseAbs().maxCoeff() / ma.cwiseAbs().maxCoeff() ==
          Approx(0.0).margin(1.0e-15));
}

TEST_CASE("per-element scaling enters the assembly linearly", "[assembly]") {
  FemModel model = make_small_plate();
  Assembler assembler(model);
  const SparseMatrix k1 = assembler.assemble_stiffness();
  const Vector scale = Vector::Constant(model.mesh().num_elements(), 0.25);
  const SparseMatrix k2 = assembler.assemble_stiffness(&scale);
  const Matrix diff = Matrix(k2) - 0.25 * Matrix(k1);
  REQUIRE(diff.cwiseAbs().maxCoeff() / Matrix(k1).cwiseAbs().maxCoeff() ==
          Approx(0.0).margin(1.0e-15));
}

TEST_CASE("assembly rejects malformed scale vectors", "[assembly][diagnostics]") {
  FemModel model = make_small_plate();
  Assembler assembler(model);

  const Vector wrong_length = Vector::Ones(3);
  REQUIRE_THROWS_AS(assembler.assemble_stiffness(&wrong_length), ModelError);

  Vector negative = Vector::Ones(model.mesh().num_elements());
  negative(0) = -1.0;
  REQUIRE_THROWS_AS(assembler.assemble_stiffness(&negative), ModelError);

  Vector nan_scale = Vector::Ones(model.mesh().num_elements());
  nan_scale(1) = std::numeric_limits<Scalar>::quiet_NaN();
  REQUIRE_THROWS_AS(assembler.assemble_stiffness(&nan_scale), ModelError);
}

TEST_CASE("the reduced stiffness matrix is positive definite once constrained",
          "[solver][verification]") {
  FemModel model = make_small_plate();
  Assembler assembler(model);
  const SparseMatrix k = assembler.assemble_stiffness();
  const SparseMatrix kff = assembler.reduce_free_free(k);

  REQUIRE(kff.rows() == model.dofs().num_free());
  Scalar pivot_ratio = 0.0;
  REQUIRE(is_positive_definite(kff, 1.0e-14, &pivot_ratio));
  REQUIRE(pivot_ratio > 0.0);

  // Every eigenvalue of the dense equivalent is strictly positive.
  const Matrix dense_kff(kff);
  Eigen::SelfAdjointEigenSolver<Matrix> es(dense_kff);
  REQUIRE(es.eigenvalues().minCoeff() > 0.0);
}

TEST_CASE("DOF partitioning round-trips vectors", "[dofs]") {
  FemModel model = make_small_plate();
  const DofManager& dofs = model.dofs();
  REQUIRE(dofs.num_free() + dofs.num_constrained() == dofs.num_dofs());

  Vector full = Vector::Zero(dofs.num_dofs());
  for (Index d = 0; d < dofs.num_dofs(); ++d) full(d) = 0.001 * (d + 1);
  const Vector reduced = dofs.restrict_to_free(full);
  const Vector expanded = dofs.expand(reduced);
  for (Index d : dofs.free_dofs()) REQUIRE(expanded(d) == Approx(full(d)));
  for (Index d : dofs.constrained_dofs()) {
    REQUIRE(expanded(d) == Approx(dofs.prescribed_value(d)));
  }

  REQUIRE_THROWS_AS(dofs.expand(Vector::Zero(3)), ModelError);
  REQUIRE_THROWS_AS(dofs.restrict_to_free(Vector::Zero(3)), ModelError);
  REQUIRE_THROWS_AS(dofs.dof(-1, 0), ModelError);
  REQUIRE_THROWS_AS(dofs.dof(0, 2), ModelError);
}

TEST_CASE("all linear solvers agree on a small model", "[solver][verification]") {
  FemModel model = make_small_plate(6, 4);
  Assembler assembler(model);

  const std::vector<LinearSolverType> types = {
      LinearSolverType::DenseLu,          LinearSolverType::SimplicialLdlt,
      LinearSolverType::SimplicialLlt,    LinearSolverType::SparseLu,
      LinearSolverType::ConjugateGradient};

  Vector reference;
  Scalar reference_compliance = 0.0;
  for (LinearSolverType type : types) {
    StaticAnalysisOptions options;
    options.linear.type = type;
    options.linear.iterative_tolerance = 1.0e-14;
    options.linear.residual_tolerance = 1.0e-7;
    StaticAnalysis analysis(model, assembler, options);
    const std::vector<StaticSolution> solutions = analysis.solve_all();
    const StaticSolution& sol = solutions.front();

    if (reference.size() == 0) {
      reference = sol.displacement;
      reference_compliance = sol.compliance;
      continue;
    }
    const Scalar rel = (sol.displacement - reference).cwiseAbs().maxCoeff() /
                       reference.cwiseAbs().maxCoeff();
    REQUIRE(rel == Approx(0.0).margin(1.0e-8));
    REQUIRE(sol.compliance == Approx(reference_compliance).epsilon(1.0e-9));
  }
}

TEST_CASE("the dense solver refuses large systems", "[solver][diagnostics]") {
  LinearSolverOptions options;
  options.type = LinearSolverType::DenseLu;
  auto solver = make_linear_solver(options);
  SparseMatrix big(5000, 5000);
  big.setIdentity();
  REQUIRE_THROWS_AS(solver->factorize(big), SolverError);
}

TEST_CASE("a singular system is reported, not silently solved",
          "[solver][diagnostics]") {
  LinearSolverOptions options;
  auto solver = make_linear_solver(options);

  // Zero matrix: all pivots vanish.
  SparseMatrix zero(4, 4);
  zero.setZero();
  zero.makeCompressed();
  REQUIRE_THROWS_AS(solver->factorize(zero), SolverError);

  // Indefinite matrix: LDLT finds a non-positive pivot.
  TripletList triplets = {Triplet(0, 0, 1.0), Triplet(1, 1, -1.0)};
  SparseMatrix indefinite(2, 2);
  indefinite.setFromTriplets(triplets.begin(), triplets.end());
  REQUIRE_THROWS_AS(solver->factorize(indefinite), SolverError);

  // Non-finite entries.
  TripletList bad = {Triplet(0, 0, std::numeric_limits<Scalar>::quiet_NaN()),
                     Triplet(1, 1, 1.0)};
  SparseMatrix nan_matrix(2, 2);
  nan_matrix.setFromTriplets(bad.begin(), bad.end());
  REQUIRE_THROWS_AS(solver->factorize(nan_matrix), SolverError);

  // Empty system.
  SparseMatrix empty(0, 0);
  REQUIRE_THROWS_AS(solver->factorize(empty), SolverError);

  REQUIRE_THROWS_AS(parse_linear_solver_type("no_such_solver"), ConfigError);
}

TEST_CASE("reactions balance the applied load exactly",
          "[reactions][verification]") {
  CantileverCase c;
  c.poisson = 0.3;
  FemModel model = make_cantilever(c, 24, 6);
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
  const std::vector<StaticSolution> solutions = analysis.solve_all();
  const StaticSolution& sol = solutions.front();

  // Force balance.
  REQUIRE(sol.equilibrium.applied_force.y() == Approx(c.tip_load));
  REQUIRE(sol.equilibrium.reaction_force.y() == Approx(-c.tip_load).epsilon(1.0e-10));
  REQUIRE(sol.equilibrium.reaction_force.x() == Approx(0.0).margin(1.0e-8));
  REQUIRE(sol.equilibrium.relative_force_error == Approx(0.0).margin(1.0e-10));

  // Moment balance about the origin: the root moment must equal P * L.
  REQUIRE(sol.equilibrium.applied_moment.z() ==
          Approx(c.tip_load * c.length).epsilon(1.0e-10));
  REQUIRE(sol.equilibrium.applied_moment.x() == 0.0);
  REQUIRE(sol.equilibrium.applied_moment.y() == 0.0);
  REQUIRE(sol.equilibrium.relative_moment_error == Approx(0.0).margin(1.0e-10));

  // Reactions vanish at free DOFs.
  for (Index d : model.dofs().free_dofs()) {
    REQUIRE(sol.reactions(d) == Approx(0.0));
  }

  // Compliance equals twice the strain energy for homogeneous Dirichlet data.
  REQUIRE(sol.compliance == Approx(2.0 * sol.strain_energy).epsilon(1.0e-10));
  REQUIRE(sol.scaled_residual < 1.0e-10);
}

TEST_CASE("non-zero prescribed displacements are handled by condensation",
          "[reactions][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 8;
  spec.ny = 4;
  spec.lx = 0.8;
  spec.ly = 0.4;

  const Scalar e_mod = 70.0e9;
  const Scalar thickness = 0.01;
  const Scalar stretch = 1.0e-4;

  FemModel model(make_structured_quad_mesh(spec),
                 IsotropicMaterial(e_mod, 0.0, 2700.0, "uniaxial"), thickness,
                 StressState::PlaneStress, IntegrationOptions());

  DisplacementConstraint left;
  left.region.name = "left";
  Selector lbox;
  lbox.kind = SelectorKind::Box;
  lbox.xmax = 0.0;
  left.region.members.push_back(lbox);
  left.fix_x = true;
  left.fix_y = true;
  model.constraints().push_back(left);

  DisplacementConstraint right;
  right.region.name = "right";
  Selector rbox;
  rbox.kind = SelectorKind::Box;
  rbox.xmin = spec.lx;
  right.region.members.push_back(rbox);
  right.fix_x = true;
  right.value_x = stretch;
  model.constraints().push_back(right);

  LoadCaseSpec load;
  load.name = "prescribed_stretch";
  PointLoadSpec zero;
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector3(spec.lx * 0.5, spec.ly * 0.5, 0.0);
  zero.region.members.push_back(nearest);
  zero.force = Vector3::Zero();
  load.point_loads.push_back(zero);
  model.load_case_specs().push_back(load);
  model.finalize();

  REQUIRE(model.dofs().has_nonzero_prescribed());

  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler, StaticAnalysisOptions());
  const std::vector<StaticSolution> solutions = analysis.solve_all();
  const StaticSolution& sol = solutions.front();

  // With nu = 0 and a fully clamped left edge the exact answer is a uniform
  // uniaxial stretch, so the reaction resultant is E * A * strain.
  const Scalar strain = stretch / spec.lx;
  const Scalar expected = e_mod * (spec.ly * thickness) * strain;
  Scalar right_reaction = 0.0;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    if (std::abs(model.mesh().node(n).x() - spec.lx) < 1.0e-12) {
      right_reaction += sol.reactions(n * 2 + 0);
    }
  }
  REQUIRE(right_reaction == Approx(expected).epsilon(1.0e-9));

  // Net force balance still holds: the two edge reactions cancel.
  REQUIRE(sol.equilibrium.reaction_force.x() == Approx(0.0).margin(1.0e-6));

  // Applied force is zero, so compliance is zero while strain energy is not.
  REQUIRE(sol.compliance == Approx(0.0).margin(1.0e-20));
  REQUIRE(sol.strain_energy > 0.0);
  REQUIRE(sol.strain_energy ==
          Approx(0.5 * expected * stretch).epsilon(1.0e-9));
}

TEST_CASE("under-constrained models are detected before solving",
          "[diagnostics][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 3;

  const auto build_with = [&](const std::vector<DisplacementConstraint>& bcs) {
    FemModel model(make_structured_quad_mesh(spec), default_material(), 0.01,
                   StressState::PlaneStress, IntegrationOptions());
    model.constraints() = bcs;
    LoadCaseSpec load;
    load.name = "corner";
    PointLoadSpec p;
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector3(1.0, 0.0, 0.0);
    p.region.members.push_back(nearest);
    p.force = Vector3(0.0, -10.0, 0.0);
    load.point_loads.push_back(p);
    model.load_case_specs().push_back(load);
    model.finalize();
    return model;
  };

  SECTION("no constraints at all") {
    FemModel model = build_with({});
    const ModelDiagnostics diag = diagnose_model(model);
    REQUIRE_FALSE(diag.well_posed());
    REQUIRE(diag.components.size() == 1);
    REQUIRE(diag.components.front().rigid_null_dimension == 3);
    REQUIRE_THROWS_AS(require_well_posed(model), ModelError);
    // The static analysis refuses to run.
    Assembler assembler(model);
    REQUIRE_THROWS_AS(StaticAnalysis(model, assembler, StaticAnalysisOptions()),
                      ModelError);
  }

  SECTION("a single pinned node leaves the rotation free") {
    DisplacementConstraint bc;
    bc.region.name = "one_node";
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector3(0.0, 0.0, 0.0);
    bc.region.members.push_back(nearest);
    bc.fix_x = true;
    bc.fix_y = true;
    FemModel model = build_with({bc});
    const ModelDiagnostics diag = diagnose_model(model);
    REQUIRE_FALSE(diag.well_posed());
    REQUIRE(diag.components.front().rigid_null_dimension == 1);
  }

  SECTION("two nodes fixed in y only leave the x translation free") {
    DisplacementConstraint bc;
    bc.region.name = "bottom";
    Selector box;
    box.kind = SelectorKind::Box;
    box.ymax = 0.0;
    bc.region.members.push_back(box);
    bc.fix_y = true;
    FemModel model = build_with({bc});
    const ModelDiagnostics diag = diagnose_model(model);
    REQUIRE_FALSE(diag.well_posed());
    REQUIRE(diag.components.front().rigid_null_dimension == 1);
  }

  SECTION("a properly clamped edge is well posed") {
    DisplacementConstraint bc;
    bc.region.name = "left";
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.0;
    bc.region.members.push_back(box);
    bc.fix_x = true;
    bc.fix_y = true;
    FemModel model = build_with({bc});
    const ModelDiagnostics diag = diagnose_model(model);
    REQUIRE(diag.well_posed());
    REQUIRE(diag.components.front().rigid_null_dimension == 0);
    REQUIRE_NOTHROW(require_well_posed(model));
  }

  SECTION("three well-placed single constraints suffice") {
    DisplacementConstraint pin;
    pin.region.name = "pin";
    Selector n1;
    n1.kind = SelectorKind::NearestNode;
    n1.point = Vector3(0.0, 0.0, 0.0);
    pin.region.members.push_back(n1);
    pin.fix_x = true;
    pin.fix_y = true;

    DisplacementConstraint roller;
    roller.region.name = "roller";
    Selector n2;
    n2.kind = SelectorKind::NearestNode;
    n2.point = Vector3(1.0, 0.0, 0.0);
    roller.region.members.push_back(n2);
    roller.fix_y = true;

    FemModel model = build_with({pin, roller});
    const ModelDiagnostics diag = diagnose_model(model);
    REQUIRE(diag.well_posed());
  }
}

TEST_CASE("a floating element group is reported separately",
          "[diagnostics][verification]") {
  // Two 1x1 blocks separated by a gap, constrained only on the left block.
  Matrix coords(2, 8);
  coords << 0, 1, 1, 0, 3, 4, 4, 3,
            0, 0, 1, 1, 0, 0, 1, 1;
  const std::vector<Index> connectivity = {0, 1, 2, 3, 4, 5, 6, 7};

  FemModel model(Mesh(coords, connectivity, ElementType::Quad4), default_material(),
                 0.01, StressState::PlaneStress, IntegrationOptions());
  DisplacementConstraint bc;
  bc.region.name = "left_block";
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  bc.region.members.push_back(box);
  bc.fix_x = true;
  bc.fix_y = true;
  model.constraints().push_back(bc);

  LoadCaseSpec load;
  load.name = "pull";
  PointLoadSpec p;
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector3(4.0, 1.0, 0.0);
  p.region.members.push_back(nearest);
  p.force = Vector3(100.0, 0.0, 0.0);
  load.point_loads.push_back(p);
  model.load_case_specs().push_back(load);
  model.finalize();

  const ModelDiagnostics diag = diagnose_model(model);
  REQUIRE(diag.components.size() == 2);
  REQUIRE_FALSE(diag.well_posed());
  // Exactly one group is unconstrained, and the message says so.
  int floating = 0;
  for (const MeshComponent& c : diag.components) {
    if (c.prescribed_dofs == 0) ++floating;
  }
  REQUIRE(floating == 1);
  REQUIRE(diag.problems.size() == 1);
  REQUIRE(diag.problems.front().find("floating region") != std::string::npos);
}
