/// \file test_element.cpp
/// \brief Element-level verification: symmetry, rigid-body modes, exact
///        integration, mass conservation and consistent edge loads (Q4).
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quad4.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

TEST_CASE("element stiffness is symmetric", "[element][verification]") {
  Quad4Element element;
  const IsotropicMaterial material = default_material();
  const Matrix3 d = material.plane_stress_matrix();
  const IntegrationOptions opts;

  for (const auto& coords : {unit_square_coords(2.0, 0.5), distorted_quad_coords()}) {
    const Matrix ke = element.stiffness(coords, d, 0.01, opts);
    REQUIRE(ke.rows() == 8);
    REQUIRE(ke.cols() == 8);
    const Scalar asymmetry = (ke - ke.transpose()).cwiseAbs().maxCoeff();
    REQUIRE(asymmetry / ke.cwiseAbs().maxCoeff() == Approx(0.0).margin(1.0e-15));
  }
}

TEST_CASE("element stiffness annihilates the three rigid-body modes",
          "[element][verification]") {
  Quad4Element element;
  const Matrix3 d = default_material().plane_stress_matrix();
  const IntegrationOptions opts;

  for (const auto& coords : {unit_square_coords(1.5, 0.8), distorted_quad_coords()}) {
    const Matrix ke = element.stiffness(coords, d, 0.02, opts);

    // Modes about the element centroid.
    Vector2 centroid = Vector2::Zero();
    for (int a = 0; a < 4; ++a) centroid += coords.col(a);
    centroid /= 4.0;

    Matrix modes = Matrix::Zero(8, 3);
    for (int a = 0; a < 4; ++a) {
      const Vector2 r = coords.col(a) - centroid;
      modes(2 * a + 0, 0) = 1.0;
      modes(2 * a + 1, 1) = 1.0;
      modes(2 * a + 0, 2) = -r.y();
      modes(2 * a + 1, 2) = r.x();
    }

    const Scalar scale = ke.cwiseAbs().maxCoeff();
    for (int m = 0; m < 3; ++m) {
      const Vector residual = ke * modes.col(m);
      REQUIRE(residual.cwiseAbs().maxCoeff() / (scale * modes.col(m).norm()) ==
              Approx(0.0).margin(1.0e-14));
      // Zero strain energy as well.
      REQUIRE(modes.col(m).dot(ke * modes.col(m)) /
                  (scale * modes.col(m).squaredNorm()) ==
              Approx(0.0).margin(1.0e-14));
    }

    // The null space is exactly three-dimensional: the other five eigenvalues
    // must be strictly positive.
    Eigen::SelfAdjointEigenSolver<Matrix> es(ke);
    const Vector eig = es.eigenvalues();
    int near_zero = 0;
    for (Eigen::Index i = 0; i < eig.size(); ++i) {
      if (std::abs(eig(i)) < 1.0e-9 * scale) ++near_zero;
    }
    REQUIRE(near_zero == 3);
    REQUIRE(eig(eig.size() - 1) > 0.0);
  }
}

TEST_CASE("element stiffness reproduces the closed form for a unit square",
          "[element][verification]") {
  // For a rectangle the 2x2 rule is exact, so the result must be independent of
  // higher quadrature orders.
  Quad4Element element;
  const Matrix3 d = default_material().plane_stress_matrix();
  const auto coords = unit_square_coords(0.3, 0.2);

  IntegrationOptions two;
  two.stiffness_points = 2;
  IntegrationOptions three;
  three.stiffness_points = 3;
  IntegrationOptions four;
  four.stiffness_points = 4;

  const Matrix k2 = element.stiffness(coords, d, 0.01, two);
  const Matrix k3 = element.stiffness(coords, d, 0.01, three);
  const Matrix k4 = element.stiffness(coords, d, 0.01, four);

  const Scalar scale = k2.cwiseAbs().maxCoeff();
  REQUIRE((k3 - k2).cwiseAbs().maxCoeff() / scale == Approx(0.0).margin(1.0e-13));
  REQUIRE((k4 - k2).cwiseAbs().maxCoeff() / scale == Approx(0.0).margin(1.0e-13));
}

TEST_CASE("element stiffness scales linearly with thickness and modulus",
          "[element]") {
  Quad4Element element;
  const IntegrationOptions opts;
  const auto coords = distorted_quad_coords();

  const Matrix3 d1 = IsotropicMaterial(1.0e9, 0.3, 1.0).plane_stress_matrix();
  const Matrix3 d2 = IsotropicMaterial(3.0e9, 0.3, 1.0).plane_stress_matrix();
  const Matrix k1 = element.stiffness(coords, d1, 0.01, opts);
  const Matrix k2 = element.stiffness(coords, d2, 0.01, opts);
  const Matrix k3 = element.stiffness(coords, d1, 0.03, opts);

  REQUIRE((k2 - 3.0 * k1).cwiseAbs().maxCoeff() / k2.cwiseAbs().maxCoeff() ==
          Approx(0.0).margin(1.0e-14));
  REQUIRE((k3 - 3.0 * k1).cwiseAbs().maxCoeff() / k3.cwiseAbs().maxCoeff() ==
          Approx(0.0).margin(1.0e-14));
}

TEST_CASE("strain operator reproduces a linear displacement field exactly",
          "[element][verification]") {
  Quad4Element element;
  const auto coords = distorted_quad_coords();

  // u = a + G x with constant strain.
  const Vector2 a(1.0e-3, -2.0e-3);
  Matrix2 g;
  g << 2.0e-3, 5.0e-4, 5.0e-4, -1.0e-3;

  Vector ue(8);
  for (int n = 0; n < 4; ++n) {
    const Vector2 u = a + g * coords.col(n);
    ue(2 * n + 0) = u.x();
    ue(2 * n + 1) = u.y();
  }
  const Vector3 expected(g(0, 0), g(1, 1), g(0, 1) + g(1, 0));

  // Exact at every parametric point, including well away from Gauss points.
  const Scalar pts[5][2] = {{0, 0}, {-0.9, 0.7}, {0.6, -0.3}, {1.0, 1.0}, {-1.0, -1.0}};
  for (const auto& p : pts) {
    NaturalPoint np;
    np.xi = p[0];
    np.eta = p[1];
    const StrainOperator op = element.strain_operator(coords, np);
    REQUIRE(op.detJ > 0.0);
    const Vector strain = op.b * ue;
    REQUIRE((strain - expected).cwiseAbs().maxCoeff() ==
            Approx(0.0).margin(1.0e-17));
  }
}

TEST_CASE("Jacobian determinant integrates to the element area", "[element]") {
  Quad4Element element;
  const auto coords = distorted_quad_coords();
  IntegrationOptions opts;
  opts.stiffness_points = 3;

  Scalar area = 0.0;
  for (const auto& gp : gauss_legendre_square(opts.stiffness_points)) {
    NaturalPoint p;
    p.xi = gp.xi;
    p.eta = gp.eta;
    area += gp.weight * element.strain_operator(coords, p).detJ;
  }

  // Shoelace area of the same quadrilateral.
  Scalar twice = 0.0;
  for (int a = 0; a < 4; ++a) {
    const Vector2 p = coords.col(a);
    const Vector2 q = coords.col((a + 1) % 4);
    twice += p.x() * q.y() - q.x() * p.y();
  }
  REQUIRE(area == Approx(0.5 * twice).epsilon(1.0e-12));
}

TEST_CASE("inverted and degenerate elements are rejected",
          "[element][diagnostics]") {
  Quad4Element element;
  const Matrix3 d = default_material().plane_stress_matrix();
  const IntegrationOptions opts;

  Eigen::Matrix<Scalar, 2, Eigen::Dynamic> clockwise(2, 4);
  clockwise << 0, 0, 1, 1,
               0, 1, 1, 0;
  REQUIRE_THROWS_AS(element.stiffness(clockwise, d, 0.01, opts), MeshError);

  Eigen::Matrix<Scalar, 2, Eigen::Dynamic> collapsed(2, 4);
  collapsed << 0, 1, 1, 0,
               0, 0, 0, 0;
  REQUIRE_THROWS_AS(element.stiffness(collapsed, d, 0.01, opts), MeshError);

  Eigen::Matrix<Scalar, 2, Eigen::Dynamic> wrong_size(2, 3);
  wrong_size << 0, 1, 1,
                0, 0, 1;
  NaturalPoint p;
  REQUIRE_THROWS_AS(element.strain_operator(wrong_size, p), MeshError);

  REQUIRE_THROWS_AS(element.stiffness(unit_square_coords(), d, 0.0, opts), ConfigError);
  REQUIRE_THROWS_AS(element.stiffness(unit_square_coords(), d, -1.0, opts), ConfigError);
}

TEST_CASE("consistent mass matrix is symmetric and conserves mass",
          "[element][modal][verification]") {
  Quad4Element element;
  const Scalar density = 2700.0;
  const Scalar thickness = 0.02;
  const IntegrationOptions opts;

  for (const auto& coords : {unit_square_coords(1.3, 0.7), distorted_quad_coords()}) {
    const Matrix me = element.consistent_mass(coords, density, thickness, opts);
    REQUIRE((me - me.transpose()).cwiseAbs().maxCoeff() / me.cwiseAbs().maxCoeff() ==
            Approx(0.0).margin(1.0e-15));

    // Positive definite: every mode has positive kinetic energy.
    Eigen::SelfAdjointEigenSolver<Matrix> es(me);
    REQUIRE(es.eigenvalues().minCoeff() > 0.0);

    // Sum of all entries = 2 * rho * t * A (one factor per translation).
    Scalar twice = 0.0;
    for (int a = 0; a < 4; ++a) {
      const Vector2 p = coords.col(a);
      const Vector2 q = coords.col((a + 1) % 4);
      twice += p.x() * q.y() - q.x() * p.y();
    }
    const Scalar mass = density * thickness * 0.5 * twice;
    REQUIRE(me.sum() == Approx(2.0 * mass).epsilon(1.0e-12));

    // A unit rigid translation carries exactly the element mass.
    Vector tx = Vector::Zero(8);
    for (int a = 0; a < 4; ++a) tx(2 * a) = 1.0;
    REQUIRE(tx.dot(me * tx) == Approx(mass).epsilon(1.0e-12));
  }
}

TEST_CASE("row-sum lumping conserves the element mass", "[element][modal]") {
  Quad4Element element;
  const auto coords = distorted_quad_coords();
  const Matrix me = element.consistent_mass(coords, 2700.0, 0.02, IntegrationOptions());

  Scalar lumped_total = 0.0;
  for (int i = 0; i < 8; ++i) lumped_total += me.row(i).sum();
  REQUIRE(lumped_total == Approx(me.sum()).epsilon(1.0e-14));
  for (int i = 0; i < 8; ++i) REQUIRE(me.row(i).sum() > 0.0);
}

TEST_CASE("edge traction produces the exact consistent nodal load",
          "[element][loads][verification]") {
  Quad4Element element;
  const auto coords = unit_square_coords(2.0, 0.5);
  const Scalar thickness = 0.01;
  const Vector3 traction(0.0, -1.0e5, 0.0);  // Pa, downward
  IntegrationOptions opts;

  // Local edge 2 runs from node 2 (2, 0.5) to node 3 (0, 0.5): length 2 m.
  const Vector fe = element.boundary_traction(coords, 2, traction, thickness, opts);
  const Scalar expected_total = traction.y() * 2.0 * thickness;
  REQUIRE(fe.sum() == Approx(expected_total));

  // A constant traction on a straight edge splits evenly between its two nodes.
  REQUIRE(fe(2 * 2 + 1) == Approx(0.5 * expected_total));
  REQUIRE(fe(2 * 3 + 1) == Approx(0.5 * expected_total));
  // No force anywhere else.
  REQUIRE(fe(0) == Approx(0.0));
  REQUIRE(fe(1) == Approx(0.0));
  REQUIRE(fe(2 * 1 + 1) == Approx(0.0));

  REQUIRE_THROWS_AS(element.boundary_traction(coords, 4, traction, thickness, opts),
                    MeshError);
  REQUIRE_THROWS_AS(element.boundary_traction(coords, -1, traction, thickness, opts),
                    MeshError);
}

TEST_CASE("element factory returns the requested topology", "[element]") {
  const auto element = make_element(ElementType::Quad4);
  REQUIRE(element != nullptr);
  REQUIRE(element->type() == ElementType::Quad4);
  REQUIRE(element->num_nodes() == 4);
  REQUIRE(element->num_dofs() == 8);
  REQUIRE(element->dim() == 2);
  REQUIRE(element->num_faces() == 4);
  REQUIRE(element->face_nodes(0)[0] == 0);
  REQUIRE(element->face_nodes(3)[1] == 0);
  REQUIRE_THROWS_AS(element->face_nodes(4), MeshError);
}
