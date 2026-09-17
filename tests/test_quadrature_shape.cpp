/// \file test_quadrature_shape.cpp
/// \brief Quadrature exactness and shape-function identities.
#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quad4.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace sparlab;
using Catch::Approx;

TEST_CASE("Gauss-Legendre line rules integrate polynomials exactly", "[quadrature]") {
  // An n-point rule is exact up to degree 2n-1.
  for (int n = 1; n <= 4; ++n) {
    const auto& rule = gauss_legendre_line(n);
    REQUIRE(static_cast<int>(rule.size()) == n);

    for (int degree = 0; degree <= 2 * n - 1; ++degree) {
      Scalar numeric = 0.0;
      for (const auto& p : rule) numeric += p.weight * std::pow(p.xi, degree);
      // exact integral of x^d over [-1,1]
      const Scalar exact = (degree % 2 == 1) ? 0.0 : 2.0 / (degree + 1);
      REQUIRE(numeric == Approx(exact).margin(1.0e-14));
    }
  }
}

TEST_CASE("Gauss-Legendre square rules are tensor products with total weight 4",
          "[quadrature]") {
  for (int n = 1; n <= 4; ++n) {
    const auto& rule = gauss_legendre_square(n);
    REQUIRE(static_cast<int>(rule.size()) == n * n);
    Scalar total = 0.0;
    for (const auto& p : rule) total += p.weight;
    REQUIRE(total == Approx(4.0).margin(1.0e-14));
  }

  // 2x2 integrates the bilinear monomial xi^3 eta^3 exactly (degree 3 each).
  const auto& rule = gauss_legendre_square(2);
  Scalar value = 0.0;
  for (const auto& p : rule) {
    value += p.weight * std::pow(p.xi, 3) * std::pow(p.eta, 3);
  }
  REQUIRE(value == Approx(0.0).margin(1.0e-14));

  Scalar quadratic = 0.0;
  for (const auto& p : rule) quadratic += p.weight * p.xi * p.xi * p.eta * p.eta;
  REQUIRE(quadratic == Approx(4.0 / 9.0).margin(1.0e-14));
}

TEST_CASE("unsupported quadrature orders are rejected", "[quadrature][diagnostics]") {
  REQUIRE_THROWS_AS(gauss_legendre_line(0), ConfigError);
  REQUIRE_THROWS_AS(gauss_legendre_line(5), ConfigError);
  REQUIRE_THROWS_AS(gauss_legendre_square(7), ConfigError);
}

TEST_CASE("Q4 shape functions satisfy the partition of unity and nodal delta",
          "[shape]") {
  const Scalar nodes[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  for (int a = 0; a < 4; ++a) {
    const Eigen::Vector4d n = quad4_shape_functions(nodes[a][0], nodes[a][1]);
    for (int b = 0; b < 4; ++b) {
      REQUIRE(n(b) == Approx(a == b ? 1.0 : 0.0).margin(1.0e-15));
    }
  }

  // Partition of unity and zero gradient sum at arbitrary points.
  const Scalar samples[5][2] = {{0, 0}, {0.3, -0.7}, {-0.9, 0.2}, {0.5, 0.5}, {-0.25, -0.4}};
  for (const auto& s : samples) {
    const Eigen::Vector4d n = quad4_shape_functions(s[0], s[1]);
    REQUIRE(n.sum() == Approx(1.0).margin(1.0e-15));
    const Eigen::Matrix<Scalar, 4, 2> g = quad4_shape_gradients_natural(s[0], s[1]);
    REQUIRE(g.col(0).sum() == Approx(0.0).margin(1.0e-15));
    REQUIRE(g.col(1).sum() == Approx(0.0).margin(1.0e-15));
  }
}

TEST_CASE("Q4 shape gradients match central differences", "[shape]") {
  const Scalar h = 1.0e-6;
  const Scalar xi = 0.23;
  const Scalar eta = -0.41;
  const Eigen::Matrix<Scalar, 4, 2> analytic = quad4_shape_gradients_natural(xi, eta);
  for (int a = 0; a < 4; ++a) {
    const Scalar dxi = (quad4_shape_functions(xi + h, eta)(a) -
                        quad4_shape_functions(xi - h, eta)(a)) /
                       (2.0 * h);
    const Scalar deta = (quad4_shape_functions(xi, eta + h)(a) -
                         quad4_shape_functions(xi, eta - h)(a)) /
                        (2.0 * h);
    REQUIRE(analytic(a, 0) == Approx(dxi).margin(1.0e-9));
    REQUIRE(analytic(a, 1) == Approx(deta).margin(1.0e-9));
  }
}
