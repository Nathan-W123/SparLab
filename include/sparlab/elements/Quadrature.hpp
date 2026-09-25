/// \file Quadrature.hpp
/// \brief Gauss-Legendre quadrature rules on the reference line [-1,1], the
///        reference square [-1,1]^2 and the reference cube [-1,1]^3.
///
/// A rule with `n` points per direction integrates polynomials of degree
/// \f$2n-1\f$ exactly. SparLab defaults to:
///   * stiffness:  2 per direction (exact for the bilinear Q4 and trilinear
///                 Hex8 on parallelepipeds),
///   * mass:       3 per direction (exact for the N^T N product on general
///                 cells, where the Jacobian is not constant),
///   * face loads: 2 per direction (exact for a constant traction times the
///                 face shape functions).
#pragma once

#include "sparlab/core/Types.hpp"

#include <vector>

namespace sparlab {

/// One quadrature point on the reference square.
struct QuadraturePoint2D {
  Scalar xi = 0.0;
  Scalar eta = 0.0;
  Scalar weight = 0.0;
};

/// One quadrature point on the reference cube.
struct QuadraturePoint3D {
  Scalar xi = 0.0;
  Scalar eta = 0.0;
  Scalar zeta = 0.0;
  Scalar weight = 0.0;
};

/// One quadrature point on the reference line.
struct QuadraturePoint1D {
  Scalar xi = 0.0;
  Scalar weight = 0.0;
};

/// Tensor-product Gauss-Legendre rule on the square with
/// `points_per_direction` points per direction (supported: 1, 2, 3, 4).
/// \throws ConfigError for unsupported orders.
const std::vector<QuadraturePoint2D>& gauss_legendre_square(int points_per_direction);

/// Tensor-product Gauss-Legendre rule on the cube (supported: 1, 2, 3, 4).
/// Points are ordered with xi varying fastest, then eta, then zeta.
/// \throws ConfigError for unsupported orders.
const std::vector<QuadraturePoint3D>& gauss_legendre_cube(int points_per_direction);

/// One-dimensional Gauss-Legendre rule (supported: 1, 2, 3, 4).
/// \throws ConfigError for unsupported orders.
const std::vector<QuadraturePoint1D>& gauss_legendre_line(int points);

/// Symmetric four-point rule on the reference tetrahedron
/// \f$\{\xi,\eta,\zeta \ge 0,\ \xi+\eta+\zeta \le 1\}\f$ (volume 1/6): the
/// points sit at barycentric coordinates \f$(\alpha,\beta,\beta,\beta)\f$ and
/// permutations with \f$\alpha = (5+3\sqrt5)/20\f$, \f$\beta = (5-\sqrt5)/20\f$,
/// each weighted 1/24. Exact for polynomials of degree 2, which is the
/// integrand \f$B^T D B\f$ of a straight-sided quadratic tetrahedron.
const std::vector<QuadraturePoint3D>& tetrahedron_rule_4();

/// Collapsed Gauss rule on the reference tetrahedron with `n` Gauss-Legendre
/// points per direction (n^3 points, supported n = 1..4): the unit cube
/// \f$(u,v,w)\f$ is mapped onto the tetrahedron by
/// \f$\xi = u,\ \eta = v(1-u),\ \zeta = w(1-u)(1-v)\f$, whose Jacobian
/// \f$(1-u)^2(1-v)\f$ is folded into the weights (Stroud's conical product).
/// Exact for polynomials of total degree \f$2n-3\f$, so it takes n >= 2
/// to integrate even a constant exactly.
const std::vector<QuadraturePoint3D>& collapsed_gauss_tetrahedron(int n);

/// Collapsed Gauss rule on the reference triangle
/// \f$\{\xi,\eta \ge 0,\ \xi+\eta \le 1\}\f$ (area 1/2): \f$\xi = u\f$,
/// \f$\eta = v(1-u)\f$ with Jacobian \f$1-u\f$; n^2 points (n = 1..4), exact
/// for polynomials of total degree \f$2n-2\f$.
const std::vector<QuadraturePoint2D>& collapsed_gauss_triangle(int n);

}  // namespace sparlab
