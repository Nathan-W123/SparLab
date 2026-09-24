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

}  // namespace sparlab
