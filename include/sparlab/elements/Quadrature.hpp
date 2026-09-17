/// \file Quadrature.hpp
/// \brief Gauss-Legendre quadrature rules on the reference square [-1,1]^2 and
///        on the reference line [-1,1].
///
/// A rule with `n` points per direction integrates polynomials of degree
/// \f$2n-1\f$ exactly. SparLab defaults to:
///   * stiffness:  2 x 2 (exact for the bilinear Q4 on parallelograms),
///   * mass:       3 x 3 (exact for the biquadratic N^T N product on general
///                 quadrilaterals, where the Jacobian is not constant),
///   * edge loads: 2 points (exact for linear traction x linear shape function).
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

/// One quadrature point on the reference line.
struct QuadraturePoint1D {
  Scalar xi = 0.0;
  Scalar weight = 0.0;
};

/// Tensor-product Gauss-Legendre rule with `points_per_direction` points per
/// direction (supported: 1, 2, 3, 4).
/// \throws ConfigError for unsupported orders.
const std::vector<QuadraturePoint2D>& gauss_legendre_square(int points_per_direction);

/// One-dimensional Gauss-Legendre rule (supported: 1, 2, 3, 4).
/// \throws ConfigError for unsupported orders.
const std::vector<QuadraturePoint1D>& gauss_legendre_line(int points);

}  // namespace sparlab
