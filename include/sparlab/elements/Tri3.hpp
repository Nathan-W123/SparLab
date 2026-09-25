/// \file Tri3.hpp
/// \brief Three-node linear triangle (constant-strain plane element).
///
/// Shape functions on the reference triangle with vertices (0,0), (1,0),
/// (0,1) and the nodes listed counter-clockwise:
/// \f[
///   N_1 = 1 - \xi - \eta, \qquad N_2 = \xi, \qquad N_3 = \eta .
/// \f]
/// The map is affine, so the Jacobian \f$ J = [x_2 - x_1,\; x_3 - x_1] \f$ is
/// constant, \f$\det J = 2A\f$, and the strain operator \f$B\f$ (3 x 6) is
/// constant over the element: every kernel below is evaluated in closed form
/// and is exact, whatever `IntegrationOptions` requests.
///
///   * stiffness   \f$ K_e = t A\, B^T D B \f$,
///   * mass        \f$ M_e = \dfrac{\rho t A}{12}
///                 \begin{bmatrix} 2&1&1\\1&2&1\\1&1&2 \end{bmatrix}
///                 \otimes I_2 \f$,
///   * edge load   a constant traction on an edge of length \f$L\f$ puts
///                 \f$ t L \bar t / 2 \f$ on each of its two nodes.
///
/// The constant-strain triangle passes the patch test exactly but is stiff in
/// bending: it needs a much finer mesh than the Q4 for the same accuracy,
/// which the verification studies quantify rather than hide.
#pragma once

#include "sparlab/elements/Element.hpp"

namespace sparlab {

/// Cartesian shape-function gradients of a triangle (row a = grad N_a) and
/// its signed area [m^2] from 2 x 3 nodal coordinates.
Eigen::Matrix<Scalar, 3, 2> tri3_shape_gradients(const Matrix& coords, Scalar* area);

class Tri3Element final : public Element {
 public:
  ElementType type() const override { return ElementType::Tri3; }
  int dim() const override { return 2; }
  int num_nodes() const override { return 3; }
  int num_faces() const override { return 3; }
  NaturalPoint reference_centroid() const override;

  Matrix stiffness(const Matrix& coords, const Matrix& d, Scalar thickness,
                   const IntegrationOptions& opts) const override;

  Matrix consistent_mass(const Matrix& coords, Scalar density, Scalar thickness,
                         const IntegrationOptions& opts) const override;

  StrainOperator strain_operator(const Matrix& coords,
                                 const NaturalPoint& point) const override;

  Vector shape_functions(const NaturalPoint& point) const override;

  std::vector<NaturalPoint> stress_evaluation_points(
      const IntegrationOptions& opts) const override;

  std::vector<IntegrationPoint> integration_rule(
      const IntegrationOptions& opts) const override;

  Vector boundary_traction(const Matrix& coords, int local_face, const Vector3& traction,
                           Scalar thickness, const IntegrationOptions& opts) const override;
};

}  // namespace sparlab
