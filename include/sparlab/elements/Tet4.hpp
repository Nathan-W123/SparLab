/// \file Tet4.hpp
/// \brief Four-node linear tetrahedron (constant-strain solid element).
///
/// Shape functions on the reference tetrahedron with vertices (0,0,0),
/// (1,0,0), (0,1,0), (0,0,1):
/// \f[
///   N_1 = 1 - \xi - \eta - \zeta, \quad N_2 = \xi, \quad N_3 = \eta,
///   \quad N_4 = \zeta .
/// \f]
/// The nodes must be positively oriented,
/// \f$ (x_2 - x_1)\times(x_3 - x_1)\cdot(x_4 - x_1) > 0 \f$ (the VTK and Gmsh
/// convention). The map is affine, \f$\det J = 6V\f$, and the 6 x 12 strain
/// operator is constant, so every kernel is evaluated in closed form and is
/// exact:
///
///   * stiffness   \f$ K_e = V\, B^T D B \f$,
///   * mass        \f$ M_e = \dfrac{\rho V}{20}
///                 \begin{bmatrix} 2&1&1&1\\1&2&1&1\\1&1&2&1\\1&1&1&2 \end{bmatrix}
///                 \otimes I_3 \f$,
///   * face load   a constant traction on a triangular face of area \f$S\f$
///                 puts \f$ S \bar t / 3 \f$ on each of its three nodes.
///
/// Strain is stored in the 6-component Voigt order of the Hex8, with
/// engineering shear strains, and the DOF ordering is node-major. Like the
/// Hex8 the element has no thickness; `thickness` must be 1.
///
/// The linear tetrahedron is what a mesh generator produces for an arbitrary
/// CAD shape, and it is the stiffest element here in bending: accuracy comes
/// from refinement, which the verification studies measure.
#pragma once

#include "sparlab/elements/Element.hpp"

namespace sparlab {

/// Cartesian shape-function gradients of a tetrahedron (row a = grad N_a) and
/// its signed volume [m^3] from 3 x 4 nodal coordinates.
Eigen::Matrix<Scalar, 4, 3> tet4_shape_gradients(const Matrix& coords, Scalar* volume);

class Tet4Element final : public Element {
 public:
  ElementType type() const override { return ElementType::Tet4; }
  int dim() const override { return 3; }
  int num_nodes() const override { return 4; }
  int num_faces() const override { return 4; }
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
