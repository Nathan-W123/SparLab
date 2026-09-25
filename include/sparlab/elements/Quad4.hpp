/// \file Quad4.hpp
/// \brief Four-node bilinear isoparametric quadrilateral (plane element).
///
/// Shape functions on the reference square \f$(\xi,\eta)\in[-1,1]^2\f$ with the
/// counter-clockwise node ordering of StructuredMesh.hpp:
/// \f[
///   N_1 = \tfrac{1}{4}(1-\xi)(1-\eta), \quad
///   N_2 = \tfrac{1}{4}(1+\xi)(1-\eta), \quad
///   N_3 = \tfrac{1}{4}(1+\xi)(1+\eta), \quad
///   N_4 = \tfrac{1}{4}(1-\xi)(1+\eta).
/// \f]
/// Isoparametric mapping \f$ x = \sum_a N_a x_a \f$, Jacobian
/// \f$ J_{ij} = \partial x_i / \partial \xi_j \f$, and
/// \f$ \partial N_a/\partial x = J^{-T} \partial N_a/\partial \xi \f$.
///
/// The element DOF ordering is node-major:
/// \f$ \{u_{x1}, u_{y1}, u_{x2}, u_{y2}, u_{x3}, u_{y3}, u_{x4}, u_{y4}\} \f$.
#pragma once

#include "sparlab/elements/Element.hpp"

namespace sparlab {

/// Bilinear shape functions of the Q4 element (free functions so they can be
/// unit-tested and reused by post-processing without constructing an element).
/// \{
Eigen::Vector4d quad4_shape_functions(Scalar xi, Scalar eta);
Eigen::Matrix<Scalar, 4, 2> quad4_shape_gradients_natural(Scalar xi, Scalar eta);
/// \}

class Quad4Element final : public Element {
 public:
  ElementType type() const override { return ElementType::Quad4; }
  int dim() const override { return 2; }
  int num_nodes() const override { return 4; }
  int num_faces() const override { return 4; }

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
