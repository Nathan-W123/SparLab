/// \file Element.hpp
/// \brief Abstract element interface.
///
/// Every element implementation supplies the four kernels the rest of the
/// library needs: a stiffness matrix, a consistent mass matrix, a
/// strain-displacement operator at a requested parametric location, and a
/// consistent nodal load vector for an edge traction. Adding a new topology
/// (Q8, T3, T6, ...) means adding one subclass and one entry to `make_element`;
/// no other component changes.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <memory>

namespace sparlab {

/// Element-local coordinate pair on the reference domain.
struct NaturalPoint {
  Scalar xi = 0.0;
  Scalar eta = 0.0;
};

/// Quadrature orders used by an element's kernels.
struct IntegrationOptions {
  int stiffness_points = 2;  ///< points per direction for K_e
  int mass_points = 3;       ///< points per direction for M_e
  int edge_points = 2;       ///< points along an edge for traction loads
};

/// Strain-displacement data evaluated at one parametric point.
struct StrainOperator {
  Eigen::Matrix<Scalar, kVoigt, Eigen::Dynamic> b;  ///< 3 x (dofs) operator [1/m]
  Scalar detJ = 0.0;                                ///< Jacobian determinant [m^2]
};

/// Abstract 2-D continuum element.
class Element {
 public:
  virtual ~Element() = default;

  virtual ElementType type() const = 0;
  virtual int num_nodes() const = 0;

  /// Degrees of freedom carried by the element (num_nodes * kDofsPerNode).
  int num_dofs() const { return num_nodes() * kDofsPerNode; }

  /// Element stiffness matrix
  /// \f$ K_e = \int_{\Omega_e} t\, B^T D B \, d\Omega \f$ [N/m].
  /// \param coords 2 x num_nodes nodal coordinates [m].
  /// \param d constitutive matrix [Pa].
  /// \param thickness out-of-plane thickness [m].
  virtual Matrix stiffness(const Eigen::Matrix<Scalar, 2, Eigen::Dynamic>& coords,
                           const Matrix3& d, Scalar thickness,
                           const IntegrationOptions& opts) const = 0;

  /// Consistent element mass matrix
  /// \f$ M_e = \int_{\Omega_e} \rho\, t\, N^T N \, d\Omega \f$ [kg].
  virtual Matrix consistent_mass(const Eigen::Matrix<Scalar, 2, Eigen::Dynamic>& coords,
                                 Scalar density, Scalar thickness,
                                 const IntegrationOptions& opts) const = 0;

  /// Strain-displacement operator and Jacobian determinant at a natural point.
  /// \throws MeshError when detJ <= 0 (inverted or degenerate element).
  virtual StrainOperator strain_operator(
      const Eigen::Matrix<Scalar, 2, Eigen::Dynamic>& coords,
      const NaturalPoint& point) const = 0;

  /// Shape function values at a natural point (size num_nodes).
  virtual Vector shape_functions(const NaturalPoint& point) const = 0;

  /// Parametric locations of the stiffness quadrature points, in the order used
  /// by `stiffness`. Needed for stress recovery at Gauss points.
  virtual std::vector<NaturalPoint> stress_evaluation_points(
      const IntegrationOptions& opts) const = 0;

  /// Consistent nodal forces for a constant traction on local edge `local_edge`.
  /// \f$ f_e = \int_{\Gamma_e} t\, N^T \bar{t} \, ds \f$ [N].
  /// \param traction traction vector [Pa].
  /// \return vector of length num_dofs [N].
  virtual Vector edge_traction(const Eigen::Matrix<Scalar, 2, Eigen::Dynamic>& coords,
                               int local_edge, const Vector2& traction, Scalar thickness,
                               const IntegrationOptions& opts) const = 0;

  /// Local node indices bounding local edge `local_edge`.
  virtual std::array<int, 2> edge_nodes(int local_edge) const = 0;

  /// Number of edges of the element.
  virtual int num_edges() const = 0;
};

/// Factory for the supported element topologies.
/// \throws ConfigError for an unsupported type.
std::unique_ptr<Element> make_element(ElementType type);

}  // namespace sparlab
