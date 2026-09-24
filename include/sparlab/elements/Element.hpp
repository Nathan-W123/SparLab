/// \file Element.hpp
/// \brief Abstract element interface.
///
/// Every element implementation supplies the four kernels the rest of the
/// library needs: a stiffness matrix, a consistent mass matrix, a
/// strain-displacement operator at a requested parametric location, and a
/// consistent nodal load vector for a traction on one boundary face. The
/// interface is dimension-generic: coordinates, constitutive matrices and the
/// strain operator are dynamic Eigen matrices sized by `dim()` and
/// `num_voigt()`, so the assembler, stress recovery and load application are
/// written once for the plane Q4 and the solid Hex8. Adding a new topology
/// (Q8, T3, Tet4, ...) means adding one subclass and one entry to
/// `make_element`; no other component changes.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <memory>
#include <vector>

namespace sparlab {

/// Element-local coordinates on the reference domain. `zeta` is ignored by
/// plane elements.
struct NaturalPoint {
  Scalar xi = 0.0;
  Scalar eta = 0.0;
  Scalar zeta = 0.0;
};

/// Quadrature orders used by an element's kernels.
struct IntegrationOptions {
  int stiffness_points = 2;  ///< points per direction for K_e
  int mass_points = 3;       ///< points per direction for M_e
  int edge_points = 2;       ///< points per direction on a boundary face for tractions
};

/// Strain-displacement data evaluated at one parametric point.
struct StrainOperator {
  Matrix b;           ///< num_voigt x num_dofs operator [1/m]
  Scalar detJ = 0.0;  ///< Jacobian determinant [m^dim]
};

/// Abstract continuum element.
class Element {
 public:
  virtual ~Element() = default;

  virtual ElementType type() const = 0;
  /// Spatial dimension (2 or 3), also the translational DOFs per node.
  virtual int dim() const = 0;
  virtual int num_nodes() const = 0;
  /// Number of boundary faces (edges in 2-D).
  virtual int num_faces() const = 0;

  /// Degrees of freedom carried by the element (num_nodes * dim).
  int num_dofs() const { return num_nodes() * dim(); }

  /// Voigt components of the element's strain operator.
  int num_voigt() const { return voigt_components(dim()); }

  /// Element stiffness matrix
  /// \f$ K_e = \int_{\Omega_e} t\, B^T D B \, d\Omega \f$ [N/m].
  /// \param coords dim x num_nodes nodal coordinates [m].
  /// \param d constitutive matrix, num_voigt x num_voigt [Pa].
  /// \param thickness out-of-plane thickness [m] (must be 1 for 3-D elements).
  virtual Matrix stiffness(const Matrix& coords, const Matrix& d, Scalar thickness,
                           const IntegrationOptions& opts) const = 0;

  /// Consistent element mass matrix
  /// \f$ M_e = \int_{\Omega_e} \rho\, t\, N^T N \, d\Omega \f$ [kg].
  virtual Matrix consistent_mass(const Matrix& coords, Scalar density, Scalar thickness,
                                 const IntegrationOptions& opts) const = 0;

  /// Strain-displacement operator and Jacobian determinant at a natural point.
  /// \throws MeshError when detJ <= 0 (inverted or degenerate element).
  virtual StrainOperator strain_operator(const Matrix& coords,
                                         const NaturalPoint& point) const = 0;

  /// Shape function values at a natural point (size num_nodes).
  virtual Vector shape_functions(const NaturalPoint& point) const = 0;

  /// Parametric locations of the stiffness quadrature points, in the order used
  /// by `stiffness`. Needed for stress recovery at Gauss points.
  virtual std::vector<NaturalPoint> stress_evaluation_points(
      const IntegrationOptions& opts) const = 0;

  /// Consistent nodal forces for a constant traction on local face
  /// `local_face` (an edge in 2-D):
  /// \f$ f_e = \int_{\Gamma_e} t\, N^T \bar{t} \, d\Gamma \f$ [N].
  /// \param traction traction vector [Pa]; the z component must be zero for a
  ///        plane element.
  /// \return vector of length num_dofs [N].
  virtual Vector boundary_traction(const Matrix& coords, int local_face,
                                   const Vector3& traction, Scalar thickness,
                                   const IntegrationOptions& opts) const = 0;

  /// Local node indices of boundary face `local_face`, from the topology's
  /// shared face table.
  /// \throws MeshError for an out-of-range face index.
  const std::vector<int>& face_nodes(int local_face) const;
};

/// Factory for the supported element topologies.
/// \throws ConfigError for an unsupported type.
std::unique_ptr<Element> make_element(ElementType type);

}  // namespace sparlab
