/// \file Beam2.hpp
/// \brief Two-node Timoshenko beam in 3-D: three translations and three
///        rotations per node, the interdependent interpolation that makes it
///        exact for nodal loads.
///
/// **Local axes.** x' runs along the element from node 0 to node 1; y' is
/// the component of the section's orientation vector normal to x' (by
/// default the one that makes z' the projection of global Z - "up" - or of
/// global X for an element within 0.1 degree of vertical); z' = x' x y'.
/// The section's centroid lies on the axis and y', z' are its principal
/// axes (BeamSection.hpp).
///
/// **Kinematics.** A point (y', z') of the section at x' moves by
/// \f$u_x = u - y'\theta_z + z'\theta_y\f$, \f$u_y = v - z'\theta_x\f$,
/// \f$u_z = w + y'\theta_x\f$ (u, v, w the axis's displacement and
/// \f$\theta\f$ the section's rotation, in local components), so the
/// generalised strains are the stretch \f$u'\f$, the shears
/// \f$\gamma_y = v' - \theta_z\f$ and \f$\gamma_z = w' + \theta_y\f$, the
/// twist rate \f$\theta_x'\f$ and the curvatures \f$\theta_y'\f$,
/// \f$\theta_z'\f$, with the stiffnesses EA, \f$k_y GA\f$, \f$k_z GA\f$, GJ,
/// \f$EI_y\f$, \f$EI_z\f$. Torsion is Saint-Venant's (no warping restraint).
///
/// **Interpolation.** u and \f$\theta_x\f$ are linear. Bending in each plane -
/// \f$(v, \theta_z)\f$ with \f$\Phi_y = 12 EI_z / (k_y G A L^2)\f$, and
/// \f$(w, -\theta_y)\f$ with \f$\Phi_z = 12 EI_y / (k_z G A L^2)\f$ - takes
/// the interdependent interpolation (Reddy 1997): the deflection cubic, the
/// rotation quadratic and the shear strain constant, which is the solution of
/// the homogeneous Timoshenko equations. With \f$\mu = 1/(1+\Phi)\f$ and
/// \f$\xi = x'/L\f$,
/// \f[
///   v = \mu[1-3\xi^2+2\xi^3+\Phi(1-\xi)]\,v_0
///     + \mu L[\xi-2\xi^2+\xi^3+\tfrac{\Phi}{2}(\xi-\xi^2)]\,\theta_0
///     + \mu[3\xi^2-2\xi^3+\Phi\xi]\,v_1
///     + \mu L[-\xi^2+\xi^3-\tfrac{\Phi}{2}(\xi-\xi^2)]\,\theta_1,
/// \f]
/// \f[
///   \theta = \tfrac{6\mu}{L}(\xi^2-\xi)\,v_0 + \mu[1-4\xi+3\xi^2+\Phi(1-\xi)]\,\theta_0
///          + \tfrac{6\mu}{L}(\xi-\xi^2)\,v_1 + \mu[3\xi^2-2\xi+\Phi\xi]\,\theta_1 .
/// \f]
/// The nodal displacements of any assembly of these elements under nodal
/// loads are those of the Timoshenko beam itself, and there is no shear
/// locking; with \f$k = 0\f$ in a direction (a section without shear
/// deformation) \f$\Phi = 0\f$ and the element is the Euler-Bernoulli beam
/// with Hermite cubics.
///
/// **Matrices.** Stiffness \f$\int B^T D B\,dx'\f$; consistent mass
/// \f$\int H^T \mathrm{diag}(\rho A, \rho A, \rho A, \rho I_p, \rho I_y,
/// \rho I_z) H\,dx'\f$ (rotary inertia included, \f$I_p = I_y + I_z\f$);
/// geometric stiffness of the axial force N,
/// \f$N\int [u'^2+v'^2+w'^2+\tfrac{I_p}{A}\theta_x'^2+\tfrac{I_y}{A}\theta_y'^2
/// +\tfrac{I_z}{A}\theta_z'^2]\,dx'\f$ - the continuum's
/// \f$\int\sigma_{xx} u_{k,x} u_{k,x}\,dV\f$ of a uniform axial stress, the
/// bending stresses' part (lateral-torsional buckling) left out. All are
/// integrated exactly (4 Gauss points) and turned to global axes.
///
/// **Geometry input.** A 14 x 2 matrix (`geometry`): the nodal coordinates
/// in rows 0-2, the orientation vector in rows 3-5 (zero for the default),
/// the section in rows 6-11 - A, I_y, I_z, J, k_y, k_z - and the material's
/// E and G in rows 12-13, repeated in both columns: the interpolation, and
/// so the mass and the loads, depend on E / G through Phi. The constitutive
/// matrix is diag(E, G) (2 x 2), and must match.
#pragma once

#include "sparlab/elements/BeamSection.hpp"
#include "sparlab/elements/Element.hpp"

namespace sparlab {

/// Local axes of a beam element, as the rows of `rotation` (local = R
/// global).
struct BeamFrame {
  Scalar length = 0.0;
  Matrix3 rotation = Matrix3::Identity();
  Vector3 x_axis() const { return rotation.row(0).transpose(); }
  Vector3 y_axis() const { return rotation.row(1).transpose(); }
  Vector3 z_axis() const { return rotation.row(2).transpose(); }
};

/// The section resultants of a beam element at its two ends, in its local
/// axes: N, Q_y, Q_z, T, M_y, M_z [N, N m] (M_y the moment about y',
/// conjugate to the curvature theta_y'; a positive M_y stretches the fibres
/// at z' > 0, a positive M_z compresses those at y' > 0).
struct BeamEndForces {
  BeamFrame frame;
  Eigen::Matrix<Scalar, 6, 1> start = Eigen::Matrix<Scalar, 6, 1>::Zero();
  Eigen::Matrix<Scalar, 6, 1> end = Eigen::Matrix<Scalar, 6, 1>::Zero();
};

class Beam2Element final : public Element {
 public:
  ElementType type() const override { return ElementType::Beam2; }
  int dim() const override { return 3; }
  int num_nodes() const override { return 2; }
  int num_faces() const override { return 0; }
  int dofs_per_node() const override { return 6; }

  /// \param d diag(E, G) [Pa]. `thickness` is not used.
  Matrix stiffness(const Matrix& geometry, const Matrix& d, Scalar thickness,
                   const IntegrationOptions& opts) const override;

  Matrix consistent_mass(const Matrix& geometry, Scalar density, Scalar thickness,
                         const IntegrationOptions& opts) const override;

  /// Lumped mass: half the element's mass \f$\rho A L\f$ on each node's
  /// translations and half its sections' rotary inertia,
  /// \f$\tfrac12 \rho L\, R^T \mathrm{diag}(I_p, I_y, I_z) R\f$, on its
  /// rotations - a 3 x 3 block, the inertia tensor about the local axes
  /// turned to global ones. (The consistent matrix's rotational rows also
  /// carry the translational inertia the interpolation ties to the
  /// rotations, which the nodal masses already hold.)
  Matrix lumped_mass(const Matrix& geometry, Scalar density, Scalar thickness,
                     const IntegrationOptions& opts) const override;

  /// The generalised strains (u', gamma_y, gamma_z, theta_x', theta_y',
  /// theta_z') at xi in [-1, 1] as an operator on the global DOFs (6 x 12),
  /// and det J = L / 2.
  StrainOperator strain_operator(const Matrix& geometry,
                                 const NaturalPoint& point) const override;

  /// The linear shape functions (1 - xi)/2, (1 + xi)/2 of the axis.
  Vector shape_functions(const NaturalPoint& point) const override;

  /// Two Gauss points on [-1, 1].
  std::vector<NaturalPoint> stress_evaluation_points(
      const IntegrationOptions& opts) const override;
  std::vector<IntegrationPoint> integration_rule(const IntegrationOptions& opts) const override;

  /// The geometric stiffness of the axial force of `ue` (class comment).
  Matrix geometric_stiffness(const Matrix& geometry, const Matrix& d, const Vector& ue,
                             Scalar stress_scale, Scalar thickness,
                             const IntegrationOptions& opts) const override;

  Vector geometric_stiffness_derivative(const Matrix& geometry, const Matrix& d,
                                        const Vector& phi, Scalar stress_scale,
                                        Scalar thickness,
                                        const IntegrationOptions& opts) const override;

  /// A beam has no faces. \throws ConfigError (beams take line loads).
  Vector boundary_traction(const Matrix& geometry, int local_face, const Vector3& traction,
                           Scalar thickness, const IntegrationOptions& opts) const override;

  /// Consistent nodal forces of a uniform force per unit length `q` [N/m]
  /// in global components: \f$\int H^T q\,dx'\f$.
  Vector line_load(const Matrix& geometry, const Vector3& q) const;

  /// The end resultants of the element displacement `ue` (global) under the
  /// uniform line load `q` [N/m] it carries: the local end forces
  /// \f$K u - f_q\f$, the start's negated.
  BeamEndForces end_forces(const Matrix& geometry, const Matrix& d, const Vector& ue,
                           const Vector3& q) const;

  /// The 14 x 2 geometry matrix of a beam between `x0` and `x1` with the
  /// resolved section `section` (its orientation included) and the moduli
  /// E and G [Pa].
  static Matrix geometry(const Vector3& x0, const Vector3& x1, const BeamSection& section,
                         Scalar youngs_modulus, Scalar shear_modulus);

  /// The local axes and length of the element.
  /// \throws MeshError for a zero length or an orientation vector along the
  ///         axis.
  static BeamFrame frame(const Matrix& geometry);
};

}  // namespace sparlab
