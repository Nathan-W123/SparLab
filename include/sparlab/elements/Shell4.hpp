/// \file Shell4.hpp
/// \brief Four-node MITC4 shell (Dvorkin and Bathe): a degenerated
///        continuum with nodal directors, global rotations and assumed
///        transverse shear strains.
///
/// **Kinematics.** A point of the shell at the natural coordinates
/// \f$(r, s)\in[-1,1]^2\f$ of its mid-surface and \f$\zeta\in[-1,1]\f$ through
/// its thickness t lies at, and moves by,
/// \f[
///   X = \sum_k N_k \left(x_k + \zeta\tfrac{t}{2} V_k\right), \qquad
///   u = \sum_k N_k \left(u_k + \zeta\tfrac{t}{2}\, \theta_k \times V_k\right),
/// \f]
/// with the bilinear \f$N_k\f$ of the Q4, the nodal directors \f$V_k\f$ (unit
/// vectors normal to the surface at the nodes) and, per node, the three
/// translations \f$u_k\f$ and the three rotations \f$\theta_k\f$ about the
/// global axes - six DOFs, in that order. Small rotations: the analysis is
/// linear.
///
/// **Strains.** The covariant strain components
/// \f$\tilde\varepsilon_{ij} = \tfrac12 (g_i\cdot u_{,j} + g_j\cdot u_{,i})\f$
/// on the base vectors \f$g_r, g_s, g_\zeta\f$ of the point; the in-plane ones
/// (rr, ss, rs) at the point, the transverse shears (MITC) interpolated from
/// the mid-edges - \f$\tilde\varepsilon_{r\zeta}\f$ from A = (0, 1) and
/// C = (0, -1), \f$\tilde\varepsilon_{s\zeta}\f$ from B = (-1, 0) and
/// D = (1, 0) - which keeps a thin shell from locking in shear. The normal
/// strain through the thickness is not used (plane stress). They are turned
/// into the local Cartesian components by
/// \f$\varepsilon_{ab} = T_{ai} T_{bj} \tilde\varepsilon_{ij}\f$,
/// \f$T = E^T G^{-T}\f$, E the local frame and G the base vectors.
///
/// **Local frame.** \f$e_3\f$ along the interpolated director; \f$e_1\f$ the
/// projection of the global x axis onto the plane normal to it (the global z
/// axis when x lies within 0.1 degree of \f$e_3\f$), \f$e_2 = e_3\times e_1\f$
/// - the frame the resultants are reported in.
///
/// **Material.** Plane stress in \f$(e_1, e_2)\f$ from the 3 x 3 plane-stress
/// matrix `d`, and transverse shear \f$k G\f$ with \f$k = 5/6\f$ and
/// \f$G = d_{33}\f$.
///
/// **Integration.** 2 x 2 points in the plane by 2 through the thickness for
/// the stiffness and the geometric stiffness; `mass_points` in the plane by
/// 3 through the thickness for the mass (the rotary inertia
/// \f$\rho t^3/12\f$ included).
///
/// **Drilling.** A rotation about the normal moves no point of the shell,
/// so it is tied to the in-plane rotation of the mid-surface by a penalty -
/// Hughes and Brezzi's drilling constraint - at the 2 x 2 points:
/// \f$\tfrac12 k_d \int (n\cdot\theta - \omega)^2 dA\f$ with
/// \f$\omega = \tfrac12\, n\cdot(a^\alpha\times u_{,\alpha})\f$ the rotation of
/// the mid-surface about its normal n (\f$a^\alpha\f$ the in-plane dual base
/// vectors) and \f$k_d = \alpha\, G t\f$ (\f$\alpha\f$ = `drilling_factor`).
/// Under a rigid rotation both terms equal \f$n\cdot\omega\f$ on any
/// geometry, so the six rigid-body motions stay free of stiffness; at a fold,
/// where a rotation about one element's normal bends its neighbour, the
/// penalty is what makes the two compatible.
///
/// **Geometry input.** The element geometry is a 6 x 4 matrix: the nodal
/// coordinates in rows 0-2 and the directors in rows 3-5 (FemModel builds
/// them, averaged over the elements at a node). A 3 x 4 matrix of coordinates
/// alone takes each node's director as the element's own normal there.
#pragma once

#include "sparlab/elements/Element.hpp"

namespace sparlab {

/// The stress resultants of a shell element at one point of its mid-surface,
/// in the local frame (e1, e2, e3) at that point.
struct ShellResultants {
  Vector3 e1 = Vector3::Zero();
  Vector3 e2 = Vector3::Zero();
  Vector3 e3 = Vector3::Zero();
  /// Membrane forces N11, N22, N12 [N/m].
  Vector3 membrane = Vector3::Zero();
  /// Bending moments M11, M22, M12 [N m / m]; M11 > 0 stretches the side the
  /// director points to (zeta > 0) along e1.
  Vector3 moment = Vector3::Zero();
  /// Transverse shear forces Q13, Q23 [N/m].
  Eigen::Vector2d shear = Eigen::Vector2d::Zero();
  /// In-plane stresses sigma11, sigma22, sigma12 on the side the director
  /// points to (zeta = 1, "top") and on the other one (zeta = -1) [Pa].
  Vector3 stress_top = Vector3::Zero();
  Vector3 stress_bottom = Vector3::Zero();
  /// von Mises stress on the two faces (plane stress, no transverse shear
  /// there) and on the mid-surface, with the membrane stresses and the
  /// transverse shear stresses at their parabolic peak 3 Q / (2 t) [Pa].
  Scalar von_mises_top = 0.0;
  Scalar von_mises_bottom = 0.0;
  Scalar von_mises_mid = 0.0;
};

class Shell4Element final : public Element {
 public:
  /// The default drilling stiffness factor (the class comment).
  static constexpr Scalar kDefaultDrillingFactor = 1.0e-3;
  /// Transverse shear correction factor.
  static constexpr Scalar kShearFactor = 5.0 / 6.0;

  explicit Shell4Element(Scalar drilling_factor = kDefaultDrillingFactor);

  ElementType type() const override { return ElementType::Shell4; }
  int dim() const override { return 3; }
  int num_nodes() const override { return 4; }
  int num_faces() const override { return 4; }
  int dofs_per_node() const override { return 6; }
  Scalar drilling_factor() const { return drilling_factor_; }

  /// \param d the 3 x 3 plane-stress constitutive matrix [Pa].
  /// \param thickness the shell thickness [m].
  Matrix stiffness(const Matrix& geometry, const Matrix& d, Scalar thickness,
                   const IntegrationOptions& opts) const override;

  Matrix consistent_mass(const Matrix& geometry, Scalar density, Scalar thickness,
                         const IntegrationOptions& opts) const override;

  /// Lumped mass: the translations' diagonal of the consistent mass scaled
  /// to the element's mass (Hinton, Rock and Zienkiewicz), and on each
  /// node's rotations the same share of the element's rotary-inertia tensor
  /// - the consistent rotational blocks summed over all node pairs, the
  /// inertia \f$\int \rho\, z^2 (|d|^2 I - d d^T)\,dV\f$ of a uniform
  /// rotation of the fibres (d the interpolated director, z the distance from
  /// the mid-surface). On a flat element that is \f$\rho t^3/12\f$ per unit
  /// area about the in-plane axes and nothing about the normal, in any frame.
  Matrix lumped_mass(const Matrix& geometry, Scalar density, Scalar thickness,
                     const IntegrationOptions& opts) const override;

  /// The local strain operator at (r, s, zeta) - six rows in the Voigt order
  /// of a solid, (11, 22, 33, 12, 23, 13) in the local frame, the 33 row
  /// zero - with the transverse shears of the MITC interpolation, and the
  /// volume Jacobian. The geometry must carry its directors (6 x 4), since
  /// the thickness enters the operator: use `strain_operator_at`.
  StrainOperator strain_operator(const Matrix& geometry,
                                 const NaturalPoint& point) const override;

  /// The same with the thickness given.
  StrainOperator strain_operator_at(const Matrix& geometry, Scalar thickness,
                                    const NaturalPoint& point) const;

  Vector shape_functions(const NaturalPoint& point) const override;

  /// 2 x 2 in the plane by 2 through the thickness.
  std::vector<NaturalPoint> stress_evaluation_points(
      const IntegrationOptions& opts) const override;
  std::vector<IntegrationPoint> integration_rule(const IntegrationOptions& opts) const override;

  /// The geometric stiffness of the in-plane stresses of `ue`,
  /// \f$K_G = \int \sigma_{ab}\, D_a^T D_b\, dV\f$ over the local in-plane
  /// directions a, b, with \f$D_a q = \partial u/\partial x_a\f$ the gradient
  /// of the shell's displacement field (rotation terms included).
  Matrix geometric_stiffness(const Matrix& geometry, const Matrix& d, const Vector& ue,
                             Scalar stress_scale, Scalar thickness,
                             const IntegrationOptions& opts) const override;

  /// \f$g\f$ with \f$\phi^T K_G(u)\phi = g^T u\f$ (the buckling-load
  /// sensitivity's adjoint load).
  Vector geometric_stiffness_derivative(const Matrix& geometry, const Matrix& d,
                                        const Vector& phi, Scalar stress_scale,
                                        Scalar thickness,
                                        const IntegrationOptions& opts) const override;

  /// Consistent nodal forces of a constant traction [Pa] on edge
  /// `local_face`, over the edge's area (its length times the thickness).
  Vector boundary_traction(const Matrix& geometry, int local_face, const Vector3& traction,
                           Scalar thickness, const IntegrationOptions& opts) const override;

  /// Consistent nodal forces of a uniform pressure on the mid-surface,
  /// acting against its normal \f$g_r\times g_s\f$ (right-handed about the
  /// node order): \f$f_k = -p\int N_k\, (g_r\times g_s)\, dr\, ds\f$.
  Vector pressure_load(const Matrix& geometry, Scalar pressure,
                       const IntegrationOptions& opts) const;

  /// The resultants of the element displacement `ue` at (r, s).
  ShellResultants resultants(const Matrix& geometry, const Matrix& d, const Vector& ue,
                             Scalar thickness, Scalar r, Scalar s) const;

  /// Mid-surface area [m^2].
  static Scalar area(const Matrix& geometry);

  /// Unit normal \f$g_r\times g_s\f$ of the mid-surface at (r, s).
  static Vector3 normal(const Matrix& geometry, Scalar r, Scalar s);

 private:
  Scalar drilling_factor_;
};

}  // namespace sparlab
