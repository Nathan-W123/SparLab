/// \file StressRecovery.hpp
/// \brief Element strain / stress recovery and derived invariants.
///
/// Strain is obtained from the element's strain-displacement operator,
/// \f$ \boldsymbol{\varepsilon} = \mathbf{B}(\xi,\eta)\, \mathbf{u}_e \f$, and
/// stress from the constitutive law, \f$ \boldsymbol{\sigma} = s_e \mathbf{D}
/// \boldsymbol{\varepsilon} \f$, where \f$ s_e \f$ is the element stiffness
/// scale factor (1 for a plain analysis, \f$E(\rho_e)/E_0\f$ for a SIMP
/// design). With \f$ s_e \f$ applied the result is the *macroscopic* stress
/// carried by the element; the unscaled value is the stress in the underlying
/// solid material and is reported alongside it.
///
/// Reported element values are the average over the stiffness quadrature
/// points. Nodal values are area-weighted averages of the adjacent element
/// values (simple averaging, not superconvergent patch recovery), and are used
/// for smooth contour plots only.
///
/// von Mises stress uses the correct out-of-plane stress for the active
/// idealisation: \f$\sigma_{zz} = 0\f$ for plane stress and
/// \f$\sigma_{zz} = \nu(\sigma_{xx}+\sigma_{yy})\f$ for plane strain, so
/// \f[
///   \sigma_{vm} = \sqrt{\tfrac{1}{2}\left[(\sigma_{xx}-\sigma_{yy})^2 +
///   (\sigma_{yy}-\sigma_{zz})^2 + (\sigma_{zz}-\sigma_{xx})^2\right] +
///   3\sigma_{xy}^2}.
/// \f]
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"

namespace sparlab {

/// Voigt stress/strain storage: 3 x N columns, one per entity.
using VoigtField = Eigen::Matrix<Scalar, kVoigt, Eigen::Dynamic>;

struct StressField {
  VoigtField element_strain;            ///< 3 x num_elements [-]
  VoigtField element_stress;            ///< 3 x num_elements, macroscopic [Pa]
  VoigtField element_solid_stress;      ///< 3 x num_elements, solid material [Pa]
  Vector element_von_mises;             ///< num_elements, from element_stress [Pa]
  Vector element_solid_von_mises;       ///< num_elements, from solid stress [Pa]
  Vector element_principal_max;         ///< num_elements, sigma_1 [Pa]
  Vector element_principal_min;         ///< num_elements, sigma_2 [Pa]
  Vector element_strain_energy;         ///< num_elements [J]

  VoigtField nodal_stress;              ///< 3 x num_nodes, averaged [Pa]
  Vector nodal_von_mises;               ///< num_nodes [Pa]
};

/// von Mises stress from a 2-D Voigt stress vector.
Scalar von_mises(const Vector3& voigt_stress, StressState state, Scalar poisson);

/// In-plane principal stresses (sigma_1 >= sigma_2).
void principal_stresses(const Vector3& voigt_stress, Scalar& s1, Scalar& s2);

/// Strain at one parametric point of one element.
Vector3 element_strain_at(const FemModel& model, Index element, const NaturalPoint& point,
                          const Vector& displacement);

/// Stress at one parametric point of one element (macroscopic, scaled).
Vector3 element_stress_at(const FemModel& model, Index element, const NaturalPoint& point,
                          const Vector& displacement, Scalar stiffness_scale = 1.0);

/// Recover all strain / stress data for a displacement field.
/// \param stiffness_scale optional per-element SIMP factors \f$E(\rho)/E_0\f$.
StressField recover_stresses(const FemModel& model, const Assembler& assembler,
                             const Vector& displacement,
                             const Vector* stiffness_scale = nullptr);

}  // namespace sparlab
