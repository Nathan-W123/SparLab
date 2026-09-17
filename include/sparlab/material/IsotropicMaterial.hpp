/// \file IsotropicMaterial.hpp
/// \brief Linear-elastic isotropic material and its 2-D constitutive matrices.
///
/// Voigt ordering used throughout SparLab is
/// \f$ \{\sigma_{xx}, \sigma_{yy}, \sigma_{xy}\} \f$ paired with the
/// *engineering* strain vector \f$ \{\varepsilon_{xx}, \varepsilon_{yy},
/// \gamma_{xy}\} \f$ where \f$\gamma_{xy} = 2\varepsilon_{xy}\f$.
///
/// Plane stress (\f$\sigma_{zz}=0\f$):
/// \f[
///   \mathbf{D} = \frac{E}{1-\nu^2}
///   \begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \tfrac{1-\nu}{2}
///   \end{bmatrix}
/// \f]
/// Plane strain (\f$\varepsilon_{zz}=0\f$):
/// \f[
///   \mathbf{D} = \frac{E}{(1+\nu)(1-2\nu)}
///   \begin{bmatrix} 1-\nu & \nu & 0 \\ \nu & 1-\nu & 0 \\ 0 & 0 &
///   \tfrac{1-2\nu}{2}\end{bmatrix}
/// \f]
/// In plane strain the out-of-plane stress is
/// \f$\sigma_{zz} = \nu(\sigma_{xx}+\sigma_{yy})\f$, which StressRecovery
/// accounts for when forming the von Mises stress.
#pragma once

#include "sparlab/core/Types.hpp"

#include <string>

namespace sparlab {

/// Isotropic linear-elastic material described by \f$(E, \nu, \rho)\f$.
class IsotropicMaterial {
 public:
  /// \param youngs_modulus E [Pa], must be > 0.
  /// \param poisson_ratio nu [-], must lie in (-1, 0.5) for plane strain and
  ///        (-1, 1) excluding +-1 for plane stress. SparLab restricts the input
  ///        to (-1, 0.5) so both idealisations remain well posed.
  /// \param density rho [kg/m^3], must be >= 0 (0 disables modal analysis).
  /// \param name optional label carried into result files.
  IsotropicMaterial(Scalar youngs_modulus, Scalar poisson_ratio, Scalar density,
                    std::string name = "material");

  Scalar youngs_modulus() const { return e_; }
  Scalar poisson_ratio() const { return nu_; }
  Scalar density() const { return rho_; }
  const std::string& name() const { return name_; }

  /// Shear modulus \f$G = E / (2(1+\nu))\f$ [Pa].
  Scalar shear_modulus() const { return e_ / (2.0 * (1.0 + nu_)); }

  /// Constitutive matrix for the requested 2-D idealisation [Pa].
  Matrix3 constitutive(StressState state) const;

  /// Plane-stress constitutive matrix [Pa].
  Matrix3 plane_stress_matrix() const;

  /// Plane-strain constitutive matrix [Pa].
  Matrix3 plane_strain_matrix() const;

  /// Return a copy of this material with a scaled Young's modulus. Used by the
  /// aerospace material-stiffness sweep.
  IsotropicMaterial with_youngs_modulus(Scalar e) const {
    return IsotropicMaterial(e, nu_, rho_, name_);
  }

 private:
  Scalar e_;
  Scalar nu_;
  Scalar rho_;
  std::string name_;
};

}  // namespace sparlab
