/// \file SimpInterpolation.hpp
/// \brief SIMP material interpolation and its derivative.
///
/// SparLab uses the *modified* SIMP law (Sigmund 2007), which keeps a stiffness
/// floor without making the derivative vanish at \f$\rho = 0\f$:
/// \f[
///   \frac{E(\rho)}{E_0} = \epsilon + (1-\epsilon)\, \rho^{p},
///   \qquad
///   \frac{1}{E_0}\frac{\mathrm{d}E}{\mathrm{d}\rho} =
///   p\,(1-\epsilon)\,\rho^{p-1},
/// \f]
/// with \f$\epsilon = E_{min}/E_0\f$ (default \f$10^{-9}\f$) and penalty
/// \f$p \ge 1\f$. Because the floor is additive rather than a lower bound on
/// \f$\rho\f$, design variables may reach exactly 0 and 1 while
/// \f$\mathbf{K}\f$ stays positive definite.
///
/// The mass of an element is interpolated separately, since penalising mass the
/// same way as stiffness is not physical and using a linear mass law with a
/// penalised stiffness law creates spurious low-frequency modes localised in
/// void regions (the classical SIMP eigenvalue artefact). The available laws
/// are documented in `MassInterpolation`.
#pragma once

#include "sparlab/core/Types.hpp"

#include <string>

namespace sparlab {

/// Mass interpolation law \f$ m(\rho) \f$, used for modal analysis of a design.
enum class MassInterpolation {
  /// \f$ m = \rho \f$. Physically correct for the mass of a graded material but
  /// pairs badly with a penalised stiffness at low density: the stiffness/mass
  /// ratio collapses and void regions produce spurious low-frequency modes.
  Linear,
  /// \f$ m = \epsilon_m + (1-\epsilon_m)\rho^{p} \f$ with the *same* penalty as
  /// the stiffness law, so \f$\omega^2 \sim E/m\f$ stays bounded in void
  /// regions and no spurious mode appears. Default for SIMP modal analysis.
  PenaltyMatched
};

std::string to_string(MassInterpolation law);
MassInterpolation parse_mass_interpolation(const std::string& text);

struct SimpOptions {
  Scalar penalty = 3.0;        ///< p [-], >= 1
  Scalar emin_ratio = 1.0e-9;  ///< epsilon = E_min / E_0 [-], in (0, 1)
  Scalar mass_floor = 1.0e-9;  ///< epsilon_m for PenaltyMatched [-]
  MassInterpolation mass_law = MassInterpolation::PenaltyMatched;

  void validate() const;
};

/// \f$E(\rho)/E_0\f$ for one density.
Scalar simp_stiffness_factor(Scalar rho, const SimpOptions& options);

/// \f$\mathrm{d}(E/E_0)/\mathrm{d}\rho\f$ for one density.
Scalar simp_stiffness_derivative(Scalar rho, const SimpOptions& options);

/// \f$m(\rho)\f$ for one density.
Scalar simp_mass_factor(Scalar rho, const SimpOptions& options);

/// Element-wise \f$E(\rho)/E_0\f$ for a density vector.
Vector simp_stiffness_factors(const Vector& rho, const SimpOptions& options);

/// Element-wise \f$\mathrm{d}(E/E_0)/\mathrm{d}\rho\f$ for a density vector.
Vector simp_stiffness_derivatives(const Vector& rho, const SimpOptions& options);

/// Element-wise \f$m(\rho)\f$ for a density vector.
Vector simp_mass_factors(const Vector& rho, const SimpOptions& options);

}  // namespace sparlab
