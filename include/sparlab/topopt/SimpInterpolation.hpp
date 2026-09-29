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
///
/// **Body loads.** Self-weight, a body force density and the centrifugal load
/// act on the material present, so an element's body load is its full-density
/// load times \f$\gamma(\rho)\f$, and physically \f$\gamma = \rho\f$: a
/// graded element weighs its volume fraction. Against the penalised stiffness
/// that ratio is unbounded as \f$\rho \to 0\f$ - load \f$\rho\f$ on stiffness
/// \f$\rho^p\f$ - and a near-void element sags under its own weight without
/// limit (the parasitic effect of Bruyneel and Duysinx, 2005). Below a
/// threshold \f$\rho_t\f$ the load therefore follows
/// \f[
///   \gamma(\rho) = \rho_t \left[ p\,x^{p} - (p-1)\,x^{p+1} \right], \qquad
///   x = \rho/\rho_t < 1 ,
/// \f]
/// and \f$\gamma = \rho\f$ above it: continuous with a continuous slope at
/// \f$\rho_t\f$ (value \f$\rho_t\f$, slope 1), increasing, and vanishing like
/// \f$\rho^p\f$, so the load of an element never exceeds
/// \f$p\,\rho_t^{1-p}\f$ times its stiffness's share of the solid's. At
/// \f$p = 1\f$ it is \f$\rho\f$ everywhere. The load of every element at or
/// above the threshold - and of a 0/1 design - is the physical one.
///
/// **Thermal loads.** The thermal load of a SIMP element is its free strain
/// \f$\alpha\,\Delta T\f$ against its own stiffness,
/// \f$\int B^T D(\rho)\,\varepsilon_0\,dV\f$, i.e. the full-density thermal
/// load times the stiffness factor \f$E(\rho)/E_0\f$: the thermal stress
/// coefficient \f$E\alpha\f$ interpolated with the stiffness, so a graded
/// element expands freely without stress, as the material it models would.
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
  /// \f$\rho_t\f$ of the body-load interpolation [-], in [0, 1); 0 keeps
  /// \f$\gamma = \rho\f$ at every density.
  Scalar body_load_threshold = 0.1;

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

/// \f$\gamma(\rho)\f$, an element's share of its full-density body load.
Scalar body_load_factor(Scalar rho, const SimpOptions& options);

/// \f$\mathrm{d}\gamma/\mathrm{d}\rho\f$.
Scalar body_load_derivative(Scalar rho, const SimpOptions& options);

}  // namespace sparlab
