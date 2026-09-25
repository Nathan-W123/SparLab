/// \file Projection.hpp
/// \brief Smoothed Heaviside projection of the filtered density, with
///        continuation of its sharpness.
///
/// A density filter makes the optimisation well posed but leaves a band of
/// intermediate density one filter radius wide around every member. SIMP
/// credits that grey material with \f$\rho^p\f$ of its stiffness, so the
/// optimiser's compliance describes a structure that does not exist; the
/// thresholded solid the run exports can be markedly stiffer or softer.
///
/// The projection of Wang, Lazarov and Sigmund (2011) maps the filtered
/// density \f$\tilde\rho\f$ to the physical density
/// \f[
///   \bar\rho = \frac{\tanh(\beta\eta) + \tanh\!\big(\beta(\tilde\rho - \eta)\big)}
///                   {\tanh(\beta\eta) + \tanh\!\big(\beta(1 - \eta)\big)},
///   \qquad
///   \frac{d\bar\rho}{d\tilde\rho} =
///   \frac{\beta\,\big(1 - \tanh^2(\beta(\tilde\rho - \eta))\big)}
///        {\tanh(\beta\eta) + \tanh\!\big(\beta(1 - \eta)\big)} ,
/// \f]
/// a smooth step at the threshold \f$\eta\f$ whose sharpness \f$\beta\f$ is
/// raised during the run. \f$\bar\rho(0) = 0\f$ and \f$\bar\rho(1) = 1\f$ for
/// every \f$\beta\f$, the map is strictly increasing, and it tends to the
/// identity as \f$\beta \to 0\f$ and to a step as \f$\beta \to \infty\f$.
///
/// The chain rule stays exact: with the density filter
/// \f$\tilde\rho = \hat H x\f$,
/// \f[
///   \frac{\partial f}{\partial x} =
///   \hat H^T\!\left(\frac{d\bar\rho}{d\tilde\rho}\odot
///   \frac{\partial f}{\partial \bar\rho}\right)
/// \f]
/// for the compliance, the volume and the stress constraint alike, and the
/// tests compare it with central differences. The volume constraint acts on
/// \f$\bar\rho\f$, the material that is actually there.
///
/// Continuation starts at a small \f$\beta\f$ (nearly the plain filter, where
/// the problem is smooth) and multiplies it by `beta_factor` every
/// `beta_interval` iterations - or earlier, once the design has stopped
/// changing at the current \f$\beta\f$ - up to `beta_max`. Convergence is only
/// accepted at `beta_max`. The projection needs the density filter (or no
/// filter): the sensitivity filter has no chain rule to extend.
///
/// **Robust formulation** (Wang, Lazarov and Sigmund 2011). One filtered
/// field projected at three thresholds gives three designs: the *eroded*
/// one (\f$\eta + \Delta\eta\f$), the *intermediate* blueprint
/// (\f$\eta\f$) and the *dilated* one (\f$\eta - \Delta\eta\f$), as if the
/// part came out uniformly thinner or thicker than drawn. The objective and
/// every constraint but the volume act on the eroded design - for compliance
/// and buckling the worst of the three - and the volume constraint on the
/// dilated one, with its target rescaled every `robust_volume_interval`
/// iterations to \f$V^*\,V_d/V_i\f$ so that the blueprint meets the volume
/// fraction. With OC, rescaling at every iteration locked the MBB beam into a
/// period-2 cycle at beta = 32 that no stopping rule accepts; rescaling every
/// 20 iterations left a cycle small enough for the objective-stall test
/// (docs/benchmarks.md, section 13). A member thinner than the erosion
/// vanishes from the eroded
/// design, so the optimiser gains nothing from it: the blueprint's members
/// and gaps keep a minimum size set by the filter radius and
/// \f$\Delta\eta\f$, which LengthScale.hpp measures after the run. The
/// blueprint is the design that is reported and exported.
#pragma once

#include "sparlab/core/Types.hpp"

namespace sparlab {

struct ProjectionOptions {
  bool enabled = false;
  Scalar eta = 0.5;          ///< threshold of the step, in (0, 1)
  Scalar beta_start = 1.0;   ///< initial sharpness (> 0)
  Scalar beta_max = 32.0;    ///< final sharpness (>= beta_start)
  Scalar beta_factor = 2.0;  ///< multiplier per continuation stage (> 1)
  int beta_interval = 50;    ///< most iterations spent at one beta (>= 1)
  /// Advance to the next beta as soon as the design change at the current
  /// one falls below the optimiser's change tolerance.
  bool advance_on_convergence = true;
  /// Robust formulation over eroded / intermediate / dilated designs.
  bool robust = false;
  Scalar robust_delta = 0.1;       ///< Delta eta: thresholds eta +- robust_delta
  int robust_volume_interval = 1;  ///< iterations between dilated-target updates
  /// Without the robust formulation: evaluate the final design's eroded and
  /// dilated variants (thresholds eta +- robust_delta) once, to show how the
  /// part would perform if it came out uniformly thinner or thicker.
  bool erosion_check = false;

  Scalar eroded_eta() const { return eta + robust_delta; }
  Scalar dilated_eta() const { return eta - robust_delta; }

  /// \throws ConfigError with the offending key.
  void validate() const;

  /// Number of beta levels from beta_start to beta_max.
  int num_stages() const;
  /// Beta of stage `k` (0-based), capped at beta_max.
  Scalar beta_of_stage(int k) const;
};

/// \f$\bar\rho\f$ of every entry of `filtered`.
Vector heaviside_project(const Vector& filtered, Scalar beta, Scalar eta);

/// \f$d\bar\rho / d\tilde\rho\f$ of every entry of `filtered`.
Vector heaviside_derivative(const Vector& filtered, Scalar beta, Scalar eta);

/// Scalar forms.
/// \{
Scalar heaviside_project(Scalar filtered, Scalar beta, Scalar eta);
Scalar heaviside_derivative(Scalar filtered, Scalar beta, Scalar eta);
/// \}

}  // namespace sparlab
