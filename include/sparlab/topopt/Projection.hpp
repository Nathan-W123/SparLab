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
