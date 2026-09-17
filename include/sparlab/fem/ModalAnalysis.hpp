/// \file ModalAnalysis.hpp
/// \brief Undamped free-vibration analysis of the constrained model.
///
/// The generalised eigenproblem on the *free* DOFs is
/// \f[
///   K_{ff}\, \phi = \lambda\, M_{ff}\, \phi, \qquad
///   \lambda = \omega^2, \qquad f = \frac{\omega}{2\pi}.
/// \f]
/// Prescribed DOFs are removed by the same partitioning used for the static
/// solve, so no constrained DOF can contribute a spurious mode and no penalty
/// stiffness distorts the spectrum.
///
/// **Algorithm.** Bathe's subspace iteration with shift-invert acceleration:
/// starting from a subspace \f$X_0 \in \mathbb{R}^{n\times q}\f$ with
/// \f$q = \min(n, \max(2m, m+8))\f$,
/// \f[
///   \bar{X}_{k+1} = (K_{ff} - \sigma M_{ff})^{-1} M_{ff} X_k,
/// \f]
/// followed by a Rayleigh-Ritz projection onto the \f$q\f$-dimensional
/// subspace, a dense symmetric generalised eigensolve of the projected pair,
/// and rotation of the basis. The single sparse Cholesky factorisation is
/// reused for every iteration and every vector. Convergence is measured on the
/// relative change of the \f$m\f$ requested eigenvalues.
///
/// **Validity checks.** Each converged pair is screened for
///   * negative \f$\lambda\f$ (non-physical: the stiffness matrix is not
///     positive definite on the free DOFs),
///   * \f$\lambda\f$ indistinguishable from zero relative to
///     \f$\|K_{ff}\|/\|M_{ff}\|\f$ (an unsuppressed rigid-body or mechanism
///     mode),
///   * a large residual \f$\|K\phi - \lambda M\phi\| / \|\lambda M \phi\|\f$.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"

#include <string>
#include <vector>

namespace sparlab {

struct ModalAnalysisOptions {
  int num_modes = 6;               ///< number of lowest modes requested
  MassType mass_type = MassType::Consistent;
  int max_iterations = 200;        ///< subspace iteration cap
  Scalar tolerance = 1.0e-10;      ///< relative eigenvalue change for convergence
  /// Shift \f$\sigma\f$ applied as a multiple of the smallest diagonal ratio
  /// \f$\min_i K_{ii}/M_{ii}\f$. A small negative shift keeps
  /// \f$K - \sigma M\f$ positive definite while accelerating convergence to the
  /// lowest modes; 0 disables shifting.
  Scalar shift_factor = 0.0;
  /// Modes whose eigenvalue is below `rigid_body_ratio` times the mean
  /// stiffness/mass scale are flagged as rigid-body / mechanism modes.
  Scalar rigid_body_ratio = 1.0e-10;
  /// Maximum accepted eigenpair residual.
  Scalar residual_tolerance = 1.0e-6;
  /// Fixed seed for the random part of the starting subspace, so runs are
  /// bit-for-bit reproducible.
  unsigned int seed = 20240917u;
};

struct ModalResult {
  Vector eigenvalues;              ///< lambda = omega^2 [1/s^2], ascending
  Vector angular_frequencies;      ///< omega [rad/s]
  Vector frequencies_hz;           ///< f [Hz]
  Matrix mode_shapes;              ///< num_dofs x num_modes, M-orthonormal [m]
  Vector modal_residuals;          ///< per mode, ||K phi - lambda M phi|| / ||lambda M phi||
  Vector modal_mass_fraction;      ///< fraction of the mode's kinetic energy
                                   ///< carried by low-density elements (topology runs)
  Scalar total_mass = 0.0;         ///< sum of the assembled mass matrix / 2 per direction [kg]
  int iterations = 0;
  int subspace_size = 0;
  bool converged = false;
  Scalar final_change = 0.0;
  std::vector<std::string> warnings;
};

/// Solve the generalised eigenproblem for the lowest `num_modes` modes.
/// \param stiffness_scale optional per-element stiffness factors.
/// \param mass_scale optional per-element mass factors.
/// \throws ModelError for an ill-posed model, SolverError for a non-positive
///         definite pair, ConvergenceError when the iteration budget is spent.
ModalResult solve_modal(const FemModel& model, const Assembler& assembler,
                        const ModalAnalysisOptions& options,
                        const Vector* stiffness_scale = nullptr,
                        const Vector* mass_scale = nullptr);

/// Reference first bending frequency of a uniform cantilever from
/// Euler-Bernoulli beam theory [Hz]:
/// \f$ f_1 = \frac{\beta_1^2}{2\pi} \sqrt{\dfrac{EI}{\rho A L^4}} \f$ with
/// \f$\beta_1 L = 1.87510407\f$.
/// \param e Young's modulus [Pa], \param rho density [kg/m^3],
/// \param length beam length [m], \param height section height [m],
/// \param thickness section width [m].
Scalar cantilever_bending_frequency(Scalar e, Scalar rho, Scalar length, Scalar height,
                                    Scalar thickness, int mode_number = 1);

/// Longitudinal (axial) natural frequency of a fixed-free uniform rod [Hz]:
/// \f$ f_n = \dfrac{2n-1}{4L}\sqrt{E/\rho} \f$.
///
/// A 2-D plane-stress cantilever carries these axial modes interleaved with its
/// bending modes, and they are a *different* analytical reference from the
/// bending series: the axial wave speed \f$\sqrt{E/\rho}\f$ is exact for a bar
/// in uniaxial stress, which the plane-stress model reproduces exactly when
/// \f$\nu = 0\f$. At non-zero Poisson ratio the lateral constraint of a finite
/// section raises the computed frequency slightly above this value.
/// \param mode_number n >= 1.
Scalar rod_axial_frequency(Scalar e, Scalar rho, Scalar length, int mode_number = 1);

}  // namespace sparlab
