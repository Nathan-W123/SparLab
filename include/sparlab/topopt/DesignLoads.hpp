/// \file DesignLoads.hpp
/// \brief The load vectors of a density design: mechanical loads fixed, body
///        and thermal loads following the material the design places.
///
/// A load case's force vector splits into three parts (LoadCaseData.hpp):
/// point loads, tractions and pressures act on the domain whatever the design;
/// self-weight, body force densities and the centrifugal load act on the
/// material present; and a temperature field loads the material through its
/// thermal strain against its own stiffness. For physical densities
/// \f$\rho\f$,
/// \f[
///   f(\rho) = f_{mech} + \sum_e \gamma(\rho_e)\, A_e^T f_e^{b}
///           + \sum_e \frac{E(\rho_e)}{E_0}\, A_e^T f_e^{th},
/// \f]
/// with \f$f_e^b\f$ and \f$f_e^{th}\f$ element \f$e\f$'s body and thermal
/// loads at full density (Loads.hpp), \f$\gamma\f$ the body-load
/// interpolation and \f$E/E_0\f$ the SIMP stiffness factor
/// (SimpInterpolation.hpp). At \f$\rho = 1\f$ everywhere it is the model's own
/// load vector.
///
/// The compliance \f$c = f^T u\f$, \f$K u = f\f$, then has the gradient
/// \f[
///   \frac{\partial c}{\partial \rho_e} = 2\,u_e^T \frac{\partial f_e}{\partial \rho_e}
///     - u_e^T \frac{\partial K_e}{\partial \rho_e} u_e ,
/// \f]
/// whose first term is positive where adding material adds more load than
/// stiffness - which is why a design-dependent load needs an update that
/// accepts gradients of either sign (MMA), and why an adjoint of any other
/// response gains the term \f$\lambda_e^T \partial f_e/\partial\rho_e\f$.
///
/// The temperature field itself is the load case's and does not depend on the
/// design: a uniform or regional field is prescribed, and a conducted one -
/// whose conduction path the design would change - is refused by the
/// optimiser (TopologyOptimizer.hpp).
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/topopt/SimpInterpolation.hpp"

#include <vector>

namespace sparlab {

class DesignLoads {
 public:
  /// Split every load case of a finalised model into its fixed and
  /// design-dependent parts.
  /// \throws ModelError for a model that is not finalised.
  explicit DesignLoads(const FemModel& model);

  std::size_t num_cases() const { return fixed_.size(); }

  /// True when any load case has a body or a thermal load.
  bool design_dependent() const { return dependent_; }
  bool has_body(std::size_t l) const { return body_[l].size() > 0; }
  bool has_thermal(std::size_t l) const { return thermal_[l].size() > 0; }

  /// \f$f(\rho)\f$ of load case `l` [N] for physical densities `rho`. A case
  /// with neither part returns the model's load vector itself.
  Vector load(std::size_t l, const Vector& rho, const SimpOptions& simp) const;

  /// \f$a_e^T\,\partial f_e/\partial\rho_e\f$ for every element, for a
  /// full-length nodal vector `a` (twice the displacement for the compliance,
  /// an adjoint for any other response); zero for a case with neither part.
  Vector contract(std::size_t l, const Vector& rho, const SimpOptions& simp,
                  const Vector& a) const;

  /// \f$u_e^T f_e^{th}\f$ for every element at full density [J]; zero for a
  /// case without a temperature.
  Vector thermal_work(std::size_t l, const Vector& u) const;

  /// \f$\tfrac12\int_{\Omega_e}\varepsilon_0^T D_0\,\varepsilon_0\,dV\f$ per
  /// element [J]; empty for a case without a temperature.
  const Vector& thermal_self_energy(std::size_t l) const { return self_energy_[l]; }

  /// Nodal temperatures of case `l` [K]; empty without a temperature field.
  const Vector& temperature(std::size_t l) const;

 private:
  const FemModel& model_;
  std::vector<Vector> fixed_;        ///< mechanical part, per case
  std::vector<Matrix> body_;         ///< edofs x ne unit loads, per case (0 x 0: none)
  std::vector<Matrix> thermal_;      ///< edofs x ne full-density loads (0 x 0: none)
  std::vector<Vector> self_energy_;  ///< per element, per thermal case
  bool dependent_ = false;

  /// Scatter `factor(e) * columns.col(e)` into a zero vector.
  Vector scatter(const Matrix& columns, const Vector& factor) const;
};

}  // namespace sparlab
