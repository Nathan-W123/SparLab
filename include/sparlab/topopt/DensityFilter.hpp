/// \file DensityFilter.hpp
/// \brief Linear-hat density and sensitivity filters.
///
/// Both filters use the same weight kernel over element centroids,
/// \f[
///   H_{ei} = \max\!\left(0,\; r_{min} - \lVert x_e - x_i \rVert\right),
/// \f]
/// assembled into a sparse matrix whose rows are normalised so that
/// \f$\sum_i \hat{H}_{ei} = 1\f$.
///
/// **Density filter** (Bruns & Tortorelli 2001, Bourdin 2001) - the default.
/// The physical density is a linear map of the design variables,
/// \f$ \tilde{\rho} = \hat{H} x \f$, so the chain rule is exact:
/// \f[
///   \frac{\partial f}{\partial x_i} =
///   \sum_e \hat{H}_{ei}\, \frac{\partial f}{\partial \tilde{\rho}_e}
///   \;=\; \left(\hat{H}^T \nabla_{\tilde\rho} f\right)_i .
/// \f]
/// Because the map is exact, analytical sensitivities can be - and are -
/// verified against central finite differences.
///
/// **Sensitivity filter** (Sigmund 1997). The design variable *is* the density
/// and only the gradient is smoothed,
/// \f[
///   \widehat{\frac{\partial c}{\partial x_e}} =
///   \frac{\sum_i H_{ei}\, x_i\, \partial c/\partial x_i}
///        {\max(\gamma, x_e)\sum_i H_{ei}},
/// \f]
/// which suppresses checkerboarding just as effectively but is a heuristic: the
/// filtered gradient is *not* the gradient of any objective, so a
/// finite-difference check of the filtered value is meaningless by
/// construction. This is why SparLab defaults to the density filter and runs
/// its sensitivity verification on that path.
///
/// A radius below one element size leaves the filter inactive and permits
/// checkerboarding; the constructor warns in that case.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <string>

namespace sparlab {

enum class FilterType {
  None,        ///< no filtering (mesh-dependent, checkerboard-prone; for study use)
  Density,     ///< linear density filter, exact chain rule (default)
  Sensitivity  ///< Sigmund sensitivity filter, heuristic
};

std::string to_string(FilterType type);
FilterType parse_filter_type(const std::string& text);

class DensityFilter {
 public:
  /// \param mesh design mesh (element centroids define the kernel support).
  /// \param type filter variant.
  /// \param radius \f$r_{min}\f$ [m]; must be > 0 unless type is None.
  DensityFilter(const Mesh& mesh, FilterType type, Scalar radius);

  FilterType type() const { return type_; }
  Scalar radius() const { return radius_; }
  Index num_elements() const { return num_elements_; }

  /// Mean number of elements inside the kernel support (diagnostic).
  Scalar average_support() const { return average_support_; }

  /// Row-normalised weight matrix \f$\hat{H}\f$ (empty for FilterType::None).
  const SparseMatrix& weights() const { return h_; }

  /// Physical density from design variables. Identity for None and Sensitivity.
  Vector to_physical(const Vector& x) const;

  /// Exact chain rule \f$\hat{H}^T v\f$. Identity for None and Sensitivity.
  Vector pull_back(const Vector& v) const;

  /// Sigmund sensitivity filter applied to a raw gradient.
  /// \param x design variables, \param dc raw \f$\partial c/\partial x\f$.
  Vector filter_sensitivity(const Vector& x, const Vector& dc) const;

  /// Apply whichever gradient transformation the active filter requires.
  Vector transform_gradient(const Vector& x, const Vector& dc_dphysical) const;

 private:
  void build(const Mesh& mesh);

  FilterType type_;
  Scalar radius_;
  Index num_elements_ = 0;
  SparseMatrix h_;             ///< row-normalised weights
  SparseMatrix h_raw_;         ///< un-normalised weights (sensitivity filter)
  Vector row_sums_;            ///< raw row sums
  Scalar average_support_ = 0.0;
};

}  // namespace sparlab
