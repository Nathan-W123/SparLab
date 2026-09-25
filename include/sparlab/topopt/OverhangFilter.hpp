/// \file OverhangFilter.hpp
/// \brief Additive-manufacturing overhang filter (Langelaar 2016, 2017) and
///        the overhang check of a thresholded design.
///
/// A part printed layer by layer can only deposit material where the layer
/// below supports it. On a structured grid with square (cubic) cells whose
/// layers are normal to the build direction, an element is supported when
/// solid material lies directly below it or diagonally below it: 3 elements
/// in 2-D, a cross of 5 in 3-D (the element below and its four face
/// neighbours in that layer). This encodes an overhang limit of
/// \f$\arctan(\Delta_{build}/\Delta_{lateral})\f$ from the build plate,
/// 45 degrees for square cells.
///
/// **Filter.** Layer by layer from the build plate, the printable density is
/// \f[
///   \xi_e = \mathrm{smin}\big(x_e,\ \Xi_e\big), \qquad
///   \Xi_e = \mathrm{smax}_{s \in S(e)} \xi_s ,
/// \f]
/// with \f$\xi_e = x_e\f$ on the first layer (the plate supports it), the
/// smooth maximum \f$\mathrm{smax}(\xi) = (\sum_s \xi_s^P)^{1/Q}\f$ with
/// \f$Q = P + \ln n_s / \ln \xi_0\f$ (so \f$n_s\f$ equal values
/// \f$\xi_0 = 0.5\f$ return \f$\xi_0\f$ exactly) and the smooth minimum
/// \f$\mathrm{smin}(x, \Xi) = \tfrac12\big(x + \Xi - \sqrt{(x-\Xi)^2+\epsilon}
/// + \sqrt\epsilon\big)\f$. An element can therefore hold at most as much
/// material as its supports. Passive elements keep their input density: a
/// passive solid pad is taken as fixture, and the check below reports it if
/// it overhangs.
///
/// **Sensitivities.** \f$\xi\f$ depends on every input below it through the
/// supports, so the chain rule runs as an adjoint recursion from the top
/// layer down: with \f$\lambda_e = \partial f/\partial\xi_e + \sum_{c:\,e\in
/// S(c)} \lambda_c\,\partial\xi_c/\partial\Xi_c\,\partial\Xi_c/\partial\xi_e\f$,
/// \f$\partial f/\partial x_e = \lambda_e\,\partial\xi_e/\partial x_e\f$ -
/// one pass over the elements, exact, and checked against central
/// differences in the tests.
///
/// In the optimiser the filter sits between the density filter and the
/// Heaviside projection: \f$x \to \tilde\rho = \hat H x \to \xi =
/// \mathrm{AM}(\tilde\rho) \to \bar\rho = H_\beta(\xi)\f$.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Build direction of an additive process: the axis and its sense.
struct BuildDirection {
  int axis = 1;       ///< 0 = x, 1 = y, 2 = z
  int sense = 1;      ///< +1: the plate is at the low end of the axis; -1: the high end
  std::string label() const;
};

/// Parse "+x", "-y", "z" (= "+z"), ...
/// \throws ConfigError for anything else.
BuildDirection parse_build_direction(const std::string& text, int dim);

struct OverhangOptions {
  bool filter = false;           ///< apply the filter in the optimisation
  BuildDirection direction;      ///< the build direction
  Scalar smax_exponent = 40.0;   ///< P of the smooth maximum
  Scalar smax_reference = 0.5;   ///< xi_0 at which smax is exact
  Scalar smin_epsilon = 1.0e-4;  ///< epsilon of the smooth minimum
  void validate() const;
};

class OverhangFilter {
 public:
  /// \param mesh a structured Q4 or Hex8 mesh (StructuredGridInfo present).
  /// \param passive per-element flags; flagged elements keep their density.
  /// \throws ConfigError for another mesh or an out-of-plane direction.
  OverhangFilter(const Mesh& mesh, OverhangOptions options,
                 std::vector<char> passive = {});

  const OverhangOptions& options() const { return options_; }
  Index num_elements() const { return static_cast<Index>(layer_of_.size()); }

  /// Printable density \f$\xi\f$ of an input density.
  Vector apply(const Vector& input) const;

  /// \f$\partial f/\partial x\f$ from \f$\partial f/\partial\xi\f$, at the
  /// input `input` (the adjoint recursion).
  Vector pull_back(const Vector& input, const Vector& df_dxi) const;

  /// Supporting elements of `e` (empty on the first layer).
  const std::vector<Index>& supports(Index e) const {
    return supports_[static_cast<std::size_t>(e)];
  }
  /// Layer of `e` along the build direction, 0 on the build plate.
  Index layer(Index e) const { return layer_of_[static_cast<std::size_t>(e)]; }

  /// Overhang angle the support stencil encodes [degrees]:
  /// atan(layer thickness / lateral cell size).
  Scalar overhang_angle_degrees() const { return angle_degrees_; }

 private:
  struct Forward {
    Vector xi;
    Vector support_max;  ///< Xi per element (0 on the first layer)
  };
  /// The smooth maximum over the supports of `e`, with the terms of its
  /// derivative: the largest support m and S = sum_j (xi_j / m)^P.
  struct SmoothMax {
    Scalar value = 0.0;
    Scalar largest = 0.0;
    Scalar scaled_sum = 0.0;
    Scalar q = 1.0;
  };
  Forward forward(const Vector& input) const;
  SmoothMax smooth_max(Index e, const Vector& xi) const;

  OverhangOptions options_;
  std::vector<std::vector<Index>> supports_;
  std::vector<Index> layer_of_;
  std::vector<Index> order_;        ///< elements sorted by layer, bottom first
  std::vector<char> passive_;
  std::vector<Scalar> q_of_count_;  ///< Q for a support count
  Scalar angle_degrees_ = 45.0;
};

/// Result of checking a thresholded design for unsupported material.
struct OverhangReport {
  std::string build_direction;
  Scalar threshold = 0.5;
  Scalar overhang_angle_degrees = 45.0;
  Index solid_elements = 0;
  Index unsupported_elements = 0;    ///< solid, off the plate, no solid support
  Scalar solid_volume = 0.0;         ///< [m^3]
  Scalar unsupported_volume = 0.0;   ///< [m^3]
  Scalar unsupported_fraction = 0.0; ///< of the solid volume
  Index lowest_unsupported_layer = -1;
};

/// Threshold `density` and count the solid elements none of whose supports is
/// solid (the first layer rests on the plate). Uses the same stencil as the
/// filter.
OverhangReport check_overhang(const OverhangFilter& stencil, const Vector& density,
                              const Vector& element_volumes, Scalar threshold);

}  // namespace sparlab
