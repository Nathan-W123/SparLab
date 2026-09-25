/// \file LengthScale.hpp
/// \brief Minimum length-scale check of a thresholded design by morphological
///        opening and closing on the element centroids.
///
/// A design is thresholded into a solid set \f$S\f$, and a ball of radius
/// \f$r\f$ probes it: the *opening* \f$(S \ominus B_r) \oplus B_r\f$ keeps
/// every solid element a ball of radius \f$r\f$ inside the solid can reach, so
/// the solid elements it removes lie in members thinner than about \f$2r\f$;
/// the *closing* \f$(S \oplus B_r) \ominus B_r\f$ fills every void element a
/// ball of radius \f$r\f$ inside the void cannot reach, i.e. gaps narrower
/// than about \f$2r\f$. Erosion and dilation act on element centroids with
/// the neighbours inside the mesh only, so the domain boundary neither erodes
/// a member lying against it nor closes a gap that opens onto it. The
/// fractions reported are by volume and resolve a length to about one
/// element: they measure what a design delivers, whereas the robust
/// formulation is what is meant to guarantee it. Opening also trims convex
/// corners, whose cells no ball of the probe's size reaches, so a probe
/// passes while it removes at most a tolerance of the volume: a pixelated
/// corner costs a few cells, a member thinner than the probe all of its.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

namespace sparlab {

struct LengthScaleReport {
  Scalar radius = 0.0;                    ///< probe ball radius r [m]
  Scalar threshold = 0.5;
  Index solid_elements = 0;
  Index solid_violations = 0;             ///< removed by the opening
  Scalar solid_violation_fraction = 0.0;  ///< of the solid volume
  Index void_elements = 0;
  Index void_violations = 0;              ///< filled by the closing
  Scalar void_violation_fraction = 0.0;   ///< of the void volume
};

/// \throws ConfigError for a non-positive radius or mismatched lengths.
LengthScaleReport check_length_scale(const Mesh& mesh, const Vector& density,
                                     const Vector& element_volumes, Scalar threshold,
                                     Scalar radius);

/// Probes of increasing radius and the minimum sizes they bound.
struct LengthScaleScan {
  std::vector<LengthScaleReport> probes;  ///< radii step, 2 step, ... up to the cap
  Scalar step = 0.0;                      ///< radius increment [m]
  Scalar tolerance = 0.02;                ///< volume fraction a probe may flag
  /// About the narrowest solid member: twice the largest probe radius whose
  /// opening removes at most `tolerance` of the solid volume [m]; zero when
  /// even the smallest probe removes more.
  Scalar solid_min_size = 0.0;
  /// The same for the void gaps, from the closing.
  Scalar void_min_size = 0.0;
  bool solid_bound_reached_cap = false;   ///< no probe up to the cap failed
  bool void_bound_reached_cap = false;
};

/// Scan probe radii step, 2 step, ..., max_radius and report the minimum
/// member and gap sizes the design keeps, each to within one step.
LengthScaleScan scan_length_scale(const Mesh& mesh, const Vector& density,
                                  const Vector& element_volumes, Scalar threshold,
                                  Scalar step, Scalar max_radius, Scalar tolerance = 0.02);

}  // namespace sparlab
