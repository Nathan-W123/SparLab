/// \file BeamSection.hpp
/// \brief Cross-sections of the Timoshenko beam (Beam2.hpp): their area,
///        second moments, torsion constant, shear coefficients and extreme
///        fibres, given directly or computed from a shape.
///
/// The section lies in the beam's local (y', z') plane, its centroid on the
/// beam's axis and y', z' its principal axes: \f$I_y = \int z'^2 dA\f$ resists
/// bending in the x'-z' plane (about y'), \f$I_z = \int y'^2 dA\f$ bending in
/// the x'-y' plane (about z'). A shape's width runs along y', its height
/// along z'.
///
/// Shear coefficients \f$k\f$ (the shear area \f$k A\f$) of the shapes are
/// Cowper's (1966), which depend on Poisson's ratio: a rectangle
/// \f$10(1+\nu)/(12+11\nu)\f$, a solid circle \f$6(1+\nu)/(7+6\nu)\f$, a tube
/// of radii \f$r_i < r_o\f$, \f$m = r_i / r_o\f$,
/// \f$6(1+\nu)(1+m^2)^2 / ((7+6\nu)(1+m^2)^2 + (20+12\nu)m^2)\f$. The torsion
/// constant is Saint-Venant's, the section's warping free: \f$\pi r^4 / 2\f$
/// for a circle, \f$\pi (r_o^4 - r_i^4)/2\f$ for a tube, and for a
/// rectangle \f$a \times b\f$, \f$a \ge b\f$, the series
/// \f$\tfrac{a b^3}{3}\left[1 - \tfrac{192 b}{\pi^5 a}\sum_{n\ \mathrm{odd}}
/// \tfrac{\tanh(n\pi a / 2b)}{n^5}\right]\f$.
#pragma once

#include "sparlab/core/Types.hpp"

#include <string>

namespace sparlab {

enum class BeamSectionShape {
  General,    ///< area, second moments, torsion constant given directly
  Rectangle,  ///< width (along y') x height (along z')
  Circle,     ///< solid circle of a radius
  Tube        ///< circular tube of an outer and an inner radius
};

std::string to_string(BeamSectionShape shape);
/// \throws ConfigError for an unknown name.
BeamSectionShape parse_beam_section_shape(const std::string& text);

/// A beam cross-section. A shape's properties follow from its dimensions
/// (`resolve`); a general section states them.
struct BeamSection {
  std::string name;
  BeamSectionShape shape = BeamSectionShape::General;
  Scalar width = 0.0;         ///< rectangle, along y' [m]
  Scalar height = 0.0;        ///< rectangle, along z' [m]
  Scalar radius = 0.0;        ///< circle; the tube's outer radius [m]
  Scalar inner_radius = 0.0;  ///< tube [m]
  Scalar area = 0.0;          ///< A [m^2]
  Scalar iy = 0.0;            ///< I_y = int z'^2 dA [m^4]
  Scalar iz = 0.0;            ///< I_z = int y'^2 dA [m^4]
  Scalar torsion = 0.0;       ///< Saint-Venant torsion constant J [m^4]
  /// Shear coefficients k_y, k_z (the shear areas k A along y' and z'); a
  /// shape's are Cowper's when left at 0. Both 0 on a general section with
  /// `shear_deformation` false: the Euler-Bernoulli beam.
  Scalar shear_y = 0.0;
  Scalar shear_z = 0.0;
  /// False: no shear deformation (the Euler-Bernoulli beam).
  bool shear_deformation = true;
  /// Largest |y'| and |z'| of the section [m], for its extreme-fibre
  /// stress; 0 when unknown (a general section may state them).
  Scalar fibre_y = 0.0;
  Scalar fibre_z = 0.0;
  /// Direction of the local y' axis, approximately (its component along
  /// the beam's axis is dropped); zero for the default (Beam2.hpp).
  Vector3 orientation = Vector3::Zero();

  /// Polar moment I_y + I_z [m^4].
  Scalar polar() const { return iy + iz; }
  /// True for a shape whose extreme fibre is a radius (circle, tube).
  bool round() const { return shape == BeamSectionShape::Circle || shape == BeamSectionShape::Tube; }
};

/// The section's properties for Poisson's ratio `nu`: a shape's area, second
/// moments, torsion constant and extreme fibres from its dimensions, and its
/// Cowper shear coefficients unless set; a general section as given, checked.
/// \throws ConfigError for a non-positive dimension or property, an inner
///         radius not below the outer one, a shear coefficient outside
///         (0, 1] - or not 0 on a section without shear deformation.
BeamSection resolve(const BeamSection& section, Scalar nu);

/// Saint-Venant torsion constant of a solid rectangle [m^4].
Scalar rectangle_torsion_constant(Scalar width, Scalar height);

/// Cowper's shear coefficients.
/// \{
Scalar cowper_rectangle(Scalar nu);
Scalar cowper_circle(Scalar nu);
Scalar cowper_tube(Scalar nu, Scalar radius_ratio);
/// \}

}  // namespace sparlab
