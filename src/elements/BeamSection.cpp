#include "sparlab/elements/BeamSection.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

constexpr Scalar kPi = 3.14159265358979323846;

void require_positive(Scalar value, const std::string& what, const std::string& section) {
  if (!(value > 0.0) || !std::isfinite(value)) {
    std::ostringstream os;
    os << "beam section '" << section << "': " << what << " must be positive (got " << value
       << ")";
    throw ConfigError(os.str());
  }
}

void require_coefficient(Scalar k, const std::string& what, const std::string& section) {
  if (!(k > 0.0) || !(k <= 1.0)) {
    std::ostringstream os;
    os << "beam section '" << section << "': the shear coefficient " << what
       << " must lie in (0, 1] (got " << k << "); a section without shear deformation "
          "sets 'shear_deformation': false instead";
    throw ConfigError(os.str());
  }
}

}  // namespace

std::string to_string(BeamSectionShape shape) {
  switch (shape) {
    case BeamSectionShape::General: return "general";
    case BeamSectionShape::Rectangle: return "rectangle";
    case BeamSectionShape::Circle: return "circle";
    case BeamSectionShape::Tube: return "tube";
  }
  return "unknown";
}

BeamSectionShape parse_beam_section_shape(const std::string& text) {
  if (text == "general") return BeamSectionShape::General;
  if (text == "rectangle") return BeamSectionShape::Rectangle;
  if (text == "circle") return BeamSectionShape::Circle;
  if (text == "tube") return BeamSectionShape::Tube;
  throw ConfigError("unknown beam section shape '" + text +
                    "' (expected general, rectangle, circle or tube)");
}

Scalar rectangle_torsion_constant(Scalar width, Scalar height) {
  const Scalar a = std::max(width, height);
  const Scalar b = std::min(width, height);
  // The series converges as n^-5: n up to 49 leaves less than 1e-12.
  Scalar sum = 0.0;
  for (int n = 1; n < 50; n += 2) {
    const Scalar nn = static_cast<Scalar>(n);
    sum += std::tanh(nn * kPi * a / (2.0 * b)) / (nn * nn * nn * nn * nn);
  }
  return a * b * b * b / 3.0 * (1.0 - 192.0 * b / (std::pow(kPi, 5) * a) * sum);
}

Scalar cowper_rectangle(Scalar nu) { return 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu); }

Scalar cowper_circle(Scalar nu) { return 6.0 * (1.0 + nu) / (7.0 + 6.0 * nu); }

Scalar cowper_tube(Scalar nu, Scalar radius_ratio) {
  const Scalar m2 = radius_ratio * radius_ratio;
  const Scalar s = (1.0 + m2) * (1.0 + m2);
  return 6.0 * (1.0 + nu) * s / ((7.0 + 6.0 * nu) * s + (20.0 + 12.0 * nu) * m2);
}

BeamSection resolve(const BeamSection& section, Scalar nu) {
  BeamSection s = section;
  const std::string& name = s.name;
  Scalar cowper = 0.0;
  switch (s.shape) {
    case BeamSectionShape::Rectangle: {
      require_positive(s.width, "the width", name);
      require_positive(s.height, "the height", name);
      s.area = s.width * s.height;
      s.iy = s.width * s.height * s.height * s.height / 12.0;
      s.iz = s.height * s.width * s.width * s.width / 12.0;
      s.torsion = rectangle_torsion_constant(s.width, s.height);
      s.fibre_y = 0.5 * s.width;
      s.fibre_z = 0.5 * s.height;
      cowper = cowper_rectangle(nu);
      break;
    }
    case BeamSectionShape::Circle: {
      require_positive(s.radius, "the radius", name);
      const Scalar r2 = s.radius * s.radius;
      s.area = kPi * r2;
      s.iy = s.iz = 0.25 * kPi * r2 * r2;
      s.torsion = 0.5 * kPi * r2 * r2;
      s.fibre_y = s.fibre_z = s.radius;
      cowper = cowper_circle(nu);
      break;
    }
    case BeamSectionShape::Tube: {
      require_positive(s.radius, "the outer radius", name);
      require_positive(s.inner_radius, "the inner radius", name);
      if (!(s.inner_radius < s.radius)) {
        std::ostringstream os;
        os << "beam section '" << name << "': the inner radius " << s.inner_radius
           << " must be smaller than the outer radius " << s.radius;
        throw ConfigError(os.str());
      }
      const Scalar ro2 = s.radius * s.radius;
      const Scalar ri2 = s.inner_radius * s.inner_radius;
      s.area = kPi * (ro2 - ri2);
      s.iy = s.iz = 0.25 * kPi * (ro2 * ro2 - ri2 * ri2);
      s.torsion = 0.5 * kPi * (ro2 * ro2 - ri2 * ri2);
      s.fibre_y = s.fibre_z = s.radius;
      cowper = cowper_tube(nu, s.inner_radius / s.radius);
      break;
    }
    case BeamSectionShape::General: {
      require_positive(s.area, "the area", name);
      require_positive(s.iy, "I_y", name);
      require_positive(s.iz, "I_z", name);
      require_positive(s.torsion, "the torsion constant", name);
      if (s.fibre_y < 0.0 || s.fibre_z < 0.0) {
        throw ConfigError("beam section '" + name + "': the extreme fibres cannot be negative");
      }
      break;
    }
  }
  if (!s.shear_deformation) {
    if (s.shear_y != 0.0 || s.shear_z != 0.0) {
      throw ConfigError("beam section '" + name +
                        "': a section without shear deformation takes no shear coefficients");
    }
    return s;
  }
  if (s.shape != BeamSectionShape::General) {
    if (s.shear_y == 0.0) s.shear_y = cowper;
    if (s.shear_z == 0.0) s.shear_z = cowper;
  } else if (s.shear_y == 0.0 || s.shear_z == 0.0) {
    throw ConfigError("beam section '" + name +
                      "': a general section states its shear coefficients (k_y, k_z), or "
                      "'shear_deformation': false for the Euler-Bernoulli beam");
  }
  require_coefficient(s.shear_y, "k_y", name);
  require_coefficient(s.shear_z, "k_z", name);
  return s;
}

}  // namespace sparlab
