#include "sparlab/elements/Quadrature.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <array>
#include <cmath>
#include <map>
#include <sstream>

namespace sparlab {
namespace {

/// Abscissae / weights of Gauss-Legendre rules on [-1, 1].
const std::vector<QuadraturePoint1D>& line_rule(int points) {
  static const std::map<int, std::vector<QuadraturePoint1D>> rules = [] {
    std::map<int, std::vector<QuadraturePoint1D>> table;
    table[1] = {{0.0, 2.0}};

    const Scalar a2 = 1.0 / std::sqrt(3.0);
    table[2] = {{-a2, 1.0}, {a2, 1.0}};

    const Scalar a3 = std::sqrt(3.0 / 5.0);
    const Scalar w3a = 5.0 / 9.0;
    const Scalar w3b = 8.0 / 9.0;
    table[3] = {{-a3, w3a}, {0.0, w3b}, {a3, w3a}};

    // 4-point rule: x = sqrt((3 -+ 2 sqrt(6/5))/7)
    const Scalar r = std::sqrt(6.0 / 5.0);
    const Scalar x1 = std::sqrt((3.0 - 2.0 * r) / 7.0);
    const Scalar x2 = std::sqrt((3.0 + 2.0 * r) / 7.0);
    const Scalar s = std::sqrt(30.0);
    const Scalar w1 = (18.0 + s) / 36.0;
    const Scalar w2 = (18.0 - s) / 36.0;
    table[4] = {{-x2, w2}, {-x1, w1}, {x1, w1}, {x2, w2}};
    return table;
  }();

  const auto it = rules.find(points);
  if (it == rules.end()) {
    std::ostringstream os;
    os << "Gauss-Legendre rule with " << points
       << " points is not available (supported: 1, 2, 3, 4)";
    throw ConfigError(os.str());
  }
  return it->second;
}

}  // namespace

const std::vector<QuadraturePoint1D>& gauss_legendre_line(int points) {
  return line_rule(points);
}

const std::vector<QuadraturePoint2D>& gauss_legendre_square(int points_per_direction) {
  static std::map<int, std::vector<QuadraturePoint2D>> cache;
  const auto cached = cache.find(points_per_direction);
  if (cached != cache.end()) return cached->second;

  const auto& rule = line_rule(points_per_direction);  // validates the order
  std::vector<QuadraturePoint2D> tensor;
  tensor.reserve(rule.size() * rule.size());
  for (const auto& gj : rule) {
    for (const auto& gi : rule) {
      QuadraturePoint2D p;
      p.xi = gi.xi;
      p.eta = gj.xi;
      p.weight = gi.weight * gj.weight;
      tensor.push_back(p);
    }
  }
  return cache.emplace(points_per_direction, std::move(tensor)).first->second;
}

const std::vector<QuadraturePoint3D>& gauss_legendre_cube(int points_per_direction) {
  static std::map<int, std::vector<QuadraturePoint3D>> cache;
  const auto cached = cache.find(points_per_direction);
  if (cached != cache.end()) return cached->second;

  const auto& rule = line_rule(points_per_direction);
  std::vector<QuadraturePoint3D> tensor;
  tensor.reserve(rule.size() * rule.size() * rule.size());
  for (const auto& gk : rule) {
    for (const auto& gj : rule) {
      for (const auto& gi : rule) {
        QuadraturePoint3D p;
        p.xi = gi.xi;
        p.eta = gj.xi;
        p.zeta = gk.xi;
        p.weight = gi.weight * gj.weight * gk.weight;
        tensor.push_back(p);
      }
    }
  }
  return cache.emplace(points_per_direction, std::move(tensor)).first->second;
}

}  // namespace sparlab
