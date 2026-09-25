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

const std::vector<QuadraturePoint3D>& tetrahedron_rule_4() {
  // Built once, on first use (a thread-safe static initialisation), because
  // the simplex kernels run inside parallel loops.
  static const std::vector<QuadraturePoint3D> rule = [] {
    const Scalar alpha = (5.0 + 3.0 * std::sqrt(5.0)) / 20.0;
    const Scalar beta = (5.0 - std::sqrt(5.0)) / 20.0;
    const Scalar w = 1.0 / 24.0;
    // (xi, eta, zeta) are the barycentric coordinates L1, L2, L3; the first
    // point has L0 = alpha.
    return std::vector<QuadraturePoint3D>{{beta, beta, beta, w},
                                          {alpha, beta, beta, w},
                                          {beta, alpha, beta, w},
                                          {beta, beta, alpha, w}};
  }();
  return rule;
}

namespace {

/// Gauss-Legendre points and weights mapped from [-1, 1] to [0, 1].
std::vector<QuadraturePoint1D> unit_interval_rule(int n) {
  std::vector<QuadraturePoint1D> out;
  for (const QuadraturePoint1D& p : line_rule(n)) {
    out.push_back({0.5 * (p.xi + 1.0), 0.5 * p.weight});
  }
  return out;
}

void check_collapsed_order(int n, const char* shape) {
  if (n < 1 || n > 4) {
    std::ostringstream os;
    os << "collapsed Gauss rule on the " << shape << " with " << n
       << " points per direction is not available (supported: 1 to 4)";
    throw ConfigError(os.str());
  }
}

}  // namespace

const std::vector<QuadraturePoint3D>& collapsed_gauss_tetrahedron(int n) {
  check_collapsed_order(n, "tetrahedron");
  static const std::array<std::vector<QuadraturePoint3D>, 5> rules = [] {
    std::array<std::vector<QuadraturePoint3D>, 5> table;
    for (int order = 1; order <= 4; ++order) {
      const std::vector<QuadraturePoint1D> g = unit_interval_rule(order);
      for (const QuadraturePoint1D& pu : g) {
        for (const QuadraturePoint1D& pv : g) {
          for (const QuadraturePoint1D& pw : g) {
            QuadraturePoint3D p;
            p.xi = pu.xi;
            p.eta = pv.xi * (1.0 - pu.xi);
            p.zeta = pw.xi * (1.0 - pu.xi) * (1.0 - pv.xi);
            p.weight = pu.weight * pv.weight * pw.weight * (1.0 - pu.xi) * (1.0 - pu.xi) *
                       (1.0 - pv.xi);
            table[static_cast<std::size_t>(order)].push_back(p);
          }
        }
      }
    }
    return table;
  }();
  return rules[static_cast<std::size_t>(n)];
}

const std::vector<QuadraturePoint2D>& collapsed_gauss_triangle(int n) {
  check_collapsed_order(n, "triangle");
  static const std::array<std::vector<QuadraturePoint2D>, 5> rules = [] {
    std::array<std::vector<QuadraturePoint2D>, 5> table;
    for (int order = 1; order <= 4; ++order) {
      const std::vector<QuadraturePoint1D> g = unit_interval_rule(order);
      for (const QuadraturePoint1D& pu : g) {
        for (const QuadraturePoint1D& pv : g) {
          QuadraturePoint2D p;
          p.xi = pu.xi;
          p.eta = pv.xi * (1.0 - pu.xi);
          p.weight = pu.weight * pv.weight * (1.0 - pu.xi);
          table[static_cast<std::size_t>(order)].push_back(p);
        }
      }
    }
    return table;
  }();
  return rules[static_cast<std::size_t>(n)];
}

}  // namespace sparlab
