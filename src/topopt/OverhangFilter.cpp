#include "sparlab/topopt/OverhangFilter.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

namespace sparlab {

std::string BuildDirection::label() const {
  static const char* axes = "xyz";
  std::string out(1, sense > 0 ? '+' : '-');
  out.push_back(axes[axis]);
  return out;
}

BuildDirection parse_build_direction(const std::string& text, int dim) {
  std::string t = text;
  BuildDirection d;
  d.sense = 1;
  if (!t.empty() && (t.front() == '+' || t.front() == '-')) {
    d.sense = t.front() == '-' ? -1 : 1;
    t.erase(t.begin());
  }
  if (t == "x") {
    d.axis = 0;
  } else if (t == "y") {
    d.axis = 1;
  } else if (t == "z") {
    d.axis = 2;
  } else {
    throw ConfigError("build direction '" + text + "' must be one of +x, -x, +y, -y" +
                      std::string(dim == 3 ? ", +z, -z" : ""));
  }
  if (d.axis >= dim) {
    throw ConfigError("build direction '" + text + "' leaves the plane of a 2-D model; use "
                      "+x, -x, +y or -y");
  }
  return d;
}

void OverhangOptions::validate() const {
  if (!(smax_exponent >= 2.0 && smax_exponent <= 200.0)) {
    throw ConfigError("topology.overhang.smax_exponent must lie in [2, 200]");
  }
  if (!(smax_reference > 0.0 && smax_reference < 1.0)) {
    throw ConfigError("topology.overhang.smax_reference must lie in (0, 1)");
  }
  if (!(smin_epsilon > 0.0 && smin_epsilon < 0.1)) {
    throw ConfigError("topology.overhang.smin_epsilon must lie in (0, 0.1)");
  }
}

OverhangFilter::OverhangFilter(const Mesh& mesh, OverhangOptions options,
                               std::vector<char> passive)
    : options_(options), passive_(std::move(passive)) {
  options_.validate();
  const auto& info = mesh.structured_info();
  if (!info.has_value() ||
      (mesh.element_type() != ElementType::Quad4 && mesh.element_type() != ElementType::Hex8)) {
    throw ConfigError("the overhang filter needs a structured Q4 or Hex8 grid "
                      "(mesh.type structured_quad or structured_hex): its supports are the "
                      "grid cells of the layer below");
  }
  const int dim = mesh.dim();
  if (options_.direction.axis >= dim) {
    throw ConfigError("the build direction " + options_.direction.label() +
                      " leaves the plane of a 2-D model");
  }
  const Index n[3] = {info->nx, info->ny, dim == 3 ? info->nz : 1};
  const Scalar h[3] = {info->lx / static_cast<Scalar>(info->nx),
                       info->ly / static_cast<Scalar>(info->ny),
                       dim == 3 ? info->lz / static_cast<Scalar>(info->nz) : 1.0};
  const Index ne = mesh.num_elements();
  if (ne != n[0] * n[1] * n[2]) throw ConfigError("structured grid size does not match the mesh");
  if (!passive_.empty() && static_cast<Index>(passive_.size()) != ne) {
    throw ConfigError("overhang filter: passive flags have the wrong length");
  }
  if (passive_.empty()) passive_.assign(static_cast<std::size_t>(ne), 0);

  const int a = options_.direction.axis;
  const auto index_of = [&](const Index c[3]) {
    return c[0] + n[0] * (c[1] + n[1] * c[2]);
  };
  supports_.assign(static_cast<std::size_t>(ne), {});
  layer_of_.assign(static_cast<std::size_t>(ne), 0);
  std::vector<int> lateral;
  for (int k = 0; k < dim; ++k) {
    if (k != a) lateral.push_back(k);
  }
  for (Index e = 0; e < ne; ++e) {
    Index c[3] = {e % n[0], (e / n[0]) % n[1], e / (n[0] * n[1])};
    const Index layer = options_.direction.sense > 0 ? c[a] : n[a] - 1 - c[a];
    layer_of_[static_cast<std::size_t>(e)] = layer;
    if (layer == 0) continue;
    Index below[3] = {c[0], c[1], c[2]};
    below[a] -= options_.direction.sense;
    std::vector<Index>& s = supports_[static_cast<std::size_t>(e)];
    s.push_back(index_of(below));
    for (int k : lateral) {
      for (int step : {-1, 1}) {
        Index side[3] = {below[0], below[1], below[2]};
        side[k] += step;
        if (side[k] < 0 || side[k] >= n[k]) continue;
        s.push_back(index_of(side));
      }
    }
  }
  order_.resize(static_cast<std::size_t>(ne));
  std::iota(order_.begin(), order_.end(), Index{0});
  std::stable_sort(order_.begin(), order_.end(), [&](Index x, Index y) {
    return layer_of_[static_cast<std::size_t>(x)] < layer_of_[static_cast<std::size_t>(y)];
  });
  q_of_count_.assign(8, options_.smax_exponent);
  std::size_t most_supports = 0;
  for (const std::vector<Index>& s : supports_) most_supports = std::max(most_supports, s.size());
  for (std::size_t count = 1; count < q_of_count_.size(); ++count) {
    q_of_count_[count] = options_.smax_exponent + std::log(static_cast<Scalar>(count)) /
                                                      std::log(options_.smax_reference);
    if (count <= most_supports && !(q_of_count_[count] > 0.0)) {
      std::ostringstream os;
      os << "topology.overhang.smax_exponent " << options_.smax_exponent
         << " is too small for " << count << " supports at smax_reference "
         << options_.smax_reference << ": the exponent Q = P + ln(n)/ln(xi_0) of the "
         << "smooth maximum must stay positive, so P must exceed "
         << std::log(static_cast<Scalar>(count)) / -std::log(options_.smax_reference);
      throw ConfigError(os.str());
    }
  }
  Scalar lateral_size = 0.0;
  for (int k : lateral) lateral_size += h[k];
  lateral_size /= static_cast<Scalar>(lateral.size());
  angle_degrees_ = std::atan(h[a] / lateral_size) * 180.0 / 3.14159265358979323846;
}

OverhangFilter::SmoothMax OverhangFilter::smooth_max(Index e, const Vector& xi) const {
  // Relative to the largest support m, with S = sum_j (xi_j / m)^P in
  // [1, n_s]: smax = m^(P/Q) S^(1/Q). Summing xi_j^P directly underflows
  // for the near-void densities of an optimisation (1e-8^40 = 1e-320), and
  // the derivative's S^(1/Q - 1) then overflows to infinity.
  const std::vector<Index>& s = supports_[static_cast<std::size_t>(e)];
  SmoothMax out;
  for (Index j : s) out.largest = std::max(out.largest, xi(j));
  if (!(out.largest > 0.0)) return out;
  const Scalar p = options_.smax_exponent;
  for (Index j : s) out.scaled_sum += std::pow(std::max(xi(j), 0.0) / out.largest, p);
  out.q = q_of_count_[s.size()];
  out.value = std::pow(out.largest, p / out.q) * std::pow(out.scaled_sum, 1.0 / out.q);
  return out;
}

OverhangFilter::Forward OverhangFilter::forward(const Vector& input) const {
  if (input.size() != num_elements()) {
    throw ConfigError("overhang filter applied to a density of the wrong length");
  }
  Forward out;
  out.xi = input;
  out.support_max = Vector::Zero(input.size());
  const Scalar eps = options_.smin_epsilon;
  const Scalar root_eps = std::sqrt(eps);
  for (Index e : order_) {
    if (layer_of_[static_cast<std::size_t>(e)] == 0 || passive_[static_cast<std::size_t>(e)]) {
      continue;
    }
    const Scalar big_xi = smooth_max(e, out.xi).value;
    out.support_max(e) = big_xi;
    const Scalar x = input(e);
    const Scalar d = x - big_xi;
    out.xi(e) = 0.5 * (x + big_xi - std::sqrt(d * d + eps) + root_eps);
  }
  return out;
}

Vector OverhangFilter::apply(const Vector& input) const { return forward(input).xi; }

Vector OverhangFilter::pull_back(const Vector& input, const Vector& df_dxi) const {
  if (df_dxi.size() != num_elements()) {
    throw ConfigError("overhang filter pull-back received a gradient of the wrong length");
  }
  const Forward fw = forward(input);
  const Scalar eps = options_.smin_epsilon;
  const Scalar p = options_.smax_exponent;
  Vector lambda = df_dxi;
  Vector out = Vector::Zero(num_elements());
  for (auto it = order_.rbegin(); it != order_.rend(); ++it) {
    const Index e = *it;
    if (layer_of_[static_cast<std::size_t>(e)] == 0 || passive_[static_cast<std::size_t>(e)]) {
      out(e) = lambda(e);
      continue;
    }
    const Scalar d = input(e) - fw.support_max(e);
    const Scalar root = std::sqrt(d * d + eps);
    out(e) = lambda(e) * 0.5 * (1.0 - d / root);
    const Scalar through = lambda(e) * 0.5 * (1.0 + d / root);
    if (through == 0.0) continue;
    // d smax / d xi_j = (P/Q) m^(P/Q - 1) S^(1/Q - 1) (xi_j / m)^(P - 1);
    // P/Q >= 1, so the power of m stays bounded as m -> 0.
    const SmoothMax sm = smooth_max(e, fw.xi);
    if (!(sm.largest > 0.0)) continue;
    const Scalar common = through * (p / sm.q) * std::pow(sm.largest, p / sm.q - 1.0) *
                          std::pow(sm.scaled_sum, 1.0 / sm.q - 1.0);
    for (Index j : supports_[static_cast<std::size_t>(e)]) {
      lambda(j) += common * std::pow(std::max(fw.xi(j), 0.0) / sm.largest, p - 1.0);
    }
  }
  return out;
}

OverhangReport check_overhang(const OverhangFilter& stencil, const Vector& density,
                              const Vector& element_volumes, Scalar threshold) {
  const Index ne = stencil.num_elements();
  if (density.size() != ne || element_volumes.size() != ne) {
    throw ConfigError("overhang check received vectors of the wrong length");
  }
  OverhangReport report;
  report.build_direction = stencil.options().direction.label();
  report.threshold = threshold;
  report.overhang_angle_degrees = stencil.overhang_angle_degrees();
  for (Index e = 0; e < ne; ++e) {
    if (density(e) < threshold) continue;
    ++report.solid_elements;
    report.solid_volume += element_volumes(e);
    if (stencil.layer(e) == 0) continue;
    bool supported = false;
    for (Index s : stencil.supports(e)) supported = supported || density(s) >= threshold;
    if (supported) continue;
    ++report.unsupported_elements;
    report.unsupported_volume += element_volumes(e);
    if (report.lowest_unsupported_layer < 0 ||
        stencil.layer(e) < report.lowest_unsupported_layer) {
      report.lowest_unsupported_layer = stencil.layer(e);
    }
  }
  report.unsupported_fraction =
      report.solid_volume > 0.0 ? report.unsupported_volume / report.solid_volume : 0.0;
  return report;
}

}  // namespace sparlab
