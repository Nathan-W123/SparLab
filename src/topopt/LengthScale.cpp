#include "sparlab/topopt/LengthScale.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace sparlab {
namespace {

/// Neighbours (the element itself included) within `radius` of each centroid,
/// found with a uniform bucket grid of cell size `radius`.
std::vector<std::vector<Index>> neighbourhoods(const Mesh& mesh, Scalar radius) {
  const Index ne = mesh.num_elements();
  std::vector<Vector3> c(static_cast<std::size_t>(ne));
  for (Index e = 0; e < ne; ++e) c[static_cast<std::size_t>(e)] = mesh.element_centroid(e);
  const auto cell = [&](Scalar v) { return static_cast<long long>(std::floor(v / radius)); };
  const auto key = [](long long i, long long j, long long k) {
    return (static_cast<std::size_t>(i) * 73856093u) ^ (static_cast<std::size_t>(j) * 19349663u) ^
           (static_cast<std::size_t>(k) * 83492791u);
  };
  std::unordered_map<std::size_t, std::vector<Index>> buckets;
  for (Index e = 0; e < ne; ++e) {
    const Vector3& p = c[static_cast<std::size_t>(e)];
    buckets[key(cell(p.x()), cell(p.y()), cell(p.z()))].push_back(e);
  }
  const Scalar r2 = radius * radius * (1.0 + 1.0e-12);
  std::vector<std::vector<Index>> out(static_cast<std::size_t>(ne));
  for (Index e = 0; e < ne; ++e) {
    const Vector3& p = c[static_cast<std::size_t>(e)];
    const long long ci = cell(p.x()), cj = cell(p.y()), ck = cell(p.z());
    std::vector<Index>& list = out[static_cast<std::size_t>(e)];
    for (long long di = -1; di <= 1; ++di) {
      for (long long dj = -1; dj <= 1; ++dj) {
        for (long long dk = -1; dk <= 1; ++dk) {
          const auto it = buckets.find(key(ci + di, cj + dj, ck + dk));
          if (it == buckets.end()) continue;
          for (Index f : it->second) {
            if ((c[static_cast<std::size_t>(f)] - p).squaredNorm() <= r2) list.push_back(f);
          }
        }
      }
    }
    std::sort(list.begin(), list.end());
    list.erase(std::unique(list.begin(), list.end()), list.end());
  }
  return out;
}

std::vector<char> erode(const std::vector<char>& set, const std::vector<std::vector<Index>>& nbr) {
  std::vector<char> out(set.size(), 0);
  for (std::size_t e = 0; e < set.size(); ++e) {
    if (!set[e]) continue;
    bool all = true;
    for (Index f : nbr[e]) all = all && set[static_cast<std::size_t>(f)];
    out[e] = all ? 1 : 0;
  }
  return out;
}

std::vector<char> dilate(const std::vector<char>& set, const std::vector<std::vector<Index>>& nbr) {
  std::vector<char> out(set.size(), 0);
  for (std::size_t e = 0; e < set.size(); ++e) {
    bool any = false;
    for (Index f : nbr[e]) any = any || set[static_cast<std::size_t>(f)];
    out[e] = any ? 1 : 0;
  }
  return out;
}

}  // namespace

LengthScaleReport check_length_scale(const Mesh& mesh, const Vector& density,
                                     const Vector& element_volumes, Scalar threshold,
                                     Scalar radius) {
  const Index ne = mesh.num_elements();
  if (density.size() != ne || element_volumes.size() != ne) {
    throw ConfigError("length-scale check received vectors of the wrong length");
  }
  if (!(radius > 0.0)) {
    std::ostringstream os;
    os << "the length-scale check needs a positive probe radius (got " << radius << " m)";
    throw ConfigError(os.str());
  }
  LengthScaleReport report;
  report.radius = radius;
  report.threshold = threshold;
  std::vector<char> solid(static_cast<std::size_t>(ne));
  for (Index e = 0; e < ne; ++e) solid[static_cast<std::size_t>(e)] = density(e) >= threshold ? 1 : 0;
  const std::vector<std::vector<Index>> nbr = neighbourhoods(mesh, radius);
  const std::vector<char> opened = dilate(erode(solid, nbr), nbr);
  const std::vector<char> closed = erode(dilate(solid, nbr), nbr);
  Scalar solid_volume = 0.0;
  Scalar void_volume = 0.0;
  Scalar solid_lost = 0.0;
  Scalar void_filled = 0.0;
  for (Index e = 0; e < ne; ++e) {
    const std::size_t i = static_cast<std::size_t>(e);
    const Scalar v = element_volumes(e);
    if (solid[i]) {
      ++report.solid_elements;
      solid_volume += v;
      if (!opened[i]) {
        ++report.solid_violations;
        solid_lost += v;
      }
    } else {
      ++report.void_elements;
      void_volume += v;
      if (closed[i]) {
        ++report.void_violations;
        void_filled += v;
      }
    }
  }
  report.solid_violation_fraction = solid_volume > 0.0 ? solid_lost / solid_volume : 0.0;
  report.void_violation_fraction = void_volume > 0.0 ? void_filled / void_volume : 0.0;
  return report;
}

LengthScaleScan scan_length_scale(const Mesh& mesh, const Vector& density,
                                  const Vector& element_volumes, Scalar threshold,
                                  Scalar step, Scalar max_radius, Scalar tolerance) {
  if (!(step > 0.0) || !(max_radius >= step)) {
    throw ConfigError("the length-scale scan needs a positive step and a cap of at least one "
                      "step");
  }
  LengthScaleScan scan;
  scan.step = step;
  scan.tolerance = tolerance;
  bool solid_open = true;
  bool void_open = true;
  for (int k = 1; step * k <= max_radius * (1.0 + 1.0e-12) && (solid_open || void_open); ++k) {
    const Scalar r = step * k;
    const LengthScaleReport probe =
        check_length_scale(mesh, density, element_volumes, threshold, r);
    scan.probes.push_back(probe);
    if (solid_open) {
      if (probe.solid_violation_fraction <= tolerance) {
        scan.solid_min_size = 2.0 * r;
      } else {
        solid_open = false;
      }
    }
    if (void_open) {
      if (probe.void_violation_fraction <= tolerance) {
        scan.void_min_size = 2.0 * r;
      } else {
        void_open = false;
      }
    }
  }
  scan.solid_bound_reached_cap = solid_open;
  scan.void_bound_reached_cap = void_open;
  return scan;
}

}  // namespace sparlab
