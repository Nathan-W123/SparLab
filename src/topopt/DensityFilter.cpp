#include "sparlab/topopt/DensityFilter.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace sparlab {
namespace {

/// Uniform spatial hash over element centroids with cell size = filter radius,
/// so each query only inspects the 3 x 3 block of neighbouring buckets.
class CentroidHash {
 public:
  CentroidHash(const Mesh& mesh, Scalar cell) : cell_(cell) {
    centroids_.resize(2, mesh.num_elements());
    for (Index e = 0; e < mesh.num_elements(); ++e) {
      centroids_.col(e) = mesh.element_centroid(e);
    }
    for (Index e = 0; e < mesh.num_elements(); ++e) {
      buckets_[key(centroids_(0, e), centroids_(1, e))].push_back(e);
    }
  }

  const Eigen::Matrix2Xd& centroids() const { return centroids_; }

  /// Append every element whose bucket can contain a point within `cell_`.
  void neighbours(Index e, std::vector<Index>& out) const {
    out.clear();
    const long ix = cell_index(centroids_(0, e));
    const long iy = cell_index(centroids_(1, e));
    for (long dy = -1; dy <= 1; ++dy) {
      for (long dx = -1; dx <= 1; ++dx) {
        const auto it = buckets_.find(hash(ix + dx, iy + dy));
        if (it == buckets_.end()) continue;
        out.insert(out.end(), it->second.begin(), it->second.end());
      }
    }
  }

 private:
  long cell_index(Scalar v) const {
    return static_cast<long>(std::floor(v / cell_));
  }
  static std::size_t hash(long ix, long iy) {
    // Cantor-style mixing of two signed integers into one bucket key.
    const std::size_t a = static_cast<std::size_t>(ix * 2654435761L);
    const std::size_t b = static_cast<std::size_t>(iy * 2246822519L);
    return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
  }
  std::size_t key(Scalar x, Scalar y) const {
    return hash(cell_index(x), cell_index(y));
  }

  Scalar cell_;
  Eigen::Matrix2Xd centroids_;
  std::unordered_map<std::size_t, std::vector<Index>> buckets_;
};

}  // namespace

std::string to_string(FilterType type) {
  switch (type) {
    case FilterType::None: return "none";
    case FilterType::Density: return "density";
    case FilterType::Sensitivity: return "sensitivity";
  }
  return "unknown";
}

FilterType parse_filter_type(const std::string& text) {
  std::string lower;
  std::transform(text.begin(), text.end(), std::back_inserter(lower),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lower == "none" || lower == "off") return FilterType::None;
  if (lower == "density") return FilterType::Density;
  if (lower == "sensitivity") return FilterType::Sensitivity;
  throw ConfigError("unknown filter type '" + text +
                    "' (expected none|density|sensitivity)");
}

DensityFilter::DensityFilter(const Mesh& mesh, FilterType type, Scalar radius)
    : type_(type), radius_(radius), num_elements_(mesh.num_elements()) {
  if (type_ == FilterType::None) {
    log::warn(
        "density/sensitivity filtering is disabled; the optimised topology will be "
        "mesh dependent and may show checkerboard patterns");
    return;
  }
  if (!(radius_ > 0.0)) {
    std::ostringstream os;
    os << "filter radius must be positive (got " << radius_
       << " m); use filter.type = none to disable filtering explicitly";
    throw ConfigError(os.str());
  }
  build(mesh);
}

void DensityFilter::build(const Mesh& mesh) {
  CentroidHash hash(mesh, radius_);
  const Eigen::Matrix2Xd& c = hash.centroids();

  TripletList raw;
  raw.reserve(static_cast<std::size_t>(num_elements_) * 12);
  row_sums_.setZero(num_elements_);

  std::vector<Index> candidates;
  std::size_t support_total = 0;
  for (Index e = 0; e < num_elements_; ++e) {
    hash.neighbours(e, candidates);
    Scalar sum = 0.0;
    std::size_t support = 0;
    for (Index i : candidates) {
      const Scalar dist = (c.col(e) - c.col(i)).norm();
      const Scalar w = radius_ - dist;
      if (w <= 0.0) continue;
      raw.emplace_back(e, i, w);
      sum += w;
      ++support;
    }
    if (!(sum > 0.0)) {
      std::ostringstream os;
      os << "filter row for element " << e
         << " is empty, which cannot happen for a positive radius; the mesh may contain "
            "non-finite centroids";
      throw ConfigError(os.str());
    }
    row_sums_(e) = sum;
    support_total += support;
  }
  average_support_ =
      static_cast<Scalar>(support_total) / static_cast<Scalar>(num_elements_);

  h_raw_ = SparseMatrix(num_elements_, num_elements_);
  h_raw_.setFromTriplets(raw.begin(), raw.end());
  h_raw_.makeCompressed();

  // Row-normalised copy.
  TripletList normalised;
  normalised.reserve(raw.size());
  for (const Triplet& t : raw) {
    normalised.emplace_back(t.row(), t.col(), t.value() / row_sums_(t.row()));
  }
  h_ = SparseMatrix(num_elements_, num_elements_);
  h_.setFromTriplets(normalised.begin(), normalised.end());
  h_.makeCompressed();

  if (average_support_ < 2.0) {
    log::warn("filter radius ", radius_,
              " m only reaches ", average_support_,
              " elements on average; a radius below one element size does not suppress "
              "checkerboarding. Increase filter.radius or filter.radius_elements");
  }
  log::debug("built ", to_string(type_), " filter: radius ", radius_, " m, ",
             h_.nonZeros(), " weights, average support ", average_support_,
             " elements");
}

Vector DensityFilter::to_physical(const Vector& x) const {
  if (x.size() != num_elements_) {
    std::ostringstream os;
    os << "design vector has length " << x.size() << " but the filter was built for "
       << num_elements_ << " elements";
    throw ConfigError(os.str());
  }
  if (type_ != FilterType::Density) return x;
  return h_ * x;
}

Vector DensityFilter::pull_back(const Vector& v) const {
  if (v.size() != num_elements_) {
    std::ostringstream os;
    os << "gradient vector has length " << v.size() << " but the filter was built for "
       << num_elements_ << " elements";
    throw ConfigError(os.str());
  }
  if (type_ != FilterType::Density) return v;
  return h_.transpose() * v;
}

Vector DensityFilter::filter_sensitivity(const Vector& x, const Vector& dc) const {
  if (type_ != FilterType::Sensitivity) return dc;
  if (x.size() != num_elements_ || dc.size() != num_elements_) {
    throw ConfigError("sensitivity filter received vectors of the wrong length");
  }
  // gamma guards against division by zero at rho = 0 (Sigmund's 1e-3).
  const Scalar gamma = 1.0e-3;
  const Vector weighted = h_raw_ * (x.cwiseProduct(dc));
  Vector out(num_elements_);
  for (Index e = 0; e < num_elements_; ++e) {
    out(e) = weighted(e) / (std::max(gamma, x(e)) * row_sums_(e));
  }
  return out;
}

Vector DensityFilter::transform_gradient(const Vector& x,
                                         const Vector& dc_dphysical) const {
  switch (type_) {
    case FilterType::None: return dc_dphysical;
    case FilterType::Density: return pull_back(dc_dphysical);
    case FilterType::Sensitivity: return filter_sensitivity(x, dc_dphysical);
  }
  throw ConfigError("unhandled filter type");
}

}  // namespace sparlab
