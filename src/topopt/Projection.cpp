#include "sparlab/topopt/Projection.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <cmath>
#include <sstream>

namespace sparlab {

void ProjectionOptions::validate() const {
  if (!enabled) return;
  std::ostringstream os;
  if (!(eta > 0.0 && eta < 1.0)) {
    os << "topology.projection.eta must lie in (0, 1), got " << eta;
  } else if (!(beta_start > 0.0) || !std::isfinite(beta_start)) {
    os << "topology.projection.beta_start must be positive, got " << beta_start;
  } else if (!(beta_max >= beta_start) || !std::isfinite(beta_max)) {
    os << "topology.projection.beta_max (" << beta_max
       << ") must be at least beta_start (" << beta_start << ")";
  } else if (beta_max > 1024.0) {
    os << "topology.projection.beta_max = " << beta_max
       << " is beyond any useful sharpness; the derivative becomes a spike that no "
          "optimiser can follow. Use at most a few hundred (32 to 64 is usual)";
  } else if (!(beta_factor > 1.0)) {
    os << "topology.projection.beta_factor must exceed 1, got " << beta_factor;
  } else if (beta_interval < 1) {
    os << "topology.projection.beta_interval must be at least 1, got " << beta_interval;
  } else {
    return;
  }
  throw ConfigError(os.str());
}

int ProjectionOptions::num_stages() const {
  int stages = 1;
  Scalar beta = beta_start;
  while (beta < beta_max && stages < 64) {
    beta *= beta_factor;
    ++stages;
  }
  return stages;
}

Scalar ProjectionOptions::beta_of_stage(int k) const {
  Scalar beta = beta_start;
  for (int i = 0; i < k && beta < beta_max; ++i) beta *= beta_factor;
  return std::min(beta, beta_max);
}

Scalar heaviside_project(Scalar filtered, Scalar beta, Scalar eta) {
  const Scalar a = std::tanh(beta * eta);
  const Scalar denom = a + std::tanh(beta * (1.0 - eta));
  return (a + std::tanh(beta * (filtered - eta))) / denom;
}

Scalar heaviside_derivative(Scalar filtered, Scalar beta, Scalar eta) {
  const Scalar denom = std::tanh(beta * eta) + std::tanh(beta * (1.0 - eta));
  const Scalar t = std::tanh(beta * (filtered - eta));
  return beta * (1.0 - t * t) / denom;
}

Vector heaviside_project(const Vector& filtered, Scalar beta, Scalar eta) {
  Vector out(filtered.size());
  for (Eigen::Index i = 0; i < filtered.size(); ++i) {
    out(i) = heaviside_project(filtered(i), beta, eta);
  }
  return out;
}

Vector heaviside_derivative(const Vector& filtered, Scalar beta, Scalar eta) {
  Vector out(filtered.size());
  for (Eigen::Index i = 0; i < filtered.size(); ++i) {
    out(i) = heaviside_derivative(filtered(i), beta, eta);
  }
  return out;
}

}  // namespace sparlab
