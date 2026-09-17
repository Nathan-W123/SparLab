#include "sparlab/topopt/SimpInterpolation.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sparlab {

std::string to_string(MassInterpolation law) {
  switch (law) {
    case MassInterpolation::Linear: return "linear";
    case MassInterpolation::PenaltyMatched: return "penalty_matched";
  }
  return "unknown";
}

MassInterpolation parse_mass_interpolation(const std::string& text) {
  std::string lower;
  std::transform(text.begin(), text.end(), std::back_inserter(lower),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lower == "linear") return MassInterpolation::Linear;
  if (lower == "penalty_matched" || lower == "penalised" || lower == "penalized") {
    return MassInterpolation::PenaltyMatched;
  }
  throw ConfigError("unknown mass interpolation '" + text +
                    "' (expected linear|penalty_matched)");
}

void SimpOptions::validate() const {
  if (!(penalty >= 1.0)) {
    std::ostringstream os;
    os << "SIMP penalty must be >= 1 (got " << penalty
       << "); p < 1 makes intermediate densities favourable and the 0/1 design "
          "unreachable";
    throw ConfigError(os.str());
  }
  if (!(emin_ratio > 0.0 && emin_ratio < 1.0)) {
    std::ostringstream os;
    os << "SIMP emin_ratio must lie in (0, 1) (got " << emin_ratio
       << "); a zero floor makes the stiffness matrix singular in void regions";
    throw ConfigError(os.str());
  }
  if (!(mass_floor > 0.0 && mass_floor < 1.0)) {
    std::ostringstream os;
    os << "SIMP mass_floor must lie in (0, 1) (got " << mass_floor << ")";
    throw ConfigError(os.str());
  }
}

Scalar simp_stiffness_factor(Scalar rho, const SimpOptions& options) {
  const Scalar r = std::clamp(rho, 0.0, 1.0);
  return options.emin_ratio + (1.0 - options.emin_ratio) * std::pow(r, options.penalty);
}

Scalar simp_stiffness_derivative(Scalar rho, const SimpOptions& options) {
  const Scalar r = std::clamp(rho, 0.0, 1.0);
  // For p > 1 the derivative vanishes at rho = 0, which is the intended
  // behaviour of the modified SIMP law.
  return options.penalty * (1.0 - options.emin_ratio) *
         std::pow(r, options.penalty - 1.0);
}

Scalar simp_mass_factor(Scalar rho, const SimpOptions& options) {
  const Scalar r = std::clamp(rho, 0.0, 1.0);
  switch (options.mass_law) {
    case MassInterpolation::Linear:
      return r;
    case MassInterpolation::PenaltyMatched:
      return options.mass_floor +
             (1.0 - options.mass_floor) * std::pow(r, options.penalty);
  }
  throw ConfigError("unhandled mass interpolation law");
}

Vector simp_stiffness_factors(const Vector& rho, const SimpOptions& options) {
  Vector out(rho.size());
  for (Eigen::Index e = 0; e < rho.size(); ++e) {
    out(e) = simp_stiffness_factor(rho(e), options);
  }
  return out;
}

Vector simp_stiffness_derivatives(const Vector& rho, const SimpOptions& options) {
  Vector out(rho.size());
  for (Eigen::Index e = 0; e < rho.size(); ++e) {
    out(e) = simp_stiffness_derivative(rho(e), options);
  }
  return out;
}

Vector simp_mass_factors(const Vector& rho, const SimpOptions& options) {
  Vector out(rho.size());
  for (Eigen::Index e = 0; e < rho.size(); ++e) {
    out(e) = simp_mass_factor(rho(e), options);
  }
  return out;
}

}  // namespace sparlab
