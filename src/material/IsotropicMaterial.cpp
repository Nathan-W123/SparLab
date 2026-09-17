#include "sparlab/material/IsotropicMaterial.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <sstream>
#include <utility>

namespace sparlab {

IsotropicMaterial::IsotropicMaterial(Scalar youngs_modulus, Scalar poisson_ratio,
                                     Scalar density, std::string name)
    : e_(youngs_modulus), nu_(poisson_ratio), rho_(density), name_(std::move(name)) {
  if (!(e_ > 0.0)) {
    std::ostringstream os;
    os << "material '" << name_ << "' needs a positive Young's modulus (got " << e_
       << " Pa)";
    throw ConfigError(os.str());
  }
  if (!(nu_ > -1.0 && nu_ < 0.5)) {
    std::ostringstream os;
    os << "material '" << name_ << "' needs a Poisson ratio in (-1, 0.5) for a "
       << "positive-definite 2-D constitutive matrix (got " << nu_ << ")";
    throw ConfigError(os.str());
  }
  if (!(rho_ >= 0.0)) {
    std::ostringstream os;
    os << "material '" << name_ << "' needs a non-negative density (got " << rho_
       << " kg/m^3)";
    throw ConfigError(os.str());
  }
}

Matrix3 IsotropicMaterial::plane_stress_matrix() const {
  Matrix3 d = Matrix3::Zero();
  const Scalar f = e_ / (1.0 - nu_ * nu_);
  d(0, 0) = f;
  d(0, 1) = f * nu_;
  d(1, 0) = f * nu_;
  d(1, 1) = f;
  d(2, 2) = f * (1.0 - nu_) / 2.0;
  return d;
}

Matrix3 IsotropicMaterial::plane_strain_matrix() const {
  Matrix3 d = Matrix3::Zero();
  const Scalar f = e_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
  d(0, 0) = f * (1.0 - nu_);
  d(0, 1) = f * nu_;
  d(1, 0) = f * nu_;
  d(1, 1) = f * (1.0 - nu_);
  d(2, 2) = f * (1.0 - 2.0 * nu_) / 2.0;
  return d;
}

Matrix3 IsotropicMaterial::constitutive(StressState state) const {
  switch (state) {
    case StressState::PlaneStress: return plane_stress_matrix();
    case StressState::PlaneStrain: return plane_strain_matrix();
  }
  throw ConfigError("unhandled 2-D stress state");
}

}  // namespace sparlab
