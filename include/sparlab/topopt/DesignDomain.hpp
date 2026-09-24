/// \file DesignDomain.hpp
/// \brief Design variables, their bounds, passive regions and volume data.
///
/// Every element carries one design variable \f$x_e \in [0,1]\f$. Passive
/// regions are imposed as *equal* lower and upper bounds on the design
/// variable, so the optimiser never moves them and the objective remains an
/// exact function of the free variables (which is what makes the
/// finite-difference verification meaningful).
///
/// The volume constraint is written on the **physical** density,
/// \f[
///   g(x) = \sum_e \tilde{\rho}_e(x)\, v_e - \nu\, V_{domain} \le 0,
/// \f]
/// with \f$V_{domain} = \sum_e v_e\f$ the volume of the *whole* design domain,
/// including passive regions. A passive solid region therefore consumes part of
/// the volume budget, which is the correct reading of a mass constraint on a
/// real part. Feasibility (passive solid volume below the target) is checked up
/// front.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/Selector.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// A region whose density is held fixed.
struct PassiveRegionSpec {
  SelectorGroup region;
  /// true: solid (x = 1); false: void (x = `void_density`).
  bool solid = true;
  /// Density imposed on a void region. Kept at 0 by default: the modified SIMP
  /// law supplies the stiffness floor, so x = 0 is safe.
  Scalar void_density = 0.0;
};

enum class PassiveTag : unsigned char { Free = 0, Solid = 1, Void = 2 };

class DesignDomain {
 public:
  /// \param model finalised FE model providing the mesh and element volumes.
  /// \param volume_fraction \f$\nu\f$ [-] in (0, 1].
  /// \param initial_density starting value for free variables; when negative,
  ///        `volume_fraction` is used (the standard uniform start).
  /// \param passive passive region specifications.
  DesignDomain(const FemModel& model, Scalar volume_fraction, Scalar initial_density,
               const std::vector<PassiveRegionSpec>& passive);

  Index num_elements() const { return num_elements_; }
  Index num_free_variables() const { return num_free_; }
  Index num_passive_solid() const { return num_solid_; }
  Index num_passive_void() const { return num_void_; }

  const Vector& lower_bounds() const { return lower_; }
  const Vector& upper_bounds() const { return upper_; }
  const Vector& initial_design() const { return initial_; }
  const Vector& element_volumes() const { return volumes_; }
  const std::vector<PassiveTag>& tags() const { return tags_; }

  bool is_free(Index e) const { return tags_[static_cast<std::size_t>(e)] == PassiveTag::Free; }

  Scalar domain_volume() const { return domain_volume_; }
  Scalar volume_fraction() const { return volume_fraction_; }
  /// Starting value of the free variables: `initial_density` when one was
  /// given, otherwise the volume fraction.
  Scalar initial_density() const { return initial_density_; }
  Scalar volume_target() const { return volume_fraction_ * domain_volume_; }
  Scalar passive_solid_volume() const { return solid_volume_; }
  Scalar passive_void_volume() const { return void_volume_; }

  /// Physical volume of a density field, \f$\sum_e \rho_e v_e\f$ [m^3].
  Scalar volume_of(const Vector& density) const { return density.dot(volumes_); }

  /// Volume fraction of a density field [-].
  Scalar fraction_of(const Vector& density) const {
    return volume_of(density) / domain_volume_;
  }

  /// Clamp a design vector to its bounds (in place).
  void clamp(Vector& x) const;

  /// Human-readable description used in run summaries.
  std::string describe() const;

 private:
  Index num_elements_ = 0;
  Index num_free_ = 0;
  Index num_solid_ = 0;
  Index num_void_ = 0;
  Scalar volume_fraction_ = 0.5;
  Scalar initial_density_ = 0.5;
  Scalar domain_volume_ = 0.0;
  Scalar solid_volume_ = 0.0;
  Scalar void_volume_ = 0.0;
  Vector volumes_;
  Vector lower_;
  Vector upper_;
  Vector initial_;
  std::vector<PassiveTag> tags_;
};

}  // namespace sparlab
