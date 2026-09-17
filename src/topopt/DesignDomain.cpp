#include "sparlab/topopt/DesignDomain.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
#include <sstream>

namespace sparlab {

DesignDomain::DesignDomain(const FemModel& model, Scalar volume_fraction,
                           Scalar initial_density,
                           const std::vector<PassiveRegionSpec>& passive)
    : num_elements_(model.mesh().num_elements()), volume_fraction_(volume_fraction) {
  if (!(volume_fraction_ > 0.0 && volume_fraction_ <= 1.0)) {
    std::ostringstream os;
    os << "volume fraction must lie in (0, 1] (got " << volume_fraction_ << ")";
    throw ConfigError(os.str());
  }

  volumes_ = model.element_volumes();
  domain_volume_ = volumes_.sum();
  if (!(domain_volume_ > 0.0)) {
    throw ConfigError("design domain has non-positive volume");
  }

  tags_.assign(static_cast<std::size_t>(num_elements_), PassiveTag::Free);
  lower_ = Vector::Zero(num_elements_);
  upper_ = Vector::Ones(num_elements_);

  const Scalar start = initial_density >= 0.0 ? initial_density : volume_fraction_;
  if (!(start >= 0.0 && start <= 1.0)) {
    std::ostringstream os;
    os << "initial density must lie in [0, 1] (got " << start << ")";
    throw ConfigError(os.str());
  }
  initial_ = Vector::Constant(num_elements_, start);

  for (const PassiveRegionSpec& spec : passive) {
    const std::vector<Index> elements = spec.region.select_elements(model.mesh());
    if (elements.empty()) {
      throw ConfigError("passive region '" + spec.region.name +
                        "' selected no elements; check its geometry against the design "
                        "domain extents");
    }
    if (!spec.solid && !(spec.void_density >= 0.0 && spec.void_density < 1.0)) {
      std::ostringstream os;
      os << "passive void region '" << spec.region.name
         << "' needs a void density in [0, 1) (got " << spec.void_density << ")";
      throw ConfigError(os.str());
    }
    const Scalar value = spec.solid ? 1.0 : spec.void_density;
    const PassiveTag tag = spec.solid ? PassiveTag::Solid : PassiveTag::Void;
    for (Index e : elements) {
      const PassiveTag existing = tags_[static_cast<std::size_t>(e)];
      if (existing != PassiveTag::Free && existing != tag) {
        std::ostringstream os;
        os << "element " << e << " is claimed by both a passive solid and a passive "
           << "void region (last one: '" << spec.region.name
           << "'); make the regions disjoint";
        throw ConfigError(os.str());
      }
      tags_[static_cast<std::size_t>(e)] = tag;
      lower_(e) = value;
      upper_(e) = value;
      initial_(e) = value;
    }
    log::debug("passive ", (spec.solid ? "solid" : "void"), " region '",
               spec.region.name, "': ", elements.size(), " elements at density ", value);
  }

  for (Index e = 0; e < num_elements_; ++e) {
    switch (tags_[static_cast<std::size_t>(e)]) {
      case PassiveTag::Free: ++num_free_; break;
      case PassiveTag::Solid:
        ++num_solid_;
        solid_volume_ += volumes_(e);
        break;
      case PassiveTag::Void:
        ++num_void_;
        void_volume_ += volumes_(e);
        break;
    }
  }

  if (num_free_ == 0) {
    throw ConfigError(
        "every element is passive; there is nothing left to optimise. Shrink the "
        "passive regions");
  }

  const Scalar target = volume_target();
  if (solid_volume_ > target) {
    std::ostringstream os;
    os << "the passive solid regions already occupy " << solid_volume_ << " m^3, more "
       << "than the volume target of " << target << " m^3 (" << volume_fraction_
       << " of the " << domain_volume_
       << " m^3 domain). Raise topology.volume_fraction or shrink the passive solid "
          "regions";
    throw ConfigError(os.str());
  }
  // The free variables must be able to reach the target on their own.
  Scalar free_volume = 0.0;
  for (Index e = 0; e < num_elements_; ++e) {
    if (is_free(e)) free_volume += volumes_(e);
  }
  if (solid_volume_ + free_volume < target) {
    std::ostringstream os;
    os << "the volume target " << target << " m^3 is unreachable: passive solid ("
       << solid_volume_ << " m^3) plus all free elements at density 1 ("
       << free_volume << " m^3) is less than the target. Lower "
       << "topology.volume_fraction or shrink the passive void regions";
    throw ConfigError(os.str());
  }
}

void DesignDomain::clamp(Vector& x) const {
  if (x.size() != num_elements_) {
    throw ConfigError("clamp() received a design vector of the wrong length");
  }
  for (Index e = 0; e < num_elements_; ++e) {
    x(e) = std::clamp(x(e), lower_(e), upper_(e));
  }
}

std::string DesignDomain::describe() const {
  std::ostringstream os;
  os << num_elements_ << " design variables (" << num_free_ << " free, " << num_solid_
     << " passive solid, " << num_void_ << " passive void), domain volume "
     << domain_volume_ << " m^3, volume target " << volume_target() << " m^3 at fraction "
     << volume_fraction_;
  return os.str();
}

}  // namespace sparlab
