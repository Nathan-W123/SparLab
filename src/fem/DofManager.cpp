#include "sparlab/fem/DofManager.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <algorithm>
#include <sstream>

namespace sparlab {

DofManager::DofManager(Index num_nodes) : num_nodes_(num_nodes) {
  if (num_nodes_ < 0) throw ModelError("DofManager requires a non-negative node count");
  constrained_.assign(static_cast<std::size_t>(num_dofs()), 0);
  values_.assign(static_cast<std::size_t>(num_dofs()), 0.0);
}

Index DofManager::dof(Index node, int component) const {
  if (node < 0 || node >= num_nodes_) {
    std::ostringstream os;
    os << "node index " << node << " is outside [0, " << num_nodes_ - 1 << "]";
    throw ModelError(os.str());
  }
  if (component < 0 || component >= kDofsPerNode) {
    std::ostringstream os;
    os << "DOF component " << component << " is outside [0, " << kDofsPerNode - 1
       << "] (0 = x, 1 = y)";
    throw ModelError(os.str());
  }
  return node * kDofsPerNode + component;
}

void DofManager::prescribe(Index node, int component, Scalar value) {
  prescribe_dof(dof(node, component), value);
}

void DofManager::prescribe_dof(Index d, Scalar value) {
  if (d < 0 || d >= num_dofs()) {
    std::ostringstream os;
    os << "DOF index " << d << " is outside [0, " << num_dofs() - 1 << "]";
    throw ModelError(os.str());
  }
  constrained_[static_cast<std::size_t>(d)] = 1;
  values_[static_cast<std::size_t>(d)] = value;
  dirty_ = true;
}

void DofManager::rebuild() const {
  free_.clear();
  fixed_.clear();
  reduced_.assign(static_cast<std::size_t>(num_dofs()), -1);
  for (Index d = 0; d < num_dofs(); ++d) {
    if (constrained_[static_cast<std::size_t>(d)]) {
      fixed_.push_back(d);
    } else {
      reduced_[static_cast<std::size_t>(d)] = static_cast<Index>(free_.size());
      free_.push_back(d);
    }
  }
  dirty_ = false;
}

const std::vector<Index>& DofManager::free_dofs() const {
  if (dirty_) rebuild();
  return free_;
}

const std::vector<Index>& DofManager::constrained_dofs() const {
  if (dirty_) rebuild();
  return fixed_;
}

Index DofManager::reduced_index(Index dof) const {
  if (dirty_) rebuild();
  if (dof < 0 || dof >= num_dofs()) {
    std::ostringstream os;
    os << "DOF index " << dof << " is outside [0, " << num_dofs() - 1 << "]";
    throw ModelError(os.str());
  }
  return reduced_[static_cast<std::size_t>(dof)];
}

bool DofManager::has_nonzero_prescribed() const {
  for (Index d = 0; d < num_dofs(); ++d) {
    if (constrained_[static_cast<std::size_t>(d)] &&
        values_[static_cast<std::size_t>(d)] != 0.0) {
      return true;
    }
  }
  return false;
}

Vector DofManager::prescribed_vector() const {
  Vector up = Vector::Zero(num_dofs());
  for (Index d = 0; d < num_dofs(); ++d) {
    if (constrained_[static_cast<std::size_t>(d)]) {
      up(d) = values_[static_cast<std::size_t>(d)];
    }
  }
  return up;
}

Vector DofManager::expand(const Vector& reduced) const {
  const auto& fdofs = free_dofs();
  if (reduced.size() != static_cast<Eigen::Index>(fdofs.size())) {
    std::ostringstream os;
    os << "cannot expand a reduced vector of length " << reduced.size() << "; the model "
       << "has " << fdofs.size() << " free DOFs";
    throw ModelError(os.str());
  }
  Vector full = prescribed_vector();
  for (std::size_t k = 0; k < fdofs.size(); ++k) {
    full(fdofs[k]) = reduced(static_cast<Eigen::Index>(k));
  }
  return full;
}

Vector DofManager::restrict_to_free(const Vector& full) const {
  const auto& fdofs = free_dofs();
  if (full.size() != num_dofs()) {
    std::ostringstream os;
    os << "cannot restrict a vector of length " << full.size() << "; the model has "
       << num_dofs() << " DOFs";
    throw ModelError(os.str());
  }
  Vector reduced(static_cast<Eigen::Index>(fdofs.size()));
  for (std::size_t k = 0; k < fdofs.size(); ++k) {
    reduced(static_cast<Eigen::Index>(k)) = full(fdofs[k]);
  }
  return reduced;
}

void DofManager::element_dofs(const Index* nodes, int nodes_per_elem, Index* out) const {
  for (int a = 0; a < nodes_per_elem; ++a) {
    out[kDofsPerNode * a + 0] = nodes[a] * kDofsPerNode + 0;
    out[kDofsPerNode * a + 1] = nodes[a] * kDofsPerNode + 1;
  }
}

}  // namespace sparlab
