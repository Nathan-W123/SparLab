#include "sparlab/fem/Assembler.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

/// A uniform structured grid has geometrically identical cells, so the element
/// matrices need to be integrated only once.
bool mesh_is_uniform(const Mesh& mesh) {
  const auto& info = mesh.structured_info();
  return info.has_value() && info->uniform;
}

void check_scale(const Vector* scale, Index num_elements, const char* what) {
  if (scale == nullptr) return;
  if (scale->size() != num_elements) {
    std::ostringstream os;
    os << what << " scale vector has length " << scale->size() << " but the mesh has "
       << num_elements << " elements";
    throw ModelError(os.str());
  }
  if (!scale->allFinite()) {
    throw ModelError(std::string(what) + " scale vector contains NaN or Inf");
  }
  if (scale->minCoeff() < 0.0) {
    std::ostringstream os;
    os << what << " scale vector has a negative entry (min = " << scale->minCoeff()
       << "); element matrices would lose positive semi-definiteness";
    throw ModelError(os.str());
  }
}

}  // namespace

Assembler::Assembler(const FemModel& model) : model_(model) {
  uniform_ = mesh_is_uniform(model_.mesh());
  build_cache();
}

void Assembler::build_cache() {
  if (uniform_) {
    cached_k_ = compute_element_stiffness(0);
    cached_m_ = compute_element_mass(0);
    log::debug("assembler: uniform structured mesh detected, reusing one ",
               cached_k_.rows(), "x", cached_k_.cols(), " element matrix pair");
  } else {
    per_element_k_.assign(static_cast<std::size_t>(model_.mesh().num_elements()), Matrix());
    per_element_m_.assign(static_cast<std::size_t>(model_.mesh().num_elements()), Matrix());
  }
}

Matrix Assembler::compute_element_stiffness(Index e) const {
  return model_.element().stiffness(model_.mesh().element_coordinates(e),
                                    model_.constitutive(), model_.thickness(),
                                    model_.integration());
}

Matrix Assembler::compute_element_mass(Index e) const {
  return model_.element().consistent_mass(model_.mesh().element_coordinates(e),
                                          model_.material().density(),
                                          model_.thickness(), model_.integration());
}

const Matrix& Assembler::element_stiffness(Index e) const {
  if (uniform_) return cached_k_;
  Matrix& slot = per_element_k_[static_cast<std::size_t>(e)];
  if (slot.size() == 0) slot = compute_element_stiffness(e);
  return slot;
}

const Matrix& Assembler::element_mass(Index e) const {
  if (uniform_) return cached_m_;
  Matrix& slot = per_element_m_[static_cast<std::size_t>(e)];
  if (slot.size() == 0) slot = compute_element_mass(e);
  return slot;
}

SparseMatrix Assembler::assemble_stiffness(const Vector* scale) const {
  const Mesh& mesh = model_.mesh();
  const Index ne = mesh.num_elements();
  check_scale(scale, ne, "stiffness");

  const int npe = mesh.nodes_per_elem();
  const int edofs = npe * kDofsPerNode;
  TripletList triplets;
  triplets.reserve(static_cast<std::size_t>(ne) * edofs * edofs);

  std::vector<Index> gdofs(static_cast<std::size_t>(edofs));
  for (Index e = 0; e < ne; ++e) {
    const Scalar s = scale ? (*scale)(e) : 1.0;
    if (s == 0.0) continue;
    const Matrix& ke = element_stiffness(e);
    model_.dofs().element_dofs(mesh.element_nodes(e), npe, gdofs.data());
    for (int i = 0; i < edofs; ++i) {
      for (int j = 0; j < edofs; ++j) {
        triplets.emplace_back(gdofs[static_cast<std::size_t>(i)],
                              gdofs[static_cast<std::size_t>(j)], s * ke(i, j));
      }
    }
  }

  SparseMatrix k(model_.dofs().num_dofs(), model_.dofs().num_dofs());
  k.setFromTriplets(triplets.begin(), triplets.end());
  k.makeCompressed();
  return k;
}

SparseMatrix Assembler::assemble_mass(MassType type, const Vector* scale) const {
  const Mesh& mesh = model_.mesh();
  const Index ne = mesh.num_elements();
  check_scale(scale, ne, "mass");
  if (model_.material().density() <= 0.0) {
    throw ModelError(
        "mass assembly requires a positive material density; set 'material.density' "
        "to run modal analysis");
  }

  const int npe = mesh.nodes_per_elem();
  const int edofs = npe * kDofsPerNode;
  TripletList triplets;
  triplets.reserve(static_cast<std::size_t>(ne) * edofs * edofs);

  std::vector<Index> gdofs(static_cast<std::size_t>(edofs));
  for (Index e = 0; e < ne; ++e) {
    const Scalar s = scale ? (*scale)(e) : 1.0;
    if (s == 0.0) continue;
    const Matrix& me = element_mass(e);
    model_.dofs().element_dofs(mesh.element_nodes(e), npe, gdofs.data());
    if (type == MassType::Consistent) {
      for (int i = 0; i < edofs; ++i) {
        for (int j = 0; j < edofs; ++j) {
          triplets.emplace_back(gdofs[static_cast<std::size_t>(i)],
                                gdofs[static_cast<std::size_t>(j)], s * me(i, j));
        }
      }
    } else {
      // Row-sum (Hinton-Rock-Zienkiewicz style) lumping. For the Q4 with a
      // consistent mass matrix the row sums reproduce the element mass exactly,
      // so total mass is conserved.
      for (int i = 0; i < edofs; ++i) {
        const Scalar lumped = me.row(i).sum();
        triplets.emplace_back(gdofs[static_cast<std::size_t>(i)],
                              gdofs[static_cast<std::size_t>(i)], s * lumped);
      }
    }
  }

  SparseMatrix m(model_.dofs().num_dofs(), model_.dofs().num_dofs());
  m.setFromTriplets(triplets.begin(), triplets.end());
  m.makeCompressed();
  return m;
}

Scalar Assembler::total_mass(const Vector* scale) const {
  const Mesh& mesh = model_.mesh();
  check_scale(scale, mesh.num_elements(), "mass");
  Scalar mass = 0.0;
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Scalar s = scale ? (*scale)(e) : 1.0;
    mass += s * model_.material().density() * mesh.element_area(e) * model_.thickness();
  }
  return mass;
}

SparseMatrix Assembler::reduce_free_free(const SparseMatrix& full) const {
  const DofManager& dofs = model_.dofs();
  const auto& free = dofs.free_dofs();
  const Index nf = static_cast<Index>(free.size());

  TripletList triplets;
  triplets.reserve(static_cast<std::size_t>(full.nonZeros()));
  for (Eigen::Index col = 0; col < full.outerSize(); ++col) {
    const Index rc = dofs.reduced_index(static_cast<Index>(col));
    if (rc < 0) continue;
    for (SparseMatrix::InnerIterator it(full, col); it; ++it) {
      const Index rr = dofs.reduced_index(static_cast<Index>(it.row()));
      if (rr < 0) continue;
      triplets.emplace_back(rr, rc, it.value());
    }
  }
  SparseMatrix reduced(nf, nf);
  reduced.setFromTriplets(triplets.begin(), triplets.end());
  reduced.makeCompressed();
  return reduced;
}

SparseMatrix Assembler::reduce_free_prescribed(const SparseMatrix& full) const {
  const DofManager& dofs = model_.dofs();
  const auto& fixed = dofs.constrained_dofs();
  std::vector<Index> fixed_index(static_cast<std::size_t>(dofs.num_dofs()), -1);
  for (std::size_t k = 0; k < fixed.size(); ++k) {
    fixed_index[static_cast<std::size_t>(fixed[k])] = static_cast<Index>(k);
  }

  TripletList triplets;
  for (Eigen::Index col = 0; col < full.outerSize(); ++col) {
    const Index fc = fixed_index[static_cast<std::size_t>(col)];
    if (fc < 0) continue;
    for (SparseMatrix::InnerIterator it(full, col); it; ++it) {
      const Index rr = dofs.reduced_index(static_cast<Index>(it.row()));
      if (rr < 0) continue;
      triplets.emplace_back(rr, fc, it.value());
    }
  }
  SparseMatrix block(dofs.num_free(), static_cast<Index>(fixed.size()));
  block.setFromTriplets(triplets.begin(), triplets.end());
  block.makeCompressed();
  return block;
}

SparseMatrix Assembler::reduce_prescribed_rows(const SparseMatrix& full) const {
  const DofManager& dofs = model_.dofs();
  const auto& fixed = dofs.constrained_dofs();
  std::vector<Index> fixed_index(static_cast<std::size_t>(dofs.num_dofs()), -1);
  for (std::size_t k = 0; k < fixed.size(); ++k) {
    fixed_index[static_cast<std::size_t>(fixed[k])] = static_cast<Index>(k);
  }

  TripletList triplets;
  for (Eigen::Index col = 0; col < full.outerSize(); ++col) {
    for (SparseMatrix::InnerIterator it(full, col); it; ++it) {
      const Index fr = fixed_index[static_cast<std::size_t>(it.row())];
      if (fr < 0) continue;
      triplets.emplace_back(fr, static_cast<Index>(col), it.value());
    }
  }
  SparseMatrix block(static_cast<Index>(fixed.size()), dofs.num_dofs());
  block.setFromTriplets(triplets.begin(), triplets.end());
  block.makeCompressed();
  return block;
}

}  // namespace sparlab
