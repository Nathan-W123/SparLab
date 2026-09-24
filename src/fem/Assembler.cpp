#include "sparlab/fem/Assembler.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
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

void Assembler::build_pattern() const {
  const Mesh& mesh = model_.mesh();
  const Index nn = mesh.num_nodes();
  const Index ne = mesh.num_elements();
  const int npe = mesh.nodes_per_elem();
  const int dim = mesh.dim();

  // Elements around each node.
  std::vector<Index> inc_ptr(static_cast<std::size_t>(nn) + 1, 0);
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 0; a < npe; ++a) ++inc_ptr[static_cast<std::size_t>(nodes[a]) + 1];
  }
  for (Index n = 0; n < nn; ++n) {
    inc_ptr[static_cast<std::size_t>(n) + 1] += inc_ptr[static_cast<std::size_t>(n)];
  }
  std::vector<Index> inc(static_cast<std::size_t>(inc_ptr.back()));
  {
    std::vector<Index> fill(inc_ptr.begin(), inc_ptr.end() - 1);
    for (Index e = 0; e < ne; ++e) {
      const Index* nodes = mesh.element_nodes(e);
      for (int a = 0; a < npe; ++a) inc[static_cast<std::size_t>(fill[static_cast<std::size_t>(nodes[a])]++)] = e;
    }
  }
  // Sorted node neighbourhoods (the node itself included).
  std::vector<Index> nbr_ptr(static_cast<std::size_t>(nn) + 1, 0);
  std::vector<Index> nbr;
  nbr.reserve(static_cast<std::size_t>(nn) * (dim == 3 ? 27 : 9));
  {
    std::vector<Index> marker(static_cast<std::size_t>(nn), -1);
    std::vector<Index> list;
    for (Index b = 0; b < nn; ++b) {
      list.clear();
      for (Index k = inc_ptr[static_cast<std::size_t>(b)]; k < inc_ptr[static_cast<std::size_t>(b) + 1]; ++k) {
        const Index* nodes = mesh.element_nodes(inc[static_cast<std::size_t>(k)]);
        for (int a = 0; a < npe; ++a) {
          if (marker[static_cast<std::size_t>(nodes[a])] != b) {
            marker[static_cast<std::size_t>(nodes[a])] = b;
            list.push_back(nodes[a]);
          }
        }
      }
      std::sort(list.begin(), list.end());
      nbr.insert(nbr.end(), list.begin(), list.end());
      nbr_ptr[static_cast<std::size_t>(b) + 1] = static_cast<Index>(nbr.size());
    }
  }
  // Column-major structure: every column of node b lists the rows of its
  // neighbours' blocks in ascending order.
  Pattern& p = pattern_;
  const Index n = nn * dim;
  p.size = n;
  p.outer.assign(static_cast<std::size_t>(n) + 1, 0);
  for (Index b = 0; b < nn; ++b) {
    const Index deg = nbr_ptr[static_cast<std::size_t>(b) + 1] - nbr_ptr[static_cast<std::size_t>(b)];
    for (int kb = 0; kb < dim; ++kb) {
      const Index c = b * dim + kb;
      p.outer[static_cast<std::size_t>(c) + 1] = p.outer[static_cast<std::size_t>(c)] + deg * dim;
    }
  }
  p.inner.resize(static_cast<std::size_t>(p.outer.back()));
  for (Index b = 0; b < nn; ++b) {
    for (int kb = 0; kb < dim; ++kb) {
      Index q = p.outer[static_cast<std::size_t>(b * dim + kb)];
      for (Index k = nbr_ptr[static_cast<std::size_t>(b)]; k < nbr_ptr[static_cast<std::size_t>(b) + 1]; ++k) {
        for (int ka = 0; ka < dim; ++ka) p.inner[static_cast<std::size_t>(q++)] = nbr[static_cast<std::size_t>(k)] * dim + ka;
      }
    }
  }
  p.block_offset.resize(static_cast<std::size_t>(ne) * npe * npe);
  for (Index e = 0; e < ne; ++e) {
    const Index* nodes = mesh.element_nodes(e);
    for (int b = 0; b < npe; ++b) {
      const Index* begin = nbr.data() + nbr_ptr[static_cast<std::size_t>(nodes[b])];
      const Index* end = nbr.data() + nbr_ptr[static_cast<std::size_t>(nodes[b]) + 1];
      for (int a = 0; a < npe; ++a) {
        const Index* hit = std::lower_bound(begin, end, nodes[a]);
        p.block_offset[(static_cast<std::size_t>(e) * npe + a) * npe + b] =
            static_cast<Index>(hit - begin) * dim;
      }
    }
  }
  p.ready = true;
  log::debug("assembler: cached the sparsity pattern (", n, " DOFs, ", p.inner.size(),
             " stored entries)");
}

bool Assembler::matches_pattern(const SparseMatrix& full) const {
  if (!pattern_.ready || !full.isCompressed()) return false;
  if (full.rows() != pattern_.size || full.cols() != pattern_.size) return false;
  if (static_cast<std::size_t>(full.nonZeros()) != pattern_.inner.size()) return false;
  return std::equal(pattern_.outer.begin(), pattern_.outer.end(), full.outerIndexPtr()) &&
         std::equal(pattern_.inner.begin(), pattern_.inner.end(), full.innerIndexPtr());
}

template <typename ElementMatrix>
SparseMatrix Assembler::scatter(const ElementMatrix& element_matrix, const Vector* scale) const {
  if (!pattern_.ready) build_pattern();
  const Mesh& mesh = model_.mesh();
  const Index ne = mesh.num_elements();
  const int npe = mesh.nodes_per_elem();
  const int dim = mesh.dim();
  const Pattern& p = pattern_;
  SparseMatrix k(p.size, p.size);
  k.resizeNonZeros(static_cast<Eigen::Index>(p.inner.size()));
  std::copy(p.outer.begin(), p.outer.end(), k.outerIndexPtr());
  std::copy(p.inner.begin(), p.inner.end(), k.innerIndexPtr());
  Scalar* values = k.valuePtr();
  std::fill(values, values + p.inner.size(), 0.0);
  const StorageIndex* outer = p.outer.data();
  for (Index e = 0; e < ne; ++e) {
    const Scalar s = scale ? (*scale)(e) : 1.0;
    const Matrix& ke = element_matrix(e);
    const Index* nodes = mesh.element_nodes(e);
    const Index* offsets = p.block_offset.data() + static_cast<std::size_t>(e) * npe * npe;
    for (int b = 0; b < npe; ++b) {
      for (int kb = 0; kb < dim; ++kb) {
        const Index column = dim * b + kb;
        const Index base = outer[nodes[b] * dim + kb];
        for (int a = 0; a < npe; ++a) {
          Scalar* target = values + base + offsets[a * npe + b];
          for (int ka = 0; ka < dim; ++ka) target[ka] += s * ke(dim * a + ka, column);
        }
      }
    }
  }
  return k;
}

SparseMatrix Assembler::from_reduction(const Reduction& map, const SparseMatrix& full) {
  SparseMatrix out(map.rows, map.cols);
  out.resizeNonZeros(static_cast<Eigen::Index>(map.inner.size()));
  std::copy(map.outer.begin(), map.outer.end(), out.outerIndexPtr());
  std::copy(map.inner.begin(), map.inner.end(), out.innerIndexPtr());
  const Scalar* source = full.valuePtr();
  Scalar* values = out.valuePtr();
  for (std::size_t k = 0; k < map.source.size(); ++k) values[k] = source[map.source[k]];
  return out;
}

SparseMatrix Assembler::assemble_stiffness(const Vector* scale) const {
  const Mesh& mesh = model_.mesh();
  const Index ne = mesh.num_elements();
  check_scale(scale, ne, "stiffness");
  if (use_pattern_ && (scale == nullptr || scale->minCoeff() > 0.0)) {
    return scatter([this](Index e) -> const Matrix& { return element_stiffness(e); }, scale);
  }

  const int npe = mesh.nodes_per_elem();
  const int edofs = npe * mesh.dim();
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

  if (use_pattern_ && type == MassType::Consistent &&
      (scale == nullptr || scale->minCoeff() > 0.0)) {
    return scatter([this](Index e) -> const Matrix& { return element_mass(e); }, scale);
  }

  const int npe = mesh.nodes_per_elem();
  const int edofs = npe * mesh.dim();
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
      // Row-sum (Hinton-Rock-Zienkiewicz style) lumping. For the Q4 and Hex8
      // with a consistent mass matrix the row sums reproduce the element mass
      // exactly, so total mass is conserved.
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
  // The thickness factor is 1 for a solid mesh, which leaves the product exact.
  const Scalar thickness = mesh.dim() == 2 ? model_.thickness() : 1.0;
  Scalar mass = 0.0;
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Scalar s = scale ? (*scale)(e) : 1.0;
    mass += s * model_.material().density() * mesh.element_measure(e) * thickness;
  }
  return mass;
}

SparseMatrix Assembler::reduce_free_free(const SparseMatrix& full) const {
  const DofManager& dofs = model_.dofs();
  const auto& free = dofs.free_dofs();
  const Index nf = static_cast<Index>(free.size());

  if (use_pattern_ && matches_pattern(full)) {
    Reduction& map = free_free_;
    if (map.key != free || map.rows != nf) {
      map.key = free;
      map.rows = nf;
      map.cols = nf;
      map.outer.assign(static_cast<std::size_t>(nf) + 1, 0);
      map.inner.clear();
      map.source.clear();
      const StorageIndex* outer = full.outerIndexPtr();
      const StorageIndex* inner = full.innerIndexPtr();
      for (Index rc = 0; rc < nf; ++rc) {
        const Index col = free[static_cast<std::size_t>(rc)];
        for (StorageIndex k = outer[col]; k < outer[col + 1]; ++k) {
          const Index rr = dofs.reduced_index(inner[k]);
          if (rr < 0) continue;
          map.inner.push_back(rr);
          map.source.push_back(k);
        }
        map.outer[static_cast<std::size_t>(rc) + 1] = static_cast<StorageIndex>(map.inner.size());
      }
    }
    return from_reduction(map, full);
  }

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

  if (use_pattern_ && matches_pattern(full)) {
    Reduction& map = free_prescribed_;
    if (map.key != fixed || map.cols != static_cast<Index>(fixed.size())) {
      map.key = fixed;
      map.rows = dofs.num_free();
      map.cols = static_cast<Index>(fixed.size());
      map.outer.assign(fixed.size() + 1, 0);
      map.inner.clear();
      map.source.clear();
      const StorageIndex* outer = full.outerIndexPtr();
      const StorageIndex* inner = full.innerIndexPtr();
      for (std::size_t fc = 0; fc < fixed.size(); ++fc) {
        const Index col = fixed[fc];
        for (StorageIndex k = outer[col]; k < outer[col + 1]; ++k) {
          const Index rr = dofs.reduced_index(inner[k]);
          if (rr < 0) continue;
          map.inner.push_back(rr);
          map.source.push_back(k);
        }
        map.outer[fc + 1] = static_cast<StorageIndex>(map.inner.size());
      }
    }
    return from_reduction(map, full);
  }

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
