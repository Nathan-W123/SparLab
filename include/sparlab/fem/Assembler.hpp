/// \file Assembler.hpp
/// \brief Sparse global assembly of stiffness and mass matrices.
///
/// Assembly builds a triplet list and compresses it into a CSC matrix:
/// \f[
///   \mathbf{K} = \sum_e s_e\, \mathbf{A}_e^T \mathbf{K}_e^0 \mathbf{A}_e,
///   \qquad
///   \mathbf{M} = \sum_e m_e\, \mathbf{A}_e^T \mathbf{M}_e^0 \mathbf{A}_e,
/// \f]
/// where \f$\mathbf{A}_e\f$ is the (implicit) DOF gather operator, and
/// \f$s_e\f$, \f$m_e\f$ are per-element scale factors. For a plain FE analysis
/// both scale factors are 1; the topology optimizer supplies the SIMP stiffness
/// factor \f$E(\rho_e)/E_0\f$ and a mass interpolation factor.
///
/// On a *uniform structured* mesh every cell has identical geometry, so
/// \f$\mathbf{K}_e^0\f$ and \f$\mathbf{M}_e^0\f$ are computed once and reused.
/// This is the dominant cost saving in a topology optimisation loop. The
/// generic per-element path is always available and the two are checked against
/// each other in the test suite.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"

#include <vector>

namespace sparlab {

class Assembler {
 public:
  explicit Assembler(const FemModel& model);

  /// Unit-density element stiffness matrix of element `e` [N/m].
  const Matrix& element_stiffness(Index e) const;

  /// Unit-density consistent element mass matrix of element `e` [kg].
  const Matrix& element_mass(Index e) const;

  /// True when a single cached element matrix is reused for the whole mesh.
  bool uses_element_cache() const { return uniform_; }

  /// Assemble the global stiffness matrix [N/m].
  /// \param scale optional per-element stiffness factors (length
  ///        num_elements); `nullptr` means all ones.
  SparseMatrix assemble_stiffness(const Vector* scale = nullptr) const;

  /// Assemble the global mass matrix [kg].
  /// \param type consistent or row-sum lumped.
  /// \param scale optional per-element mass factors; `nullptr` means all ones.
  SparseMatrix assemble_mass(MassType type, const Vector* scale = nullptr) const;

  /// Extract the free-free block \f$K_{ff}\f$ of a full-size matrix.
  SparseMatrix reduce_free_free(const SparseMatrix& full) const;

  /// Extract the free-prescribed block \f$K_{fp}\f$ of a full-size matrix.
  SparseMatrix reduce_free_prescribed(const SparseMatrix& full) const;

  /// Extract the prescribed-all block \f$[K_{pf}\ K_{pp}]\f$ applied to a full
  /// displacement vector, i.e. rows of `full` at prescribed DOFs.
  SparseMatrix reduce_prescribed_rows(const SparseMatrix& full) const;

  /// Total structural mass of the model [kg] for the given element factors.
  Scalar total_mass(const Vector* scale = nullptr) const;

 private:
  void build_cache();
  Matrix compute_element_stiffness(Index e) const;
  Matrix compute_element_mass(Index e) const;

  const FemModel& model_;
  bool uniform_ = false;
  Matrix cached_k_;
  Matrix cached_m_;
  mutable std::vector<Matrix> per_element_k_;
  mutable std::vector<Matrix> per_element_m_;
};

}  // namespace sparlab
