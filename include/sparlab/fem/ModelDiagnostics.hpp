/// \file ModelDiagnostics.hpp
/// \brief Pre-solve detection of ill-posed models.
///
/// Two independent checks are performed, both before any factorisation:
///
/// 1. **Per-component rigid-body check.** The element graph is split into
///    connected components (elements sharing a node). For a component with
///    node set \f$\mathcal{N}\f$ the space of infinitesimal rigid-body motions
///    is spanned in 2-D by
///    \f$ r_1 = (1,0),\ r_2 = (0,1),\ r_3 = (-(y-y_c),\ x-x_c) \f$
///    and in 3-D by the three translations plus the three rotations
///    \f$ e_j \times (x - x_c) \f$. Because the Q4 and Hex8 spaces contain all
///    linear fields exactly, each \f$r_i\f$ produces exactly zero strain, hence
///    zero energy. A component is properly constrained only if no non-trivial
///    combination \f$\sum_i c_i r_i\f$ vanishes on every prescribed DOF of that
///    component; equivalently the matrix of rigid modes evaluated at the
///    component's prescribed DOFs must have full rank (3 in 2-D, 6 in 3-D).
///
/// 2. **Floating-component check.** A component with no prescribed DOF at all
///    is reported explicitly, since that is the most common authoring mistake.
///
/// Both checks are exact for rigid-body singularity. Internal mechanisms
/// (for example two blocks joined at a single node) are *not* detected
/// analytically; they surface as a non-positive Cholesky pivot, which
/// LinearSolver reports as a solver error with the same hint text.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// One connected group of elements and its constraint status.
struct MeshComponent {
  std::vector<Index> elements;
  std::vector<Index> nodes;
  Index prescribed_dofs = 0;
  /// Dimension of the rigid-body null space left by the constraints (0 = ok,
  /// 3 in 2-D / 6 in 3-D = completely free).
  int rigid_null_dimension = 0;
};

struct ModelDiagnostics {
  Index num_dofs = 0;
  Index num_free_dofs = 0;
  Index num_prescribed_dofs = 0;
  std::vector<MeshComponent> components;
  /// Human-readable, actionable descriptions; empty when the model is well
  /// posed for a static solve.
  std::vector<std::string> problems;

  bool well_posed() const { return problems.empty(); }
};

/// Number of rigid-body modes of an unconstrained body: 3 in 2-D, 6 in 3-D.
int rigid_body_mode_count(int dim);

/// Run the checks above. Never throws for an ill-posed model - the caller
/// decides what to do.
ModelDiagnostics diagnose_model(const FemModel& model);

/// Run the checks and raise `ModelError` listing every problem found.
void require_well_posed(const FemModel& model);

/// Connected components of the element graph (elements sharing >= 1 node).
std::vector<std::vector<Index>> element_connected_components(const Mesh& mesh);

}  // namespace sparlab
