/// \file BoundaryConditions.hpp
/// \brief Displacement constraint and load specifications, and their
///        application to a DofManager / global force vector.
///
/// Sign convention: forces and displacements are positive along the global
/// +x / +y axes. Tractions are given as a stress vector \f$\bar{t}\f$ [Pa] in
/// global components and are integrated over the selected boundary edges,
/// producing consistent nodal forces (not lumped ones).
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/elements/Element.hpp"
#include "sparlab/fem/DofManager.hpp"
#include "sparlab/fem/Selector.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Prescribed displacement on the x and/or y component of a node region.
struct DisplacementConstraint {
  SelectorGroup region;
  bool fix_x = false;
  bool fix_y = false;
  Scalar value_x = 0.0;  ///< prescribed u_x [m]
  Scalar value_y = 0.0;  ///< prescribed u_y [m]
};

/// Concentrated nodal force applied to a node region.
struct PointLoadSpec {
  SelectorGroup region;
  Vector2 force = Vector2::Zero();  ///< [N]
  /// When true, `force` is the *resultant* and is divided equally among the
  /// selected nodes (mesh-independent total load). When false, `force` is
  /// applied to each selected node.
  bool distribute_total = true;
};

/// Constant traction applied to the mesh boundary edges whose two end nodes
/// both lie inside the region.
struct TractionLoadSpec {
  SelectorGroup region;
  Vector2 traction = Vector2::Zero();  ///< [Pa]
};

/// One load case. Multiple cases are combined by the topology optimizer with
/// the given weights; the FE solver solves each independently.
struct LoadCaseSpec {
  std::string name = "load_case";
  Scalar weight = 1.0;  ///< weight in the multi-load compliance objective [-]
  std::vector<PointLoadSpec> point_loads;
  std::vector<TractionLoadSpec> tractions;
  /// Set when the case carries no applied force and is driven purely by the
  /// prescribed displacements of the model (enforced-deflection studies, the
  /// constant-strain patch test). Declaring it explicitly keeps an accidentally
  /// empty load case an error rather than a silently trivial solution.
  bool prescribed_displacement_only = false;
};

/// Apply every constraint to `dofs`.
/// \return the number of distinct DOFs constrained.
/// \throws ConfigError when a region selects no nodes (almost always an
///         authoring mistake, and silently ignoring it produces a singular
///         system later).
Index apply_constraints(const Mesh& mesh,
                        const std::vector<DisplacementConstraint>& constraints,
                        DofManager& dofs);

/// Assemble the global force vector [N] of one load case.
/// \throws ConfigError for empty regions or when a traction region matches no
///         boundary edge.
Vector assemble_load_vector(const Mesh& mesh, const Element& element,
                            const LoadCaseSpec& load_case, Scalar thickness,
                            const IntegrationOptions& integration);

}  // namespace sparlab
