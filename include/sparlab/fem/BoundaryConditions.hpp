/// \file BoundaryConditions.hpp
/// \brief Displacement constraint and load specifications, and their
///        application to a DofManager / global force vector.
///
/// Sign convention: forces and displacements are positive along the global
/// +x / +y / +z axes. Tractions are given as a stress vector \f$\bar{t}\f$ [Pa]
/// in global components and are integrated over the selected boundary faces
/// (edges in 2-D), producing consistent nodal forces (not lumped ones).
///
/// Vectors are three-component throughout; on a 2-D model the z entry of a
/// force or traction must be zero and `fix_z` is an error, so a deck cannot
/// silently lose a component it thought it had applied.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/elements/Element.hpp"
#include "sparlab/fem/DofManager.hpp"
#include "sparlab/fem/Selector.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Prescribed displacement on the x, y and/or z component of a node region.
struct DisplacementConstraint {
  SelectorGroup region;
  bool fix_x = false;
  bool fix_y = false;
  bool fix_z = false;
  Scalar value_x = 0.0;  ///< prescribed u_x [m]
  Scalar value_y = 0.0;  ///< prescribed u_y [m]
  Scalar value_z = 0.0;  ///< prescribed u_z [m]

  bool fixes(int component) const {
    return component == 0 ? fix_x : component == 1 ? fix_y : fix_z;
  }
  Scalar value(int component) const {
    return component == 0 ? value_x : component == 1 ? value_y : value_z;
  }
};

/// Concentrated nodal force applied to a node region.
struct PointLoadSpec {
  SelectorGroup region;
  Vector3 force = Vector3::Zero();  ///< [N]; z must be 0 on a 2-D model
  /// When true, `force` is the *resultant* and is divided equally among the
  /// selected nodes (mesh-independent total load). When false, `force` is
  /// applied to each selected node.
  bool distribute_total = true;
};

/// Constant traction applied to the mesh boundary faces (edges in 2-D) whose
/// nodes all lie inside the region.
struct TractionLoadSpec {
  SelectorGroup region;
  Vector3 traction = Vector3::Zero();  ///< [Pa]; z must be 0 on a 2-D model
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
///         system later), or when a 2-D model is asked to fix z.
Index apply_constraints(const Mesh& mesh,
                        const std::vector<DisplacementConstraint>& constraints,
                        DofManager& dofs);

/// Assemble the global force vector [N] of one load case.
/// \throws ConfigError for empty regions, a traction region that matches no
///         boundary face, or an out-of-plane component on a 2-D model.
Vector assemble_load_vector(const Mesh& mesh, const Element& element,
                            const LoadCaseSpec& load_case, Scalar thickness,
                            const IntegrationOptions& integration);

}  // namespace sparlab
