/// \file Selector.hpp
/// \brief Geometric region selectors for boundary conditions, loads and passive
///        topology regions.
///
/// Selectors are purely geometric so that the *same* configuration block can be
/// applied to meshes of different resolution, and so that a topology extracted
/// from a density field can be re-constrained without re-authoring indices.
///
/// A `SelectorGroup` is the union of its members, optionally complemented.
/// Tolerances are relative to the mesh diagonal so that degenerate boxes
/// (`xmin == xmax`) reliably capture a line of nodes.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <limits>
#include <string>
#include <vector>

namespace sparlab {

enum class SelectorKind {
  All,          ///< every entity
  Box,          ///< axis-aligned box, inclusive within tolerance
  Circle,       ///< closed disc of radius `radius` about `center`
  Annulus,      ///< closed ring with `inner_radius` <= r <= `radius`
  NodeIds,      ///< explicit node indices (node selection only)
  ElementIds,   ///< explicit element indices (element selection only)
  NearestNode   ///< the single node closest to `point`
};

/// One geometric primitive.
struct Selector {
  SelectorKind kind = SelectorKind::All;

  // Box bounds [m]; unspecified sides default to +-infinity.
  Scalar xmin = -std::numeric_limits<Scalar>::infinity();
  Scalar xmax = std::numeric_limits<Scalar>::infinity();
  Scalar ymin = -std::numeric_limits<Scalar>::infinity();
  Scalar ymax = std::numeric_limits<Scalar>::infinity();

  Vector2 center = Vector2::Zero();  ///< circle / annulus centre [m]
  Scalar radius = 0.0;               ///< outer radius [m]
  Scalar inner_radius = 0.0;         ///< annulus inner radius [m]

  Vector2 point = Vector2::Zero();   ///< NearestNode target [m]
  std::vector<Index> ids;            ///< explicit indices

  /// Absolute geometric tolerance [m]. When <= 0 a tolerance of
  /// `1e-9 * mesh diagonal` is used.
  Scalar tolerance = 0.0;

  bool contains(const Vector2& x, Scalar tol) const;
};

/// Union of selectors with an optional complement, plus a label for diagnostics.
struct SelectorGroup {
  std::string name = "region";
  std::vector<Selector> members;
  bool invert = false;  ///< select the complement of the union

  /// Node indices matched by this group (ascending, unique).
  std::vector<Index> select_nodes(const Mesh& mesh) const;

  /// Element indices matched by this group, tested at element centroids
  /// (ascending, unique).
  std::vector<Index> select_elements(const Mesh& mesh) const;

  bool empty() const { return members.empty(); }
};

}  // namespace sparlab
