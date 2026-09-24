/// \file SubMesh.hpp
/// \brief Extraction of a sub-mesh from an element subset.
///
/// A SIMP density field is a *relative material distribution*, not a solid
/// body. To analyse "the structure the topology represents" - for example to
/// compute its natural frequencies - the field has to be interpreted. SparLab
/// does this explicitly: elements with \f$\rho \ge \rho_{cut}\f$ are kept, the
/// largest face-connected group of kept elements is retained (edge-connected
/// in 2-D), and a new mesh is built from it. The interpretation (threshold
/// value, connectivity rule, discarded material fraction) is recorded in the
/// result so that nothing is presented as geometry without saying how it was
/// obtained.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <vector>

namespace sparlab {

struct SubMeshResult {
  Mesh mesh;
  std::vector<Index> element_map;    ///< new element index -> original index
  std::vector<Index> node_map;       ///< new node index -> original index
  std::vector<Index> old_to_new_node;///< original node -> new index, -1 if dropped
};

/// Build a mesh from the given element subset, renumbering nodes compactly.
/// \throws MeshError when the subset is empty.
SubMeshResult extract_element_subset(const Mesh& mesh, const std::vector<Index>& elements);

/// Connected components of the element graph using *shared faces* (two shared
/// nodes in 2-D, four in 3-D), which is the correct connectivity for load
/// transfer: elements that only touch at a corner or an edge form a hinge, not
/// a connection.
std::vector<std::vector<Index>> element_components_by_face(const Mesh& mesh);

/// Report of a density-field interpretation.
struct TopologyInterpretation {
  Scalar threshold = 0.5;
  Index elements_above_threshold = 0;
  Index elements_retained = 0;
  Index components_above_threshold = 0;
  Scalar volume_above_threshold = 0.0;  ///< [m^3]
  Scalar volume_retained = 0.0;         ///< [m^3]
  /// Material volume discarded because it sat in a smaller disconnected group.
  Scalar volume_discarded_as_islands = 0.0;
  SubMeshResult sub;
};

/// Threshold `density` at `threshold`, keep the largest face-connected group of
/// elements, and return both the extracted mesh and the bookkeeping needed to
/// report the interpretation honestly.
/// \param element_volumes per-element volume [m^3] used for the reported sums.
/// \throws MeshError if no element survives the threshold.
TopologyInterpretation interpret_density_as_solid(const Mesh& mesh, const Vector& density,
                                                  const Vector& element_volumes,
                                                  Scalar threshold,
                                                  bool largest_component_only = true);

}  // namespace sparlab
