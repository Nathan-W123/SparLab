#include "sparlab/mesh/SubMesh.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <sstream>
#include <utility>

namespace sparlab {
namespace {

class DisjointSet {
 public:
  explicit DisjointSet(std::size_t n) : parent_(n) {
    std::iota(parent_.begin(), parent_.end(), static_cast<Index>(0));
  }
  Index find(Index a) {
    while (parent_[static_cast<std::size_t>(a)] != a) {
      parent_[static_cast<std::size_t>(a)] =
          parent_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(a)])];
      a = parent_[static_cast<std::size_t>(a)];
    }
    return a;
  }
  void unite(Index a, Index b) {
    a = find(a);
    b = find(b);
    if (a != b) parent_[static_cast<std::size_t>(b)] = a;
  }

 private:
  std::vector<Index> parent_;
};

/// Components of an element subset connected through shared faces (edges in
/// 2-D). A face is identified by its sorted node list.
std::vector<std::vector<Index>> components_by_face(const Mesh& mesh,
                                                   const std::vector<Index>& subset) {
  const std::vector<std::vector<int>>& local = element_local_faces(mesh.element_type());
  std::map<std::vector<Index>, Index> first_owner;  // face -> local subset index
  DisjointSet ds(subset.size());

  std::vector<Index> key;
  for (std::size_t s = 0; s < subset.size(); ++s) {
    const Index* nodes = mesh.element_nodes(subset[s]);
    for (const std::vector<int>& face : local) {
      key.clear();
      for (int a : face) key.push_back(nodes[a]);
      std::sort(key.begin(), key.end());
      const auto it = first_owner.find(key);
      if (it == first_owner.end()) {
        first_owner.emplace(key, static_cast<Index>(s));
      } else {
        ds.unite(it->second, static_cast<Index>(s));
      }
    }
  }

  std::map<Index, Index> root_to_component;
  std::vector<std::vector<Index>> components;
  for (std::size_t s = 0; s < subset.size(); ++s) {
    const Index root = ds.find(static_cast<Index>(s));
    auto it = root_to_component.find(root);
    if (it == root_to_component.end()) {
      it = root_to_component.emplace(root, static_cast<Index>(components.size())).first;
      components.emplace_back();
    }
    components[static_cast<std::size_t>(it->second)].push_back(subset[s]);
  }
  return components;
}

}  // namespace

std::vector<std::vector<Index>> element_components_by_face(const Mesh& mesh) {
  std::vector<Index> all(static_cast<std::size_t>(mesh.num_elements()));
  std::iota(all.begin(), all.end(), static_cast<Index>(0));
  return components_by_face(mesh, all);
}

SubMeshResult extract_element_subset(const Mesh& mesh,
                                     const std::vector<Index>& elements) {
  if (elements.empty()) {
    throw MeshError("cannot extract a sub-mesh from an empty element subset");
  }
  const int npe = mesh.nodes_per_elem();

  SubMeshResult out;
  out.old_to_new_node.assign(static_cast<std::size_t>(mesh.num_nodes()), -1);
  out.element_map = elements;

  std::vector<Index> connectivity;
  connectivity.reserve(elements.size() * static_cast<std::size_t>(npe));
  for (Index e : elements) {
    if (e < 0 || e >= mesh.num_elements()) {
      std::ostringstream os;
      os << "sub-mesh element index " << e << " is outside [0, " << mesh.num_elements() - 1
         << "]";
      throw MeshError(os.str());
    }
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 0; a < npe; ++a) {
      Index& mapped = out.old_to_new_node[static_cast<std::size_t>(nodes[a])];
      if (mapped < 0) {
        mapped = static_cast<Index>(out.node_map.size());
        out.node_map.push_back(nodes[a]);
      }
      connectivity.push_back(mapped);
    }
  }

  Matrix coords(mesh.dim(), static_cast<Eigen::Index>(out.node_map.size()));
  for (std::size_t n = 0; n < out.node_map.size(); ++n) {
    coords.col(static_cast<Eigen::Index>(n)) = mesh.coordinates().col(out.node_map[n]);
  }

  out.mesh = Mesh(std::move(coords), std::move(connectivity), mesh.element_type());
  out.mesh.validate();
  return out;
}

TopologyInterpretation interpret_density_as_solid(const Mesh& mesh, const Vector& density,
                                                  const Vector& element_volumes,
                                                  Scalar threshold,
                                                  bool largest_component_only) {
  if (density.size() != mesh.num_elements()) {
    std::ostringstream os;
    os << "density vector has length " << density.size() << " but the mesh has "
       << mesh.num_elements() << " elements";
    throw MeshError(os.str());
  }
  if (element_volumes.size() != mesh.num_elements()) {
    throw MeshError("element volume vector length does not match the mesh");
  }
  if (!(threshold > 0.0 && threshold < 1.0)) {
    std::ostringstream os;
    os << "density threshold must lie strictly between 0 and 1 (got " << threshold << ")";
    throw ConfigError(os.str());
  }

  TopologyInterpretation report;
  report.threshold = threshold;

  std::vector<Index> kept;
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    if (density(e) >= threshold) {
      kept.push_back(e);
      report.volume_above_threshold += element_volumes(e);
    }
  }
  report.elements_above_threshold = static_cast<Index>(kept.size());
  if (kept.empty()) {
    std::ostringstream os;
    os << "no element reaches the density threshold " << threshold
       << " (max density is " << density.maxCoeff()
       << "); lower the threshold or increase the volume fraction";
    throw MeshError(os.str());
  }

  const auto components = components_by_face(mesh, kept);
  report.components_above_threshold = static_cast<Index>(components.size());

  std::vector<Index> retained = kept;
  if (largest_component_only && components.size() > 1) {
    std::size_t best = 0;
    Scalar best_volume = -1.0;
    for (std::size_t c = 0; c < components.size(); ++c) {
      Scalar v = 0.0;
      for (Index e : components[c]) v += element_volumes(e);
      if (v > best_volume) {
        best_volume = v;
        best = c;
      }
    }
    retained = components[best];
    std::sort(retained.begin(), retained.end());
  }

  for (Index e : retained) report.volume_retained += element_volumes(e);
  report.elements_retained = static_cast<Index>(retained.size());
  report.volume_discarded_as_islands =
      report.volume_above_threshold - report.volume_retained;

  report.sub = extract_element_subset(mesh, retained);

  log::info("density interpretation at rho >= ", threshold, ": ",
            report.elements_above_threshold, " of ", mesh.num_elements(),
            " elements above threshold in ", report.components_above_threshold,
            " group(s); retained ", report.elements_retained, " elements (",
            report.volume_retained, " m^3), discarded ",
            report.volume_discarded_as_islands, " m^3 as disconnected islands");

  // A thresholded design that falls into several pieces is not a structure, and
  // the pieces that were dropped are material the volume constraint was still
  // charged for. It is usually the signature of an unfiltered or
  // under-resolved design, so say so at warning level rather than leaving it
  // to be noticed in the summary.
  if (report.components_above_threshold > 1) {
    const Scalar share =
        report.volume_above_threshold > 0.0
            ? report.volume_discarded_as_islands / report.volume_above_threshold
            : 0.0;
    log::warn("the density field thresholded at ", threshold, " falls into ",
              report.components_above_threshold,
              " disconnected groups; keeping the largest discards ",
              report.volume_discarded_as_islands, " m^3 (",
              100.0 * share,
              " % of the material above the threshold). Check the filter's "
              "average support first - a support near one element means no "
              "filtering and a checkerboard - and otherwise expect members thin "
              "enough to break up at this threshold, which a larger filter "
              "radius or a higher volume fraction would thicken");
  }
  return report;
}

}  // namespace sparlab
