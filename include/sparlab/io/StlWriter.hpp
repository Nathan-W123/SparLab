/// \file StlWriter.hpp
/// \brief Triangulated boundary surfaces of a mesh and binary STL export.
///
/// A topology run produces two geometries worth keeping as solids: the design
/// domain the optimiser started from ("before") and the structure obtained by
/// thresholding the final density field and keeping its largest connected
/// group ("after"). Both are written as watertight triangle surfaces:
///
///   * a 3-D (Hex8) mesh contributes its boundary faces, each split into two
///     triangles wound so the right-hand normal points out of the material;
///   * a 2-D (Q4) mesh is extruded through the model thickness along +z: the
///     element polygons become the two caps and the boundary edges become the
///     side walls.
///
/// The surface is indexed (shared vertices), which makes closure checkable
/// exactly: a surface is closed and consistently oriented when every directed
/// edge is matched by exactly one edge in the opposite direction. Its area and
/// the volume it encloses (divergence theorem) are integrated as well. On a
/// mesh with planar faces the enclosed volume equals the cell volume; on a
/// distorted hex mesh a cut face is a bilinear patch that two flat triangles
/// only approximate, so the two volumes then differ by the geometric error of
/// the triangulation, which is reported rather than hidden.
///
/// The STL carries geometry only. What it represents is stated where it is
/// written: a density field interpreted at a threshold, not a design that has
/// been checked for stress, manufacturability or anything else.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <array>
#include <string>
#include <vector>

namespace sparlab {

/// One oriented triangle [m].
struct Triangle {
  Vector3 a = Vector3::Zero();
  Vector3 b = Vector3::Zero();
  Vector3 c = Vector3::Zero();

  /// Unit normal by the right-hand rule (zero for a degenerate triangle).
  Vector3 normal() const;
  Scalar area() const;
};

/// Indexed triangle surface with shared vertices.
struct TriangleSurface {
  std::vector<Vector3> points;                 ///< [m]
  std::vector<std::array<Index, 3>> faces;     ///< counter-clockwise seen from outside

  Index num_triangles() const { return static_cast<Index>(faces.size()); }
  Triangle triangle(Index i) const;
  std::vector<Triangle> triangles() const;
};

/// Integrated properties and the closure check of a triangle surface.
struct SurfaceStats {
  Index num_triangles = 0;
  Index num_points = 0;
  Scalar area = 0.0;             ///< [m^2]
  Scalar enclosed_volume = 0.0;  ///< [m^3], divergence theorem; negative if inverted
  BoundingBox bounds;
  /// True when every directed edge is paired with exactly one reversed edge,
  /// i.e. the surface is closed and consistently oriented.
  bool closed = false;
  Index unmatched_edges = 0;     ///< directed edges without a reversed partner
};

/// Boundary surface of the mesh, extruded by `thickness` along +z when the
/// mesh is two-dimensional (ignored for a solid mesh).
TriangleSurface boundary_surface(const Mesh& mesh, Scalar thickness = 1.0);

/// Area, enclosed volume, bounds and the closure check.
SurfaceStats surface_stats(const TriangleSurface& surface);

/// Write a binary STL (80-byte header, uint32 count, 50 bytes per facet).
/// \param header at most 59 characters of description stored in the file.
/// \throws IoError on failure.
void write_stl(const std::string& path, const TriangleSurface& surface,
               const std::string& header = "SparLab");

/// Read a binary STL back as a flat triangle list (used by the tests).
/// \throws IoError for a malformed file.
std::vector<Triangle> read_stl(const std::string& path);

}  // namespace sparlab
