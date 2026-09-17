/// \file VtkWriter.hpp
/// \brief Legacy ASCII VTK (`UNSTRUCTURED_GRID`) export for ParaView.
///
/// The legacy format is used deliberately: it is a single self-describing text
/// file with no XML or compression dependencies, and ParaView, VisIt and
/// pyvista all read it. Q4 elements are written as VTK_QUAD (cell type 9) with
/// z = 0.
///
/// Field naming convention in the output file:
///   * point data: `displacement` (3-vector, z = 0), `displacement_magnitude`,
///     `nodal_von_mises`, and `mode_<k>` for mode shapes;
///   * cell data: `density`, `stiffness_factor`, `sigma_xx`, `sigma_yy`,
///     `sigma_xy`, `von_mises`, `principal_max`, `principal_min`,
///     `strain_energy`.
///
/// A density field written to VTK is the SIMP *material distribution*, not a
/// solid body; docs/limitations.md states how it should be read.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <string>
#include <utility>
#include <vector>

namespace sparlab {

/// Accumulates named point/cell arrays and writes one legacy VTK file.
class VtkWriter {
 public:
  explicit VtkWriter(const Mesh& mesh, std::string title = "SparLab field data");

  /// Add a scalar array over nodes (length num_nodes).
  VtkWriter& add_point_scalars(const std::string& name, const Vector& values);

  /// Add a 2-D vector array over nodes; written as a 3-vector with z = 0.
  /// \param values length 2*num_nodes in node-major DOF ordering.
  VtkWriter& add_point_vectors(const std::string& name, const Vector& values);

  /// Add a scalar array over elements (length num_elements).
  VtkWriter& add_cell_scalars(const std::string& name, const Vector& values);

  /// Write the file.
  /// \throws IoError on failure.
  void write(const std::string& path) const;

 private:
  const Mesh& mesh_;
  std::string title_;
  std::vector<std::pair<std::string, Vector>> point_scalars_;
  std::vector<std::pair<std::string, Vector>> point_vectors_;
  std::vector<std::pair<std::string, Vector>> cell_scalars_;
};

}  // namespace sparlab
