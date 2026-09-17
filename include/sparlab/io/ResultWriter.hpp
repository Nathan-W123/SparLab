/// \file ResultWriter.hpp
/// \brief Writes every run artefact into one output directory.
///
/// Layout produced for a case named `<case>` in `<out>`:
/// \code
///   <out>/summary.json              scalars, tolerances, timings, provenance
///   <out>/config.json               verbatim echo of the input deck
///   <out>/mesh.json                 nodes, connectivity, BCs, applied loads
///   <out>/displacement_<lc>.csv     nodal displacements per load case
///   <out>/stress_<lc>.csv           element strains/stresses per load case
///   <out>/reactions_<lc>.csv        support reactions per load case
///   <out>/fields_<lc>.vtk           ParaView fields per load case
///   <out>/modes.csv                 eigenvalues, frequencies, residuals
///   <out>/mode_shapes.csv           M-orthonormal mode shapes
///   <out>/mode_<k>.vtk              one file per mode shape
///   <out>/history.csv               optimisation iteration history
///   <out>/density_final.csv         final design and physical density
///   <out>/density_history.csv       density snapshots for the animation
///   <out>/topology_<lc>.vtk         optimised design fields
/// \endcode
/// Every file is plain text. `summary.json` is the single source of truth for
/// the numbers quoted in the documentation.
#pragma once

#include "sparlab/core/Timer.hpp"
#include "sparlab/core/Types.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/mesh/SubMesh.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <string>
#include <vector>

namespace sparlab {

class ResultWriter {
 public:
  /// Create (or reuse) the output directory.
  ResultWriter(std::string directory, const Configuration& config);

  const std::string& directory() const { return directory_; }

  /// Echo the input deck so a result directory is self-describing.
  void write_config() const;

  /// Nodes, connectivity, prescribed DOFs and the applied load vectors.
  void write_mesh(const FemModel& model) const;

  /// Nodal displacements of one load case.
  void write_displacement(const Mesh& mesh, const std::string& load_case,
                          const Vector& displacement) const;

  /// Element strain/stress data of one load case.
  void write_stress(const Mesh& mesh, const std::string& load_case,
                    const StressField& field, const Vector* density = nullptr) const;

  /// Non-zero support reactions of one load case.
  void write_reactions(const Mesh& mesh, const DofManager& dofs,
                       const std::string& load_case, const Vector& reactions) const;

  /// Combined VTK file for one load case.
  void write_static_vtk(const Mesh& mesh, const std::string& load_case,
                        const Vector& displacement, const StressField& field,
                        const Vector* density, const Vector* stiffness_factor) const;

  /// Eigenvalues, frequencies, residuals and (optionally) mode shapes.
  void write_modal(const Mesh& mesh, const ModalResult& modal,
                   const std::string& tag = "") const;

  /// Optimisation iteration history.
  void write_history(const TopologyOptimizationResult& result) const;

  /// Final design, physical density and per-element energy.
  void write_density(const Mesh& mesh, const DesignDomain& domain,
                     const TopologyOptimizationResult& result) const;

  /// Density snapshots used by the evolution animation.
  /// \param max_frames snapshots are strided down to at most this many frames.
  void write_density_history(const TopologyOptimizationResult& result,
                             int max_frames = 60) const;

  /// Write an arbitrary JSON document into the directory.
  void write_json(const std::string& file_name, const json::Value& value) const;

  /// Path of a file inside the output directory.
  std::string file(const std::string& name) const;

 private:
  std::string directory_;
  const Configuration& config_;
};

/// Build the JSON summary of a static (plus optional modal) analysis.
json::Value make_static_summary(const Configuration& config, const FemModel& model,
                                const ModelDiagnostics& diagnostics,
                                const std::vector<StaticSolution>& solutions,
                                const std::vector<StressField>& stresses,
                                const ModalResult* modal, const TimingLedger& timings);

/// Build the JSON summary of a topology-optimisation run.
json::Value make_topology_summary(const Configuration& config, const FemModel& model,
                                  const DesignDomain& domain,
                                  const DensityFilter& filter,
                                  const TopologyOptimizationResult& result,
                                  const TopologyInterpretation* interpretation,
                                  const ModalResult* modal_initial,
                                  const ModalResult* modal_optimised,
                                  const ModalResult* modal_mass_matched,
                                  const TimingLedger& timings);

/// Common provenance block: version, build type, timestamp, tolerances.
json::Value make_provenance(const Configuration& config);

}  // namespace sparlab
