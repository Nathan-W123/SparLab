/// \file Config.hpp
/// \brief Declarative JSON input deck: schema, parsing and model construction.
///
/// An input deck fully determines a run. The schema is documented in
/// docs/configuration.md; the section names below mirror it one to one:
/// \code
/// {
///   "name": "cantilever_beam",
///   "mesh":     { "type": "structured_quad", "nx":.., "ny":.., "lx":.., "ly":.. },
///               // or { "type": "structured_hex", "nx","ny","nz", "lx","ly","lz" }
///   "material": { "youngs_modulus":.., "poisson_ratio":.., "density":.. },
///   "model":    { "thickness":.., "stress_state": "plane_stress" },
///   "boundary_conditions": [ { "fix": ["x","y"], "region": {..} } ],
///   "load_cases":          [ { "name":.., "weight":.., "point_loads": [..] } ],
///   "solver":   { "linear": {..}, "equilibrium_tolerance":.. },
///   "modal":    { "enabled": true, "num_modes":.. },
///   "topology": { "enabled": true, "volume_fraction":.., "filter": {..} },
///   "output":   { "vtk": true, "csv": true }
/// }
/// \endcode
///
/// The mesh type fixes the spatial dimension of the whole deck: a
/// `structured_hex` mesh takes three-entry vectors, may fix `z`, and uses the
/// three-dimensional constitutive law; a `structured_quad` mesh takes
/// two-entry vectors and a plane stress state. Mixing the two is an error
/// rather than a silently truncated component.
///
/// Every numeric field is in SI units. Unknown keys are reported (not ignored),
/// because a misspelled key that silently takes its default is a direct route
/// to a wrong published number.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/LinearSolver.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/DesignDomain.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <optional>
#include <string>
#include <vector>

namespace sparlab {

/// Mesh generators a deck can name.
enum class MeshKind {
  StructuredQuad,  ///< "structured_quad": 2-D Q4 grid
  StructuredHex    ///< "structured_hex": 3-D Hex8 grid
};

std::string to_string(MeshKind kind);

struct ModalConfig {
  bool enabled = false;
  ModalAnalysisOptions options;
  /// Also run modal analysis on the density-thresholded topology (topopt runs).
  bool analyse_optimised_topology = true;
  /// Compare against a uniformly thinned domain of the same mass (2-D only:
  /// the argument rests on thickness scaling and has no 3-D counterpart).
  bool compare_mass_matched_baseline = true;
};

struct TopologyConfig {
  bool enabled = false;
  Scalar volume_fraction = 0.4;
  /// Starting density for free variables; negative means "use volume_fraction".
  Scalar initial_density = -1.0;
  FilterType filter_type = FilterType::Density;
  /// Filter radius in metres. When `filter_radius_elements` > 0 the radius is
  /// computed as that multiple of the mean element size instead.
  Scalar filter_radius = 0.0;
  Scalar filter_radius_elements = 1.5;
  TopologyOptimizerOptions optimizer;
  std::vector<PassiveRegionSpec> passive_regions;
};

struct OutputConfig {
  bool write_csv = true;
  bool write_vtk = true;
  bool write_density_history = true;
  bool write_mode_shapes = true;
};

/// A complete, validated run specification.
class Configuration {
 public:
  std::string name = "case";
  std::string description;
  std::string source_path;

  MeshKind mesh_kind = MeshKind::StructuredQuad;
  StructuredMeshSpec mesh_spec;
  Scalar thickness = 1.0;
  StressState stress_state = StressState::PlaneStress;
  IntegrationOptions integration;

  std::vector<DisplacementConstraint> constraints;
  std::vector<LoadCaseSpec> load_cases;

  StaticAnalysisOptions analysis;
  ModalConfig modal;
  TopologyConfig topology;
  OutputConfig output;

  /// The parsed document, retained so run summaries can record exactly what
  /// was requested.
  json::Value document;

  /// Spatial dimension implied by the mesh type (2 or 3).
  int dim() const { return mesh_kind == MeshKind::StructuredHex ? 3 : 2; }

  const IsotropicMaterial& material() const;
  void set_material(IsotropicMaterial material);

  /// Filter radius resolved against a concrete mesh [m].
  Scalar resolved_filter_radius(const Mesh& mesh) const;

 private:
  std::optional<IsotropicMaterial> material_;
};

/// Parse a configuration from an already-parsed document.
/// \param source name used in error messages.
/// \param strict when true, unknown keys raise ConfigError instead of a warning.
Configuration parse_configuration(const json::Value& document, const std::string& source,
                                  bool strict = false);

/// Read and parse a configuration file.
Configuration load_configuration(const std::string& path, bool strict = false);

/// Build the mesh described by a configuration (without the model).
Mesh build_mesh(const Configuration& config);

/// Build the mesh, material, boundary conditions and load cases described by a
/// configuration and finalise the model.
FemModel build_model(const Configuration& config);

/// Build the design domain for a topology optimisation run.
DesignDomain build_design_domain(const Configuration& config, const FemModel& model);

/// Parse a region specification (used by tests and by the config parser).
/// \param dim mesh dimension (2 or 3), which decides how many entries the
///        vector-valued fields must carry; 0 accepts either form.
SelectorGroup parse_region(const ConfigNode& node, const std::string& default_name,
                           int dim = 0);

}  // namespace sparlab
