/// \file FemModel.hpp
/// \brief Aggregation of everything that defines one discrete structural model.
///
/// The model is a value type: it owns a mesh, the material, the idealisation
/// (plane stress, plane strain or full 3-D elasticity), integration options,
/// the DOF partitioning and the load cases. It deliberately performs no
/// numerics — assembly, solution, stress recovery and optimisation all consume
/// a `const FemModel&`.
///
/// The spatial dimension is the mesh's. A 2-D mesh takes a plane stress state
/// and a positive thickness; a solid mesh takes `StressState::ThreeDimensional`
/// and no thickness (the value must be 1); a shell mesh takes
/// `StressState::Shell` and a positive thickness, which regions may override
/// (`assign_thickness`). Any mismatch is rejected in the constructor, so a
/// deck cannot pair a solid mesh with a plane idealisation.
///
/// **Shell directors.** The director of a shell element at a node is the
/// surface's exact normal there when the mesh carries one (a generated
/// plate, cylinder or sphere, Mesh::node_normals), else the average of the
/// normals of the elements there that lie within the fold angle of one
/// another (a smooth surface), turned to face each element's own side; an
/// element across a fold - a T-junction, the corner of a box - keeps its own
/// normal there, and the global rotations at the node couple the two.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/elements/Element.hpp"
#include "sparlab/fem/BoundaryConditions.hpp"
#include "sparlab/fem/DofManager.hpp"
#include "sparlab/fem/LoadCaseData.hpp"
#include "sparlab/material/IsotropicMaterial.hpp"
#include "sparlab/mesh/Mesh.hpp"

#include <memory>
#include <string>
#include <vector>

namespace sparlab {

struct LinearSolverOptions;

/// Mass-matrix formulation.
enum class MassType {
  Consistent,  ///< M_e = rho t int N^T N dOmega (default)
  Lumped       ///< row-sum lumping of the consistent matrix (diagonal)
};

std::string to_string(MassType type);

/// Settings of a shell model.
struct ShellOptions {
  /// Drilling penalty relative to G t (Shell4.hpp).
  Scalar drilling_factor = 1.0e-3;
  /// Elements whose normals at a node differ by more than this share no
  /// director there [degrees].
  Scalar fold_angle_deg = 20.0;
};

/// A fully specified linear-elastic model plus its load cases.
class FemModel {
 public:
  FemModel(Mesh mesh, IsotropicMaterial material, Scalar thickness,
           StressState stress_state, IntegrationOptions integration);

  const Mesh& mesh() const { return mesh_; }
  /// Spatial dimension of the model (2 or 3): the translations per node.
  int dim() const { return mesh_.dim(); }
  /// Degrees of freedom per node: `dim()` for a continuum model, six (three
  /// translations, three rotations) for a shell or beam model.
  int dofs_per_node() const { return element_->dofs_per_node(); }
  /// The primary material (the deck's `material`); every element uses it
  /// unless `assign_material` gave it another.
  const IsotropicMaterial& material() const { return materials_.front(); }
  /// Out-of-plane thickness [m] of a 2-D model, the default thickness of a
  /// shell model; 1 for a solid model.
  Scalar thickness() const { return thickness_; }
  /// Thickness of element `e` [m]: a shell element's section, else
  /// `thickness()`.
  Scalar thickness_of(Index e) const {
    return element_thickness_.empty() ? thickness_
                                      : element_thickness_[static_cast<std::size_t>(e)];
  }
  /// True for a shell model (StressState::Shell, Shell4 elements).
  bool is_shell() const { return stress_state_ == StressState::Shell; }
  /// Give the listed shell elements the thickness `t` [m] (a later
  /// assignment overrides an earlier one).
  /// \throws ModelError on a model that is not a shell, for an element out
  ///         of range or a non-positive thickness.
  void assign_thickness(Scalar t, const std::vector<Index>& elements);
  /// Replace the shell settings (drilling factor, fold angle), rebuilding the
  /// element and the directors. \throws ModelError on a non-shell model.
  void set_shell_options(const ShellOptions& options);
  const ShellOptions& shell_options() const { return shell_options_; }
  /// What the element kernels take for element `e`: its nodal coordinates
  /// (dim x nodes), and for a shell its directors beneath them (6 x 4).
  Matrix element_geometry(Index e) const;
  StressState stress_state() const { return stress_state_; }
  const IntegrationOptions& integration() const { return integration_; }
  const Element& element() const { return *element_; }

  DofManager& dofs() { return dofs_; }
  const DofManager& dofs() const { return dofs_; }

  /// Constitutive matrix of the primary material [Pa] (3 x 3 or 6 x 6).
  const Matrix& constitutive() const { return d_.front(); }

  /// Replace the primary material (used by the material-stiffness sweep).
  void set_material(const IsotropicMaterial& material);

  /// Give the listed elements `material`, appended to the model's materials.
  /// A later assignment overrides an earlier one for the same elements.
  /// \throws ModelError for an element index out of range.
  void assign_material(const IsotropicMaterial& material, const std::vector<Index>& elements);

  const std::vector<IsotropicMaterial>& materials() const { return materials_; }
  int num_materials() const { return static_cast<int>(materials_.size()); }
  /// True when every element uses the primary material.
  bool single_material() const { return element_material_.empty(); }
  /// Index into `materials()` of element `e`'s material.
  int element_material(Index e) const {
    return element_material_.empty() ? 0 : element_material_[static_cast<std::size_t>(e)];
  }
  const IsotropicMaterial& material_of(Index e) const {
    return materials_[static_cast<std::size_t>(element_material(e))];
  }
  const Matrix& constitutive_of(Index e) const {
    return d_[static_cast<std::size_t>(element_material(e))];
  }

  std::vector<LoadCaseSpec>& load_case_specs() { return load_case_specs_; }
  const std::vector<LoadCaseSpec>& load_case_specs() const { return load_case_specs_; }

  std::vector<DisplacementConstraint>& constraints() { return constraints_; }
  const std::vector<DisplacementConstraint>& constraints() const { return constraints_; }

  /// Apply `constraints()` to the DOF manager and assemble one force vector per
  /// load case. Must be called after the constraints and load cases are set and
  /// before any assembly. Idempotent.
  /// \param require_load_cases when false, a model with no load case is
  ///        accepted. This is used for free-vibration analysis of a structure
  ///        extracted from a density field, where the original load regions may
  ///        no longer intersect the retained material.
  void finalize(bool require_load_cases = true);

  bool finalized() const { return finalized_; }

  /// Global force vectors [N], one per load case (index matches
  /// `load_case_specs()`): the mechanical, body and thermal parts together.
  const std::vector<Vector>& load_vectors() const;

  /// The parts of load case `l` and its temperature field (LoadCaseData.hpp).
  const LoadCaseData& load_case_data(std::size_t l) const;

  /// Linear-solver settings of the conduction solve that `finalize` runs for
  /// a load case with a conducted temperature field.
  void set_conduction_solver(const LinearSolverOptions& options);

  /// Load-case weights normalised to sum to one (used by the objective).
  std::vector<Scalar> normalised_weights() const;

  /// Total solid-material volume of the design domain [m^3].
  Scalar domain_volume() const;

  /// Per-element volume [m^3]: area * thickness in 2-D and for a shell, cell
  /// volume for a solid.
  Vector element_volumes() const;

 private:
  Mesh mesh_;
  std::vector<IsotropicMaterial> materials_;
  std::vector<int> element_material_;  ///< empty: every element uses materials_[0]
  Scalar thickness_;
  StressState stress_state_;
  IntegrationOptions integration_;
  std::unique_ptr<Element> element_;
  std::vector<Matrix> d_;              ///< one constitutive matrix per material
  DofManager dofs_;
  std::vector<DisplacementConstraint> constraints_;
  std::vector<LoadCaseSpec> load_case_specs_;
  std::vector<Vector> load_vectors_;
  std::vector<LoadCaseData> load_data_;
  std::shared_ptr<LinearSolverOptions> conduction_solver_;
  std::vector<Scalar> element_thickness_;  ///< empty: every element takes thickness_
  ShellOptions shell_options_;
  /// Shell directors per element, one column per node; empty otherwise.
  std::vector<Eigen::Matrix<Scalar, 3, 4>> directors_;
  bool finalized_ = false;

  void compute_directors();
};

}  // namespace sparlab
