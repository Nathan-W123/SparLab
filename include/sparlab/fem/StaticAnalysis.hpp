/// \file StaticAnalysis.hpp
/// \brief Linear static solution, reaction recovery and equilibrium checks.
///
/// The reduced system is
/// \f[ K_{ff}\, u_f = f_f - K_{fp}\, u_p, \f]
/// and reactions follow from the *full* residual
/// \f[ r = K u - f, \f]
/// whose entries at prescribed DOFs are the support reactions and whose entries
/// at free DOFs are the (tiny) solver residual. This gives exact reactions with
/// no separate assembly.
///
/// Energy definitions used consistently across the code base:
///   * strain energy \f$ U = \tfrac{1}{2} u^T K u \f$ [J],
///   * compliance \f$ C = f^T u \f$ [J] with \f$f\f$ the *applied* load vector.
/// For homogeneous Dirichlet data \f$C = 2U\f$; the two are reported separately
/// so the identity can be checked (it is, in the test suite).
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/LinearSolver.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Global force / moment balance of a solved load case. Moments are taken
/// about the origin; on a 2-D model only the z component is non-zero.
struct EquilibriumCheck {
  Vector3 applied_force = Vector3::Zero();    ///< sum of applied nodal forces [N]
  Vector3 reaction_force = Vector3::Zero();   ///< sum of support reactions [N]
  Vector3 force_residual = Vector3::Zero();   ///< applied + reaction [N]
  Vector3 applied_moment = Vector3::Zero();   ///< about the origin [N m]
  Vector3 reaction_moment = Vector3::Zero();  ///< about the origin [N m]
  Vector3 moment_residual = Vector3::Zero();  ///< applied + reaction [N m]
  Scalar relative_force_error = 0.0;          ///< |residual| / max(|applied|, tiny)
  Scalar relative_moment_error = 0.0;         ///< normalised by |applied moment| scale
};

/// Result of one solved load case.
struct StaticSolution {
  std::string load_case_name;
  Scalar weight = 1.0;
  Vector displacement;                ///< full length, [m]
  Vector reactions;                   ///< full length, non-zero at prescribed DOFs [N]
  Scalar compliance = 0.0;            ///< f^T u [J]
  Scalar strain_energy = 0.0;         ///< 0.5 u^T K u [J]
  Scalar max_displacement_magnitude = 0.0;  ///< [m]
  Index max_displacement_node = -1;
  Scalar scaled_residual = 0.0;       ///< ||K_ff u_f - rhs|| / ||rhs||
  int solver_iterations = 0;
  EquilibriumCheck equilibrium;
};

struct StaticAnalysisOptions {
  LinearSolverOptions linear;
  /// Run `require_well_posed` before assembling. Disabled internally by the
  /// topology optimizer, whose SIMP stiffness floor keeps K non-singular even
  /// in void regions.
  bool check_model = true;
  /// Relative force-balance tolerance; exceeding it raises SolverError.
  Scalar equilibrium_tolerance = 1.0e-6;
};

/// Linear static solver. One instance reuses a single factorisation of
/// \f$K_{ff}\f$ across all load cases, which is what makes multi-load-case
/// topology optimisation affordable.
class StaticAnalysis {
 public:
  StaticAnalysis(const FemModel& model, const Assembler& assembler,
                 StaticAnalysisOptions options = StaticAnalysisOptions());

  /// Factorise \f$K_{ff}\f$ for the given per-element stiffness factors.
  /// \param stiffness_scale optional \f$E(\rho_e)/E_0\f$ (length num_elements).
  void prepare(const Vector* stiffness_scale = nullptr);

  /// Solve every load case of the model. `prepare` is called automatically if
  /// it has not been called since the last change of `stiffness_scale`.
  std::vector<StaticSolution> solve_all(const Vector* stiffness_scale = nullptr);

  /// Solve one right-hand side with the current factorisation, returning the
  /// full displacement vector.
  Vector solve_load_vector(const Vector& applied_force);

  /// Solve \f$K_{ff}\,\lambda_f = r_f\f$ with the current factorisation and
  /// zeros at the prescribed DOFs: the adjoint problem of any functional
  /// evaluated on this model, whatever the prescribed displacements.
  Vector solve_homogeneous(const Vector& rhs);

  /// Weighted compliance \f$\sum_l w_l\, f_l^T u_l\f$ with normalised weights.
  static Scalar weighted_compliance(const std::vector<StaticSolution>& solutions,
                                    const std::vector<Scalar>& normalised_weights);

  /// The assembled global stiffness matrix of the last `prepare` call.
  const SparseMatrix& stiffness() const { return k_full_; }

  const StaticAnalysisOptions& options() const { return options_; }

 private:
  StaticSolution build_solution(const std::string& name, Scalar weight,
                                const Vector& applied_force);

  const FemModel& model_;
  const Assembler& assembler_;
  StaticAnalysisOptions options_;
  std::unique_ptr<LinearSolver> solver_;
  SparseMatrix k_full_;
  SparseMatrix k_ff_;
  SparseMatrix k_fp_;
  Vector prescribed_;
  bool prepared_ = false;
};

}  // namespace sparlab
