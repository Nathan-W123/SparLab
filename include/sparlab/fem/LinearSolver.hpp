/// \file LinearSolver.hpp
/// \brief Sparse and dense linear solvers with explicit failure reporting.
///
/// Every solver reports failure through `SolverError`; none of them returns a
/// silently wrong answer. After each solve the caller can (and by default does)
/// verify the scaled residual
/// \f$ \|A x - b\|_2 / \max(\|b\|_2, \epsilon) \f$ against a recorded
/// tolerance. Recorded tolerances appear in every run summary.
#pragma once

#include "sparlab/core/Types.hpp"

#include <memory>
#include <string>

namespace sparlab {

enum class LinearSolverType {
  SimplicialLdlt,     ///< sparse Cholesky (LDL^T) with AMD ordering - default
  SimplicialLlt,      ///< sparse Cholesky (LL^T), requires strict SPD
  SparseLu,           ///< sparse LU, works for indefinite systems
  ConjugateGradient,  ///< diagonal-preconditioned CG (matrix-free friendly)
  DenseLu             ///< dense partial-pivot LU; verification of small models
};

std::string to_string(LinearSolverType type);
LinearSolverType parse_linear_solver_type(const std::string& text);

struct LinearSolverOptions {
  LinearSolverType type = LinearSolverType::SimplicialLdlt;
  /// Relative residual tolerance requested from iterative solvers.
  Scalar iterative_tolerance = 1.0e-12;
  /// Iteration cap for iterative solvers; 0 selects Eigen's default (2n).
  int max_iterations = 0;
  /// Scaled-residual tolerance enforced after every solve. Exceeding it raises
  /// SolverError rather than returning a dubious displacement field.
  Scalar residual_tolerance = 1.0e-8;
  /// Smallest acceptable |pivot| / max|pivot| ratio in the Cholesky factor.
  /// Below this the system is reported as numerically singular.
  Scalar pivot_tolerance = 1.0e-14;
};

/// Polymorphic wrapper over the Eigen solvers.
class LinearSolver {
 public:
  virtual ~LinearSolver() = default;

  /// Factorise (or store, for iterative solvers) the system matrix.
  /// \throws SolverError on numerical failure, with a diagnosis of the likely
  ///         modelling cause.
  virtual void factorize(const SparseMatrix& a) = 0;

  /// Solve for one right-hand side. `factorize` must have been called.
  virtual Vector solve(const Vector& b) = 0;

  virtual std::string name() const = 0;

  /// Iterations consumed by the last solve (0 for direct solvers).
  virtual int last_iterations() const { return 0; }

  /// Solver-reported error estimate of the last solve (0 for direct solvers).
  virtual Scalar last_error() const { return 0.0; }
};

std::unique_ptr<LinearSolver> make_linear_solver(const LinearSolverOptions& options);

/// Scaled residual \f$\|Ax-b\| / \max(\|b\|, \text{tiny})\f$.
Scalar scaled_residual(const SparseMatrix& a, const Vector& x, const Vector& b);

/// Factorise, solve and verify the residual in one call.
/// \throws SolverError when the residual exceeds `options.residual_tolerance`
///         or the solution is not finite.
Vector solve_and_verify(const SparseMatrix& a, const Vector& b,
                        const LinearSolverOptions& options,
                        Scalar* out_residual = nullptr);

/// Test whether a symmetric sparse matrix is numerically positive definite by
/// inspecting the LDL^T pivots.
/// \param min_pivot_ratio smallest accepted pivot / largest pivot.
bool is_positive_definite(const SparseMatrix& a, Scalar min_pivot_ratio = 1.0e-14,
                          Scalar* out_min_pivot_ratio = nullptr);

}  // namespace sparlab
