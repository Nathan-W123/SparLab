/// \file Mma.hpp
/// \brief Method of Moving Asymptotes (Svanberg 1987) with the primal-dual
///        interior-point subproblem solver of Svanberg's reference code.
///
/// MMA solves
/// \f[
///   \min_x f_0(x) \quad \text{s.t.}\quad f_i(x) \le 0,\ i = 1..m,\quad
///   x_{min} \le x \le x_{max}
/// \f]
/// by replacing every function with a separable convex approximation built
/// from its gradient and a pair of moving asymptotes \f$L_j < x_j < U_j\f$:
/// \f[
///   \tilde f_i(x) = r_i + \sum_j \left( \frac{p_{ij}}{U_j - x_j} +
///   \frac{q_{ij}}{x_j - L_j} \right),
/// \f]
/// with \f$p_{ij}\f$ carrying the positive and \f$q_{ij}\f$ the negative part
/// of \f$\partial f_i/\partial x_j\f$ (plus a small convexifying term). The
/// asymptotes move outward while a variable keeps its direction and inward
/// when it oscillates, which is what damps the iteration. The subproblem is
/// solved to its KKT conditions with Svanberg's primal-dual Newton method on a
/// sequence of decreasing barrier parameters (`subsolv`), including the
/// artificial variables \f$y_i \ge 0\f$ and \f$z \ge 0\f$ that keep every
/// subproblem feasible; with the default \f$c_i = 1000\f$ the constraints
/// are met to the subproblem tolerance whenever that is possible.
///
/// Unlike the optimality-criteria update, MMA handles any number of
/// constraints, which is what a stress-constrained problem needs. For the
/// single volume constraint it is the more expensive method and is not the
/// default.
///
/// Reference: K. Svanberg, "The method of moving asymptotes - a new method
/// for structural optimization", Int. J. Numer. Meth. Engng 24 (1987) 359-373,
/// and the freely described `mmasub` / `subsolv` algorithms.
#pragma once

#include "sparlab/core/Types.hpp"

namespace sparlab {

struct MmaOptions {
  /// Move limit per iteration as a fraction of \f$x_{max} - x_{min}\f$.
  Scalar move_limit = 0.2;
  /// Initial asymptote distance as a fraction of \f$x_{max} - x_{min}\f$.
  Scalar asymptote_init = 0.5;
  /// Asymptote expansion when a variable keeps moving the same way.
  Scalar asymptote_increase = 1.2;
  /// Asymptote contraction when a variable oscillates.
  Scalar asymptote_decrease = 0.7;
  /// Fraction of the asymptote distance kept as a margin (Svanberg's albefa).
  Scalar albefa = 0.1;
  /// Minimum convexifying curvature (Svanberg's raa0).
  Scalar raa0 = 1.0e-5;
  /// Weight of the artificial variable z in the objective (a0).
  Scalar a0 = 1.0;
  /// Linear penalty on the constraint slacks y_i; large values make the
  /// constraints effectively hard.
  Scalar c = 1000.0;
  /// Quadratic penalty on the constraint slacks y_i.
  Scalar d = 1.0;
  /// Final barrier parameter of the subproblem solver.
  Scalar epsimin = 1.0e-7;
  /// A constraint whose value, or whose gradient times the move limit,
  /// exceeds this magnitude is scaled down to it for the subproblem (and its
  /// multiplier scaled back), so a badly violated constraint cannot make the
  /// absolute residual target of the Newton solver unreachable. Well-scaled
  /// problems never trigger it.
  Scalar constraint_scale_cap = 10.0;
  /// Newton iterations allowed per barrier level.
  int max_inner_iterations = 200;

  void validate() const;
};

/// Outcome of one MMA update.
struct MmaStep {
  Vector x;               ///< new design variables (n)
  Vector lambda;          ///< constraint multipliers (m), >= 0, in the caller's scaling
  Vector y;               ///< constraint slacks (m); > 0 marks an unmet constraint
  Vector constraint_scale;  ///< internal scale applied per constraint (1 = untouched)
  Scalar z = 0.0;         ///< artificial objective variable
  Scalar max_change = 0.0;      ///< \f$\|x_{new} - x\|_\infty\f$
  int subproblem_iterations = 0;
  Scalar subproblem_residual = 0.0;
  bool subproblem_converged = false;
};

/// Stateful MMA driver: keeps the asymptotes and the two previous designs.
class MmaOptimizer {
 public:
  /// \param num_variables n, \param num_constraints m (>= 1).
  /// \param lower / upper fixed variable bounds (length n).
  MmaOptimizer(Index num_variables, Index num_constraints, const Vector& lower,
               const Vector& upper, MmaOptions options = MmaOptions());

  /// One outer MMA iteration.
  /// \param x current design (n).
  /// \param f0 objective value (only its gradient enters the approximation).
  /// \param df0dx objective gradient (n).
  /// \param fval constraint values (m), feasible when <= 0.
  /// \param dfdx constraint gradients, m x n.
  /// \throws ConvergenceError when the subproblem solver does not reach its
  ///         tolerance (never silently accepted).
  MmaStep update(const Vector& x, Scalar f0, const Vector& df0dx, const Vector& fval,
                 const Matrix& dfdx);

  Index num_variables() const { return n_; }
  Index num_constraints() const { return m_; }
  int iteration() const { return iteration_; }
  const Vector& lower_asymptotes() const { return low_; }
  const Vector& upper_asymptotes() const { return upp_; }
  const MmaOptions& options() const { return options_; }

 private:
  MmaStep solve_subproblem(const Vector& alfa, const Vector& beta, const Vector& p0,
                           const Vector& q0, const Matrix& p, const Matrix& q,
                           const Vector& b) const;

  Index n_;
  Index m_;
  Vector xmin_;
  Vector xmax_;
  MmaOptions options_;
  int iteration_ = 0;
  Vector xold1_;
  Vector xold2_;
  Vector low_;
  Vector upp_;
};

}  // namespace sparlab
