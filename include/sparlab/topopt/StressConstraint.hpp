/// \file StressConstraint.hpp
/// \brief Aggregated von Mises stress constraint for SIMP with the qp relaxation
///        and its adjoint sensitivities.
///
/// **Relaxed stress.** The stress of the solid material in element e is
/// \f$\sigma_e = D B_e u_e\f$ (evaluated at the element centre). A void
/// element with the SIMP stiffness floor still carries a finite solid-material
/// stress, so a constraint on \f$\sigma_e\f$ directly would forbid removing
/// material (the singularity problem). The qp relaxation (Bruggi 2008; Le et al.
/// 2010) constrains instead
/// \f[
///   \sigma^{rel}_e = \tilde\rho_e^{\,q}\, \sigma_{vm}(\sigma_e),\qquad q < p,
/// \f]
/// which vanishes with the density while staying the true stress at
/// \f$\tilde\rho_e = 1\f$.
///
/// **Aggregation.** One constraint per load case replaces the element-wise
/// ones with a p-norm (Le et al. 2010):
/// \f[
///   s_e = \frac{\sigma^{rel}_e}{\sigma_{lim}},\qquad
///   g_{PN} = \Big(\sum_e s_e^{P}\Big)^{1/P},\qquad
///   g = c\, g_{PN} - 1 \le 0,
/// \f]
/// with the scale \f$c\f$ updated each iteration so that \f$c\,g_{PN}\f$
/// tracks the actual maximum, \f$c_k = \alpha\, s_{max}/g_{PN} + (1-\alpha)
/// c_{k-1}\f$, evaluated on the previous iterate. The p-norm bounds the
/// maximum from above, \f$s_{max} \le g_{PN} \le n^{1/P} s_{max}\f$, so the
/// scale converges to a value below one.
///
/// **Sensitivities.** With \f$K u = f\f$ and the von Mises matrix \f$V\f$
/// (\f$\sigma_{vm}^2 = \sigma^T V \sigma\f$),
/// \f[
///   \frac{\partial g_{PN}}{\partial\tilde\rho_e} =
///   \Big(\frac{s_e}{g_{PN}}\Big)^{P-1}
///   \frac{q\,\tilde\rho_e^{\,q-1}\sigma_{vm,e}}{\sigma_{lim}}
///   - \lambda_e^T \frac{\partial K_e}{\partial\tilde\rho_e} u_e,
///   \qquad
///   K\lambda = \sum_j A_j^T\,
///   \Big(\frac{s_j}{g_{PN}}\Big)^{P-1}\frac{\tilde\rho_j^{\,q}}{\sigma_{lim}}
///   \frac{(V\sigma_j)^T D B_j}{\sigma_{vm,j}},
/// \f]
/// one adjoint solve per load case with the factorisation already in hand.
/// The design-variable gradient follows from the density filter's exact chain
/// rule, which is why the constraint is only offered with that filter (or
/// none). Every gradient is verified against central differences in the test
/// suite and the verification app.
///
/// **What the constraint does and does not guarantee.** It bounds the
/// relaxed, aggregated stress of the SIMP model. The reported check that
/// matters for a design is the re-solve of the thresholded structure with
/// full material, whose maximum von Mises stress is compared with the limit
/// in the run summary; the two differ, and the summary says by how much.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/topopt/DensityFilter.hpp"
#include "sparlab/topopt/Sensitivity.hpp"
#include "sparlab/topopt/SimpInterpolation.hpp"

#include <vector>

namespace sparlab {

struct StressConstraintOptions {
  bool enabled = false;
  Scalar limit = 0.0;            ///< allowable von Mises stress \f$\sigma_{lim}\f$ [Pa]
  Scalar p_norm = 8.0;           ///< aggregation exponent P
  Scalar relaxation = 0.5;       ///< qp exponent q (must stay below the SIMP penalty)
  Scalar scaling_blend = 0.5;    ///< alpha of the adaptive scale update
  /// Constraint value above which the design counts as infeasible.
  Scalar feasibility_tolerance = 1.0e-3;

  void validate(const SimpOptions& simp) const;
};

/// Stress data of one load case at one design.
struct StressEvaluation {
  Vector solid_von_mises;        ///< sigma_vm of the solid material per element [Pa]
  Vector relaxed_von_mises;      ///< rho^q sigma_vm per element [Pa]
  Scalar max_relaxed_ratio = 0.0;   ///< max_e sigma_rel / sigma_lim
  Index max_element = -1;
  Scalar p_norm_ratio = 0.0;     ///< g_PN (dimensionless)
  Scalar scale = 1.0;            ///< c used for this evaluation
  Scalar constraint = 0.0;       ///< c g_PN - 1
  Vector dg_dphysical;           ///< d(constraint)/d rho (when requested)
  Vector dg_dx;                  ///< d(constraint)/d x through the filter
};

/// von Mises matrix V for the model's stress state, sized nv x nv, so that
/// \f$\sigma_{vm}^2 = \sigma^T V \sigma\f$ for a Voigt stress vector.
Matrix von_mises_matrix(StressState state, Scalar poisson);

class StressConstraint {
 public:
  StressConstraint(const FemModel& model, const Assembler& assembler,
                   const DensityFilter& filter, StressConstraintOptions options);

  const StressConstraintOptions& options() const { return options_; }

  /// Evaluate the constraint of load case `load_case` at the design behind
  /// `eval` (which must come from `objective.evaluate` on the same design).
  /// \param scale the scale c to apply.
  /// \param need_gradients compute the adjoint sensitivities.
  StressEvaluation evaluate(ComplianceObjective& objective, const ObjectiveEvaluation& eval,
                            std::size_t load_case, Scalar scale, bool need_gradients) const;

  /// Adaptive scale for the next iteration from the current maximum and p-norm.
  Scalar next_scale(Scalar previous_scale, const StressEvaluation& current) const;

 private:
  const FemModel& model_;
  const Assembler& assembler_;
  const DensityFilter& filter_;
  StressConstraintOptions options_;
  Matrix v_;                      ///< von Mises matrix
  Matrix db_uniform_;             ///< D B at the centre, cached for uniform meshes
  bool uniform_ = false;

  Matrix db_at_centre(Index e) const;
};

}  // namespace sparlab
