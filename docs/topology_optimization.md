# Topology optimisation

## 1. The problem

Minimum compliance under a volume constraint, with per-element density design
variables, on a plane (Q4) or a solid (Hex8) mesh:

```
  min over x      c(x) = sum_l  w_l  f_l^T u_l
  subject to      K(rho_tilde(x)) u_l = f_l          for every load case l
                  g(x) = rho_tilde^T v - nu V  <=  0
                  0 <= x_e <= 1                      (tightened on passive elements)

  optionally      g_l(x) = c_l  sigma_PN,l(x) / sigma_lim  -  1  <=  0   per load case
```

The optional last line is the aggregated stress constraint of section 5c; it
needs the MMA update of section 5b, since optimality criteria can carry only
the volume constraint.

where

| symbol | meaning |
|--------|---------|
| `x` | design variables, one per element |
| `rho_tilde` | *physical* density, the filtered design variables |
| `w_l` | load-case weights, normalised to sum to 1 |
| `v_e` | element volume [m^3] |
| `V = sum_e v_e` | design-domain volume [m^3] |
| `nu` | volume-fraction target |

The constraint is written on the **physical** density, not on the design
variables: that is the volume the structure actually has.

## 2. SIMP material interpolation

SparLab uses the *modified* SIMP law (`topopt/SimpInterpolation.cpp`):

```
  E(rho) / E_0  =  epsilon + (1 - epsilon) rho^p

  d(E/E_0)/drho =  p (1 - epsilon) rho^(p-1)
```

with `epsilon = E_min / E_0` (default `1e-9`) and penalty `p >= 1`.

Why the modified form rather than `E = E_min + rho^p (E_0 - E_min)` with a lower
bound `rho >= rho_min`:

* the floor is **additive**, so a design variable may reach exactly 0 and 1
  while `K` stays positive definite. No artificial lower bound on the density is
  needed, and no design variable has to be held away from its true bound;
* the derivative is well defined on the whole closed interval;
* `E(0) = epsilon E_0` keeps the stiffness matrix invertible in void regions.
  With `epsilon = 1e-9` the smallest LDL^T pivot ratio is about `1e-9`, far above
  the `1e-14` singularity threshold, so a direct solver handles it comfortably.

`p = 1` makes the problem convex but leaves a grey design, because intermediate
density is then as efficient as a mixture of solid and void. `p = 3` is the usual
choice; the design study sweeps `p` from 1.5 to 5 and shows what each end costs.

### Mass interpolation

Mass is interpolated *separately* from stiffness, because penalising mass the
same way as stiffness is not physical:

```
  linear:           m(rho) = rho
  penalty_matched:  m(rho) = epsilon_m + (1 - epsilon_m) rho^p      (default)
```

The default matters for modal analysis of a SIMP design. With a linear mass law
and a penalised stiffness law, a void element has stiffness `~1e-9` but mass
`~1e-3`, so `omega^2 ~ E/m ~ 1e-6` and the void region produces spurious
low-frequency modes that have nothing to do with the structure - the classical
SIMP eigenvalue artefact. Matching the mass exponent to the stiffness exponent
keeps `E/m` bounded (of order 1) at every density, so no spurious mode appears.

This is why SparLab additionally reports, per mode, the fraction of modal
kinetic energy carried by elements below density 0.3: it is the diagnostic for
exactly this artefact. `docs/limitations.md` explains why the *headline* modal
comparison uses the explicitly extracted solid sub-mesh instead.

## 3. Filtering

Both filters use the same linear-hat kernel over element centroids
(`topopt/DensityFilter.cpp`):

```
  H_ei = max(0, r_min - || x_e - x_i ||)
```

assembled sparsely, then row-normalised so `sum_i Hhat_ei = 1`. Neighbour
search uses a uniform spatial hash with cell size `r_min`, so building the
operator is `O(n)` and works on unstructured meshes too.

### Density filter (default)

```
  rho_tilde = Hhat x
```

The physical density is a **linear** map of the design variables, so the chain
rule is exact:

```
  df/dx_i = sum_e Hhat_ei df/drho_tilde_e = ( Hhat^T grad_rho f )_i
```

Because the map is exact, the analytical gradient can be - and is - verified
against central differences to 2e-8 relative. `pull_back` is the exact adjoint
of `to_physical`, which the test suite checks directly via
`<H x, v> == <x, H^T v>`.

Row normalisation means a constant field passes through unchanged and the
filtered field is a convex average, so it never leaves `[min x, max x]`.

### Sensitivity filter (Sigmund 1997)

```
  filtered dc/dx_e = ( sum_i H_ei x_i dc/dx_i ) / ( max(gamma, x_e) sum_i H_ei )
```

with `gamma = 1e-3` guarding the division. Here the design variable *is* the
density and only the gradient is smoothed. It suppresses checkerboarding just as
effectively, and on the MBB benchmark it reaches a *lower* compliance (203.2 J
versus 218.8 J for the density filter at the same 60x20-equivalent radius),
because it does not impose a genuine minimum length scale - it merely smooths
the search direction.

It is a heuristic: the filtered gradient is not the gradient of any objective,
so a finite-difference check of the filtered value is meaningless by
construction. That is why SparLab defaults to the density filter and runs its
sensitivity verification on that path.

### No filter

Available for the study. Without filtering the design is mesh dependent and
checkerboards: on the MBB benchmark it reaches 203.2 J with a grey level of
0.012, i.e. an almost perfectly 0/1 design - but the 0/1 pattern *is* the
checkerboard, which is a numerical artefact of the Q4 displacement field, not a
structure. A low grey level is therefore not by itself evidence of a good
design.

### Choosing the radius

`r_min` is the minimum-length-scale knob: the smallest member the design can
produce is roughly `2 r_min`. It is specified either directly in metres
(`filter.radius`) or as a multiple of the mean element size
(`filter.radius_elements`). Specifying it in **metres** is what makes the
optimum mesh convergent; specifying it in cells makes the length scale shrink
with the mesh, which is the classical mesh-dependence pathology. The design
study runs both arms side by side to show the difference.

## 4. Sensitivity analysis

Compliance is self-adjoint. Differentiating `K u = f` with `f` independent of the
design,

```
  dK/drho_e u + K du/drho_e = 0    =>    du/drho_e = -K^{-1} (dK/drho_e) u
```

and therefore

```
  dc/drho_e = sum_l w_l f_l^T du_l/drho_e
            = -sum_l w_l u_l^T (dK/drho_e) u_l
            = -sum_l w_l  d(E(rho_e)/E_0)/drho_e  u_{l,e}^T K_e^0 u_{l,e}
```

Three consequences:

* **no adjoint solve is needed.** One linear solve per load case yields both the
  objective and the exact gradient. That is what makes multi-load-case
  optimisation affordable;
* `dc/drho_e <= 0` always, since `u_e^T K_e^0 u_e >= 0` and the SIMP derivative
  is non-negative: adding material never increases compliance. The code warns if
  a positive entry ever appears;
* with the density filter, `grad_x c = Hhat^T grad_rho c`, and
  `grad_x g = Hhat^T v`.

Implemented in `topopt/Sensitivity.cpp`; verified in `sparlab_verify --study
sensitivity` (and `sensitivity-3d` on a Hex8 mesh) and regression-tested in
`tests/test_topopt.cpp`. The same class keeps the factorisation of the last
evaluation, so a non-self-adjoint quantity - the stress aggregate of section
5c - gets its adjoint solve from `ComplianceObjective::solve_adjoint` at the
cost of one back-substitution per load case.

## 5. Optimality criteria

With a single inequality constraint the KKT stationarity condition is separable.
For every variable strictly inside its bounds,

```
  B_e  =  ( -dc/dx_e ) / ( lambda dg/dx_e )  =  1
```

which gives the classical fixed-point update (`topopt/OptimalityCriteria.cpp`):

```
  x_e^{k+1} = clip(  x_e^k B_e^eta ,
                     max( x_lower_e , x_e^k - m ) ,
                     min( x_upper_e , x_e^k + m )  )
```

with damping `eta = 1/2` and move limit `m = 0.2` by default. The multiplier
`lambda` is found by bisection: the mapped volume is monotonically decreasing in
`lambda`, so bisection is globally convergent. SparLab bisects on
`sqrt(lo * hi)` rather than the arithmetic mean because `lambda` spans many
decades, and it widens the initial bracket automatically before bisecting. The
achieved volume is then within `1e-10` relative of the target - every benchmark
reports a volume-constraint violation below `1e-10`.

If the bracket cannot be made to straddle the target, `ConvergenceError` names
the reachable volume range and the two causes (move limit too small to reach the
target from the current design, or passive regions making the target
unreachable).

**When OC is the right choice.** One constraint, a separable objective in the
SIMP sense, and convergence in a few hundred iterations with no tuning. It is
the default (`optimizer.method: "oc"`) and every compliance-only benchmark
uses it. Its limitation is exactly one inequality constraint, which is what
the next section removes.

## 5b. Method of moving asymptotes

`optimizer.method: "mma"` replaces the OC update with Svanberg's method of
moving asymptotes (`topopt/Mma.cpp`), which handles any number of
constraints - the volume plus one stress constraint per load case here.

At iterate `x^k`, each function `f_i` (objective `i = 0`, constraints
`i = 1..m`) is replaced by the separable convex approximation

```
  f_i(x) ~ r_i + sum_j [ p_ij / (U_j - x_j) + q_ij / (x_j - L_j) ]

  p_ij = (U_j - x_j^k)^2 [ max(df_i/dx_j, 0) + 1e-3 |df_i/dx_j| + raa0 / (xmax_j - xmin_j) ]
  q_ij = (x_j^k - L_j)^2 [ max(-df_i/dx_j, 0) + 1e-3 |df_i/dx_j| + raa0 / (xmax_j - xmin_j) ]
```

built from the asymptotes `L_j < x_j^k < U_j`. The asymptotes start at
`x^k -+ 0.5 (xmax - xmin)` and move on the history of each variable: widened
by 1.2 when a variable has moved the same way twice, tightened by 0.7 when it
has oscillated, and never closer than `0.01` nor further than `10` times the
range. The subproblem is solved on the box
`max(xmin, L + 0.1 (x^k - L), x^k - m) <= x <= min(xmax, U - 0.1 (U - x^k), x^k + m)`,
with `m` the move limit, by the primal-dual interior-point method of
Svanberg's reference implementation (`subsolv`): Newton on the perturbed KKT
conditions with the barrier parameter driven from 1 down to `1e-7`, line
search on the residual, and the constraint multipliers `lambda_i` returned
with the step.

Two details are SparLab's rather than the textbook's:

* **scaling.** The objective is divided by its value at the first iteration
  and the volume constraint is written as `V/V_target - 1`, so every function
  is O(1). A constraint that is nevertheless violated by orders of magnitude
  (a stress limit far below what any feasible design can reach) would make
  the subproblem's absolute residual target unreachable, so such a row is
  scaled down to a cap of 10 and its multiplier scaled back afterwards; the
  scale is recorded per iteration;
* **the returned design.** Convergence is judged on the iterate whose
  functions were just evaluated and requires it to be feasible
  (`max_i g_i <= constraint_tolerance`, default `1e-4`) as well as either
  stationary (`max |dx| <= change_tolerance`) or stalled in the objective.
  That evaluated iterate is what the run returns. The unevaluated update can
  differ from it by up to the move limit under the objective-stall rule and
  has no feasibility guarantee, and returning it would have reported numbers
  the check never saw.

Only free variables enter the subproblem; passive elements keep their pinned
value. The solver throws `ConvergenceError` naming the residual if the
subproblem does not converge within its Newton budget, rather than continuing
with a bad step.

**Verification.** `tests/test_mma.cpp` solves problems with known optima: a
separable problem whose solution is `x_j = 1/n` with multiplier `2/n`, a
two-constraint problem with both constraints active at `(0.3, 0.7)` and
multipliers `(0.6, 0.8)`, and the scale-cap case; and it checks that MMA and
OC reach compliances within 5 % of each other on the same compliance-only
cantilever (they are different algorithms with different move-limit
semantics, so bit-equality is not the expectation).

**Move limit.** MMA with a moving constraint surface tolerates a smaller
move limit than OC: the stress-constrained deck uses `0.1` against the `0.2`
of the compliance decks, because at `0.2` a handful of corner elements kept
oscillating between their bounds and the run hit its iteration cap without
meeting either criterion. That behaviour is reported, not hidden, which is
how it was found.

## 5c. Aggregated stress constraint

`topology.stress` adds one constraint per load case (`topopt/StressConstraint.cpp`):

```
  sigma_e        = rho_e^q  sigma_vm( D_0 B_e(centre) u_{l,e} )      relaxed element stress
  sigma_PN       = ( sum_e  sigma_e^P )^(1/P)                        p-norm aggregate
  g_l(x)         = c_l  sigma_PN / sigma_lim  -  1  <=  0
```

with `q = 0.5` (`relaxation`), `P = 8` (`p_norm`) and `sigma_lim` the limit.
Three facts about this formulation, each of which the reader needs to hold
onto:

* **the relaxation is what makes the problem solvable.** The stress in a
  SIMP element is the solid-material stress `D_0 B u`, which does not vanish
  when the density does: a void element with a nonzero strain would violate
  any limit, and the optimiser could never remove material near a
  concentration. Multiplying by `rho^q` with `q < p` lets a vanishing element
  satisfy the constraint (the "qp relaxation", `q = 0.5` against `p = 3`).
  The price is that intermediate densities carry a stress the material would
  not: the constraint bounds the **relaxed** stress of the SIMP model, not
  the stress of a part;
* **the p-norm underestimates the maximum**, by up to `n^(1/P)` on a mesh
  of `n` elements, and the ratio changes as the design changes. The scale
  `c_l` is re-fitted every iteration to the ratio of the true relaxed maximum
  to the aggregate, blended with the previous value
  (`c_k = alpha s_max/g_PN + (1 - alpha) c_{k-1}`, `alpha = 0.5`), so the
  constraint tracks the maximum it stands for. The summary records the
  aggregate, the scale, the true maximum and the element that carries it;
* **the number that answers "does the structure meet the limit" is the
  re-solve.** Every stress-constrained run thresholds its design, analyses
  the extracted structure with full material and reports
  `interpreted_solid_analysis.max_von_mises_over_limit`. On the L-bracket
  that is 0.80 against a relaxed maximum of 0.98: the interpreted part sits
  comfortably inside a limit the relaxed model only just meets, because
  thresholding promotes the corner's intermediate densities to solid material
  and the corner stress drops. The opposite can happen on another design,
  which is why both numbers are printed.

**Gradient.** With `Psi_l = d sigma_PN / d u_l` assembled from the element
centres, one adjoint solve per load case on the cached factorisation,
`K lambda_l = Psi_l`, gives

```
  d sigma_PN / d rho_e  =  (explicit: d rho^q / d rho at fixed u)
                          -  lambda_{l,e}^T  (dE(rho_e)/d rho / E_0)  K_e^0  u_{l,e}
```

pulled back through the density filter by the exact chain rule. This is a
genuine adjoint - unlike the compliance, the stress aggregate is not
self-adjoint - and it is verified against central differences on a 2-D
problem with two load cases and a passive pad, and on a 3-D Hex8 problem
(`tests/test_mma.cpp`). The sensitivity filter is refused for this constraint,
because its filtered gradient is not the gradient of anything; so are
non-zero prescribed displacements, which the adjoint above does not include.

**What the constraint does not claim.** It bounds an aggregated, relaxed,
element-centre von Mises stress of a density model. It is not a stress
analysis of a part, not a fatigue or yield substantiation, and it says
nothing about stress concentrations the mesh does not resolve - a re-entrant
corner in a density design is mesh sensitive whether or not the constraint
is on. It is the standard research formulation (Le, Norato, Bruns, Ha and
Tortorelli 2010), implemented so the design it produces can be checked by
the re-solve.

## 6. Passive regions

A passive region is imposed as **equal lower and upper bounds** on the design
variable (`topopt/DesignDomain.cpp`):

* `solid`: `x_e = 1`, standing in for attachment bosses and bearing collars;
* `void`: `x_e = 0` (configurable), standing in for holes and system cutouts.

Pinning the *design variable* rather than the physical density keeps the
objective an exact function of the free variables, which is what makes the
finite-difference verification meaningful. The consequence, stated plainly: the
density filter blurs across a passive boundary, so the physical density adjacent
to a passive-solid/passive-void interface is intermediate rather than 0 or 1 -
the hole edge is smeared over roughly `r_min`. The extracted solid
interpretation shows what survives the threshold.

Feasibility is checked up front, with the numbers in the message:

* passive solid volume above the target -> the constraint is infeasible;
* passive solid plus all free elements at density 1 below the target -> the
  target is unreachable;
* solid and void regions claiming the same element -> the regions overlap;
* a region selecting no elements at all -> almost certainly an authoring error.

## 7. Continuation

The SIMP penalty can be ramped from `penalty_start` to `simp.penalty` in
`continuation_steps` stages of `continuation_iterations` iterations each. A low
starting penalty makes the early iterations nearly convex and reduces the
dependence of the final topology on the starting design. The penalty used at
each iteration is recorded in `history.csv`, and only the final penalty stage is
allowed to terminate the loop - otherwise a run would stop before the penalty
had been ramped.

Both aerospace decks use continuation (`p: 1.5 -> 3.0` in 3 stages of 30
iterations). The beam benchmarks do not, so the two behaviours are both
exercised.

## 8. Convergence and monotonicity

Two independent stopping indicators, either of which ends the loop:

| criterion | default | comment |
|-----------|---------|---------|
| `max \|dx\|` between iterations | `1e-2` | the classical SIMP/OC criterion |
| relative compliance change over `objective_window` iterations | `5e-5` over 20 | the backstop that terminates fine-mesh runs |

The second criterion is not a convenience. On a fine mesh the design change
stalls at a value well above `1e-2`, because thin members *migrate* one cell at
a time long after the objective has settled: a measured example on the
cantilever benchmark shows the compliance improving by only 0.35% between
iteration 250 and iteration 800 while `max |dx|` still oscillates between 0.004
and 0.12. A design-change-only rule would run to the iteration cap and report
non-convergence on a design that is, for engineering purposes, converged.

Every result records which criterion fired (`stop_reason`), the final value of
both indicators, and - when neither was met - a warning stating that the
returned design is the last iterate rather than a converged optimum.

Under MMA a third condition applies: the iterate must also be feasible
(`max_i g_i <= constraint_tolerance`). A stalled objective on a design that
violates a constraint is not convergence and is not reported as such; the run
continues to the cap and the summary then carries `feasible = false`. The
stress-constrained L-bracket stops on the objective-stall criterion after 240
iterations with a design change of 0.020 - the same fine-mesh migration
described above, now with a constraint surface that moves as the p-norm scale
is re-fitted - at a largest constraint value of `-3.8e-5`.

**Compliance is not guaranteed to fall monotonically**, and SparLab does not
pretend otherwise:

* OC is a fixed-point update with a move limit, not a line-search method, so a
  step can overshoot;
* each continuation stage raises the penalty, which raises the compliance of a
  fixed density field.

The full history is written to `history.csv` and plotted unsmoothed. In practice
the benchmark histories are monotone after the first few iterations and after
each continuation step; the test suite asserts monotonicity only over the final
third of a run, where it is a property the method actually has.

## 9. Diagnostics

The optimiser and its supporting pieces report, never hide:

| condition | what happens |
|-----------|--------------|
| positive entry in `dc/dx` | warning naming the value; the OC update clamps it |
| volume bisection cannot bracket the target | `ConvergenceError` with the reachable range |
| volume bisection bracket collapses | warning with the achieved volume and the relative gap |
| iteration cap reached | `converged = false`, `stop_reason = iteration_cap`, warning with both indicator values (and, for MMA, the largest constraint value) |
| final volume above target by more than `1e-6` relative (OC) or `constraint_tolerance` (MMA) | warning with both volumes |
| MMA subproblem fails to converge | `ConvergenceError` with the residual reached; the run stops rather than taking a bad step |
| stress constraint requested with OC, the sensitivity filter or non-zero prescribed displacements | `ConfigError` at deck load, naming the incompatibility |
| stress constraint violated at the returned design | `feasible = false`, warning with the largest constraint value and the relaxed stress ratio; a stalled objective on an infeasible design is never reported as convergence |
| interpreted structure's re-solve exceeds the stress limit | `interpreted_solid_analysis.meets_stress_limit = false` in the summary, with the ratio |
| non-finite objective | `SolverError` |
| singular or near-singular `K_ff` | `SolverError` naming the three usual modelling causes |
| density field with no element above the threshold | `MeshError` with the maximum density |
| thresholded design in several disconnected groups | reported in the interpretation block, with the discarded island volume |
| thresholded design not analysable (its largest group reaches no support, say) | warning, and `interpreted_solid_analysis.analysis_failed` with the reason. The optimisation result is still written: a design whose *interpretation* is invalid is a finding, not a reason to discard a completed run |

## 10. Measures reported

| quantity | definition | how to read it |
|----------|------------|----------------|
| compliance | `sum_l w_l f_l^T u_l` [J] | lower is stiffer |
| volume fraction | `rho_tilde^T v / V` | must equal the target |
| grey level | `M_nd = (4/n) sum_e rho_e (1 - rho_e)` | 0 = pure 0/1, 1 = every element at 0.5 |
| stiffness gain | equal-mass uniform plate compliance / optimised compliance | above 1 means the optimisation paid off |
| interpretation | threshold, retained elements, connected groups, discarded island volume | how the field was read as geometry |
| interpreted compliance | weighted compliance of the *thresholded* structure, re-solved with real material | what the design would actually deliver |
| `compliance_vs_simp_ratio` | interpreted compliance / SIMP compliance | the gap between the relaxed model and a real structure; either side of 1 (see below) |
| `constraint_violation`, `feasible` (MMA) | largest constraint value at the returned design, and whether it is within `constraint_tolerance` | a converged MMA run is feasible by definition; an iteration-capped one may not be |
| `stress.max_relaxed_stress_ratio` | true maximum of `rho^q sigma_vm` over the limit, per load case | what the constraint bounds; hovers about 1 while active |
| `stress.p_norm_ratio`, `stress.scale` | the aggregate over the limit, and the scale that maps it to the maximum | how far the aggregate sits from the maximum it stands for |
| `interpreted_solid_analysis.max_von_mises_over_limit` | peak von Mises of the re-solved thresholded structure over the limit | the number that says whether the *structure* meets the limit |
| `geometry_export` | triangles, closure, non-manifold edges, enclosed vs cell volume of `structure_after.stl` | whether the exported surface is a usable solid, and what it inherits from the interpretation |

The grey level needs care. A density filter of radius `r_min` leaves a genuinely
intermediate boundary layer roughly `r_min` wide, and on a coarse mesh that layer
is a large fraction of the domain: the MBB benchmark settles near 0.28 with the
density filter and near 0.17 with the sensitivity filter, while the *unfiltered*
run reaches 0.012 - and the unfiltered design is a checkerboard. Read the grey
level together with the topology, never alone.

`compliance_vs_simp_ratio` needs the same care, and in the opposite direction
to the obvious guess. Thresholding does two things at once: it promotes every
element above the threshold to **full** material, which the penalty had been
crediting with only `rho^p` of its stiffness, and it deletes everything below
the threshold. The first makes the interpreted structure stiffer and slightly
heavier, the second softer and lighter, and which wins depends on the design.
On the bracket the promotion dominates - the interpreted structure comes out
around 14 % stiffer than the compliance the optimiser reported, for 0.4 % more
mass - so the reported objective is *pessimistic* there rather than optimistic.
Either way it is not the compliance of a part, which is why the interpreted
value is computed and reported on every run instead of being inferred.
