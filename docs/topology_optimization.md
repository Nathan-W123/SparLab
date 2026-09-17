# Topology optimisation

## 1. The problem

Minimum compliance under a volume constraint, with per-element density design
variables:

```
  min over x      c(x) = sum_l  w_l  f_l^T u_l
  subject to      K(rho_tilde(x)) u_l = f_l          for every load case l
                  g(x) = rho_tilde^T v - nu V  <=  0
                  0 <= x_e <= 1                      (tightened on passive elements)
```

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
sensitivity` and regression-tested in `tests/test_topopt.cpp`.

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

**Why OC and not MMA.** One constraint, a separable objective in the SIMP sense,
and convergence in a few hundred iterations with no tuning. Its limitation -
exactly one inequality constraint - is stated in `docs/limitations.md`; a stress
or frequency constraint would need MMA or an augmented-Lagrangian method.

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
| iteration cap reached | `converged = false`, `stop_reason = iteration_cap`, warning with both indicator values |
| final volume above target by more than `1e-6` relative | warning with both volumes |
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
