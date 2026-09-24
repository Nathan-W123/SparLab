# Verification and validation

Two different questions, kept separate throughout:

**Verification** - *does the code solve the equations it claims to solve?*
Judged against exact answers for the discrete problem. These must pass to tight
tolerances; a failure is a bug.

**Validation** - *are the modelling assumptions appropriate?* Judged against an
independent theory whose assumptions differ from the model's. A finite gap is
**expected**, and the useful result is its size and sign, not agreement. Where
the gap has a known cause, the test asserts the *sign* of the discrepancy, which
is a stronger statement than asserting agreement.

Reproduce everything below with:

```bash
make test              # the Catch2 suite: 138 cases, 6 199 assertions
make verify            # the studies, which exit non-zero if any tolerance is missed
make cross-validation  # the same problems in CalculiX and scikit-fem, node by node
```

All numbers in this document come from `results/verification/summary.json`,
`results/cross_validation/summary.json` and the test run;
`docs/results/README.md` carries the machine-generated tables.

## Study results

| Study | Kind | Metric | Value | Tolerance | Result |
|-------|------|--------|-------|-----------|--------|
| Patch test, distorted mesh | verification | max relative error in `u`, strain and stress | `4.07e-15` | `1e-10` | PASS |
| Solver agreement | verification | max relative displacement difference vs dense LU | `7.04e-12` | `1e-8` | PASS |
| Topology sensitivity | verification | min over steps of the max relative gradient error | `2.18e-08` | `1e-5` | PASS |
| Mesh convergence | verification + validation | relative tip-deflection error vs Timoshenko, finest mesh, `nu = 0` | `4.78e-04` | `0.02` | PASS |
| Modal frequencies | validation | `f1` relative error vs Euler-Bernoulli, finest mesh | `5.52e-03` | `0.02` | PASS |
| Patch test 3-D, distorted Hex8 mesh | verification | max relative error in `u`, strain and stress | `1.63e-15` | `1e-10` | PASS |
| Topology sensitivity 3-D, Hex8 | verification | min over steps of the max relative gradient error | `1.56e-08` | `1e-5` | PASS |
| Mesh convergence 3-D, Hex8 cantilever | verification + validation | relative tip-deflection error vs Timoshenko, finest mesh | `5.42e-03` | `0.03` | PASS |
| Modal frequencies 3-D, Hex8 cantilever | validation | `f1` (weak axis) relative error vs Euler-Bernoulli, finest mesh | `1.06e-02` | `0.03` | PASS |

Supporting measurements from the same runs:

| Quantity | Value |
|----------|-------|
| Observed convergence order, tip deflection | `2.31` at `nu = 0`, `2.28` at `nu = 0.3` (Q4); `2.63` (Hex8) |
| Mass conservation, `sum(M)/dim` vs `rho V` | `<= 4.84e-14` relative (Q4), `<= 4.89e-14` (Hex8) |
| Axial mode vs fixed-free rod theory | `4.02e-06` relative |
| Elements excluded from the FD check at active bounds | 6 of 72 (Q4), 2 of 36 (Hex8) |

And from the cross-validation against two independent codes (section 11):

| Problem | Reference | Max relative nodal-displacement difference | Tolerance | Result |
|---------|-----------|-------------------------------------------:|----------:|--------|
| Plane cantilever, 1 440 Q4 | scikit-fem 12.0.2, `ElementQuad1` | `1.47e-10` | `1e-7` | PASS |
| Plane cantilever, 1 440 Q4 | CalculiX 2.21, `CPS4` | `8.75e-07` | `1e-5` | PASS |
| Solid block, tip load, 1 280 Hex8 | scikit-fem, `ElementHex1` | `5.59e-12` | `1e-7` | PASS |
| Solid block, tip load, 1 280 Hex8 | CalculiX, `C3D8` | `3.39e-06` | `1e-5` | PASS |
| Solid block, top pressure | scikit-fem, `ElementHex1` | `1.42e-11` | `1e-7` | PASS |
| Solid block, top pressure | CalculiX, `C3D8` | `2.43e-06` | `1e-5` | PASS |

## 1. Element-level verification

Every item is a unit test in `tests/test_element.cpp`.

| Property | How it is checked | Result |
|----------|-------------------|--------|
| `K_e` symmetry | `max\|K - K^T\| / max\|K\|` on a uniform and a distorted element | `< 1e-15` |
| Rigid-body modes | `K_e r_i` and `r_i^T K_e r_i` for the two translations and the rotation | `< 1e-14` relative |
| Null-space dimension | eigenvalues of `K_e`: exactly three below `1e-9 max\|K\|`, the rest positive | exact |
| Quadrature sufficiency | 2x2, 3x3 and 4x4 stiffness on a rectangle | agree to `1e-13` |
| Linear-field reproduction | `B u_e` for `u = a + G x` at five parametric points including the corners, on a distorted element | `< 1e-17` absolute |
| Jacobian integration | `sum_g w_g detJ` vs the shoelace area of a distorted quadrilateral | `1e-12` relative |
| Thickness / modulus scaling | `K_e(3E) = 3 K_e(E)`, `K_e(3t) = 3 K_e(t)` | `< 1e-14` |
| `M_e` symmetry and definiteness | symmetry plus the smallest eigenvalue | symmetric to `1e-15`, all eigenvalues `> 0` |
| Mass conservation | `sum(M_e) = 2 rho t A`, and `t_x^T M_e t_x = rho t A` | `1e-12` relative |
| Lumping conserves mass | row sums of `M_e` add to `sum(M_e)` | `1e-14` |
| Consistent edge load | total equals `traction x length x thickness`, split evenly on a straight edge | exact |
| Inverted / collapsed elements | clockwise and zero-height elements | `MeshError` raised |

## 2. Assembly, partitioning and solvers

`tests/test_assembly_solver.cpp`.

* **Global symmetry and rigid-body modes.** The assembled `K` is symmetric to
  `1e-14` relative and annihilates all three rigid-body vectors to `1e-13`. An
  unconstrained `K` is correctly reported as *not* positive definite.
* **Element-cache equivalence.** The same uniform mesh assembled through the
  cached path and through the generic per-element path gives matrices agreeing
  to `1e-15` relative, for both `K` and `M`. The optimisation loop therefore
  cannot silently diverge from the general case.
* **Linear scaling of the SIMP factor.** `assemble_stiffness(0.25)` equals
  `0.25 K` to `1e-15`.
* **Positive definiteness after valid constraints.** `K_ff` on a properly
  clamped plate: every LDL^T pivot positive, and the smallest eigenvalue of the
  dense equivalent strictly positive.
* **DOF round-trip.** `expand(restrict_to_free(u))` reproduces `u` at free DOFs
  and the prescribed values at constrained ones.
* **Solver agreement.** All five backends - `SimplicialLDLT`, `SimplicialLLT`,
  `SparseLU`, diagonally preconditioned CG, and dense `PartialPivLU` - on one
  6x4 plate agree within `1e-8` relative on the displacement field and `1e-9` on
  the compliance. The dedicated study on a 12x4 plate measures `7.04e-12`.

### Reaction-force equilibrium

On a 24x6 cantilever with a 1 kN tip resultant:

| Check | Measured |
|-------|----------|
| Vertical reaction vs applied load | `1e-10` relative |
| Horizontal reaction (should vanish) | `< 1e-8` N |
| Global force-balance residual | `< 1e-10` relative |
| Root moment vs `P L` | `1e-10` relative |
| Moment-balance residual | `< 1e-10` relative |
| Reactions at free DOFs (should vanish) | `0` |
| `C = 2U` | `1e-10` relative |
| Scaled solve residual | `< 1e-10` |

With a **non-zero prescribed displacement** (a `1e-4` m enforced stretch on a
`nu = 0` bar), the reaction resultant matches the closed form `E A eps` to `1e-9`
relative, the two edge reactions cancel to `1e-6` N, the compliance is exactly
zero (the applied load vector is zero) and the strain energy matches
`1/2 F delta` to `1e-9`. That last pair is the reason compliance and strain
energy are reported separately rather than assuming `C = 2U`.

## 3. Detection of ill-posed models

`ModelDiagnostics` runs before any factorisation. `tests/test_assembly_solver.cpp`
covers each case:

| Model | Detected | Reported null dimension |
|-------|----------|-------------------------|
| No constraints at all | yes | 3 |
| One node pinned in x and y | yes | 1 (rotation) |
| Bottom edge fixed in y only | yes | 1 (x translation) |
| Clamped edge | well posed | 0 |
| Pin plus a roller at another node | well posed | 0 |
| Two blocks, only one constrained | yes | the floating group is named |

`StaticAnalysis` refuses to construct on an ill-posed model, so the failure
appears before any solve. Where the analytic check cannot see the problem - an
internal mechanism, two blocks joined at a single node - it surfaces as a
non-positive Cholesky pivot, which `LinearSolver` reports as a singular system
with the three usual modelling causes named.

Explicitly exercised solver failures: a zero matrix, an indefinite matrix, a
matrix containing NaN, an empty system, and the dense solver's size refusal.
Each raises `SolverError` rather than returning a plausible-looking answer.

## 4. Patch test

The strongest single verification here. A constant-strain field

```
  u = (1e-4, -2e-4) + [[3e-4, 1e-4], [1e-4, -2e-4]] x
```

is prescribed on every boundary node of a 4x3 mesh, and the interior solution,
the recovered strain and the recovered stress are compared with the exact field.
Run at four interior-node perturbations, so the isoparametric mapping is
genuinely exercised: on a *uniform* grid a Q4 reproduces linear fields for
trivial reasons.

| Perturbation (fraction of cell size) | Max relative error in `u`, strain and stress |
|--------------------------------------|----------------------------------------------|
| 0.00 | `< 1e-11` |
| 0.15 | `< 1e-11` |
| 0.30 | `< 1e-11` |
| 0.40 | `< 1e-11` |

The dedicated study measures `4.07e-15`, i.e. round-off. This case is also why
`prescribed_displacement_only` exists on a load case: the test has no applied
force at all, and declaring that explicitly keeps an accidentally empty load
case an error.

## 5. Stress recovery against a known field

A linear displacement field is imposed directly on a distorted 3x2 mesh,
bypassing the solve, and every recovered quantity is compared with its closed
form:

| Quantity | Agreement |
|----------|-----------|
| Element strain | `< 1e-12` relative |
| Element stress | `< 1e-12` relative |
| von Mises | `1e-11` relative |
| Principal-stress sum vs `sxx + syy` | `1e-11` relative |
| Nodal averaged stress (constant field) | `< 1e-11` relative |
| Total strain energy vs `1/2 eps^T sigma V` | `1e-10` relative |
| SIMP scaling: macroscopic stress scales, solid stress and strain do not | `1e-10` |

Separately, the element field is reproduced from independent pointwise
`element_stress_at` calls to `1e-12`, confirming the documented averaging rule.

von Mises is checked against four closed forms: uniaxial tension gives the axial
stress, pure shear gives `sqrt(3) tau`, equal biaxial tension gives the same
value in plane stress, and plane strain gives the smaller value that including
`szz = nu(sxx + syy)` implies. The principal stresses are checked against both
invariants of the 2-D stress tensor.

## 6. Mesh convergence and beam theory

Cantilever: `L = 1 m`, `h = 0.1 m` (`L/h = 10`), `t = 10 mm`, `E = 70 GPa`, 1 kN
tip resultant spread over the tip edge. Six meshes from 10x2 to 320x64, at two
Poisson ratios.

References:

```
  Euler-Bernoulli:  delta = P L^3 / (3 E I)                        (bending only)
  Timoshenko:       delta = P L^3 / (3 E I) + P L / (k G A),  k = 5/6
```

**Verification part** - self-convergence against the finest mesh isolates the
discretisation error from the beam-theory modelling gap. The observed order is
`2.31` at `nu = 0` and `2.28` at `nu = 0.3`, consistent with the second-order
convergence expected of the Q4 displacement, and the error shrinks monotonically
under refinement (asserted in the test suite).

**Validation part** - on the finest mesh at `nu = 0`:

| Reference | Relative difference |
|-----------|--------------------|
| Timoshenko (bending + shear) | `4.78e-04` |
| Euler-Bernoulli (bending only) | about `6e-03` |

The FEM deflection is larger than Euler-Bernoulli in magnitude, which is the
right sign: beam theory omits the shear deformation the 2-D model includes. At
`nu = 0.3` a further gap appears, from the anticlastic (Poisson) coupling that
1-D beam theory also omits - a **modelling** difference, not a code error, which
is why the studies run both Poisson ratios side by side.

Two more properties are asserted rather than merely measured:

* compliance converges **from below** under refinement, because the
  displacement-based Q4 is a lower bound on the true compliance;
* the distributed-traction resultant is mesh independent to `1e-12` at every
  resolution, and the compliance changes by only a few percent between 8x2 and
  32x8.

## 7. Modal verification and validation

**Mass matrix** (`tests/test_modal.cpp`): symmetric to `1e-15`, positive
definite on the free DOFs, and `sum(M) = 2 rho V` to `1e-11` for both the
consistent and the lumped form, on a *distorted* mesh. The dedicated study
measures `4.84e-14`. A rigid translation carries exactly the total mass. A zero
density raises `ModelError` rather than producing a singular pair.

**Eigen solver.** On a 40x8 cantilever (698 free DOFs, above the 400-DOF
threshold, so the iterative path runs) the subspace iteration matches a dense
`GeneralizedSelfAdjointEigenSolver` to `1e-8` relative on all five requested
eigenvalues. Eigenvalues come out ascending, all eigenpair residuals are below
`1e-8`, the modes satisfy `phi^T M phi = 1` to `1e-8`, and every mode vanishes
at prescribed DOFs.

**Thickness invariance.** Scaling the thickness by 0.37 leaves all three
frequencies unchanged to `1e-9` while the mass scales by exactly 0.37. This is
the property that makes "equal-mass uniform plate" a well-defined baseline in
the topology study, and it is measured rather than assumed - the topology runs
repeat the check and record the maximum difference.

**Lumped vs consistent.** Row-sum lumping gives strictly lower frequencies (mass
moved toward the diagonal) and agrees with the consistent matrix within 5% on a
converged mesh. Both the sign and the size are asserted.

**Bending validation.** Cantilever bending frequencies at `nu = 0` against
Euler-Bernoulli:

| Mode | FEM (120x12) | Theory | Difference |
|------|--------------|--------|-----------|
| 1 | within 1% | `beta_1 L = 1.87510407` | below theory |
| 2 | within 6% | `beta_2 L = 4.69409113` | below theory |

Both are asserted to be *below* theory, since the 2-D model includes shear
deformation and rotary inertia that Euler-Bernoulli omits, and the gap must grow
with mode number. The dedicated study measures `5.52e-03` for `f1` on the finest
mesh.

**Axial validation.** A plane-stress cantilever also carries axial modes,
interleaved with the bending series. Against fixed-free rod theory
`f_n = (2n-1)/(4L) sqrt(E/rho)` - which is *exact* for a bar in uniaxial stress,
and plane stress at `nu = 0` is exactly that - the measured agreement is
`4.02e-06` relative. The mode is identified by its motion being more than 90%
axial, not by its index, since where it falls in the spectrum depends on the
aspect ratio. In the `cantilever_analysis` deck it appears as mode 3 at
1272.9 Hz against a theoretical 1272.9 Hz.

**Rigid-body screening.** An unconstrained plate either produces warnings
naming the near-zero eigenvalues, with `lambda_0 < 1e-6 lambda_3`, or raises
`SolverError`; it never returns rigid-body modes silently as physical
frequencies. Requesting more modes than there are free DOFs warns and truncates.

## 8. Topology sensitivity against finite differences

The required regression test, run both as a study
(`sparlab_verify --study sensitivity`) and as a test case
(`tests/test_topopt.cpp`).

**Setup.** A 12x6 cantilever (72 elements), density filter of radius 1.5 cells,
`p = 3`, `emin_ratio = 1e-9`, direct solver with a `1e-9` residual tolerance. The
design point is deliberately non-uniform,
`x = 0.35 + 0.30 sin(7x) cos(5y)` on free elements: a uniform field makes every
sensitivity nearly identical and would hide an indexing mistake. One passive
solid patch and one passive void disc are included so the exclusion logic is
exercised. **Every** element is tested, at eight step sizes.

**Comparison.** `dc/dx_e` from the analytical chain rule against
`[c(x + h e_e) - c(x - h e_e)] / (2h)`, where `c` is the *filtered* SIMP
compliance - i.e. the whole objective, filter included.

| Step `h` | Max relative error | Behaviour |
|----------|--------------------|-----------|
| `1e-1` ... `1e-3` | falls as `h^2` | truncation-error dominated |
| `1e-5` | `2.18e-08` (best) | balance point |
| `1e-7`, `1e-8` | rises again | round-off dominated: differencing two compliances that agree to ~10 digits |

The `h^2` fall and the subsequent rise is the signature of a correct gradient;
a wrong gradient shows a plateau at its own error level instead. Three error
measures are reported at every step: the maximum and RMS over elements, and the
relative error of the directional derivative along the gradient, which is the
most sensitive single scalar.

**Exclusions.** 6 of 72 elements are excluded, each with a recorded reason:

* passive variables (lower bound equals upper bound) - no admissible
  perturbation exists;
* free variables within `h` of 0 or 1 - a central difference would leave the
  feasible box, so the one-sided value is not a valid central difference.

**Pass criterion.** Maximum relative error `<= 1e-5` at the best step. Measured
`2.18e-08`, i.e. three orders of magnitude of margin. The test-suite version
additionally requires the RMS error below `1e-6`, the directional error below
`1e-5`, and asserts that passive elements *were* excluded - so a regression that
quietly stopped excluding them would fail.

The check is also run **without** a filter, where `dc/dx` must equal
`dc/drho` exactly (verified to `1e-14`), and on a **two-load-case** problem with
unequal weights, where the weighted objective's gradient is verified the same
way.

## 9. Topology optimiser verification

`tests/test_topopt.cpp`.

| Property | Check | Result |
|----------|-------|--------|
| SIMP derivative | central differences at five densities | `1e-6` relative |
| SIMP clamping | `rho` outside `[0,1]` | clamped, no NaN from `pow` |
| Mass interpolation | `E/m` at `rho = 1e-3` under the matched law | order 1, and above the linear law's value |
| Filter row sums | every row of `Hhat` | `1` to `1e-12` |
| Filter preserves a constant | `Hhat c = c` | `1e-14` |
| Filter adjoint | `<H x, v> == <x, H^T v>` | `1e-12` |
| Filter range | a 0/1 checkerboard stays in `[0,1]` and its range shrinks | yes |
| Volume constraint per OC step | achieved vs target volume over five steps | `1e-8` relative, `volume_converged` true |
| Bounds and move limit | every variable after each step | respected to `1e-12` |
| Passive regions | solid pinned at 1, void at 0, after a full run | exact |
| Gradient signs | `dc/dx <= 0`, `dv/dx > 0` | yes |
| Adding material lowers compliance | `c(x + 0.1) < c(x)` | yes |
| Strain-energy identity | `sum_e U_e = c/2` for one load case | `1e-10` |
| Multi-load weighting | `c = sum_l w_l c_l` with normalised weights | `1e-12` |
| Continuation | penalty non-decreasing, ends at the target | exact |
| Non-convergence reporting | 3-iteration cap with an unreachable tolerance | `converged = false`, warning names the iteration cap |
| Filter radius sets the length scale | total variation of the density field at 1.2 vs 3.0 cells | larger radius gives a smoother design |

**End-to-end** on a 30x15 cantilever at 40% volume: converges in 123
iterations, the volume fraction equals the target to `1e-6` at *every* recorded
iteration, compliance improves by a factor of 5.9 over the uniform start, the
grey level settles below 0.40, and the compliance history is non-increasing over
the final third of the run to within `1e-3` relative.

**KKT consistency** (reported in `docs/benchmarks.md` rather than asserted): at
the converged MBB design, `B_e = -(dc/dx_e)/(dv/dx_e)` across elements strictly
inside their bounds has a 5-95 percentile spread of about 9%, which is the
residual expected of OC with a move limit and a finite change tolerance.

## 10. Configuration and I/O

`tests/test_io.cpp` - 30+ cases. Worth naming:

* the JSON parser handles the full value grammar, comments, trailing commas and
  `\u` escapes, and reports the **line and column** of a syntax error (checked
  on a multi-line document);
* duplicate keys, unterminated strings and comments, trailing content and
  malformed numbers are all errors;
* serialisation round-trips, including a non-finite number, which becomes
  `null` rather than invalid JSON;
* `ConfigNode` type errors name the dotted path (`a.c must be a number, got
  string`);
* **unknown keys are reported.** A deck with `"thicknes"` instead of
  `"thickness"` is an error under `--strict-config` and a warning otherwise -
  and the test asserts the misspelling really did cause the default to be taken,
  which is exactly the failure mode the check exists to catch;
* a region selecting no nodes, a region naming two primitives at once, a
  negative load-case weight, modal analysis without a density, and an empty
  load case are all rejected with specific messages;
* the CSV writer enforces its declared column count; the VTK writer's header,
  cell block and data blocks are checked by parsing the file back.

## 11. The solid (Hex8) path

Every study above has a 3-D counterpart, run by `sparlab_verify --study
patch-test-3d | sensitivity-3d | mesh-convergence-3d | modal-3d` and covered
by `tests/test_solid.cpp` (27 cases).

**Element level.** Partition of unity and the nodal delta property of the
trilinear shape functions; the cube Gauss rules as tensor products with total
weight 8; `K_e` symmetric with exactly six zero eigenvalues (three
translations, three rotations `e_j x r`) on a unit cube and on a distorted
cell; the 2x2x2 rule exact on a box (agrees with 3x3x3 to round-off); the
strain operator reproducing a linear field with six independent constant
strains on a distorted cell; `sum_g w_g detJ` against the closed-form volume
of a box and of a sheared box; inverted, degenerate and mis-sized cells
rejected; the consistent mass symmetric, positive definite and summing to
`3 rho V`; a constant face traction split evenly over the four corners with
the exact resultant.

**Patch test 3-D.** `u = (1e-4, -2e-4, 0.5e-4) + G x` with a symmetric `G`
carrying all six constant strain components, prescribed on the boundary of a
3x3x3 mesh whose interior nodes are perturbed; the interior displacement,
strain and stress agree with the exact field to `1.63e-15` - round-off, as in
2-D.

**Diagnostics 3-D.** An unconstrained solid reports a null dimension of six
(and the assembled `K` annihilates all six rigid-body vectors to `1e-13`), a
single fixed node leaves three rotations, a face fixed in one component
leaves the in-face translations and the rotation about the normal, a clamped
face is well posed, and so is a 3-2-1 set of three non-collinear point
supports; a 2-D model rejects any out-of-plane input with a message naming
the key.

**Mesh convergence 3-D.** The solid cantilever `L = 1 m, h = 0.1 m, b = 0.05 m`
at `nu = 0` on four meshes from 16x2x1 to 96x12x6 (306 to 26 481 DOFs), tip
resultant 1 kN over the tip face. The fully integrated Hex8 is stiff in
bending on coarse meshes - the coarsest is 16 % below Timoshenko - and
approaches the reference from below as the section is refined: the finest
mesh is within `5.42e-03` of Timoshenko, and the self-convergence order
against it is `2.63`. The deflection approaching from *below* is the sign
shear locking must have, and the test suite asserts it.

**Modal 3-D.** The same bar at `rho = 2700 kg/m^3` has two first bending
modes, about the weak (0.05 m) and the strong (0.1 m) axis. Both are compared
with Euler-Bernoulli on three meshes: the strong-axis mode is within 2.1 %
on the coarsest mesh, the weak-axis one needs the finest mesh to reach
`1.06e-02` because it bends through only a few elements. Mass is conserved to
`4.89e-14` (`sum(M) = 3 rho V`).

**Sensitivity 3-D.** The compliance gradient of a 6x3x2 Hex8 cantilever with
a density filter, one passive solid and one passive void cell, at a
non-uniform design point, against central differences at eight steps: best
step `1e-5`, maximum relative error `1.56e-08`, the same `h^2` fall and
round-off rise as in 2-D.

**Solver, faces and I/O.** Face tractions give a mesh-independent resultant
and the reactions balance it to `1e-10`; 3-D selectors (box with z limits,
sphere, cylinder about any axis), sub-meshes and face connectivity are
exercised; a `structured_hex` deck parses into a 3-D model; dimension
mismatches in a deck are rejected with a reason; the writers produce z
columns, six stress components and VTK hexahedra.

**End to end.** A short 3-D optimisation (OC) runs, converges on its
criteria, thresholds to a single connected structure and re-solves it.

## 12. MMA and the stress constraint

`tests/test_mma.cpp` (12 cases):

| Property | Check | Result |
|----------|-------|--------|
| Separable problem with a known optimum | `min sum x_j^2 s.t. sum x_j >= 1` converges to `x_j = 1/n`, multiplier `2/n` | `1e-6` |
| Two active constraints | optimum `(0.3, 0.7)` with multipliers `(0.6, 0.8)` | `1e-5` |
| Badly scaled constraint | a constraint violated by `1e5` is scaled to the cap and its multiplier scaled back | subproblem converges; `lambda ~ 1e-5` |
| Malformed input | wrong sizes, non-finite entries, bounds crossed | rejected |
| von Mises as a quadratic form | `sigma^T V sigma` against `von_mises()` in plane stress, plane strain and 3-D | `1e-12` |
| p-norm aggregate | bounded by the maximum from below and by `n^(1/P)` times it from above; centre stress on a uniform mesh equals the recovered element average; the adaptive scale maps it to the maximum | exact / `1e-12` |
| Stress-constraint gradient, 2-D | adjoint gradient of the aggregate, two load cases with unequal weights and a passive pad, vs central differences | relative error `< 1e-5` per entry (scale-aware) |
| Stress-constraint gradient, 3-D | the same on a Hex8 cantilever | `< 1e-5` |
| No filter | `dg/dx == dg/drho` | exact |
| Sensitivity filter | refused with a `ConfigError` | yes |
| MMA vs OC | the same 24x12 compliance problem | compliances within 5 % |
| Stress-limited L-bracket | 32x32, limit at 70 % of the unconstrained design's relaxed peak, which the test locates at the re-entrant corner and not at the load pad | feasible; relaxed ratio within 5 % of 1; the peak at least a fifth below the unconstrained one; compliance higher |
| Infeasible limit | limit `1e3` Pa, far below anything reachable | `feasible = false`, warning; the run is reported, not hidden |
| Deck parsing | every `optimizer.mma.*` and `topology.stress.*` key, plus the two incompatibilities | round-trip and rejection |

## 13. Geometry export

`tests/test_geometry.cpp`: a plane mesh extrudes to a closed surface of
volume `A t` and area `2A + P t` to `1e-12`; a hex mesh's boundary is closed,
outward (checked against the centroid) and encloses exactly its volume; a
thresholded L of cells on a distorted mesh is watertight and agrees with its
cell volume to the flat-triangle error of its bilinear cut faces (`1.5e-4`,
and *not* to round-off - the test asserts both); on the uniform grid the same
cut is exact to `1e-12`; dropping one triangle is reported as three unmatched
edges; two cells touching along an edge are reported as closed with one
non-manifold edge; a binary STL round-trips through `read_stl` with the
standard 84-byte header layout.

## 14. Cross-validation against independent codes

The one item earlier versions of this document listed as absent. The same
discrete problem - nodes, connectivity, supports, consistent nodal loads - is
solved by two independent codes and the nodal displacements compared node by
node (`python/scripts/cross_validate.py`, `make cross-validation`):

* **scikit-fem 12.0.2** rebuilds the problem with `ElementQuad1` /
  `ElementHex1`, the same Lame constants (`lambda* = 2 lambda G/(lambda+2G)`
  for plane stress) and the same integration order. It is the *same element
  formulation* in an independent implementation, so the only expected
  difference is linear-solver round-off - and that is what is measured:
  `1.5e-10` (Q4) and `5.6e-12`, `1.4e-11` (Hex8) relative;
* **CalculiX 2.21** (`ccx`) runs the exported `.inp` decks. `C3D8` is the same
  trilinear element as SparLab's Hex8; the measured differences of `3.4e-06`
  and `2.4e-06` are within the six-significant-digit rounding of its `.frd`
  result file (floor `5e-6`), i.e. as close as the file format allows one to
  see. `CPS4` is a plane element CalculiX expands through the thickness - a
  different discretisation of the plane problem - and it agrees to
  `8.7e-07`.

Tolerances are `1e-7` for scikit-fem and `1e-5` for CalculiX, both recorded
in the summary with the `.frd` floor. The comparison exits non-zero if any
pair exceeds its tolerance; CI runs the scikit-fem half on every push.

## What is not covered

Stated plainly, since the absence matters as much as the presence:

* **no comparison against experiment**;
* the cross-validation covers linear static displacements on two problems.
  Stresses, natural frequencies and the optimised designs are not compared
  with another code;
* plane strain is unit-tested but no verification *study* runs in it;
* the sensitivity checks run on 72- and 36-element meshes (they need two
  extra solves per element per step); the gradients are not FD-verified at
  benchmark resolution, although they are the same code path;
* the stress-constraint gradient is FD-verified; the *constraint's* claim
  - that the relaxed aggregate bounds the stress of a part - is checked only
  by the re-solve of the thresholded structure, which is a consistency check,
  not a proof;
* the Hex8 mesh-convergence study reaches 26 481 DOFs; the direct solver's
  `O(n^2.3)` factorisation cost is what stopped it there;
* no convergence study of the *optimised topology* against mesh size in the
  verification suite - that lives in the design study, where the
  `mesh_fixed_r` and `mesh_fixed_cells` arms address it directly.
