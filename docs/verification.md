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
| Solver agreement (7 solvers, multigrid and `auto` included) | verification | max relative displacement difference vs dense LU | `8.97e-12` | `1e-8` | PASS |
| Topology sensitivity | verification | min over steps of the max relative gradient error | `2.18e-08` | `1e-5` | PASS |
| Mesh convergence | verification + validation | relative tip-deflection error vs Timoshenko, finest mesh, `nu = 0` | `4.78e-04` | `0.02` | PASS |
| Modal frequencies | validation | `f1` relative error vs Euler-Bernoulli, finest mesh | `5.52e-03` | `0.02` | PASS |
| Patch test 3-D, distorted Hex8 mesh | verification | max relative error in `u`, strain and stress | `1.63e-15` | `1e-10` | PASS |
| Topology sensitivity 3-D, Hex8 | verification | min over steps of the max relative gradient error | `1.56e-08` | `1e-5` | PASS |
| Mesh convergence 3-D, Hex8 cantilever | verification + validation | relative tip-deflection error vs Timoshenko, finest mesh | `5.42e-03` | `0.03` | PASS |
| Modal frequencies 3-D, Hex8 cantilever | validation | `f1` (weak axis) relative error vs Euler-Bernoulli, finest mesh | `1.06e-02` | `0.03` | PASS |
| Patch test, distorted Tri3 and Tet4 meshes | verification | max relative error in `u`, strain and stress | `9.17e-15` | `1e-10` | PASS |
| Mesh convergence, Tri3 and Tet4 cantilevers | verification + validation | larger of the two finest-mesh errors vs Timoshenko | `1.66e-02` | `0.03` | PASS |
| Multigrid CG vs Cholesky, Hex8 and Tet4 | verification | max relative displacement difference vs LDL^T; iteration growth under refinement at most `1.6x` | `4.37e-12` | `1e-8` | PASS |
| Sensitivity through the Heaviside projection, Q4 and Tet4 | verification | worst over `beta = 2, 8, 32` of the best scaled-entry or directional error | `9.65e-08` | `1e-5` | PASS |

Supporting measurements from the same runs:

| Quantity | Value |
|----------|-------|
| Observed convergence order, tip deflection | `2.31` at `nu = 0`, `2.28` at `nu = 0.3` (Q4); `2.63` (Hex8); `2.27` (Tri3); `3.22` (Tet4, still pre-asymptotic - section 15) |
| Multigrid CG iterations, coarsest to finest multi-level mesh | `14 -> 16` (Hex8), `16 -> 20` (Tet4); Jacobi CG `238 -> 699` and `481 -> 729` on the same meshes |
| Mass conservation, `sum(M)/dim` vs `rho V` | `<= 4.84e-14` relative (Q4), `<= 4.89e-14` (Hex8) |
| Axial mode vs fixed-free rod theory | `4.02e-06` relative |
| Elements excluded from the FD check at active bounds | 6 of 72 (Q4), 2 of 36 (Hex8) |

And from the cross-validation against two independent codes (section 14):

| Problem | Reference | Max relative nodal-displacement difference | Tolerance | Result |
|---------|-----------|-------------------------------------------:|----------:|--------|
| Plane cantilever, 1 440 Q4 | scikit-fem 12.0.2, `ElementQuad1` | `1.47e-10` | `1e-7` | PASS |
| Plane cantilever, 1 440 Q4 | CalculiX 2.21, `CPS4` | `8.75e-07` | `1e-5` | PASS |
| Solid block, tip load, 1 280 Hex8 | scikit-fem, `ElementHex1` | `5.59e-12` | `1e-7` | PASS |
| Solid block, tip load, 1 280 Hex8 | CalculiX, `C3D8` | `3.39e-06` | `1e-5` | PASS |
| Solid block, top pressure | scikit-fem, `ElementHex1` | `1.42e-11` | `1e-7` | PASS |
| Solid block, top pressure | CalculiX, `C3D8` | `2.43e-06` | `1e-5` | PASS |
| Plane cantilever, 2 880 Tri3, `nu = 0` | scikit-fem, `ElementTriP1` | `1.66e-11` | `1e-7` | PASS |
| Plane cantilever, 2 880 Tri3, `nu = 0` | CalculiX, `CPS3` | `8.86e-07` | `1e-5` | PASS |
| Solid block, 7 680 Tet4, tip load / top pressure | scikit-fem, `ElementTetP1` | `1.72e-12` / `9.46e-12` | `1e-7` | PASS |
| Solid block, 7 680 Tet4, tip load / top pressure | CalculiX, `C3D4` | `3.64e-06` / `2.97e-06` | `1e-5` | PASS |
| Gmsh lug bracket, 20 336 Tri3, `nu = 0`, two load cases | scikit-fem, `ElementTriP1` | `6.09e-13` / `7.52e-14` | `1e-7` | PASS |
| Gmsh lug bracket, 20 336 Tri3, `nu = 0`, two load cases | CalculiX, `CPS3` | `2.76e-06` / `4.29e-06` | `1e-5` | PASS |
| Gmsh lug bracket, `nu = 0.33` (the benchmark material) | scikit-fem, `ElementTriP1` | `5.04e-13` / `2.42e-13` | `1e-7` | PASS |
| Gmsh lug bracket, `nu = 0.33` | CalculiX, `CPS3` | `5.63e-04` / `1.14e-03` | - | INFO |
| Gmsh engine mount, 39 936 Tet4, two load cases | scikit-fem, `ElementTetP1` | `1.48e-12` / `5.28e-13` | `1e-7` | PASS |
| Gmsh engine mount, 39 936 Tet4, two load cases | CalculiX, `C3D4` | `1.71e-06` / `1.98e-06` | `1e-5` | PASS |

The two `INFO` rows are not a disagreement between codes but between
idealisations: CalculiX expands its plane elements into a layer of solid
elements, which reproduces plane stress only at `nu = 0` - the same mesh at
`nu = 0` agrees to the `.frd` rounding floor (section 14).

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
* **Solver agreement.** The direct backends - `SimplicialLDLT`,
  `SimplicialLLT`, `SparseLU`, dense `PartialPivLU` - and the iterative ones -
  diagonally preconditioned CG, multigrid CG, and `auto` - agree on one 6x4
  plate within `1e-8` relative on the displacement field and `1e-9` on the
  compliance. The dedicated study on a 12x4 plate (with the multigrid coarse
  size lowered to 16 unknowns, so it builds a real hierarchy on 130 DOFs)
  measures at most `8.97e-12`, from the multigrid solver in 41 iterations.
* **Cached-pattern assembly.** The assembly that scatters into the recorded
  sparsity pattern is bitwise identical to the triplet assembly for `K` and
  `M`, on Q4, Tri3, Hex8 and Tet4 meshes, with and without SIMP scale
  factors (`tests/test_multigrid.cpp`).

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

* **scikit-fem 12.0.2** rebuilds the problem with `ElementQuad1`,
  `ElementTriP1`, `ElementHex1` or `ElementTetP1`, the same Lame constants
  (`lambda* = 2 lambda G/(lambda+2G)` for plane stress) and the same
  integration order. It is the *same element formulation* in an independent
  implementation, so the only expected difference is linear-solver round-off
  - and that is what is measured, from `7.5e-14` to `1.5e-10` relative over
  the seven problems;
* **CalculiX 2.21** (`ccx`) runs the exported `.inp` decks. `C3D8` and `C3D4`
  are the same trilinear hexahedron and linear tetrahedron as SparLab's; the
  differences, `1.7e-06` to `3.6e-06`, are within the six-significant-digit
  rounding of its `.frd` result file (floor `5e-6`), i.e. as close as the
  file format allows one to see - including the 39 936-tetrahedron engine
  mount read from a Gmsh file, which SparLab solves with the multigrid
  solver.
* **CalculiX's plane elements are a different idealisation.** `CPS4` and
  `CPS3` are expanded internally into a layer of solid elements with the
  plane-stress condition imposed on it. That reproduces plane stress only
  for `nu = 0`, and the measurement says so: the Gmsh lug bracket at `nu = 0`
  agrees to `2.8e-06` and `4.3e-06`, within the `.frd` floor, while the same
  mesh at its real `nu = 0.33` differs by `5.6e-04` and `1.1e-03`. The plane
  decks judged against CalculiX therefore run at `nu = 0` (the cantilevers
  and `lug_bracket_nu0_analysis`); the `nu = 0.33` comparison is recorded as
  informational (`passed: null`, `INFO` in the tables), and scikit-fem - the
  same plane element - is the verification for it.

Tolerances are `1e-7` for scikit-fem and `1e-5` for CalculiX, both recorded
in the summary with the `.frd` floor. The comparison exits non-zero if any
judged pair exceeds its tolerance; CI runs the scikit-fem half on all seven
problems on every push.

## 15. The linear simplices (Tri3, Tet4)

`tests/test_unstructured.cpp` and `sparlab_verify --study patch-test-simplex |
mesh-convergence-simplex`.

**Element level.** `K_e` symmetric, exactly three (Tri3) or six (Tet4) zero
eigenvalues and the rest positive, a linear field reproduced exactly, the
closed-form consistent mass summing to `dim rho V`, a constant edge or face
traction split equally with the exact resultant, and inverted or flat cells
rejected with a `MeshError`.

**Patch test.** The constant-strain fields of sections 4 and 11 on triangle and tetrahedron
meshes split from perturbed Q4 and Hex8 grids (interior nodes moved by up to
40 % and 30 % of the cell size, and undistorted for comparison):
interior displacement, strain and stress agree with the exact field to
`9.17e-15` - round-off, as for the bilinear and trilinear elements.

**Mesh convergence.** The same two cantilevers as sections 6 and 11, meshed
by splitting the structured cells, so every simplex mesh has exactly the
nodes of a Q4 or Hex8 mesh:

| Element | DOFs | Error vs Timoshenko, coarsest to finest | Q4 / Hex8 on the same beam |
|---------|------|-----------------------------------------|-----------------------------|
| Tri3 (`nu = 0.3`) | 126 to 21 186 | `37.7 %, 13.4 %, 3.87 %, 1.14 %, 0.42 %` | `0.85 %` at 2 754 DOFs, `0.22 %` at 41 730 |
| Tet4 (`nu = 0`) | 306 to 59 211 | `50.2 %, 20.7 %, 6.25 %, 2.90 %, 1.66 %` | `1.21 %` at 8 775 DOFs, `0.54 %` at 26 481 |

The constant-strain elements lock in bending exactly as expected: they
approach the beam value from below, and on the *same nodes* the Tet4 error is
about five times the Hex8 error (`6.25 %` against `1.21 %` at 8 775 DOFs,
`2.90 %` against `0.54 %` at 26 481). The Tri3 needs roughly three times the
unknowns of the Q4 for the same error (interpolating the Q4 curve: 2.6x to
2.8x over the measured range). The order of the error against
Timoshenko rises towards 2 under refinement - `1.28, 1.72, 1.89, 1.95` for
the Tet4 - and the Tri3's last step (`1.44`) is already limited by the
modelling gap between the beam and plane elasticity, the same gap the Q4
curve levels off at. The self-convergence order the study also records
(`2.27` Tri3, `3.22` Tet4) is measured against the finest mesh, which for
the Tet4 is itself still `1.7 %` from the beam value, so it overstates the
rate; the error against the independent reference is the number to read.
The study passes at `1.66e-02` against a tolerance of `0.03`.

## 16. Meshes read from files

`tests/test_unstructured.cpp` (the reader cases) and the two Gmsh parts:

* Gmsh 2.2 and 4.1 files, written by the test itself, round-trip triangle
  and tetrahedron meshes node for node, with their physical groups as node
  and element sets;
* a clockwise triangle and a mirrored tetrahedron are re-ordered and
  counted; unused nodes are dropped and counted; a mesh in millimetres read
  with `scale = 0.001` has the right extent; coincident nodes are reported,
  and merged on request;
* malformed and unsupported files fail with a message naming the cause and,
  where there is one, the Gmsh option that fixes it: mixed triangles and
  quadrilaterals (`Mesh.RecombineAll`), second-order cells
  (`Mesh.ElementOrder = 1`), a binary file or an unsupported version, a
  plane mesh off `z = 0`, a cell referencing a missing node, a file holding
  only boundary elements (`Mesh.SaveAll`), a truncated section, a folded
  cell, a missing file and an unknown extension;
* SparLab's own CalculiX decks read back node for node through the Abaqus
  reader, which also handles `*INCLUDE`, `*NSET` / `*ELSET` with `GENERATE`,
  continuation lines, and refuses `*INSTANCE` translations and multi-part
  assemblies it cannot place;
* a deck addresses a group by name for supports, loads and passive regions,
  resolves a relative mesh path against its own directory, and a group that
  does not exist, or a node set used as an element region, is an error
  listing the sets the mesh has.

The two committed parts (`configs/meshes/`) are then solved and
cross-validated like the structured meshes (section 14): the 20 336-triangle
lug bracket and the 39 936-tetrahedron engine mount agree with scikit-fem to
`1e-12` and better, and with CalculiX to its `.frd` rounding wherever the two
codes solve the same problem.

## 17. The multigrid solver

`tests/test_multigrid.cpp` and `sparlab_verify --study multigrid`.

**The hierarchy's defining identities**, on a distorted Q4 and a Tet4
mesh: a near-null space of dimension 3 and 6; the tentative prolongator
reproduces it exactly (`P_hat B_c = B` to `1e-12` relative) with
orthonormal columns; every coarse operator equals `P^T A P`; the levels
shrink; the operator complexity lies between 1 and 3. On a distorted Hex8
mesh the V-cycle is symmetric (`<M^-1 x, y> = <x, M^-1 y>` to `1e-11`) and
positive definite, and reduces the error, for both smoothers - which is
what plain CG needs of a preconditioner.

**Agreement with the direct solver** on every element type, including a
SIMP-contrast model (densities down to the `1e-9` floor): the displacement
differs from LDL^T by less than the requested tolerance allows. A reused
hierarchy - aggregates kept, numbers refreshed - is bitwise equal to one
built from scratch for the new matrix. Warm starts reach the same answer in
fewer iterations (none, from the exact solution), and one thread and many
give bitwise identical results. An unsupported model is reported as
under-constrained by the coarse factorisation, and malformed input is
rejected.

**Iterations under refinement.** The study solves a cantilever block
(`nx x nx/2 x nx/4` cells, tip load) to a relative residual of `1e-10` with
multigrid CG, Jacobi CG and LDL^T:

| Element | DOFs | Levels | MG iterations | Jacobi iterations | Difference MG vs LDL^T |
|---------|------|-------:|--------------:|------------------:|-----------------------:|
| Hex8 | 2 295 | 2 | 14 | 238 | `7.2e-13` |
| Hex8 | 15 147 | 2 | 14 | 469 | `2.2e-12` |
| Hex8 | 47 775 | 3 | 16 | 699 | `3.7e-12` |
| Tet4 | 6 825 | 2 | 16 | 481 | `2.8e-12` |
| Tet4 | 21 090 | 3 | 20 | 729 | `4.4e-12` |

Jacobi's count grows with `1/h`, as `sqrt(kappa)` must; the multigrid
count barely moves, which is the property the solver exists for. The
smallest meshes of the study (up to 1 092 DOFs) sit below the coarse-grid
size, where the "hierarchy" is the direct coarse solve and converges in one
iteration; the growth limit is judged over the multi-level meshes only
(worst `1.25x`, against a limit of `1.6x`). Timings are in
`results/verification/multigrid_scaling.csv` and the larger comparison in
`docs/benchmarks.md`: for one solve at these sizes Jacobi CG is about as fast
as multigrid, because its cheap iterations cost about what the multigrid
setup does; the direct solver is 45 times slower at 47 775 Hex8 DOFs
(40.9 s against 0.90 s).

## 18. The Heaviside projection

`tests/test_projection.cpp` and `sparlab_verify --study
sensitivity-projection`.

* The map keeps `rho_bar(0) = 0` and `rho_bar(1) = 1` for every `beta`
  from 0.5 to 64, is monotone, is within `1e-6` of the identity at
  `beta = 1e-3` and a step at `beta = 256`, and its derivative matches a
  central difference of the map itself.
* The compliance and volume gradients through filter and projection match
  central differences on a distorted Q4, a Hex8 and a Tet4 mesh, and the
  stress-constraint gradient on a Q4 and a Hex8 mesh.
* The `beta` schedule steps on its interval, early on convergence when asked,
  caps at `beta_max`, and a deck switches the projection on with its
  schedule; the sensitivity filter is refused.
* End to end, a projected run finishes at `beta_max` with a much lower grey
  level than the same run unprojected.

The study checks the gradient at `beta = 2, 8, 32` on a 12x6 Q4 and a
5x2x2-cell Tet4 mesh (72 and 120 design variables), at three steps each:

| Element | `beta = 2` | `beta = 8` | `beta = 32` |
|---------|-----------:|-----------:|------------:|
| Q4 | `6.3e-08` | `8.4e-09` | `9.7e-08` |
| Tet4 | `9.4e-09` | `1.2e-08` | `5.0e-08` |

Each entry is the best step's worse of the scaled entry error and the
directional-derivative error. At `beta = 32` the plain relative error of the
worst entry is between `1.37` and `1.76` on the Q4 mesh, depending on the
step: far from the threshold the projection's derivative is tiny, so the
analytical entry is near zero and its central difference is round-off.
Judging each entry against
`max(|analytical|, |FD|, 1e-3 ||gradient||_inf)` separates that from a real
error, and the directional derivative along the full gradient, which no
entry can hide in, agrees to `1e-9`.

## What is not covered

Stated plainly, since the absence matters as much as the presence:

* **no comparison against experiment**;
* the cross-validation covers linear static displacements on seven
  problems, two of them read from mesh files. Stresses, natural frequencies
  and the optimised designs are not compared with another code;
* plane strain is unit-tested but no verification *study* runs in it;
* the sensitivity checks run on 72- and 36-element meshes (they need two
  extra solves per element per step); the gradients are not FD-verified at
  benchmark resolution, although they are the same code path;
* the stress-constraint gradient is FD-verified; the *constraint's* claim
  - that the relaxed aggregate bounds the stress of a part - is checked only
  by the re-solve of the thresholded structure, which is a consistency check,
  not a proof;
* the Hex8 mesh-convergence study still stops at 26 481 DOFs and runs on
  the direct solver, so its numbers stay comparable with earlier runs; the
  Tet4 study reaches 59 211 DOFs, where `auto` has switched to multigrid;
* the multigrid solver is verified against the direct solver up to 47 775
  DOFs in the study and at the benchmark sizes in `docs/benchmarks.md`;
  beyond the direct solver's reach (the 356 475- and 830 115-DOF runs) there
  is no second solution to compare with, only the residual checks;
* the mesh readers are tested on files written by the tests and by Gmsh
  4.15.2; files from other generators (Abaqus/CAE, HyperMesh, Salome) use the
  same keywords but have not been tried;
* no convergence study of the *optimised topology* against mesh size in the
  verification suite - that lives in the design study, where the
  `mesh_fixed_r` and `mesh_fixed_cells` arms address it directly.
