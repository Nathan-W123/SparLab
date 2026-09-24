# Architecture

## Layering

Each box is one directory. An arrow means "depends on"; there are no cycles and
no upward dependencies, so any layer can be replaced without touching the ones
below it.

```
                          +-------------------------+
                          |  apps/                  |
                          |  sparlab_solve          |   command-line drivers
                          |  sparlab_topopt         |   argument parsing,
                          |  sparlab_verify         |   orchestration, console
                          |  sparlab_bench          |   reports
                          +------------+------------+
                                       |
        +------------------------------+------------------------------+
        |                              |                              |
        v                              v                              v
+---------------+          +-------------------------+       +-----------------+
|  io/          |          |  topopt/                |       |  fem/           |
|  Json         |          |  DesignDomain           |       |  ModalAnalysis  |
|  Config       |--------->|  SimpInterpolation      |------>|  StaticAnalysis |
|  CsvWriter    |          |  DensityFilter          |       |  StressRecovery |
|  VtkWriter    |          |  Sensitivity            |       |  ModelDiagnostics|
|  StlWriter    |          |  StressConstraint       |       |  LinearSolver   |
|  CalculixWriter|         |  OptimalityCriteria     |       |  Assembler      |
|  ResultWriter |          |  Mma                    |       |  FemModel       |
+---------------+          |  TopologyOptimizer      |       |  BoundaryConds  |
                           +-------------------------+       |  DofManager     |
                                                             |  Selector       |
                                                             +--------+--------+
                                                                      |
                                    +---------------------------------+
                                    |                 |               |
                                    v                 v               v
                           +----------------+  +-------------+  +------------+
                           |  elements/     |  |  material/  |  |  mesh/     |
                           |  Element (abc) |  |  Isotropic  |  |  Mesh      |
                           |  Quad4  Hex8   |  |  Material   |  |  Structured|
                           |  Quadrature    |  |             |  |  SubMesh   |
                           +-------+--------+  +------+------+  +-----+------+
                                   |                  |               |
                                   +--------+---------+---------------+
                                            v
                                   +--------------------+
                                   |  core/             |
                                   |  Types  Exceptions |
                                   |  Logging  Timer    |
                                   |  Version           |
                                   +--------------------+

   python/sparlab_viz/  reads only the files io/ writes:
       loaders -> style -> fields / solid -> plots / plots3d / studies
   python/scripts/cross_validate.py  drives CalculiX and scikit-fem on the
       exported decks and compares nodal displacements
```

The core is dimension-generic at run time: `Mesh::dim()` is 2 or 3, the
element kernels are written for their own dimension (`Quad4` on a 2-D,
`Hex8` on a 3-D mesh) and everything above the element layer - DOF
numbering, assembly, partitioning, stress recovery, modal analysis, the
optimiser, the writers - is written against `dim`, `nodes_per_element` and
the Voigt length rather than against a fixed 2. The Q4 kernels keep
fixed-size internal matrices, so the plane path produces bit-for-bit what it
did before the solid path existed - checked by re-running every benchmark
deck and comparing each summary, CSV and VTK file with the earlier output.

## What each component owns

| Component | Owns | Deliberately does *not* |
|-----------|------|--------------------------|
| `core/` | scalar and matrix aliases, the exception hierarchy, the logger, timing, build provenance | anything numerical |
| `mesh/Mesh` | nodal coordinates (2-D or 3-D), flat connectivity with an explicit stride, validation, boundary edge / face extraction | materials, DOFs, physics |
| `mesh/StructuredMesh` | the rectangular Q4 and the box Hex8 generators, each with a node-perturbing variant for patch tests | anything unstructured |
| `mesh/SubMesh` | extracting an element subset with compact renumbering; edge- (face-) connected components; the density-to-solid interpretation | deciding *whether* to interpret |
| `material/IsotropicMaterial` | `(E, nu, rho)`, the plane-stress, plane-strain and three-dimensional matrices, input validation | element integration |
| `elements/Element` | the abstract kernel interface: stiffness, consistent mass, strain operator, edge or face traction, local face tables | the element's own geometry |
| `elements/Quad4`, `elements/Hex8` | the bilinear quadrilateral and the trilinear hexahedron | any other topology |
| `elements/Quadrature` | Gauss-Legendre rules on the line, the square and the cube | where they are used |
| `fem/DofManager` | DOF numbering, prescribed values, the free/prescribed partition, gather and scatter | assembly |
| `fem/Selector` | purely geometric region selection (box, circle, annulus, ids, nearest node), unions and complements | what a region is *for* |
| `fem/BoundaryConditions` | constraint and load specifications, and turning them into prescribed DOFs and a global force vector | solving |
| `fem/FemModel` | the complete discrete model: mesh, material, idealisation, integration orders, DOFs, load cases | any numerics |
| `fem/Assembler` | sparse assembly of `K` and `M`, block extraction, the uniform-mesh element cache | choosing scale factors |
| `fem/LinearSolver` | five solver backends behind one interface, pivot inspection, residual verification | model semantics |
| `fem/ModelDiagnostics` | per-component rigid-body and floating-region detection before any factorisation | fixing the model |
| `fem/StaticAnalysis` | the reduced solve, reaction recovery, global equilibrium checks, one factorisation shared across load cases | stress |
| `fem/StressRecovery` | strain, stress, von Mises, principal stresses, element strain energy, nodal averaging | plotting |
| `fem/ModalAnalysis` | the generalised eigenproblem, subspace iteration, validity screening, the analytical references | design variables |
| `topopt/DesignDomain` | design variables, bounds, passive tags, element volumes, feasibility checks | the objective |
| `topopt/SimpInterpolation` | `E(rho)`, its derivative, the mass laws | assembly |
| `topopt/DensityFilter` | the filter operator, its exact adjoint, the Sigmund variant | the objective |
| `topopt/Sensitivity` | the weighted-compliance objective, the exact gradients, the cached factorisation an adjoint solve reuses, and the finite-difference verifier | the update rule |
| `topopt/StressConstraint` | the relaxed, p-norm-aggregated von Mises constraint per load case, its adaptive scale and its adjoint gradient through the filter | the optimiser |
| `topopt/OptimalityCriteria` | one OC step with multiplier bisection | the loop |
| `topopt/Mma` | one MMA step: asymptotes, the convex subproblem and its primal-dual interior-point solve | the objective and constraints |
| `topopt/TopologyOptimizer` | the loop for either method, continuation, convergence, history, snapshots | I/O |
| `io/Json` | a self-contained JSON reader/writer and a path-aware, typo-catching config reader | the schema |
| `io/Config` | the input-deck schema, validation, and model construction | numerics |
| `io/CsvWriter`, `io/VtkWriter` | plain-text export with explicit precision | what to export |
| `io/StlWriter` | the boundary surface of a mesh (extruded for a plane one) as an indexed triangle surface, its closure and manifold checks, binary STL in and out | choosing what to export |
| `io/CalculixWriter` | one CalculiX input deck per load case for the same discrete problem | running CalculiX |
| `io/ResultWriter` | the result-directory layout, the geometry export and the summary documents | computing anything |
| `apps/` | argument parsing, orchestration, console reports, exit codes | physics |
| `python/sparlab_viz` | reading result files and drawing figures (plane fields in `fields`/`plots`, solid surfaces in `solid`/`plots3d`) | recomputing physics |
| `python/scripts/cross_validate.py` | driving CalculiX and scikit-fem on the exported problems and comparing nodal displacements | judging which code is right |

## Design decisions worth naming

**Dirichlet conditions by partitioning, not penalty.** The reduced system
`K_ff u_f = f_f - K_fp u_p` keeps `K_ff` symmetric positive definite, gives
exact reactions from the full residual `r = K u - f`, and leaves the eigenvalue
spectrum undistorted. A penalty method would have made the modal analysis wrong
in a way that is hard to notice.

**One factorisation per optimiser iteration.** `StaticAnalysis::prepare()`
factorises `K_ff` once and `solve_load_vector()` reuses it, so a three-load-case
design costs one factorisation and three back-substitutions per iteration rather
than three factorisations. Combined with self-adjointness (no adjoint solve for
the compliance gradient), that is what makes multi-load-case optimisation on a
38 000-element mesh a 5-minute job rather than an hour.

**Element-matrix cache for uniform grids.** On a uniform structured mesh every
cell is geometrically identical, so `K_e^0` and `M_e^0` are integrated once. The
generic per-element path is always compiled in and the test suite asserts the
two agree, so the optimisation cannot silently diverge from the general case.

**Purely geometric region selectors.** Boundary conditions, loads and passive
regions are all specified by geometry, never by node or element indices (though
indices are available). The same deck therefore applies unchanged at any mesh
resolution, and a topology extracted from a density field can be re-constrained
without re-authoring anything - which is what makes the static and modal
analysis of the interpreted structure possible at all. That re-analysis runs on
every topology job, because the compliance of the *thresholded* structure is
what the design would really deliver; if the extracted structure turns out not
to be analysable, the reason is recorded in the summary and warned about rather
than allowed to destroy a completed optimisation.

**Passive regions as equal bounds.** Pinning the design variable rather than the
physical density keeps the objective an exact function of the free variables,
which is what makes the finite-difference verification of the gradient
meaningful.

**A hand-written JSON parser.** The build needs nothing but Eigen and a C++17
compiler. The parser accepts `//` and `/* */` comments and a trailing comma,
because an annotated input deck is far easier to maintain, and it reports the
line and column of a syntax error. The `ConfigNode` wrapper records every key it
reads so unknown keys can be reported: a misspelled key that silently takes its
default is one of the easiest ways to publish a wrong number.

**The Python layer never computes physics.** It reads `mesh.json`,
`summary.json` and the CSV tables. Derived values in figures are limited to
ratios and least-squares slopes, and each is labelled as derived. That means a
number in a figure can always be traced to a solver run.

**One dimension-generic core, not two solvers.** Rather than a 3-D copy of
the 2-D code, the dimension is a run-time property of the mesh. What that
bought: every verification study, the optimiser, the writers and the
diagnostics run on the solid path unchanged, and the plane path is protected
by the fixed-size Q4 kernels and the byte-for-byte reproduction check. What it
cost: a handful of `dim`-dependent branches in the constitutive law, the
von Mises formula, the rigid-body-mode table and the face tables.

**MMA returns the iterate it judged.** Convergence is decided on the iterate
whose objective and constraints were just evaluated, and that iterate is the
design returned - not the unevaluated update, which under the objective-stall
criterion can still differ by up to the move limit and need not be feasible.
A stalled objective on an infeasible design is reported as non-convergence,
never as a result.

**Geometry export says what it is.** The STL surface is the boundary of the
thresholded cells, so it inherits the mesh's staircase and every
interpretation choice. The writer checks closure (every directed edge has its
reverse) and manifoldness (no edge used twice each way), compares the enclosed
volume with the cell volume, and records all of it in the summary; a surface
that is closed but non-manifold because cells touch along an edge is reported
as such rather than passed off as a solid.

## Extension points

The interfaces were shaped so the following are additions, not rewrites.

**A new element topology.** Implement `Element` (the kernels plus the local
face table), add an `ElementType` enumerator, and register it in
`make_element()`. That is exactly how `Hex8` was added: the layers above it
are written against `dim()`, `num_nodes()`, `num_faces()` and `num_voigt()`
and needed no change. `Mesh` already stores connectivity as a flat array with
an explicit stride; mixed topologies need that stride replaced by a
per-element offset table, which is the one change outside the new subclass.

**Plane strain.** Already implemented and unit-tested end to end: set
`model.stress_state` to `plane_strain`. The constitutive matrix, the von Mises
formula's out-of-plane term and the test coverage are all in place. Plane stress
remains the default and the idealisation the verification studies exercise, so
plane strain is documented as available rather than as validated to the same
depth.

**A new material model.** `FemModel` holds one `IsotropicMaterial` and hands its
constitutive matrix to the element kernels. An orthotropic or
temperature-dependent material replaces that with a small interface returning
`D` per element; nothing above the element layer depends on the material being
isotropic.

**An unstructured mesh.** `Mesh` accepts arbitrary coordinates and connectivity
today, and everything except `StructuredMesh`, the structured index helpers and
the `imshow` fast path in the Python layer already works on one. A reader for an
external mesh format slots in beside `make_structured_quad_mesh`.

**A different optimiser.** `ComplianceObjective::evaluate` returns the objective,
the volume constraint and both gradients; `optimality_criteria_update` and
`MmaOptimizer::update` are the two update rules, selected by
`optimizer.method`. An augmented-Lagrangian or a trust-region step would be
a third, again without touching the objective.

**A new constraint.** `StressConstraint` is the template: it takes the model,
the assembler and the filter, evaluates its value and gradient from an
`ObjectiveEvaluation` (reusing the cached factorisation for its adjoint solve
through `ComplianceObjective::solve_adjoint`) and hands one row of `fval` /
`dfdx` to MMA. A frequency constraint, a displacement bound or a second volume
budget on a region fits the same shape.

**A new objective.** The self-adjoint shortcut in `Sensitivity.cpp` is specific
to compliance. A non-self-adjoint objective (a displacement at a point, a
frequency) needs an adjoint solve, which the stress constraint already
exercises: the factorisation is cached and reused for the adjoint right-hand
side.

**A new external code for cross-validation.** `CalculixWriter` shows the
shape: write the same nodes, connectivity, supports and nodal loads in the
target's format, then add a `solve_with_<code>` function and a `compare`
entry in `python/scripts/cross_validate.py`, which already handles the
node-by-node comparison, the tolerances and the summary.

## File map

```
include/sparlab/          public headers, one per component, documented
src/                      implementations, mirroring include/
apps/                     four command-line drivers plus shared CLI support
tests/                    twelve Catch2 translation units plus shared fixtures
configs/benchmarks/       the six benchmark decks (four plane compliance cases,
                          the stress-constrained L-bracket, the solid bracket)
configs/studies/          the design-study baseline deck
configs/verification/     the plane and the solid static+modal decks
python/sparlab_viz/       loaders, style, field artists (plane and solid), figures
python/scripts/           figure and table drivers, the cross-validation driver
scripts/                  run scripts (benchmarks, verification, cross-validation,
                          study, scaling)
docs/                     this documentation
docs/figures/             committed figures and the animations
docs/results/             committed result tables, regenerated from summaries
results/                  generated output (git-ignored)
```
