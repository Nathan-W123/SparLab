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
|  ResultWriter |          |  OptimalityCriteria     |       |  LinearSolver   |
+---------------+          |  TopologyOptimizer      |       |  Assembler      |
                           +-------------------------+       |  FemModel       |
                                                             |  BoundaryConds  |
                                                             |  DofManager     |
                                                             |  Selector       |
                                                             +--------+--------+
                                                                      |
                                    +---------------------------------+
                                    |                 |               |
                                    v                 v               v
                           +----------------+  +-------------+  +------------+
                           |  elements/     |  |  material/  |  |  mesh/     |
                           |  Element (abc) |  |  Isotropic  |  |  Mesh      |
                           |  Quad4         |  |  Material   |  |  Structured|
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
       loaders -> style -> fields -> plots / studies
```

## What each component owns

| Component | Owns | Deliberately does *not* |
|-----------|------|--------------------------|
| `core/` | scalar and matrix aliases, the exception hierarchy, the logger, timing, build provenance | anything numerical |
| `mesh/Mesh` | nodal coordinates, flat connectivity with an explicit stride, validation, boundary-edge extraction | materials, DOFs, physics |
| `mesh/StructuredMesh` | the rectangular Q4 generator and a node-perturbing variant for patch tests | anything unstructured |
| `mesh/SubMesh` | extracting an element subset with compact renumbering; edge-connected components; the density-to-solid interpretation | deciding *whether* to interpret |
| `material/IsotropicMaterial` | `(E, nu, rho)`, the plane-stress and plane-strain matrices, input validation | element integration |
| `elements/Element` | the abstract kernel interface: stiffness, consistent mass, strain operator, edge traction | the element's own geometry |
| `elements/Quad4` | the bilinear quadrilateral | any other topology |
| `elements/Quadrature` | Gauss-Legendre rules on the square and the line | where they are used |
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
| `topopt/Sensitivity` | the weighted-compliance objective, the exact gradients, and the finite-difference verifier | the update rule |
| `topopt/OptimalityCriteria` | one OC step with multiplier bisection | the loop |
| `topopt/TopologyOptimizer` | the loop, continuation, convergence, history, snapshots | I/O |
| `io/Json` | a self-contained JSON reader/writer and a path-aware, typo-catching config reader | the schema |
| `io/Config` | the input-deck schema, validation, and model construction | numerics |
| `io/CsvWriter`, `io/VtkWriter` | plain-text export with explicit precision | what to export |
| `io/ResultWriter` | the result-directory layout and the summary documents | computing anything |
| `apps/` | argument parsing, orchestration, console reports, exit codes | physics |
| `python/sparlab_viz` | reading result files and drawing figures | recomputing physics |

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

## Extension points

The interfaces were shaped so the following are additions, not rewrites.

**A new element topology.** Implement `Element` (five kernels plus edge
metadata), add an `ElementType` enumerator, and register it in
`make_element()`. `Mesh` already stores connectivity as a flat array with an
explicit stride; mixed topologies need that stride replaced by a per-element
offset table, which is the one change outside the new subclass. `Assembler`,
`DofManager`, `StaticAnalysis`, `StressRecovery` and `ModalAnalysis` are written
against `num_nodes()` / `num_dofs()` and need no change.

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
the constraint and both gradients; `optimality_criteria_update` is a free
function taking those plus the bounds. An MMA or augmented-Lagrangian step
replaces the latter without touching the former, which is the change a stress or
frequency constraint would need.

**A new objective.** The self-adjoint shortcut in `Sensitivity.cpp` is specific
to compliance. A non-self-adjoint objective (a displacement at a point, a
frequency) needs an adjoint solve, which fits the same structure: the
factorisation is already cached and reusable for the adjoint right-hand side.

## File map

```
include/sparlab/          public headers, one per component, documented
src/                      implementations, mirroring include/
apps/                     four command-line drivers plus shared CLI support
tests/                    nine Catch2 translation units plus shared fixtures
configs/benchmarks/       the four benchmark decks
configs/studies/          the design-study baseline deck
configs/verification/     the small static+modal deck used by the README example
python/sparlab_viz/       loaders, style, field artists, figures
python/scripts/           figure and table drivers
scripts/                  run scripts (benchmarks, verification, study, scaling)
docs/                     this documentation
docs/figures/             committed figures and the animation
docs/results/             committed result tables, regenerated from summaries
results/                  generated output (git-ignored)
```
