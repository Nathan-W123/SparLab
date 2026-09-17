# SparLab

**A 2-D finite-element structural solver with density-based topology
optimisation, built for lightweight aerospace structures.**

Modern C++17 numerical core (Eigen sparse), Python visualisation, JSON input
decks, and a verification suite that compares every claim against either an
exact answer or an independent theory.

<p align="center">
  <img src="docs/figures/aerospace_bracket_evolution.gif"
       alt="A two-bolt aerospace bracket evolving from a solid rectangular design domain into an optimised truss topology over 375 SIMP iterations"
       width="620">
</p>

<p align="center">
  <em>A two-bolt mounting bracket with a pin lug, optimised for minimum weighted
  compliance under three load cases at a 35&nbsp;% volume constraint. 38&nbsp;400
  elements, 375 iterations, 1&nbsp;128 linear solves, 422&nbsp;s. Black is
  material, the page colour is void.</em>
</p>

The result: **6.5 % stiffer and 41 % higher in fundamental frequency than a
uniform plate of the same mass**, with the volume constraint met to 7×10⁻¹¹
relative.

---

## What it does

| Capability | Detail |
|------------|--------|
| **Finite elements** | Plane stress and plane strain, four-node bilinear quadrilaterals, Gauss-Legendre quadrature, sparse assembly, exact Dirichlet partitioning |
| **Loads** | Point loads and consistently integrated edge tractions, multiple load cases with weights |
| **Recovery** | Displacements, exact support reactions, element and nodal strain/stress, von Mises, principal stresses, element strain energy, compliance |
| **Modal analysis** | Consistent or lumped mass, generalised eigenproblem by shift-invert subspace iteration, validity screening for negative and rigid-body eigenvalues |
| **Topology optimisation** | SIMP with penalty continuation, density and sensitivity filters, analytical sensitivities, optimality criteria with multiplier bisection, passive solid/void regions, multi-load-case objective |
| **Verification** | Patch test, rigid-body modes, positive definiteness, reaction equilibrium, sparse-vs-dense agreement, finite-difference gradient check, mass conservation, beam and rod theory, mesh convergence |
| **Diagnostics** | Pre-solve detection of rigid-body under-constraint and floating regions, singular-matrix reporting with the likely modelling cause, explicit non-convergence reporting |
| **Output** | `summary.json`, CSV tables, legacy VTK for ParaView, publication-quality figures and animations |

## Quick start

```bash
# 1. Dependencies (Debian/Ubuntu; see scripts/setup_deps.sh for other platforms)
./scripts/setup_deps.sh

# 2. Build and test          (~1 min build, ~10 s tests)
make build
make test

# 3. Solve one case          (~0.1 s)
make benchmark CASE=cantilever_analysis

# 4. Verification studies    (~2 s, exits non-zero if any tolerance is missed)
make verify

# 5. Everything: all benchmarks, the design study, scaling, figures, tables
make all
```

Only Eigen 3.3+ and a C++17 compiler are required to build. Catch2 is used for
the tests and is fetched automatically if the system package is absent. The
Python layer needs `numpy`, `pandas`, `matplotlib` and `pillow`.

Every target is a thin wrapper over a script or a binary, so anything can also
be typed by hand:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/bin/sparlab_solve  --config configs/verification/cantilever_analysis.json --output results/x
./build/bin/sparlab_topopt --config configs/benchmarks/aerospace_bracket.json    --output results/y
./build/bin/sparlab_verify --study all --output results/verification
./build/bin/sparlab_bench  --sizes 20,40,80,160,320 --repeats 3 --output results/benchmark
```

Run any executable with `--help` for its full flag list.

## Commands

| Task | Command |
|------|---------|
| Configure and build | `make build` |
| Run all tests | `make test` |
| Run the verification / validation studies | `make verify` |
| Solve one benchmark | `make benchmark CASE=aerospace_bracket` |
| Run all benchmarks | `make benchmarks` |
| Run the aerospace design study | `make study` |
| Runtime scaling benchmark | `make scaling` |
| Regenerate all figures and animations | `make figures` |
| Refresh the committed result tables | `make results` |
| Everything, in order | `make all` |
| Clean | `make clean`, `make distclean` |

## Results

All numbers below are read from the `summary.json` of the run named beside them.
`docs/results/README.md` is the machine-generated version and is refreshed by
`make results`, so the documentation cannot drift from the solver output.

### Verification and validation

| Study | Kind | Metric | Value | Tolerance |
|-------|------|--------|------:|----------:|
| Patch test, distorted mesh | verification | max relative error in `u`, strain, stress | `4.07e-15` | `1e-10` |
| Solver agreement (5 backends) | verification | max relative difference vs dense LU | `7.04e-12` | `1e-8` |
| **Topology sensitivity vs central differences** | verification | max relative gradient error, best step | `2.18e-08` | `1e-5` |
| Mass conservation, `sum(M)/2 = rho V` | verification | relative error | `4.84e-14` | - |
| Mesh convergence | verification | observed order, tip deflection | `2.31` | - |
| Cantilever tip deflection vs Timoshenko | validation | relative difference, finest mesh | `4.78e-04` | `0.02` |
| Cantilever `f1` vs Euler-Bernoulli | validation | relative difference | `5.52e-03` | `0.02` |
| **Axial mode vs fixed-free rod theory** | validation | relative difference | `4.02e-06` | - |

`make verify` exits non-zero if any tolerance is missed, so it is a usable
numerical regression gate. Full detail, including what is *not* covered, in
[`docs/verification.md`](docs/verification.md).

### Benchmarks

| Case | Elements | `nu` | Iterations | Compliance [J] | Equal-mass plate [J] | Stiffness gain | Grey | Runtime [s] |
|------|---------:|-----:|-----------:|---------------:|---------------------:|---------------:|-----:|------------:|
| Cantilever beam | 12 800 | 0.40 | 362 | 1.11252 | 1.37813 | **1.239** | 0.095 | 44.4 |
| MBB beam | 10 800 | 0.50 | 434 | 221.551 | 259.521 | **1.171** | 0.279 | 50.2 |
| Aerospace bracket | 38 400 | 0.35 | 375 | 3.00015 | 3.19665 | **1.065** | 0.116 | 422.5 |
| Wing rib | 12 500 | 0.40 | 530 | 1.35608 | 1.21194 | **0.894** | 0.221 | 56.0 |

The volume constraint is satisfied to between `1.1e-11` and `7.2e-11` relative
in every case, and at every recorded iteration.

"Stiffness gain" compares against a *uniform plate of the same mass*, not
against the solid domain. In this 2-D idealisation `K` and `M` both scale
linearly with thickness, so such a plate has compliance `C_solid / nu` and
exactly the full-solid frequencies - the runs verify that frequency invariance
numerically (to `2.5e-12`) rather than assuming it.

The wing rib's gain below 1 is a real result, not a solver failure: diagnostic
runs show the optimiser beats uniform thinning by 4.6 % with full design
freedom, and that the *mandated attachment features* - spar pads and skin
flanges, forced solid and charged against the same volume budget - cost 16 %
between them. Measured on its **interpreted** structure instead of on the
objective, the rib comes out level with the plate rather than 10.6 % short of
it (1.006 after correcting for the 1.3 % mass thresholding adds) - the reported
objective is the pessimistic one, for the reason the design study takes apart.
See [`docs/benchmarks.md`](docs/benchmarks.md).

### Modal: what optimisation does to the spectrum

| Case | Structure | Mass [kg] | `f1` [Hz] | `f2` [Hz] |
|------|-----------|----------:|----------:|----------:|
| Cantilever | full solid / equal-mass plate | 0.899 / 0.360 | 872.1 | 3173.0 |
| Cantilever | optimised topology | 0.360 | **1016.4** | 1685.7 |
| Bracket | full solid / equal-mass plate | 1.349 / 0.472 | 1437.2 | 4252.3 |
| Bracket | optimised topology | 0.473 | **2027.4** | 3917.9 |
| Wing rib | full solid / equal-mass plate | 0.278 / 0.111 | 1488.7 | 2623.1 |
| Wing rib | optimised topology | 0.113 | **386.6** | 487.1 |

Compliance optimisation raises the fundamental frequency of the cantilever
(+16 %) and the bracket (+41 %) at equal mass, and **lowers** it sharply for the
wing rib (-74 %). In every case the *higher* modes fall, because a truss has
local member modes a continuous plate does not. A compliance objective sees only
the static load path; a design with a vibration requirement needs that
requirement in the optimisation, which this single-constraint optimiser cannot
carry. Discussion in [`docs/benchmarks.md`](docs/benchmarks.md) and
[`docs/aerospace_study.md`](docs/aerospace_study.md).

### Design study: 41 runs, 7 arms

The bracket is swept over volume fraction, load-case weighting, mesh resolution
(twice), SIMP penalty, filter radius and Young's modulus. All 41 points
converge on the stated criteria and meet the volume constraint. The findings
that changed how the rest of this project is set up:

| Finding | Measurement |
|---------|-------------|
| **The reported objective ranks the numerical settings backwards** | Reported stiffness gain spans **0.72-1.22** and crosses parity with a uniform plate eight times; the *interpreted* structures of the 38 sound points span **1.01-1.27** and never cross it. The three best points on the reported objective are the three worst structures in the study |
| A mandated feature can dominate the answer | The passive attachment collars are 3.9 % of the domain but **25.8 %** of the volume budget at `nu` = 0.15, which is most of why the payoff falls from 26.6 % to 1.0 % as the budget shrinks |
| Quote the filter radius in metres, never in cells | Radius fixed in metres, five meshes agree to **1.2 %**; radius fixed at 1.5 cells, the same meshes spread **28 %** and the topology chases the mesh |
| Check the filter is resolved and the penalty high enough | Three points produced designs that fall apart when thresholded - two from an identity filter (average support **1.000** element, up to **494** disconnected groups) and one from `p` = 1.5, whose interpretation is **6.75x** softer than the field it came from |
| A minimum length scale is not free | **36.9 %** more reported compliance between a 3.75 mm and a 15 mm filter radius |
| The local optimum costs something real | Two runs of the same problem reached designs **1.3 %** apart on the same objective |
| One load case was doing nothing | `up_reversal` is exactly `-0.5` times `down_limit`, so its compliance is `0.25 C_down` to ten digits and the compliance objective cannot see it |
| The runs are reproducible to the last bit | The study was run twice on two separately compiled binaries; all 41 points reproduced with a relative difference of **exactly zero** |

Full tables, figures and reasoning in
[`docs/aerospace_study.md`](docs/aerospace_study.md).

### Runtime scaling

Fitted slopes of `log(time)` against `log(DOFs)` over the largest three of eight
sizes, up to 194 922 DOFs:

| Phase | Slope | Expected |
|-------|------:|----------|
| Assemble `K` | 1.15 | `O(n)` plus cache effects |
| Sparse Cholesky factorisation | **1.57** | the `O(n^1.5)` of a good 2-D ordering |
| One back-substitution | 1.31 | - |
| Objective + gradient (one optimiser iteration) | 1.55 | factorisation-dominated |

One optimiser iteration costs essentially one factorisation, which is the
consequence of two design decisions: compliance is self-adjoint so no adjoint
solve is needed, and one factorisation of `K_ff` is shared across every load
case.

## Figures

| | |
|---|---|
| ![Bracket mesh and boundary conditions](docs/figures/aerospace_bracket_mesh_bcs.png) | ![Bracket optimised topology](docs/figures/aerospace_bracket_topology.png) |
| Design domain, rigid bolt attachments, three load cases, and the passive solid/void regions | The SIMP density field and its explicit interpretation as solid material |
| ![Sensitivity verification](docs/figures/verify_sensitivity.png) | ![Mesh convergence](docs/figures/verify_mesh_convergence.png) |
| Analytical gradients against central differences over eight step sizes, with every element tested | Tip deflection against Euler-Bernoulli and Timoshenko theory, plus self-convergence order |
| ![Bracket convergence history](docs/figures/aerospace_bracket_convergence.png) | ![Runtime scaling](docs/figures/runtime_scaling.png) |
| Compliance, volume constraint, design change and grey level against iteration | Wall-clock cost of each solver phase against problem size |
| ![Bracket modal comparison](docs/figures/aerospace_bracket_modal_comparison.png) | ![Mass-stiffness trade](docs/figures/study_pareto.png) |
| Natural frequencies before and after optimisation, at equal mass | Mass-stiffness trade and the vibration metric across the volume sweep |

Further figures in [`docs/figures/`](docs/figures/): deformed shapes,
displacement and stress fields, reaction forces, mode shapes, density evolution
montages, mesh-dependence and settings studies, and one animation per benchmark.

## Formulation

The full derivation is in [`docs/formulation.md`](docs/formulation.md); the
essentials:

**Element.** Bilinear Q4 shape functions on `(xi, eta) in [-1,1]^2`,

```
  N_1 = 1/4 (1-xi)(1-eta)   N_2 = 1/4 (1+xi)(1-eta)
  N_3 = 1/4 (1+xi)(1+eta)   N_4 = 1/4 (1-xi)(1+eta)
```

isoparametric mapping `x = sum_a N_a x_a`, Jacobian `J_ij = dx_i/dxi_j`, and

```
  K_e = sum_g w_g t detJ B^T D B        (2 x 2 Gauss, exact for the Q4)
  M_e = sum_g w_g rho t detJ N^T N      (3 x 3 Gauss)
  f_e = sum_g w_g t (L/2) N^T t_bar     (2-point, along an edge)
```

**Boundary conditions by partitioning, not penalty.** The reduced system is
`K_ff u_f = f_f - K_fp u_p` and the reactions come from the full residual
`r = K u - f`. This keeps `K_ff` symmetric positive definite, gives exact
reactions with no second assembly, and leaves the eigenvalue spectrum
undistorted.

**SIMP.** The modified law, so a design variable may reach exactly 0 or 1 while
`K` stays invertible:

```
  E(rho)/E_0     = eps + (1 - eps) rho^p ,        eps = E_min/E_0 = 1e-9
  d(E/E_0)/drho  = p (1 - eps) rho^(p-1)
```

Mass is interpolated separately, `m(rho) = eps_m + (1-eps_m) rho^p`, matching
the stiffness exponent so `E/m` stays bounded in void regions and no spurious
low-frequency mode appears.

**Filtering.** `rho_tilde = Hhat x` with row-normalised linear-hat weights
`H_ei = max(0, r_min - |x_e - x_i|)`. Because the map is linear the chain rule is
exact, which is what makes the finite-difference verification meaningful.

**Sensitivities.** Compliance is self-adjoint, so no adjoint solve is needed:

```
  dc/drho_e = -sum_l w_l  d(E(rho_e)/E_0)/drho_e  u_{l,e}^T K_e^0 u_{l,e}  <= 0
  grad_x c  = Hhat^T grad_rho c ,      grad_x g = Hhat^T v
```

**Optimality criteria.** The separable KKT condition
`B_e = (-dc/dx_e)/(lambda dg/dx_e) = 1` gives

```
  x_e^{k+1} = clip( x_e^k B_e^eta , max(x_lo, x_e^k - m) , min(x_hi, x_e^k + m) )
```

with `eta = 1/2`, move limit `m = 0.2`, and `lambda` found by geometric bisection
on the monotone volume map.

See [`docs/topology_optimization.md`](docs/topology_optimization.md) for the
filters, continuation, convergence criteria and diagnostics.

## Conventions

Strict SI everywhere: metres, newtons, pascals, kilograms, hertz, joules. No
unit-conversion helpers and no implicit scaling, so a value in a result file is
already in the unit its column header states.

Right-handed axes with the model in the x-y plane. Structured node numbering
`node(i,j) = j(nx+1)+i` and element numbering `elem(i,j) = j·nx+i`, both with x
fastest. Element nodes counter-clockwise from the lower-left corner. Two
translational DOFs per node, `dof(n,c) = 2n + c`. Forces and displacements
positive along `+x`/`+y`; tensile stress positive; tractions given as a global
stress vector, not a normal pressure. Voigt ordering
`{sxx, syy, sxy}` paired with *engineering* shear strain.

Full table, including every tolerance and its default, in
[`docs/conventions.md`](docs/conventions.md).

## Architecture

```
  apps/            sparlab_solve  sparlab_topopt  sparlab_verify  sparlab_bench
                            |
     +----------------------+----------------------+
     |                      |                      |
  io/                   topopt/                  fem/
  Json  Config          DesignDomain             StaticAnalysis  ModalAnalysis
  CsvWriter VtkWriter   SimpInterpolation        StressRecovery  ModelDiagnostics
  ResultWriter          DensityFilter            LinearSolver    Assembler
                        Sensitivity              FemModel  DofManager  Selector
                        OptimalityCriteria       BoundaryConditions
                        TopologyOptimizer                |
                                     +-------------------+-------------------+
                                     |                   |                   |
                                 elements/          material/             mesh/
                                 Element (abc)      IsotropicMaterial     Mesh
                                 Quad4  Quadrature                        Structured
                                     |                   |                SubMesh
                                     +---------+---------+-------------------+
                                               |
                                            core/  Types Exceptions Logging Timer

  python/sparlab_viz/   loaders -> style -> fields -> plots / studies
                        (reads only what io/ writes; never recomputes physics)
```

No cycles, no upward dependencies. Each component's responsibilities, the
design decisions behind them, and what it takes to add a new element type, a new
material model, an unstructured mesh or a different optimiser are in
[`docs/architecture.md`](docs/architecture.md).

## Input decks

A run is fully determined by one JSON file. Comments and trailing commas are
accepted; syntax errors report their line and column; **unknown keys are
reported**, and `--strict-config` makes them an error - a misspelled key that
silently takes its default is one of the easiest ways to publish a wrong number.

```json
{
  "name": "cantilever_beam",
  "mesh":     { "type": "structured_quad", "nx": 160, "ny": 80,
                "lx": 0.40, "ly": 0.20 },
  "material": { "name": "Al 7075-T6", "youngs_modulus": 71.7e9,
                "poisson_ratio": 0.33, "density": 2810.0 },
  "model":    { "thickness": 0.004, "stress_state": "plane_stress" },
  "boundary_conditions": [
    { "name": "clamped_root", "fix": ["x", "y"],
      "region": { "box": { "xmax": 0.0 } } }
  ],
  "load_cases": [
    { "name": "tip_down", "weight": 1.0,
      "point_loads": [
        { "force": [0.0, -2000.0], "distribution": "total",
          "region": { "box": { "xmin": 0.40, "ymin": 0.0975, "ymax": 0.1025 } } }
      ] }
  ],
  "modal":    { "enabled": true, "num_modes": 6 },
  "topology": { "enabled": true, "volume_fraction": 0.4,
                "filter": { "type": "density", "radius_elements": 2.0 } }
}
```

Regions are purely geometric - box, circle, annulus, nearest node, explicit ids,
unions and complements - so the same deck applies unchanged at any mesh
resolution. Every field is documented in
[`docs/configuration.md`](docs/configuration.md).

Shipped decks: [`configs/benchmarks/`](configs/benchmarks/) (the four benchmark
cases), [`configs/studies/`](configs/studies/) (the design-study baseline),
[`configs/verification/`](configs/verification/) (a small static + modal deck).

## Output

Each run writes a self-describing directory:

```
results/<case>/
  summary.json              every scalar, tolerance, timing and provenance field
  config.json               verbatim echo of the input deck
  mesh.json                 nodes, connectivity, prescribed DOFs, applied loads
  displacement_<lc>.csv     nodal displacements per load case
  stress_<lc>.csv           element strains, stresses, von Mises, energy
  reactions_<lc>.csv        support reactions
  fields_<lc>.vtk           ParaView fields (point and cell data)
  modes*.csv, mode_*.vtk    eigenvalues, frequencies, residuals, mode shapes
  history.csv               optimisation iteration history
  density_final.csv         final design, physical density, passive tags
  density_history.csv       density snapshots for the animation
```

`summary.json` is the single source of truth for every number quoted in the
documentation.

## Engineering integrity

Things this project deliberately does, because the opposite is easy and wrong:

* **a density field is never presented as geometry without saying how it was
  converted.** Every topology summary records the threshold, the elements
  retained, the number of disconnected groups and the material discarded as
  islands. Every run then re-solves that *extracted* structure with real
  material and reports its own compliance, mass and peak stress, because the
  objective the optimiser minimised is not the compliance of the part a
  threshold would produce - on the bracket the extracted structure comes out
  9-29 % stiffer. The modal analysis of an "optimised structure" runs on the
  same sub-mesh;
* **no manufacturability is claimed.** A filter radius is a minimum length
  scale, not a manufacturing constraint. No draw direction, tool access, wall
  thickness, overhang angle or fillet is modelled;
* **no comparison value is invented.** The MBB compliance is reported as
  computed, with the density-, sensitivity- and no-filter variants side by side,
  and without claiming agreement with a literature number that was not
  re-derived under identical settings;
* **failures are reported, not hidden.** Singular matrices, non-positive
  Cholesky pivots, residuals above tolerance, non-convergence, disconnected
  topologies and unknown configuration keys all produce a specific, actionable
  message - and the optimiser records *which* convergence criterion fired and
  the final value of both;
* **every tolerance is recorded** in the run summary;
* **runs are deterministic.** The only randomness is a seeded filler in the
  subspace-iteration starting basis and the node-perturbation generator, both
  with explicit default seeds.

## Assumptions and limitations

The short version: 2-D only, linear and small-strain, static plus undamped free
vibration, no body forces, Q4 elements only, attachment features are
*representations* and not joint models, load cases are illustrative and not
flight loads, the objective is weighted compliance and not a load envelope, one
inequality constraint only, no stress constraint, no manufacturability
modelling, no comparison against a commercial code or experiment.

The long version, with what it would take to lift each item, is in
[`docs/limitations.md`](docs/limitations.md). It is worth reading before
treating any number here as a design answer.

## Documentation

| Document | Contents |
|----------|----------|
| [`docs/formulation.md`](docs/formulation.md) | continuum problem, shape functions, quadrature, assembly, solvers, stress recovery, modal algorithm |
| [`docs/topology_optimization.md`](docs/topology_optimization.md) | SIMP, filters, sensitivity derivation, optimality criteria, passive regions, continuation, convergence, diagnostics |
| [`docs/conventions.md`](docs/conventions.md) | units, coordinates, numbering, signs, energy definitions, tolerances, determinism |
| [`docs/architecture.md`](docs/architecture.md) | layering, component responsibilities, design decisions, extension points |
| [`docs/configuration.md`](docs/configuration.md) | complete input-deck reference and the command-line overrides |
| [`docs/verification.md`](docs/verification.md) | every verification and validation check, with measured values and what is not covered |
| [`docs/benchmarks.md`](docs/benchmarks.md) | the four benchmark cases in detail, convergence behaviour, runtime and scaling |
| [`docs/aerospace_study.md`](docs/aerospace_study.md) | the parametric design study: mass-stiffness trade, load weighting, mesh dependence, penalty, filter radius, material stiffness |
| [`docs/limitations.md`](docs/limitations.md) | assumptions and scope boundaries |
| [`docs/results/README.md`](docs/results/README.md) | machine-generated result tables |

## Continuous integration

[`.github/workflows/ci.yml`](.github/workflows/ci.yml) runs three jobs:

* **build and test** on GCC and Clang, in Release and Debug, with
  `-DSPARLAB_WARNINGS_AS_ERRORS=ON`. The Debug build enables Eigen's own
  assertions, which is the configuration most likely to catch an indexing
  mistake;
* **verification** runs all the studies and fails the build if any documented
  tolerance is missed, uploading the summary either way;
* **benchmark** runs a static+modal analysis, a reduced-mesh topology
  optimisation including the modal comparison of the interpreted topology, the
  scaling benchmark, and regenerates the figures - so a break in the whole
  pipeline, not just the library, is caught.

## Licence

MIT. See [`LICENSE`](LICENSE).
