# Assumptions and limitations

This is what SparLab does **not** do, and what its results do **not** claim.
Nothing here is a bug list; every item is a deliberate scope boundary, and each
says what would be needed to lift it.

## Physics and idealisation

**Plane or solid continuum, nothing in between.** Plane stress (default),
plane strain, and three-dimensional solids on hexahedral meshes. There are no
plate, shell or beam elements, so a thin-walled part is either a plane
load-path study or a solid mesh with enough elements through the wall to
resolve its bending - which the solid element needs several of (below). No
model here sees buckling: plate buckling in particular can govern a thin
lightened web, and neither the plane nor the solid linear analysis can see
it.

**Linear, small strain, small displacement.** One factorisation, one solve, no
load stepping. Geometric non-linearity, plasticity, contact, creep and thermal
strain are all absent. The deformation figures are exaggerated by a stated
factor purely for visibility; the analysis behind them is linear.

**Static and undamped free vibration only.** No transient response, no damping,
no forced response, no fatigue. A natural frequency here is the undamped
eigenvalue of the constrained model.

**No body forces.** Loads enter as boundary tractions and concentrated nodal
forces. Self-weight and inertia relief are not implemented; the wing-rib deck's
"fuel inertia" case is a *representative edge pressure*, not a body-force
calculation.

**Q4 and Hex8 elements only.** The four-node bilinear quadrilateral and the
eight-node trilinear hexahedron are good workhorses but are stiff in bending
(shear locking): they need several elements through a bending depth to
resolve the stress field. The mesh-convergence studies quantify this - the
coarsest solid cantilever is 16 % too stiff, the finest within 0.5 %. A Q8 /
Hex20, an enhanced-assumed-strain or a B-bar element would do better per
DOF, and `docs/architecture.md` describes what adding one involves.

**Structured meshes only.** Both generators produce axis-aligned grids of
identical cells (optionally perturbed for the patch tests). Everything above
the mesh layer accepts arbitrary connectivity, but there is no reader for an
external mesh, so a curved boundary is a staircase of cells and a hole is a
passive void region.

**Point loads are singular.** A concentrated force on a node produces a stress
that does not converge under refinement - the peak grows without bound. Reported
peak stresses near a point load are mesh-dependent artefacts. The benchmark
decks spread their resultants over a small region for this reason, and the
distributed-traction path is the honest way to apply a pressure.

**Plane strain is implemented but less exercised.** The constitutive matrix, the
out-of-plane stress in the von Mises formula and the unit tests are all in
place, and the test suite confirms plane strain is stiffer than plane stress for
the same beam. But every verification and validation study runs in plane stress,
so plane strain is documented as *available*, not as validated to the same
depth.

## Attachment and joint representation

The bolt holes, bearing collars and pin lug in the aerospace bracket are
**representations**, not joint models:

* a bolt shank is idealised as perfectly rigid - every node inside a bolt hole
  is fully fixed. Bearing compliance, fastener flexibility, hole clearance and
  bolt preload are absent;
* the lug load is applied to the ring of nodes around the hole as if bearing
  pressure were uniform. A real pin joint has a cosine-like bearing distribution
  and local yielding;
* passive solid collars stand in for machined bosses. Their size was chosen to
  look like plausible hardware, not derived from a bearing-stress allowable;
* no fretting, no bushing, no through-thickness load transfer.

A real joint analysis needs 3-D contact and a fastener model. Nothing here
should be read as a joint substantiation.

## Loads

**The load cases are illustrative.** The bracket's 9 kN down-load, its 4.5 kN
reversal and its 4.5 kN lateral case, and the wing rib's skin pressures, are
plausible magnitudes chosen to exercise multi-load-case optimisation. They are
**not** derived from a flight-loads analysis, not factored to any airworthiness
requirement, and carry no limit/ultimate distinction.

**The objective is weighted compliance, not a load envelope.** Minimising
`sum_l w_l C_l` makes the design stiff on *average* across the weighted cases.
It does not guarantee any single case meets a deflection limit, and it is not
the same as optimising for the worst case. A worst-case (min-max) formulation
needs a different objective and a bound-constrained optimiser.

## Topology optimisation

**A density field is not geometry.** The SIMP result is a relative material
distribution over a fixed mesh. Elements at intermediate density represent
neither solid nor void; they are an artefact of relaxing a 0/1 problem into a
continuous one. SparLab never presents a density field as a solid body without
saying how it was converted: the "interpretation" block in every topology
summary records the threshold, how many elements survived, how many
disconnected groups they formed and how much material was discarded as islands.
The modal comparison of an "optimised structure" is run on that explicitly
extracted sub-mesh.

**Nothing here is a manufacturability claim.** The filter radius imposes a
minimum *length scale*, which is not the same as a manufacturing constraint. The
code models no draw direction, no milling-tool access, no minimum wall
thickness in the third dimension, no overhang angle for additive manufacture, no
support-structure cost, no fillet radii, and no machining setup. The designs
should be read as load-path guidance for a designer, not as parts.

**The optimum is local.** SIMP with `p > 1` is non-convex; a different starting
design, penalty schedule or mesh can converge to a different local optimum. The
design study shows exactly that: the mesh-refinement arm with the filter radius
fixed in cells produces visibly different topologies at different resolutions.
Continuation reduces the dependence but does not remove it.

**Grey boundaries are intrinsic to density filtering.** A density filter of
radius `r_min` leaves an intermediate band roughly `r_min` wide at every
material boundary. On a coarse mesh that band is a large fraction of the domain,
which is why the MBB benchmark settles near a grey level of 0.28. A Heaviside
projection would sharpen it; it is not implemented.

**Two constraint types.** Optimality criteria handles the volume constraint
alone; MMA handles the volume plus one aggregated stress constraint per load
case. Per-mode frequency constraints, displacement bounds, multiple volume
budgets on separate regions and manufacturing constraints are not written -
`docs/architecture.md` says what each would take, and `StressConstraint` is
the template - and the OC path remains the default because it converges with
no tuning where MMA needs a move limit chosen per problem.

**The stress constraint bounds a relaxed aggregate, not a part's stress.**
The constrained quantity is `rho^q sigma_vm` at each element centre,
aggregated with a p-norm whose scale is re-fitted every iteration. Three
consequences: intermediate-density elements carry a stress the material would
not; the aggregate underestimates the true maximum between re-fits, so the
relaxed peak hovers about the limit rather than sitting under it; and a
stress concentration the mesh does not resolve is not bounded at all - a
re-entrant corner in a density design is mesh sensitive whether or not the
constraint is on. Every constrained run therefore re-solves its thresholded
structure with full material and reports that peak against the limit, which
is a consistency check, not a yield or fatigue substantiation. The reported
macroscopic stress in a low-density element is still an average over a
near-void cell and not a material stress, and the figures still mask elements
below the interpretation threshold.

**MMA needs tuning where OC does not.** With a moving constraint surface the
default move limit of 0.2 left the stress-constrained L-bracket oscillating
at its iteration cap; 0.1 converged. Asymptote parameters and the p-norm
exponent are further knobs, each recorded in the summary. A run that hits the
cap says so, but it does not say which knob to turn.

**The design change stalls on fine meshes.** OC's `max |dx|` measure does not
reach `1e-2` on the finer benchmark meshes because thin members migrate one cell
at a time long after the objective has settled. A measured example: on the
cantilever benchmark the compliance improves by only 0.35% between iteration 250
and iteration 800 while `max |dx|` keeps oscillating between 0.004 and 0.12.
This is why a relative-objective-change criterion is the backstop, and why every
run records both indicators and which one fired.

**Compliance is not monotone by construction.** OC is a fixed-point update with
a move limit, not a line-search method, and each continuation stage raises the
penalty, which raises the compliance of a fixed density field. Histories are
recorded and plotted unsmoothed.

**The exported geometry is the interpretation, staircase and all.**
`structure_after.stl` is the boundary of the cells at or above the threshold
- a voxel surface, not a smoothed or fitted one - and inherits every choice
the interpretation made. Two retained cells that touch only along an edge
leave a non-manifold edge in it (32 on the 3-D bracket); the writer reports
the count and a slicer may split the surface into shells there. A plane
design is extruded by its thickness. None of this is a printable or
machinable part: no smoothing, no minimum feature size beyond the filter
radius, no fillets, no draft.

## Numerics

**Direct solvers by default, which bounds the solid problem size.**
`SimplicialLDLT` is excellent for 2-D problems up to a few hundred thousand
DOFs and is what makes the optimisation loop cheap (one factorisation shared
by every load case, and for compliance no adjoint solve). In 3-D its fill-in
grows as measured `n^2.3`: 9.5 s per factorisation at 28 413 DOFs against
0.09 s for a plane mesh of the same size, so the solid benchmark stops at
4 096 cells and 2.8 s per iteration. A finer solid study needs an iterative
solver with a multigrid preconditioner, which is not implemented; the
Jacobi-preconditioned CG that is available is a poor fit for the
`1e-9`-stiffness-floor conditioning of a SIMP system.

**Single-threaded assembly.** Eigen's own operations may use OpenMP where it is
available, but the element loop is serial: a parallel assembly needs either
per-thread triplet buffers or a colouring, and the extra complexity was not worth
it at these problem sizes. The reported runtimes are therefore single-core
numbers.

**The SIMP stiffness floor limits conditioning.** With `emin_ratio = 1e-9` the
smallest LDL^T pivot ratio sits around `1e-9`, comfortably above the `1e-14`
singularity threshold, but the condition number of `K_ff` is correspondingly
large. That is why the direct solver is the default for topology runs; the CG
solver would need many more iterations there.

**Mechanism modes are not detected analytically.** The pre-solve diagnostics
detect rigid-body under-constraint exactly, per connected element group, and
detect floating regions. An internal mechanism - two blocks joined at a single
node, say - is not caught analytically; it surfaces as a non-positive Cholesky
pivot, which the solver reports as a singular system with the same hint text.

**Nodal stress averaging is simple, not superconvergent.** Nodal values are
area-weighted averages of adjacent element values, used for smooth contours
only. A superconvergent patch recovery would give better nodal accuracy; element
values, which are what the CSV tables and the reported peaks use, are unaffected.

## Verification and validation

**Verification is strong; validation is narrow.** The implementation is verified
against exact answers in both dimensions - the patch tests to 4e-15 and
2e-15, solver agreement to 7e-12, gradient finite differences to 2e-8 and
2e-8, mass conservation to 5e-14 - and its linear static displacements are
cross-validated node by node against two independent codes: scikit-fem
agrees to solver round-off (1e-10 and 1e-11 relative) and CalculiX to the
rounding of its own result file (3e-6). Validation against independent
theory covers exactly three references: Euler-Bernoulli and Timoshenko
cantilever deflection, Euler-Bernoulli bending frequencies, and fixed-free
rod axial frequencies. The cross-validation covers displacements on two
problems only; stresses, frequencies and optimised designs are not compared
with another code, and there is **no comparison against experiment**.

**The MBB compliance is not compared to a published number.** SparLab reports
what it computes (218.8 J with the density filter, 203.2 J with the sensitivity
filter, at the 60x20-equivalent radius on a 180x60 mesh). Which value a
published result corresponds to depends on the filter variant, the radius
convention, the convergence criterion and the penalty schedule, so quoting
agreement without re-deriving the reference under the same settings would be
guesswork. The internal consistency checks - the verified gradient, the exact
volume constraint, the KKT spread across interior elements - are what support
the result.

## Reporting

**The mass-stiffness comparison is 2-D specific.** In the plane idealisation
`K` and `M` both scale linearly with thickness, so a uniformly thinned plate
has exactly the full-solid natural frequencies while its compliance rises as
`1/thickness`. That is what makes "equal-mass uniform plate" a well-defined
baseline, and the runs verify the frequency invariance numerically rather
than assuming it. A solid has no thickness to thin, so the 3-D bracket
reports only the full solid domain as its reference and no "stiffness gain";
its `mass_stiffness_comparison.note` says so. The baseline is also not fair
where passive voids remove part of the domain from the design but not from
the plate, which is why the L-bracket quotes its unconstrained run instead.

**Solid density designs are greyer, and the re-solve matters more.** A
1.5-cell filter averages 17 hexahedra rather than 9 quadrilaterals, and the
intermediate boundary layer is a surface rather than a curve, so the 3-D
bracket settles at a grey level of 0.34 and its thresholded structure is 39 %
stiffer than the SIMP field. The reported objective of a solid density design
is a poor estimate of what the interpreted part delivers; the re-solved value
is the one to plan with, and a Heaviside projection (not implemented) would
narrow the gap.

**Runtimes are machine-specific.** Every reported time comes from one machine
and one build, both recorded in the result summaries. Absolute values will
differ elsewhere; the scaling exponents are the portable part, and even those
carry cache effects at these sizes.
