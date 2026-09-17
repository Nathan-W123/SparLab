# Assumptions and limitations

This is what SparLab does **not** do, and what its results do **not** claim.
Nothing here is a bug list; every item is a deliberate scope boundary, and each
says what would be needed to lift it.

## Physics and idealisation

**Two dimensions only.** Plane stress (default) and plane strain. Out-of-plane
bending, warping, buckling and shell behaviour are absent. A bracket or rib
optimised here is a *planar load-path study*, not a substitute for a 3-D
analysis of the real part. Plate buckling in particular can govern a thin
lightened web, and this model cannot see it.

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

**Q4 elements only.** The four-node bilinear quadrilateral is a good workhorse
but is stiff in bending: it needs several elements through the thickness to
resolve a bending stress field. The mesh-convergence study quantifies this. A
Q8 or an enhanced-assumed-strain element would do better per DOF, and
`docs/architecture.md` describes what adding one involves.

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

**One inequality constraint only.** Optimality criteria handles exactly one -
here the volume. Stress constraints, per-mode frequency constraints, multiple
volume constraints on separate regions and manufacturing constraints all need a
different optimiser (MMA or augmented Lagrangian). The interfaces were shaped so
that swapping the update step is a local change, but the step itself is not
written.

**Stress is not constrained, and SIMP stress is delicate.** The reported
macroscopic stress in a low-density element is an average over a near-void cell
and is not a material stress. The figures therefore mask elements below the
interpretation threshold. Stress-constrained topology optimisation additionally
needs a stress-relaxation scheme (qp or similar) that is not implemented.

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

## Numerics

**Direct solvers by default.** `SimplicialLDLT` is excellent for 2-D problems up
to a few hundred thousand DOFs and is what makes the optimisation loop cheap
(one factorisation shared by every load case, and no adjoint solve). It is not
the right tool for a 3-D problem of the same element count, where a multigrid or
domain-decomposition preconditioner would be needed.

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
against exact answers - the patch test to 4e-15, solver agreement to 7e-12,
gradient finite differences to 2e-8, mass conservation to 5e-14. Validation
against independent theory covers exactly three references: Euler-Bernoulli and
Timoshenko cantilever deflection, Euler-Bernoulli bending frequencies, and
fixed-free rod axial frequencies. There is **no comparison against a commercial
FE code and no comparison against experiment.** Those are the two things that
would raise confidence further, and neither is here.

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

**The mass-stiffness comparison is 2-D specific.** In this idealisation `K` and
`M` both scale linearly with thickness, so a uniformly thinned plate has exactly
the full-solid natural frequencies while its compliance rises as `1/thickness`.
That is what makes "equal-mass uniform plate" a well-defined baseline, and the
runs verify the frequency invariance numerically rather than assuming it. The
argument does not carry over to a 3-D part, where thinning changes the bending
stiffness and the frequencies.

**Runtimes are machine-specific.** Every reported time comes from one machine
and one build, both recorded in the result summaries. Absolute values will
differ elsewhere; the scaling exponents are the portable part, and even those
carry cache effects at these sizes.
