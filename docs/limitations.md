# Assumptions and limitations

This is what SparLab does **not** do, and what its results do **not** claim.
Nothing here is a bug list; every item is a deliberate scope boundary, and each
says what would be needed to lift it.

## Physics and idealisation

**Plane or solid continuum, nothing in between.** Plane stress (default),
plane strain, and three-dimensional solids on hexahedral or tetrahedral
meshes. There are no
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

**Linear elements only.** The four-node quadrilateral, the three-node
triangle, the eight-node hexahedron and the four-node tetrahedron. All four
are stiff in bending and need several elements through a bending depth; the
constant-strain simplices are much stiffer than the bilinear and trilinear
elements. The mesh-convergence studies quantify it: the coarsest solid
cantilever is 16 % too stiff on Hex8 and 50 % on Tet4, and on the same nodes
the Tet4 error is about five times the Hex8 error. A tetrahedral mesh of a
bending-dominated part therefore needs to be fine, and its compliance reads
low until it is. Quadratic elements (Q8, Tri6, Hex20, Tet10), enhanced
assumed strain or B-bar would do better per DOF; `docs/architecture.md` says
what adding one involves.

**Meshes: one linear cell type, straight-sided.** The structured generators
make boxes of quadrilaterals, triangles, hexahedra or tetrahedra. Real
geometry comes in through Gmsh (MSH 2.2 / 4.1, ASCII) and Abaqus / CalculiX
`.inp` files, with these limits: one cell type per mesh (a quad-dominant or
a hex-dominant mesh with prisms is refused, with the re-export that fixes
it); linear cells only, so a curved boundary is a polygon of cell edges;
Abaqus parts must be flat (no instances with translations, one part);
supports, loads and materials in the file are ignored in favour of the deck;
and SparLab does not mesh - the two committed parts come from
`python/scripts/make_meshes.py`, which needs Gmsh. There is no adaptive
refinement, and on a graded mesh a filter radius given in cells
(`radius_elements`) multiplies the *mean* cell size, so it is a different
physical radius in the fine and the coarse regions than a reader might
assume; a radius in metres avoids that.

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

**The projection sharpens the design; it does not give it a length scale.**
A density filter alone leaves an intermediate band about one radius wide at
every material boundary (the MBB benchmark settles at a grey level of 0.28,
the solid bracket at 0.34). The Heaviside projection removes most of it - the
projected solid bracket ends at 0.013 - and with it the gap between the
objective and the thresholded structure. What it does not do is guarantee a
minimum member size: a single projection at `eta = 0.5` lets features
thinner than the filter radius survive, and the robust formulation that does
guarantee one (eroded, intermediate and dilated designs optimised together)
is not implemented. At a large `beta` the optimisation also becomes
strongly non-convex: single elements at the solid-void interface flip by the
full move limit while the compliance stands still, which is why the
projected decks stop on the objective-stall criterion rather than on the
design change. The `beta` schedule is a further setting the result depends
on, and it is recorded with every run.

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
at its iteration cap. At 0.1 it still spends long stretches in a cycle of
period 5 at the move limit and needs about 500 iterations to converge on the
design change; one of its subproblems needed more Newton iterations than the
original budget of 200 per barrier level. Asymptote parameters, the p-norm
exponent and the subproblem budget are further knobs, each recorded in the
summary. A run that hits a cap says so, but it does not say which knob to
turn.

**The design change stalls on fine meshes.** OC's `max |dx|` measure does not
reach `1e-2` on the finer benchmark meshes because thin members migrate one cell
at a time long after the objective has settled. A measured example: on the
cantilever benchmark the compliance improves by only 0.35% between iteration 250
and iteration 800, while `max |dx|` stays between 0.004 and 0.09 in nine
iterations out of ten and reaches the move limit of 0.2 around iteration 260.
This is why an objective-stall criterion is the backstop - the relative spread
of the compliance over a window, which an oscillating design cannot satisfy
mid-cycle - and why every run records both indicators and which one fired.

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

**The automatic solver choice is a size rule.** `auto`, the default,
factorises up to 50 000 free unknowns in 2-D and 10 000 in 3-D and uses
multigrid CG above that. The limits sit above the measured single-solve
crossover on purpose (`docs/benchmarks.md`, *Runtime and scaling*): a
factorisation serves every load case and every modal solve, and it has no
iteration count that can go wrong. They are not tuned per problem. A deck
that names a solver gets that solver, and `solver.linear.auto_direct_limit`
moves the limits.

**The direct solver bounds the solid problem size.** `SimplicialLDLT` with
AMD ordering is what makes a small optimisation loop cheap, with one
factorisation shared by every load case and no adjoint solve for
compliance. In 3-D its measured cost grows like `n^2.35`: 43 s per
factorisation at 47 775 DOFs, with 789 factor entries per DOF. No run here
factorises a solid system above that size. A nested-dissection ordering
(METIS) or a supernodal factorisation (CHOLMOD) would move the limit, and
neither is a dependency.

**Multigrid is tested on what the benchmarks need, not in general.** The
smoothed-aggregation hierarchy uses the rigid-body modes as its near-null
space, and its strength threshold keeps solid and SIMP void in separate
aggregates. On the problems tested the CG iteration count barely moves
with the mesh: 14 to 17 on Q4 and Hex8 blocks up to 830 115 unknowns, 16
to 21 on Tet4, and 18 to 24 per solve inside the optimisation runs with
their `1e-9` stiffness floor. Not tested: nearly incompressible material
(`nu` near 0.5), where the rigid-body modes are not a sufficient near-null
space; badly shaped or strongly graded meshes; and plane strain. A solve
that does not converge is reported with the residual reached and the keys
to change, never returned as an answer. The setup is not free either. For
one solve from scratch, Jacobi-preconditioned CG is as fast as multigrid up
to 28 413 unknowns on the Hex8 block, and a modal analysis with CG repeats
the whole solve for every vector of the subspace. The two modal analyses
of the 356 475-DOF bracket take a fifth of its run.

**Element loops are serial; only the iterative solvers are threaded.**
Assembly, sensitivities, filtering and stress recovery loop over elements
on one core, and the sparse Cholesky factorisation is serial too. The
multigrid kernels use OpenMP, and Eigen threads the matrix-vector product
of its Jacobi-preconditioned CG. Both give the same bits on one thread as
on many; a test asserts it for multigrid. Direct-solver runtimes are
therefore single-core numbers and iterative-solver runtimes are
four-thread numbers. Parallel element loops would need per-thread buffers
or a colouring to keep the results deterministic.

**The SIMP stiffness floor limits conditioning.** With `emin_ratio = 1e-9`
the smallest LDL^T pivot ratio sits around `1e-9`, comfortably above the
`1e-14` singularity threshold, but the condition number of `K_ff` is
correspondingly large. The multigrid solver copes because void and solid do
not share aggregates. Jacobi-preconditioned CG is kept as a baseline and is
measured on uniform material only; no benchmark uses it.

**Mechanism modes are not detected analytically.** The pre-solve
diagnostics detect rigid-body under-constraint exactly, per connected
element group, and detect floating regions. An internal mechanism, such as
two blocks joined at a single node, is not caught analytically. With the
direct solver it surfaces as a non-positive Cholesky pivot, which the
solver reports as a singular system with the same hint text. With the
multigrid solver it should surface as a singular coarsest level or a solve
that does not converge. A rigid-body mode reaching the coarsest level is
tested, an internal mechanism is not, so a model suspected of one is better
checked with the direct solver.

**Nodal stress averaging is simple, not superconvergent.** Nodal values are
area-weighted averages of adjacent element values, used for smooth contours
only. A superconvergent patch recovery would give better nodal accuracy; element
values, which are what the CSV tables and the reported peaks use, are unaffected.

## Verification and validation

**Verification is strong; validation is narrow.** The implementation is
verified against exact answers on all four element types. The patch tests
pass to 4e-15 (Q4), 2e-15 (Hex8) and 9e-15 (Tri3 and Tet4). The sparse and
dense solvers agree to 9e-12, and multigrid CG agrees with Cholesky to
4e-12. The compliance gradient matches central differences to 2e-8 on Q4
and Hex8, and to 1e-7 through the Heaviside projection on Q4 and Tet4.
Mass is conserved to 5e-14. Linear static displacements are
cross-validated node by node against two independent codes on seven
problems with twelve load cases, covering all four element types and both
mesh-file formats. scikit-fem agrees to solver round-off, 1.5e-10 or
better. CalculiX agrees to the rounding of its own result file, 4.3e-6 or
better, wherever the two codes solve the same discrete problem; its
plane-stress comparisons at `nu != 0` are recorded but not judged.
Validation against independent theory covers exactly three references:
Euler-Bernoulli and Timoshenko cantilever deflection, Euler-Bernoulli
bending frequencies, and fixed-free rod axial frequencies. Stresses,
frequencies and optimised designs are not compared with another code, and
there is **no comparison against experiment**.

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
than assuming it. A solid has no thickness to thin, so the solid cases
report only the full solid domain as their reference and no "stiffness
gain"; their `mass_stiffness_comparison.note` says so. The baseline is also
not fair where passive voids remove part of the domain from the design but
not from the plate, which is why the L-bracket quotes its unconstrained run
instead. The lug bracket's holes and solid rims are in both the part and the
plate, so there the comparison is like for like.

**Without a projection, solid density designs are greyer and the re-solve
matters more.** A 1.5-cell filter averages 17 hexahedra rather than 9
quadrilaterals, and the intermediate boundary layer is a surface rather
than a curve. The 3-D bracket without a projection therefore settles at a
grey level of 0.34, and its thresholded structure is 39 % stiffer than the
SIMP field. The Heaviside projection closes most of that gap: the projected
bracket ends at a grey level of 0.013 with its thresholded structure within
2 % of the objective, and the 356 475-DOF bracket within 0.03 %. The engine
mount, whose projection stops at `beta = 16` on a tetrahedral mesh, still
differs by 4.5 %. The re-solved value is the one to plan with, because it
is the analysis of the part that is actually exported.

**A projected run reports the projected density.** With the projection on,
the compliance, the volume, the grey level and every figure refer to the
projected density, the one the finite-element model sees. The grey level
before projection is recorded beside it (`filtered_grey_level`) and is much
higher: 0.12 to 0.35 on the benchmarks. The volume constraint holds at every
design update, but two kinds of iterate are analysed before an update has
put them on the target: the uniform start seen through the projection, and
the first design after each increase of `beta`. Their recorded volume is off
by up to 16 % at the start (the lug bracket) and up to 1.8 % after a `beta`
step (the projected solid bracket). The convergence figure says so beside
those points.

**Runtimes are machine-specific.** Every reported time comes from one machine
and one build, both recorded in the result summaries. Absolute values will
differ elsewhere; the scaling exponents are the portable part, and even those
carry cache effects at these sizes.
