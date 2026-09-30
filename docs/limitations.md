# Assumptions and limitations

This is what SparLab does **not** do, and what its results do **not** claim.
Nothing here is a bug list; every item is a deliberate scope boundary, and each
says what would be needed to lift it.

## Physics and idealisation

**Plane, solid, shell and beam models, one kind at a time.** Plane stress
(default), plane strain, three-dimensional solids on hexahedral or
tetrahedral meshes, four-node MITC4 shells for linear statics, natural
frequencies, the harmonic response and linear buckling, and two-node
Timoshenko beams in space for the same and the transient analysis (below).
No model mixes element types, so a stiffened panel is either all shell -
the stiffeners as shell walls - or a solid mesh, and a frame cannot carry a
shell deck or be tied to a solid. A plane model cannot
buckle out of its plane, so plate buckling of a thin lightened web - which
can govern it - stays invisible to the 2-D decks; a shell model of the web
shows it, and a solid mesh does if it is fine enough to bend.

**Elastic and linear unless asked otherwise.** The static, modal and
buckling analyses and the topology optimisation are linear and elastic:
small strain, small displacement, one factorisation - a material's
`plasticity` block is ignored there (the solve and optimise apps say so). The
deformation figures of those runs are exaggerated by a stated factor purely
for visibility. `sparlab_solve` adds a non-linear static analysis
(`nonlinear` in the deck; `docs/formulation.md`, sections 7c and 7d): large
displacement and rotation, with the Saint Venant-Kirchhoff law (small strain)
or a compressible neo-Hookean law (large strain), follower pressures, a
centrifugal load at the deformed position, the finite-strain thermal split,
J2 plasticity with isotropic (linear and Voce) and kinematic (Prager)
hardening in small strain or in the Green-Lagrange strain at large rotation,
load paths that unload and reverse, load control that stops at limit and
bifurcation points and at a plastic collapse and says so, and the arc-length
method through limit points. What it does not do:

* plasticity is rate-independent J2 at small strain: no finite-strain
  plasticity (with `finite` kinematics the J2 return in the Green-Lagrange
  strain is sound for large rotation with small strains, and the run warns
  beyond a strain of 0.05), no creep, viscoplasticity, damage, fracture,
  non-associative or pressure-dependent yield, nonlinear (Armstrong-Frederick)
  kinematic hardening, anisotropy or temperature-dependent properties; the
  elastic part of the neo-Hookean law cannot be combined with plasticity;
* constant-strain elements (Tri3, Tet4) can lock under the isochoric flow of
  a fully plastic state in plane strain and 3-D, depending on the mesh
  pattern, and fully integrated Q4 and Hex8 lock without the mean
  dilatation that is their default - a locked collapse load comes out high,
  which is unconservative;
* the static path is quasi-static, and a snap-through that the arc-length
  method follows is a sequence of equilibria, not the dynamic jump a real
  structure would make - the non-linear transient (below) integrates that
  jump, with its own limits; contact has limits of its own (below);
* at a bifurcation of a perfect structure there is no branch switching: load
  control stops there, and the arc-length method stays on the fundamental
  path. A post-buckling analysis needs an imperfection built into the mesh
  (a perturbed geometry), which the deck does not generate;
* the Saint Venant-Kirchhoff law is valid for small strain only - the summary
  warns above a Green-Lagrange strain of 0.05, and under strong compression
  (a stretch below `1/sqrt(3)`) it softens unphysically. The neo-Hookean law
  has no plane-stress form and no thermal strain; its constants are the small-
  strain `E` and `nu`, not a fit to test data of a real elastomer;
* a follower pressure over free edges makes the tangent non-symmetric; LU
  factorises it, which reports no inertia, so the stability of such a state
  is not assessed (the summary says so);
* body forces and self-weight stay dead loads per reference volume, a
  prescribed displacement scales with the load factor, and every load of a
  case shares one load factor - a `load_path` scales them all together
  (load, unload, reverse), but there are no load sequences of different
  loads (a preload followed by a service load) within a case;
* the topology optimisation, the modal analysis and the buckling analysis do
  not use the non-linear state (`sparlab_topopt` refuses the block).

**Buckling is linear bifurcation.** The buckling check is the eigenvalue
problem `(K + lambda K_G(u)) phi = 0` of the linear static state: the
bifurcation load of the perfect geometry, which is an upper bound on the
collapse load of a real part with its imperfections, residual stresses and
plasticity. No imperfection sensitivity, no follower loads, no knock-down
factor. A load factor of 6 is not a safety factor of 6 against collapse; for
a shell-like or imperfection-sensitive structure the difference can be
large. The non-linear analysis gives the limit load of the geometry as
meshed - an imperfect mesh included - and brackets the first loss of
stability under load control, but not a post-buckling branch that starts at
a bifurcation of the perfect geometry.

**Dynamics: direct integration with a constant step, and the direct harmonic
response.** The transient analysis integrates `M a + C v + K u = A(t) f` (or
its non-linear counterpart) by HHT-alpha with a constant step, and the
harmonic response solves the complex dynamic stiffness at every frequency.
What they do not do:

* one amplitude scales every load and prescribed displacement of a case: no
  loads that follow different time histories within a case, no initial
  velocity field, and no load applied at a moving position;
* the step is constant and chosen by the user: no automatic step control by
  an error estimate (a failing Newton step is halved, and the halves grow
  back), and no explicit (central-difference) integration for short
  impact-like events, which the implicit method can only follow with a small
  step at a high cost per step;
* damping is Rayleigh's (and, for the harmonic response, a structural loss
  factor): no modal damping ratios, no discrete dampers, no frequency-dependent
  or viscoelastic material, and in a non-linear run the Rayleigh stiffness
  term uses the linear elastic stiffness, not the tangent;
* the non-linear transient starts at rest (a preloaded start is a static
  analysis of its own) and has the non-linear statics' limits (section
  above); the trapezoidal rule conserves the energy of a non-linear system
  only to `O(dt^2)`, and the run warns when the balance shows energy created
  beyond 1 % - it does not bound the error of the motion itself;
* no modal superposition (neither for the transient nor for the harmonic
  response), no response spectrum, random vibration or power spectral
  densities, no rotor dynamics (Coriolis and gyroscopic terms), no
  fluid-structure coupling, no wave-absorbing boundaries, and no fatigue;
* the harmonic response is linear: the load case's linear stiffness, no
  prestress from a static preload (the geometric stiffness of a preloaded
  structure is not added), and no harmonic balance for a non-linear
  structure;
* a natural frequency of the modal analysis is the undamped eigenvalue of the
  constrained model; the dynamic analyses are not coupled to the topology
  optimisation (`sparlab_topopt` refuses their blocks).

**Contact: small sliding, in statics.** The non-linear static analysis
models unilateral contact between a surface and a rigid plane, cylinder or
sphere, or between two surfaces of the model (dual mortar), frictionless or
with Coulomb friction (`docs/formulation.md`, section 7f). What it does not
do:

* the contact geometry is that of the reference configuration and the gap is
  linear in the displacement: small displacements and small sliding. Surfaces
  that slide by more than a fraction of an element, rotate, or come into
  contact after a large motion are outside the model, and contact refuses
  `finite` kinematics and the arc-length method (load control only, so a
  contact problem with a limit point cannot be followed past it);
* linear elements only (Q4, Tri3, Hex8, Tet4): the dual basis does not exist
  on a Tet10 face, so a quadratic mesh has no contact;
* statics only: the linear static, modal, buckling, transient and
  frequency-response analyses ignore contact - no impact, no contact in a
  transient, no modes or buckling of a model held by its contact - and the
  topology optimisation does not take it;
* a body that only its contact holds must touch its support at the start. A
  gap under a load that nothing else resists is not closed (the stiffness is
  singular until it closes); a gap is closed by prescribed displacements or
  a moving obstacle, and a body held only by frictionless contact that can
  slide stops with the reason. There is no stabilisation (no artificial
  damping or soft springs) to carry such a body into contact;
* friction is isotropic Coulomb with one constant coefficient: no static and
  kinetic coefficients, no dependence on velocity, pressure or temperature,
  no anisotropy. The slip is measured over each load step, so with friction
  the load steps are part of the answer, as they are physically - a
  frictional path is not reversible;
* the constraint is enforced on the slave side. A curved master surface
  meshed coarser than its slave adds the error of its chords, first order in
  its element size (a centre-pressure error of `3.6e-2` at `a / h = 10.6`
  for a master twice as coarse, in the Hertz study) - mesh the master as finely
  as the slave or finer. Slave nodes at the edge of the master surface, which
  it covers only in part, are left out of contact with a warning, as are
  nodes with nothing opposite them within the search distance;
* rigid obstacles are analytic planes, cylinders and spheres that translate
  with the load factor: no rotation, no other shapes and no meshed rigid
  bodies. Against a curved obstacle the gap is taken at the nodes (exact on
  a plane), and its force acts along its own normal: the Hertz study shows
  the `a / R` difference to the half-space theory this makes;
* no automatic contact search: every pair is declared, a node takes one
  constraint (it cannot be on two slave surfaces, or on a slave and a master
  one), and there is no self-contact detection. No tied (bonded) interfaces,
  adhesion, cohesive zones, wear or thermal contact conductance, and no
  contact of shells or beams (both refused);
* two mortar surfaces that start a gap apart transmit friction across it, as
  the small-sliding model on the reference geometry does in every code: the
  couple of that force pair (the gap times the tangential force) remains in
  the moment balance.

**Shells: MITC4, linear, on four-node quadrilaterals.** The shell is the
degenerated continuum with Dvorkin and Bathe's transverse-shear tying and a
drilling penalty (`docs/formulation.md`, section 7g), verified against
exact plate, shell and ring solutions and the MacNeal-Harder benchmarks
(`docs/verification.md`, section 27). What it does not do:

* linear only: static, modal, harmonic response and linear buckling. The
  non-linear analysis (large rotations, plasticity), contact and the
  transient analysis refuse a shell model, the last because the rotation
  about the director moves no mass - the mass matrix is singular - so the
  initial accelerations are not determined. Temperatures and the
  centrifugal load are refused too. Topology optimisation takes a shell -
  each cell's density scales its whole stiffness, a perforated sheet of fixed
  thickness, not a thickness distribution - with the buckling constraint on
  its out-of-plane buckling, but not the stress constraint (its stress varies
  through the thickness), the overhang filter or the manufacturability
  checks;
* one element: the four-node quadrilateral. No triangles (S3), no quadratic
  shells (S8R, S9) - a mesh file with them is refused with the reason - so
  a curved surface is a mesh of flat-ish facets, and the geometry error of
  those facets is part of the discretisation error (the thick-ring study
  converges to the circle at second order);
* MITC4 is free of shear locking on meshes of parallelograms, not on
  distorted ones: a distorted 4 x 4 mesh of a thin clamped plate reaches
  0.15 of the deflection at `t / a = 1e-3`, while from 8 x 8 on the
  deflection no longer falls with the thickness. It is not free of membrane
  locking: the pinched cylinder reaches 0.38, 0.75, 0.93 and 0.99 of its
  reference on 4 to 32 cells a side. Coarse meshes of bending-dominated
  curved shells are too stiff;
* the rotation about the director has no stiffness of its own: a penalty
  `alpha G t` (`alpha = 1e-3`) ties it to the in-plane rotation of the
  displacement field. A moment about a shell's normal therefore acts only
  through the penalty, and at a fold, where one wall's drilling rotation is
  its neighbour's bending rotation, the penalty couples the two walls. The
  box beam moves by at most `1.7e-4` for `alpha` from `1e-6` to `1e-2`; a
  penalty of `1e-1` begins to stiffen the walls;
* one isotropic material per element and a thickness constant over each
  element (`model.shell.sections`): no composite layups, no offsets of the
  reference surface from the mid-surface, no thickness varying within an
  element, no stress through the thickness (`sigma_33 = 0`), and the
  transverse shear stress is recovered as the parabolic `1.5 Q / t` at the
  mid-surface. Stresses and resultants are those at each element's centre;
* the directors are the mesh's exact normals when it has them (the
  generated plate, cylinder and sphere), or averaged over the elements that
  meet at a node within the fold angle (20 degrees by default): a coarse
  mesh of a curved surface read from a file has directors that are not its
  true normals. Self-weight and mass are those of the solid the directors
  sweep, which on the flat facets of a curved surface falls short of `t A`
  by a discretisation error of order `h^2` (1.8 % on 5 x 4 cells of a
  sphere zone, 0.45 % on 10 x 8);
* a point load on a shell is singular twice over: its bending stress and,
  in a shear-deformable shell, its deflection, whose transverse-shear part
  grows like `log(1 / h)` under the load. The pinched cylinder passes its
  thin-shell reference on the finest mesh for that reason; spread a load
  over an area to get a deflection that converges;
* shells cannot be tied to solids or beams (no shell-to-solid coupling, no
  rigid links or multi-point constraints), and a model is all shells or
  none;
* the frequency response's field files carry the translations only, and a
  monitor reads a translation component;
* CalculiX's `S4`, which the exporter writes, is not the same element: it
  expands each shell into incompatible-mode solids over normals it averages
  itself, and imposes held rotations and nodal moments through rigid knots,
  which stiffen its model where rotations are held (tenfold on a quarter of
  the Scordelis-Lo roof). The export is a starting point for a CalculiX
  model, not the same discrete problem.

**Beams: two-node Timoshenko, linear, Saint-Venant torsion.** The beam
element uses the interdependent interpolation, exact for nodal and uniform
loads on straight members and free of shear locking (`docs/formulation.md`,
section 7h), verified against the exact solutions of its model
(`docs/verification.md`, section 28). What it does not do:

* linear only: static, modal, harmonic response, transient and linear
  buckling. The non-linear analysis (large rotations, plasticity), contact,
  temperatures, the centrifugal load, tractions and pressures (a beam has no
  faces) and topology optimisation are refused, each with its reason;
* torsion is Saint-Venant's: the section warps freely, so an open
  thin-walled section (a channel, an I) twists too easily near a restrained
  end, where its warping stiffness carries part of the torque; the torsion
  constant of a rectangle comes from its series, of a general section from
  the deck;
* the geometric stiffness is that of the axial force alone
  (`N int [u'^2 + v'^2 + w'^2 + (I/A) theta'^2]`, the element's mean axial
  force): the bending stresses' part is left out, so there is no
  lateral-torsional buckling, and every twist buckles at `G J A / I_p`
  (without warping stiffness, the classical torsional buckling load);
* the section's axes are principal and its centroid, shear centre and
  reference axis coincide: no products of inertia, no offsets or
  eccentric connections, no end releases (hinges) - every joint is rigid -
  no tapered members, and a curved member is a polygon of chords (the
  quarter-circle study converges to the arc at second order);
* shapes are the rectangle, the solid circle and the circular tube, with
  Cowper's shear coefficients; anything else is a general section whose
  properties the deck states;
* the normal stress is `|N|/A + |M_y| c_z / I_y + |M_z| c_y / I_z` at the
  extreme fibres of each element end - exact for a rectangle, a bound for a
  section inside that box - or with `sqrt(M_y^2 + M_z^2)` for a round
  section; shear stresses are not reported, and under a uniform load the
  moment's extreme may lie inside an element, not at its ends;
* frequencies and buckling loads converge at `O(h^2)` once the elements
  are shorter than the section is deep: model a stocky member with as many
  elements as a slender one needs for the same accuracy;
* CalculiX's `B31`, which the exporter writes for rectangular sections, is
  a different model: it expands each element into bricks over its
  rectangle and joins members at a node through a rigid knot. Its
  displacements come out smaller than the centre-line beam's, and more so
  once members meet: 0.5 % on a lone clamped column, 3.2 % on a portal of
  two columns and a beam (1.8 % against the Euler-Bernoulli limit of the
  same portal, so not the shear deformation), 5.4 % on the space frame of
  the cross-validation. Its `U1` beam is
  the Timoshenko beam in statics, but version 2.21's shear term stiffens it
  (the deflection falls below Euler-Bernoulli's as the shear coefficient
  falls), so it is compared in the Euler-Bernoulli limit only.

**Body loads: self-weight, force densities and steady rotation.** Gravity,
uniform body force densities on element regions and the centrifugal load of a
steady rotation are integrated exactly from the consistent mass. Rotation is
the static centrifugal load only: no Coriolis or gyroscopic terms. In the
linear analyses it acts at the undeformed position; the non-linear analysis
applies it at the deformed position, with the spin-softening stiffness that
follows. There is no inertia relief (a
free-flying body balanced by its own acceleration), so a body load must be
reacted by supports. The wing-rib deck's "fuel inertia" case is still a
*representative edge pressure*, not a body-force calculation.

**Thermal strain is linear thermoelasticity.** A temperature field - uniform,
given on regions, or the solution of steady conduction with fixed
temperatures, surface fluxes, convection and heat generation - enters as the
free strain `alpha (T - T_ref)` with temperature-independent `E`, `nu`,
`alpha` and `k`. The conduction is steady (no transient heat transfer, heat
capacity or time), has no radiation, no contact conductance between parts and
no heat from deformation; the temperature is solved first and does not depend
on the displacement. Materials may differ by element region
(`material_regions`), each with its own reference temperature, but the
CalculiX export needs one common reference temperature (CalculiX measures
thermal strain from the initial nodal temperature). First-order hexahedra and
quadrilaterals integrate the interpolated temperature (the consistent load),
which cannot represent the free expansion of a linear temperature gradient
exactly - that needs a quadratic displacement - so a free Hex8 or Q4 part in
a gradient shows a small spurious stress that vanishes under refinement; the
Tet10 represents it exactly. In the non-linear analysis the thermal stretch
`1 + alpha dT` splits off multiplicatively (Saint Venant-Kirchhoff only),
with the same temperature-independent constants. Topology optimisation takes
body loads and uniform or regional temperatures (next section), but not a
conducted temperature field, whose conduction path the design would change,
nor several materials (`sparlab_topopt` refuses both).

**Linear elements, and one quadratic element.** The four-node
quadrilateral, the three-node triangle, the eight-node hexahedron and the
four-node tetrahedron are stiff in bending and need several elements through
a bending depth; the constant-strain simplices are much stiffer than the
bilinear and trilinear elements. The mesh-convergence studies quantify it:
the coarsest solid cantilever is 16 % too stiff on Hex8 and 50 % on Tet4. The
ten-node tetrahedron removes most of this for solids: on a 10 x 2 x 1-cell
cantilever, where the Hex8 is 33 % and the Tet4 66 % too stiff, it is within
0.06 % of beam theory, and on the engine mount meshed from CAD the Tet4 compliance
is 19 % below the finest Tet10 run on the benchmark's own mesh
(`docs/benchmarks.md`, section 11). There is no quadratic plane element
(Q8, Tri6) and no Hex20, so a 2-D model or a hexahedral mesh still needs
refinement to bend. On a Tet10 cell with curved edges the 4-point stiffness
rule is not exact - the integrand is rational - as in every code that uses
C3D10; the study's curved and straight-sided meshes of the same cells agree
to 0.13 % at 3 mm. Enhanced assumed strain or B-bar would help the linear
elements per DOF; `docs/architecture.md` says what adding an element
involves.

**Meshes: one cell type per mesh.** The structured generators make boxes of
quadrilaterals, triangles, hexahedra or tetrahedra, linear or (tetrahedra,
`mesh.order: 2`) quadratic. Real geometry comes in through Gmsh (MSH 2.2 /
4.1, ASCII) and Abaqus / CalculiX `.inp` files, with these limits: one cell
type per mesh (a quad-dominant or a hex-dominant mesh with prisms is
refused, with the re-export that fixes it); the only second-order cell is
the ten-node tetrahedron (C3D10, Gmsh type 11), whose edge nodes may follow
a curved surface - every other boundary is a polygon of cell edges, and
elevating a Tet4 mesh keeps its facets;
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

**Plane strain is verified on curved sections, not on beams.** The thick
cylinder under pressure, the conducted thermal cylinder, the rotating
cylinder and the finite-strain tube run in plane strain against exact
solutions, the Hex8 and Tet10 sections held at `u_z = 0` reproduce the Q4
answers, and CalculiX's plane-strain expansion (`CPE4`) agrees node by node.
The beam and frequency studies against Timoshenko and Euler-Bernoulli run
in plane stress only.

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

A real joint analysis needs contact at the bolt and bearing faces - which
the non-linear analysis now provides for small sliding, and the benchmark
decks do not use - and a fastener model. Nothing here should be read as a
joint substantiation.

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

**Manufacturability: two rules modelled, nothing more.** The overhang
filter models one additive-manufacturing rule - no material without support
directly or diagonally below it, a 45-degree limit on square or cubic cells -
and the robust formulation models one uniform manufacturing error, the whole
part coming out thinner or thicker by the erosion. A design that passes the
overhang check is printable *under that rule*: support removal, residual
stress and distortion, surface finish, anisotropic material, a minimum wall
and the build plate's own limits are not modelled. Nothing models a draw
direction, milling-tool access, fillet radii or a machining setup. The
filter radius and the robust formulation impose a minimum *length scale*,
which is a geometric property, not a process. Without these options the
designs are load-path guidance for a designer, not parts; with them they
are still not certified parts.

**Loads that follow the design are interpolated, and the compliance still
decides.** Self-weight, body forces and the centrifugal load scale with
`gamma(rho)`, which is the element's volume fraction at and above
`topology.simp.body_load_threshold` (0.1) and falls like `rho^p` below it, so
that near-void material cannot sag without bound under its own weight; below
the threshold the load is therefore smaller than a graded material's would
be. A 0/1 design, and every element at or above the threshold, carries its
physical load. Thermal loads scale with the stiffness factor `E(rho)/E_0` -
the SIMP material expands with the solid's coefficient - and the temperature
field is fixed: a conducted field is refused, because the design would
change it. The objective is still the compliance `f^T u` of the whole load.
Under a temperature field the thermal compliance of a graded design falls
with its stiffness, so the minimum-compliance design can use less than its
volume allowance and stay grey, down to no element at the interpretation
threshold and so no part to export. The run warns of both; the projection sharpens such a design, and a stress constraint bounds
what the heating does to the material that remains, but compliance is not a
strength criterion for a thermal load. The body-load interpolation and the
gradients are verified (`docs/verification.md`, section 29); no optimised
design under these loads is compared with another code or an experiment.

**The optimum is local.** SIMP with `p > 1` is non-convex; a different starting
design, penalty schedule or mesh can converge to a different local optimum. The
design study shows exactly that: the mesh-refinement arm with the filter radius
fixed in cells produces visibly different topologies at different resolutions.
Continuation reduces the dependence but does not remove it.

**The projection sharpens the design; alone it does not give it a length
scale.** A density filter alone leaves an intermediate band about one radius
wide at every material boundary (the MBB benchmark settles at a grey level
of 0.28, the solid bracket at 0.34). The Heaviside projection removes most
of it - the projected solid bracket ends at 0.013 - and with it the gap
between the objective and the thresholded structure. A single projection at
`eta = 0.5` does not secure a minimum member size: features thinner than the
filter radius survive, and on the projected MBB beam two diagonals vanish
under a 0.1 erosion and cost 22 % of the stiffness. The robust formulation
(`topology.projection.robust`) is the tool for one, and the length-scale
check measures what a design delivers. The size follows from the filter
radius and `robust_delta`; it is not a stated number the solver enforces.
At `robust_delta = 0.1` the MBB's smallest member grew only from 3 to 4
cells, and at 0.05 the robust column still has members about one cell
thick. At a large `beta` the optimisation also becomes
strongly non-convex: single elements at the solid-void interface flip by the
full move limit while the compliance stands still, which is why the
projected decks stop on the objective-stall criterion rather than on the
design change. The `beta` schedule is a further setting the result depends
on, and it is recorded with every run.

**Three constraint types.** Optimality criteria handles the volume
constraint alone (of the dilated design in a robust run); MMA handles the
volume plus one aggregated stress constraint and one aggregated buckling
constraint per load case. Per-mode frequency constraints, displacement
bounds and multiple volume budgets on separate regions are not written -
`docs/architecture.md` says what each would take, and `StressConstraint` and
`BucklingConstraint` are the templates - and the OC path remains the default
because it converges with no tuning where MMA needs a move limit chosen per
problem.

**The buckling constraint acts on the SIMP model, not on the part.** The
constrained load factors are those of the density field, in which
intermediate densities keep `rho^p` of their stress stiffness. The exported
part - thresholded, full material, grey members gone - can buckle earlier:
on the column benchmark the constraint holds `lambda >= 6` on the SIMP model
while the part buckles at 3.65 with a plain projection and 5.81 with the
robust formulation. The buckling check of the exported part is therefore
the number to judge, and a design should be re-run with a margin when it
falls short. The KS aggregate overestimates the largest ratio by at most
`ln(m)/P` and is conservative; the sensitivity formula holds for simple
eigenvalues, and at a repeated eigenvalue only the aggregate stays
differentiable. Void material can still produce pseudo modes; the
energy-fraction diagnostic names them but does not remove them. The cost is
one eigensolve and one adjoint per aggregated mode every iteration - 242
linear solves per iteration on the column against one without the
constraint - and the benchmark needed a move limit of 0.05 and 3000 Newton
iterations for MMA's subproblem.

**The overhang filter needs a structured grid, and OC struggles with it.**
The supports of an element are grid cells of the layer below, so the filter
and its check run only on `structured_quad` and `structured_hex` meshes,
with the build direction along a grid axis; the overhang angle is
`atan(layer thickness / cell width)`, 45 degrees on square cells. A
45-degree staircase touches its supports along an edge, which the
interpretation's face connectivity does not count as a joint, so a filtered
3-D design can leave small islands (0.7 % of the bracket's material). With
OC, the filtered MBB and bracket runs locked into period-2 cycles at the
final `beta`; both overhang decks use MMA.

**The robust formulation needs care with OC.** Each iteration projects the
filtered field twice more, for the dilated and the blueprint designs,
without an extra solve: only the eroded design is analysed. The dilated
design's volume target is rescaled from the blueprint's, so the blueprint
meets the volume fraction only as closely as the ratio of the two volumes
holds between rescalings. On the robust MBB, which rescales every 20
iterations, interface elements flip at `beta = 32`, the blueprint's fraction
alternates between about 0.498 and 0.502, and the returned design is 0.41 %
below it. Rescaled at every iteration, OC locked into a period-2 cycle at
`beta = 32` that met no stopping criterion.

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

**The non-linear analysis factorises its tangents.** Without contact every
Newton step factorises the tangent directly - LDL^T, whose inertia the
stability checks need, or LU for a non-symmetric one - whatever
`solver.linear` says, so a large solid non-linear model meets the direct
solver's limit above. With contact the symmetric steps follow
`solver.linear` - multigrid CG above the `auto` limits, which carries the
frictionless 3-D Hertz study to 62 544 unknowns, with LDL^T taking over
where it does not converge - but a step in which a node slips under friction
is factorised by sparse LU.

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
the smallest pivot relative to its diagonal entry stays many decades above
the `1e-14` threshold (it is logged at debug verbosity) - the singularity
check judges each pivot against its own diagonal, so neither the floor nor a
shell's rotations beside its translations trip it - but the condition number
of `K_ff` is correspondingly large. The multigrid solver copes because void and solid do
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
verified against exact answers on all seven element types. The patch tests
pass to 4e-15 (Q4), 2e-15 (Hex8), 9e-15 (Tri3 and Tet4) and 1.3e-12 (the
MITC4 shell, membrane and bending states on distorted meshes in a turned
plane), and the Tet10 passes the quadratic (pure-bending) patch test to
1.2e-14. The shell converges at second order to the exact Reissner-Mindlin
deflection, frequencies, harmonic response and buckling loads of a plate,
without shear locking from `t / a = 1e-1` to `1e-4`, and to the exact state
of a thick ring. The Timoshenko beam is exact for end loads and uniform loads
on straight members, one element per member included (2.2e-12), and
converges at second order to the exact frequencies, harmonic response and
buckling loads of its model and, through its chords, to a curved member. The
sparse and dense solvers agree to 9e-12, and multigrid CG agrees with Cholesky to
4e-12. The compliance gradient matches central differences to 2e-8 on Q4
and Hex8, and to 1e-7 through the Heaviside projection on Q4 and Tet4; the
buckling constraint's gradient to 1.7e-6 and the overhang filter's to
2.8e-6. Mass is conserved to 5e-14. A homogeneous large deformation on
distorted meshes is reproduced to 1e-9 m, and the tube at finite strain
converges at every element's order to its exact solution. A transient run
equals the exact solution of its discrete equations (every mode integrated
exactly, in extended precision) to 3.3e-10, a rod's harmonic response its
exact discrete solution to 2.9e-10, and the rod's harmonic and transient
responses and a one-element oscillator's finite-strain elastic and
elastoplastic motion converge at second order to their exact solutions.
Linear static displacements are cross-validated node by node against two
independent codes on thirty-eight problems with sixty-three load cases,
covering all five continuum element types and both mesh-file formats; so
are buckling load factors on three columns, large-deflection states on four
decks, elastoplastic states on ten, transient histories on five, harmonic
responses on two and contact states on six; and the shell's displacements,
rotations, frequencies and buckling factors on four decks against an
independent MITC4 written in NumPy, to 6.2e-9 or better but for the
rotations of the worst-conditioned deck, 2.6e-7, below its round-off
scale; and the beam's displacements, rotations, end resultants,
frequencies, buckling factors and harmonic response on three decks against
an independent Timoshenko frame written in NumPy, to 1.3e-10 or better.
scikit-fem agrees to solver round-off, displacements and
load factors alike: every linear difference lies below the round-off scale
of its system (its condition number times eps), 1.5e-10 or better except on
the three worst-conditioned systems - a plane-strain strip and two slender
Tet10 cantilevers, 1.4e-9 to 2.2e-8; scikit-fem's own total Lagrangian
solver agrees with SparLab's non-linear states to 1e-14, and its own J2
solver with the elastoplastic states, small strain and finite, E-bar
included, to 7e-12; an HHT-alpha integration and a direct complex solve on
its matrices agree with the transient histories to 2.2e-8 (the final
acceleration of an elastoplastic run; its displacements to 1.2e-12) and
with the harmonic responses to 1.0e-10, and an independent contact solve
with the contact states to 8.6e-13. CalculiX's displacements agree to
the rounding of its own result file, 4.3e-6 or better, linear, `NLGEOM`,
`*PLASTIC` and `*DYNAMIC` alike, wherever the two codes solve the same
discrete problem; its plane-stress
comparisons at `nu != 0` and its finite-strain plasticity are recorded but
not judged, its kinematic hardening is not compared (it does not
reproduce Prager's rule), and neither are its plane elements in dynamics
(their `*DYNAMIC` response contradicts CalculiX's own `*FREQUENCY`); its
buckling factors differ by up to 8.3e-5 for a reason not identified, and
its shell (`S4`, solids over its own normals, rotations held through rigid
knots) and its linear beam (`B31`, a brick over the rectangle) are
different discretisations, recorded and not judged; its `U1` beam, in its
Euler-Bernoulli limit, agrees with the Euler-Bernoulli frame to 9.1e-8, at
the seven digits it prints.
Validation against independent theory covers exactly these references:
Euler-Bernoulli and Timoshenko cantilever deflection, Euler-Bernoulli
bending frequencies, fixed-free rod axial frequencies, the Euler-Engesser
buckling load of a clamped column, Euler's elastica, Hertz's line and point
contact, the thin-shell references of three MacNeal-Harder shell
benchmarks, and beam theory and Bredt's torsion for a box beam. Stresses,
the continuum elements' frequencies, non-linear static load paths and
optimised designs are not compared with another code, and there is **no
comparison against experiment**.

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
