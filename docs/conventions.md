# Conventions

Everything in SparLab is strict SI. There are no unit-conversion helpers and no
implicit scaling anywhere in the code, so a value read from a result file is
already in the unit its column header states.

## Units

| Quantity | Unit | Symbol in code / files |
|----------|------|------------------------|
| Length, coordinate, displacement | metre | `m` |
| Out-of-plane thickness | metre | `thickness_m` |
| Force, reaction | newton | `N` |
| Traction, stress, Young's modulus | pascal | `Pa` |
| Mass density | kilogram per cubic metre | `kg/m^3` |
| Mass | kilogram | `kg` |
| Area | square metre | `m^2` |
| Volume (area x thickness) | cubic metre | `m^3` |
| Strain, density, volume fraction | dimensionless | `[-]` |
| Compliance, strain energy | joule | `J` (= N m) |
| Eigenvalue `lambda` | 1 / second squared | `1/s^2` |
| Angular frequency `omega` | radian per second | `rad/s` |
| Frequency `f` | hertz | `Hz` |
| Stiffness entry | newton per metre | `N/m` |
| Time | second | `s` |

The MBB benchmark deliberately uses `E = 1 Pa`, a 1 m cell and `F = 1 N` so that
its compliance can be read against the dimensionless values usually quoted for
that problem. It is still SI - just a particular choice of scale.

## Coordinate system

Right-handed Cartesian, with the model in the x-y plane and z out of plane:

```
      y
      ^
      |            thickness t extends +- t/2 in z
      |
      +-------> x
     /
    z  (out of the page, towards the reader)
```

* The structured mesh generator places the domain with its lower-left corner at
  `(x0, y0)`, default `(0, 0)`.
* Node numbering on a structured grid is `node(i, j) = j * (nx + 1) + i`, so the
  x index runs fastest.
* Element numbering is `elem(i, j) = j * nx + i`, again x fastest. This is why a
  C-order reshape of an element vector to `(ny, nx)` produces an image whose
  row 0 is the *bottom* of the domain.
* Element nodes are stored counter-clockwise from the lower-left corner:

```
  3 --------- 2          eta
  |           |           ^
  |    (e)    |           |
  |           |           +---> xi
  0 --------- 1
```

Counter-clockwise ordering makes the Jacobian determinant positive; a clockwise
or self-intersecting element is rejected by `Mesh::validate()` with a message
naming the element.

## Degrees of freedom

Two translations per node, numbered node-major:

```
  dof(node n, component c) = 2 n + c,     c = 0 -> u_x,  c = 1 -> u_y
```

An element's DOF vector follows the same ordering:
`{u_x1, u_y1, u_x2, u_y2, u_x3, u_y3, u_x4, u_y4}`.

There are no rotational degrees of freedom: these are continuum elements, not
beam or shell elements.

## Signs

* Displacements and forces are positive along `+x` and `+y`. A downward load is
  a negative `y` force.
* Tractions are given as a stress vector in **global** components, not as
  pressure normal to a surface. A downward traction on a horizontal top edge is
  `[0, -p]`; to push *into* a surface you must supply the correct sign yourself.
  This is deliberate: a "pressure" convention needs an outward normal, which is
  ambiguous on a corner.
* Tensile stress is positive. `sigma_xx > 0` means the material is in tension
  along x.
* Reactions are reported as the force the *supports apply to the structure*, so
  applied load plus reaction sums to zero. Every static result records that sum
  and its relative error.
* Positive moment is counter-clockwise about `+z`, i.e. `M = x F_y - y F_x`,
  taken about the origin.

## Voigt notation

Stress and strain are stored as three-component vectors with the *engineering*
shear strain:

```
  sigma = { sigma_xx, sigma_yy, sigma_xy }
  eps   = { eps_xx,   eps_yy,   gamma_xy },    gamma_xy = 2 eps_xy
```

This pairing makes `sigma^T eps` the strain-energy density, which is why the
constitutive matrix carries `G` rather than `2G` in its (3,3) entry.

## Energy definitions

| Quantity | Definition | Note |
|----------|------------|------|
| Strain energy | `U = 1/2 u^T K u` | always the elastic energy stored |
| Compliance | `C = f^T u` | `f` is the **applied** load vector, excluding reactions |

For homogeneous Dirichlet data these satisfy `C = 2U` exactly, and every static
result reports the ratio `C / 2U` so the identity can be checked. With a
non-zero prescribed displacement the two differ: the patch test, for example,
has `C = 0` and `U > 0`, because all the work is done by the supports.

The topology objective is the weight-normalised sum of the per-case
compliances,

```
  c(x) = sum_l w_l f_l^T u_l ,     sum_l w_l = 1
```

so a run's objective value does not change if every weight is scaled by the same
factor.

## Density and volume fraction

* A design variable `x_e` and a physical density `rho_e` both live in `[0, 1]`
  and are dimensionless.
* The volume fraction is `sum_e rho_e v_e / sum_e v_e`, with `v_e` the element
  volume in `m^3`. It is computed on the **physical** (filtered) density, which
  is the quantity the constraint is written on.
* Passive solid regions count towards the volume budget. A mass constraint on a
  real part includes its attachment bosses, so the budget must too.
* Mass is `rho_material * sum_e rho_e v_e`, i.e. the SIMP density multiplies the
  material density linearly.

## Interpreting a density field

A SIMP density field is a *relative material distribution*, not a solid body.
SparLab never presents one as geometry without saying how it was converted:

* every figure of a density field states that black is `rho = 1`;
* the "interpretation" block in a topology summary records the threshold, how
  many elements survived it, how many disconnected groups they formed, and how
  much material was discarded as islands;
* the modal analysis of an "optimised structure" is run on that explicitly
  extracted sub-mesh, not on the density field.

See `docs/limitations.md` for what this interpretation does and does not claim.

## Tolerances

Every tolerance is configurable and every run records the value it used in
`summary.json` under `tolerances`. The defaults are:

| Tolerance | Default | What it bounds |
|-----------|---------|----------------|
| `solver.linear.residual_tolerance` | `1e-8` | scaled residual `\|\|Ku-f\|\| / \|\|f\|\|` after each solve |
| `solver.linear.pivot_tolerance` | `1e-14` | smallest / largest LDL^T pivot before the system is called singular |
| `solver.linear.iterative_tolerance` | `1e-12` | requested relative residual for the CG solver |
| `solver.equilibrium_tolerance` | `1e-6` | relative global force-balance error |
| `modal.tolerance` | `1e-10` | relative eigenvalue change between subspace iterations |
| `modal.residual_tolerance` | `1e-6` | relative eigenpair residual, also part of the stopping rule |
| `topology.optimizer.change_tolerance` | `1e-2` | `max \|dx\|` between iterations |
| `topology.optimizer.objective_tolerance` | `5e-5` | relative compliance change over `objective_window` iterations |
| `topology.optimizer.volume_tolerance` | `1e-10` | relative volume error of the multiplier bisection |

Exceeding a solver tolerance raises an exception with a diagnosis; exceeding an
optimiser tolerance is reported as non-convergence in both the console output
and `summary.json`. Nothing is silently accepted.

## Determinism

The only randomness in the code is the filler part of the subspace-iteration
starting basis and the node perturbation of `make_perturbed_quad_mesh`. Both
take an explicit seed (`modal.seed`, default `20240917`) and both default to
deterministic values, so repeated runs of the same deck on the same build give
bit-identical results.
