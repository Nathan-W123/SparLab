# Mathematical formulation

This document states exactly what SparLab solves. Notation follows
`docs/conventions.md`; every equation below is implemented in the file named
beside it.

## 1. Continuum problem

Two-dimensional small-strain linear elasticity on a domain `Omega` of constant
out-of-plane thickness `t`, with displacement boundary `Gamma_u` and traction
boundary `Gamma_t`:

```
  div sigma + b = 0            in Omega          (equilibrium)
  eps = 1/2 (grad u + grad u^T)                  (kinematics)
  sigma = D : eps                                (constitutive)
  u = u_bar                    on Gamma_u
  sigma . n = t_bar            on Gamma_t
```

Body forces `b` are not implemented; loads enter as boundary tractions and
concentrated nodal forces. The weak form, for all admissible `v`, is

```
  integral_Omega  t eps(v)^T D eps(u) dOmega
      =  integral_Gamma_t  t v^T t_bar dGamma  +  sum_k v(x_k)^T F_k
```

with the second term the point loads.

### Plane idealisations

With Voigt ordering `{sigma_xx, sigma_yy, sigma_xy}` and
`{eps_xx, eps_yy, gamma_xy}` (`gamma_xy = 2 eps_xy`), an isotropic material
gives (`material/IsotropicMaterial.cpp`):

**Plane stress** (`sigma_zz = 0`; a thin sheet loaded in its plane)

```
           E     | 1    nu       0     |
  D  =  -------- | nu    1       0     |
        1 - nu^2 | 0     0  (1-nu)/2   |
```

**Plane strain** (`eps_zz = 0`; a long prismatic body, or a plane of symmetry)

```
                 E          | 1-nu   nu        0      |
  D  =  ------------------   | nu    1-nu       0      |
        (1+nu)(1-2nu)        | 0      0    (1-2nu)/2   |
```

In plane strain the out-of-plane stress is `sigma_zz = nu (sigma_xx + sigma_yy)`
and the von Mises formula must use it; `fem/StressRecovery.cpp` does.

The `(3,3)` entry of both matrices equals the shear modulus
`G = E / (2(1+nu))`, which is a consequence of pairing Voigt stress with
*engineering* shear strain and is checked in the test suite.

Both idealisations are implemented and unit-tested. Plane stress is the default
and is the path the benchmark and verification studies exercise.

## 2. Element formulation

### Shape functions

The four-node bilinear isoparametric quadrilateral (`elements/Quad4.cpp`), on
the reference square `(xi, eta) in [-1,1]^2`:

```
  N_1 = 1/4 (1 - xi)(1 - eta)        N_3 = 1/4 (1 + xi)(1 + eta)
  N_2 = 1/4 (1 + xi)(1 - eta)        N_4 = 1/4 (1 - xi)(1 + eta)
```

with natural derivatives

```
  dN_1/dxi  = -1/4 (1 - eta)     dN_1/deta = -1/4 (1 - xi)
  dN_2/dxi  = +1/4 (1 - eta)     dN_2/deta = -1/4 (1 + xi)
  dN_3/dxi  = +1/4 (1 + eta)     dN_3/deta = +1/4 (1 + xi)
  dN_4/dxi  = -1/4 (1 + eta)     dN_4/deta = +1/4 (1 - xi)
```

These satisfy `sum_a N_a = 1` everywhere, `N_a = delta_ab` at node `b`, and
`sum_a grad N_a = 0` - all three are unit-tested.

### Isoparametric mapping

Geometry and displacement use the same basis:

```
  x(xi, eta) = sum_a N_a(xi, eta) x_a        u(xi, eta) = sum_a N_a(xi, eta) u_a
```

The Jacobian and its determinant are

```
  J_ij = dx_i / dxi_j = sum_a x_a,i  dN_a/dxi_j            (a 2 x 2 matrix)
  detJ = J_11 J_22 - J_12 J_21
```

`detJ <= 0` at any quadrature point means the element is inverted, collapsed or
wound clockwise, and raises `MeshError` naming the element and the parametric
point. Physical shape-function gradients follow from
`dN_a/dx = (dN_a/dxi) J^{-1}` (the 2x2 inverse is written out explicitly; the
assembly loop runs it millions of times).

### Strain-displacement operator

```
          | dN_1/dx    0     dN_2/dx    0     ...  |
  B  =    |   0     dN_1/dy     0    dN_2/dy  ...  |      (3 x 8)
          | dN_1/dy dN_1/dx  dN_2/dy dN_2/dx  ...  |
```

so `eps = B u_e` with `u_e` in node-major order. Because the Q4 space contains
all linear fields exactly, `B u_e` reproduces a constant strain field exactly at
*every* parametric point, on distorted meshes as well as uniform ones. That is
what the patch test checks, and it holds to round-off.

### Quadrature

Gauss-Legendre tensor-product rules on the reference square
(`elements/Quadrature.cpp`). An `n`-point rule per direction integrates
polynomials of degree `2n - 1` exactly.

| Kernel | Default rule | Why |
|--------|--------------|-----|
| Stiffness `K_e` | 2 x 2 | exact for the Q4 on a parallelogram; full integration, no reduced-integration hourglass modes |
| Mass `M_e` | 3 x 3 | `N^T N` is biquadratic, and on a general quadrilateral `detJ` is not constant |
| Edge traction | 2 points | exact for a linear traction times a linear shape function |

All are configurable per run via `model.integration`. On a rectangle the 2x2,
3x3 and 4x4 stiffness matrices agree to round-off, which is unit-tested.

### Element matrices

```
  K_e = integral_Omega_e  t B^T D B  dOmega
      = sum_g  w_g  t  detJ(xi_g)  B(xi_g)^T D B(xi_g)              [N/m]

  M_e = integral_Omega_e  rho t N^T N dOmega
      = sum_g  w_g  rho t detJ(xi_g)  N(xi_g)^T N(xi_g)             [kg]
```

where `N` is the `2 x 8` shape-function matrix. Both are symmetrised after
integration, since any asymmetry is pure round-off.

`K_e` has exactly three zero eigenvalues - the two translations and the
infinitesimal rotation - and five strictly positive ones. Both parts of that
statement are unit-tested on a uniform and on a distorted element.

Row-sum lumping of `M_e` gives the optional diagonal mass matrix. For the Q4 the
row sums add up to the exact element mass, so total mass is conserved either
way.

### Consistent edge loads

For an edge with end nodes `a`, `b` and a constant traction `t_bar`,

```
  f_e = integral_Gamma_e  t N^T t_bar ds
```

parametrised by `s in [-1,1]` with the constant Jacobian `|x_b - x_a| / 2`. For
a constant traction on a straight edge this splits the resultant evenly between
the two nodes, and the total is exactly `t_bar * (edge length) * t`. The test
suite checks the total is mesh independent under refinement.

## 3. Assembly

`fem/Assembler.cpp` builds a triplet list and compresses it to CSC:

```
  K = sum_e  s_e  A_e^T K_e^0 A_e            M = sum_e  m_e  A_e^T M_e^0 A_e
```

`A_e` is the (implicit) DOF gather operator, and `s_e`, `m_e` are per-element
scale factors: both are 1 for a plain analysis, and the topology optimizer
supplies `s_e = E(rho_e)/E_0` and a mass interpolation factor `m_e`.

On a **uniform structured mesh** all cells are geometrically identical, so
`K_e^0` and `M_e^0` are integrated once and reused; this is the dominant saving
in an optimisation loop. The generic per-element path is always available, and
the test suite asserts the two produce bit-comparable matrices.

## 4. Boundary conditions and the reduced system

Dirichlet conditions are applied by **partitioning**, not by penalty terms
(`fem/DofManager.cpp`, `fem/StaticAnalysis.cpp`). Splitting the DOFs into free
`f` and prescribed `p`,

```
  | K_ff  K_fp | | u_f |     | f_f     |
  |            | |     |  =  |         |
  | K_pf  K_pp | | u_p |     | f_p + r |
```

the first block row gives the system actually solved,

```
  K_ff u_f = f_f - K_fp u_p
```

and the reactions come from the *full* residual

```
  r = K u - f
```

whose entries at prescribed DOFs are the support reactions and whose entries at
free DOFs are the (tiny) solver residual. This gives exact reactions with no
second assembly, and keeps `K_ff` symmetric positive definite - no artificial
large diagonal entries to wreck the conditioning or the eigenvalue spectrum.

Because `K t = 0` for a rigid translation `t`, summing the residual over all
DOFs in one direction gives `sum(r) + sum(f) = 0` identically. That is the
global force balance every static result reports, and it is a genuine check on
the assembly and the solve rather than an identity of the post-processing.

## 5. Linear solvers

| Name | Eigen type | Use |
|------|------------|-----|
| `simplicial_ldlt` (default) | `SimplicialLDLT` with AMD ordering | SPD systems; reports non-positive or tiny pivots as a singular system |
| `simplicial_llt` | `SimplicialLLT` with AMD | strictly SPD |
| `sparse_lu` | `SparseLU` with COLAMD | indefinite systems |
| `conjugate_gradient` | `ConjugateGradient`, Jacobi preconditioner | large systems, iterative |
| `dense_lu` | `PartialPivLU` | verification of small models; refuses above 4000 DOFs |

After every solve the scaled residual `||A x - b|| / max(||b||, tiny)` is
compared with `solver.linear.residual_tolerance`. Exceeding it raises
`SolverError` with the measured value; a non-finite solution does the same. The
LDL^T path additionally inspects its pivots: a non-positive pivot, or a
`min/max` pivot ratio below `pivot_tolerance`, is reported as a singular system
together with the three modelling causes that usually produce it.

## 6. Stress recovery

```
  eps_e = B(xi_g) u_e                  sigma_e = s_e D eps_e
```

Element values are the average over the stiffness quadrature points; nodal
values are area-weighted averages of the adjacent element values (simple
averaging, used for smooth contours only - not a superconvergent patch
recovery). With `s_e` applied the result is the *macroscopic* stress carried by
the element, and the unscaled solid-material stress is reported alongside it.

Derived quantities:

```
  sigma_vm = sqrt( 1/2 [ (sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 ] + 3 sxy^2 )
  szz = 0                            (plane stress)
  szz = nu (sxx + syy)               (plane strain)

  sigma_{1,2} = (sxx+syy)/2  +-  sqrt( ((sxx-syy)/2)^2 + sxy^2 )

  U_e = 1/2 s_e u_e^T K_e^0 u_e      (element strain energy, exact)
```

## 7. Modal analysis

Undamped free vibration on the free DOFs (`fem/ModalAnalysis.cpp`):

```
  K_ff phi = lambda M_ff phi ,    lambda = omega^2 ,    f = omega / (2 pi)
```

Prescribed DOFs are removed by the same partitioning as the static solve, so no
constrained DOF can contribute a spurious mode.

**Algorithm.** Bathe subspace iteration with optional shift-invert. With
`q = min(n, max(2m, m+8))` trial vectors,

```
  Xbar_{k+1} = (K_ff - sigma M_ff)^{-1} M_ff X_k
  K_r = Xbar^T K_ff Xbar ,   M_r = Xbar^T M_ff Xbar          (q x q, dense)
  solve  K_r Q = mu M_r Q ,   X_{k+1} = Xbar Q
```

One sparse Cholesky factorisation serves every iteration and every vector. The
starting basis is deterministic: the mass diagonal, then unit vectors at the
DOFs with the largest `m_ii / k_ii` ratio, then a seeded PRNG if more are
needed. Because `Q^T M_r Q = I`, the rotated basis is automatically
M-orthonormal, and the modes are renormalised at the end so `phi^T M phi = 1`
holds to round-off.

For systems with 400 or fewer free DOFs a dense
`GeneralizedSelfAdjointEigenSolver` is used instead: it is faster at that size
*and* gives the test suite an independent reference for the iterative path. The
two agree to 1e-8 relative on a 40x8 cantilever.

**Stopping rule.** Both the relative eigenvalue change and every eigenpair
residual `||K phi - lambda M phi|| / ||lambda M phi||` must be below their
tolerances. Including the residual matters because the highest requested mode
converges last: its eigenvalue can look settled while its vector is not.

**Validity screening.** Each converged pair is checked for a non-finite
eigenvalue (error), a negative eigenvalue beyond round-off (error - `K_ff` is
not positive definite, which is impossible for a properly constrained
linear-elastic model), and an eigenvalue indistinguishable from zero at the
problem's stiffness/mass scale (warning - an unsuppressed rigid-body or
mechanism mode).

**Mass conservation.** Summing every entry of the assembled mass matrix gives
`sum(M) = 2 * rho V`, one factor per translation direction. This identity is
exact for both the consistent and the lumped matrix and is used as a
verification check (measured error `<= 5e-14`).

### Analytical references

Two independent closed forms are used as validation references:

```
  bending, fixed-free beam:   f_n = (beta_n L)^2 / (2 pi L^2) sqrt(E h^2 / (12 rho))
                              beta_1 L = 1.87510407, beta_2 L = 4.69409113, ...

  axial, fixed-free rod:      f_n = (2n - 1) / (4 L) sqrt(E / rho)
```

The axial reference is *exact* for a bar in uniaxial stress, which plane stress
at `nu = 0` reproduces exactly; the measured agreement is 4e-6 relative. The
bending reference is not exact for a 2-D model, because Euler-Bernoulli theory
omits the shear deformation and rotary inertia the 2-D model includes - so the
computed frequencies must sit slightly *below* theory, and the gap must grow
with mode number. Both behaviours are asserted in the test suite, which is
stronger than asserting agreement.

## 8. Topology optimisation

See `docs/topology_optimization.md` for the SIMP interpolation, the filters, the
sensitivity derivation and the optimality-criteria update.
