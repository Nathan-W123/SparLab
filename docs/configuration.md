# Input-deck reference

A run is fully determined by one JSON deck. Every numeric field is SI
(`docs/conventions.md`). Unknown keys are **reported**, and with
`--strict-config` they are an error - a misspelled key that silently takes its
default is one of the easiest ways to publish a wrong number.

Two conveniences beyond strict JSON: `//` and `/* */` comments, and a trailing
comma before `}` or `]`. Everything else follows RFC 8259, and a syntax error
reports its line and column.

## Top level

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | `"case"` | case name; also the default output directory |
| `description` | string | `""` | free text, echoed into `summary.json` |
| `units` | string | `"SI"` | informational only |
| `mesh` | object | required | see below |
| `material` | object | required | see below |
| `model` | object | `{}` | idealisation and integration |
| `boundary_conditions` | array | required, non-empty | displacement constraints |
| `load_cases` | array | required, non-empty | loading |
| `solver` | object | `{}` | linear solver and equilibrium tolerances |
| `modal` | object | `{}` | free-vibration analysis |
| `topology` | object | `{}` | topology optimisation |
| `output` | object | `{}` | which artefacts to write |

## `mesh`

```json
"mesh": { "type": "structured_quad", "nx": 240, "ny": 160,
          "lx": 0.30, "ly": 0.20, "x0": 0.0, "y0": 0.0 }
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `type` | string | `structured_quad` | the only generator implemented |
| `nx`, `ny` | integer | required | elements per direction, `>= 1` |
| `lx`, `ly` | number | required | domain extents [m], `> 0` |
| `x0`, `y0` | number | `0` | lower-left corner [m] |

## `material`

```json
"material": { "name": "Al 7075-T6", "youngs_modulus": 71.7e9,
              "poisson_ratio": 0.33, "density": 2810.0 }
```

| Key | Type | Default | Constraint |
|-----|------|---------|------------|
| `name` | string | `"material"` | carried into result files |
| `youngs_modulus` | number [Pa] | required | `> 0` |
| `poisson_ratio` | number | required | in `(-1, 0.5)` so both idealisations stay positive definite |
| `density` | number [kg/m^3] | `0` | `>= 0`; a positive value is required for modal analysis |

## `model`

```json
"model": { "thickness": 0.008, "stress_state": "plane_stress",
           "integration": { "stiffness_points": 2, "mass_points": 3,
                            "edge_points": 2 } }
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `thickness` | number [m] | `1.0` | out-of-plane thickness, `> 0` |
| `stress_state` | string | `plane_stress` | or `plane_strain` |
| `integration.stiffness_points` | integer | `2` | Gauss points per direction for `K_e`, 1-4 |
| `integration.mass_points` | integer | `3` | Gauss points per direction for `M_e`, 1-4 |
| `integration.edge_points` | integer | `2` | Gauss points along an edge, 1-4 |

## Regions

Boundary conditions, loads and passive regions all take a `region` object. A
region is the union of one or more primitives, optionally complemented.

```json
"region": { "name": "bolt_attachments", "invert": false, "tolerance": 0.0,
            "any_of": [
              { "circle": { "center": [0.030, 0.045], "radius": 0.0105 } },
              { "circle": { "center": [0.030, 0.155], "radius": 0.0105 } }
            ] }
```

A single primitive may also be written directly in the `region` object, without
`any_of`.

| Primitive | Form | Selects |
|-----------|------|---------|
| `all` | `"all": true` | everything |
| `box` | `{ "xmin":, "xmax":, "ymin":, "ymax": }`, any subset | points inside the axis-aligned box; omitted sides are unbounded |
| `circle` | `{ "center": [x, y], "radius": r }` | the closed disc |
| `annulus` | `{ "center": [x, y], "inner_radius": r0, "radius": r1 }` | the closed ring, `r0 < r1` |
| `node_ids` | `[...]` | explicit nodes (node regions only) |
| `element_ids` | `[...]` | explicit elements (element regions only) |
| `nearest_node` | `[x, y]` | the single node closest to the point (node regions only) |

| Region key | Type | Default | Meaning |
|------------|------|---------|---------|
| `name` | string | auto | used in diagnostics |
| `invert` | bool | `false` | select the complement of the union |
| `tolerance` | number [m] | `0` | absolute geometric tolerance; `0` means `1e-9` times the mesh diagonal |

A degenerate box such as `{"xmax": 0.0}` reliably captures the line of nodes at
`x = 0` because of that relative tolerance. Element regions are tested at the
element centroid. Naming two primitives in one region object without `any_of` is
an error, as is a region that selects nothing.

## `boundary_conditions`

```json
"boundary_conditions": [
  { "name": "clamped_root", "fix": ["x", "y"], "value": [0.0, 0.0],
    "region": { "box": { "xmax": 0.0 } } }
]
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | auto | diagnostics |
| `fix` | array of `"x"` / `"y"` | required, non-empty | which components to prescribe |
| `value` | `[u_x, u_y]` [m] | `[0, 0]` | prescribed displacement |
| `region` | object | required | node region |

Non-zero values are handled by static condensation, not by modifying the load
vector, so the reactions stay exact.

## `load_cases`

```json
"load_cases": [
  { "name": "down_limit", "weight": 1.0,
    "point_loads": [
      { "name": "lug", "force": [0.0, -9000.0], "distribution": "total",
        "region": { "annulus": { "center": [0.270, 0.100],
                                 "inner_radius": 0.0075, "radius": 0.0105 } } }
    ],
    "tractions": [
      { "name": "upper_skin", "traction": [0.0, -1.2e5],
        "region": { "box": { "ymin": 0.10 } } }
    ]
  }
]
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | auto | used in file names and reports |
| `weight` | number | `1.0` | weight in the multi-load objective, `>= 0`; weights are normalised to sum to 1 |
| `point_loads` | array | `[]` | concentrated nodal forces |
| `tractions` | array | `[]` | distributed edge loads |
| `prescribed_displacement_only` | bool | `false` | declares a case with no applied force, driven by the prescribed displacements (the patch test, enforced-deflection studies) |

**Point load**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `force` | `[F_x, F_y]` [N] | required | see `distribution` |
| `distribution` | `"total"` / `"per_node"` | `"total"` | `total` divides the resultant among the selected nodes, so the total force is mesh independent; `per_node` applies `force` to each node |
| `region` | object | required | node region |

**Traction**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `traction` | `[t_x, t_y]` [Pa] | required | stress vector in **global** components, not a normal pressure |
| `region` | object | required | node region; a boundary edge is loaded when *both* its end nodes lie inside |

Tractions are integrated to consistent nodal forces, so the resultant is exactly
`traction * (loaded length) * thickness` at any mesh resolution. A traction
region matching no boundary edge is an error.

A load case with neither `point_loads` nor `tractions` is an error unless
`prescribed_displacement_only` is set.

## `solver`

```json
"solver": {
  "linear": { "type": "simplicial_ldlt", "residual_tolerance": 1e-8,
              "iterative_tolerance": 1e-12, "max_iterations": 0,
              "pivot_tolerance": 1e-14 },
  "equilibrium_tolerance": 1e-6,
  "check_model": true
}
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `linear.type` | string | `simplicial_ldlt` | `simplicial_ldlt`, `simplicial_llt`, `sparse_lu`, `conjugate_gradient`, `dense_lu` |
| `linear.residual_tolerance` | number | `1e-8` | scaled residual accepted after each solve |
| `linear.iterative_tolerance` | number | `1e-12` | requested relative residual for CG |
| `linear.max_iterations` | integer | `0` | CG iteration cap; `0` uses Eigen's default of `2n` |
| `linear.pivot_tolerance` | number | `1e-14` | smallest accepted `min/max` LDL^T pivot ratio |
| `equilibrium_tolerance` | number | `1e-6` | relative force-balance error before `SolverError` |
| `check_model` | bool | `true` | run the rigid-body and floating-region diagnostics before assembling |

## `modal`

```json
"modal": { "enabled": true, "num_modes": 6, "mass_type": "consistent",
           "tolerance": 1e-10, "residual_tolerance": 1e-6,
           "max_iterations": 200, "shift_factor": 0.0,
           "rigid_body_ratio": 1e-10, "seed": 20240917,
           "analyse_optimised_topology": true,
           "compare_mass_matched_baseline": true }
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | run free-vibration analysis |
| `num_modes` | integer | `6` | lowest modes requested, `>= 1` |
| `mass_type` | string | `consistent` | or `lumped` (row-sum) |
| `tolerance` | number | `1e-10` | relative eigenvalue change for convergence |
| `residual_tolerance` | number | `1e-6` | eigenpair residual, also part of the stopping rule |
| `max_iterations` | integer | `200` | subspace-iteration cap |
| `shift_factor` | number | `0` | shift as a multiple of `min_i K_ii/M_ii`; `0` disables |
| `rigid_body_ratio` | number | `1e-10` | below this times the stiffness/mass scale, a mode is flagged as rigid-body |
| `seed` | integer | `20240917` | seed for the filler part of the starting basis |
| `analyse_optimised_topology` | bool | `true` | also compute *mode shapes* of the thresholded solid interpretation |
| `compare_mass_matched_baseline` | bool | `true` | topology runs also analyse an equal-mass uniform plate |

## `topology`

```json
"topology": {
  "enabled": true,
  "volume_fraction": 0.35,
  "initial_density": -1,
  "simp": { "penalty": 3.0, "emin_ratio": 1e-9, "mass_floor": 1e-9,
            "mass_interpolation": "penalty_matched" },
  "filter": { "type": "density", "radius_elements": 3.0 },
  "optimizer": { "max_iterations": 600, "change_tolerance": 1e-2,
                 "objective_tolerance": 5e-5, "objective_window": 20,
                 "move_limit": 0.2, "damping": 0.5,
                 "continuation_steps": 3, "penalty_start": 1.5,
                 "continuation_iterations": 30,
                 "volume_tolerance": 1e-10, "max_bisections": 200,
                 "history_stride": 1, "interpretation_threshold": 0.5 },
  "passive_regions": [
    { "name": "bolt_holes", "type": "void",
      "region": { "circle": { "center": [0.030, 0.045], "radius": 0.010 } } }
  ]
}
```

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | `false` | `sparlab_topopt` requires `true` |
| `volume_fraction` | number | `0.4` | target `nu`, in `(0, 1]`; measured on the physical density over the **whole** domain, passive regions included |
| `initial_density` | number | `-1` | starting value for free variables; negative means "use `volume_fraction`" |

**`simp`**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `penalty` | number | `3.0` | `p >= 1` |
| `emin_ratio` | number | `1e-9` | `E_min/E_0`, in `(0, 1)` |
| `mass_floor` | number | `1e-9` | mass floor for `penalty_matched` |
| `mass_interpolation` | string | `penalty_matched` | or `linear` |

**`filter`**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `type` | string | `density` | `density`, `sensitivity` or `none` |
| `radius` | number [m] | - | radius in metres; takes precedence when `> 0` |
| `radius_elements` | number | `1.5` | radius as a multiple of the mean element size |

Specify the radius in **metres** when comparing meshes: that holds the minimum
length scale fixed, which is what makes the optimum mesh convergent.

**`optimizer`**

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `max_iterations` | integer | `200` | iteration cap |
| `change_tolerance` | number | `1e-2` | stop when `max |dx|` falls below this |
| `objective_tolerance` | number | `0` | stop when the relative compliance change over `objective_window` iterations falls below this; `0` disables |
| `objective_window` | integer | `5` | window for the above |
| `move_limit` | number | `0.2` | maximum per-iteration density change, in `(0, 1]` |
| `damping` | number | `0.5` | OC exponent `eta`, in `(0, 1]` |
| `continuation_steps` | integer | `1` | `> 1` enables penalty continuation |
| `penalty_start` | number | `1.0` | starting penalty, `>= 1` and `<= simp.penalty` |
| `continuation_iterations` | integer | `25` | iterations per continuation stage |
| `volume_tolerance` | number | `1e-10` | relative volume error of the multiplier bisection |
| `max_bisections` | integer | `200` | bisection cap |
| `history_stride` | integer | `1` | record a density snapshot every N iterations; `0` disables |
| `interpretation_threshold` | number | `0.5` | density threshold used when reading the field as geometry, in `(0, 1)` |

**`passive_regions`** - an array of

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | auto | diagnostics |
| `type` | string | `solid` | `solid` (`x = 1`) or `void` (`x = density`) |
| `density` | number | `0` | density imposed on a void region, in `[0, 1)` |
| `region` | object | required | element region |

## `output`

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `csv` | bool | `true` | per-node and per-element CSV tables |
| `vtk` | bool | `true` | legacy VTK field files for ParaView |
| `density_history` | bool | `true` | density snapshots for the animation |
| `mode_shapes` | bool | `true` | mode-shape CSV and per-mode VTK |

`summary.json` and `mesh.json` are always written: a result directory should be
self-describing.

## Command-line overrides

`sparlab_topopt` accepts overrides so a parametric study never has to edit or
duplicate a deck. Each corresponds to one deck field:

```
--volume-fraction        topology.volume_fraction
--penalty                topology.simp.penalty
--filter-radius          topology.filter.radius            (metres)
--filter-radius-elements topology.filter.radius_elements   (cells)
--filter-type            topology.filter.type
--max-iterations         topology.optimizer.max_iterations
--nx --ny                mesh.nx, mesh.ny
--youngs-modulus         material.youngs_modulus
--load-weights w1,w2,... one weight per load case, in deck order
--modes N                modal.enabled = true, modal.num_modes = N
--tag NAME               appended to the case name in summary.json
```

Run any app with `--help` for its full flag list.
