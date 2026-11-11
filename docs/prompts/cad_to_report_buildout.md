# Build-out prompt: CAD file + load spec → structural analysis PDF report

## Goal

Add a single end-to-end entry point to SparLab that takes **a CAD file and a
declarative load specification** and produces **a complete structural analysis
report as a PDF**, with no manual meshing step and no hand-written solver deck.

Target usage:

```bash
python3 python/scripts/analyze_part.py --spec specs/my_bracket.yaml --output results/my_bracket
# -> results/my_bracket/report.pdf
```

Everything between those two things is the deliverable.

## What already exists (do not rebuild)

Read these before writing code. The solver, the mesh importer, the load
machinery and the plotting layer are all already in the repository; this task
is a front-end and a reporting back-end around them.

- `src/io/MeshReader.cpp` / `include/sparlab/io/MeshReader.hpp` — reads Gmsh
  `.msh` (2.2, 4.1) and Abaqus/CalculiX `.inp`. Handles `scale` (0.001 for mm
  CAD), orientation repair, duplicate-node merge, quality refusal of inverted
  and folded cells. **Named sets**: Gmsh physical groups and Abaqus
  `*NSET`/`*ELSET` become node/element sets addressable from the deck.
- `src/io/Config.cpp` — the JSON deck schema. `mesh.type: "file"` with `path`,
  `scale`, `format`; `boundary_conditions[].region.group`; `load_cases[]` with
  `tractions`, `pressures`, point loads, self-weight, rotation, thermal.
- `configs/verification/engine_mount_tet10_analysis.json` — a worked deck that
  reads an imported mesh, clamps a named hole group, and applies two
  pin-hole traction load cases. **This is the template the generator should
  emit.** Match its structure.
- `apps/sparlab_solve.cpp` — the analysis binary. Flags include
  `--config`, `--output`, `--buckling N`, `--export-calculix`.
- `src/io/ResultWriter.cpp` — writes `summary.json` (equilibrium check, modal
  results, per-load-case norms), CSV tables, legacy VTK.
- `python/scripts/make_meshes.py` — a worked, reproducible Gmsh driver:
  pinned `General.NumThreads`, `Mesh.ElementOrder`, `Mesh.MeshSizeMin/Max`,
  `Mesh.MeshSizeFromCurvature`, `Mesh.Algorithm`, `Mesh.Algorithm3D`, and a
  companion `.json` recording the Gmsh version and every option.
  **Reuse this pattern exactly** — the new CAD mesher must be as reproducible.
- `python/sparlab_viz/` — `loaders.py` (`load_case`, `load_mesh`, `load_csv`,
  `CaseResults`, legacy-VTK reader), `plots.py`, `plots3d.py`, `fields.py`,
  `style.py`. Figure generation for a result directory already exists.
- `python/scripts/plot_case.py` — generates every applicable figure for one
  result directory. The report should call into the same `sparlab_viz`
  functions rather than drawing fields from scratch.

## Stage 1 — CAD import and meshing

New module: `python/sparlab_viz/cad.py` (or `python/scripts/cad_to_mesh.py`
plus a library module — your call, but the geometric classification must be
importable and unit-testable, not buried in a script).

- Load CAD via the Gmsh Python API with the OpenCASCADE kernel:
  `gmsh.model.occ.importShapes(path)` for STEP (`.step`, `.stp`), IGES
  (`.iges`, `.igs`) and BREP (`.brep`). Gmsh is already an optional dependency
  (`pip install gmsh`); follow the existing `make_meshes.py` import guard that
  exits with an actionable message.
- Mesh with linear or quadratic tetrahedra. Default to **Tet10**
  (`Mesh.ElementOrder = 2`) for curved CAD solids: the reader supports C3D10
  and Gmsh second-order tets, curved edges land on the real hole surfaces, and
  it avoids the shear/volumetric stiffness of Tet4. Allow `element_order: 1`
  in the spec.
- Pin every Gmsh option and record the Gmsh version plus the full option dict
  in a companion `<part>_mesh.json`, as `make_meshes.py` does. Same version
  and same spec must give the same mesh.
- Sizing: curvature-driven by default (`Mesh.MeshSizeFromCurvature`) so holes
  and fillets refine without the user sizing them by hand. Honour
  `mesh.size`, `mesh.size_min`, `mesh.curvature` from the spec.
- **Hard constraint from the reader**: it refuses mixed cell types, prisms and
  pyramids. So set `Mesh.RecombineAll = 0` and emit all-tetrahedral meshes.
  If the user asks for hexes, fail with a clear message rather than writing a
  mesh the reader will reject.
- Write the mesh where the reader can take it: Gmsh `.msh` 4.1, or `.inp`.
  Prefer `.msh` 4.1 (ASCII, `Mesh.Binary = 0`).
- Report after meshing: node and cell counts, bounding box in both CAD units
  and metres, min/max element size, a worst-element quality figure, and the
  mesh volume against the CAD solid volume from
  `gmsh.model.occ.getMass(3, tag)` — a volume mismatch beyond a tolerance is
  the single best early warning that the mesh does not represent the part.

## Stage 2 — Face naming (the part that needs real design thought)

A STEP file's faces are numbered, not named, but the whole deck addresses
geometry **by group name**. This is the crux of the task. Support three
mechanisms, in this precedence order:

1. **Explicit map** (always wins, always available as an escape hatch). The
   spec may name CAD entities directly:
   ```yaml
   groups:
     bolt_holes: {surfaces: [21, 24, 27, 30]}
     pin_hole:   {surfaces: [18]}
   ```
   Mirrors the `groups` block already in `configs/meshes/engine_mount_3d.json`.

2. **Geometric selectors** — the automatic path, and what makes this usable.
   Query OCC for each face's type and parameters
   (`gmsh.model.getType`, `getBoundingBox`, `occ.getMass`, surface normals
   at sample points) and select by predicate:
   ```yaml
   groups:
     bolt_holes:
       cylindrical: true
       radius: {min: 4.0, max: 6.0}     # CAD units
       axis: z
     base:
       planar: true
       normal: [0, 0, -1]
       at: {z: min}                      # extreme face along an axis
     pin_hole:
       cylindrical: true
       radius: {min: 9.0, max: 11.0}
       axis: x
   ```
   Implement at minimum: `cylindrical` / `planar` / `conical` / `spherical`,
   `radius` range, `axis` direction (with tolerance), `normal` direction,
   `at: {axis: min|max}` for extreme faces, `area` range, and a
   `box` / `inside` bounding-box filter. Each selector resolves to a list of
   CAD face tags, which then become a Gmsh physical group of dimension 2.

3. **Discovery mode** — `--list-faces` dumps every CAD face with its type,
   area, radius (if analytic), axis and centroid, plus a rendered image with
   face tags labelled if that is cheap to do. This is how a user writes the
   explicit map for geometry the selectors cannot express. Do not skip this:
   it is the difference between a tool that works on the test part and a tool
   that works on the next part.

Requirements on naming:
- Every selector must report **what it matched**: the group name, the face
  tags, and each face's area and type. Silent matching is a correctness trap.
- A selector matching **zero** faces is a hard error naming the selector, not
  a warning. A deck that clamps nothing solves a floating body; the solver's
  rigid-body diagnostic will catch it, but the error belongs here where the
  cause is legible.
- A face matched by **two** groups is an error unless the spec sets
  `allow_overlap: true`.
- Also create volume groups (`design`, `rings`, passive regions) the same way,
  so this front-end can feed `sparlab_topopt` later, not just `sparlab_solve`.

## Stage 3 — The input specification

One file carries everything: the CAD path, units, mesh controls, material,
groups, supports, load cases, and which analyses to run. YAML for authoring
(`pyyaml`), but **accept JSON too** and treat YAML as optional — if `pyyaml`
is absent, JSON must still work, so the tool has no new hard dependency.

Sketch (fill in against the real `Config.cpp` schema — do not invent keys the
deck does not accept):

```yaml
name: my_bracket
cad:
  path: parts/bracket.step
  units: mm              # -> mesh.scale 0.001

mesh:
  element_order: 2
  size: 4.0
  curvature: 12

material:
  name: Al 7075-T6
  youngs_modulus: 71.7e9
  poisson_ratio: 0.33
  density: 2810.0
  yield_strength: 503e6   # for the safety factor in the report

groups:
  bolt_holes: {cylindrical: true, radius: {min: 4, max: 6}, axis: z}
  pin_hole:   {cylindrical: true, radius: {min: 9, max: 11}, axis: x}

supports:
  - {group: bolt_holes, fix: [x, y, z]}

load_cases:
  - name: vertical
    pressures:  [{group: pin_hole, resultant: [0, 0, -12000]}]   # N
  - name: lateral
    tractions:  [{group: pin_hole, resultant: [5000, 0, 0]}]
  - name: vertical_plus_gravity
    pressures:  [{group: pin_hole, resultant: [0, 0, -12000]}]
    self_weight: {direction: [0, 0, -1], g: 9.81}

analyses:
  static: true
  modal:    {enabled: true, num_modes: 6}
  buckling: {enabled: true, num_modes: 4}
```

**Resultant-to-traction conversion is a required convenience.** The existing
deck has a hand-computed comment — `12 kN / (2 pi 10 mm x 60 mm)` — turned into
a magic number `3.183098861837907e6`. Users should give a force in newtons and
the tool should divide by the **actual meshed face area** (sum of the group's
face areas from OCC, or the integrated mesh facet area) and emit the traction
vector, printing both the area and the resulting traction so the arithmetic is
auditable. Keep raw `traction:` as an option for anyone who wants it.

Validate the spec with clear per-field errors before meshing. A typo in a
material property should not surface as a solver failure ten minutes later.

## Stage 4 — Deck generation and solve

- Emit a deck matching `configs/verification/engine_mount_tet10_analysis.json`
  in structure, written into the output directory so the run is reproducible
  and hand-editable afterwards. Include generated comments recording where
  each number came from (the group, the face area, the resultant).
- Invoke `./build/bin/sparlab_solve --config <deck> --output <dir>` with
  `--buckling N` when requested. Locate the binary relative to the repository
  root as the existing scripts do; fail with a message pointing at
  `make build` if it is absent.
- Surface solver diagnostics verbatim. SparLab already detects rigid-body
  under-constraint, floating regions, invalid cells and singular matrices with
  the element and node named — that output is the most useful thing a user
  gets when their model is wrong. Do not swallow it.
- Non-zero exit from the solver is a non-zero exit here, with the solver's own
  message, and no half-written PDF.

## Stage 5 — The PDF report

**Use `matplotlib.backends.backend_pdf.PdfPages`.** matplotlib is already a
dependency and `python/sparlab_viz/style.py` already defines the house figure
style, so this adds no new dependency and the report matches the existing
figures. (If you find you need flowing text and real tables badly enough,
`reportlab` may be added as an *optional* extra with a matplotlib fallback —
but try PdfPages first and keep it dependency-free.)

New module: `python/sparlab_viz/report.py`, driven by
`python/scripts/analyze_part.py`.

Contents, in order:

1. **Cover** — part name, CAD file, date, SparLab version/commit, a rendered
   view of the meshed part.
2. **Model summary** — material table; mesh statistics (nodes, cells, element
   type, size range, worst quality); **mesh volume vs CAD volume**; the
   resolved group table (name → face tags → total area); supports; load cases
   with their resultants and derived tractions.
3. **Model verification** — this is what makes the report trustworthy rather
   than decorative. Pull from `summary.json`:
   - the **equilibrium check** (`applied_force_N`, `reaction_force_N`,
     `force_residual_N`, `relative_force_error`, and the moment equivalents).
     State it as a pass/fail against a tolerance. A model whose reactions do
     not balance its loads has not been solved, whatever the stress plot shows.
   - solver type, iteration count and final residual.
   - any mesh-quality or orientation repairs the reader reported.
4. **Per load case**, one or two pages each:
   - deformed shape (exaggerated) over the undeformed outline
   - displacement magnitude field, and the peak displacement with its node and
     location
   - von Mises field, and the peak with its element and location
   - principal stresses where they clarify the state
   - support reactions vs applied load
   - **safety factor** = `yield_strength / peak_von_mises`, stated plainly,
     with the caveat below.
5. **Modal** — frequency table, modal mass fractions, lowest mode shapes.
6. **Buckling** — load factors per load case and the lowest mode shapes; flag
   any factor below 1 prominently.
7. **Summary and caveats** — a short table of the governing numbers across
   load cases (peak von Mises, peak displacement, min safety factor, lowest
   frequency, lowest buckling factor), then the limitations. State them, do
   not bury them:
   - linear elastic analysis: the safety factor is meaningful only below yield,
     and peak stress at a re-entrant corner or a sharp CAD edge is a
     **mesh-dependent singularity**, not a real stress
   - no mesh convergence study was performed unless the user ran one
   - contact, preload and bolt stiffness are not modelled unless specified
   - single-material, room-temperature unless a thermal case was given

Reuse `sparlab_viz.plots` / `plots3d` / `fields` for every field plot. If a
needed plot does not exist, add it there, where `plot_case.py` can use it too.

## Stage 6 — Convergence (do this if the budget allows; say so if you skip it)

`--convergence` re-meshes at 2–3 element sizes, re-solves, and adds a page
plotting peak von Mises and peak displacement against element count. This is
the difference between a number and an engineering result, and the repository
already values this kind of check everywhere else. If you leave it out, leave
a clearly-scoped TODO and say so in the summary — do not let its absence be
discovered from the code.

## Testing

- Unit-test the geometric selectors against a small STEP file committed under
  `tests/data/` (or generated by OCC primitives in the test itself, which
  avoids committing a binary): a block with known cylindrical holes, asserting
  that a radius+axis selector picks exactly the intended face tags.
- Unit-test resultant→traction conversion against a hand-computed case. The
  existing `12 kN / (2 pi 10 mm x 60 mm)` number in
  `engine_mount_tet10_analysis.json` is a ready-made check.
- Spec validation tests: missing field, unknown field, zero-match selector,
  overlapping groups.
- One end-to-end test behind a marker that skips without `gmsh`: a simple
  cantilever-like CAD solid under a tip load, asserting the tip deflection
  against beam theory within the discretisation error, and asserting
  `report.pdf` exists and is non-trivially sized.
- Follow the repository's existing test conventions (Catch2 for C++; match
  whatever `tests/` and the Python scripts already do for Python).

## Documentation

- A `README.md` section: the one-command usage, a complete worked spec, and
  the face-naming story including `--list-faces`.
- A worked example committed under `specs/` or `configs/parts/` that runs
  against a committed CAD file, so `make`-level invocation is possible.
- A `make` target in the spirit of the existing thin wrappers, e.g.
  `make analyze SPEC=specs/my_bracket.yaml`.

## Constraints

- **No new hard dependencies.** `gmsh` stays an optional import with an
  actionable error, as in `make_meshes.py`. The PDF must build with the
  existing `numpy`/`pandas`/`matplotlib`/`pillow` layer. `pyyaml` optional,
  JSON always works.
- **Do not modify the C++ solver core** unless something is genuinely
  missing; this is a front-end plus a reporting layer. If you do find a real
  gap in `MeshReader` or `Config`, flag it rather than working around it with
  a fragile shim.
- Match the existing code's style: module docstrings that explain the physics
  and the reproducibility story, not just the API.
- Reproducibility is a requirement, not a nicety: same spec plus same Gmsh
  version gives the same mesh, the same deck and the same numbers.

## Definition of done

`python3 python/scripts/analyze_part.py --spec <worked example> --output results/x`
runs from a clean build on a committed CAD file and produces a `report.pdf`
whose equilibrium check passes, whose peak stress and displacement are
reported with their locations, and whose caveats are stated. Tests pass.
README documents it. Work is committed on `claude/brave-bardeen-2h386b`.

## Order of work

1. Stage 1 + 2 (`cad.py`, selectors, `--list-faces`) with their unit tests —
   this is the risky part; get it right before building on it.
2. Stage 3 + 4 (spec, deck emission, solve) — prove the pipeline end to end
   on the worked example and check the numbers against the existing
   `engine_mount_tet10_analysis` results.
3. Stage 5 (the PDF) once there are real results to report.
4. Stage 6, docs, `make` target.

Report honestly at each stage: what works, what was verified against what, and
what was left out.
