# Twisted tower: a tower skin optimised by SparLab

`sparlab_lattice_tower.mp4`: 12.0 s, 1080 × 1080, 30 fps, H.264 (High), 1.9 MB.
It is small enough that no 720p copy is needed.

The material loss comes from a SparLab topology optimisation, not an
effect.

- **Design space:** the skin of a 100 × 100 × 250 mm square tower: four
  walls and a cap, 11 264 MITC4 shells, 6 mm of Ti-6Al-4V.
- **Supports and loads:** the base is clamped. The top is pressed down
  with 1 kN and twisted with 500 N along each top edge, in two load cases,
  one for each sense of twist.
- **Solid regions:** the cap and the top two rows, where the twist goes
  in, and the bottom row are kept solid.
- **Optimisation:** keeps 35 % of the skin. SIMP with p = 3, a density
  filter of 2.5 cells and no projection.
- **Result:** converged on its objective-stall test after 176 iterations.
  Thresholded at 0.5, the part is 4 072 of the 11 264 shells in one
  connected group with no islands.
- **Engine:** SparLab 1.0.0, Release build (gcc 13.3.0, Eigen 3.4.0), as
  recorded in the run's `summary.json`.

The first frame is the whole skin and the last is the engine's final
physical density. The frames in between show the optimiser's density
snapshots, saved every 2 iterations. Each is blended with its neighbour
and contoured at the level that keeps a smoothly falling share of the
skin, so the material it rates lowest goes first. They are a visualisation
of the run: the optimiser holds the volume at 35 % from its first
iteration and moves material around, rather than carving it away.

## Look

- **Camera path:** a static perspective camera, 30° field of view.
  - It sits on the direction (3, 1.6, 4), about 18° above the horizon.
  - Its distance is set so the tower's bounding sphere (half-diagonal
    0.1436 m, drawn at 10 units per metre) fills the frame with an 8 %
    margin: 5.99 units.
  - It looks at the tower's centre.
  - The tower turns about its vertical axis, linearly from −45° to +45°
    over the 360 frames.
- **Lighting:** image-based only. three.js `RoomEnvironment` goes through
  a `PMREMGenerator` (blur 0.04) as the scene environment, with no direct
  lights. three.js r169 defaults: no tone mapping, sRGB output.
- **Colour map:** none. One `MeshStandardMaterial`, colour `#c9ced6`,
  metalness 0.7, roughness 0.3. Normals come from the mesh
  (`computeVertexNormals`).
- **Background:** solid `#0d1117`.
- **Typography:** none. The video has no text.
- **Timing / beats** (30 fps, 360 frames):
  - Frames 0–44 (1.5 s): the closed tower (mesh 0).
  - Frames 45–254 (7 s): the material goes. Mesh index = round((frame − 45) / 209 × 210).
  - Frames 255–359 (3.5 s): holds on the final design (mesh 210).
- **Drawing the skin:** each shell's density is painted onto a 32 × 32 × 80
  voxel grid, cut to a fixed 6.25 mm skin, and contoured by marching cubes
  (twice upsampled). Taubin smoothing (40 passes) keeps face vertices on
  their faces.
- **Encode:** libx264, preset slow, CRF 18, High profile, yuv420p,
  `+faststart`.

## Files (exactly as run)

| File | Role |
|---|---|
| `skin_mesh.py` | Writes the shell mesh `skin.inp`, `skin_cells.json` and the base deck `config.json` |
| `config_simp.json` | The deck the engine ran (`config.json` with projection off and 220 iterations at most) |
| `skin_frames.py` | Turns the run's density history into the 211 frame meshes |
| `box_surface.py`, `meshio_bin.py` | Helper modules `skin_frames.py` imports |
| `box.html` | three.js page that draws one frame (timeline, camera, material) |
| `capture.mjs` | Renders `box.html` frame by frame in headless Chromium (Playwright) |
| `encode.sh` | ffmpeg encode of the PNG frames |

The scripts keep the absolute paths of the machine they ran on: the
`sys.path` entries in `skin_frames.py`, Playwright's location in
`capture.mjs` and the ffmpeg binary in `encode.sh`. Adjust them to rerun
elsewhere. Run from this folder, `skin_frames.py` finds its helpers here,
since Python puts a script's own directory first on its path.

Tools: Python 3.11 with numpy 2.4.6, scipy 1.17.1, scikit-image 0.26.0 and
pandas 3.0.5; Node 22.22.2 with Playwright 1.56.1 (Chromium with
SwiftShader WebGL); three.js 0.169.0 in `node_modules` next to
`box.html`; ffmpeg 7.0.2 (imageio-ffmpeg build).

## Commands (as run)

```bash
O=/tmp/claude-0/-home-user-SparLab/d56879da-faaf-5259-b306-3c582599460d/scratchpad/options
S="$O/skin"
LD=/tmp/claude-0/-home-user-SparLab/d56879da-faaf-5259-b306-3c582599460d/scratchpad/linkedin

# 1. Shell mesh and base deck (defaults: 32 x 80 cells per wall, 35 %, filter 2.5 cells)
python3 "$S/skin_mesh.py" "$S"

# 2. The deck that was run: projection off, at most 220 iterations
python3 - "$S" "$O" <<'EOF'
import json, sys
S, O = sys.argv[1], sys.argv[2]
c = json.load(open(f"{S}/config.json"))
c["name"] = "tower_skin_simp"
c["topology"]["projection"] = {"enabled": False}
c["topology"]["optimizer"]["max_iterations"] = 220
json.dump(c, open(f"{S}/config_simp.json", "w"), indent=1)
EOF

# 3. The optimisation (from the repository root, in the background, beside another run)
cd /home/user/SparLab && nohup bash -c "OMP_NUM_THREADS=2 ./build/bin/sparlab_topopt --config '$S/config_simp.json' --output '$S/run_simp' > '$S/run_simp.log' 2>&1; echo \"exit \$?\" >> '$S/run_simp.log'" > /dev/null 2>&1 &

# 4. Frame meshes, then the turn the page reads from index.json
python3 "$S/skin_frames.py" "$S/run_simp" "$S" "$O/render/tower_meshes" 211 10
python3 - "$O/render/tower_meshes/index.json" <<'EOF'
import json, sys
p = sys.argv[1]; d = json.load(open(p)); d["spin0"] = -45.0; d["spin"] = 90.0; json.dump(d, open(p, "w"))
EOF

# 5. Frames: box.html served from $O/render (node_modules holds three 0.169.0);
#    the second argument is unused when QUERY is set
cd "$O/render" && PAGE=box.html QUERY='meshes=tower_meshes' node capture.mjs "$O/render" x "$O/render/tower_frames" all

# 6. Encode
"$LD/encode.sh" "$O/render/tower_frames" "$O/sparlab_lattice_tower.mp4"
```
