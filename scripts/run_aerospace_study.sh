#!/usr/bin/env bash
# Aerospace bracket parametric design study.
#
# usage: scripts/run_aerospace_study.sh [output-dir]
#
# Every arm changes exactly one thing relative to the shared baseline in
# configs/studies/aerospace_bracket_study.json, so the effect of each setting is
# isolated. All overrides go through sparlab_topopt command-line flags; the deck
# itself is never edited, which is what makes the study reproducible.
#
# Arms
#   volume_fraction  mass-stiffness trade: the Pareto-style curve. Modal
#                    analysis is enabled here so the first natural frequency of
#                    the interpreted topology is recorded against mass.
#   load_weighting   weight of the lateral load case, from 0 (pure vertical) to
#                    2 (lateral dominant).
#   mesh_fixed_r     mesh refinement with the filter radius fixed in METRES, so
#                    the minimum length scale is held constant. This is the
#                    correct way to look for mesh dependence.
#   mesh_fixed_cells mesh refinement with the filter radius fixed in CELLS, so
#                    the length scale shrinks with the mesh. Included to show
#                    the classical mesh-dependence pathology that filtering in
#                    physical units avoids.
#   penalty          SIMP penalty p from 1.5 (nearly convex, grey) to 5.
#   filter_radius    filter radius in cells, i.e. the minimum member size.
#   youngs_modulus   material stiffness from a magnesium-like to a steel-like
#                    modulus.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

OUT="${1:-$SPARLAB_ROOT/results/study}"
CONFIG="$SPARLAB_ROOT/configs/studies/aerospace_bracket_study.json"
require_binaries sparlab_topopt
mkdir -p "$OUT"

run_point() {
  local arm="$1"; shift
  local tag="$1"; shift
  local dir="$OUT/$arm/$tag"
  mkdir -p "$dir"
  printf '  %-16s %-22s' "$arm" "$tag"
  local start
  start=$(date +%s.%N)
  # Warn level, not error: the per-iteration log would be 40 000 lines of
  # noise, but a warning is the only record of something like an unresolved
  # filter radius or a non-convergent point, and each point's stderr is the
  # place to keep it.
  if SPARLAB_LOG_LEVEL=warn "$BIN_DIR/sparlab_topopt" \
        --config "$CONFIG" --output "$dir" --tag "${arm}_${tag}" \
        "$@" >"$dir/stdout.txt" 2>"$dir/stderr.txt"; then
    local elapsed
    elapsed=$(echo "$(date +%s.%N) - $start" | bc)
    printf 'ok   %6.1f s\n' "$elapsed"
  else
    printf 'FAILED\n'
    echo "    see $dir/stderr.txt" >&2
    tail -3 "$dir/stderr.txt" >&2 || true
    return 1
  fi
}

banner "aerospace bracket design study"
echo "baseline deck: $CONFIG"
echo "output:        $OUT"
echo

# --- 1. volume fraction (mass-stiffness trade, with a vibration metric) ------
for vf in 0.15 0.20 0.25 0.30 0.35 0.40 0.50 0.60; do
  run_point volume_fraction "vf_$vf" --volume-fraction "$vf" --modes 4
done

# --- 2. load-case weighting --------------------------------------------------
# weights are (down_limit, up_reversal, lateral)
for w in "1,0.5,0" "1,0.5,0.25" "1,0.5,0.5" "1,0.5,1" "1,0.5,2" "1,0,0"; do
  tag="w_$(echo "$w" | tr ',' '_')"
  run_point load_weighting "$tag" --load-weights "$w"
done

# --- 3. mesh refinement, filter radius fixed in metres -----------------------
# The baseline deck uses radius_elements = 1.5 on a 120 x 80 mesh of a 0.30 m
# domain, i.e. 1.5 * 0.0025 = 0.00375 m. Holding that in metres isolates the
# discretisation effect from the length-scale effect.
for n in 60 90 120 180 240; do
  ny=$((n * 2 / 3))
  run_point mesh_fixed_r "n_${n}x${ny}" --nx "$n" --ny "$ny" --filter-radius 0.00375
done

# --- 4. mesh refinement, filter radius fixed in cells ------------------------
for n in 60 90 120 180 240; do
  ny=$((n * 2 / 3))
  run_point mesh_fixed_cells "n_${n}x${ny}" --nx "$n" --ny "$ny" \
            --filter-radius-elements 1.5
done

# --- 5. SIMP penalty ---------------------------------------------------------
for p in 1.5 2.0 2.5 3.0 3.5 4.0 5.0; do
  run_point penalty "p_$p" --penalty "$p"
done

# --- 6. filter radius --------------------------------------------------------
for r in 1.0 1.5 2.0 3.0 4.0 6.0; do
  run_point filter_radius "r_$r" --filter-radius-elements "$r"
done

# --- 7. material stiffness ---------------------------------------------------
for e in 45e9 71.7e9 110e9 200e9; do
  run_point youngs_modulus "E_$e" --youngs-modulus "$e" --modes 4
done

banner "study complete"
echo "aggregate it with:"
echo "  python3 python/scripts/plot_study.py --study $OUT --figures docs/figures"
