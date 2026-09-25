#!/usr/bin/env bash
# Cross-validate the static solver against CalculiX and scikit-fem.
#
# usage: scripts/run_cross_validation.sh [results-dir]
#
# Solves the reference analysis decks - the 2-D cantilever and the 3-D block
# on Q4 / Hex8 and on Tri3 / Tet4 elements, the two Gmsh parts (the lug
# bracket at its real Poisson ratio and at nu = 0, the engine mount on Tet4
# and on curved Tet10 cells), and an axially loaded column on Hex8, Tet4 and
# Tet10 with its linear buckling check - with sparlab_solve (exporting
# CalculiX decks), then compares the nodal displacements node by node, and
# the buckling load factors mode by mode, with CalculiX (ccx, static and
# *BUCKLE) and scikit-fem. Exits non-zero if any comparison exceeds its
# documented tolerance.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RESULTS="${1:-$SPARLAB_ROOT/results}"
PYTHON="${PYTHON:-python3}"
require_binaries sparlab_solve

banner "cross-validation: solving the reference decks"
for case in cantilever_analysis block_3d_analysis cantilever_tri_analysis \
            block_tet_analysis lug_bracket_nu0_analysis engine_mount_tet10_analysis \
            column_hex_buckling_analysis column_tet4_buckling_analysis \
            column_tet10_buckling_analysis; do
  "$BIN_DIR/sparlab_solve" --config "$SPARLAB_ROOT/configs/verification/$case.json" \
                           --output "$RESULTS/$case" --export-calculix
done
# The Gmsh parts as static analyses (their topology sections are ignored).
for case in lug_bracket_2d engine_mount_3d; do
  "$BIN_DIR/sparlab_solve" --config "$SPARLAB_ROOT/configs/benchmarks/$case.json" \
                           --output "$RESULTS/${case}_analysis" --export-calculix
done

banner "cross-validation: CalculiX and scikit-fem"
cd "$SPARLAB_ROOT"
CASES=()
for case in cantilever_analysis block_3d_analysis cantilever_tri_analysis \
            block_tet_analysis lug_bracket_nu0_analysis lug_bracket_2d_analysis \
            engine_mount_3d_analysis engine_mount_tet10_analysis \
            column_hex_buckling_analysis column_tet4_buckling_analysis \
            column_tet10_buckling_analysis; do
  CASES+=(--case "$RESULTS/$case")
done
"$PYTHON" python/scripts/cross_validate.py "${CASES[@]}" --output "$RESULTS/cross_validation"
