#!/usr/bin/env bash
# Cross-validate the static solver against CalculiX and scikit-fem.
#
# usage: scripts/run_cross_validation.sh [results-dir]
#
# Solves the 2-D cantilever and the 3-D block analysis decks with
# sparlab_solve (exporting CalculiX decks), then compares the nodal
# displacements node by node with CalculiX (ccx) and scikit-fem. Exits
# non-zero if any comparison exceeds its documented tolerance.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RESULTS="${1:-$SPARLAB_ROOT/results}"
PYTHON="${PYTHON:-python3}"
require_binaries sparlab_solve

banner "cross-validation: solving the reference decks"
for case in cantilever_analysis block_3d_analysis; do
  "$BIN_DIR/sparlab_solve" --config "$SPARLAB_ROOT/configs/verification/$case.json" \
                           --output "$RESULTS/$case" --export-calculix
done

banner "cross-validation: CalculiX and scikit-fem"
cd "$SPARLAB_ROOT"
"$PYTHON" python/scripts/cross_validate.py \
  --case "$RESULTS/cantilever_analysis" --case "$RESULTS/block_3d_analysis" \
  --output "$RESULTS/cross_validation"
