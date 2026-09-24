#!/usr/bin/env bash
# Run one benchmark configuration.
#
# usage: scripts/run_benchmark.sh <case-name> [results-dir]
#
# <case-name> is the stem of a file in configs/benchmarks/ (cantilever_beam,
# mbb_beam, aerospace_bracket, wing_rib, l_bracket_stress, bracket_3d) or in
# configs/verification/ (cantilever_analysis, block_3d_analysis). A deck with a
# topology section is run through sparlab_topopt; one without is run through
# sparlab_solve.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

CASE="${1:-cantilever_beam}"
RESULTS="${2:-$SPARLAB_ROOT/results}"
CONFIG="$SPARLAB_ROOT/configs/benchmarks/$CASE.json"

if [[ ! -f "$CONFIG" ]]; then
  CONFIG="$SPARLAB_ROOT/configs/verification/$CASE.json"
fi
if [[ ! -f "$CONFIG" ]]; then
  echo "error: no configuration found for case '$CASE'." >&2
  echo "available:" >&2
  ls "$SPARLAB_ROOT"/configs/benchmarks/*.json "$SPARLAB_ROOT"/configs/verification/*.json >&2
  exit 1
fi

require_binaries sparlab_solve sparlab_topopt

OUT="$RESULTS/$CASE"
mkdir -p "$OUT"

if grep -q '"enabled": *true' <<<"$(sed -n '/"topology"/,/^  }/p' "$CONFIG")"; then
  banner "topology optimization: $CASE"
  "$BIN_DIR/sparlab_topopt" --config "$CONFIG" --output "$OUT"
else
  banner "static + modal analysis: $CASE"
  "$BIN_DIR/sparlab_solve" --config "$CONFIG" --output "$OUT"
fi
