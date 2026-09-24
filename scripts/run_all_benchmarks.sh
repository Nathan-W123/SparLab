#!/usr/bin/env bash
# Run every benchmark configuration in configs/benchmarks/, plus the small
# static+modal analysis decks (2-D and 3-D), into <results-dir>/<case>/.
#
# The four compliance benchmarks take seconds to a minute each; the
# stress-constrained L-bracket a few tens of seconds; the 3-D bracket
# (15k DOFs, 150 iterations) around ten minutes.
#
# usage: scripts/run_all_benchmarks.sh [results-dir]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RESULTS="${1:-$SPARLAB_ROOT/results}"
require_binaries sparlab_solve sparlab_topopt

"$SPARLAB_ROOT/scripts/run_benchmark.sh" cantilever_analysis "$RESULTS"
"$SPARLAB_ROOT/scripts/run_benchmark.sh" block_3d_analysis "$RESULTS"
for case in cantilever_beam mbb_beam aerospace_bracket wing_rib l_bracket_stress bracket_3d; do
  "$SPARLAB_ROOT/scripts/run_benchmark.sh" "$case" "$RESULTS"
done

banner "all benchmarks complete"
echo "results under: $RESULTS"
