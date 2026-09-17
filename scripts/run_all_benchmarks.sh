#!/usr/bin/env bash
# Run every benchmark configuration in configs/benchmarks/, plus the small
# static+modal analysis deck, into <results-dir>/<case>/.
#
# usage: scripts/run_all_benchmarks.sh [results-dir]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RESULTS="${1:-$SPARLAB_ROOT/results}"
require_binaries sparlab_solve sparlab_topopt

"$SPARLAB_ROOT/scripts/run_benchmark.sh" cantilever_analysis "$RESULTS"
for case in cantilever_beam mbb_beam aerospace_bracket wing_rib; do
  "$SPARLAB_ROOT/scripts/run_benchmark.sh" "$case" "$RESULTS"
done

banner "all benchmarks complete"
echo "results under: $RESULTS"
