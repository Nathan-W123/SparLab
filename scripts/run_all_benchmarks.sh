#!/usr/bin/env bash
# Run every benchmark configuration in configs/benchmarks/, plus the small
# static+modal analysis decks (2-D and 3-D), into <results-dir>/<case>/.
#
# The four compliance benchmarks take seconds to a minute each; the
# stress-constrained L-bracket a few tens of seconds; the 3-D bracket
# (15k DOFs, 150 iterations) around ten minutes with the direct solver; the
# Gmsh parts a few minutes; the 356k-DOF bracket_3d_large (multigrid CG)
# roughly half an hour.
#
# usage: scripts/run_all_benchmarks.sh [results-dir]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RESULTS="${1:-$SPARLAB_ROOT/results}"
require_binaries sparlab_solve sparlab_topopt

"$SPARLAB_ROOT/scripts/run_benchmark.sh" cantilever_analysis "$RESULTS"
"$SPARLAB_ROOT/scripts/run_benchmark.sh" block_3d_analysis "$RESULTS"
for case in cantilever_beam mbb_beam mbb_beam_projected aerospace_bracket wing_rib \
            l_bracket_stress bracket_3d bracket_3d_projected lug_bracket_2d engine_mount_3d \
            bracket_3d_large; do
  "$SPARLAB_ROOT/scripts/run_benchmark.sh" "$case" "$RESULTS"
done

# The stress-constrained deck once more with the constraint switched off: the
# reference the stress table and the comparison figure set beside it.
banner "topology optimization: l_bracket_stress with the stress constraint off"
"$BIN_DIR/sparlab_topopt" --config "$SPARLAB_ROOT/configs/benchmarks/l_bracket_stress.json" \
                          --no-stress --output "$RESULTS/l_bracket_unconstrained"

banner "all benchmarks complete"
echo "results under: $RESULTS"
