#!/usr/bin/env bash
# Run every benchmark configuration in configs/benchmarks/, plus the small
# static+modal analysis decks (2-D and 3-D), into <results-dir>/<case>/, and
# the comparison runs the documentation sets beside them (the same deck with
# one feature switched off).
#
# The four compliance benchmarks take seconds to a minute each; the
# stress-constrained L-bracket a few tens of seconds; the 3-D bracket
# (15k DOFs, 150 iterations) around ten minutes with the direct solver; the
# Gmsh parts a few minutes; the 356k-DOF bracket_3d_large (multigrid CG)
# roughly half an hour. The buckling-constrained column takes under three
# minutes (without the robust formulation about two); the robust and
# overhang MBB beams about half a minute each, the overhang bracket under
# two minutes.
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
            bracket_3d_large column_buckling mbb_beam_robust mbb_beam_overhang \
            bracket_3d_overhang; do
  "$SPARLAB_ROOT/scripts/run_benchmark.sh" "$case" "$RESULTS"
done

# The stress-constrained deck once more with the constraint switched off: the
# reference the stress table and the comparison figure set beside it.
banner "topology optimization: l_bracket_stress with the stress constraint off"
"$BIN_DIR/sparlab_topopt" --config "$SPARLAB_ROOT/configs/benchmarks/l_bracket_stress.json" \
                          --no-stress --output "$RESULTS/l_bracket_unconstrained"

# The column without the robust formulation, and without the buckling
# constraint: the part the constraint alone delivers, and the compliance
# optimum it moves away from.
comparison() {
  local deck="$1" out="$2"; shift 2
  banner "topology optimization: $deck $* -> $out"
  "$BIN_DIR/sparlab_topopt" --config "$SPARLAB_ROOT/configs/benchmarks/$deck.json" \
                            --output "$RESULTS/$out" "$@"
}
comparison column_buckling column_buckling_nonrobust --no-robust --tag nonrobust
comparison column_buckling column_buckling_unconstrained --no-buckling-constraint \
           --tag unconstrained
# The robust MBB without the robust formulation: mbb_beam_projected plus the
# length-scale and erosion checks.
comparison mbb_beam_robust mbb_beam_robust_off --no-robust --tag robust_off
# The overhang decks without the filter (still checked), and the MBB built
# from its top edge down.
comparison mbb_beam_overhang mbb_beam_overhang_off --no-overhang-filter --tag overhang_off
comparison mbb_beam_overhang mbb_beam_overhang_down --overhang -y --tag build_down
comparison bracket_3d_overhang bracket_3d_overhang_off --no-overhang-filter --tag overhang_off

banner "all benchmarks complete"
echo "results under: $RESULTS"
