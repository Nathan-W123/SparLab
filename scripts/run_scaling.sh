#!/usr/bin/env bash
# Runtime scaling benchmarks: the Q4 plate (runtime_scaling.*) and the Hex8
# block (runtime_scaling_3d.*).
#
# usage: scripts/run_scaling.sh [output-dir] [repeats]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

OUT="${1:-$SPARLAB_ROOT/results/benchmark}"
REPEATS="${2:-3}"
require_binaries sparlab_bench
banner "runtime scaling benchmark (2-D, Q4)"
"$BIN_DIR/sparlab_bench" --sizes 20,40,80,120,160,240,320,440 \
                         --repeats "$REPEATS" --output "$OUT"
banner "runtime scaling benchmark (3-D, Hex8)"
"$BIN_DIR/sparlab_bench" --dim 3 --sizes 8,16,24,32,40 \
                         --repeats "$REPEATS" --output "$OUT"
