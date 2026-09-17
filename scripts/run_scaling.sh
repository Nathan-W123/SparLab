#!/usr/bin/env bash
# Runtime scaling benchmark.
#
# usage: scripts/run_scaling.sh [output-dir] [repeats]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

OUT="${1:-$SPARLAB_ROOT/results/benchmark}"
REPEATS="${2:-3}"
require_binaries sparlab_bench
banner "runtime scaling benchmark"
"$BIN_DIR/sparlab_bench" --sizes 20,40,80,120,160,240,320,440 \
                         --repeats "$REPEATS" --output "$OUT"
