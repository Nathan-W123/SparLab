#!/usr/bin/env bash
# Regenerate every figure and the animation from the results on disk.
#
# usage: scripts/make_figures.sh [results-dir] [figures-dir]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RESULTS="${1:-$SPARLAB_ROOT/results}"
FIGURES="${2:-$SPARLAB_ROOT/docs/figures}"
PYTHON="${PYTHON:-python3}"

banner "regenerating figures"
cd "$SPARLAB_ROOT"
"$PYTHON" python/scripts/make_all_figures.py --results "$RESULTS" --figures "$FIGURES"
