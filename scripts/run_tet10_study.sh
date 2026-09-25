#!/usr/bin/env bash
# Linear against quadratic tetrahedra on the engine mount meshed from CAD.
#
# usage: scripts/run_tet10_study.sh [results-dir]
#
# Meshes the part with Gmsh at 8, 6, 4 and 3 mm with linear and with curved
# quadratic tetrahedra (plus 2 and 1.5 mm linear meshes), solves each with
# sparlab_solve - Tet4, the linear mesh elevated to straight-sided Tet10, and
# the curved Tet10 mesh - and tabulates the compliances against the finest
# curved Tet10 run and a Richardson estimate. Needs the gmsh Python package;
# about four minutes on four cores.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RESULTS="${1:-$SPARLAB_ROOT/results}"
PYTHON="${PYTHON:-python3}"
require_binaries sparlab_solve
banner "Tet4 / Tet10 study on the engine mount"
cd "$SPARLAB_ROOT"
"$PYTHON" python/scripts/tet10_part_study.py --output "$RESULTS/tet10_part_study" \
                                             --bin "$BIN_DIR"
