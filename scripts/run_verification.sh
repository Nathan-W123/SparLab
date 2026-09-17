#!/usr/bin/env bash
# Run every verification / validation study.
#
# usage: scripts/run_verification.sh [output-dir]
#
# Exits non-zero if any study fails its documented tolerance.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

OUT="${1:-$SPARLAB_ROOT/results/verification}"
require_binaries sparlab_verify
banner "verification and validation studies"
"$BIN_DIR/sparlab_verify" --study all --output "$OUT"
