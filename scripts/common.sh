#!/usr/bin/env bash
# Shared helpers for the SparLab run scripts.
set -euo pipefail

SPARLAB_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$SPARLAB_ROOT/build}"
BIN_DIR="$BUILD_DIR/bin"

require_binaries() {
  local missing=0
  for exe in "$@"; do
    if [[ ! -x "$BIN_DIR/$exe" ]]; then
      echo "error: $BIN_DIR/$exe not found." >&2
      missing=1
    fi
  done
  if [[ $missing -ne 0 ]]; then
    echo "Build first:  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j" >&2
    exit 1
  fi
}

banner() {
  printf '\n=== %s ===\n' "$*"
}
