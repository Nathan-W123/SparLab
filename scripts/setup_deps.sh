#!/usr/bin/env bash
# Install SparLab's dependencies.
#
# usage: scripts/setup_deps.sh [--python-only] [--system-only]
#
# System packages (Debian/Ubuntu; needs sudo):
#   build-essential cmake ninja-build libeigen3-dev catch2
#
# Python packages for the visualisation layer:
#   numpy pandas matplotlib pillow
#
# On other platforms install the equivalents:
#   Fedora/RHEL   dnf install gcc-c++ cmake ninja-build eigen3-devel catch2-devel
#   macOS         brew install cmake ninja eigen catch2
#   Windows       vcpkg install eigen3 catch2
#
# Catch2 is optional: if the system package is absent, CMake fetches v3.5.2 with
# FetchContent, so a machine with network access needs only Eigen.
set -euo pipefail

PYTHON_ONLY=0
SYSTEM_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --python-only) PYTHON_ONLY=1 ;;
    --system-only) SYSTEM_ONLY=1 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option '$arg'" >&2; exit 2 ;;
  esac
done

SUDO=""
if [[ $EUID -ne 0 ]] && command -v sudo >/dev/null 2>&1; then SUDO="sudo"; fi

if [[ $PYTHON_ONLY -eq 0 ]]; then
  if command -v apt-get >/dev/null 2>&1; then
    echo "--- installing system packages with apt-get"
    $SUDO apt-get update -qq
    $SUDO apt-get install -y build-essential cmake ninja-build libeigen3-dev catch2
  else
    echo "apt-get not found; install the equivalents for your platform:" >&2
    sed -n '/On other platforms/,/vcpkg/p' "$0" >&2
    [[ $SYSTEM_ONLY -eq 1 ]] && exit 1
  fi
fi

if [[ $SYSTEM_ONLY -eq 0 ]]; then
  echo "--- installing Python packages"
  PY="${PYTHON:-python3}"
  "$PY" -m pip install --upgrade pip >/dev/null
  "$PY" -m pip install numpy pandas matplotlib pillow
fi

echo "--- done. Next:  make build && make test"
