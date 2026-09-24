#!/usr/bin/env bash
# Runtime scaling benchmarks.
#
#   runtime_scaling.*              Q4 plate, sparse Cholesky (LDL^T)
#   runtime_scaling_amg_cg.*       Q4 plate, multigrid-preconditioned CG
#   runtime_scaling_3d.*           Hex8 block, sparse Cholesky
#   runtime_scaling_3d_amg_cg.*    Hex8 block, multigrid CG, to 830k DOFs
#   runtime_scaling_3d_conjugate_gradient.*  Hex8 block, Jacobi CG
#   runtime_scaling_3d_tet.*       Tet4 block, sparse Cholesky
#   runtime_scaling_3d_tet_amg_cg.*  Tet4 block, multigrid CG
#
# The iterative solvers run to a relative residual of 1e-10. The Jacobi CG
# series skips the objective + gradient timing (--no-objective): it is there
# to show how its iteration count grows, not as an optimisation solver.
#
# usage: scripts/run_scaling.sh [output-dir] [repeats]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

OUT="${1:-$SPARLAB_ROOT/results/benchmark}"
REPEATS="${2:-3}"
require_binaries sparlab_bench
B="$BIN_DIR/sparlab_bench"

banner "runtime scaling (2-D, Q4): sparse Cholesky"
"$B" --sizes 20,40,80,120,160,240,320,440 --repeats "$REPEATS" --output "$OUT"
banner "runtime scaling (2-D, Q4): multigrid CG"
"$B" --sizes 20,40,80,160,320,440,640,1000 --repeats 2 --solver amg_cg --output "$OUT"

banner "runtime scaling (3-D, Hex8): sparse Cholesky"
"$B" --dim 3 --sizes 8,16,24,32,40,48 --repeats 2 --output "$OUT"
banner "runtime scaling (3-D, Hex8): multigrid CG"
"$B" --dim 3 --sizes 8,16,24,32,40,48,64,80,96,128 --repeats 2 --solver amg_cg \
     --output "$OUT"
banner "runtime scaling (3-D, Hex8): Jacobi-preconditioned CG"
"$B" --dim 3 --sizes 8,16,24,32,40,48,64 --repeats 2 --solver conjugate_gradient \
     --no-objective --output "$OUT"

banner "runtime scaling (3-D, Tet4): sparse Cholesky"
"$B" --dim 3 --element tet --sizes 8,12,16,24,32,40 --repeats 2 --output "$OUT"
banner "runtime scaling (3-D, Tet4): multigrid CG"
"$B" --dim 3 --element tet --sizes 8,12,16,24,32,48,64 --repeats 2 --solver amg_cg \
     --output "$OUT"
