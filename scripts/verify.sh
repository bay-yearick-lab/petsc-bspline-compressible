#!/usr/bin/env bash
# Build the solver and run the manufactured-solution checks on 4 MPI ranks.
# Each check prints the max-norm error between the discrete Galerkin result
# and the exact operator applied to the manufactured solution.
#
# Usage: PETSC_DIR=... [PETSC_ARCH=...] scripts/verify.sh [nprocs]
#        (inside the conda env: PETSC_DIR=$CONDA_PREFIX scripts/verify.sh)
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
np="${1:-4}"
build="${root}/build"
export PETSC_DIR="${PETSC_DIR:-${CONDA_PREFIX:-}}"

cmake -S "${root}" -B "${build}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${build}" -j

checks=(linearsolve_x linearsolve_y linearsolve_z
        nonlinear_sum_BBB divergence viscous)
for t in "${checks[@]}"; do
   mpirun -n "${np}" "${build}/bin/test_operators" -test "${t}"
done
