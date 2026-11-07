#!/usr/bin/env bash
# Strong-scaling timing of the channel-flow driver from PETSc -log_view.
# Runs a fixed number of RK3 steps on the 80 x 120 x 80 channel for each rank
# count and records wall time of KSPSolve (wall-normal line solves), MatMult
# (FFTW transforms and Galerkin operators) and the whole run.
#
# Usage: benchmarks/run_line_solves.sh [nsteps] [ranks...]
#        e.g. benchmarks/run_line_solves.sh 20 1 2 4 8 16
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exe="${root}/build/bin/channel_flow"
nsteps="${1:-20}"; shift || true
ranks=("${@:-1 2 4 8}")
ranks=(${ranks[@]})
out="${root}/benchmarks/results.csv"
work="$(mktemp -d)"

echo "nprocs,nsteps,ksp_solve_s,mat_mult_s,total_s" > "${out}"
for np in "${ranks[@]}"; do
   log="${work}/log_${np}.txt"
   (cd "${work}" && mpirun -n "${np}" "${exe}" -nsteps "${nsteps}" \
        -output_every 100000000 -log_view > "${log}")
   ksp=$(awk '$1=="KSPSolve"{print $4}' "${log}")
   mm=$(awk '$1=="MatMult"{print $4}' "${log}")
   tot=$(awk '/^Time \(sec\):/{print $3}' "${log}")
   echo "${np},${nsteps},${ksp},${mm},${tot}" | tee -a "${out}"
done
echo "wrote ${out}"
