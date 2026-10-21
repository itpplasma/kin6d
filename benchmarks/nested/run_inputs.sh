#!/bin/bash
# Run inside a Slurm allocation, following CLUSTER_BUILD_RECIPE.md.
set -euo pipefail
test -n "${SLURM_JOB_ID:-}" || { echo 'Run this ladder in Slurm' >&2; exit 2; }
unset DISPLAY
export MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
build=$(realpath "$1")
out=$(realpath -m "$2")
repo=$(cd "$(dirname "$0")/../.." && pwd)
for case in landreman-integer issan-integer; do
    for m in 8 12 16; do
        dest="$out/$case/m$m"
        mkdir -p "$dest"
        sed "/^\//i\ mmax=$m, nmax=$m, angles=64, quadrature=8, turns=1" \
            "$repo/examples/nested-oracles/$case.nml" > "$dest/input.nml"
        (cd "$dest"; /usr/bin/time -v "$build/nested_oracle_export" \
            --config input.nml > export.log 2> time.txt)
    done
done
