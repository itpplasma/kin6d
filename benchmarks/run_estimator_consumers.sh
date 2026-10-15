#!/usr/bin/env bash
# Build Release with OpenMP MUMPS first. Results and transient files stay in OUT.
set -euo pipefail
build=${1:?build directory}
out=${2:?output directory}
mkdir -p "$out"
export OPENBLAS_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1
export OMP_PROC_BIND=spread OMP_PLACES=cores
for consumer in nn1 nn2 ampere d3; do
    case "$consumer" in
        nn1) sizes='16 32 64';;
        nn2) sizes='8 16 24';;
        ampere) sizes='4 8 12';;
        d3) sizes='48 96 144';;
    esac
    for size in $sizes; do
        for threads in 1 2 4 8 16; do
            OMP_NUM_THREADS=$threads timeout 60 \
                "$build/benchmark_estimator_consumers" "$consumer" "$size" 2 \
                > "$out/$consumer-$size-$threads.log" 2>&1
            echo "$consumer $size $threads complete"
        done
    done
done
