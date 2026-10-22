#!/usr/bin/env bash
# Focused NN2 profiling. Keep logs outside Git. No GUI or private inputs.
# Usage: study_nonnested_scaling.sh BUILD RESULTS [local|large|krylov|threads]
set -euo pipefail
build=$(realpath "${1:?build directory}")
results=$(realpath -m "${2:?results directory}")
scope=${3:-local}
mkdir -p "$results"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg KIN6D_PROFILE=1
unset DISPLAY
run() {
    local tag=$1
    shift
    /usr/bin/time -v "$build/test_nonnested_island" "$@" \
        > "$results/$tag.log" 2>&1
}
case "$scope" in
    local|large)
        for k in 2 3; do
            for nr in 8 16; do
                for n in 1 2 4; do
                    # Large all-pair assembly goes to a guarded Slurm job.
                    if [[ $scope == local && $nr == 16 && $n != 1 ]]; then
                        continue
                    fi
                    tag="dense-p$k-r$nr-n$n"
                    status=0
                    run "$tag" profile "$k" "$nr" "$n" 1e2 dense || status=$?
                    printf '%s\t%d\n' "$tag" "$status" >> "$results/status.tsv"
                done
            done
        done
        ;;
    krylov)
        for n in 1 2 4 8 16 20; do
            for ratio in 1e2 1e3 1e4 1e5 1e6; do
                # A failed solve is evidence, retained with its exit status.
                tag="matrixfree-p2-r8-n$n-chi$ratio"
                status=0
                run "$tag" profile 2 8 "$n" "$ratio" matrixfree 30 100 || status=$?
                printf '%s\t%d\n' "$tag" "$status" >> "$results/status.tsv"
            done
        done
        ;;
    threads)
        for threads in 1 2 4 8; do
            export OMP_NUM_THREADS=$threads
            tag="matrixfree-p2-r8-n2-t$threads"
            status=0
            run "$tag" profile 2 8 2 1e2 matrixfree || status=$?
            printf '%s\t%d\n' "$tag" "$status" >> "$results/status.tsv"
        done
        ;;
    *) echo "unknown scope: $scope" >&2; exit 2 ;;
esac
