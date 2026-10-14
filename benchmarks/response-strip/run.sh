#!/usr/bin/env bash
# Run after the native Release oracle gates. Results stay outside source history.
set -euo pipefail
strip_evidence=${1:?usage: run.sh ABSOLUTE_EVIDENCE_DIRECTORY}
mkdir -p "$strip_evidence"
cmake --build --preset cpu --target test_response_strip_preconditioner -j 6
ctest --test-dir build -R '^response_strip_preconditioner_.*$' --output-on-failure
lscpu > "$strip_evidence/hardware.txt"
for strip_threads in 1 2 4; do
    for strip_trial in 1 2 3; do
        OMP_NUM_THREADS="$strip_threads" OMP_DYNAMIC=false OMP_WAIT_POLICY=passive \
            OPENBLAS_NUM_THREADS=1 /usr/bin/time -v \
            build/test_response_strip_preconditioner scaling \
            > "$strip_evidence/threads_${strip_threads}_trial_${strip_trial}.log" \
            2> "$strip_evidence/threads_${strip_threads}_trial_${strip_trial}.resources"
    done
done
