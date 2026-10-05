#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
lake build
lake env lean tests/Accepted.lean
if [[ $# -gt 0 ]]; then
    mkdir -p ../output
    "$1" ../output/ComputedKrook.lean
    lake env lean ../output/ComputedKrook.lean
fi
for rejected in BadResidual BadRadius; do
    output=$(mktemp)
    if lake env lean "tests/$rejected.lean" >"$output" 2>&1; then
        cat -- "$output"
        rm -f -- "$output"
        echo "ERROR: corrupt certificate $rejected was accepted" >&2
        exit 1
    fi
    if ! rg -q 'unsolved goals' "$output"; then
        cat -- "$output"
        rm -f -- "$output"
        echo "ERROR: $rejected failed for an unexpected reason" >&2
        exit 1
    fi
    rm -f -- "$output"
    echo "Rejected corrupt certificate: $rejected"
done
