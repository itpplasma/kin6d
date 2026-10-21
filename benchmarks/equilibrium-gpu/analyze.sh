#!/usr/bin/env bash
set -euo pipefail
if (($# != 1)); then
    printf 'usage: %s BENCHMARK_LOG\n' "$0" >&2
    exit 1
fi
printf '[V] observed medians; CPU ratios are sensitive to host contention\n'
printf 'action precision CPU/GPU_kernel GPU_kernel_ms copy+kernel_ms scaled_error\n'
awk '$1 ~ /^[01]$/ && NF == 13 {
    printf "[V] %s %s %.3f %.6f %.6f %.6g\n", \
        $1, $2, $4/$7, 1000*$7, 1000*($9+$7+$10), $13
}' "$1"
