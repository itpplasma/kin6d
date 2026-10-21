#!/usr/bin/env bash
# Finite CPU benchmark run; registered CTest comparisons remain authoritative.
set -euo pipefail
build=$(cd "${1:?build directory}" && pwd)
mkdir -p "${2:?output directory}"
out=$(cd "$2" && pwd)
threads=${3:-'1 2 4 8 16'}
size=${4:-131072}
repeats=${5:-15}
trap 'printf "[V] failed at %s:%s (exit %s)\n" "${BASH_SOURCE[0]}" "$LINENO" "$?" > "$out/failure.txt"' ERR
unset DISPLAY
export MPLBACKEND=Agg OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=1 OMP_STACKSIZE=64M
if [[ $(uname -s) == Darwin ]]; then
    # GNU OpenMP on macOS reports that CPU affinity is unsupported.
    unset OMP_PROC_BIND OMP_PLACES
    timer=(/usr/bin/time -l)
else
    export OMP_PROC_BIND=spread OMP_PLACES=cores
    timer=(/usr/bin/time -v)
fi
ctest --test-dir "$build" -j 2 -R '^cpu_precision_' --output-on-failure \
    > "$out/tests.log" 2>&1
for t in $threads; do
    for kernel in gs map; do
        OMP_NUM_THREADS=$t "${timer[@]}" "$build/benchmark_cpu_precision" \
            "$kernel" "$size" "$repeats" > "$out/$kernel-$t.log" 2>&1
    done
done
{
    uname -a
    if command -v sha256sum >/dev/null; then
        sha256sum "$build/benchmark_cpu_precision"
    else
        shasum -a 256 "$build/benchmark_cpu_precision"
    fi
    rg 'CMAKE_Fortran_(COMPILER:|FLAGS:|FLAGS_RELEASE:)' "$build/CMakeCache.txt" || true
} > "$out/provenance.txt"
