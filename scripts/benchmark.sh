#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p output/benchmarks
cmake --preset cpu
cmake --build --preset cpu --target benchmark_fp
git rev-parse HEAD > output/benchmarks/source-base-commit.txt
git diff --binary HEAD | sha256sum > output/benchmarks/source-worktree-patch.sha256
sha256sum build/benchmark_fp > output/benchmarks/executable.sha256
uname -a > output/benchmarks/machine.txt
if command -v lscpu >/dev/null; then
    lscpu >> output/benchmarks/machine.txt
fi
rg 'CMAKE_(Fortran_COMPILER|Fortran_FLAGS[A-Z_]*|BUILD_TYPE):' build/CMakeCache.txt \
    > output/benchmarks/compiler-cache.txt
ninja -C build -t commands CMakeFiles/benchmark_fp.dir/app/benchmark_fp.f90.o \
    > output/benchmarks/compile-commands.txt
"$(sed -n 's/^CMAKE_Fortran_COMPILER:FILEPATH=//p' build/CMakeCache.txt)" --version \
    > output/benchmarks/compiler-version.txt
# Explicit audit. Numerical validation belongs to the executable, not timings.
case_sizes=("$@")
if (( ${#case_sizes[@]} == 0 )); then
    case_sizes=(16 32 64 256 4096)
fi
for cells in "${case_sizes[@]}"; do
    if [[ ! "$cells" =~ ^[0-9]+$ ]]; then
        echo 'cell counts must be nonnegative decimal integers' >&2
        exit 2
    fi
    /usr/bin/time -v -o "output/benchmarks/fp-${cells}-resources.txt" \
        ./build/benchmark_fp "$cells" > "output/benchmarks/fp-${cells}.txt"
    cat "output/benchmarks/fp-${cells}.txt"
done
