#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p output
git rev-parse HEAD > output/source-base-commit.txt
git diff --binary HEAD | sha256sum > output/source-worktree-patch.sha256
git ls-files -z -- '*.f90' 'proof/*' 'CMakeLists.txt' 'CMakePresets.json' \
    'scripts/*' 'benchmarks/atlas.json' | xargs -0 sha256sum > output/source-manifest.sha256
cmake --preset cpu
cmake --build --preset cpu
ctest --preset cpu
./build/kin6d_demo "$@"
./build/steady_laboratory_demo
./build/transport_demo
