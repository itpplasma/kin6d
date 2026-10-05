#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
cmake --preset cpu
cmake --build --preset cpu
ctest --preset cpu
./build/kin6d_demo "$@"
