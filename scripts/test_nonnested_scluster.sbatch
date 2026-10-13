#!/bin/bash
#SBATCH -p compute
#SBATCH -c 32
#SBATCH -t 03:00:00
#SBATCH -o nn1-%j.log
# Run on sCluster from a private checkout of eq-nonnested-b. This mirrors
# ~/kt2.sbatch: conda gfortran, ~/fs1 MUMPS, FortSym, Release, full CTest.
# Usage: sbatch scripts/test_nonnested_scluster.sbatch /path/to/checkout
set -euo pipefail
cd "${1:-$HOME/kin6d-NN1}"
test "$(git branch --show-current)" = eq-nonnested-b
git rev-parse HEAD > nn1-revision.txt
FS1=$HOME/fs1
E=$FS1/env
export PATH="$E/bin:$PATH"
export FC="$E/bin/x86_64-conda-linux-gnu-gfortran"
export CC="$E/bin/x86_64-conda-linux-gnu-gcc"
export LD_LIBRARY_PATH="$E/lib:$HOME/fsym-env/lib:$HOME/fortsym-diag/build-ci/lib:${LD_LIBRARY_PATH:-}"
cmake --preset cpu -DCMAKE_PREFIX_PATH="$FS1/mumps;$E" \
    -DKIN6D_FORTSYM_BUILD="$HOME/fortsym-diag/build-ci" > cfg.log 2>&1
cmake --build build --target gen_nonnested gen_equilibrium_core -j32 \
    > gen-build.log 2>&1
./build/gen_nonnested
./build/gen_equilibrium_core
cmake --build build -j32 > build.log 2>&1
cd build
ctest -j8 --output-on-failure --timeout 900 > ../ctest.log 2>&1
