#!/bin/bash
#SBATCH -p compute
#SBATCH -c 32
#SBATCH -t 03:00:00
#SBATCH -o /home/ert/nn2-ctest-%j.log
# NN2 full CTest on sCluster from private checkouts. Mirrors the NN1 recipe
# (~/kt2.sbatch: conda gfortran, ~/fs1 MUMPS, pinned FortSym, Release).
# The FortFEM pin may be an unpushed commit, so its checkout is passed.
# Usage: sbatch scripts/test_nonnested_nn2_scluster.sbatch KIN6D FORTFEM
set -euo pipefail
K=${1:?kin6d checkout}
F=${2:?fortfem checkout}
cd "$K"
test "$(git -C "$F" rev-parse HEAD)" = \
    "$(sed -n 's/^set(KIN6D_FORTFEM_REVISION "\(.*\)")/\1/p' CMakeLists.txt)"
git rev-parse HEAD > nn2-revision.txt
FS1=$HOME/fs1
E=$FS1/env
export PATH="$E/bin:$PATH"
export FC="$E/bin/x86_64-conda-linux-gnu-gfortran"
export CC="$E/bin/x86_64-conda-linux-gnu-gcc"
export LD_LIBRARY_PATH="$E/lib:$HOME/fsym-env/lib:$HOME/fortsym-diag/build-ci/lib:${LD_LIBRARY_PATH:-}"
cmake --preset cpu -DCMAKE_PREFIX_PATH="$FS1/mumps;$E" \
    -DKIN6D_FORTSYM_BUILD="$HOME/fortsym-diag/build-ci" \
    -DFETCHCONTENT_SOURCE_DIR_FORTFEM="$F" > cfg.log 2>&1
cmake --build build --target gen_nonnested -j32 > gen-build.log 2>&1
./build/gen_nonnested
git diff --exit-code src/generated > generated-drift.log || echo DRIFT >> generated-drift.log
cmake --build build -j32 > build.log 2>&1
cd build
ctest -j8 --output-on-failure --timeout 900 > ../ctest.log 2>&1
