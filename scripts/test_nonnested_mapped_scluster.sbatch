#!/bin/bash
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 05:00:00
#SBATCH --mem=48G
#SBATCH -o /home/ert/nn1map-%j.log
# NN1 mapped 3D: large ellipse ladders, then the full CTest, on sCluster.
# Recipe as scripts/test_nonnested_nn2_scluster.sbatch. The FortFEM pin is passed
# as a checkout. Usage: cd $HOME; sbatch .../test_nonnested_mapped_scluster.sbatch KIN6D FORTFEM
set -uo pipefail
K=${1:?kin6d checkout}
F=${2:?fortfem checkout}
cd "$K"
test "$(git -C "$F" rev-parse HEAD)" = \
    "$(sed -n 's/^set(KIN6D_FORTFEM_REVISION "\(.*\)")/\1/p' CMakeLists.txt)" || exit 2
git rev-parse HEAD > nn1map-revision.txt
FS1=$HOME/fs1
E=$FS1/env
export PATH="$E/bin:$PATH"
export FC="$E/bin/x86_64-conda-linux-gnu-gfortran"
export CC="$E/bin/x86_64-conda-linux-gnu-gcc"
export LD_LIBRARY_PATH="$E/lib:$HOME/fsym-env/lib:$HOME/fortsym-diag/build-ci/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=spread OMP_DYNAMIC=false
cmake --preset cpu -DCMAKE_PREFIX_PATH="$FS1/mumps;$E" \
    -DKIN6D_FORTSYM_BUILD="$HOME/fortsym-diag/build-ci" \
    -DFETCHCONTENT_SOURCE_DIR_FORTFEM="$F" > cfg.log 2>&1 || exit 3
cmake --build build -j16 > build.log 2>&1 || exit 4
T=./build/test_nonnested_mapped
/usr/bin/time -v $T hladder 6 2 4 8 16 > ladder-h.log 2> ladder-h.time
echo "hladder exit $?" >> ladder-h.log
/usr/bin/time -v $T modes 8 1 2 3 4 6 8 > ladder-modes.log 2> ladder-modes.time
echo "modes exit $?" >> ladder-modes.log
cd build
ctest -j4 --output-on-failure --timeout 900 > ../ctest.log 2>&1
