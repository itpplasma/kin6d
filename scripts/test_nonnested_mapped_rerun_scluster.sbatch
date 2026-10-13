#!/bin/bash
#SBATCH -J nn1qa_rerun_20261013
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 04:00:00
#SBATCH --mem=48G
#SBATCH -o /home/ert/nn1qa-rerun-%j.log
# Ellipse nr = 16 mode ladder (N = 4 6 8 10) and the full CTest of the branch.
# Recipe as test_nonnested_qa_scluster.sbatch. Usage: cd $HOME; sbatch THIS KIN6D FORTFEM
set -uo pipefail
K=${1:?kin6d checkout}
F=${2:?fortfem checkout}
cd "$K"
test "$(git -C "$F" rev-parse HEAD)" = \
    "$(sed -n 's/^set(KIN6D_FORTFEM_REVISION "\(.*\)")/\1/p' CMakeLists.txt)" || exit 2
cat nn1qa-revision.txt > /dev/null
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
/usr/bin/time -v ./build/test_nonnested_mapped modes 16 4 6 8 10 > ell-modes.log 2> ell-modes.time
echo "ell modes exit $?" >> ell-modes.log
cd build
ctest -j2 --output-on-failure --timeout 1800 > ../ctest.log 2>&1
echo "ctest exit $?" >> ../ctest.log
