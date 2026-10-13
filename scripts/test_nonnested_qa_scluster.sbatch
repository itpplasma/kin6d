#!/bin/bash
#SBATCH -J nn1qa_20261013
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 02:00:00
#SBATCH --mem=32G
#SBATCH -o /home/ert/nn1qa-%j.log
# Mapped NN1 precise-QA vacuum ladders (docs/nn1-qa-benchmark.md) plus the larger
# ellipse ladders with N tied to the mesh. Recipe as test_nonnested_mapped_scluster.sbatch.
# Usage: cd $HOME; sbatch .../test_nonnested_qa_scluster.sbatch KIN6D FORTFEM
set -uo pipefail
K=${1:?kin6d checkout}
F=${2:?fortfem checkout}
cd "$K"
test "$(git -C "$F" rev-parse HEAD)" = \
    "$(sed -n 's/^set(KIN6D_FORTFEM_REVISION "\(.*\)")/\1/p' CMakeLists.txt)" || exit 2
git rev-parse HEAD > nn1qa-revision.txt
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
cmake --build build --target test_nonnested_qa test_nonnested_mapped -j16 > build.log 2>&1 || exit 4
Q="./build/test_nonnested_qa benchmarks/nn1_qa"
M=./build/test_nonnested_mapped
export KIN6D_MAPPED_MIN_ORDER=0
$Q check > qa-check.log 2>&1; echo "check exit $?" >> qa-check.log
/usr/bin/time -v $Q modes 4 4 6 8 10 12 > qa-modes.log 2> qa-modes.time
echo "modes exit $?" >> qa-modes.log
/usr/bin/time -v $Q hladder 10 2 4 8 16 > qa-h.log 2> qa-h.time
echo "hladder exit $?" >> qa-h.log
$Q modes 16 4 6 8 10 > qa-modes16.log 2>&1; echo "modes16 exit $?" >> qa-modes16.log
$Q case 8 10 60 > qa-nphi.log 2>&1
$Q case 8 10 > qa-nphi-ref.log 2>&1
$M hladder 10 4 8 16 > ell-h.log 2>&1; echo "ell hladder exit $?" >> ell-h.log
$M modes 16 4 6 8 10 > ell-modes.log 2>&1; echo "ell modes exit $?" >> ell-modes.log
