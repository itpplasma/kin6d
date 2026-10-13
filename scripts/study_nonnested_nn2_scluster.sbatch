#!/bin/bash
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 00:45:00
#SBATCH -o /home/ert/nn2-study-%j.log
# NN2 chi-ratio convergence study and 3D manufactured ladder, run on the
# build of test_nonnested_nn2_scluster.sbatch (submit with afterany).
# The coupled Solov'ev ladders (increment 2) run at three chi/g ratios.
# Usage: sbatch scripts/study_nonnested_nn2_scluster.sbatch KIN6D [all|coupled]
set -euo pipefail
cd "${1:?kin6d checkout}/build"
E=$HOME/fs1/env
export LD_LIBRARY_PATH="$E/lib:$HOME/fsym-env/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=16
git rev-parse HEAD
for k in 2 3; do
    for ratio in 1e4 1e6 1e8; do
        /usr/bin/time -v ./test_nonnested_coupled solovev $k $ratio \
            > ../nn2-coupled-p$k-$ratio.log 2> ../nn2-coupled-p$k-$ratio.time || true
    done
done
grep -h -E 'field/p|orders|PASS|STOP' ../nn2-coupled-*.log
test "${2:-all}" = coupled && exit 0
for k in 2 3 4; do
    /usr/bin/time -v ./test_nonnested_transport sovinec $k \
        > ../nn2-sovinec-p$k.log 2> ../nn2-sovinec-p$k.time
done
/usr/bin/time -v ./test_nonnested_transport manufactured \
    > ../nn2-manufactured.log 2> ../nn2-manufactured.time
grep -h -E 'best|observed|PASS' ../nn2-*.log
