#!/bin/bash
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 00:45:00
#SBATCH -o /home/ert/nn2-study-%j.log
# NN2 chi-ratio convergence study and 3D manufactured ladder, run on the
# build of test_nonnested_nn2_scluster.sbatch (submit with afterany).
# Usage: sbatch scripts/study_nonnested_nn2_scluster.sbatch KIN6D
set -euo pipefail
cd "${1:?kin6d checkout}/build"
E=$HOME/fs1/env
export LD_LIBRARY_PATH="$E/lib:$HOME/fsym-env/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=16
git rev-parse HEAD
for k in 2 3 4; do
    /usr/bin/time -v ./test_nonnested_transport sovinec $k \
        > ../nn2-sovinec-p$k.log 2> ../nn2-sovinec-p$k.time
done
/usr/bin/time -v ./test_nonnested_transport manufactured \
    > ../nn2-manufactured.log 2> ../nn2-manufactured.time
grep -h -E 'best|observed|PASS' ../nn2-*.log
