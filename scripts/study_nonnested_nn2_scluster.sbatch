#!/bin/bash
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH -o /home/ert/nn2-study-%j.log
# NN2 chi-ratio convergence study and 3D manufactured ladder, run on the
# build of test_nonnested_nn2_scluster.sbatch (submit with afterany).
# The coupled Solov'ev ladders (increment 2) run at three chi/g ratios.
# The island mode (increment 3) runs the finite-pressure island ratio scans
# and D3 ladders on the NN1 helical-coil annulus.
# The landreman mode runs the exact 3D family-1 ladders (per degree K, ratios 1/1e3/1e6,
# N = 4/6 on nr 4/8/16, continuation in the ratio).
# Usage: sbatch scripts/study_nonnested_nn2_scluster.sbatch KIN6D [all|coupled|island|landreman [K...]]
set -euo pipefail
cd "${1:?kin6d checkout}/build"
E=$HOME/fs1/env
export LD_LIBRARY_PATH="$E/lib:$HOME/fsym-env/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=16
git rev-parse HEAD
if test "${2:-all}" = island; then
    run() { /usr/bin/time -v ./test_nonnested_island "$@" \
        > "../nn2-island-$(echo "$*" | tr ' ' _).log" 2>&1 || true; }
    run ratio 3 8 1 1e8
    run ratio 3 16 1 1e8
    run ratio 2 16 1 1e6
    run ratio 3 8 2 1e6
    run ladder 2 4 1 1e4 4
    run ladder 3 4 1 1e6 3
    grep -h -E 'unknowns|flattening|edges|D3 ladder|D3 order|PASS|STOP' \
        ../nn2-island-*.log
    exit 0
fi
if test "${2:-all}" = landreman; then
    # Further arguments: polynomial degrees (default 2 3). Ratios 1, 1e3, 1e6
    # (ratio > 1 by continuation in log ratio, 1 decade per stage, with the
    # 2-decade path check on nr = 4); N = 4 and 6 on nr = 4/8/16.
    shift 2
    run() { /usr/bin/time -v ./test_nonnested_landreman "$@" \
        > "../nn2-landreman-$(echo "$*" | tr ' ' _).log" 2>&1 || true; }
    for K in ${@:-2 3}; do
        for r in 1 1e3 1e6; do
            run "$K" "$r" 3 4
            run "$K" "$r" 3 6
        done
    done
    grep -H -E 'stages|path check|orders|truncation|PASS|STOP|Maximum resident|Elapsed' \
        ../nn2-landreman-*.log
    exit 0
fi
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
