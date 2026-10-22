#!/bin/bash
# IH2 chart pilot/sweep on sCluster. Usage (from a clone at the recorded revision):
#   RUN=$HOME/runs/kin6d-c1-20261022 sbatch --export=ALL,CASES=cases.txt,TAG=pilot \
#       benchmarks/c1-charts/c1.sbatch
# CASES lines: CHART P N NPHI CHI_RATIO [TOL]; MAP is the fitted nested map.
# BUILD=1 (re)builds and tests first; sweep jobs depend on one such build job.
#SBATCH --job-name=kin6d-c1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=04:00:00
set -euo pipefail
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-16} OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
export OMP_PROC_BIND=close OMP_PLACES=cores
unset DISPLAY
src=$(git -C "$SLURM_SUBMIT_DIR" rev-parse --show-toplevel)/benchmarks/c1-charts
build=$RUN/build
out=$RUN/$TAG-$SLURM_JOB_ID.csv
printf 'HOST=%s JOB=%s KIN6D=%s THREADS=%s\n' "$(hostname)" "$SLURM_JOB_ID" \
    "$(git -C "$src" rev-parse HEAD)" "$OMP_NUM_THREADS"
if [ "${BUILD:-0}" = 1 ]; then
    cmake -S "$src" -B "$build" -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_Fortran_FLAGS='-O3 -march=native' -DKIN6D_C1_MAP="$MAP"
    cmake --build "$build" -j "$SLURM_CPUS_PER_TASK"
    OMP_NUM_THREADS=4 ctest --test-dir "$build" -L c1 --output-on-failure
fi
"$build/c1_charts" --header > "$out"
while read -r chart p n nphi chi tol; do
    [ -z "$chart" ] && continue
    case "$chart" in \#*) continue ;; esac
    /usr/bin/time -f "rss_kib=%M wall=%e" "$build/c1_charts" "$MAP" "$chart" "$p" "$n" \
        "$nphi" "$chi" "${tol:-1e-12}" >> "$out" 2> "$out.$chart-$p-$n-$nphi-$chi.time" \
        || echo "$chart,$p,$n,$nphi,NA,NA,$chi,${tol:-1e-12},FAILED" >> "$out"
done < "$CASES"
echo "done: $out"
