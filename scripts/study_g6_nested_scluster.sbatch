#!/bin/bash
#SBATCH --job-name=g6nested
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
# Submit pilot first, inspect its result/peak, then submit study separately.
# Usage: sbatch -o ROOT/pilot-%j.log SCRIPT ROOT pilot
# ROOT contains source/, deps/{fortnum,fortfem,fortsparse,fortio,fortplot,
# fortgen_runtime}/ and the cluster skill's mem_guard.sh.
set -Eeuo pipefail
root=${1:?run root}
mode=${2:?pilot or study}
source /home/ert/opt/kin6d-env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_DYNAMIC=false
export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
unset DISPLAY
guard=("$root/mem_guard.sh" 47.2)
kin6d_toolchain_report "$root/toolchain-$SLURM_JOB_ID.txt"
case $mode in
    pilot)
        args=()
        for dep in fortnum fortfem fortsparse fortio fortplot fortgen_runtime; do
            args+=("-DFETCHCONTENT_SOURCE_DIR_${dep^^}=$root/deps/$dep")
        done
        "${guard[@]}" cmake --preset cpu -S "$root/source" -B "$root/build" \
            "${KIN6D_CMAKE_ARGS[@]}" "${args[@]}" -DFORTSPARSE_ENABLE_AMD=ON
        "${guard[@]}" cmake --build "$root/build" --target test_nonnested_coupled \
            benchmark_g6_solovev_input -j "$SLURM_CPUS_PER_TASK"
        "${guard[@]}" ctest --test-dir "$root/build" \
            -R '^nonnested_g6_solovev_readback_p3$' --output-on-failure -j4
        "$root/build/benchmark_g6_solovev_input" "$root/input.g6-solovev"
        sha256sum "$root/build/test_nonnested_coupled" > "$root/pilot-binary.sha256"
        ;;
    study)
        sha256sum -c "$root/pilot-binary.sha256"
        # Each command records a failed gate as a failure and stops the study.
        # Three spatial levels and three Fourier levels are separate ladders.
        for nr in 8 16 32; do
            "${guard[@]}" /usr/bin/time -v "$root/build/test_nonnested_coupled" \
                g6 3 1e6 "$nr" 0 > "$root/h-nr$nr.log" 2>&1
        done
        for modes in 0 1 2; do
            "${guard[@]}" /usr/bin/time -v "$root/build/test_nonnested_coupled" \
                g6 3 1e6 8 "$modes" > "$root/m-N$modes.log" 2>&1
        done
        ;;
    *) echo "unknown mode: $mode" >&2; exit 2 ;;
esac
