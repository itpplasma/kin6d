#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -c 16
#SBATCH -t 08:00:00
#SBATCH -J nn2-v3scale
# Usage: sbatch --output=RUN/slurm-%j.log THIS SOURCE RUN MEM_GUARD [all|krylov]
# 16 cores reserve 188.8 GiB under sCluster CR_CORE. Run on eight threads,
# leaving the extra cores idle to reserve memory for dense assembly/fill.
set -euo pipefail
source_dir=$(realpath "${1:?source export or worktree}")
run_root=$(realpath -m "${2:?run root}")
guard=$(realpath "${3:?cluster skill mem_guard.sh}")
scope=${4:-all}
cpus=${SLURM_CPUS_PER_TASK:-16}
guard_limit=$(awk -v cores="$cpus" 'BEGIN {printf "%.1f", cores*11.8}')
mkdir -p "$run_root"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
unset DISPLAY
if [[ -x $HOME/fs1/env/bin/x86_64-conda-linux-gnu-gfortran ]]; then
    export FC=$HOME/fs1/env/bin/x86_64-conda-linux-gnu-gfortran
    export CC=/usr/bin/gcc
    export PATH=$HOME/fs1/env/bin:$PATH
    export LD_LIBRARY_PATH=$HOME/fs1/env/lib:${LD_LIBRARY_PATH:-}
fi
cmake --preset dev-fast -S "$source_dir" -B "$run_root/build" \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DKIN6D_SHARED_DEPS_DIR="$run_root/deps" > "$run_root/configure.log" 2>&1
"$guard" "$guard_limit" cmake --build "$run_root/build" --target \
    test_nonnested_island test_nonnested_coupled test_nonnested_leaves -j8 > "$run_root/build.log" 2>&1
# One job, with a completed small-case probe before the larger cases.
"$guard" "$guard_limit" /usr/bin/time -v "$run_root/build/test_nonnested_island" \
    profile 3 8 1 1e2 dense > "$run_root/probe.log" 2>&1
if [[ $scope == all ]]; then
    "$guard" "$guard_limit" "$source_dir/scripts/study_nonnested_scaling.sh" \
        "$run_root/build" "$run_root/large" large
fi
"$guard" "$guard_limit" "$source_dir/scripts/study_nonnested_scaling.sh" \
    "$run_root/build" "$run_root/krylov" krylov

# High-anisotropy equivalence is a gate, never silently relaxed.
mkdir -p "$run_root/equivalence"
for k in 2 3; do
    for n in 1 2; do
        status=0
        "$guard" "$guard_limit" /usr/bin/time -v "$run_root/build/test_nonnested_island" \
            compare "$k" 8 "$n" 1e6 \
            > "$run_root/equivalence/p$k-n$n.log" 2>&1 || status=$?
        printf 'p%d-n%d\t%d\n' "$k" "$n" "$status" \
            >> "$run_root/equivalence/status.tsv"
    done
done
