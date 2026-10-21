#!/usr/bin/env bash
# Submit from a dedicated unpacked snapshot; resource requests belong to sbatch.
set -euo pipefail
trap 'printf "[V] failed: %s:%s: %s\n" "${BASH_SOURCE[0]}" "$LINENO" "$BASH_COMMAND" >&2' ERR
unset DISPLAY
export MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_WAIT_POLICY=passive
prototype=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$prototype/../.." && pwd)
run_root=${SLURM_SUBMIT_DIR:?submit from the dedicated snapshot directory}
printf '[V] host=%s job=%s cores=%s\n' "$(hostname)" "$SLURM_JOB_ID" "$SLURM_CPUS_PER_TASK"
cat "$root/REVISION"
cd "$root"
sha256sum --check SHA256SUMS
nvidia-smi
nvcc --version
compiler_args=()
if command -v g++-11 >/dev/null; then
    compiler_args+=(-DCMAKE_CUDA_HOST_COMPILER="$(command -v g++-11)")
fi
cmake -S "$prototype" -B "$run_root/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=75 \
    "${compiler_args[@]}"
cmake --build "$run_root/build" --parallel "$SLURM_CPUS_PER_TASK"
sha256sum "$run_root/build/benchmark_equilibrium_gpu"
ldd "$run_root/build/benchmark_equilibrium_gpu"
ctest --test-dir "$run_root/build" -j "$SLURM_CPUS_PER_TASK" \
    -R '^gpu_equilibrium_leaf_oracles$' --output-on-failure
for threads in 1 2 4 8; do
    for points in 32768 262144; do
        OMP_NUM_THREADS=$threads "$run_root/build/benchmark_equilibrium_gpu" \
            "$points" 21 > "$run_root/gpu-t${threads}-n${points}.log"
    done
done
