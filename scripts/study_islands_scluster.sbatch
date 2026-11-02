#!/bin/bash
#SBATCH -p compute
#SBATCH -c 8
#SBATCH -t 03:00:00
# Supply --chdir and --output in a dedicated run root containing mem_guard.sh
# from the cluster skill. Build once, then submit studies with afterok:BUILD.
# Usage: sbatch ... SCRIPT CHECKOUT BUILD_ROOT build|vacuum [DEGREE]
set -euo pipefail
checkout=${1:?checkout}
build_root=${2:?build root}
mode=${3:?build or vacuum}
run_root=$PWD
guard="$run_root/mem_guard.sh"
test -x "$guard"
unset DISPLAY
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 OMP_WAIT_POLICY=passive
env_root=$HOME/fs1/env
export PATH="$env_root/bin:$PATH"
export FC="$env_root/bin/x86_64-conda-linux-gnu-gfortran"
export CC="$env_root/bin/x86_64-conda-linux-gnu-gcc"
export LD_LIBRARY_PATH="$env_root/lib:${LD_LIBRARY_PATH:-}"
git -C "$checkout" rev-parse HEAD
git -C "$checkout" status --short
case "$mode" in
    build)
        "$guard" 94.4 cmake -S "$checkout" -B "$build_root" -G Ninja \
            -DCMAKE_BUILD_TYPE=Release -DKIN6D_BUILD_EXTRAS=OFF \
            -DKIN6D_CHECK_LEAN=OFF -DCMAKE_PREFIX_PATH="$HOME/fs1/mumps;$env_root" \
            -DCMAKE_BUILD_RPATH="$env_root/lib" \
            -DCMAKE_EXE_LINKER_FLAGS=-Wl,--enable-new-dtags \
            -DFORTSPARSE_ENABLE_AMD=ON \
            -DFORTSPARSE_AMD_INCLUDE_DIR=/usr/include/suitesparse \
            -DFORTSPARSE_AMD_LIBRARY=/usr/lib/x86_64-linux-gnu/libamd.so
        "$guard" 94.4 cmake --build "$build_root" -j 8 --target \
            test_island_diagnostics test_island_vacuum_oracle test_nonnested_island
        "$guard" 94.4 ctest --test-dir "$build_root" -j 8 --output-on-failure \
            -R '^island_(diagnostics_pendulum|coil_quadrature|vacuum_seed_p3)$'
        ;;
    vacuum)
        degree=${4:?degree}
        sha256sum "$build_root/test_island_vacuum_oracle"
        "$guard" 94.4 "$build_root/test_island_vacuum_oracle" ladder "$degree" 4 4
        ;;
    *) exit 2 ;;
esac
