#!/bin/bash
# The same native test selection runs in a build container or a heavy Slurm tier.
set -euo pipefail
src=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$src"
if [ -n "${SLURM_JOB_ID:-}" ]; then
    root=$(dirname "$src")
    tier=$(cat "$root/tier.txt")
    . "$root/task.env"
    . /home/ert/opt/kin6d-env.sh
    export PATH=/home/ert/.local/bin:$PATH
else
    tier=${1:?usage: ci/run.sh TIER}
    root=$(mktemp -d "${RUNNER_TEMP:-/tmp}/native-ci.XXXXXX")
    export FC=/usr/bin/gfortran CC=/usr/bin/gcc CXX=/usr/bin/g++
    export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
fi
unset DISPLAY
export CMAKE_TOOLCHAIN_FILE=$src/ci/cluster-toolchain.cmake
if [ -d "$root/mirrors" ]; then
    export GIT_CONFIG_GLOBAL=$root/gitconfig
    : > "$GIT_CONFIG_GLOBAL"
    git config --global protocol.file.allow always
    shopt -s nullglob
    for mirror in "$root"/mirrors/*/*.git; do
        slug=${mirror#"$root/mirrors/"}; slug=${slug%.git}
        git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug.git"
        git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug"
        git config --global --add "url.file://$mirror.insteadOf" "git@github.com:$slug.git"
    done
fi
export GITHUB_WORKSPACE=$src
export GITHUB_ENV=$root/github-env GITHUB_PATH=$root/github-path GITHUB_OUTPUT=$root/github-output
: > "$GITHUB_ENV"; : > "$GITHUB_PATH"; : > "$GITHUB_OUTPUT"
export SCI_CPUS=${SCI_CPUS:-8}
export CMAKE_BUILD_PARALLEL_LEVEL=$SCI_CPUS CTEST_PARALLEL_LEVEL=$SCI_CPUS
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export RUNNER_TEMP=$root/tmp
mkdir -p "$RUNNER_TEMP" ci-out
{ hostname; /usr/bin/gfortran --version | head -1; git rev-parse HEAD; } | tee ci-out/toolchain.txt
load_job_environment() {
    while IFS='=' read -r key value; do
        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] && export "$key=$value"
    done < "$GITHUB_ENV"
    while IFS= read -r item; do [ -z "$item" ] || export PATH="$item:$PATH"; done < "$GITHUB_PATH"
}
if [ -n "${SLURM_JOB_ID:-}" ]; then
mkdir -p "$root/deps"
for lib in fortplot fortio fortnum fortfem fortsparse; do
    variable=KIN6D_$(echo "$lib" | tr 'a-z' 'A-Z')_REVISION
    revision=$(sed -nE "s/^set\\($variable \"?([0-9a-f]{40})\"?.*/\\1/p" CMakeLists.txt)
    test -n "$revision"
    mkdir -p "$root/deps/$lib"
    git --git-dir="$root/mirrors/lazy-fortran/$lib.git" archive "$revision" | tar -x -C "$root/deps/$lib"
done
else
    kin6d_configure() {
        local source=$1 build=$2; shift 3
        cmake -S "$source" -B "$build" -G Ninja -DBLA_VENDOR=OpenBLAS "$@"
    }
fi
case "$tier" in
  cpu-precision-comparison)
    echo 'Build the native CPU benchmark'
    (
cmake -S benchmarks/cpu-precision -B build-cpu-precision -G Ninja -DCMAKE_BUILD_TYPE=Release

    )
    load_job_environment
    echo step-3
    (
cmake --build build-cpu-precision --target benchmark_cpu_precision -j 2

    )
    load_job_environment
    echo 'Compare generated variants and independent analytic values'
    (
ctest --test-dir build-cpu-precision -j 2 -R '^cpu_precision_' --output-on-failure

    )
    load_job_environment
    ;;
  equilibrium-equilibrium)
    export OMP_NUM_THREADS="4"
    export OPENBLAS_NUM_THREADS="1"
    export OMP_WAIT_POLICY="passive"
    export MPLBACKEND="Agg"
    echo Configure
    (
kin6d_configure "$src" build "$root/deps" -DCMAKE_BUILD_TYPE=Release -DKIN6D_BUILD_EXTRAS=ON -DKIN6D_CHECK_LEAN=OFF

    )
    load_job_environment
    echo 'Build affected consumers'
    (
cmake --build build -j 4 --target test_ideal_mhd_coordinates test_ideal_xi test_map_equilibrium test_nested_oracles benchmark_nested_oracles test_nested_boundary_input nested_oracle_export test_curve_chart_solver test_curve_chart test_ideal_map_axisymmetric test_ideal_free_boundary test_landreman_free_boundary test_equilibrium_core test_equilibrium_boundary test_equilibrium_geometry test_gs_readback_free gs_phase1 test_gs_native_adapter test_gs_fb_jacobian test_gs_direct test_mhd_ampere_error test_nonnested_coupled test_nonnested_transport_error test_gs_prescribed_q test_grad_shafranov test_gs_fb_manufactured test_gs_profile_law_readback response_qualify test_response_interfaces test_ideal_seed_3d

    )
    load_job_environment
    echo 'Run affected native tests'
    (
ctest --test-dir build -j 4 --output-on-failure -R '^(ideal_mhd_coordinate_analytic_limits|ideal_xi_hessian_and_analytic_limits|equilibrium_nested_.*|ideal_map_curve_chart_.*|equilibrium_curve_chart_.*|ideal_map_derivatives_and_rotation_null|ideal_map_reduced_omega_derivatives|ideal_map_n3_leakage|ideal_fb_derivatives|ideal_fb_external_kink|ideal_fb_landreman_axisymmetric|equilibrium_core_.*|gs_k6_mini_regression|grad_shafranov_profile_law_readback|response_qualification_external_adapter|response_toroidal_strip_two_ring_traces|response_cut_total_pressure_balance|response_toroidal_strip_edge_merged_domains|grad_shafranov_p[23]_readback_free|gs_native_adapter_contract|gs_free_boundary_t[37]_.*|grad_shafranov_csr_direct_newton_oracles|grad_shafranov_analytic_and_derivative_oracles|ideal_fb_vertical_newton|mhd_ampere_n_error_status|nonnested_coupled_solovev_p2|nonnested_transport_equilibrated_n0|kin6d_layering|ideal_seed_3d_boundary_profiles)$'

    )
    load_job_environment
    ;;
  equilibrium-m6-d8)
    export OMP_NUM_THREADS="4"
    export OPENBLAS_NUM_THREADS="1"
    export OMP_WAIT_POLICY="passive"
    export MPLBACKEND="Agg"
    echo Configure
    (
kin6d_configure "$src" build "$root/deps" -DCMAKE_BUILD_TYPE=Release -DKIN6D_BUILD_EXTRAS=OFF -DKIN6D_CHECK_LEAN=OFF

    )
    load_job_environment
    echo 'Build affected consumers'
    (
cmake --build build -j 4 --target test_equilibrium_harness test_gs_fb_manufactured test_ideal_map_axisymmetric

    )
    load_job_environment
    echo 'Check registry and layering'
    (
ctest --test-dir build -j 4 --output-on-failure -R '^(equilibrium_d8_registry_contract|kin6d_layering)$'

    )
    load_job_environment
    echo 'Check psi ladders and offsets'
    if [ "${SCI_MATRIX_LANE}" = psi ]; then (
ctest --test-dir build -j 4 --output-on-failure -R '^(equilibrium_d8_psi_offset_.*|gs_free_boundary_t3_manufactured_newton)$'

    ); fi
    load_job_environment
    echo 'Check map ladders and offset returns'
    if [ "${SCI_MATRIX_LANE}" = map ]; then (
ctest --test-dir build -j 4 --output-on-failure -R '^(ideal_seed_boundary_label|ideal_map_n0_(solovev_a3_p[23]|shaped_solovev_p[23]|seed_return_(a3_p[23]|shaped_p2)|seed_nested_(a3_p[23]|shaped_p3)))$'

    ); fi
    load_job_environment
    ;;
  *) echo "unknown CI tier: $tier" >&2; exit 2 ;;
esac
