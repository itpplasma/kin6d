#!/usr/bin/env bash
# Coupled/refitted gs_free fixture. OUT may be removed after saving numeric evidence.
set -euo pipefail
build=${1:?absolute build directory}
out=${2:?absolute output directory}
root=$(cd "$(dirname "$0")/.." && pwd)
export OPENBLAS_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1
export OMP_PROC_BIND=spread OMP_PLACES=cores
mkdir -p "$out"
for size in 48 72 96; do
    for threads in 1 2 4 8 16; do
        for repeat in 1 2; do
            dir="$out/gs_free-$size-$threads-$repeat"
            mkdir -p "$dir"
            cp "$root"/test/fixtures/gs_free/* "$dir/"
            cat > "$dir/input.nml" <<NML
&gs_free
boundary_file='boundary.dat', boundary_nodes=$size, degree=2,
coil_file='coils.txt', profile_file='profile.dat',
plasma_ellipse=3.0,0.7,1.0, chain_nodes=$((size/2)), lcfs_mode='limiter',
points_file='limiter.dat', axis_guess=3.0,0.0, search_radius=0.5,
seed_amplitude=0.5, max_refits=3, currents=-10000.0,0.0
/
NML
            (cd "$dir" && OMP_NUM_THREADS=$threads timeout 60 \
                "$build/gs_free" --config input.nml > run.log 2>&1)
            ncdump -v estimator_s "$dir/output/gs_free.nc" \
                | sed -n '/data:/,$p' > "$dir/estimator.dat"
            echo "gs_free $size $threads $repeat complete"
        done
    done
done
