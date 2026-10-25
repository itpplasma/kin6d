#!/bin/bash
#SBATCH --job-name=k6-vmec-ladder
#SBATCH --cpus-per-task=96
#SBATCH --exclusive
#SBATCH --nodelist=node13
#SBATCH --time=01:00:00
#SBATCH --output=/home/ert/runs/kin6d-gap-baseline-20261025-01/logs/job-%j.out
# VMEC++ 0.9.1 reference ladder only (decision 84). Exclusive EPYC node13, as job 3367225.
# Everything timed runs from node-local /tmp (BeeGFS rule): interpreter, venv, binary, input, outputs.
# Results are copied to $R/out-$SLURM_JOB_ID at the end. Peak RSS < 0.3 GiB per run (exclusive is for timing).
set -uo pipefail
source /home/ert/opt/kin6d-env.sh
R=/home/ert/runs/kin6d-gap-baseline-20261025-01; mkdir -p "$R/logs"
OUT=$R/out-$SLURM_JOB_ID; L=/tmp/k6-vmec-$SLURM_JOB_ID; mkdir -p "$OUT" "$L/out"
OLD=/home/ert/k6bench-20261021/venv; PYHOME=/home/ert/.local/share/uv/python/cpython-3.12-linux-x86_64-gnu
cp -a "$PYHOME" "$L/py"; cp -a "$OLD" "$L/venv"
ln -sfn "$L/py/bin/python3.12" "$L/venv/bin/python"; sed -i "s#^home = .*#home = $L/py/bin#" "$L/venv/pyvenv.cfg"
cp "$R/build-main/benchmark_nested_reference" "$R/vmec_run.py" "$R/helical_mode.py" "$L/"
cp /home/ert/k6bench-20261021/corrected-inputs/landreman-integer/m12/output/input.nested_equilibrium "$L/input"
VPY=$L/venv/bin/python; "$VPY" -c 'import vmecpp,simsopt;print("vmecpp",vmecpp.__file__)' > "$OUT/staging.txt" 2>&1
sha256sum "$L/benchmark_nested_reference" "$L/input" >> "$OUT/staging.txt"; uname -r >> "$OUT/staging.txt"
cpus_for() { python3 -c 'import os,sys; print(",".join(map(str,sorted(os.sched_getaffinity(0))[:int(sys.argv[1])])))' "$1"; }
export OPENBLAS_NUM_THREADS=1 OMP_WAIT_POLICY=passive OMP_PROC_BIND=close OMP_PLACES=cores
# tag mpol NS_ARRAY FTOL_ARRAY NITER_ARRAY
POINTS=(
 "A m5 9,13 1e-6,1e-8 3000,5000"
 "B m5 13,25 1e-6,1e-8 3000,5000"
 "C m8 13,25 1e-6,1e-10 3000,20000"
 "D m8 13,25,49 1e-6,1e-8,1e-10 3000,5000,20000"
 "E m12 9,25 1e-6,1e-10 3000,20000"
 "F m12 13,25 1e-6,1e-10 3000,20000"
 "G m12 13,25,49 1e-6,1e-8,1e-10 3000,5000,20000"
 "H m16 13,25,49 1e-6,1e-8,1e-10 3000,5000,20000"
)
for r in 1 2 3; do for p in "${POINTS[@]}"; do set -- $p; tag=$1; m=${2#m}
  for t in 1 8 16; do
    d=$L/out/$tag-t$t-r$r; mkdir -p "$d"
    ( cd "$L"; OMP_NUM_THREADS=$t taskset -c "$(cpus_for $t)" timeout 900 "$VPY" "$L/vmec_run.py" vmecpp "$L/input" $m "$d" \
        --clockwise --ns-array $3 --ftol-array $4 --niter-array $5 > "$d/solve.log" 2> "$d/err" ) || echo "FAIL vmec $tag t$t r$r" >> "$OUT/failures.txt"
    if [ -f "$d/samples.txt" ]; then
      ( cd "$L"; ./benchmark_nested_reference "$d/samples.txt" 0.5 0 0.015625 "$d/rows.txt" > "$d/compare.log" 2>&1 ); fi
  done; done; done
cp -a "$L/out" "$OUT/vmec"; echo "$(printf '%s\n' "${POINTS[@]}")" > "$OUT/points.txt"; touch "$OUT/DONE"; rm -rf "$L"
