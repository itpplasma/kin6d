#!/bin/bash
#SBATCH -p compute
#SBATCH -J hwscan
#SBATCH -t 16:00:00
#SBATCH -o hwscan_%A_%a.log
# One Hasegawa-Wakatani case of cases.tsv per array task. Submit one array per
# thread count (column 14): sbatch -c THREADS --array=IDS scan.sbatch, with
# HW_ROOT holding cases.tsv and the executable hasegawa_wakatani.
set -euo pipefail
. "${MILCA_ENV:-$HOME/opt/milca-env.sh}"
root=${HW_ROOT:?set HW_ROOT to the run directory}
line=$(grep -v '^#' $root/cases.tsv | awk -F'\t' -v id=$SLURM_ARRAY_TASK_ID '$1==id')
IFS=$'\t' read -r id name alpha kappa d order modified n dt tend tavg ss seed threads <<< "$line"
dir=$root/runs/$name; mkdir -p $dir; cd $dir
cat > case.nml <<EOF
&hasegawa_wakatani
 alpha = $alpha, kappa = $kappa, d = $d, order = $order, modified = $modified,
 nx = $n, ny = $n, k0 = 0.15, dt = $dt, time_end = $tend, average_from = $tavg,
 sample_steps = $ss, segment = 512, modes = 16, amplitude = 1.0e-2, seed = $seed /
EOF
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=1 OMP_WAIT_POLICY=passive
export OMP_PROC_BIND=close OMP_PLACES=cores
hostname > node.txt
/usr/bin/time -f "%e s wall, %M KB peak" $root/hasegawa_wakatani --config case.nml > run.log 2> run.err
cat run.log; tail -2 run.err
