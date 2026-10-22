#!/usr/bin/env python3
"""Summarize NN2 phase samples, preserving inclusive phase names and failures."""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('results', type=Path)
args = parser.parse_args()
out = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
out.writerow(['case', 'phase', 'calls', 'seconds', 'max_end_rss_KiB',
              'max_process_hwm_KiB', 'solve_status', 'krylov_iterations',
              'krylov_status', 'process_exit'])
for path in sorted(args.results.glob('*.log')):
    phases = defaultdict(lambda: [0, 0.0, 0, 0])
    solve_status = 'unavailable'
    krylov = ''
    krylov_status = ''
    process_exit = 'unavailable'
    elapsed = 0.0
    peak = 0
    for line in path.read_text(errors='replace').splitlines():
        match = re.match(r'PROF (nn2-[\w-]+)\s+([\d.]+)$', line)
        if match and match[1] != 'nn2-jacobian-nnz':
            entry = phases[match[1]]
            entry[0] += 1
            entry[1] += float(match[2])
        match = re.match(r'MEM (\S+) rss_KiB=(\d+) hwm_KiB=(\d+)', line)
        if match:
            entry = phases[match[1]]
            entry[2] = max(entry[2], int(match[2]))
            entry[3] = max(entry[3], int(match[3]))
        match = re.match(r'solve status\s+(\d+)', line)
        if match:
            solve_status = match[1]
        if line.startswith('Krylov status'):
            krylov_status = ','.join(line.split()[2:])
        if 'Maximum resident set size' in line:
            peak = int(line.split()[-1])
        if 'Elapsed (wall clock)' in line:
            parts = line.rsplit(': ', 1)[1].split(':')
            elapsed = sum(float(x)*60**i for i, x in enumerate(reversed(parts)))
        if 'Exit status:' in line:
            process_exit = line.split()[-1]
        if line.startswith('Krylov iterations'):

            krylov = ','.join(line.split()[2:])
    phases['process-total'] = [1, elapsed, 0, peak]
    for phase, entry in sorted(phases.items()):
        out.writerow([path.stem, phase, entry[0], f'{entry[1]:.3f}',
                      entry[2], entry[3], solve_status, krylov, krylov_status, process_exit])
