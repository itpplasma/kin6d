#!/usr/bin/env python3
"""Report observed cross-code errors and refinement differences, not bounds."""
import argparse
import json
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('samples', type=Path)
a = p.parse_args()

def read(name):
    return np.loadtxt(a.samples / name, skiprows=1)

def difference(x, y):
    if x.shape != y.shape or np.max(np.abs(x[:, :3] - y[:, :3])) > 1e-10:
        raise ValueError('Refinement samples are not at the same physical points')
    return float(np.linalg.norm(x[:, 3:] - y[:, 3:]) / np.linalg.norm(y[:, 3:]))

out = {}
for case in ('axis', 'helical'):
    reference = read(f'{case}-M8.samples')
    result = dict(spec_delta_4_to_6=difference(read(f'{case}-M4.samples'),
                                             read(f'{case}-M6.samples')),
                  spec_delta_6_to_8=difference(read(f'{case}-M6.samples'), reference))
    modes = 0 if case == 'axis' else 2
    for rings in (2, 4, 8):
        candidate = read(f'candidate_{case}_{rings}_{modes}.dat')
        result[f'vector_error_r{rings}'] = difference(candidate, reference)
        scale = np.linalg.norm(reference[:, 3:])
        result[f'modB_error_r{rings}'] = float(np.linalg.norm(
            np.linalg.norm(candidate[:, 3:], axis=1) -
            np.linalg.norm(reference[:, 3:], axis=1)) / scale)
    result['mesh_delta_4_to_8'] = difference(
        read(f'candidate_{case}_4_{modes}.dat'), read(f'candidate_{case}_8_{modes}.dat'))
    if case == 'helical':
        for low, high in ((2, 4), (4, 6)):
            result[f'mode_delta_{low}_to_{high}_r4'] = difference(
                read(f'candidate_helical_4_{low}.dat'),
                read(f'candidate_helical_4_{high}.dat'))
    out[case] = result
print(json.dumps(out, indent=2))
