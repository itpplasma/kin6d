#!/usr/bin/env python3
"""Run synthetic native SPEC cases and extract common physical sample points."""
import argparse
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

p = argparse.ArgumentParser()
p.add_argument('xspec', type=Path)
p.add_argument('output', type=Path)
a = p.parse_args()
xspec, output = a.xspec.resolve(), a.output.resolve()
output.mkdir(parents=True, exist_ok=True)
recipe = Path(__file__).resolve().parent
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
env.pop('DISPLAY', None)
for case in ('axis', 'helical'):
    for poloidal, radial in ((4, 8), (6, 10), (8, 12)):
        run = output / f'{case}-M{poloidal}-L{radial}'
        run.mkdir(exist_ok=True)
        text = (recipe / f'{case}.sp').read_text()
        text = re.sub(r'Mpol=\d+', f'Mpol={poloidal}', text)
        text = re.sub(r'Lrad=\d+', f'Lrad={radial}', text)
        # SPEC preset.f90:1292: Nt=4*Ndiscrete*Mpol; keep Nt=192.
        text = re.sub(r'Ndiscrete=\d+', f'Ndiscrete={48//poloidal}', text)
        if case == 'helical':
            text = re.sub(r'Ntor=\d+', f'Ntor={poloidal//2}', text)
        (run / f'{case}.sp').write_text(text)
        with (run / 'run.log').open('w') as log:
            subprocess.run([str(xspec), case], cwd=run, env=env,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        assert 'completion' in (run / 'run.log').read_text()
        sample = output / f'{case}-M{poloidal}.samples'
        subprocess.run([sys.executable, str(recipe / 'extract_spec_grid.py'),
                        str(run / f'{case}.sp.h5'), str(sample)], check=True)
    shutil.copyfile(output / f'{case}-M8.samples', output / f'{case}.samples')
(output / 'spec-binary.sha256').write_text(
    hashlib.sha256(xspec.read_bytes()).hexdigest() + '  ' + str(xspec) + '\n')
