#!/usr/bin/env python3
"""Extract native SPEC grid samples; this optional adapter does not solve a field."""
import argparse
from pathlib import Path
import h5py
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('h5')
p.add_argument('samples')
a = p.parse_args()
with h5py.File(a.h5) as f:
    g = f['grid']
    nt, nz = int(g['Nt'][0]), int(g['Nz'][0])
    period = float(g['pi2nfp'][0])
    r, z = np.asarray(g['Rij']), np.asarray(g['Zij'])
    br, bphi, bz = (np.asarray(g[k]) for k in ('BR', 'Bp', 'BZ'))
    assert nt == 192 and nz in (1, 96), (nt, nz)
    # HDF5 reverses the Fortran (radial, angular) dimensions.
    rows = []
    for radial in (3, 7, 11):
        for plane in ((0,) if nz == 1 else (1, 11, 25, 37, 49, 61, 77, 89)):
            for theta in (5, 23, 47, 71, 97, 121, 149, 173):
                angular = (theta % nt) + plane*nt
                rr = r[angular, radial]
                # SPEC sphdf5.f90:662-666 writes B^phi, so physical B_phi=R B^phi.
                rows.append([rr, z[angular, radial], plane*period/nz,
                             br[angular, radial], rr*bphi[angular, radial],
                             bz[angular, radial]])
np.savetxt(a.samples, rows, header=str(len(rows)), comments='')
print(f'{len(rows)} native SPEC samples: {Path(a.samples)}')
