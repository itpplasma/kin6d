"""b(theta, zeta) on nt x nz points (theta fastest) of one surface from the Boozer spectrum."""
import json, sys
import numpy as np
src, idx, nt, nz, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
r = json.load(open(src))[idx]
th = 2*np.pi*np.arange(nt)/nt
ze = 2*np.pi*np.arange(nz)/(nz*r['nfp'])
T, Z = np.meshgrid(th, ze, indexing='ij')
B = np.zeros_like(T)
for m, n, c, s in zip(r['xm'], r['xn'], r['bmnc'], r['bmns']):
    B += c*np.cos(m*T - n*Z) + s*np.sin(m*T - n*Z)
b0 = r['bmnc'][0]
with open(out, 'w') as f:
    f.write(f"{nt} {nz} {r['nfp']}\n{r['s']!r} {r['iota']!r} {r['G']!r} {r['I']!r} {b0!r}\n")
    for v in (B/b0).flatten(order='F'):
        f.write(f"{float(v)!r}\n")
