"""Boozer spectra of the half-grid surfaces nearest to the requested s.

Usage: python boozer_surfaces.py WOUT OUT.json S1,S2,... MBOZ NBOZ
(booz_xform 0.1.0; LASYM kept; nboz > 0 keeps the toroidal harmonics).
"""
import sys, json
import numpy as np
import booz_xform
from netCDF4 import Dataset
wout = sys.argv[1]; out = sys.argv[2]; svals = [float(x) for x in sys.argv[3].split(',')]
mb, nb = int(sys.argv[4]), int(sys.argv[5])
with Dataset(wout) as w:
    ns = int(w['ns'][:])
idx = sorted({int(round(s*(ns - 1) - 0.5)) + 1 for s in svals})  # half-grid index (1-based)
bx = booz_xform.Booz_xform()
bx.read_wout(wout, flux=True)
bx.mboz, bx.nboz = mb, nb
bx.compute_surfs = [i - 1 for i in idx]
bx.verbose = 0
bx.run()
res = []
for k, js in enumerate(bx.compute_surfs):
    s = (js + 1 - 0.5)/(ns - 1)
    xm = np.array(bx.xm_b); xn = np.array(bx.xn_b)
    bmnc = np.array(bx.bmnc_b)[:, k]
    bmns = np.array(bx.bmns_b)[:, k] if bx.asym else 0*bmnc
    res.append(dict(s=s, iota=float(bx.iota[js]), G=float(bx.Boozer_G[k]),
                    I=float(bx.Boozer_I[k]), nfp=int(bx.nfp), xm=xm.tolist(), xn=xn.tolist(),
                    bmnc=bmnc.tolist(), bmns=bmns.tolist()))
json.dump(res, open(out, 'w'))
for r in res:
    print(f"s={r['s']:.4f} iota={r['iota']:.5f} G={r['G']:.4f} I={r['I']:.5f} B00={r['bmnc'][0]:.4f}")
