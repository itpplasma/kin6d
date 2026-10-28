"""Compare a KIN6D Gamma_p (<prefix>.lcfs, <prefix>.json) with the VMEC++ LCFS.

For every KIN6D boundary point X the distance to the VMEC++ LCFS surface is
the minimum over a fine (theta, phi) sample refined by Newton on the squared
distance; |B| is compared at the foot point (VMEC++ bmnc on the outermost half
grid surface, extrapolated half a cell to s = 1). Prints one JSON line.
"""
import argparse
import json

import netCDF4
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("kin6d_prefix")
ap.add_argument("wout")
a = ap.parse_args()

w = netCDF4.Dataset(a.wout)
xm, xn = w["xm"][:], w["xn"][:]
rmnc, zmns = w["rmnc"][:][-1], w["zmns"][:][-1]
xmb, xnb = w["xm_nyq"][:], w["xn_nyq"][:]
bmnc = w["bmnc"][:]
b_edge = 1.5 * bmnc[-1] - 0.5 * bmnc[-2]


def surface(theta, phi):
    ang = np.multiply.outer(theta, xm) - np.multiply.outer(phi, xn)
    r = (rmnc * np.cos(ang)).sum(-1)
    z = (zmns * np.sin(ang)).sum(-1)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], -1)


def modb(theta, phi):
    ang = np.multiply.outer(theta, xmb) - np.multiply.outer(phi, xnb)
    return (b_edge * np.cos(ang)).sum(-1)


th, ph = np.meshgrid(np.linspace(0, 2 * np.pi, 361)[:-1],
                     np.linspace(0, 2 * np.pi, 721)[:-1], indexing="ij")
pts = surface(th, ph).reshape(-1, 3)
tree = cKDTree(pts)
k = np.loadtxt(a.kin6d_prefix + ".lcfs")
x, bk = k[:, 1:4], k[:, 4]
_, idx = tree.query(x)
dist, berr = [], []
for p, i in zip(x, idx):
    t0, p0 = th.ravel()[i], ph.ravel()[i]
    f = lambda v: np.sum((surface(np.array(v[0]), np.array(v[1])) - p) ** 2)
    r = minimize(f, [t0, p0], method="Nelder-Mead",
                 options=dict(xatol=1e-12, fatol=1e-20, maxiter=400))
    dist.append(np.sqrt(r.fun))
    berr.append(modb(np.array(r.x[0]), np.array(r.x[1])))
dist = np.array(dist)
bv = np.array(berr)
rel_b = np.abs(bk - bv) / np.abs(bv)
info = json.load(open(a.kin6d_prefix + ".json"))
out = dict(points=len(dist), max_distance_m=float(dist.max()),
           rms_distance_m=float(np.sqrt(np.mean(dist ** 2))),
           max_rel_modB=float(rel_b.max()), rms_rel_modB=float(np.sqrt(np.mean(rel_b ** 2))),
           aminor_vmec=float(w["Aminor_p"][:]), kin6d=info)
print(json.dumps(out))
