"""Independent beta = 0 oracle: the vacuum flux surface of the coil field.

At zero pressure and zero plasma current both codes should return the coil
field's flux surface that encloses the toroidal flux phiedge. Field lines are
traced with direct Biot-Savart of the MAKEGRID polygons (make_coils.py) from
outboard points (R, Z = 0) of the phi = 0 symmetry plane. The surface whose
phi = 0 crossings enclose phiedge (bisection on R) is the reference. The phi = 0
sections of KIN6D (<prefix>.lcfs, periodic spline in phi-hat) and VMEC++
(wout) are compared with it by the distance of their points to the closed
crossing curve. Prints one JSON line.
"""
import argparse
import json
import sys
import pathlib

import netCDF4
import numpy as np
from matplotlib.path import Path
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import make_coils as mc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("phiedge", type=float)
ap.add_argument("--kin6d", nargs="*", default=[])
ap.add_argument("--wout", nargs="*", default=[])
ap.add_argument("--turns", type=int, default=60)
a = ap.parse_args()

polys, current = mc.coils()
A = np.concatenate(polys)
Bq = np.concatenate([np.roll(p, -1, 0) for p in polys])


def bfield(x):
    r1, r2 = x - A, x - Bq
    n1, n2 = np.linalg.norm(r1, axis=1), np.linalg.norm(r2, axis=1)
    f = (n1 + n2) / (n1 * n2 * (n1 * n2 + np.einsum("ij,ij->i", r1, r2)))
    return 1e-7 * current * (np.cross(r1, r2) * f[:, None]).sum(0)


def rhs(p, y):
    r, z = y
    b = bfield(np.array([r * np.cos(p), r * np.sin(p), z]))
    br = b[0] * np.cos(p) + b[1] * np.sin(p)
    bp = -b[0] * np.sin(p) + b[1] * np.cos(p)
    return [r * br / bp, r * b[2] / bp]


def crossings(r0):
    s = solve_ivp(rhs, (0, 2 * np.pi * a.turns), [r0, 0.0],
                  t_eval=2 * np.pi * np.arange(a.turns + 1), rtol=1e-11, atol=1e-13)
    c = s.y.T
    centre = c.mean(0)
    order = np.argsort(np.arctan2(c[:, 1] - centre[1], c[:, 0] - centre[0]))
    return c[order]


def flux(curve, h=0.002):
    path = Path(curve)
    rr, zz = np.meshgrid(np.arange(curve[:, 0].min(), curve[:, 0].max(), h) + h / 2,
                         np.arange(curve[:, 1].min(), curve[:, 1].max(), h) + h / 2)
    pts = np.c_[rr.ravel(), zz.ravel()]
    pts = pts[path.contains_points(pts)]
    return sum(bfield(np.array([r, 0.0, z]))[1] for r, z in pts) * h * h


def dist(points, curve):
    # distance to the closed polygon through the (dense) crossings
    seg_a, seg_b = curve, np.roll(curve, -1, 0)
    d = seg_b - seg_a
    out = []
    for p in points:
        t = np.clip(((p - seg_a) * d).sum(1) / (d * d).sum(1), 0, 1)
        out.append(np.min(np.linalg.norm(seg_a + t[:, None] * d - p, axis=1)))
    return np.array(out)


lo, hi = 1.16, 1.32
for _ in range(14):
    mid = 0.5 * (lo + hi)
    if flux(crossings(mid)) < a.phiedge:
        lo = mid
    else:
        hi = mid
ref = crossings(0.5 * (lo + hi))
res = dict(reference_outboard_R=0.5 * (lo + hi), reference_flux=flux(ref),
           reference_R=[float(ref[:, 0].min()), float(ref[:, 0].max())],
           reference_Z=[float(ref[:, 1].min()), float(ref[:, 1].max())])
for prefix in a.kin6d:
    k = np.loadtxt(prefix + ".lcfs")
    sec = []
    for node in range(int(np.sum(k[:, 0] == k[0, 0]))):
        rows = k[node::int(np.sum(k[:, 0] == k[0, 0]))]
        ph = np.r_[rows[:, 0], 2 * np.pi]
        xyz = np.r_[rows[:, 1:4], rows[:1, 1:4]]
        spl = CubicSpline(ph, xyz, bc_type="periodic")
        g = np.linspace(0, 2 * np.pi, 4001)
        phys = np.arctan2(spl(g)[:, 1], spl(g)[:, 0])
        j = np.where((phys[:-1] < 0) & (phys[1:] >= 0) & (np.abs(phys[:-1]) < 1))[0]
        for jj in j:
            t = -phys[jj] / (phys[jj + 1] - phys[jj])
            p = spl(g[jj] + t * (g[jj + 1] - g[jj]))
            sec.append([np.hypot(p[0], p[1]), p[2]])
    sec = np.array(sec)
    d = dist(sec, ref)
    res[prefix] = dict(points=len(sec), max_m=float(d.max()), rms_m=float(np.sqrt((d ** 2).mean())))
for wout in a.wout:
    w = netCDF4.Dataset(wout)
    th = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    xm = w["xm"][:]
    r = (w["rmnc"][:][-1] * np.cos(np.outer(th, xm))).sum(1)
    z = (w["zmns"][:][-1] * np.sin(np.outer(th, xm))).sum(1)
    d = dist(np.c_[r, z], ref)
    res[wout] = dict(points=len(th), max_m=float(d.max()), rms_m=float(np.sqrt((d ** 2).mean())))
print(json.dumps(res))
