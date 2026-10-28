"""Multi-turn straight-field-line oracle for the iota sign in the chart frame.

KIN6D's frozen field (app/fb3d_synthetic.f90, iota_field_t) has
chi_gradient = sign * b3 * iota * y on the reference point y = (q1, q2) of the
curve chart (chart_frame.py). With b_pol = (-d chi/dq2, d chi/dq1) and b3 the
toroidal density (src/mhd/ideal_map.f90) a field line obeys
dq/dphi-hat = sign * iota * (-q2, q1): it turns counter-clockwise in (q1, q2)
exactly when sign * iota > 0.

The oracle traces the coil vacuum field in chart coordinates for many turns
and takes the rotation number about the vacuum axis (fixed point of the
period map) from the unwrapped poloidal angle. Its sign is the sign the frozen
field needs for a positive iota table. The rotation in the laboratory (R, Z)
plane at phi = 0 (one crossing per turn) is reported next to the VMEC++ iota
and the handedness of its theta (sign of zmns(m = 1, n = 0)). One JSON line.
"""
import argparse
import json
import pathlib
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import chart_frame as cf  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--turns", type=int, default=40)
ap.add_argument("--radii", type=float, nargs="*", default=[0.02, 0.05, 0.08])
ap.add_argument("--wout", default=None)
a = ap.parse_args()


def lab_rhs(p, y):
    r, z = y
    b = cf.bfield(np.array([r * np.cos(p), r * np.sin(p), z]))
    br = b[0] * np.cos(p) + b[1] * np.sin(p)
    bp = -b[0] * np.sin(p) + b[1] * np.cos(p)
    return [r * br / bp, r * b[2] / bp]


def lab_sections(y0, turns):
    ts = 2 * np.pi * np.arange(turns + 1)
    return solve_ivp(lab_rhs, (0, ts[-1]), y0, method="DOP853", rtol=1e-10, atol=1e-12,
                     t_eval=ts).y.T


axis = cf.chart_axis()
lab_axis = cf.to_lab(axis[0], axis[1], 0.0)
lab_axis = np.array([np.hypot(lab_axis[0], lab_axis[1]), lab_axis[2]])
res = dict(turns=a.turns, chart_axis_q=axis.tolist(), lab_axis_RZ=lab_axis.tolist(),
           chart=[], lab=[])
for r in a.radii:
    t = np.linspace(0, 2 * np.pi * a.turns, 64 * a.turns + 1)
    path = cf.trace(axis + [0.0, r], t[-1], t_eval=t).y.T
    ang = np.unwrap(np.arctan2(path[:, 1] - axis[1], path[:, 0] - axis[0]))
    res["chart"].append(dict(start_offset_q2=r,
                             iota_ccw_q1q2=float((ang[-1] - ang[0]) / t[-1]),
                             min_dist_axis=float(np.min(np.linalg.norm(path - axis, axis=1)))))
    p0 = cf.to_lab(axis[0], axis[1] + r, 0.0)
    sec = lab_sections([np.hypot(p0[0], p0[1]), p0[2]], a.turns)
    inc = np.diff(np.unwrap(np.arctan2(sec[:, 1] - lab_axis[1], sec[:, 0] - lab_axis[0])))
    res["lab"].append(dict(start_R=float(np.hypot(p0[0], p0[1])),
                           iota_ccw_RZ=float(np.mean(inc) / (2 * np.pi)),
                           spread=float(np.std(inc) / (2 * np.pi))))
signs = {float(np.sign(c["iota_ccw_q1q2"])) for c in res["chart"]}
res["required_sign_for_positive_iota_dat"] = signs.pop() if len(signs) == 1 else None
if a.wout:
    import netCDF4
    w = netCDF4.Dataset(a.wout)
    xm, xn = w["xm"][:], w["xn"][:]
    k = int(np.where((xm == 1) & (xn == 0))[0][0])
    zmns = float(w["zmns"][:][-1][k])
    res["vmec"] = dict(zmns_m1n0_edge=zmns, iota_axis=float(w["iotaf"][:][0]),
                       iota_edge=float(w["iotaf"][:][-1]), theta_ccw_in_RZ=zmns > 0)
print(json.dumps(res))
