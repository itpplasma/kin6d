"""VMEC++ free-boundary reference for the synthetic nfp = 3 case.

beta = 0, zero net toroidal current (ncurr = 1, curtor = 0), phiedge from
--phiedge. The external field is VMEC++'s makegrid response table computed
from the same MAKEGRID coils file that KIN6D reads (one circuit, EXTCUR = 1).
Writes <out>/wout.nc, <out>/lcfs.dat (R, Z on phi cuts), <out>/iota.dat and
<out>/timing.json (wall time of the makegrid step and of the run).
"""
import argparse
import json
import pathlib
import time

import numpy as np
import vmecpp

ap = argparse.ArgumentParser()
ap.add_argument("coils")
ap.add_argument("out")
ap.add_argument("--phiedge", type=float, default=0.0314159)
ap.add_argument("--mpol", type=int, default=6)
ap.add_argument("--ntor", type=int, default=4)
ap.add_argument("--ns", type=int, nargs="+", default=[16, 32])
ap.add_argument("--ftol", type=float, default=1e-12)
ap.add_argument("--nzeta", type=int, default=36)
ap.add_argument("--grid", type=int, nargs=3, default=[121, 121, 36])
a = ap.parse_args()
out = pathlib.Path(a.out)
out.mkdir(parents=True, exist_ok=True)
R0, EPS, A = 1.0, 0.15, 0.1

t0 = time.perf_counter()
grid = vmecpp.MakegridParameters(
    normalize_by_currents=False, assume_stellarator_symmetry=True,
    number_of_field_periods=3, r_grid_minimum=0.70, r_grid_maximum=1.30,
    number_of_r_grid_points=a.grid[0], z_grid_minimum=-0.30, z_grid_maximum=0.30,
    number_of_z_grid_points=a.grid[1], number_of_phi_grid_points=a.grid[2])
table = vmecpp.MagneticFieldResponseTable.from_coils_file(a.coils, grid)
t_grid = time.perf_counter() - t0

ns = a.ns
d = dict(
    lasym=False, nfp=3, mpol=a.mpol, ntor=a.ntor, ntheta=0, nzeta=a.nzeta,
    ns_array=ns, ftol_array=[a.ftol] * len(ns), niter_array=[20000] * len(ns),
    delt=0.7, tcon0=1.0, aphi=[1.0], phiedge=a.phiedge, nstep=200,
    pmass_type="power_series", am=[0.0], pres_scale=0.0, gamma=0.0, spres_ped=1.0,
    ncurr=1, pcurr_type="power_series", ac=[0.0], curtor=0.0, bloat=1.0,
    lfreeb=True, mgrid_file="", extcur=[1.0], nvacskip=6, lforbal=False,
    raxis_c=[R0, EPS] + [0.0] * (a.ntor - 1), zaxis_s=[0.0, -EPS] + [0.0] * (a.ntor - 1),
    rbc=[{"n": 0, "m": 0, "value": R0}, {"n": 1, "m": 0, "value": EPS},
         {"n": 0, "m": 1, "value": A}],
    zbs=[{"n": 1, "m": 0, "value": -EPS}, {"n": 0, "m": 1, "value": A}],
)
path = out / "input.json"
path.write_text(json.dumps(d, indent=1))
vin = vmecpp.VmecInput.from_file(path)
t0 = time.perf_counter()
res = vmecpp.run(vin, magnetic_field=table, verbose=False)
t_run = time.perf_counter() - t0
w = res.wout
w.save(str(out / "wout.nc"))
# LCFS on phi cuts (per period phi = 0, pi/6, pi/3 of the full turn) and |B|.
theta = np.linspace(0, 2 * np.pi, 257)[:-1]
rows = []
for phi in [0.0, np.pi / 9, np.pi / 6, 2 * np.pi / 9]:
    ang = np.outer(theta, w.xm) - phi * w.xn[None, :]
    r = (w.rmnc[:, -1][None, :] * np.cos(ang)).sum(1)
    z = (w.zmns[:, -1][None, :] * np.sin(ang)).sum(1)
    angb = np.outer(theta, w.xm_nyq) - phi * w.xn_nyq[None, :]
    b = (w.bmnc[:, -1][None, :] * np.cos(angb)).sum(1)
    for t_, r_, z_, b_ in zip(theta, r, z, b):
        rows.append((phi, t_, r_, z_, b_))
np.savetxt(out / "lcfs.dat", np.array(rows), header="phi theta R Z |B|(half-grid edge)")
s = np.linspace(0, 1, w.ns)
np.savetxt(out / "iota.dat", np.c_[s, w.iotaf], header="s iota_full")
info = dict(makegrid_s=t_grid, run_s=t_run, fsqr=float(w.fsqr), fsqz=float(w.fsqz),
            fsql=float(w.fsql), ns=int(w.ns), mpol=a.mpol, ntor=a.ntor,
            iterations=int(w.itfsq), volume=float(w.volume_p),
            iota_axis=float(w.iotaf[0]), iota_edge=float(w.iotaf[-1]),
            aminor=float(w.Aminor_p), rmajor=float(w.Rmajor_p), phiedge=a.phiedge)
(out / "timing.json").write_text(json.dumps(info, indent=1))
print(json.dumps(info))
