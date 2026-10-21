#!/usr/bin/env python3
"""Retain public FreeGS 01-freeboundary.py data; never needed by CTest.

Run with the pinned environment in requirements.txt. All quantities are SI;
psi is Wb/rad, zero on the symmetry axis. No FreeGS implementation is copied.
"""

import argparse
from importlib.metadata import version
from pathlib import Path
import time
import platform
from functools import partial
from unittest.mock import patch

import freegs
import numpy as np
from freegs import critical
from freegs.gradshafranov import mu0
from scipy.integrate import romb
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import root

COMMIT = "8b838c1df162ca770a6937ac79d8c73e8f10a53b"
SOURCE = f"https://github.com/freegs-plasma/freegs/blob/{COMMIT}"
LEVELS = np.array([0.01, 0.5, 0.9, 0.95])


def solve(n, boundary):
    """Published physical case, selected vacuum algorithm, tighter Picard."""
    machine = freegs.machine.TestTokamak()
    eq = freegs.Equilibrium(
        tokamak=machine, Rmin=0.1, Rmax=2.0, Zmin=-1.0, Zmax=1.0,
        nx=n, ny=n, order=4, boundary=getattr(freegs.boundary, boundary),
    )
    profile = freegs.jtor.ConstrainPaxisIp(
        eq, 1e3, 2e5, 2.0, alpha_m=1.0, alpha_n=2.0, Raxis=1.0,
    )
    control = freegs.control.constrain(
        xpoints=[(1.1, -0.6), (1.1, 0.8)],
        isoflux=[(1.1, -0.6, 1.1, 0.6)], gamma=1e-12,
    )
    start = time.monotonic()
    freegs.solve(eq, profile, control, rtol=1e-9, atol=1e-12, maxits=500,
                 show=False)
    print(f"grid {n}: solve {time.monotonic() - start:.3f} s", flush=True)
    return eq, profile


def safety(eq, ntheta, nray):
    # FreeGS exposes the ray budget on find_psisurface, but find_safety does
    # not forward it. Bind that existing argument without changing its maths.
    ray = partial(critical.find_psisurface, n=nray)
    with patch.object(critical, "find_psisurface", ray):
        return critical.find_safety(eq, psinorm=LEVELS, ntheta=ntheta)


def observables(eq):
    op, xp = critical.find_critical(eq.R, eq.Z, eq.psi())
    # Match X-points by position, not their flux-dependent sort order.
    xp = sorted(xp, key=lambda x: x[1])
    return np.r_[
        eq.psi_axis, eq.psi_bndry, op[0][:2], xp[0][:2], xp[-1][:2],
        safety(eq, 1024, 800), eq.poloidalBeta(), eq.internalInductance2(),
        eq.plasmaCurrent(), eq.plasmaVolume(),
    ]


def critical_defects(eq):
    """Measure the fixed 1 mT stopping floor in FreeGS find_critical."""
    op, xp = critical.find_critical(eq.R, eq.Z, eq.psi())
    xp = sorted(xp, key=lambda x: x[1])
    spline = RectBivariateSpline(eq.R[:, 0], eq.Z[0, :], eq.psi())
    polished = []
    for point in [op[0], xp[0], xp[-1]]:
        gradient = lambda v: [spline.ev(*v, dx=1), spline.ev(*v, dy=1)]
        result = root(gradient, point[:2], tol=1e-11)
        if not result.success or np.linalg.norm(gradient(result.x)) > 1e-10:
            raise RuntimeError("critical-point diagnostic failed")
        polished.append(np.r_[result.x, spline.ev(*result.x)])
    return np.r_[polished[0][2] - op[0][2], polished[1][2] - xp[0][2],
                 polished[0][:2] - op[0][:2], polished[1][:2] - xp[0][:2],
                 polished[2][:2] - xp[-1][:2]]


def volume_readback(eq, n):
    """Refine masked volume quadrature on a fixed FreeGS field, no PDE solve.

    The bounding window contains the LCFS and both X-points. Use FreeGS's
    public field, pressure and core-mask routines, and its Romberg convention.
    """
    r = np.linspace(0.5, 1.8, n)
    z = np.linspace(-0.75, 0.9, n)
    R, Z = np.meshgrid(r, z, indexing="ij")
    psi = eq.psiRZ(R, Z)
    op, xp = critical.find_critical(eq.R, eq.Z, eq.psi())
    mask = critical.core_mask(R, Z, psi, op, xp, eq.psi_bndry)
    if (np.any(mask[[0, -1], :]) or np.any(mask[:, [0, -1]])):
        raise RuntimeError("volume readback window cuts the plasma")
    pn = (psi - eq.psi_axis) / (eq.psi_bndry - eq.psi_axis)
    pressure = np.zeros_like(pn)
    pressure[mask > 0] = eq.pressure(pn[mask > 0])
    dv = 2*np.pi*R*(r[1] - r[0])*(z[1] - z[0])*mask
    b2 = eq.Br(R, Z)**2 + eq.Bz(R, Z)**2
    p_integral = romb(romb(pressure*dv))
    b_integral = romb(romb(b2*dv))
    return [2*mu0*p_integral/b_integral,
            2*b_integral/(mu0**2*op[0][0]*eq.plasmaCurrent()**2),
            romb(romb(dv))]


def write_data(out, solutions, boundary="freeBoundary"):
    out.mkdir(parents=True, exist_ok=True)
    header = (
        f"FreeGS {version('freegs')}, tag v0.8.2, commit {COMMIT}\n"
        f"{SOURCE}/01-freeboundary.py; machine.py:TestTokamak\n"
        f"numpy {np.__version__}; scipy {version('scipy')}; "
        f"python {platform.python_version()} environment: requirements.txt\n"
        "ConstrainPaxisIp: paxis=1000 Pa, Ip=200000 A, fvac=2 T m, "
        "alpha_m=1, alpha_n=2, Raxis=1 m\n"
        "X targets (1.1,-0.6),(1.1,0.8); isoflux (1.1,-0.6)-(1.1,0.6)\n"
        f"{boundary}; order=4; domain R=[0.1,2], Z=[-1,1]; "
        "rtol=1e-9, atol=1e-12; control gamma=1e-12\n"
        "q: ntheta=1024, ray points=800; li2 uses magnetic-axis R\n"
        "candidate reference; continuous bounds unavailable; "
        "refinement differences are indicators\n"
    )
    rows = []
    for n, eq, profile in solutions:
        currents = [c.current for _, c in eq.tokamak.coils]
        rows.append(np.r_[n, currents, profile.pprime(0), profile.ffprime(0),
                          observables(eq)])
    np.savetxt(
        out / "reference.dat", rows, fmt="%.17e", header=header +
        "columns: n I_P1L I_P1U I_P2L I_P2U pp0 ff0 psi_axis psi_boundary "
        "axis_R axis_Z Xlower_R Xlower_Z Xupper_R Xupper_Z "
        "q(.01) q(.5) q(.9) q(.95) beta_p li2 Ip volume",
    )
    q_rows = []
    for n, eq, _ in solutions:
        for nt, nr in [(128, 100), (256, 200), (512, 400), (1024, 800)]:
            q_rows.append(np.r_[n, nt, nr, safety(eq, nt, nr)])
    np.savetxt(out / "q-convergence.dat", q_rows, fmt="%.17e",
               header=header + "columns: n_grid ntheta n_ray q(.01) q(.5) q(.9) q(.95)")
    critical_rows = [np.r_[n, critical_defects(e)] for n, e, _ in solutions]
    np.savetxt(out / "critical-refinement.dat", critical_rows, fmt="%.17e",
               header=header + "polished minus raw: n psi_axis psi_boundary "
               "axis_R axis_Z Xlower_R Xlower_Z Xupper_R Xupper_Z")
    _, eq, profile = solutions[-1]
    # Export the actual quadrature filaments, including their turn weights:
    # ShapedCoil uses polygon quadrature, not a MultiCoil rectangular grid.
    lines = [f"# {line}" for line in header.splitlines()]
    lines.append("kin6d-coil-table 1")
    for name, coil in eq.tokamak.coils:
        points = getattr(coil, "_points", [(coil.R, coil.Z, 1.0)])
        for r, z, weight in points:
            lines.append(f"coil {name} filament {r:.17e} {z:.17e} "
                         f"{weight * coil.turns:.17e}")
    for name, _ in eq.tokamak.coils:
        lines.append(f"circuit {name} -1e7 1e7 {name}")
    (out / "coils.dat").write_text("\n".join(lines) + "\n")
    # Small branch-selection seed on a regular grid, NOT boundary conditions.
    # KIN6D must remove its interpolation/offset error by solving its equations.
    r = np.linspace(0.2, 1.65, 25)
    z = np.linspace(-0.95, 0.95, 33)
    values = eq.psiRZ(r[:, None], z[None, :])
    np.savetxt(out / "seed.dat", values, fmt="%.17e", header=header +
               "25 R rows [0.2,1.65], 33 Z columns [-.95,.95]; branch seed only")
    ray = partial(critical.find_psisurface, n=800)
    with patch.object(critical, "find_psisurface", ray):
        lcfs = eq.separatrix(npoints=64)
    np.savetxt(out / "lcfs.dat", lcfs, fmt="%.17e",
               header=header + "LCFS: R Z (64 angles, FreeGS separatrix routine)")
    gaps = [[n, np.max(np.abs(e.psiRZ(lcfs[:, 0], lcfs[:, 1]) - e.psi_bndry))]
            for n, e, _ in solutions]
    np.savetxt(out / "lcfs-convergence.dat", gaps, fmt="%.17e",
               header=header + "columns: n_grid max_flux_gap_at_finest_LCFS (Wb/rad)")
    volume_rows = [np.r_[n, volume_readback(eq, n)] for n in [257, 513, 1025]]
    np.savetxt(out / "volume-convergence.dat", volume_rows, fmt="%.17e",
               header=header + "fixed finest field, window R=[.5,1.8], Z=[-.75,.9]; "
               "columns: n_quadrature beta_p li2 volume")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    parser.add_argument("--boundary", choices=["freeBoundary", "freeBoundaryHagenow"],
                        default="freeBoundary")
    args = parser.parse_args()
    if version("freegs") != "0.8.2":
        raise RuntimeError("reference requires FreeGS 0.8.2")
    write_data(args.output, [(n, *solve(n, args.boundary)) for n in [65, 129, 257]],
               args.boundary)
