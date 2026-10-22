#!/usr/bin/env python3
"""Fit a VMEC/VMEC++ wout to the axis-regular nested map read by nested_chart.f90.

Each stellarator-symmetric mode becomes rho^m * P(s) with a degree nk-1
polynomial P in s = rho^2 (least squares over the wout radial grid; lambda on
the half grid). The fitted map, not the wout, is the test geometry; the script
reports how far it is from the wout and checks the B-direction convention
(iota - lambda_zeta) / (1 + lambda_theta) = B^theta / B^zeta against the wout
Nyquist B^u, B^v on half-grid surfaces.
"""
import argparse
import sys

import netCDF4
import numpy as np


def fit(values, s, m, nk):
    rho = np.sqrt(s)
    basis = rho[:, None] ** m * s[:, None] ** np.arange(nk)[None, :]
    coef, *_ = np.linalg.lstsq(basis, values, rcond=None)
    return coef, np.max(np.abs(basis @ coef - values))


def series(coef_c, coef_s, xm, xn, s, theta, zeta, nk):
    """Sum of the fitted modes (cos part from coef_c, sin part from coef_s)."""
    rho = np.sqrt(s)
    out = np.zeros_like(theta)
    for j, (m, n) in enumerate(zip(xm, xn)):
        radial = rho**m * np.polyval(coef_c[j][::-1], s) if coef_c is not None else 0
        radial_s = rho**m * np.polyval(coef_s[j][::-1], s) if coef_s is not None else 0
        out += radial * np.cos(m * theta - n * zeta) + radial_s * np.sin(m * theta - n * zeta)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wout")
    parser.add_argument("output")
    parser.add_argument("--nk", type=int, default=6)
    args = parser.parse_args()
    data = netCDF4.Dataset(args.wout)
    if int(data.variables["lasym__logical__"][...]) != 0:
        sys.exit("only stellarator-symmetric wout files are supported")
    ns = int(data.variables["ns"][...])
    nfp = int(data.variables["nfp"][...])
    xm = np.asarray(data.variables["xm"][:], dtype=int)
    xn = np.asarray(data.variables["xn"][:], dtype=int)
    rmnc = np.asarray(data.variables["rmnc"][:])
    zmns = np.asarray(data.variables["zmns"][:])
    lmns = np.asarray(data.variables["lmns"][:])
    iotaf = np.asarray(data.variables["iotaf"][:])
    aminor = float(data.variables["Aminor_p"][...])
    s_full = np.linspace(0.0, 1.0, ns)
    s_half = (np.arange(1, ns) - 0.5) / (ns - 1)
    nk = args.nk
    rc, zs, ls = [], [], []
    err_r = err_z = err_l = 0.0
    for j, m in enumerate(xm):
        c, e = fit(rmnc[:, j], s_full, m, nk)
        rc.append(c)
        err_r = max(err_r, e)
        c, e = fit(zmns[:, j], s_full, m, nk)
        zs.append(c)
        err_z = max(err_z, e)
        c, e = fit(lmns[1:, j], s_half, m, nk)
        ls.append(c)
        err_l = max(err_l, e)
    iota, err_i = fit(iotaf, s_full, 0, nk)
    with open(args.output, "w") as out:
        out.write(f"# kin6d nested map v1 from {args.wout} nk={nk}\n")
        out.write(f"{nfp} {len(xm)} {nk}\n")
        out.write(" ".join(f"{v:.17e}" for v in iota) + "\n")
        for j in range(len(xm)):
            values = list(rc[j]) + list(zs[j]) + list(ls[j])
            out.write(f"{xm[j]} {xn[j]} " + " ".join(f"{v:.17e}" for v in values) + "\n")
    print(f"nfp={nfp} modes={len(xm)} nk={nk} ns={ns} aminor={aminor:.6g}")
    print(f"max coefficient residual / aminor: R {err_r / aminor:.3e} Z {err_z / aminor:.3e}")
    print(f"max lambda residual {err_l:.3e} iota residual {err_i:.3e}")
    # B-direction convention check on half-grid surfaces.
    xm_nyq = np.asarray(data.variables["xm_nyq"][:], dtype=int)
    xn_nyq = np.asarray(data.variables["xn_nyq"][:], dtype=int)
    bu = np.asarray(data.variables["bsupumnc"][:])
    bv = np.asarray(data.variables["bsupvmnc"][:])
    iotas = np.asarray(data.variables["iotas"][:])
    theta, zeta = np.meshgrid(np.linspace(0, 2 * np.pi, 37)[:-1],
                              np.linspace(0, 2 * np.pi / nfp, 19)[:-1])
    worst = 0.0
    for j in (ns // 4, ns // 2, ns - 2):
        lt = sum(m * lmns[j, k] * np.cos(m * theta - n * zeta) for k, (m, n) in enumerate(zip(xm, xn)))
        lz = sum(-n * lmns[j, k] * np.cos(m * theta - n * zeta) for k, (m, n) in enumerate(zip(xm, xn)))
        ratio = (iotas[j] - lz) / (1 + lt)
        bus = sum(bu[j, k] * np.cos(m * theta - n * zeta) for k, (m, n) in enumerate(zip(xm_nyq, xn_nyq)))
        bvs = sum(bv[j, k] * np.cos(m * theta - n * zeta) for k, (m, n) in enumerate(zip(xm_nyq, xn_nyq)))
        worst = max(worst, np.max(np.abs(ratio - bus / bvs)) / np.max(np.abs(bus / bvs)))
    print(f"B^theta/B^zeta convention check (wout lambda vs bsup): max rel diff {worst:.3e}")


if __name__ == "__main__":
    main()
