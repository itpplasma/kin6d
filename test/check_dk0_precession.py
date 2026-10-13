"""DK0 gate: bounce-averaged toroidal precession of the local (theta, alpha) ansatz.

Independent side: v_d.grad(alpha) from the curl form v_d = (v_par/Omega) curl(v_par b)
(fixed energy and mu), alpha = phi - q(s) theta, with the secular shear term
-q'(s) theta v_d.grad(s) kept (drift-kinetics-dk0.md, Eqs. K3-K5).
Oracle: NEO-RT orbit.f90:timestep integrand Om_tB/v^2 (commit 8411537), with and
without its dqds shear term (`noshear`).

Field: axisymmetric circular large-aspect Boozer data on a surface neighbourhood,
B = B00(s) (1 - eps(s) cos theta), covariant B_theta = I(s), B_phi = G(s), B_s = 0,
psi_pr = dpsi_tor/ds, q(s); right-handed (s, theta, phi) so sign_theta = +1.
Units m = e = c = v = 1 (Omega = B).  The curl uses only covariant components and the
Boozer Jacobian, no metric and no NEO-RT expression.
"""
import json

import numpy as np
import sympy as sp

s, th, eta, sig = sp.symbols("s theta eta sigma", real=True)
R0, a, B0, s0 = 6.2, 2.0, 5.3, 0.09
eps = a * sp.sqrt(s) / R0
B00 = B0 * (1 + 0.15 * (s - s0))  # mean |B| rising with s (shape offset term)
Bmod = B00 * (1 - eps * sp.cos(th))
psi_pr = B0 * a**2 / 2


def model(qshear):
    q = 1.6 + qshear * (s - s0)
    G = B0 * R0 * (1 - 0.02 * (s - s0))  # B_phi covariant (poloidal current)
    I = B0 * a**2 * s / (q * R0)  # B_theta covariant (toroidal current)
    iota = 1 / q
    J = psi_pr * (G + iota * I) / Bmod**2  # Boozer Jacobian, sign_theta = +1
    vpar = sig * sp.sqrt(1 - eta * Bmod)
    A_th, A_ph = vpar * I / Bmod, vpar * G / Bmod  # covariant v_par b, A_s = 0
    dq = sp.diff(q, s)
    # contravariant curl: (s, theta, phi) = (d_th A_ph, -d_s A_ph, d_s A_th)/J
    # alpha = phi - q theta: covariant grad alpha = (-q' theta, -q, 1)
    curl_s = sp.diff(A_ph, th) / J
    curl_th = -sp.diff(A_ph, s) / J
    curl_ph = sp.diff(A_th, s) / J
    drift = (vpar / Bmod) * (curl_th * (-q) + curl_ph)
    shear = (vpar / Bmod) * curl_s * (-dq * th)
    # NEO-RT orbit.f90 timestep, Om_tB/v^2, hctrvr(3) = B^theta/B
    hth = iota * psi_pr / (J * Bmod)
    neort = lambda st: q / (2 * psi_pr * Bmod) * (
        -(2 - eta * Bmod) * sp.diff(Bmod, s)
        + 2 * (1 - eta * Bmod) * hth * (sp.diff(I, s) + q * sp.diff(G, s) + st * G * dq))
    dtdth = J * Bmod / (sp.Abs(vpar) * iota * psi_pr)  # dt/dtheta along a leg
    sub = {s: s0, sig: 1}
    f = lambda e: sp.lambdify((th, eta), e.subs(sub), "numpy")
    return dict(drift=f(drift), shear=f(shear), neort_on=f(neort(1)),
                neort_off=f(neort(0)), dt=f(dtdth))


def bounce_average(m, key, et, n=320):
    # sin(theta/2) = k sin(chi) removes the turning-point singularity uniformly in k
    e0, b00 = float(eps.subs(s, s0)), float(B00.subs(s, s0))
    k = np.sqrt((1 - et * b00 * (1 - e0)) / (2 * et * b00 * e0))
    x, w = np.polynomial.legendre.leggauss(n)
    chi = 0.5 * np.pi * x
    t = 2 * np.arcsin(k * np.sin(chi))
    jac = 0.5 * np.pi * w * 2 * k * np.cos(chi) / np.sqrt(1 - (k * np.sin(chi))**2)
    dt = m["dt"](t, et) * jac
    val = sum(np.asarray(m[name](t, et), dtype=float) for name in key)
    return np.sum(val * dt) / np.sum(dt)


Bmin = float(B00.subs(s, s0) * (1 - eps.subs(s, s0)))
Bmax = float(B00.subs(s, s0) * (1 + eps.subs(s, s0)))
rows, worst = [], 0.0
scale = 0.03  # |Om_tB/v^2| of deeply trapped particles on this surface
for label, qs in [("shear_on", 1.3), ("shear_off", 0.0)]:
    m = model(qs)
    for k2 in [0.02, 0.3, 0.6, 0.83, 0.95, 0.995]:
        et = 1.0 / (Bmin + k2 * (Bmax - Bmin))  # exact kappa^2 of this field
        ours = bounce_average(m, ["drift", "shear"], et)
        oracle = bounce_average(m, ["neort_on"], et)
        conv = abs(bounce_average(m, ["drift", "shear"], et, n=640) - ours)
        rel = abs(ours - oracle) / scale  # scaled: Om_tB crosses zero near kappa^2 0.8
        worst = max(worst, rel)
        shear_part = bounce_average(m, ["shear"], et)
        rows.append(dict(case=label, kappa2=k2, Om_tB_over_v2=ours, neort=oracle,
                         rel_err=rel, quad_change=conv, shear_part=shear_part))
    if qs == 0.0:  # noshear switch of NEO-RT equals q' = 0 of the ansatz
        m_on = model(1.3)
        et = 1.0 / (Bmin + 0.5 * (Bmax - Bmin))
        d = abs(bounce_average(m_on, ["drift"], et) - bounce_average(m_on, ["neort_off"], et))
        worst = max(worst, d / scale)

# shear term must matter: on/off differ by more than 1 % for some kappa^2
on = [r["Om_tB_over_v2"] for r in rows if r["case"] == "shear_on"]
off = [r["Om_tB_over_v2"] for r in rows if r["case"] == "shear_off"]
shear_effect = max(abs(x - y) / abs(y) for x, y in zip(on, off))

# deeply trapped limit (theta = 0, eta B = 1, curvature weight 1 - eta B -> 0):
# Om_tB/v^2 -> q (B00 eps' - B00' (1 - eps))/(2 psi_pr Bmin); with B00' = 0 and
# w = v^2/2 this is the stage-1 omega_gradB_scale w eps'/chi', chi' = psi_pr/q
eps_p = float(sp.diff(eps, s).subs(s, s0))
e0, b00, b00p = (float(x) for x in (eps.subs(s, s0), B00.subs(s, s0), sp.diff(B00, s).subs(s, s0)))
limit = 1.6 * (b00 * eps_p - b00p * (1 - e0)) / (2 * psi_pr * Bmin)
deep = [r for r in rows if r["case"] == "shear_off" and r["kappa2"] == 0.02][0]
deep_rel = abs(deep["Om_tB_over_v2"] / limit - 1)

# K8: (7) -> (5) under the |B|-only replacement kappa.xi -> dB_L/B (a model switch, not an
# identity): mu dB_L + m v_par^2 dB_L/B = (m v^2/2)(2 - lambda B) dB_L/B
mm, vv, lam, BB, dBL = sp.symbols("m v lambda B dB_L", positive=True)
k8 = sp.simplify(mm * vv**2 * lam / 2 * dBL + mm * vv**2 * (1 - lam * BB) * dBL / BB
                 - mm * vv**2 / 2 * (2 - lam * BB) * dBL / BB)
assert k8 == 0

out = dict(rows=rows, worst_rel_err=worst, shear_effect=shear_effect,
           deep_trapped_limit=limit, deep_trapped_rel=deep_rel)
print(json.dumps(out, indent=1))
assert worst < 1e-9, worst  # observed ~1e-11 (rounding near turning points)
assert max(r["quad_change"] for r in rows) / scale < 1e-9
assert shear_effect > 1e-2, shear_effect
assert deep_rel < 0.05, deep_rel  # O(kappa^2) curvature/shape corrections at kappa^2 = 0.02
