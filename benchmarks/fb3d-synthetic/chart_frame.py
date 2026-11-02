"""Coil vacuum field and the helical curve chart of app/fb3d_synthetic.f90.

The chart is X = a(t) + q1 n(t) + q2 b(t) with a(t) the helical axis
R = 1 + eps cos 3t, Z = eps sin 3t, t the unit tangent, b = t x e_z/|t x e_z|
and n = b x t (projected-guide frame with guide e_z, see
src/equilibrium/equilibrium_chart.f90 and the generated curve_compose).
Field lines are traced in chart coordinates (q1, q2) against phi-hat = t.
"""
import pathlib
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import make_coils as mc  # noqa: E402

POLYS, CURRENT = mc.coils()
A = np.concatenate(POLYS)
BQ = np.concatenate([np.roll(p, -1, 0) for p in POLYS])
EZ = np.array([0.0, 0.0, 1.0])
PERIOD = 2 * np.pi / mc.NFP


def bfield(x):
    r1, r2 = x - A, x - BQ
    n1, n2 = np.linalg.norm(r1, axis=1), np.linalg.norm(r2, axis=1)
    f = (n1 + n2) / (n1 * n2 * (n1 * n2 + np.einsum("ij,ij->i", r1, r2)))
    return 1e-7 * CURRENT * (np.cross(r1, r2) * f[:, None]).sum(0)


def curve(t):
    r = mc.R0 + mc.EPS * np.cos(mc.NFP * t)
    return np.array([r * np.cos(t), r * np.sin(t), mc.EPS * np.sin(mc.NFP * t)])


def dcurve(t):
    r = mc.R0 + mc.EPS * np.cos(mc.NFP * t)
    dr = -mc.NFP * mc.EPS * np.sin(mc.NFP * t)
    return np.array([dr * np.cos(t) - r * np.sin(t), dr * np.sin(t) + r * np.cos(t),
                     mc.NFP * mc.EPS * np.cos(mc.NFP * t)])


def _unit(v):
    # analytic in t (no abs), so that a complex step differentiates it
    return v / np.sqrt(np.sum(v * v))


def frame(t):
    tan = _unit(dcurve(t))
    b = _unit(np.cross(tan, EZ))
    return curve(t), tan, np.cross(b, tan), b


def to_lab(q1, q2, t):
    x, _, n, b = frame(t)
    return x + q1 * n + q2 * b


def chart_rhs(t, y):
    # d X/d phi-hat by a complex step: exact to roundoff, so the ODE is smooth
    h = 1e-30
    x = to_lab(y[0], y[1], t)
    dt = np.imag(to_lab(y[0], y[1], t + 1j * h)) / h
    _, _, n, b = frame(t)
    c = np.linalg.solve(np.column_stack([n, b, dt]), bfield(x))
    return [c[0] / c[2], c[1] / c[2]]


def trace(y0, t_end, t_eval=None, t0=0.0):
    s = solve_ivp(chart_rhs, (t0, t_end), y0, method="DOP853", rtol=1e-10, atol=1e-12,
                  t_eval=t_eval, dense_output=t_eval is None)
    return s


def chart_axis(guess=(0.0, 0.0)):
    """Fixed point of the one-period map in chart coordinates (Newton, FD)."""
    y = np.array(guess, float)
    for _ in range(20):
        f = trace(y, PERIOD).sol(PERIOD) - y
        if np.linalg.norm(f) < 1e-11:
            break
        jac = np.empty((2, 2))
        for k in range(2):
            e = np.zeros(2)
            e[k] = 1e-6
            jac[:, k] = (trace(y + e, PERIOD).sol(PERIOD) - y - e - f) / 1e-6
        y = y - np.linalg.solve(jac, f)
    return y


def lab_point_to_chart(r, z, phi=0.0):
    """Chart coordinates (q1, q2, t) of the lab point (R, phi, Z) (Newton)."""
    p = np.array([r * np.cos(phi), r * np.sin(phi), z])
    q = np.array([0.0, 0.0, phi])
    for _ in range(50):
        f = to_lab(*q) - p
        if np.linalg.norm(f) < 1e-14:
            break
        h = 1e-7
        jac = np.column_stack([(to_lab(*(q + h * e)) - to_lab(*(q - h * e))) / (2 * h)
                               for e in np.eye(3)])
        q = q - np.linalg.solve(jac, f)
    return q
