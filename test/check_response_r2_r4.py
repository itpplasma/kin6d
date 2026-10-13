"""T8 gate: linearised rotating MHD response (R2)-(R4) of response-r2-r4.md.

Exact identities are checked on concrete smooth fields in Cartesian coordinates
(x, y, z), with R = sqrt(x^2 + y^2), exp(i n phi) = ((x + i y)/R)^n and an
axisymmetric B0 = sigma grad psi x grad phi + F grad phi, V0 = Omega(psi) R^2 grad phi.
No cylindrical vector formula is used; every residual is evaluated at random points.
Layer checks: constant-psi slab layers (single fluid, drift-MHD) against closed forms
(FKR, Cole and Fitzpatrick 2006 Eqs. 41/47) and against SLAYER's layer ODE
(GPEC slayer/delta.f, w_der) in the limit pr, c_beta, ds, pe -> 0.
"""
import json

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp

I = sp.I
x, y, z = sp.symbols("x y z", real=True)
X = (x, y, z)
n = 2
om = sp.Rational(3, 10) + sp.Rational(1, 20) * I  # complex mode frequency
R = sp.sqrt(x**2 + y**2)
ph = ((x + I * y) / R) ** n
eR = sp.Matrix([x / R, y / R, 0])
ePh = sp.Matrix([-y / R, x / R, 0])
eZ = sp.Matrix([0, 0, 1])
grad_phi = ePh / R


def grad(f):
    return sp.Matrix([sp.diff(f, v) for v in X])


def div(a):
    return sum(sp.diff(a[i], X[i]) for i in range(3))


def curl(a):
    return sp.Matrix([sp.diff(a[2], y) - sp.diff(a[1], z), sp.diff(a[0], z) - sp.diff(a[2], x),
                      sp.diff(a[1], x) - sp.diff(a[0], y)])


def ddir(a, v):  # (a . grad) v
    return sp.Matrix([a.dot(grad(v[i])) for i in range(3)])


rng = np.random.default_rng(7)
PTS = [dict(zip(X, (float(rng.uniform(0.8, 1.6)) * c, float(rng.uniform(0.8, 1.6)) * s,
                    float(rng.uniform(-0.4, 0.4)))))
       for c, s in [(1, 0.3), (-0.4, 1), (0.7, -0.9)]]


def resid(vec, scale):
    """max |vec| over the points, relative to max |scale|."""
    v = sp.Matrix([vec]) if not isinstance(vec, sp.MatrixBase) else vec
    s = sp.Matrix([scale]) if not isinstance(scale, sp.MatrixBase) else scale
    num = max(abs(complex(sp.N(c.subs(p)))) for p in PTS for c in v)
    den = max(abs(complex(sp.N(c.subs(p)))) for p in PTS for c in s)
    return num / den


out = {}
for sigma in (1, -1):
    psi = R**2 * (1 + sp.Rational(3, 10) * z) + sp.Rational(1, 5) * z**2
    F = sp.Rational(17, 10) + sp.Rational(1, 5) * psi
    s_ = sp.Symbol("s")
    Om_of = lambda u: sp.Rational(2, 5) + sp.Rational(3, 10) * u + sp.Rational(1, 10) * u**2
    Om, Omp = Om_of(psi), sp.diff(Om_of(s_), s_).subs(s_, psi)
    B0 = sigma * grad(psi).cross(grad_phi) + F * grad_phi
    V0 = Om * R**2 * grad_phi
    w = om - n * Om
    poly = lambda a, b, c: (a + b * R + c * z + I * (b - c * R * z)) * ph
    xi = poly(0.3, 0.2, -0.1) * eR + poly(-0.2, 0.4, 0.3) * ePh + poly(0.5, -0.1, 0.2) * eZ
    A1 = poly(0.1, -0.3, 0.2) * eR + poly(0.4, 0.1, -0.2) * ePh + poly(-0.3, 0.2, 0.1) * eZ
    Phi1 = poly(0.2, 0.3, -0.4)
    dt = lambda v: -I * om * v
    xgp = xi.dot(grad(psi))
    r = {}
    # equilibrium: E0 = -grad Phi0, Phi0' = sigma Omega, satisfies E0 + V0 x B0 = 0
    r["equilibrium_E0"] = resid(-sigma * Om * grad(psi) + V0.cross(B0), V0.cross(B0))
    # C1: Eulerian and Lagrangian velocity of the displacement (Lie bracket form)
    v1 = dt(xi) + ddir(V0, xi) - ddir(xi, V0)
    r["v1_doppler_shear"] = resid(v1 - (-I * w * xi - Omp * xgp * R**2 * grad_phi), v1)
    dv = dt(xi) + ddir(V0, xi)
    r["Dv_coriolis"] = resid(dv - (-I * w * xi + Om * eZ.cross(xi)), dv)
    # C5: Frieman-Rotenberg inertia rho (d_t + V0.grad)^2 xi = rho (-i w + Omega z x)^2 xi
    lhs = dt(dv) + ddir(V0, dv)
    u = -I * w * xi + Om * eZ.cross(xi)
    r["FR_inertia"] = resid(lhs - (-I * w * u + Om * eZ.cross(u)), lhs)
    # C2: co-moving Ohm identity for arbitrary A1, Phi1, xi
    E1 = I * om * A1 - grad(Phi1)
    B1 = curl(A1)
    ohm = E1 + v1.cross(B0) + V0.cross(B1)
    G = A1 - xi.cross(B0)
    Gphi = (R**2 * grad_phi).dot(G)
    rhs = I * w * G - Omp * Gphi * grad(psi) - grad(Phi1 - V0.dot(A1))
    r["comoving_ohm"] = resid(ohm - rhs, ohm)
    # ideal limit, gauge A = xi x B0, Phi1 = -sigma Omega xi.grad psi
    Ai, Pi = xi.cross(B0), -sigma * Om * xgp
    ohm_i = I * om * Ai - grad(Pi) + v1.cross(B0) + V0.cross(curl(Ai))
    r["ideal_xi_x_B0"] = resid(ohm_i, v1.cross(B0))
    # ideal limit in the code gauge R A_phi = 0 (mhd_ampere_n.f90:12-13)
    Ac = Ai - I * sigma * grad(xgp) / n
    Pc = sigma * w / n * xgp
    ohm_c = I * om * Ac - grad(Pc) + v1.cross(B0) + V0.cross(curl(Ac))
    r["ideal_code_gauge"] = resid(ohm_c, v1.cross(B0))
    r["code_gauge_Aphi"] = resid((R**2 * grad_phi).dot(Ac), xgp)
    # C3: ideal radial field and the MEPHIT-form current closure identity
    r["Bpsi_mde"] = resid(grad(psi).dot(curl(xi.cross(B0))) - B0.dot(grad(xgp)), B0.dot(grad(xgp)))
    J0 = curl(sp.Matrix([z * R * 0.3, x * z + 0.2 * y, 0.5 * x * y]))  # div J0 = 0
    Bt = B1  # any divergence-free per-n field
    p1 = poly(0.3, -0.2, 0.4)
    b2 = B0.dot(B0)
    Jperp = B0.cross(grad(p1) - J0.cross(Bt)) / b2
    sig0 = J0.dot(B0) / b2
    code = -grad(p1).dot(curl(B0 / b2)) + J0.dot(grad(B0.dot(Bt) / b2)) - Bt.dot(grad(sig0))
    r["closure_divJperp"] = resid(-div(Jperp) - code, code)
    # C4: Lagrangian |B| and the gyrotropic pressure-tensor contraction
    Bm = sp.sqrt(b2)
    b = B0 / Bm
    kap = ddir(b, b)
    dBL = b.dot(curl(xi.cross(B0))) + xi.dot(grad(Bm))
    xpar = xi.dot(b)
    form = Bm * (-div(xi) + b.dot(grad(xpar)) - kap.dot(xi))
    r["deltaB_L"] = resid(dBL - form, dBL)
    pp, pq = sp.Rational(7, 10) + I / 3, sp.Rational(-2, 5) + I / 7  # p_perp1, p_par1
    gx = sp.Matrix(3, 3, lambda i, j: sp.diff(xi[j], X[i]))  # (grad xi)_{ij} = d_i xi_j
    P1 = pp * sp.eye(3) + (pq - pp) * b * b.T
    contr = sum(gx[i, j] * P1[i, j] for i in range(3) for j in range(3))
    r["gradxi_P1"] = resid(contr - (-pp * dBL / Bm + pq * (b.dot(grad(xpar)) - kap.dot(xi))), contr)
    out[f"sigma={sigma}"] = r
    assert max(r.values()) < 1e-12, (sigma, r)


# ---- slab layers ----------------------------------------------------------------
X_, d_, g, gE, gi, ge, eta, rho, kB, nu, Dp = sp.symbols(
    "X delta gamma gamma_E gamma_i gamma_e eta rho kB nu Delta_p", positive=True)
xx = sp.Symbol("x", positive=True)
Yf = sp.Function("Y")
layers = {}
# constant-psi drift-MHD layer (psi = 1): Ohm  ge (1 - x xi) = (eta/mu0) psi'',
# vorticity  mu0 rho gE gi xi'' - mu0 rho nu gi xi'''' = -kB^2 x psi''  (mu0 = 1)
for label, inertia, visc, dpow in [("RI", 1, 0, 4), ("VR", 0, 1, 6)]:
    xi_x = Yf(xx / d_) / d_  # xi = Y(X)/delta, x = delta X
    psi2 = ge / eta * (1 - xx * xi_x)
    eq = (inertia * rho * gE * gi * sp.diff(xi_x, xx, 2) - visc * rho * nu * gi * sp.diff(xi_x, xx, 4)
          + kB**2 * xx * psi2)
    dsol = sp.solve(sp.Eq(d_**dpow, (rho * gE * gi * inertia + rho * nu * gi * visc) * eta
                          / (kB**2 * ge)), d_)[0]
    eqX = sp.simplify((eq.subs(d_, dsol).subs(xx, X_ * dsol)) / (kB**2 * ge * dsol / eta))
    target = (sp.diff(Yf(X_), X_, 2) if inertia else -sp.diff(Yf(X_), X_, 4)) + X_ * (1 - X_ * Yf(X_))
    layers[label] = dict(delta=str(dsol), parameter_free=bool(sp.simplify(eqX - target) == 0))
    # Delta = int psi'' dx = (ge/eta) delta int (1 - X Y) dX
    Dl = sp.powsimp(ge / eta * dsol, force=True)
    layers[label]["Delta_over_hat"] = str(Dl)
    assert layers[label]["parameter_free"], label
# exponents: RI  Delta ~ ge^(3/4) (gE gi)^(1/4) eta^(-3/4); VR  ge^(5/6) gi^(1/6) nu^(1/6) eta^(-5/6)
DRI = sp.powsimp(ge / eta * (rho * gE * gi * eta / (kB**2 * ge)) ** sp.Rational(1, 4), force=True)
assert sp.simplify(sp.diff(sp.log(DRI), ge) * ge - sp.Rational(3, 4)) == 0
assert sp.simplify(sp.diff(sp.log(DRI), gE) * gE - sp.Rational(1, 4)) == 0
assert sp.simplify(sp.diff(sp.log(DRI), gi) * gi - sp.Rational(1, 4)) == 0
assert sp.simplify(sp.diff(sp.log(DRI), eta) * eta + sp.Rational(3, 4)) == 0
DVR = sp.powsimp(ge / eta * (rho * nu * gi * eta / (kB**2 * ge)) ** sp.Rational(1, 6), force=True)
assert sp.simplify(sp.diff(sp.log(DVR), ge) * ge - sp.Rational(5, 6)) == 0
assert sp.simplify(sp.diff(sp.log(DVR), gi) * gi - sp.Rational(1, 6)) == 0
assert sp.simplify(sp.diff(sp.log(DVR), nu) * nu - sp.Rational(1, 6)) == 0

# single fluid: ge = gE = gi = gamma; FKR: Delta_layer(gamma) = Delta' with eta ~ 1/S
DRI1 = DRI.subs({ge: g, gE: g, gi: g})
assert sp.simplify(sp.diff(sp.log(DRI1), g) * g - sp.Rational(5, 4)) == 0
S = sp.Symbol("S", positive=True)
gam = sp.solve(sp.Eq(DRI1.subs(eta, 1 / S), Dp), g)[0]
dl = (rho * g / (S * kB**2)) ** sp.Rational(1, 4)
fkr = dict(gamma_S=sp.simplify(sp.diff(sp.log(gam), S) * S),
           gamma_Dp=sp.simplify(sp.diff(sp.log(gam), Dp) * Dp),
           width_S=sp.simplify(sp.diff(sp.log(dl.subs(g, gam)), S) * S))
assert fkr == dict(gamma_S=sp.Rational(-3, 5), gamma_Dp=sp.Rational(4, 5), width_S=sp.Rational(-2, 5))
# constant-psi validity delta * Delta ~ 1 at gamma ~ S^(-1/3), width ~ S^(-1/3) there
cross = sp.solve(sp.Eq(dl * DRI1.subs(eta, 1 / S), 1), g)[0]
cpsi = dict(gamma_S=sp.simplify(sp.diff(sp.log(cross), S) * S),
            width_S=sp.simplify(sp.diff(sp.log(dl.subs(g, cross)), S) * S))
assert cpsi == dict(gamma_S=sp.Rational(-1, 3), width_S=sp.Rational(-1, 3))
# visco-resistive width delta^6 ~ rho nu eta/kB^2: Pm^(1/6) S^(-1/3) with nu = Pm eta
Pm = sp.Symbol("P_m", positive=True)
dvr = (rho * Pm / S * (1 / S) / kB**2) ** sp.Rational(1, 6)
assert sp.simplify(sp.diff(sp.log(dvr), S) * S + sp.Rational(1, 3)) == 0


# numerical layer constants: BVP on [0, L], odd Y, Y ~ 1/X (ideal outer solution)
def delta_hat(order, L=12.0):
    if order == 2:  # Y'' = X^2 Y - X
        f = lambda t, u: np.vstack([u[1], t**2 * u[0] - t])
        bc = lambda a, b: np.array([a[0], b[0] - (1 / L + 2 / L**5)])
        u0 = np.vstack([np.tanh(np.linspace(0, L, 400)), np.zeros(400)])
    else:  # Y'''' = X - X^2 Y
        f = lambda t, u: np.vstack([u[1], u[2], u[3], t - t**2 * u[0]])
        bc = lambda a, b: np.array([a[0], a[2], b[0] - (1 / L - 24 / L**7), b[1] - (-1 / L**2 + 168 / L**8)])
        u0 = np.zeros((4, 400))
        u0[0] = np.tanh(np.linspace(0, L, 400))
    t = np.linspace(0, L, 400)
    sol = solve_bvp(f, bc, t, u0, tol=1e-10, max_nodes=200000)
    assert sol.success
    tt = np.linspace(0, L, 200001)
    val = 2 * np.trapezoid(1 - tt * sol.sol(tt)[0], tt)
    # int_L^inf of 1 - XY: Y = 1/X + 2/X^5 (RI), Y = 1/X - 24/X^7 (VR)
    tail = 2 * (-2 / (3 * L**3) if order == 2 else 24 / (5 * L**5))
    return val + tail


ri_exact = float(2 * mp.pi * mp.gamma(0.75) / mp.gamma(0.25))
vr_exact = float(6 ** (2 / 3) * mp.pi * mp.gamma(5 / 6) / mp.gamma(1 / 6))
ri_num, vr_num = delta_hat(2), delta_hat(4)
# integral representation Y = (X/2) int_0^1 (1 - t^2)^(-1/4) exp(-X^2 t/2) dt
Yrep = lambda Xv: Xv / 2 * mp.quad(lambda t: (1 - t**2) ** (-0.25) * mp.exp(-Xv**2 * t / 2), [0, 1])
rep_res = max(abs(mp.diff(Yrep, Xv, 2) - Xv**2 * Yrep(Xv) + Xv) for Xv in (0.3, 1.1, 2.5))
layers["constants"] = dict(RI_num=ri_num, RI_exact=ri_exact, VR_num=vr_num, VR_exact=vr_exact,
                           Yrep_residual=float(rep_res))
assert abs(ri_num / ri_exact - 1) < 1e-6 and abs(vr_num / vr_exact - 1) < 1e-6
assert rep_res < 1e-10

# SLAYER (GPEC slayer/delta.f w_der, PeOhmOnly): (C3 Y')' = G Y in Fourier variable p.
p, Q, Qe, Qi, pr, cb, ds, tau = sp.symbols("p Q Q_e Q_i pr c_beta ds tau")
G_sl = ((cb**2 * pr * p**4 - Q * (Q - Qi) + I * (cb**2 + pr) * p**2 * (Q - Qi))
        / (ds**2 * pr * (1 + tau) * p**4 + I * (Q - Qe) + p**2 * (cb**2 + I * ds**2 * (Q - Qi)))) * p**2
C3_sl = p**2 / (p**2 + I * (Q - Qe))
G0 = sp.simplify(G_sl.subs({pr: 0, cb: 0, ds: 0}))
# ours, Fourier transform x -> i d/dp, d/dx -> i p of the drift-MHD layer (all scales 1):
#   ge (psi - x xi) = psi'',  gE gi xi'' = -x psi''
xh, ph_ = sp.Function("xh")(p), sp.Function("ph")(p)
ohm_f = sp.Eq(ge * (ph_ - I * sp.diff(xh, p)), -p**2 * ph_)
psol = sp.solve(ohm_f, ph_)[0]
vort = gE * gi * (-p**2 * xh) + I * sp.diff(-p**2 * psol, p)  # = 0
ours = sp.simplify(vort / (-1))  # (p^2/(ge + p^2) xh')' - (gE gi/ge) p^2 xh (times ge)
ref = sp.diff(p**2 / (ge + p**2) * sp.diff(xh, p), p) - gE * gi / ge * p**2 * xh
ratio = sp.simplify(ours / ref)
conv = {ge: I * (Q - Qe), gE: I * Q, gi: I * (Q - Qi)}
slayer = dict(our_ode_over_ref=str(ratio),
              C3_match=bool(sp.simplify(C3_sl - (p**2 / (ge + p**2)).subs(conv)) == 0),
              G_match=bool(sp.simplify(G0 - (gE * gi / ge * p**2).subs(conv)) == 0))
assert ratio.free_symbols <= {ge} and slayer["C3_match"] and slayer["G_match"]
# penetration: Delta_RI -> 0 as ge -> 0 (electron fluid at rest in the frame of the drive)
assert sp.limit(DRI, ge, 0) == 0

print(json.dumps(dict(identities=out, layers=layers, fkr={k: str(v) for k, v in fkr.items()},
                      constant_psi_limit={k: str(v) for k, v in cpsi.items()}, slayer=slayer),
                 indent=1, default=str))
