# DK0: local drift-kinetic equation with finite shear and its limits

- Scope: the radially local, linearised drift-kinetic equation of
  [drift-kinetics.md](../../docs/drift-kinetics.md) on one closed flux
  surface, its reduction to the bounce-averaged toroidal precession of NEO-RT
  (shear on and off), and the quasilinear option (a) as a switch set on (b).
- Assume nested surfaces near psi_0, static (or slowly rotating) fields,
  subsonic rotation, Phi_0 = Phi_0(psi), thin orbits (rho* expansion to first
  order in the orbit width), gyro-averaged f.
- Units in the formulas are Gaussian-like with m, Ze, c explicit only where
  they matter; orientation signs are fixed by COCOS 3 in the implementation,
  not here. Status keys follow [theory-map.md](../../docs/theory-map.md).

## K1. Equation in (theta, alpha, lambda, x)

Field-line label from the D7 geometry, no straight-field-line angle needed:

    alpha = phi - nu(psi, theta),  nu = int_0^theta (B.grad phi)/(B.grad theta) dtheta',
    nu(psi, theta + 2 pi) = nu(psi, theta) + 2 pi q(psi).

Then b.grad alpha = 0, and with lambda = v_perp^2/(v^2 B), x = v/v_th
(invariants of the unperturbed motion at the order kept), the steady
linearised equation for f = sum_n f_n(theta, lambda, x, sigma) exp(i n alpha)
is

    sigma v |xi| b.grad theta d_theta f_n + i (n omega_alpha - omega) f_n - C(f_n) = S_n,
    omega_alpha = v_E.grad alpha + v_d.grad alpha,   |xi| = sqrt(1 - lambda B).      (K1)

Dropped: d_psi f at fixed (theta, alpha) (locality), v_d.grad theta d_theta
(O(rho*) against streaming), the radial drift except in the drive. The
twisted periodicity f_n(theta + 2 pi) = exp(2 pi i n q) f_n(theta) follows
from periodicity in (theta, phi). Status: stated.

## K2. Tangential drift with the shear term

    v_d.grad alpha = v_d.grad phi - (d_theta nu) v_d.grad theta - (d_psi nu) v_d.grad psi,
    d_psi nu = q'(psi) theta + (2 pi-periodic in theta).                                 (K2)

The last term is secular in theta. It is the finite-shear term: it is the
same order as the grad-B and curvature precession because theta = O(1).
Status: stated (algebra of a coordinate change).

## K3. Bounce average for trapped particles

Along the exact guiding-centre orbit, alpha-dot = v_d.grad alpha exactly
(b.grad alpha = 0 at the instantaneous position). For a trapped orbit theta
and psi are periodic in time, so nu(psi(t), theta(t)) is periodic and
<alpha-dot>_b = <phi-dot>_b: the bounce average of (K2) is the toroidal
precession. Evaluating (K2) on psi_0 is the first-order (thin-orbit) term.

Label invariance: alpha' = alpha + g(psi) changes v_d.grad alpha by
g' v_d.grad psi, whose trapped bounce average vanishes (<psi-dot>_b = 0).
The precession is therefore independent of the poloidal angle and of the
label origin; Boozer angles may be used to compare with oracles. Status:
stated.

## K4. Shear term by parts

In axisymmetry p_phi conservation gives, to the order kept,

    v_d.grad psi = v_par b.grad u,   u = I v_par / Omega,   I = R B_phi.          (K4)

For a closed trapped orbit, integrating the secular term by parts (u and
theta periodic, no boundary term):

    < -(d_psi nu) v_d.grad psi >_b = < u v_par b.grad theta  d_theta d_psi nu >_b
                                   = q' < u v_par b.grad theta >_b  (Boozer/SFL theta).  (K5)

d_theta d_psi nu = d_psi (B.grad phi / B.grad theta) is the radial derivative
of the local field-line pitch, defined without flux coordinates. Physical
reading: psi* = psi - u is conserved, the two legs sit at different radii
and see different q. The expansion parameter is n q' Delta_psi_b (banana
width); beyond it the problem is global. Status: hand derivation;
limit-checked numerically through K6.

## K5. Reduction to NEO-RT Om_tB

Boozer data (s, theta, phi), covariant B = I(s) grad theta + G(s) grad phi
(B_s = 0 as in NEO-RT, `hcovar(1) = 0`), Jacobian J = psi_tor'(G + iota I)/B^2,
alpha = phi - q theta. With v_d = (v_par/Omega) curl(v_par b) (the
perpendicular part equals the grad-B plus curvature drift; the parallel part
drops from v_d.grad alpha) and A = v_par b, A_s = 0:

    v_d.grad alpha = (v_par / (Omega J)) [ q d_s A_phi + d_s A_theta - q' theta d_theta A_phi ].   (K6)

NEO-RT `orbit.f90:timestep` (commit 8411537) integrates

    Om_tB / v^2 = m c q / (2 Z e psi_tor' B) [ -(2 - eta B) d_s B
                  + 2 (1 - eta B) h^theta (d_s B_theta + q d_s B_phi + B_phi q'_s) ],    (K7)

with `noshear` removing the q'_s term. The shear term of (K7) equals the
by-parts form (K5) of the secular term of (K6) identically (both reduce to
(q'_s q / psi_tor') times the bounce integral of v_par G/B dtheta). The
non-shear parts agree after bounce averaging.

Numerical check: `test/check_dk0_precession.py` builds (K6) with
SymPy derivatives of the covariant components (no metric, no NEO-RT
expression) and bounce-averages it in the secular form (not the by-parts
form) against (K7). Field: B = B00(s)(1 - eps(s) cos theta), eps = a sqrt(s)/R0,
R0 = 6.2, a = 2, s = 0.09, q = 1.6, I(s), G(s) and B00(s) with nonzero
s-derivatives. Gauss-Legendre in sin(theta/2) = k sin(chi), 320 nodes.

| case | kappa^2 range | max abs error / Om_tB(deep) | quadrature change 320 -> 640 |
|---|---|---:|---:|
| shear on, q'_s = 1.3 | 0.02-0.995 | 1e-11 | 2e-10 |
| shear off, q'_s = 0 (= NEO-RT `noshear`) | 0.02-0.995 | 2e-15 | 2e-10 |

The quadrature change is rounding near the turning points (v_par -> 0), not
truncation; the script gates both columns at 1e-9.

The shear term changes the precession by O(1) near the reversal point
(relative change up to 1.2; Om_tB crosses zero near kappa^2 = 0.8 here). The
deeply trapped value matches the stage-1 `omega_gradB_scale` limit
w eps'/chi' (plus the B00' shape term) to 2.7 % at kappa^2 = 0.02, as
expected from the O(kappa^2) curvature weight. Status: limit-checked
(trapped). The ansatz is adopted for trapped particles.

Not covered: passing particles. For a transit the by-parts step leaves the
boundary term -2 pi q' u(theta_0), which depends on the cut theta_0. It is
the gauge phase exp(i n q' theta u) of the twisted periodicity evaluated at
psi* instead of psi, and NEO-RT's passing frequency Om_th/iota + Om_tB must
be reproduced with the same cut convention before the passing branch is
adopted. Status: open (`needs review` below).

## K8. Drive of option (b) and the quasilinear switch set (a)

Linear response to the perturbed guiding-centre Hamiltonian (ideal gauge
delta A = xi x B0, pull-back along unperturbed motion):

    delta H_n = mu delta B_L,n + m v_par^2 kappa.xi_n + Z e delta Phi_L,n,
    S_n = (i n c/(Z e)) delta H_n d_psi F_M - i omega delta H_n d_E F_M,            (K8)

delta B_L = b.delta B + xi.grad B0. Status: stated from background knowledge
[K]; the sign and the omega term are fixed in DK1 against the stage-1
harness, not assumed.

Option (a) is (b) with four switches:

1. |B|-only drive: kappa.xi -> delta B_L/B. Then (K8) becomes
   (m v^2/2)(2 - lambda B) delta B_L/B, the NEO-RT H_n weight; this algebra
   is checked in the script. The replacement itself is a model change, not
   an identity: kappa.xi = xi.grad_perp ln B + mu0 xi.grad p/B^2, so it
   equals delta B_L/B only when the Eulerian b.delta B and the pressure term
   vanish. Its defect is measured in DK4 (MARS KNTVEUL evidence).
2. delta Phi = 0.
3. Static perturbation, omega = 0 (NEO-2-QL has no mode frequency).
4. Unperturbed surfaces and metric: delta B enters only through |B| on
   psi_0 in Boozer form with B_s dropped.

Common to (a) and (b) are the validity assumptions of the quasilinear
torque: amplitude small enough that unperturbed orbits are the basis, and
collisional broadening wider than nonlinear trapping. The torque is the
bilinear functional

    T_phi,n = -(n/2) Im < int d^3v delta H_n^* f_n >   (sign per COCOS 3),

[K], equivalent to the kinetic delta W_k form of
[Park 2011](https://doi.org/10.1063/1.3662039) (abstract-level evidence).
Status: stated; verified only through the DK1 harness.

## Open questions

- needs review: passing-particle shear term and cut convention (K5).
- needs review: K4 beyond axisymmetry (3D, option (c)) is not derived.
- needs review: FortSym leaf for (K6) with JVP/VJP; the SymPy script is an
  oracle, not the production source.
- needs review: sign and omega term of (K8) against the stage-1 harness.

Chris&AI
