# Circular-torus flux-pumping case and its toroidal limits

This case is the 3D part of RM-B8 and a physics case of RM7, after the straight
cylinder, in [the MHD design](nonlinear-resistive-mhd.md). It asks whether
plasma pressure helps a saturated m = 1, n = 1 mode keep the central safety
factor `q0` near one (flux pumping) at Lundquist numbers `S = 1e6` to `1e7`.
The helically symmetric cylinder of [vrmhd-cylinder.md](vrmhd-cylinder.md)
cannot answer this, because the pressure effects that matter near `q0 = 1` are
toroidal: Bussac stabilisation of the ideal internal kink and the
quasi-interchange come from the coupling of m = 1 to m = 2 through the
toroidal metric and the Shafranov shift. In single-helicity reduced MHD of the
cylinder case, flux pumping ends near `S = 1.3e6` at magnetic Prandtl number
10. JOREK runs of this torus at `R0/a = 10`, branched from pumped states at
`S = 1e6` to `S = 3.2e6` and `1e7`, show repeated relaxations of `q0` (to 0.95-0.98)
instead of flux pumping at `beta_p1 = 0.005` and 0.23, up to `1e5 tau_A` after the
branch; they are not finished. The setup below is that of these JOREK runs,
so that their numbers can be compared.

The page also collects the analytic limits that a toroidal evolution has to
meet before its pressure scans are trusted. The Bussac, Shafranov,
geodesic-acoustic and Suydam reference values are computed by the test-side
references of [linear-mhd-oracles.md](linear-mhd-oracles.md)
(`test/support/rm_mhd_oracles.f90`, leaves in `gen/gen_rm_oracles.f90`) and
checked against an independent evaluation in the tests `rm_oracle_bussac`,
`rm_oracle_shafranov` and `rm_oracle_torus_closed`. The Grad-Shafranov column
of the Shafranov table comes from that independent evaluation only. The
proposed checks on an evolution are not registered yet.

## Geometry, units and symmetry

Circular fixed boundary of minor radius `a = 1` about the major radius `R0`,
with `R0/a = 3`, 5 and 10. The vacuum toroidal field is `B0 = 1` at `R0`,
`mu0 = mi = 1`, and the mass density is uniform, so time is in Alfven times
`tau_A = a sqrt(mu0 rho)/B0`. The symmetry is the rotation about the major
axis of `gs_symmetry_t`, and the 3D state is the meridional section times the
toroidal harmonics `n = 0, 1` (optionally 2) of
[the section contract](vrmhd-section.md#fourier-harmonics-about-an-invariant-state).
The wall is an ideal no-slip conductor with fixed temperature. The loop
voltage drives the current through the rotation harmonic form
`h = grad(phi)/(2 pi)` of [the driven-state conventions](vrmhd-section.md#driven-state-conventions).

## Profiles

The profile recipe is that of the cylinder case (`vrmhd_cylinder_profiles`):
uniform perpendicular conductivity, heat source
`Q0 (1 - tanh((psi_N - 0.3)/0.1))/2` with `psi_N` the normalised poloidal flux and
`Q0` fixed by the recipe, Ohmic current `j ~ T^(3/2)`, edge
temperature `0.05` of the central-to-edge difference and `q0 = 0.9`. As a
function of the minor radius of that cylinder state, `q(a) = 3.352`, `q = 1` at
`r1 = 0.264` with shear `s1 = r1 q'(r1) = 0.216`, and `q = 2` at `r2 = 0.729`.
The pressure is `p = T` with `Te = Ti = T/2`, scaled to the global poloidal beta
`beta_p = 2 <p>/B_theta(a)^2` (`<.>` the area average). The poloidal beta inside
`q = 1`,

```text
beta_p1 = (2/B_theta(r1)^2) int_0^r1 (r/r1)^2 (-dp/dr) dr,
```

is `0.234 beta_p` for this profile. The central beta `beta_0 = 2 p(0)/B0^2` is

| beta_p | beta_p1 | beta_0, R0/a = 3 | 5 | 10 |
|---|---|---|---|---|
| 0.02 | 0.005 | 0.05 % | 0.02 % | 0.005 % |
| 0.3 | 0.070 | 0.80 % | 0.29 % | 0.07 % |
| 1 | 0.234 | 2.66 % | 0.96 % | 0.24 % |
| 3 | 0.702 | 7.97 % | 2.87 % | 0.72 % |

In the torus the profiles are flux functions of the axisymmetric state:
`p(psi_N)` and `FF'(psi_N)` of the cylinder state, with `F = R B_phi`, define a
fixed-boundary Grad-Shafranov problem whose solution is the initial `n = 0`
state. Its `q0` is reported, not rescaled to 0.9.

## Closures and sources

As in [the cylinder case](vrmhd-cylinder.md#closures-and-sources), with

| quantity | value |
|---|---|
| resistivity | isotropic, `eta0 (Te/T0)^(-3/2)`, `T0 = Te(0)`, `eta0 = 1/S`, `S = 1e6`, `3e6`, `1e7` |
| viscosity | constant, `Pm/S`, `Pm = 10` |
| perpendicular conduction | `chi = 2.4e-5` in `d_t T = chi lap T` |
| parallel conduction | `chi_par = 24` (`chi_par/chi_perp = 1e6`), also 240 and 2400 (`1e7`, `1e8`) |

The heat source and the loop voltage are built on the flux surfaces of the
Grad-Shafranov state. A source built on the cylinder radius instead leaves a
temperature gradient along the shifted surfaces that `chi_par/chi_perp = 1e6`
removes in a transient, which moves `eta`.

## Shafranov shift

The flux surfaces of the large-aspect-ratio equilibrium are circles of radius
`r` about `R0 + Delta(r)`, `Delta(a) = 0` (Wesson, *Tokamaks*, Secs. 3.6 and
3.7, Eq. 3.6.7). One integration from the axis gives

```text
dDelta/dr = -(r/R0) [beta_p(r) + l_i(r)/2],
beta_p(r) = 2 [<p>_r - p(r)]/B_theta(r)^2,
l_i(r)    = 2 int_0^r B_theta^2 r' dr' / (r^2 B_theta(r)^2),
```

with `<p>_r` the area average inside `r`. The axis shift `Delta(0)` of this
formula was compared with the magnetic axis of a fixed-boundary
Grad-Shafranov solve with the same `p(psi_N)` and `FF'(psi_N)`:

| R0/a | beta_p | Delta(0)/a, formula | Grad-Shafranov | ratio | max abs dDelta/dr |
|---|---|---|---|---|---|
| 3 | 0.02 | 0.0744 | 0.0755 | 1.016 | 0.22 |
| 3 | 0.3 | 0.0984 | 0.0993 | 1.009 | 0.30 |
| 3 | 1 | 0.159 | 0.157 | 0.988 | 0.51 |
| 3 | 3 | (0.331) | 0.289 | | 1.09 |
| 5 | 0.02 | 0.0446 | 0.0449 | 1.006 | 0.13 |
| 5 | 1 | 0.0952 | 0.0948 | 0.996 | 0.30 |
| 5 | 3 | 0.198 | 0.188 | 0.949 | 0.65 |
| 10 | 0.02 | 0.0223 | 0.0223 | 1.002 | 0.07 |
| 10 | 1 | 0.0476 | 0.0476 | 0.999 | 0.15 |
| 10 | 3 | 0.0992 | 0.0979 | 0.987 | 0.33 |

At `R0/a = 3`, `beta_p = 3` the shifted circles of the formula intersect
(`|dDelta/dr| > 1`), so only the Grad-Shafranov value applies. Even at
`beta_p = 0.02` the shift is not small, because `l_i(a)/2` drives it.
Proposed check: the `n = 0` state has its axis at the formula value within
3 % where `max |dDelta/dr| <= 0.3`, and at the Grad-Shafranov value within
1 % elsewhere.

## Steady Ohmic state

The flux-surface average of the parallel Ohm's law in steady state
(Wesson Eqs. 4.4.6-4.4.9) is `eta <J.B> = V_psi F <R^-2>`, with `V_psi` the
uniform rate of the poloidal flux and `<.>` here the flux-surface average. Its
cylinder limit `eta J_z = E_z` is the condition of the cylinder case.

- At `p = 0` a static state exists exactly, with `FF'(psi) = V_psi/eta(psi)`
  and `R J_phi = V_psi/eta` uniform on each surface. A prescribed current
  source `S_j e_phi` keeps this only if `R S_j` is a flux function. With the
  heat source built on the surfaces, this state has to stay static to roundoff.
- At `p != 0` no static Ohmic state exists: `R^2 p' + FF'` cannot equal the
  `R`-independent `V_psi/eta(psi)` on a surface. Pfirsch-Schlueter currents
  and a flow appear. With isotropic `eta`, the surface-averaged radial flow is
  `<v R>/R0 = -(1 + 2 q^2) eta p'/B0^2 - E_phi B_theta/B0^2` (Wesson Eq. 4.5.5,
  `E_phi` the toroidal electric field), of order `1e-7 a/tau_A` at
  `beta_p = 1`, `R0/a = 3`, `S = 1e6`. Proposed check: over `1e4 tau_A` the
  axisymmetric state drifts by less than `1e-3` in `q0` and `1e-3 a` in
  `Delta(0)`, and its flow is at most ten times this estimate.

## Bussac limit of the ideal internal kink

Bussac, Pellat, Edery and Soule, Phys. Rev. Lett. 35, 1638 (1975), Eqs. (3)
and (4), give the energy of a rigid shift `xi0` inside `r1` in a circular
large-aspect-ratio torus as
`delta W = 2 pi^2 R0 B0^2 eps1^4 |xi0|^2 delta W_T` with `eps1 = r1/R0` and

```text
delta W_T = [8 s_B(1+b-c) + 9b(1-c) - 24bc(beta_p1+s_B) - 16c(1+b)(beta_p1+s_B)^2] / [16(1+b-c)],
s_B = int_0^r1 (dr/r1) (r/r1)^3 (q^-2 - 1),
b = (r/4) d/dr ln(u1 r^-2) at r1,   c = (r/4) d/dr ln(u2 r^2) at r1,
(1/r)(r u')' - 4u/r^2 + 2 [d/dr ln(1/q - 1/2)] (u' - u/r) = 0,
```

where `u1` is regular at the axis and `u2` regular at `q = 2` (`u' = u/r`
there) or zero at the wall when `q(a) < 2`. To this order the cylindrical
contribution vanishes for `n = 1`, so in a straight cylinder the ideal m = 1
mode is marginal at order `eps1^2` and grows resistively. The numerical
evaluation reproduces Bussac's small-`Delta q` formula
`3 nu Delta q [13/(48(nu+4)) - beta_p1^2 int_r1^r2 (dr/r1)(r/r1)^(nu-5)]`
for `q = 1 - Delta q (1 - (r/r1)^nu)` to `1e-4` at `nu = 2`,
`Delta q = 1e-4`, `r1 = 1e-3`.

Above the critical `beta_p1` the mode is ideally unstable. Minimising the
Lagrangian of the layer at `r1`, where `k_par = -q' x/R0`, against `delta W`
gives

```text
gamma tau_A = -2 pi eps1^2 delta W_T / (R0 s1 sqrt(f_in)).
```

The inertia factor `f_in = 1 + 2 q^2 = 3` at `q = 1` (Bussac, p. 1641, writes
`1 + 2/n^2`) comes from the toroidal displacement `xi_phi = 2 i xi_r/n` that
makes the shift divergence-free. Because the shift is divergence-free, the
compressional energy `Gamma p (div xi)^2` (`Gamma = 5/3` the adiabatic index)
does not enter `delta W_T`, and compressible MHD with parallel flow has
`f_in = 3`. This rate agrees with Porcelli, Boucher and Rosenbluth, Plasma
Phys. Control. Fusion 38, 2163 (1996), Eqs. (1) and (8). Bussac's printed rate
is half of it (open question below), so only the threshold and the scaling are
used as checks.

Ohmic profile of this case (`s_B = 0.0188`, `b = 0.0378`, `c = 0.324`):
`beta_p1,crit = 0.206`, a global `beta_p` of 0.88.

| beta_p | beta_p1 | delta W_T | gamma tau_A, R0/a = 3 | 5 | 10 |
|---|---|---|---|---|---|
| 0.02 | 0.005 | +0.0287 | stable | stable | stable |
| 0.3 | 0.070 | +0.0235 | stable | stable | stable |
| 1 | 0.234 | -0.0071 | 3.1e-4 | 6.6e-5 | 8.3e-6 |
| 3 | 0.702 | -0.234 | 1.0e-2 | 2.2e-3 | 2.7e-4 |

Parabolic family `q = 0.9 + 0.1 (r/r1)^2` (`s1 = 0.200`, `s_B = 0.0180`,
`b = 0.0363`), conducting wall at `a`:

| r1 | q(a) | c | beta_p1,crit | delta W_T at beta_p1 = 0 / 0.3 / 0.6 | gamma tau_A at 0.6, R0/a = 3 / 10 |
|---|---|---|---|---|---|
| 0.30 | 2.01 | 0.276 | 0.232 | +0.0280 / -0.0158 / -0.127 | 7.7e-3 / 2.1e-4 |
| 0.35 | 1.72 | 0.234 | 0.264 | +0.0281 / -0.0072 / -0.097 | 8.0e-3 / 2.2e-4 |
| 0.40 | 1.53 | 0.202 | 0.294 | +0.0282 / -0.0010 / -0.076 | 8.1e-3 / 2.2e-4 |

The ideal rate exceeds the resistive m = 1 rate
`(s1/(r1 R0))^(2/3) eta(r1)^(1/3)` (Ara et al., Ann. Phys. 112, 443, 1978) in
the target range only at high `beta_p` and small `R0/a`. Elsewhere the sign of
`delta W_T` shows in the resistive growth: tearing-like scaling `S^(-3/5)` below
the threshold and resistive-kink scaling `S^(-1/3)` above it. Proposed checks:

- threshold from linear runs at `R0/a = 10` with an `S` scan: the `beta_p1` at
  which the `S -> infinity` limit of `gamma` departs from zero, equivalently
  where the scaling changes, at `beta_p1,crit` within 0.03 for the parabolic
  case `r1 = 0.3` and for the Ohmic profile;
- scaling, independent of the prefactor: `gamma R0^3/r1^2` constant within 10 %
  between `R0/a = 5` and 10 at fixed `beta_p1`, and
  `gamma(0.35)/gamma(0.30) = 1.04`, `gamma(0.40)/gamma(0.30) = 1.06` within 10 %
  in the parabolic family at `beta_p1 = 0.6`;
- plateau at `R0/a = 3`: the growth rate is independent of `S` between `1e8`
  and `1e9` within 5 %.

## Geodesic acoustic mode

An `m = 0`, `n = 0` poloidal flow compresses the plasma through the geodesic
curvature into an `m = 1` pressure perturbation that drives it back.
With the parallel sound flow included, the frequency at minor radius `r0` is

```text
omega_GAM^2 = (2 Gamma p / (rho R0^2)) (1 + 1/(2 q^2))
```

at `p`, `rho` and `q` of `r0` (Winsor, Johnson and Dawson, Phys. Fluids 11,
2448, 1968). Proposed check: `n = 0`, ideal, `beta_0 = 1e-3`, a small `m = 0`
flow at `t = 0`; the frequency of the local flow at `r0` over the first periods
agrees within 5 % for `r0/R0 <= 0.05`, with the difference decreasing in
`r0/R0`. The 5 % is a choice that covers finite-aspect-ratio corrections and
phase mixing, not a derived bound.

## Pressure check in the straight cylinder

The pressure terms that carry the toroidal limits above can be checked
without toroidal geometry. In a straight cylinder, localized interchange modes
at a resonant surface are stable only if `D_S = -2 p' q^2/(r B_z^2 q'^2) < 1/4`
(Suydam 1958, in the form of Freidberg, Rev. Mod. Phys. 54, 801, 1982,
Eq. 5.99). The criterion is necessary, not sufficient, and modes near the
threshold grow very weakly (Gupta et al., Phys. Plasmas 9, 3395, 2002), so a
check brackets it. Proposed check, translation symmetry with Fourier harmonics
in `z`: weak shear `r q'/q = 0.1` at a single resonant surface `q_s = 1.5`,
`r_s = 0.5`, for `m/n = 3k/2k`, `k = 1` to 4. At `D_S = 0.15` no harmonic grows
faster than `1e-5 per tau_A` over `1e4 tau_A` (a growth there is first checked
against Newcomb's criterion), and at `D_S = 1` the harmonics grow. These values
give `|p'| = 0.01 D_S` at `r_s`, a beta near 1 %.

## JOREK baseline

The matching JOREK run at `R0/a = 3`, `S = 3e6`, `beta_p = 1`, with toroidal
harmonics `n = 0, 1`, 60 radial and 48 poloidal elements and a step of
`5 tau_A`, takes about 3 s per step on 24 cores of one cluster node, or about
`1.4e4` core-seconds per `1000 tau_A`. At `R0/a = 10` with `n = 0..2` and a step
of `40 tau_A` it takes about `3.5e3` core-seconds per `1000 tau_A` on the same
24 cores. Steps of 20 and `40 tau_A` keep the saturated `R0/a = 10` state of
JOREK unchanged to four digits in `q0` over `5e4 tau_A`. RM-B9 compares with
JOREK on the same node as the MILCA run.

The step matters here. On a polar grid about the geometric centre, an explicit
step for the ideal terms at `R0/a = 3`, `beta_p = 1` was measured at
`0.021 tau_A`: the Shafranov shift puts the magnetic axis about `0.15 a` off the
grid origin, where the poloidal field is about 0.05 and the poloidal spacing of
that grid about a third of its radial spacing, so the poloidal Alfven rate is
about 22 per `tau_A`.

## Open questions

- The factor 2 in the Bussac growth rate. If Bussac's `2 pi^2 R0` prefactor
  belongs to the complex-amplitude convention, the rates above are twice too
  large. The layer analysis of Rosenbluth, Dagazian and Rutherford, Phys.
  Fluids 16, 1894 (1973) should settle it. Until then the absolute rates are
  not a check.
- The geodesic-acoustic frequency is taken in the form in which Winsor et al.
  are commonly quoted, not from the paper itself.
