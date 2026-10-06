# Q7: static screening hierarchy and Wigner -> Vlasov

Campaign case Q7 (docs/hierarchy-campaign.md). Two independent slices share
the species record, the finite-T Fermi EOS and the free-streaming leaf.
Status labels follow AGENTS.md: *bound* is rigorous, *indicator* is not.

## Part 1: Debye -> finite-T Fermi -> Lindhard (`src/screening_hierarchy.f90`)

Model: homogeneous ideal 3D spin-two electron gas with a neutralizing
background, static linear response (ideal-gas RPA). A static local-field
correction, correlations and the dynamic response are out of scope.
Inputs `(n, T)` in the caller's reduced units map through `kin6d_species`
to `Theta, r_s, Gamma, lambda_D`; `kin6d_fermi_gas` supplies `mu` and
`k_s^2 = (e^2/eps0) dn/dmu` (charge `sqrt(g)`, `eps0 = 1/(4 pi)`).
Units inside: `kF = 1`, `EF = 1`, `x = q/(2kF)`,
`G(x) = -chi0(q)/N(0)`, `eps = 1 + kappa2 G(q/2)/q^2`,
`kappa2 = (kTF0/kF)^2 = 4/(pi kF aB)`.

- T = 0: closed-form Lindhard `F(x)` from four generated charts (series for
  `x <= 0.1` and `x >= 10`, two log charts without `abs`).
- Route A (momentum): `G = (1/2x) int r f(r) log|(r+x)/(r-x)| dr`.
- Route B (Maldague): `G = int w(e) sqrt(e) F(x/sqrt e) de`, `w = -df/de`,
  `e = s^2`, over the exact T=0 response.
- Both use composite Gauss-Legendre (12 and 20 nodes per panel) with
  geometric grading toward the log singularity (A) or the 2kF kink (B).
  Grading stops 4096 ulp short of the singular point. `|Q20 - Q12|` is an
  indicator.
- Truncation bounds (emitted as outward interval leaves):
  route A, `R >= 2x`: `f <= exp((mu-r^2)/theta)` and `log <= 4x/r` give
  `theta exp((mu-R^2)/theta)/R`. Route B, `smax^2 > mu`:
  `w <= exp(-(e-mu)/theta)/theta`, `F <= 1` and the tangent bound of
  `sqrt(e)` give `exp((mu-smax^2)/theta)(smax + theta/(2 smax))`.
- Screened test charge: `phi/phi_bare = 1 - (2/pi) int sin(qr) kappa2 G/(q(q^2+kappa2 G)) dq`.
  The q-truncation bound `4 kappa2/(3 pi qmax^4)` uses `G(x) <= 2/(3x^2)` for
  every theta. Lemma: for `y >= sqrt 2`, `F = sum y^(-2n-2)/((2n+1)(2n+3)) <=
  1/(3(y^2-1)) <= 2/(3y^2)`. On `[sqrt(2/3), sqrt 2]`, F decreases, and a
  400-cell partition check in the test gives `F(y_a) <= 2/(3 y_b^2)`. Below
  that interval the bound is at least 1, which already dominates F. The
  Maldague average and `int w e^(3/2) = (3/2) int f sqrt(e) = 1` (fixed
  density) carry the lemma to finite theta. The G indicator propagates to
  first order.
- Derivatives: the generated leaves carry JVP/VJP. `d(k_s)/d theta` comes
  from `fermi_gas_jvp`. A theta-derivative of the full Lindhard
  `G(x, theta)` has no consumer yet and is not provided.

Evidence (`test_screening`, `app/screening_demo.f90`, r_s = 1):

| check | result |
|---|---|
| route A at T=0 vs closed form | 4.7e-15 |
| route A vs route B, theta 0.01-10, x 0.05-4 | 8.7e-15 (indicators < 1e-11, tails < 1e-15) |
| G(0) vs Sommerfeld at theta = 0.01 | 2.1e-8 (O(theta^4)) |
| theta^2 approach to T = 0 (deviation ratio, theta 0.02/0.01) | 4.007 |
| theta = 400 vs Boltzmann Dawson response D(y)/y | 3.3e-5 |
| constant-G inversion vs exact Yukawa | 5.0e-7, below its own q-tail 5.6e-6 |
| Friedel period, 20 <= r kF <= 60 | 3.1382/kF (pi/kF = 3.1416, 0.1 %) |
| theta = 10, phi/phi_bare vs exp(-k_s r), 1 <= k_s r <= 4 | 7.2e-4 relative |
| G(0) inside the certified fermi_gas `k_s` enclosure | all theta |

Finding (resolved): the `kin6d_fermi_gas` *candidate* `k_s/kTF0` was off by
up to 1.8e-3 (relative) at `theta ~ 0.025`, inside its certified 1% enclosure.
The cause was one Gauss rule across the Fermi window; the graded mesh split
at the Fermi surface (docs/fermi-gas.md) brings `G(0)/k_s - 1` at
`theta = 0.01` from 8.3e-6 to 3e-10.

## Part 2: 1X1V Wigner-Poisson (`src/wigner_poisson.f90`)

Electrons (`q = -1`, `m = 1`) with a neutralizing background, periodic in x;
units `omega_p`, `v_ref`, `lambda = v_ref/omega_p`,
`hbar = hbar omega_p/(m v_ref^2)` (= `hbar omega_p/T` for `v_ref = v_th`).
The grid is Fourier in x and periodic Fourier in v (eta). Strang step:

1. stream dt/2 exactly in (k, v) through the shared
   `kin6d_streaming_leaves:streaming_rotation` leaf;
2. Poisson: `Vhat_k = nhat_k/k^2` (`V = -phi`, `phi'' = n - 1`);
3. kick dt exactly in (x, eta): `fhat *= exp(i dt sum_k s_k(eta) Vhat_k e^{ikx})`,
   `s_k = 2 sin(hbar k eta/2)/hbar`, which is the shift `v +- hbar k/2`. The
   Vlasov limit uses `s_k = k eta` on the same grid;
4. stream dt/2.

Nyquist modes are filtered, and nothing else is. The Wigner function is
never clipped; its negative mass is a diagnostic. The forward tangent
`d/d hbar` propagates through every substep (`wigner_kick_tangent`,
`wigner_phase_coefficient_jvp`), so the discrete JVP matches finite
differences.

The independent oracle (`src/wigner_dispersion.f90`) shares no formula with
the solver: `eps = 1 - (1/(k^3 hbar)) sum_b w_b/(sqrt2 sigma_b)
[Z(xi+) - Z(xi-)]`, `xi+- = (omega/k +- hbar k/2 - u_b)/(sqrt2 sigma_b)`.
Z is integrated from `Z' = -2(1 + zeta Z)`, `Z(0) = i sqrt(pi)`, which gives
the Landau continuation. Newton finds the root, and the hbar-derivative is
the equation-level implicit `-eps_hbar/eps_omega`. A linear-prediction fit
extracts frequencies from the runs after ballistic phase mixing.

| check | result |
|---|---|
| classical Landau k=0.5 root vs literature 1.4157-0.1534i | < 1e-4 |
| Vlasov / Wigner (hbar=2) runs vs roots, k = 0.5 | 7e-5 |
| cool quantum Langmuir (sigma 0.1, k 0.4-1.4, hbar 0-2) runs vs roots | <= 1.0e-4 |
| Bohm-Gross `1 + 3k^2 sigma^2 + hbar^2 k^4/4` vs root | asymptotic, <= 3.3 % at k = 1.4, hbar = 2 |
| exact cold limit (sigma 0.02) vs Bohm-Gross | < 2e-4 |
| two-stream growth (beams +-1, sigma 0.2), 26 unstable (k, hbar) | <= 1.7e-4 |
| hbar=1 band extends to k = 1.4; hbar=2 stabilizes k >= 0.9 | cold quartic agrees |
| `||f_W - f_V||_2` vs hbar (0.8 ... 0.05) | rate 1.991 (two-stream t=20), 1.999 (Landau t=10) |
| d gamma/d hbar (k=0.6, hbar=1): tangent / implicit / FD | -0.03708 / -0.03737 / -0.03737 |
| discrete tangent vs central FD of the run | < 1e-6 relative |
| Strang order (dt 0.1 -> 0.05) | 2.07 |
| mass drift (movie run, 128 x 512, dt 0.039, t = 45) | 1e-14; momentum 5e-7 (W) / 6e-4 (V, unresolved after t ~ 24.5); dE/E 3.2e-6 (W) / 1.6e-5 (V) |
| negative mass at t = 20, hbar 1 / 0.5 / 0 (128 x 512) | 0.262 / 0.085 / 1.3e-13 |

The tangent growth slope uses the finite window [24, 40], so it keeps a
residual transient (0.8 %). The root derivative is the oracle. The Vlasov
negative mass is pure numerics: it stays at roundoff until the nonlinear
filamentation outruns the grid (`t ~ 22` at 64 x 256, `t ~ 24.5` at the
movie grid 128 x 512; refining v alone does not move it), so it serves as
the resolution indicator. Momentum is conserved only up to the
Nyquist/periodic-v truncation once the dynamics is nonlinear. All errors
in this part are indicators; no finite-basis residual bound exists yet.

## Outputs and cost

`screening_demo` and `wigner_demo` accept `--config FILE` (namelists
`&screening rs, qmax, theta_eps, theta_fr /` and
`&wigner dt, t_movie, hbar_movie /`). They write `output/<app>.nc` and CSV
files, plus `output/screening-hierarchy.png`,
`output/wigner-dispersion.png`, `output/wigner-two-stream.png` and
`output/wigner-two-stream.mp4`. The movie runs 13 s, 312 frames at 24 fps,
to t = 45 (linear growth, saturation near t ~ 28, trapped-particle vortex),
on a 128 x 512 grid with split steps <= 0.05. It shows f_W | f_V on a
coolwarm map centred at zero, the negative part of f_W on its own colour
scale, and |n_1| plus negative-mass histories with the time at which the
Vlasov indicator exceeds 1e-6 (grid under-resolved after that).

Cost (Apple M-series, gfortran 15, Release, one core). `screening_demo`
takes 13 s: k_s scan 0.6 s, epsilon 0.1 s, potentials 9.1 s, Friedel 3.3 s.
`wigner_demo` takes 113 s: Langmuir scan 23 s, two-stream scan 14 s,
hbar-limit 6 s, error ledger 11 s, movie runs 13 s, rendering the rest.
One Wigner step with diagnostics on 64 x 256 costs about 1 ms. Memory is the
`nx x nv` state plus four same-size work arrays (phase table, shift symbol
and their tangents). CTest: `test_screening` takes 3.5 s and `test_wigner`
2.2 s in Release; in Debug they take 7 s each.

## Open

- Generic pieces belong upstream: the plasma dispersion function Z (FortNum
  special functions); graded composite Gauss rules (FortNum quadrature).
- The screening figure uses per-panel log axes (fortplot 5380eb6). Panel 3
  floors |phi/phi_bare| at 1e-8 and panel 4 clips to its y range for display
  only; the CSV/NetCDF data are unclipped.
- No rigorous enclosure for the Lindhard quadrature (singular integrands) or
  for the Wigner run; a static local-field correction and an interacting
  RPA/HF comparison (QG7 next step) are not implemented.

## Error-control cost

Optional `screening_q_table` profiling separates low-rule estimation and
strict tail evaluation from the high-rule candidate. On the measured theta=1
fixture the diagnostic work exceeds candidate cost; it misses the always-on
10–20% target. Cache immutable response tables or schedule reference
checks at checkpoints. The quadrature and linear propagation terms remain
indicators; strict tail arithmetic alone does not certify the full response.
See [the measured strategy](error-control-strategy.md).
