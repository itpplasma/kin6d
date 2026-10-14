# Plasma response: one operator family from ideal to kinetic (R)

Status: design proposal for review, 2026-10-13; no code change. Extends
[mhd-linear-response.md](mhd-linear-response.md) (L4, L6; L5 on branch
`gs-3d-resistive`) under [equilibrium-core.md](equilibrium-core.md) D5-D7,
with the same principles as [drift-kinetics.md](drift-kinetics.md). Evidence
from the iter_tc24 response ladder (rungs a-c, L5 TC24 array, 2026-10-13).

## Principles

1. Not a collection of replacements. GPEC/DCON/rdcon, MARS-F, MARS-K, ARES,
   MEPHIT and SLAYER stay external oracles; none is reimplemented. Their
   physics models are all reachable as settings of one KIN6D formulation.
2. Same physics models, better numerics. The claim is robustness and accuracy
   per cost at rational surfaces (no fit rings, no delta extrapolation, no
   Mercier-exponent bookkeeping), not new physics.
3. One operator, closures per region. One linearised per-n problem; ideal,
   resistive, viscous, two-fluid and kinetic are switches on the perturbed
   current closure, chosen per region (outer volume, layer strips), with
   measured model defects (as the DK switches (a) in drift-kinetics.md).
4. One outer operator (D5): FEEC Ampere, ideal volume current, exact exterior
   S_n (D6), factored once per mesh and n, shared by every model.
5. Code-vs-code only within one model lane (ideal static, resistive static,
   resistive rotating, kinetic). A cross-lane difference is a model result.

## Formulation

Per toroidal mode n, perturbation X(R,Z) e^{i(n phi - omega t)} of an
axisymmetric state from D7, plasma rotating toroidally with Omega(psi)
(rigid per surface), Doppler frequency w = omega - n Omega (static coil drive:
omega = 0, w = -n Omega). Unknowns: A (Nedelec, gauge R A_phi = 0), the
displacement xi where a region needs it, and the closure current J_pl.

    curl_n curl_n A = mu0 (J_coil + J_pl),   A outside: exact S_n (D6)          (R1)

    rho (-i w + Omega e_Z x)^2 xi = J_pl x B0 + J0 x curl_n A - div P1
                                    + div(rho nu grad v1) + (centrifugal)       (R2)

    i w G - Omega' G_phi grad psi - grad phi_c = eta J_pl
        + s_2f (J_pl x B0 + J0 x curl_n A - grad p_e1 - (n1/n0) grad p_i0)/(n_e e),
    G = A - xi x B0,  G_phi = R^2 grad phi . G                                  (R3)

    P1 = p1 I,  p1 = -xi.grad p0 - Gamma p0 div xi    (fluid)
    P1 = sum_a int d^3v (m v v) f_a1[xi, A, w]           (kinetic, DK)          (R4)

(R3) is Faraday plus Ohm in the co-moving frame (B1 = curl A); its ideal
limit eta = 0, s_2f = 0 is the frozen-in law G = grad chi, phi_c = i w chi:
A = xi x B0 with Phi = -sigma Omega xi.grad psi in the DK gauge
(drift-kinetics.md eq. 2; delta Phi_L = 0 there), A = xi x B0 -
(i sigma/n) grad(xi.grad psi) in the code gauge R A_phi = 0, where the
co-moving potential phi_c is an unknown of its own. The J0 x curl_n A term
freezes the flux into the electron fluid (Omega_e). v1 = -i w xi -
Omega' (xi.grad psi) R^2 grad phi; for the static drive w = -n Omega the
Coriolis and centrifugal terms are of the order of rho w^2 and are kept with
it; the s_2f row needs the ion gyroviscous stress in (R2) [K]. Signs follow
the per-n convention of `kin6d_coils_3d`. Derived in
[response-r2-r4.md](../theory/derivations/response-r2-r4.md) (T8), with
`test/check_response_r2_r4.py`; open points listed there.

Static ideal outer region (rho w^2, nu, eta, s_2f -> 0, fluid P1): (R2)
reduces to J_pl x B0 + J0 x B1 = grad p1 with B0.grad p1 = -p0' B1.grad psi and
div J_pl = 0. This is exactly the MEPHIT-form closure already in
`kin6d_mhd_ideal_n` (MDEs solved in the Fourier basis), so xi need not be
carried there; off rationals p1 = -xi.grad p0 and div xi = 0, so Gamma does
not enter. The MDE denominators (m - n q) are singular at every rational:
each rational needs a layer treatment (below), never a raw MDE.

### Switch table

| Setting (outer / layer) | rho w^2 | eta | nu | s_2f, omega_* | P1 | Class reproduced |
|---|---|---|---|---|---|---|
| vacuum + shielding sheets | - | - | - | - | - | L4 `shield` |
| ideal static / ideal sheet (A: L6 Frobenius) | 0 | 0 | 0 | 0 | fluid | GPEC/DCON singflx, MEPHIT ideal |
| ideal static / Delta' data (A) | 0 | 0 | 0 | 0 | fluid | rdcon, PEST-3, STRIDE matching data |
| ideal static / asymptotic resistive layer (A: L5) | 0 / on | 0 / on | 0 / on | 0 | fluid | rdcon + SLAYER-type, constant-psi FKR |
| ideal static / resolved resistive strip (B) | 0 / on | 0 / on | 0 / on | 0 / opt | fluid | MARS-F, ARES in the layer region |
| resistive everywhere (B with strip = plasma) | on | on | on | 0 | fluid | MARS-F, ARES (check switch) |
| drift-MHD layer (A or B) | on | on | on | on | fluid | SLAYER drift-MHD, Cole-Fitzpatrick |
| kinetic outer delta W_k | 0 | 0 | 0 | 0 | DK | MARS-K perturbative or self-consistent |
| kinetic layer (B) | on | from DK | - | in DK | DK | KiLCA, MEPHIT hybrid kinetic-MHD |

A "GPEC-type", "MARS-F-type" or "MARS-K-type" result is a row of this table
evaluated with one Ampere operator, one coil drive, one readback; there is no
code-named mode in the input, only the switches (AGENTS rule).

## Rational surfaces: two interchangeable treatments

Every rational q = m/n inside the plasma (all helicities coupled by
toroidicity, not only the listed ones) gets a layer object with one interface:
`layer_closure_t` returning its contribution to the outer system as
constraints on the layer-edge sheets (A) or as a Schur complement on two
interface rings (B).

### (A) Asymptotic matching closures (exists: L4, L6, L5)

Outer ideal solution near the rational, x = psin - psin_s,

    Phi_m = A |x|^(1/2 - mu) + B |x|^(1/2 + mu) + C x + ...,  mu = sqrt(-D_I)   (A1)

(Glasser, Greene, Johnson 1975). L6 sets A_+- = 0 (ideal); L5 imposes the
layer relation [Phi']/Phi_s = Delta_layer(w) on the same edge sheets, with the
FKR/Fitzpatrick scaled slab layer: constant-psi resistive-inertial
Dhat = 2 pi Gamma(3/4)/Gamma(1/4) g^(5/4)/eps, numerically otherwise.
Strengths: cheap, interpretable (Delta', D_I, P = -S/(r_s Delta')), exact
ideal limit. Measured failure modes (2026-10-13):

1. Fit-ring dependence. TC24 full pressure, q = 3: layer-current spread
   between ring grading 1.5 and 2.0 of 24-195 %; q = 10/3: 69-78 %, and
   14-105 % even at zero pressure.
2. mu != 1/2 with pressure (mu = 0.553 at q = 3, 0.747 at q = 10/3): S(delta)
   is not linear in delta; a delta^(2 eps) or (x_s/delta)^eps term (eps =
   mu - 1/2) is invisible to a two-point extrapolation. GPEC's singflx
   depends on `sing_spot` for the same reason.
3. Delta extrapolation itself: the large-solution residual |a1/Phi| outside
   is 17-620, i.e. the fits extract a small coefficient under a large one.
4. Missing sheets on sideband rationals (q = 3 near the edge with only q = 2
   listed): an unshielded 1/(u - u_s) pole of the strip current made the
   q = 2 result depend on the nearest half ring, up to +14 % (fixed on the
   sheet closure by `rat_auto`; not yet for the Frobenius closure).
5. Layers near the cut or near crowded rationals (q = 10/3 at psin 0.990,
   m +- 1 neighbours unlisted) do not converge at all.
6. The layer model is the cylinder-equivalent slab: no Glasser curvature
   (D_R), no toroidal sideband coupling inside the layer, no inertia
   enhancement beyond a caller factor.

Zero-pressure q = 3 converges (spread < 7 %), and the cylinder Delta' and
edge-current oracles pass, so (A) is correct where its assumptions hold.

### (B) Resolved layer strip (proposed)

The layer is solved, not matched. Around each rational a strip
|x| < X_s (X_s = 20-50 layer widths, capped at half the distance to the next
rational; overlapping strips merge) carries (R1)-(R4) with all layer terms
on, on a locally refined radial mesh graded geometrically from h_min ~
delta/10 at x = 0 to the outer ring spacing. Outside the strips the ideal
static outer region of (R1)-(R2) is unchanged. Coupling is by continuity of
the tangential A trace and normal flux on the two interface rings (mortar
when the bases differ), i.e. a Schur complement of the strip on the outer
system. No coefficient is fitted, mu and D_I enter nowhere, there is no
delta to extrapolate, and toroidal coupling of m +- k inside the strip is
resolved because the strip carries the poloidal resolution of the ring.

How (B) removes the failure modes: 1-3 (fits, mu, delta) do not exist;
4 disappears because strips are created for every rational from q(psin)
(the `rat_auto` list), so the outer MDE is never evaluated near any
m - nq = 0; 5 becomes a resolution question with an error estimate instead
of a non-converging fit; 6 is replaced by the full toroidal operator.
Remaining risks: interface placement error (the outer region drops
rho w^2, eta, nu, which are O((delta/X_s)^k) there; measured by moving X_s),
and conditioning at large S (scaled strip variables as in L5).

Ideal limit. As eta -> 0 at fixed w, (B) must reproduce the L6 sheet
current; this is an asymptotic check of (B) by (A), not a production route
for ideal runs (the width vanishes; (A) with eps = 0 is exact there).

### Layer widths and cost at ITER parameters

Widths in units of minor radius a with S = tau_R/tau_A, P_m = nu/(eta/mu0),
s magnetic shear (Fitzpatrick 1998; Cole and Fitzpatrick 2006):

| Regime | width scaling | ITER-like value (S ~ 1e8-1e9, w tau_A ~ 1e-3-1e-2) |
|---|---|---|
| static resistive (Alfven-resonance) | (n s)^(-1/3) S^(-1/3) | 0.5-1.5 mm (MEPHIT `resistive_layer_width` on TC24, m = 4-11) |
| resistive-inertial, rotating | (w tau_A)^(1/4) (n s)^(-1/2) S^(-1/4) | ~ 1e-3 a, i.e. 1-3 mm [estimate] |
| visco-resistive | P_m^(1/6) (n s)^(-1/3) S^(-1/3) | 1-2 x static [estimate] |
| ion gyroradius rho_i (5 keV D, 5 T) | - | ~ 3 mm |

Widths omit the geometric factor (R r_s/a^2)^(1/2). The S^(-1/3) width is
the w-independent constant-psi boundary: constant-psi RI needs w tau_A <<
(n s)^(2/3) S^(-1/3) ~ 1e-3, so at ITER w tau_A the layer is
non-constant-psi (T8, R6). The ITER fluid layer is at or below rho_i: single-fluid resistive layer
results are a model, and the two-fluid and kinetic rows are needed for
quantitative penetration. This is an order-of-magnitude statement from the
scalings, to be replaced by the R-B2 profiles.

Cost of (B) with geometric grading ratio r: rings per side ~ ln(h_out/h_min)
/ ln r ~ 25 at r = 1.2, h_out/h_min = 100, i.e. ~ 60-80 rings per strip. Cost
grows like ln S, not like a uniform-mesh S^(1/3). In the strip the poloidal
structure is dominated by the resonant m and its sidebands; with N_m = 20-40
straight-field-line harmonics and ~ 8 complex fields the strip is
block-tridiagonal with blocks of ~ 300, i.e. ~ 2e4 unknowns and ~ 1e9 flop
per factorisation per rational. For comparison MARS-F solves the whole plasma
with this structure (35 s on 1 CPU for the large-aspect-ratio rung c0). A
strip in the FEEC ring basis (512-2048 vertices per edge ring) would instead
cost ~ 1e6 unknowns per rational with dense ring blocks; hence decision 2.

### Solver and preconditioning (D5)

Outer: Ampere `K + C^H S^-1 C` (Hermitian positive definite) factored once
per mesh and n, as now. The full problem (outer closure + strips) is solved by
right-preconditioned complex FGMRES (FortNum, shared with DK-pre-3) with the
ideal operator as preconditioner: the ideal-static system (outer closure plus
L6 sheets, or delta W_n from the K-E1 leaves at the reference) times
block-diagonal strip LUs. At w, eta -> 0 the preconditioned operator tends to
the identity; near ideal or tearing marginality (W_free -> 0, Delta' -> 0-,
seen at q0 2.2, A <= 7) a few outlying eigenvalues are deflated (Arnoldi, as
MEPHIT does for its fixed-point iteration), not damped. Many right-hand sides
(coil phasings, w scans) reuse the factorisations.

### When to use which

- Ideal static response (GPEC lane, NTV inputs from ideal xi): (A) L6, with
  `rat_auto` and fit diagnostics; (B) not used.
- Delta', D_I, marginality and interpretation: (A).
- Resistive or drift-MHD response with rotation (MARS-F lane), rationals
  near the edge or crowded, finite pressure with mu != 1/2: (B) as
  production; (A) L5 as its constant-psi/asymptotic check.
- Kinetic layer: (B) only.
- Parameter scans over S, w: (A) for screening, (B) at selected points.

## Interfaces

- Equilibrium (D7): B, p, J, metric, q, q', <B.grad theta>, <|grad psin|^2>,
  SFL angles (PEST for oracle comparison only); the kinetic geometry view
  (DK-pre-1, in progress) also supplies the strip with |B|, kappa and
  flux-surface averages. Profiles n_e, T_e, T_i, Omega, Z_eff for eta (Spitzer
  with the same NRL ln Lambda as MEPHIT) and rho come with the case, not from
  the gfile. The response lane owns no geometry.
- L7 export (DK-pre-2, in progress): per n, Eulerian delta B, xi,
  Lagrangian delta B_L = b.delta B + xi.grad B0, kappa.xi, optional
  delta Phi and w. Outside the strips xi follows from the outer closure (MDE
  inverse); inside, from the strip unknowns, so the export is continuous
  through the layer and carries no 1/x pole. Distance to rationals and strip
  widths are exported; the local DK keeps refusing surfaces inside a strip.
- DK solver: (i) outer delta W_k: P1 of (R4) from the per-surface DK
  solution with the vector drive (DK4), entering (R2) as div P1 (returned:
  delta p_par,n and delta p_perp,n per surface and theta, the delta W_k
  density and its xi_par term; T8 R7); perturbative
  (one DK pass on the ideal xi, MARS-K perturbative, Liu et al. 2008) first,
  self-consistent by fixed point with the D5 preconditioner later. (ii)
  Kinetic layer: a strip DK with x as a resolved coordinate, k_par = k_par' x,
  supplying the parallel conductivity that replaces eta and s_2f in (R3)
  (KiLCA-class, Heyn et al. 2006). It reuses the DK kernels in (theta,
  lambda, x_v) and the collision interface; it is a new radial coupling, not
  the local operator.
- NTV: torque density T_phi,n from the DK Hermitian form on the L7 export
  (drift-kinetics.md), plus the electromagnetic layer torque of each strip,
  T_EM ~ Im(Delta_layer) |Phi_s|^2 (Fitzpatrick 1998), computed from the
  resolved strip fields, not from the asymptotic formula.
- Torque balance and rotation feedback: Omega(psi) enters only through w.
  A later outer Newton on Omega with NTV + T_EM against a viscous profile
  model uses tangents of the response (same factorisations, adjoint
  A^dagger); penetration bifurcations by continuation. Interface only here.

## Verification ladder

Rungs as in the iter_tc24 response plan; gates per model lane.

| Rung | Problem | Oracles | Applies to |
|---|---|---|---|
| a | vacuum fields of coils and sheets | Biot-Savart, I_m/K_m cylinder | all (done: L3 = GPEC to ~4 digits) |
| b-ideal | cylinder, one rational, external harmonic: S/(2m) | `cyl_response.py` shooting; GPEC LAR (0.2-0.5 %) | L6, (B) ideal limit |
| b-resistive | Delta' (step profile closed form), P = -S/(r_s Delta'), constant-psi Dhat, non-constant-psi and viscous scaled layer | closed forms; FKR (Furth et al. 1963); Fitzpatrick 1998 constants; rdcon Delta' (0.1-1.5 %) | L5, R-B1 |
| c | large-A circular torus, A = 3-20 | GPEC singflx and MARS-F ideal (J_MARS/J_GPEC 0.99-1.02 at A 4-20, rung c1); MARS-F/ARES resistive at equal eta, ROTE (ARES = MARS to 1e-5) | L6, R-B2 |
| d | shaped Solov'ev, several rationals, edge current and cut ladder | GPEC, MARS-F, MEPHIT per lane | all |
| e | TC24 at cuts 0.995/0.997/full | lane-matched, converged quantities only | all |
| NTV | stage-1 frozen field and its successor fed by L7 | iter_tc24 `ntv/` harness, `decision.py`; NEO-2 via `.bc` | L7 + DK |

Internal identities: Hess W = L4_n (D5); (B) -> (A) as S -> infinity in the
constant-psi regime; (B) with strip = whole plasma equals (B) with narrow
strips to the interface error; reciprocity (M Hermitian at w = 0, eta = 0);
energy balance (dissipated power of the strips = Re of the coil work).
Kinked or near-marginal cases (DCON W_free < 0.1) are report-only in every
gate.

## Increments and gates

- R-A: consolidate L4/L5/L6 behind `layer_closure_t`; automatic all-rational
  layers from q(psin) for every closure (Frobenius included; the fix
  552161a covers only the sheet closure); fit diagnostics (spread over
  `layer_ratio`, |a1/Phi|) reported per rational; three-point delta fit with
  the delta^(2 eps) term as the last change to the fit machinery. Gates: all
  `mhd_*` ctests unchanged; c1_kin6d parity with GPEC within the c1 gate;
  TC24 zero-pressure q = 3 spread < 7 %; full pressure reported, not gated.
- R-B1: resolved resistive strip in the periodic-cylinder chart. Gates:
  constant-psi Dhat and the L5 numerical layer (non-constant psi, viscous) to
  1 %; P = -S/(r_s Delta') on the step profile to 1 %; ideal limit to the L6
  sheet current to 1 %; interface error decreasing with X_s; FGMRES
  iterations flat for S = 1e6-1e10; rings per strip growing like ln S.
- R-B2: torus. Gates: large-A limit returns R-B1 within the rung-c O(eps^2)
  trend; MARS-F and ARES resistive at equal eta and ROTE on rung c (same
  extractor `mars_shield.py`, flux observable in one chart); ideal limit to
  GPEC singflx on rung c; TC24 full pressure q = 3 and q = 10/3 converged in
  strip resolution and X_s without any fit; two-fluid switch s_2f against
  SLAYER drift-MHD (Park 2022) in the cylinder.

  Increment 1 (2026-10-14, `test_mhd_closure strip`, ctest
  `response_strip_torus_coupling`):
  This gate couples a single-helicity, zero-beta cylindrical strip to the
  toroidal outer operator. It does not resolve toroidal sidebands or pressure
  inside the strip, and it does not expose a strip closure in `gs_response_3d`.
  The R-B1 strip is coupled to the FEEC torus by trace continuity on the two
  edge rings. Case: tanh step (q0 1.2, r1 0.5, width 0.02, zero beta),
  n = 1, 2/1 strip with edges at psin_s -+ xd, 3/1 and 4/1 L4 sheets.
  Rules declared before the first run: (a) eta = 0, gamma = 1e-6: strip net
  edge current / L4 sheet current, |ratio - 1| <= 1e-3 at R0/a 40, xd 0.02;
  (b) torus growth rate / R-B1 sharp-step cylinder rate at S = 1e6, 1e7,
  1e8 within 2 %. (b) was re-declared on 2026-10-14 (controller): the
  sharp-step cylinder (r_s Delta' = 5.000) is a different problem from the
  tanh profile solved (ODE r_s Delta' = 5.247, +3.9 % in gamma), so it does
  not test the coupling. (b) now compares with the cylinder rate of the same
  tanh profile (same strip coupling, cylinder tearing-ODE outer map), limit
  unchanged at 2 %, and requires the gap to shrink with nrad. Checks kept:
  coupled closure = Schur form (I - S M)^-1 S P0 (5e-4); the same coupling
  with the analytic sharp-step map returns the R-B1 rate (<= 1e-3, measured
  2e-4). The same-profile map first integrated the ODE through r_s for the
  outer edge; it now continues current-free across the strip, as the FEEC
  interior map does (same/step 3.9 %, matching (5.247/5)^(4/5)).

  nrad ladder, R0/a 40, xd 0.02, GMRES 3e-8, MUMPS (sCluster job 3225264,
  kin6d 8021456; times 11-29 s on 4 threads):

  | nrad | rings | (a) rel | (b) S 1e6 | S 1e7 | S 1e8 |
  |---:|---:|---:|---:|---:|---:|
  | 32 | 42 | 1.02e-3 | 9.46e-2 | 9.59e-2 | 9.63e-2 |
  | 48 | 56 | 1.02e-3 | 2.86e-2 | 2.90e-2 | 2.91e-2 |
  | 64 | 72 | 1.02e-3 | 1.60e-2 | 1.62e-2 | 1.63e-2 |
  | 80 | 86 | 1.01e-3 | 9.60e-3 | 9.72e-3 | 9.75e-3 |
  | 96 | 103 | 1.01e-3 | 6.36e-3 | 6.45e-3 | 6.47e-3 |
  | 112 | 120 | 1.01e-3 | 4.47e-3 | 4.53e-3 | 4.54e-3 |
  | 128 | 133 | 1.01e-3 | 3.32e-3 | 3.37e-3 | 3.38e-3 |
  | 160 | 163 | 1.01e-3 | 2.06e-3 | 2.09e-3 | 2.09e-3 |

  (b) converges to the same-profile cylinder at about second order in
  1/nrad (64 -> 128: 4.8x; 80 -> 160: 4.6x) and passes from nrad 64, which
  is the gated resolution. R0/a 20 / 80 at nrad 96: (b) 9.0e-3 / 5.7e-3;
  xd 0.01 / 0.04 at nrad 96: (b) 6.6e-3 / 6.2e-3; xd 0.01 at nrad 64:
  1.65e-2.
  (a) does not close under refinement: 1.02e-3 -> 1.01e-3 over nrad 32-160,
  independent of the GMRES tolerance (1.04e-3 at 1e-7, 1.02e-3 at 3e-8). It
  scales with the strip width: 1.55e-4 (xd 0.01, nrad 64-128), 1.01e-3
  (0.02), 4.40e-3 (0.04), i.e. about xd^2-xd^2.7, with a weak R0/a trend
  (6.2e-4, 1.01e-3, 1.11e-3 at R0/a 20, 40, 80). It is a finite-strip
  truncation of the comparison with a sheet on the rational (presumably the
  linear-F local equations of the strip; not isolated), not a mesh or
  coupling error. Under the first declaration (xd 0.02) gate (a) failed by
  1 %. (a) was re-declared on 2026-10-14 (controller) because its error is
  the finite strip width relative to the sheet limit, which is what (a)
  tests: at xd 0.01, |ratio - 1| <= 1e-3 (limit unchanged), and the order
  of that error between xd 0.02 and 0.01 at the same nrad must be >= 1.8
  (measured 1.01e-3 / 1.55e-4: order 2.7). The gated run is therefore at
  xd 0.01, nrad 64, plus an (a)-only run at xd 0.02 for the order. (b) at
  xd 0.01 (job 3225264): nrad 64 1.65e-2 / 1.63e-2 / 1.64e-2, nrad 96
  6.6e-3, nrad 128 3.4e-3, so nrad 64 stays the gated resolution.
- R-K: kinetic. K1 outer delta W_k (perturbative, then self-consistent)
  after DK4 passes, against MARS-K perturbative on rung c; K2 kinetic strip
  against KiLCA and MEPHIT hybrid kinetic-MHD (Lainer et al. 2026) in the
  cylinder. Stop rule: if DK1 stops, R-K stops; R-B results stand alone.

Not done: no xi-form Newcomb/DCON shooting code, no Fourier-radial
whole-plasma MARS-F clone, no SLAYER or PENTRC reimplementation (oracles
only), no resistive wall or RWM eigenvalue solver (operator kept compatible,
D6 S_n can carry a thin wall later), no nonlinear islands or time-dependent
penetration, no 3D (stellarator) equilibria, no MEPHIT-style artificial
k_par damping (it is a layer width without a layer model).

## Full-pressure strip residual

`kin6d_response_fluid` supplies the non-rotating, single-fluid point residual
of (R2)-(R4), in Cartesian physical components and SI units. With
`gamma = -i omega`, the supplied spatial jets give

    p1 = -xi.grad p0 - Gamma p0 div xi,
    momentum = rho gamma^2 xi - J1 x B0 - J0 x B1 + grad p1,
    Ohm = -gamma A1 - grad Phi1 + gamma xi x B0 - eta J1.

The pressure gradient keeps both the equilibrium pressure Hessian and
`grad(div xi)`. This admits parallel compression and pressure forces in curved
fields without an inserted Mercier coefficient or fitted layer parameter.
`J1 = curl B1/mu0` and `B1 = curl A1` remain the magnetic discretisation's
responsibility. Electrical resistivity here is in Ohm m, unlike the normalised
`eta = 1/S` of the cylindrical strip. The input jets must include chart and
basis derivatives before conversion to Cartesian components.

The source is `gen/gen_response_fluid.f90`, which uses FortSym to emit the
pressure and force leaves and their real JVP/VJP products. Regenerate from the
repository root with `gen_response_fluid` linked to FortSym; the generator
applies the project formatter. The complex residual is linear in perturbation
jets at fixed equilibrium and gamma, so its perturbation JVP is the same
action. Parameter sensitivities of a future strip solve require differentiating
the discretised equation and its interfaces; point products alone do not supply
those sensitivities. The batch action uses OpenMP with fixed scratch arrays.

CTest `response_fluid_full_pressure_oracles` checks fast and slow magnetosonic
dispersion at nonzero pressure, resistive Alfven damping, translation of a
curved pressure-balanced equilibrium, pressure gradients against independently
differenced polynomials, batch linearity, and pressure/force JVP/VJP products.
These are point-operator checks. No continuous or discrete error certificate
is available, and no TC24 convergence or Glasser-layer result is claimed.

`kin6d_response_toroidal_strip` adds a candidate mapped-annulus discretisation
at fixed nonzero toroidal n, in the outer operator's `A_phi = 0` gauge. Its six
unknowns are `(A_R, A_Z, xi_R, xi_phi, xi_Z, Phi)`. A complete odd-size poloidal
Fourier grid retains both signs of every represented harmonic; equilibrium and
metric coefficients couple these harmonics. Five-node local polynomial radial
derivatives come from FortNum, and the Fourier transform comes from FortFEM.
The supplied inverse map transforms derivatives to physical R/Z coordinates;
cylindrical basis derivatives are included before the Cartesian point kernel.
The caller must supply the equilibrium pressure gradient and Cartesian Hessian.
There are no rational denominators, fitted exponents or single-helicity kernels.

The isolated strip prescribes tangential A and normal displacement on both
rings, retaining the other four equations there. FortNum GMRES solves the
resulting linear equation with explicit caller-supplied equation scales and
rechecks its residual. The perturbation JVP is the same linear action; the
storage VJP composes the generated point VJPs with transposed spatial actions.
Its pairing is the stored-grid Hermitian pairing. A physical weighted adjoint
needs the volume-quadrature mass factors. RHS tangents solve the equation again.
Equilibrium/mesh/parameter derivatives and a validated adjoint inverse are
unavailable. No continuous or discrete error certificate is supplied. Nodal
`div B` and `div J` defects are reported in SI units because this collocation
action does not have FEEC's exact commuting-complex identities.

Registered CTests check analytic forcing on a pressure-balanced curved toroidal
equilibrium, all represented poloidal modes at n = +/-3, the forced state and
drive tangent, JVP/VJP pairings with and without trace rows, zero-divergence
analytic identities, and radial refinement against an independently
differentiated degree-six electrostatic field. These establish an isolated
candidate strip, not a TC24 layer-current or Glasser result. Representative
small-eta/gamma conditioning, a scalable preconditioner, TC24 equilibrium-jet
sampling and coupled-interface refinement are still required.

### Two-ring traces and cut balance (candidate increment)

`toroidal_strip_traces` samples both strip rings with one normal
`N = grad u / |grad u|` and poloidal tangent `T = (-N_Z, 0, N_R)`.
Its ten entries are `(A.T, A_phi, B1.N, E1.T, E1_phi, xi.N, traction(3), Ptot1)`.
Here `E1 = -gamma A1 - grad Phi1`, `Ptot1 = p1 + B0.B1/mu0`, and

    traction = (B0 (B1.N) + B1 (B0.N))/mu0 - Ptot1 N.

The pressure uses the retained full `p0` in the adiabatic leaf. The point stress
and cut-balance leaves and their real JVP/VJP are generated by
`gen/gen_response_fluid.f90`. At fixed background/chart/gamma the complex trace
JVP is the linear trace action; `toroidal_strip_traces_vjp` composes the generated
pressure/stress VJPs with transposed spatial actions. The stored-grid Hermitian
pairing applies. Physical surface weights multiply the seed, and a physical
volume adjoint additionally needs the volume Riesz map. Hot loops use OpenMP
and existing strip workspace. No full trace Jacobian is formed.

An outer trace must use the same normal before subtracting the two sides.
Outward-normal tractions instead add. These fixed artificial-ring traces assume
no imposed sheet current or traction. They supply the physical quantities for
coupling; they do not yet impose their continuity in the FEEC outer solve.

For a physical flux-surface cut (`B0.N = 0`), `response_cut_balance` evaluates
the Lagrangian normal-pressure jump

    Ptot1_inner - Ptot1_outer + xi_n d_n(Ptot0_inner - Ptot0_outer).

Its fixed-gradient transpose is available. The background traction balance
must also hold, or an external support must be declared. Continuous `B0` with
vacuum outside and retained positive `p_edge` has a background pressure defect
of `p_edge`; the operator does not remove that pressure or make the imported
cut an equilibrium. Full moving-surface geometry and support modelling remain
the outer boundary consumer's responsibility.

`toroidal_strip_domains` clips supplied rational-centered extents to the physical
domain and merges overlapping or touching intervals, returning ownership for
every supplied rational. Near-edge and exactly-on-cut centers are retained;
centers outside the domain are rejected rather than silently omitted. Topology
changes have no smooth derivative. The setup API does not discover rationals,
sample equilibrium jets or create boundary conditions for a strip at a cut.

Three CTests check independent analytic two-ring traces for both signs of n,
trace JVP/VJP, the full nonlinear Maxwell-stress derivative, a finite-pressure
moving-boundary pressure oracle, and domain coverage against direct membership.
The full-cut q=11/3 center at psin=0.999036564 is covered by the merged-edge
domain oracle. These are candidate operator tests with unavailable enclosures;
they do not qualify a coupled TC24 response or calibrate GPEC singflx.

The full-pressure R-B2 application still needs enforcement of electric and
mechanical interface balance and the two-ring magnetic Schur/mortar coupling.
`gs_response_3d` currently has only sheet and asymptotic closures.
Its cut-edge pressure jump is also explicitly outside the ideal outer model;
TC24 pressure-retaining cuts need consistent total-pressure boundary balance.
The increment-1 test and the new point residual cannot substitute for these
missing pieces.

## Decisions for Chris

1. Production treatment for resistive and kinetic response: (B) strips;
   (A) kept as ideal production and asymptotic oracle. Recommended. The
   alternative (patch (A) further: more fit rings, nonanalytic
   extrapolation, curvature terms) chases fit artefacts per rational.
2. Strip poloidal basis: straight-field-line Fourier (N_m 20-40, mortar to
   the FEEC rings) versus refined FEEC rings. Recommended: Fourier in the
   strip, decided by measured cost in R-B2; the outer region stays FEEC.
3. Two-fluid switch in R-B2 (s_2f, omega_*), before kinetics, because the
   ITER layer width is at or below rho_i. Recommended.
4. Strip = whole plasma (MARS-F-like) only as a verification switch, not a
   production mode. Recommended.
5. Rotation feedback and torque balance deferred until R-B2 passes.
   Recommended.

## References

- Glasser, Greene, Johnson, Phys. Fluids 18, 875 (1975), https://doi.org/10.1063/1.861224
- Furth, Killeen, Rosenbluth, Phys. Fluids 6, 459 (1963), https://doi.org/10.1063/1.1706761
- Pletzer, Bondeson, Dewar, J. Comput. Phys. 115, 530 (1994) (PEST-3 outer
  matching data), https://doi.org/10.1006/jcph.1994.1215
- Glasser, Phys. Plasmas 23, 072505 (2016) (DCON), https://doi.org/10.1063/1.4958328
- Glasser, Wang, Park, Phys. Plasmas 23, 112506 (2016) (resistive DCON,
  matched asymptotics), https://doi.org/10.1063/1.4967862
- Glasser, Kolemen, Glasser, Phys. Plasmas 25, 032507 (2018) (Riccati ideal
  response; STRIDE lineage to be confirmed in its source),
  https://doi.org/10.1063/1.5007042
- Park, Logan, Phys. Plasmas 24, 032505 (2017) (GPEC), https://doi.org/10.1063/1.4977898
- Park, Phys. Plasmas 29, 072506 (2022) (SLAYER layer regimes), https://doi.org/10.1063/5.0093079
- Liu, Bondeson, Fransson et al., Phys. Plasmas 7, 3681 (2000) (MARS-F),
  https://doi.org/10.1063/1.1287744
- Liu, Chu, Chapman et al., Phys. Plasmas 15, 112503 (2008) (MARS-K),
  https://doi.org/10.1063/1.3008045
- Fitzpatrick, Phys. Plasmas 5, 3325 (1998), https://doi.org/10.1063/1.873000
- Cole, Fitzpatrick, Phys. Plasmas 13, 032503 (2006), https://doi.org/10.1063/1.2178167
- Albert, Lainer, Biro et al., Comput. Phys. Commun. 277, 108401 (2022)
  (MEPHIT FEEC form), https://doi.org/10.1016/j.cpc.2022.108401
- Markl, Ulbl, Albert et al., Nucl. Fusion 63, 126007 (2023), https://doi.org/10.1088/1741-4326/acf20c
- Lainer, Zenz, Markl et al., Plasma Phys. Control. Fusion 68, 055037 (2026) (hybrid
  kinetic-MHD), https://doi.org/10.1088/1361-6587/ae681c
- Heyn, Ivanov, Kasilov et al., Nucl. Fusion 46, S159 (2006) (KiLCA),
  https://doi.org/10.1088/0029-5515/46/4/S07

DOIs, volumes and pages checked against Crossref (2026-10-13); width and cost
numbers marked [estimate] are scalings, not runs.
