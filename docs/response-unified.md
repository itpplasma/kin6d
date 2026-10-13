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

    -rho w^2 xi = J_pl x B0 + J0 x curl_n A - div P1 - i w div(rho nu grad xi)  (R2)

    i w (A - xi x B0) = eta J_pl - s_2f (J_pl x B0 - grad p_e1)/(n_e e)        (R3)

    P1 = p1 I,  p1 = -xi.grad p0 - Gamma p0 div xi    (fluid)
    P1 = sum_a int d^3v (m v v) f_a1[xi, A, w]           (kinetic, DK)          (R4)

(R3) is Faraday plus Ohm in the co-moving frame (B1 = curl A); its ideal
limit eta = 0, s_2f = 0 is the frozen-in law A = xi x B0, the ideal gauge used
for the DK drive (drift-kinetics.md eq. 2). Coriolis, centrifugal and
flow-shear terms of MARS-F (Liu et al. 2000) are additional explicit terms in
(R2), stated here from background knowledge [K] and added only with a gate.
Signs follow the per-n convention of `kin6d_coils_3d`; (R2)-(R4) are not yet
derived in theory/ and need the T-ledger check of
[theory-map.md](theory-map.md) before R-B1 code.

Static ideal outer region (rho w^2, nu, eta, s_2f -> 0, fluid P1): (R2)
reduces to J_pl x B0 + J0 x B1 = grad p1 with B0.grad p1 = -p0' B1.grad psi and
div J_pl = 0. This is exactly the MEPHIT-form closure already in
`kin6d_mhd_ideal_n` (MDEs solved in the Fourier basis), so xi need not be
carried there. The MDE denominators (m - n q) are singular at every rational:
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

The ITER fluid layer is at or below rho_i: single-fluid resistive layer
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
  solution with the vector drive (DK4), entering (R2) as div P1; perturbative
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
