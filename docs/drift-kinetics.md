# Drift kinetics: local neoclassical transport and NTV (DK)

Status: design decision, 2026-10-13, approved by Chris. No solver code yet.
DK0 (theory gate for the local ansatz with shear) is done for trapped
particles; see [the DK0 derivation](../theory/derivations/drift-kinetics-dk0.md).
The programme is gated on a numerics claim (DK1, below): if it fails, the
programme stops before the full operator.

## Principles

1. Not a collection of replacements of existing codes. NEO-RT, POTATO,
   MARS-K, PENTRC and the Shaing formulas stay external oracles; none is
   reimplemented. NEO-2's physics is kept; what changes is the numerics.
2. Same physics models, better numerics. The claim to be shown is accuracy
   per cost and robustness in the regimes where existing codes are known to
   fail (superbanana plateau, omega_E -> 0), not new physics.
3. One operator, modular collision operators. One local DK operator
   L[geometry, n, omega, collision model]; reductions (Krook, |B|-only,
   quasilinear, no shear, Lorentz) are switches with measured model defects,
   not named-code modes ([AGENTS.md](../AGENTS.md)). Collision operators come
   from the [collisions.md](collisions.md) ladder through one interface.
4. Geometry through D7 readback, not flux coordinates. The operator reads B,
   |B|, grad B, the metric and a field-line label from one kinetic geometry
   view of [equilibrium-core.md](equilibrium-core.md) D7. Boozer data are an
   export for oracle codes only.
5. Verification through the stage-1 NTV harness (iter_tc24 `ntv/`: frozen
   field, full-spectrum bounce reference, `decision.py`) plus Onsager and
   entropy-production checks ([adjoints-structure.md](adjoints-structure.md)).
   Derivations follow [theory-map.md](theory-map.md): reproduce, do not rely.

## Problem

Steady, linearised drift-kinetic equation for each species a,
f_a = F_M,a + f_a1, with F_M a local Maxwellian on psi and Phi_0(psi):

    v_par b.grad f1 + (v_E + v_d).grad f1 - sum_b C_ab^lin(f1)
        = -(v_d.grad psi)(A1 + x^2 A2) F_M + (Z e v_par B / T) A3 F_M,            (1)

A1 = dln p/dpsi + Z e Phi_0'/T - (3/2) dln T/dpsi, A2 = dln T/dpsi,
A3 = <E_par B>/<B^2> (NEO-2 drive convention), x = v/v_th. Ordering: slow
rotation, v_E ~ v_d ~ rho* v_th.

In scope: multispecies local transport (D_ij, bootstrap, Spitzer/Ware,
parallel flows) as the n = 0 problem; NTV (D_ij^NA, torque density per
species and n, offset rotation, ambipolar E_r) across 1/nu, nu-sqrt(nu),
superbanana plateau and bounce/transit resonances, including omega_E -> 0;
tangents and adjoints of the outputs.

Out of scope in this programme: finite orbit width beyond the first-order
shear term (POTATO remains the FOW oracle); islands and resonant layers at
q = m/n (detected and refused); non-nested or stochastic regions (refused);
sonic rotation (Mach defect reported); kinetic-MHD feedback into the L-lane
(kept compatible, not built); turbulence and time dependence.

## Formulation

(b) Default: full linear DK response of a 2D (axisymmetric) equilibrium to
the vector perturbation (delta B, xi, delta Phi) of the L-lane. The drive is
the perturbed guiding-centre Hamiltonian pulled back along unperturbed
motion in the ideal gauge delta A = xi x B0:

    delta H_n = mu delta B_L,n + m v_par^2 kappa.xi_n + Z e delta Phi_L,n,
    delta B_L = b.delta B + xi.grad B0,                                           (2)

and the torque density is the bilinear functional
T_phi,n = -(n/2) Im < int d^3v delta H_n^* f_n > (sign per COCOS 3), the
kinetic delta W_k form of [Park 2011](https://doi.org/10.1063/1.3662039).
Equation (2) and the torque form are stated from background knowledge and
are not yet verified; DK1 fixes signs against the stage-1 harness. The
NEO-2 two-stage structure is kept: f_10 (n = 0) first, then each n.

(a) NEO-2-style quasilinear: an explicit switch set on (b), not a second
code ([Kernbichler et al. 2016](https://doi.org/10.1088/0741-3335/58/10/104001),
[Martitsch et al. 2016](https://doi.org/10.1088/0741-3335/58/7/074007)). The
four reduction assumptions:

1. |B|-only drive: kappa.xi -> delta B_L/B, so delta H_n =
   (m v^2/2)(2 - lambda B) delta B_L,n/B (the NEO-RT H_n weight). This is a
   model change, not an identity (it needs b.delta B = 0 and no pressure
   term in kappa); its defect is measured.
2. delta Phi = 0.
3. Static perturbation, omega = 0.
4. Unperturbed surfaces and metric: delta B enters only through |B| on the
   unperturbed surface (Boozer, covariant B_s dropped).

Common to (a) and (b): amplitude small enough that unperturbed orbits are
the basis, and collisional broadening wider than nonlinear trapping. Their
violation is the superbanana/delta-B-trapping regime.

(c) 3D DK on 3D equilibria (DKES/SFINCS/MONKES class, map-lane states, finite
amplitude): the design stays compatible (same kernels on (theta, alpha,
lambda, x), same collision interface, same D7 view), but no increment is
planned for it.

## Local ansatz with finite shear

f_n lives on one surface; d/dpsi of f is dropped at fixed (theta, alpha) with
the geometry-agnostic label

    alpha = phi - nu(psi, theta),  nu = int_0^theta (B.grad phi)/(B.grad theta) dtheta',   (3)

so v_d.grad alpha contains the secular finite-shear term -(d_psi nu)
v_d.grad psi, d_psi nu = q' theta + periodic. For trapped particles its
bounce average equals q' <(I v_par/Omega) v_par b.grad theta>_b: the two legs
of a banana sit at different radii (psi* = psi - I v_par/Omega is conserved)
and see different q. "Local with finite shear" therefore means local in psi*
to first order in the orbit width, valid while n q' Delta_psi_b << 1. NEO-2's
`isw_mag_shear` is this iota' term in the tangential drift and NEO-RT's
`noshear` removes the same term from Om_tB
([Albert et al. 2016](https://doi.org/10.1063/1.4961084), Sec. VI).

DK0 gate (passed for trapped particles): the bounce-averaged precession
built from the curl form of v_d with the secular term reproduces NEO-RT's
`orbit.f90` Om_tB on a circular large-aspect Boozer field, shear on and off,
to 1e-11 relative to the deep-trapped precession
(`test/check_dk0_precession.py`). The ansatz is adopted for trapped
particles. The passing branch (cut convention of the gauge phase) is open
and gates the passing part of DK1.

## Variables and operator

- Space: surface psi_0; theta any poloidal angle provided by D7; alpha from
  (3). Per n: f_n(theta, lambda, x, sigma) exp(i n alpha) with twisted
  periodicity f_n(theta + 2 pi) = exp(2 pi i n q) f_n(theta).
- Velocity: x = v/v_th; lambda = v_perp^2/(v^2 B) (NEO-2's eta);
  sigma = sign v_par. lambda is conserved along the line, so every layer to
  be resolved (trapped-passing boundary lambda_c = 1/B_max, precession
  reversal lambda_r(x), bounce resonances) is a fixed lambda-line. The price
  is the curved accessible domain {lambda B(theta) <= 1}.

Operator per harmonic, with the mode frequency carried from the start:

    L_n f = sigma v |xi| b.grad theta d_theta f + i (n omega_alpha - omega) f - C(f),
    omega_alpha = omega_E + v_d.grad alpha,   |xi| = sqrt(1 - lambda B).        (4)

omega_E = v_E.grad alpha is rigid in psi. For n = 0 the second term
vanishes and (4) is the standard neoclassical problem. The Lorentz operator
in lambda, with the measure B dlambda/|xi| per sigma, is

    C_L f = (2 nu_D |xi|/B) d_lambda (lambda |xi| d_lambda f),
    int g C_L f (B/|xi|) dlambda = -2 nu_D int lambda |xi| d_lambda g d_lambda f dlambda,   (5)

symmetric and negative semi-definite, degenerate at lambda = 0 and at the
reflection curve. The weak form (5) is what is discretised: null space and
dissipation hold by construction (bracket form, [collisions.md](collisions.md)).
Status of (4)-(5): stated, V0 checks below.

Drifts: v_d.grad psi only in the drive; v_d.grad alpha (grad-B, finite-beta
curvature, shear term) local in the operator, not bounce-averaged (as
NEO-2-QL `isw_calc_magdrift`). All drift leaves come from FortSym in KIN6D
`gen/`, shared by operator, frequency placer and tracing.

## Discretisation

Pitch and theta: discontinuous Galerkin of order p = 2-4 on the curved
domain {lambda B(theta) <= 1}, upwind flux in theta per sigma, SIPG for the
lambda diffusion of (5). Anisotropic layers (thin in lambda, long in theta)
are placed a priori at lambda_c (width ~ sqrt(nu)), at lambda_r(x) and at
bounce-resonance lines predicted by bounce integrals of the same drift
leaves, and refined a posteriori by dual-weighted residuals for the torque
and D_ij functionals. NEO-2 refines only at collisional boundaries, never at
the reversal layer; that is the measured cause of its superbanana-plateau
failure (below).

Energy: nodal GLL elements in x on [0, x_max] (pattern `src/radial_dvr.f90`),
shared with the collision laboratory speed grid. The drift term is then
diagonal in x. Sonine/Laguerre stays as the field-part projection basis.

Justification required: DG/SIPG and GLL are adopted only if DK1 shows they
beat a simpler scheme at equal accuracy. DK1 therefore also runs one simpler
reference discretisation on the same problem, either conservative finite
volumes on lambda levels (NEO-2-like) or Legendre in xi = v_par/v
(SFINCS/MONKES-like, [Landreman et al. 2014](https://doi.org/10.1063/1.4870077),
[Escoto et al. 2024](https://doi.org/10.1088/1741-4326/ad3fc9)). If the
simpler scheme reaches the gate at comparable cost, it is kept and DG is
dropped.

Collision ladder (one `apply`/`jvp`/`diag_block` interface per spatial
point; field parts as low-rank updates inside the Krylov matvec, never by
Richardson iteration):

| Increment | Operator | Role |
|---|---|---|
| DK1-DK2 | Lorentz (pitch-angle scattering) | numerics gate, regime ladder |
| DK3 | Sugama model operator ([Sugama, Watanabe, Nunami 2009](https://doi.org/10.1063/1.3257907), metadata not yet verified) | production default, momentum and energy conservation |
| only if NEO-2 matching needs it | full linearised Landau | reference |
| parity runs only | Krook (incl. PENTRC's harmonic rate) | MARS-K/PENTRC matching, never a default |

Conservation defects are measured with the solver's moment functionals;
Lorentz rows are reservoirs and are not corrected.

## Solver

The DK problem is linear. Complex FGMRES with right preconditioning:

- Level 1: streaming, drift and Lorentz are block diagonal in nodal x. One
  sparse direct factorisation of the 2D (theta, lambda, sigma) system per
  energy node, batched and reused across right-hand sides (3 drives x
  species x n), instead of one 3D LU. Krylov carries energy scattering and
  the field part.
- Level 2: bounce-averaged coarse level. In the 1/nu limit the trapped
  solution is bounce-constant at leading order; the coarse space is the
  bounce-averaged equation per well (1D in lambda, Lorentz plus
  i n Omega_alpha). Target: iteration counts uniform in nu and in omega_E
  down to 0. Bounce-averaged models (Shaing, PENTRC-type) appear only here,
  as a preconditioner, not as a solver path.
- Null spaces by explicit constraints: the n = 0 density and energy null
  directions by bordering (<int f> = 0, <int x^2 f> = 0 with source
  multipliers), not by regularisation. For n != 0 at omega_E -> 0 the
  trapped near-null space (bounce-constant functions) is controlled by the
  coarse level with the passing-boundary condition. The claim that the
  bounce-averaged Lorentz problem with that condition is invertible is
  background knowledge [K] and must be demonstrated in DK1.

Newton enters only for the ambipolar E_r root and torque-balance coupling,
with exact tangents through the same factorisations; adjoints for torque,
D_ij and bootstrap reuse them (A^dagger = -T + D).

## Interfaces

- Geometry: one batched D7 kinetic geometry view (below, DK-pre-1). The DK
  code holds no private geometry. Refusal: no label (non-nested state) or a
  surface that is not closed means no local DK solve.
- Perturbation: in memory from the L-lane per n, evaluated at the DK
  quadrature points (DK-pre-2): Eulerian delta B_n, xi_n, delta B_L,n,
  kappa.xi_n, optional delta Phi_n and the mode frequency omega. The run
  refuses surfaces within the estimated layer width of a rational q = m/n.
  The `.bc` Boozer export exists only to feed oracle codes (NEO-2, NEO-RT).
- Rotation: subsonic only; omega in (4) and in the drive from the start.
- Outputs per surface, species, n: fluxes, <u_par B>, j_bs, L_ij, D_ij^NA,
  T_phi(psi), offset rotation, ambipolar E_r; torque per bounce harmonic for
  oracle comparison; adjoint sensitivities (torque response to delta B_n).
  Self-consistent kinetic feedback into the L-lane uses the same Hermitian
  torque form but is not part of this programme.

## Prerequisites (DK-pre)

These serve other lanes as well and come first.

| Item | Content | Owner |
|---|---|---|
| DK-pre-1 | D7 kinetic geometry view: batched B, \|B\|, grad\|B\|, kappa, metric, grad psi, field-line label nu and d_psi nu, q, q', flux-surface averages (coarea, with JVP/VJP) | KIN6D (D7), using FortFEM `flux_surface_average` |
| DK-pre-2 | L7 export: xi_n, Lagrangian delta B_L,n, kappa.xi_n in memory; distance to rationals; `.bc` for oracles | KIN6D (L-lane) |
| DK-pre-3 | Preconditioned complex (F)GMRES with right preconditioner and many right-hand sides | FortNum |
| DK-pre-4 | DG advection plus SIPG assembly on curved 2D domains | FortFEM |

FortSparse batched factorisations are used as they exist; gaps found in DK1
are repaired there.

## Programme and gates

- DK0 theory (done for trapped particles): local ansatz with shear,
  reproduction of NEO-RT Om_tB with shear on and off; (2) -> (a) drive
  algebra under switch 1. Open: passing branch.
- DK1 numerics gate. Lorentz collisions only.
  1. Axisymmetric neoclassical (n = 0, monoenergetic then Maxwellian
     weighted): bordered null space; D_11, D_31, D_33 against NEO-2 and the
     analytic limits; V0 operator checks.
  2. Stage-1 frozen-field NTV on the circular torus (iter_tc24 `ntv/`):
     1/nu and superbanana plateau against the full-spectrum bounce
     reference, scored by `decision.py` (|ratio - 1| <= max(0.10, 3 eps),
     trend intercept |a - 1| <= max(0.05, 2 sigma_a), >= 4 gated rows with
     eps < 0.05).
  3. The claim: beat NEO-2 where it failed. NEO-2's superbanana-plateau error
     grows like eps^-3 at fixed eta levels because the drift-reversal layer
     is not resolved, and the eta200 run costs 5.6e6 unknowns and 42 min of
     factorisation; at omega_E -> 0 with the drift off its flux grows like
     omega_E^-0.8 instead of reaching the 1/nu plateau (ratio 3.51 at
     omega_E = 1e-3; test-particle variant singular). DK1 passes only if the
     adaptive pitch FEM plus preconditioner (i) pass the SBP gate on the
     eps ladder where NEO-2 fails, (ii) give a finite torque that plateaus
     as omega_E -> 0 with nu-uniform and omega_E-uniform iteration counts,
     and (iii) do so at lower unknowns and wall time than NEO-2 eta200 for
     the same accuracy, with the simpler reference discretisation reported
     alongside.
  4. Stop rule: without a clear advantage in (i)-(iii), the programme stops
     here, before energy scattering, field parts and the vector drive. The
     result is reported as a negative finding.
- DK2 energy grid and test-particle operator, Krook parity switch,
  frequency-based layer placement, DWR estimators. Gates: full stage-1
  regime ladder; PENTRC Krook parity; NEO-RT collisionless limit.
- DK3 Sugama operator, multispecies, n = 0 standard neoclassics. Gates:
  NEO-2 D_ij/bootstrap/Spitzer; Onsager-Casimir L_ij(B) = eps_i eps_j L_ji(-B)
  from independent solves; L_sym >= 0 and entropy production; NTV with
  momentum conservation against NEO-2-QL. Full linearised Landau only if
  these comparisons require it.
- DK4 vector drive from the L-lane (needs DK-pre-2, L4/L5). Gates: Eulerian
  vs Lagrangian and compressive vs E x B fields from stage 1 (MARS KNTVEUL);
  the (a) switch set against NEO-2-QL fed by the same `.bc`; measured
  defects of the four (a) switches; rational-layer refusal.
- No increment for option (c), global DK, islands or kinetic feedback.

## Verification

- V0 operator: manufactured solutions for DG streaming on the curved domain
  with twisted periodicity (order p + 1); Lorentz eigenfunctions P_l(xi)
  with eigenvalue -nu_D l(l+1)/2; null-space dimension at n = 0;
  conservation defects; FortSym JVP/VJP pairing of the drift leaves.
- V1 stage-1 harness unchanged: a KIN6D DK run is one more root next to
  NEO-RT, PENTRC, NEO-2 and ARES.
- V2 collisionless limits: NEO-RT superbanana plateau (stage-1 ratios
  1.02-1.12, PASS), per bounce harmonic, shear on/off; POTATO reported as
  the finite-orbit-width bound of the local defect.
- V3 structure: Onsager-Casimir symmetry and nonnegative entropy production
  ([adjoints-structure.md](adjoints-structure.md)); momentum conservation
  on/off effect on 1/nu transport and flows.
- Result states follow [certification.md](certification.md): candidate
  through DK3; no certified claim is promised.

## Oracles and comparability

| Code | Model (as read) | Use |
|---|---|---|
| NEO-2 / NEO-2-QL ([Kernbichler et al. 2008](https://doi.org/10.1585/pfr.3.S1061), [Kasilov et al. 2014](https://doi.org/10.1063/1.4894479)) | full linearised collisions, field-line integration, quasilinear in delta\|B\| | n = 0 transport; (a) switch set where its eta ladder converges |
| NEO-RT ([Albert et al. 2016](https://doi.org/10.1063/1.4961084)) | collisionless resonant transport in action-angle form, Om_tB with shear | precession (DK0), SBP and resonance limits |
| POTATO | finite-orbit-width orbits | bound on the local defect |
| PENTRC ([Logan et al. 2013](https://doi.org/10.1063/1.4849395)) | bounce-averaged, Krook | Krook parity |
| MARS-K ([Liu et al. 2008](https://doi.org/10.1063/1.3008045)) | bounce-averaged kinetic delta W, Krook | Lagrangian/Eulerian drive, parity |
| Shaing ([Shaing 2003](https://doi.org/10.1063/1.1567285); [Shaing et al. 2009](https://doi.org/10.1088/0741-3335/51/3/035009)) | asymptotic regime formulas | leading-order limits only |

## Risks and limits

- The ITER benchmark payoff needs DK-pre-2 and the L-lane (L5, L7) first;
  until then DK runs only on prescribed and frozen fields.
- Finite orbit width beyond the first-order shear term and islands are out
  of scope; near rationals and at large n q' Delta_psi_b the local answer is
  refused, not approximated.
- The numerics claim may fail: if DG layers plus the coarse level do not
  beat NEO-2 on the SBP and omega_E -> 0 gates, the stop rule applies.
- NEO-RT's residual O(eps) SBP slope (2.05 eps in stage 1) limits it as a
  gate tighter than 10 % until explained.
- DG conditioning at the degenerate reflection boundary for p >= 3 is
  untested (NN2's anisotropy floor is the precedent).

## Open questions

- needs review: passing-particle shear term and cut convention (DK0).
- needs review: sign and omega term of the drive (2) and the torque form.
- needs review: invertibility of the bounce-averaged Lorentz problem with
  the passing-boundary condition, and whether the coarse level stays
  nu-uniform once the field part is on.
- needs review: Sugama 2009 operator metadata and equations from the paper
  body before DK3.

Chris&AI
