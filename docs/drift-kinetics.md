# Drift kinetics: local neoclassical transport and NTV (DK)

Status: design decision, 2026-10-13, approved by Chris; local 3D and global
steady extensions approved 2026-11-03; full NEO-2 feature parity approved
2026-11-06 (tracker #198). No solver code yet.
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

Finite orbit width beyond the first-order shear term is admitted only in
the K2 extension below; the local ansatz keeps its ordering and refusal
conditions. POTATO remains an external finite-orbit-width (FOW) oracle.
Out of scope: islands and resonant layers of the perturbative NTV problem at
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
lambda, x), same collision interface, same D7 view). It is now the K1
extension below; the axisymmetric harmonic operator (4) alone is insufficient
for this two-angle problem.

Extensions (adaptive velocity mesh, finite orbit width, nonlinear resonances), their switches
and benchmarks: [drift-kinetics-extensions.md](drift-kinetics-extensions.md).

## Ambipolar transport on 3D equilibria (K1, K2)

Decision, 2026-11-03: extend this steady operator in two stages for the
ambipolar radial electric field E_r in W7-X. K1 is local three-dimensional
transport; K2 retains radial orbit coupling and determines the potential
jointly with the species distributions. Both use the collision ladder, D7
kinetic geometry view and equation-level derivatives below. This is planned
scope, not an implemented or validated solver. “Full orbits” here means
complete guiding-centre orbits within the drift-kinetic ordering, not
resolution of the gyromotion.

### K1: local transport and the ambipolar root

Activate option (c) as the equilibrium neoclassical response (n = 0 in the
perturbation labelling). The distribution still depends on both surface
angles: the local transport residual contains dot(theta) partial_theta f1
and dot(alpha) partial_alpha f1, with streaming, E x B and tangential
magnetic drift, plus the velocity-space terms of the selected trajectory
model and collisions. Radial magnetic drift acts on the Maxwellian drive;
its radial derivative on the unknown is dropped. Setting n = 0 in (4)
removes that equation's drift term and does not implement K1.

Before implementation, derive this local projection, its phase-space measure,
velocity-coordinate transport, null spaces and force/flux pairing. Distinguish
the incompressible monoenergetic ICNTS model from the full E x B/tangential-drift
model used for SFINCS comparisons. Every reduction has a measured model defect.

Extend the DG pitch/surface and nodal GLL energy discretisation, using FGMRES
with sparse direct blocks per energy as the starting preconditioner. A
three-dimensional field enlarges the surface block; its cost and robustness
are not established. The existing DK1 simpler-discretisation comparison and
stop rule still apply. Generic missing numerics belong in the owning library.

First gate D11, D31 and D33 against the International Collaboration on
Neoclassical Transport in Stellarators (ICNTS,
[Beidler et al. 2011](https://doi.org/10.1088/0029-5515/51/7/076001)) and matched
SFINCS/MONKES/DKES coefficients on published fields. Check the 1/nu and
sqrt(nu) limits, where nu is the deflection collision frequency. Intrinsic
ambipolarity in an axisymmetric limit requires momentum-conserving collisions;
Lorentz-only rows do not establish that gate.

Then compute the particle flux Gamma_a of each species with charge q_a and
solve sum_a q_a Gamma_a(E_r) = 0 by Newton with implicit tangents. Bracket and
continue roots; report folds and multiple solutions rather than infer
uniqueness from convergence. Check the root solver on an independently solved
reduced flux model and compare all bracketed kinetic roots with matched
SFINCS flux curves, with declared numerical and reference errors.

### K2: global steady distributions and quasineutrality

Retain v_d.grad(psi) partial_psi f on the unknown over a radial domain, where
v_d is the magnetic drift and psi labels nested flux surfaces. Use full-f or
a Maxwellian-subtracted representation of the same steady guiding-centre
transport/collision equation. Specify the particle and energy sources,
absorbing outer boundary (no incoming particles), axis regularity or inner
inflow data, and potential gauge before solving. Reservoir exchanges and
boundary losses enter the particle, charge and energy ledgers.

Determine Phi(psi) jointly with the distributions from quasineutrality by
Newton–Krylov iteration; E_r = -dPhi/dr for the declared radius r. The first
K2 derivation gate must supply the quasineutrality residual from the chosen
field action, any retained polarisation term, its charge/source compatibility
and null-space constraints. This scope decision does not assert existence or
uniqueness of a steady branch. Time-dependent E_r is outside this programme.

The local operator becomes a switch of the global one: remove radial
advection of the unknown while retaining the Maxwellian drive. Geometry,
collision and discretisation interfaces remain shared. Per-energy local
solves are the starting preconditioner, not a demonstrated global solution.
Differentiate the coupled equations for JVP/VJP, never the Newton iterations.
Use a consistent D7 psi across nested surfaces; refuse islands and non-nested
regions. Resolve numerical diffusion below the radial transport signal and
report separate numerical, model and input error; no certificate is promised
without the required stability and reconstruction bounds.

K2 promotion requires independent checks:

- a small-orbit-width scan recovering K1 and matched SFINCS/MONKES
  coefficients, with radial and velocity refinement;
- the closed momentum-conserving axisymmetric limit and analytic regime/root
  limits in their orderings; source or boundary currents are separate terms;
- particle and energy balance with all source, boundary and field-work terms;
- comparison with DK1's simpler finite-volume or Legendre reference on
  problems both representations solve;
- ion transport at fixed potential against public GORILLA. Its public
  monoenergetic applet can check D11 at fixed normalised E x B velocity;
  this is partial coverage. A Maxwellian flux comparison at arbitrary fixed
  potential remains a required open check until that reference is available.
  This programme does not modify GORILLA or its SCEF applet.

Use only published or synthetic cases in this public repository. Agreement
between deterministic and particle methods alone does not replace the
analytic and independent reference-code gates.

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
Status of (4)-(5): stated, V0 checks below. Implemented for one point of
the field line (`src/continuum/dk_pitch.f90`): elements bounded by lambda
levels carry continuous bases polynomial in t = sqrt(1/B - lambda) for each
sigma, joined at the reflection node, so both parities of P_l(xi) with
l <= p are exact to round-off and the null space is the constants; the first
mode outside the space converges like h^(2p) (test `dk_pitch_lorentz_operator`).
A basis polynomial in lambda instead reaches only first order for odd l,
because sqrt(1 - lambda B) is odd in t. The theta direction is next.

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

Energy: nodal GLL elements in x on [0, x_max] (pattern `src/pairs/radial_dvr.f90`),
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

Status: the Legendre reference is implemented for the axisymmetric n = 0
Lorentz problem (`src/continuum/dk_mono.f90`): Legendre polynomials in xi,
Fourier collocation in theta, block elimination in l and the density null
space bordered by a source multiplier. Against limits derived independently
by averaging the same equation (test `dk_mono_*`): collisionless D33 and D31
reproduce the passing and trapped fractions to 5e-6, banana D11/nu to 5e-6,
the Pfirsch-Schlueter D11/nu to 6e-10 and the plateau value to 1e-6;
D31 = -D13 and the border multiplier vanish to rounding. The
finite-collisionality correction of the banana limits scales like
sqrt(nu*) (-4.4e-2 at nu* = 1e-4). The cost of the low-collisionality end is
in theta, not in xi: at nu* = 2.5e-5 the coefficients converge to 1e-5 with
257 theta points and 800 Legendre polynomials, and theta resolution grows as
nu* falls, because the trapped-passing layer moves in theta. That
measurement is the case DG on lambda levels has to beat.

Against NEO-2-PAR in Lorentz mode on libneo's tok02 file (s = 0.50781,
nu* = nu R/(iota v) from 4e-4 to 4e1; NEO-2's collpar per cm is 2 nu_D/v):
NEO-2's conductivity ratio gamma_E = (pi/4) 3 nu D33/<b^2> and the
normalisation-free ratio D13 D31/(D11 D33) agree to 2e-4 and 1.1e-3, and the
gap shrinks as NEO-2 is refined (test `dk_mono_neo2_tok02`). The KIN6D values
are converged to 1e-6 at 65 theta points and 200 Legendre polynomials in
0.011 s of one core; NEO-2 approaches 1e-3 only at four times its default
resolution, in 1.9 to 4.9 s on two MPI ranks with 1 to 3.6 GB.

The solver is factorised once (`mono_factor`: LU of each Schur block S_l and
Y_l = S_l^-1 Lo_l, and the bordered l = 0 block) and every further
right-hand side reuses it (`mono_solve`). Derivatives of D_ij with respect
to nu, kappa and the field shape b(theta) follow the implicit rule
A dg = ds - dA g through the same factors, with the leaf tangents from the
FortSym JVP of `gen_dk_mono`; they agree with centred differences to the
differences' truncation error, 1e-10 to 4e-8 (test `dk_mono_jvp`). Each
coefficient comes with two indicators in the sense of
[certification.md](certification.md), not bounds: the discretisation
indicator compares with doubled theta and Legendre resolution (b resampled
by its trigonometric interpolant), effectivity 0.96 against a reference at
four times the resolution, and the algebraic indicator is the functional of
one solve with the true residual, at rounding level (test `dk_mono_estimate`).

Data layout for later batched or GPU execution: the blocks are dense,
contiguous n_theta x n_theta arrays indexed by l, and the factorisation is a
sequence of small dense LU and triangular solves per l. The batch dimension
is outside this routine: energy nodes, species and toroidal harmonics give
independent block systems of the same shape, solved as one batch. A port
waits for a profile of the full operator and must beat the CPU at matched
accuracy.

Orbit-class spectral elements (`src/continuum/dk_hp.f90`) solve the same
axisymmetric n = 0 equation on a mesh of the (theta, xi) cylinder whose rows
follow the conserved pitch lambda = (1 - xi^2)/b: passing rows are levels of
lambda, trapped rows are scaled separatrices zeta xi_c(theta), and the
separatrix lambda_c = 1/b_max itself is a mesh line. Columns start and end at
the X point theta_x of b_max, where the trapped elements collapse and their
nodes form one unknown. Rows are graded geometrically towards the separatrix,
columns towards theta_x. Elements carry tensor Lagrange polynomials of degree
p on Gauss-Lobatto nodes, continuous across elements; reflection points need
no treatment in (theta, xi). In the measure dtheta dxi/(b b.grad theta)
streaming is the bracket with the stream function -lambda/2; with its nodal
interpolant the bracket form is polynomial and integrated exactly, so it is
skew and constants are exact left and right null vectors. The system is
factorised by FortSparse's supernodal LU; the density null space is removed
by pinning one unknown and projecting, since a dense border row defeats the
fill-reducing ordering. A dual-weighted-residual indicator of a coefficient
weights the residual of the solution in degree p + 1 with the adjoint of that
degree (an indicator, not a bound) and drives bisection of rows and columns.

Checks (tests `dk_hp_*`): against `dk_mono` at nu = 1e-3 to 1e-1 the four
coefficients agree to 4e-9 with p = 6 and convergence is exponential in p;
D31 + D13 vanishes to 1e-10 of sqrt(D11 D33) from nu = 1e-6 to 1e2. From
nu = 1e-6 and 2.5e-7, outside the reach of the Legendre reference, Richardson
extrapolation reproduces the collisionless fractions f_c and 1 - f_c to
1.1e-6 and the banana D11/nu to 1.1e-6; the Pfirsch-Schlueter D11/nu at
nu = 1e3 agrees to 2e-8 (at nu = 1e4 rounding, amplified by the conditioning
of the collision-dominated system, limits it to 1e-6). On a coarse mesh at nu = 1e-4 the indicator has
effectivity 1.6, and one marking step reduces the D31 error twentyfold. After
refinement the effectivity fell to 0.05 in that check: the enrichment by one
degree misses error in layers that neither degree resolves.

Time to accuracy for D31 (eps = 0.1, eight cores of one AMD EPYC node per run;
reference: the mean of p = 7 and p = 8 element runs, which differ by at most
5e-8 relative): at nu = 1e-4 the Legendre reference is faster, 3e-10 in 0.5 s
against 8e-9 in 1.7 s. At nu = 1e-5 the elements give 1e-8 in 0.5 s, the
reference 4e-9 in 14 s and 4 GB. At nu = 1e-6 the elements give 2e-8 in 0.6 s
and 0.2 GB; the reference at 513 theta points and 3200 polynomials, 72 s and
13 GB, is 7e-5 off. At nu = 1e-7 the elements give 3e-8 in 1.2 s, while the
reference is 8e-4 off at its largest resolution. The adaptive loop, which
repeats the degree p + 1 solve at every step, is slower than an a priori
graded mesh of the same accuracy.

DK1 step 2 with the same reference discretisation (`src/continuum/dk_ntv.f90`,
#183): one toroidal harmonic n of the local Lorentz equation with the
streaming shift -i n q, the grad-B precession Omega = omega_E +
pre lambda (cos theta - offset) and the drive lambda S, whose bounce average
is the stage-1 bounce reference. (1 - xi^2) couples l to l +- 2, so pairs
(2j, 2j + 1) form blocks of size 2 n_theta, eliminated as for n = 0;
the nu, omega_E and amplitude sensitivities share one coefficient-space
torque adjoint through the same factors, so their cost no longer grows
with the number of amplitudes; the dense Schur products use BLAS. The
repeated spectrum-size timings behind that design change are retained
in `benchmarks/ntv`, a kernel comparison without a full-collision or
device claim. Derivatives in nu and omega_E agree with centred
differences to 4e-8 (test `dk_ntv_jvp`). So do the derivatives with respect
to every perturbation amplitude c_m and along (eps, k, offset, c_J), with the
operator and drive tangents from the FortSym JVP of the pointwise leaves
(test `dk_ntv_tangents`, 4e-8). On the stage-1 frozen circular
field (n = 3, full |B| spectrum), scored against the full-phase reference
(all trapped bounce harmonics plus passing), at nu_eff/omega_gradB = 0.006:

| eps | ratio | resolution spread | tolerance |
|---|---|---|---|
| 0.010 | 0.957 | 5e-4 | 0.10 |
| 0.020 | 0.961 | 1.2e-3 | 0.10 |
| 0.030 | 0.954 | 4e-4 | 0.10 |
| 0.040 | 0.962 | 8e-4 | 0.12 |
| 0.045 | 0.968 | 4e-4 | 0.13 |
| 0.060 | 0.982 | 3e-4 | 0.18 |

Every row passes, but the linear fit has intercept 0.948 +- 0.006, so
|a - 1| = 0.052 exceeds 0.05 and the SBP gate fails as scored. The spread is
the change over 193 to 321 theta points and 1600 to 6400 Legendre
polynomials. At eps = 0.03 the deficit falls with collisionality, 0.047,
0.042 and 0.037 at nu, nu/2 and nu/4, so finite collisionality explains part
of it; that it vanishes as nu -> 0 is not shown. Cost of a converged row:
308800 unknowns per energy node, 3.7 GiB and 5 min for 16 nodes on four
cores; resolution does not grow as eps falls from 0.06 to 0.01. In the 1/nu
regime with the drift off (eps = 0.03) the torque agrees with the reference
to 0.6 % and stays finite from omega_E = 0.012 rad/s down to omega_E = 0,
rising by 0.68 % like the reference (0.69 %), solved directly. DG in (theta, lambda) has not been needed for
resolution on this ladder.

Rescoring against NEO-2-QL (#183). NEO-2-QL (main 98ae645) was run on Boozer
files written from the same frozen field (n q = 6.5), for the same plasma and
collisionality, with pitch-angle scattering only: the test-particle operator
without energy scattering and field part (`isw_energy = isw_integral = 0`),
i.e. the Lorentz operator with the same nu_D(x). Its equation differs from
the stage-1 reference in three places, each a switch of `dk_ntv` (#238):
drift and drive weight (1 + xi^2)/b instead of lambda, the curvature drift
c_J xi^2 of the file's toroidal current (c_J = 2 B0 I'/(G q eps') = 0.045 to
0.27), and b.grad theta = k b with the flux measure of Boozer angles. With
all three on (KIN6D at 193 theta points and 1600 to 3200 Legendre
polynomials, NEO-2 at eta_part_trapped = 140, lag = 5):

| eps | KIN6D / Pa | NEO-2 / Pa | NEO-2/KIN6D | NEO-2 lag 5 -> 9 | NEO-2 eta 70 -> 140 |
|---|---|---|---|---|---|
| 0.010 | 2.834e-10 | 2.884e-10 | 1.018 | -0.2 % | -2.8 % |
| 0.020 | 2.446e-7 | 2.399e-7 | 0.981 | +0.5 % | -0.0 % |
| 0.030 | 1.178e-5 | 1.173e-5 | 0.996 | -0.1 % | +0.2 % |
| 0.040 | 1.691e-4 | 1.677e-4 | 0.992 | +0.4 % | +0.1 % |
| 0.045 | 4.354e-4 | 4.302e-4 | 0.988 | +0.7 % | +0.0 % |
| 0.060 | 2.828e-3 | 2.774e-3 | 0.981 | +0.5 % | +0.3 % |

The torque is the ion NTV torque density. NEO-2's lag 9 -> 13 change is
0.09 % (eps = 0.02); at eps = 0.01 its pitch grid is not converged and
eta_part_trapped = 200 diverges in its iteration. KIN6D changes by at most
1.3e-3 between 193 and 321 theta points and 1600 to 6400 polynomials, and by
+0.13 % (eps = 0.02) and +0.42 % (eps = 0.03) from 16 to 24 energy nodes
(24 -> 32: 0.04 %). With these drift and collision choices matched, the two
codes agree to 0.4 % to 1.9 %, NEO-2 below KIN6D at every eps >= 0.02;
at eps = 0.02 and
0.06 the gap exceeds both codes' numerical spread. NEO-2-QL also solves two source
terms of Kasilov et al. (2014), Eq. (61): the rest of the perturbed radial
drift beyond the toroidal-derivative drive, and the perturbed mirror and
collision operators acting on the axisymmetric solution f_10. With both
switched on in `dk_ntv` (`source_dpar`, `source_dlf10`, f_10 from `ntv_f10`,
same rows and resolution) the drift source changes the torque by -6.4 % to
-7.8 % and the f_10 source by +4.1 % to +6.7 %; together they cancel to

| eps | 0.010 | 0.020 | 0.030 | 0.040 | 0.045 | 0.060 |
|---|---|---|---|---|---|---|
| both sources | -3.67 % | -2.58 % | -0.96 % | -0.22 % | -0.09 % | -0.01 % |
| NEO-2/KIN6D | 1.056 | 1.007 | 1.005 | 0.994 | 0.989 | 0.981 |

They close the gap at eps = 0.02 and 0.03 and leave it unchanged for
eps >= 0.04; the 1.1 % and 1.9 % at eps = 0.045 and 0.06 are not explained
by them, and at eps = 0.01 NEO-2's pitch grid is not converged. The torques
are energy sums with the node weights wx (F1 + (x - 5/2) F2), which
reproduce the reported KIN6D column to 2e-3. The
switch defects on these rows are -0.3 % to +10 % (all three together, growing
with eps), so on some rows the stage-1 model and NEO-2's differ by more than
the codes do.
At eps = 0.03 both codes rise as nu falls: KIN6D by 0.10 % and 0.95 % at nu/2
and nu/4, NEO-2 by 0.22 % and 0.42 % (not pitch-converged at nu/4).

Interpretation of #183: an independent discretisation (NEO-2's finite elements
on eta levels along the field line) reproduces KIN6D's torque to 2 % with the
matched terms above. This does not establish the cause of the 4 to 5 % deficit
against the collisionless full-phase reference. The collisionality trend
supports finite-collisionality effects as a candidate explanation; that the
deficit vanishes as nu -> 0 is still not shown, and the source differences
above remain unresolved. With 24 energy nodes the stage-1
ratios rise by 0.1 to 0.4 % (0.962 and 0.957 at eps = 0.02 and 0.03). The
SBP gate stays failed as scored against the reference. Development proceeds
without treating this gate as passed.

Timing at the reported resolution errors, both codes on one AMD EPYC 7452
node with eight
threads each (KIN6D: the stage-1 model, 16 energy nodes solved one after the
other; NEO-2: test-particle pitch-angle operator). Error is the change to each
code's own finest run (KIN6D 193 x 3200 or 6400; NEO-2 eta 140, lag 9
where available); unknowns of KIN6D are summed over the 16 energy nodes:

| eps | code, resolution | error | unknowns | wall | peak memory |
|---|---|---|---|---|---|
| 0.045 | NEO-2, eta 70, lag 5 | 7e-3 (lag) | 0.57 M | 175 s | 18.2 GiB |
| 0.045 | NEO-2, eta 70, lag 9 | ~1e-3 | 0.95 M | 624 s | 60.7 GiB |
| 0.045 | KIN6D, 129 x 1600 | 1.0e-3 | 3.3 M | 170 s | 1.6 GiB |
| 0.045 | KIN6D, 193 x 1600 | 1e-5 | 4.9 M | 551 s | 3.6 GiB |
| 0.030 | NEO-2, eta 70, lag 5 | 2.3e-3 | 0.52 M | 166 s | 16.5 GiB |
| 0.030 | NEO-2, eta 140, lag 5 | 8e-4 (lag) | 1.32 M | 421 s | 44.4 GiB |
| 0.030 | KIN6D, 161 x 1200 | 2.0e-3 | 3.1 M | 240 s | 1.9 GiB |
| 0.030 | KIN6D, 193 x 1600 | 2e-5 | 4.9 M | 546 s | 3.6 GiB |

The errors in the table are those of the theta/Legendre and pitch/energy
resolution only. KIN6D's 16-node energy rule adds 0.13 % (eps = 0.02) and
0.42 % (eps = 0.03); 24 nodes reduce it and cost 1.5 times more. Comparing
KIN6D's 129 x 1600 row with NEO-2's eta 70, lag 9 row at eps = 0.045,
KIN6D is 3.7 times faster, uses 38 times less memory and has 3.5 times more
unknowns. Comparing KIN6D's 193 x 1600 row with NEO-2's eta 140, lag 5 row
at eps = 0.03, NEO-2 is 1.3 times faster; KIN6D uses 12 times less memory
and has 3.7 times more unknowns. Scaling the KIN6D timings by the stated
1.5 energy-node cost gives estimated speed factors of 2.4 in KIN6D's favour
and 2 in NEO-2's favour, respectively; these are not matched-total-accuracy
measurements. For the two coarser eps = 0.03 rows, NEO-2 is 1.4 times faster;
KIN6D uses 8.7 times less memory and has 6 times more unknowns. NEO-2 cannot
refine its pitch grid past
eta 140 at eps = 0.01, where it diverges at 200, so that row has no matched
accuracy. Not measured: the 16 energy nodes solved in parallel, which
KIN6D's layout allows, and other node types.

Energy grid and test-particle operator (`src/continuum/dk_collop.f90`, #184):
each species has the speed x = v/v_ta on continuous GLL elements of degree
at least 2 on [0, x_max], and the perturbation is f_a = F_a h_a P_l(xi) with
h_a in that space, so 1, x and x^2 are represented exactly. The operator is
stored per Legendre degree as the Galerkin matrix of int phi_i P_l C[f] d^3v
(rows: species and node; columns: partner species and node); pitch-angle
scattering, energy scattering and, for unequal temperatures, the thermal
exchange term come from FortSym leaves (`gen/gen_dk_collop.f90`, with a
series for the Chandrasekhar function below argument 1), integrated by
composite Gauss-Legendre rules converged to rounding. Streaming stays
diagonal in the nodes through the lumped GLL mass. Checks (tests
`dk_collop_*`): at equal temperatures the Maxwellian is a null vector, density
is conserved, the matrix is symmetric and its largest eigenvalue is at most
1e-14 of its norm, on one quadratic element as on 31 nodes; the deflection,
parallel-diffusion and flux rates agree with Rosenbluth potentials of the
Maxwellian computed by quadrature (`rosenbluth_reference`) to 6e-15, 3e-11 and
3e-11 (the oracle's own quadrature error); the momentum and heat taken from a
reservoir agree with the closed Maxwellian friction and temperature-exchange
rates to 2e-14 and 4e-13 for electrons and deuterons with T_i/T_e = 0.7.
Tangents with respect to temperatures and densities come from the leaf JVPs
and agree with centred differences to 4e-10. Krook (-nu_K f) is a parity
switch.

Field-particle part of self-collisions (same module, model `COLLOP_FULL`,
#185): since (grad grad G).v = grad(v G' - G), the flux of C(F_a, f_b) is F_a
times the gradient of a potential Psi P_l built from the Rosenbluth
potentials of f_b, and the operator is assembled in that weak form,
(g', Psi') + l(l+1) g Psi/v^2, with the potentials from the Legendre
expansions of 1/|v - v'| and |v - v'| (kernels from `gen_dk_collop`) by
quadrature split at the kink and graded towards it and towards the origin.
Every basis function is integrated over [0, oo) with the last element
continued polynomially, so that the invariants hold independently of
x_max; density is conserved identically. For electrons and for deuterons
alone (test `dk_collop_conservation`), density, momentum and energy are left
null vectors to 2e-14 of the cancelling terms on one quadratic element and to
1.3e-15 on 26 nodes, and the null spaces of l = 0 and l = 1 have dimensions
2 and 1. The self-collision operator is symmetric to 1e-16 and its largest
eigenvalue is below 1e-16 of its norm for l <= 6 (`dk_collop_self_adjoint`);
tangents in temperature and density agree with centred differences to
7e-10. Between species (`dk_collop_conservation_pair`, electrons and
deuterons on the same grids, T_i/T_e = 1 and 0.7): the density of each
species, total momentum and total energy are left null vectors to 1e-14,
the null spaces of l = 0 and l = 1 have dimensions 3 and 1, and at equal
temperatures the operator is symmetric to 5e-15 and dissipative. Friction
and linearised heat exchange between drifting or heated Maxwellians agree
with the closed rates to 3e-14 and 1e-12 (`dk_collop_exchange_pair`). The
temperature tangent of the field-particle part includes the boundary term of
the moving split point y = v/v_tb, where the derivative of the kernel jumps.
Against NEO-2's matrix elements for electrons and deuterons of the AUG 30835
benchmark run (`collop.h5`, T_i/T_e = 0.70, quadratic B-splines on knots
0:0.8:4, which lie exactly in cubic elements on the same knots; test
`dk_collop_neo2_matrix`, data in `test/data`): Lorentz and energy-scattering
elements agree to 2e-12, the field part to 2e-6 for self-collisions and to
8e-4 of the largest element for the electron-on-deuteron block at l = 0,
whose elements are 1e-4 of the Lorentz ones; that block is where the two
quadratures differ most, and which of them is closer is not established
beyond the KIN6D conservation and exchange-rate checks. Intrinsic
ambipolarity and the kinetic solve with energy scattering follow.

The primal field-particle assembly evaluates two scaled moment sweeps, below
and above each target speed, instead of repeating the entire source integral
at every target. The moments use the two separable powers of the same retained
Rosenbluth kernel; FortSym derives their coefficients and leaves for values,
Jacobian-vector products and vector-Jacobian products.
Quadrature and basis tables are shared between independent Legendre degrees,
which run in parallel with OpenMP. This removes repeated source integration;
the stored Galerkin matrices still have quadratic size in the energy nodes.
Derivative directions that leave the thermal-speed ratio unchanged reuse the
moment sweeps. Directions that change that ratio retain the direct kernel
integration and its moving-boundary term. A zero direction returns a zero tangent.

Test-particle quadrature sums accumulate in `real128` before conversion to the
`real64` operator. This reduces rounding error when the large self-collision
contributions cancel in the much smaller electron-ion heat exchange.
The conservation, exchange-rate and NEO-2 matrix checks remain diagnostic
oracles; this assembly does not yet supply a rigorous continuous error bound.

For a reproducible single-deuteron assembly measurement with 17 energy nodes
and degrees 0 through 128, run:

```sh
OMP_NUM_THREADS=1 build/test_dk_collop assembly_profile 128
```

Change the thread count to measure the parallel sweep. This reports assembly
time, not kinetic-solve time.

Collision ladder (one `apply`/`jvp`/`diag_block` interface per spatial
point; field parts as low-rank updates inside the Krylov matvec, never by
Richardson iteration):

| Increment | Operator | Role |
|---|---|---|
| DK1-DK2 | Lorentz (pitch-angle scattering) | numerics gate, regime ladder |
| DK3 | full linearised Landau, several species | NEO-2 parity target and production default; invariants by construction |
| DK3 | Sugama model operator ([Sugama, Watanabe, Nunami 2009](https://doi.org/10.1063/1.3257907), metadata not yet verified) | cheaper reduction with measured defect |
| parity | linearised relativistic operator (Braams–Karney and Beliaev–Budker levels, as NEO-2 `isw_relativistic`) | ECCD and hot electrons; nonrelativistic Landau is its limit |
| parity runs only | Krook (incl. PENTRC's harmonic rate) | MARS-K/PENTRC matching, never a default |

The energy space must contain $1$ and $x^2$ (GLL degree at least 2). NEO-2's
order-2 B-spline (hat) basis does not, and its discrete operator then conserves
energy only to $3.4\times10^{-2}$ at four basis functions (quadratic B-splines:
$5\times10^{-7}$); constraints imposed on such a basis hide the error instead
of removing it. Every operator is cross-checked against the richer model above
it (#188). Conservation defects are measured with the solver's moment functionals;
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

## Energy-coupled NTV kernel

`ntv_energy_coefficients` solves the NTV kernel with the linearised collision
operator of `dk_collop` instead of monoenergetic pitch-angle scattering. The
speed direction uses the element basis of `dk_collop` without its x = 0 node
(h = 0 there, exact for l >= 1 and of zero weight for l = 0). Streaming,
rotation, precession and pitch-angle scattering are lumped at the nodes, as
the mass of a spectral-element method, while energy scattering and the
field-particle part keep their Galerkin matrices, which couple the nodes and,
through the field-particle part, the species. f_10 comes from the same
operator at n = 0, with the theta average of each node's degree-0 part
constrained to zero, and enters through the mirror and collision sources of
the kernel. The coupled system is solved by GMRES, right-preconditioned with
one dense factorisation of all nodes and species for the two lowest Legendre
pair blocks and pair-block factorisations per node and species above them.
`ntv_profile` selects it with a collision model on its fourth input line.

## Axisymmetric Boozer oracle adapter

`kin6d_ntv_geometry` adapts an axisymmetric Boozer surface to the existing
local Lorentz kernel through `ntv_line_t`. The circular provider reproduces
the original kernel, including collision-frequency, electric-precession and
perturbation-amplitude tangents and resolution indicators. An imported line
supplies parallel streaming, mirror force, the two magnetic-drift weights,
drive and flux measure from the retained FortSym definition in
`gen/gen_ntv_boozer.f90`. The adapter retains the low-beta curvature and
omitted radial covariant-component assumptions of the quasilinear model.

This adapter supports comparisons with external codes; it does not replace
the D7 geometry interface or implement three-dimensional transport. Nonzero
toroidal geometry harmonics, absent coefficients and singular reference flux,
rotational transform or Jacobian are refused. Geometry-direction tangents of
an imported line are unavailable and return `NTV_BAD_INPUT`; the generated
pointwise JVP and VJP remain available. Synthetic analytic coefficients,
pointwise derivative checks, circular-provider equivalence and pinned leaf
regeneration verify this slice. Its resolution differences remain indicators,
not certified error bounds; no new device comparison or physical validation
is claimed.

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

## NEO-2 parity (tracker #198)

KIN6D anticipates NEO-2's complete feature set, as one operator instead of the
QL/PAR pair: monoenergetic and multispecies transport, bootstrap current,
conductivity and Ware pinch (#189), the generalized Spitzer function for current
drive with the relativistic operator (#190, #187), a fast-ion source (#191;
NEO-2's own is an unfinished stub), quasilinear NTV as a switch set (#192), the
radial electric field (#193), local 3D transport (#194), surface scans and the
exports downstream codes read (#196), and a time-to-accuracy comparison (#197).
Implementation starts after the generators move to the reworked FortSym API
(#179). The DK1 stop rule below still applies.

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
- K1 activates option (c) after its derivation and local coefficient/root
  gates above. K2 adds global steady radial transport only after K1 passes,
  with its own derivation, source/boundary and independent comparison gates.
  Islands and kinetic feedback remain outside this programme.

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
- The local ansatz still refuses large n q' Delta_psi_b and perturbative
  resonant layers. K2 admits finite radial orbit width through the global
  equation, not by using a local answer outside its ordering. Islands and
  non-nested regions remain outside this programme.
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
