# Nonlinear visco-resistive MHD: general 3D and symmetry-reduced 2D

Status: design decided 2026-11-02; implementation state 2026-11-09.
Implemented: RM0 (symmetry-reduced GS), RM1 (generated weak-form leaves, the
invariant section residual, Jacobian and mass in rotation, translation and
screw symmetry with manufactured, energy and spectral gates), the 2D
stationary Newton wrapper with Coulomb gauge and current or loop-voltage
drive, BDF2 evolution with supernodal direct factors, and, ahead of RM3, the
linear and nonlinear Fourier operators about an invariant state including the
closed helical coframe ([section contract](vrmhd-section.md)). Not yet
implemented: the RM-B5 island rates, the RM-B8 published-code and RM-B10 slow-manifold
gates, the 3D boundary solves and time evolution, the field-split
preconditioner (RM4), and the milestones RM9 to RM12 and gates RM-B11 to
RM-B13 of the model ladder, which are planned only. This document owns the
nonlinear visco-resistive part of [R12](https://github.com/itpplasma/milca/blob/94caec95e474bd23640278c5e03b50204943f369/ROADMAP.md#r12--equilibrium-mhdkinetic-mhd-vacuum-and-circuits)
item 12; Hall/XMHD stays with R12. Scope: tokamak and stellarator core MHD on
resistive time scales. The [model ladder](#model-ladder-on-the-resistive-time-scale)
(plan of 2026-11-13) orders the models from a surface-averaged background to
general 3D.

## Why

Nonlinear resistive MHD of the tokamak core on the resistive time scale is
expensive in established codes. The published 3D JOREK simulation of flux
pumping in ASDEX Upgrade (Zhang et al., Nucl. Fusion 65 066001, 2025) used
about 0.5 M core-hours on Marconi for 700 ms, about 6000 core-seconds per time
step.

Hypothesis: the cost is dominated by
preconditioning. The inter-harmonic coupling through parallel heat conduction
along the perturbed field grows with the perturbation amplitude and is
discarded by per-harmonic preconditioners. The design aims at preconditioners
whose iteration count does not grow with anisotropy, amplitude or time step,
tested by RM-B4.

## Target regime

The model, discretisation and solver choices are made for the parameters
between present medium-size tokamaks and a DEMO-class reactor (core values,
order-of-magnitude estimates).

| quantity | ASDEX Upgrade | ITER | DEMO |
|---|---|---|---|
| Lundquist number S | ~3e9 | ~1e10-1e11 | ~1e11 |
| tau_R/tau_A | >1e6 | >1e7 (Zhang et al. 2025: 10x AUG) | 40x AUG (Zhang et al. 2025) |
| chi_par/chi_perp (electrons), collisional | ~1e11 | ~1e13 | ~2e13 |
| chi_par/chi_perp (electrons), flux-limited | ~1e9-1e10 | ~1e10 | ~1e10 |
| magnetic Prandtl number | ~1e3 | ~1e2-1e3 | ~1e2-1e3 |
| relevant toroidal harmonics | n <= 2-4 (core 1/1 mode) | n <= 10-20 | n <= 10-20 or more |
| beta | low | low | moderate, high bootstrap fraction |

The collisional (Spitzer) parallel conductivity overstates the parallel heat
flux in these cores, where the electron mean free path (kilometres) exceeds the
connection length qR (10-60 m); the effective value is set by free streaming,
chi_par ~ v_te L_par. Resistive layers thin as S^(-1/3), so resonant surfaces
need local refinement rather than uniformly finer meshes.

Common structure: k_par << k_perp, nested surfaces with islands, slow
evolution on top of fast waves. Edge, X-point and high-n (ELM) applications are
not targets of this path; its architecture must not exclude them. Concretely:
the discretisation is not aligned to flux surfaces, so X-points and separatrices
are regular points of the mesh; wall, target and sheath conditions enter as
boundary terms of the weak form; the FCI parallel operator of the
preconditioner (U3) handles field lines that end on a boundary (intersection
with the boundary value, or volume penalisation as in GRILLIX, Stegmeir et al.,
Phys. Plasmas 26 052517, 2019), so the parallel block on open field lines
is a two-point problem along each line. RM-B3 and RM-B4 include an
open-field-line and a single-null X-point case to test this.

## Model

One residual `M(U) dU/dt + R(U; theta) = 0` (DESIGN.md), Cartesian physical
vectors, and no global toroidal angle in the residual or gauge:

- mass density rho with parallel and perpendicular diffusion and a source;
- velocity v with isotropic viscosity nu(T), optional parallel viscosity,
  Lorentz force J x B and pressure gradient;
- vector potential A and auxiliary current J in H(curl), B = curl A,
  with J defined by weak Ampère's law (J = curl B/mu0 for smooth fields),
  Ohm's law
  E = -v x B + eta(T)(J - J_s); the gauge condition is chart-free (the
  [NN2](equilibrium-nonnested-nn2.md) gauge, which is written through phi,
  is generalised before reuse);
- electron and ion temperatures with conduction tensor
  K = chi_par(T) b b + chi_perp (I - b b), Ohmic heating, optional
  equipartition and heat sources.

| symbol | meaning |
|---|---|
| U, theta | state vector, parameters |
| J_s | prescribed non-inductive current source |
| b | B/sqrt(B.B + epsilon_B²), with model parameter epsilon_B > 0 |
| eta, nu | resistivity (Ohm m, SI) and dynamic viscosity |
| tau_A, tau_R, tau_H | Alfven, resistive and poloidal Alfven times |
| S | tau_R/tau_A |
| F' | radial derivative of k.B at the resonant surface |
| Delta', W | tearing stability index, island width |
| N_fp | number of field periods |
| L_xi | Lie derivative along the Killing field xi |

This is the model class of JOREK model 711 (Pamela et al., Phys. Plasmas 27
102510, 2020; code overview Hoelzl et al., Nucl. Fusion 61 065001, 2021) and of
Zhang et al. 2025. Reduced MHD is a declared reduction of the same residual
(ordering, lift and model-defect estimate), not a second model.

Every pointwise leaf (weak-form integrands for each field, transport tensors,
eta(T), nu(T), chi(T)) has one FortSym definition. Primal, JVP and VJP come
from it. The Jacobian is applied, not assembled, except for the blocks a
preconditioner factorises.

The mixed-current volume integrands, units, source convention, mass terms
and energy-exchange identities are specified in
[the implemented weak forms](vrmhd-weak-form.md).

## Representations

### Symmetry-reduced 2D

A continuous isometry of Euclidean space generated by a Killing field xi
reduces the problem to the 2D orbit space for fields with L_xi U = 0:

| xi | reduction |
|---|---|
| rotation about a fixed axis | axisymmetric |
| translation | straight periodic cylinder or slab |
| screw motion (rotation plus translation, pitch h) | helical symmetry |

The reduced residual is the 3D residual restricted to invariant fields. FortSym
derives it from the same leaves by substituting the invariant ansatz, so there
is no second physics source. The symmetry is case input (a Killing field and an
orbit-space chart), not a coordinate built into the code, as
[geometry.md](geometry.md) requires.

First use: the axisymmetric Grad-Shafranov (GS) solver hard-codes the R dR dZ
measure (`src/equilibrium/grad_shafranov.f90`). It becomes the rotation
instance of a symmetry-reduced GS; its existing tests are the regression gate.
The new oracle is the helically symmetric GS equation (Johnson, Frieman,
Kulsrud and Oberman, Phys. Fluids 1 281, 1958).

Scope: exact for straight periodic cylinders and helices. A torus has no
helical isometry. Near-helical toroidal states use a few harmonics of the
symmetry-breaking coordinate about the reduced state; their model defect is
measured against the 3D representation, never assumed small. Stellarators have
no continuous isometry and use the 3D representation.

### General 3D

Phase 1 (correctness): the existing meridional FEEC complex times a real
Fourier series in a case-supplied symmetry angle, assembled on the CPU. It
covers axisymmetric devices with perturbations and the straight cylinder. The
Fourier chart is an optimisation for those geometries, as geometry.md permits;
the residual and gauge stay chart-free, and Phase 1 is not the stellarator
(Figure-8) path.

Phase 2 (speed): tensor-product isoparametric quads/hexes with
Gauss-Lobatto-Legendre (GLL) elements of degree 6-8, matrix-free
sum-factorised operators and p-multigrid. In the FortFEM discretisation
laboratory this was the fastest family at matched L1 error for a scalar
diffusion problem (FortFEM `doc/design/discretisation_laboratory.md`,
`benchmark/lab/RESULTS.md`); whether that carries over to MHD is tested by
RM-B9. Requires compatible GLL edge/face bases (U2). Stellarators enter here
through mapped tensor patches on the curve-following chart that MILCA already
has (`src/equilibrium/equilibrium_chart.f90`: X = Gamma(zeta) + q1 N(zeta) +
q2 B(zeta) about a closed, non-planar axis curve, as in GVEC's G-Frame). The
Figure-8 stellarator, which has no usable cylindrical toroidal angle, is the
acceptance case: its equilibrium is the R12 gate `g_frame_figure8_equilibrium`
(against GVEC on the same frame and boundary), and RM8 runs the nonlinear
model on it.

Discrete symmetries are used exactly. With N_fp field periods, the
linearisation about a symmetric state couples only harmonics n and n' with
n' = +-n mod N_fp; for real fields this gives floor(N_fp/2)+1 independent
preconditioner families.

Flux-coordinate-independent (FCI) parallel operators are used inside the
preconditioner. Whether they also improve accuracy as the parallel
discretisation is decided by RM-B3 in 3D (RM5); FCI is not the fundamental
geometry (DESIGN.md).

## Solvers

Stationary states, the main speed lever: Newton-Krylov with Eisenstat-Walker
forcing, pseudo-transient and pseudo-arclength continuation (FortFEM),
continuation in chi_par/chi_perp, S and drive. Driven steady states such as
flux pumping, saturated islands and the response to resonant perturbations are
computed directly instead of over 1e5-1e6 tau_A of evolution. Their stability
is checked by the [S0 consumer](https://github.com/itpplasma/milca/blob/94caec95e474bd23640278c5e03b50204943f369/ROADMAP.md#s0--stabilityresponse-consumer-cross-cutting-starts-with-r1r3).
Oscillating states need shooting or harmonic balance; that is a later step.

Time evolution: BDF2 first, then IMEX and
exponential integrators with Krylov phi-functions (Einkemmer, Tokman and
Loffeld, J. Comput. Phys. 330 550, 2017, studied them for resistive MHD), with
error-controlled steps. Transient adjoints follow
[differentiability.md](differentiability.md).

Preconditioning, the core of the design, splits the linearised operator into:

1. an anisotropic transport block: parallel solve along field lines on all
   harmonics at once, with the linearised FCI parallel operator (U3), plus
   perpendicular multigrid; this block carries the chi_par coupling between
   harmonics;
2. an MHD block: a physics-based Schur complement that parabolises the
   velocity equation (Chacon, Phys. Plasmas 15 056103, 2008), with
   per-harmonic or per-family blocks only as inner solves (U5);
3. an outer FGMRES on the full operator.

Both large parameters of the target regime are used, not only survived.

- Anisotropy, epsilon = chi_perp/chi_par ~ 1e-10: as epsilon -> 0 the
  temperature becomes constant along field lines (a flux function on nested
  surfaces, constant across an island, averaged over a stochastic region), and
  only the field-line-averaged perpendicular transport remains. The transport
  block uses this structure, either as an asymptotic-preserving micro-macro
  decomposition (field-line average plus an O(epsilon) correction; Degond,
  Deluzet, Narski and Negulescu and later work) or as a two-level
  preconditioner whose coarse space is the near-null space of b.grad, the
  field-line-constant functions. Conditioning and iteration count then do not
  depend on epsilon by construction (RM-B4).
- Time scales, tau_A/tau_R ~ 1e-6-1e-7: besides implicit stepping, the
  Schur-complement treatment of the stiff wave part and direct stationary
  solves, slow resistive evolution with a quasi-stationary helical or island
  state is computed on the slow manifold: a sequence of driven stationary 3D
  states parametrised by the slowly evolving profiles (adiabatic elimination of
  the fast dynamics), with a model-defect estimate against full time evolution
  and a local switch to it where the separation fails (mode onset, crash), as
  the model-switching contract of AGENTS.md requires (RM-B10).

Closures are separate FortSym leaves with one interface each, so a closure can
be replaced without touching the residual: eta(T), nu(T), chi_par(T) in
collisional (Spitzer) form, a flux-limited chi_par (parallel heat flux bounded
by a multiple of n T v_te) for the nearly collisionless core, and the
non-inductive current J_s as a prescribed source. Comparisons with published
results of other codes use the coefficients stated in those publications;
physics runs state the closure they use.

Closure requirements, from the first closure on:

- Asymptotic consistency: every closure reproduces the classical limits in
  their regimes, each a registered gate with the classical formula as oracle:
  Spitzer and Braginskii parallel and perpendicular transport for
  lambda_mfp << L_par and Omega tau >> 1; free streaming and the
  Hammett-Perkins Landau-damping response (Phys. Rev. Lett. 64 3019, 1990)
  for lambda_mfp >> L_par; ideal and kinetic MHD for rho_i/L -> 0;
  gyroviscous (FLR) MHD at first order in k_perp rho_i (Roberts and Taylor,
  Phys. Rev. Lett. 8 197, 1962). A blended closure that is good in between
  must still meet each limit to its stated order.
- Validity monitors: every run computes, locally in space and time, rho_i and
  d_i over the smallest resolved scale (layer, island, pedestal width), banana
  and fast-ion orbit widths over that scale, the Knudsen number
  lambda_mfp/L_par, the Landau parameter omega/(k_par v_th), omega/Omega_i,
  beta and S, and reports them with the result. A result computed where the
  closure in use is outside its validity range is flagged, never reported
  silently; the model-defect estimate drives promotion to a richer closure.

Physics-run rule: a transport coefficient that differs from its physical value
(an enhanced resistivity or viscosity, for example) is stated with its reason
and with a scan in S that shows its effect. Physical (neoclassical) resistivity
for current diffusion is kept separate from any declared layer model
(two-fluid terms, hyper-resistivity, electron viscosity), and physics results
are quoted only where the scaling regime (Sweet-Parker or plasmoid, layer
width against rho_i) is identified. The sensitivity of a result to layer
physics is itself reported.

Every reduction carries a model-defect estimate and a switch to the fuller
model where the estimate exceeds its budget: symmetry-reduced 2D against 3D,
near-helical toroidal states, reduced MHD, flux-limited closure and the slow
manifold. Scale separation fails locally (resistive layers, crashes,
stochastic regions), so the switch is local in space and time.

Research option, not on the critical path: matched asymptotics, an outer ideal
solution coupled through Delta' to an asymptotic resistive layer, as coarse
model or control variate for the layer in linear and weakly nonlinear
response.

Research option, not on the critical path: algebraic decoupling of the
frozen-coefficient linearised operator over Weyl-algebra coefficients
(Stafford, J. London Math. Soc. 18 429, 1978; algorithms of Quadrat and
Robertz, J. Symbolic Comput. 42 1113, 2007), generated with FortSym, to derive
block-triangular preconditioners. Admitted only through the same gates.

## Model ladder on the resistive time scale

Plan of 2026-11-13. This section adds no implementation; the state column
says what runs today.

Two questions drive this path: how the safety-factor profile q of a driven
tokamak core evolves on the resistive time scale, and which stationary states
it reaches with and without a saturated m/n = 1/1 mode. They are answered on a
ladder of models that are restrictions or declared reductions of the one
residual above. A cheap rung isolates a mechanism and can be converged and
scanned; a fuller rung confirms it. Wherever it is meaningful, a rung offers
three solves of the same discrete residual: time evolution, a direct Newton
solve for the stationary state, and the spectrum of the linearisation about
that state (the analysis is owned by [stability.md](stability.md)).

| rung | model | dropped relative to the 3D residual | state | milestone |
|---|---|---|---|---|
| L0 | surface-averaged transport of q, pressure and density with a Grad-Shafranov constraint (1.5D) | inertia, viscosity, variation of density and temperature on a flux surface, all n > 0 | designed | RM10 |
| L1 | helically symmetric straight cylinder (2D), [vrmhd-cylinder.md](vrmhd-cylinder.md) | toroidal curvature; every helicity except multiples of m/n = 1/1 | stationary solve, spectrum and BDF2 run; trajectories are candidates | RM2 |
| L2 | axisymmetric torus (2D), [vrmhd-torus.md](vrmhd-torus.md) | all n > 0 | Grad-Shafranov initial state registered; the 2D stationary and BDF2 solves exist for the rotation symmetry but are not yet applied to this case | RM9 |
| L3 | torus with a few (m, n) harmonics on a radial mesh (proposal) | all other harmonics | proposal, to be tested against L4 | RM11 |
| L4 | torus with toroidal harmonics n = 0, 1, [vrmhd-torus.md](vrmhd-torus.md) | n >= 2, including the 2/2 harmonic of a saturated 1/1 mode | volume operators run; boundary solves and evolution open | RM5 |
| L5 | more toroidal harmonics; general 3D ([Representations](#representations)) | nothing | designed | RM5, RM6, RM8 |
| LM | L0 with a mean-field closure for the mode (hypothesis) | the mode itself | hypothesis | RM12 |

Symbols: m and n poloidal and toroidal mode numbers, `<.>` flux-surface
average, `Phi` toroidal flux inside a surface, `Vol' = dVol/dPhi` with `Vol`
the enclosed volume, `iota = 1/q`, `phi` the toroidal angle, `F = R B_phi`,
`V` the applied loop voltage and `V_L(Phi)` the loop voltage of one surface,
`L_z` the axial period of the cylinder, `Gamma` particle flux, `rho` mass
density with source `Srho`, `M = a v_A rho/nu` the viscous Lundquist number
(`S/Pm` in Alfven units) and `H = sqrt(M S) = a B/sqrt(eta nu)` the Hartmann
number.

### Checks common to all rungs

Tolerances are fixed here, before any measurement. The default is the
observable radius of AGENTS.md, 1 % of a declared physical scale.

1. Reduction against parent. In the limit both models share, named
   observables (q at the axis and at fixed flux labels, `V_L` per surface,
   mode energies, growth rates) agree within the larger of 1 % and the sum of
   the two discretisation estimates, each taken from two resolutions. Where
   the reduction has a small parameter (`a/R0`; `1/S`, `chi_perp/chi_par` and
   `1/H^2` for the inertia, surface variation and viscosity that L0 drops),
   the difference is also measured at three values of it at least,
   and its order is reported; a difference that does not shrink is a model
   defect, not a tolerance to widen.
2. Newton state against evolution. Both solve the same discrete residual, so
   they agree to solver tolerance, not to discretisation error: an evolution
   started from the Newton state with a perturbation of relative size 1e-3
   returns to it within ten times the Newton tolerance if the spectrum is
   stable, and its slowest decay rate equals the least damped eigenvalue
   within 1 % once the time step is halved without change at that level
   (BDF2 is second order). A helical state is compared up to a shift of its
   phase, and the neutral eigenvalue of that shift is excluded.
3. Eigenvalue against growth. For an unstable state the growth rate of the
   perturbation in the evolution equals the leading eigenvalue under the same
   1 % rule.
4. In the checks of this ladder, published results of closed codes are
   context for a matched-model comparison, not a pass criterion
   ([literature-policy.md](literature-policy.md)), and comparisons use
   published results and published input data only (AGENTS.md).

### L0: surface-averaged background

Model (Jardin 2010, Ch. 6, in the single-fluid form of his Sec. 6.4, with the
Ohm's law `E + v x B = eta (J - Js)` of this residual):

```text
V_L   = 2 pi eta <(J - Js).B> / <B.grad phi>          (Jardin 6.38, 6.52, 6.64)
d iota/dt = d V_L/d Phi   at fixed Phi                (6.48)
d(rho Vol')/dt + d(Vol' Gamma)/d Phi = Vol' <Srho>    (6.45)
Gamma = 2 pi q rho [<R^2 grad phi.Rv> - <B.Rv>/<B.grad phi>],  Rv = eta (J - Js)   (6.40)
```

with the energy equation (6.46) closed by the perpendicular conduction and
sources of this residual, density diffusion added to `Gamma` as in Jardin's
diagonal model (Sec. 6.2.4 C), and the geometry from a Grad-Shafranov solve at
fixed q(Phi) and entropy (6.58; Grad and Hogan 1970). On this time scale only
the surface-averaged parallel Ohm's law moves q; the perpendicular part fixes
the velocity relative to the toroidal-flux surfaces (6.34) and hence `Gamma`,
which for a scalar resistivity is Pfirsch-Schlueter and classical diffusion
plus an inward pinch (6.68, 6.69). A stationary q profile needs `V_L` uniform
(6.48). `Js` occupies the current-drive slot (6.52). Force balance replaces
the momentum equation with an error of second order in `1/S` (6.10, 6.11).
The two conditions per surface that close the magnetostatic equations for a
slowly diffusing plasma go back to Kruskal and Kulsrud (1958); Hinton and
Hazeltine (1976, Sec. II.D) give the two-fluid form of the fluxes, which this
single-fluid residual does not contain.

- Derivation: FortSym forms the integrands of `V_L` and `Gamma` from the Ohm
  and continuity leaves of the rotation instance and proves the pointwise
  identity behind (6.34); surface averages use the existing equilibrium
  geometry ([equilibrium-core.md](equilibrium-core.md)), the constraint the
  prescribed-q solve ([gs-inverse.md](gs-inverse.md)). L0 is the first 1.5D
  instance of [transport-multiscale.md](transport-multiscale.md) and follows
  its closure contract.
- Internal checks (RM-B12): q(Phi, t) over one resistive time against L2 from
  the same initial state and closures; the stationary states of L0 and L2
  (profiles and the flow field of 6.34).
- External: in a cylinder with uniform resistivity and a fixed wall electric
  field, (6.66) gives an axial-current perturbation `J0(j01 r/a)` that decays
  at the rate `(eta/mu0) (j01/a)^2`, `j01 = 2.405` the first zero of the
  Bessel function J0 (our evaluation of 6.66); the ordering of the skin,
  particle and heating times of Jardin Sec. 6.4.
- Acceptance: the decay rate to 1e-4 and at design order, since the mode is
  exact and the problem one-dimensional; rule 1 against L2.
- Open: the prescribed-q solve takes q against poloidal flux and a pressure
  law, not the adiabatic variables of (6.58); the free-boundary circuit
  equations (6.56) are not planned here.

### L1: helically symmetric cylinder

The full residual on fields that depend on `r` and `theta - z/R0`. It is the
cheap mechanism model: saturation of the 1/1 mode, its electromotive force,
`eta(T)` and anisotropic conduction, at the cost of a 2D problem. It cannot
contain the Shafranov shift, Pfirsch-Schlueter flows, the coupling of m = 1
to m = 0 and 2 or the pressure drive of the internal kink, which are toroidal
([vrmhd-torus.md](vrmhd-torus.md)).

- Derivation: screw-symmetric ansatz substituted into the same leaves (RM1).
- Internal checks: RM-B7 against the Fourier representation; rules 2 and 3 on
  the stationary helical state (today a probe, not a gate); the loop-voltage
  ledger RM-B13: at every radius of a stationary helical state the average
  over the helical angle obeys
  `V/L_z = <eta> <J_z> + <eta~ J_z~> - <(v x B)_z>`, `~` marking the deviation
  from the average and `J` standing for `J - Js` when a source is present.
  The last term holds the electromotive force of the mode and the part
  `<v_r> <B_theta>` of the mean radial flow, the second the correlation of
  resistivity and current. The ledger is exact, because the electrostatic
  part of `E_z` averages to zero.
- External: linear rates of RM-B2 (Coppi et al. 1976; Ara et al. 1978) and
  island rates of RM-B5.
- Acceptance: ledger residual below 1e-3 of `V/L_z` and falling at design
  order, because a residual `dV` changes q at the relative rate `dV/V` per
  resistive time (6.48); rules 2 and 3.
- Open: no published single-helicity stationary state with matched inputs is
  recorded here as an external reference. The case holds its poloidal current
  with a source `Js = J_theta e_theta`, the cylinder form of the switch
  discussed for L2; without it the stationary cylinder has the classical
  diffusion and pinch flow of (6.69) in its cylinder limit and needs the same
  particle closure. Whether L1 moves to the source-free state is an open
  decision, and until then a comparison of L4 with L1 has to state which
  convention each side uses.

### L2: axisymmetric torus

The rotation instance of the section residual. Target: the stationary state
driven by the loop voltage alone, with uniform `V_L`, Pfirsch-Schlueter and
pinch flows and a stated particle closure. With a resistivity that is
constant on flux surfaces it has no static counterpart: the argument is on
the torus page and is consistent with Montgomery, Bates and Lewis (1997) and
Bates and Montgomery (1998). The case and the current-source switch are in
[vrmhd-torus.md](vrmhd-torus.md#steady-ohmic-state).

- Derivation: rotation ansatz substituted into the same leaves (RM1); the
  lift of the Grad-Shafranov state into the section.
- Internal checks (RM-B11): `V_L` per surface against the applied voltage;
  stationary profiles and flow against L0; rules 2 and 3.
- External: the Shafranov shift (registered, `vrmhd_torus_shafranov`); the
  large-aspect-ratio flux of Wesson Eq. 4.5.5, which is Jardin's (6.69) in
  that limit; the explicit low-Hartmann-number flow of Bates and Montgomery
  (1998) in a circular torus with no-slip wall and uniform coefficients.
  Context: low-Reynolds-number flows in a rectangular section with
  stress-free wall (Montgomery, Bates and Li 1997; Kamp, Montgomery and Bates
  1998), numerical steady states at Hartmann numbers from far below to far
  above one in the same geometry (Kamp and Montgomery 2003) and
  two-fluid M3D-C1 steady states whose radial flows follow the
  Pfirsch-Schlueter form (Ferraro and Jardin 2009, Eq. 39).
- Acceptance: `V_L` uniform within 1e-3 of the applied voltage and improving
  at design order (the identity is exact, the reason for 1e-3 as in L1);
  flow against L0 by rule 1, outside a wall layer of a few `a/H`; the
  Bates-Montgomery flow within 1 % at `H^2 <= 1e-3` and a squared viscous
  Lundquist number `<= 1e-3`, their expansion parameters, with the difference
  shown to fall linearly in both.
- Open: the stress-free wall of the three rectangular-section papers is not
  in the weak form; the compressible residual reaches the incompressible
  limit of these papers only at high background pressure, which the check has
  to demonstrate. The Bates-Montgomery check tests the operator at `H << 1`,
  while the case has `H` near `3e5` (`S = 1e6`, `Pm = 10`); in the target
  regime the flow is checked only against L0.

### L3: few harmonics on a radial mesh (proposal)

Proposal, not a result: a radial mesh times a small set of (m, n) harmonics,
for example 0/0, 1/0, 1/1 and 2/1, with the exact toroidal metric,
products evaluated at sampled angles as in the nonlinear Fourier operator, and
the same implicit solvers. It would give toroidal coupling at the cost of a
few 1D problems. The metric is not expanded in `a/R0`: Bussac's energy is of
fourth order in `r1/R0`, two orders beyond the order at which the cylinder
mode is marginal
([vrmhd-torus.md](vrmhd-torus.md#bussac-limit-of-the-ideal-internal-kink)),
so second-order terms of the metric are of the size sought.

The harmonics of the symmetry-breaking coordinate about a helical state,
foreseen under [Symmetry-reduced 2D](#symmetry-reduced-2d), are a different
truncation: a mode (m, n) has harmonic number m - n there, so numbers 0 and
+-1 hold all m = n and m = n +- 1. On the section mesh this costs three
blocks, as L4 does, so it is no cheaper than L4; it differs in what it keeps.

- Derivation: FortSym substitutes the double-Fourier ansatz into the same
  leaves; axis regularity of each m is part of the derivation.
- Checks: rule 1 against L4 on the internal-kink growth rate across the
  Bussac threshold, the saturated mode energy and the stationary q profile,
  at `R0/a = 3` and 10, with the harmonic set enlarged until the difference
  stops changing.
- Acceptance: admitted only if it meets rule 1 against L4 and costs at most a
  tenth of L4 at that accuracy; otherwise it is recorded as failed and
  dropped. The factor is a choice: below it the extra model does not pay for
  its verification.
- Open: everything. A radial-Fourier discretisation of the nonlinear residual
  needs admission under [numerical-methods.md](numerical-methods.md).

### L4: torus with n = 0, 1

The case of [vrmhd-torus.md](vrmhd-torus.md): toroidal coupling of m = 1 to
m = 0 and 2, the pressure drive of the internal kink (Bussac et al. 1975) and
a stationary helical core in a torus.

- Derivation: Fourier factor in the rotation angle about the L2 state, orbit
  leaves from `gen_vrmhd_section` (RM3).
- Internal checks: n = 0 block against L2; rules 2 and 3 on the stationary
  helical state; the ledger RM-B13 for the n = 0 toroidal Ohm's law, the
  toroidal analogue of the cylinder form (Krebs et al. 2017, Eq. 4); L1 with
  axial period `2 pi R0`, as `R0/a` grows at fixed q profile, `S` and
  closures, with the order of the difference measured
  and not assumed: the ideal energy of the cylinder mode vanishes at second
  order in `r1/R0` and Bussac's toroidal energy is of fourth order, so the two
  models may differ at leading order.
- External: the Bussac threshold and ideal rate through MILCA's ideal
  stability, and the geodesic acoustic frequency, as on the torus page.
  Context: M3D-C1 single-fluid runs with published equations and input
  tables,
  where sawteeth give way to a stationary 1/1 state above a threshold in
  beta (Krebs et al. 2017, appendix Tables 1 and 2; Jardin, Ferraro and Krebs
  2015, not read for this plan), and JOREK runs of ASDEX Upgrade with states
  ordered by Hartmann number (Zhang et al. 2025, 2026), whose device inputs
  are not public.
- Acceptance: rule 1 for the reductions; the proposed checks of the torus
  page; the model defect of n <= 1 is the change of the observables when
  n = 2 is added (L5), reported with every result.
- Open: boundary solves and evolution of the Fourier representation (RM5);
  whether the Krebs case can be matched from the published tables alone.

### L5: more harmonics and general 3D

RM5, RM6 and RM8 as planned. On this ladder its task is the convergence of
L4: observables under added toroidal harmonics, by rule 1.

### LM: mean-field closure (hypothesis)

Hypothesis: on the resistive time scale a saturated 1/1 mode acts on the
background only through an extra loop voltage in L0,

```text
V_L = 2 pi eta <(J - Js).B>/<B.grad phi> + V_mode(Phi; profiles),
```

where `V_mode` collects the electromotive and resistivity-correlation terms
of the ledger RM-B13 and is a function of the axisymmetric profiles. The
split of the n = 0 Ohm's law into these terms is published (Krebs et al.
2017, Eqs. 4 to 6); that they are determined by the profiles alone is not
established. The distinct states found at different dissipation, and the
transition from marginal flux pumping to sawteeth in long runs (Zhang et al.
2026, abstract), indicate that a closure for stationary states cannot cover
every regime. The closure is a provider under the
[closure contract](transport-multiscale.md#closure-provider-contract) and
enters only after an independent check: calibrated on stationary states of
L1 and L4, it must predict q on states not used for calibration (other `S`,
`beta_p` and heat-source widths) and reproduce the slow-manifold sequence of
RM-B10 within rule 1, and it reports where it is outside its calibrated
range. Until then no result is quoted from it.

### Sources of the ladder

- S. C. Jardin, *Computational Methods in Plasma Physics*, CRC Press (2010),
  Ch. 6, [doi:10.1201/EBK1439810958](https://doi.org/10.1201/EBK1439810958).
- H. Grad and J. Hogan, Phys. Rev. Lett. 24, 1337 (1970),
  [doi:10.1103/PhysRevLett.24.1337](https://doi.org/10.1103/PhysRevLett.24.1337).
- M. D. Kruskal and R. M. Kulsrud, Phys. Fluids 1, 265 (1958),
  [doi:10.1063/1.1705884](https://doi.org/10.1063/1.1705884): conditions for
  a plasma diffusing slowly across toroidal surfaces (Secs. J to L).
- F. L. Hinton and R. D. Hazeltine, Rev. Mod. Phys. 48, 239 (1976),
  [doi:10.1103/RevModPhys.48.239](https://doi.org/10.1103/RevModPhys.48.239),
  Sec. II.D: Pfirsch-Schlueter flux in the two-fluid description.
- D. Montgomery, J. W. Bates and H. R. Lewis, Phys. Plasmas 4, 1080 (1997),
  [doi:10.1063/1.872195](https://doi.org/10.1063/1.872195).
- D. Montgomery, J. W. Bates and S. Li, Phys. Fluids 9, 1188 (1997),
  [doi:10.1063/1.869206](https://doi.org/10.1063/1.869206).
- J. W. Bates and D. C. Montgomery, Phys. Plasmas 5, 2649 (1998),
  [doi:10.1063/1.872952](https://doi.org/10.1063/1.872952).
- L. P. Kamp, D. C. Montgomery and J. W. Bates, Phys. Fluids 10, 1757 (1998),
  [doi:10.1063/1.869692](https://doi.org/10.1063/1.869692).
- L. P. Kamp and D. C. Montgomery, Phys. Plasmas 10, 157 (2003),
  [doi:10.1063/1.1524629](https://doi.org/10.1063/1.1524629).
- N. M. Ferraro and S. C. Jardin, J. Comput. Phys. 228, 7742 (2009),
  [doi:10.1016/j.jcp.2009.07.015](https://doi.org/10.1016/j.jcp.2009.07.015).
- M. N. Bussac, R. Pellat, D. Edery and J. L. Soule, Phys. Rev. Lett. 35,
  1638 (1975),
  [doi:10.1103/PhysRevLett.35.1638](https://doi.org/10.1103/PhysRevLett.35.1638).
- A. A. Ware, Phys. Rev. Lett. 25, 15 (1970),
  [doi:10.1103/PhysRevLett.25.15](https://doi.org/10.1103/PhysRevLett.25.15):
  the pinch of trapped particles. It is kinetic and outside this ladder; the
  residual contains only the fluid pinch `E_phi B_theta/B^2` (his Eq. 1).
- I. Krebs et al., Phys. Plasmas 24, 102511 (2017),
  [doi:10.1063/1.4990704](https://doi.org/10.1063/1.4990704).
- S. C. Jardin, N. Ferraro and I. Krebs, Phys. Rev. Lett. 115, 215001 (2015),
  [doi:10.1103/PhysRevLett.115.215001](https://doi.org/10.1103/PhysRevLett.115.215001)
  (not read for this plan).
- H. Zhang et al., Nucl. Fusion 65, 066001 (2025),
  [doi:10.1088/1741-4326/adcfb4](https://doi.org/10.1088/1741-4326/adcfb4);
  Nucl. Fusion 66, 076027 (2026),
  [doi:10.1088/1741-4326/ae7613](https://doi.org/10.1088/1741-4326/ae7613).
- J. Wesson, *Tokamaks*, Secs. 4.4 and 4.5, as cited on the torus page.

No open-source repository was used for this plan. The Hinton-Hazeltine and
Wesson entries are secondary; the formulas above were checked in Jardin.

## Benchmarks and gates

Each gate is a registered CTest with an analytic or manufactured oracle first
and an external code second.

The test-side [linear cylinder and island references](linear-mhd-oracles.md)
state the dimensional conventions, analytic checks and validity limits for
RM-B2 and RM-B5. They are independent references, not completed solver gates.

| ID | gate | oracle |
|---|---|---|
| RM-B1 | manufactured solution of the full residual, all fields | design order of convergence |
| RM-B2 | linear waves and instabilities | Alfven and magnetosonic dispersion in a uniform slab; resistive internal kink gamma ~ (eta F'^2)^(1/3) ~ tau_H^(-2/3) tau_R^(-1/3) (Coppi et al., Sov. J. Plasma Phys. 2 533, 1976); Furth-Killeen-Rosenbluth tearing (Phys. Fluids 6 459, 1963); existing `response_strip_cylinder_layer` values |
| RM-B3 | anisotropic transport | Sovinec chi-ratio ladder (Sovinec et al., J. Comput. Phys. 195 355, 2004; tests `nonnested_transport_sovinec_chi_ratio_p*`) extended to 1e13; closed-field-line benchmark of Guenter et al., J. Comput. Phys. 209 354, 2005; anisotropic diffusion with field lines ending on plates and in a single-null X-point geometry (manufactured solutions); overlapping islands with an ergodic layer, where the field-line average runs over the stochastic region (Rechester and Rosenbluth, Phys. Rev. Lett. 40 38, 1978) |
| RM-B4 | preconditioner independence | FGMRES iterations bounded independent of chi_par/chi_perp (1e6-1e13), Lundquist number (up to 1e12), perturbation amplitude and dt, on closed and open field lines; linear amplitudes in RM4, saturated states in RM5 |
| RM-B5 | nonlinear islands | Rutherford growth dW/dt = 1.22 (eta/mu0) Delta'(W) (Rutherford, Phys. Fluids 16 1903, 1973; coefficient from White et al., Phys. Fluids 20 800, 1977); saturated width (Escande and Ottaviani, Phys. Lett. A 323 278, 2004; Militello and Porcelli, Phys. Plasmas 11 L13, 2004); plasmoid onset of a Sweet-Parker sheet above S ~ 1e4 and the S-independent reconnection rate ~ 0.01 beyond it (Loureiro, Schekochihin and Cowley, Phys. Plasmas 14 100703, 2007; Bhattacharjee, Huang, Yang and Rogers, Phys. Plasmas 16 112102, 2009), which also fixes the resolution needed in crashes |
| RM-B6 | ledgers | energy and mass; helicity in the ideal limit |
| RM-B7 | 2D against 3D | identical equations on a straight periodic cylinder agree to discretisation error; the 3D solver preserves helical symmetry |
| RM-B8 | published codes | published nonlinear results of other codes for tearing and internal kink in circular geometry and for the flux-pumping states of the straight cylinder and the circular torus with its toroidal limits ([vrmhd-torus.md](vrmhd-torus.md)), with the model, inputs and resolution stated in the publication; 2D part in RM2, 3D part in RM5 |
| RM-B9 | performance | time to matched observable accuracy against published costs of comparable simulations (Zhang et al. 2025), strong scaling, memory, CPU and GPU |
| RM-B10 | slow manifold | quasi-stationary sequence of driven stationary states against full BDF2 evolution on the straight periodic cylinder: q profile and island width to the stated tolerance, defect estimate bounds the observed difference, switch to time evolution at mode onset |
| RM-B11 | stationary Ohmic torus without current source ([ladder](#l2-axisymmetric-torus)) | loop voltage uniform over the flux surfaces; profiles and flow of the surface-averaged model L0; large-aspect-ratio particle flux (Wesson Eq. 4.5.5; Jardin 2010, Eq. 6.69); low-Hartmann-number flow of Bates and Montgomery, Phys. Plasmas 5 2649, 1998 |
| RM-B12 | surface-averaged background L0 ([ladder](#l0-surface-averaged-background)) | decay of an axial-current eigenmode in a uniform cylinder (Jardin 2010, Eq. 6.66); q(Phi, t) and stationary state of the axisymmetric torus L2 |
| RM-B13 | loop-voltage ledger | averaged Ohm's law of a stationary helical state, at each radius of the cylinder and for the n = 0 toroidal component in the torus: applied voltage = resistive + resistivity-correlation + electromotive part |

Admission (AGENTS.md):

- Certification: the symmetry-reduced GS and stationary residuals target
  certified-discrete through the existing GS certificate path. Transient and
  nonlinear outputs are `candidate` until a stability estimate for the
  linearisation is available.
- Differentiability: stationary states carry equation-level JVP/VJP from the
  first stationary solve; transient adjoints use checkpointing.
- Memory: Krylov and preconditioner memory are budgeted per gate and reported
  with RM-B9.

Speed targets (to be measured, not claims): a 2D steady state within a minute
on one core; a 3D steady state within an hour on one GPU node; 3D time
evolution at least 10x faster than the published cost of comparable JOREK
simulations (Zhang et al. 2025) at matched accuracy. Physics
applications start only after RM-B1 to RM-B9 pass on the representation used.

Non-public device data (experimental profiles, partner restart files) never
enter this repository. Device comparisons beyond published tables live in the
owning private projects.

## Upstream work

Generic pieces are added to the owning libraries on their main branches and
pinned here. FortGen is FortSym's kernel intermediate representation and
emitter; `emit_tensor_fortran` and `emit_tensor_cuda` are its sum-factorised
tensor-product schedules.

| ID | library | work |
|---|---|---|
| U1 | FortSym, FortGen | multi-field, vector-valued tensor-product leaf with Piola channels and shared transforms; one array packing for CPU and CUDA |
| U2 | FortFEM | compatible GLL edge/face bases (H(curl), H(div)) on curved quads/hexes; generated matrix-free apply through `emit_tensor_*` |
| U3 | FortFEM | linearised FCI parallel operator in the field split; FCI plane multigrid; smoothers for nonsymmetric, indefinite blocks |
| U4 | FortNum | BDF/IMEX steppers with error control, Krylov phi-functions, device-resident Krylov |
| U5 | FortSparse | nonsymmetric AMG or block smoothers for the MHD Schur blocks |

Generated CUDA kernels sit behind the ISO_C_BINDING C ABI of
[tooling.md](tooling.md) and are verified against the CPU reference path.

## Milestones

| ID | content | depends on |
|---|---|---|
| RM0 | symmetry-reduced chart; GS generalised (axisymmetric regression, helical GS oracle) | - |
| RM1 | FortSym leaves of the full residual; symmetry-reduced 2D residual; RM-B1, RM-B2 in 2D | RM0 |
| RM2 | 2D stationary and BDF2 solves; RM-B3, RM-B5, RM-B6, RM-B8 in 2D; helically symmetric straight-cylinder flux-pumping state by Newton and RM-B10 in 2D | RM1 |
| RM3 | 3D residual on the meridional complex x Fourier; RM-B1, RM-B2, RM-B7 | RM1 |
| RM4 | field-split preconditioner; RM-B4 at linear amplitude | RM3, U3, U5 |
| RM5 | 3D nonlinear stationary and evolution; RM-B3 (FCI decision), RM-B4 saturated, RM-B5, RM-B6, RM-B8 in 3D; slow-manifold path, RM-B10 | RM4 |
| RM6 | GLL matrix-free path, CPU and generated CUDA; exponential integrators; RM-B1, RM-B2, RM-B4, RM-B7 rerun on it; RM-B9 | RM5, U1, U2, U4 |
| RM7 | physics: straight-cylinder flux pumping, sawtoothing versus flux-pumping states, then published tokamak cases | RM2, RM6 |
| RM8 | Figure-8 stellarator on the curve-following chart: RM-B1, RM-B6 and a resistive island case, equilibrium checked against GVEC | RM6, `g_frame_figure8_equilibrium` |
| RM9 | ladder rung L2: axisymmetric torus evolution and stationary state without current source; RM-B11 | RM2 |
| RM10 | ladder rung L0: surface-averaged background; RM-B12 | RM9, [transport-multiscale.md](transport-multiscale.md) |
| RM11 | ladder rung L3 (proposal): few harmonics on a radial mesh, admitted or dropped by its defect and cost against L4 | RM5 |
| RM12 | ladder rung LM (hypothesis): RM-B13 ledger on stationary states of L1 and L4; mean-field closure of L0 from it, checked on states not used for calibration and against RM-B10 | RM10, RM2, RM5 |

## Later: extended MHD and edge (not scheduled)

- Diamagnetic terms (1/1-mode stabilisation, sawtooth onset), then the Hall and
  electron-pressure terms of a two-fluid Ohm's law (crash times). The Hall
  term does no work in the mixed-current form ((J x B).J = 0), so the energy
  ledger stays exact. Rotating states are stationary in a co-moving frame and
  are solved with the frequency as an extra unknown. Gates: RM-B4 rerun with
  these terms and linear drift-tearing and whistler dispersion oracles.
- ELM cycles need, in addition, X-point geometry with scrape-off layer and
  sheath conditions, a free boundary with wall currents, neoclassical edge
  flows and n up to 20-40, as in the JOREK ASDEX Upgrade ELM-cycle simulations
  (Cathey et al., Nucl. Fusion 60 124007, 2020).
- The L-H transition is turbulence physics and outside this model; it needs
  global edge turbulence (gyrokinetic or drift-reduced fluid), which may later
  reuse this path's FCI and X-point infrastructure.

## Later: kinetic coupling (not scheduled)

Target: a closure valid from the core to the wall. The reference is a
drift-kinetic parallel closure with collisions computed by MILCA's steady
drift-kinetic solver, valid at any collisionality; cheaper closures
(Braginskii in the scrape-off layer and divertor, flux-limited or Landau-fluid
closures in the core, gyrofluid closures where k_perp rho_i ~ 1) are used where
their defect estimate against the reference is within budget, switched locally
in space and time with conservative interfaces, under the consistency and
monitor requirements above. Where perpendicular MHD itself fails (structures
of width ~ rho_i), the promotion target is gyrokinetics. The rigorous
counterpart (fluid limits uniform in collisionality) is pursued in a separate
private research project.


Kinetic effects enter later through the closure interface, after RM7, and are
not part of RM0-RM8:

- neoclassical resistivity and bootstrap current from MILCA's steady
  drift-kinetic solver (see [drift-kinetics.md](drift-kinetics.md)), checked
  against published neoclassical limits; the bootstrap current enters as J_s;
- a kinetic or Landau-fluid parallel heat-flux closure replacing the flux
  limiter;
- fast ions (beam ions, alphas) through hybrid kinetic-MHD coupling. Whether
  they are represented by MILCA's own particle, Eulerian or semi-Lagrangian
  carriers ([representations.md](representations.md)) or supplied by an external
  Monte Carlo orbit code such as SIMPLE or GORILLA through a provider interface
  is decided when the work starts. MILCA does not add a specialised
  guiding-centre orbit-following Monte Carlo solver of its own.
