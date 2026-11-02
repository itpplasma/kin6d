# Nonlinear visco-resistive MHD: general 3D and symmetry-reduced 2D

Status: design, decided 2026-11-02. Implemented: the RM0 fixed-boundary
symmetry-reduced GS and the RM1a pointwise weak-form leaves; assembled nonlinear
residuals and solvers remain milestones. This document owns the
nonlinear visco-resistive part of [R12](../ROADMAP.md#r12--equilibrium-mhdkinetic-mhd-vacuum-and-circuits)
item 12; Hall/XMHD stays with R12. Scope: tokamak and stellarator core MHD on
resistive time scales.

## Why

Nonlinear resistive MHD of the tokamak core on the resistive time scale is
expensive in established codes. The published 3D JOREK simulation of flux
pumping in ASDEX Upgrade (Zhang et al., Nucl. Fusion 65 066001, 2025) used
about 0.5 M core-hours on Marconi for 700 ms, about 6000 core-seconds per time
step. One reproduction on a 48-core node (in a private project, not in this
repository) observed:

- with chi_par/chi_perp ~ 1e10, GMRES with a preconditioner that factorises
  each toroidal harmonic separately rose from ~8 to its 300-iteration limit
  within 35 steps, when the n = 1 magnetic energy was ~4e-7 of the n = 0
  energy;
- splitting the harmonics into two blocks, a looser tolerance or a smaller
  time step did not restore convergence at acceptable cost;
- factorising all harmonics together converged, but the factorisation went
  stale within ~10 steps while the mode grew.

Hypothesis, from that one reproduction: the cost is dominated by
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
through mapped tensor patches on the curve-following chart that KIN6D already
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
is checked by the [S0 consumer](../ROADMAP.md#s0--stabilityresponse-consumer-cross-cutting-starts-with-r1r3).
Oscillating states need shooting or harmonic balance; that is a later step.

Time evolution: BDF2 first (the scheme of the JOREK comparisons), then IMEX and
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
non-inductive current J_s as a prescribed source. Benchmarks against JOREK use
JOREK's prescribed coefficients; physics runs state the closure they use.

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
| RM-B8 | JOREK | public JOREK regression cases (tearing and kink in circular geometry) and a straight-cylinder flux-pumping case, matched model and resolution; 2D part in RM2, 3D part in RM5 |
| RM-B9 | performance | time to matched observable accuracy against JOREK on the same node, strong scaling, memory, CPU and GPU |
| RM-B10 | slow manifold | quasi-stationary sequence of driven stationary states against full BDF2 evolution on the straight periodic cylinder: q profile and island width to the stated tolerance, defect estimate bounds the observed difference, switch to time evolution at mode onset |

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
evolution at least 10x faster than JOREK at matched accuracy. Physics
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
| RM2 | 2D stationary and BDF2 solves; RM-B3, RM-B5, RM-B6, RM-B8 in 2D; helically symmetric straight-cylinder flux-pumping state by Newton and RM-B10 in 2D, against the JOREK cylinder runs | RM1 |
| RM3 | 3D residual on the meridional complex x Fourier; RM-B1, RM-B2, RM-B7 | RM1 |
| RM4 | field-split preconditioner; RM-B4 at linear amplitude | RM3, U3, U5 |
| RM5 | 3D nonlinear stationary and evolution; RM-B3 (FCI decision), RM-B4 saturated, RM-B5, RM-B6, RM-B8 in 3D; slow-manifold path, RM-B10 | RM4 |
| RM6 | GLL matrix-free path, CPU and generated CUDA; exponential integrators; RM-B1, RM-B2, RM-B4, RM-B7 rerun on it; RM-B9 | RM5, U1, U2, U4 |
| RM7 | physics: straight-cylinder flux pumping, sawtoothing versus flux-pumping states, then published tokamak cases | RM2, RM6 |
| RM8 | Figure-8 stellarator on the curve-following chart: RM-B1, RM-B6 and a resistive island case, equilibrium checked against GVEC | RM6, `g_frame_figure8_equilibrium` |

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
drift-kinetic parallel closure with collisions computed by KIN6D's steady
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

- neoclassical resistivity and bootstrap current from KIN6D's steady
  drift-kinetic solver (see [drift-kinetics.md](drift-kinetics.md)), checked
  against published neoclassical limits; the bootstrap current enters as J_s;
- a kinetic or Landau-fluid parallel heat-flux closure replacing the flux
  limiter;
- fast ions (beam ions, alphas) through hybrid kinetic-MHD coupling. Whether
  they are represented by KIN6D's own particle, Eulerian or semi-Lagrangian
  carriers ([representations.md](representations.md)) or supplied by an external
  Monte Carlo orbit code such as SIMPLE or GORILLA through a provider interface
  is decided when the work starts. KIN6D does not add a specialised
  guiding-centre orbit-following Monte Carlo solver of its own.
