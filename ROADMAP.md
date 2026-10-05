# kin6d roadmap

Status: staged delivery plan, 2026-10-05.

The roadmap is ordered by scientific leverage and verification value, not by
the final breadth of the code. The steady kinetic programme is the main spine.
Student dynamics, turbulence and alternative representations are introduced
early where they exercise the same operators cheaply.

## Rules for every milestone

A capability is promoted only when it has:

1. a precisely named mathematical model and domain;
2. an analytic/manufactured/structural oracle;
3. an independent external or cross-representation oracle when available;
4. primal/JVP/VJP behavior for continuous production operators;
5. bounds/status for reported states, fields, observables and sensitivities;
6. a small correctness case and a representative performance/memory case;
7. preserved earlier limiting cases.

An optimization or representation change is not allowed to obscure the
underlying mathematical operator.

Certification and differentiation start with the operator, rather than after
the solver is complete. Report unavailable bounds explicitly. Pointwise plots
need pointwise/cellwise bounds; weak particle observables and statistical
confidence statements retain their distinct meanings.

The practical workhorse target is about 1% error on a named observable with a
declared scale, including an absolute scale near zero. Use a validated coarse
rounding budget when sufficient; tighter arithmetic audits are optional while
deterministic containment for certified results remains required. Optimize
measured total setup/solve/certificate/requested-derivative time at matched
model, accuracy and claim strength. No speed ranking is asserted in advance.

## Immediate delivery order

1. Establish reusable benchmark support in the owning upstream tools: fixed
   inputs, warm/repeated timings, measurement variability, memory, candidate,
   certificate and derivative costs. KIN6D supplies the physical cases.
   Initial CPU fixtures are delivered in FortNum/FortSym, with independent
   arithmetic/derivative oracles and raw samples; KIN6D has a resolution sweep.
   Production arithmetic selection and external speed comparisons remain open.
2. Generalize prescribed-field damped transport to independently supplied
   sources. Bound the continuous reconstructed residual and source error using
   a proved inverse/stability estimate, rather than a known-solution formula.
   Delivered for exact supplied Fourier/P1 sources, with whole-cell state and
   tangent bounds, generated nodal JVP/VJP, a cosine observable, independent
   tests and a native demo/audit. A continuous physical goal adjoint and
   external analytic-source reconstruction bounds remain separate gates.
3. Extend FP to independently supplied forcing with declared boundary flux,
   compatibility and normalization. Certify flux/reconstruction defects through
   the applicable stability estimate, including boundary and source errors.
4. Introduce the simple representation comparisons below alongside this steady
   spine. General nonlinear validation follows once its stability, nullspace
   and unresolved-complement hypotheses can be checked.

Retain all analytical laboratories as independent regression oracles. A small
algebraic residual or agreement between representations does not replace a
continuous certificate.

## R0 — retain and harden the first 1X1V slice

**Implemented baseline:** prescribed-field free streaming and Krook balance,
single-mode Poisson, a coupled source-balanced Poisson--FP family, and bounded
homogeneous OU/FP equilibrium. Generated derivative leaves and reviewed
continuous bounds cover the declared analytical cases. The
[general transport case](docs/general-transport.md) additionally certifies
unknown solutions for independently supplied Fourier/P1 sources through the
continuous residual and inverse bound.

Near-term work:

- keep the existing exact characteristic and manufactured steady oracles green;
- retain the general transport certificate and extend FP to independently
  supplied forcing with flux/stability, derivative and observable bounds;
- make the approximately 1% observable target a useful default;
- improve bound/status packaging without requiring Lean;
- retain strict/reference and optimized actions;
- keep the demo understandable enough for a student to trace end to end.

Do not generalize the current Fourier state layout into a universal 6D type.

## R1 — self-consistent 1X1V Vlasov--Poisson laboratory

The implemented single Fourier mode and manufactured coupled family retain
exact neutrality compatibility and potential gauge. Extend these operators to
general self-consistent 1X1V states. Required additional cases:

- linear plasma oscillation;
- Landau damping;
- two-stream instability;
- a genuinely self-consistent BGK equilibrium.

The steady residual is primary. Transients use the same field/transport
operators. Distinguish exact linear dispersion roots from finite-time nonlinear
solutions.

General coupled states require the implicit-state adjoint:
\(R_U^\dagger\lambda=Q_U^\dagger\).

## R2 — representation laboratory: one physics problem, several carriers

Begin before the full R1 case list is complete:

1. Exact free streaming: continuum versus deterministic characteristic
   particles, with analytical density modes and parameter sensitivities.
2. Neutral single-mode Poisson: continuum projection, direct PIF and simple
   local-basis PIC deposit/solve/gather, with analytical field and energy.
   This is a coupling test, without a self-consistent trajectory claim.
3. Short self-consistent Vlasov--Poisson: the same initial distribution,
   boundaries and observation times across continuum, PIF and PIC. State the
   actual trajectory-bound status; validated field solves alone are insufficient.

Use deterministic sampling first. Share physical transport, field laws and
pairings, while exposing quadrature, basis/shape and time errors separately.
Carry generated JVP/VJP from the start and declare PIC shape/event domains.
Compare total time, memory, bounds and requested sensitivities at matched
observable accuracy. Preserve the continuum reference.

Conservative semi-Lagrangian paths, geometric/current-conserving PIC,
blobs/remapping, delta-f/control variates and randomized QMC/MC follow their
concrete use cases. Early electrostatic PIC does not establish Maxwell current
continuity or geometric long-time structure. Defer NUFFT, adaptive particles,
stochastic collisions and GPU/MPI until these small cases require them.

PIF/PIC are particle-in-basis specializations, not new physics modes. A generic
particle runtime is extracted upstream only after at least two concrete KIN6D
consumers demonstrate the stable abstraction.

## R3 — early turbulence gallery

In parallel with R1/R2, build cheap structure-preserving demonstrations:

- Hasegawa--Mima, starting with one linear mode and energy/enstrophy ledger;
- Hasegawa--Wakatani;
- interchange/blob dynamics;
- minimal ITG gyrofluid;
- electromagnetic drift-Alfvén.

The nonlinear bracket comes from one FortSym definition with verified
antisymmetry/pairing behavior and alias-controlled projection. A movie is not a
certificate; analytic dispersion/growth and discrete free-energy/source/sink
ledgers are mandatory.

This branch is a technology and teaching demonstrator, not the main research
spine.

## R4 — 0X3V collision laboratory

Introduce the future velocity representation in the smallest domain. Progress:

\[
{\rm Lorentz}
\rightarrow
{\rm LBD/Dougherty}
\rightarrow
{\rm Landau\ test}
\rightarrow
{\rm full\ linearized\ Landau}
\rightarrow
{\rm nonlinear\ Landau}
\rightarrow
{\rm Beliaev\!-\!Budker}.
\]

Use spherical harmonics and species-specific radial momentum discretization
where they remain the best reference. Validate number/momentum/energy
invariants, entropy production, nullspaces, weighted self-adjointness and
analytic harmonic/radial limits.

Factor radial elliptic operators once and benchmark many-RHS batched solves.
Stochastic binary-collision variants may be added as independent particle
realizations after the deterministic operator is established.

## R5 — 1X3V steady transport

Combine spatial transport with the collision laboratory. Targets include:

- parallel conductivity/Spitzer--Härm limits;
- mirror/trapped-particle test cases;
- slowing-down/source-balanced states;
- classical transport coefficients;
- simple sheath/inflow cases where the boundary model is explicit.

This is the first nontrivial stationary kinetic transport publication rung.
Develop goal-oriented adjoint bounds for selected fluxes.

## R6 — real-device DK neoclassical transport

Move to imported real tokamak/stellarator equilibria and a drift-kinetic model.

Required comparisons:

- NEO-2;
- SFINCS;
- DKES/KNOSOS/MONKES where scope matches;
- rabe/asymptotic limits where applicable.

Required structure:

- independently computed reciprocal coefficients;
- Onsager--Casimir transformation under field/time reversal;
- nonnegative entropy production in the correct force/flux pairing;
- collision invariant and adjoint tests;
- explicit geometry/input uncertainty.

This is the first major production target.

## R7 — gyrofluid, Landau-fluid and gyromoment hierarchy

Promote the early turbulence branch into a systematic hierarchy derived from a
named gyrokinetic parent where available.

Sequence:

1. low-order gyrofluid;
2. Landau-fluid closures;
3. Hermite--Laguerre gyromoments;
4. controlled increase of retained moments;
5. comparison/convergence toward GK observables.

Reference Scott/GEM free-energy structure and Kendl TIFF/T3FF cases only where
their exact model assumptions match. Moment count alone is not an accuracy
claim.

## R8 — gyrokinetics

Add DK/GK species mixtures and a consistent reduced field action. Progress from
small stationary/linear cases to transient turbulence and global/full-f cases.

Required:

- gyrocenter phase-space measure and push-forward;
- polarization/magnetization consistency;
- moment/GK reduction comparisons;
- local analytic wave/instability limits;
- cross-code comparison with stella/GS2/GENE/GX and later GENE-X/XGC where
  model scope matches.

Delta-f, full-f, continuum and particle GK representations may coexist; do not
duplicate the physical model for each representation.

## R9 — global FK and nonlinear Coulomb reference

Build toward multispecies FK on general mapped geometry:

- nonlinear Landau;
- relativistic Beliaev--Budker;
- Poisson/quasineutral/Darwin/Maxwell field choices;
- fixed-ion/electron-only and energetic-particle configurations;
- direct steady states before unnecessarily long transients.

Particle, blob and continuum FK are alternative representations. Use whichever
reaches the requested observable accuracy cheapest; retain cross-representation
validation.

## R10 — edge, neutrals, reactions and walls

Extend the same engine through:

- open field lines/X-points/FCI;
- prescribed/closure/resolved cross-field transport;
- neutral fluid, deterministic kinetic and MC kinetic backends;
- charge exchange, ionization/recombination, molecular channels and CR models;
- sheath hierarchy through locally resolved FK;
- reflection, recycling, sputtering, pumping and wall state;
- radiation losses.

SOLPS-ITER/EIRENE, EMC3-EIRENE, SOLEDGE3X-EIRENE, GBS, GRILLIX and GENE-X are
comparison points, not modes.

## R11 — heating, energetic species and nuclear sources

Add when a concrete experiment/design case requires them:

- prescribed quasilinear diffusion;
- ray/QL and later full-wave kinetic response;
- ECRH/ECCD, ICRH and LHCD;
- NBI neutral attenuation, birth distribution and active fast ions;
- D-T/D-D/D-3He rates and alpha source distributions;
- runaway/relativistic radiation-reaction pieces where required.

KAMEL/KiLCA/KIM and appropriate external RF/FP codes provide independent
oracles.

## R12 — MHD, kinetic-MHD, vacuum and circuits

MHD grows as a parallel bulk-model branch using the same geometry/field
foundation:

1. imported equilibrium candidates;
2. native axisymmetric Grad--Shafranov;
3. ideal linear response/eigenmodes;
4. resistive response;
5. kinetic-MHD;
6. Hall/XMHD/multi-fluid;
7. nonlinear extended-MHD cases driven by an actual scientific target;
8. conducting structures, vacuum/BEM, coils and circuits.

ARES/MEPHIT/MARS/GPEC and later JOREK/NIMROD-type cases are validation
landmarks. Do not reimplement their entire feature sets before a KIN6D use case
requires them.

## R13 — Solar-System plasma configurations

Reuse the same hierarchy with open domains and additional forces:

- Parker wind/spiral;
- Harris/GEM reconnection;
- CGL and Landau-fluid closures;
- hybrid Vlasov/PIC configurations;
- collisionless shocks and pickup ions;
- kinetic charge exchange/ENA diagnostics;
- photoionization, mass loading and surface charging.

Vlasiator, Gkeyll, Smilei, WarpX and related codes provide cross-validation.

## R14 — diagnostics, inference, design and adaptive fidelity

Synthetic diagnostic operators enter once their forward physics exists.
Implicit adjoints support parameter inference, equilibrium reconstruction,
stellarator/coil/divertor design and control.

The long-range adaptive objective is

\[
\arg\min_{\rm model,representation,resolution} T
\]

subject to a requested bound on a physical observable. Adjoint-weighted defects
rank which numerical or physical refinement buys the most observable accuracy
per cost.

Automatic model switching requires actual reduction/interface remainder
bounds; until then it is an experimental heuristic with explicit status.

## Permanent numerical workstreams

### Symbolic mathematics

FortSym remains the source of project-owned local mathematics. Extend it
upstream only on concrete need. Keep mathematical DAG, indexed reduction IR and
hardware schedule separate.

### Matrix-free/tiled execution

Develop sum-factorized FortFEM tensor actions, batched many-RHS operators and
KeOps/FlashAttention-style tiled reductions. Never materialize pairwise
interaction matrices or global kinetic Jacobians when a streamed reduction is
the actual mathematics.

### Performance

For each rung retain:

- small correctness benchmark;
- representative memory/performance benchmark;
- best-known hardware/compiler Pareto baseline;
- solver and certifier profiles;
- forward/JVP/VJP timing where derivatives are material.

Optimize the owning dependency rather than accumulating KIN6D workarounds.

### Certification

Prioritize useful Fortran interval/ball certificates. Develop Lean only for
generic theorem/checker value once the numerical certificate shape is stable.
SPD/coercive pieces should exploit residual-energy, equilibrated-flux and
adjoint QoI bounds. Non-normal/hypocoercive pieces need their own stability
theorems.

### Literature and validation

Every new rung enters the analytical benchmark atlas and has a frozen
cross-code/literature reference where possible. Existing itpplasma codes
NEO-2, NEO-RT, SIMPLE, GORILLA, KAMEL/KiLCA/KIM, MEPHIT, ARES and rabe are
valuable retained oracles as their production roles are subsumed.

## Out of initial scope

General relativity, strong-field QED/pair creation, full reactor neutronics,
detailed solid-state material microphysics and general radiation
hydrodynamics remain external.

Chris&AI
