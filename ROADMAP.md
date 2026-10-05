# kin6d roadmap

Status: staged delivery plan, 2026-10-05.

The roadmap is ordered by scientific leverage and verification value, not by
the final breadth of the code. The steady kinetic programme is the main spine.
Student dynamics, turbulence and alternative representations are introduced
early where they exercise the same operators cheaply.

The long-range product goal is a **one-stop shop for classical fusion and
Solar-System plasma physics**. Breadth must come from reusable model,
representation, analysis, transport, boundary and coupling primitives rather
than domain-specific solver silos. A new application area should normally add
a thin configuration plus one or two genuinely reusable capabilities to the
common engine.

## Rules for every milestone

A capability is promoted only when it has:

1. a precisely named mathematical model and domain;
2. an analytic/manufactured/structural oracle;
3. an independent external or cross-representation oracle when available;
4. primal/JVP/VJP behavior for continuous production operators;
5. a stated rigorous-certificate status;
6. a small correctness case and a representative performance/memory case;
7. preserved earlier limiting cases.

An optimization or representation change is not allowed to obscure the
underlying mathematical operator.

## R0 — retain and harden the first 1X1V slice

**Implemented baseline:** prescribed-field free streaming and a direct
source-balanced Krook steady state with generated derivative leaves and a
continuous snapshot interval bound.

Near-term work:

- keep the existing exact characteristic and manufactured steady oracles green;
- improve certificate packaging without requiring Lean;
- retain strict/reference and optimized actions;
- keep the demo understandable enough for a student to trace end to end.

Do not generalize the current Fourier state layout into a universal 6D type.

## R1 — self-consistent 1X1V Vlasov--Poisson laboratory

Add periodic electrostatic Poisson with explicit neutrality compatibility and
gauge. Required cases:

- single Fourier Poisson mode;
- linear plasma oscillation;
- Landau damping;
- two-stream instability;
- a genuinely self-consistent BGK equilibrium.

The steady residual is primary. Transients use the same field/transport
operators. Distinguish exact linear dispersion roots from finite-time nonlinear
solutions.

This milestone creates the first implicit-state adjoint:
\(R_U^\dagger\lambda=Q_U^\dagger\).

## R2 — representation laboratory: one physics problem, several carriers

Use the same 1X1V Vlasov--Poisson cases to add representation alternatives one
at a time:

1. Eulerian/reference continuum;
2. conservative/structure-preserving semi-Lagrangian characteristic path;
3. deterministic particle quadrature with PIF;
4. geometric/charge-conserving PIC;
5. finite blobs/remapping;
6. delta-f/control-variate variants;
7. QMC/MC sampling as explicitly statistical variants.

The acceptance goal is not to crown one method universally. Compare accuracy,
structure, memory, time and noise at equal physical observables. Preserve the
continuum reference.

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

## S0 — stability/response consumer (cross-cutting, starts with R1/R3)

Do not postpone stability to the late MHD branch. Establish the generic
matrix-free analysis consumer as soon as a nontrivial JVP exists.

Early sequence:

1. 1X1V plasma-wave eigenvalue/dispersion case;
2. HM/drift-wave linear spectrum before nonlinear turbulence;
3. simple ideal-MHD/Alfvén eigenproblem on an analytic equilibrium;
4. forced response/resolvent;
5. nonmodal/optimal-growth case using JVP+VJP;
6. continuation to a marginal-stability point.

Later the same consumer covers ideal kink/ballooning, tearing/RWM, Alfvén
eigenmodes, NTM reduced/nonlinear evolution, and DK/GK/FK microinstabilities.
See [stability.md](docs/stability.md).

## TR0 — macro transport core and cheap integrated modelling (cross-cutting, starts before R6)

Add a conservative transport state early rather than waiting for global
kinetics. Begin with exact 1D diffusion/advection/source cases but keep the
FortFEM state capable of 2D/3D extension.

The first useful integrated chain is

\[
\text{profiles/equilibrium}
\rightarrow
\text{prescribed or neoclassical/QL closure}
\rightarrow
\text{1D transport}
\rightarrow
\text{updated profiles/equilibrium}.
\]

Reproduce the KAMEL/KiLCA/KIM -> QL-Balance pattern in memory as an internal
oracle. Then generalize the closure interface beyond scalar coefficients.

Before or alongside the real-device DK milestone, add:

- neoclassical Onsager coefficients as a transport provider;
- linear-stability/QL closure provider;
- closure Jacobian/JVP for implicit profile solves;
- first 1.5D equilibrium--transport loop.

Next grow to 2D/3D macro transport and nonlocal global-orbit kernels.
A GORILLA-like characteristic engine may provide transition/response kernels
without assuming local radial diffusion. Then add a prescribed electromagnetic
turbulence cross-spectrum along global orbits before coupling a live
microturbulence provider.

The multiscale ladder is offline table -> in-memory loose coupling -> iterated
closure -> separated-timescale local/global microproblems -> sparse expensive
microphysics sampling -> fully coupled global evolution. See
[transport-multiscale.md](docs/transport-multiscale.md).

## UX0 — configuration, reproducibility and operator UX (cross-cutting)

Establish the user-facing contract before the model space becomes large:

1. a typed Fortran configuration model with Fortran namelist as the canonical
   editable syntax;
2. validate and fully resolved configuration output so defaults/automatic
   choices are explicit and recoverable;
3. a machine-readable schema generated from the same Fortran-owned semantics;
4. NetCDF-4-compatible HDF5 result files through FortIO, always embedding the
   original and resolved input plus build/input provenance;
5. explicit GNU, Intel LLVM, NVIDIA and LLVM Flang build/test lanes, with
   certifier floating-point flags isolated from fast candidate builds;
6. simple CLI workflows for run/inspect/extract-input before any GUI;
7. optional Python notebook/server and TypeScript browser layers only after the
   schema/result contract is stable.

The normal interaction should expose physical model, observable and requested
accuracy before low-level numerical knobs. Experts retain an advanced numerical
view. Browser use over an SSH tunnel or an existing cluster portal is preferred
to a separate desktop application; scheduler metadata remains outside the
physics namelist.

See [interface-reproducibility.md](docs/interface-reproducibility.md).

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

Move to frozen real tokamak/stellarator field/equilibrium candidates and a
drift-kinetic model. Import is sufficient for this rung, but imported data are
not trusted physics: normalize conventions, retain provenance/input
uncertainty, and evaluate the field/geometry checks required by the DK claim.
Native equilibrium re-solves are added in R12 and should later be used to
separate imported interpolation/reconstruction error from DK transport error.

This is also the first production libneo dependency gate. Reuse existing
libneo readers/field utilities behind a KIN6D adapter rather than cloning them.
Before promoting each exercised path, audit and repair generic defects upstream
in libneo with independent tests. Add reusable spatial derivatives or
representation-level error diagnostics there when the need is generic; keep
KIN6D's DK residual/adjoint/certificate mathematics in KIN6D. Existing libneo
consumers should remain compatible and inherit the fixes.

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

## R12 — equilibrium, MHD/kinetic-MHD, vacuum and circuits

MHD grows as a parallel bulk-model branch using the same physical-space
geometry/field foundation. Equilibrium is a stationary solve, not a VMEC/EQDSK
mode. The staged path is:

1. vacuum magnetostatics from explicit coils/currents and conducting boundaries;
2. native axisymmetric Grad--Shafranov, first fixed boundary and then free
   boundary with explicit external-coil/vacuum data;
3. EQDSK/gEQDSK import as a candidate plus residual checking, native re-solve
   and certificate comparison;
4. VMEC/DESC/GVEC import as nested-surface candidates and cheap
   representations/preconditioners, followed by physical-space residual checks;
5. a general physical-space three-dimensional ideal-MHD stationary solver that
   does not require global nested surfaces and can represent islands/stochastic
   regions when the equations permit them;
6. explicit comparison with HINT/PIES-like nonnested full-pressure relaxation
   and, for a different relaxation model, SPEC/SPECTRE/MRxMHD cases;
7. resistive-MHD and Hall/XMHD stationary/evolution paths where reconnection and
   history are part of the selected physics;
8. kinetic-MHD and kinetic continuation where pressure/current closure is no
   longer adequately fluid;
9. ideal/resistive linear response, eigenmodes and nonlinear stability through
   the common S0 consumer;
10. conducting structures, vacuum/BEM, coils and external circuits.

For static smooth ideal MHD, retain the actual constraint
\(B\cdot\nabla p=0\): removing the nested-surface representation does not
license arbitrary smooth pressure gradients in ergodic volumes. Different
stationary solutions are stored as different equilibria with their own initial
state/continuation provenance; no general probabilistic branch machinery belongs
in KIN6D.

ARES/MEPHIT/MARS/GPEC and later JOREK/NIMROD-type cases remain response and
evolution landmarks. VMEC/DESC/GVEC, HINT/PIES and SPEC/SPECTRE are equilibrium
comparators with different model assumptions. Do not reimplement their entire
feature sets before a KIN6D use case requires a shared capability.

See [equilibria-inference.md](docs/equilibria-inference.md).

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

## R13A — one-stop-shop coverage completion (cross-cutting, demand-driven)

After the core fusion and Solar-System branches are usable, close remaining
classical-plasma coverage gaps by adding **generic capabilities**, not a new
named solver for each community. Prefer an abstraction that serves both fusion
and space whenever the underlying mathematics is shared.

### Strongly driven thermal/current evolution

Support rapid coupled changes of temperature, current, impurities, fields and
conducting structures using the existing transport, MHD/XMHD, reaction,
radiation, energetic-particle and wall/circuit machinery.

This generic capability should cover strongly driven events without hard-coding
one device scenario. The reusable pieces are:

- stiff source/reaction/current evolution;
- rapid equilibrium/field change;
- impurity and radiation coupling;
- relativistic energetic-particle production/transport where needed;
- conductor/wall response and load diagnostics.

### Edge/pedestal stability and finite-amplitude transition

Combine the common equilibrium, stability, transport and nonlinear-evolution
consumers so edge-gradient systems can be studied from linear margin through
reduced finite-amplitude evolution to nonlinear saturation.

Do not create a special pedestal solver. The added value is one consistent
path:

\[
\text{equilibrium}
\rightarrow
\text{linear/nonmodal stability}
\rightarrow
\text{reduced nonlinear branch}
\rightarrow
\text{full nonlinear model}
\rightarrow
\text{transport feedback}.
\]

### Partially ionized plasma

Make partial ionization a first-class multispecies/reaction regime using shared
charged/neutral populations rather than a separate atmosphere/edge package.
Needed reusable physics includes:

- ion-neutral drag and charge exchange;
- ionization/recombination;
- anisotropic/Hall/Pedersen-like conductivity closures;
- neutral flow and momentum/energy exchange;
- transition from one-/multi-fluid closures to kinetic charged and neutral
  species.

This serves fusion edge/divertor, chromospheric, ionospheric, cometary and
planetary applications.

### Conducting-boundary and ionosphere/object coupling

Generalize the existing wall/sheath/circuit boundary machinery to a finite
conducting environment coupled to field-aligned currents and plasma fluxes.

The reusable abstraction is a surface/boundary state with:

- conductance/capacitance/current closure;
- sheath/floating-potential behavior;
- photo/secondary emission where requested;
- imposed or self-consistent tangential electric response;
- coupling to an external neutral/fluid layer.

This should cover fusion conductors, ionospheric electrodynamics, moons and
spacecraft without separate numerical infrastructures.

### Focused energetic-particle transport

Add a reduced energetic-particle transport model between orbit-following and
full kinetic phase-space evolution. It may evolve pitch-angle/energy/action
coordinates with streaming, focusing, diffusion and source/loss terms.

The same transport operator should serve fast ions, radiation-belt/SEP-like
populations, pickup ions and shock-accelerated particles when their assumptions
match. It is another reduced species model/transport closure, not a special
space-physics code.

### Wave-kinetic and statistical-spectrum evolution

Promote the spectral closure interface from prescribed turbulence statistics
to an evolving wave/turbulence state, schematically

\[
\partial_t N
+
\nabla_k\omega\cdot\nabla_x N
-
\nabla_x\omega\cdot\nabla_k N
=
S_{\rm wave}+C_{\rm wave}-D_{\rm wave}.
\]

Couple \(N(x,k,t)\), cross-spectra and particle distributions through the same
quasilinear transport contract already used by RF, turbulence and energetic
particles. This provides an intermediate level between fixed spectra and
fully nonlinear field evolution.

### Magnetic-topology diagnostics as a service

Add a reusable field-topology analysis layer on top of the common field-line
map/geometry service. Candidate diagnostics include resonances, island widths,
connection lengths, stochasticity/overlap, field-line diffusion,
finite-time-Lyapunov/squashing measures and topology/helicity observables.

These are diagnostics and closure inputs for transport/stability; they must not
become a second magnetic-geometry implementation.

### Charged grains and embedded bodies

Generalize the particle/surface machinery to finite charged macroparticles or
bodies with mass, charge state and optional size/material state. Reuse the same
external-force, charging, collision, sheath and wall interaction infrastructure.

This covers dusty/plume/ring/spacecraft-like applications without perturbing
the main fusion path. Implement only when a concrete case requires it.

### Rotating/stratified bulk plasmas

Extend bulk-fluid/MHD models with reusable gravity, rotation, buoyancy,
thermodynamic/EOS and anisotropic-conduction terms when needed. This enables
stratified convection/dynamo and rotating space-plasma applications while
remaining the same bulk-model infrastructure used elsewhere.

### Admission criterion for coverage extensions

A coverage extension belongs in KIN6D when it can be expressed mostly through
existing common state, transport, stability, field, reaction, surface or
particle interfaces and adds a reusable primitive with at least two plausible
consumers. Domain-specific data preparation or diagnostics may remain thin
application adapters.

Do **not** internalize large adjacent disciplines merely to claim completeness.
In particular, full non-LTE radiation transfer, reactor neutronics,
solid-state material microphysics, strong-field QED and GR remain external
coupled systems unless a later scientific target justifies a narrow reusable
interface.

## R14 — diagnostics, TIAGO inference integration, design and adaptive fidelity

Synthetic diagnostic forward operators enter KIN6D once their underlying
physics exists, but reconstruction/inference is owned by TIAGO. KIN6D exposes
states/observables together with matrix-free JVP/VJP and implicit-equilibrium
adjoint actions; TIAGO owns measurements, likelihood/covariance models,
diagnostic/calibration nuisance parameters, priors and posterior uncertainty.

The first integration should keep TIAGO's current VMEC++ reconstruction as an
independent reference and add a provider-neutral interface through which KIN6D
can supply Grad--Shafranov and later general MHD/DK/GK/FK stationary states.
The default UQ remains a local Laplace/Gauss--Newton approximation. When the
same data support several distinct KIN6D equilibria, fit them separately and
optionally report a finite mixture of local Laplace approximations; do not
average distinct equilibria into one state. More elaborate sampling is
deferred until a concrete within-solution posterior is demonstrably
non-Gaussian.

Report reconstruction/measurement/calibration uncertainty separately from
KIN6D numerical certificates and model-reduction error bands. Combine them only
under an explicitly justified probabilistic/mathematical model.

Implicit adjoints also support stellarator/coil/divertor design and control.
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

### Interchange, configuration and result contract

Use libneo as the shared adapter/convention layer for established machine and
code formats where it already owns the capability. KIN6D consumes normalized
physical data and records the source format/conventions; it does not turn VMEC,
EQDSK, SPECTRE, JOREK or another source into a fundamental internal model.

Treat every new KIN6D dependency on a libneo code path as a boy-scout
modernization opportunity: fix generic bugs and convention errors upstream,
add missing reusable tests/derivatives/error diagnostics there, and keep
existing downstream consumers working. The KIN6D adapter remains narrow so a
backend can still be replaced later; swappability is not a reason to duplicate
or abandon reusable libneo mathematics.

Keep Fortran namelist configuration and NetCDF-4-compatible HDF5 results stable
across physics milestones. Every result retains original/resolved input and
provenance. Optional Python/web UX is layered above this contract.

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
