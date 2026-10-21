# kin6d roadmap

Status: staged delivery plan, 2026-10-06.

The roadmap is ordered by scientific leverage and verification value, not by
the final breadth of the code. The steady kinetic programme is the main spine.
Student dynamics, turbulence and alternative representations are introduced
early where they exercise the same operators cheaply.

The long-range product goal is a **one-stop shop for classical plasma
physics**: fusion, Solar-System/natural, industrial/low-temperature, beam/
accelerator, atmospheric, dusty/non-neutral and classical high-energy-density
regimes. Breadth must come from reusable model, representation, analysis,
transport, reaction, radiation, boundary and coupling primitives rather than
domain-specific solver silos.

The complete design envelope lives in
[the capability catalogue](docs/capability-catalogue.md) and
[the numerical-method catalogue](docs/numerical-methods.md). Those catalogues
are not backlogs. The roadmap maintains only a small active path chosen by
scientific value and leverage.

## Controller continuation checkpoint: October 21, 2026

The user requested that agents conclude and preserve this work for the next
controller. Recheck forge and live-job state before acting. Keep queued branches
immutable. The controller owns integration and promotion. Do not start further
layers, jobs, agents or campaigns during this conclusion checkpoint.

[CI migration #52](https://github.com/itpplasma/kin6d/pull/52) merged through the
normal queue as 5068dfcb7f7255e0dc3c0fd0c9c71a0108f865e0. All required native gates
passed on owned hardware. All PR CI skips, including full-ci-labelled PRs;
actual required execution occurs in the merge queue. Preserve the workflows,
native commands, compiler constraints and selectors. Hosted checks are advisory.
The earlier queue lost runner communication despite native map output passing
11/11; retain the infrastructure failure without fabricating a green status.

- [#39 S0 volume spectrum](https://github.com/itpplasma/kin6d/pull/39) is queued
  at 6500ddccc6611f75a0a04891127e2a457deb6943 under controller approval. Native
  21817433 passes 6/6 and four exact identities. All 488 sealed native inputs
  match qualified 4451be6a except test print prefixes. Strict workflows are
  byte-identical; only S0 build/test selectors were added. Await actual merge.
- [#35 3D readback](https://github.com/itpplasma/kin6d/pull/35) is published at
  558f1cee39c697b7422562077c2497d57bdd1a86, awaiting controller final review.
  Native 21817109 passes 14/14; nine inverse identities are exact. All 488 native
  inputs match qualified 616ab0b. Workflows remain byte-identical, with qualified
  selectors added. Do not self-promote.
- [#41 local continuum/rigid subspace](https://github.com/itpplasma/kin6d/pull/41)
  remains published at 1c9cfce9d360c04e3098a6d06d7e2a99eb0f878d. Continuation
  branch completion/kin6d41-20261021 holds
  f50308e2cefcc14c32e479894679ed7ec1103714. Native 21817477 passes 9/9, with four
  S0 and five toroidal exact identities per kind. All 489 sealed inputs retain
  production/generated/dependency bytes; two test print prefixes changed. After
  #39 actually squashes, transfer only the upper layer using prepared lower
  4451be6ad1de121b37d7f9c78611c31d59d25b8c, retarget, port selectors and compare
  seals. Controller approval permits normal queue after that verified transfer
  and a final-head receipt.
- [#44 wall-preserving prolongation](https://github.com/itpplasma/kin6d/pull/44)
  remains published at 99e66f91cbaef82033e6e85e0d8bdfb5b0465848. Continuation
  branch completion/kin6d44-20261021 holds
  db30648cd3b1765e36502a2c7be53eeb75d16f0d. Native 21817548 passes 14/14 in
  33.33 seconds; total job 2:16, guard peak 183 MiB. CAS 3319449 seed44 proves
  closed-flux/iota/position/poloidal-phase identities. Generated leaves, support
  and symbolic statements match; only a generator print prefix changed. Its
  prepared lower is 558f1cee39c697b7422562077c2497d57bdd1a86; original published
  lower was d5386131048ab9e38b68cc3c17e7cc2ed87ecf9f. After #35 actually squashes,
  transfer the prepared layer using 558f1ce, compare seals and return to the
  controller. Do not retarget the published PR before the lower squash.
- [#46 active axis columns](https://github.com/itpplasma/kin6d/pull/46) remains
  unchanged draft f580c39678cba9ed0a168eaec85e24f65d3cb04f. Its old lower is
  99e66f91cbaef82033e6e85e0d8bdfb5b0465848. After #44 actually squashes, transfer
  its independent 12-line production change and Hessian-action oracle, qualify
  affected native tests, and return evidence to the controller. No new #46 job
  was started during conclusion.
- [#47 swept BEM coupling](https://github.com/itpplasma/kin6d/pull/47) remains
  draft ea255dbe66b763bcdce01313fae467b406db9a56. Continuation branch
  completion/kin6d47-performance-20261021 holds
  d34d080e2f4125b6429b311365ffea6839d3dff9. Original 21817424 passes 6/7, with
  Landreman alone timing out at the unchanged 1800-second limit. Diagnostic
  21817533 locates adaptive geometry assembly on eight active cores. The bounded
  consumer repair passes 8/8 direct equivalence/native controls in 21817540.
  Original slow gates 21817541 remain running. HOLD promotion until those exact
  gates pass; zero external MATMUL relocations alone are not acceptance.

Never replay a whole lower stack after squash. Preserve upper-layer changes
and compare sealed native inputs before reusing evidence. All prepared candidates
above were pushed to their named continuation branches without changing the
published upper PR bases or queued #39.

BEM47 retains merged [FortFEM](https://github.com/lazy-fortran/fortfem) pin
77c6b8c2862009bd6c7dfd969abfad5da54b6467. The controller reviewed the consumer
repair: nondeferred first-jet fallback, the same cylindrical generated
position/first leaves, preserved custom/subclass dispatch, and callback-only
tiny contractions. Geometry, edge basis, quadrature, resolutions, tolerances,
case deadlines and global flags are unchanged. No upstream cache was started;
that proposal remains held until a concrete remaining need is measured.

Durable user service kin47-perf-d34d080 collects 21817541 through final sync,
report and checksums. Service kin6d-native-21817548 has collected #44's terminal
native output. Neither starts an agent or promotes code. The local continuation
index ~/data/handoffs/20261021_completion/kin6d.md links source seals, transfer
receipts, raw outputs, reports and manifests outside Git. Follow it before
repeating tests. The original unmanaged CLI was fenced and retired without
cancelling scientific jobs, runners or branches.

Scientific limits remain explicit: three KS1d inertia pilots exhausted six
hours without a count. Failed dependent chains were removed. The c256 equilibrium
is admitted but its eigenpair fails the original residual gate. No ITER stability
sign, lowest branch, marginality or production 3D convergence is established.
Six admitted near-zero probes and four DCON grids retain unresolved physical-map
and normalization differences. The extracted large rest-block fixture remains
private and outside Git, with no admitted completed factor/count/solve. Do not
restart the six-hour campaign because a small-case repair passed. S0 evidence
covers volume/analytic contracts, local continuum and projected rigid subspaces,
not a full global TAE/kink or continuous certificate.

The next controller should collect #47's live result, preserve #39's normal queue,
review #35, and then transfer approved upper layers in order.

## Rules for every milestone

The non-negotiable scientific principles in
[DESIGN.md](DESIGN.md#non-negotiable-scientific-principles) are **promotion
gates, not feature scope**. YAGNI may reduce how much of a capability is built;
it may not remove differentiability/error semantics/structure/provenance and
leave them as unspecified future cleanup.

A capability is promoted only when it has:

1. a precisely named mathematical model and domain;
2. an analytic/manufactured/structural oracle;
3. an independent external or cross-representation oracle when available;
4. primal/JVP/VJP behavior for continuous production operators;
5. bounds/status for reported states, fields, observables and sensitivities;
6. a small correctness case and a representative performance/memory case,
   including extra time and multipliers relative to primal-only and requested
   derivative solves under matched setup/reuse conditions;
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
## Prioritization instead of feature accumulation

Do not implement the catalogue in document order. At each planning cycle rank
candidate vertical slices by:

1. scientific value for an actual research/teaching question;
2. cross-domain reuse of the primitive introduced;
3. how many later capabilities the work unlocks;
4. strength and independence of available validation oracles;
5. expected implementation/maintenance cost and architectural risk;
6. achievable performance and memory envelope;
7. certifiability/error-control leverage.

Prefer a small end-to-end slice that creates a reusable primitive over a broad
framework with no scientific consumer. A capability remains documented but
dormant until it wins this prioritization.

The active spine below is the current best path, not a permanent chronological
contract. Cross-cutting tracks may move earlier when they unlock more valuable
science.


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

The first rung has an [animated particle/Eulerian comparison](docs/streaming-comparison.md)
for supplied Fourier/P1 free-streaming data, with continuous state and density
bounds and generated smooth leaf derivatives. Certified continuous sensitivity
bounds remain unavailable. The next representation gate is the neutral Poisson
deposit/solve/gather comparison below; this does not replace the forced-FP
work on the primary steady spine.

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

## MB0 — microscopic relativistic many-body reference (cross-cutting)

Add a small-box microscopic branch beneath the continuum kinetic hierarchy.
Its endpoint is the multispecies relativistic Maxwell--Lorentz many-particle
system; Vlasov/FK remains the highest continuum mean-field reference. The
[gallery](docs/chapter-gallery.md) supplies the teaching comparisons; the
slab and overdamped planar shielding equations there have different
boundary conditions and assumptions from a periodic microscopic plasma.

Delivered on the shared foundation (species record, one symbolic Coulomb and
Ewald source, general-N KDK with encounter-resolved steps, namelist/NetCDF
provenance): the two-body Rutherford encounter, the periodic pair quench, the
isolated cloud, the periodic OCP plasma oscillation and emergent Debye
shielding (candidates, final evidence pending), and the classical H1a
correlation map. H1a found that Kirkwood superposition fails at
Gamma >~ 10 (32% energy error at Gamma = 10) while HNC stays within 1%.

Remaining ladder, in order:

1. C3/C4 final evidence and \(\Lambda=n\lambda_D^3=100,300,1000,3000,10000\)
   scans, with \(L/\lambda_D\) varied separately and mean-response,
   fluctuation and collision-rate convergence reported separately;
2. self-consistent BGY--Kirkwood and the second YBG equation; \(S(k,\omega)\)
   and dielectric response;
3. microscopic slowing-down, pitch-angle diffusion, momentum/temperature
   equilibration and Green--Kubo transport compared with continuum collision
   operators and the quantum-corrected encounter kernels;
4. weak-to-strong-coupling scans testing Landau, Beliaev--Budker,
   Balescu--Lenard and mean-force closures;
5. a prescribed magnetic-field few-particle case.

The first field-mediated electromagnetic slice is a teaching case: a 1D
Yee/Mur Maxwell field coupled to relativistic Boris macroparticle sheets
(Vlasov--Maxwell PIC, not physical particles) for wave reflection and
transmission at a plasma edge ([em-slab-reflection.md](docs/em-slab-reflection.md),
candidate). Relativistic Coulomb, Darwin, full Maxwell and the tiny direct retarded
Lienard--Wiechert comparison follow when their regime/use case requires them;
they are not prerequisites for the early Coulomb gallery. The classical parent
remains full Maxwell--Lorentz with an explicitly chosen self-field prescription.
The \(\Lambda=1000,L=4\lambda_D\) candidate contains 64000 physical
particles per species; choose its force algorithm only after a measured
small-box gate. It is a planned collective-limit benchmark, not a microscopic
tokamak equivalent or an established cheap direct-pair run.

Keep the periodic screening reference separate from the isolated Yukawa
overlay. Radial averaging of 3D particles does not turn the dynamics into a
one-dimensional line or a spherical-shell model. Spherical/slab reductions
need their own derived force law and model-error status. Full trajectory
certification and statistical estimator intervals are different gates.

Numerically retain direct \(O(N^2)\) pair summation as the gold-standard tiny
oracle, then add FMM/tree or Ewald/spectral-Ewald/P3M for electrostatic/periodic
large-\(N\) interaction. For full electrodynamics prefer local Maxwell-field
evolution for large \(N\); direct retarded Lienard--Wiechert interaction is a
small-\(N\) oracle, with space-time/hierarchical retarded kernels a research
optimization only if measured useful.

Do not reduce particle count and call the result the same microscopic plasma.
PIC/PIF superparticles and finite blobs change finite-\(N\) fluctuations and
close-collision correlations unless every marker is a physical particle or the
modified model is explicitly declared.

The microscopic branch is primarily a closure/correlation validator and
teaching/reference laboratory. It need not delay the steady kinetic spine until
its validation leverage justifies prioritization.

See [microscopic-many-body.md](docs/microscopic-many-body.md).

## QF0 — quantum plasma fundamentals and dense-matter bridge (cross-cutting)

Add a quantum-fundamentals branch that makes the limits of the classical
hierarchy explicit without turning every plasma calculation into a quantum
many-body problem.

Delivered: the certified-discrete periodic pair quench (Q0), the
partial-wave two-electron encounter with singlet/triplet/classical comparison
and a Rydberg packet (Q1/Q1b), Gaussian packet models measured against it
(Q2), finite-N ideal fermions (Q3), the interacting jellium cell with RDM
cumulants (Q4), Hartree/HF/Mermin/TDHF against that cell (Q5), Debye to
Lindhard screening and Wigner to Vlasov (Q7), the relativistic Fermi EOS
(Q8) and hybrid quantum encounters in a classical plasma (QG6). Findings that
set the next steps: Gaussian packets fail for eta >= 1 and classical
ensembles beat them for eta >= 3; open-shell HF breaks translation symmetry;
the N = 3 three-body cumulant reaches 0.24 of the three-body RDM norm at
r_s = 8; the quantum
Coulomb-logarithm correction is ln eta - Re psi(1 + i eta); visible
encounter fringes come from exchange interference.

Active path, in order: closures beyond the first YBG equation, quantum BBGKY with measured cumulants, thermodynamic-limit
references (PIMC/literature), Lindhard with local-field corrections, exchange
in the hybrid collision Monte Carlo, scalable packet/orbital methods only
where validated, and adaptive correlation-order and classical/quantum
switching. Claims, findings and gates are in
[microscopic quantum validation](docs/quantum-microscopic-validation.md);
task state is in [the hierarchy campaign](docs/hierarchy-campaign.md).
Model-level derivations from the microscopic classical and quantum parents
to every model used in research and teaching (for example a NEO-2-type
linearized collision operator traced to quantum mechanics) are planned in
[the theory map](docs/theory-map.md); literature proofs are reproduced, not
relied on.
Stellar structure remains outside this sequence.

Physics ladder:

1. finite-temperature ideal Fermi gas and regime coordinates
   \(\Theta=k_BT/E_F\), \(\Gamma\) and \(x_F=p_F/(m_ec)\);
2. classical Debye \(\rightarrow\) finite-\(T\) Fermi/Thomas--Fermi screening;
3. Lindhard/RPA response and Friedel oscillations;
4. quantum hydrodynamic/Fermi-EOS closures where derived;
5. Wigner--Poisson and later Wigner--Maxwell with
   \(\hbar\rightarrow0\) recovery of Vlasov;
6. Hartree/Hartree--Fock and DFT/TDDFT-like microscopic reference providers;
7. quantum BBGKY/reduced-density-matrix and cumulant/low-rank correlation
   models;
8. quantum collision/dielectric corrections: diffraction, exchange, Pauli
   blocking, Uehling--Uhlenbeck/quantum Lenard--Balescu-like closures;
9. equilibrium exact-diagonalization/QMC/PIMC reference cases where practical;
10. stellar/dense-plasma closure examples from classical ideal gas through
    degenerate and relativistically degenerate electron pressure.

Items 1--3, 5, the Hartree/HF part of 6 and the ideal-EOS part of 10 exist
as teaching references; items 7--9 are on the active path above.

Numerical representation remains independent:

- exact \(\Psi\) or \(\rho_N\): tiny-\(N\) gold-standard references;
- Hartree/HF/DFT/TDDFT: ordinary 3D orbital/density discretizations using
  FEM/spectral/plane-wave/finite-difference families;
- QMC/PIMC: stochastic microscopic equilibrium reference;
- Wigner: Eulerian/semi-Lagrangian/DG/spectral/tensor phase-space methods;
- quantum-fluid/Fermi-EOS closures: ordinary macroscopic FEM/FV/DG.

Planning scales are method dependent, not support claims: exact full
wavefunctions are generally only a few 3D particles; QMC/PIMC microscopic
boxes are typically tens--thousands; orbital methods may reach
\(10^2\)--\(10^4\)-scale particle/orbital cells; Wigner inherits the
high-dimensional FK problem; quantum-fluid/EOS closures can cover whole-device
or whole-star domains.

The fundamental native boundary is **QED**. Pair creation/annihilation,
quantized-radiation dynamics, strong-field QED cascades and radiative QED
corrections remain external provider territory. Relativistic Fermi EOS and
controlled semirelativistic/Dirac-Coulomb references remain admissible below
that boundary when QED effects are negligible.

See [quantum-plasma.md](docs/quantum-plasma.md).

## AM0 — adaptive multiscale, multifidelity and UQ controller (cross-cutting)

The long-range controller chooses the cheapest **local-in-space-and-time**
combination of numerical and physical fidelity that closes the requested
observable error budget.

The candidate refinement set eventually includes

\[
\{h,p,\Delta t,\text{velocity/harmonic order},N_p,
\text{solver/precision},
\text{MHD/GK/FK},
\text{collision/correlation order},
\text{classical/quantum fidelity},
\text{neutral/radiation/chemistry fidelity},
\text{microprovider sample},
\text{surrogate refresh}\}.
\]

A rich model \(R_1\) evaluated on a lifted cheap-model state \(U_0\) supplies a
model defect \(d=R_1(LU_0)\). Adjoint weighting estimates its effect on the
requested observable,

\[
\eta\simeq-\lambda^\dagger d.
\]

Analytical asymptotic remainder bounds, validated model-defect estimates and
sparse richer-model computations are preferred over inventing a random model
error. Probabilistic discrepancy is used only with an explicit stochastic
interpretation.

Adaptation is not merely spatial AMR. It may be local in space, time,
phase-space, species, energy, correlation order, frequency/angle, reaction
state and uncertainty/sample space.

Examples include:

- MHD globally with a transient local GK/FK region;
- \(f_1+g_2\) globally with \(g_3\) activated only where its observable defect
  matters;
- classical Maxwell--Lorentz dynamics with short quantum-diffraction
  microproblems for rare close encounters;
- fluid neutrals with a kinetic rarefied boundary layer;
- grey radiation with transient multigroup/angular refinement;
- sparse expensive turbulence/FK/N-body/Wigner/QMC microproblems requested by
  macro adjoint indicators.

Quantum pair promotion is appropriate for local diffraction/scattering only
when its assumptions hold. Degeneracy/exchange generally requires promotion of
a region/species/many-body model, not just one close pair. Demotion requires an
explicit information-loss/decoherence/coarse-graining rule.

Micro/macro coupling uses lift/projection maps, conservative exchange,
interface-defect estimates and local temporal subcycling/multirate integration.
If locality or scale separation fails, enlarge the rich region or promote the
global target model.

The controller's conceptual greedy score is

\[
k^\star=\arg\max_k \frac{\eta_k}{C_k},
\]

where \(\eta_k\) is expected observable-error reduction and \(C_k\) incremental
cost. Recompute after each significant promotion because refinement effects are
not independent.

UQ is multifidelity as well: cheap levels receive many samples and expensive
corrections few samples using control-variate/MLMC/multifidelity patterns, while
deterministic error bounds remain distinct from statistical variance.

Development sequence:

1. \(h/p/\Delta t\) observable-driven refinement;
2. one cheap/rich physical-model defect with a validated remainder;
3. one sparse microprovider invoked by an adjoint indicator;
4. one local-in-time promotion/demotion case;
5. one spatial model patch with conservative transfer;
6. one \(g_2\to g_3\) correlation-order indicator;
7. one classical close encounter promoted to a quantum microproblem;
8. one multifidelity UQ/control-variate example;
9. only then extract a reusable generic scheduler.

See [adaptive-multiscale-uq.md](docs/adaptive-multiscale-uq.md).

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
6. weak/current-sheet ideal MHD with the function-space contract of
   [equilibria-inference.md](docs/equilibria-inference.md#solution-space-contract)
   (broken/interface FEEC, B in H(div), interface-aware J = curl B);
7. native MRxMHD, with SPEC/SPECTRE as comparators;
8. a native nested-surface fast path: a spectral/surface representation of the
   same physics residual, not only VMEC/DESC/GVEC imports;
9. flowing/rotating stationary MHD equilibria;
10. anisotropic/tensor-pressure equilibria (CGL, then kinetic-MHD; ANIMEC as
    comparator);
11. explicit comparison with HINT/PIES-like nonnested full-pressure relaxation
    and MRX/Zhang--Zhu--Hegna-type non-nested solvers, within matched models;
12. resistive-MHD and Hall/XMHD stationary/evolution paths where reconnection and
    history are part of the selected physics;
13. kinetic-MHD and kinetic continuation where pressure/current closure is no
    longer adequately fluid;
14. ideal/resistive linear response, eigenmodes and nonlinear stability through
    the common S0 consumer;
15. conducting structures, vacuum/BEM, coils and external circuits.

Milestone: general 3D free-boundary physical-space equilibrium with
islands/stochastic regions.

Ambition: the broadest framework, with near-best-in-class performance per
regime through specialised representations, not one discretisation for
everything. MRX is the main performance target for the non-nested branch.
Dispatch picks the cheapest valid representation:

| Physics/topology | Representation |
|---|---|
| smooth nested | surface/spectral fast path |
| unusual nested geometry | mapped surface |
| smooth islands / non-nested ideal | physical-space FEEC |
| current sheets | broken/interface FEEC |
| chaotic relaxed regions | MRxMHD |
| reconnection | resistive |
| Hall | Hall/XMHD |
| anisotropic pressure | CGL/kinetic-MHD |
| strong flow | flowing MHD |
| kinetic closure | DK/GK/FK |
| free boundary | plasma + BEM vacuum + coils/circuits |

For static smooth ideal MHD, retain the actual constraint
\(B\cdot\nabla p=0\): removing the nested-surface representation does not
license arbitrary smooth pressure gradients in ergodic volumes. Different
stationary solutions are stored as different equilibria with their own initial
state/continuation provenance; no general probabilistic branch machinery belongs
in KIN6D.

The nested-surface 3D equilibrium lane (energy minimisation over a map from a
reference domain) keeps its reference chart behind an abstract interface: the
energy, quadrature, generated leaves and solver see only the map Jacobian. The
first chart is cylindrical (TC24 RMP). A later step adds a curve-following
frame, X = Gamma(zeta) + q1 N(zeta) + q2 B(zeta) with a closed axis curve and a
Frenet or user-defined frame, as GVEC's G-Frame does, so that non-planar-axis
and Figure-8 stellarators can be solved. Its gate is the atlas case
`g_frame_figure8_equilibrium`: first an analytic curve-frame case, then a
frozen public Figure-8 equilibrium compared with GVEC on the same frame and
boundary. Not scheduled before the TC24 lanes close.

ARES/MEPHIT/MARS/GPEC and later JOREK/NIMROD-type cases remain response and
evolution landmarks. VMEC/DESC/GVEC, MRX, HINT/PIES and SPEC/SPECTRE are equilibrium
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

## Coverage tracks — selected by priority, not by domain

The remaining classical-plasma envelope is maintained as reusable capability
tracks rather than a sequence of domain-specific solver milestones. Pull one
forward when an application makes it the highest-value next vertical slice.

### Reacting and partially ionized matter

Add the shared low-temperature/reacting hierarchy:

\[
\text{local-field/electron-energy}
\rightarrow
\text{Boltzmann/EEDF}
\rightarrow
\text{FK electrons},
\]

together with heavy-species reaction-diffusion, compressible multispecies flow
and rarefied kinetic/DSMC-like neutrals. Large sparse internal-state chemistry,
ion-neutral drag, charge exchange and multi-temperature energy exchange use the
same species/reaction graph as fusion edge and space applications.

Drift-diffusion/Poisson and front-capable models supply the cheap end for
ionization-front/streamer-like regimes. See
[reacting-partially-ionized.md](docs/reacting-partially-ionized.md).

### General surfaces, dielectrics and adjacent media

Generalize the existing sheath/wall/circuit contract to conducting electrodes,
dielectrics with evolving surface charge, secondary/thermionic/field/
photo-emission, reactive coverages, deposition/etching and finite conducting
bodies. Detailed solid/liquid physics remains an external provider behind a
conservative particle/current/heat/reaction interface.

### Radiation hierarchy

Make radiation native through the useful reduced hierarchy

\[
\text{optically thin}
\rightarrow
\text{opacity/escape closure}
\rightarrow
\text{grey diffusion}
\rightarrow
\text{multigroup/moment transport}
\rightarrow
\text{angular transport},
\]

with high-fidelity line-by-line angle/frequency non-LTE transfer available
through an external provider until a concrete case justifies internalizing more.
See [radiation.md](docs/radiation.md).

### Compressible, shock/front and classical HEDP regimes

Reuse the fluid/MHD/reacting/radiation infrastructure for strongly compressible
and shock-dominated plasma. Add nonideal EOS/opacity as versioned closure data,
laser/beam deposition and stiff matter-radiation coupling only when needed.

This track covers classical portions of arcs, propulsion, atmospheric entry,
Z-pinch/magnetized-target/ICF-like and laser-produced plasmas without creating a
separate HEDP architecture. Quantum dense-matter capabilities follow QF0
and require their own validation; QED remains an external-provider boundary.

### Microscopic many-body and strong-correlation reference

Use the MB0 Maxwell--Lorentz hierarchy for both weakly coupled closure
validation and strongly coupled reference problems. Explicit BBGKY/cumulant
models, low-rank pair correlations and mean-force closures sit between the
finite-(N) microscopic system and continuum FK. Quantum-many-body numerical
ideas may inspire compression, but any reduced correlation representation
carries an explicit truncation defect.

### Wave/statistical and energetic-particle transport

Retain the planned focused energetic-particle hierarchy,
wave-kinetic/statistical-spectrum evolution and global-orbit/turbulence
cross-spectrum coupling. These are shared fusion/space/accelerator primitives,
not domain-specific additions.

### Topology, embedded bodies and rotating/stratified fluids

Retain magnetic-topology diagnostics, charged grains/finite bodies, gravity,
rotation, buoyancy and anisotropic conduction as common services. Reconnection,
shocks, edge transitions, non-neutral traps and beam/wakefield applications
should be expressed through existing model/analysis consumers rather than new
top-level solver families.

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
The long-range form is local and temporal rather than one global model choice:
the final result records the space-time/phase-space model map and the dominant
remaining error/UQ contributions. See
[adaptive-multiscale-uq.md](docs/adaptive-multiscale-uq.md).

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

### Method diversity and adaptivity

The common engine must admit the numerical family appropriate to the selected
regime. In addition to the primary tensor/FEEC/spectral path, retain deliberate
routes for conservative finite-volume/entropy-stable DG and WENO-like
shock/front methods, block AMR/hp adaptivity, semi-Lagrangian kinetic methods,
PIC/PIF/blobs, DSMC/MCC, photon Monte Carlo and small direct N-body/MD
references.

See [numerical-methods.md](docs/numerical-methods.md). Do not implement a method
until a model requires it, but do not let today's preferred discretization
become an architectural prohibition.

### Time-scale and asymptotic numerics

Add IMEX, exponential, splitting, multirate/subcycling and
asymptotic-preserving/uniformly accurate formulations when scale separation
makes them materially cheaper than resolving irrelevant fast physics. Verify
both the parent equation and the intended limiting model.

### Reduced-order and uncertainty acceleration

Low-rank/tensor compression, sparse grids, ROMs and surrogate closures are
optional accelerators with explicit reconstruction/error status. MC/QMC and
stochastic-collocation/ensemble tools remain uncertainty or sampling methods,
not substitutes for deterministic certificates.


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

### Build structure and fast development loop

Decision 2026-10-15 (Chris + Claude). KIN6D stays one repository and one CMake/CTest
project, but stops being one monolithic `kin6d` library: today about 106 modules,
145 executables and about 200 unlabelled tests. Separate repositories may follow
later along the same seams. The build stays on CMake; no conversion to fpm.

Vertical layers (a layer uses only layers below it):

| Layer | Target | Content |
|---|---|---|
| L0 | `kin6d_core` | intervals, ODE, species, sort keys, series statistics, laboratory, identity helpers |
| L1 | `kin6d_fem2d` | FortFEM slice, triangle/curved bases, CSR and direct solvers, GS options, timers, load batches |
| L1 | `kin6d_pairs` | Ewald, near encounter, PME, N-body, pair laws, radial DVR, quantum collision kernels |
| L2 | model classes | `kin6d_equilibrium` (GS, equilibrium, 3D coils, equilibrated bounds); `kin6d_mhd` (MHD, ideal, response, vacuum surfaces, nonnested; may use equilibrium); `kin6d_classical`; `kin6d_quantum`; `kin6d_continuum` |
| L3 | apps | applications, demos, animation export; only L3 links plotting and graphics IO |

Rules:

- `kin6d` remains an umbrella INTERFACE target, so existing consumers and tests
  keep working.
- A layering CTest checks the Fortran `use` graph against the table; an upward
  `use` fails it.
- Every test carries CTest `LABELS`: its model class (or layer) plus a tier
  `fast`, `slow`, `regression` or `performance`.
- Each class has a `kin6d_check_<class>` build target that builds its tests'
  executables, and a CTest preset per class.
- A `dev-fast` configure preset turns plotting and IO extras off and shares the
  FetchContent cache across worktrees.

Migration:

1. CMake-only partition: assign the existing sources to layer targets, add labels,
   per-class targets, presets and the layering test. No files move. This starts
   after the open estimator branches are merged.
   Implemented 2026-10-15 (`cmake/kin6d_layers.cmake`, `kin6d_layering` test,
   `cmake/kin6d_test_labels.cmake`, `dev-fast` and per-class test presets). The
   use graph required `interaction_laboratory` in `kin6d_pairs` (the quantum
   encounter uses it), Wigner-Poisson/dispersion in `kin6d_continuum` (shared
   streaming leaves) and `equilibrium_nonnested` in `kin6d_mhd`; `kin6d_io`
   (result IO, animation export) stays an L3 library because only applications
   use it. Tiers are `fast`/`slow` (over 60 s). Executables link the layers
   their sources use instead of the umbrella. The dev-fast preset shares the
   MakeAvailable checkouts (fortplot, fortio, fortsparse) by pin; fortnum and
   fortfem stay per build tree.
2. Implemented: handwritten library sources live in `src/<layer>/`, matching
   their target ownership. Generated leaves stay in `src/generated/`, and
   generators retain their canonical output paths.

Development loop (to go into AGENTS.md with step 1): edit, then
`cmake --build build --target kin6d_check_<class>` and `ctest -L <class> -LE slow`;
never run the full suite, including at controller integration. With `fo`, the same
loop is `fo test -L <class> -LE slow`, building only the selected tests' targets, and a Gremlin lane
scoped by label. That `fo` support is implemented in `fo` as a
backend-neutral project model (fo doc/CMAKE_GREMLIN_PLAN.md, 2026-10-15) that reads
targets from the CMake file API and labels from CTest; KIN6D needs no fo-specific
files. Baseline for the speed-up: an incremental rebuild after touching one source
takes 2-6 s, a single target 2.9 s, a clean build about 800 CPU-s; single tests run
up to 1600 s.

### Literature and validation

Every new rung enters the analytical benchmark atlas and has a frozen
cross-code/literature reference where possible. Existing itpplasma codes
NEO-2, NEO-RT, SIMPLE, GORILLA, KAMEL/KiLCA/KIM, MEPHIT, ARES and rabe are
valuable retained oracles as their production roles are subsumed.

## Deliberate outer boundary

The fundamental quantum boundary is **QED**, not "all quantum matter".
KIN6D may natively include quantum degeneracy, exchange/diffraction,
Wigner/quantum kinetic models, quantum BBGKY/cumulants, orbital/QMC
microscopic references and quantum-fluid/EOS closures. It may consume
QED-derived rates as data.

KIN6D does not initially own general quantized-radiation dynamics, pair
creation/annihilation, strong-field QED cascades or radiative QED corrections.
It likewise does not own general-relativistic plasma/gravity, neutron-star/
nuclear-density matter, full reactor neutronics/nuclear engineering, detailed
bulk solid mechanics/material microstructure, comprehensive liquid/solid
chemistry, or the highest-fidelity line-by-line angle/frequency non-LTE
radiation-transfer discipline. Reduced radiation hydrodynamics through
grey/multigroup/moment models remains inside scope.

Chris&AI
