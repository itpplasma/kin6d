# kin6d agent instructions

## Read before changing the project

Read README.md, DESIGN.md and ROADMAP.md first. Then read the specification
that owns the proposed change. In particular:

- docs/representations.md for Eulerian/semi-Lagrangian/particle/blob,
  full-f/delta-f, PIC/PIF and deterministic/QMC/MC choices;
- docs/student-laboratory.md for early priorities;
- docs/stability.md before adding eigenmodes, forced response, nonmodal growth,
  continuation, bifurcation or nonlinear stability;
- docs/transport-multiscale.md before adding macro transport, closure providers,
  integrated modelling or multiscale coupling;
- docs/model-family.md and docs/species-models.md for MHD/reduced/hybrid models;
- docs/moment-models.md for gyrofluid/Landau-fluid/gyromoments;
- docs/collisions.md for collision operators and thermodynamic structure;
- docs/numerics.md and docs/performance.md before hot-path work;
- docs/differentiability.md and docs/adjoints-structure.md before derivative,
  sensitivity, optimization or goal-error work;
- docs/certification.md before a rigorous-result claim;
- docs/tooling.md before implementation or dependency changes.

The first prescribed-field 1X1V slice already runs. Read
docs/first-example.md before changing it or claiming broader evidence.

## Architectural rules

- Keep **one mathematical engine**. Do not add named-code modes such as
  neo2_mode, gene_mode, jorek_mode, pic_code or separate application
  implementations for physics already expressible through shared operators.
- Physics model, distribution form, numerical carrier, particle-field coupling,
  sampling strategy, collision realization, field model and active/test/
  reservoir role are separate axes. Reject unsupported combinations explicitly.
- Direct steady kinetic science is the primary research spine. Early transients
  and turbulence are required demonstrations and verification branches, not the
  architectural center.
- Stability/response is an analysis-consumer axis available to every compatible
  MHD/fluid/gyrofluid/DK/GK/FK model. Do not create separate stability physics.
- Macro transport dimensionality (0D/1D/1.5D/2D/3D) and closure fidelity are
  independent axes. A transport closure may return coefficients, fluxes,
  nonlocal kernels or spectral/cross-spectral data; do not force every closure
  into scalar D/chi/V.
- Structure preservation is the first principle. Semi-Lagrangian,
  characteristic, Eulerian, PIC/PIF and blob methods are realizations selected
  by the problem, not project identities.
- Deterministic continuum is the preferred high-accuracy reference where
  competitive, but Monte Carlo, randomized QMC, real PIC and stochastic
  collisions/neutrals are valid first-class methods when justified. Never call
  a statistical interval a deterministic certificate.
- Full-f, delta-f and control-variate formulations are first-class. A reduced
  model may serve as a preconditioner or control variate without changing the
  target model solved.
- Physical vectors are Cartesian. Charts and momentum frames are independent.
  Do not require global toroidal angle, nested surfaces or one geometry class.
- Bulk MHD aggregates explicitly owned populations. It is not a per-species
  fluid switch. Prevent double counting of charge/current/pressure/energy.

## Mathematics, derivatives and certification

- Project-owned continuous numerical leaves come from retained FortSym
  definitions. Generate primal/JVP/VJP from the same symbolic source; add HVP
  only for a concrete consumer. Keep symbolic mathematics, indexed operator
  structure and hardware schedules separate.
- FortAD, finite/complex-step checks and analytical identities are derivative
  oracles. Production correctness must not require LLVM/Enzyme.
- Differentiate implicit equations, not solver iterations: nonlinear states,
  linear solves, roots and eigenproblems use their equation-level tangent and
  adjoint systems.
- Declare physical pairings, signs, nullspaces and time-reversal parity.
  Self-adjointness, skew-adjointness, positivity, entropy production and
  Onsager--Casimir relations are explicit contracts, not informal expectations.
- The fast candidate solver is untrusted. The primary rigorous path is a
  Fortran interval/ball certificate. Lean is optional downstream assurance and
  is not an early milestone requirement.
- Distinguish candidate, certified-discrete, certified-continuous,
  optionally formally-certified-continuous, probabilistically-certified
  and certification-failed. State the exact claim.
- A numerical enclosure proves the named model, not physical reality. Model
  reduction, physical-data uncertainty and stochastic sampling need separate
  bounds/status.
- Goal-oriented indicators are not rigorous remainders until the stability,
  lifting and nonlinear remainder hypotheses are checked.

## Numerical/performance rules

- Review mathematical complexity, separability, unavoidable bytes and solver
  vector count before implementation.
- Prefer tensor-product mapped patches, sum factorization, partial assembly,
  matrix-free actions and many-RHS batching where the model permits them.
- Never materialize a global kinetic Jacobian or logical pairwise matrix merely
  because the formula contains two indices. Use tiled map/reduce or sparse
  interaction structure.
- Hot production kernels allocate no memory per call. Avoid hot-loop dynamic
  dispatch, hidden array temporaries and duplicate full-state layouts.
- Keep one canonical state layout plus bounded tile scratch. Adaptivity should
  use block/resolution classes rather than arbitrary per-cell raggedness unless
  evidence demands otherwise.
- Configuration-space ownership is the default distributed decomposition so
  local collision/velocity work remains local. Change only on measured need.
- Candidate mixed precision, fast math and nondeterministic GPU reductions are
  allowed only when strict residual/certificate evaluation remains independent.
- Maintain simple reference, optimized CPU and optimized accelerator actions
  for hot operators when those paths exist. Verify them against one another and
  independent physics/mathematical oracles.
- Performance evidence is multi-level: microkernel, operator, residual, solve
  and scientific case. Track time, peak/attributed memory, iterations,
  communication/reductions and derivative overhead. Compare against best-known
  fixed-input baselines, not only the previous commit.
- Fix generic bottlenecks in the dependency that owns them. Do not duplicate a
  generic optimized kernel downstream.

## Physics admission rules

A new production physics feature needs:

1. a concrete scientific or teaching use case;
2. a complete mathematical model and coupling/energy/source ledger;
3. an analytic/manufactured or structural oracle;
4. an independent cross-code/literature oracle where available;
5. a differentiability statement;
6. a certification target/status;
7. a performance/memory target.

Implement the smallest complete case. Do not build speculative framework
hierarchies in advance.

Important long-range branches include collisions, gyrofluid/gyromoments,
stability/bifurcation, 0D--3D integrated transport and true multiscale
micro/macro coupling, MHD/kinetic-MHD, PIC/PIF/blobs, real-device DK/GK/FK,
edge/neutrals/walls, RF/NBI/reactions, space plasmas, diagnostics and
optimization. Planned scope is
not implemented functionality.

## Repository and dependency policy

- Scientific/application code, generators and rigorous checker are Fortran.
  No mandatory Python layer.
- FortFEM owns generic spatial mathematics; FortNum generic numerics and
  verified arithmetic; FortSparse sparse kernels; FortSym symbolic generation;
  FortAD derivative tooling; FortIO/FortPlot/FortArray their respective generic
  domains. KIN6D owns plasma-specific mathematics and composition.
- Native CUDA is optional behind ISO_C_BINDING and a thin C ABI. The CPU path
  remains a reference.
- When an owned upstream dependency has a concrete defect or missing generic
  feature, add a focused independent reproducer there, repair it there, verify
  it, push the owning main as allowed, pin/recheck KIN6D. Do not keep permanent
  consumer workarounds.
- Preserve unrelated edits. One controller owns integration/promotion.
- Never touch protected faepop* or faepcr* hosts without an explicit user
  request naming the host.
- Keep generated images, movies, PDFs and large benchmark data out of source
  history. Retain recipes, numeric evidence and provenance.

## Communication

Scientific and technical communication is English. Messages drafted or sent on
Chris's behalf to other people end with Chris&AI on its own final line unless
explicitly instructed otherwise.

Chris&AI
