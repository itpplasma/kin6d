# kin6d

KIN6D is a **general plasma multiphysics engine**. Its deepest classical
reference is the multispecies **relativistic Maxwell--Lorentz many-particle
system**; relativistic 3X3P Vlasov--Maxwell/FK is the highest continuum
mean-field reference. A quantum-fundamentals branch covers degeneracy,
exchange/diffraction, quantum kinetic/statistical reductions and dense-plasma
reference models up to, but not including, general QED. Its practical purpose is broader: assemble the
cheapest mathematically controlled plasma model and numerical representation
that answer the requested question at the requested accuracy.

The intended application envelope spans fusion, Solar-System/natural plasmas,
industrial and low-temperature discharges, beam/accelerator and laser plasmas,
atmospheric plasmas, dusty/non-neutral systems, and classical high-energy-density
or strongly coupled regimes where an appropriate classical/nonideal model
remains valid.

The project is built around one common physical and numerical architecture.
Equilibrium, stability/response, MHD, gyrofluid, gyromoments, drift/gyro/full
kinetics, macro transport, continuum and particle representations, collisions,
turbulence, edge physics, neutrals, walls, heating, space plasmas, diagnostics
and optimization are configurations, consumers or reductions of shared
operators rather than separate named-code modes.

## Non-negotiable principles

KIN6D's broad capability catalogue is optional scope. The following are not:
**one mathematical source of truth, differentiability, controlled/certified
error, structure preservation, independent verification, reproducibility and
performance at controlled accuracy**.

Every production vertical slice must preserve these contracts at the level
appropriate to that slice. In particular, differentiability and error control
are not later features to bolt on after the physics works. Smooth production
operators are designed with primal/JVP/VJP semantics from one definition, and
numerical answers carry explicit error/status information with numerical,
model, input/reconstruction and stochastic contributions kept distinct.

See the canonical statement in [DESIGN.md](DESIGN.md#non-negotiable-scientific-principles).

## Current status

The CPU laboratory includes prescribed-field free streaming and Krook balance,
a single-mode Poisson solve, a self-consistent source-balanced 1X1V
Poisson--Fokker--Planck family, and a bounded homogeneous OU/FP steady solve.
It retains FortSym-generated primal/JVP/VJP leaves, Fortran plots and reviewed
continuous bounds for the declared cases. See
[the first example](docs/first-example.md) and
[the steady laboratory](docs/steady-laboratory.md) for their exact scope.

The coupled solver resolves a specified Fourier/polynomial source family;
certification of arbitrary sources beyond the damped-transport model and
nonlinear branch uniqueness remain future work. PIC/PIF, the collision hierarchy, gyrofluid
turbulence and real-device DK are also planned capabilities.

Upstream CPU benchmark infrastructure and a general damped-transport
residual/inverse certificate are implemented. The latter accepts independently
supplied Fourier/P1 sources and certifies full-domain state and directional
derivative errors without an exact solution; see
[general transport](docs/general-transport.md). Forced-FP flux/stability bounds
are next. Simple continuum/particle/PIF/PIC comparisons enter early, before hard
nonlinear or particle cases. The practical default targets about 1% accuracy
for a named observable with a declared scale; tighter accuracy remains a
request rather than a prerequisite for useful results.

The [animated representation comparison](docs/streaming-comparison.md) shows
weighted characteristic particles alongside an Eulerian Fourier/P1 solution
of prescribed-field free streaming, with certified distribution and density
bands. It is the first transport comparison; self-consistent particle-field
coupling remains planned.

## Scientific spine

The main research progression is deliberately steady-state first:

\[
1X1V\ {\rm steady}
\rightarrow
0X3V\ {\rm collisions}
\rightarrow
1X3V\ {\rm transport}
\rightarrow
{\rm real\ device\ DK}
\rightarrow
{\rm gyromoments/GK}
\rightarrow
{\rm global\ FK}.
\]

Tiny transient and turbulence examples are developed early in parallel because
they are excellent teaching, verification, structure-preservation and
performance tests. Hasegawa--Mima, Hasegawa--Wakatani, interchange/blob and
gyrofluid examples must use the same operator machinery as the steady spine.

Two additional axes are first-class from the beginning. **Analysis type**
covers linear spectra, forced response, nonmodal growth, continuation,
bifurcation and nonlinear stability for every physical model that supports
them. **Transport** spans 0D/1D/1.5D/2D/3D macro evolution with closure
providers ranging from prescribed coefficients through neoclassical and
quasilinear response to nonlocal orbit kernels, nonlinear microturbulence and
fully coupled global evolution. See [stability](docs/stability.md) and
[transport/multiscale coupling](docs/transport-multiscale.md).

A third cross-cutting axis is **adaptive fidelity and uncertainty**. Eventually
KIN6D should refine or coarsen locally in space, time, phase space, species,
correlation order and physical model, spending expensive FK, microscopic,
quantum, radiation or detailed-chemistry work only where it materially changes
the requested observable. Analytical asymptotic remainders, rich-model defects,
adjoints, sparse enriched solves and UQ all contribute to one explicit error
ledger. See [adaptive multiscale and uncertainty](docs/adaptive-multiscale-uq.md).

## Model hierarchy

The intended charged-plasma hierarchy is

\[
{\rm MHD}
\rightarrow
{\rm Hall/XMHD,\ multi-fluid}
\rightarrow
{\rm gyrofluid/Landau\ fluid}
\rightarrow
{\rm gyromoments}
\rightarrow
{\rm DK/GK}
\rightarrow
{\rm FK}.
\]

Species may use different models in one calculation. Bulk MHD is an aggregate
model with explicit population ownership, not merely another per-species flag.
Spatially hybrid models are allowed only when their coupling and interface
balances are derived.

The field hierarchy is independent:

\[
{\rm prescribed}
\rightarrow
{\rm electrostatic}
\rightarrow
{\rm quasineutral\ electromagnetic}
\rightarrow
{\rm Darwin}
\rightarrow
{\rm full\ Maxwell}.
\]

Background, externally imposed drive and induced plasma response are separate
components.

## Representation is a separate axis

KIN6D does not identify kinetic physics with one numerical representation.
For a distribution, select independently:

- **form:** full-f, delta-f, or a declared control-variate split;
- **carrier:** Eulerian basis/grid, semi-Lagrangian continuum, particles,
  finite blobs or a hybrid;
- **particle-field coupling:** PIC, PIF, particle-in-FE or general
  particle-in-basis;
- **sampling:** deterministic quadrature, QMC/randomized QMC or Monte Carlo;
- **collision realization:** deterministic operator or justified stochastic
  process.

Particles are not synonymous with PIC, and PIC is not synonymous with Monte
Carlo. Deterministic continuum remains the preferred high-accuracy reference
where it is competitive, but the project principle is **no uncontrolled
numerical noise**, not a ban on stochastic methods.

See [docs/representations.md](docs/representations.md).

## Foundational structure

Where the physics permits it, reversible dynamics is written in Hamiltonian or
Poisson form and irreversible dynamics in metric/metriplectic form,

\[
\dot U =
J(U)\frac{\delta H}{\delta U}
+
M(U)\frac{\delta S}{\delta U}
+
P_{\rm open}(U),
\]

with the appropriate discrete antisymmetry, positivity, degeneracy,
conservation and source/sink ledgers. Semi-Lagrangian/characteristic methods are
important realizations of transport, but **structure preservation is the first
principle**, not one particular transport algorithm.

Stationary states are primary objects:

\[
R(U;\theta)=0.
\]

Time evolution uses the same residual/brackets rather than defining a separate
physics implementation. An equilibrium is therefore one stationary solution of
a named physical model, not a privileged file format and not necessarily a
unique function of a small parameter set. Imported EQDSK, VMEC, SPECTRE or
other states are candidates, initial guesses and comparison data. KIN6D must be
able to evaluate their residuals, recompute them with a native model where that
model exists, and certify the resulting numerical claim. Nested flux surfaces
are an optional exploitable structure, never a core requirement. See
[equilibria and inference](docs/equilibria-inference.md).

## Geometry and spatial numerics

Physical vectors are Cartesian. Computational coordinates are arbitrary smooth
maps \(x=X(\xi)\); no core model may require a global toroidal angle or nested
flux surfaces. Cylindrical/Fourier geometry is an optimized fusion chart.
Curve/G-frame maps, spherical/open domains, multipatch geometries, X-points and
the Figure-8 stellarator are part of the intended geometry envelope.

FortFEM owns generic spatial mathematics: FEEC, tensor-product B-splines/NURBS,
DG/HDG, traces, mortars, FCI and boundary operators. KIN6D owns plasma laws.
Mapped tensor patches and matrix-free tensor actions are preferred where they
preserve both performance and geometry; unstructured/cut methods are escape
paths, not the default everywhere.

## Collisions and thermodynamic structure

The collision hierarchy includes none, Krook/BGK, Lorentz, conservative
LBD/Dougherty, reduced linearized FP, Landau test/field-particle pieces, full
linearized Landau, nonlinear Landau, relativistic Beliaev--Budker and later
screened/matched research models.

Collision physics, discrete invariant enforcement and certificate strength are
separate contracts. Linearized Coulomb models should expose the correct
weighted adjoint, invariant nullspace, entropy-production form and
Onsager--Casimir structure. NEO-2-style independently computed reciprocal
transport coefficients remain permanent end-to-end verification tests.

## Physics envelope

Long-range scope includes:

- imported and internal axisymmetric/3D equilibria;
- ideal/resistive/extended MHD, kinetic-MHD and response/eigenvalue problems;
- gyrofluid, Landau-fluid and Hermite--Laguerre gyromoment hierarchies;
- DK, GK and relativistic FK;
- PIC, PIF, variational particle-FE, blobs, remapping and delta-f methods;
- deterministic and stochastic Coulomb/reaction/neutral realizations;
- core, pedestal, SOL, divertor, sheath and wall physics;
- atomic/molecular/internal-state chemistry, partial ionization and reacting
  neutral gas from fluid to rarefied kinetic descriptions;
- drift-diffusion/front models, dielectric/electrode charging and reactive
  surface/interface physics;
- native radiation from optically thin losses through grey/multigroup/moment
  transport, with high-fidelity non-LTE spectral transfer available as an
  external provider;
- RF/quasilinear/full-wave and wave-kinetic physics, NBI and energetic/fusion
  product sources;
- compressible/shock-capable classical plasma and nonideal EOS/opacity
  closures;
- a full relativistic microscopic Maxwell--Lorentz many-body reference,
  including direct pair, FMM/Ewald/P3M and small retarded Lienard--Wiechert
  validation paths for correlations/collision closures;
- quantum plasma fundamentals: finite-temperature Fermi statistics,
  Thomas--Fermi/Lindhard screening, Wigner kinetics, quantum
  BBGKY/cumulants, quantum collision corrections, orbital/QMC reference
  methods and quantum-fluid/EOS closures for dense/stellar plasmas;
- coils, vacuum, conducting walls, dielectrics and external circuits;
- Solar-System, atmospheric, industrial and other classical-plasma
  configurations;
- synthetic diagnostic forward operators, adjoints, optimization and control;
- integration with TIAGO for equilibrium/state reconstruction and inference;
- local-in-space-and-time adaptive fidelity across numerical resolution,
  MHD/GK/FK/microscopic/quantum/correlation/radiation/reaction levels, with
  propagated deterministic, epistemic and stochastic uncertainty.

The complete intended physics/numerics envelope is catalogued separately from
implementation priority. See [the capability catalogue](docs/capability-catalogue.md)
and [the numerical-method catalogue](docs/numerical-methods.md). The fundamental
quantum boundary is **QED**: KIN6D may own nonrelativistic/semirelativistic
quantum plasma, degeneracy, exchange/correlation, Wigner/quantum kinetic and
quantum-matter reference methods, but pair creation, quantized-radiation
dynamics, strong-field QED and radiative QED corrections remain external.
General relativity, full reactor neutronics/nuclear engineering, detailed bulk
material mechanics and the highest-fidelity line-by-line non-LTE
radiation-transfer discipline likewise remain external coupled systems unless a
concrete reusable interface is justified.

## Differentiability

Continuous production operators are defined once symbolically where practical.
FortSym is the mathematical source of project-owned numerical leaves and emits
primal/JVP/VJP families, HVPs when needed, strict/reference arithmetic,
optimized CPU/CUDA leaves and rigorous enclosure leaves where supported.

Global nonlinear, linear, root and eigen solves are differentiated
**implicitly**. Do not differentiate through Newton/Krylov iterations when an
equation-level derivative exists. FortAD and finite/complex-step checks are
independent derivative oracles; production differentiability does not require
LLVM/Enzyme.

Hardware scheduling is separate from symbolic mathematics. Indexed
map/reduce/tensor operators may be tiled in KeOps/FlashAttention style so
logical pairwise or Jacobian matrices are never materialized.

## Certification

The fast solver is untrusted. The primary rigorous path is

\[
{\rm fast\ solve}
\rightarrow
{\rm Fortran\ interval/ball\ certifier}
\rightarrow
{\rm checked\ enclosure}.
\]

Lean is optional downstream assurance for generic theorems or especially
important claims; it is not required to obtain rigorous numerical intervals and
is not on the early critical path.

A result distinguishes numerical error, model-reduction error, physical-input
uncertainty and, for Monte Carlo, statistical uncertainty. Reconstruction and
diagnostic uncertainty are owned by TIAGO and remain a separate probabilistic
layer conditional on the selected KIN6D model. Deterministic certificates and
probabilistic confidence statements are never conflated. Adjoints are used for
sensitivities and for goal-oriented numerical/model-error estimation.

Every reported state, field, observable and sensitivity carries an applicable
bound and its claim status, or explicitly states that the bound is unavailable.
Plots show the full applicable intervals. Pointwise bands require pointwise or
whole-cell bounds; an L2 radius cannot be drawn as a pointwise band. Particle
Dirac measures require weak-observable or specified reconstruction bounds.

Certification and differentiation accompany operators from their first use.
The workhorse certificate may use a validated coarse rounding budget when it
meets the requested accuracy. Deterministic containment remains mandatory for
a deterministic certified result; a tighter or higher-precision audit is
optional. General nonlinear certificates follow the linear/stable cases.

## Performance

Performance is a permanent scientific workstream. The design principles are:

1. exploit mathematical separability before parallelism;
2. exploit regular tensor structure before generality;
3. store states and compact factors, not global operators;
4. avoid duplicate full-state layouts;
5. tile so reductions consume data while resident in cache/shared memory;
6. reuse factorizations and many-RHS structure;
7. make reduced models serve as preconditioners/control variates;
8. profile and optimize the actual owning dependency;
9. preserve a simple reference action for every optimized path;
10. measure time and memory to a requested certified observable accuracy.

The goal is the fastest measured total time at matched model, observable,
accuracy and claim strength. Include setup, candidate solve, certification and
requested derivatives, and report memory and measurement variability. This is
an optimization objective, not an established speed ranking. Generic benchmark
support and reproducible baselines belong upstream before new hot paths.

The controlled Fortran stack is FortSym, FortFEM, FortNum, FortSparse, FortAD,
FortIO, FortPlot, FortArray and Fo, with PETSc and standard accelerator
libraries used behind replaceable interfaces where useful.

## Documentation map

Read [DESIGN.md](DESIGN.md) for the canonical architecture and
[ROADMAP.md](ROADMAP.md) for staged delivery. Important companion documents are:

- [Equilibria, imported candidates and TIAGO inference](docs/equilibria-inference.md)
- [Configuration, UX, I/O and reproducibility](docs/interface-reproducibility.md)
- [Complete plasma-physics capability catalogue](docs/capability-catalogue.md)
- [Microscopic relativistic many-body reference](docs/microscopic-many-body.md)
- [Quantum plasma fundamentals and dense-matter bridge](docs/quantum-plasma.md)
- [Adaptive multiscale, multifidelity and uncertainty architecture](docs/adaptive-multiscale-uq.md)
- [Numerical-method catalogue](docs/numerical-methods.md)
- [Reacting and partially ionized plasmas](docs/reacting-partially-ionized.md)
- [Radiation hierarchy](docs/radiation.md)
- [Representations](docs/representations.md)
- [Stability, response and bifurcation](docs/stability.md)
- [Transport and multiscale coupling](docs/transport-multiscale.md)
- [Student laboratory](docs/student-laboratory.md)
- [Model family](docs/model-family.md)
- [Species models](docs/species-models.md)
- [Moment models](docs/moment-models.md)
- [Collisions](docs/collisions.md)
- [Numerics](docs/numerics.md)
- [Differentiability](docs/differentiability.md)
- [Adjoints and structure](docs/adjoints-structure.md)
- [Certification](docs/certification.md)
- [General transport demonstration](docs/general-transport.md)
- [Transport certificate argument](docs/transport-certificate.md)
- [Performance](docs/performance.md)
- [Geometry](docs/geometry.md)
- [Coupling](docs/coupling.md)
- [Edge physics](docs/edge-physics.md)
- [Space plasmas](docs/space-plasmas.md)
- [Extensions](docs/extensions.md)
- [Analytical benchmark atlas](docs/benchmarks.md)
- [Tooling](docs/tooling.md)

Scientific/application code and the certifier are Fortran. No Python runtime is
required to configure, run, certify or inspect KIN6D. Optional Python and web
layers may provide notebook, cluster and browser convenience on top of the same
Fortran-owned configuration/schema and result contract. Generated graphics and
large benchmark products are not source artifacts.

Chris&AI
