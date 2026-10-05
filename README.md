# kin6d

KIN6D is a plasma-physics engine whose highest-fidelity reference is full
relativistic 3X3P kinetic dynamics, but whose practical purpose is broader:
solve the cheapest mathematically controlled plasma model that answers the
requested question at the requested accuracy.

The project is built around one common physical and numerical architecture.
Equilibrium, MHD, gyrofluid, gyromoments, drift/gyro/full kinetics, continuum
and particle representations, collisions, turbulence, edge physics, neutrals,
walls, heating, space plasmas, diagnostics and optimization are configurations
or reductions of shared operators rather than separate named-code modes.

## Current status

The first executable slice is implemented: a CPU 1X1V prescribed-field
laboratory with a direct source-balanced steady solve, exact free-streaming
dynamics, FortSym-generated primal/JVP/VJP leaves, Fortran plotting and a
reviewed continuous snapshot bound evaluated with interval arithmetic.
See [docs/first-example.md](docs/first-example.md).

Self-consistent Poisson, PIC/PIF, the collision hierarchy, gyrofluid turbulence,
real-device DK and general nonlinear certification remain future milestones.
The repository must not describe planned capabilities as implemented.

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
physics implementation.

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
- atomic/molecular charge-state networks and radiation;
- RF/quasilinear/full-wave heating, NBI and fusion-product sources;
- coils, vacuum, conducting walls and external circuits;
- Solar-System plasma configurations;
- synthetic diagnostics, inverse problems, adjoints, optimization and control.

General relativity, strong-field QED, complete reactor neutronics, detailed
solid-state material evolution and general radiation hydrodynamics are outside
the initial scope.

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
uncertainty and, for Monte Carlo, statistical uncertainty. Deterministic
certificates and probabilistic confidence statements are never conflated.
Adjoints are used for sensitivities and for goal-oriented numerical/model-error
estimation.

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

The controlled Fortran stack is FortSym, FortFEM, FortNum, FortSparse, FortAD,
FortIO, FortPlot, FortArray and Fo, with PETSc and standard accelerator
libraries used behind replaceable interfaces where useful.

## Documentation map

Read [DESIGN.md](DESIGN.md) for the canonical architecture and
[ROADMAP.md](ROADMAP.md) for staged delivery. Important companion documents are:

- [Representations](docs/representations.md)
- [Student laboratory](docs/student-laboratory.md)
- [Model family](docs/model-family.md)
- [Species models](docs/species-models.md)
- [Moment models](docs/moment-models.md)
- [Collisions](docs/collisions.md)
- [Numerics](docs/numerics.md)
- [Differentiability](docs/differentiability.md)
- [Adjoints and structure](docs/adjoints-structure.md)
- [Certification](docs/certification.md)
- [Performance](docs/performance.md)
- [Geometry](docs/geometry.md)
- [Coupling](docs/coupling.md)
- [Edge physics](docs/edge-physics.md)
- [Space plasmas](docs/space-plasmas.md)
- [Extensions](docs/extensions.md)
- [Analytical benchmark atlas](docs/benchmarks.md)
- [Tooling](docs/tooling.md)

Scientific/application code and the certifier are Fortran. No Python runtime is
required. Generated graphics and large benchmark products are not source
artifacts.

Chris&AI
