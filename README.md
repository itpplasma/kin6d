# kin6d

Design for a deterministic, stationary, full-distribution plasma kinetic
solver with rigorous solution enclosures.

The full-kinetic representation includes three spatial and three momentum
dimensions from its first version. Physical space and vectors are Cartesian;
computational charts are selectable. The initial fusion chart uses cylindrical
elements and toroidal Fourier modes, with the compatible momentum frame so
axisymmetric runs activate only the zero mode. Relativistic kinematics, species-specific momentum
grids, spherical harmonics, and a direct steady solve form the numerical
spine. A separate certifier must establish bounds for the stated continuous
problem before a result can be accepted as certified.

All-species full kinetics is the reference. The planned family also supports
per-species reduced and hybrid models, including kinetic electrons with fixed
ions. Its physical mission spans TCV, AUG, LHD, W7-X, SPARC, ITER, ARC,
DEMO-like designs, Helios, and Stellaris. Generic 3D fields and domain geometry
are the physics boundary; symmetry is a case choice. Production promotion
requires physical-parameter device benchmarks and certified observables.

The project is one mathematical engine with selectable species, bulk-plasma,
collision, field, geometry and time configurations. MHD and kinetic–MHD are
explicit long-range branches. Established-code names describe comparison
cases, never separate solver architectures.

**Repository status:** design only. There is no solver, certifier, or Lean
proof implementation yet. Full nonlinear kinetic certification is a research
objective; this repository does not claim that such bounds already exist.

Read the [design document](DESIGN.md) first. Its companion specifications are:

- [Certification](docs/certification.md): theorem obligations, result states,
  error budgets, observable intervals, and the first continuous certificate.
- [Roadmap](ROADMAP.md): the next deliverable, milestone gates,
  independent validators, and resource estimates.
- [Time-dependent extension](docs/transients.md): shared operators, physical
  timescales, temporal methods, and trajectory certificates.
- [References](docs/references.md): checked sources and their actual scope.
- [Performance engineering](docs/performance.md): calibrated benchmarks,
  best-known time–memory–accuracy frontiers, profiling and upstream ownership.
- [Numerical foundation](docs/numerics.md): FortFEM ownership, FEEC/IGA/DG,
  tensor actions, batching, layouts, derivatives and solver memory.
- [Tooling and kernel policy](docs/tooling.md): the Fortran-native stack,
  mandatory fortsym generation, CUDA boundary, and upstream repair rules.
- [Collision hierarchy](docs/collisions.md): shared pair primitives,
  differential/integral choices, conservation treatment, and model intervals.
- [Runnable demonstrations](docs/demonstrations.md): plots from the first
  executable increment, certified bands, and retained limiting-case overlays.
- [Species models](docs/species-models.md): FK/GK/DK and cheaper closures,
  consistent mixed-species coupling, fixed reservoirs, and reduction bounds.
- [Device envelope](docs/devices.md): machine/design missions, generic 3D
  geometry, layered CI gates, physical parameters, and memory planning.
- [Model family](docs/model-family.md): the chameleon philosophy, bulk MHD,
  equilibrium candidates, species/spatial hybrids, and response certification.
- [Analytical Benchmark Atlas](docs/benchmarks.md): exact/asymptotic references,
  simple circular/helical geometry, and the [planned catalogue](benchmarks/atlas.json).
- [Geometry](docs/geometry.md): Cartesian physical vectors, mapped charts,
  optional FCI, and the Figure-8 geometry gate.
- [Coupling](docs/coupling.md): background, external drive and induced response;
  active, test and reservoir populations.
- [Edge/reaction/surface physics](docs/edge-physics.md): transport, neutral,
  atomic, radiation, sheath and wall axes, with uncertainty accounting.
- [Solar System plasmas](docs/space-plasmas.md): open domains, forces,
  anisotropic closures, wave scattering and space benchmarks.
- [Further extensions](docs/extensions.md): field reductions, heating, fusion
  reactions, circuits, diagnostics, adaptation, units and provenance.
- [Literature policy](docs/literature-policy.md): external corpus, lawful OA
  resolution, Zotero import, provenance, and reviewed citation expansion.

Scientific and application code is Fortran. Every project-owned numerical
kernel is emitted by fortsym; native CUDA is an optional accelerator behind a
thin C ABI. The rigorous certifier is Fortran and the trusted proof layer is
Lean. Python is not a required build, runtime, plotting, or analysis layer.
Every executable milestone includes an end-to-end runnable demonstration and
Fortran postprocessing. Plot source is tracked; generated figures are ignored
and regenerated on demand.

The defining rule is **every accepted scientific solution carries a checked
enclosure**. Exploratory candidates remain available when certification fails,
with the missing bound recorded explicitly. Refinement is driven by proved
error contributions; mesh agreement is a debugging aid.

The bounds describe solutions of named mathematical models. Bounds for
operator differences and for reduction from many-particle physics are
separate claims, each requiring its own theorem.

Prepared on 2026-10-05 for the itpplasma group.

Chris&AI
