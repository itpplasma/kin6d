# kin6d roadmap

Status: proposed progression, 2026-10-05. No milestone is implemented yet.
The first product is a certified steady solver; the eventual product also
supports certified evolution through the same physics operators.

## Destination

Deliver stationary, deterministic, full-distribution, relativistic 3X3P
kinetics for electrons and several ion species, self-consistent Maxwell
fields, and nonlinear Beliaev–Budker collisions in general physical 3D
geometry. Cartesian vectors and arbitrary maps are fundamental; the initial
fusion chart is cylindrical/Fourier with a compatible momentum frame. Every accepted solution has a continuous-problem enclosure and certified requested
observables. Lean verification is developed from the first reference case.

Deliver three connected scientific tiers: certified full-kinetic reference,
certified reduced/hybrid production, and certified comparisons between them.
Full kinetics serves local patches, surfaces, global axisymmetric steady
states, and affordable 3D response. Reduced models target larger domains and
longer evolution, including whole-device gyrokinetic turbulence. Per-species
choices extend through FK/GK/DK, parallel kinetic, orbit averaged, fluid,
adiabatic, and fixed; only mathematically supported combinations are enabled.

The physical mission spans TCV/AUG/LHD/W7-X, then SPARC/ITER/ARC/DEMO-like
designs and Helios/Stellaris. [Device envelope](docs/devices.md) specifies
published design points, physical-parameter gates, and generic geometry.
For each frozen equilibrium progress through geometry smoke tests, small
kinetic cases, and production-resolution certificates. This ladder is part
of the project, not an optional final validation campaign.

Later deliver transient trajectory enclosures and controlled comparisons with
more detailed collision/correlation models. The mathematical target is a
named continuum model, with additional model-reduction guarantees only when
their theorems are established. Full global BBGKY or real-particle simulation
is not a production milestone.

## Decision rule

Every new abstraction must remove existing duplication, enable the immediately
next capability, or address a measured performance or accuracy problem.
Decide at the latest point that still permits the next complete deliverable.
Keep the old problem executable when a new capability arrives. A configuration
change should activate additional toroidal modes in the cylindrical chart
rather than a new solver. General charts do not require a toroidal angle.

The first executable commit is runnable end to end and produces meaningful
plots through Fortran postprocessing. Every later increment preserves that
working demo, shows candidates with complete certified interval bands, and
overlays its general model's limiting configuration with the retained simpler
case. Track plot sources and regenerate outputs on demand; never commit plot
images, PDFs, SVGs, or animation files.

Reserve only the durable mathematical boundaries now: shared evolution and
constraint actions, species and basis definitions, reconstruction, moments,
and a separate certificate path. Concrete arrays, field unknowns, radial
bases, preconditioners, and proof representations may evolve behind them.
Avoid speculative backend registries and a universal phase-space framework.

Keep external-field/domain evaluation independent of machine names and
species evolution. Do not implement every species reduction now: introduce
its smallest state/reconstruction/field/collision seam when the next case
requires it. A shared operator interface does not make FK/GK/DK equations
identical. Explicitly derived model actions use the same steady and transient
consumers; fixed species have no evolution block.

Fortran, mandatory fortsym numerical-kernel emission, and the Lean proof
boundary are settled policy. Native CUDA remains an optional performance
backend; CMake/Ninja and fo own the build workflow. Flexibility applies to
algorithms and the smallest interfaces required by each milestone, while
preserving these choices. Generic dependency gaps are repaired at their
owning repositories and promoted to main promptly after focused verification;
kin6d never adopts a permanent private numerical workaround.

Each row below introduces one principal difficulty. Small substeps within a
row must still preserve an executable candidate and an independent oracle.
The next deliverable is specific; later rows state outcomes and dependencies,
rather than pretending that unresolved mathematics has an implementation date.

## Analytical Benchmark Atlas

The [atlas specification](docs/benchmarks.md) and
[machine-readable inventory](benchmarks/atlas.json) are maintained before
implementation. Activate the smallest references with each milestone: R0
bounded FP; R1 exact characteristics; R2/P1 analytic fields and manufactured
3D; B/C collision identities and rigorous closed-form/potential references;
G0 circular tokamak and straight/toroidal helical constructions; later
large-aspect-ratio, classical/neoclassical, GK and MHD limits.

Each capability reports its applicable E0/E1/A0/M/I/X reference matrix and
actual proof level. Unimplemented references stay planned; leading-order
asymptotes without remainders remain diagnostics. Generate candidate/reference
interval overlays on demand, preserve earlier cases, and repair required
generic special-function/enclosure support upstream. Build no universal test
framework or entire special-function library before the next case needs it.

## Next deliverable

Build the bounded homogeneous 1P Fokker–Planck reference case specified in
[certification](docs/certification.md#first-continuous-certificate). It has
positive known coefficients, prescribed Dirichlet data, an analytic
manufactured solution, a coercivity-based continuous error bound, and one
certified moment. Include the one-dimensional pointwise-band theorem and a
Fortran postprocessor producing candidate/exact/enclosure plots, two-resolution
overlays, and the moment interval. The default CPU demo needs no Python or GPU.

Its gate requires a numerical candidate, a separately evaluated rigorous
residual, a proved stability constant, a Lean-checked instantiation of the
continuous enclosure theorem, and rejection of an invalid certificate.
At two resolutions, both enclosures must contain the independent exact
solution. Report actual widths and error contributions. This first bounded
problem does not need Landau, magnetic geometry, parallelism, or an unbounded
momentum-tail theorem.

Use Fortran, fortsym-generated numerical leaves, fortnum verified arithmetic,
CMake/CTest through fo, and Lean for this vertical slice. Pin the dependency
revisions and record runtime, peak memory, certificate time, and theorem
assumptions. Select the smallest required numerical interfaces and artifact
format. Add unsupported generic enclosure/codegen operations upstream when
R0 requires them; do not write a handwritten substitute in kin6d.

## Steady development sequence

| Milestone | New difficulty | Retained or independent oracle | Continuous acceptance gate |
|---|---|---|---|
| R0: first reference | Rigorous residual-to-error bridge in bounded 0X1P | Analytic driven FP solution | Coercivity theorem, checked arithmetic, one moment interval, Lean application |
| R1: 1X1P reference | Characteristic transport and inflow/periodic boundaries | Analytic transport/source cases | Transport and boundary defects with a verified stability bound |
| R2: coupled reference | Electrostatic field feedback and direct nonlinear solve | Manufactured Poisson-coupled case; VlaPy where equations and BCs overlap | Continuous nonlinear field/kinetic enclosure on a constrained space |
| P0: production skeleton | Real toroidal/harmonic representation and cylindrical momentum frame | Cartesian free-flight and rotationally symmetric analytic cases | Reconstruction, coordinate terms, and projected/complement estimates for simple physics |
| P1: early 3D smoke test | `n=0,+/-1`, mode derivatives, and products | Analytic three-dimensional manufactured solution | Fourier normalization, padded product, and continuous tail checks |
| B0: relaxation baseline | None/Krook selection and fixed reservoir | Analytic relaxation/source problem | Stated reservoir and continuous enclosure; small model metadata |
| B1: angular model | Lorentz pitch-angle scattering | Exact spherical-harmonic decay spectrum | Angular/frequency bounds and appropriate reservoir balances |
| B2: conservative relaxation | Moment-dependent BGK target | Exact matched-target equilibrium | Target moment/derivative bounds and declared conservation/entropy theorem |
| B3: simple FP model | LBD drag and isotropic diffusion | Manufactured FP and equilibrium cases | Coefficient/moment closure and continuous FP bound |
| B4: discrete correction | Reusable constrained conservation enforcement | Exact invariant functionals and pair balances | Solvable correction, preservation claims, and consistency/derivative defect bounds |
| C0: Coulomb test action | Differential/test-particle Landau in fixed background | Analytic coefficient/harmonic limits | Continuous background-potential and reservoir-model bounds |
| C1: Coulomb linear action | Integral/field response composed with C0 | Bilinear identity, AD derivative, and agreed linearized external case | Reference assumptions, full linear stability, and potential complement |
| C2: nonlinear collision | Single-species Landau quadratic correction | Maxwellian equilibrium and pinned PETSc comparison | Nonlinear derivatives, collision tails, and nullspace constraints |
| C3: interspecies collision | Pair exchange on unequal momentum grids | Two-species equilibrium and independent PETSc cases | Particle invariants, paired exchange, and cross-grid error bounds |
| C4: relativistic test action | BB pair potentials and reservoir coefficients | Maxwell–Jüttner background and nonrelativistic limit | Relativistic radial potential and test-action bounds |
| C5: relativistic linear action | BB field-particle response | Shared pair/AD identities | Reference-dependent relativistic linear stability |
| C6: relativistic nonlinear action | BB quadratic correction and coupled invariants | Maxwell–Jüttner equilibrium and consistent limiting cases | Nonlinear pair/derivative/complement bounds and relativistic balances |
| K0: kinetic assembly | Marry orbit transport and nonlinear collisions in a simple prescribed-field case | Reference transport and collision cases separately | Coupled continuous residual and local inverse/stability estimate |
| E0: electron experiment | Electron kinetics in 3D fields with fixed ions | Preserved reservoir/collision/orbit cases and manufactured electron problem | Frozen ion distribution, charge/current, and certified electron observables |
| K1: several species | Electrons, D/T, ash, and impurity distributions | Preserved one/two-species cases | All pair bounds and global balances with species-specific scales |
| F0: electrostatic production | Self-consistent charge response in production representation | R2 and manufactured harmonic Poisson problems | Gauge, field reconstruction, input and boundary bounds |
| E1: electron feedback | Self-consistent fields with ions still fixed | E0 prescribed-field limit and F0 field oracle | Coupled electron/field enclosure and explicit reservoir exchanges |
| F1: electromagnetic production | Self-consistent magnetic response and compatible Maxwell geometry | Analytic field/charge-continuity cases | Curl/divergence, exterior treatment, constraints, and coupled certificate |
| G0: idealized torus | Full-orbit axisymmetric geometry and sources/walls | Symmetry slices and exact prescribed-orbit cases | Geometry, boundary traces, and orbit-stability obligations |
| G1: real tokamak | Frozen AUG/TCV external field and discharge inputs | Independent equilibrium evaluations; named kinetic codes in agreed limits | Physical-parameter device scope, field/input enclosures, and certified observables |
| G2: nonlinear 3D response | Several toroidal modes with realistic mode coupling | `N=0` reproduces G1; KAMEL in a specified linear-response limit | Nonlinear toroidal complement and observable bounds |
| G3: stellarator geometry | W7-X/LHD field periods and a genuinely 3D boundary | Independent field/period maps and retained G0/G2 cases | Geometry transforms, symmetry scope, and physical-parameter 3D observables |
| G4: reactor envelope | High-field/large-size design inputs | Frozen SPARC/ITER/ARC/DEMO and Helios/Stellaris revisions | Same numerical claims at physical masses, fields, profiles, and requested tolerance |

P0 can start at `n=0`, low `L`, one species, and a tiny meridional mesh.
Low angular order is a manufactured/projection test, not a claim that an
arbitrary full distribution is accurately resolved. P1 happens before
complex collisions or geometry can embed accidental axisymmetric assumptions.
Afterwards use `n=0` for most physics development and retain P1 as a regression.

Relativistic momentum kinematics are present in P0. Nonrelativistic collision
benchmarks explicitly select their consistent limiting kinematics. C4–C6
follow Landau directly because the BB potential formulation fits the basis;
their completion still depends on analytical certification and species coupling.
Test, field, linear, and nonlinear actions are compositions of each family's
shared generated pair primitive. [Collision specification](docs/collisions.md)
defines model metadata, conservation modes, and certified model differences.
Field-only actions are retained as diagnostics. Reduced linearized models
are optional comparison branches, not prerequisites for production BB.

The homogeneous collision gate initially certifies equilibrium or a clearly
specified steady forced problem. Non-equilibrium relaxation histories need
the transient certificate below, even if time evolution is used as a
numerical comparison. Invariant diagnostics and entropy checks are independent
behavioral evidence, not a substitute for a continuous solution enclosure.

E0 can branch from P1 and B0/B1 with the smallest physically stated electron
model; it need not wait for C6 or evolving ions. Promote its collision model
along the retained collision ladder as those certificates become available.
E1 depends on the selected self-consistent field block. Geometry smoke tests
for all available frozen device inputs can begin before the kinetic device
gates. Implement G3 as soon as the first 3D boundary case requires it rather
than imposing a tokamak-only development order.

## Reduced and hybrid development sequence

These branches begin when a frozen scientific case cannot meet its cost or
scope target with the reference model. They do not delay R0 or require a
universal species hierarchy at P0. Each row needs its own continuous-model
certificate, runnable plots, and the nearest affordable FK comparison.

| Branch | First useful configuration | New obligation |
|---|---|---|
| H0: prescribed species | Kinetic electrons with fixed ions | Reservoir balances and independent field selection, retained as E0/E1 |
| H1: simplest closure | Adiabatic or fluid electrons with kinetic ions | Declared closure, consistent field response, and observable comparison |
| H2: DK/FK hybrid | DK electrons with FK ions | Guiding-center action, Jacobians, polarization/magnetization, transformed pair exchanges |
| H3: GK production | DK/GK electrons and GK thermal ions | Ordered gyroaveraging/FLR, compatible fields and collisions, full-surface case |
| H4: energetic hybrid | DK electrons, GK D/T, FK alphas | Supported cross-representation pairs and physical energy/charge balances |
| H5: cheaper kinetic closures | Parallel kinetic, bounce/orbit averaged, or moments when demanded | Perpendicular/orbit/moment closure, boundary and exchange assumptions |
| H6: reduction bridge | Matching full and reduced observables | Independent interval differences first; rigorous reduction-defect/stability theorem later |
| H7: adaptive fidelity | Cheapest model meeting the requested total bound | Proved numerical plus model bounds and certified state transfer when switching |

GK, DK, and cheaper closures discard different physics; promotion order is
case-dependent, not a universal accuracy ladder. The model catalogue is an
intended feature set, not an assertion that every species/pair/field product
is supported. [Species models](docs/species-models.md) defines compatibility.
An interval difference between two certified solutions is useful before a
general asymptotic remainder theorem exists, but does not replace that theorem.

## Bulk MHD and response branches

[Model-family specification](docs/model-family.md) defines bulk population
ownership, equilibrium imports, MHD closures, and species/spatial hybrids.
These branches are scientific options, not prerequisites for R0 or FK.
Names below describe mathematical capabilities, not separate executables.

| Branch | First useful case | Added obligation |
|---|---|---|
| EQ0: imported equilibrium | Next frozen device candidate | Format/sign/geometry oracle; input uncertainty distinct from equilibrium certificate |
| EQ1: axisymmetric equilibrium | Fixed-boundary Grad–Shafranov | Exact Solov'ev/shaped reference, continuous field/domain enclosure |
| EQ2: free boundary | Coils and external structures | Exterior, conductor and moving-boundary bounds; X-points as a later increment |
| EQ3: external 3D candidate | VMEC/DESC/GVEC or separately SPEC | Correct topology/model, constraints, complement and stability theorem |
| MHD0: ideal static | Simple circular ideal equilibrium | Force balance and bulk closure, atlas and retained field geometry |
| MHD1: linear forced | Ideal response on enclosed equilibrium | Equilibrium uncertainty, constraints, forcing and local resolvent bounds; ARES/GPEC overlap only |
| MHD2: dynamic linear | Analytic Alfvén/sound or selected eigenmode | Spectral isolation/normalization or certified forced trajectory |
| MHD3: resistive linear | Simple layer/induction response | Ohm/thermal law and wall/exterior conditions before tearing/RWM claims |
| MHD4: kinetic feedback | MHD thermal bulk plus DK/FK population | No double counting, derived exchange/field action, frozen versus self-consistent distinction |
| MHD5: extended/multi-fluid | First Hall/two-fluid case requested | Complete derived Ohm/energy closure and reduced-limit certificate |
| MHD6: spatial hybrid | Bulk MHD with a kinetic resonant layer | Interface/matching, location uncertainty, flux and model-defect bounds |
| MHD7: nonlinear XMHD | A concrete nonlinear device question | Nonlinear trajectory certificate, existing analytic gates and scoped device promotion |

Imported candidates remain available while internal equilibrium solvers
mature. A 3D residual norm alone is not an equilibrium certificate. Introduce
new state/closure/linearization seams as the first branch needs them. Internal
3D iteration, inverse equilibrium, walls and nonlinear disruption machinery
remain deferred until required; do not reproduce all mature codes in advance.

## Time-dependent development sequence

The operator seam is present at P0. Implement a transient consumer only when
its smallest useful problem is requested; it can branch from B0/B1/C2 or K0
without waiting for realistic tokamak geometry.

| Milestone | New difficulty | Oracle and acceptance gate |
|---|---|---|
| T0: short trajectory | Temporal reconstruction and propagated uncertainty | Analytic orbit or linear relaxation; continuous tube and one observable |
| T1: nonlinear relaxation | Stiff collision evolution | Homogeneous collision invariants, independent relaxation benchmark, tube with momentum complements |
| T2: kinetic waves | Dynamic Maxwell with consistent initial constraints | Analytic linear-wave limit and charge continuity; field/kinetic/temporal certificate |
| T3: multiscale evolution | Large scale separation | Resolved short-interval oracle; case-specific growth estimate and rigorous fast-component error |
| T4: slow evolution | Useful long horizons and bounded output cost | Previous slab certificates propagated without reset; stated full-state or observable-only guarantee |

[Time-dependent extension](docs/transients.md) defines the timescales and
trajectory theorem obligations. Small physical periods range from roughly
picoseconds to tens of nanoseconds in the example parameters. Transport can
evolve over milliseconds to seconds. Coverage of this range is an eventual
scope, not a claim that one coarse timestep or one stable integrator resolves
all of it. Tight long-time bounds remain conditional on verified stability.

Select implicit, exponential, multirate, or space–time methods only when the
preceding transient exposes the need. Preserve the steady consumer. Steady
existence and transient stability remain separate certificates.

## Higher-fidelity research sequence

After the certified Landau/Beliaev–Budker spine works, introduce one reference
model at a time in homogeneous or small geometries:

1. Dynamically screened Balescu–Lenard, with specified dielectric stability
   and short-distance treatment.
2. A matched soft/strong-collision construction, with derived overlap
   subtraction, pair invariants, and demonstrated domain of validity.
3. Mean-force comparisons when equilibrium pair data and their applicability
   are available; no automatic ordering above every nonequilibrium closure.
4. Small explicit pair-correlation experiments with a specified higher-order
   closure. Treat truncated BBGKY as another model until its remainder is bound.
5. Certified operator-difference intervals or validated parameter continuation.
   Rich-model defect certification may avoid a full rich-model solve when
   all its derivative and complement bounds can be established.
6. Closure remainder theorems connecting kinetic predictions to many-particle
   marginals in verified regimes. This is a separate research program with
no assumed completion date.

Hard collisions, radiation, ionization, recombination, and reaction physics
enter when a target regime requires them, with explicit balances and matching.
No extension inherits a certificate merely by sharing a software interface.

## Spatial numerical foundation and performance gates

[FortFEM and production numerics](docs/numerics.md) defines ownership and the
checked upstream snapshot. Reuse its mathematical spatial layer rather than
copying generic kernels into kin6d. FEEC is the compatible structure for
Maxwell/magnetic blocks; fields and distributions can use distinct spaces.
Mapped tensor patches are preferred where suitable, with mortar/cut/unstructured
escape paths when needed. Keep R0's reference small and executable.

Before the first spatial client, select a pinned supported FortFEM interface
and test a minimal independent map/space/operator case. Before large kinetic
spatial work, measure span-local basis evaluation, sum-factorized partial
assembly and many-RHS batching. Implement missing generic slices upstream with
focused independent oracles and prompt verified-main promotion. Retain assembled
reference actions; add matrix-free CPU first and CUDA only after profiling.

Production stores compact factors and a canonical kinetic layout, with tile
scratch and explicit work-vector budgets. Matrix-free generated JVPs, reused
radial factors and low-memory physics preconditioning precede GPU/MPI complexity.
AoSoA tiling, resolution classes, convolution/FFT and Gaunt/angular-transform
crossovers are measured choices. No global kinetic Jacobian, permanent duplicate
full-state layouts or hidden large-GMRES workspace. Exact mathematical
transformations still need rounding/conditioning bounds in the certifier.

Measure time and peak memory to certified observable accuracy. Include solver,
certifier, padding, collision potentials, communication and I/O. Structure
identities and equivalent optimized actions do not alone prove continuous
stability or tails. Missing certificate or throughput evidence stays explicit.

## Long-range application branches

These branches are demand-driven and do not change the specific next R0
deliverable. Each row is split into runnable increments with retained limit
plots and an applicable continuous certificate before production promotion.
Shared libraries may gain thin application drivers as consumers arise.

| Branch | Progression and admission evidence |
|---|---|
| Geometry | Cartesian identities → mapped elements → cylindrical fusion chart → Figure-8 G-frame geometry gate; independent map/field checks and rigorous Jacobian/inverse bounds |
| Coupling | Prescribed test field → fixed bulk with induced response → coupled bulk evolution; separate population roles, complete source and reservoir ownership |
| Edge/neutrals | Exact neutral streaming/reactions → slab recycling → charge exchange/sheath → optional EIRENE coupling → deterministic neutral certification → physical exhaust benchmark |
| Atomic/molecular | Effective rates → CR states → non-Maxwellian cross-section integrals → differential products; data validity, provenance and exchange bounds at each increment |
| Solar System | Parker/Harris/linear waves → open-boundary reduced case → selected hybrid/kinetic application; gravity/rotation/wave operators only when the case requires them |
| Fields | First Poisson/Maxwell subset → a specifically derived quasineutral or Darwin model when cost requires it; constraint/energy identities and full-model comparison status |
| Heating/sources | Prescribed source/QL benchmark → NBI birth/attenuation → ray/QL or full-wave coupling; power/momentum ledger and resonance/spectrum bounds |
| Fusion | Validated rate quadrature → fuel depletion and product sources → alpha/neutron diagnostics; identical-pair factors, mass-defect and product kinematics |
| Conductors/circuits | Vacuum/coil inputs → conducting-wall/circuit response → free-boundary evolution/control candidate; gauge, power balance and compatible induced fields |
| Diagnostics | One requested bounded observable → instrument mapping → physical-device measurement comparison; observation noise separated from numerical intervals |
| Inference/optimization | Correct forward diagnostic → independently checked adjoint → constrained inference/design; global optimality/identifiability never inferred from convergence |
| Adaptation | Certified numerical error localization → one needed mesh/mode adaptation → validated species/region model transfer; reduction and interface bounds required |

See [geometry](docs/geometry.md), [coupling](docs/coupling.md),
[edge physics](docs/edge-physics.md), [space scope](docs/space-plasmas.md) and
[extensions](docs/extensions.md) for mathematical limits. A Monte Carlo
neutral/noise estimate has statistical status; it cannot discharge a
deterministic continuous certificate without additional deterministic bounds.
A norm enclosure does not automatically bound loss/hitting observables.

Units/normalization and full result provenance are obligations from the first
case. Report total time and peak memory to reach a certified observable target,
including solve, certifier and postprocessing costs. General relativity, QED,
pair creation, full neutron transport and comprehensive material/radiation
hydrodynamics remain outside initial scope. Planned features are not
implemented switches or committed performance promises.

## Deferred decisions and triggers

| Decision | Decide when | Required evidence |
|---|---|---|
| Dependency interfaces and build profiles | First vertical slice | Working Fortran/CMake/fo/Lean integration without Python |
| Missing generic dependency capabilities | First concrete blocker | Independent upstream reproducer, verified main fix, consumer recheck |
| Final continuous norm and inverse construction | Before each new certification gate | Applicable theorem and computable constants |
| Infinite-momentum basis or tail treatment | First whole-space claim | Verified tail, origin, and potential-boundary estimates |
| Compatible Maxwell discretization | F1 preparation | FortFEM FEEC complex, independent curl/divergence identities, boundary/gauge and source transfer |
| Spatial batch/tensor optimization | First measured many-channel bottleneck | Upstream span/sum-factorization/batch action with assembled and analytical oracles |
| Preconditioner and matrix-free actions | Krylov/assembly cost becomes material | Profile and retained behavioral oracle |
| GPU collision batching | Collision work dominates | End-to-end profile, including certification |
| MPI decomposition | Workstation capacity or time is insufficient | Memory/runtime model and measured communication need |
| Local magnetic frame or angular adaptivity | Full harmonic cost dominates | Frame derivation and certified complement control |
| AMR, compression, mixed precision | Measured state or arithmetic bottleneck | Error contribution included in the acceptance theorem |
| Temporal method beyond first driver | A transient fails accuracy/cost target | Resolved oracle and certified temporal defect |
| First reduced species seam | A frozen case requires a cheaper model | Derived field/pair compatibility and retained reference comparison |
| Device input adapter or mapped boundary | Next device geometry gate | Pinned input, independent geometry oracle, and rigorous transformation bounds |
| AP formulation/preconditioner | Physical scale separation prevents a useful full-model solve | Original PDE certificate, limiting oracle, and measured accuracy/cost |

## Resource planning

For four species with `128` radial coefficients and full `L=8`, there are
`81` angular coefficients per radial node. Real-equivalent FP64 kinetic
state storage is `8 N_x N_mode sum_s(N_p,s (L_s+1)^2)` bytes.

| Configuration | One kinetic state |
|---|---:|
| Homogeneous, `n=0` | `0.316 MiB` |
| 256 configuration degrees, `n=0` | `81 MiB` |
| `64 x 128` configuration degrees, `n=0` | `2.53 GiB` |
| Same, `n=-4,...,4` | `22.78 GiB` |

Here configuration counts mean actual scalar degrees of freedom, including
element polynomial order. They are not counts of high-order elements.
Toroidal storage assumes exploiting reality; complex storage conventions can
increase the total. Fields, radial potentials, dealiased angular work arrays,
Newton/Krylov vectors, preconditioners, and interval data are additional.

A 64 GiB workstation is a plausible early platform, not a guaranteed fit for
the last two examples. Twenty Krylov vectors for a 2.53 GiB state already cost
about 51 GiB. Measure memory before adding hardware or abstractions; keep
Krylov storage bounded and stream certificate work where possible.

[Device memory envelope](docs/devices.md#memory-envelope) adds 48-radial-node
machine-scale examples and contrasts the initial 81-coefficient basis with a
later 39-coefficient optimization. It also states the eight-device planning
scenario without treating state-vector storage as a full capacity estimate.

## Milestone completion

Freeze the model, source revision, and candidate bytes. Run the independent
behavioral oracle and the required continuous checker. Record norm, assumptions,
radius, observable intervals, residual/stability budget, and resource use.
Run the end-to-end demo and postprocessor, inspect real interval bands and
earlier limiting-case overlays, and verify that no generated plot artifacts
are tracked. [Demonstration policy](docs/demonstrations.md) defines the figure
contract. Runnable behavior and meaningful certified plots are part of every
executable milestone, including R0.
Preserve earlier cases, including the early 3D smoke test. A failed continuous
gate leaves the capability experimental and records the missing theorem or
bound; it is never replaced by apparent mesh convergence.

Production promotion additionally passes the scoped physical-parameter
AUG-like and W7-X-like gates. Expand coverage across the full device envelope
through geometry, small kinetic, and production levels. Name the model,
equilibrium revision, observable, physical parameters, and certified width in
each device report. Component reference completion remains possible before
these machine gates; it must not be advertised as production readiness.

Chris&AI
