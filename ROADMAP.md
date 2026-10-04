# kin6d roadmap

Status: proposed progression, 2026-10-05. No milestone is implemented yet.
The first product is a certified steady solver; the eventual product also
supports certified evolution through the same physics operators.

## Destination

Deliver stationary, deterministic, full-distribution, relativistic 3X3P
kinetics for electrons and several ion species, self-consistent Maxwell
fields, and nonlinear Beliaev–Budker collisions in tokamak geometry. Every
accepted solution has a continuous-problem enclosure and certified requested
observables. Lean verification is developed from the first reference case.

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
change should activate additional toroidal modes rather than a new solver.

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
| K1: several species | Electrons, D/T, ash, and impurity distributions | Preserved one/two-species cases | All pair bounds and global balances with species-specific scales |
| F0: electrostatic production | Self-consistent charge response in production representation | R2 and manufactured harmonic Poisson problems | Gauge, field reconstruction, input and boundary bounds |
| F1: electromagnetic production | Self-consistent magnetic response and compatible Maxwell geometry | Analytic field/charge-continuity cases | Curl/divergence, exterior treatment, constraints, and coupled certificate |
| G0: idealized torus | Full-orbit axisymmetric geometry and sources/walls | Symmetry slices and exact prescribed-orbit cases | Geometry, boundary traces, and orbit-stability obligations |
| G1: realistic tokamak | AUG/ITER-like external field and discharge inputs | NEO-2/CQL3D and other named codes only in agreed overlapping limits | Field/input enclosures and certified observables for a stated case |
| G2: nonlinear 3D response | Several toroidal modes with realistic mode coupling | `N=0` reproduces G1; KAMEL in a specified linear-response limit | Nonlinear toroidal complement and observable bounds |

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

## Deferred decisions and triggers

| Decision | Decide when | Required evidence |
|---|---|---|
| Dependency interfaces and build profiles | First vertical slice | Working Fortran/CMake/fo/Lean integration without Python |
| Missing generic dependency capabilities | First concrete blocker | Independent upstream reproducer, verified main fix, consumer recheck |
| Final continuous norm and inverse construction | Before each new certification gate | Applicable theorem and computable constants |
| Infinite-momentum basis or tail treatment | First whole-space claim | Verified tail, origin, and potential-boundary estimates |
| Compatible Maxwell discretization | F1 preparation | Curl/divergence identities and boundary/gauge design |
| Preconditioner and matrix-free actions | Krylov/assembly cost becomes material | Profile and retained behavioral oracle |
| GPU collision batching | Collision work dominates | End-to-end profile, including certification |
| MPI decomposition | Workstation capacity or time is insufficient | Memory/runtime model and measured communication need |
| Local magnetic frame or angular adaptivity | Full harmonic cost dominates | Frame derivation and certified complement control |
| AMR, compression, mixed precision | Measured state or arithmetic bottleneck | Error contribution included in the acceptance theorem |
| Temporal method beyond first driver | A transient fails accuracy/cost target | Resolved oracle and certified temporal defect |

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

Chris&AI
