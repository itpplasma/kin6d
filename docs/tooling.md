# Tooling and kernel policy

Status: required design, checked against upstream source on 2026-10-05.
The repository contains documentation only; the integration gates below
remain implementation work.

## Stack

| Responsibility | Tool | Boundary |
|---|---|---|
| Scientific/application code | Modern Fortran | Physics orchestration, concrete arrays, model configuration |
| Numerical leaves | fortsym | Mandatory symbolic definition and final kernel emission |
| Spatial mathematics | FortFEM | Generic mapped spaces, FEEC/IGA/DG, traces, interfaces, FCI and weak actions; plasma laws remain external |
| Sparse blocks | FortSparse/PETSc | Assembled references/coarse/field blocks, not a full kinetic Jacobian |
| Numerics | fortnum, BLAS/LAPACK, PETSc | Existing generic algorithms; PETSc Fortran interfaces for SNES/KSP/PC/TS |
| Automatic differentiation | fortad | Supported JVP/VJP/HVP transformations; final numerical leaves still emitted by fortsym |
| Rigorous evaluation | fortnum verified modules | Interval/ball arithmetic, rigorous FFT bounds, verified linear algebra |
| Data files | fortio | Its supported HDF5/NetCDF subset |
| Plots and diagnostics | fortplot | Native PNG/PDF/text output and supported plot/animation facilities |
| Labeled analysis arrays | fortarray | Profiles, field slices, and diagnostics; outside hot kinetic storage |
| Build and registered tests | CMake/Ninja, CTest | CPU initially; optional native CUDA language when needed |
| Development workflow | fo | Focused test/exec commands and resident Gremlin feedback |
| Accelerator kernels | Generated native CUDA | Thin C ABI and ISO_C_BINDING; no compiler-specific CUDA Fortran |
| Formal proof/checking | Lean | Certificate soundness theorem and checked numerical premises |

No required Python orchestration, generation, plotting, postprocessing, or
runtime layer belongs in kin6d. Use a Fortran generator driver against
fortsym's native or approved compiled CAS engines. Disable implicit fallback
to Python/SymPy adapters. A missing native operation is an upstream feature
request to implement, not a reason to introduce a hidden Python dependency.

Libraries may expose optional Python bindings or external reference tests;
those are not required kin6d interfaces. The first executable gate must
exercise the selected build, generation, tests, certifier, and plot path with
Python unavailable. Do not drop required scientific tests just to satisfy
that gate; repair a mandatory dependency instead.

## Kernel generation policy

**Every project-owned numerical kernel must be emitted by fortsym.** This
covers orbit/geometry coefficients, basis and moment factors, collision
leaves, source formulas, field coefficients, reconstruction leaves, rigorous
enclosures, and derivative leaves. Generated CUDA kernels follow the same
rule. Handwritten code may manage arrays, assemble or schedule calls, perform
I/O, and implement backend ownership or launch glue. It must not conceal a
new scalar numerical formula inside that glue.

Existing generic algorithms are supplied by upstream numerical libraries;
the rule does not require rewriting all of BLAS, PETSc, or fortnum. New
generic algorithms or enclosure primitives belong in their owning library.
Any application-specific numerical leaf they call is generated under this
policy. A toy or reference implementation is not exempt.

Retain a Fortran generator in `gen/` with the symbolic DAG, exact constants,
domain assumptions, requested targets, generator/dependency revisions, and
independent analytical oracle. Generated outputs carry provenance and are
never edited by hand. Pin the generator and regeneration settings so a
clean rebuild can reproduce the committed artifacts. Source equality is
a provenance check; behavioral evidence is a separate gate.

Use one mathematical definition to emit fast Fortran and rigorous leaves.
Where supported, derive the derivative expressions and CUDA leaves from the
same definition. Source AD through fortad can supply transformed expressions
or orchestration, but final numerical leaves must pass through fortsym
emission. If an AD-to-DAG bridge is missing, implement that generic bridge
upstream or use fortsym differentiation for the immediate kernel. Do not
waive the emission rule for fortad output.

Independence of an oracle comes from an exact analytical identity, checked
Lean result, or external algorithm/model with matching scope. A second
emission of the same formula is not an independent numerical oracle.
Application-owned oracle kernels also obey the emission rule; exact expected
constants and theorem statements are allowed as fixtures.

Keep plasma-specific symbolic derivations in kin6d's `gen/`. Generic CAS,
AD, quadrature, sparse/FEM, and verified-arithmetic capabilities stay in
their respective dependencies. Do not put a kin6d collision derivation in
`fortsym/app/` merely to make generation convenient.

## FortFEM ownership and performance

[Spatial policy](numerics.md) records FortFEM main
`2d7baa19fe1bef02ce279fa407830c8731bdc663` and source-inspection limits.
FortSparse main was observed at `7ad20738b3f888a3bbec37324eff3117c630221d`;
these observations are not a tested kin6d integration lock. Pin the selected
compatible dependency graph when a spatial client exists. High-level weak
forms compile to specialized generated leaves rather than hot-loop expression
interpretation. Existing generic numerical algorithms remain upstream; all
new application leaves and required generic numerical-leaf additions use
fortsym. Keep assembled reference and optimized CPU/CUDA paths with independent
behavioral oracles. No unmeasured whole-library performance rewrite is required.

## Checked capability boundaries

The upstream snapshot confirms fortnum interval/ball modules, radix-2
rigorous FFT/convolution bounds, and norm-wise verified linear algebra.
Its rigorous FFT uses binary128 interval twiddle construction; the first
integration must confirm that the selected compiler/profile supports the
required kinds and arithmetic. A cuFFT result is not automatically covered
by the CPU radix-2 algorithm's error theorem.

fortsym has a CUDA numerical emission target and paired Fortran float/
enclosure emission. The checked rigorous emitter supports algebraic
operations, integer powers, and square roots; its documented runtime
interface does not yet emit general `sin` or `exp` expressions. Its CUDA
target emits numerical leaves, not scheduling or device-memory management.
Do not advertise every expression/target combination as already supported.

R0's Gaussian and cosine require rigorous transcendental handling. fortnum
already exposes interval transcendental operations. When implementing R0,
add the missing generic fortsym enclosure interface/emission support with
a focused behavioral oracle, push that upstream fix, then recheck generation
in kin6d. Another exact manufactured solution may simplify the first case,
but cannot justify handwritten application kernels.

fortad supports a documented subset of Fortran and rejects unsupported
activity, aliasing, lifetime, and opaque-call cases. Design small supported
kernel interfaces; extend generic semantics upstream when a real case
requires them. Code generation and AD tests do not by themselves prove
the derivative or enclosure bridge in Lean.

fortio implements a subset of HDF5/NetCDF, rather than every API/layout.
Confirm the needed checkpoint and certificate layouts before fixing a file
contract. fortarray currently owns dense real64 arrays through rank four;
the kinetic state uses compact contiguous arrays with explicit mode/species
maps. fortarray remains useful for low-dimensional real-valued diagnostics.

## CUDA boundary

The CPU Fortran path remains an executable reference. Add CUDA only for a
profiled expensive kernel, beginning with a single batch operation. The
boundary is `Fortran -> ISO_C_BINDING -> extern C wrapper -> generated CUDA
leaf / CUDA libraries`. Plasma equations remain in their symbolic definitions.

Keep opaque contexts/device pointers and explicit dimensions, layout,
precision, species scales, and mode conventions in the interface. Specify
allocation ownership, host/device residency, stream ordering, synchronization,
and returned error status. A kernel launch wrapper must report launch and
asynchronous completion failures. CPU builds never configure a CUDA compiler
or require CUDA runtime libraries.

cuFFT, cuBLAS, cuSPARSE, CUDA graphs, and later CUDA-aware MPI are optional
implementation choices behind this boundary. Ordinary CUDA math is part of
the untrusted candidate path. The Fortran certifier independently establishes
valid bounds; extending verified arithmetic to CUDA needs its own proved
rounding/error contract and is not an initial requirement.

## Performance workflow

[Performance policy](performance.md) reuses native upstream benchmarks and fo
with distinct fast sentinels, dedicated audits and physical-science cases.
Current optional Python benchmark collectors/gates do not become mandatory
kin6d tooling; add a needed native generic gate upstream instead. Pin inputs,
compiler/hardware and exact executable, validate outputs independently, retain
raw repeated timing/RSS evidence and never silently rewrite a slower baseline.

## Build and feedback

CMake/Ninja owns configuration and mixed-language builds; CTest owns
registered tests and their fixtures, working directories, resources, and
scientific verdicts. Keep project presets/profiles authoritative. Enable
Fortran first and CUDA only in an explicit later accelerator profile.
Lean's Lake build is separate and its required certificate checks are part
of the milestone's declared verification contract.

fo delegates CMake projects to CMake/CTest in the checked version. Use its
public focused test/exec commands and one resident Gremlin development lane
when executable development starts. Public help currently exposes bare
`fo` as resident start/attach and `fo verify` as a broader gate; verify the
selected driver/profile before relying on either behavior. Standalone native
CMake compatibility remains a delivery target, not a completed universal
capability. Documentation-only work does not start an empty resident lane.

Pin the actual fo executable and SHA256 for an implementation session. Keep
source revision, compiler, flags, dependencies, and current-generation test
evidence with each promoted increment. A cache hit or an empty affected-test
selection is not an independent correctness result.

Verified arithmetic needs a strict compilation profile consistent with its
error analysis: no unsafe reassociation, unaccounted contraction, or flushing
subnormals. Check the library's exact requirements, including contraction
settings, for the selected compiler. Fast solver flags may differ because
its result remains an untrusted candidate. Plotting or exporting a bound
must preserve its outward enclosure, rather than rounding interval endpoints
inward for display.

## Upstream main policy

The controller has user authorization to fix defects and add missing generic
features in owned upstream dependencies and push locally verified changes
to their main branches promptly. This includes the owned lazy-fortran stack
and any other owned dependency needed by the next kin6d capability.

For each encountered defect or concrete missing capability:

1. Capture the consumer failure or required behavior with exact revisions,
   command, input, and an independent small behavioral oracle.
2. Locate the owning dependency and coordinate with its active controller.
   Preserve unrelated edits; use an isolated worktree when needed.
3. Make the smallest complete generic repair there. Retain the focused
   reproducer as an upstream regression test and run required local gates.
4. The controller integrates and pushes the verified repair to upstream main
   promptly. Do not defer it behind an unrelated redesign or wait for a
   broad CI drain when the applicable promotion rules allow local evidence.
5. Pin the resulting upstream commit in kin6d, regenerate affected artifacts,
   and recheck the original consumer behavior and certificate path.

Workers return evidence and never promote authoritative state. Stage explicit
paths and preserve other owners' work. For dependencies outside our control,
contribute through their actual upstream access/promotion mechanism; keep
the same reproducer and consumer recheck requirement.

Do not add permanent consumer workarounds, private numerical forks, or copies
of generic kernels. If an upstream repair is blocked, record the reproducer,
missing feature or defect, verification limit, and concrete next action.
Report the consumer capability as blocked or experimental until the needed
fix and proof obligations are met. Protected `faepop*`/`faepcr*` hosts remain
untouched without an explicit request naming the host.

Chris&AI
