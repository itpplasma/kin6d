# Tooling and kernel policy

Status: canonical implementation/tooling policy.

## Verified CMake profile repair

Fo `60c8a6b2e92396605f06bfe67217e1b7fb946b4d` repairs checked Debug flags
persisting in global `CMAKE_Fortran_FLAGS` after a public Fo call. Named profiles
select their configuration, explicit configurations retain priority, and flags
are scoped to that configuration. The independent public CLI fixture preserves
project macros, checks runtime bounds failures, and verifies a later native
Release preset in the same build directory. The previous driver fails; the
candidate passes that fixture and 73 backend assertions.

The verified task driver is `/var/tmp/fo-cmake-profile-driver-20261005/fo`, SHA256
`6bb8bb26cf7374d8206f506ac36e36009ae8b26e195357648771938c435c9d02`.
Native FPM was used for diagnostic bootstrap after an incomplete duplicate
cold Fo build. The original kin6d transport gates pass 2/2 and its cached
benchmark executes through this driver in Debug. Returning to native Release
passes 5/5, retains empty global flags and reproduces the exact final benchmark
executable digest. Final performance measurements use that clean Release.
This historical evidence concerns delegated profile behavior. The current
focused resident evidence is recorded below.

Fo `17db6d81e0cf86b2325da159aa62964755eb92e6` additionally repairs public
execution from the selected external CMake build root. The quantum slice uses
`/var/tmp/fo-cmake-exec-root-20261005/build/gfortran_B54CC5252B24789C/app/fo`,
SHA256 `23f14fb2e36f0e999d4e6db9a421f57743fe477921c786f882be945d4a6c7aa6`.
Its independent stale-binary fixture, old/new profile fixture, and kin6d
Release/Debug evolution gates pass. Existing resident lanes are preserved.

## Focused development workflow

Fo main `5bdc786` supports focused CMake target builds, label/layer selection and
resident capture of KIN6D's populated dependency sources. On KIN6D `a93e74b`,
the verified driver SHA256 is
`2e6392ba836aabf9489c85dec6971f77247429536521121e9b69b47532e37333`.
The named `real_value_identity_edge_cases` case passes through both ordinary
Fo and a resident lane. Its ordinary plan selects `kin6d_core`, `kin6d_fft`,
`kin6d_intervals`, `kin6d_ode` and `test_real_identity`, without widening.

```sh
fo test real_value_identity_edge_cases -L '^core$' --layer kin6d_core
fo gremlin start --lane development --target real_value_identity_edge_cases \
    --random 0 -L '^core$' --layer kin6d_core
```

The pre-move edit to `src/real_identity.f90` (now
`src/core/real_identity.f90`) is selected by `--only-changed` with that
source as its reason. A deliberately incorrect identity is rejected by the
registered test, a syntax error is rejected at build time, and restoring the
source passes. These edit/reject/restore checks use ordinary Fo. The resident
repeat was blocked by an unrelated ARES lane holding the host admission slots;
it does not establish resident edit/reject/restore or complete scope coverage.
Run only the affected named cases. Never run the full test suite.

## Stack and ownership

| Responsibility | Owner |
|---|---|
| Plasma application/model composition | KIN6D, modern Fortran |
| Symbolic local mathematics | FortSym |
| Spatial FE/FEEC/IGA/DG/FCI/BEM mathematics | FortFEM |
| Generic shock/front FV/AMR, rarefied/DSMC/MD or radiation-transport primitives | Existing generic owner where suitable; otherwise establish an upstream owner only when a concrete reusable implementation is required |
| Generic numerics and rigorous arithmetic | FortNum |
| Sparse algebra/factors | FortSparse |
| Independent/source AD tooling | FortAD |
| Scientific I/O and NetCDF-4/HDF5 subset | FortIO |
| Experimental reconstruction/inference | TIAGO (external consumer) |
| Plotting | FortPlot |
| Low-dimensional labeled analysis arrays | FortArray |
| Nonlinear/linear solver orchestration | PETSc initially, replaceable |
| Build/test | CMake/Ninja/CTest |
| Development feedback | Fo |
| Accelerator leaves | native CUDA behind thin C ABI |
| Optional formal theorem/checker layer | Lean |

No Python runtime, generator, plotting or postprocessing layer is required for
the scientific executable or certifier. Optional Python may provide notebook or
web-server convenience and TypeScript may provide a browser UI, but both sit
above the Fortran-owned configuration/schema and result contract. External
Python tools may remain comparison/oracle software.

## FortSym is the mathematical source of project-owned leaves

Project-owned scalar/tensor numerical formulas should be retained symbolically
and emitted by FortSym. This applies to:

- primal local physics;
- JVP/VJP leaves;
- strict/reference arithmetic;
- CUDA leaves;
- interval/enclosure leaves where supported;
- generated constants and small contractions.

Handwritten Fortran owns state, schedules, solver composition, memory,
communication, I/O and calls to generic libraries. Do not hide new numerical
formulas in orchestration glue merely to bypass generation.

Generated code is never edited by hand. Retain the generator, assumptions,
target, dependency revision and an independent behavioral oracle.

The rule does not require rewriting BLAS, PETSc, FortNum, FortFEM or other
generic algorithms in FortSym. Generic missing capabilities belong upstream.
Do not force every new numerical family into FortFEM merely because it is the
current spatial foundation: conservative finite-volume/AMR, DSMC/MD or
radiation-transport infrastructure should live at the generic owner that
emerges from real reusable consumers.

## Three layers of generated execution

Keep these distinct:

1. **symbolic mathematical IR** — expressions, tensors, exact derivatives and
   identities;
2. **indexed operator IR** — map/reduce, gather/scatter, sparse interaction and
   tensor contraction structure;
3. **schedule** — tiling, SIMD/warp mapping, shared-memory staging, fusion,
   batching and data layout.

Changing a GPU tile must not change the mathematical definition.

For operations such as

\[
y_i=\sum_jF(x_i,y_j,p),
\]

the indexed form is authoritative and the \(i,j\) matrix is virtual. CPU/CUDA
schedules consume tiles and reduce them in fast memory.

Ahead-of-time specialization is preferred initially. Generate a bounded family
of schedules/orders, benchmark on the target hardware and select the measured
winner. Runtime JIT is not required.

## Derivatives

Every continuous production operator provides primal/JVP/VJP from one
mathematical definition. HVP is added when an actual consumer exists.

FortAD is an independent source-transformation oracle/fallback. Finite
differences, complex step and analytical identities remain independent checks.
LLVM/Enzyme is optional and never a required production dependency.

Implicit operations are differentiated mathematically:

\[
R_U\,dU=-R_\theta\,d\theta,
\qquad
R_U^\dagger\lambda=Q_U^\dagger.
\]

Do not tape Newton/Krylov iterations when the implicit derivative exists.

For reverse kernels, recompute cheap tile-local intermediates instead of
keeping a global AD tape.

## Rigorous arithmetic and Lean

The rigorous production path is Fortran/FortNum verified arithmetic.
Workhorse evaluation uses ordinary centers with proved aggregate arithmetic
and approximation remainders where supported. Strict interval/ball evaluation
remains the reference and fallback. Both require conservative containment
under their execution assumptions, not exact or last-bit-tight endpoints.
A certificate may be mathematically rigorous without Lean.

Lean is an optional final checker for reusable theorems and especially
high-assurance results. It is not a dependency of ordinary solver builds or
early continuous-certification milestones.

Strict certifier compilation follows the verified arithmetic requirements:
no unsafe reassociation or unaccounted contraction. Fast candidate builds may
use different flags because their results are untrusted until checked.

Each fast certified primitive declares supported domains, schedule, rounding
mode, contraction/reassociation, subnormal handling and overflow disposition.
Ordinary libm accuracy or an empirical percentage inflation is not a proof.
Unsupported domains fall back or return explicit certification failure.
Arithmetic budget is checked after stability/adjoint amplification; precision
can rise without changing the mathematical model or output claim.

Upstream comes first: maintain native FortNum reduction and FortSym generated
primal/derivative/enclosure fixtures with independent oracles and raw timing
samples before optimizing their schedules. Use existing interfaces rather than
building a universal ball compiler. Add flux/reconstruction support to FortFEM
only when the forced-FP certificate provides its concrete consumer.

## FortFEM boundary

FortFEM owns generic mapped spaces, tensor B-splines/NURBS, FEEC, DG/HDG,
traces, mortars, FCI and boundary operators. KIN6D supplies plasma constitutive
laws and model composition.

For production hot paths, add generic matrix-free/sum-factorized/many-RHS
actions upstream when needed. Retain assembled sparse paths as small reference,
coarse-solver and debugging/certification tools.

Do not fork or copy a generic FortFEM kernel into KIN6D.

## Particle-runtime boundary

PIC/PIF/blobs initially live as KIN6D representation implementations so the
correct abstraction is discovered from real physics cases.

Generic reusable primitives may later move upstream, for example:

- AoSoA particle storage;
- deterministic/QMC initialization;
- sorting/binning;
- migration;
- deposition/gather;
- particle-in-basis/NUFFT primitives;
- blob/remapping infrastructure;
- wall-intersection acceleration.

Do not create a speculative separate particle package before multiple
independent consumers require the same contracts.

## CUDA boundary

The CPU path remains an executable reference.

Native acceleration is

\[
{\rm Fortran}
\rightarrow
{\rm ISO\_C\_BINDING}
\rightarrow
{\rm thin\ extern\ C\ wrapper}
\rightarrow
{\rm generated\ CUDA/library}.
\]

CUDA-specific types do not enter plasma physics modules.

The interface declares layout, dimensions, precision, ownership, residency,
stream ordering and error status. Candidate state stays on device across
operator calls. Add CUDA only for a measured hot path.

cuFFT/cuBLAS/cuSPARSE, CUDA graphs, asynchronous copies and CUDA-aware MPI are
implementation tools, not architectural requirements.

## PETSc boundary

PETSc may initially own SNES/KSP/PC/TS orchestration. KIN6D owns

\[
R(U),\qquad Jv,\qquad J^\dagger v,\qquad P^{-1}v,\qquad U.
\]

Do not let physics arrays or model definitions depend on PETSc internals.
This keeps replacement possible if a bespoke solver later becomes worthwhile.

## Build and verification

CMake/Ninja is authoritative for mixed Fortran/CUDA builds and CTest for
registered tests. Fo provides focused/resident development feedback but does
not redefine the scientific test contract.

The `cpu` and runtime-checked `debug` presets run the numerical gates without
Lean. The explicit `formal` preset also checks the retained theorem package and
Fortran-produced rational certificate. This changes dependency selection, not
the claimed strength of a numerical result.

MUMPS (through FortSparse, threaded with OpenMP) is required, and direct
solves use FortSparse's AUTO backend: SuperLU below 40000 stored entries,
where its factor + solve is faster, MUMPS above (measured crossover in the
FortSparse README). Configure fails unless MUMPS is found through
`MUMPS_DIR` or `CMAKE_PREFIX_PATH` (an OpenMP scivision/mumps install) or as a
distribution package, and warns when the found library has no OpenMP. The
`superlu` preset (`KIN6D_ENABLE_MUMPS=OFF`) builds the serial SuperLU path for
every direct solve; the `mumps` preset adds a focused direct-solver test
filter. A build tree configured before MUMPS became the default may still cache
`KIN6D_ENABLE_MUMPS=OFF`; reconfigure it with `-DKIN6D_ENABLE_MUMPS=ON`. System
MPI builds are supported: KIN6D initializes MPI with serialized thread support
on first factorization if needed, each MUMPS instance runs on
`MPI_COMM_SELF`, and an already finalized or insufficiently thread-capable
runtime is rejected through the solver's ordinary failure status. Applications
that manage MPI must free all solver handles before finalizing it. Sequential
MUMPS builds retain their mpiseq runtime and need no Fortran MPI module.

FortSparse's default fetched SuperLU shares the selected system BLAS/LAPACK
with MUMPS. Existing build trees that cached `FORTSPARSE_USE_SYSTEM_DEPS=ON`
must unset it (`cmake -U FORTSPARSE_USE_SYSTEM_DEPS ...`) to use this default.
An explicit system-dependency choice remains supported, but prebuilt solvers
must share the selected BLAS to avoid mixing implementations in one process.

The focused resident workflow above captures and builds the real CMake-only
consumer, including populated dependency sources; the historical capture
blocker was not reproduced. Native CMake/CTest remain authoritative. Keep one
owned resident lane, preserve other owners' lanes and respect host admission.
Concurrent builds in another build directory can produce a large watcher-event
backlog; keep that independent build outside the watched checkout.

Every generated/operator family keeps:

- independent behavior test;
- derivative pairing/directional tests;
- strict/reference comparison;
- performance case when hot;
- certificate test where a rigorous path exists.

Regeneration equality is provenance evidence, not an independent correctness
test.

## Upstream repair policy

When profiling or correctness exposes a generic defect in a controlled
dependency:

1. retain a minimal consumer reproducer;
2. identify the owning library;
3. add the generic upstream test;
4. implement the smallest complete fix upstream;
5. verify and promote it according to that repository's policy;
6. pin/recheck KIN6D.

Do not accumulate permanent downstream workarounds or duplicate optimized
kernels.

### libneo: replaceable boundary, actively improved implementation

KIN6D uses libneo behind a narrow adapter/canonical-data boundary so historical
module structure, global state and compiler-specific implementation details do
not become KIN6D architecture. That boundary is for replaceability and
containment, not to freeze libneo in its current form.

The first production use of each libneo capability by KIN6D is an upstream
hardening gate. Audit the exercised path and, when the issue is generic, leave
libneo better than it was:

- fix convention/sign/unit and interpolation defects upstream;
- remove or contain global state when it blocks reentrancy, threading,
  differentiability or independent checking;
- expose clean generic spatial derivatives/JVPs when several consumers need
  them;
- add accuracy/residual/interpolation diagnostics and reusable error/enclosure
  support when the mathematics belongs to the imported representation;
- retain simple value-only compatibility for old consumers where practical;
- add analytic/manufactured/cross-code oracles before replacing historical
  behavior.

Do not require all historical field implementations to implement one huge new
interface. Add optional/capability-specific derivative or enclosure interfaces
only when concrete consumers need them.

If making an old evaluator differentiable or rigorous would be disproportionate,
use libneo to parse/normalize the external data and reconstruct a canonical
KIN6D basis/FEM/spline representation. KIN6D then differentiates/certifies that
representation. This is preferred to copying the old evaluator into KIN6D.

The ownership split is:

- libneo: generic file/convention transforms, reusable field/geometry
  evaluation, generic spatial derivatives and representation-level error
  diagnostics/enclosures;
- KIN6D: physical-model residuals, implicit equilibrium/steady-state tangents
  and adjoints, validated solution existence/error, and physical model error;
- TIAGO: diagnostic likelihood/posterior UQ.

Upstream libneo improvements are intentionally ecosystem improvements: existing
SIMPLE, MEPHIT, NEO-2, TIAGO and other consumers should be kept compatible and
benefit from corrected generic behavior as they update their pin/version.

## Configuration and user-facing contract

Fortran owns the typed configuration model, compatibility checks, defaults and
automatic numerical choices. Fortran namelist is the canonical editable text
format. Preserve both the user's literal namelist and a fully resolved namelist
with all defaults and automatic choices materialized. Generate machine-readable
schema information from the same definitions rather than maintaining a second
Python/JavaScript schema by hand.

The CLI is the reference operator interface and should grow small orthogonal
commands such as validate, resolve, run, inspect, extract-input and doctor as
real consumers appear. A browser UI may synchronize generated forms with a raw
namelist editor, submit jobs through local scheduler adapters and display
results/certificates. It must not acquire independent scientific semantics.

Scheduler/resource configuration is not plasma physics and stays outside the
physics namelist.

## Scientific result format

Production scientific output uses NetCDF-4 semantics on the HDF5 subset
supported by FortIO. Do not use arbitrary HDF5 features that make the file
unreadable to NetCDF tooling. Named dimensions, coordinate variables, units and
descriptive attributes are preferred.

Every result embeds the literal input namelist, the fully resolved namelist,
KIN6D/dependency revisions and build configuration, compiler
identity/version/flags and strict-vs-fast mode, material hardware/runtime
information, and normalized metadata/hashes for external inputs.

Large external equilibria/databases need not be copied into every result when a
content hash and sufficient convention/provenance metadata identify them. The
invariant is that the effective KIN6D configuration can always be recovered
from the result.

## Compiler portability

Supported production compiler families are GNU gfortran, Intel LLVM ifx,
NVIDIA nvfortran and LLVM Flang. CMake must select or probe family-specific
warning/optimization/strict-FP flags; GNU-only flags are never global project
policy. Fast candidate and strict/reference/certifier targets may deliberately
use different flags.

CI/support claims mean the relevant tests actually compile and run with the
named compiler. FetchContent remains useful for development, but release/offline
cluster workflows should not require outbound network access merely to build a
pinned release.

## Reproducibility/provenance

Record exact source/generator/dependency revisions, compiler and flags,
hardware/driver, mathematical model, normalization and external-data hashes.
Fast and strict builds are separate named configurations. Result-file
provenance is the primary run record; source history need not contain generated
run products.

Generated images, animations and large performance datasets are regenerated
artifacts and do not belong in source history.

Chris&AI
