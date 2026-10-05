# Tooling and kernel policy

Status: canonical implementation/tooling policy.

## Stack and ownership

| Responsibility | Owner |
|---|---|
| Plasma application/model composition | KIN6D, modern Fortran |
| Symbolic local mathematics | FortSym |
| Spatial FE/FEEC/IGA/DG/FCI/BEM mathematics | FortFEM |
| Generic numerics and rigorous arithmetic | FortNum |
| Sparse algebra/factors | FortSparse |
| Independent/source AD tooling | FortAD |
| Scientific I/O | FortIO |
| Plotting | FortPlot |
| Low-dimensional labeled analysis arrays | FortArray |
| Nonlinear/linear solver orchestration | PETSc initially, replaceable |
| Build/test | CMake/Ninja/CTest |
| Development feedback | Fo |
| Accelerator leaves | native CUDA behind thin C ABI |
| Optional formal theorem/checker layer | Lean |

No Python runtime, generator, plotting or postprocessing layer is required.
External Python tools may remain comparison/oracle software.

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

The existing KIN6D resident lane has a documented CMake-only input-capture
defect; native CMake/CTest evidence remains authoritative while that upstream
gap is open. Do not reset another owner's lane or claim complete captured-input
evidence from it. The exact selected Fo driver and current limitation are
recorded in the steady laboratory evidence.

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

## Reproducibility/provenance

Record exact source/generator/dependency revisions, compiler and flags,
hardware/driver, mathematical model, normalization and external-data hashes.
Fast and strict builds are separate named configurations.

Generated images, animations and large performance datasets are regenerated
artifacts and do not belong in source history.

Chris&AI
