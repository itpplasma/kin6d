# Spatial foundation and performance design

Status: proposed numerical policy. No kin6d numerical backend is implemented.
The initial reference case remains small; production policies below take
effect when their corresponding operators are introduced.

## Ownership and checked FortFEM snapshot

FortFEM is the preferred spatial mathematical foundation. kin6d owns plasma
models, species/bulk ownership, collisions, reactions, closures, physical
boundary laws, machine-data adapters and certification claims. FortFEM owns
generic spatial maps, spaces, weak operators, orientations, traces, interfaces
and discretization actions. FortSparse/PETSc supply sparse/large solves;
fortnum/fortsym supply generic numerical/spectral algebra and generated leaves.
FortFEM must not become a plasma application.

The checked main revision on 2026-10-05 is
`2d7baa19fe1bef02ce279fa407830c8731bdc663`. Its
[scope roadmap](https://github.com/lazy-fortran/fortfem/blob/2d7baa19fe1bef02ce279fa407830c8731bdc663/ROADMAP.md)
separates reusable numerical ingredients from application closures. Source
inspection confirms B-spline de Rham complexes, NURBS surface/volume maps and
geometry derivatives in
[bspline_feec.f90](https://github.com/lazy-fortran/fortfem/blob/2d7baa19fe1bef02ce279fa407830c8731bdc663/src/elements/bspline_feec.f90).
The [coverage ledger](https://github.com/lazy-fortran/fortfem/blob/2d7baa19fe1bef02ce279fa407830c8731bdc663/doc/roadmap_mephit_feec_bem.md)
distinguishes implemented paths from remaining work.

The checked [3D spline assembly](https://github.com/lazy-fortran/fortfem/blob/2d7baa19fe1bef02ce279fa407830c8731bdc663/src/assembly/assembly_bspline_3d.f90)
constructs CSC/local assembly data and repeatedly calls a basis evaluator
whose checked source allocates whole-basis work/output arrays. This motivates
a span-local, sum-factorized batched path. This inspection is not a FortFEM
build/test or throughput result; it does not prove that every other upstream
path lacks such optimization. Verify the exact API and current upstream
owner before implementing a missing generic slice.

## Governing rules

Exploit separability before parallelism. Use regular mapped patches where
suitable, while retaining a general-geometry escape path. Store states and
compact factors rather than a full phase-space operator. Move data only when
an operator/ownership change requires it. Every approximation must have an
explicit defect in the eventual certificate; mathematically equivalent
accelerations still have floating-point and implementation obligations.

These rules constrain the production path without prescribing day-one AoSoA,
MPI, CUDA, AMR, or a universal execution-plan framework. Keep the smallest
runnable layout until a measured need justifies the next change.

## Compatible spaces and geometry

Adopt the discrete de Rham structure for Maxwell/compatible magnetic fields:
`H1 -> H(curl) -> H(div) -> L2`, with electric fields in H(curl) and magnetic
flux in H(div) where the chosen formulation requires them. Validate incidence,
orientation, Piola maps, commuting projections, boundary conditions and
charge/current transfer. `div_h curl_h=0` is an algebraic/topological identity;
it does not bound continuous field error or guarantee a correctly initialized
Gauss law. Floating-point residuals and gauge/harmonic nullspaces remain.
A different potential formulation needs an equivalent compatible complex.

Fields, fluid variables and distributions need not use identical spaces.
Smooth equilibria/fields may benefit from spline/IGA continuity; kinetic
transport near sharp gradients and walls may prefer broken spline/DG/SEM
spaces. Knot multiplicity can vary continuity within a spline family, but
broken spaces still need their own fluxes, traces and compatible transfers.
Shared physical geometry does not make these transfers automatic.

Prefer mapped tensor-product patches where suitable, including cylindrical
fusion charts, G-frame Figure-8 patches, Cartesian/spherical space blocks.
Mortar, fitted/cut or unstructured pieces handle actual geometric/topological
needs. Do not mandate tetrahedra everywhere or assume one tensor patch can
cover every wall/X-point domain. NURBS maps and other embeddings must pass
[geometry](geometry.md)'s regularity/invertibility certificate obligations.

FortFEM supplies FCI interpolation/support actions and generic trace/interface
machinery. The application supplies the magnetic law, tracer choice, map
validity and terminal physical boundary law, reusing generic upstream tracing
where appropriate. FCI is not the physical coordinate system. BEM/DtN/PML,
IGA multipatches, cut enrichment and HDG are later options for a concrete
vacuum/open-domain/interface need, each with its own error treatment.

## Tensor actions and batches

For suitable charts, the full state has spatial, species-specific radial
momentum, angular-harmonic and optional periodic factors. This factorization
of representation is not a low-rank approximation to the distribution.
Never form the full six-dimensional tensor-product matrix or global kinetic
Jacobian. Apply low-dimensional factors and local variable-coefficient
couplings; not every nonlinear/curved operator has an exact finite Kronecker
sum. Retain compact factors, quadrature data and necessary sparse coarse/
field blocks, rather than interpreting “store states” as banning all matrices.

For tensor degree `k` and quadrature proportional to `k`, dense 3D local
application costs roughly O(k^6), while sum factorization can give O(k^4)
work with O(k^3) state. Curvature, vector components, interfaces and nonlinear
coefficient evaluation affect the constants and applicability. Benchmark
end-to-end certified accuracy, not this asymptotic estimate alone.

The first needed generic optimization belongs upstream in FortFEM: a
span-local evaluator with caller-owned scratch, precomputed 1D basis/
derivative tables, successive tensor contractions, and partial assembly of
compact quadrature geometry/Hodge data. No global spline basis or element
matrix is built in the hot application path. Cache only what improves the
measured memory/time tradeoff and invalidate geometry-dependent data when
geometry changes. Preserve the existing assembled path as a reference.

The plasma workload applies the same spatial action to many momentum/angular
channels. Introduce a generic batched many-RHS action upstream when this
workload exists. Reuse basis/geometry data across channels and resolve runtime
configuration outside point/element hot loops. Mathematical forms compile
to specialized generated kernels; they are not dynamically interpreted at
every quadrature point. All new numerical leaves obey the fortsym policy,
including upstream additions and CUDA/derivative/enclosure targets.

A DG spatial space does not independently select a second transport backend.
Choose one coherent characteristic/weak transport scheme first, with its
reconstruction, mass, flux and boundary contracts. Add alternate schemes only
for a demonstrated need and retain an independent analytical oracle.

## Layout, modes and collision locality

Use one canonical packed kinetic layout, with candidate/residual/update and
necessary solver vectors accounted separately. When profiling requires it,
use tiled AoSoA/batched storage and tile scratch to serve spatial actions,
velocity-local collisions and transforms. Avoid permanent full-state copies
for separate collision/transport/FFT orientations. Species blocks may differ
in radial grid, scale and harmonic degree; they do not share an imposed
common electron/ion/alpha momentum grid.

Store only independent real-state Fourier coefficients where valid, handling
normalization, conjugacy and Nyquist modes explicitly. Complex perturbation/
eigenvectors may require a different representation. In cylindrical fusion,
`n=0` is axisymmetry only with the compatible rotating momentum/vector frame;
a generalized periodic chart's zero mode is not necessarily a symmetry.
Sparse convolution for tiny mode sets versus padded FFT for larger ones is
a later measured crossover, not two initial implementations. Both need
consistent projection/alias bounds, including nonpolynomial operations.

Spherical harmonics give sparse rotation and low-order moment extraction:
number in degree zero, Cartesian momentum/current in degree one, stress in
degrees zero/two, with appropriate relativistic radial weights. There is no
need to reconstruct all angular samples to compute these FK moments.
Reduced-model polarization/magnetization and collision potential products
still require their own actions. Low-degree Gaunt contractions versus angular
transforms are similarly measured alternatives with rigorous quadrature/tails.

For fixed grid, boundary condition and pair coefficients, harmonic radial
potential operators can reuse one factorization per species/grid/degree for
many spatial RHS. Moving grids, changed tails/boundaries or state-dependent
coefficients invalidate relevant caches. Landau and Braams–Karney equations
use their own consistent coordinate/scaling conventions; a velocity-space
Rosenbluth radial formula cannot silently be reused in relativistic momentum.
This batch opportunity does not remove singular/tail certification work or
cross-species grid transfers.

## Solves, derivatives and memory

Use PETSc nonlinear/Krylov infrastructure through its Fortran interface.
Production residuals and JVPs are matrix-free; exact mathematical derivatives
come from the operator algebra/fortsym and supported fortad paths. The computed
JVP has ordinary floating error and requires independently checked actions.
A weighted VJP uses its stated coefficient/physical inner product, not an
unqualified array transpose. Finite differences are diagnostics or explicitly
bounded fallback experiments, not the default production derivative contract.

Assembled tiny/coarse/field blocks are useful in preconditioners. Reduced DK/GK,
fluid/MHD, local collision, degree-block and mode-block actions may accelerate
an FK solve without changing its converged residual's requested physics.
Track transfers and avoid accidentally solving the preconditioner's equations
as the final model. Strong preconditioning has to be measured, not assumed.

Budget global vector count explicitly. A candidate, residual, update and
roughly a handful of Krylov/work vectors is a goal, not a universal guarantee.
Small restarted FGMRES, short-recurrence variants or nonlinear globalization
are decisions based on actual stability/convergence. Fifty 100-GB vectors
are not an acceptable unreported default. Tiny full-GMRES diagnostics remain
useful. Include preconditioner factors, collisional auxiliaries, padded mode
scratch, certifier and checkpoint memory in the budget.

Static condensation/HDG can remove local interiors when the chosen discrete
system admits the elimination. Exact algebraic elimination preserves that
system, but conditioning, local solve error, interface norms and reconstruction
still enter the certificate. It is not automatically easier to certify.

## Execution and retained evidence

When distribution is necessary, decompose primarily in configuration space,
keeping complete local velocity and preferably periodic mode sets together.
Collisions are spatially local; transport needs halos and fields need global
coupling. Fourier physical-angle nonlinear evaluation also needs its mode
set/transform communication. Move to velocity/mode decomposition only after
capacity or timing demonstrates its need. Regular patch resolution classes
can retain batching under adaptivity; transfer and class-selection error
must be controlled.

Overlap independent local collisions with transport halo exchange when valid.
CUDA batches use the thin ISO_C boundary. Keep state/work vectors resident
where a measured GPU path benefits, reduce diagnostics in situ and checkpoint
deliberately. Certification may run on CPU but data-transfer/staging costs
must be counted. No routine full-distribution output at every iteration.

Mixed precision, approximate preconditioners and nondeterministic candidate
reductions are allowed because the candidate is untrusted. The independent
certifier evaluates the exact reconstruction under its own rigorous arithmetic
contract. A failed certificate triggers refinement/re-solve or an explicit
candidate-only result; it never converts an uncontrolled approximation into
an accepted interval.

Retain assembled reference, matrix-free CPU and eventual generated CUDA
realizations for small common operators. Cross-check their actions, then use
independent polynomial/invariant/manufactured/analytical expectations; agreement
between three implementations sharing a mistaken formula is insufficient.
Keep derivative duality, orientations, boundary balances and multiple RHS
in the verification cases. Reference implementations also use generated leaves.

Local residual, map, quadrature and tensor-factor bounds can make certification
tractable without constructing an interval global matrix. Continuous stability,
nullspaces and unresolved spatial/momentum/harmonic tails remain the hard
proof obligations; matrix-free JVP access alone cannot establish an inverse
norm bound. Every optimization retains or supplies its needed error analysis.

Measure total time and peak memory to certified observable tolerance, with
solver, certification, communication, transfer and postprocessing breakdowns.
Add each generic performance slice at its owning upstream repository, retain
a focused independent reproducer, push a locally verified main fix promptly,
and recheck the consumer. This document authorizes no speculative library
redesign before a concrete client needs it.

Chris&AI
