# Symbolic operators, derivatives and bounded memory

Status: required design policy and staged upstream work. No kin6d kernel,
adjoint, lazy-reduction backend or sensitivity certificate is implemented.
Readers implementing the next operator should preserve its mathematical
definition while choosing only the execution machinery that case needs.

## Mathematical source and execution plans

FortSym owns the mathematical definition of every project-owned numerical
leaf. Each continuous production operator provides primal, JVP and VJP from
that same definition; HVP/second variations are added for a concrete consumer.
The definition includes assumptions, units, constraints, quadrature/weights,
boundary actions and parameter dependencies. Handwritten Fortran orchestrates
arrays, nonlinear iterations, tiling, MPI, device placement and I/O.
Execution control is not part of the symbolic physics definition.

The target pipeline separates three responsibilities:

| Layer | Meaning | Owner |
|---|---|---|
| Mathematical expression DAG | Exact expressions, tensor algebra, domain assumptions and differentiation | FortSym; plasma-specific definitions remain in kin6d generators |
| Indexed operator IR | Maps, reductions, segmented reductions, gathers/scatters and tensor contractions, with explicit dependencies | Generic FortSym generation; FortFEM supplies spatial operators |
| Execution schedule | Layout, tiles, SIMD/warps, batching, fusion and memory movement | Generic backend schedules upstream; application composition in kin6d |

Generate ordinary Fortran CPU/reference and optional native CUDA leaves,
with strict and rigorous evaluators from the shared definition. Compiler AD
and runtime JIT are not requirements. A bounded ahead-of-time variant family
can be compiled and tuned by hardware/compiler/operator/shape. Store the
selected executable and generation/schedule provenance.

Local FortSym source at `a66f3756585752fb3e3aa9cf4d27646919656734`
contains a [typed kernel IR](https://github.com/lazy-fortran/fortsym/blob/a66f3756585752fb3e3aa9cf4d27646919656734/src/codegen/fortsym_kernel_ir.f90)
lowered from an expression DAG. Its [roadmap](https://github.com/lazy-fortran/fortsym/blob/a66f3756585752fb3e3aa9cf4d27646919656734/ROADMAP.md)
still lists broad indexed/implicit/tensor differentiation as incomplete.
This inspection does not establish the full pipeline above. Verify each
required target and implement a missing generic slice at its owning upstream
repository, with an independent behavioral test and consumer recheck.

FortAD supplies independent derivative generation/checks and supported
parser-level fallback. Any production numerical leaf still passes through
FortSym emission. Enzyme may be an optional additional LLVM comparison.
Neither defines the production mathematics or becomes a required compiler.
The [FortFEM roadmap](https://github.com/lazy-fortran/fortfem/blob/2d7baa19fe1bef02ce279fa407830c8731bdc663/ROADMAP.md)
provides a compatible residual/JVP/VJP direction; planned contracts are not
evidence that all upstream physics/derivative paths already exist.

## Lazy indexed reductions

For `a_i = sum_j w_j F(x_i,y_j,theta)`, generate a tiled evaluation that
consumes each pair contribution immediately. Keep accumulators and reusable
tiles in registers/cache/shared memory; never persist the complete logical
pairwise matrix or its Jacobian. Inputs, outputs, bounded scratch and any
partial reductions remain explicit memory costs. Tiling reduces intermediate
storage, but does not reduce an intrinsically O(M N) arithmetic workload.

[KeOps map-reduce](https://www.kernel-operations.io/keops/engine/map_reduce_schemes.html)
is the concrete reference for evaluating/reducing pairwise formulas without
a device-memory kernel matrix. Its different schedules can require partial
buffers; the contract bounds those buffers rather than promising zero scratch.
[FlashAttention](https://arxiv.org/abs/2205.14135) motivates measuring movement
between accelerator memory levels; [FlashAttention-2](https://arxiv.org/abs/2307.08691)
illustrates the effect of work partitioning on the same mathematical task.
These are architectural examples, not plasma speedup claims or dependencies.

Use this pattern for appropriate nonlocal/BEM, kinetic-reaction, angular
coupling and quadrature reductions. Prefer superior physical algorithms:
harmonic Rosenbluth radial solves remain the proposed Landau route rather
than replacing them with a quadratic direct integral. A direct tiled kernel
may still be an independently derived small reference. Singular kernels need
regularization/self-term and quadrature contracts before tuning.

The symbolic family also generates reverse reductions. Cotangents for `x_i`
reduce over `j`; those for `y_j` reduce over `i`. Weights, parameter reductions,
normalizations and complex conjugation follow the declared pairing. Gather
transposes become scatter/reductions with an explicit collision/ownership
strategy. Differentiating a reduction does not license an unbounded tape,
hidden dense Jacobian, racing writes or omitted MPI reverse communication.

## Recompute within tiles

Default to no global derivative tape for local generated operators. Recompute
cheap forward intermediates inside VJP tiles; checkpoint an expensive reusable
quantity only when measured arithmetic/traffic costs justify its storage.
Recomputation is a schedule choice with numerical roundoff obligations, not
a change in the physical expression. Record forward and reverse peak scratch,
register pressure/spills, recomputation work and bytes transferred.

Fuse tensor interpolation, local physics and quadrature accumulation where
measured reuse warrants it. Small 1D basis tables and radial factors may be
stored; full tensor-basis or pairwise intermediates need an explicit reuse
justification. Compare legal fusion/CSE strategies: fewer FLOPs can require
more live registers and cause spills. Generic generation eventually records
operation counts, temporary lifetimes and candidate FMA/live-scalar estimates.

The [CUTLASS overview](https://docs.nvidia.com/cutlass/latest/overview.html)
provides a model for separating operation/layout from hardware decomposition.
Asynchronous tile pipelines, persistent blocks and graphs are later schedule
options for a profiled bottleneck. They do not justify a required C++ template
framework or device-specific formulas in the plasma model.

## Implicit derivatives of solved equations

For a stationary branch `R(U,theta)=0`, differentiate the equation:
`R_U dU = -R_theta dtheta`. For a real scalar objective, solve
`R_U^T lambda = Q_U^T` and obtain
`dQ/dtheta = Q_theta - lambda^T R_theta` in coefficient dual coordinates.
A physical weighted adjoint uses the corresponding Riesz/weight maps;
complex formulations require the declared conjugate pairing.

This requires a differentiable branch and an invertible constrained/gauge-fixed
Jacobian. Constraints and parameter-dependent domains/boundaries belong in
the differentiated residual. Near bifurcations, nullspaces or nonunique roots,
report failure of these assumptions rather than a fictitious sensitivity.
Primal/adjoint residuals and conditioning determine derivative accuracy.
Do not differentiate Newton/Krylov histories as the production definition.

Linear solves use `A dx = db - dA x` and a transpose solve in reverse mode;
factor derivatives need not be taped. Simple resonance roots use implicit
root derivatives only when the root is isolated and its derivative nonzero.
Eigenvalue derivatives require simplicity, normalization and phase/gauge
conditions; degeneracy may require a subspace formulation. Differentiating
the equations avoids dependence on accidental iteration counts and tolerances.

Transient adjoints still require trajectory access. Use a declared checkpoint/
recompute schedule, deterministic replay requirements and memory/I/O budget.
The local no-tape policy does not mean a long nonlinear trajectory's adjoint
can run without saved states or recomputation. Distinguish a discrete time-step
gradient from a validated continuous-time sensitivity.

## Validity and independent gates

Each differentiable operator declares its smooth domain, fixed-topology
assumptions and event margins. Wall-hit/loss changes, mesh/cut topology,
mode activation, resonance collisions and limiter branches can be nonsmooth.
Freeze the relevant topology for an ordinary derivative or provide explicit
event/saltation/one-sided semantics with their hypotheses. Otherwise return
an unavailable derivative. Optimizer surrogates are separately labeled.

Compare optimized/reference primal actions and derivative actions, then test
against an independent analytical derivative or FortAD route. Directional
finite differences need a step/roundoff study; complex step requires a
holomorphic implementation and cannot validate arbitrary abs/min/events.
Check adjoint duality with the actual weights, orientations and parallel
ownership. Agreement or duality alone can preserve the same erroneous
formula; it is necessary evidence, not the only derivative oracle.

A rigorous sensitivity enclosure requires a certified primal branch,
derivative enclosures, verified constrained forward/adjoint inverse bounds,
parameter and discretization/tail control, and objective regularity. The
shared DAG makes consistency manageable but does not prove its derivation,
generated compiler semantics or continuous-model validity automatically.
Lean-checkable sensitivity claims remain a later explicit theorem milestone.

## Physical state capacity

Keep one canonical packed distribution, plus accounted solver vectors and
bounded tile scratch. Avoid separate global transport/collision/FFT layouts.
Lazy operators remove auxiliary matrices; they do not remove the O(N_state)
physical distribution, residual or necessary work vectors.

Prefer contiguous owned blocks per resolution class. Block pools are useful
for activated patches, adaptive harmonic/radial classes or rare populations;
they need explicit ownership, stable IDs, transfer and fragmentation budgets.
[PagedAttention](https://arxiv.org/abs/2309.06180) supplies an allocation example;
[vAttention](https://arxiv.org/abs/2405.04437) motivates comparing the layout
overhead with contiguous alternatives. Their serving results do not predict
kinetic-solver performance.

When state capacity exceeds one accelerator, distribute configuration-space
ownership with full local velocity state before defaulting to host/NVMe
streaming. Other decomposition is permitted when a measured capacity/traffic
constraint demands it. Any out-of-core path declares working-set, bandwidth,
prefetch, fault and checkpoint/restart behavior. Count transfers end to end;
no automatic unified-memory paging is assumed to make an infeasible run cheap.

## Staged adoption and performance

Start with mandatory generated primal/JVP/VJP for the next continuous leaf.
Add a generic indexed reduction only when a client needs it, then its CPU
and optional CUDA schedules, independent reference and bounded-memory tests.
KeOps comparisons can be optional external audits of matching formulas and
accuracy; they introduce no Python runtime into kin6d. Apply the same generic
machinery to FortFEM tensor/quadrature and collision reductions where useful.
Strict/rigorous emission and derivative domains accompany each supported
operator; unsupported enclosure/AD targets fail explicitly.

Add implicit solve sensitivities when an actual diagnostic/inverse/design
question needs them. HVPs, certified gradients, persistent schedules and
out-of-core pools follow their concrete consumer and theorem obligations.
This sequence does not postpone R0 behind an entire new symbolic compiler.

Extend the [performance objective](performance.md) with time to a requested
gradient tolerance, stating validated versus rigorously certified status.
Record forward/JVP/VJP time and memory, adjoint solve/checkpoint costs and
the resulting gradient error. Tune a Pareto vector over these returned
products; a faster primal with a much slower VJP need not win an inverse task.

Chris&AI
