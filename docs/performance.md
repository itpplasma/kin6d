# Continuous performance engineering

Status: required workstream, with no kin6d measurements yet. Optimize total
time and memory to the requested certified observable accuracy from the first
executable milestone. Track the best known frontier rather than accepting
successive small regressions relative to the previous commit.

## Objective and review before implementation

Minimize candidate plus certification time subject to the requested interval
width, available memory and the stated model/proof level. Report postprocessing,
I/O and transfers too. A faster kernel that increases nonlinear iterations or
widens the accepted enclosure can make the scientific calculation slower.
Different reduced models are comparable only when their model-error claims
and target observable are compatible.

Review model choice, representation, algorithmic complexity, data movement,
parallel decomposition and microarchitecture in that order. Before a hot
feature, identify minimum required state, separable dimensions, diagonal/sparse
algebra, unavoidable work/traffic/communication and reusable quantities.
Then implement a correct generated reference, supply an independent oracle,
measure a baseline, profile one bottleneck, change one thing, remeasure and
retain only a validated time/memory/accuracy improvement.

Each milestone has a small correctness case and a representative performance
case large enough to expose its actual polynomial/harmonic orders and cache/
memory behavior. R0 can use a modest CPU case; it does not require a GPU,
cluster or implementation of every future calibration kernel.

## Reuse benchmark infrastructure

Use existing numerical-library harnesses and fo orchestration rather than
building another general framework. FortNum's
[benchmark guide](https://github.com/lazy-fortran/fortnum/blob/main/benchmark/README.md)
documents machine-readable measurements, RSS evidence and controlled reference
gates. FortAD's
[vector-JVP design](https://github.com/lazy-fortran/fortad/blob/main/docs/design/vector-jvp-benchmark.md)
and FortIO's [performance guide](https://github.com/lazy-fortran/fortio/blob/main/docs/performance.md)
are starting points to review/pin when their clients are measured.
These references are not kin6d integration or throughput evidence.

Some upstream benchmark collection/gate scripts currently use Python. Do not
make them mandatory in kin6d's build, run, tests, plots or performance workflow.
Reuse native measurements and metadata where supported; if a needed native
gate is missing, implement the smallest reusable support at its owning library
or fo. Do not clone the same benchmark/gate logic privately downstream.
Generated calibration/reference numerical leaves follow the fortsym policy;
instrumentation and orchestration are not new plasma numerical kernels.

## Calibration and operator contract

Calibrate actual sustained behavior on each authorized target machine:
CPU/cache copy and arithmetic, small dense and batched banded solves, relevant
FFT sizes, and later GPU HBM/FP64/reductions/copies and transfer links.
Distributed cases additionally measure message latency/bandwidth, reductions
and representative halo exchanges. Record working-set size, affinity, NUMA
placement, clocks/thermal conditions, software and timing method. Vendor peak
is not a measured roofline. Do not automatically touch protected hosts or
schedule/install calibration jobs; runner access remains explicitly scoped.

An operator record contains input/output/scratch bytes, estimated operations,
minimum traffic under stated reuse assumptions, allocation count, synchronization
count, communication volume and arithmetic intensity. Production local hot
actions should allocate no dynamic storage and introduce no internal global
synchronization; record unavoidable setup/event exceptions separately.
FortSym work/temporary estimates can later supplement measured counters.

For calibrated bandwidth B and arithmetic throughput P, a roofline estimate
is `t_est=max(bytes_min/B,flops/P)`. Record `t_est/t_measured` as a diagnostic
with all assumptions. This is not a rigorous universal lower bound: latency,
cache levels, dependencies, measured workload mismatch and overlap matter.
A low ratio prompts investigation; a fixed ratio threshold does not prove
that an implementation is optimal or incorrect.

## Five measurement levels

| Level | Workload | Required record |
|---|---|---|
| 0 | Microkernel | Observable output, shape/order, latency, scratch/allocations, independent oracle |
| 1 | Complete operator | Representative arrays, bytes/work, throughput and peak memory |
| 2 | Residual | Inclusive/exclusive phase times, reconstruction and coupling costs |
| 3 | Nonlinear/linear solve | Residual/JVP/preconditioner calls, iterations, rebuilds, vector count and peak memory |
| 4 | Scientific calculation | Model/inputs, proof level, interval width, solve/certifier/I/O times and total memory |

Use cheap nested timers/counters for transport, collisions, reactions, fields,
boundaries, derivatives, preconditioning, linear solve, certificate and I/O.
Avoid double-counting overlapping/inclusive phases. Emit machine-readable
records from normal runs; heavier profiles are opt-in.

Find where time goes before investigating why. CPU compiler reports and
`perf`/available hardware counters follow timers; detailed sampling or code
inspection targets only important kernels. On NVIDIA use a system timeline
before detailed kernel counters; record registers, spills, memory throughput
and occupancy without optimizing occupancy as an end in itself. MPI profiles
include waits, reductions, message sizes, overlap and imbalance. Profiler and
instrumentation overhead are reported rather than treated as application time.

## Statistical evidence and baselines

Warm up, repeat and retain raw measurements, median and dispersion. Alternate
A/B order or randomize controlled ordering when thermal/drift effects matter.
Make workloads observable so the compiler cannot eliminate them. Use separate
processes where candidate-specific peak RSS must be measured. Noise larger
than the claimed gain means inconclusive evidence, not a speedup.

Maintain a best-known Pareto set of time, memory and achieved accuracy for
fixed benchmark/input, model/proof level, hardware, compiler/flags and shape.
Do not silently reset a baseline to a slower parent. Only the controller
promotes a new validated baseline with its exact revision and evidence.
A faster but much larger-memory candidate may be unsuitable. Changes in model,
solver tolerance, returned derivative product or output obligations require
a new labelled comparison, not an apparent optimization.

Pinned authorized runners provide hard timing gates; ordinary hosted runners
can provide correctness/gross signals but not reliable fine regressions.
No new nightly job is created by this document. Expensive audits run explicitly
or under a separately authorized schedule. Stable material regressions can
trigger a focused dependency/application bisect with frozen inputs.

## Implementation priorities and measured choices

[Spatial numerics](numerics.md) defines tensorization, partial assembly,
batched RHS, FEEC and one canonical tiled state. Additional production rules
are measured contracts, not a reason to implement every specialization now.

| Area | Priority and qualification |
|---|---|
| CPU hot loops | Caller scratch, contiguous ownership, simple strides, no per-point polymorphism; inspect array temporaries/vectorization |
| Small algebra | FortSym-generated contractions/specializations at measured common orders; retain generic/reference path |
| Generation | Compare contraction/CSE choices using FLOPs, live temporaries, code size and register pressure; fewer FLOPs may spill more |
| Geometry | Reuse fixed-map/basis tables, uniform spans and useful metrics; compare recompute cost with added bytes |
| Fourier | Reuse plans, batch transforms, independent real spectrum where valid; measured direct-triad/FFT crossover |
| Angular algebra | Packed valid modes, sparse coupling graph; measured Gaunt/transform crossover, low-harmonic moments |
| Collisions | Cache valid radial factors, batch RHS and species pairs by spatial block; conserving correction must preserve its stated algebra |
| Reactions | Parse/normalize once, sparse reaction graph, batched generated interpolation; pruning needs omitted-population bounds |
| Boundaries | Generic BVH/acceleration structure at its owner, batch intersections, rigorous event bounds near grazing |
| CUDA | Resident state, large regular batches, launch/sync accounting; fuse only while measured benefit exceeds register/spill cost |
| CUDA graphs/shared memory | Add only after stable measured DAG/reuse; neither is a universal improvement |
| MPI | Spatial ownership preserving collision locality, overlap valid local work with halos, count/reduce global reductions |
| CPU scale | Explicit NUMA first touch/affinity, regular scheduling, bandwidth saturation; huge pages only after evidence |
| Metadata | Checked narrow local indices where safe, sufficiently wide global identities/counts; reject overflow |
| I/O | In-situ moments, deliberate authoritative checkpoints, large writes; asynchronous/chunk/compression choices measured |

Retain the simple reference and optimized CPU path; add the generated GPU
realization only when a profiled workload justifies it. Backend fusion and
vectorization do not change the model. Verify independent physical/analytical
expectations as well as equivalent-path agreement; shared mistakes can make
all implementations agree. Compiler vectorization evidence is relevant on
a declared compiler/shape, not a portable universal correctness theorem.

Bounded compile-time autotuning may select tile/batch/block sizes, fusion,
specializations and transform crossover by hardware/compiler/kernel/shape.
Store the measured selection and exact executable; no runtime JIT is required.
Compiler tournaments, PGO/LTO, assembly/PTX/SASS inspection and TLB tuning are
secondary tools for an identified bottleneck. Keep physics free of hard-coded
SIMD widths or one accelerator's occupancy assumptions.

## Solver, certifier and arithmetic

Profile primal, JVP and VJP separately with matching returned products; ratios
are diagnostics, not universal constant limits. Batched directions are useful
when the consumer actually has several. Independent derivative/duality oracles
remain mandatory, and numerical leaves still pass through fortsym.

Count every state/work vector and factor. Attribute state, residual/update,
Krylov, geometry, preconditioner, collision, scratch, padding, host/device
staging and certification memory. Strong reduced-physics preconditioning,
inexact Newton/continuation and carefully selected low-memory solvers often
matter more than kernel throughput. Record failures and iterations; never
silently change the final model to obtain a faster solve.

Reuse static factors and proof data only while their exact dependency hashes
remain valid. Approximate/lagged factors may be candidate preconditioners;
they are not automatically valid certificate inverse bounds. Reuse across
parameter scans includes branch checks and changed-data invalidation.
Checkpoint frequency trades measured write/rebuild cost against lost work;
reconstructible caches need not become authoritative checkpoint data.

Optimize the certifier too: local interval residuals, strict arithmetic,
verified linear algebra, norms, tails and reusable geometry/basis constants.
Do not loosen proof obligations to improve time. The strict certified path's
rounding, contraction/subnormal and transcendental contracts are independent
of fast-candidate compiler flags.

Class I transformations are mathematically equivalent factorization/fusion/
elimination/storage changes; their floating-point conditioning/roundoff still
needs coverage. Class II candidate arithmetic includes mixed precision,
reassociation and fast approximations, allowed only as untrusted search.
Class III changes the represented mathematical action/state, such as quadrature,
compression, dropped modes, FCI interpolation or physics reduction: include an
explicit consistent defect/remainder. A strict residual alone does not bound
a missing continuous complement. Compression remains optional, with a proved
truncation norm and an affordable fallback when rank grows.

## Cadence and upstream ownership

Fast sentinels accompany ordinary development: independent correctness,
allocation/vectorization evidence where stable, and gross microkernel changes.
A dedicated audit covers calibrated operators, derivatives, memory and scaling.
Release/science runs cover actual physical-device/space inputs and certified
observables. No commit needs to wait for unrelated hours-long audits, but its
applicable scientific and promotion gates still have to pass.

Repair a bottleneck at the lowest actual owner: FortFEM spatial actions,
FortNum generic numerics/verified arithmetic, FortSparse factors/sparse algebra,
FortSym generation, FortAD derivative machinery, FortIO I/O, and libneo field
interpolation where used. kin6d owns layout, physical composition/scheduling,
preconditioning and distributed execution. PETSc/MPI/CUDA bindings remain
adapters; their object layouts do not define the parent model/state contract.
No duplicate downstream optimization, private numerical fork or unverified
main promotion. Preserve owners' unrelated work, retain a focused oracle,
push verified owned-upstream main repairs promptly, pin and recheck the client.

Every major science result reports candidate time, certificate time and peak
memory with the resulting interval. Future reports generate time/memory/accuracy
frontier plots and earlier-case overlays through Fortran/fortplot; images are
regenerated, never committed. Being near the best known frontier is a continuous
engineering aim demonstrated by evidence, not an unmeasured optimality claim.

Chris&AI
