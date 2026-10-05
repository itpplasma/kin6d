# Performance engineering

Status: permanent workstream. Performance is designed and measured from the
first executable, not postponed until feature completion.

## Objective

The useful final metric is not raw FLOP/s. It is approximately

\[
T_{\rm result}(\epsilon_Q)
=
T_{\rm candidate}
+
T_{\rm certificate}
\]

subject to the requested observable uncertainty and available memory.

The project goal is to beat the best available matched codes in **total time
to the requested result**, including all reported error bands and requested
derivatives. This is a target, not an established speed claim. Record setup,
candidate, reconstruction/residual, stability check, bound evaluation,
derivative solves and postprocessing separately. Report cold single-use and
amortized repeated-use costs; reuse only data with checked validity keys.

The default requested observable radius is 1% of a declared physical scale
with an absolute floor near zero. Bound-evaluation roundoff initially receives
1% of that radius after inverse/adjoint amplification. Measure achieved width
as well as time; a tighter reference run is not a matched 1% comparison.

Every reported output keeps its appropriate bound/status. Selected tight
goal bounds may guide refinement, while other outputs retain valid wider
enclosures. Do not save time by dropping bands or relabeling norm errors as
pointwise errors.

For stochastic cases, include the sampling/confidence target explicitly.

A faster kernel that increases nonlinear iterations, memory pressure or
certificate width may be an end-to-end regression.

## Optimization order

Always optimize in this order:

1. mathematical model/representation;
2. asymptotic algorithm;
3. data movement and state count;
4. parallel decomposition;
5. kernel schedule;
6. instruction-level microoptimization.

Do not attempt to rescue an intrinsically expensive architecture with CUDA.

Before implementing a hot feature, answer:

- What state is mathematically necessary?
- What is the lowest plausible work complexity?
- Which dimensions are separable?
- Which transforms make operators diagonal/sparse?
- What bytes must cross each memory level?
- What can be reused across species, modes, timesteps or Newton iterations?
- How many full-state work vectors are required?

## Machine calibration

For each representative machine record sustained, not advertised, quantities:

- DRAM/HBM bandwidth;
- cache/shared-memory bandwidth where measurable;
- FP64/FMA throughput;
- small dense/banded solve performance;
- FFT throughput at actual sizes;
- device/host and inter-GPU transfer;
- MPI latency/bandwidth/allreduce.

For a kernel with operation count \(F\) and unavoidable traffic \(B_{\min}\),

\[
t_{\rm lower}
=
\max(F/P_{\rm sustained},B_{\min}/BW_{\rm sustained}).
\]

Track efficiency relative to this bound where the model is meaningful.

## Five benchmark levels

### Microkernel

Basis evaluation, 3x3 transforms, angular coupling, radial solve, deposition,
PIF factor, reduction.

### Operator

Transport, collision, Maxwell block, Poisson solve, PIC deposit/gather.

### Residual/JVP/VJP

One complete application with realistic channel counts.

### Solve

Record nonlinear iterations, Krylov iterations, preconditioner applications,
work-vector count and rebuild costs.

### Scientific case

Time and memory to a physical observable/certificate target.

Never use a microkernel improvement as proof of an application improvement.

## Performance contracts

Hot operators should have a machine-readable conceptual contract containing:

- input/output bytes;
- scratch bytes;
- dynamic allocations;
- operation count where tractable;
- unavoidable/estimated traffic;
- synchronizations;
- communication volume;
- arithmetic intensity;
- primal/JVP/VJP timings.

For production hot paths, target zero per-call dynamic allocation and zero
avoidable synchronization.

## Baselines

Maintain a best-known Pareto frontier for fixed input, compiler and hardware:

\[
(\text{time},\text{peak memory},\text{accuracy/certificate width}).
\]

Compare to the best known result, not merely the parent commit, to prevent slow
performance erosion.

Noisy benchmark machines do not gate percent-level regressions. Dedicated
performance runners use warmup, repeated measurements, median/dispersion and
alternating A/B order where thermal drift matters.

External comparisons fix the model, boundaries, normalization, observables,
requested error type/width, precision, compiler, hardware and thread count.
Publish candidate-only and certified totals side by side. If a competitor
cannot supply a comparable certificate, report that gap and its independently
validated accuracy; do not infer a matched certified speed ranking. Include
requested forward/adjoint results in optimization workload comparisons.

## Certification benchmarks before optimization

The current bounded-OU laboratory checks against an explicitly known Gaussian
profile. It is an independent reference/performance case, not the general
production checker. Its 4096-cell strict certificate previously took about
425 times the direct solve on the local Ryzen 5950X. Those cells achieved a
radius near 3.37e-7, far tighter than a 1% request; the ratio cannot predict
general residual-based certification cost.

Before expanding physics, retain two upstream native fixtures:

1. FortNum: strict reductions versus FP centers with proved aggregate error
   bounds, including cancellation, gradual underflow and explicit overflow
   failure; preserve exact independent value oracles and raw timer samples.
2. FortSym: one shared definition emitted as primal/JVP/VJP and enclosure
   leaves, independent value/derivative oracles, and separately timed paths.

KIN6D's `scripts/benchmark.sh` sweeps the FP and general damped-transport cases
at several resolutions
and records requested/achieved radius, solver/certificate/acceptance costs,
state bytes and provenance. Correctness validation is outside timed sections.
Run these as explicit audits, not noisy timing gates. The implemented transport
certificate uses arbitrary represented-source residuals and the continuous
inverse bound, including implicit tangents, cell bands and goal intervals.
Its [measured ledger](general-transport.md) identifies interval reconstruction
as the first dominant cost. Paired interval sine/cosine at FortNum, optional
paired emission at FortSym and a private cached Fourier basis now share work
across state/tangent bands. Cache setup and owned bytes are reported alongside
direct and cached complete-pipeline samples; warm reuse is a separate result.
Reference-profile evaluation stays in the benchmark suite.

Profile strict exponential/reduction costs before optimizing them. At the
owner, consider batched range-reduced evaluation, proved approximation and
rounding remainders, reusable validated tables and controlled reduction trees.
Choose schedules by measured total cost and accepted width. A proposed GPU
schedule needs its own reduction/rounding theorem and hardware measurements;
CPU speedups do not establish GPU support.

## Continuous profiling

The application always records a cheap hierarchical time ledger, for example:

- residual;
  - spatial transport;
  - velocity transport;
  - collisions;
  - fields;
  - reactions;
  - boundaries;
- JVP/VJP;
- preconditioner;
- Krylov;
- communication/reductions;
- certificate;
- I/O.

Also record attributed memory:

\[
M_{\rm state}+M_{\rm residual}+M_{\rm Krylov}+M_{\rm geometry}
+M_{\rm preconditioner}+M_{\rm scratch}.
\]

External profiling escalates only as needed:

1. internal ledger;
2. compiler vectorization reports/perf/Nsight Systems;
3. counters/Nsight Compute;
4. assembly/PTX/SASS inspection for the true top kernels.

## CPU microoptimization

Hot CPU kernels should favor:

- caller-owned or fixed scratch;
- contiguous arrays and simple strides;
- compile-time/specialized common orders;
- explicit small contractions;
- SIMD-friendly tiles;
- precomputed invariants when reuse beats traffic;
- no polymorphic dispatch in inner loops.

Fortran array syntax is not automatically preferred when it creates hidden
temporaries.

Specialize common spline/order/harmonic cases but retain a generic fallback.
FortSym may emit several CSE/unrolling candidates because minimum FLOPs may
increase register pressure.

## Tensor/spline execution

For tensor-product bases:

- span-local nonzero basis evaluation;
- precomputed 1D tables;
- sum factorization;
- partial assembly;
- many-RHS batching.

Do not store a 3D basis matrix or global sparse operator when compact 1D factors
and geometry coefficients suffice.

Geometry itself is subject to a recompute-vs-store roofline decision. A cheap
metric expression may be faster to recompute than to load repeatedly.

## Fourier and angular execution

Use exact Fourier multipliers for derivatives. Nonlinear products choose
between sparse/direct convolution and FFT/pseudospectral evaluation at a
measured crossover.

Spherical-harmonic couplings use sparse Gaunt/angular graphs when sparse and
transforms when dense enough. Store only allowed coefficients.

## Collision execution

Rosenbluth/radial operators that depend only on species/grid/\(\ell\) are
factored once and reused over spatial/mode RHS.

Benchmark specialized/banded batched solves before generic sparse machinery.
Fuse low-order moment/conservation reductions only when doing so reduces total
traffic without causing register spill.

## Particle/PIC/PIF execution

Particle kernels use SoA/AoSoA storage and large batches.

PIC considers:

- sorting/binning amortization;
- deposition contention;
- charge-conserving current kernels;
- gather locality;
- migration buffers.

PIF/particle-in-basis uses tiled particle-mode map/reduce. Compare direct
\(N_pN_k\) execution with NUFFT only at real sizes; do not assume a crossover.

Blobs/remapping add reconstruction cost and must justify it through lower
particle count/noise or better accuracy.

Monte Carlo/QMC comparisons report statistical efficiency, not just
particles/s.

## GPU rules

Keep state resident on device.

Measure:

- HBM throughput;
- registers/thread;
- spills/local memory;
- occupancy;
- achieved compute;
- launch count;
- synchronization;
- transfer volume.

Do not optimize occupancy as an isolated metric. Shared memory and fusion are
used only when measured reuse pays.

Once the kernel graph stabilizes, CUDA Graphs or persistent execution may
reduce launch overhead. Asynchronous staged loads may overlap movement and
compute.

## Distributed execution

Default decomposition is configuration space with complete local velocity state.

Overlap spatial halo exchange with local collision/reaction work.

Count and reduce MPI collectives aggressively: fuse dot products/norms or use
pipelined methods when they improve actual solve time.

NUMA placement, first touch, thread/rank binding and huge-page/TLB effects are
production-scale concerns once state sizes justify them.

## Solver performance

Total steady cost behaves roughly as

\[
T\sim N_{\rm Newton}N_{\rm Krylov}T_{\rm JVP}
\]

plus preconditioning and field work.

Therefore solver/preconditioner quality normally matters more than a small
kernel speedup.

Reduced physics is encouraged inside preconditioners:

\[
P_{\rm FK}^{-1}
\approx
\text{DK/GK transport}
+
\text{local collisions}
+
\text{reduced field solve}.
\]

This does not change the target FK residual.

Krylov memory is explicitly budgeted. Large default GMRES restart dimensions
are forbidden when they multiply a huge kinetic state beyond practical memory.

## Mixed precision

The candidate path may use lower precision, tensor-core arithmetic or
approximate factors when convergence and final strict residuals remain good.

Certification arithmetic is separate and follows its rounding theorem.

Measure end-to-end benefit; a faster low-precision preconditioner that causes
many more iterations is not a win.

## Derivative performance

Record

\[
T_{\rm JVP}/T_{\rm primal},\qquad
T_{\rm VJP}/T_{\rm primal}
\]

and derivative scratch/state overhead.

Batched JVPs are useful when several tangent directions share primal data.
VJPs should recompute cheap tile-local intermediates rather than maintain
large tapes.

For optimization workloads, schedule selection considers the primal and
adjoint together.

## I/O

Routine diagnostics reduce in situ. Full 6D distributions are deliberate
checkpoints, not frequent plots.

Chunk layouts follow state/tile structure. Use large buffered/asynchronous
writes when supported. Do not checkpoint reconstructible caches/factors unless
restart time justifies them.

FortIO owns generic throughput improvements.

## Upstream ownership

If a controlled dependency is material in the profile, optimize it at the
owning layer:

- FortFEM: tensor/matrix-free spatial actions;
- FortNum: FFT/special functions/banded/verified numerics;
- FortSparse: sparse/factor/preconditioner kernels;
- FortSym: generated algebra/schedules;
- FortAD: derivative transformations;
- FortIO: I/O;
- KIN6D: state layout, physics composition, global scheduling.

Retain a focused reproducer and recheck the original KIN6D case after upstream
promotion.

## Performance CI

Use three levels:

1. **fast sentinels** — seconds; allocation/vectorization/gross regressions;
2. **dedicated audits** — pinned hardware; realistic operators/solves;
3. **science benchmarks** — physical device cases and certificate target.

Normal development must remain fast; continuous performance awareness does not
mean every commit waits for a large device run.

Chris&AI
