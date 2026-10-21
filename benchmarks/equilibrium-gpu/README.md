# Equilibrium CUDA precision prototype

[V] This experiment measures local Hessians and Hessian-vector products (HVPs)
from the frozen-in MHD energy defined in `gen/gen_ideal_map.f90`. It does not
replace production assembly, Newton solves, residuals or estimators. Its claim
is **candidate**, with **enclosure unavailable**. Every deformation used for
comparison is a dimensionless, well-conditioned sample, not a solved equilibrium.

[V] The generator emits CUDA FP64, FP32 and mixed leaves from the production
Hessian roots; HVP is the symbolic JVP of the same gradient. Mixed means FP32
inputs/arithmetic and FP64 outputs. It adds no FP64 accumulation or refinement.
The scalar Euclidean Hessian is symmetric, so its transpose action is the same
HVP. Chart chains, constraints and coefficient-space pairings remain outside it.

[V] [FortSym](https://github.com/lazy-fortran/fortsym) revision
7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6 delegates scalar emission to
[FortGen](https://github.com/lazy-fortran/fortgen) revision
42a4c6cce9272a6bc92e0539b21f3ffc0c97a068. Both already support FP64/FP32/mixed
CUDA leaves; OpenMP target and OpenACC are alternative Fortran decorations.
They provide no launch schedule, transfers, sparse solver or automatic
mixed-precision iterative refinement. Offload runtime performance was not tested.
No new generic emitter feature was needed.

[V] The previous KIN6D pin emitted overlong Fortran statements and general CUDA
`pow` calls for negative integer powers in this consumer. Existing upstream
fixes supply continued statements and typed reciprocal products. An isolated
build of those fixes was required: the pre-existing FortSym build's linked
FortGen was 488335bf491bfdf99c257cf6da822e5e95f5d2d6 despite newer source HEAD.
The consumer was regenerated with the rebuilt libraries; the production energy,
gradient and Hessian changed only in provenance comments.

## CPU profile and limits

[V] CPU baseline a25fa0173459e3c29d564fac7fae19de99cfaa42, Ryzen 5950X,
GNU Fortran 16, Release, OpenBLAS one thread. The actual map acceptance ladder
was run at 1 and 8 OpenMP threads. For its finest P2/N=0 case (18,567 free DOFs):

| Status | Threads | Total s | Energy/gradient s | Hessian s | Factor/solve s |
|---|---:|---:|---:|---:|---:|
| [V] | 1 | 12.15 | 2.28 | 4.87 | 4.90 |
| [V] | 8 | 4.47 | 0.34 | 0.97 | 3.09 |

[V] Observed total scaling is 2.72x and Hessian scaling 5.02x. Reported endpoint
energy, displacement norm and oracle errors match at the printed precision;
the finest exact-solution defect is 2.1522105769793812e-6. These timings are
shared-host observations, not isolated speed rankings or certified error bounds.

[V] `perf` of the real map binary attributes 12.87% of cycles to the element
Hessian loop, 10.17% to `local_derivatives`, 7.03% to CSR block assembly,
5.74% to `point_jet`, and only 2.35% to the scalar Hessian leaf. See
`src/mhd/ideal_map.f90` (`local_derivatives`, `ideal_map_hessian`) and
FortFEM's `csr_assemble_blocks`. Eliminating just that leaf could improve this
sampled CPU execution by at most 1/(1-0.0235) = 1.024x with everything else fixed.
The prototype therefore establishes a building block, not dominant-phase or
whole-solver acceleration. GPU element contractions and assembly remain unported.

[V] Main ead3cca switched sparse solves to native FortSparse CPU factors and was
merged into this lane. A bounded resample of that real consumer was stopped
after about 150 s without a completed ladder timing. A five-second sample
attributes 39.55% to `fortsparse_native_ldlt.f90:629` (`row_put`), 29.32% to
`fortsparse_status.f90:42` (`status_set`), and 23.10% to
`fortsparse_native_ldlt.f90:270` (`numeric`). Upstream main was checked and
remains the pinned 09fc2b4 revision. This is a factorization hotspot observation,
not a completed solve or an established cause of every slowdown. The earlier
AUTO-backend totals above must not be attributed to the new native backend.
The timed ladder is the N=0 restriction of the map lane; a completed
non-axisymmetric solve profile remains unavailable.

[V] The current element schedule forms dense local `nl` by `nl` matrices with
`nl = 3*nmode*nloc` and retains per-cell reduced blocks. Its contraction work
is quadratic in `nmode*nloc` per quadrature point; the scalar HVP alone does
not remove that work. [H] A complete matrix-free element action can replace
those quadratic coefficient contractions by gather/action/scatter work linear
in local coefficient count, subject to chart chains and gauge constraints.

[V] K6-mini GS (3,667 DOFs) took 61.26 s at one thread and 12.16 s at eight
(5.04x, the eight-thread run includes sampling overhead). The source-cell timer
was 61.00/11.70 s. Sampling attributes 31.40% to `source_products`, 20.12% to
`triangle_basis_jet`, 7.69% to profile evaluation and 7.34% to cut products.
The relevant path is `src/equilibrium/gs_nonlinear.f90:799` onward. The cell timer
includes current normalization (`app/gs_phase1.f90:248` and
`gs_normalize_current`), whereas `GS_T_ASSEMBLY` counts Newton assembly calls;
those overlapping totals must not be added. The estimator was 0.203/0.0644 s
and its majorant agreed to 9.7e-14 relative (status 4, residual verified,
approximate frozen-source estimate, continuous certification unavailable).

## Device measurements

[V] Local device 1 was an RTX 5060 Ti with about 1.2 GiB free beside an existing
service; no service was stopped. CUDA 13.4, sm_120, 128 threads/block, CUDA FMA
enabled, strict FP64 Fortran reference, OpenMP passive wait, 21 repetitions.
Each sample is independent; all allocations occur outside timed repeated calls.
The ABI is Fortran ISO_C_BINDING to a C CUDA wrapper; layout is
`(sample, component)` in Fortran, with explicit dimensions, precision and status.

[V] Medians for 32,768 samples and eight CPU threads follow. The full ranges,
transfer costs and large 262,144-sample case are in `measurements.tsv`.
CPU variability and shared GPU use limit generalization of the observed ratios.

| Status | Action | Precision | CPU ms | GPU resident ms | Scaled error |
|---|---|---|---:|---:|---:|
| [V] | Hessian | FP64 | 2.6995 | 0.150816 | 2.17e-15 |
| [V] | Hessian | FP32 | 2.6995 | 0.014048 | 7.92e-7 |
| [V] | Hessian | mixed | 2.6995 | 0.018624 | 7.92e-7 |
| [V] | HVP | FP64 | 1.9318 | 0.165056 | 1.39e-15 |
| [V] | HVP | FP32 | 1.9318 | 0.008416 | 6.21e-7 |
| [V] | HVP | mixed | 1.9318 | 0.010624 | 6.21e-7 |

[V] Observed FP32 resident-kernel speedups versus the eight-thread reference
are 192x (Hessian) and 230x (HVP); versus GPU FP64 they are 10.7x and 19.6x.
For the Hessian, H2D + kernel + D2H is 1.465 ms in FP32 versus 0.014 ms resident.
The first FP64 call including context/setup/conversion/transfers took 188 ms;
subsequent precision variants start with a live CUDA context. Timings exclude
production geometry, assembly, solves and error estimation. Copy timing is one
observation per variant, not a distribution. HVP emits nine scalars instead of
45, reducing output traffic/storage by five; it is not an element HVP yet.

[V] At 262,144 samples FP32 Hessian/HVP medians are 0.178464/0.037216 ms,
with maximum scaled errors 8.88e-7/6.22e-7. Explicit device storage peaks at
144,703,488 bytes across the tested variants. The compiler reports no spills;
Hessian FP64/FP32/mixed use 252/168/168 registers, HVP 96/56/56. These are
compiler resource reports, not measured occupancy or memory-throughput counters.

[V] Error is max over samples of max-component difference divided by
max(1, max-component reference magnitude). The unit floor applies only to this
dimensionless fixture. Every output is checked for finiteness. Independent
scalar-energy directional curvature agrees to 5.27e-8. Generated FP64 HVP is
also checked against the production packed Hessian times the direction.
Tolerance is 2e-12 for FP64 and 3e-5 for FP32/mixed, below the requested few-percent
candidate target on these samples; neither tolerance is an enclosure theorem.

[V] Before native-backend integration, derivative/rotation, normalized map
residual, D3 upper-bound, layering and K6-mini gates passed. After integration,
four gates pass, including K6-mini with its estimator scope/status and reference
values. The normalized-map residual gate was explicitly stopped at the bounded
task limit after 286.66 s; Fo records a signal failure, not acceptance. Its
expensive solve is at `test/test_ideal_map_axisymmetric.f90:427`, on the native
LDLT path identified above. That current-backend gate remains unverified.
No GPU-produced equilibrium has been passed to the
FP64 D3/Repin pipeline: preservation of estimates for an accelerated nonlinear
solve, global mixed-precision convergence and ill-conditioned geometry are
**unverified**. The optional prototype changes no production result path.

## Next measured steps

- [H] Resolve the measured native LDLT insertion/status-update cost in
  FortSparse before extrapolating full-solver GPU gains. Retain a focused
  factor/solve reproducer, unchanged residuals and an independent dense oracle;
  do not duplicate the sparse algorithm in KIN6D. No upstream factorization
  repair or indefinite GPU factor was implemented in this bounded prototype.
- [H] First port the complete map element action: geometry jets, chart chain,
  generated local HVP, test-function contractions and element gather/reduce.
  Keep matrix-free actions resident; preallocate per-thread/per-block scratch.
  Existing `ideal_map_hessian` allocates scratch per element at line 1019;
  move that schedule allocation out of the loop before treating it as a target
  CPU performance reference. Generic basis/Piola kernels belong in FortFEM;
  physics leaves stay in KIN6D `gen/`.
- [H] For GS, batch accepted cut-quadrature points and generated basis/source
  products across cells. Retain FP64 cut topology, acceptance quadrature and
  final residuals until a separate error analysis admits a lower precision.
  Accelerating only `gs_source` will miss the measured callbacks and basis work.
- [V] The pinned FortSparse already has optional native CUDA SPD panel,
  Schur-update and solve kernels (`cmake/FortSparseNativeCUDA.cmake`), using
  FP64 and an explicit device option. The profiled baseline selected CPU AUTO
  (SuperLU/MUMPS). Integrated main ead3cca selects native CPU factors; map
  gauge-constrained Newton systems can be indefinite.
  Existing SPD CUDA capability is not evidence for the required indefinite or
  FP32 mixed solver. It was inspected, not benchmarked here.
- [H] Reuse and test that upstream SPD path for eligible GS blocks; implement
  FP32 preconditioner/inner iterations and FP64 residual/refinement in the owning
  FortSparse/FortNum libraries. Use a gauge-aware flexible Krylov method for
  varying inexact preconditioners, with bounded vector memory and a robust
  FP64 fallback. Verify equation-level tangent/adjoint solves and D5 identities.
- [H] Compare total setup + Newton + residual + estimator time at matched
  1–3% observable accuracy and unchanged claim strength. Keep strict FP64
  residuals, stability/error estimates and certification arithmetic independent;
  retain FP64 GPU kernels for GPUs where that precision is competitive.
  A full GPU solve needs residency through assembly/operator, Krylov vectors,
  preconditioner and reductions. No full-GPU speedup is established here.

## Reproduction

[V] Regenerate from a KIN6D CMake tree configured with the pinned FortSym build:

```sh
mkdir -p benchmarks/equilibrium-gpu/generated
cmake --build build-dev-fast --target gen_ideal_map
build-dev-fast/gen_ideal_map --gpu
cmake -S benchmarks/equilibrium-gpu -B /tmp/kin6d-gpu -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=120
cmake --build /tmp/kin6d-gpu -j 8
CUDA_VISIBLE_DEVICES=1 ctest --test-dir /tmp/kin6d-gpu -j 4 \
  -R '^gpu_equilibrium_leaf_oracles$' --output-on-failure
OMP_NUM_THREADS=8 OMP_WAIT_POLICY=passive CUDA_VISIBLE_DEVICES=1 \
  /tmp/kin6d-gpu/benchmark_equilibrium_gpu 32768 21 > gpu.log
benchmarks/equilibrium-gpu/analyze.sh gpu.log
```

[V] `run_slurm.sh` rebuilds a source snapshot on the allocated GPU, runs the
same CTest and two sizes at 1/2/4/8 CPU threads. The T4 job uses sm_75 and
GCC 11 with CUDA 11.8. Prototype build needs only CMake, Fortran/OpenMP, CUDA
and the C++ runtime; no production backend configuration was designed.
The controller still owns build/runtime configuration choices.

[V] Full logs, profiles, compiler/resource reports and binary/library hashes:
`/home/ert/data/handoffs/20261021_cluster_triage/kin6d_G/`.
The original Fo driver failed on an MPI module in its default build tree, then
on unexpanded `__LINE__` at `test/test_curve_chart.f90:60` in native compilation.
The latter has an existing owning Fo fix, 18df703f78dd276429abac4ba3dd4f341e1dfd53.
Focused checks used its documented `FO_CMAKE_DELEGATE_BUILD=1` CMake route;
no consumer compiler workaround was added. The task-owned Gremlin was stopped;
no detached watcher remains.

Chris&AI
