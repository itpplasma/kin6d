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
`fortsparse_native_ldlt.f90:270` (`numeric`). At that measurement, upstream main was checked and
matched the 09fc2b4 pin. This is a factorization hotspot observation,
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

## Native CPU follow-up, 2026-10-22

[V] KIN6D 072d226 uses native supernodal LU with AMD for the map Newton/KKT
matrix (`src/mhd/ideal_map.f90`, `factor_masked`). The FortSparse pin now includes
[the default scalar LDLT row-merge fix](https://github.com/lazy-fortran/fortsparse/pull/28),
faa7b5348992dbb67eb2c8a2b60c33b17d211759. The historical stalled LDLT sample
above and the current production LU measurements are different solver paths.

[V] Default LDLT previously repeated binary searches and a 120-character status
write for each Schur entry. Upstream #26 had fixed blocked recovery, while default
scalar execution still missed the merge workspace. The new default schedules
independent row merges through OpenMP, preserves scalar pivot choices and
contraction order, accounts merge workspace under strict caps, and retains
entry-wise updates under LOW_MEMORY. Its no-fill regression checks a known
solution, the original-matrix residual and a relative timing bound against LOW.

[V] sCluster EPYC 9455, GCC 16.2 Release, OpenBLAS 0.3.34 at one thread, AMD 2.4.6,
OpenMP passive wait, close binding within an eight-core Slurm allocation.
Five fresh processes per thread count; medians below. Runs used identical source
snapshots on node4/node18. `Map` times map construction plus Newton at the finest
rung. `Ladder` is the complete three-rung process, including meshing, independent
2D solves, defects and D3 checks. CPU is full-process user plus system time;
RSS is the largest observed individual-process peak, in KiB. Phase medians need
not add exactly to the median map wall time.

| Free DOFs | Threads | Map s | Energy/gradient s | Hessian s | Factor/solve s | Ladder s | CPU s | Peak RSS KiB |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 18,567 | 1 | 5.46 | 0.82 | 2.21 | 2.39 | 8.45 | 8.40 | 95,932 |
| 18,567 | 2 | 3.91 | 0.43 | 1.30 | 2.11 | 6.40 | 8.80 | 100,580 |
| 18,567 | 4 | 3.16 | 0.24 | 0.84 | 2.05 | 5.34 | 9.02 | 105,596 |
| 18,567 | 8 | 2.80 | 0.14 | 0.62 | 2.01 | 4.83 | 9.56 | 120,392 |
| 41,615 | 1 | 5.67 | 0.68 | 2.66 | 2.25 | 14.93 | 14.87 | 173,048 |
| 41,615 | 2 | 3.99 | 0.35 | 1.64 | 1.94 | 11.37 | 14.95 | 180,616 |
| 41,615 | 4 | 3.18 | 0.20 | 1.16 | 1.77 | 10.00 | 15.65 | 183,088 |
| 41,615 | 8 | 2.81 | 0.12 | 0.88 | 1.76 | 9.26 | 16.92 | 192,952 |

[V] Matched old-pin controls (4f716f6, five processes per count, same compiler,
libraries, inputs and binding policy; node4) are below. The measured native LU
source is unchanged between pins. New-pin ladder wall medians are 2.7–10.5%
lower in these observations, with overlapping timing ranges and different Slurm
allocations, so no LU speedup is attributed to the pin change. Peak RSS changes
range from -1.5% to +4.2%. All printed endpoints agree exactly between pins.
Control jobs 3355225/3355226 passed after pilot 3355205.

| Free DOFs | Threads | Old map s | Old factor/solve s | Old ladder s | Old CPU s | Old peak RSS KiB |
|---:|---:|---:|---:|---:|---:|---:|
| 18,567 | 1 | 5.91 | 2.60 | 9.10 | 9.05 | 97,356 |
| 18,567 | 2 | 4.04 | 2.21 | 6.58 | 9.00 | 101,536 |
| 18,567 | 4 | 3.35 | 2.18 | 5.65 | 9.25 | 102,720 |
| 18,567 | 8 | 3.09 | 2.21 | 5.31 | 10.24 | 115,548 |
| 41,615 | 1 | 6.05 | 2.50 | 15.90 | 15.85 | 172,944 |
| 41,615 | 2 | 4.48 | 2.27 | 12.71 | 16.71 | 180,636 |
| 41,615 | 4 | 3.50 | 1.99 | 11.00 | 17.38 | 184,068 |
| 41,615 | 8 | 3.01 | 1.92 | 10.17 | 18.45 | 192,560 |

[V] Map scaling is 1.95x/2.02x, while factor/solve scaling is 1.19x/1.28x.
The finest cases take 21/10 Newton iterations respectively, so their raw times
are not a fixed-iteration size-scaling comparison. Every run passes the exact
flux and independent 2D convergence-order gates. Endpoints are identical across
all five trials and thread counts within each ladder: normalized exact-flux
defect 2.1522071047864953e-6 / 5.931685014394489e-7, energy
9361.367520435855 / 9361.36784912228. They remain candidate solutions with
continuous enclosure unavailable.

[V] Generic scalar LDLT strong scaling is measured separately on a 768-by-768
signed-diagonal/rank-one fixture. Median factor times at 1/2/4/8 threads are
0.983715/0.577523/0.330505/0.215278 s (4.57x, 57% efficiency at eight threads).
Actual maximum teams are 1/2/4/8; maximum absolute original-matrix residuals are 2.89e-14 and
known-solution errors stay below 2e-13. Full-process wall medians are
1.02/0.63/0.38/0.26 s, CPU 0.99/1.12/1.14/1.15 s, peak RSS maxima
37,392/39,312/39,828/40,552 KiB. This longer diagnostic is outside FortSparse's
sub-second admission corpus. It does not establish scaling of the production LU.

[V] A separate five-process replay of that initial Newton matrix (18,824
stored rows including fixed identity rows, 18,567 free DOFs) factors in
4.729101/4.307087/4.109780/3.989423 s at 1/2/4/8 threads, only 1.19x scaling.
The normalized original-matrix residual is 1.26364e-15 and maximum known-solution
error 2.11457e-10 at every count. Full-process wall medians are
4.79/4.37/4.18/4.15 s, CPU 4.76/4.93/5.03/4.94 s and peak RSS maxima
57,632/58,156/63,972/61,112 KiB. Replay job 3355283 passed. This sparse consumer
matrix establishes a narrower scaling result than the dense clique.

[V] Reproduce the normal ladder with `test_ideal_map_axisymmetric solovev 2`.
The larger ladder changes only `base = 32` to `base = 48` in that test driver
(boundary-node ladders 32/64/128 and 48/96/192). Use `OMP_NUM_THREADS=1,2,4,8`
in separate processes and `OPENBLAS_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`,
`OMP_PROC_BIND=close`, `OMP_PLACES=cores`, `OMP_WAIT_POLICY=PASSIVE`.
Raw logs, full timing ranges, binary hashes, the larger-driver source and build
recipes are retained outside Git at `/home/ert/runs/p1ldlt-20261022` on sCluster
and `/mnt/storage/scratch/kin6d-islands/P1` locally. Measurement jobs:
3355146/3355168 (normal), 3355147/3355169 (larger), 3355148 (generic LDLT).

[V] A separate blocked-LDLT trial was rejected for the map. Its first frozen
18,824-row Newton system factors quickly (1.03 s on a shared workstation at four threads,
normalized residual 7.48e-16), but that sample does not represent later Newton
matrices. The complete 18,567-DOF one-thread ladder timed out after 300 s on
sCluster at 1,755,232 KiB peak RSS. The 41,615-DOF one-thread finest solve took
40.51 s, including 37.00 s in factor/solve, with complete-ladder wall/CPU
59.95/59.90 s and peak RSS 224,932 KiB. Those single observations support keeping
production supernodal LU; no consumer factor-family change is included here.
The outstanding LDLT limitation concerns full Newton matrix sequences, beyond
the entry-update regression and clique scaling fixture fixed upstream.

[V] These are the N=0 map restriction. They establish no matched-accuracy timing
comparison with VMEC++ CPU/GPU, no whole-equilibrium GPU speedup and no general
non-axisymmetric solver performance claim. Saturated-workstation production
samples at 1/8 threads both timed out at 60 s; those censored observations are
retained separately from the completed cluster measurements.

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
The top-level build offers `KIN6D_BUILD_CUDA_PRECISION_BENCHMARK=ON` with an
explicit `CMAKE_CUDA_ARCHITECTURES`; the option defaults to OFF. The standalone
build above also requires an explicit architecture. Production FP64 defaults
remain unchanged. Submit `run_slurm.sh` from the dedicated snapshot root;
Slurm copies the script to its spool directory before execution.

[V] Full logs, profiles, compiler/resource reports and binary/library hashes:
`/home/ert/data/handoffs/20261021_cluster_triage/kin6d_G/`.
The original Fo driver failed on an MPI module in its default build tree, then
on unexpanded `__LINE__` at `test/test_curve_chart.f90:60` in native compilation.
The latter has an existing owning Fo fix, 18df703f78dd276429abac4ba3dd4f341e1dfd53.
Focused checks used its documented `FO_CMAKE_DELEGATE_BUILD=1` CMake route;
no consumer compiler workaround was added. The task-owned Gremlin was stopped;
no detached watcher remains.

Chris&AI

## Generated batch entry points

[V] The handwritten `evaluate` wrapper has been replaced by FortSym's
`emit_cuda_batch_ir`. The same Hessian/HVP roots now produce both scalar leaves
and the coalesced SoA global kernels. Launch geometry, allocations and transfers
remain in the consumer; CUDA remains opt-in with an explicit architecture list.

[V] RTX 5060 Ti, sm_120, 65536 points, 11 repetitions, 8-thread CPU reference:
Hessian/HVP FP32 resident medians are 14.82/8.03 microseconds, versus
271.42/290.66 microseconds in GPU FP64. FP32 maximum scaled defects are
8.14e-7/6.21e-7. Transfer-inclusive Hessian/HVP costs are approximately
1.65/0.83 ms, compared with CPU medians 2.29/0.85 ms. Cold setup is separate;
these shared-host kernel numbers establish no full-solver speedup.
