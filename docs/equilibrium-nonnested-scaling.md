# NN2 Fourier solve scaling prototype

[V] This opt-in prototype belongs to the resistive coupled model in
[equilibrium-nonnested-nn2.md](equilibrium-nonnested-nn2.md). The assembled
Newton path remains the default. Its states retain the `candidate` claim.
Neither Krylov convergence nor agreement with the assembled path certifies an
MHD equilibrium or a continuous discretization error.

## Operator and storage

[V] `nc_problem_t%matrix_free = .true.` selects FortNum flexible GMRES with a
right preconditioner. The Jacobian action curls the perturbation using
FortFEM, synthesizes B, grad p and grad Phi on the residual's phi grid using
FortNum FFTs, evaluates the existing generated directional leaves, and
projects the current and flux against the real Fourier test functions.
Projection uses phi means, so the cosine and sine test coefficients do not
carry the factor two used for Fourier amplitudes. Toroidal derivatives use
the physical n. Dirichlet and floating-conductor row/column reductions are
the same reductions as in the assembled solve.

The preconditioner averages the **pointwise leaf derivatives**, rather than
averaging B and differentiating the averaged field. It contracts those
coefficients with each n's test/trial operators, retaining the two real
sine/cosine blocks within n. The loop-voltage unknown belongs to n=0 only.
FortSparse factors each reduced n block independently. OpenMP distributes
spatial quadrature/phi work over elements and factors/solves over n blocks.
There is no private Krylov implementation and no fallback to assembling the
full Jacobian when GMRES fails.

[V] The mean coefficients use O(elements × quadrature points) memory. Per-mode
matrices, factors and bounded-restart Krylov vectors use O(N) memory at a
fixed mesh and restart. Element preconditioner arrays contain at most two
real blocks. These are storage-complexity statements, not a measured
wall-time law. The current residual evaluation and source loading retain
the reference trigonometric loops.

[V] `nfp` defaults to one. For nfp>1 the field layout still labels physical
modes, but only n=0,nfp,2nfp,... are active. Omitted harmonics of the seed
must vanish, the datum hooks must have the declared periodicity, and the
maximum physical n must be divisible by nfp. The nonlinear rule samples
one field period with default 4(N/nfp)+1 points. FFT synthesis/projection
compresses the active harmonics while retaining physical n in derivatives.
Inactive n blocks are not factored. This symmetry restriction is unsuitable
for a perturbation that breaks the declared field period.

[V] The matrix-free path uses the mean-block absolute action for the Newton
backward-error and line-search scales. It does not claim to reproduce the
assembled absolute Jacobian action. The comparison tests therefore check
**every accepted Newton iterate**, their count and the final state against
the unchanged assembled solver at a relative tolerance of 1e-10. The island
gate is applied separately to A, p, Phi and V, and the vanishing multiplier
is normalized by the Ampere scale A/mu0. Both solves run to 1e-13 backward
error. A failed gate remains a failure.

## Reproduction and phase attribution

Configure the native `dev-fast` preset with `RelWithDebInfo`, build
`test_nonnested_island`, `test_nonnested_coupled` and
`test_nonnested_leaves`, then run the affected registered CTests:

```sh
ctest --test-dir BUILD -j8 --output-on-failure -R '^nonnested_matrixfree_'
ctest --test-dir BUILD -j8 --output-on-failure \
  -R '^(nonnested_symbolic_derivatives_and_d3_identity|nonnested_coupled_fourier_reduction_p[23]|kin6d_layering)$'
OMP_NUM_THREADS=8 scripts/study_nonnested_scaling.sh BUILD RESULTS local
OMP_NUM_THREADS=8 scripts/study_nonnested_scaling.sh BUILD RESULTS threads
```

Set `KIN6D_PROFILE=1` to print wall time and Linux resident-memory samples
at phase boundaries. `nn2-dense-leaves` is the pointwise/Jacobian contraction
inside `element`, `nn2-triplet-scatter` is `scatter_jacobian`,
`nn2-compress-B` and `nn2-compress-J` are FortSparse triplet compression,
and `nn2-chain-curl` is the sparse product with the curl. `nn2-factor`
measures the assembled reduced LU. The matrix-free phases are
`nn2-mean-leaves`, `nn2-mode-assembly` and `nn2-mode-factor`.
`nn2-linear-solve` includes the direct/Krylov solve and any factorization.
`nn2-line-search` includes scaling and trial residuals. Inclusive phases
must not be added to their children.

`MEM` reports current RSS and process-wide high-water RSS in KiB, sampled
at the end of each phase. It does **not** claim an independently reset peak
for each phase. `/usr/bin/time -v` supplies the whole-process maximum.
The benchmark driver retains the exit code of every Krylov scan. That scan
varies chi_parallel/chi_perp at fixed g_parallel/g_perp=100, separating
transport anisotropy from the conductivity scan in the original island test.
The `profile` arguments are degree, nr, maximum physical n, chi ratio,
`dense|matrixfree`, optional Newton-step limit, optional conductivity ratio.

Larger cases use `study_nonnested_scaling_scluster.sbatch`. Its sixteen-core
reservation provides a nominal 188.8 GiB share under the cluster's CR_CORE
policy, with eight solver threads and the cluster skill's memory guard.
It completes a small probe before larger solves, runs only the focused
executables, and retains individual raw logs outside Git.

## Attribution to source

| Measured phase | Implementation |
|---|---|
| Constant operators | `src/mhd/nonnested_coupled.f90:528` and FortFEM `assembly_meridian_fourier_curl.f90` |
| Fused leaf evaluation and local contraction | `src/mhd/nonnested_coupled.f90:797`, generated JVP calls at 897/905, contractions at 938/939 |
| Jacobian triplets and compression | `src/mhd/nonnested_coupled.f90:391`, scatter at 982, FortSparse `fortsparse_csc.f90:368` |
| Curl chain | `src/mhd/nonnested_coupled.f90:445`, FortSparse `fortsparse_csc.f90:274` |
| Direct factorization | `src/mhd/nonnested_coupled.f90:1141` and `factor_reduced` at 1200 |
| Line search | `src/mhd/nonnested_coupled.f90:598` |
| FFT action | `src/mhd/nonnested_coupled_matrix_free.f90:156` |
| Mean leaf and per-n assembly/factorization | `src/mhd/nonnested_coupled_matrix_free.f90:260` and 503 |

The fused local phase includes the generated leaves and their dense local
contractions. Its timer does not isolate individual generated functions.


## Measured dense costs [V]

All tables below are single runs, not confidence intervals. The synthetic
island fixture is the existing periodic helical coil case, with
chi_parallel/chi_perp = g_parallel/g_perp = 100 unless stated otherwise.
Eight solver threads were used. Times are seconds summed over the solve,
including its final Jacobian rebuild. The leaf column is part of Jacobian
assembly and must not be added again. Peak is whole-process RSS in GiB.
Per-phase end RSS and cumulative process HWM are in the accompanying TSVs.

| Degree | nr | N | Leaf/contraction | Jacobian | Factor | Line search | Total | Peak GiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P2 | 8 | 1 | 0.59 | 1.16 | 1.30 | 0.020 | 3.58 | 0.61 |
| P2 | 8 | 2 | 2.78 | 4.51 | 2.33 | 0.046 | 9.30 | 1.57 |
| P2 | 8 | 4 | 15.71 | 21.76 | 7.71 | 0.138 | 35.62 | 4.84 |
| P2 | 16 | 1 | 2.33 | 5.15 | 6.48 | 0.074 | 16.10 | 2.45 |
| P2 | 16 | 2 | 11.21 | 19.44 | 11.73 | 0.179 | 41.45 | 6.62 |
| P2 | 16 | 4 | 62.86 | 90.62 | 41.11 | 0.588 | 159.09 | 21.24 |
| P3 | 8 | 1 | 3.04 | 5.02 | 3.52 | 0.053 | 11.17 | 1.83 |
| P3 | 8 | 2 | 14.28 | 19.98 | 9.88 | 0.128 | 35.99 | 4.75 |
| P3 | 8 | 4 | 88.87 | 109.28 | 35.16 | 0.393 | 161.70 | 15.26 |

The dense contraction loops form all mode pairs at each phi point, so their
work includes an additional nphi = 4N+1 factor beyond quadratic storage.
At chi = g = 1e6 the nine available cells also solve. Peak memory changes
little: P2/nr=16/N=4 reaches 21.40 GiB, with total time 254.88 s and four
Newton steps, versus two at ratio 100. The high-ratio phase table is retained
in `benchmarks/nn2_scaling/scluster-dense-highchi.tsv`.

All P3/nr=16/N=1,2,4 cells were attempted and fail before the coupled solve.
The existing NN1 vacuum seed rejects its n=1 backward error
1.1594441745099641e-9 against 1e-9 at
`src/mhd/equilibrium_nonnested.f90:173`. Projection and the n=0 solve pass.
This pre-existing limitation also affects the unchanged dense path. Its
setup RSS is not a coupled-solve measurement. A short reproduction is
`test_nonnested_island profile 3 16 1 1e2 dense`.

### Sampled CPU attribution [V]

A separate local P2/nr=8/N=2 dense run used `gprofng collect app -p 1`.
Its 86.898 sampled CPU seconds include all threads and profiler overhead;
they are not the wall times above. Exclusive samples attribute 41.775 s
(48.07%) to `element::add_block`, specifically the two `matmul` updates at
`src/mhd/nonnested_coupled.f90:938-939`. The generated current and flux JVPs
account for 0.385 s (0.44%) and 0.143 s (0.16%), respectively, at
`src/generated/kin6d_nonnested_leaves.f90:402` and `:324`. Thus the dense
local contractions dominate the generated derivative arithmetic in this
case. The coil fixture's segment evaluation accounts for 18.001 s (20.72%)
across setup and data hooks. Samples and phase timers measure different
quantities and are retained separately.

## Mode and thread scaling [V]

The original generic identity-tensor operator construction also allocated
all mode pairs. The prototype therefore uses the exact sparse identity
forms added in [FortFEM PR 106](https://github.com/lazy-fortran/fortfem/pull/106),
pinned at 7e3e2d281887081c7280b5c4fafc050be7c633b1. The old general API is
unchanged. Its P2/P3, N=0..2 quadrature and curl-grad checks pass 69/69 in
Debug. This dependency is necessary for linear construction storage.

At fixed P2/nr=8, ratio 100, the sCluster matrix-free scan gives:

| N | Peak GiB | Total seconds | GMRES iterations per Newton step |
|---:|---:|---:|---|
| 1 | 0.304 | 5.61 | 8, 9 |
| 2 | 0.572 | 6.38 | 9, 9 |
| 4 | 1.095 | 11.16 | 9, 9 |
| 8 | 1.990 | 15.08 | 9, 10 |
| 16 | 3.328 | 34.43 | 9, 10 |
| 20 | 3.997 | 42.44 | 9, 10 |

N=20 has 572033 unknowns. RSS/(2N+1) stays about 0.10--0.12 GiB,
consistent with the linear storage design, with fixed restart 40.
This is a synthetic-fixture result, not a W7-X or LHD equilibrium test.
At fixed physical N, increasing nfp removes inactive factor blocks but the
public state layout still allocates all physical modes. The nfp=5 gate
checks active physical n=5 against the dense operator.

Local strong scaling uses P2/nr=8/N=2, two Newton steps and 9,9 GMRES
iterations throughout, on an AMD Ryzen 9 5950X:

| Threads | Total seconds | Speedup | Peak GiB |
|---:|---:|---:|---:|
| 1 | 36.38 | 1.00 | 0.442 |
| 2 | 23.54 | 1.55 | 0.517 |
| 4 | 16.95 | 2.15 | 0.518 |
| 8 | 12.72 | 2.86 | 0.580 |

A separate registered check compares full final states for 1,2,4,8 threads.
The largest normalized block difference in its final verbose run is
1.5001e-13, below 1e-10. A single wall-time sample per thread count does not
establish a statistical scaling limit.

## Anisotropy and qualification

[V] With g_parallel/g_perp held at 100, increasing chi causes true linear
residual stagnation. At N=20 the successful ratio-100 solve takes 9,10
iterations and ratio 1e3 takes 13,13,13. The first solve at larger ratios is:

| chi ratio | GMRES iterations | True normalized residual | Result |
|---:|---:|---:|---|
| 1e4 | 25 | 1.664e-13 | stagnation, info=3 |
| 1e5 | 36 | 9.204e-13 | stagnation, info=3 |
| 1e6 | 41 | 3.529e-12 | stagnation, info=3 |

[V] The pattern persists at smaller N: for chi=1e6 and N=1,2,4 the first
solve takes 37,39,41 iterations and stagnates at approximately 3.2--3.5e-12.
The original joint chi=g=1e6 island gate is worse: P2/N=1,2 take 63,69
iterations and stop at 3.61e-11,4.27e-11; P3/N=1,2 take 79,117 iterations
and stop at 2.40e-10,2.54e-10. These are explicit failures at the unchanged
linear target 1e-13. No dense fallback or relaxed acceptance is used.
They occur at the vacuum seed p=0, where the pressure-to-field feedback
term (derivative of transport tensor times grad p) vanishes.

[H] Averaging the anisotropic tensor loses the cancellation between Fourier
modes that represents a field-line-constant pressure. In particular,
v.Kbar.v = chi_perp |v|^2 + (chi_par-chi_perp) mean((b.v)^2) need not be
small for a direction participating in that cancellation. This explains a
plausible loss of the pressure near-kernel. The fixed-g experiment isolates
a pressure-transport limitation; the stronger joint-ratio failures also
suggest a conductivity/Phi near-kernel problem. Neither experiment proves
an irreducible floating-point floor.

[H] Retain phi coupling in a pressure/Phi Schur block, using an approximate
field-line solve and a coarse space for line-constant components, while
keeping the per-n A factors. At nonzero pressure include the transport/field
feedback in that Schur approximation. FortNum's existing inner solver can
be used with the current outer flexible GMRES. Independently test variable
and row equilibration and cancellation in leaf accumulation. A better
preconditioner alone has not been shown to remove the measured stagnation.

[V] The reviewed focused native suite passes 13/14 checks: all four
reduction comparisons, three of four island comparisons, field-period and
thread-invariance checks, the existing P2/P3 reduction tests, generated-leaf
derivatives and layering. P3/N=1 fails the strict island trace gate:
Phi relative difference 1.1461e-10 exceeds 1e-10. Its A, normalized multiplier,
p and V differences are 4.0524e-13, 2.4119e-13, 9.9066e-14 and 1.0503e-12.
The high-anisotropy gates also fail as above. A later focused smoke set and
thread-invariance rerun pass after formatting and scope-only edits. The
prototype must remain a draft and must not be queued. There is no green Fo
Gremlin generation; native checks are the evidence reported here.

[V] Reproduction scripts, compact phase summaries and source/log provenance
are in [benchmarks/nn2_scaling](../benchmarks/nn2_scaling/README.md).
Raw experiments and profiler data remain outside Git. The next acceptance
step is to improve solve accuracy and pass all
unchanged dense-equivalence gates before considering promotion.

## I34 derivative and conditioning checks

[V] `nc_vjp` applies the transpose of `nc_jvp` in the Euclidean pairing of
coefficient arrays, using the existing generated current/flux VJP leaves.
It reverses the floating-row reduction before the element pass and retains
boundary input columns. Its FFT testing, spatial contractions and transposed
curl use O(M) storage at fixed mesh, with OpenMP over elements. The reduction
fixture checks both the JVP/VJP pairing and the complete assembled transpose,
and the island fixture independently checks the floating-conductor pairing.
The finite-difference directional gate is now 1e-7, matching I3.

[V] Both solver paths now use FortFEM's exact constant identity forms.
The assembled nonlinear element Jacobian remains the independent comparator.
This removes quadrature-only off-mode entries in the constant curl/gauge
operators, rather than compensating for them by relaxing the iterate gate.
FortFEM is pinned to merged #107 (21f236191c39688bebe25f98b185c32b8209ec70),
which contains the identity-form work from superseded draft #106.

[V] A small failing cell was examined independently: P2, nr=4, N=1,
chi_parallel/chi_perp = g_parallel/g_perp = 1e6, vacuum seed. The pressure
Dirichlet block has 1,344 free coefficients. Generalized symmetric eigenvalues
of the full pressure block against the phi-averaged block range from
0.4054201588 to 1.5755062. The relative symmetry defects are below 7.5e-17.
Thus the pressure near-kernel-loss hypothesis is unsupported at this cell.
The existing averaged preconditioner already retains p/Phi coupling within
each n block. This pressure-only spectrum does not exclude coupled
nonnormality or failures on finer meshes and larger mode sets.

| P2/nr4/N1, chi=g=1e6 | Krylov iterations | normalized residual | result |
|---|---:|---:|---|
| V3 path | 31 | 4.984e-12 | fails 1e-13 |
| Trial update-based equilibration | 28 | 5.306e-13 | fails, discarded |
| Trial identity plus mean-defect action | 45 | 2.375e-10 | fails, discarded |

[V] On the failed V3 increment, the assembled and FFT actions differ by
2.1401e-11 in the row-scaled, normalized right-hand-side metric. Their
residuals are 2.1730e-11 and 4.9837e-12 respectively. This is evidence of
arithmetic sensitivity at the requested linear tolerance. It does not prove
an irreducible precision floor or identify the complete high-chi root cause.
Neither unsuccessful trial is retained, and the linear tolerance, nonlinear
and iterate gates remain unchanged. No new generic numerical primitive was
justified by these experiments.

[H] A successful repair still needs a stable coupled action or solve and a
measured high-chi comparison against the assembled states. The matrix-free
path remains opt-in. The older memory, iteration and strong-scaling tables
above are V3 baseline evidence, not measurements of the I34 revision.

[V] The four-thread sCluster pilot passes 10/10 registered selections:
all four reduction comparisons, all four island comparisons, nfp=5 and the
symbolic-leaf check. P3/nr8/N1's Phi iterate difference is 8.5292e-11 against
the unchanged 1e-10 gate. This is close to the threshold, and the independent
eight-thread repeat remains required. On the final local code, the P2/N1
VJP/assembled-transpose error is 6.0773e-16 and the floating-conductor pairing
error is 2.7322e-14. The selected reduction, floating island, leaf and layering
checks pass. These results do not close the high-chi gates.
