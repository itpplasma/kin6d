# Nested-equilibrium comparison

[V] The authoritative exact fields, surfaces, flux and pressure are retained in
`gen/gen_nested_oracles.f90`. `nested_oracle_export` emits matched integer-family
inputs for VMEC++, DESC and GVEC. Sheared-family exports remain geometry only.
Use `examples/nested-oracles/{landreman,issan}-integer.nml`; the cases have
eps=0.5, delta=1/64 and tau=0 or 0.4. All reported state bounds remain unavailable.

[V] `input.landreman-paul-qa` is the fixed-boundary, zero-pressure, zero-current
Landreman–Paul QA input from
[stellarator_multiobjective](https://github.com/mishapadidar/stellarator_multiobjective/blob/be13ff521f59f3f6430d2900f1a5ee399ad36946/vmec_input_files/input.LandremanPaul2021_QA).
It is a cross-code case without an exact equilibrium oracle. Runtime settings
in the adapter replace the original very tight iteration tolerances.

[V] Build the opt-in `benchmark_nested_reference` target. The optional Python
adapter calls external codes and their readers; it does not implement KIN6D's
physics. In a Slurm allocation with the selected package installed:

```sh
python benchmarks/nested/run_reference.py vmecpp INPUT 8 OUT
build/benchmark_nested_reference OUT/samples.txt 0.5 0.4 0.015625
```

[V] Use tau=0 for Landreman, 0.4 for Issan. The adapter also accepts `desc` and
`gvec`. Add `--clockwise` for the integer-family inputs, including when using
VMEC: it records iota in the input convention after VMEC's internal theta flip.
For GVEC this reverses theta as well as zeta, preserving positive Jacobians.
Omit it for the supplied QA input. Never apply the integer-family comparator
to the QA case. VMEC uses
radial splines from SIMSOPT; no DESC fit or secondary solve is substituted.
Samples exclude the axis and edge and use three normalized toroidal-flux
surfaces. Vector components are Cartesian. The verifier compares B and pressure
at the candidate's physical positions, the exact flux label and the strong
force defect. These are equal-weight sample diagnostics, not volume norms,
surface distances, a D3 dual norm or continuous error bounds.

[V] The DESC adapter translates VMEC's all-zero axis sentinel into DESC's
native automatic-axis option. Treating it literally puts the initial axis at
R=0 and gave an unusable QA state despite an optimizer `xtol` success. Optimizer
termination is recorded separately from physical accuracy. Higher-resolution
GVEC cold starts can fail its positive-Jacobian gate; retain that failure rather
than reporting an equilibrium error or a performance success.

[V] Each point records the resolved input, state, solve and total time, process
peak RSS, thread count, external version and code-specific convergence status.
Compare physical errors before interpreting timings. `FTOL=1e-8` in VMEC is a
squared force stopping criterion; its square root is 1e-4 in that code's norm,
not an error bound on B. DESC/GVEC use 1e-6 stopping tolerances. No regression
gate between code versions is inferred from these settings.

[V] `boundary-inputs.tsv` records the admitted Fourier-input ladder, measured
on sCluster with gfortran 12.2, eight threads and the compiler-floor recipe.
Timing includes strict point-evaluation intervals and file output. It excludes
the memory guard's polling delay. Each point has one measurement; these are
input-preparation costs and cannot rank equilibrium solvers. The source was
acdcb99923c84bdac1a435ccb7f8758f7bb015b1. All listed points passed the unchanged
1% boundary-radius gate. Inside Slurm, reproduce the ladder with
`bash benchmarks/nested/run_inputs.sh BUILD OUT`. Plot the retained ledger
with `python benchmarks/nested/plot_inputs.py OUT`.

## Accuracy and cost study, 2026-10-22

[V] In both integer families, a coarse P3 point passes the sampled-fields gate with
2,295 coefficients and 1,098 free unknowns. Its coefficient count is within 4% of the
least measured VMEC point, although its native Newton phases cost much more. P3 reduces
the measured count by 38% relative to P2, while finer P3 meshes can worsen pressure
accuracy.

[V] The study uses KIN6D main 072d226 and its pinned FortSparse 4f716f6, plus a second
build with the P1 FortSparse patch d835267 (the same source tree as published faa7b53 in
[FortSparse #28](https://github.com/lazy-fortran/fortsparse/pull/28)). Main already
selects native supernodal pivoted LU with AMD in `src/mhd/ideal_map.f90`, so these are not
timings of the older scalar-LDLT map path. The native boundary/profile driver and
physical-point diagnostics are taken from pending [PR
#35](https://github.com/itpplasma/kin6d/pull/35), pinned at 558f1ce, and linked against
the listed main solver. This benchmark does not integrate that PR or qualify a 3D
common-core adapter on main.

[V] All reported timings run sequentially on one exclusive sCluster node, AMD EPYC 9455,
with affinity restricted to either one or eight physical cores. VMEC++ caps its native
team at `min(requested, (ns-1)/2)` for fixed boundaries. The low-radial-count points
therefore use four or six native workers within the eight-core budget, and the table marks
that cap. Other points request one/eight OpenMP threads. OpenBLAS and MKL use one
thread. KIN6D uses GNU Fortran 12.2, Release, OpenMP, AMD ordering and the
compiler-floor recipe. External packages are VMEC++ 0.9.1, GVEC 1.5.0 and DESC 0.17.3
with JAX 0.9.2, all on CPU. The existing DESC environment was reused. The VMEC++
checkout 0bd147f also built and solved the Landreman input out of tree, and the
fixed-boundary solver sources are unchanged from its v0.9.1 tag used for cluster timing.
The table uses one completed process per point, with superseded instrumentation attempts
retained separately, so small timing differences do not establish a regression or
speedup. No complete GPU equilibrium path is provided by this KIN6D prototype or the
inspected VMEC++ build. CUDA leaf timings in `benchmarks/equilibrium-gpu` cannot supply
a solver speedup.

[V] C counts the code-specific spatial coefficients: R/Z/lambda in VMEC++ and DESC,
X1/X2/lambda in GVEC, and three map-coordinate fields in KIN6D. VMEC++ follows the
requested `ns * mnmax * 3` count, doubled for a nonsymmetric run. DESC uses its actual
R, Z and lambda basis sizes, and GVEC its native radial degree-five B-spline basis size
times the actual Fourier mode count for each field. These counts precede axis, boundary
and gauge constraints. KIN6D also reports the free reduced unknown count after its
constraints. They are distinct quantities, and an exact cross-code comparison of
independent free unknowns is not available. Only P2 and P3 are admitted by
`ideal_seed_disk_mesh`. KIN6D has meridional triangular Pk elements and `2*N+1` toroidal
modes, and its poloidal resolution is the FE mesh, rather than a second Fourier basis.

[V] The 1% sampled-fields gate requires all of B RMS, pressure RMS and maximum
normalized flux-label error to be at most 0.01, with a successful native solve at both
CPU budgets. Integer cases use the exact Cartesian oracle, which is stronger than a
best-resolved numerical reference. QA uses the unmodified original input schedule as its
reference: final ns=201, MPOL=16, NTOR=12 and FTOL=2e-17. Its eight-thread native solve
took 40.90 s (69.19 s workflow, 661.0 MiB peak RSS). Comparisons are at identical
physical points, with native SIMSOPT splines and checked inverse coordinates. QA
requires B and flux-label differences below 1%, iota-profile differences below 1% of the
reference maximum, and zero absolute pressure defect (the vacuum pressure is fixed to
zero in all three codes). This is not a continuous error bound, an independent
total-flux integral or a force-residual match. “DOFs at 1%” means the least measured C
on this ladder, not a proven minimum over discretizations.

[V] In the table, `pK-nN-bB` means degree K, N toroidal harmonics and B boundary mesh
nodes, and the native boundary fit retains poloidal modes through 8. Error columns take the
larger value at one/eight threads. S marks the native timed-phase sum. The complete
per-point ledger retains both thread results and failed attempts.

| Code | Case | Resolution | C / free | B / p / flux-label (%) | Strong force | Δiota | Workflow s, 1/8 | Native call or S s, 1/8 | RSS MiB, 1/8 | C at 1% fields |
|---|---|---|---:|---|---|---|---|---|---|---:|
|KIN6D main|Landreman|p2-n4-b16|3,699 / 1,890|0.342 / 0.52 / 0.941|—|0 (input)|3.86/1.37|3.81/1.35 S|50.6/60.6|3,699|
|KIN6D main|Landreman|p3-n4-b8|2,295 / 1,098|0.162 / 0.396 / 0.748|—|0 (input)|4.21/1.96|4.17/1.94 S|49.6/62.9|2,295|
|KIN6D P1 pin|Landreman|p2-n4-b16|3,699 / 1,890|0.342 / 0.52 / 0.941|—|0 (input)|3.84/1.38|3.79/1.36 S|50.8/60.8|3,699|
|VMEC++|Landreman|m5-ns13 (worker cap 6)|2,379 / —|0.079 / 0.646 / 0.807|0.0832|4.44e-16|5.19/5.16|0.028/0.0233|441/455|2,379|
|GVEC|Landreman|m8 (Jacobian abort)|11,258 / —|— / — / —|—|—|≥69.9/≥18.4|—/—|475/499|—|
|DESC|Landreman|m4|443 / —|0.168 / 7.08 / 7.66|0.00644|0|47/26.8|31.6/14.8|2.49e+03/2.51e+03|—|
|DESC|Landreman|m8|2,877 / —|0.165 / 8.15 / 10.9|0.000323|0|275/154|259/142|1.05e+04/1.04e+04|—|
|KIN6D main|Issan|p2-n4-b16|3,699 / 1,890|0.257 / 0.435 / 0.711|—|0 (input)|5.41/1.93|5.36/1.91 S|50.4/66.6|3,699|
|KIN6D main|Issan|p3-n4-b8|2,295 / 1,098|0.115 / 0.335 / 0.55|—|0 (input)|3.87/1.85|3.83/1.83 S|49.9/68.8|2,295|
|KIN6D P1 pin|Issan|p2-n4-b16|3,699 / 1,890|0.257 / 0.435 / 0.711|—|0 (input)|5.23/1.93|5.19/1.9 S|52.4/60.8|3,699|
|VMEC++|Issan|m4-ns9 (worker cap 4)|2,214 / —|0.0835 / 0.702 / 0.843|0.188|4.44e-16|5.1/5.08|0.0361/0.026|437/449|2,214|
|GVEC|Issan|m8 (Jacobian abort)|11,258 / —|— / — / —|—|—|≥56.8/≥15.9|—/—|473/499|—|
|DESC|Issan|m4|891 / —|0.155 / 6.25 / 6.14|0.00181|0|245/105|231/94.2|5.86e+03/5.91e+03|—|
|DESC|Issan|m8 (600 s limit)|5,763 / —|— / — / —|—|—|≥600/≥600|—/—|3.63e+04/3.61e+04|—|
| KIN6D main | QA | current-driven seed unavailable | — | — | — | — | — | — | — | — |
|VMEC++|QA|m4-ftol12|2,091 / —|0.0418 / 0 / 0.269|0.0237|0.00023|5.4/5.44|0.0903/0.187|435/447|2,091|
|GVEC|QA|m4|1,089 / —|0.0172 / 0 / 0.234|0.0106|0.000326|19/16.8|17.1/14.9|243/251|1,089|
|DESC|QA|m4|443 / —|0.073 / 0 / 1.17|0.00243|0.00074|61/33.7|41.6/18.3|3.72e+03/3.71e+03|2,877|
|DESC|QA|m8|2,877 / —|0.00167 / 0 / 0.0541|2.23e-06|5.57e-06|408/207|388/191|1.71e+04/1.69e+04|2,877|

[V] QA pressure defect is absolute in mu0*p (T²), and QA strong force is absolute in
T²/m because the reference pressure scale is zero. Integer strong force is the existing
relative sampled norm. A missing force value is unavailable, and no algebraic-gradient value
is substituted.

[V] Landreman: the least measured passing native point uses 2,295 coefficients and 1,098
free unknowns, versus 2,379 VMEC coefficients (3.5% fewer by the stated coefficient
convention). P3 needs 2,295 coefficients, compared with 3,699 for P2 (38.0% fewer).

[V] Issan: the least measured passing native point uses 2,295 coefficients and 1,098
free unknowns, versus 2,214 VMEC coefficients (3.7% more by the stated coefficient
convention). P3 needs 2,295 coefficients, compared with 3,699 for P2 (38.0% fewer).

[V] The full matched-accuracy gate remains open. The measured main's common-core map
adapter rejects nonzero toroidal modes (`src/mhd/map_equilibrium.f90:62`), so no
comparable imported/native 3D D3 residual was measured here. The relative free Newton
gradient is an algebraic stopping test and does not replace D3 or bound field/pressure
errors. The integer comparator's strong force norm also omits interelement face
currents. Iota profiles are reported in the input orientation: DESC's positive-Jacobian
theta reversal and GVEC's phi=-zeta chart must be accounted for. The fixed iota and
total toroidal-flux constraints are not independent D7 readback integrals or field-line
traces. The QA current-driven nested initializer is still unavailable in KIN6D. The NN1
coil-vacuum QA study is a different boundary problem and was not substituted. G6's
non-nested reduction study is separate.

[V] GVEC cold-start positive-Jacobian failures are retained with their cost and the
basis of the last saved stage, while the resolved input retains the requested final
basis. The failure originates in `src/functionals/mhd3d/mhd3d_minimize.F90:471`. The
integer Landreman export also contains zero RBS/ZBC arrays, which GVEC's input
translator treats as permission for a nonsymmetric basis. Its actual expanded basis is
counted. No cold-start failure is presented as a converged equilibrium or a physical
error. P3 can reduce B error while worsening pressure/label error at the same mesh and
angular resolution, so B alone cannot select a matched point.

[V] The adapter records separate native-call and complete workflow timers, actual basis
sizes, affinity, solver status and parent/child peak RSS. Workflow time includes
imports, conversion, readback and output, while KIN6D's workflow also includes boundary
preparation. External native-call timers exclude those costs. For KIN6D, S is the sum of
disjoint energy/gradient, Hessian and factor/solve timers, a lower bound on the complete
Newton call. Whole-workflow totals include Python startup and external readback, which
can dominate short VMEC++ solves. RSS is the measured process maximum rather than the
requested allocation. The final adapter uses VMEC's `resize` API so its saved input has
consistent axis and boundary array dimensions. Native consumer samples before and after
this correction agree exactly for all three input families at m=4. Failed Python calls
retain an incomplete metadata record instead of disappearing.

[V] The QA stopping criterion matters independently of DOFs. At unchanged m=8 and ns=33,
tightening VMEC FTOL from 1e-8 to 1e-12 removes most of the flux-label discrepancy
against the original-input reference. The loose m=12 run is therefore not used as the
final QA reference. The m=4/m=8 tighter-tolerance points and independent radial ladder
are retained alongside the default settings. Those controls are also provided by the
separate [G6 adapter extension](https://github.com/itpplasma/kin6d/pull/72).

[V] The P1 comparison covers N=4 and b=16 at P2/P3. It is a comparison at fixed
discretization, rather than an independent search for P1's minimum DOF count. The
recorded oracle diagnostics agree with the main-pin runs, and one process per point cannot
resolve small timing changes.
