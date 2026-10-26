# Landreman profile harness (measurement branch only, never a PR)

Overlay of the PF2 boundary/profile adapter (scratch kin6d-islands/PF2/adapter) on a
KIN6D head: `landreman_prof TAU MPOL NTOR BOUNDARY_NODES DEGREE`, built by
`include(benchmarks/landreman-prof/CMakeLists.txt OPTIONAL)`. Set
`KIN6D_PROF_SPECTRUM=1` for the dense reduced-Hessian spectrum/PCG probe
(`map_spectrum_probe.f90`). `gll/spectrum_gll.f90` is the same probe for the #82
GLL spike operator (add it to benchmarks/gll/CMakeLists.txt of an #82 tree).

Solver knobs (environment): `KIN6D_PROF_KSTEP=1` Krylov step, `KIN6D_PROF_LAGGED=1`
lagged Hessian, `KIN6D_PROF_FPOLICY` factor policy (default LDL^T in the library),
`KIN6D_PROF_FTHREADS` factor team. `KIN6D_PROF_CONT='ntor degree boundary_nodes tol'`
sets the continuation ladder. Without it the harness solves the same ntor and mesh
at degree - 1 to 1e-2 first (decision 126; '4 2 12 1e-2' on b12), and
`KIN6D_PROF_CONT=off` solves without a ladder.
