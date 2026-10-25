# Landreman profile harness (measurement branch only, never a PR)

Overlay of the PF2 boundary/profile adapter (scratch kin6d-islands/PF2/adapter) on a
KIN6D head: `landreman_prof TAU MPOL NTOR BOUNDARY_NODES DEGREE`, built by
`include(benchmarks/landreman-prof/CMakeLists.txt OPTIONAL)`. Set
`KIN6D_PROF_SPECTRUM=1` for the dense reduced-Hessian spectrum/PCG probe
(`map_spectrum_probe.f90`). `gll/spectrum_gll.f90` is the same probe for the #82
GLL spike operator (add it to benchmarks/gll/CMakeLists.txt of an #82 tree).
