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
`gvec`. Never apply the integer-family comparator to the QA case. VMEC uses
radial splines from SIMSOPT; no DESC fit or secondary solve is substituted.
Samples exclude the axis and edge and use three normalized toroidal-flux
surfaces. Vector components are Cartesian. The verifier compares B and pressure
at the candidate's physical positions, the exact flux label and the strong
force defect. These are equal-weight sample diagnostics, not volume norms,
surface distances, a D3 dual norm or continuous error bounds.

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
c0ffe910938bd7fdec051b0a9eccfb629e0c9470. Four modes failed the unchanged 1%
boundary-radius gate; 8, 12 and 16 modes passed for both cases.
