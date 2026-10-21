# CPU precision experiment

[V] This opt-in benchmark calls GS source/tangent and ideal-map Hessian leaves
emitted from the existing FortSym definitions. It compares fp32 arithmetic and
portable fp16/bf16 storage with fp32 tile arithmetic against fp64. It leaves
production candidates, final residuals and estimators unchanged. The storage
codec is the shared FortGen fwprec implementation, pinned on its feature branch
in `cmake/kin6d_cpu_precision.cmake`; integration to main remains pending.

Configure the ordinary KIN6D tree with
`-DKIN6D_BUILD_CPU_PRECISION_BENCHMARK=ON`, or build only these same microkernels:

```sh
cmake -S benchmarks/cpu-precision -B build-cpu-precision -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS_RELEASE='-O3 -g -march=native -fopt-info-vec-all'
cmake --build build-cpu-precision --target benchmark_cpu_precision -j 4
benchmarks/cpu-precision/run.sh build-cpu-precision /absolute/result/path
```

On Apple ARM use `-mcpu=native` instead of `-march=native`. The runner accepts
`BUILD OUT 'THREAD COUNTS' SIZE REPEATS`. It runs finite work, records the binary
identity, invokes registered CTest comparisons and captures wall time/RSS.
Regenerate with the pinned FortSym build and `gen_grad_shafranov` / `gen_ideal_map`.

[V] The fixtures are well-conditioned tables near R=3 and the identity map,
not quadrature tables from a solved equilibrium. Independent analytic checks
cover the GS source/tangents and three identity-map Hessian entries. Scaled
maximum output defects must be below 1% for fp32 and 3% for half storage.
These comparisons are empirical indicators; a continuous certificate and a
rounding majorant are unavailable. Half arithmetic is not claimed: the portable
shared codec's unpack and repack cost is included in each timed tile. Initial
packing/casting is reported separately. Trials alternate fp64/fp32 order; times
are minima over repetitions, with no significance claim on shared machines.

[V] With GNU Fortran, the benchmark raises the inlining limit so the large
Hessian leaf can enter the SIMD loop. Verify the actual compiler report on
each architecture. A branch-only long-line option supports the batch emitter's
remaining long call sites; scalar-leaf wrapping already comes from the pinned
upstream fix. No private storage or precision framework is implemented here.

[V] The ordinary KIN6D build also supplies `benchmark_equilibrium_cpu NR NTOR
REPEATS`. It profiles energy/gradient, Hessian assembly and the fp64 D3 residual
on the existing Solovev fixture's initial candidate. It does not solve that
candidate or claim an equilibrium error bound. `gs_load_bench` profiles the real
2D load; use a quadrature tolerance such as 1e-3 and compare its output arrays
with a tighter reference. K6-mini remains the native solved-state and estimator
regression. [V] The branch-only `benchmark_cpu_solve NR REPEATS` consumes the unchanged
fwprec precision/GMRES-IR sources when the ordinary build is configured with
`-DKIN6D_CPU_PRECISION_FRAMEWORK_ROOT=/path/to/fwprec/checkouts`. That root must
contain the owning `fortnum` and `fortsparse` checkouts. It solves the actual GS
CSR matrix, checks the result through the production exchange-form operator,
and compares the state to the native GS solution at a 1% target. It reports the
actual final storage/compute policy, including any escalation. This is a
prototype source seam until those library APIs are published; record their
source hashes with the binary. Production mixed-precision integration remains
pending.
