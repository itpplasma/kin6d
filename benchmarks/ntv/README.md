# NTV torque and gradient benchmark

`lorentz-spectrum64.in` is a synthetic circular-field case with 64 perturbation
harmonics, 129 poloidal points and 128 Legendre polynomials. It uses the full
magnetic-drift weight, Boozer streaming and pitch-angle scattering. It is a
kernel benchmark, not a full-collision or device validation case.

Build and run the existing native driver:

```sh
cmake --preset cpu
cmake --build --preset cpu --target test_dk_ntv -j 32
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    /usr/bin/time -v ./build/test_dk_ntv kernel benchmarks/ntv/lorentz-spectrum64.in
```

The driver computes torque and all perturbation-amplitude, collision-rate and
rotation sensitivities. Its printed columns are the torque kernel, angle and
pitch refinement differences, and the collision-rate and rotation derivatives.
The two refinement columns are zero in this timing input because refinement
solves are explicitly disabled. They are not zero error bounds. Change the
third entry on the resolution line from 0 to 1 to include those solves.

Compare revisions on the same machine, compiler, dependencies, thread counts
and input. Warm each executable, alternate the revisions, retain repeated wall
times and peak memory, and check the numerical outputs before comparing cost.
The registered `dk_ntv_*` tests independently check the operator and derivatives;
timing does not replace them. A device time-to-accuracy claim additionally
requires convergence in angle, pitch and energy and the full collision model.

`timings.csv` retains the repeated local measurements; `measurement.json`
records hardware, compiler, flags, baseline revision and scope. The `before`
variant is the baseline, `after` replaces parameter-wise tangent solves by one
torque adjoint, and `blas` also uses BLAS for the dense Schur products.
`equivalence.json` compares each scalar output separately and the amplitude
gradient in its maximum norm, so quantities with different units are not mixed.
These arithmetic comparisons accompany the independent derivative tests; they
are not discretization-error estimates.

The spectrum inputs cover 1, 16 and 64 harmonics. To reproduce the local-stencil
rows, change the final entry of the resolution line from 0 to 6. Keep the other
entries fixed when comparing the two implementations of that discretization.
