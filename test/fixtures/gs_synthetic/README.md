# Synthetic Grad–Shafranov regression

These inputs are constructed analytically and contain no machine data.
The historical CTest name `gs_k6_mini_regression` is retained. Its cubic-profile
Newton, P3 direct solve, current normalization, readback, field-query and majorant
path is exercised by a Solov'ev solution with R0=3 m, a=0.6 m, C=0.25 and F=6 T m:

    psi = C*((R^2-R0^2)^2 + 4*R^2*Z^2) - 4*C*R0^2*a^2
    Delta*psi = 16*C*R^2; pprime = -16*C/mu0; FFprime = 0
    R(t) = sqrt(R0^2 + 2*R0*a*cos(t)); Z(t) = R0*a*sin(t)/R(t)

The existing FortSym identity test `grad_shafranov_symbolic_identities` retains
the exact operator/field identities for this polynomial. The independent
`synthetic_analytic` check in `test_equilibrium_core` compares every queried
field component with its analytic value. Its 0.1% physical-scale gate includes
coarse-mesh and spline-boundary error; it is not a continuous certificate.
The separate version-regression gate is 1e-12, matching the requested nonlinear
weak residual. Flux-estimator differences use the solution flux scale, and
the symmetry-zero axis height uses the major radius as its absolute scale.

Generate the inputs by compiling `test/support/generate_gs_synthetic.f90` with
gfortran and running it in this directory. The source has 513 normalized knots;
its constant derivative law is independent of the solved axis. The prescribed
current is an independent periodic quadrature of the exact boundary integral.
Query points lie on the half-radius analytic surface. Physical source data for
K6 production runs remain in the private TC24 repository.

The reference files were captured from the unchanged production solver at
6ac7bd4 on sCluster node4, job 3309757, with gfortran 12.2.0, OMP_NUM_THREADS=1
and OPENBLAS_NUM_THREADS=1. The gs_phase1 executable SHA256 was
262d0f2e407ec1dd206abdbdc53d51e4c715338c01569b8e6a423fe77b40b7b4.
The analytic scaled field error was 5.664e-7. To regenerate references, run
`gs_phase1 --config input.nml`, then
`gs_phase1 --evaluate mesh.dat 6 0 profile-law.dat < points.dat > field.dat`.
Copy field.dat/profiles.dat to reference_field.dat/reference_profiles.dat and
retain the physical/nondimensional keys already listed in reference_summary.dat;
exclude timings, provenance strings and solver diagnostic status strings.
