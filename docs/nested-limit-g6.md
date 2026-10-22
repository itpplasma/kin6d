# Nested-limit agreement: G6

[V] **G6 remains open.** This regression verifies the axisymmetric Solov'ev
limit of the relaxed A+D solve against an independent analytic oracle. It
does not establish agreement with VMEC++ or GVEC, or a solved public
stellarator case. Errors are observed diagnostics, with continuous bounds
and certification unavailable. The thresholds in
[island-equilibrium.md](island-equilibrium.md#verification-ladder-before-machines)
are unchanged.

## Matched axisymmetric problem

[V] The existing synthetic analytic fixture
`test/support/nn_solovev_fixture.f90` has R0=6.2 m, a=0.62 m, F=32.86 T m,
c=F/(12 R0^3), mu0=1 and p'=-16 c. Its wall is
R=sqrt(R0^2+2 R0 a u cos(theta)), Z=R0 a u sin(theta)/R at u=1.
The equilibrium has psi=4 c R0^2 a^2 (u^2-1), p=p' psi and zero boundary
pressure. Newton starts with an interior displacement of 0.1 m and V=0.

[V] The frozen A+D data are S=-chi_perp p' laplacian(psi),
g_parallel=g0 (R/R0)^2, I_phi=p' integral R dR dZ,
chi_parallel/chi_perp=1e6 and g_parallel/g_perp=1e6. They admit the same
analytic ideal state at finite ratio, with Phi=0 and V=2 pi R0^2 p'/g0.
This special nested cancellation does not establish the ideal limit for
general sources or current profiles. The existing current equations and
source fixture are reused without changing the closure.

[V] `benchmark_g6_solovev_input` exports this wall, toroidal flux, pressure
and transform in normalized toroidal flux s. Pressure is divided by
4 pi 1e-7 for VMEC's SI convention. NCURR=0 prescribes matched iota(s),
without adding an independent current constraint. The wall uses m<=16 and
512-point Fourier projection. Its staggered-grid defect is 2.15e-14 minor
radii. The 65 profile knots use cubic-spline interchange, whose interpolation
error is not certified. The independent flux integral is 6.5236848 Wb and
the on-axis iota is 2/3.

## Native regression and diagnostics

[V] `nonnested_g6_solovev_readback_p3` runs P3, eight rings and N=0 through
`nc_solve`, `free_b_field_t%sample` and `nn_residual`. It enforces sampled
B,p errors <=1e-3, relative toroidal-flux error <=1e-4, absolute iota error
<=1e-3 and relative unintended B,p coefficients <=1e-10. Existing Newton
and algebraic-divergence checks remain active. This uses D7's point
interface. The free-B type does not implement the common labelled-field
`read_equilibrium` contract.

[V] B,p errors use the existing 128 frozen points. B uses analytic
poloidal-field RMS as denominator, which is stricter here than total-field
RMS. Pressure uses analytic RMS. Toroidal flux integrates sampled B_phi
over the physical meridional mesh. Iota is reconstructed from B on numerical
psi contours corresponding to analytic u=0.30, 0.55 and 0.80. Bisection
locates each contour and checks the root. A 512-point angular integral of
dphi/dtheta=B_phi/(R B.grad theta) gives q and iota=1/q. This is an
axisymmetric contour diagnostic, not a general 3D surface finder.

[V] Local Release, GNU Fortran 16, four OpenMP threads, eight rings:

| diagnostic | observed |
|---|---:|
| relative B RMS | 5.010e-5 |
| relative p RMS | 3.334e-6 |
| relative toroidal flux error | 3.077e-8 |
| maximum absolute iota error | 6.092e-5 |
| normalized common-core D3 | 4.042e-6 |

[V] The G6 CTest and `kin6d_layering` pass. The existing registered P2
Fourier-reduction test also passes: unintended N=1 coefficient content is
9.16e-17 and D3 is unchanged from N=0. This is a separate reduction check,
not P3 mode-convergence evidence.

[V] sCluster job 3355111 used GNU Fortran 12.2, four threads and the native
FortSparse/AMD backend, after successful build pilots 3355093 and 3355104.
The following spatial ladder is P3, N=0 and ratio 1e6. Iota denotes the
largest absolute error on the three contours, not an error enclosure.

| rings | B RMS | p RMS | relative flux | iota | D3 |
|---:|---:|---:|---:|---:|---:|
| 8 | 5.010e-5 | 3.334e-6 | 3.077e-8 | 6.092e-5 | 4.042e-6 |
| 16 | 6.681e-6 | 1.211e-7 | 1.924e-9 | 3.123e-6 | 3.780e-7 |
| 32 | 1.055e-6 | 4.494e-9 | 1.203e-10 | 6.435e-7 | 3.206e-8 |

[V] Observed endpoint orders are 2.78 for B, 4.77 for sampled p, 4.00 for
flux and 3.49 for D3. These are fixture measurements, not certified rates.
The independent mode ladder at eight rings gives:

| maximum n | unintended B coefficients | unintended p coefficients |
|---:|---:|---:|
| 0 | 0 | 0 |
| 1 | 8.571e-16 | 2.130e-16 |
| 2 | 9.781e-16 | 2.814e-16 |

[V] B,p and D3 agree to the printed digits across these mode levels. The
largest iota error changes by 2.4e-11 absolute from N=0 to N=1. Adding
toroidal modes to axisymmetric data produces no accuracy gain. Timed
commands took 4.48/9.63/34.84 s for the spatial ladder and
4.55/16.27/34.01 s for N=0/1/2; the largest per-command peak RSS was
1.86 GiB. These are analytic-oracle differences, not VMEC/GVEC differences.

## External comparison and public stellarator admission

[V] The optional `benchmarks/nested/run_reference.py` adapter accepts
VMEC-only `--axisymmetric`, `--ns` and `--ftol`, separating radial and mode
resolution. It retains the resolved input, state, native residuals, iota,
toroidal flux and Cartesian samples. A sixth argument to the native `g6`
command reads these samples and checks total-vector B and mu0 p relative
RMS <=1e-3 at the external solver's physical points. Empty, malformed,
nonfinite or out-of-domain samples fail. This alone cannot close G6:
native convergence, matched flux/iota and common-core residuals for both
candidates must also be checked. VMEC's force norm is not D3.

[V] The public stellarator input is
[`input.landreman-paul-qa`](../benchmarks/nested/input.landreman-paul-qa),
from the [pinned public input](https://github.com/mishapadidar/stellarator_multiobjective/blob/be13ff521f59f3f6430d2900f1a5ee399ad36946/vmec_input_files/input.LandremanPaul2021_QA).
Its SHA256 is
e8fd222daf6e37e5dad7d9d38646ac47bb99c9d7a9a6cf262063b8a6e153b9cb.
It is fixed boundary, NFP=2, p=0, zero net current and Phi_edge=0.08385727554 Wb.
It needs no private machine data.

[V] Current A+D has no `nn_map_t` argument: `nc_problem_t` in
`src/mhd/nonnested_coupled.f90` supplies sources and conductivities, while
`element` samples cylindrical meridian bases and F/R. The mapped QA wall
and metric live in the separate NN1 vacuum path `nn_mapped_solve` in
`src/mhd/nonnested_mapped.f90`. Its published-wall vacuum result is not an
A+D result. Its coil field also has nonzero normal field on the published
wall, as measured in [nn1-qa-benchmark.md](nn1-qa-benchmark.md), so it cannot
replace the fixed-boundary VMEC equilibrium.

[V] GVEC is restricted to a validated exact-start state for this comparison.
The existing `run_reference.py gvec` branch calls `gvec.run(parameters)`
without a restart or exact-state initializer and is therefore not used for
G6. No admitted exact-start state for these cases is supplied here.

[H] Closing G6 needs mapped A+D with the same QA wall and flux, a validated
exact-start GVEC adapter, and common D3/readback support for imported
candidates. Then run at least three independent h and mode levels per
solver, reporting differences and quadrature/profile/trace sensitivity
separately and retaining failed thresholds without retuning.

## Reproduction

[V] Configure the CPU preset with AMD enabled and an out-of-tree build:

```sh
cmake --preset cpu -B /mnt/storage/kin6d-g6nested-build -DFORTSPARSE_ENABLE_AMD=ON
cmake --build /mnt/storage/kin6d-g6nested-build --target \
  test_nonnested_coupled benchmark_g6_solovev_input -j4
ctest --test-dir /mnt/storage/kin6d-g6nested-build \
  -R '^(nonnested_g6_solovev_readback_p3|nonnested_coupled_fourier_reduction_p2|kin6d_layering)$' \
  --output-on-failure -j4
/mnt/storage/kin6d-g6nested-build/benchmark_g6_solovev_input /tmp/input.g6-solovev
```

[V] In a Slurm allocation with VMEC++ and SIMSOPT installed:

```sh
python benchmarks/nested/run_reference.py vmecpp /tmp/input.g6-solovev 16 OUT \
  --axisymmetric --ns 65 --ftol 1e-12
BUILD/test_nonnested_coupled g6 3 1e6 8 0 OUT/samples.txt
```

[V] `scripts/study_g6_nested_scluster.sbatch ROOT pilot` configures AMD and
builds/checks the actual binary. After verifying that pilot and its peak,
submit `ROOT study` for P3 nr=8/16/32 at N=0 and N=0/1/2 at nr=8. Each task
reserves four cores and uses the cluster memory guard at 47.2 GiB. Failed
fixture gates stop the study. Environments and outputs stay outside Git.
