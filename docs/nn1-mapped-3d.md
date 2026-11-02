# NN1 on a phi-dependent domain (mapped 3D)

Status: candidate. Source `src/mhd/nonnested_mapped.f90`
(`kin6d_nonnested_mapped`), tests `test/test_nonnested_mapped.f90`. NN1 itself is
described in [equilibrium-nonnested.md](equilibrium-nonnested.md).

## Purpose

NN1 solves the vacuum field `B = F grad(phi) + curl A` on an axisymmetric
reference domain with a meridian FEEC complex times a Fourier basis in `phi`.
This step lets the physical domain depend on `phi`, for example a rotating
ellipse, without changing the complex.

## Map

The physical point is `x = F(xi)` with reference `xi = (R, Z, phi)`. The map keeps
`phi`; only `(R, Z)` are displaced:

    (R, Z)_phys = (R, Z)_ref + d(R, Z, phi).

`nn_map_t` stores `d` as Lagrange V0 coefficients times real Fourier blocks in
`phi` (the N1 `q1, q2` representation on the FEEC mesh). `nn_map_build`
interpolates an analytic map. `nn_map_jet` returns `DF` in orthonormal cylindrical
frames and `J = det DF`.

## Pulled-back metric

Magnetic 2-forms transform with the Piola map. The complex on the reference
domain is unchanged, so `div B = 0` holds algebraically for any map. Only the
Hodge (mass) tensor changes:

    K = DF^T DF / J.

`K` enters upstream `meridian_fourier_curl_curl` as a `meridian_fourier_tensor_t`
(FortFEM a27e21f7d498). The field is `B = curl A` on the reference domain; the
physical field is `DF B / J`.

## Solve

With `h = F/R e_phi` (the vacuum harmonic) the discrete problem is

    S A = -C^T M2(K) h,    S = C^T M2(K) C,

with the pulled-back potential `DF^T A + F ln R e_Z` as data (`nn_mapped_data`): an
L2 projection followed by the exact discrete gauge `R A_phi = 0` for `n > 0`
(`A - grad_F lambda`). `S` is singular on the gradient kernel; the right
preconditioner adds `eps M1` on the `n = 0` Nedelec block.

Solver: real FGMRES (fortnum) with a right preconditioner built from the
phi-averaged metric `Kbar`. That operator is diagonal in the Fourier mode `n`, is
factored once (`nn_factor_t`, the former `nn_linear_solve` split) and applied
per mode. The preconditioned operator is spectrally equivalent with ratio
`sup K / Kbar`, so the iteration count should be independent of `h` and `N`.

## Mode coupling

`K` is a rational function of the map, so in principle all `n, n'` couple. With a
map of Fourier degree `d`, the coupling is dominated by `|n - n'| <= 2d`, and with
`nfp` field periods only `n = n' mod nfp` couples. The `phi` quadrature needs
`nphi >= 2N + 1` for polynomial integrands; the tests use `4N + 4`, which is not
exact for rational `K`. The quadrature error enters the verification ladder and
must be separated from the discretisation error when a ladder stalls.

## Verification

Oracle: `B = F grad(phi) + B_v e_Z + sum_{m=1..3} grad(c_m R0 (R/R0)^m cos m phi)`
on the disk `R0 = 3, a = 1`; rotating ellipse with semi-axes 1.25 / 0.8 turning
with `phi` (map degree 2). Ctests (`ctest -R nonnested_mapped`):

1. `identity_regression`: identity map reproduces NN1 (nr=2, P2, N=3) to
   `max|b - b_NN1|/max|b| = 5.3e-14`, 1 FGMRES iteration, `div = 7e-15`; field
   error 1.479e-2 vs the oracle, identical to NN1.
2. `ellipse_h_ladder` (`hladder 10 4 8 16`) and 3. `ellipse_mode_ladder`
   (`modes 2 1 3`).

Measured (local, 4 threads, 2026-10-13):

| nr | N | ndof | field err | energy rel | FGMRES | residual | div |
|----|---|------|-----------|------------|--------|----------|-----|
| 2 | 5 | 1633 | 6.27e-2 | 5.5e-4 | 16 | 3.3e-11 | 5e-15 |
| 4 | 5 | 6001 | 2.48e-2 | 1.7e-4 | 17 | 8.4e-11 | 5.5e-15 |
| 8 | 5 | 22033 | 9.22e-3 | 4.6e-5 | 17 | 5.0e-10 | 4.7e-15 |

Observed field order 1.34 (nr 2->4) and 1.43 (4->8), rising but not yet the
expected 2 for P2. The mode ladder at nr=2 falls 1.09 (N=1), 0.105 (N=3), 6.27e-2
(N=5), 6.25e-2 (N=7): it saturates at the spatial error from N=5. Iterations are
flat in `h` and `N` as predicted; `div B` stays at roundoff.

The CTest h ladder runs at N = 10 (`hladder 10 4 8 16`) with the program gate 1.7 on the
finest pair (`KIN6D_MAPPED_MIN_ORDER` default): measured 1.64 (nr 4->8) and 1.80 (nr 8->16).
The order loss at N = 5-6 was mode truncation (docs/nn1-qa-benchmark.md).

## Open items

- Matrix-free apply. Upstream `M2(K)` is element-dense in modes, `O((nl(2N+1))^2)`
  per element (about 1.5 GB of triplets at nr=16, N=6). Production needs a
  matrix-free apply by FFT in `phi` (grid values of `K`, FFT per quadrature
  point, OpenMP over elements). It belongs in FortFEM, not here.
- `csc_matvec` is serial. This is a defect (serial hot path) in the sparse layer.
  The repair belongs in fortsparse/fortfem; it is not part of this step and no
  private copy is made.

## Consequences for NN2, NN3 and SPEC/SPECTRE

- NN2: the same `K` enters its mass and transport forms; the map interface is
  reusable unchanged.
- NN3 takes `X` from the N4 shell or from an N1 map.
- A fixed-boundary rotating-ellipse case for the SPEC/SPECTRE comparison becomes
  possible: the boundary is the image of the reference disk under the map, and
  the vacuum field here is the reference for the vacuum region.

## Fixed full Cartesian chart operators

`nn_map_set_chart(chart,map)` owns a copy of the existing `ideal_chart_t`
contract `H(q)`, with `q=(R,Z,phi)`. It is an alternative to the retained
Fourier displacement, not an admission certificate. At every evaluated point,
`nn_map_point` returns the physical Cartesian position and covariant/contravariant
Piola matrices acting on reference orthonormal cylindrical components ordered
`(R,phi,Z)`. This order differs from the left-handed coordinate tuple `q`.
For `R>0`, define

    D = [H_R, H_phi/R, H_Z],  J = det(D) = -det(dH/dq)/R,
    P1 = D^(-T),  P2 = D/J,
    K1 = J P1^T P1,  K2 = J P2^T P2.

The existing generated FortFEM Piola implementation owns the inverse and
local mapping arithmetic. `map_metric_t%space=1` selects K1 and its default
`space=2` selects K2 through the existing phi-dependent tensor API. The
cylindrical chart has D equal to the Cartesian cylindrical rotation, J=1,
and both Hodge tensors equal to the identity. Nonfinite or nonpositive
orientation is rejected, including the provider's existing singularity guard.
No absolute determinant repairs an incorrectly oriented chart. Positivity at
sampled points does not prove global injectivity, embedding or regularity.

`nn_mapped_cartesian_data` takes physical Cartesian positions and returns
Cartesian potential components from its callback. `nn_mapped_cartesian_error`
uses the same explicit convention for its exact-field callback. The original
`nn_mapped_data` and `nn_mapped_error` keep physical `(R,Z,phi)` samples and
cylindrical vector components and reject a full-chart map. The common
`equilibrium_state` physical sampling/export contract is unchanged; generic
chart inversion and a whole mapped NN2 state are not provided by this slice.

The new data entry reuses the existing V1 projection and discrete gauge,
with pullback `D^T A` and the same reference harmonic compensation. The
existing `nn_mapped_solve` then solves the metric vacuum equation with K2,
including its harmonic/cohomology load. A Piola pushforward of `F/R e_phi`
without that solve is only a transported cohomology representative and is not
claimed to be a physical curl-free vacuum field.

Both error entries use the toroidal mean measure
`J R dR dZ dphi/(2*pi)`, not an extra full-turn or field-period multiplier.
The reported error is the absolute physical L2 norm and the two energies are
`mean_integral(|B|^2/(2*mu0))`. Coordinates are metres, field components tesla,
potential components tesla metres and J/P1/P2 are dimensionless when H and
reference lengths have the same unit. No rigorous remainder is attached to
these empirical errors or quadrature comparisons.

This is a fixed-geometry operator composition. Field transforms are linear:
tangents use the same P1/P2 and cotangents their transposes. Existing exact
incidence/curl and generated point leaves remain their authorities. Chart
parameter optimization, an NN2 nonlinear mapped equilibrium, moving physical
interfaces, full exterior BEM and GPU execution are separate unqualified
consumers. The assembled NN1 modal cost is unchanged: element storage scales
quadratically with the number of Fourier blocks. Point mapping adds constant
3-by-3 scratch and no per-point heap allocation; this is not a performance
claim or a new private tensor backend.
