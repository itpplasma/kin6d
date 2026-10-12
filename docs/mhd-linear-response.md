# Linear per-n MHD response (L4)

Status: first stages of L4 of the 3D design (`kin6d_3d_tc24/DESIGN.md`, §2, §8.2).
The ideal response is computed in the FEEC form of Albert, Lainer and Biro,
Comput. Phys. Commun. 277 (2022) 108401, as used by MEPHIT, not in the
ξ-form of the design: the plasma current is obtained from linear force balance
on flux surfaces and the field from the per-n Ampère operator.

## Model

A perturbation is `X(R, Z) e^{i n φ}` (physical field `X e^{i n φ} + c.c.`,
the convention of `kin6d_coils_3d`), `n ≠ 0`, static (`ω = 0`), on an
axisymmetric G-EQDSK equilibrium.

- **Coil field**: closed-form Biot–Savart of the 3D coils (L3), projected on
  `e^{i n φ}`.
- **Response field** of the plasma currents: `B = curl_n A`, `A` in the
  Nédélec space in the gauge `R A_φ = 0` (exact for `n ≠ 0`); outside the
  plasma `B = grad_n u` with the exact exterior Steklov–Poincaré map `S_n`
  (`fortfem_toroidal_mode_boundary_operators_2d`). The operator
  `K + C^H S^{-1} C` (`C` the normal-flux pairing,
  `fortfem_meridian_magnetic_coupling`) is Hermitian positive definite and
  factored once per mesh and `n` (`kin6d_mhd_ampere_n`).
- **Shielding sheets** on the rational surfaces `q = m/n`: surface currents
  `n̂ × grad_n Λ`, `Λ = λ e^{i(nζ − mθ_s)}`, with `λ` fixed by zero total
  resonant flux (`kin6d_mhd_response_n`).
- **Volume current** (`model = 'ideal'`, `kin6d_mhd_ideal_n`):
  `B0·∇p = −p0' B·∇ψ`, `J⊥ = B0 × (∇p − J0 × B)/B0²`,
  `B0·∇(J∥/B0) = −∇·J⊥`, with the magnetic differential equations solved
  exactly in the geometric-angle Fourier basis. The current leaving through
  the plasma boundary continues in a field-aligned boundary sheet. Surfaces
  with `|q| < qlow` carry no volume response (the DCON `qlow` start).
  The fixed point `a = K⁻¹ L(B(a))` is solved by GMRES (fortnum) with the
  shielding sheets eliminated.

Not modelled: the pressure jump of a cut equilibrium at the boundary (vacuum
outside with `p0_edge > 0`), rotation, resistivity, kinetic terms.

## Discretisation

Flux-aligned mesh (`kin6d_mhd_mesh`): rings at `u = sqrt(ψ_N)`, uniform plus
every rational surface; vertices per ring double outward; isoparametric
geometry of order ≤ 2 whose ring-edge nodes lie on the flux surfaces. FEEC
order `k` (P_k Nédélec). The volume current of the strip between two rings is
that of its mid surface, sampled at `np_sample` uniform angles.

## Verification (ctest `mhd_*`)

| Test | Oracle | Result |
|---|---|---|
| `mhd_ampere_n_exterior_manufactured` | current `curl_n W` of a compact bump; field from toroidal-angle quadrature of `1/(4π|x−y|)` | field order k (k = 1, 2), exterior trace order ≥ 2 |
| `mhd_shielding_cylinder_limit_reciprocity` | screened response of a periodic cylinder (modified Bessel functions); sheet reciprocity | 1e-5 at R0/a = 40, second order in a/R0; M Hermitian to 1e-13 |
| `mhd_ideal_n_mde_force_free_identity` | MDE inverse; force-free `μ0 J = σ B` against a direct FEEC solve; charge conservation with `p' ≠ 0`, `F' ≠ 0` | 5e-10; second order in the ring spacing (8.5e-5); 3e-5 |

`gs_response_3d` also recomputes the coil resonant flux through the L3
surface path (`surface_harmonics`), agreeing to 1e-5.

## Usage

`gs_response_3d input.nml`, namelist `/l4/`: `eqdsk`, `coils`, `currents`
(one table per phasing), `n`, `rat_num`/`rat_den`, `nrad`, `nth_max`,
`order`, `model` (`'shield'` or `'ideal'`), `np_sample`, `qlow`,
`boundary_sheet`. Output `prefix_res.dat`: per phasing and rational surface
the vacuum resonant flux `Φ_vac` and `b_vac = Φ_vac/area`, the sheet amplitude
`λ` and `b_shield = −M_ss λ/area`, the resonant flux of the shielding current
alone per area (GPEC's `singflx`); `m = n q` with the signed `q` of
`kin6d_vacuum_surfaces` (GPEC labels it `−m`).
