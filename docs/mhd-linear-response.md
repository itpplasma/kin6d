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

Cut boundary. With `p0_edge > 0`, vacuum outside and `B0` continuous (the
cut G-EQDSK), the linearised force balance has a `p0_edge ξ_n δ'` term that
no finite field perturbation balances: the cut is not an equilibrium
plasma–vacuum interface, and DCON/GPEC drop this term as well. The traction
`ξ_n [[∂_n(p0 + B0²/2μ0)]]` vanishes identically for continuous `B0`. The
remaining boundary terms are the field-aligned closure sheet and, with
`edge_moved_current`, the moved current `ξ_n J0` of the displaced boundary
(DESIGN §2.5). TC24 results depend strongly on this choice.
Not modelled: rotation, resistivity, kinetic terms.

## Discretisation

Flux-aligned mesh (`kin6d_mhd_mesh`): rings at `u = sqrt(ψ_N)`, uniform plus
every rational surface; vertices per ring double outward; isoparametric
geometry of order ≤ 2 whose ring-edge nodes lie on the flux surfaces. FEEC
order `k` (P_k Nédélec). The volume current of the strip between two rings is
that of its mid surface, sampled at `np_sample` uniform angles.

## Error diagnostics

Every `ampere_n_solve` stores `amp%error`. Its `discrete_energy_error` is
`(r^H (K + C^H S^-1 C)^-1 r)^(1/2)`, evaluated with the retained coupled
operator and one correction solve using the existing factorization. It measures
algebraic error in the interior-plus-exterior energy norm with the same
`R dR dZ` pairing and toroidal amplitude convention as the solve. Ordinary
arithmetic and factorization errors are not enclosed, so its claim is
`indicator`, including when the value is near zero.

Passing `fs`, an enriched `fsd` on the same mesh and quadrature, and a
thread-safe `current` callback to `ampere_n_solve` also computes FortFEM's
per-mode equilibrated reconstruction. The callback returns the already
`mu0`-scaled smooth volume current in T/m. The result is available as
`amp%error%interior`, including per-cell indicators. The same diagnostic for an
arbitrary candidate is exposed by `ampere_n_equilibrated_error`.

This reconstruction estimates the **interior problem with the candidate's
tangential trace fixed**. It does not include surface-current loads or the
exterior Steklov remainder. Consequently `interior_claim` is
`indicator-fixed-tangential-trace`, `continuous_claim` remains
`unavailable-exterior-remainder`, and no finite exterior radius is reported.
`ampere_n_coupled_error` instead imposes the sheet jump and the boundary
dual trace in a broken scalar potential `sigma = R H_phi`. It bounds the
coupled problem with the **Galerkin** exterior and interpolated sheet data.
Its zero Steklov term is an algebraic identity for that exterior operator;
it does not bound its difference from the exact exterior map.

An optional thread-safe `sheet_derivative` callback returning `dG/dtau`
adds `sheet_remainder` to `coupled_bound`, with `sheet_status=practical-lift`.
For each oriented sheet edge, set `r = G - I_k G`. On its listed adjacent
triangle the potential `sigma_r = -(lambda_a + lambda_b)^4 r(tau)`,
`tau = lambda_b/(lambda_a + lambda_b)`, vanishes on the other two edges.
It has jump `-r`. Summing these lifts gives
`H_r = (grad(sigma_r)/(i n), sigma_r/R)`, satisfying
`(curl_n v,H_r)_R = L_r(v)` by cellwise integration by parts. Testing the
correction equation with its solution gives `|||A_G-A_IG||| <= ||H_r||_R`,
also with a nonnegative exterior energy. This is a continuous inequality
with constant one. The numerical norm uses ordinary sampled quadrature,
so the returned correction remains practical. Missing callbacks and
inconsistent shared-node jumps have explicit unavailable statuses. Without
the callback the original bound still concerns interpolated sheet data.

The same integration-by-parts identity gives constant one for the interior
Prager--Synge reconstruction: `||curl_n e||^2 + ||H-B||^2 = ||H-B_h||^2`.
No unknown Friedrichs constant is needed. If equilibration has a volume
defect `d`, the gauge implies `||v/sqrt(R)|| <= ||curl_n v||_R/|n|`, hence
the additional term is `||R^(3/2) d||/|n|` (or `Rmax ||d||_R/|n|`).
The nonnegative exterior term only strengthens this estimate. These
constants do not enclose quadrature or floating-point evaluation.

In the small toroidal ring fixture the new sheet term costs 17--41% of a
reused Ampere solve (including its residual diagnostic), or below 0.05% of
boundary assembly, factorization and solve. Three repeated measurements
were taken on the shared host. The independent cylinder test's eight loaded
cells give about 2x speedup from one to four threads, with unchanged radii
to `1e-12` relative. This is evidence for the added term only; the existing
global potential reconstruction and dense BEM path are not O(N) algorithms.

The exact-exterior remainder still requires a computable continuous trace
residual norm and its reliability constant; the assembled BEM matrices
alone do not provide them. A small coupled algebraic residual does not close
that gap. The derivative enclosure also remains
explicitly unavailable; the primal diagnostic does not certify response
sensitivities. Mode zero still requires a separate gauge and is rejected.

`ideal_response(..., error=report)` audits its returned GMRES vector with the
same production current and sheet elimination used by the solve. Write
`F(a) = K^-1 L(B(a))`, with `K` the coupled discrete Ampere operator, `T` the
homogeneous response map, and `d = F(a)-a`. The report's `frozen_current`
contains the algebraic Ampere diagnostic for the load frozen at this final
state, and `defect_energy = ||d||_K` up to solve/rounding error. It replaces
the diagnostic of the last internal solve in `amp%error` when requested.

The theorem `||a*-a||_K <= ||d||_K/(1-q)` requires a genuine upper bound
`||T||_K <= q < 1`. Instead, this implementation measures only the direction
`q_probe = ||T d||_K/||d||_K`; that ratio can miss nonnormal amplification
and is not an operator-norm upper bound. With `q_probe < 1`,
`algebraic_estimate` uses the geometric-tail formula with claim
`practical-directional-contraction`. A noncontractive probe returns
`unavailable-noncontractive-probe` and no finite tail estimate. The exact
zero-current/zero-response case returns zero without dividing by zero.
This is an equation-level diagnostic, independent of GMRES's stopping norm;
it can also describe a candidate whose GMRES iteration limit was reached.

This norm refers to the volume-response coefficients, with the discrete
exterior energy included. It does not bound the error of the combined volume
and resonant-sheet field, sheet amplitudes, or any pointwise observable.
`discretization_claim = unavailable-current-reconstruction`,
`sampling_claim = unavailable-current-sampling`, and
`continuous_claim = unavailable-ideal-stability` keep the missing parts
separate. In particular a continuous curl_n reconstruction for the frozen
load alone would not control the feedback of the ideal current. No rigorous
algebraic or continuous ideal-MHD certificate is claimed. Smooth scalar
ratio/tail leaves and their JVP/VJP are generated by
`gen/gen_ideal_response_error.f90`; the derivative of the complete report
remains unavailable. The diagnostic reuses the parallel current and load
assembly; field sampling over rings is OpenMP-parallel. Additional retained
state is linear in coefficient and sample counts, but inherited factored
solves prevent an O(N) time claim.

Measured on 2026-10-15, force-free n = 2, sigma = 0.5, P2, 128 angular
samples, three trials per mesh on the shared host:

| Volume DOFs | Discrete algebraic index | Diagnostic seconds (4 threads) | Diagnostic / solve |
|---|---|---|---|
| 592 | 1.087 | 0.00198 | 0.148 |
| 2064 | 1.181 | 0.00734 | 0.123 |
| 7184 | 1.087 | 0.0366 | 0.126 |

The indices use loose tolerances 1e-2/1e-3/1e-4 against a 1e-12 discrete
reference. Timings and paired ratios are medians at tolerance 1e-10
(14 GMRES iterations),
excluding mesh/operator construction and factorization. At 7184 DOFs the
one-thread diagnostic median is 0.0339 s: 1-to-4 speedup is 0.93, with
individual-trial ratios 0.79--1.24. This test establishes no useful strong
scaling; all printed error values agree across thread counts. Continuous
efficiency and an O(N) runtime bound remain unavailable.

## Verification (ctest `mhd_*`)

| Test | Oracle | Result |
|---|---|---|
| `mhd_ampere_n_error_status` | affine/curved meshes, both mode signs, independent coupled energy, nonzero exactly representable field, perturbed candidates | energy diagnostic and explicit exterior/derivative status |
| `mhd_ampere_n_sheet_interpolation_remainder` | continuous cosine-series solution for a cylinder sheet, cubic exact recovery, mode parity, generated JVP/VJP | index 2.23--2.39; order 3.47--3.49; practical |
| `mhd_ampere_n_sheet_coupled_equilibrated` | toroidal ring sheets, perturbed candidates, missing loads, optional interpolation correction | discrete-exterior reconstruction plus explicit sheet remainder/status |
| `mhd_ampere_n_exterior_manufactured` | current `curl_n W` of a compact bump; field from toroidal-angle quadrature of `1/(4π|x−y|)` | field order k (k = 1, 2), exterior trace order ≥ 2 |
| `mhd_shielding_cylinder_limit_reciprocity` | screened response of a periodic cylinder (modified Bessel functions); sheet reciprocity | 1e-5 at R0/a = 40, second order in a/R0; M Hermitian to 1e-13 |
| `mhd_newcomb_screw_pinch_limit` | Newcomb ODE (fortnum Cash–Karp), shape of the response inside the plasma | β=0: second order in a/R0, 6e-3 at R0/a=40; β_p≈2 converges to the same limit (pressure-sensitive test needs a periodic-cylinder chart) |
| `mhd_ideal_n_mde_force_free_identity` | MDE inverse; force-free `μ0 J = σ B` against a direct FEEC solve; charge conservation with `p' ≠ 0`, `F' ≠ 0` | 5e-10; second order in the ring spacing (8.5e-5); 3e-5 |
| `mhd_ideal_n_final_state_error` | zero-current exact recovery; final frozen load rebuilt separately; force-free discrete reference on three meshes; ratio/tail derivative identities; noncontraction rejection | practical algebraic estimate only; continuous error unavailable |
| `gs_response_3d_ideal_error_contract` | production response with resonant and boundary sheets | final-state diagnostic and separate unavailable contributions reach the result file |
| `mhd_closure_tearing_delta_prime` | V5: step-profile closed form `r_sΔ' = −2m(m−1−nq0)/(m−1−nq0+(nq0/m)^m)` for the independent Newcomb-type ODE; FEEC closure Δ' against the ODE (Wesson `j=(1−r²)²`, m/n = 2/1, extrapolated in δ; tanh step); edge current jump (uniform current, response `1+1/(m−nq0−1)`); `D_I = −1/4` at `p' = 0` | ODE: O(w²) to the step limit (3.7e-3 at w = 0.0025); Δ' 1.3e-2 (R0/a = 40), 2.1e-2 (10), tanh step 2.4e-2, at nrad 32; edge 3.2e-3 at R0/a = 40, wrong by 4x without the boundary sheet |

`gs_response_3d` also recomputes the coil resonant flux through the L3
surface path (`surface_harmonics`), agreeing to 1e-5.

## Ideal closure (L6)

`kin6d_mhd_closure_n`, `closure = 'frobenius'` in `gs_response_3d`. Near a
rational `q = m/n` at `psin_s`, `x = psin − psin_s`, the resonant flux of an
ideal outer solution is `Φ_m = A|x|^(1/2−μ) + B|x|^(1/2+μ) + Cx + …` on each
side (Glasser–Greene–Johnson; `μ = sqrt(−D_I)`, Mercier index from the DCON
formula). The layer `|x| < δ` carries no volume current; two resonant sheets
on its edge rings (the shared sheet machinery of `kin6d_mhd_response_n`,
`sheet_closure_t`) are fixed by vanishing large coefficients `A±`, read from
fits of `Φ_m` on three graded rings outside the layer on each side (basis
`t^−ε, t, (t^(1+ε)−t)/ε`, `ε = μ − 1/2`, which tends to `1, t, t ln t` at zero
beta). Prescribed `A±` give Δ' data. The net layer current (sum of the edge
sheets) is the jump of `dΦ_m/dpsin` across the layer; carried to GPEC's
`sing_spot` with `(spot/δ)^(μ−1/2)` it is GPEC's `singcur`, and its self flux
per area GPEC's `singflx`. For `μ ≠ 1/2` GPEC's value depends on `sing_spot`.
Rings are graded geometrically (`layer_ratio`) from `±δ` up to the uniform
spacing or 0.45 of the distance to the next rational and smoothed
(`ring_labels_smooth`). The Ampère operator, coil drive and the MEPHIT-form
volume current are unchanged.

Verification: `mhd_closure_tearing_delta_prime` (above). With nrad 96 and
grading 1.25 the δ-extrapolated Δ' is within 0.8 % of the ODE (R0/a = 20, 40;
1.4 % at R0/a = 10), the tanh step within 5e-4. `D_I` of the TC24 cut-0.995
gfile agrees with DCON's `di0` to 3e-4 at q = 4/3 … 8/3, 1 % at q = 3 and
0.15 % at q = 10/3; `dq/dpsin` to 0.2 %.

Edge model. DCON/GPEC take the plasma energy with continuous `B0`, zero
traction and no `p_edge` term; its natural boundary condition is linearised
total-pressure balance. In the current picture this is the field-aligned
boundary sheet plus the moved current `ξ_n J0` (its perpendicular part is the
diamagnetic sheet `δp/B0` of pressure balance; the parallel sheet is fixed by
charge conservation), i.e. `boundary_sheet = .true.`,
`edge_moved_current = .true.`. This is the exact ideal response of the
pressure-offset cut equilibrium (`p_edge = 0`, same `p'`, `FF'`); the cut
gfile itself has `p_edge = 10.5 kPa` with vacuum outside and is not an
equilibrium at the cut. At zero beta the moved current reduces to the
boundary sheet, checked by the edge-current-jump oracle. On TC24 the moved
current changes `singflx` by 4–35x (nrad 96).

## Resistive closure with rotation (L5)

`kin6d_mhd_layer_n`, `closure = 'resistive'` in `gs_response_3d`. A separate
current model on the shared operator (equilibrium-core D5): the outer region
is the static ideal response above (Ampère operator, coil drive, MEPHIT-form
volume current, unchanged); each rational keeps the L6 layer and its two edge
sheets, but their amplitudes follow from matching to the resistive inner
layer of a rotating plasma instead of ideal shielding:

- continuous large coefficient across the layer, `a1(in) = a1(out) = A`
  (the reconnected flux of the outer solution);
- net layer current `J = junit (λ_in + λ_out) = Δ_layer A` (`J` the jump of
  `dΦ_m/dψ_N`; the sheet amplitudes enter the constraint directly,
  `sheet_closure_t%wl`).

Layer response. Slab resistive-inertial(-viscous) layer of the MARS-F
single-fluid model (Doppler rate `γ = −i(ω − nΩ)`, `γ = inΩ` for the static
coil field), `γ b = iF'x v + (η/μ0) b''`, `γ v'' − ν v'''' = i(F'x/μ0ρ) b''`,
scaled with `L⁴ = |γ|ρη/F'²` to two parameters `ε = Q^(−3/2)` and `P`
(viscosity) and solved numerically (finite differences on a sinh grid,
Richardson): `Δ = Δ̂(γ/|γ|, ε, P)/L` and the flux at the rational relative to
the matching value, `b(0)/b0` (1 in the constant-ψ regime, small in the
inertial one). No constant-ψ assumption. Toroidal coefficients are the
cylinder-equivalent ones: `F' = n|q'|⟨B·∇θ⟩G`, `⟨B·∇θ⟩ = 4π²ψ_o/V'`,
`G = ⟨|∇ψ_N|²⟩^½`, `Δ_ψN = Δ/G`. Not included: curvature (Glasser) terms
of the layer (TC24 has `D_I ≈ −0.28`, `μ − 1/2 ≈ 0.03`; the large coefficient
is referred to the layer edge `δ`), the toroidal inertia enhancement (option
`inertia_q2`, ρ → ρ(1 + 2q²)), MARS-F's parallel sound-wave damping and
compressibility, and rotation in the outer region (`γτ_A ≪ 1`).

Limits. η → 0 at fixed Ω does not give L6: it gives the ideal-inertial
layer, with a finite `Δ` of order `1/x_ii`, `x_ii = |γ|(μ0ρ)^½/F'`
(dissipated at the Alfvén resonances `|x| = x_ii`); the flux at the rational
`b(0) → 0`, but the layer current differs from the static ideal one by
`O(x_ii Δ')`. L6 is recovered as `x_ii → 0`.

Verification (`ctest mhd_resistive_layer_rotation`, `test_mhd_closure
resistive`):

| Oracle | Result |
|---|---|
| constant-ψ RI limit `Δ̂ = 2πΓ(3/4)/Γ(1/4) g^(5/4)/ε`, `g = 1, i, e^(−iπ/4)` | 8e-5 at ε = 1e4, order 1 in 1/ε |
| layer resolution (Q = 1), `g → ḡ` | 7e-9; conjugate to roundoff |
| viscous-resistive `Δ̂ ∝ g P^(1/6)/ε` | constant 2.107, 2.116 at P = 1e4, 1e6 (literature 2.104) |
| cylinder, Wesson `q = 2/1`, static drive, rotating plasma: matched (ODE Δ', drive jump E, layer Δ and `b(0)/b0`) against the full linear resistive reduced-MHD BVP (no matching) | `ψ(r_s)` within 2.6e-3 … 5.4e-3 for Q = 0.18 (constant-ψ), 1.8, 4.6, 11.7 (inertial); Ω → −Ω conjugates |
| FEEC resistive closure, torus R0/a = 40, against the cylinder BVP | 3.3e-3 (Q = 1.8), 6.2e-3 (Q = 11.7); outer A against the matched cylinder 4.9e-3 |
| ideal limit: η/μ0 = 8.8e-11 at fixed γ, ρ → ρ/100, ρ/10⁴, layer current against L6 | 4.4e-3, 4.4e-4, 4.4e-5 (order 1/2 in ρ, i.e. 1 in x_ii) |

Status (2026-10-13). The finite-mu (eps > 0) asymptotic matching was
falsified by its pre-declared test: the net layer current exponent and the
J_inf spread at q = 3 and 10/3 did not meet the declared criteria. The
closure is kept as the ideal/asymptotic check (the oracles above). Resistive
and finite-pressure layers go to resolved strips
([response-unified.md](response-unified.md), R-B1 on main).

## Usage

`gs_response_3d input.nml`, namelist `/l4/`: `eqdsk`, `coils`, `currents`
(one table per phasing), `n`, `rat_num`/`rat_den`, `nrad`, `nth_max`,
`order`, `model` (`'shield'` or `'ideal'`), `np_sample`, `qlow`,
`boundary_sheet`, `rat_pad`/`rat_grade` (rings graded geometrically toward
each rational), `edge_moved_current`, `closure` (`'sheet'` or `'frobenius'`),
`layer_delta`, `layer_ratio`, `layer_nfit`, `spot`, `error_estimate`
(default true: equilibrated coupled bound of every unit-sheet solve in
`prefix_sheet_error.dat`, plus the final-state ideal algebraic diagnostic in
comments in `prefix_res.dat`); with `closure =
'resistive'` also `profiles` (table ψ_N, η [Ω m], ρ [kg/m³], Ω [rad/s]),
`omega_drive`, `layer_nu`, `eta_scale`, `inertia_q2`, and the output
`prefix_resist.dat` (Δ_layer, L, Q, A, A/area, layer jump at δ and at spot,
|Φ_vac|, `b(0)/b0` and the flux at the rational per area `A b(0)/b0/area`). With MUMPS the Ampère operator is
factored as the real SPD form (1M unknowns at nrad 128). Output `prefix_res.dat`: per phasing and rational surface
the vacuum resonant flux `Φ_vac` and `b_vac = Φ_vac/area`, the sheet amplitude
`λ` and `b_shield = −M_ss λ/area`, the resonant flux of the shielding current
alone per area (GPEC's `singflx`); `m = n q` with the signed `q` of
`kin6d_vacuum_surfaces` (GPEC labels it `−m`). The file also carries a
`# ampere_discrete_energy_error_max` comment: the worst discrete residual-energy
indicator of the app's unit-sheet Ampère solves (claim `indicator`; the
continuous claim stays `unavailable-exterior-remainder`).

## Perturbation export (L7 core)

`kin6d_mhd_perturbation_n` (`src/mhd/mhd_perturbation_n.f90`), test
`mhd_perturbation_l7_xi_lagrangian_charts`. From the total Eulerian per-n field
`dB_E` (coil + plasma + sheets, coefficient of `e^{i n phi}`, physical field
`dB_E e^{i n phi} + c.c.`, cylindrical `(R, phi, Z)`, SI) on a flux surface it
reconstructs the ideal displacement `xi = xi_n nhat + xi_b ebin`
(`nhat = sgn grad psi/|grad psi|` outward, `ebin = b x nhat`, `xi_par = 0`) from
`dB.grad psi = B0.grad(xi.grad psi)` and the binormal component of
`dB = curl(xi x B0)` (both first-order ODEs along the field line, solved by
periodic Fourier collocation in the geometric angle; `div xi` drops out), and
returns `dB_L = dB_E + (xi.grad) B0`, `|B|_L = b.dB_E + xi.grad|B|`,
`kappa = (b.grad) b` and `kappa.xi` (`perturbation_t`, in memory).
`gs_response_3d` writes it for `export_psin(:)`, `export_np`, `export_mmax` to
`prefix_l7_ph<k>.dat` (spectrum quantity `phi_n` = flux-density harmonics of the normal
Eulerian field in T m^2 per rad^2, one-sided: the peak-amplitude `phi_mn` of
`surface_harmonics`, `res.dat` and GPEC/MARS is twice it) (format `kin6d-l7-v1`: header with n, n-sign, declared
COCOS 3, `sgn_dpsi`, `s_pol`, `ip` sign, chart definitions; point table; PEST and
Hamada spectra) with a `sha256sum`-compatible sidecar `.sha256`.

Limits. The singular set is `n q = m`: surfaces within `dist_rational < 1e-3`
are refused (status 2); `xi` near a rational layer is the ideal-equivalent
displacement of the sheet-closed field, not a resolved layer solution.
`xi_par` is a gauge that shifts `|B|_L` by `xi_par b.grad|B|`; it is set to
zero. `grad B0` comes from central differences of the spline equilibrium.
Boozer angles are not provided (no Boozer map in the D7 layer yet), so no
`.bc` file is written; NEO-2/NEO-RT, PENTRC and MARS (`BPLASMA`/`XPLASMA`)
conversions consume this file or the in-memory type.

The [production scaling ledger](benchmarks.md#production-estimator-scaling)
measures the same unit-sheet Ampere APIs used by `gs_response_3d`, including
both fresh setup and a reused per-n factor. The sheet reconstruction currently
requires FortFEM's global potential solve; its one-shot complex factorization
selects serial SuperLU. Timings do not establish O(N) complexity or a
small-fraction cost relative to a reused sheet solve.
