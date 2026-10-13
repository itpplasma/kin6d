# Free-boundary Grad-Shafranov: refit loop and `gs_free`

Status: P8 (outer refit loop) and P11 (application) of the free-boundary design.
The formulation (FEM-BEM coupling on a fixed curve Gamma_c, cut-cell source with
exact interface Jacobian, Newton on the coupled operator) is in the P6/P7 modules
`kin6d_gs_coupling` and `kin6d_gs_free_boundary`; the coil-current inverse is
`kin6d_gs_coil_inverse` (P9). This note covers what P8 and P11 add.

## Refit loop (`gs_fb_refit_solve`, `src/gs_free_boundary.f90`)

1. Newton on the given mesh (`gs_fb_newton`).
2. Distance of the Gamma_p chain from the computed boundary,
   `max_j |psi - psi_b| / |grad psi|` at the chain points, in units of the chain
   spacing h (`gs_fb_chain_movement`). Stop when it is below `tolerance` (1e-3).
3. Otherwise `gs_fb_refit_mesh`: trace {psi = psi_b} around the axis along rays,
   resample at equal arc length, project onto the curve (Newton along grad psi),
   rebuild the two-loop mesh (same Gamma_c, order and curvature, the chain as
   curved P2/P3 edges), interpolate psi onto the new nodes, re-solve.
4. The chain cannot follow the boundary closer than the flux discretisation error.
   If a refit does not reduce the distance below `contraction` (0.7) times the
   previous one, the loop ends and reports `stagnated`. Exhausting `max_refits`
   while still contracting returns `GS_FB_REFIT_NOT_CONVERGED` with a valid solution.

The boundary is traced along rays from the axis, so the plasma must be star-shaped
about the axis. A mesh-dependent psi-independent load (`sys%extra_load`, used by
manufactured tests) is rebuilt through the `extra_cb` callback after each mesh change.

### Faces on the boundary and slivers

On a fitted mesh the nodes of a chain edge have psi = psi_b up to rounding. In
`gs_fb_evaluate`, a face of a plasma cell whose nodes all satisfy
`|psi - psi_b| <= face_tolerance |psi_a - psi_b|` (default 1e-10) is snapped to
exactly psi_b in the level-cut input (`cell_psi` of `gs_cubic_load`). This removes
the rounding-level slivers on both sides. The cell is then assigned by the sign of
psi - psi_b at its centroid: on the plasma side it keeps the face as its boundary
and the Jacobian receives the face interface term
`sigma int_face f_edge v_i v_j / |grad psi| dl` (gradient of the plasma cell); on
the vacuum side the cell leaves the plasma. A cell that is zero at its centroid
too counts as vacuum, so the old status 4 ("cell constant on a kink") cannot occur
for it. Faces beyond the tolerance take the ordinary cut path. `test_gs_fb_refit`
checks: the source is continuous in the chain-node shift (1e-13, 1e-9, 1e-7), and
the source derivative equals a one-sided finite difference to 4e-4 (P2) and 1e-5
(P3), against 0.6 and 0.9 without the face term.

## Order with and without refit

`test_gs_fb_refit study`: strong-jump source (edge value 70 % of the axis value) plus
a smooth load, limiter mode, coarse-to-fine against a refit P3 reference with 144
boundary edges. Errors are nodal flux at six probes, psi_b, plasma current and the
outboard boundary radius. Finest-pair orders (72 vs 36 boundary edges, one halving):

| mesh  | psi probes | psi_b | Ip   | R_lcfs |
|-------|-----------:|------:|-----:|-------:|
| P2 fixed | 2.2 | 0.9 | 0.8 | -2.0 |
| P2 refit | 1.3 | 2.0 | 3.8 | 2.8 |
| P3 fixed | 4.1 | -1.0 | 2.2 | 0.8 |
| P3 refit | 3.9 | 3.2 | 3.0 | 4.0 |

The fixed mesh stalls on the functionals that depend on the cut (psi_b, Ip, boundary
position); the refit recovers the order. The probe flux is smooth away from the
boundary and keeps its order on both meshes (the P2 probe pair is limited by a
probe near the boundary). The manufactured case T3 (`test_gs_fb_manufactured`)
cannot discriminate: its nodal load subtracts the source of the interpolated exact
flux, which cancels the jump on any mesh. With the refit it gives nodal orders 3.7
(P2) and 3.8 (P3) against 2.7 and 3.3 fixed (`test_gs_fb_refit sweep`).

## Error estimates (`src/gs_majorant.f90`, DESIGN section 4, P10)

Energy norm |||v|||^2 = int |grad v|^2/R dA + <S v, v>. For the frozen problem
a(psi,v) + <S psi,v> = (f,v) + <g,v> and every continuous vector P2/P3 flux y,

    |||psi - psi_h||| <= D + B_T + N,
    D   = ||sqrt(R)(y - grad psi_h/R)||,
    B_T = C_F ||Pi0 r|| + (sum_T (h_T/pi)^2 R_max,T ||r - Pi0 r||_T^2)^(1/2),
    N   = ||y.n + S psi_h - g||_{S^-1},   r = f + div y,

from the residual identity a(e,v) + <S e,v> = (y - grad psi_h/R, grad v) +
(f + div y, v) + <g - S psi_h - y.n, v>, checked exactly by FortSym
(`test_gs_coupled_majorant_identities`, with omitted-defect and flipped-normal
controls). B_T splits the balance residual into cell means (coupled Friedrichs
constant C_F) and mean-free parts (Payne-Weinberger h_T/pi); it is at most C_F
||r||, the form in DESIGN, and makes B_T of higher order than D. y minimises
D^2 + sum of the B_T squares + N^2 (PCG preconditioned with the sparse LU of the
volume part). The estimate is practical, not rigorous: S and the dual norm of N
come from a reference BEM with three P2 panels per mesh edge
(`gs_coupled_estimator_t`), and C_F^2 is the power-iteration estimate of the
largest eigenvalue of (K + S_ref)^-1 M, approached from below.

Nonlinear free boundary (`gs_fb_majorant`): M0 is the frozen-source majorant at
psi_h (law in the plasma cells where s < 1, plus a pointwise forcing that must
match `sys%extra_load`); beta is the smallest singular value of J in the energy
pairing (`gs_fb_inf_sup`, power iteration on J^-1 A J^-T A, which approaches
1/beta^2 from below); the bound |||e||| <~ M0/beta is first order. Goals
(`gs_coupled_goal`): Q(psi) = int q psi + <q_b, psi> lies in Q(psi_h) + rho +-
M_psi M_z for the frozen problem, M_z the adjoint majorant; the Ip goal is
Ampere on Gamma_c (`gs_coupled_ip_moments`, z = 1/mu0 is exact). For the
nonlinear solution `gs_fb_goal` returns |||z_h||| M0, with J^T z_h = Q' (first
order). `gs_coil_current_bounds` propagates row bounds b_i to |dI_k| <=
sum_i |(H^-1 D^T)_ki| b_i on the free set (attained by the sign vector).

`test_gs_coupled_majorant` (exact psi* = coils + C^5 bump, OMP 4):

| case | edges | energy error | majorant | I_eff |
|------|------:|--------:|---------:|------:|
| frozen P2 | 24 / 48 / 96 | 2.78e-1 / 9.09e-2 / 2.24e-2 | 5.64e-1 / 1.48e-1 / 2.93e-2 | 2.03 / 1.63 / 1.31 |
| frozen P3 | 24 / 48 | 7.38e-2 / 1.53e-2 | 1.79e-1 / 2.44e-2 | 2.42 / 1.60 |
| nonlinear P2, M0/beta | 24 / 48 / 96 | 2.78e-1 / 9.09e-2 / 2.24e-2 | 5.89e-1 / 1.57e-1 / 3.16e-2 | 2.12 / 1.73 / 1.41 |

The efficiency index decreases towards 1 (majorant order 2.34 against 2.02 for
P2, 2.87 against 2.27 for P3); beta (0.956) equals the dense SVD to 1e-5; the
shape and Ip enclosures contain the exact values.

## Large aspect ratio (T5, `test_gs_fb_large_aspect`, P13)

Circular plasma, a = 0.5 m, R0 = a/eps, uniform vertical field from a ring pair at
(3 R0, +-1.5 R0); the Shafranov field mu0 Ip/(4 pi R_g)(ln(8 R_g/a_g) + beta_p +
l_i/2 - 3/2), from the measured boundary centre, half width, beta_p and l_i,
against the coil field at (R_g, 0):

| eps | 1/4 | 1/8 | 1/16 | 1/32 |
|-----|----:|----:|-----:|-----:|
| relative mismatch | -2.48e-2 | -6.57e-3 | -1.71e-3 | -4.08e-4 |

observed order 1.9 to 2.1. A quadrupole with current of the sign of Ip raises the
elongation (1.024 to 1.092 at eps = 1/8), the opposite sign lowers it (0.962).

## Wall time of `gs_free` (P13)

`gs_free` forward solve with up to three refits and readback (its `wall_s`) on the
app fixture, limiter mode, sCluster exclusive node (AMD EPYC 9455), conda gfortran,
seconds:

| degree | boundary nodes | nodes | 1 thread | 4 | 16 | 32 |
|-------:|---------------:|------:|---------:|--:|---:|---:|
| 2 | 48 | 1 401 | 2.27 | 1.28 | 1.05 | 1.08 |
| 2 | 96 | 5 497 | 5.02 | 3.15 | 2.67 | 2.64 |
| 2 | 192 | 21 797 | 15.1 | 10.5 | 9.33 | 9.32 |
| 2 | 384 | 86 693 | 52.1 | 38.7 | 36.0 | 35.3 |
| 3 | 48 | 3 115 | 6.26 | 2.49 | 1.63 | 1.74 |
| 3 | 96 | 12 295 | 14.6 | 5.80 | 3.63 | 3.31 |
| 3 | 192 | 48 898 | 51.9 | 23.7 | 17.4 | 15.8 |
| 3 | 384 | 194 770 | 311 | 202 | 176 | 167 |

Strong scaling saturates near 1.5x (P2) and 2-4x (P3). A profile (P2, 192 boundary
nodes) puts about 40 % in the Gamma_c BEM assembly, which `gs_fb_refit_solve`
repeats in `gs_fb_setup` on every refit although Gamma_c does not change, and a
few per cent each in the serial SuperLU factorisation and the cut-cell source.
The coil flux at points (`gs_circuit_flux_vectors`, about 19 %) is threaded over
points (32-thread time at 384 nodes: 46.8 s before, 35.3 s after for P2; 189 s
and 167 s for P3).

## Application `gs_free`

`gs_free --config input.nml` reads `&gs_free`; `gs_free --evaluate` reuses a result.
Units: m, A, Wb/rad. Flux is absolute with psi = 0 on the symmetry axis.

| group | variables |
|---|---|
| mesh | `boundary_file` (count, then R Z rows: closed curve Gamma_c, as `gs_phase1`), `boundary_nodes`, `max_area`, `degree` (2, 3), `curved`, `plasma_ellipse` = R0, a, b (initial Gamma_p chain), `chain_nodes` |
| coils | `coil_file` (`gs_coils` text table), `currents(:)` (circuit currents, table order) |
| profiles | `profile_file` (normalised law file, as `gs_phase1`), `mu0` |
| LCFS | `lcfs_mode` = `xpoint`, `limiter`, `target_mean`; `cut_fraction`; `axis_guess`, `xpoint_guess`, `search_radius`; `points_file` (limiter or target points); `axis_sign` |
| solve | `max_refits`, `refit_tolerance`, `newton_tolerance`, `quadrature_tolerance`, `max_newton_iterations`, `seed_amplitude` |
| inverse | `fit` = `none`, `unconstrained`, `bounded`, `fixed_at_limit`; `constraints_file`; `lambda`; `reference_currents`, `n_reference`; `tie_a`, `tie_b`; `emulate_freegs`, `relaxation`, `max_fit_iterations`, `fit_step_tolerance` |
| output | `mesh_file`, `profile_points` |

Constraint rows, one per line, `#` comments: `KIND R Z R2 Z2 TARGET SCALE`, with KIND
`isoflux_pair`, `xpoint_br`, `xpoint_bz`, `psi_point`, `isoflux_lcfs`, `psi_bound`
(definitions in `gs_coil_inverse`; the circuit limits come from the coil table).
With `fit /= none` each pass solves, fits the currents, and refits the mesh onto the
new boundary until the chain is within `refit_tolerance` h (the model is rebuilt on
every mesh). `tie_a(k) = tie_b(k)` ties two circuit currents to be equal.

Outputs, in the working directory: `output/gs_free.nc` (NetCDF with the provenance
attributes of `kin6d_result_io`: command line, literal and resolved namelist,
commit, build, plus `lcfs_mode`, `fit`, coil and profile file names, thread count);
`mesh.dat` (native mesh, flux measured from psi_b), `profile-law.dat` (executed
law, axis = psi_a - psi_b), `summary.dat`. NetCDF variables:

- equilibrium: `points`, `triangles`, `psi`, `psi_axis`, `psi_bound`, `axis`,
  `xpoint`, `plasma_current`, `volume`, `lcfs` (256 points), `chain`,
  `boundary_faces`, `refits`, `refit_converged`, `refit_movement`;
- profile readback of the executed law on the plasma flux range: `profile_psi`,
  `profile_s`, `profile_pressure`, `profile_f`, `profile_pprime`,
  `profile_ffprime`;
- coils: `circuit_currents`, `coil_currents`, `circuit_i_min`, `circuit_i_max`;
- inverse (when `fit /= none`): `fit_status`, `fit_iterations`, `fit_phi`,
  `fit_phi_data`, `fit_phi_reg`, `fit_residual`, `fit_gradient`, `fit_pinned`, and
  the GSVD report `gsvd_values` (ascending), `gsvd_modes`, `gsvd_factor_x`,
  `gsvd_null_mode`, `gsvd_condition`.

`gs_free --evaluate mesh.dat F_edge FFprime [profile-law.dat [psi_bound]]` reads
points from stdin (count, then R Z rows) and writes psi, B_R, B_Z, B_T per point,
as `gs_phase1 --evaluate`; psi is absolute when `psi_bound` is given. B_T = F / R
with F from the law inside the plasma and `F_edge` outside it. Points outside
Gamma_c give NaN. `--evaluate-stream` serves repeated requests until EOF.

Limits. The flux-surface readback of `gs_phase1` (q, toroidal flux, volume by the
whole-mesh sums) assumes psi = 0 on the mesh boundary and is not used; the volume
is 2 pi times the plasma current of a unit-dp/dpsi law on the same cut cells, and no
q profile is written yet. The coil table, profile law and constraints for ITER
(P12) stay in the iter_tc24 producer.

Tests: `gs_free_app_contract` (CTest; forward solve with refit, NetCDF values and
provenance, `--evaluate` against the axis flux and stream parity, the three fit
variants on twin data, rejected input) and `gs_free_boundary_p8_refit_loop`.
