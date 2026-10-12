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

## Application `gs_free`

`gs_free --config input.nml` reads `&gs_free`; `gs_free --evaluate` reuses a result.
Units: m, A, Wb/rad. Flux is absolute with psi = 0 on the symmetry axis.

| group | variables |
|---|---|
| mesh | `boundary_file` (count, then R Z rows: closed curve Gamma_c, as `gs_phase1`), `boundary_nodes`, `max_area`, `degree` (2, 3), `curved`, `plasma_ellipse` = R0, a, b (initial Gamma_p chain), `chain_nodes` |
| coils | `coil_file` (`gs_coils` text table), `currents(:)` (circuit currents, table order) |
| profiles | `profile_file` (normalised law file, as `gs_phase1`), `mu0` |
| LCFS | `lcfs_mode` = `xpoint`, `limiter`, `target_mean`; `cut_fraction`; `axis_guess`, `xpoint_guess`, `search_radius`; `points_file` (limiter or target points); `axis_sign` |
| solve | `max_refits`, `refit_tolerance`, `newton_tolerance`, `quadrature_tolerance`, `max_newton_iterations`, `seed_amplitude`, `seed_file` (initial flux map: `nr nz`, R, Z, psi with R fastest; replaces the vacuum-plus-Gaussian seed) |
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
- flux surfaces (`read_gs_equilibrium` with `psi_bound` on the plasma cells):
  `q_status` (0 when written), `q_psi`, `q_s` (normalised poloidal flux, axis 0,
  boundary 1), `q`, `s_tor` (normalised toroidal flux), `toroidal_flux` (Wb);
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

Flux surfaces. `read_gs_equilibrium` (`kin6d_gs_readback`) takes the optional
`psi_bound` and `plasma_cell`: it reads back on the sub-mesh of the plasma cells
with the flux measured from psi_b, so the same contour code gives q and the
toroidal flux (contour integral of F |dpsi| / (R |grad psi|)) up to the cut level.
The whole-mesh volume and flux sums of the fixed-boundary path are not used there;
the volume is 2 pi times the plasma current of a unit-dp/dpsi law on the same cut
cells. `test_gs_readback_free` checks the free path against the fixed path (gauge
shift, and an interior surface taken as the boundary) and the exact Solov'ev
surfaces. The coil table, profile law and constraints for ITER
(P12) stay in the iter_tc24 producer.

Tests: `gs_free_app_contract` (CTest; forward solve with refit, NetCDF values and
provenance, `--evaluate` against the axis flux and stream parity, the three fit
variants on twin data, rejected input) and `gs_free_boundary_p8_refit_loop`.
