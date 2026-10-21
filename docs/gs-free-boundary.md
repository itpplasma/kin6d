# Free-boundary Grad-Shafranov: refit loop and `gs_free`

Status: P8 (outer refit loop) and P11 (application) of the free-boundary design.
The formulation (FEM-BEM coupling on a fixed curve Gamma_c, cut-cell source with
exact interface Jacobian, Newton on the coupled operator) is in the P6/P7 modules
`kin6d_gs_coupling` and `kin6d_gs_free_boundary`; the coil-current inverse is
`kin6d_gs_coil_inverse` (P9). This note covers what P8 and P11 add.

## Refit loop (`gs_fb_refit_solve`, `src/equilibrium/gs_free_boundary.f90`)

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

## Error estimates (`src/equilibrium/gs_majorant.f90`, DESIGN section 4, P10)

Energy norm `|||v|||^2 = int |grad v|^2/R dA + <S v,v>`. For the frozen
problem `a(psi,v)+<S psi,v>=(f,v)+<g,v>` and every conforming H(div) flux sigma,

    |||psi-psi_h||| <= D + B + N,
    D = ||sqrt(R)(sigma-grad psi_h/R)||,
    B = source_oscillation + C_F ||div sigma+Pi_p f||,
    N = ||sigma.n+S psi_h-g||_{S^-1}.

The continuous residual identity is checked by the retained FortSym identities.
Recovery now uses FortFEM's mapped RT vertex patches. Boundary flux moments
remain free reconstruction unknowns; their Piola-transformed traces enter the
same reference-BEM mismatch N. This preserves arbitrary boundary loads without
asserting homogeneous Dirichlet data. The global vector recovery matrix and
factorization have been removed. Affine cells use local mean-free source
bounds; curved cells use the coupled Friedrichs constant.

The estimate remains practical: S and its inverse norm come from a refined
reference BEM, the coupled Friedrichs constant is estimated, and integrals use
ordinary quadrature. Patch recovery does not validate those remaining premises.
In addition, N tests the residual only against the finite reference trace
space. A nonzero full-trace residual can annihilate that space, so increasing
the reference panel count supplies no upper bound on the omitted part.
`gs_coupled_majorant_t%exterior_status` explicitly retains
`unavailable-bem-remainder`, independently of the nonlinear stability status.

FortFEM's `bem_residual_factor` and `bem_residual_energy` quantify flux
enrichment on a common curve: for single-layer matrix V_h and residual
`r_h=b_h-V_h lambda_H`, the diagnostic is `sqrt(r_h^T V_h^-1 r_h)`.
It equals the V_h energy of the enriched Galerkin flux correction. With exact
continuous operators it is a lower bound on the omitted flux error; it is not
added to M0 as a purported continuous upper bound. The fast
`gs_bem_single_layer_reference_residual` test checks the correction identity
on the GS boundary operator and exact recovery in the reference space. The
independent FortFEM Fourier test exhibits a zero restricted residual with
nonzero error. Neither test supplies the missing continuous BEM remainder.

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

A guaranteed continuous lower stability bound is still unavailable for this
normalized free-boundary formulation. Its Jacobian contains the point-value
functionals `g_a` and `g_b` (`gs_free_boundary.f90`, module contract). Point
evaluation is unbounded in the two-dimensional H1 energy space: a logarithmic
spike, constant inside radius epsilon and proportional to
`log(1/r)/sqrt(log(1/epsilon))` outside, has bounded gradient norm but an
unbounded center value. Thus an enclosed eigenvalue of the discrete Jacobian
alone would not establish the required continuous operator bound. A stronger
function-space/regularity argument or bounded flux-level functionals, a
continuous FEM-BEM remainder and a nonlinear neighborhood theorem are needed
before replacing beta or asserting nonlinear product goal bands. The verified
fixed-Dirichlet absolute-flux coercivity path does not apply to this model.

The short `gs_free_boundary_equilibrated_quick` gate uses the retained
exact coil-plus-bump solution on curved meshes (two threads):

| Case | Boundary edges | Energy error | Majorant |
|---|---:|---:|---:|
| Frozen P2 | 12 | 0.9031 | 4.679 |
| Frozen P2 | 24 | 0.2783 | 0.5208 |
| Frozen P3 | 12 | 0.4724 | 2.788 |
| Nonlinear P2, M0/beta | 24 | 0.2783 | 0.5444 |

It also checks the shape/current goals and a perturbed boundary-load ledger.
At the nonlinear rung beta is about 0.95557 and agrees with the independent
dense SVD. These small cases establish their stated containment/refinement
checks, without a general efficiency or rigorous-certification claim.

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

### Parallel free-boundary path (K8)

The P13 stall had three causes, all fixed on `gs-free-parallel`:

- The Gamma_c Steklov operator was reassembled by every `gs_fb_setup`. It is
  now cached per curve (`gs_coupling_steklov`, keyed on the bitwise BEM trace nodes;
  a refit that enters the same curve at another panel is relabelled, not
  reassembled) and shared with the majorant estimator.
- `kin6d_spatial`, which holds the upstream BEM assembly, was compiled without
  OpenMP, so its panel-pair colouring ran serially. It now links OpenMP.
- The indefinite J0 = A - M (symmetric, checked at run time) uses native
  FortSparse pivoted LDL^T with AMD ordering. An asymmetric J0 uses native
  general LU. `KIN6D_FB_SOLVER=lu` selects general LU for diagnostic comparisons;
  the old value `superlu` remains an alias for this native path.
  `KIN6D_FB_ORDERING=nd` selects nested dissection. Circuit sensitivities use
  one retained block solve, including on the LU path: a factor owns mutable
  solve/refinement workspaces and cannot be shared by concurrent calls.

The majorant estimator factors its reference energy operator with native
Cholesky, uses independent local RT recovery solves and threads the Friedrichs
mass assembly. The timing baselines below predate the RT recovery replacement. The axis and X-point cell
scans and the CSC pattern build are threaded. `KIN6D_PROFILE=1` prints per-phase
wall times (`PROF label seconds`) on stderr.

`gs_free` wall time on the same fixture and node, seconds, before / after
(1 / 4 / 16 / 32 threads, sCluster exclusive node):

| degree | boundary nodes | before | after | speedup at 32 | 1 to 32 threads |
|-------:|---------------:|--------|-------|--------------:|----------------:|
| 2 | 96 | 5.0 / 3.1 / 2.7 / 2.6 | 3.43 / 1.04 / 0.45 / 0.38 | 6.8x | 8.9x (was 1.9x) |
| 2 | 192 | 14.9 / 10.5 / 9.4 / 9.3 | 9.00 / 2.80 / 1.37 / 1.22 | 7.6x | 7.4x (was 1.6x) |
| 2 | 384 | 51.2 / 38.9 / 35.6 / 35.6 | 26.84 / 8.43 / 4.02 / 3.64 | 9.8x | 7.4x (was 1.4x) |
| 3 | 96 | 14.7 / 5.8 / 3.6 / 3.3 | 13.06 / 3.66 / 1.35 / 1.02 | 3.3x | 12.8x (was 4.4x) |
| 3 | 192 | 50.4 / 23.4 / 17.0 / 16.0 | 41.02 / 11.69 / 4.24 / 3.14 | 5.1x | 13.1x (was 3.1x) |
| 3 | 384 | 305 / 199 / 172 / 168 | 164.9 / 49.3 / 20.2 / 15.6 | 10.7x | 10.6x (was 1.8x) |

Results agree across all thread counts and with the previous code to 2e-6
relative (psi_axis, psi_bound, plasma current). Part of that spread is old: the
refit loop flips between two nearly equal end states with roundoff for P2. The
remaining time at 32 threads is the MUMPS factorisation (about 30 %) and the
per-element source quadrature of the cut cells.

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
| output | `mesh_file`, `profile_points`, `error_estimate` (default true) |

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
- error estimate of the final psi (`gs_fb_majorant` through `gs_fb_certifier_t`):
  `majorant_status` (0 available, -1 not requested, otherwise the lane status),
  `majorant` (M0/beta), `majorant_frozen` (M0), `majorant_stability` (beta),
  `majorant_flux` (D), `majorant_balance` (B), `majorant_interface` (N),
  `estimator_s`, with the attributes `majorant_claim` (`practical estimate`, or
  `unavailable`/`unavailable-not-requested`), `majorant_method`,
  `majorant_quadrature_status` (`approximate`) and `majorant_scope`. The claim
  stays practical for the reasons listed under Error estimates; `wall_s`
  excludes the estimator;
- interface scope/status: `majorant_interface_scope` states that N is a dual
  norm restricted to the reference-BEM trace space; `majorant_exterior_status`
  is `unavailable-bem-remainder` for an available practical estimate and the
  estimator's unavailable status otherwise. No finite continuous exterior
  radius is supplied;
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

The [production scaling ledger](benchmarks.md#production-estimator-scaling)
records the coupled/refitted `gs_free` fixture at three boundary resolutions.
`estimator_s` includes recovery, reference BEM and practical stability work;
`wall_s` includes mesh setup, the nonlinear solve, refits and readback.
Consequently their ratio describes the application workload, not one Newton
factorization. Refinement may change the final triangulation between threads.
