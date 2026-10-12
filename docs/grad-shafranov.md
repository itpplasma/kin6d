# Fixed-boundary Grad–Shafranov

- Scope: 2D axisymmetric fixed-boundary equilibrium in physical `(R,Z)`, `R>0`.
  Circular, shaped and supplied polygonal boundaries use the same engine.
- SI; physical counter-clockwise `phi`; flux per radian `psi`;
  `B_R=-psi_Z/R`, `B_Z=psi_R/R`, `B_phi=F/R`.
- Equation: `Delta*psi = -mu0 R² pprime - FFprime`.
  The weak form uses `grad(psi)·grad(v)/R`, FortFEM meshing/generated basis
  and FortNum preconditioned CG. `element_order=1` uses three-point quadrature;
  `element_order=2` uses continuous six-node triangles with shared edge DOFs
  and positive seven-point degree-five stiffness quadrature. Fixed-source P3
  uses ten-node triangles and positive degree-six quadrature. Profile sources
  use FortFEM integration split at the native P1/P2 flux level cuts, with
  positive retained quadrature nodes, paired error estimates and work limits.
  The rational `1/R` stiffness weight still requires quadrature refinement.
- `q_seed` sets the source strength in built-in circular models. Actual `q`
  is an output in forward mode. The [prescribed-q mode](gs-inverse.md) instead
  accepts q(s_pol), constant pprime and signed total toroidal flux, including
  curved P2/P3 meshes.

The [forward energy-majorant derivation](../theory/derivations/grad-shafranov-forward-majorant.md)
proves a continuous weighted residual bound under explicit domain, trace,
flux-conformity and stability assumptions. Native FortSym tests check its
algebraic identities and reject omitted jump, boundary and integration premises.
The numerical fixed-source estimator below is implemented independently of
these retained identities.

## Run and compare

```sh
cmake --preset cpu
cmake --build --preset cpu --target grad_shafranov_demo test_grad_shafranov
build/grad_shafranov_demo --config examples/grad_shafranov_circle.in
ctest --preset cpu -R '^grad_shafranov_'
```

- Namelist `/grad_shafranov/`: `r0,a,kappa,delta,b0,q_seed,max_area`,
  `boundary_nodes,element_order,nonlinear_tolerance,pressure_mode,mesh_file,`
  `boundary_file,profile_file`. Order defaults to 1; only 1 and 2 are supported.
- `pressure_mode='zero_beta'`: `pprime=0`, `FFprime=-2*b0/q_seed`.
  `finite_beta`: `pprime=-2*b0/(q_seed*mu0*r0²)`, `FFprime=0`.
  Boundary `psi=0`; edge pressure zero; edge `F=r0*b0`.
- Boundary file: number of points, then `R Z` rows, no repeated closing point.
- Profile file: number of rows, then ascending `psi pprime FFprime p F`.
  Sources are piecewise linear with constant derivative endpoint extension.
- Supplied profiles define a **derivative-controlled reconstruction**:
  `p=p_edge+integral(pprime dpsi)` and
  `F=sign(F_edge)*sqrt(F_edge²+2*integral(FFprime dpsi))`.
  Integrate each segment exactly; reject nonpositive `F²` or zero edge branch.
  Received `p,F` columns remain reference data, never silently overwritten.
- Nonlinear profiles use Picard relaxation 0.3, maximum 1000 iterations;
  failure returns status 2.
  Successful discrete weak residual is at most `nonlinear_tolerance`, default
  `1e-8`. The native `solve_gs_profiles` API accepts the same optional relative
  tolerance. Source quadrature uses the smaller of `1e-11` and one percent of
  the requested weak tolerance; unresolved integration returns failure.
  This algebraic residual and quadrature estimates are not observable error
  bounds or rigorous enclosures.
- Mesh export: comment, node/triangle counts, `R Z psi boundary` node rows,
  then one-based triangle connectivity. P1 retains `kin6d-gs-mesh-v1` with
  three indices. P2 emits `kin6d-gs-mesh-v2`, with vertex indices 1–3 followed
  by edge midpoint indices for edges (1,2), (2,3), (3,1). The default geometry remains
  affine; the explicit curved API below exports mesh-v3. Supplied-profile runs also export
  `mesh_file.profiles.dat` with the executed five-column law.
- `output/grad_shafranov.nc` retains literal/resolved configuration, source
  profiles, physical fields, units and provenance. Copy it per run before
  repeating; the default filename is reused.
- Physical fields, element-current samples and interior face-sheet samples
  use the same represented flux. `gs_element_jet` returns its physical value,
  gradient and Hessian inside a named element; `state%gradient` contains only
  affine geometry gradients and must not be used as a P2 flux gradient.
  Physical gradient/Hessian contraction subtracts the first local nodal flux
  before applying derivatives, using partition of unity to remove a constant
  gauge before cancellation. A constant nodal tangent likewise gives zero
  current/axis-q JVP. Binary64 rounding of nodal values after a large gauge
  shift remains an input error; the polynomial oracle bounds that amplification.
- `gs_classical_current` evaluates `j_phi` inside each element.
  `gs_current_jump` returns the distributional sheet coefficient in A/m on
  a shared face. Its unit normal points from the left element to the right.
  Sheet samples are not integrated total currents or current-error bounds.
- `gs_axis_candidate` minimises the represented flux over every mesh entity:
  cell interior stationary points, the stationary dip along each edge, and
  every node, since a coarse mesh puts the minimum exactly on a shared entity.
  An entity shared by neighbouring cells is one axis and is reported once, keyed
  by a node of that entity, so the adjacent cell cannot register the same
  physical point as a second minimum. Status 0 requires that minimum to lie
  strictly inside the domain boundary, and its flux is never above the nodal
  minimum; status 3 rejects P1, 4 says the flux stays lowest on the boundary or
  no enclosed minimum was found, and 5 rejects two distinct entities at the same
  flux. At a shared-entity minimum the reported element is the lowest-index
  incident cell and its Hessian is that cell's. A P2 element order alone does
  not ensure axis regularity. `gs_axis_safety_factor` uses the signed F at that
  axis and its positive Hessian; it is a smooth-axis candidate, without a
  certified bound.

The fixed-source solver targets a relative weak residual of `1e-10`. The
operator applies symmetric nodal differences, preserving the constant-flux
nullspace before Dirichlet elimination. This avoids cancellation of absolute
fluxes in the action. CG is limited to `max(500,10*ceil(sqrt(N)))` iterations;
failure remains explicit. Equation responses accept a requested tolerance, and
nonlinear solves use one percent of their outer tolerance (at most `1e-11`).
`grad_shafranov_shaped_solver_work` checks the nullspace, true residual and
iteration growth on two shaped meshes; it has no wall-time assertion.

## Curved P2 and P3 geometry

`initialize_gs(state, contour, max_area, status, element_order=2, curved=.true.)`
constructs six-node isoparametric triangles; `element_order=3` constructs
ten-node triangles for fixed-source and cubic-profile solves. A periodic C2 cubic spline passes
through the supplied contour samples at uniform parameters `(j-1)/n`.
Alternatively supply `curve=callback`, where `callback(t, point)` evaluates the
closed curve for `0 <= t < 1`; its values at those parameters must match the
contour. P2 moves boundary midpoints. P3 places two nodes on each boundary
edge and extends its leading quadratic curvature to the interior node. Cells
away from the boundary remain affine. Neighbors share every edge DOF. Boundary vertices are retained
without subdivision, so callers refine both the contour and the interior mesh.

The map is `X(xi)=sum_j N_j(xi)*X_j`, using FortFEM's generated basis products.
Assembly, `gs_cell_data`, physical point location and `gs_element_jet` use the
same geometry. Physical Hessians include the map-curvature connection term.
The original `solve_gs`, `gs_linear_response`, `sample_gs` and axis APIs remain
available. `read_gs_mesh` restores an evaluation-only field for `sample_gs`;
it does not assemble a new solve operator. The contour and global readback is
specified in [gs-coarea.md](gs-coarea.md). `sample_gs` selects the first containing cell on shared edges;
its gradient is a one-sided FE value there. Mesh-v3 exports six P2 geometry
nodes. Mesh-v4 exports ten P3 nodes: vertices, two nodes per directed edge
(1,2), (2,3), (3,1), then the interior node. Mesh-v1/v2 remain affine.

Fixed-geometry coefficient derivatives use the same symmetric assembled
operator and equation-level linear response. Shape derivatives remain unavailable.
The legacy piecewise-linear profile API returns status 20 on curved/P3 meshes;
the cubic-profile API below provides their fixed-geometry state derivatives.
Fixed-source integration uses positive degree-five (P2) or degree-six (P3)
triangle rules. These are ordinary floating numerical results.

The generic maps, connectivity and periodic curve remain in KIN6D. FortFEM
supplies the generated P1/P2 basis, cardinal/factor products and scalar-jet
transform. Its arbitrary-order triangle API lacks Hessians and allocates per
evaluation, so the P3 wrapper here reuses its generated cardinal products and
restricts its generated tetrahedral Hessian action to the triangle. Moving
these missing generic pieces upstream is a controller decision.

`grad_shafranov_curved_geometry` checks an independent polynomial map, inverse,
signed physical derivatives, Hessian connection and periodic spline order.
`grad_shafranov_curved_convergence` solves the exact polynomial Solovev problem
on its physical LCFS with 48, 96 and 192 boundary segments and proportional
interior refinement. It checks psi L2 order 3 and Bpol L2 order 2 within 0.3,
and the degraded affine-boundary rates. It prints DOF and total mesh/assembly/
solve wall time. Error integrals use higher-order independent quadrature over
each represented domain; they are unnormalized meridional L2 norms.

`grad_shafranov_p3_curved_convergence` checks the exact A3 Solovev LCFS on
96, 192 and 384 boundary segments: psi L2 order 4 and Bpol L2 order 3, each
within 0.3 on both refinement pairs. A3 reaches this range earlier than the
thin-torus quartic case, whose boundary-dominated B error initially drops at
order 3.5. Independent P3 tests cover polynomial values/gradients/Hessians,
curved-map inversion and connection terms, cubic contour intersections,
implicit tangent/adjoint pairing, signed q and native export/readback.
Contour sensitivities remain unsupported. Nonlinear cubic profiles use the
Newton path below; the recovered-flux estimator also supports P2 and P3.

## Phase 1 executable and producer

Build `gs_phase1` with the CPU preset. Run `gs_phase1 --config input.nml` with
`&gs_phase1` entries `boundary_file`, `boundary_nodes`, `max_area`, `pprime`,
`ffprime`, `f_edge`, `mu0`, `curved`, `degree` (2 or 3), `radial_points` and `angular_points`.
The boundary file contains a count followed by uniformly parameterized R,Z
samples, without a closing duplicate. The app interpolates that input with a
periodic cubic spline and samples the requested boundary mesh. Nonpositive
`max_area` selects half the squared mean boundary-segment length. The optional
`interior_points_file` seeds the mesh with supplied interior SI `(R,Z)` points:
its first line is the point count, followed by one pair per line. The library
accepts the same coordinates as `initialize_gs(..., interior_points=points)`.
Absent or empty input keeps the original mesh path; boundary segments and
`max_area` retain their meaning. Source
coefficients and F_edge are supplied in SI; no exact-solution values enter
the solve. Defaults are degree 2 and 96 boundary nodes. Zero radial/angular
settings resolve to `max(17,ceil(sqrt(Ncells))+1)` radial points and angular
budget `16*(degree+2)`. Explicit positive budgets remain available. The boundary
gauge is psi=0.

For cost comparisons, use the automatic readback budgets first and refine them
against the requested q accuracy. A fixed dense profile can dominate a coarse
P3 solve. Contour roots use FortNum's safeguarded Brent interpolation; the pinned
version tests convergence within bounded iteration counts, including one secant
step for a linear residual. The recovered-flux majorant and its quadrature are
independent of this root-solver optimization.

The app reports free DOF, CG iterations and total mesh/assembly/solve wall time,
and writes `mesh.dat`, `profiles.dat`, `summary.dat` and `output/gs_phase1.nc`.
The NetCDF file embeds the original boundary samples, literal/resolved input,
state, observables and standard build/runtime provenance. Readback cost is
separate from primal solve time; `summary.dat` and NetCDF record contour and
area-integration times, and the resolved namelist retains the sampling budgets.
A batch query is
`gs_phase1 --evaluate mesh.dat F_edge FFprime`; stdin supplies a point count
and R,Z pairs, stdout returns psi, BR, BZ and Bphi (NaN outside the mesh).
`--evaluate-stream` accepts repeated batches until EOF, retaining the loaded
mesh; the producer owns and closes that reader. `--q profiles.dat` accepts
a count and s_tor values and returns the native cubic q interpolation.

`equilibrium/phase1/producers/kin6d.py` implements
`run(case, resolution, workdir) -> Result` and the shared Phase 1 attributes.
It uses this native executable for both solving and physical field evaluation;
q(s_tor) calls the native cubic profile reader (`--q profiles.dat`). Set
`KIN6D_BINARY` for an installed executable. `source.kind` is `solovev`, with
SI `pprime` and `ffprime`
provided directly or under `source.coefficients`. Supported resolution keys
are `boundary_nodes`, `max_area`, `radial_points`, `angular_points`, `curved`, `degree`.
The producer forces one thread and a disk scratch directory. Its sign map is
the identity under the TC24 case contract's COCOS 3 ledger. The common driver checks the producer against the exact physical E0 and
Cerfon LCFSs, including fields loaded back from its native mesh export.

The exported axis location uses `gs_magnetic_axis`: a quartic least-squares
recovery of nodal flux near the nodal minimum, followed by Newton on the
recovered gradient. The patch radius is six local element lengths and shrinks
with refinement; compact radial weights vanish at the patch edge. QR uses
translated/scaled coordinates and gauge differences.
This location diagnostic is distinct from the represented element minimum
used internally for contour levels and axis-q. Recovery sensitivities are
unavailable. `grad_shafranov_axis_refinement` checks at least field-order
convergence against the exact Solovev axis on three independently meshed LCFSs.

## Fixed-source error estimator

`kin6d_gs_majorant::gs_functional_majorant` computes a Repin functional
majorant for `-div(grad(psi)/R)=mu0*pprime*R+FFprime/R`, with homogeneous
Dirichlet data. It recovers a continuous vector P2 flux `y` by minimizing
`D^2+C_F^2*B^2`, where `D=||sqrt(R)*(y-grad(psi)/R)||` and
`B=||source+div(y)||`. The returned estimate is `D+C_F*B` for the weighted
energy error. The positive-R Bezier hull bounds the represented domain by a
rectangle; zero extension supplies its Poincare constant, multiplied by
`sqrt(Rmax)`. All integrals and the recovery solve use ordinary floating
arithmetic. This is an a posteriori numerical estimator, not a rigorous
certificate, and does not include the physical-boundary approximation error.

The curved convergence CTest checks effectivity in `[1,3]` on every mesh
against the independent exact weighted energy error and checks the zero-source
limit. It prints the flux and balance contributions separately. Recovery uses
O(DOF) storage and a separate symmetric positive-definite solve; its cost is
not included in the reported primal mesh/assembly/solve timing.

## Cubic-profile Newton solve on curved P2/P3

`solve_gs_cubic_profiles` uses the same stiffness, geometry and field as
`solve_gs`. `gs_profile_law_t` stores continuous piecewise cubic **physical**
`pprime=dp/dpsi` and `FFprime=d(F²/2)/dpsi`, either in absolute psi or in
`s=1-psi/psi_axis`, with boundary psi=0. A normalized table does not prescribe
the flux span: the represented stationary extremum is recomputed during Newton.
Both axis minima and maxima are supported. A zero or nonunique axis, invalid
geometry, exhausted quadrature, or failed damped Newton solve returns failure.
The initial state is the constant edge-source solve (mid-profile source if both
edge sources vanish); a supplied nonzero FE state may select another branch.

The table is the executed polynomial law, with endpoint-polynomial extension.
Supplying explicit coefficients avoids changing a pinned spline through a
second fit. Cubic splines, their quadratic derivatives, and continuous linear
limits all fit the format. `gs_phase1` accepts `profile_file` and
`nonlinear_tolerance` (default `1e-9`) and the independent source
`quadrature_tolerance` (default `1e-8`). Early Newton steps use looser integration
and tighten to this target before convergence is accepted. Algebraic residual
and quadrature error are distinct numerical quantities. The text file is:

```text
nknots coordinate_mode
p_edge F_edge psi_axis_seed
knot_1 ... knot_n
pp_c0 pp_c1 pp_c2 pp_c3 ff_c0 ff_c1 ff_c2 ff_c3
... one coefficient row per interval ...
```

For a normalized derivative table only, `gs_phase1 total_current` imposes a
signed SI toroidal-current constraint. The exact scaling symmetry multiplies
psi, pprime and FFprime by the same positive factor; normalized flux is
unchanged, p-p_edge and F²-F_edge² scale by its square, and F_edge is fixed.
The input therefore prescribes derivative **shapes**, not fixed pressure values.
The current is the volume-source integral before boundary elimination, divided
by mu0. Zero means no constraint; opposite-polarity or absolute-law requests
fail. This does not implement a prescribed Istar profile or prescribed p(s).

Mode 0 means absolute psi; mode 1 means normalized s, whose endpoints must be
0 and 1. Coefficients use ascending powers of coordinate minus the interval's
left knot. Units of the **values** remain Pa/(Wb/rad) and T² m²/(Wb/rad);
normalization does not rescale these derivatives. The seed supplies the initial
profile coordinate if the edge sources vanish; it is not an additional forward
constraint. Normalized laws solve their own axis flux. Pressure and F come from exact polynomial primitives anchored at the declared edge values.
Nonpositive F² and negative pressure fail readback. `profile-law.dat` retains
the executed table and solved normalization; NetCDF embeds its full content.
The pinned TC24 law uses mode 0, and must never be rescaled to the solved axis.

Source integration uses the physical Jacobian of each mapped element. Cells
without cuts use paired Duffy rules with orders tied to the element degree.
C1 laws first try positive Duffy rules with 8, 12, 24 and 48 points per
coordinate, caching basis values and gradients for each rule across cells.
Acceptance uses separate load and Jacobian component error estimates. The
Jacobian used by Newton may have tolerance `1e-3`; the load must meet the
requested quadrature tolerance before convergence. Near convergence Newton
reuses its Jacobian for trial loads. Neither approximation changes the equation
or the derivative products returned by the default `gs_cubic_load` API.

If paired integration fails, cut P2 cells use FortFEM level integration.
For P3, the inner coordinate is partitioned
at **all cubic roots** of every relevant profile knot. Outer panels include
edge intersections and vertical tangencies (including closed interior cuts);
LAPACK roots of the source-cut discriminant locate the latter. Paired positive
rules on these strips fall back to FortFEM adaptive integration when needed.
C1 laws without possible closed cuts integrate the inner slices together;
general laws retain separate bands and topology events. Exact edge extrema
detect relevant levels outside the edge range, retaining the topology-aware
path for interior islands. C0 and C1 cut-island tests exercise regions with no
edge intersection. Bernstein hulls discard only irrelevant knots. These are
floating quadrature estimates with explicit work limits, not bounds.

`gs_cubic_load` returns a local source mass matrix and a rank-one axis term.
For normalized laws, `ds=-dpsi/axis+psi*daxis/axis²`, with `daxis` equal to the
FE value perturbation at the stationary point. No Newton iterations are
differentiated. `apply_gs_profile_jacobian` and its transpose use coefficient
and coefficient-dual pairings, with fixed geometry, fixed knots/coefficient
law and fixed Dirichlet data. `gs_profile_response` solves the tangent or
adjoint equation. A symmetric absolute-flux Jacobian first uses diagonal
preconditioned CG, with an explicit residual check. Breakdown or failure falls
back to GMRES using the stiffness solve as preconditioner; the nonsymmetric
normalized-axis problem uses that general path. Shape/axis-branch
switch derivatives and continuous derivative error bounds remain unavailable.
C0 source laws have a derivative only away from a whole cell on a genuine kink;
smooth spline laws are the intended Newton input.

`read_gs_equilibrium(..., law=law)` returns signed q and full toroidal Phi for
both axis polarities and both nonzero F branches. `profiles.dat` adds p, F,
pprime and FFprime after the original psi/s_tor/q columns. Native batch field
readback accepts a fifth `profile-law.dat` argument; `--source profile-law.dat`
reads a count and psi values and returns p/F/pprime/FFprime. The iter_tc24
producer/EQDSK extension consumes that executed law, including positive-axis
q signs and edge pressure. Exterior EQDSK continuation is still not a vacuum
solution.

### Nonlinear recovered-flux estimator assumptions

`gs_functional_majorant(..., law=law)` recovers y against the actual nonlinear
source evaluated at the candidate. Its default return is `M0=D+C_F B`, with
status **4: nonlinear stability unavailable**. M0 controls the residual in the
weighted energy dual norm on the represented domain; it is not, by itself, an
error bound for the nonlinear solution. Geometry, integration error, profile
input uncertainty and model uncertainty remain separate.

If the nonlinear residual is strongly monotone on a set containing the exact
and candidate solutions, with energy constant alpha>0, then the conditional
energy estimate is M0/alpha. Supply `stability_lower_bound=alpha` to request
that value. For an absolute-flux law with nonpositive source derivative,
alpha=1. More generally, a one-sided source Lipschitz bound L gives
alpha=1-C_F² L when positive. A normalized law also requires control of the
axis functional and its rank-one source contribution in the chosen norm;
point evaluation is not bounded on general 2D H1. Do not reuse the frozen-axis
constant as a bound for that coupled problem. The caller must justify branch
existence and the stability assumption; this API does not validate them.

The opt-in `estimate_stability=.true.` path for absolute-flux laws returns
**status 5: practical linearized estimate** when its numerical stability test
succeeds. At the converged candidate form `M_+` from the positive part of the
source derivative, so the linearized operator satisfies `J >= A-M_+` in the
energy pairing. Power iteration on `A^-1 M_+` estimates its largest eigenvalue
theta. The reported `stability_estimate` is `alpha_hat=1-theta-rho`, where rho
is the computed energy norm of the eigenpair residual. Iteration stops at
`rho <= 1e-5*max(1,abs(theta))`, with at most 60 iterations. A positive result
gives `M0/alpha_hat`; otherwise the return stays at status 4. The app requests
this path, and writes the estimate and status to both text and NetCDF.

This iteration does not prove it found the largest eigenvalue, its residual
does not enclose that eigenvalue, and the represented-space test does not
establish continuous coercivity. Source-mass and recovery quadrature errors,
nonlinear remainder and branch existence are also unbounded. Consequently
status 5 is a practical estimate, not a nonlinear bound. The recovered residual
already includes the algebraic Newton residual; adding it again would double
count it. A caller-supplied `stability_lower_bound` retains the conditional
bound path and takes precedence. Normalized-axis laws remain status 4 without
such a supplied bound. Recovery convergence and iteration count are reported
separately: an unfinished minimization can give a very conservative estimate.

The manufactured nonlinear CTests check P2 orders 3/2 and P3 orders 4/3,
bounded mesh-independent Newton steps, and monotone-law majorant effectivity
against an independent exact error. The profile oracle checks the Solovev
limit, pressure/F primitives, both polarity reversals, physical source-cut
integrals, normalized finite differences, transpose and implicit pairing.
Additional exact-error checks cover positive-source amplification with a
nonzero algebraic residual. Other checks cover rejection of an unstable
comparison, normalized-axis unavailability, and the source evaluation count for
one quadratic expressed using 513 knots (208 evaluations versus 41,808 before
the C1 path on a P3 cell).

### Dense spline knots and convergence

A C2 cubic pressure or F-squared primitive has a C1 piecewise quadratic source.
Across a regular flux surface at a knot, that source has bounded piecewise
second derivatives and no delta term in its weak second derivative. Under the
usual smooth-domain elliptic regularity and adjoint regularity assumptions,
this permits the H4 solution regularity needed for P3 energy order 3 and L2
order 4. It does not inherently limit psi convergence to order 2. Large jumps
in source curvature, small knot spacing and boundary geometry can delay that
asymptotic regime. These regularity assumptions must be checked separately
for singular boundaries or source laws with lower continuity.

`grad_shafranov_p3_spline_knot_convergence` uses a 513-knot C1 quadratic source,
an independently summed manufactured forcing and curved P3 geometry. Its
three meshes give psi orders 4.11/4.08 and Bpol orders 3.36/3.29. On the pinned
TC24 data, n=48/96/192/384 boundary-node comparisons give successive psi orders
1.94 then 3.62 and Bpol orders 1.20 then 2.49. These results support a
preasymptotic explanation; the available Bpol sequence has not yet established
order 3, so its reduced rate remains an open defect candidate.
Its 128-mode Fourier boundary and sharply varying source need substantially
more resolution than the constant-source manufactured problem. Timing,
common-domain field comparisons, quadrature sensitivity and estimator versus
Richardson evidence are retained by the iter_tc24 experiment, not a rigorous
continuous error claim.

## Profile-source residual and derivatives

- `gs_profile_load` assembles the declared piecewise-linear pressure and FFprime
  laws on a fixed represented flux, including Dirichlet elimination. Internal
  profile cuts must be integrated even when they cross a single element.
- `relative_tolerance` defaults to `1e-11`; `max_evaluations` limits each cell.
  `error_estimate` reports the norm of accumulated absolute load-error estimates
  on free rows. Failure clears the returned load. These are numerical estimates,
  with no rigorous continuous error enclosure.
- `gs_profile_load_jvp` and `gs_profile_load_vjp` act on nodal flux, mu0, profile
  values and profile-knot positions at fixed mesh geometry. Generated GS source
  products and FortNum interpolation products define their local derivatives.
  KIN6D integrates those products on the same native cut geometry, with separate
  component estimates; it differentiates the continuous weak source rather
  than the adaptive quadrature's branch decisions.
- Moving-cut boundary terms cancel because the source law is continuous at
  each knot. A whole cell lying on a genuine profile kink returns status 4
  for unavailable derivatives; there is no unique linear derivative there.
  Geometry/shape derivatives and rigorous derivative bounds remain unavailable.
- Flux and knot coordinates share a local gauge shift before profile evaluation.
  This preserves representable constant-gauge invariance without subtracting
  large reconstructed flux values after interpolation.

## Evidence and limits

### Exact cylinder with an independent axial period

The Gold–Hoyle field is an exact zero-pressure equilibrium on a straight
periodic cylinder, with `B_z=B_axis/(1+tau^2*r^2)` and `B_theta=tau*r*B_z`.
Specify the axial period `L` independently: the right-handed winding is
`q_c=2*pi*r*B_z/(L*B_theta)=2*pi/(L*tau)`. The existing cylinder evaluator's
`r0` parameter sets `L=2*pi*r0`; it is a period scale, not a claim of finite-R
toroidal geometry. Full axial and meridional fluxes satisfy
`Phi'=2*pi*r*B_z`, `Psi_mer'=L*B_theta`, and `dPhi/dPsi_mer=q_c`.

For comparison with a thin torus of major radius `Rmajor`, hold `L` and
physical pitch fixed. Right-handed `q_t=L*q_c/(2*pi*Rmajor)` therefore varies
with major radius. In the local embedding `x=R-Rmajor`, `y=-Z`,
`z=Rmajor*phi`, the native GS flux has the leading-order map
`psi=Rmajor*A_z`, `A_z'=-B_theta`; native signed q is `-q_t`.
These coordinate identities do not make the cylinder an exact finite-R GS
solution. Perturbations with axial integer n require a separate integer
toroidal-mode map; no perturbation solve is covered here.

For later KIM equilibrium-only use, its CGS scalar `u=B^2` obeys
`u'=-2*r*u/(q_c^2*R0^2+r^2)-8*pi*p'`, with `R0=L/(2*pi)` in cm.
At a positive first radius, pass the actual field's
`u0_seed=B_z(r_min)^2+B_theta(r_min)^2` in G^2. Substituting axis amplitude
for `B_z(r_min)` in the default seed is only an axis approximation.
The reader/seed path remains separately qualified; this proof test does not
execute KIM.

The optional existing pinned FortSym build registers
`grad_shafranov_periodic_cylinder_exact_identities`: seventeen exact
curl, radial-force, explicit-period, flux and leading-order map identities,
without numerical probes. Its punctured-domain algebra assumes nonzero
L/tau/major radius and a positive logarithm argument; axis regularity needs
the independent analytic series/field oracle. No production kernel changes
or new runtime CAS dependency are introduced. The one-dimensional force
balance and winding definitions follow J. P. Freidberg,
[*Ideal MHD* (2014), Chapter 5, pp. 85–122](https://doi.org/10.1017/CBO9780511795046.006),
and [MIT 22.615 Lecture 5, pp. 1, 4, 7](https://ocw.mit.edu/courses/22-615-mhd-theory-of-fusion-systems-spring-2007/a5c962e43b645b825b13715a131166e2_lecture5.pdf).
Exact book equation numbers are not asserted.

- [Derivation](../theory/derivations/grad-shafranov.md): signs, weak form,
  Solovev/logarithmic families, Gold–Hoyle cylinder and flux/safety-factor map.
- `gold_hoyle_demo --config cylinder.in` exports the independent straight
  cylinder to `output/gold_hoyle.nc`; namelist `/gold_hoyle/` accepts
  `r0,a,b0,q,samples`. Generated C0 kernel and physical force-balance,
  constant-q and flux-derivative oracles keep it separate from finite-R GS.
  Example configuration: `examples/gold_hoyle.in`.
- `test_grad_shafranov`: independent polynomial refinement; nonlinear
  `psi=R²` cancellation oracle; linear implicit tangent/adjoint pairing;
  local generated field JVP/VJP and analytical profile primitive.
- `test_gs_p2`: independent quadratic value/gradient/Hessian and signed field,
  classical/sheet-current, direct elliptic coarea axis-q, boundary-only weak
  solves, nonlinear cancellation, equation derivative and cusp-rejection oracles.
  `test_gs_p2_refinement` checks off-grid physical field and classical-current
  refinement against the exact quartic GS solution, including P1 at the same mesh.
- The P1 reference diagnostic accepts only mesh-v1 and rejects P2.
  Known-solution helpers and their generated leaves live under `test/`,
  outside the production library.
  Its cellwise certificate is not valid for a quadratic representation.
- Optional pinned FortSym build registers `grad_shafranov_symbolic_identities`:
  eight exact identities, no numerical probes used as proof.
- `gen/gen_grad_shafranov.f90` generates source/field primal, JVP and
  VJP leaves. Fixed-geometry constant-source linear equation responses are
  supported. Cubic-profile state tangents/adjoints use the equation Jacobian;
  branch switching and boundary-shape derivatives are unavailable.
- Every export is a **candidate**. Polygonal boundary, P1/P2 field jumps,
  discretization and profile-model errors need independent refinement.
  A small residual does not certify the continuous equilibrium.
- Free boundary is a library path only (`src/gs_free_boundary.f90`: FEM-BEM
  Newton, limiter / X-point / target-mean flux levels, cut source with the
  interface Jacobian; no application, mesh refit or inverse problem yet).
  No 3D MHD, kinetic transport or perturbation model.
  Flux charts and GEQDSK conversion are downstream adapters, not GS inputs.
- Current ITER comparison inputs/results and disagreement classification live
  in the TC24 repository's `equilibrium/`; received profiles are not exact
  oracles. Do not use cross-code agreement alone as validation.
