# Fixed-boundary Grad–Shafranov

- Scope: 2D axisymmetric forward equilibrium in physical `(R,Z)`, `R>0`.
  Circular, shaped and supplied polygonal boundaries use the same engine.
- SI; physical counter-clockwise `phi`; flux per radian `psi`;
  `B_R=-psi_Z/R`, `B_Z=psi_R/R`, `B_phi=F/R`.
- Equation: `Delta*psi = -mu0 R² pprime - FFprime`.
  The weak form uses `grad(psi)·grad(v)/R`, FortFEM meshing/generated basis
  and FortNum preconditioned CG. `element_order=1` uses three-point quadrature;
  `element_order=2` uses continuous six-node triangles with shared edge DOFs
  and positive seven-point degree-five stiffness quadrature. Profile sources
  use FortFEM integration split at the native P1/P2 flux level cuts, with
  positive retained quadrature nodes, paired error estimates and work limits.
  The rational `1/R` stiffness weight still requires quadrature refinement.
- `q_seed` sets the source strength in built-in circular models. Actual `q`
  is an output; constant-q inversion is unavailable.

The [forward energy-majorant derivation](../theory/derivations/grad-shafranov-forward-majorant.md)
proves a continuous weighted residual bound under explicit domain, trace,
flux-conformity and stability assumptions. Native FortSym tests check its
algebraic identities and reject omitted jump, boundary and integration premises.
An unknown-equilibrium estimator API is not implemented by these proofs.

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
- `gs_axis_candidate` searches quadratic stationary minima and all edge
  extrema. Status 0 requires a strict interior global candidate minimum;
  status 3 rejects P1, 4 leaves vertex/edge minima unqualified, and 5 rejects
  nonunique interior minima. A P2 element order alone does not ensure axis
  regularity. `gs_axis_safety_factor` uses the signed F at that axis and its
  positive Hessian; it is a smooth-axis candidate, without a certified bound.

## Curved P2 geometry

`initialize_gs(state, contour, max_area, status, element_order=2, curved=.true.)`
constructs six-node isoparametric triangles. A periodic C2 cubic spline passes
through the supplied contour samples at uniform parameters `(j-1)/n`.
Alternatively supply `curve=callback`, where `callback(t, point)` evaluates the
closed curve for `0 <= t < 1`; its values at those parameters must match the
contour. Only boundary edge midpoints move. Interior elements remain affine,
and neighboring elements share every edge DOF. Boundary vertices are retained
without subdivision, so callers refine both the contour and the interior mesh.

The map is `X(xi)=sum_j N_j(xi)*X_j`, using FortFEM's generated P2 basis.
Assembly, `gs_cell_data`, physical point location and `gs_element_jet` use the
same geometry. Physical Hessians include the map-curvature connection term.
The original `solve_gs`, `gs_linear_response`, `sample_gs` and axis APIs remain
available. `sample_gs` selects the first containing cell on shared edges;
its gradient is a one-sided FE value there. Mesh-v3 exports the six physical
geometry nodes explicitly. Mesh-v1/v2 retain their affine interpretation.

Fixed-geometry coefficient derivatives use the same symmetric assembled
operator and equation-level linear response. Shape derivatives and curved
nonlinear profile-load derivatives are unavailable; their APIs return status
20. Integration uses the positive seven-point triangle rule. These are ordinary
floating numerical results, without interval or rigorous certification.

The small generic map and periodic curve implementation currently live in
KIN6D because the pinned FortFEM supplies the P2 basis and scalar-jet transform
but not the six-node map. Moving this layer upstream is a controller decision.
P3 geometry is not implemented.

`grad_shafranov_curved_geometry` checks an independent polynomial map, inverse,
signed physical derivatives, Hessian connection and periodic spline order.
`grad_shafranov_curved_convergence` solves the exact polynomial Solovev problem
on its physical LCFS with 48, 96 and 192 boundary segments and proportional
interior refinement. It checks psi L2 order 3 and Bpol L2 order 2 within 0.3,
and the degraded affine-boundary rates. It prints DOF and total mesh/assembly/
solve wall time. Error integrals use higher-order independent quadrature over
each represented domain; they are unnormalized meridional L2 norms.

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
- The P1 reference enclosure executable accepts only mesh-v1 and rejects P2.
  Its cellwise certificate is not valid for a quadratic representation.
- Optional pinned FortSym build registers `grad_shafranov_symbolic_identities`:
  eight exact identities, no numerical probes used as proof.
- `gen/gen_grad_shafranov.f90` generates source/field/Solovev primal, JVP and
  VJP leaves. Fixed-geometry constant-source linear equation responses are
  supported; nonlinear branch and boundary-shape derivatives are unavailable.
- Every export is a **candidate**. Polygonal boundary, P1/P2 field jumps,
  discretization and profile-model errors need independent refinement.
  A small residual does not certify the continuous equilibrium.
- No free-boundary solver, 3D MHD, kinetic transport or perturbation model.
  Flux charts and GEQDSK conversion are downstream adapters, not GS inputs.
- Current ITER comparison inputs/results and disagreement classification live
  in the TC24 repository's `equilibrium/`; received profiles are not exact
  oracles. Do not use cross-code agreement alone as validation.
