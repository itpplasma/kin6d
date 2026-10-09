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
  by edge midpoint indices for edges (1,2), (2,3), (3,1). The geometry remains
  affine and the boundary remains polygonal. Supplied-profile runs also export
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

## Smooth-domain geometry primitive

`kin6d_solovev_geometry` supplies a **candidate geometry callback only** for
`H(R,Z)=(R^2-r0^2)^2+4 R^2 Z^2=level`, centered at `(r0,0)`. Its explicit
inputs are `r0`, positive `level<r0^4`, `0<tau<1`, and a CCW convex polygon.
They contain no source, pressure, F, flux amplitude or reference-field values.
The caller must independently establish the reviewed positive-R convex-domain,
polygon inclusion, inner homothetic inclusion and positive radial-derivative
premises. Constructor finite/vertex/facet checks do not prove those premises.

`initialize_solovev_geometry` constructs a reusable geometry record.
`sample_solovev_geometry` returns `geometry_sample_t`: physical `value(2)`,
`jacobian(a,i)` and `hessian(a,i,j)`, the active polygon `facet`, an `identity`
flag, interface flags and `status`. Its inputs are polygon physical coordinates;
the output is the radial inverse G mapping that polygon to the smooth domain.
Inside the inner homothetic domain G is identity, including at the axis.
On an annular smooth piece, a FortNum bracketed root solves the quartic boundary
equation and the stable smaller quadratic root gives G's radial scale.
FortSym-generated leaves evaluate that geometry and its equation-level spatial
jets; bisection iterations are not differentiated. `solovev_geometry_jvp` and
`solovev_geometry_vjp` act on point directions and Euclidean Cartesian weights
at fixed geometry. Parameter/facet/topology derivatives are not exposed.

Status0 means successful **floating candidate evaluation**. Status1 rejects
invalid/uninitialized input, status2 reports an inner or sector derivative cut,
and status3 reports failed root/arithmetic/positivity checks. A native guard band
around a cut also reports status2. No one-sided derivative is selected silently;
nonzero status leaves value, J, H and derivative products zero, with only status,
facet and interface metadata available. Boundary and near-interface samples may
be rejected by arithmetic guards; success on all valid polygon points is not
claimed. `certified` remains false; `root_residual` on success is a diagnostic,
not a root-error enclosure or a positive-Jacobian certificate.
The 64-epsilon bracket padding and scale-range allowance are candidate
evaluation recovery guards, not rigorous roundoff estimates or admission bounds.

The map is piecewise smooth and is not globally C2. Its Hessian tensor order is
physical-component first; FortFEM's generic sampled tensor uses `(i,j,a)`.
Actual certified mapped panels, root/Jacobian arithmetic budgets and common
native geometry/assembly/
readback wiring remain open. This primitive changes no mesh format, assembled
GS equation, existing affine solve, coarea or classical-current capability and
admits no new scientific PDE run. The CTest
`grad_shafranov_solovev_geometry_primitive` checks noncircular independent formal
series jets, known radial forward/inverse and curvature identities, nonsymmetric
Cartesian directions, pairing, axis behavior and explicit cut rejection.

`kin6d_solovev_cell_jet:solovev_cell_scalar_jet` connects this sampled geometry
to FortFEM's generic scalar-jet transform. Inputs are a CCW affine triangle,
point `xi` inside the reference triangle, and the gradient and Hessian of the
scalar pullback `u o T` in reference coordinates. FortFEM's P1 geometry supplies
`A(xi)` and its constant Jacobian `M`.
For `T=G o A`, the adapter forms `JT=JG M` and
`HT(a,:,:)=transpose(M) HG(a,:,:) M`, then explicitly reorders `(a,i,j)` to
FortFEM's `(i,j,a)`. FortFEM computes the physical gradient and Hessian,
including `Href-sum_a gradient(a)*HT(:,:,a)` before inverse transformation.

The returned geometry contains the composed point, `JT(a,i)` and `HT(a,i,j)`;
the inverse Jacobian has order `(i,a)`. Facet and interface tags retain their
callback meaning. `identity` describes G, so an identity sample still returns
the affine cell Jacobian M. The root residual remains G's diagnostic.
Finite reference jets and vertices, a point inside the reference triangle,
positive computed affine/composed determinants and finite resulting jets are
required. KIN interface status2 is propagated only from G; FortFEM determinant
or arithmetic failures become KIN arithmetic status3. All numerical outputs
remain zero on failure; diagnostic metadata may remain. `certified` is false.

CTest `grad_shafranov_solovev_cell_scalar_jet` uses an independent noncircular
formal-series map jet and a Cartesian quadratic to recover physical derivatives
on a nonsymmetric cell. Omitting curvature or permuting tensor components fails
that oracle. Signed cylindrical field/current diagnostics, the identity branch,
cuts and invalid/arithmetic failure clearing are also checked. This candidate
adapter supplies no global geometry or roundoff certificate and changes no
assembly, mesh, PDE solve or scientific admission.
