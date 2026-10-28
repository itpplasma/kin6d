# Boundary and toroidal-profile initializer

`kin6d_ideal_seed_3d` initializes the existing frozen-in map solver from a
Fourier boundary and polynomial profiles. It does not use an equilibrium's
interior fields. `ideal_seed_3d_init` accepts all four R/Z sine/cosine channels
in the convention `m*theta - n*nfp*phi`, with general positive `nfp`. The
coefficient array has shape `(4,2*ntor+1,0:mpol-1)`; its second dimension is
ordered from `-ntor` to `ntor`, and the first is `(Rcos,Rsin,Zcos,Zsin)`.

Profiles use normalized **toroidal** flux `s`, with ascending monomial iota
and pressure coefficients. These differ from the poloidal-flux profiles of
`ideal_seed`. In reference-disk coordinates `(u,v)=(Rref-2,Zref)`, set
`s=u^2+v^2`, `b3=Phi_edge/pi`, and
`chi=Phi_edge/(2*pi)*integral_1^s iota(t) dt`. The frozen density is
`(-chi_v,chi_u,b3)`. Its divergence vanishes and its field-line rotation
number is the supplied iota. The generator proves these two algebraic
identities exactly, including the chain rule for an arbitrary smooth iota.

The boundary is extended through disk harmonics `(u+i*v)^m`, which are
regular at the axis. The initializer samples the represented Jacobian for
uniform orientation and rejects folding. Clockwise and counterclockwise
poloidal angles are admitted; the prescribed iota refers to that input angle. Without the
field angle below, the pinned wall and the frozen form force the input angle to be the
straight-field-line angle of the target on the boundary. It fixes all boundary omega increments so
Newton cannot move the prescribed nonaxisymmetric wall. The native Newton solver is retained. This cylindrical positive-toroidal-flux
slice uses its existing gauge that fixes the toroidal coordinate and allows
both R/Z increments. Fixed reference-disk radial directions can become tangent
to seeded surfaces, even for a clockwise circle. FortSym verifies the
flux-form invariance of the corresponding angle relabelling, conditional on a
positive relabelling Jacobian. The target must admit B_phi>0 in this chart;
toroidal-field reversal is outside this initializer's contract.

## Boundary field-line angle

The wall nodes are pinned at the input angle, `X(rho=1, theta, phi) =
wall(theta, phi)`, and the frozen form `b = (F/2 pi) ds^dtheta - chi'(s) ds^dphi`
(`F = Phi_edge`) has field lines that are straight in the reference angle
`theta`. A pinned map therefore makes the field lines on the boundary straight in
the input angle, which means that VMEC's `lambda` must vanish there. This is an
extra tangential constraint, because the fixed-boundary problem constrains only
`xi.n`. With a non-straight-field-line input angle the error does not converge
away. For example, the axisymmetric sheared oracle in its chart angle keeps
2-4 % B errors, and AUG with `|lambda| ~ 1.1 rad` on the LCFS shows a rigid
interior shift.

`ideal_seed_3d_init(..., angle_mpol=M, angle_ntor=N)` keeps the wall pinned and
relabels the frozen field instead. The reference is relabelled along its
surfaces by `psi(s, theta, phi) = (s, theta + eta, phi)`, and the frozen form
becomes its pull-back

    b' = psi^* b = (F/2 pi) ds^d(theta + eta) - chi'(s) ds^dphi
       = b + (F/2 pi) ds^d eta,

which is linear in `eta`. In the seed channels this means `chi_(u,v) += -(F/pi)
eta_phi (u, v)` and `b3 += (F/pi) eta_theta`. The field lines of `b'` are
straight in `theta + eta` and have the same iota.

- **Unknowns.** The coefficients of `eta = sum Re/Im c_mn (u + i v)^m
  exp(-i n nfp phi)` for `m = 0, 1 <= n <= N` and `1 <= m <= M, |n| <= N`, with
  two real coefficients per pair. The extension `(u+iv)^m` is harmonic and
  regular at the axis, so `b'` is smooth and fixed by the boundary trace
  `eta(rho=1)`. The extension does not add interior freedom, because the map
  already owns that.
- **Same energy and same problem.** `b'` is the pull-back of `b` under a
  diffeomorphism that preserves every reference surface (when `1 + eta_theta > 0`,
  which is checked at every quadrature point). It is therefore exactly closed,
  tangent to the surfaces, and carries the same toroidal flux (`oint d eta = 0`),
  poloidal flux and iota on every surface, while the pressure `p(s)` does not
  change. A state `(X, eta)` has the physical field `X_* b' = (X o psi^-1)_* b`,
  which is the field of the map `Y = X o psi^-1`. Its energy is the ideal
  energy `W[Y]` of the same frozen-in problem, and `Y` maps the reference boundary
  onto the same wall with the poloidal reparametrisation `psi^-1` on the
  boundary. Any orientation-preserving poloidal reparametrisation of the wall at
  fixed physical `phi` is of this form, so the union over `eta` of the pinned
  problems is the fixed-boundary problem with tangential sliding. The wall
  geometry never moves, which keeps the discretisation conforming: every
  discrete state is an admissible ideal state on the fixed (interpolated)
  domain, so the discrete energy stays an upper bound up to quadrature. The
  sliding-node variant (kin6d#119) lacks this property, because there the
  minimiser bulges the edges outwards.
- **Gauge.** The `m = n = 0` constant is omitted. A rigid shift of `theta` leaves
  `b` unchanged, so it is an exact null. The toroidal labels stay physical
  (position gauge, `omega = 0`), so the relabelling is purely poloidal. Interior
  relabellings are owned by the map and not by `eta`.
- **Reachability.** The exact equilibrium has straight field lines in
  `theta* = theta_in + lambda(theta_in, phi)` on the boundary. Hence
  `eta(rho=1) = lambda` (up to a constant) together with `X = X_PEST o psi`
  reproduces it for any input angle with geometric `phi`, provided that `lambda`
  is resolved by the `M, N` series (spectrally for an analytic boundary).

The pointwise leaves (the mode `eta, eta_theta, eta_phi`, the form increment,
and `dw/db`, `d2w/db2`, `d2w/dv db` through the chart jets) are generated by
`gen/gen_ideal_field_angle.f90` and `gen/gen_ideal_map.f90`. The generators
prove that the mode derivatives are the exact angle derivatives, that the
increment is closed, that the field lines are straight in `theta + eta`, and
that the energy is homogeneous of degree 2 in `b`. `kin6d_ideal_map` appends the
coefficients after the node and axis dofs and adds the dense angle rows to every
cell block (`hessian_angle`). The same Newton, LDL^T and trust-region path then
applies, with the angle step capped like the rotation.

Choose `M` at most about half the number of boundary nodes. On the axisymmetric
sheared oracle (P3, chart angle) `M = nb/2` gives B errors of 8.7e-4, 3.2e-4
and 1.3e-4 at nb = 8, 12, 16. The straight-field-line input gives 1.4e-3,
7.1e-4 and 2.0e-4. A larger `M` lowers the energy but uses the high harmonics
to compensate discretisation error (b8 M8: 1.2e-3, b16 M16: 3.0e-4). The
option is off by default. On 3D references that are ideal saddles (the integer
family with constant iota 2, the sheared m8n4 case) a toroidal `N > 0` opens
the saddle, and Newton leaves for a lower-energy helical state, as the
sliding-node variant did. CTest `ideal_map_field_angle` checks the leaves, the
map gradient and Hessian, and the chart-angle regression.

`ideal_seed_disk_mesh` builds the common P2/P3 reference mesh. Generated
leaves include primal/JVP/VJP; the initializer parallelizes independent nodal
coefficient assembly. FP64 is the current reference path. The seed is a
candidate: sampled admission is not a whole-cell enclosure, and no continuous
state/derivative bound or successful equilibrium solve follows from it.
Iota must have one nonzero sign on the sampled profile-admission grid. This
is a diagnostic restriction of the current map gauge, not a certified
polynomial positivity test. Arbitrary nonconvex boundaries may fail harmonic
extension; there is no claim that every Fourier boundary is admitted.

Build `test_ideal_seed_3d` and run CTest `ideal_seed_3d_boundary_profiles` for
independent off-grid boundary reconstruction, toroidal flux, sheared profile
integration, JVP finite differences and VJP pairing, fixed-wall constraints,
thread parity and folded-boundary rejection. With the pinned FortSym build,
`ideal_seed_3d_regeneration` proves the identities and compares regenerated
leaves against the retained source.
