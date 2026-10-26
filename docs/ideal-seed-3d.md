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
poloidal angles are admitted; the prescribed iota refers to that input angle.
By default every wall node is held at its input angle (both R/Z increments of
all modes). Because the frozen density has straight field lines in the
reference angle, the pinned wall is exact only when the input angle is the
straight-field-line angle of the target on the boundary. Any other input angle
over-constrains the wall: the sheared family in its chart angle keeps a (1,0)
error of 2-4 % in B, which disappears when the same boundary is supplied in
its straight-field-line angle.

`sliding_wall=.true.` (experimental, `ideal_map_wall_slide`) constrains only
the normal displacement: wall node a sits at `wall(theta_a+g_a(phi),phi)` with
`g_a` a Fourier series on the map's toroidal modes, placed by the discrete
projection on the planes. The terms of the wall series and their first two
theta-derivatives are the generated `seed_wall_mode` leaf, and the Newton solve
uses the wall Jacobian and the curvature term of the reduced Hessian, so the
iteration converges quadratically. The rigid rotation of the reference angle
is the only null of the frozen density, and the mean slide of the wall node
nearest angle 0 is held to fix it. The result no longer depends on the input
angle, but the discretisation does not conform. The P_k edges between sliding
nodes are not the wall, and the energy minimiser enlarges the enclosed area
through them. On the axisymmetric sheared case at P3 the area grows by
3.4e-3 at 8 boundary nodes and 2.5e-4 at 20, the energy falls below the exact
one, and B converges only as 2.9/1.5/0.73/0.37 % at 8/12/16/20 nodes. On the
integer family (a saddle on fine meshes) the slide leaves the reference state
already at 12 nodes, with a label error of 27 %. The pinned default is therefore unchanged. The native Newton solver
is retained. This cylindrical positive-toroidal-flux
slice uses its existing gauge that fixes the toroidal coordinate and allows
both R/Z increments. Fixed reference-disk radial directions can become tangent
to seeded surfaces, even for a clockwise circle. FortSym verifies the
flux-form invariance of the corresponding angle relabelling, conditional on a
positive relabelling Jacobian. The target must admit B_phi>0 in this chart;
toroidal-field reversal is outside this initializer's contract.

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
thread parity and folded-boundary rejection. `ideal_map_wall_slide_derivatives`
checks the wall leaf and the reduced gradient and Hessian of a sliding 3D wall
against differences of the energy along the retraction of the solve. With the pinned FortSym build,
`ideal_seed_3d_regeneration` proves the identities and compares regenerated
leaves against the retained source.
