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
positivity and rejects folding. It fixes all boundary omega increments so
Newton cannot move the prescribed nonaxisymmetric wall. The existing
transversal/axis increment space and Newton solver are retained.

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
