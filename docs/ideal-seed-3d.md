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
poloidal angles are admitted; the prescribed iota refers to that input angle. It fixes all boundary omega increments so
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

## Native D8 consumers

The case registry admits `landreman_3d_nested` and `issan_3d_nonsymmetric`
through `run_boundary_profile_equilibrium`. The caller supplies one
`equilibrium_3d_driver_t` containing the Fourier boundary, toroidal profiles,
flux and independent oracle parameters. `nested_integer_driver` prepares
matched boundary/flux/profiles for the integer-transform family only. It
uses the retained generated exact formulas; it does not import interior
oracle fields into the initializer. Matched toroidal-profile preparation for
the sheared family is unavailable. The Landreman-Paul QA fixture prescribes
current and cannot yet be reproduced by this iota-driven consumer.

`map_3d_equilibrium` retains the map and frozen form in the common field
interface, evaluates cylindrical B, pressure and the poloidal-flux label,
and solves the geometric inverse to sample physical points. It rejects a
frozen form inconsistent with the candidate's quadrature data. Generated
field/inverse leaves supply JVP/VJP, and FortSym proves all nine cofactor
inverse identities. Sampling and reference-consistency checks use OpenMP.
It is a candidate evaluator with unavailable continuous bounds.

The D8 consumer checks the attained free reduced-gradient norm against the
requested Newton tolerance even if the native solver reports success after
stalling. This coefficient norm is named explicitly; it is not the D3 dual
norm or an a posteriori state estimator. Reports additionally expose B/p
RMS errors and (for the integer family) maximum label error on the retained
3-by-24-by-24 grid in candidate coordinates. These are sampled diagnostics,
not volume norms. Force and estimator effectivity remain unavailable in this
consumer. The retained external benchmark recipe computes separate local
force diagnostics with a stencil-halving check; it omits cell-face jumps.

Run `test_equilibrium_3d TAU MPOL NTOR BOUNDARY_NODES [DEGREE]` to reproduce
one case (`TAU=0` for Landreman, `0.4` for Issan). A physical-accuracy ladder
must refine both the boundary/unknown Fourier resolution and the P2/P3 mesh;
passing the coarse CI smoke does not establish 1% accuracy.
