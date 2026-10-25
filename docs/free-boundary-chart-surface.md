# Fixed chart-backed coupling surface

`ideal_fb_surface` factors the existing PR47 ideal-map boundary constructor.
It returns a physical `bem_surface_t`, closed oriented connectivity and the
restriction from node-fast real Fourier coefficients to continuous P2 trace
nodes. Both geometry and restriction use the same counter-clockwise reference
edges and the existing generated chart jets and triangle edge basis. The
reference edge polynomial remains P2 or P3 when its trace is subdivided.

The optional positive `subdivisions` count splits each reference edge for the
BEM trace, independently of the physical polynomial geometry and pair
quadrature. Its default is one, preserving the cylindrical chart's original
connectivity and restriction. The map's admitted orientation `plasma%sigma`
determines a common permutation of physical face IDs and parameter vertices.
This preserves outward normals when a positive-orientation curve chart
replaces the negative-orientation cylindrical chart.

The caller supplies an admitted ideal-map chart and a counter-clockwise closed
reference chain. Constructor validation rejects malformed trace connectivity,
nonfinite coordinates and sampled degenerate tangents. These finite checks
are not a proof of global embedding, reach or constant map orientation.

The focused `ideal_fb_chart_surface_*` tests use the existing circular axis
and closed toroidal helix, together with the cylindrical default, at P2 and P3.
They compare the restriction with an independent Lagrange/Fourier evaluation,
surface tangents with finite differences and physical outward normals with
the independently pulled-back conormal. For the signed contravariant Piola
field `B = J*b/det(J)`, outward physical flux pairs with
`sigma*b.n_ref`; dropping the orientation factor would reverse this check.
The full-turn period of the pulled-back one-form `dphi` is `2*pi`, without
an additional field-period rescaling. These are numerical operator oracles,
not rigorous error bounds or an equilibrium certificate.

`ideal_fb_use_bem3d` uses this same constructor before the existing generic
FortFEM assembly. Its recovered R2 attachment keeps full-turn exterior energy
divided by `nfp` exactly once. The current provider pins and default pair
policy remain unchanged. The constructor tests perform no BEM assembly.
Geometry is in metres and the dimensionless restriction preserves magnetic
scalar potential in tesla metres. The full-surface DtN matrix has units of
metres, its normal-flux load has units of tesla square metres and the energy
is in joules for SI `mu0`. Real Fourier modes are constant, cosine and sine
at integer multiples of `nfp*phi`; they introduce no extra full-turn factor.

This is a D6 chart-backed ideal-map geometry slice. `ideal_fb_init` still uses
the cylindrical shell initializer, and NN2 still needs coherent mapped volume,
trace and source metrics before a nonaxisymmetric control surface can qualify
non-nested G8. Full assembly accuracy, independent interface error, geometric
return and control-surface variation remain separate unqualified gates.
