# Native P2/P3 coarea candidate

`kin6d_gs_coarea::gs_level_coarea` evaluates

\[
K(\psi_s)=\oint_{\psi=\psi_s}\frac{dl}{R|\nabla\psi|},
\qquad q=\frac{FK}{2\pi}.
\]

The flux convention has a strict interior minimum at the magnetic axis and
increases outwards. `R` is positive. The fixed geometry is an affine or
isoparametric triangular mesh with six P2 or ten P3 nodal coefficients per
cell, shared globally on faces. Gradients and Hessians come from the same native polynomial as the contour; there is no
smooth interpolant layered over a P1 or P2 state.

The candidate filters edges using precomputed polynomial ranges, locates contour
intersections with native mesh faces, sorts their polar angles using FortNum,
and integrates each smooth cell arc with FortNum Gauss–Legendre quadrature. `angular_points` is a budget: the
per-arc Gauss order is its integer quotient by 16. Doubling a budget that is a
multiple of 16 doubles this order, including at high orders. A contour wholly
inside the axis cell uses sixteen regular angular sectors. A constant-Dirichlet
outer contour is integrated directly along its boundary faces. On curved
cells, roots are found in the affine chart of the reference polynomial;
every quadrature weight is multiplied by `det(DX)/det(DA)*R_chart/R_physical`.
This is the change of variables for the coarea measure, including the boundary
contour. Physical volume and toroidal flux use the same isoparametric map.

The nodal/level JVP differentiates the converged contour equation. It includes
axis displacement and the Leibniz terms at moving mesh-face intersections.
The event retains its native face intersection point; endpoint weights use
that point directly rather than re-solving a rounded radial endpoint.
These terms contain the difference of the two one-sided coarea weights;
subtracting them silently would differentiate a different map. The generator
`gen/gen_gs_coarea.f90` retains the mathematical leaves and emits primal, JVP
and VJP leaves through the pinned FortSym engine. The assembled coarea API
currently supports a P2 affine-geometry JVP under Euclidean nodal/level pairings;
its assembled VJP and shape derivative are unavailable. Geometry is inactive.

Supported contours must be closed, nested, regular and star-shaped about the
qualified native axis, with an outward-monotone polynomial trace in each cell
used by an arc. A contour through a mesh vertex, a tangent face intersection,
a missing/open contour, an axis cusp, or a boundary endpoint with zero
one-sided gradient returns failure. In particular, no finite boundary coarea
is claimed for a P2 corner cell with two constant-Dirichlet edges and a
vanishing endpoint gradient. These are representation/domain qualifications,
not claims about the corresponding smooth Grad–Shafranov equilibrium.

Returned values and derivatives are numerical candidates. Comparing successive
quadrature orders supplies an indicator, not an interval remainder or a
whole-contour certificate. Nodal interpolation, mesh faces and the straight
boundary polygon contribute separate errors. This increment does not certify
field/current/force error or supply a coupled prescribed-q equilibrium solver.

The registered native test independently checks exact finite-radius circular
quadratics, inverse amplitude, translated-axis and gauge derivatives. A
separate C0 P2 fixture has genuine gradient jumps and checks the assembled JVP
against central differences of separately located contours. Deliberately
omitting the moving-face terms fails that fixture. Existing P1/P2 physical,
current, gauge, symbolic, variational and output contracts remain gates.

`kin6d_gs_readback::read_gs_equilibrium` combines these contours with the
constant FFprime law and returns signed q, the physical axis, volume and full
edge toroidal flux. Zero radial/angular settings select automatic budgets:
`max(17,ceil(sqrt(Ncells))+1)` radial samples and `16*(degree+2)` angular budget.
Thus the radial spacing follows the element scale, while each cell arc gets
Gauss order `degree+2` and becomes shorter with mesh refinement.

The radial grid is uniform in `rho=sqrt(s_pol)`. Integrating the local cubic
interpolant of `F*K*dpsi/drho` with two-point Gauss quadrature gives cumulative
Phi. The rho factor gives the axis value zero measure. This replaces the
mesh-independent trapezoidal floor. Independent element-area quadrature gives
Phi_edge and volume; `flux_consistency` reports its relative difference from
the contour integral. The area rule has Duffy order `degree+2` and reuses basis
values, without reconstructing unused field Hessians. Normalized s_tor and
Phi_edge together give Phi(s). Nonmonotone cumulative profiles fail explicitly.

`gs_readback_q` evaluates a local four-node cubic in s_tor, including one-sided
stencils at endpoints; queries outside [0,1] return NaN. The native `--q` batch
command applies this same contract to `profiles.dat`, so external consumers
need not choose another interpolant. Returned values remain numerical
candidates; radial/angular quadrature and interpolation are not enclosures.

The readback-order CTests solve the exact curved A10 Solovev problem on 48,
96 and 192 boundary nodes. They require the three-mesh fitted q rate within
0.3 of degree 2 or 3; individual meshes may fluctuate. Separate off-grid q and
toroidal-label errors must retain at least that field order. Volume/Phi geometry
rates, F reversal and native export/readback remain tested. The work CTest
counts edge polynomial solves and cell-location trials per contour and rejects
area-like growth when the mesh spacing halves. No wall-clock threshold is used.

P3 ray evaluation caches exact cell polynomial coefficients contracted from
FortFEM's generated basis jets. Its FortSym-generated value/gradient schedule
has JVP/VJP leaves and is compared with the retained direct-basis schedule
(`reference_cubic=.true.`). The cache is local to one readback call, not a new
field fit or a separate physical representation.

`sample_gs_many` uses FortFEM's walking locator for batches of physical field
queries. It retains the scalar physical inverse-map/first-cell behavior on
faces and outside the affine boundary polygon, including curved boundary cells.
The producer keeps one native batch reader open while a result is used, avoiding
repeated mesh parsing for the common driver's field/force-balance queries.

P3 edge and ray roots use FortNum Lagrange interpolation and Brent solves,
split at the polynomial's stationary points so every crossing is considered.
The same arc quadrature, physical measure factor and signed field convention
apply to both degrees. P3 contour derivatives return unavailable status 20.
The exported axis-location diagnostic uses local quartic recovery; the contour
levels and axis-q use the represented element minimum.
