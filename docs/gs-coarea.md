# Native P2/P3 coarea candidate

`kin6d_gs_coarea::gs_level_coarea` evaluates

\[
K(\psi_s)=\oint_{\psi=\psi_s}\frac{dl}{R|\nabla\psi|},
\qquad q=\frac{FK}{2\pi}.
\]

The flux convention has a strict interior minimum at the magnetic axis and
increases outwards. `R` is positive. The fixed geometry is an affine or
isoparametric triangular mesh with six globally shared P2 nodal coefficients
per cell. Gradients and
Hessians come from the same native polynomial as the contour; there is no
smooth interpolant layered over a P1 or P2 state.

The candidate finds polynomial contour intersections with native mesh faces,
sorts their polar angles using FortNum, and integrates each smooth cell arc
with FortNum Gauss–Legendre quadrature. `angular_points` is a budget: the
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
edge toroidal flux. Its radial grid is uniform in sqrt(s_pol); cumulative
trapezoidal integration of `F*K` supplies normalized s_tor. The independent
area integral supplies Phi_edge. `flux_consistency` reports their relative
difference and can be reduced by radial refinement. Neither radial nor angular
quadrature is a rigorous enclosure. The curved readback CTest measures q,
volume, Phi and flux-label convergence against the exact physical-LCFS case,
checks reversal of F, and tests native mesh export/readback round trips.

P3 edge and ray roots use FortNum Lagrange interpolation and Brent solves,
split at the polynomial's stationary points so every crossing is considered.
The same arc quadrature, physical measure factor and signed field convention
apply to both degrees. P3 contour derivatives return unavailable status 20.
The exported axis-location diagnostic uses local quartic recovery; the contour
levels and axis-q use the represented element minimum.
