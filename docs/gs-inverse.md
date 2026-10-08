# Constant-q zero-beta inverse candidate

The public inverse solves fixed-boundary, axisymmetric Grad–Shafranov with
zero pressure, constant signed canonical q, and signed total toroidal flux.
It uses the native affine P2 state, contour coarea K, and F=q·2π/K. The
continuous FF′ profile is piecewise linear in normalized poloidal flux;
F² follows its exact integral. A coupled profile Newton solve differentiates
the converged weak GS equation and the native cell-arc contour equation.

Positive zero-beta homogeneity leaves q unchanged. The normalized-flux mean
FF′ is fixed to −1 during the profile solve, then the requested physical Phi
selects the positive amplitude. Fixing a single edge coefficient instead is
poorly conditioned on a refined profile. The exact piecewise-linear mean
eliminates the last coefficient; its primal/JVP/VJP algebra comes from
`gen/gen_gs_inverse.f90`. Finite pressure does not permit this normalization.

The F² domain check covers the complete closed axis-to-edge interval,
including every quadratic minimum where linear FF′ crosses from negative
to positive. Positive knot values alone are insufficient. Signed q must be
finite and nonzero. The native q/Phi signs must agree; reversal changes F
and toroidal flux while retaining the poloidal state and FF′.

Convergence requires a weak source residual, independent off-knot q probes,
coarea quadrature refinement, and physical Phi quadrature refinement. The
default q gate remains 0.1%; development fixtures with a looser declared
threshold do not establish comparison accuracy. A converged source with a
failed q gate returns status 9 and retains the physically normalized candidate.
Newton or contour failures retain their raw unnormalized states.

This candidate supports only zero pressure and constant q. Finite-pressure
and variable-q coupling, full inverse-input JVP/VJP, shape derivatives, and
rigorous whole-contour or field/current/force certificates are unavailable.
The native representation and contour domain restrictions in
`gs-coarea.md` apply. Native derivative/oracle tests are separate from
registered physical comparison runs.

The continuous piecewise-linear FFprime law is integrated on native quadratic
flux cuts. The implicit response caches the source-slope weak mass and the
stationary-axis flux response on the same cut geometry before GMRES. Its
matrix products reuse those caches; they do not repeat adaptive integration.
Numerical paired quadrature estimates and per-cell work limits fail closed.
They are not rigorous source, derivative or physical equilibrium bounds.
The independently integrated source residual is required before any
constant-q candidate is accepted. Historical candidates solved with the
seven-point profile-source rule require a new solve and cannot inherit the
corrected source residual.

Each profile Jacobian prepares one immutable response context and reuses it
for its coefficient directions. The context owns the numerical state, weak
operator, preconditioner and FFprime knots; any changed input is rejected.
It is local to that Jacobian and is prepared again for a new state. Sharing
this context changes neither integration tolerances nor response equations.

`gs_profile_response` is a primal Picard map with the affine fixed-Dirichlet
load. `gs_profile_implicit_response` differentiates coefficient forcing at
fixed boundary data: that boundary offset has zero tangent and is omitted.
The coefficient forcing follows the native source JVP with zero state,
profile-knot, pressure and mu0 directions.

The same response context integrates every nodal FFprime derivative column
in one vector traversal of the fixed native quadratic cuts. Each column
retains the original load tolerance and global estimator check. A coefficient
direction contracts those columns; if its conservative combined estimator
exceeds the original relative budget through cancellation, the original
directional source JVP reintegrates it. Fixed Dirichlet data has zero tangent.
The hat basis and source composition use the existing interpolation adjoint
and generated GS source definition.
