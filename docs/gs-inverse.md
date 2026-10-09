# Prescribed-q fixed-boundary equilibrium

`solve_gs_prescribed_q` solves the current Phase 3 contract on curved P2/P3:
constant physical pprime, a signed q(s_pol) table, and signed total toroidal
Phi. The boundary gauge is psi=0. q and Phi have the same nonzero sign;
the supported branch has an interior minimum of psi and nonnegative pressure.
Reversing both q and Phi reverses F and Bphi while preserving psi and FFprime.
The input F_edge is not an inverse constraint. The solver recovers it.

The current cubic-profile forward engine supplies the weak GS equation.
Continuous cubic B-splines represent FFprime(s_pol); their exact polynomial
primitive defines F squared, anchored at the solved F_edge. FortNum supplies
the spline basis and derivatives. The outer equations set the native area
integral of F/R to Phi and the weighted moments of (q_native/q_input - 1) to
zero. The weights are B_i(s)*s*(1-s): axis Hessians and one-sided boundary
values receive zero measure. The residual is sampled uniformly in sqrt(s).
This avoids feeding pointwise Hessian error into the recovered source law.
The full off-moment profile is checked separately, including the endpoints.

FortNum's numerical root solver proposes the profile coefficients; its central
secant Jacobian is a solver approximation, not a sensitivity product. The
forward state's fixed-law tangent/adjoint remains the equation-level cubic
Jacobian. Full inverse-input and curved/P3 contour sensitivities are explicitly
unavailable; the legacy affine-P2 response oracles below remain supported.
No derivative of an iteration sequence is reported.

`gs_phase1` selects this mode with `q_file` and `phi_target`, along with the
ordinary `pprime`, `mu0`, boundary and mesh options. `q_file` contains a row
count followed by ascending `s_pol q` pairs, at least four rows, spanning 0..1.
Local four-point cubic interpolation defines the input q law. `profile_file`
and `q_file` are mutually exclusive. `inverse_intervals` (default 2) is the
number of cubic FFprime spline intervals; `nonlinear_tolerance` controls the
outer residual (default 1e-9), and `q_tolerance` controls the independent
relative q defect (default 1e-3). Radial/angular budgets have the same meaning
as forward readback. Refine both mesh and profile budgets when needed.

After convergence, an interleaved radial grid and doubled angular budget
check q independently of the fitted moments. Status 9 retains the physical
candidate but rejects its q accuracy. Source, contour, boundary geometry and
profile truncation errors remain numerical errors measured by refinement.
The nonlinear recovered-flux diagnostic has status 4: it controls the frozen
source defect, without an inverse-equilibrium stability bound.

`profile-law.dat` contains the recovered pprime/FFprime law and signed F_edge.
The native field reader consumes that file; `profiles.dat` contains the actual
q, p, F, pprime and FFprime. NetCDF retains the literal/resolved namelist,
input q table, Phi constraint, solved law and provenance. The iter_tc24
producer extension uses the Phase 3 `prescribed_q` schema without using the
case's forward FFprime or nominal F_edge as constraints.

The behavioral CTests cover an independent Phase 1 Solovev q integral, P3
psi/field/source convergence, constant-q zero/finite-pressure self-convergence,
signed mesh/profile round trips, off-profile rejection, and the CLI/NetCDF
contract. They contain no frozen solution snapshots.

## Retained affine-P2 constant-q entry point

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
That reintegration accumulates directional quadrature estimates separately
from primal estimates and checks them against the directional load norm.
Cancellation can therefore return failure even with zero base profiles;
the exact zero-load pressure fixture exercises this case. These estimates
remain numerical error indicators, not rigorous enclosures.
The source JVP differentiates the same locally shifted flux and profile knots
as the primal. A common gauge direction therefore produces zero source
activity before quadrature; fixed boundary derivatives retain their original
direction. The existing gauge and all-active pairing controls cover this map.
The hat basis and source composition use the existing interpolation adjoint
and generated GS source definition.
