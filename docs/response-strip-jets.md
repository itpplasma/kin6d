# Imported strip jets and radial Fourier preconditioning

These modules extend the isolated, candidate full-pressure strip. They do not
couple it to the ideal outer response or impose a physical cut boundary.

`kin6d_response_strip_jets` samples the existing `psi_grid_t` representation.
Its radial coordinate is `u=sqrt(psin)`; the poloidal coordinate is the geometric
angle about the imported magnetic axis. The first ray crossing defines the map.
Its Jacobian and second derivatives follow by implicit differentiation of
`psi(X(u,theta))=psi_axis+u^2*(psi_edge-psi_axis)`. Root iterations are inactive.
The inverse map is `d(u,theta)/d(R,Z)`. The axis chart is singular and rejected;
only regular, star-shaped annuli inside the represented flux domain are supported.

The returned magnetic derivative and pressure Hessian use physical Cartesian
components `(R,phi,Z)` at `phi=0`, including derivatives of the cylindrical basis.
The full pressure spline, including its nonzero edge value, is retained. Second
profile derivatives are evaluated by the existing cardinal spline definition in
`profiles_of_psin`; no differencing step or alternate pressure fit is introduced.
FortSym source `gen/gen_response_strip_jets.f90` emits the map/field jet leaves
and their JVP/VJP. Input coefficients, axis-location sensitivities and a smooth
import/root-selection derivative are not supplied by the sampling wrapper.
The generated leaves expose derivatives of their explicit jet inputs.

The fluid background current follows the existing outer `equilibrium_point`
force-balance definition. The independently derived `curl_current` need not
agree with it for an imported spline candidate. `ampere_defect` reports their
absolute difference in A/m^2; `label_defect` reports the solved level residual.
These are indicators, with no representation-error or physical-model enclosure.
They must not be silently promoted to a stationary-equilibrium certificate.

`kin6d_response_strip_preconditioner` changes solver orchestration, retaining the
same six-field strip action and its prescribed isolated ring rows. Fields and
interior equations are rotated into local normal/tangent components. Nine radial
colors recover five-node first-derivative and composed second-derivative bands.
One-sided endpoint rows reach six radial nodes, so each complex factor has scalar
half-bandwidth 41. No whole-strip Jacobian or full theta-pair matrix is stored.

Three low-mode probes determine the quadratic continuous angular symbol, and
six probes at the represented Fourier edges retain dominant collocation aliasing.
Smaller grids use every mode. Omitted off-diagonal harmonic blocks and alias tails
belong only to the inverse approximation. FGMRES still applies the original full
operator. There are at most 486 streamed strip actions in setup, independent of
poloidal resolution; band storage and factor/solve work grow linearly with grid
size at fixed bandwidth. Setup also includes the existing quadratic-cost angular
collocation action. No mesh-independent iteration count is promised.

Two-sided row/column equilibration addresses mixed SI units and small rate
scales, without a diagonal shift. Krylov RHS/state normalization reduces
finite-precision cancellation; the original RHS is reconstructed for the
independent physical residual check, including for small-amplitude drives.
The caller supplies the same positive equation scales used to measure the full
residual. Initialization must receive the same
geometry and dimensions as the strip; rebuild after changing geometry,
equilibrium, resistivity or gamma. Each solver instance owns its mutable scratch.
FortNum's shared read-only FFT plan and per-transform reusable workspace keep
mixed-radix/Bluestein scratch outside the hot call. Independent mode
factorizations and radial factor solves run in parallel, with no shared LAPACK
status or temporary band matrix in the factorization loop.
The required upstream FFT change is retained on local branch
`work/sol-rb2-jets-fft`, pending controller integration. Development CMake uses its
source override; no dependency hash pin is introduced by this task.

The inverse approximation has a Hermitian storage transpose. Forward solves
use `A/scale`; adjoint solves use `(A/scale)^*`. A physical adjoint additionally
requires the caller's volume/surface Riesz weights. Drive tangents solve the same
linear equation. Equilibrium/mesh/rate derivatives of the complete solve remain
caller-owned through `A dx=db-dA x`; differentiating the preconditioner is not
required. Krylov uses half the requested residual tolerance to allow for
reconstruction and the final physical action. Stagnation at that tighter
internal target is accepted only when the independently recomputed physical
residual meets the original requested tolerance. Other Krylov failures remain
failures. Residuals are re-evaluated after Krylov termination. Returned states
remain candidates; neither an inverse bound nor a discrete/continuous certificate
is available. Gamma=0 is rejected rather than regularizing a static nullspace.

Registered tests check analytic polynomial forcing on a curved pressure-balanced
annulus, both signs of n, retained pressure forces, inverse/transposed inverse,
generated jet products, and an independent elliptical-map/pressure oracle.
The finite rate ladder covers eta/mu0 from 0.2 to 2e-7 and eta=0 at fixed gamma,
then gamma from `0.7+0.2i` to `0.007+0.002i` with eta=0. These bounded manufactured
limits do not establish a TC24 ideal-sheet or Glasser limit. The resolution
ladder uses gamma=`2+0.2i` on 9x9, 17x17 and 33x33 grids with a fixed state-error
and residual target. The finer grid at gamma=`0.7+0.2i` did not meet that target;
its failed exploration remains in external evidence, not a qualified limit.

Configure `KIN6D_STRIP_EQDSK_TEST` to register the imported TC24 jet checks. They
compare B/map/pressure jets with independent spatial differences and retain the
actual pressure. Both columns of the implicit map Hessian and pressure Hessian
are checked. An additional imported-background inverse oracle sets `A=xi=0`,
`Phi=exp(i*n*phi)` and independently supplies the toroidal `-grad Phi` drive.
On a 7x25 annulus spanning psin=0.94..0.96, it covers eta=1e-7, 1e-9 and zero
Ohm m at gamma=`2000+500i`, then gamma=`600+150i` and `200+50i` at eta=0,
with rho=1e-7 kg/m^3. A bounded restart of 600 is needed; restart 240 failed
at the first imported rate. This verifies sampling and an isolated inverse,
not a coupled response solve or the pre-declared TC24 sheet-limit ladders.
`benchmarks/response-strip/run.sh` records matched-input one/two/four-thread
trials, setup/solve/warm-apply time, iterations, state error, attributed array
storage and process peak RSS. Outputs must be written outside source history.
