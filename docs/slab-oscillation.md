# Cold planar electron slab oscillation

The Chapter 3 teaching model is a cold, unmagnetized, collisionless electron slab
with fixed singly charged ions and infinite planar transverse geometry. The
small collective displacement produces boundary sheet charges and a restoring
field. With `tau=omega_pe*t`, `q=xi/xi0` and `p=dq/dtau`, the declared reduced
model is `q'=p`, `p'=-q`, with independently supplied finite initial `q0,p0`.
It is not a three-dimensional periodic uniform-translation restoring OCP or a
microscopic N-body simulation. Pressure, collisions, magnetic forces and the
finite-amplitude correction to the reduced slab model are excluded.

`slab_oscillation_demo` shows electron/ion teaching glyphs, the Eulerian view of
that same collective state, its time trace, and a magnified numerical error
panel. Position is `x/a`, where `a` is the displayed slab half-width and
`xi0/a=0.1`. Density is `n/n0`, displacement is `xi/xi0`, time is `omega_pe*t`,
and normalized field/force are `E/E0` and `F/F0`; all plotted units are `[1]`.
Species lanes are display positions. Glyphs move with the shared collective
coordinate; they are neither independent physical particles nor a PIC method.
The box density is a view of the rigid-slab reduced model, not an independently
integrated Eulerian kinetic PDE or a proof of finite-amplitude slab dynamics.

The untrusted candidate uses the pinned upstream FortNum `cash_karp_step`, a
Cash–Karp RK5(4) step, with a generated oscillator RHS. The accepted time mesh is
fixed and every reported solve reuses preallocated trajectory and stage buffers.
The embedded RK error is not used as a rigorous bound. No analytical solution
is supplied to the production checker. FortSym generates the RHS, Hermite
reconstruction, residual, derivative and enclosure leaves from retained Fortran
source. Consumer builds use the committed emitted leaves without requiring CAS.

## Continuous residual certificate

The certificate is the generic FortNum module
`fortnum_ode_residual_certificate` (pinned FortNum `7f6a0ba`). On each
strictly ordered time cell `[tl,tr]` with `h=tr-tl`, each component is its own
cubic Hermite: `q_h` with endpoint values `ql,qr` and slopes `pl,pr`, and
`p_h` with endpoint values `pl,pr` and slopes `-ql,-qr` (exact negation; the
slopes are the generated RHS at the nodes). The residual
`r=(q_h'/h-p_h, p_h'/h+q_h)` (in the cell variable `s`) is a cubic vector
polynomial; FortNum encloses its four Bernstein control vectors in outward
interval arithmetic and bounds `sup_s ||r(s)||_2` by the largest control
norm. The interval action of `J` is the generated interval RHS leaf.

Bernstein basis functions are nonnegative and sum to one, so this is a
continuous-cell claim, not a sampled residual. Both reconstructions are
globally C1 and hence absolutely continuous, so the following
variation-of-constants argument applies across all cell boundaries.

Until 2026-10 the slab used a second-order reconstruction (`p_h=q_h'`,
scalar residual `q_h''+q_h`). That residual is only `O(h^2)` because
`q_h''` is affine per cell; the componentwise residual is `O(h^3)`. On the
128-step two-period animation case the final radius fell from `6.61e-3` to
`3.30e-4`; on 40/80 steps over `tau=8` from `1.77e-2/4.30e-3` to
`1.78e-3/2.22e-4`.

For the exact reduced-model state and its error `e`,
`e'=J*e-r`, `J=[[0,1],[-1,0]]`. Since `J` is skew-symmetric, its evolution
preserves the Euclidean norm. Duhamel gives
`||e(t)||2 <= ||e(0)||2 + integral_0^t ||r(s)||2 ds`.
The generated initial discrepancy bound uses `abs(qh0-q0)+abs(ph0-p0)`, a safe
upper bound on the Euclidean norm; it avoids a strict square-root domain test
on an outward sum-of-squares interval with a slightly negative rounded lower
endpoint. Every full-cell product and accumulated radius uses outward interval
arithmetic. Both displacement and velocity are therefore enclosed by the
same radius. Strict interval reconstruction also covers evaluation rounding.
The checker accepts arbitrary finite candidate traces; corruption must enlarge
the resulting radius rather than masquerade as a small certified error.

Point queries use the radius at the containing cell's right endpoint. Whole
time-interval queries split at the original cells, enclose each cubic/derivative
and take their interval hull. Invalid/nonfinite inputs, nonmonotone meshes,
nonfinite enclosures and stale candidate-only packets fail explicitly.

The connected time-trace bands use strict endpoint packets widened by the
uniform final residual radius and a chord remainder `M*(tb-ta)^2/8`, where `M`
bounds the candidate's second derivative on every intersected original cell.
Although `q_h''` can jump at a knot, `q_h'` is continuous and Lipschitz: `q_h` is
in `W^{2,infinity}`. The standard chord remainder follows by integrating its
bounded weak second derivative, so it applies across reconstruction knots.
The error panel also adds `abs(q0)+abs(p0)` to the curvature bound for the
analytical comparator. Its full bounds are preserved in the CSV and plot;
no smoothing or width clipping is applied.

## Density edges and derivatives

The numerical density view translates a box by the numerical collective
coordinate. The displacement interval becomes intervals for the two sharp
boundaries. Density uncertainty inside either boundary strip is the full
possible `[0,n0]`; a small displacement radius is not presented as a small
pointwise density error. Outside these strips the box density is determined by
the common shift. Both edge intervals share the same correlated displacement.

RHS and Hermite primal/JVP/VJP leaves share one FortSym definition. Initial-state
sensitivity is the linear variational IVP on the same fixed mesh, and its own
residual certificate encloses the continuous directional state. The terminal
VJP is a numerical continuous-adjoint candidate using generated RHS VJP and
the upstream RK step on the reversed trace, with Euclidean initial/terminal
pairings. For this constant linear RHS its real-arithmetic stability polynomial
has the transpose identity; no bitwise floating-point identity is implied.
Continuous VJP error bounds are unavailable. Mesh/count changes and density
edge crossings are not smooth derivative events. HVP is deferred.

## Verification and reproduction

`slab_oscillation_oracles` checks arbitrary initial data against independent
binary128 harmonic-oscillator evaluation, whole-time interval inclusion,
residual-radius refinement, zero state, a hand-derived two-node polynomial
Bernstein control bound and sampled residual,
corrupted traces, invalid meshes and generated rational derivative fixtures.
It also checks the continuous tangent, the adjoint against its independent
analytical flow, a tolerance-based pairing, and connected bands across knots.

`benchmark_slab_oscillation` records 15 raw CPU samples of 160 repetitions each:
solve alone versus the same solve plus strict certificate, with validation
outside timed loops and setup/rendering excluded. Its 128-step, two-period
case requires the final state radius to be at most 0.01. A cost factor includes
certificate work and is not a claim about other integrators or tolerances.
The default animation also checks that every connected displacement half-band,
including the chord margin, is at most 0.01 in the declared displacement scale.
The solver API permits finer meshes for stricter targets.

`slab_oscillation_demo [--snapshot [frame]]` writes
`output/plasma-oscillation.csv`, `.png` and, by default, `.mp4`. The main trace
retains the analytical comparison as an independent overlay; the checker does
not call it. The gallery renderer converts the movie to GIF. Plot artifacts
are not committed. The palette follows `docs/figure-style.md`: numerical
curves blue with translucent bands, analytical references dashed black,
electrons vermillion, and species marker shape as a redundant encoding. Since
2026-10-06 the default movie has 193 frames over two periods (8 s at 24 fps)
plus a 1 s hold; the finer frame spacing narrows the chord margin, so the
connected half-band is now below 0.002.
