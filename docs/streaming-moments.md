# Chapter 1: particles, distribution and moments

The `streaming_moments_demo` animation uses the existing certified free-streaming
model, with prescribed zero electric field, position period `2*pi` and velocity
in `[-1,1]`. Its exact binary64 Fourier/P1 input represents

\[
f_0(x,v)=(1+\tfrac12\cos x)(1+\tfrac12v).
\]

The distribution is positive, with minimum `1/4`. This bounded-velocity teaching
profile replaces the original textbook Maxwellian-plus-beam illustration; it is
not a Maxwellian, an unbounded velocity approximation or an interacting plasma.
All plotted coordinates and moments are nondimensional.
For physical period `L`, velocity scale `v_ref` and density scale `n0`, the
plotted variables are `X=2*pi*x/L`, `V=v/v_ref`, `T=2*pi*v_ref*t/L` and
`F=2*v_ref*f/n0`, each with unit `[1]`. The plotted density and bulk velocity
are `n/n0` and `u/v_ref`, respectively. The equations below use these
dimensionless variables.

The four panels show 1536 weighted characteristic markers, the Eulerian
candidate distribution, normalized density and bulk velocity. Marker area is
proportional to the fixed positive quadrature weight. Particle wrapping changes
only the display chart. The markers and Eulerian coefficients evolve from the
same initial distribution; moment profiles integrate the Eulerian reconstruction,
rather than claiming estimates from the plotted sample.

## Moment contract

The normalized moments are

\[
n=\tfrac12\int_{-1}^1 f\,dv,\qquad
j=\tfrac12\int_{-1}^1 vf\,dv,\qquad u=j/n.
\]

Density uses the existing strict free-streaming moment packet. For an affine
coefficient with endpoint values `a_l,a_r`, the generated current leaf integrates
it exactly in real arithmetic:

\[
\tfrac12\int_{v_l}^{v_r}v I(a)\,dv
=\frac{v_r-v_l}{12}
  [(2v_l+v_r)a_l+(v_l+2v_r)a_r].
\]

FortNum evaluates this identity outward on the enclosed nodal rotations. If
`r` bounds the continuous distribution's velocity interpolation error uniformly,
then the current remainder is at most
`(r/2)*integral[-1,1]|v|dv = r/2`. Nodal rotation, moment accumulation and query
arithmetic remain in the strict coefficient packets; `r/2` alone is not their
complete numerical error budget.

A bulk-velocity interval divides the current interval by the density interval
only when the latter has a strictly positive lower endpoint. Otherwise the
bulk query explicitly fails. Dependencies between these intervals can make the
quotient conservative, but cannot invalidate containment.

Every display x-cell is queried as an interval. Its full density/current/bulk
ranges widen both adjacent vertices. Connecting the resulting endpoints
linearly encloses the entire cell: both lower vertices lie below its minimum,
and both upper vertices lie above its maximum. This avoids assuming a
second-derivative bound for the true quotient. The wider band includes spatial
variation within the display cell as well as numerical error.

The distribution heatmap shows candidate centers with a fixed color scale.
The existing strict distribution point queries are retained in
`output/particles-to-moments-distribution.csv`: each displayed node has its
frame, time, coordinates, candidate value and full pointwise lower/upper bounds.
These are the already computed sampling bounds, with no extra numerical work.
CSV writing follows the postprocessing timer. Heatmap pixels
are not a separate certified observable. No continuous sensitivity enclosure
is claimed. The generated moment primal/JVP/VJP leaves have Euclidean
input/output pairings; a checked leaf derivative is not a certified physical
moment derivative through the complete evolution.

## Rendering and evidence

Run `streaming_moments_demo --snapshot 1`, `--snapshot 97` or `--snapshot 193`
for the first, middle and last frames. The full run produces 193 frames over
`t=0..8` at 24 frames per second, followed by a one-second hold of the final
frame (81 frames at 16 fps before 2026-10-06). It writes `particles-to-moments.png`,
`particles-to-moments.mp4`, full displayed moment bands and CPU timings in
`output/`. Generated media remain untracked.

Timings compare the matching 64-cell primal nodal evolution with evolution plus
the strict distribution/moment packet, using 20 repetitions per frame.
Postprocessing is measured separately; rendering and encoding are excluded.
The factor is against the primal nodal solver, not a fully matched uncertified
moment-reporting pipeline. Cold setup and certified sensitivities are not timed.

Independent tests integrate the translated teaching profile analytically in
binary128 and check density, current, bulk and whole-cell containment. An exact
Simpson integral supplies an independent current-cell oracle; its complex-step
derivative checks the generated JVP. The VJP has a separate transpose pairing
check. Invalid time, stale preparation and unproved positive density fail.
These checks support reviewed continuous arithmetic evidence, not a new formal
proof of FortNum primitives or a self-consistent kinetic solver.

Suggested teaching caption: **Weighted particles, an evolved distribution and
its moments.** A positive bounded-velocity profile free-streams with zero field.
The lower panels retain density and bulk flow, while the distribution preserves
velocity-dependent structure. Shaded connected bands enclose every displayed
cell. This illustrates moment reduction without implying a fluid closure;
position has period `2*pi`, velocity lies in `[-1,1]`, and all axes are normalized.

Local worker evidence (2026-10-05): the independent moment test passed with
checked GNU Fortran in an isolated scratch build. Snapshots at frames 1, 41
and 81 were rendered; the first/middle color views and final grayscale view
were inspected. The full movie was validated with `ffprobe`: 1200 by 800,
81 frames, 16 fps, approximately 5.06 seconds. The accompanying CSVs contain 81 unique
frame records and 97 displayed moment vertices per frame. Native registered
CTest/Fo integration evidence belongs to the controller's final receipt.
Encoded-frame inspection exposed retained raster pixels after FortPlot's
figure `clear()`. The gallery pins the owning-repository repair
`9949d818b2ba5f43d582d6aff58c8032cdeec3de`; the controller rerenders and
checks encoded frames with that repair before publication.
