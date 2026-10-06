# Animated particles and Eulerian free streaming

The first representation comparison solves

\[
\partial_t f+v\partial_x f=0,\qquad
x\in\mathbb R/(2\pi\mathbb Z),\quad v\in[-1,1].
\]

Velocity is constant, the electric field is prescribed zero, and no velocity
boundary condition is required. The demo starts from
\(f_0(x,v)=1+\tfrac12\cos x\). All quantities are nondimensional;
the plotted spatial density is \(n(x)=\tfrac12\int_{-1}^1 f(x,v)\,dv\).
This is the first transport rung of the representation laboratory. It does
not yet include a particle deposit/field solve/gather cycle or field feedback.

## Numerical carriers and derivatives

The Eulerian state stores finite Fourier coefficients in position and
piecewise-linear coefficients on a velocity mesh. The linear free-streaming
operator is diagonal in each Fourier mode and velocity node. The solver
applies its rotation propagator to the supplied initial coefficients;
it does not evaluate the demo's analytical distribution to construct the
candidate. This exact subflow avoids time-stepping error for this operator.
Velocity reconstruction and floating-point arithmetic still introduce errors.
Arbitrary finite Fourier/P1 initial data are supported within this model;
error in projecting an external distribution into that model is not included.

The particles carry deterministic tensor midpoint quadrature weights from the
same initial data and move along unwrapped characteristics. Wrapping position
onto the plot is a display operation. Marker area represents the fixed weight;
tracked cohorts make individual motion visible. These are characteristic
particles, with no Monte Carlo confidence claim.

FortSym retains the definitions and emits the characteristic and Fourier
propagator primal/JVP/VJP families. Unwrapped motion and the propagator are
smooth in their declared finite inputs. A wrapped coordinate has a chart jump
at a period crossing; derivatives of that display jump are not a production
derivative. Particle loading and weak-observable derivatives are documented
with their actual implementation status rather than inferred from the pusher.
The emitted VJPs use Euclidean input/output pairings. Mesh, integer Fourier
degree, quadrature counts and cohort selection are fixed discrete choices.
Certified continuous density derivatives and parameter-derivative enclosures
are currently unavailable; a tested nodal derivative is not such an enclosure.

## Continuous distribution enclosure

The checker encloses the nodal rotations using FortNum outward arithmetic.
For each initial velocity cell and spatial harmonic, define

\[
P_k(x,v,t)=a_k(v)\cos(k(x-vt))+b_k(v)\sin(k(x-vt)).
\]

The initial coefficients are linear on that cell. With \(h\) its width,

\[
|\partial_v^2P_k|\le M_k(t)
=2|kt|(|a_k'|+|b_k'|)
+(kt)^2\max_{\text{cell endpoints}}(|a_k|+|b_k|).
\]

Thus the difference between this physical harmonic and its nodal linear
interpolant is at most \(M_kh^2/8\). Sum the mode remainders and include the
enclosed nodal arithmetic. This bound applies directly to the physical
harmonic: \(M_k\) does not bound the sum of the two rotated coefficient
curvatures. Density enclosures integrate the pointwise bound over velocity.
No unknown exact solution enters the checker; the admitted free-streaming
propagator supplies the a priori derivative estimate.

The density panel compares the Eulerian result with Fourier moments of the
weighted particles. These are weak observables; the particle Dirac measure
has no pointwise distribution enclosure. Particle quadrature, arithmetic and
continuum velocity interpolation have distinct budgets. Connected plot bands
include an additional chord remainder over their entire displayed segments.
For a retained degree \(K\), the particle curve is the Fourier polynomial
with mean \(\sum_p w_p/(4\pi)\), cosine coefficients
\(\sum_p w_p\cos(kx_p)/(2\pi)\) and analogous sine coefficients.
The ideal position midpoint count must exceed \(2K\); velocity quadrature
subdivisions must stay inside each original P1 cell. The admitted continuum
retains its finite spatial degree under free streaming, so an enclosed
quadrature error for that projection can enclose its true density.
Arbitrary supplied initial coefficients may give signed quadrature weights.
The demo's initial distribution is at least \(1/2\), so its markers have
positive weights and an area interpretation; positivity is not inferred for
arbitrary packets.
Heatmaps display candidate centers; the companion velocity slice displays the
full distribution enclosure. Rendering pixels are not a separately proved
physical observable.

## Verification and cost

Independent higher-precision tests evaluate translated arbitrary Fourier/P1
initial data off mesh nodes, including nonuniform meshes and several modes.
Derivative tests check independent directional derivatives and adjoint
pairings. The analytical density of the demo,
\(1+\tfrac12(\sin t/t)\cos x\), with its continuous value at \(t=0\),
is a benchmark oracle, not a candidate or checker input.

Measure numerical evolution, interval preparation and requested observable
checking separately from plotting and movie encoding. Record the extra time
and total multiplier against the matching solve alone. This small diagonal
operator is a reference case, not evidence of certification overhead or
performance leadership for self-consistent Vlasov--Poisson.

## Run and retained evidence

```sh
./scripts/demo.sh --snapshot
# The complete existing gallery, without movie encoding.
./build/streaming_comparison_demo
# output/particles-vlasov.mp4: 193 frames (t=0..12) plus a 1 s hold, 24 fps.
./build/streaming_comparison_demo --snapshot 193
# Select one frame after recording the full numerical evidence.
./build/benchmark_streaming_comparison
```

The comparison uses 128 Eulerian velocity cells, 97 spatial display nodes,
and 32 by 48 deterministic particles. The constant initial velocity profile
is represented identically on both carriers, although the particle source
packet needs only the two endpoints. Density and slice CSVs retain the actual
connected-band endpoints. `streaming-particles.csv` retains marker data and
rigorous ideal-quadrature geometry/weight packets; `streaming-comparison.csv`
retains every frame's numerical costs and statuses. CSV evidence follows the
current demo workflow; the planned canonical UX0 scientific-result schema is
not implemented by this gallery.

The two particle multipliers have different denominators: `particle_factor`
uses characteristic advance alone, while `particle_observable_factor` uses
advance plus Fourier projection. Absolute extra times against both baselines
are retained alongside them. Eulerian complete time includes the prepared
interval state/density packet; band postprocessing is recorded separately.
Rendering and encoding are excluded. These are correctness-reference timings,
not a claim of an optimized workhorse or comparative performance leadership.

On the Ryzen 9 5950X, CPU 4, GNU Fortran Release, five warm trials of 50
repetitions at the final time gave:

| Evolution stage | Median | Observed range |
|---|---:|---:|
| Primal nodal evolution | 3.06 microseconds | 3.06–3.38 microseconds |
| Same evolution plus interval preparation | 1.26826 milliseconds | 1.26498–1.28606 milliseconds |

The extra preparation costs 1.26520 milliseconds, or 414.46 times primal
evolution including the primal work. The velocity interpolation radius is
0.00219727; the sampling packets include additional rounding. This measurement
excludes observable sampling, plotting, derivatives and cold setup; their
costs must not be inferred from this factor. Raw trials are emitted by the
benchmark. Repeated interval trigonometry in sampling remains a measured
optimization target; cached basis/geometry and aggregate rounding budgets
must preserve the same enclosures when optimized.

The animation's full certified sampling currently costs about 0.115 seconds
per frame, beyond the evolution/packet times above. It dominates numerical
postprocessing for this small diagonal problem. The per-frame ledger records
it; a matching uncertified postprocessing baseline is still unavailable.
Cache shared Fourier factors and geometry before claiming a fast production
workhorse. This is adjacent performance work, not a reason to omit bands.

Since 2026-10-06 the movie samples 193 frames over the same interval and its
fourth panel plots both density bands relative to the Eulerian density, where
they are wider than a pixel; the velocity-slice bands remain in
`streaming-slice-bands.csv`. The earlier gallery's maximum connected
density-band half-widths over all 121
frames and 97 nodes were 0.00220920 (Eulerian) and 0.00738337 (particle
projection). Both meet the declared nominal radius 0.01000001 at physical
scale 1. This is an absolute scale-based target, not a relative-error claim
at zeros. Sampling arithmetic and spatial chord padding are included.

Verification on base `8120dcf874739d8feed2158281d07200840a0290`:

- Native Release and checked Debug CTest: 8/8 pass in each profile.
- Pinned Fo delegated Release gates: both new registered tests pass. The
  existing resident CMake capture blocker remains recorded in Fo's plan;
  this is not resident-green evidence.
- `scripts/demo.sh --snapshot`: the complete default gallery passes.
- Independent mathematical review: final source patch SHA256
  `4936f701ac8b6500f29c2f1266ded5188f830ead138606f08c18c8706a70a066`
  passes the stated composition; two identified bridges were repaired.
- Color/grayscale inspection of the snapshot and encoded initial/middle/final
  frames passes. MP4: H.264, 1320 by 880, 121 frames, 24 fps, 5.041667 seconds.

This is reviewed continuous mathematics conditional on FortNum's inherited
arithmetic contracts, with independent numerical oracles. No new formal Lean,
self-consistent-field, interacting-particle, CUDA or sensitivity certificate
is claimed. The microscopic Coulomb and Debye-volume scan additions are plans.
