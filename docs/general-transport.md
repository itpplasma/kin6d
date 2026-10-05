# General damped transport laboratory

The direct steady solver accepts an independently supplied source and certifies
its reconstructed candidate from a continuous residual. It does not evaluate
the true solution to obtain an error radius. This is the first general linear
certificate beyond the retained analytical laboratories; general nonlinear
self-consistent kinetics remains outside its scope.

## Model and input

The dimensionless equation is

\[
v\partial_x f+\nu f=g,\qquad \nu>0,
\quad x\in\mathbb R/(2\pi\mathbb Z),\quad v\in[-1,1].
\]

The prescribed electric field is zero. Velocity is a parameter of this
operator, so neither a velocity boundary condition nor artificial velocity
periodicity is imposed. The source and damping form an open balance; closed
collision invariants are not asserted.

Input consists of finite Fourier coefficients in x, each continuous
piecewise-linear in v. Ordered velocity nodes and all coefficient bytes define
the exact represented source. The zero-mode sine coefficient must vanish.
Invalid parameters, shapes, nonfinite values and incompatible grids fail
explicitly. No source tail or whole-velocity-line physics is silently appended.
Approximating a different physical source needs a separately bounded load error.

The demonstration supplies

\[
g(x,v)=1+v/8+(1-v)\cos(x)/4+v\sin(2x)/8,\qquad \nu=1.
\]

These are source definitions, not coefficients derived from a prescribed
solution. Dyadic nodes and coefficients make the 16- and 64-cell runs represent
exactly the same continuous source. The candidate is obtained by direct nodal
mode solves and reconstructed with the supplied finite Fourier/P1 layout.
All project-owned numerical expressions, derivatives, moment integrals and
enclosure leaves are emitted by the retained FortSym generator.

## Continuous state and tangent bands

For a reconstructed candidate \(\bar f\), the checker encloses

\[
r=g-(\nu+v\partial_x)\bar f.
\]

On each velocity cell the residual Fourier coefficients are quadratic
polynomials. The default evaluator uses their Bernstein coefficient hulls,
then bounds Fourier factors over the entire spatial period. The independently
structured strict audit evaluates whole-cell interval expressions directly;
losing correlations may make its bounds wider. Both currently use outward
interval arithmetic. Neither path substitutes an analytical solution.

The continuous characteristic inverse has norm at most \(1/\nu\), giving

\[
\|f-\bar f\|_\infty\le R/\nu,\qquad R\ge\|r\|_\infty.
\]

The bound covers all velocity cells and all x, including v=0. Its premises,
independent behavioral oracles and precise theorem status are documented in
[the transport certificate](transport-certificate.md). It is a numerical
continuous bound for this named model, not a general nonlinear existence or
formal whole-program claim.

For the nu tangent the source is held fixed. The implicit tangent equation is
\((\nu+v\partial_x)h=-f\). The checker bounds the tangent residual using
\(\bar f\), then includes the already certified state error:

\[
\|h-\bar h\|_\infty\le(R_h+r_f)/\nu.
\]

The implementation also accepts general source/nu directions with the matching
state-error amplification. It does not differentiate solver iterations.
The nodal VJP uses a declared Euclidean coefficient pairing; a continuous
adjoint goal gradient is not supplied or advertised by this demonstration.

## Observable and accuracy

The nontrivial density-mode observable is

\[
Q(f)=\frac{1}{2\pi}\int_{-1}^{1}\int_0^{2\pi}
\cos(x)f(x,v)\,dx\,dv
=\frac12\int_{-1}^{1}a_1(v)\,dv.
\]

The generated trapezoid expression integrates each P1 candidate cell exactly;
interval evaluation and reduction enclose its arithmetic. The coarse functional
norm bound 2 converts the state radius to a goal error, and the same argument
encloses the fixed-source nu derivative of Q. Both complete intervals include
integration and endpoint arithmetic.

The declared physical scale is 0.25 for Q and its nu derivative. Default
acceptance requires interval radius at most `1e-8 + 0.01*0.25`; an interval
evaluation of that target supplies a conservative comparison endpoint. Radius
is the conservative FortNum distance from the rounded interval midpoint to
either endpoint. A wider valid enclosure
is reported with `requested_radius_met=F`; it is not relabeled as a failed
containment theorem. No percentage is divided by a potentially zero candidate.

## Runnable figures

After the normal project configure/build:

```sh
./build/transport_demo
./build/benchmark_transport sweep
```

`transport_demo` writes:

- `output/general-transport.png`: supplied source, two state resolutions and
  fixed-source nu tangents at v=5/16, in dimensionless units;
- `output/general-transport-krook-limit.png`: the separately declared zero-mode
  slice at v=0, with source 3/4 and nu=1, compared to the earlier polynomial
  Krook reference at the same point and source;
- `output/transport-cell-bands.csv`: complete source/state/tangent intervals for
  each plotted x cell. Decimal endpoints encode binary64 values for round trips.

The main slice lies inside a coarse velocity cell, exposing reconstruction
error. The plotting schedule evaluates the candidate over every full x cell
and surrounds that range by its global certified radius. Step envelopes cover
the entire displayed cells; connected endpoint bands are not used to certify
the trigonometric curves between samples. The source ranges include interval
evaluation. Candidate lines are visual guides to the declared reconstruction.

The blue solid and orange dashed resolution curves use the color-vision-safe
Okabe--Ito blue/vermillion pair, with distinct line styles in grayscale. All
bands retain visible boundary lines. The earlier Krook reference is black and
dashed, with its own rigorous range. Its full velocity source differs from the
constant-source limiting run: the overlay asserts only the matching v=0
subproblem. Reference evaluation occurs solely in that comparison and never
in the production checker.

Prepared state/tangent packets are immutable during sampling. Certificate
preparation is global; batched point/cell postprocessing reuses it without a
global residual audit per sample. Changes to sources, candidates, parameters
or tangent arrays require preparation again. Plotting uses native Fortran
FortPlot; no Python runtime or postprocessor is required.

## Benchmark scope

`benchmark_transport [cells|sweep]` uses two modes, nu=1 and the same source.
The sweep uses 32, 64, 256 and 4096 velocity cells. A single cell count must be
a power of two to retain the exactly matched source preset; this restriction
belongs to the benchmark, not the general library.

Each evaluator reports achieved state/tangent radii, Q/nu-derivative intervals,
target status, attributed real-buffer bytes and five raw CPU timer trials.
The phases are primal solve, standalone primal preparation, standalone
recomputing acceptance, tangent solve, combined state/tangent preparation,
96-cell band/goal postprocessing, and cheap prepared target acceptance.
The combined total counts one primal solve, one tangent solve, one combined
certificate, postprocessing and prepared acceptance. The standalone diagnostic
preparation/recheck columns are excluded from that total. A separate first-call
measurement times the complete combined pipeline. Setup-plus-first and
setup-plus-warm totals add measured initialization; these are CPU costs and
do not claim cold-cache or whole-process startup measurements. Repeated warm
totals retain the initialized buffers.

Validation and output occur outside raw timed phases. The numerical APIs'
required domain checks remain part of their measured cost. Allocation occurs
at initialization; numerical solve, certificate and postprocessing loops reuse
caller-owned buffers. Attributed bytes describe those real buffers, not process
RSS or the plotting library. Benchmark postprocessing produces band data and
goal intervals; PNG rendering is a separate demo activity.

The Bernstein and strict-box paths are compared by both cost and achieved
width. A faster path whose bounds miss the requested target does not establish
matched-accuracy superiority. These CPU measurements establish local baselines;
they do not rank external codes, GPUs or distributed execution. The objective
remains the fastest measured total time at matched model, requested outputs,
accuracy and claim strength, including the requested derivatives and bands.

## Local execution evidence

Native Release and Debug app builds and the registered
`general_transport_and_tangent_oracles` test passed. The Fortran demo ran in
both profiles, including a public Fo Debug execution. Both generated figures
were inspected in color and grayscale: full-cell envelopes, the separate
Krook slice, labels and legends remained readable. The intentionally constant
limiting curves produce FortPlot horizontal-line warnings; their explicit
vertical scale preserves the intended comparison.

The Release sweep used GNU Fortran 16.2.1 with `-O3` on an AMD Ryzen 9 5950X,
pinned to CPU 4, with five raw trials per evaluator/resolution. Other host work
was not isolated. The retained command reproduces the workload; timings are
local evidence, not portable promises.

| Cells | Bernstein combined warm median | Strict-box combined warm median | Both goal targets met, Bernstein / strict |
|---|---|---|---|
| 32 | 4.134 ms | 4.120 ms | no / no |
| 64 | 4.460 ms | 4.398 ms | yes / no |
| 256 | 6.317 ms | 6.003 ms | yes / no |
| 4096 | 43.363 ms | 38.036 ms | yes / yes |

For the 64-cell Bernstein run, the whole-domain state radius was
`3.8984020059891483e-4` and the nu-tangent radius was
`9.2504017592987236e-4`. The complete Q interval was approximately
`[0.195559688, 0.197119049]`; its derivative interval was approximately
`[-0.126839909, -0.123139748]`. Rounded values here are summaries; the executable
prints round-trip binary64 endpoints and conservatively evaluates target
acceptance. The 16-cell demo comparison retains valid wider bands and reports
both requested goal targets as unmet.

Chris&AI
