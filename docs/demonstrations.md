# Runnable demonstrations and plotting

Status: the first CPU prescribed-field 1X1V slice is implemented.
See [the runnable example and evidence](first-example.md). The following
contracts govern completion of the first formal milestone and later increments.

## Runnable from the first code

The first executable increment is a complete CPU vertical slice: build the
Fortran generator, emit the numerical leaves, solve the smallest reference
case, generate/check its certificate, and run a Fortran postprocessor that
produces meaningful plots through fortplot. Document one command for this
demo and one for regenerating plots from retained run/certificate data.
The default demo needs no GPU, cluster, Python, network account, or private
discharge data. Pin and document its build/proof prerequisites.

Keep that command working as features grow. Every integrated numerical
increment preserves runnable earlier demonstrations and adds or extends a
meaningful case for its new behavior. A broken experimental branch cannot
replace the working default. Resident feedback can retain a previous working
executable, but its source/generation identity must be visible; it is not
evidence that the edited candidate passed.

The initial demo must generate real candidate and interval data. Certification
failure remains explicit, leaves diagnostics runnable, and prevents presenting
those diagnostics as a certified result. At least the retained default demo
must still complete its required certificate and interval plots.

## Figure contract

Before implementing a figure, define its purpose, data, axes/units, scale,
normalization, labels, output formats, and observable theorem. Each figure
shows the numerical candidate and the entire certified lower/upper interval
for its plotted quantity. Use a filled band with visible boundary lines,
or interval rectangles/error bars where the domain is discrete. The interval
is a deterministic enclosure, not a statistical confidence band.

| First-demo figure | Message | Representation |
|---|---|---|
| Distribution against momentum | The certified continuous solution is enclosed around the computed state | `p/P` versus explicitly normalized `f`; candidate line, pointwise/cellwise certified band, analytic exact curve |
| Resolution comparison | Certificate widths change as approximation improves | Two resolutions on the same axes; candidate lines and both interval envelopes |
| Moment interval | The requested observable lies in the reported interval | Named moment and units; estimate with full interval and exact reference value |

The norm-bound-to-pointwise bridge for the first reference is specified in
[certification](certification.md#first-continuous-certificate). For later
problems, a generic norm radius cannot be drawn directly around `f(x,p)`.
Derive the appropriate pointwise/trace bound or plot a certified bounded
observable, such as a cell-average profile, with that definition in the axis
label. A distribution slice advertised as pointwise requires a pointwise
theorem. Define suitable plotted observables when introducing each capability,
so its interval figure is part of the acceptance gate.

Plot enclosures must hold over the entire represented horizontal cells or
segments, not just at sampled nodes. Bound reconstruction and interpolation
between nodes or use rigorously enclosed cellwise step/rectangle envelopes.
Include input, quadrature, reconstruction, compression, and display rounding
where they affect the reported endpoints. Do not narrow a band by smoothing,
clipping negative values, omitting a tail, or rounding inward.

Use documented color-vision-safe colors, with redundant line styles/markers
and direct labels so comparison remains readable in grayscale. Distinguish
candidate, exact/reference curve, and interval edges. Keep legends off the
data; name model, resolution, norm/observable, proof level, and the interval's
meaning in the caption or adjacent run summary. Use physical units or state
the nondimensional scale. Use perceptually uniform maps for field images.

fortplot's checked source contains fill-between and error-bar facilities.
Confirm the exact data/segment semantics and rendering in the first demo.
Missing generic plotting behavior is repaired in fortplot main with a
focused oracle, rather than implemented as a private plotting kernel or
Python postprocessor.

## Limiting-case overlays

Use the [Analytical Benchmark Atlas](benchmarks.md) as the permanent overlay
inventory. Circular tokamak, simple helical field and large-aspect-ratio
families remain plotted beside their later geometry/transport cases. Rigorous
reference quadrature/root/special-function uncertainty is displayed as part
of the reference band. Leading asymptotes without remainder bounds are labeled
diagnostic curves, never certified exact-solution bands. MHD/hybrid comparisons
use the same physical observable with bulk/kinetic accounting stated.

Every new general capability produces a comparison figure overlaying its
limiting configuration with the retained earlier simpler implementation.
Use matching units, physical inputs, observable definitions, and reference
conventions; show both candidates and both full certified intervals. Identify
any remaining physical reduction error separately.

Examples include simple transport against its reference case, a homogeneous
collision reduction, nonlinear Landau approaching a fixed-reference linear
response, BB approaching the nonrelativistic limit, zero toroidal modes
against the retained axisymmetric case, and a constant transient against its
steady solution. A retained linear-response comparison must specify the
background and scaling rather than silently dropping reference residuals.

Also plot model-difference intervals when comparing closures. Overlaying
bands visually is useful evidence, but overlap alone does not prove equality
or common branch identity. An asymptotic reduction comparison requires its
own remainder theorem to turn it into a guaranteed reduction-error claim.

## Species and device comparisons

For each frozen device equilibrium retain its geometry plots and successive
small/production kinetic demonstrations. Label geometry-only figures as
geometry evidence; they have no plasma-solution certificate. Kinetic figures
identify the device/design revision, species/pair/field models, physical
parameters, symmetry scope, observable, and checked interval status.

Overlay fixed-ion electron cases with the corresponding prescribed-field
limit, then with evolved-ion response as it becomes available. Plot reduced
and full models using the same physical observable, including Jacobians and
polarization/magnetization where required. Display both certified bands and
their difference interval. Never identify band overlap with a proved model
reduction. Resolution and species-model refinement remain separate axes.
[Device ladder](devices.md#production-gate) and
[species comparisons](species-models.md#certified-reductions)
define these acceptance claims.

## Generated output policy

Never commit plot outputs: PNG, PDF, SVG, EPS, PS, image frames, video, or
animation artifacts. Generate them on demand under ignored build/output
directories. Commit plot source, generator source, model/plot configuration,
and small approved independent numerical fixtures instead. Keep run and
certificate metadata sufficient to reconstruct the figures from their data.
Large run data and scientific proof artifacts have a separate retention
policy; a screenshot is never a substitute for the underlying certificate.

The repository ignore rules cover image/video formats and build/run/output
directories. The implementation gate also checks tracked files so forced
adds cannot bypass the policy. CI may retain generated figures as artifacts,
with numerical checks and selected rendered inspection, without committing
them. Do not treat pixel equality as the scientific oracle.

Render and inspect the actual figure at each milestone: labels, units, band
containment, scale, reference overlays, legends, and grayscale distinctions.
Numerically check certified interval containment over the represented domain
against the analytic oracle where available. The postprocessor displays the
checked certificate status and cannot manufacture a successful one.

## Broader-model plot rules

Show background, external drive and induced response separately where needed
for interpretation, and label active/test/reservoir populations. For reaction,
wall, source and circuit examples include the relevant exchange/power ledger.
Synthetic diagnostics use physical instrument signals with their complete
numerical enclosure and separately labelled observation uncertainty.

Statistical neutral/noise confidence intervals use a distinct visual style and
legend from deterministic certified bands. Data/model intervals are labelled
with their own assumptions. No band is drawn from an L2 radius for a pointwise
quantity or from an assumed continuity of a boundary-hitting event.

Atlas cases include non-cylindrical Figure-8/map tests, open space boundaries,
reactions, pickup populations and diagnostic integrals. Geometry or field
slices can provide meaningful early plots without claiming that the full
machine solution has been certified. Retain the earlier exact/simple limit
as soon as a general implementation exists, with matched units and observables.

Chris&AI
