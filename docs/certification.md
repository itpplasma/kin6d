# Certification specification

Status: design and proof obligations. No certificate implementation exists
yet. This specification applies to the [kin6d design](../DESIGN.md).

## Result contract

The [model family](model-family.md) extends this contract to enclosed
equilibria, linear forced response, spectra and hybrid trajectories. Imported
equilibrium uncertainty is propagated into response claims. Singular ideal
operators and non-normal kinetic/resistive operators need their own constrained
inverse, spectral isolation or resolvent/growth theorem. A finite eigenvalue
interval is not a steady nonlinear root or a forced-response stability bound.
The [atlas](benchmarks.md) independently certifies reference evaluation; exact
identities and asymptotic checks support but do not replace this contract.

A steady certificate names a continuous operator `F`, its boundary data,
constraints, input parameters, and function spaces. It supplies a reconstructed
candidate `ubar` and proves, under listed assumptions,

$$
\exists U^*\in\overline B_X(\bar U,r): F(U^*)=0.
$$

Uniqueness is a separate claim, usually within this ball and on the constrained
space. It is neither global uniqueness nor dynamical stability. Positivity,
regularity, and physical admissibility need their own verified premises.
An `L²` error ball alone does not establish pointwise positivity or pointwise
error bounds.

The report must state the norm and its dimensional scales. For example, a
distribution norm must define momentum weights, species weights, spatial
measure, and any derivatives. Field components need explicit normalization.
Choose the first norm for the first theorem; the final kinetic norm remains
an analytical decision. A finite Euclidean coefficient norm is not silently
identified with a continuous function norm.

| Result state | Supported statement | Accepted scientific solution |
|---|---|---|
| `candidate` | Floating-point approximation and diagnostics | No |
| `certified-discrete` | Enclosure of a finite algebraic problem | No |
| `certified-continuous` | Continuous theorem with reviewed analysis and checked rigorous arithmetic | Yes, with that proof level stated |
| `formally-certified-continuous` | Continuous theorem and instantiated certificate checked in Lean | Yes |
| `certification-failed` | A specific bound, assumption, or checker obligation was not met | No |

A case declares the required proof level and observable widths before a run.
No automatic fallback can satisfy a request for formal certification by
returning a reviewed analytical argument. The first closed reference case
targets the formal level. Larger models may have useful candidates for a
long time before their continuous certification gate is passed.

## Candidate, generator, and checker

The numerical solver is untrusted. It can use ordinary FP64, approximate
preconditioners, GPUs, mixed precision, and nondeterministic reductions when
those capabilities become worthwhile. The certificate generator is also
untrusted: its output must be checked against the mathematical definitions.
Both are Fortran applications; the later fast path may invoke generated CUDA
leaves. The rigorous evaluator uses fortnum interval/ball arithmetic,
rigorous transforms, and verified linear algebra within their checked scope.
Project-owned numerical leaves, including enclosure and derivative leaves,
must be emitted by fortsym under the [kernel policy](tooling.md#kernel-generation-policy).

Sharing a symbolic DAG reduces formula drift but does not establish that the
DAG expresses the intended continuous physics. Retain independent analytical
oracles and a checked bridge from the mathematical definition to the kernel
and certificate data. AD results and generic numerical libraries are not
formal proofs merely because their output is generated or called verified.

The checker verifies directed interval or ball arithmetic, reconstruction,
residual bounds, stability bounds, and the theorem's premises. Exact rational
or dyadic endpoints make numerical data reproducible. The Lean layer proves
soundness of the relevant calculations and applies the continuous theorem.
[LeanCert](https://github.com/alerad/leancert) is a candidate arithmetic tool;
its existing root and integration machinery is not a ready-made kinetic PDE
certifier.

The formal trust base includes Lean's kernel, admitted foundational axioms,
the modeled equations and input interpretation, and the checked bridge from
serialized data to mathematical objects. Pin Lean and mathlib versions,
record theorem axioms, and reject `sorry`, unsafe proof shortcuts, and
unproved analytical assumptions masquerading as established conclusions.
Hashes bind artifacts to a theorem application; hashes themselves prove no
numerical inequality.

## Continuous contraction theorem

Let `F:D subset X -> Y` be continuously Fréchet differentiable on a neighborhood
of a closed ball in a Banach space. Construct a bounded linear approximate
inverse `A:Y -> X` and `T(U)=U-A F(U)`. Verify that the ball stays inside the
domain where all operators, boundaries, and assumptions are defined.

Compute rigorous bounds

$$
\|AF(\bar U)\|_X\le Y_0,\qquad
\sup_{U\in\overline B_X(\bar U,r)}
\|I-A DF(U)\|_{X\to X}\le Z(r)<1.
$$

If `Y_0 + Z(r) r <= r`, the map takes the ball into itself and is a contraction.
Completeness gives a unique fixed point there. To conclude `F(U*)=0`, prove
that `A` is injective on the relevant residual space, or provide an equivalent
argument eliminating the residual kernel. A finite projected inverse alone
does not meet this requirement.

For example, a bound `Z(r)=Z_0+Z_1 r` reduces the self-mapping test to
`Y_0 + (Z_0-1)r + Z_1 r² <= 0`, together with `Z(r)<1`. All constants and the
chosen radius require rigorous arithmetic. This describes an intended
mechanism, rather than a proof that the final plasma operator satisfies its
premises. [RadiiPolynomial.jl](https://olivierhnt.github.io/RadiiPolynomial.jl/stable/)
provides examples of this style of computer-assisted argument.

Collision invariants, potential gauges, translations, and equilibrium families
can make the derivative singular. Prescribe compatible densities, total
momentum, energy, gauges, or phase conditions as appropriate, then establish
the theorem on the resulting constrained space. Do not hide a nullspace in
an empirical inverse estimate.

## Resolved state and infinite complement

Embed the numerical reconstruction in the continuous space using `P_h`.
The operator, approximate inverse, and derivative bounds must control both
the retained block and `(I-P_h)X`, including their coupling. Tail control is
an analytical estimate, not extrapolation from the last few coefficients.

The following defects may be needed:

- Spatial reconstruction, DG interface jumps, metric factors, and boundary
  traces, including grazing orbit effects.
- Radial momentum approximation, regularity at zero, omitted momentum domain,
  and collision-potential far-field conditions.
- Angular harmonics above `L_s`, toroidal modes outside the retained set, and
  nonlinear coupling to both complements.
- Collision-potential solves, cross-species projection, singular-kernel
  treatment, quadrature, and nonlinear evaluation.
- Characteristic integration, interpolation, and any split-step consistency
  defect relative to the stationary differential equation.
- Field reconstruction, exterior boundary treatment, input interpolation,
  divergence constraints, and source/boundary compatibility.
- Nonlinear solve residual and arithmetic rounding.

For a discontinuous numerical state, either construct a conforming candidate
with a rigorously bounded reconstruction difference or use a theorem that
directly accounts for jumps in an appropriate weak norm. A cellwise strong
residual omitting interface distributions is insufficient.

## Error budget and refinement

Begin with rigorous residual contributions in a common norm. If
`F(ubar)=sum_i d_i` and `||A d_i|| <= Y_i`, then `Y_0 <= sum_i Y_i` follows by
the triangle inequality. Record derivative and stability contributions too.
Do not label a list of raw residual norms as additive solution-error bounds.

The contraction theorem converts this combined defect and stability margin
into a radius. Assigning separate solution-error contributions requires an
additional justified inequality. Include cross-couplings and model-parameter
variation where the theorem uses them.

Refinement acts on spatial resolution, radial extent/order, `L_s`, toroidal
bandwidth, quadrature, nonlinear tolerance, or interval precision according
to the dominant checked contribution. Recompute the certificate after each
change. A larger approximation space does not guarantee a smaller radius:
the stability estimate or interval dependency may worsen. Report the actual
certified widths and keep the previous successful result available.

Compare refinements in the same continuous norm and under the same problem
definition. Intersect observable intervals only when a theorem identifies
the enclosed solution as the same branch. Two locally unique roots can still
be different physical steady states.

## Certified observables

For a bounded linear functional `Q`,

$$
Q(U^*)\in Q(\bar U)+[-\|Q\|_{X^*}r,\|Q\|_{X^*}r].
$$

Rigorous quadrature and arithmetic additionally enclose `Q(ubar)`. Nonlinear
observables use derivative bounds over the ball or an independent validated
evaluation. Ratios require a denominator bounded away from zero. Wall fluxes
need trace regularity; entropy needs control near zero; reaction rates need
a specified cross-section model and its input uncertainty. None follows from
a generic norm-ball label alone.

Report numerical intervals separately from physical-input intervals and
operator differences. If two certified models give intervals `I_A,I_B`,
their observable difference lies in `I_B-I_A`. This does not establish that
either model encloses the real plasma.

An expensive model can sometimes be certified around a cheap-model candidate
by checking `A F_rich(ubar)` and full derivative/tail bounds for the rich
operator. One correction evaluation without those surrounding bounds is
insufficient. Validated continuation can identify a common branch when its
parameter-dependent existence and uniqueness obligations are proved.

## First continuous certificate

Start with a driven, spatially homogeneous, bounded one-dimensional momentum
Fokker–Planck problem. This is a reference problem, separate from production.
For `p in [-P,P]`, choose positive smooth `M(p)` and `D(p)` and solve

$$
-\partial_p\bigl(DM\partial_p(f/M)\bigr)=S(p),\qquad
f(\pm P)=M(\pm P).
$$

For example, use a normalized Gaussian `M` and the manufactured solution
`f_exact=M(1+a cos(pi p/(2P)))` with `0<a<1`; derive `S` analytically.
The boundary condition fixes the homogeneous nullspace. This driven model
does not claim closed-system particle or energy conservation.

Let `ebar=f-fbar`, with `fbar` a reconstruction satisfying the boundaries.
Define `||e||² = integral e²/M dp`, and the continuous residual
`rbar=S+partial_p(DM partial_p(fbar/M))`. Dirichlet Poincaré and bounds on
the positive coefficients give a coercivity constant

$$
\lambda_{\rm lb}=
D_{\min}\frac{M_{\min}}{M_{\max}}\left(\frac{\pi}{2P}\right)^2,
\qquad
\|f-fbar\|_{L^2(M^{-1})}
\le \|rbar\|_{L^2(M^{-1})}/\lambda_{\rm lb}.
$$

Indeed, writing `g=ebar/M` gives energy `integral DM (g')²`, bounded below
by `lambda_lb integral M g²`; testing the error equation and applying
Cauchy–Schwarz yields the enclosure. Existence comes from the coercive weak
problem. This is a planned proof to formalize, not an existing Lean theorem.
Require enough reconstruction regularity for the stated `L²` residual, or
replace it with a proved weak-residual estimate.

The first demo also needs a pointwise band. Let `Rbar` rigorously bound the
weighted residual norm and use `g=(f-fbar)/M`, with zero endpoint traces.
The energy inequality implies
`||g'||_L2 <= Rbar/sqrt(lambda_lb D_min M_min)`. Integration from either
endpoint and Cauchy–Schwarz then give, for every `p in [-P,P]`,

$$
|f(p)-fbar(p)|\le
M(p)\sqrt{\min(p+P,P-p)}\,
\frac{Rbar}{\sqrt{\lambda_{\rm lb}D_{\min}M_{\min}}}.
$$

Formalize this one-dimensional trace argument with the coercivity theorem
and evaluate its envelope with rigorous generated arithmetic. Enclose the
reconstructed candidate and envelope over whole plotting cells; interpolation
between sampled endpoints must not narrow the true band. This pointwise
bridge is specific to the reference problem and does not transfer to six
dimensions without appropriate regularity estimates.

Check coefficient bounds, residual integration, and the Poincaré constant
rigorously. Certify one bounded moment of `f`. Repeat at two resolutions and
verify both contain the independent analytic solution; report their actual
widths. Generate distribution bands, resolution overlays, and a moment-interval
plot through the Fortran postprocessor. Perturb a residual certificate or
coefficient bound and require the
checker to reject it. No finite-domain result is advertised as a solution on
the whole momentum line.

## Species reductions and device scope

Certify each selected species/field/pair system in its own continuous spaces,
with its measures, transforms, closures, and boundaries. A FK certificate
does not become a DK/GK certificate by dropping coefficients. Model identity
includes fixed reservoirs and their nonzero momentum/energy exchanges.

Matching observables from independently certified full and reduced problems
give a rigorous interval difference. A reduction-error guarantee for an
uncomputed parent problem additionally requires a justified lift, continuous
reduction defect, initial/boundary layers, and verified parent stability or
another applicable comparison theorem. State transfer and fidelity switching
carry explicit error into subsequent certificates. See
[species reduction bounds](species-models.md#certified-reductions).

Device certificates include frozen input/geometry enclosures and the actual
physical parameters. A theorem in a symmetry-restricted subspace declares
that scope. The [device gates](devices.md#production-gate) distinguish
component correctness from production maturity. Neither a geometry smoke
test nor softened scale separation discharges the device solution theorem.

## Certificate artifact

Version the format with the first executable certificate. Each accepted
artifact records the model and input hashes; basis, mesh, frame, and
normalizations; exact candidate bytes and reconstruction; norm and theorem;
checked constants, tails, boundary obligations, and radius; observable
definitions and intervals; proof level and assumptions; checker/dependency
versions; source commit; and verification result.

Keep diagnostic estimates separate from proof premises. A checker exits with
failure for missing bounds, unsupported models, malformed data, violated
assumptions, or numerical uncertainty too large to establish the theorem.
Archive the candidate and precise failed obligation so refinement or new
analysis can resume from evidence.

## Coupled models and uncertainty scope

[Coupling](coupling.md) separates frozen background, prescribed drive and
induced response. Certifying response around an imported background must
include its residual/uncertainty or explicitly state a conditional problem
with that background prescribed. Species and bulk ownership prevent duplicate
source terms; reservoir transfers belong to the stated balance.

Neutral Monte Carlo and externally sampled fluctuations have statistical
status. Fixed seeds do not give deterministic enclosures. A full coupled
certificate needs deterministic bounds on every neutral/source action used
in its premises; otherwise certify only the explicitly conditional model.
Atomic/table/surface inputs are part of the mathematical problem. A nominal
input certificate is not a physical-data uncertainty bound. To cover uncertain
inputs prove the premises uniformly over a stated input set.

Combined model/numerical intervals require a proved comparison chain to the
same target observable. Do not add unrelated paired-model differences as
remainders to physical reality. Spatial model transfers, wave partitions,
field reductions and AP numerical formulations each retain their distinct
obligations. [Extensions](extensions.md) specifies this accounting.

Boundary hits, rare losses and shock solutions require dedicated event/trace
or weak-solution estimates. A generic norm radius need not imply a narrow
interval for a discontinuous loss indicator. Unavailable bounds remain
explicitly unavailable; plots may show the candidate without a fabricated band.

Chris&AI
