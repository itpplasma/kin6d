# Certification specification

Status: general certification policy and the next implementation gates. The
prescribed-field snapshot, bounded OU/FP and manufactured self-consistent
VPFP laboratory have reviewed continuous reference enclosures; the separate
polynomial Krook case also has a checked Lean application. Their exact-profile
comparisons are reference evidence, not an arbitrary-source production
certifier. [Local evidence](steady-laboratory.md) records their boundaries,
tests and proof scopes. General residual/stability certification is now
implemented for independently supplied Fourier/P1 damped transport, including
directional sensitivities and full-cell bands; see
[the certificate argument](transport-certificate.md). Forced FP and general
nonlinear self-consistent certification remain pending.

## Result states

A KIN6D result must state exactly what has been established.

| State | Meaning |
|---|---|
| candidate | floating-point result and diagnostics only |
| certified-discrete | rigorous enclosure of a finite algebraic problem |
| certified-continuous | rigorous enclosure for the named continuous model |
| formally-certified-continuous | continuous certificate additionally checked by a proof assistant |
| probabilistically-certified | statistical confidence statement with declared probability assumptions |
| certification-failed | one or more required premises did not close |

The project priority is **useful certified-continuous results**. Lean or another
proof assistant may later strengthen generic theorem/checker assurance, but it
is not required for interval arithmetic or for early scientific promotion.

A statistical interval from Monte Carlo is not a deterministic inclusion.
When deterministic discretization bounds and statistical uncertainty coexist,
report both.

## Certified model and outputs

A certificate names:

- the continuous equations/operator;
- phase-space/domain and boundary conditions;
- constraints/gauges/nullspaces;
- all physical parameters and input data;
- the function space/norm;
- the reconstructed numerical candidate;
- the observable;
- the exact stability/existence theorem used.

Numerical certification concerns the named model. It does not automatically
cover model reduction, uncertain atomic data, experimental calibration,
posterior/reconstruction uncertainty or physical adequacy.

## Reconstruction and inference uncertainty

KIN6D certification and experimental reconstruction answer different
questions. A TIAGO posterior/covariance describes uncertainty in physical and
nuisance parameters conditional on a selected forward model and diagnostic
model. A KIN6D numerical certificate encloses the solution/observable of that
forward model. A KIN6D model-reduction bound compares physical models. These
must be reported as separate components.

Do not combine a Gaussian/Laplace posterior standard deviation with a
deterministic interval radius or a model-reduction enclosure by root-sum-square
unless an explicit probabilistic dependence model justifies that operation.

If the same measurements support several distinct equilibria, certify each
KIN6D solution locally and let TIAGO fit each separately. TIAGO may report
separate fits or a finite mixture of local Laplace approximations. KIN6D does
not claim global uniqueness merely because one local radii/Krawczyk argument
closes, and TIAGO should not average physically/topologically distinct states
into one equilibrium.

Every reported scientific result carries its own bound and status: states,
fields, cell quantities, observables and requested gradients. Cover the entire
reported/requested set rather than displaying only the quantities whose bounds
closed. A pointwise band needs a pointwise or whole-cell theorem; an L2 or
energy enclosure can instead bound an appropriate cell average or weak
observable. Statistical bands state their probability assumptions. When a
requested bound is unavailable, report that quantity as a candidate or failed
certificate and name the missing premise. Do not manufacture a finite band.

Plots preserve full interval widths and reconstruction semantics. A sampled
endpoint band does not enclose a nonlinear curve between endpoints without a
remainder; whole-cell ranges do. The declared reconstruction, not a smoothed
display curve, is the object covered by the theorem.

## Accuracy targets and arithmetic budgets

For each observable or gradient component Q, declare its units, scientific
scale s_Q, relative tolerance rtol_Q and absolute floor atol_Q. Accept a
certified interval only when its radius is at most
`tau_Q = atol_Q + rtol_Q*s_Q`. The default requested relative tolerance is
1%. A near-zero or signed quantity needs a meaningful scale and absolute
floor; relative error divided by a possibly zero candidate is undefined.

Initially allocate at most 1% of tau_Q to arithmetic after stability and
observable amplification. With a negligible absolute floor this is about
`1e-4*s_Q`. This allowance is a stopping/refinement budget, not an assertion
that floating-point operations have that error. Bound residual evaluation,
inverse/stability constants, goal correction and final observable arithmetic,
then propagate them through the theorem. Discretization, source approximation,
boundary/complement and nonlinear remainder terms must also fit the total
radius. Parameter/model/statistical uncertainties retain separate budgets.

## Solver and certifier separation

The production solver is untrusted and may use:

- ordinary FP64;
- GPUs;
- approximate preconditioners;
- mixed precision;
- relaxed/nondeterministic reductions;
- heuristic continuation.

It produces a candidate \(\bar U\).

The rigorous Fortran certifier independently evaluates the quantities required
by the theorem using FortNum interval/ball and verified numerical primitives.
Its workhorse path uses FP64 centers with proved aggregate rounding and
analytical remainder bounds. A reduction bound, for example, needs its actual
operation count, magnitude bounds and rounding assumptions; a fixed percentage
inflation proves nothing. Enclose transcendental errors, conditioning,
reductions, underflow/overflow and conversion to reported endpoints. Unsupported
arithmetic conditions fail explicitly or select a supported evaluator.

Strict mode uses outward interval operations and, where needed, higher
precision for fallback and audit. Workhorse and strict modes satisfy the same
containment contract; they may differ in cost and interval tightness. The
workhorse error-bound accumulation and emitted endpoints themselves require
rigorous arithmetic. Strict-mode agreement is useful verification, not a
replacement for the workhorse bound's proof. Project-owned numerical leaves
come from retained FortSym definitions, with rigorous targets where supported.

Sharing a symbolic DAG prevents formula drift but is not itself a proof that
the formula is the intended physics. Independent mathematical/physical oracles
remain mandatory.

Prepare a certificate once for immutable candidate/source/parameter bytes and
reuse its bounds in batched postprocessing. Point/cell rendering must not
repeat a global normalization, inverse or residual audit per sample. Record
solve, certificate preparation and postprocessing time/memory separately, and
recompute the certificate when any bound input changes.

## First certificates for unknown solutions

The first general gate is damped periodic transport on a finite velocity
interval, with prescribed E=0, constant nu>0 and arbitrary declared smooth
periodic source `g=nu*B+S`:

```
L f = (nu + v*d_x) f = g,
L^-1 g(x,v) = integral[0,infinity] exp(-nu*t)*g(x-v*t,v) dt.
```

The characteristic resolvent establishes existence and `||L^-1||_infinity <=
1/nu`, including v=0. The production checker evaluates only the continuous
residual `r=g-L fbar` and concludes `||f-fbar||_infinity <= ||r||_infinity/nu`.
It never evaluates the solution formula. Require a periodic reconstruction
differentiable in x, rigorous whole-domain source/residual evaluation, all
omitted source modes and velocity-cell remainders, and a nontrivial bounded
observable. For `Q(f)=integral q*f`, the coarse goal bound is
`||q||_L1*||r||_infinity/nu`, with observable arithmetic/quadrature added.
An approximate adjoint can tighten it under the matching continuous theorem.
The continuous inverse controls every spatial mode; a finite residual alone
does not discharge that scope. No whole-velocity-line claim is made.

Next admit arbitrary compatible sources in the bounded reflecting FP problem.
Use its positive equilibrium weight to write a symmetric diffusion operator,
fix the mass nullspace, prove the continuous weighted stability bound, and
certify a weak equilibrated-flux defect. Piecewise-linear candidates have
derivative jumps; cellwise strong residuals omitting them are invalid. Account
for load integration, flux balance, boundary and mass defects, and establish
any requested pointwise bridge separately. See
[the flux contract](adjoints-structure.md#continuous-flux-and-variational-bounds).

General nonlinear self-consistent VPFP follows only after a validated local
inverse/existence argument controls kinetic/field coupling, charge
compatibility, potential gauge, boundaries and continuous complements. It does
not inherit the prescribed-field `1/nu` inverse bound. All three gates retain
independent analytical benchmarks, corrupted-data rejection, two-resolution
evidence and timing/memory at the requested observable/gradient accuracy.
Exact solutions are benchmark oracles, not production error estimators.

## Nonlinear steady states

For

\[
F(U)=0
\]

a preferred general mechanism is a Newton--Kantorovich/Krawczyk/radii-polynomial
style theorem.

With an approximate inverse \(A\) and

\[
T(U)=U-AF(U),
\]

obtain rigorous bounds such as

\[
\|AF(\bar U)\|\le Y,
\qquad
\sup_{U\in B(\bar U,r)}\|I-ADF(U)\|\le Z(r).
\]

A suitable self-mapping/contraction inequality proves a nearby solution and,
when the theorem provides it, local uniqueness.

Nullspaces from conservation, gauge, translations or equilibrium families must
be constrained explicitly. Hiding them in a pseudoinverse is not a proof.

The finite resolved space and infinite/unresolved complement must both be
controlled for a continuous certificate.

## Coercive/SPD problems

Positive/coercive operators are high-value early certificate targets.

For

\[
Au=b,\qquad A=A^\dagger\succeq\alpha I,\quad \alpha>0,
\]

and residual

\[
r=b-A\bar u,
\]

\[
\|u-\bar u\|\le \frac{\|r\|}{\alpha}.
\]

In the energy norm,

\[
\|u-\bar u\|_A^2=r^\dagger A^{-1}r.
\]

If a preconditioner \(P\) has verified spectral equivalence

\[
c_1P\preceq A\preceq c_2P,
\]

then \(r^\dagger P^{-1}r\) gives inexpensive upper/lower energy estimates.
Thus the same mathematical preconditioner can accelerate the solver and
support its error estimator.

CG/Lanczos recurrence data may provide additional energy-error bounds when
their spectral assumptions are verified.

## Continuous PDE discretization error

An algebraic residual bound does not alone certify the continuous PDE.

For elliptic/coercive FEM problems, prefer equilibrated-flux and
Prager--Synge/functional a posteriori estimates. FortFEM's compatible
\(H^1\) and \(H({\rm div})\) spaces are especially suitable.

A typical pattern constructs a conforming approximate state \(v\) and an
equilibrated flux \(y\) satisfying the continuous balance, then bounds the
energy error by a computable constitutive defect.

DG/broken states need a conforming reconstruction or a theorem that explicitly
accounts for jump distributions.

## Kinetic/spectral complements

Continuous kinetic certification may need bounds for:

- omitted spatial modes/elements;
- momentum-domain tails;
- radial basis truncation;
- spherical harmonics above \(L\);
- Fourier modes outside the retained band;
- collision-potential solves;
- nonlinear mode coupling into omitted complements;
- characteristic/interpolation defects;
- field reconstruction and exterior boundaries;
- arithmetic rounding.

Tail estimates are analytical/numerical bounds, not extrapolation from the
last few coefficients.

## Adjoint quantity-of-interest bounds

For a linear problem

\[
Au=b,\qquad Q(u)=q^\dagger u,
\]

the adjoint solves

\[
A^\dagger z=q.
\]

With approximate primal/adjoint states,

\[
r=b-A\bar u,\qquad r_z=q-A^\dagger\bar z,
\]

the observable correction is centered by

\[
\bar z^\dagger r
\]

and the remaining error is bounded using primal and adjoint dual norms. In the
SPD case,

\[
\left|
Q(u)-Q(\bar u)-\bar z^\dagger r
\right|
\le
\|r_z\|_{A^{-1}}\|r\|_{A^{-1}}.
\]

This is the preferred route to tight certificates on actual scientific
observables instead of unnecessarily tight full-state norms.

Nonlinear problems use the corresponding linearized adjoint plus a rigorous
nonlinear remainder/stability argument.

## Goal-oriented refinement and model defects

If a residual defect decomposes as

\[
d=\sum_k d_k,
\]

then

\[
-\lambda^\dagger d_k
\]

is a natural first-order estimate of each contribution to a scalar observable.
This is useful for deciding whether to refine spatial resolution, velocity
tails, collision fidelity or a reduced model.

It is only an indicator until the appropriate remainder/stability theorem is
established.

For a cheap model \(F_0\) and richer model

\[
F_1=F_0+\Delta F,
\]

one may evaluate the rich-model defect at a cheap-model candidate and attempt a
validated nearby-solution theorem. If it closes, a full expensive nonlinear
rich solve may be unnecessary. If it does not, no model-error certificate is
claimed.

## Onsager, positivity and structural checks

For collisional transport, retain scoped structural evidence for:

- invariant nullspaces;
- weighted self-adjointness or Onsager--Casimir transformation;
- nonnegative entropy production;
- constrained coercivity/spectral gap;
- reciprocal transport coefficients.

A discrete PSD matrix is not automatically a continuous collision spectral
gap. The function space and complement must be part of the theorem.

## Rare-event and wall observables

Loss fractions, hitting times and wall loads are often discontinuous
functionals of trajectories. A generic \(L^2\) ball does not directly bound
them.

Certification may require:

- event-location enclosures;
- transversality/grazing margins;
- boundary trace bounds;
- measure of an ambiguous set;
- validated quadrature/QMC discrepancy;
- probabilistic bounds for MC estimators.

The target alpha-loss certificate is an example:

\[
L_-(T)\le L(T)\le L_+(T),
\]

with the unresolved/grazing initial-condition measure explicitly bounded.

## Monte Carlo and QMC

For random estimators the certificate names the sampling law, independence or
dependence assumptions, estimator, variance/concentration method and
confidence level.

A statement such as

\[
\Pr(Q\in I)\ge1-\delta
\]

has different semantics from deterministic inclusion.

Deterministic QMC may instead use a discrepancy/variation bound when those
quantities can be rigorously established. Randomized QMC carries both
low-discrepancy structure and statistical uncertainty.

## Physical/model uncertainty

Keep separate fields in result metadata for:

1. numerical error within the model;
2. model-reduction/comparison error;
3. physical-data/input uncertainty;
4. stochastic sampling uncertainty.

Intervals can be combined only when a proved mapping shows they refer to the
same physical quantity and their dependence has been handled correctly.

## Formal verification

Lean is optional and downstream.

Good eventual formal targets are small generic theorems reused across many
runs:

- contraction/Krawczyk/radii-polynomial result;
- SPD residual/energy bounds;
- adjoint QoI identity/bound;
- Prager--Synge/functional estimate;
- selected model-defect theorem;
- exact finite-dimensional conservation/Onsager identities.

Do not formalize the whole HPC solver. The trusted numerical core should stay
small even without Lean: exact model definition, rigorous arithmetic,
certificate data and a small checker.

## Refinement policy

Refinement is driven by checked error contributions and observable sensitivity,
not by visual mesh convergence alone.

Possible controls include spatial \(h/p\), momentum extent/order, harmonic
cutoff, Fourier bandwidth, collision potential tolerance, nonlinear tolerance,
interval precision and model fidelity.

A larger approximation space can worsen conditioning or interval dependency, so
promotion uses measured certificate width, not resolution as a proxy.

Chris&AI
