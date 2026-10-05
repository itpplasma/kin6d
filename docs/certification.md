# Certification specification

Status: general certification policy. The first 1X1V example already has a
reviewed continuous snapshot enclosure evaluated with interval arithmetic.
General nonlinear/kinetic certificates remain research work.

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

## What is certified

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
cover model reduction, uncertain atomic data, experimental calibration or
physical adequacy.

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
Project-owned rigorous numerical leaves come from the same retained FortSym
mathematical definitions as the ordinary leaves where supported.

Sharing a symbolic DAG prevents formula drift but is not itself a proof that
the formula is the intended physics. Independent mathematical/physical oracles
remain mandatory.

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

For collisional transport, certify or rigorously test where possible:

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
