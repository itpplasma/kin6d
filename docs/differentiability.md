# Differentiability and adjoint policy

Status: canonical derivative policy. The [general damped-transport case](general-transport.md)
has generated primal/JVP/VJP leaves, implicit source/damping tangents,
continuous tangent bands and a certified cosine-density directional derivative.
Its nodal VJP uses a Euclidean coefficient pairing. Continuous physical goal
adjoints and general nonlinear certified gradients remain implementation gates.
Earlier analytical checks are retained in [the laboratory evidence](steady-laboratory.md).

## Principle

KIN6D differentiates mathematical operators and implicit equations, not the
incidental algorithms used to solve them.

Every continuous production operator should expose, where its domain is smooth,

\[
F,\qquad DF\,v,\qquad DF^\dagger w.
\]

Hessian-vector products are added only for a concrete consumer.

FortSym is the primary source of project-owned primal/JVP/VJP expressions.
FortAD, analytical formulas, finite differences and complex-step checks are
independent verification tools. LLVM/Enzyme is optional and must never limit the
supported Fortran compiler set.

## Symbolic source and generated families

A local symbolic definition should be able to emit a family such as

\[
\{K,\ K_{\rm JVP},\ K_{\rm VJP},\ K_{\rm HVP},
  \ K_{\rm strict},\ K_{\rm enclosure},\ K_{\rm CUDA}\}.
\]

Not every member must exist from day one, but no handwritten derivative fork is
introduced merely because an emitter is missing. Add the generic capability
upstream or leave the derivative target explicitly unavailable.

Symbolic common-subexpression elimination is cost-aware. Fewer FLOPs can lose
to additional live temporaries/register spills. Hardware scheduling therefore
remains separate from mathematical simplification.

## Indexed reductions and bounded memory

Derivative generation must preserve lazy/tiled operator structure.

For

\[
a_i=\sum_jF(x_i,y_j),
\]

the JVP and VJP are themselves reductions. Do not materialize the pairwise
matrix or Jacobian.

The VJP should reverse the reduction direction as needed and recompute cheap
tile-local intermediates. The default is **no global derivative tape**.

This policy applies to:

- pairwise kernels;
- PIF/particle-in-basis;
- FEM quadrature;
- sparse angular couplings;
- reaction graphs;
- collision reductions.

## Implicit function theorem

For a steady/equilibrium state

\[
R(U,\theta)=0,
\]

the tangent solves

\[
R_U\,\delta U=-R_\theta\,\delta\theta.
\]

For a scalar objective \(Q(U,\theta)\), solve the adjoint

\[
R_U^\dagger\lambda=Q_U^\dagger,
\]

then

\[
\frac{dQ}{d\theta}
=
Q_\theta-\lambda^\dagger R_\theta.
\]

This is the default for equilibrium-constrained optimization and mirrors the
useful DESC strategy while remaining compiler-independent and matrix-free.

When the number of parameter directions is small, forward tangents/JVPs are
natural. When the number of scalar objectives is small relative to parameters,
the adjoint/VJP formulation is preferred.

Requested gradients are scientific outputs with their own units, scales,
absolute floors, interval radii and status. Apply the observable accuracy and
arithmetic budgets in [certification](certification.md#accuracy-targets-and-arithmetic-budgets)
to each requested component. Report components whose stability or event
premises fail; do not omit them or present an uncertified direction as a
certified gradient. Statistical sensitivities retain separate probability and
sampling-error statements.
## External inference consumers

TIAGO is the reconstruction/inference consumer for experimental data. KIN6D
must expose enough of the implicit stationary map that TIAGO can compute
diagnostic derivatives without differentiating through the solver algorithm.
The preferred contract is matrix-free: evaluate the forward state/observable,
apply JVPs, and apply VJPs/adjoints in the declared physical pairings.

For a diagnostic map H(U,eta), TIAGO may compose its own reverse-mode
diagnostic derivative with KIN6D's implicit adjoint. KIN6D does not own the
likelihood, measurement covariance, calibration nuisance parameters, priors or
posterior covariance. The existing TIAGO/VMEC++ exact-adjoint path remains an
independent reference implementation while a future KIN6D provider is added.

If several distinct equilibria are supplied to TIAGO, derivatives/covariances
are local to each one. A finite mixture of local Laplace approximations is an
inference-layer construction; no derivative of an artificial average
equilibrium is required.

## Linear solves, roots and eigenproblems

Do not differentiate through iterative solvers when the exact implicit rule is
known.

For

\[
Ax=b,
\]

\[
A\,dx=db-dA\,x.
\]

The reverse solve uses \(A^\dagger\).

For a simple root

\[
g(z,\theta)=0,
\]

\[
dz=-g_z^{-1}g_\theta\,d\theta.
\]

A derivative certificate requires a transversality bound
\(|g_z|\ge\gamma>0\). At root mergers/grazing events, the exact derivative may
be undefined or unbounded; do not hide this with regularization.

Eigenvalue/eigenvector derivatives require spectral isolation, normalization
and treatment of degeneracies.

## Physical adjoints

The adjoint depends on the physical/discrete pairing. If

\[
\langle x,y\rangle_M=x^\dagger My,
\]

then

\[
A^\dagger=M^{-1}A^*M.
\]

Do not use raw Euclidean transpose as a universal definition.

Operator metadata should declare:

- natural pairing;
- self-adjoint or skew-adjoint structure;
- nullspace/complement;
- time-reversal parity;
- positivity/coercivity status.

This connects derivative correctness to energy principles and
Onsager--Casimir reciprocity.

## DESC-style equilibrium projection

When a design parameter changes an equilibrium constraint \(F(x,c)=0\), the
design derivative contains

\[
G_c-G_xF_x^{-1}F_c.
\]

KIN6D should realize this with matrix-free JVP/VJP and linearized solves rather
than forming full Jacobians or sensitivity matrices.

A re-solve/projection after each design step is a solver policy. The derivative
of the equilibrium manifold remains the IFT relation above.

## Delta-f and control variates

For

\[
f=F_0(\theta)+\delta f,
\]

differentiate both pieces consistently. A changing reference state contributes
to the tangent/adjoint even if only \(\delta f\) is represented dynamically.

Control variates used only for variance/preconditioning do not change the
target derivative. Record whether \(F_0\) is a physical model component or an
estimator aid.

## Particles, PIC, PIF and blobs

Deterministic particle and blob trajectories are ordinary differentiable
finite-dimensional systems away from events.

PIF is particularly convenient because

\[
\rho_k=\sum_pq_pw_p\widehat S(k)e^{-ikX_p}
\]

has analytic coordinate/weight derivatives that FortSym can emit directly.

For PIC, deposition/gather and field solves must use mutually adjoint
definitions in the declared pairing if that property is claimed.

Particle sorting/binning changes execution order but need not change the
mathematical map. Discrete cell-crossing/event logic may be nonsmooth and must
declare its derivative domain.

## Stochastic methods

Random algorithms have separate derivative semantics. Possible choices are:

- pathwise/reparameterized derivatives with fixed base randomness;
- common-random-number directional comparisons;
- score/likelihood-ratio estimators;
- finite-difference statistical estimators.

These carry variance and may be poor optimization directions. When a
deterministic/control-variate model exists, use its implicit adjoint for design
and treat stochastic gradients/validation separately unless the scientific
target itself is statistical.

## Transient adjoints

For long trajectories, reverse sensitivities require a deliberate
checkpoint/recompute strategy. Do not store the entire trajectory by default.

The plan records:

- state checkpoint memory;
- recomputation factor;
- event/topology handling;
- adjoint stability;
- observable time support.

A chaotic/turbulent long-time average may need shadowing/statistical response
methods rather than a naive terminal adjoint.

## Topology and events

Examples of derivative-event surfaces include:

- particle-wall impact;
- trapped/passing transition;
- separatrix crossing;
- magnetic topology change;
- adaptive mesh/model switch;
- collision/reaction event selection.

A differentiability report states a margin to these events or declares the
gradient unavailable/non-smooth.

Optimization may use a deliberately smoothed surrogate, but the surrogate is a
different model and must be labeled.

## Verification

Every derivative family has at least one independent check:

- analytic derivative;
- centered directional difference;
- complex step where valid;
- JVP/VJP pairing identity;
- FortAD/Enzyme comparison where available.

For implicit sensitivities also verify the linearized residual equation.

Check the weighted JVP/VJP identity using the declared physical pairings and
independently evaluated metric/normalization factors. Euclidean leaf duals,
continuous adjoints and the discrete Riesz maps have different interfaces.
Verify source, boundary, constraint and metric derivatives as well as the
interior operator. Integer modes, branch switches and adaptive events are
declared discrete choices; their algebraic partial derivatives do not establish
a physical derivative across the event.

Performance records include primal, JVP and VJP time/memory for material
operators. A fast primal with an unusable VJP is not automatically the winning
schedule for optimization workloads.

For an implicit consumer, also measure tangent/adjoint solves, rigorous
derivative evaluation, certificate preparation and postprocessing. Compare
total time and peak memory to the requested gradient radius, including primal
error amplification and necessary checkpoint/recompute work. Reuse immutable
certificate data across requested components when their theorem inputs permit.

## Certified sensitivities

A generated derivative is not automatically a certified derivative of the
continuous solution.

A rigorous sensitivity/gradient claim additionally needs:

- continuous differentiability on the verified branch;
- inverse/stability bound for \(R_U\);
- rigorous \(R_\theta,Q_U,Q_\theta\) evaluation;
- adjoint/tangent solve error;
- event/topology margin;
- discretization/tail bounds.

The first general consumer uses the arbitrary-source transport equation
`L f=g`, `L=nu+v*d_x`, nu>0, on a fixed periodic spatial and bounded velocity
domain. For source/damping parameter directions its tangent obeys
`L df=dg-dnu*f`. If the certified primal radius is r_f and a reconstructed
tangent has continuous residual `r_t=dg-dnu*fbar-L dfbar`, then

```
||df-dfbar||_infinity <= (||r_t||_infinity + |dnu|*r_f)/nu.
```

Rigorous source/parameter derivative evaluation errors enter r_t. This bound
uses the same proved continuous inverse as the primal; it does not
differentiate the numerical solver. With the unweighted physical pairing
`integral dx dv`, the adjoint is `nu-v*d_x` under periodic boundaries. Its
continuous residual and goal-evaluation errors need their own enclosures.
Changing the domain, chart, velocity measure or pairing adds the corresponding
derivative terms.

For nonlinear states, primal uncertainty also changes the tangent/adjoint
operator and parameter/goal derivatives. Bound those variations on the
certified branch, together with solve error and nonlinear remainders; a small
linearized residual at the candidate alone is insufficient. Record event
margins and explicitly unavailable derivatives. Interval gradients for the
supported smooth cases are required outputs when requested; Lean assurance
is optional and does not replace the Fortran containment argument.

Chris&AI
