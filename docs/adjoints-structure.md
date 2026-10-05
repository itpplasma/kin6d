# Adjoint goals and structural certificates

Status: required contracts and planned general certification routes. The
[restored laboratory](steady-laboratory.md) retains scoped continuous reference
bounds and independent local derivative checks. The [general damped-transport
case](general-transport.md) additionally certifies arbitrary represented-source
residuals, implicit tangents and a cosine-density directional derivative.
Equilibrated-flux estimators, continuous physical goal adjoints and general
Onsager certificates remain unimplemented.
This extends [differentiability](differentiability.md) and the
[continuous certificate contract](certification.md).

## Tangent and adjoint consumers

For constrained stationary equations `F(U,theta)=0`, tangent directions solve
`F_U dU = -F_theta dtheta`. A scalar goal uses an adjoint solve and the implicit
gradient described in the differentiability policy. Forward solves scale with
the number of requested independent parameter directions; adjoints scale with
the number of scalar goals. Choose by the actual products and conditioning,
with batching where useful. Neither consumer stores a full sensitivity matrix
or differentiates Newton/Krylov/continuation control flow.

[DESC Part II](https://doi.org/10.1017/S0022377823000399) derives perturbation
and continuation from equilibrium force balance and develops higher-order
perturbations. It is a precedent for equation-level derivatives, not evidence
of kin6d's matrix-free or rigorous sensitivity capabilities. Reuse primal
operator/preconditioner infrastructure where valid; an adjoint preconditioner
needs the appropriate transpose/weighted action rather than an assumed reuse
of a nonsymmetric forward action.

## Pairings and evidence-bearing metadata

An operator declares domain/range spaces, their natural pairings, derivative
domain/event margins, constraints/nullspaces and parameter dependencies.
Structural metadata records self/skew adjointness, energy/free-energy and
entropy forms, invariant transfers, positivity, time-reversal parities,
Casimir transformation, and any coercivity/inf-sup bound. Each entry states
its assumptions, norm, continuous/discrete scope and evidence/proof status.
A boolean declaration or numerical spot check is not a theorem premise.

For coefficient metrics `M_X,M_Y`, the physical adjoint of `A:X -> Y` is
`A^dagger = M_X^{-1} A^* M_Y`, where `*` is conjugate transpose. Raw VJPs act
on coefficient duals; Riesz maps turn them into physical-space adjoints.
Mass, species, momentum, phase-space and boundary weights are explicit.
Metric dependence on parameters also enters sensitivities. Flux transpose,
FCI gather/scatter, Fourier normalization and MPI reduction ownership must
respect the same pairing.

Classify only after establishing the applicable structure:

| Class | Needed stability evidence |
|---|---|
| Self-adjoint positive/coercive | Verified positive lower bound in the declared norm |
| Positive semidefinite with nullspace | Compatibility and constrained-space bound; removing a kernel alone does not prove a gap |
| Skew plus coercive dissipation | Domain/boundary adjoint identity and coercive estimate, with solvability |
| Hypocoercive/inf-sup stable | Modified norm or inf-sup and adjoint/surjectivity control |
| Indefinite/non-normal | Verified inverse/resolvent bound on the required constrained branch |

For `A=S+K`, skew `K` contributes zero real energy only with the correct
domain and boundary conditions. Open boundaries add flux terms. A coercive
real-part bound gives an inverse estimate when existence/surjectivity holds;
injectivity alone is insufficient for an arbitrary unbounded operator.

The first general case is periodic `L=nu+v*d_x`, nu>0, on a declared finite
velocity interval. In the `integral dx dv` pairing, `L^dagger=nu-v*d_x`.
The damped-characteristic resolvent supplies existence and the continuous
sup-norm inverse bound `1/nu` for both operators, including v=0. The same
constant must not be transferred to self-consistent electric coupling or FP
diffusion without a theorem for their actual domains and pairings.

## Residual and energy bounds

For a self-adjoint positive `A`, `Au=b`, and `r=b-A ubar`, the error satisfies
`Ae=r` and exactly `||e||_A^2 = <r,A^{-1}r>`. If `A >= alpha I` with a proved
`alpha>0`, then `||e|| <= ||r||/alpha` and
`||e||_A <= ||r||/sqrt(alpha)`. Specify whether these are finite algebraic
or continuous-space statements; a discrete eigenvalue bound does not control
unresolved PDE modes, quadrature, geometry or velocity tails.

For SPD `P` with proved `c1 P <= A <= c2 P`, the energy error squared lies
between `<r,P^{-1}r>/c2` and `<r,P^{-1}r>/c1`. A useful preconditioner becomes
a rigorous estimator only after proving this equivalence and enclosing its
action. Fast approximate/FP32/lagged candidate preconditioning supplies no
such constants automatically.

[Meurant–Tichý](https://doi.org/10.1007/s11075-012-9591-9) develops quadrature
error bounds from CG information. A later generic FortNum/FortSparse route
may exploit them for SPD subsolves. Certified stopping requires valid spectral
endpoints and control of finite-precision recurrence/orthogonality errors;
floating Gauss–Radau estimates are not automatically guaranteed bounds.

For a strongly convex differentiable energy, a proved Hessian lower bound
on a convex region containing the candidate and an existing minimizer yields
state and energy-gap bounds from its gradient. Constrained minimizers require
the corresponding KKT/projected residual. Local positivity without the
region, feasibility and existence proof is insufficient. Gauss–Newton
`F_U^* F_U >= 0` is optimization structure; it establishes neither physical
coercivity of `F_U` nor dynamical stability. Normal equations can also square
conditioning and are not the default kinetic solve.

## Continuous flux and variational bounds

For symmetric uniformly positive diffusion tensor `D`, a conforming `v`
and equilibrated `y in H(div)` satisfy a hypercircle estimate when
`-div y=f` and the appropriate Neumann flux is matched. With
`sigma=D grad u`, energy error plus complementary-flux error equals
`||D grad v-y||_{D^{-1}}^2`, subject to matching essential boundary data.
This controls the continuous PDE rather than only the algebraic FEM solve.

FortFEM's H1 and RT/BDM H(div) spaces are natural ingredients. Approximate
equilibration or projected load data requires a residual/data-oscillation
term, with certified domain/trace constants; do not assert exact balance of
the original continuous load from finite quadrature alone. Nonconforming
approximations need reconstruction/interface terms. Generic equilibrated
flux, majorant/minorant and feasibility machinery belongs upstream.

[Repin–Sauter](https://doi.org/10.1016/j.crma.2006.06.024) supplies functional
reaction-diffusion bounds, and [Repin–Sauter–Smolianski](https://doi.org/10.1137/050641533)
addresses mixed elliptic formulations. Their hypotheses must be instantiated
for each kin6d operator; they are not general full-kinetic error theorems.
Primal/complementary-dual feasible energies can give a certified gap.
The variational dual flux is distinct from a goal sensitivity adjoint.

The next bounded FP case uses prescribed constant theta>0 and positive smooth
weight w on `[-P,P]`:

```
A f = -d_v(theta*w*d_v(f/w)) = S,
theta*w*d_v(f/w) = 0 at both endpoints,
integral f dv = prescribed mass,  integral S dv = 0.
```

For the OU operator, w is proportional to `exp(-v^2/(2*theta))`; it is a
coefficient of the weighted formulation, not an evaluated forced solution.
The physical state pairing is `integral f*h/w dv`. The nullspace is the span
of w, so compatibility and the mass constraint must be checked explicitly.

Write `gbar=fbar/w` with gbar in H1 and construct a flux y satisfying
`-d_v y=S` and zero endpoint flux. On the mass-constrained space the energy
error is bounded by the constitutive defect

```
eta^2 = integral |y-theta*w*gbar'|^2/(theta*w) dv.
```

Testing the weak error equation and applying weighted Cauchy--Schwarz proves
this estimate. It covers derivative jumps of a continuous piecewise-linear
candidate without interpreting its cellwise strong residual as a global
L2 function. If y balances only a projected load or approximate endpoint
flux, add independently bounded data/boundary defects. Reconstruction mass
error requires a separate nullspace correction or bound.
An arbitrary compatible signed source does not guarantee a nonnegative
distribution; physical admissibility requires an additional checked premise.

A conservative continuous weighted Poincare lower bound is
`alpha = theta*(w_min/w_max)*(pi/(2*P))^2` on `integral w*g dv=0`.
To justify it, compare the weighted variance with the variance about the
ordinary mean, apply the interval Neumann Poincare inequality, and use
`w_min<=w<=w_max`. Enclose these constants and discharge weak existence before
applying `||f-fbar||_L2(1/w) <= eta/sqrt(alpha)`. This analytical route still
requires implementation and coefficient/compatibility verification; a
discrete eigenvalue or removal of the mass mode alone supplies no such bound.

For an exactly mass-constrained error in one dimension, weighted mean zero
gives a zero of its continuous H1 representative. Integration and
Cauchy--Schwarz then give the pointwise bridge
`|f-fbar|(v) <= w(v)*sqrt(integral 1/(theta*w) dv)*eta`.
Include mass/reconstruction defects before using this as a plotting band.
Higher-dimensional pointwise claims need their own regularity theorem.

## Certified goal correction

For bounded linear goal `Q(u)=<q,u>` and adjoint `A^dagger z=q`,
`Q(u)-Q(ubar)=<z,r>`. For SPD `A`, an approximate adjoint `zbar` with residual
`r_z=q-A zbar` gives
`Q(u)-Q(ubar)=<zbar,r>+<r_z,A^{-1}r>`.
If certified primal/adjoint dual-energy errors are bounded by `eta_u,eta_z`,
then `Q(u)` lies in the rigorous enclosure of
`Q(ubar)+<zbar,r>+[-eta_u eta_z,eta_u eta_z]`.
Enclose the correction, observable evaluation and all continuous error terms.

Apply these enclosures to every requested observable and gradient, with the
per-quantity scales and budgets in [certification](certification.md).
The workhorse FP64 center plus proved aggregate arithmetic radius and the
strict outward/high-precision evaluator have the same containment obligation.
Propagate rounding through inverse, adjoint and goal norms before comparing it
with its allowance; presumed percentages or an unbounded adjoint correction
cannot close a goal certificate. Cache primal/adjoint certificate preparation
for immutable inputs and report its cost separately from rendering.

Nonsymmetric/indefinite systems need corresponding primal/adjoint norms and
stability; the SPD formula's norm bound cannot be transferred by analogy.
Nonlinear goals/residuals require a controlled remainder. Point hits, loss
indicators and unbounded point evaluations need regularity/event analysis.
[Becker–Rannacher](https://doi.org/10.1017/S0962492901000010) and
[goal-oriented FEM](https://doi.org/10.1016/S0898-1221(00)00317-5) motivate
adjoint-weighted refinement. Estimated weights and omitted remainders remain
indicators until rigorous upper/lower bounds have been established.

## Model defects and adaptive fidelity

A higher-model defect is evaluated in a common space. Derive a lifting `E`
of the cheap state and observable/boundary correspondence, then evaluate
`d=F_high(E U_low)`. FK and GK, or fluid and kinetic neutrals, do not have
the same state variables; the raw expression `F_FK-F_GK` is undefined without
this bridge. Compatible collision actions may share a space but still need
consistent backgrounds, perturbations, units and invariants.

An adjoint-weighted defect gives a first-order correction under a justified
linearization. Turning it into a higher-model interval requires a verified
higher-model inverse/stability bound near the lifted state, nonlinear
remainder, lifting defect and goal correspondence. Stability of the cheap
Jacobian alone does not prove a nearby higher-model solution. The cheap solve
may already supply a usable higher-model candidate if the certificate closes;
otherwise refine the model or perform the expensive solve.

Separate algebraic, discretization/tail, model-reduction and uncertain-input
contributions, retaining the proved comparison chain. Local signed adjoint
weights may guide h/p, radial/harmonic/Fourier, collision, neutral or atomic
refinement. Rank predicted goal-error reduction by measured added cost, then
recheck after the change. Cancellation-prone estimates and a successful
controller choice do not by themselves certify the final observable.

## Collision and transport reciprocity

Use `C_lin` for a dissipative collision generator and `D=-C_lin >= 0`.
A steady dissipative transport form is commonly `A=T-C_lin=T+D`, with
`T^dagger=-T` only under its phase-space/boundary contract. Then
`A^dagger=-T+D`. Keep signs explicit; the physical collision generator is
not a positive diffusion matrix under the same convention.

Full linearized Landau at common equilibrium has weighted entropy/adjoint
structure and collision invariants. It is not legitimate to infer a uniform
unweighted L2 spectral gap merely by removing invariants for the Coulomb
whole-velocity problem. Weighted coercivity/tails or hypocoercive analysis
is needed; a finite radial truncation does not establish the infinite-domain
bound. See the [Coulomb Landau analysis](https://www.ceremade.dauphine.fr/~mischler/expo/2016Mittag.pdf).

For declared conjugate forces/fluxes, entropy production depends on
`L_sym=(L+L^T)/2`. Certify its nonnegativity in the stated real force space
using verified eigenvalue/LDL or constrained semidefinite evidence. An
antisymmetric part does not produce entropy and can be physically legitimate.
Reciprocity instead compares `L_ij(B)=epsilon_i epsilon_j L_ji(-B)` under the
complete time-reversal/background transformation. Derive signs and force/flux
normalizations; do not impose `L=L^T` at fixed B without a matching theorem.

[Sugama–Horton](https://doi.org/10.1063/1.871856), published in 1996,
connects neoclassical entropy production and Onsager symmetry in general
toroidal multispecies plasmas. [NEO-2](https://doi.org/10.1585/pfr.3.S1061)
provides independent transport-coefficient and collision-model comparisons.
Retain separate source/response solves for reciprocal coefficients rather
than filling half the matrix by symmetry or projecting a failed result.
Use a scaled absolute/relative defect safe at zero coefficients; a single
relative ratio with a potentially zero denominator is not a universal gate.

Krook's fixed-target derivative `-nu I` with prescribed real nonnegative
`nu` can be self-adjoint/dissipative; this alone implies no physical pair
conservation or neoclassical reciprocity. Nonlinear BGK/LBD/Landau/BB require
their actual equilibrium linearization/bracket, not an unconditional linear
adjoint flag. Test-particle operators exchange with reservoirs. Gyroaveraging,
finite domains, boundary conditions and discrete corrections may alter these
properties. Record the verified construction for each collision configuration.

## Acceptance gates

Retain independent analytical/invariant and cross-code tests, weighted JVP/VJP
duality, invariant complements, dissipation, independent reciprocal solves,
analytic transport limits and NEO-2/SFINCS comparisons. Structural diagnostics
do not replace a continuous interval or prove adequacy of the collision model.
Poisson/diffusion/mass/radial blocks may admit coercive certification after
their constraints; Grad–Shafranov Newton and stable ideal-MHD blocks require
actual positivity evidence. Near marginality an inverse bound may diverge;
certificate failure alone does not prove a physical instability.

Chris&AI
