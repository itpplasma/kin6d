# Adaptive multiscale, multifidelity and uncertainty architecture

Status: canonical long-range vision. This document defines the common
error-driven adaptation policy across numerical resolution, physical model,
representation, correlation order, quantum/classical fidelity, multiphysics
coupling and uncertainty. It is not an instruction to implement a universal
controller before concrete vertical slices exist.

## Objective

KIN6D should eventually solve only as much physics, correlation, quantum
mechanics, numerical resolution and stochastic detail as the requested
observable actually requires.

For observable \(Q\) with requested tolerance \(\varepsilon_Q\), the long-range
objective is

\[
\boxed{
\min_{\mathcal C} \; {\rm cost}(\mathcal C)
\quad\text{subject to}\quad
{\rm error}_Q(\mathcal C)\le\varepsilon_Q ,
}
\]

where configuration \(\mathcal C\) includes **local and temporal** choices of

- physical model;
- representation;
- spatial/velocity/angular resolution;
- time resolution/integrator;
- collision/correlation fidelity;
- classical/quantum fidelity;
- radiation/reaction fidelity;
- micro/macro coupling level;
- deterministic/stochastic estimator;
- arithmetic precision and nonlinear/linear solver tolerances.

The smallest physical scale must not automatically dictate the global mesh,
global timestep or global model.

## Adaptation is local in more than physical space

An adaptive region may be a subset of

- configuration space \(x\);
- phase space \((x,p)\);
- species;
- energy/pitch/action space;
- correlation cluster/order;
- radiation frequency/angle;
- reaction/internal-state graph;
- time interval/event window;
- parameter/uncertainty space.

A refinement can therefore be local in **space, time, phase space, species,
physics and stochastic sample space simultaneously**.

Examples include:

- MHD globally with a transient GK/FK resonant patch;
- fluid neutrals except for a rarefied kinetic boundary layer;
- classical Maxwell--Lorentz particles except for isolated quantum-diffraction
  encounters;
- \(f_1+g_2\) globally with \(g_3\) activated only where its observable defect
  matters;
- grey radiation except for selected multigroup/angle-resolved regions and
  times;
- sparse expensive turbulence microproblems only where the macro adjoint
  indicates high leverage.

## Error ledger

Do not collapse all uncertainty into one Gaussian standard deviation.

For an observable, maintain conceptually distinct contributions such as

\[
\Delta Q =
\Delta Q_{\rm discretization}
+\Delta Q_{\rm time}
+\Delta Q_{\rm solver}
+\Delta Q_{\rm arithmetic}
+\Delta Q_{\rm model}
+\Delta Q_{\rm interface}
+\Delta Q_{\rm correlation}
+\Delta Q_{\rm quantum}
+\Delta Q_{\rm surrogate}
+\Delta Q_{\rm input}
+\Delta Q_{\rm reconstruction}
+\Delta Q_{\rm stochastic}.
\]

This sum is bookkeeping, not a rule for arithmetic combination. Each component
retains its semantics and dependencies.

### Deterministic bounded uncertainty

Use rigorous intervals/balls, residual estimates, asymptotic remainder bounds,
validated interpolation/transfer bounds and model-defect theorems where
available.

### Probabilistic uncertainty

Use probability only when there is a declared stochastic law: MC/RQMC sampling,
experimental/posterior uncertainty, random inputs or an explicitly chosen
probabilistic discrepancy model.

A statement

\[
\Pr(Q\in I)\ge1-\delta
\]

is not interchangeable with deterministic inclusion.

### Epistemic/model uncertainty

If no rigorous model-error theorem exists, prefer in this order:

1. analytical asymptotic/remainder bound;
2. a posteriori rich-model defect plus stability/remainder estimate;
3. verified/enclosed difference against sparse richer-model calculations;
4. empirical cross-model envelope with explicit non-rigorous status;
5. probabilistic discrepancy model only when it is a useful explicit model of
   ignorance.

Do not manufacture randomness merely because a deterministic missing term is
unknown.

## Rich-model defect and adjoint weighting

For a cheap model

\[
R_0(U_0)=0
\]

and a richer model \(R_1\), evaluate the rich residual at a lifted cheap
solution:

\[
d_{0\to1}=R_1(LU_0).
\]

For observable \(Q\), solve the richer-model adjoint approximately or
rigorously enough for the intended claim,

\[
R_{1,U}^{\dagger}\lambda=Q_U^\dagger.
\]

Then

\[
\eta_{0\to1}
=
-\lambda^\dagger d_{0\to1}
\]

is the leading observable-specific model correction. With a checked stability
constant and nonlinear remainder this can become a certified model-error bound.

This pattern should apply to

- MHD \(\leftrightarrow\) kinetic-MHD/GK/FK;
- DK/GK \(\leftrightarrow\) FK;
- quasineutral/Darwin \(\leftrightarrow\) full Maxwell;
- Landau/BB \(\leftrightarrow\) Balescu/mean-force/correlation models;
- \(g_2\) truncation \(\leftrightarrow\) \(g_3\);
- classical \(\leftrightarrow\) quantum correction;
- fluid \(\leftrightarrow\) kinetic neutrals;
- diffusion/moment radiation \(\leftrightarrow\) angular transport;
- reduced \(\leftrightarrow\) state-resolved chemistry.

A cheap model may therefore remain the state carrier while a richer model is
used only as an error oracle.

## Analytical asymptotics and numerical enrichment

When a derivation provides

\[
R_\epsilon
=
R^{(0)}
+\epsilon R^{(1)}
+\cdots+
\epsilon^p R^{(p)}
+\mathcal R_{p+1},
\]

retain the actual asymptotic parameter and, where possible, a bound

\[
\|\mathcal R_{p+1}\|
\le C\epsilon^{p+1}.
\]

Combine analytical estimates with sparse numerical enriched solves. Neither is
universally superior:

- analysis can prove scaling and safe limits far outside directly affordable
  simulations;
- enriched numerics can estimate constants, detect breakdown and cover
  geometries/closures where a theorem is unavailable.

The controller may use both.

## Correlation-order adaptivity

The classical and quantum BBGKY/cumulant hierarchies are natural adaptive axes.

With

\[
f_2=f_1f_1+g_2
\]

and higher cumulants \(g_3,\ldots\), solve a retained hierarchy through order
\(k\) and evaluate the residual of the next omitted level.

For example, if \(g_3=0\) in the retained model, define an omitted-correlation
defect \(d_{g_3}\). Its observable impact is estimated by

\[
\eta_{g_3}=|\lambda^\dagger d_{g_3}|.
\]

If this is negligible relative to the requested tolerance, do not solve
\(g_3\). If it dominates, refine the correlation model locally in
space/phase-space/time or launch an embedded microscopic many-particle
calculation.

Quantum reduced-density-matrix/cumulant hierarchies use the same logic, with
fermionic symmetry and quantum observables respected.

## Classical--quantum local promotion

A hot nondegenerate plasma may be classical almost everywhere while a rare
close encounter needs quantum diffraction.

A local encounter indicator can depend on quantities such as

\[
\chi_{\rm diff}\sim\lambda_{\rm dB}/b
\]

or an equivalent scattering criterion. When a close pair enters the quantum
regime, an embedded two-body quantum scattering problem may replace the
classical close interaction and return an asymptotic outgoing state or
transition kernel.

This is a **local event refinement**, not a declaration that the whole plasma
is quantum.

Degeneracy/exchange is different. When

\[
\Theta=T/T_F\lesssim O(1),
\]

Pauli blocking and exchange depend on the occupied many-body state and cannot
generally be handled by promoting only one close pair. The promoted object may
need to be a region/species/model such as Fermi/Wigner/HF/quantum kinetic.

Demotion back to a classical representation requires a declared coarse-graining
or decoherence/phase-information criterion. Do not discard persistent quantum
coherence by fiat.

## Microscopic near/far adaptivity

For the Maxwell--Lorentz branch, split interactions into controlled near and
far pieces,

\[
G=G_{\rm near}+G_{\rm far}.
\]

Possible hierarchy:

- far smooth interaction: FMM/tree, Ewald/P3M or Maxwell field solve;
- near classical interaction: direct accurate pair dynamics;
- near quantum interaction: quantum scattering microproblem when triggered.

Overlap subtraction/matching must avoid double counting and carry an interface
error estimate.

This is how sub-Debye or quantum encounter scales can remain local rather than
forcing the global electromagnetic mesh to resolve them.

## Spatial and temporal model promotion

Every fidelity axis should eventually support **promotion and demotion in
space and time** where the governing assumptions permit it.

A local model map may look schematically like

\[
{\rm MHD}
\rightarrow
{\rm GK}
\rightarrow
{\rm FK}
\rightarrow
{\rm Maxwell\!-\!Lorentz}
\rightarrow
{\rm quantum\ microproblem},
\]

but not every neighboring pair is valid without a derived interface.

Temporal adaptivity is equally important:

- a region may require FK only during a resonance or crash;
- a quantum encounter may last femtoseconds inside a much longer classical
  trajectory;
- a radiation burst may temporarily require higher spectral fidelity;
- a fast microproblem may be subcycled while macro transport takes long steps.

Use local timestepping, subcycling, multirate or event-driven integration when
their coupling error is controlled.

Promotion/demotion should use hysteresis or comparable safeguards to prevent
model thrashing.

## State transfer and interface obligations

A model switch is not just a change of array type.

Every promotion/demotion interface declares

- lift/projection/reconstruction maps;
- conserved particle/nuclei/charge/momentum/energy quantities;
- electromagnetic constraints/gauges;
- phase-space measure;
- entropy/free-energy information when meaningful;
- lost unresolved information;
- transfer/interface defect;
- derivative/JVP/VJP semantics;
- certificate/UQ semantics.

If a reduced state does not determine a unique richer state, promotion needs a
declared maximum-entropy, conditional, sampled or otherwise modeled
reconstruction. That reconstruction uncertainty remains explicit.

## Heterogeneous multiscale method pattern

A macro state can request only the microscopic information it needs:

\[
U_{\rm macro}
\longmapsto
\mathcal M_{\rm micro}[U_{\rm macro}]
\longmapsto
\{\Gamma,Q,\Pi,\epsilon,\ldots\}.
\]

The microprovider may be

- neoclassical DK;
- local/global GK/FK;
- Maxwell--Lorentz \(N\)-body;
- Wigner/HF/QMC;
- detailed radiation transport;
- state-resolved chemistry.

Evaluate expensive microproblems only at selected space-time points. Reuse and
interpolate previous microresults when a verified validity key and error
envelope allow it.

If locality or scale separation fails, enlarge the promoted region or move
toward fully coupled evolution. The controller must be allowed to conclude
that the cheap multiscale assumption is invalid.

## Surrogates and caches

A surrogate is an accelerator, never an unlabelled replacement for the target
model.

Every surrogate/cache entry records

- training/support domain;
- state/parameter coordinates;
- model and numerical fidelity;
- derivative support;
- error/validity estimate;
- provenance/version;
- whether extrapolation is occurring.

An adjoint or model-defect indicator can request a fresh expensive sample when
surrogate uncertainty matters to \(Q\).

## Multilevel/multifidelity UQ

Use the hierarchy itself to reduce UQ cost.

For nested fidelities,

\[
Q_L
=
Q_0
+
(Q_1-Q_0)
+\cdots+
(Q_L-Q_{L-1}),
\]

sample cheap levels frequently and expensive corrections sparsely. This
includes control-variate, multilevel Monte Carlo and multifidelity/QMC
patterns.

Examples include

- coarse/fine mesh;
- MHD/GK/FK;
- continuum/\(N\)-body;
- classical/quantum correction;
- low/high correlation order;
- low/high radiation fidelity.

Allocation of samples and model calls is driven by variance/error reduction per
cost, while deterministic model and discretization bounds remain separate from
sampling variance.

## Unified refinement candidates

At a point in the solve, the candidate set may include

\[
\mathcal K=
\{
h,p,\Delta t,
N_p,L_{\rm harmonic},
p_{\rm velocity},
\text{solver tolerance},
\text{arithmetic precision},
\text{collision model},
\text{correlation order},
\text{classical/quantum promotion},
\text{MHD/GK/FK promotion},
\text{neutral model},
\text{radiation model},
\text{chemistry fidelity},
\text{microproblem sample},
\text{surrogate refresh}
\}.
\]

For each candidate \(k\), estimate an observable benefit \(\eta_k\) and
incremental cost \(C_k\). A generic greedy policy is

\[
k^\star
=
\arg\max_k
\frac{\eta_k}{C_k},
\]

subject to compatibility, stability, memory and coupling constraints.

This score is a planning/controller abstraction, not a theorem. Competing
refinements can interact nonlinearly, so periodically recompute indicators
after promotion.

## Termination criterion

The adaptive calculation stops refining when the requested observable error
budget closes with its declared semantics.

A result should record not only \(Q\) and its uncertainty but also the final
space-time/phase-space model map and dominant residual uncertainty
contributions.

If no affordable configuration closes the tolerance, return the best candidate
with an explicit unmet-tolerance/certification-failed status; do not silently
relax the user's requested accuracy.

## Initial development strategy

Do **not** begin by implementing the universal controller.

Grow the architecture through concrete vertical slices:

1. ordinary \(h/p/\Delta t\) refinement with certified/estimated observable
   contribution;
2. one cheap/rich model defect, e.g. simple collision or field fidelity;
3. one sparse microprovider invoked only where an adjoint indicator asks for it;
4. one local-in-time promotion/demotion example;
5. one local spatial model-patch example with conservative transfer;
6. one correlation-order \(g_2\to g_3\) indicator;
7. one classical close-encounter \(\to\) quantum microproblem;
8. one multilevel/multifidelity UQ example.

Only then extract a generic refinement scheduler from repeated real patterns.

Chris&AI
