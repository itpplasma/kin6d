# Error control and performance

Status: measured policy and implemented Q0 slice, 2026-10-06. The audience is
KIN6D implementers choosing the cheapest defensible statement for a requested
observable. The [measurement ledger](../benchmarks/error-control/baseline.json)
and [raw samples](../benchmarks/error-control/baseline.csv) cover existing
case paths. The [instrumented follow-up](../benchmarks/error-control/measured.json)
covers 23 paths; [mode records](../benchmarks/error-control/modes.json) cover
the new modes. Planned tiers in the table
below do not change the current claim of an unimplemented checker.

The fast solver remains untrusted. Every reported observable or requested
derivative needs a named error statement or explicit unavailable/failed status.
Numerical, model, input/reconstruction and stochastic contributions remain
separate. An inexpensive indicator can drive refinement while a final check
decides whether the requested claim closes.

## Two independent relaxations

**Loosening** retains deterministic containment, `e_true <= eta`. A coarse
enclosure, crude stability constant or goal-only check can reduce cost while
overestimating the error. Report the efficiency index `eta/e_true` against an
independent reference where available; about ten is a useful target, not a
proof premise or universal cap. An efficiency index is unavailable on a
production problem with unknown true error. Zero true error gives no finite
index; report that explicitly. Validation cases must not hide unfavorable
times, norms or parameters by reporting only their best index.

**Error in the bound** allows underestimation: `e_true <= (1+epsilon_b)*eta`.
Report `approximate-bound`, the declared epsilon_b (default 0.05), its
justification, and worst observed `e_true/eta` and audit-radius/eta ratio.
Evaluation error includes omitted rounding, integration remainders, sampled
suprema, linearization or surrogate residual terms. A percentage chosen
without analysis or scoped audit evidence is an indicator, not Tier A.
An approximate bound is never displayed as an uninflated certificate.
Empirical validation applies to its documented domain and cannot prove a
large-system theorem by calibration on tiny systems.

These ingredients can coexist: a pessimistic mathematical bound may also be
evaluated approximately. Report its loosening and epsilon_b separately.
Probability is a third independent ingredient, not a way to hide unknown
systematic error.

For each observable use `tau = atol + 0.01*scale` by default. Tier D/S requires
`eta <= tau`; Tier A requires `(1+epsilon_b)*eta <= tau`; Tier P requires its
confidence radius to fit tau. A tenfold loose bound that exceeds tau fails the
request. Initially reserve 1% of tau for bound-evaluation error after stability
and observable amplification. Precision and refinement follow this budget.

## Tiers and claim strength

| Tier | Work and target | Reported status | Acceptance and escalation |
|---|---|---|---|
| E | Already computed RK defects, floating residuals, embedded quadrature, cutoff or goal indicators; incremental overhead target 10–20% or less | `indicator`, with norm, time/cell and missing terms | Drive adaptivity; if this target fails, sample at checkpoints or replace the estimator; cannot close a certificate |
| P | Randomized a posteriori norm/trace/integral estimator with a proved tail inequality | `probabilistically-certified`, confidence `1-delta`, named quantity and algorithmic law | Include numerical evaluation and stability/integration terms; budget repeated/adaptive failure probabilities; escalate if interval is too wide |
| D | Deterministic workhorse: FP centers plus proved aggregate budgets, coarse enclosures, final/checkpoint residual audits, reported goals only where theorem permits | `certified-discrete` or `certified-continuous`, exact model and norm | Target total/solve about 2–5 for large evolution runs; arithmetic allowance and total radius must close; use S or refine when premises fail |
| A | Bound formula with controlled approximate evaluation; default epsilon_b=5% | `approximate-bound`, epsilon_b, analytical/empirical support and validation domain | Fits inflated observable budget; failed epsilon_b gate means `certification-failed`; use D/S for production promotion requiring a certificate |
| S | Outward arithmetic throughout checker and strict residual action | Same deterministic status vocabulary as D; `strict audit` is an arithmetic mode | Tiny references/audits can tolerate about 100×; substantially larger multipliers need an explicit cost record; avoid as routine large-run schedule |

Targets are scheduling/design goals, not measured universal guarantees. The
status names the established mathematical claim; the tier names its
evaluation strategy. `candidate`, `certification-failed`, reference comparison
and statistical interval remain useful distinct states. H1 batch-means/CLT
intervals are statistical with mixing assumptions, not deterministic bounds or
automatically rigorous finite-sample Tier P statements.

## Probabilistic estimators for deterministic computations

The probability space can be the estimator's own randomness while the
Hamiltonian, input and candidate are fixed. Gaussian sketches can use a
chi-square lower-tail bound for a norm; Rademacher quadratic forms can use
Hanson–Wright bounds; Hutchinson trace estimators need their actual matrix
assumptions and concentration inequality. Randomized QMC needs a declared
randomization and an applicable concentration/error theorem. These are policy
options, not claims that all implementations exist in KIN6D.

The implemented FortNum sampler draws iid uniform coordinates of a fixed
vector, given a proved envelope `M >= max_i |r_i|`. With `m` samples,
`X_j=(r_(I_j)/M)^2` lies in [0,1]. A one-sided
[Hoeffding bound](https://doi.org/10.1080/01621459.1963.10500830) gives

```
||r||_2 <= M sqrt(N min(1, mean(X) + sqrt(log(1/delta)/(2m))))
```

with probability at least `1-delta`. Callback values, scalar reductions and
tail arithmetic are enclosed. Coordinate rejection avoids integer modulo
bias. Threefry seeds and draw counts are recorded; the confidence theorem
assumes the ideal independent uniform-bit law. A fixed pseudorandom seed is a
replayable realization, not a proof of statistical independence.

Repeated checkpoints use fresh conditional randomness and a total failure
budget, such as `sum delta_k <= delta_run`. Adaptive selection cannot reuse
confidence statements silently. Input/model randomness, MC ensemble error and
algorithmic probability get separate ledger fields. A small failure
probability can still matter in many runs; both delta_run and the consequence
of a false acceptance belong to the decision.

Q0 now exposes a sampled residual **norm at one time**. It does not certify
the integrated residual, trajectory or observable. The global row-sum
envelope is rigorous but cancellation-blind: measured norm upper bounds of
about 103 and 368 on the 125/343-mode inputs miss the residual scale budget.
Its roughly 1% overhead offers no useful accuracy benefit here. Keep it
opt-in; a useful future large-state route needs a tighter analytic envelope
or sketches with controlled action/evaluation error.

## Case choices

The multipliers below are fixed-input baseline medians; they include the
named error path, exclude rendering and generally exclude immutable setup.
Submicrosecond solves have poor ratio resolution; absolute times matter.
Embedded Q1/Q7 diagnostics need their instrumented stage notes. A large
multiplier on a tiny analytic solve does not imply an expensive production
run. See the CSV for inputs, scales, widths, timing method and all repetitions.

| Case class | Current evidence | Default schedule and next escalation |
|---|---|---|
| Elliptic production consumers | Repeated toroidal timings; practical bounds; [scaling ledger](benchmarks.md#production-estimator-scaling) | NN1/NN2/D3 remain above the small-fraction cost target; sheet recovery retains a global serial solve; keep setup and reused-solve denominators distinct |
| Streaming snapshot / comparison | strict analytical/reconstruction checks; snapshot about 111× | E during evolution, final D where implemented; retain S small references; prepare immutable source/cell packets once |
| Damped transport | D residual/inverse check about 44× on 64 cells | E residual/adaptation, final D goal and tangent checks; tighten goal inverse before full-state over-refinement |
| Bounded OU/FP | strict reference check about 370×, submillisecond checker | E flux diagnostic, final S reference packet; general forced-FP D remains pending; reuse a prepared check for immutable data |
| Poisson / manufactured VPFP | tiny interval coefficient/source-balanced oracle checks; coupled about 56× | final D/S; general nonlinear existence and complement are separate missing premises |
| Per-n ideal-MHD response | final frozen-current algebraic defect and directional contraction indicator | E at final state; discretization, current sampling and ideal stability remain unavailable |
| Cold slab | strict Hermite residual about 9× | E embedded RK indicator, final S; larger trajectories can use D action budgets when measured need appears |
| Planar shielding | strict step residual about 102× | E during candidate evolution, checkpoint/final D redesign next; preserve all segment residuals, not just endpoint agreement |
| Debye sphere | strict analytic point evaluation about 371×, about 9 microseconds | prepare/reuse S for small references; D transcendental tables belong upstream when many queries dominate |
| C0 interaction tube | strict Picard continuation about 22× | E structure/step diagnostics; S short encounters; local tube existence margins require online continuation, so deferral is not automatic |
| Q0 / quantum evolution | strict about 26×; D now about 1.4–1.5× | E RK indicator; final/checkpoint D; S audit, A opt-in, P only as a named norm diagnostic |
| Q1 encounter | residual-rate and projection formulas with unenclosed upstream quantities | E with direct estimator timing; D needs enclosed eigen/residual/orthogonality and continuum complement; do not upgrade the formula wrapper alone |
| Q2 packet models | sampled Duhamel about 380×; radius 2.5% of scale | E only at sparse reference checkpoints; estimator exceeds always-on cost target and can be vacuous; goal-specific stability is next |
| Q3 ideal fermions | rigorous orbital-tail interval recurrence about 13× | final S for current small canonical sums; D batching/aggregate recurrence budgets only after conditioning audit |
| Q8 NR / relativistic Fermi EOS | D quadrature/root enclosures about 120× / 338×; widths 0.94% / 0.62% of scale | reuse prepared EOS packets; D at queried thermodynamic states, S regression audit; quadrature/node and observable budgets govern reuse |
| Q4 jellium kick | strict finite-state band about 5.6× | E Krylov/RK residuals, final S for current small cells; matrix-free D requires its own action proof, not dense Q0 conversion |
| Q5 HF/Temple | parent residual plus candidate-gap input about 1.25× | E by default; D requires a verified excited-state/gap premise; interval arithmetic around candidate gap is insufficient |
| Q7 screening | enclosed tails plus quadrature indicators; instrumented total about 2.6× with estimator work larger than candidate | checkpoint/cache E quadrature differences and separate D/S tails; full response remains candidate until integration/input remainders close |
| C3/C4 PME | RMS force estimate about 2.1× / 3.6× for one force action | cache parameter estimate outside hot steps; E only, since positional law and omitted mesh/self terms do not establish an all-configuration bound |
| NN1 axisymmetric P2/P3 | practical fixed-trace RT bound; P3 index 1.179/1.076 on 8/16 rings | local mapped source oscillation; quadrature/arithmetic and physical trace error remain unenclosed |
| H1a/H1b | batch interval summarization negligible beside sampling | statistical intervals with autocorrelation/mixing checks; algorithmic P needs an explicit finite-sample law; increase independent batches/realizations when width misses budget |

## Implemented Q0 arithmetic and approximate mode

The finite Hermitian Hamiltonian is immutable binary64 input. FortNum stores
outward row magnitude sums once and performs point matrix actions with
`gamma_(2n)*sum |H_ij|*max |v_j|`, interval-input widths and absolute underflow
allowances. KIN6D composes real/imaginary actions and hbar division through
the retained symbolic RHS/scaling definitions. Hermite coefficients,
Bernstein residual reduction, norm, radius accumulation and endpoints remain
outward. No dense interval multiplication is required in D; S retains it.
The [upstream argument](/tmp/fortnum-error-strategy/docs/design/error-control-primitives.md)
is retained on local FortNum branch `error-strategy` at `a431cf2`.

`certify_quantum_evolution` can audit a previously solved trace once. It still
integrates every cell's residual; reporting at selected frames never permits
skipping the intervening residual. Partial-cell queries now use the
appropriate prefix plus fraction of the current-cell residual allowance.

A evaluates the dense action at point coefficient endpoints without outward
action rounding. A separate generated allowance covers at most `8n+4`
rounded operations per path, input coefficient widths, hbar scaling and
underflow. The sum of coefficient-action allowances bounds the residual
polynomial's omitted contribution. Each prefix must pass
`evaluation_allowance <= epsilon_b*eta`; otherwise the error query reports
`certification-failed`. Thus `e_true <= eta + allowance <= (1+epsilon_b)*eta`
has analytical support, in addition to strict audits. A's certification flag
stays false, and `(1+epsilon_b)*eta` must fit the observable budget.

The independent binary128 diagonal-unitary reference reports worst
`e_true/eta = 0.220`, worst strict-radius/eta minus one `2.06e-11`, and D
midpoint efficiency index at most 5.34 over three step sizes. Across all
seven sampled fractions the index reaches 49.3 near a right endpoint. This
misses the approximately tenfold efficiency target while still fitting the
1% radius budget; stricter arithmetic alone cannot remove this residual
triangle-bound pessimism. [Raw reference samples](../benchmarks/error-control/reference.csv)
retain all times and modes. These are
scoped validation numbers, not a universal tenfold trajectory-efficiency
claim. A costs about 1.5× solve-only, slightly more than D, so D wins here.
Independent complex two-level, tangent and adjoint tests also exercise D.
Generic sampler coverage uses 2000 streams on each of constant,
heterogeneous and sparse vectors, and Q0 checks 200 streams against an
independent residual. Coverage experiments supplement the tail theorem.

## Elliptic discretization error: equilibrated residuals

The production-consumer scaling ledger is in
[benchmarks](benchmarks.md#production-estimator-scaling). It records repeated
1/2/4/8/16-thread measurements on toroidal NN1/NN2, sheet Ampere, D3 and
coupled/refitted `gs_free` paths. These timings do not upgrade their practical
claims. In particular, setup amortization matters for the sheet solves.

For every elliptic or coercive spatial model, the standard numerical-error
metric is the patchwise equilibrated-residual (Prager--Synge) energy bound.
FortFEM owns one implementation per operator family; KIN6D wraps it with the
physical coefficient, source and guaranteed constants
([GS estimator](grad-shafranov.md#fixed-source-error-estimator)). Local
patch problems are O(N) and parallel; any H(div)/H(curl) reconstruction gives
a valid bound, so local solve errors only loosen it, and the measured
divergence defect absorbs algebraic and load-quadrature inconsistency. The
sampled paths retain approximate status because their non-polynomial
quadrature is unchecked. For affine P1/P2 constant-profile GS, a separate
strict audit (FortFEM `verified_bernstein_majorant`) turns the sampled RTp
flux into an exact H(div) Bernstein field, bounds the 1/R weights by chords
and integrates the resulting polynomials exactly in outward interval
arithmetic. Its certified radius is reported in place of the sampled one; it
is not a check of the sampled quadrature. Other source/geometry families
remain approximate.

Observable bands use `Q_h + r(z_h) +- M_psi*M_z`, with both majorants
equilibrated and one adjoint solve per goal. Nonlinear and free-boundary
problems add a stability constant. Affine fixed-Dirichlet GS with a single
absolute-flux cubic/affine profile now has an outward-rounded global
monotonicity lower bound, with no eigenvalue solve. Its sampled residual still
keeps the complete nonlinear result practical. Free-boundary GS retains the
power-iteration estimate: its flux-level point evaluations are unbounded in
the current H1 energy space, so a discrete eigenvalue enclosure alone cannot
validate the continuous stability premise or nonlinear goal bands.

Production fixed-boundary GS now selects RT equilibration on affine and mapped
P1--P3 meshes, with explicit legacy fallback. Absolute-flux profile laws have a
strong-monotonicity test from a bounded one-sided source slope; unsupported
normalization or noncoercivity retains frozen-residual status. Fixed-source
linear goals use an independently solved adjoint and the product majorant.
Coupled/free-boundary GS uses RT recovery with its existing measured BEM
interface defect and practical stability estimate. The `gs_free` application
writes this practical estimate, M0/beta with its terms, for its final state.
The interface norm tests only the finite reference-BEM trace space; it does
not control its continuous complement. Results therefore separately report
`majorant_exterior_status=unavailable-bem-remainder`. FortFEM's finite
single-layer residual energy diagnostic measures an enriched flux-space
correction, not an upper bound on that complement. A proved full-trace
Neumann residual bound and a continuous single/double-layer remainder are
still required before combining it with the interior equilibrated bound.
The NN1 axisymmetric restriction uses isoparametric RT recovery. Against the
exact Solov'ev flux on 8/16 rings, its efficiency index is 1.014/1.011 at P2
and 1.179/1.076 at P3; the discrete Dirichlet trace is fixed. The former P3
index 1.83/1.82 came from global-C_F source oscillation on curved cells, not
the RT3 flux defect. FortFEM now uses a local mapped-cell Poincare bound with
Bernstein Jacobian/determinant hulls and the directly integrated reference
source residual, retaining the global fallback. Quadrature and trace controls
are recorded in [NN1](equilibrium-nonnested.md). The `grad_shafranov_demo` constant-q
inverse state is estimated through the same fixed-source path with its executed
zero-pressure FF' table as a profile law (claim `practical estimate`, method
`equilibrated-rt`, quadrature `approximate`). The D3 Riesz solve
(`equilibrium_residual` with `error_estimate`) adds a practical upper bound
`sqrt(dual^2 + M^2)` of the continuous dual norm to its finite-test-space lower
bound: `M` is the equilibrated bound of the weighted vector-Laplace Riesz error
with discrete flux `A grad r_h - g`, and the estimator quadrature is not
enclosed. On a coarse Solov'ev psi-lane state the lower bound of a finer test
space lies below it. The NN2 transport solve
(`nt_solve(..., estimate_error=.true.)`) bounds its N = 0 block with tensor RT
recovery (A = R K_pol, A_min = R_min chi_perp) when K does not depend on phi;
a phi-dependent K with modes > 0 couples the modes and reports
`unavailable-mode-coupling`. The bound holds for any anisotropy, but the
hat-function localisation is not robust in it: on the Sovinec test the index
is 1.05 at chi ratio 1, 3--6 at 1e3, and grows with h at 1e6--1e9 (P2
solution polluted). Opt-in FortFEM `correction_sweeps` adds parallel,
damped divergence-free vertex corrections; it leaves the production default
unchanged. Eight sweeps reduce the P2 Sovinec index at ratio 1e9 from 27.5
to 4.5 at nr = 8 and from 389 to 48 at nr = 32. This is an improvement,
not an anisotropy-robust localisation: its contraction degrades with mesh
refinement and its cost exceeds the solve in this small study. The independent
directional-Gram/Cholesky and oversampled-load checks identify the dominant
P2 pollution as a discretization effect
([NN2 diagnosis](equilibrium-nonnested-nn2.md#anisotropy-estimator-diagnostic)).
Per-mode Ampere solves report a coupled discrete
residual-energy indicator (the algebraic residual energy of the discrete
system, also written by `gs_response_3d` for its sheet solves), the optional
fixed-trace interior reconstruction and the coupled reconstruction `ampere_n_coupled_error`: sheet loads enter as
a jump of R H_phi and the discrete Steklov-Poincare trace as a prescribed
boundary trace, so the bound holds for the problem with the Galerkin exterior
(Steklov term zero by construction). For a ring sheet the Galerkin bound
converges at order 1.6--1.8 towards 2 (P2), with index near 1 for large
perturbations. `gs_response_3d` writes it per unit-sheet solve
(`prefix_sheet_error.dat`). The exact-exterior and derivative remainders, the
quadrature and the sheet-data interpolation remain unavailable. The ideal
response of `mhd_ideal_n` now optionally recomputes the frozen-current
algebraic defect at the actual GMRES final state, including sheet elimination,
and reports a practical directional-contraction estimate in the coupled
discrete Ampere energy norm. The sampled ratio is not a contraction upper
bound. Current-sampling and discretization remainders and a continuous ideal
stability bound remain explicitly unavailable; the frozen-current defect is
not a final continuous error bound. These results do not certify all MHD
operators.

Priorities: (1) the equilibrated bound for each existing operator family (H1
affine, curved and Neumann, tetrahedral, mixed and SIPG; planar H(curl) and
the per-n `curl_n` operator); (2) guaranteed stability lower bounds; (3) goal
bands, later two-sided parallelogram bounds; (4) explicit quadrature
remainders.

| Approach | Use | Limitation |
|---|---|---|
| Equilibrated flux, patchwise | Default energy bound and refinement indicator | Needs a reconstruction per operator family |
| NN2 tensor divergence-free correction, opt-in | Reduces weighted flux defect; practical sampled bound | Uniform tensor robustness and small cost remain unavailable |
| Global functional majorant | Superseded fallback | O(N^1.5) global recovery; efficiency about 2 |
| `\|\|\|z_h\|\|\|*M` goal radius | None for linear goals | First order: 10^2--10^6 too wide |
| Product goal band `M_psi*M_z` | Default goal band | One adjoint and estimator per goal |
| Outward global GS coercivity `1-L*C_F^2` | Affine fixed-Dirichlet, one absolute-flux polynomial | Refuses noncoercivity, unbounded positive slope, curved/normalized/multiple-piece laws; residual quadrature remains separate |
| DWR, enriched adjoint | Indicator; nonlinear primal only | Equals a P(p+1) re-solve for linear problems; misses energy goals; radius does not shrink |
| p-, h-refined re-solve, Richardson | Reference checks | Not a bound; Richardson unreliable on non-nested meshes |
| Outward Bernstein audit of the RTp flux | Supported continuous energy certificate | Affine P1/P2, polynomial source numerator, Dirichlet only |
| Free-boundary FEM-BEM interface | Practical finite trace-space indicator | Continuous exterior remainder `unavailable-bem-remainder`; enriched flux residual norms are not upper bounds |

The fixed-boundary measurements behind the table are on a manufactured GS box at 4
threads and up to 7.6e5 interior DOF. The sampled equilibrated bound has
efficiency 1.044--1.08 (P1) and 1.017--1.04 (P2); the certified audit has
1.045--1.12 (P1) and 1.018--1.08 (P2), against 7 (P1) and 140--7000 (P2) for
the former RT0 audit. On the constant-profile production path the certified
radius exceeds the sampled one by at most 0.05%. Estimator and audit cost
grow with DOF exponents 0.97--1.06. At 7.6e5 DOF they cost 20% and 15% (P1)
and 15% and 4% (P2) of mesh, assembly and the Jacobi-CG solve, whose own
exponent is 1.55--1.64. From 1 to 16 threads at 1.9e5 DOF the estimator
speeds up 7.3--8.1 times and the audit 10 times, on a shared host. Goal bands with
equilibrated majorants are 3--8 times the true error. An enriched DWR
correction reproduced the P(p+1) primal result.

Local formulas of these estimators (integrands, geometry factors, reference
jets) are FortSym-generated like every other project leaf
([tooling](tooling.md#fortsym-is-the-mathematical-source-of-project-owned-leaves)).

## Derivatives and long-time dynamics

An equation-level tangent/adjoint has its own norm, scale, stability,
reconstruction, arithmetic and model ledger. A primal certificate does not
certify a gradient. Q0 D/S tangents keep the forcing and primal-radius term;
A tangent certification is unavailable. The initial-state adjoint enclosure
remains separate from the unavailable Hamiltonian-gradient remainder. No
Tier E/P indicator automatically upgrades these claims. Hold random draws
fixed when comparing smooth candidates; seed, sample count, mesh and model
switches are inactive/discrete controls with explicit derivative semantics.

Chaotic C3/C4/H1 trajectories amplify local error and tangent norms until
trajectory bounds are scientifically useless. A slower, tighter tube cannot
restore a useful long-time observable bound by itself. Short-time deterministic
checks still validate kernels; long-time science needs justified ensemble,
invariant-measure or shadowing statements and mixing/finite-sample hypotheses.
A confidence interval is not a deterministic trajectory certificate, and a
shadowing orbit must name the changed initial condition/model and theorem.

Large runs normally perform E continuously and D/P at final/checkpoint goals.
Escalate when the radius budget fails, epsilon_b support fails, probability
assumptions fail, or an invariant/model/interface defect changes applicability.
Record total solve-plus-error cost, setup/reuse, attributed memory and requested
derivatives at matched model and observable width. Costly audit work is explicit
reference evidence, never the hidden default for a large run.

Chris&AI
