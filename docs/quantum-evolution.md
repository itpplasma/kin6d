# General finite-basis quantum evolution

Status: implemented Fortran candidate and rigorous continuous-time residual
certificate for a supplied constant Hermitian matrix. This certificate names
the finite Hilbert space and exact represented input parameters. It does not
certify the continuum Coulomb Hamiltonian, a basis truncation, or a many-body
model reduction. See [the quantum scope](quantum-plasma.md) and
[the certificate policy](certification.md).

The admitted model is

\[
i\hbar\dot\psi=H\psi,\qquad \psi(0)=\psi_0\in\mathbb C^D,
\quad H=H^*,\quad \hbar>0.
\]

`Hre`, `Him`, and the initial real/imaginary amplitudes are independently
supplied finite binary64 arrays, interpreted as exact real parameters.
Hermiticity is checked exactly, including a zero imaginary diagonal. Negative
energies, signed or complex amplitudes, a zero state, and arbitrary dimension
are admitted; no Gaussian or physical-density positivity assumption enters.
The assembled Hamiltonian action is a general dense reference, independent of
the Fock/CAR matrix-free parent action. Constant complex Hermitian input is
supported; time-dependent or non-Hermitian generators are rejected by scope.

## Candidate and continuous certificate

`initialize_quantum_evolution` allocates the trajectory and reusable stage
buffers once. The untrusted candidate uses FortNum's existing six-stage
Cash-Karp RK5(4) method. KIN6D supplies generated local `-i H psi / hbar`
contractions; it does not implement another Runge-Kutta method. The embedded
error estimate is not the certificate. All project-owned scalar action,
sampling, enclosure, and derivative leaves come from
`gen/gen_quantum_evolution.f90` and its shared FortSym definitions. The
generic Hermite enclosure, Bernstein residual bound and radius accumulation
are FortNum's `fortnum_ode_residual_certificate` (pinned `a431cf2`); KIN6D
supplies the action `-i H x/hbar` (and `-i dH x/hbar` as tangent forcing).
The default `QE_WORKHORSE` uses point matrix actions with proved aggregate
rounding, coefficient-width and underflow allowances from FortNum.
`QE_STRICT` retains the dense outward action for audits. Both enclose the
same mathematical finite model. The historical strict cost of about 15 ms
per 125-mode cell motivated this change; the measured default now takes
about 1.4–1.5 times solve-only total, versus 25–27 times for strict audit.
See [the error-control strategy](error-control-strategy.md) for assumptions,
matched fixtures and efficiency indices.

Each time cell stores its binary64 candidate endpoint amplitudes and RHS
slopes. These bits define an exact cubic Hermite reconstruction `p(s)`,
`0 <= s <= 1`; the exact reconstruction is continuous at every node. Strict
interval coefficients enclose its exact polynomial coefficients. The checker
forms, independently of any exact solution,

\[
r(s)=p'(s)/h+iHp(s)/\hbar.
\]

This vector polynomial has degree at most three. Its power coefficients are
converted with exact rational weights into four outward interval Bernstein
control vectors. For a control box, the squared absolute component maxima
bound its Hilbert norm; the maximum of the four norm upper endpoints bounds
`sup_s ||r(s)||_2`. A Bernstein polynomial is a convex combination of its
controls, so this covers the complete cell, including its interior.

The exact generator `A=-iH/hbar` is skew-adjoint. Its propagator is unitary,
and Duhamel's formula gives

\[
\|\psi(t)-p(t)\|_2\le
\|\psi(0)-p(0)\|_2+\int_0^t\|r(\tau)\|_2\,d\tau.
\]

The initial defect is zero. Accumulating `h * residual_bound` with outward
arithmetic gives endpoint bounds; a cell's endpoint radius conservatively
covers every time inside that cell. The nonnegative domain of an analytically
known sum of squares is intersected before strict square root, after checking
finite endpoints. No overflow is discarded by this tightening.

Sampling returns a floating candidate and a Hilbert radius that additionally
covers strict reconstruction versus floating evaluation. Component boxes may
be widened by the same radius, because a coordinate projection has norm one.
This is a finite-vector statement; it is not a pointwise continuum
wavefunction, spatial density, or pair-density band.

`solve_quantum_evolution(...,certify=.false.)` followed by
`certify_quantum_evolution(...,error_mode=QE_WORKHORSE)` supports final or
checkpoint checks of a frozen trace. Every intervening cell is audited.
Queries at fraction s use the prefix radius plus s times the current-cell
residual budget, and enclose sampling error separately.

`QE_APPROXIMATE` is opt-in. It reports `approximate-bound` with epsilon_b
(default 5%) through `quantum_evolution_error`; its `bounded` flag is false.
A separate generated allowance must prove omitted evaluation error at most
epsilon_b times the displayed radius at each prefix, or the mode fails.
The inflated radius must fit the observable budget. It currently costs
slightly more than the workhorse, so offers no speed advantage on these cases.

`quantum_evolution_indicator` returns the already computed maximum embedded
RK norm with status `indicator`. `quantum_evolution_sampled_residual` returns
a law-declared probabilistic upper bound for the residual norm at one chosen
time. It does not bound residual integration or the state trajectory; its
current global envelope is too loose on the benchmark inputs. Neither query
upgrades the primal, continuum or derivative status.

## Derivatives and observables

For a Hermitian direction `dH` and arbitrary complex `dpsi0`, the equation-level
tangent satisfies `z'=Az+dA*psi`, with fixed `hbar` and physical time. Its
candidate uses the retained primal polynomial as forcing. The checker bounds
the full-cell tangent polynomial residual `rt=q'-Aq-dA*p`. With primal radius
`rho`, unitarity gives the incremental tangent radius

\[
h\bigl(\|r_t\|_{\infty,2}+\|dA\|_F\rho\bigr).
\]

The Frobenius norm is a valid upper bound for the induced norm; the cell's
endpoint primal radius covers the forcing discrepancy for its whole cell.
Both initial data are exact, so neither initial defect is omitted.

`quantum_evolution_vjp` uses the real complex Hilbert pairing
`Re(lambda* psi)`. For a terminal covector it solves the backward adjoint
equation and returns an initial-state VJP with its own certified Hilbert
radius. The full real/imaginary matrix VJP is an equation-level adjoint
quadrature of `Re(lambda* (-i*dH*psi/hbar))`, restricted to Hermitian directions.
It is a candidate: `h_gradient_bounded=false` explicitly reports the missing
continuous adjoint-quadrature goal certificate. Neither adaptive decisions nor
Cash-Karp iteration steps are differentiated. Evolution derivatives with
respect to `hbar` or observation time are not exposed in this increment.
The generated scalar square-root norm derivative has domain `value > 0`;
its primal and the state evolution/tangent certificate also admit a zero state.
Discrete dimensions, step counts and cell selection are scheduling choices.

For a bounded operator `B`, an unnormalized expectation discrepancy is bounded
by `||B|| * rho * (2*||p|| + rho)`, with strict contraction arithmetic added
separately. Generated `qe_goal_error_interval` emits this perturbation factor.
Normalized observables need a positive denominator enclosure; normalization
and physical RDM representability are not inferred from arbitrary amplitudes.

## API and result status

The module is `kin6d_quantum_evolution`, with private `quantum_evolution_t`.

- Initialize with `(state,Hre,Him,psi_re,psi_im,steps,h,hbar,ok[,budget])`.
- Solve with `(state,ok[,certify,error_mode])`; candidate-only calls clear bounded status.
- Sample with `(state,cell,fraction,re,im,radius,ok,bounded)`; the interval-box
  variant returns already widened true component boxes.
- The exact model time is `(cell-1+fraction)*h`, with represented fraction in
  `[0,1]`; generated `qe_time_interval` encloses that time for reporting.
- Tangents take `(state,dHre,dHim,dpsi_re,dpsi_im,ok[,certify])` and have a
  separate sampler/radius/status. Certifying a tangent also prepares any
  missing primal certificate once.

Nonfinite data, wrong shapes, non-Hermitian input, invalid time/fraction,
allocation failure, and failed interval arithmetic reject their corresponding
operation. A certificate is reusable only for the private immutable inputs and
retained trajectory; rebuilding the primal clears derivative results.
`bounded` means a valid continuous-time finite-basis enclosure. The separate
`quantum_evolution_accepted` flag checks the reconstruction radius against the
declared accuracy budget. Its default is `1e-12 + .01 * initial_norm_upper`;
an unmet target never becomes an accepted result. Consumers also check the
actual sampled/observable radius against each output's declared target.

Hamiltonian assembly rounding, continuum-basis complement, physical data and
model errors are separate ledgers and are unavailable through this API.
Passing the represented-matrix certificate does not silently certify those
errors. A supplied operator-norm matrix defect could support a separate
unitary perturbation bound, but that input ledger is not implemented here.

## Verification and costs

The independent oracle uses a binary128 two-level matrix exponential only in
the test, plus a signed five-dimensional diagonal unitary case and the zero
generator. It checks interior Hilbert/component containment, candidate order,
tangent finite differences against the independent exponential, a certified
initial adjoint, terminal matrix/initial pairing, and invalid/failure states.
No analytical trajectory enters the production checker.

Work is `O(steps*D^2)` for the dense action and strict residual audit. Retaining
nodes, slopes and cubic coefficients costs `O(steps*D)`; tangent storage is
allocated once on first request. The VJP additionally retains an adjoint
trajectory and `O(D^2)` quadrature scratch. There is no global Jacobian or
iteration tape and no per-stage allocation. The aggregate-rounding evaluator
is implemented; matrix-free large-basis evolution and checkpoint recomputation
remain later measured optimizations, not current support claims.

`test_quantum_evolution --benchmark` retains 20 raw warm CPU trials. Each trial
averages 100 complete pipelines for the same two-level inputs and four time
cells, with a 1% state norm-scale target and 1% tangent scale1 target. It records
candidate-only, default complete, candidate-plus-tangent and default combined
costs, achieved radii, absolute extra times and multipliers. Validation and
sampling occur outside timers; primal certification is not counted twice in
the combined path. This is warm compute evidence, excluding initialization,
sampling, rendering and external ranking. Cold setup/first-result and the
larger gallery workload require their own receipts.
