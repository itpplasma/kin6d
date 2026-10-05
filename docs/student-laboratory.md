# Student laboratory and verification ladder

Status: prescribed-field free streaming/Krook balance, single-mode Poisson,
a source-balanced coupled Poisson--FP family and bounded OU/FP equilibrium are
implemented. Their reviewed continuous bounds cover the declared analytical
cases. [Damped transport](general-transport.md) now also accepts arbitrary
represented sources and provides continuous state/tangent bands. Forced FP
and general nonlinear certification remain planned.
The laboratory uses the same engine that later supports real-device DK/GK/FK.

## Purpose

The first examples should let a student understand the equation, derive the
analytic prediction, run the code, inspect the conserved/dissipated quantities,
and see exactly what is and is not certified.

Two doors are visible from the beginning:

1. **kinetic/steady laboratory** — direct equilibria, phase mixing, fields and
   collisions;
2. **turbulence laboratory** — drift waves, zonal structures, interchange and
   gyrofluid dynamics.

The steady branch determines the research architecture. The turbulence branch
reuses the same operators and provides fast nonlinear structure/performance
tests.

## Shared 1X1V problem

The basic electrostatic kinetic model is

\[
\partial_t f_s
+v\,\partial_x f_s
+\frac{q_s}{m_s}E\,\partial_v f_s
=
C_s[f]+S_s,
\]

with optional Poisson coupling

\[
-\partial_x^2\phi=\rho/\epsilon_0,\qquad
E=-\partial_x\phi.
\]

The stationary residual is primary:

\[
R_f=v\partial_x f+(q/m)E\partial_vf-C[f]-S,\qquad
R_\phi=-\partial_x^2\phi-\rho/\epsilon_0.
\]

The first implemented case keeps \(E=0\) prescribed and uses a manufactured
source-balanced Krook state; see [first-example.md](first-example.md).
The [steady laboratory](steady-laboratory.md) adds a genuinely computed Poisson
field within a specified source family and a bounded FP equilibrium. These
cases retain analytical oracles while the next solvers accept independently
supplied sources.

The next steady exercises are general damped transport with a continuous
residual/inverse bound, then forced FP with flux/reconstruction and stability
bounds. State source regularity, velocity boundaries, compatibility and any
normalization constraints. Their derivatives and certificate inputs accompany
the operators from the first implementation. General nonlinear certificates
follow these stable linear cases.

## Microscopic many-particle gallery

A parallel microscopic gallery makes the continuum reductions tangible:

1. two nonrelativistic repulsive Coulomb particles scatter against a Rutherford
   oracle; animate their real 3D trajectories;
2. a few interacting 3D3V particles, then the prescribed magnetic-field case;
3. many physical particles produce a plasma oscillation without a Vlasov
   equation;
4. a weak test charge produces **emergent Debye shielding**. Plot the noisy
   microscopic density/field, its ensemble/time average and the linear
   Debye--Huckel prediction with finite-box and finite-\(N_D\) uncertainty;
5. scan \(\Lambda=n\lambda_D^3=100,300,1000,3000,10000\), with separate
   domain, estimator, fluctuation and collisionality gates;
6. measure \(g(r)\), \(S(k)\), \(S(k,\omega)\) and simple relaxation rates;
7. add relativistic/Darwin/Maxwell and tiny retarded-field comparisons when
   the example's regime needs them.

These are reference/teaching cases, not a requirement to simulate a whole
tokamak with individual physical particles.
Use an electron OCP with an explicit background for the first collective
microscopic cases, preventing attractive classical point-particle collapse.
The [canonical microscopic plan](microscopic-many-body.md) specifies periodic
forces, responsive-species Debye lengths, finite-k wave and estimator caveats.

## Quantum plasma fundamentals gallery

Keep the quantum branch visually and mathematically connected to the classical
examples:

1. Maxwell--Boltzmann versus Fermi--Dirac occupation and the Fermi sphere;
2. map \((\Theta,\Gamma,x_F)\) and place tokamak, solar-core,
   warm-dense and white-dwarf-like regimes;
3. compare classical thermal pressure with nonrelativistic and relativistic
   degeneracy pressure;
4. continue the microscopic shielding example through
   Debye \(\rightarrow\) finite-\(T\) Fermi/Thomas--Fermi
   \(\rightarrow\) Lindhard response and Friedel oscillations;
5. evolve a simple Wigner problem and demonstrate convergence to Vlasov as
   \(\hbar\to0\);
6. show Pauli blocking/quantum diffraction changing a simple collision rate;
7. compare two-electron distinguishable versus antisymmetrized/exchange
   behavior in a tiny reference problem;
8. solve a reduced stellar-structure/EOS example that crosses from thermal to
   degeneracy support without resolving microscopic quantum wavelengths.

The point is not a separate quantum course hidden in KIN6D. Each example
explains exactly where a classical plasma assumption succeeds or fails and
retains the same provenance/error/verification conventions as the classical
gallery.

## Kinetic teaching sequence

Retain the following named examples as the permanent bottom of the verification
ladder.

| # | Example | Main analytic/structural oracle |
|---|---|---|
| 01 | Free streaming | \(f(x,v,t)=f_0(x-vt,v)\) |
| 02 | Constant acceleration | exact characteristic map |
| 03 | Harmonic phase-space rotation | exact symplectic rotation |
| 04 | Poisson Fourier mode | \(\phi_k=\rho_k/(\epsilon_0k^2)\) |
| 05 | Plasma oscillation | linear dispersion/frequency |
| 06 | Landau damping | plasma-dispersion root plus finite-time caveats |
| 07 | Two-stream | linear growth before nonlinear trapping |
| 08 | BGK equilibrium | direct nonlinear steady residual |
| 09 | Krook balance/relaxation | exact steady and exponential limits |
| 10 | OU/Fokker--Planck | bounded zero-flux equilibrium; Gaussian formulas in the matching unbounded model |
| 11 | Lorentz scattering | spherical-harmonic decay rates |
| 12 | Landau collision box | invariants, entropy and independent rates |

Every example states domain, normalization, boundary conditions, exact scope of
its analytical formula, reconstruction error and certificate status. Do not
periodize a velocity or harmonic-oscillator domain merely to reuse a periodic
implementation.

Landau damping is an initial-value problem with phase mixing; a fitted complex
root is not the complete nonlinear finite-time solution. BGK examples must be
self-consistent solutions, not the statement that any arbitrary \(F(H)\) is a
plasma equilibrium.

## Representation laboratory

Introduce the first comparisons before nonlinear BGK, Landau damping or the
full particle infrastructure:

| Case | Carriers/coupling | Independent oracle |
|---|---|---|
| Exact free streaming | continuum and deterministic characteristic particles | translated distribution, density mode and parameter derivatives |
| Neutral single-mode Poisson | continuum projection, direct PIF and simple local-basis PIC | field amplitude and electric energy |
| Short self-consistent Vlasov--Poisson | continuum, PIF and PIC with common initial data | earlier limiting cases, independent refinement and balance tests |

Free-streaming particles alone do not exercise PIC/PIF field coupling. The
first rung now has an [animated comparison](streaming-comparison.md) with
deterministic weighted particles and a Fourier/P1 Eulerian solution. The
Poisson comparison exercises deposit/solve/gather without claiming a
self-consistent trajectory. The third case supplies field feedback and reports
the trajectory's actual bound status. Share the physical laws, parameters and
boundaries; expose each representation's quadrature, shape/basis and time
errors separately.

Use deterministic sampling first. Conservative semi-Lagrangian methods,
geometric/current-conserving PIC, blobs/remapping, delta-f/control variates and
randomized variants enter when a concrete case needs them. Adaptive particles,
NUFFT, stochastic collisions and GPU/MPI are later steps. Early electrostatic
PIC establishes only the conservation and derivative properties actually tested.

Students should be able to ask: "same Vlasov equation, why do these numerical
representations behave differently?"

Compare:

- field history and growth/damping;
- phase-space reconstruction;
- energy/momentum/charge;
- noise;
- memory;
- runtime;
- derivative behavior.

Optimize measured total time at matched model, observable accuracy and bound
status, including certification and requested derivatives. Continuum remains
the low-noise reference where practical; PIC/PIF/MC are first-class methods.

## Turbulence sequence

The first turbulence primitive is the two-dimensional Poisson bracket

\[
[a,b]=\partial_xa\,\partial_yb-\partial_ya\,\partial_xb,
\]

defined once in FortSym with generated primal/JVP/VJP and an explicitly declared
discrete pairing/projection.

### Hasegawa--Mima

Start from a periodic normalized HM model such as

\[
q=\nabla_\perp^2\phi-\phi,\qquad
\partial_tq+[\phi,q]-\kappa\partial_y\phi=0.
\]

Verify first the linear drift-wave frequency for one Fourier mode, then the
nonlinear bracket. In the undriven ideal case monitor the correct energy and
generalized enstrophy. Antisymmetry alone does not prove Jacobi or a full-time
conservation theorem.

### Hasegawa--Wakatani

Add resistive parallel electron response. State the zonal convention and the
precise large-response limit rather than assuming every HW normalization tends
to the same HM equation.

### Interchange/blob

Add curvature/interchange drive with an analytic linear-growth limit and
explicit particle/energy source and sink ledger. This is the first visually
direct edge-transport example.

### Minimal ITG gyrofluid

Add temperature/FLR physics only after the chosen low-moment model and its
free-energy structure are written down. Do not call an unrelated isothermal
gyrofluid code an ITG oracle.

### Electromagnetic drift-Alfvén

Introduce \(A_\parallel\), parallel response and an Alfvénic linear benchmark
before nonlinear turbulence.

The turbulence gallery is intentionally small and physically transparent. Its
operators later become pieces of the gyrofluid/gyromoment and reduced-MHD
branches.

## Structure from the first lecture

The laboratory should expose the common project structure rather than hide it.

For reversible transport, show characteristics and the relevant Hamiltonian or
Poisson invariant. For dissipative terms, show the appropriate entropy/free
energy or open-reservoir exchange. For driven/open systems, show

\[
\frac{dE}{dt}=P_{\rm source}-P_{\rm diss}-P_{\rm boundary}
\]

with all internal exchanges cancelling where the model says they should.

This gives a direct conceptual progression:

\[
{\rm characteristics}
\rightarrow
{\rm Vlasov\ self-consistency}
\rightarrow
{\rm collisions}
\rightarrow
{\rm nonlinear\ turbulence}
\rightarrow
{\rm reduced/kinetic\ hierarchies}.
\]

## Differentiability

Simple cases carry exact parameter derivatives from the start. Compare
FortSym-generated JVP/VJP with analytical derivatives and independent
directional differences.

Particle derivatives include weights, positions, deposit/gather and field
response. Declare shape regularity and cell/event surfaces; test the stated
smooth domain. Randomized derivatives state whether random variables are held
fixed, and retain statistical uncertainty separately from deterministic error.

Self-consistent Poisson/BGK introduces the first implicit-state relation

\[
R_UU_\theta=-R_\theta
\]

and adjoint

\[
R_U^\dagger\lambda=Q_U^\dagger.
\]

This is the same formalism later used for equilibrium and stellarator
optimization.

## Certification

Each reported state, field, observable and sensitivity carries its applicable
bound/status, or an explicit statement that a bound is unavailable. The current
analytical laboratories have reviewed continuous enclosures in their declared
domains. General transport and forced-FP certificates must bound the candidate
without substituting a known exact solution. Nonlinear BGK needs a validated
local root/stability theorem and unresolved-complement control when available.

The practical default is approximately 1% accuracy for a named observable with
a declared scale; use an absolute scale when its value can vanish. A validated
coarse rounding budget is sufficient if the full bound meets that target.
Deterministic containment is required for a certified result, while tighter or
higher-precision audits are optional.

Plot the full applicable intervals. Distribution/field bands require pointwise
or whole-cell bounds, including the plotted reconstruction; an L2 radius is
not a pointwise radius. Atomic particle measures use weak-observable bounds or
bounds for an explicitly defined smooth reconstruction. Histories identify
whether bounds cover observation times or also the interpolation between them.

For turbulence, invariant ledgers and convergence are not by themselves
continuous trajectory certificates. A case remains a candidate until its
actual temporal/spatial theorem is established.

Monte Carlo examples report statistical confidence separately from
deterministic discretization bounds.

Lean is not required for laboratory promotion. An optional formal checker may
later validate generic certificate theorems.

## Performance

Every example has a tiny laptop case and, once relevant, a realistic
performance shape. The laboratory is also the first benchmark arena for:

- allocation-free generated leaves;
- vectorization;
- matrix-free/tensor actions;
- particle tiling;
- PIC/PIF crossover;
- direct-vs-FFT reductions;
- forward/JVP/VJP overhead.

A colorful movie is not an excuse to bypass accuracy or performance evidence.

Reusable upstream benchmark support comes before new hot-path implementation.
Retain fixed inputs and repeated measurements with variability, memory and
separate setup, solve, certificate and derivative costs. Compare total time at
the same requested accuracy and claim strength; the fastest result is a
measured goal, not a claim made from one kernel timing.

## Progression to research cases

The teaching ladder feeds directly into:

\[
0X3V\ {\rm collisions}
\rightarrow
1X3V\ {\rm steady\ transport}
\rightarrow
{\rm NEO\!-\!2/SFINCS\ level\ DK}
\rightarrow
{\rm gyrofluid/gyromoments/GK}
\rightarrow
{\rm global\ FK}.
\]

Earlier examples remain permanent regressions and limiting-case plots.

Chris&AI
