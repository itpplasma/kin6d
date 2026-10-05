# Student laboratory and verification ladder

Status: the first prescribed-field 1X1V steady/free-streaming slice is
implemented. The laboratory is the public teaching face of the same engine that
later supports real-device DK/GK/FK; it is not a disposable toy code.

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
Self-consistent Poisson is the next field milestone.

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
| 10 | OU/Fokker--Planck | exact Gaussian mean/variance/equilibrium |
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

Once self-consistent Vlasov--Poisson exists, the same physical cases become the
first representation comparison:

- Eulerian continuum;
- conservative/structure-preserving semi-Lagrangian continuum;
- deterministic PIF;
- charge-conserving PIC;
- finite blobs/remapped particles;
- full-f and delta-f/control-variate variants;
- deterministic/QMC/MC initialization where meaningful.

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

The purpose is not to prove one representation universally superior. Continuum
remains the low-noise reference where practical; PIC/PIF/MC are genuine
first-class methods.

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

Certification is introduced immediately but the exact level follows the
mathematics.

The existing free-streaming example has a reviewed continuous snapshot
enclosure. Poisson should provide coercive/energy-norm bounds. Nonlinear BGK
should use a validated local root theorem when available.

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
