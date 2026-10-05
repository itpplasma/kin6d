# Transport and multiscale coupling

Status: canonical transport/multiscale specification. KAMEL/QL-Balance and
GORILLA provide important internal heritage/oracles, but no generic KIN6D
transport engine is implemented yet.

## Principle

Transport is a **macro evolution problem** whose closure can range from an
empirical coefficient to the direct flux of a global nonlinear kinetic solve.
KIN6D should support this entire ladder without forcing all microphysics into a
scalar diffusion coefficient.

For transported macro variables \(u_a\),

\[
\partial_t u_a + \nabla\cdot\Gamma_a = S_a.
\]

The transport state may contain, as demanded by the chosen model,

\[
n_s,\ T_s,\ \text{rotation/momentum},\ j_\parallel,\ \psi,\
n_Z,\ \text{fast-particle moments},\ldots
\]

and their 1D, 2D or 3D spatial dependence.

## Transport dimensionality

Transport dimensionality is independent of the microphysics model:

| Level | Meaning |
|---|---|
| 0D | global inventory/burn/circuit-like balance |
| 1D | flux-surface or other single-coordinate transport |
| 1.5D | 1D transport with repeatedly solved 2D equilibrium |
| 2D | transport retaining a second physical/angular coordinate |
| 3D | full configuration-space macro transport |
| 5D/6D | direct GK/FK evolution with no separate macro closure |

The 2D coordinates are not hard-coded. They may be \((R,Z)\),
\((\rho,\theta)\), or another mapped pair. Three-dimensional transport must not
assume flux surfaces and is important for stellarators, RMPs, islands,
stochastic fields and edge/divertor geometry.

## Closure-provider contract

The core abstraction is not a scalar diffusion-coefficient getter. A
microphysics provider maps the macro state and geometry to whatever closure
data it can justify:

\[
\mathcal M[U]
\mapsto
\{
\Gamma,\ Q,\ \Pi,\ J,\
L_{ij},\
K(x,x'),\
C_{\alpha\beta}(x,k,\omega),\
D_U\Gamma,\
\epsilon
\}.
\]

A provider need only supply the supported subset. Examples:

- prescribed \(D,\chi,V\);
- neoclassical Onsager matrix;
- quasilinear mode spectrum and fluxes;
- nonlinear local flux-tube fluxes/Jacobians;
- nonlocal orbit/global transport kernel;
- direct global gyrofluid/GK/FK flux fields;
- electromagnetic cross-spectral statistics.

The transport solver consumes a mathematically declared closure, not a
code-specific file format.

## Closure-fidelity ladder

The intended progression is

\[
\text{empirical/prescribed}
\rightarrow
\text{neoclassical/local}
\rightarrow
\text{linear/QL}
\rightarrow
\text{nonlinear local}
\rightarrow
\text{nonlocal orbit/global kernel}
\rightarrow
\text{global nonlinear microphysics}
\rightarrow
\text{fully coupled evolution}.
\]

A local Fick-like law

\[
\Gamma=-D\nabla n + Vn
\]

is the cheapest member, not the universal interface.

A nonlocal model may instead use

\[
\Gamma(x)=\int K(x,x';U)X(x')\,dx',
\]

where \(X\) denotes thermodynamic forces or profile perturbations.

## Traditional integrated modelling as an early mode

KIN6D should support the historically successful equilibrium/transport
separation directly:

\[
\text{equilibrium}
\rightarrow
\text{transport step}
\rightarrow
\text{updated profiles/current}
\rightarrow
\text{new equilibrium}.
\]

Fidelity can then increase continuously:

\[
\text{frozen equilibrium}
\rightarrow
\text{periodically re-solved equilibrium}
\rightarrow
\text{implicitly coupled quasi-static equilibrium+transport}
\rightarrow
\text{dynamic MHD/kinetic evolution}.
\]

This cheap mode is important for scenarios, scans, control and teaching and
should not be delayed until global GK/FK exists.

ASTRA-like 1D/1.5D transport is therefore a valid KIN6D configuration, not an
external workflow concept.

## KAMEL/QL-Balance heritage

The current KAMEL repository already provides a concrete internal example:

\[
\text{profiles}
\rightarrow
\text{KiLCA/KIM linear kinetic response}
\rightarrow
\text{quasilinear transport tensor}
\rightarrow
\text{QL-Balance profile evolution}.
\]

KIN6D should first reproduce this pattern in memory through common state and
closure interfaces, preserving QL-Balance/KAMEL as independent regression
oracles.

Then generalize beyond one-dimensional profiles and beyond one response model.

## Stability-to-transport coupling

A linear microstability solve can provide

\[
\{\gamma_k,\omega_k,\xi_k\}
\]

plus a saturation rule and cross-phases to produce transport.

The closure should retain richer information than \(D,\chi\) whenever useful.
At minimum, support complex mode structures and cross-spectral quantities
needed for particle, heat, momentum and electromagnetic transport.

This connects [stability](stability.md) directly to transport and lets the same
linearized operator serve both science and integrated modelling.

## Electromagnetic spectrum/cross-spectrum closure

A turbulence provider may return spectral statistics such as

\[
C_{\alpha\beta}(x,k,\omega)
=
\left\langle
\delta A_\alpha(x,k,\omega)
\delta A_\beta^*(x,k,\omega)
\right\rangle,
\]

with \(A_\alpha\) drawn from the selected fluctuating field/distribution
variables.

Cross-spectra and phases matter. A scalar power spectrum is insufficient for
fluxes such as \(\langle\tilde n\tilde v_E\rangle\).

Possible providers include:

- linear modes + saturation model;
- nonlinear local flux tube;
- sparse samples from a global simulation;
- a global nonlinear gyrofluid/GK/FK state.

The consumer sees one closure contract with declared statistical/temporal
averaging and uncertainty.

## Global orbit transport

Guiding-center/full-orbit transport provides an important intermediate regime
between local 1D coefficients and global 5D/6D turbulence.

A GORILLA-like orbit engine can accumulate transition probabilities, response
kernels or direct particle/energy/momentum transfer between macro cells:

\[
K(x,x';E,\mu,\ldots)
\]

or an equivalent conservative transport operator.

This is especially important for 2D/3D stellarator, RMP, island, stochastic and
finite-orbit-width transport where a radial local coefficient is inadequate.

The public GORILLA code is currently an orbit/characteristic engine; KIN6D
should treat richer global-transport/spectrum coupling as a planned extension,
not claim it as existing public GORILLA functionality.

## Turbulence spectrum along global orbits

A particularly valuable multiscale mode is

\[
\text{local/global turbulence statistics}
\rightarrow
\text{orbit-resolved resonant interaction}
\rightarrow
\text{2D/3D macro transport}.
\]

Schematically, for action variables \(J\),

\[
D^{ij}(J)
=
\sum_k\int d\omega\,
\mathcal R^{ij}(J;k,\omega)\,
C(k,\omega),
\]

where \(\mathcal R\) is computed using the actual global orbit/geometry.

This preserves geometry and resonance information that is lost when
microturbulence is collapsed immediately into a radial scalar diffusivity.

The first implementation should start with a simple prescribed spectrum and an
analytic orbit/resonance limit before consuming turbulence generated elsewhere.

## Coupling-strength ladder

KIN6D should support progressively tighter coupling:

### C0 — offline/table closure
Precomputed coefficients/tables feed transport.

### C1 — in-memory loose coupling
At macro step \(n\),

\[
U_n\rightarrow \mathcal M[U_n]\rightarrow U_{n+1}.
\]

### C2 — iterated self-consistent closure
Macro state and micro closure are iterated to convergence at each macro step or
steady transport solve.

### C3 — separated-timescale multiscale
Fast local/global microproblems are embedded in slow macro transport, using
spatial/time scale separation.

### C4 — sparse global-micro sampling
Global nonlinear microphysics is sampled only where/when the macro solver or
adjoint error indicator requires it; interpolation/surrogates carry explicit
validity/error status.

### C5 — fully coupled evolution
When scale separation/locality fails, evolve the global nonlinear model
directly.

TRINITY- and Tango-like strategies are conceptual precedents for C3/C4.
GENE-X-like full-f global evolution is a reference for the C5 end of the
spectrum. These are comparison architectures, not code modes.

## Locality and timescale are separate questions

True multiscale modelling must distinguish:

1. whether the fast and slow timescales separate;
2. whether the microphysics is spatially local.

A transport calculation may exploit

\[
\tau_{\rm micro}\ll\tau_{\rm transport}
\]

while still using **global** gyrokinetic/orbit microphysics rather than local
flux tubes.

Thus the hierarchy is two-dimensional:

\[
\text{local}\leftrightarrow\text{nonlocal/global},
\qquad
\text{loose}\leftrightarrow\text{fully coupled}.
\]

Do not hard-wire multiscale to mean local flux tubes.

## Differentiable coupled transport

A closure provider should expose, when available,

\[
D_U\Gamma\,v
\quad\text{and/or}\quad
D_U\Gamma^\dagger w.
\]

Then the macro transport balance can be solved implicitly:

\[
R_{\rm tr}(U,\mathcal M[U])=0
\]

rather than by fragile explicit profile iteration.

For cheap symbolic closures these derivatives are generated directly. For
expensive nonlinear microphysics they may come from tangent/adjoint solves,
validated secants or a surrogate with an explicit error/validity envelope.

Adjoints also determine which closure evaluations matter most to a requested
observable and therefore where expensive microphysics should be recomputed.

## Verification and conservation

Every transport model declares conservation of:

- particles/nuclei;
- charge;
- momentum;
- energy;
- current/flux where relevant.

Internal micro-to-macro transfers cancel in the total ledger. Sources, sinks,
radiation and walls are explicit ports.

Verification starts with exact diffusion/advection/source problems in 1D,
then separable 2D/3D manufactured cases, Onsager/neoclassical limits, and
cross-code integrated cases.

The first 2D/3D transport implementation should use the same FortFEM mapped
spaces as the rest of KIN6D rather than introduce an independent transport mesh
stack.

## Certification

Transport certificates separate:

- macro discretization/time error;
- closure-provider numerical error;
- closure/model-reduction error;
- interpolation/surrogate error;
- statistical error for sampled microphysics;
- equilibrium/geometry uncertainty.

A quasilinear coefficient with a tiny numerical error is not a rigorous
turbulence-model error bound.

Goal-oriented adjoints should eventually drive where to refine the macro grid,
recompute microphysics, increase spectrum resolution or replace a cheap closure
by a richer one.

Chris&AI
