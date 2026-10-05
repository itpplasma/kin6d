# Stability, response and bifurcation

Status: canonical analysis-consumer specification. No general stability engine is
implemented yet; the current 1X1V derivative machinery is the first prerequisite.

## Principle

Stability is an **analysis type**, not a separate plasma model. Every physical
model that defines a differentiable residual/evolution law should be able, when
mathematically meaningful, to expose steady, linear-response, eigenvalue,
nonmodal and nonlinear-stability questions through the same operator family.

For a state \(U_0\) satisfying

\[
R(U_0;\mu)=0,
\]

let

\[
A = R_U(U_0;\mu).
\]

For an evolutionary model

\[
M(U)\dot U + R(U)=0,
\]

the linear spectrum is a generalized operator problem

\[
A\,\delta U = \lambda M\,\delta U.
\]

KIN6D should consume matrix-free actions \(Av\) and \(A^\dagger v\) generated
from the same model definition. A full assembled Jacobian is not a prerequisite.

## Analysis hierarchy

The intended consumers are:

1. **linear initial-value response**;
2. **linear eigenvalue/spectrum**;
3. **forced response/resolvent**;
4. **transient/nonmodal amplification**;
5. **parameter continuation and marginality**;
6. **bifurcation detection and branch following**;
7. **nonlinear saturation/finite-amplitude thresholds**;
8. **energy/free-energy or Lyapunov stability where available**;
9. **fully nonlinear statistical state**.

These are distinct mathematical claims. A stable eigenvalue spectrum does not
exclude transient amplification in a non-normal operator, and a linearly stable
state need not be nonlinearly stable.

## Matrix-free spectral analysis

The first implementation should use only

\[
v\mapsto Av,\qquad
w\mapsto A^\dagger w,
\]

plus mass/preconditioner actions. Suitable algorithms include Krylov--Arnoldi,
Lanczos for self-adjoint cases, Jacobi--Davidson and shift-invert methods when a
verified/preconditioned inverse is available.

The same reduced-physics preconditioners used for steady Newton--Krylov should
be reusable for spectral work. Shift-invert must not silently promote an
approximate preconditioner into an exact inverse.

For certified eigenvalues, specify spectral isolation, normalization, nullspace
constraints and continuous-discretization error. A finite matrix eigenvalue is
not automatically a continuous plasma stability theorem.

## Non-normal and resolvent analysis

For non-normal operators, support resolvent and optimal-growth questions such as

\[
\|(i\omega M+A)^{-1}\|,
\qquad
\max_{\delta U(0)}
\frac{\|\delta U(t)\|}{\|\delta U(0)\|}.
\]

These consumers require both primal and adjoint actions and therefore fit the
FortSym JVP/VJP architecture directly.

Pseudospectra or transient-growth diagnostics are especially relevant for
resistive, kinetic and open-boundary problems where eigenvalues alone can be
misleading.

## Continuation and bifurcation

Treat a parameterized equilibrium/steady branch as

\[
R(U,\mu)=0.
\]

Use pseudo-arclength or another robust continuation method rather than assuming
\(U(\mu)\) is always single-valued. Detect and characterize:

- folds/turning points;
- pitchfork/symmetry-breaking branches;
- Hopf points;
- mode crossings;
- marginal-stability surfaces;
- disappearance of a constrained inverse/coercivity gap.

Continuation is an equation-level consumer. Its derivatives come from the same
IFT/JVP/VJP infrastructure as optimization. It is also provenance: two
continuation histories may reach two distinct stationary solutions at the same
external parameters. KIN6D keeps those as separate equilibria and analyzes the
stability of each separately. It does not average them or assign posterior
weights; TIAGO may later fit the separate states and form an optional finite
mixture of local Laplace approximations.

## MHD stability ladder

The MHD branch should eventually cover, at appropriate model fidelity:

- ideal internal/external kink;
- interchange and ballooning;
- Alfvén continua and discrete Alfvén eigenmodes;
- TAE/EAE/BAE/RSAE and energetic-particle drive;
- resistive tearing;
- resistive-wall modes;
- vertical/free-boundary modes;
- Hall/XMHD/two-fluid modifications;
- nonlinear island and mode saturation.

Tearing and NTM illustrate the distinction between linear and nonlinear
stability. A useful hierarchy is

\[
\text{linear tearing}
\rightarrow
\text{Rutherford island evolution}
\rightarrow
\text{modified Rutherford/NTM}
\rightarrow
\text{full nonlinear resistive/XMHD/kinetic evolution}.
\]

An NTM model must include its declared bootstrap deficit, seed threshold,
polarization/current-drive terms and equilibrium assumptions. It is not merely
a tearing eigenvalue with a renamed label.

ARES, MEPHIT, MARS/GPEC and later JOREK/NIMROD-type cases are comparison
landmarks only in their overlapping models.

## Microstability ladder

For gyrofluid/gyromoment/DK/GK/FK models, stability consumers should cover the
modes supported by the selected physics, including where appropriate:

- drift waves;
- ITG/TEM/ETG;
- kinetic ballooning;
- microtearing;
- drift-Alfvén;
- zonal/GAM response;
- mirror/firehose/whistler;
- two-stream/Weibel;
- energetic-particle resonances.

Every nonlinear turbulence benchmark begins with a matching linear
dispersion/growth calculation from the same model whenever such a linear limit
exists.

## From stability to transport

Linear stability is also a transport-closure provider. A cheap quasilinear
path may map

\[
U
\rightarrow
\{\gamma_k,\omega_k,\xi_k\}
\rightarrow
\text{saturation model}
\rightarrow
\{\Gamma,Q,\Pi,\ldots\}.
\]

Do not hard-code this as scalar diffusivity. The closure may retain mode
structure, cross-phases, spectral tensors and derivatives with respect to the
macro state.

The same stability result should therefore be consumable by both a researcher
studying modes and a transport solver asking for a closure.

## Nonlinear stability and structure

When a model possesses a usable energy/free-energy principle, support stronger
stability statements than eigenvalues alone.

Examples include:

- ideal-MHD energy principle;
- energy-Casimir methods;
- gyrofluid/GK free-energy bounds;
- strongly convex/coercive dissipative blocks;
- Lyapunov or trapping-region arguments.

A verified lower bound on a stability/coercivity margin is scientifically
valuable. As the bound approaches zero, growing sensitivities and certificate
widths are physically expected near marginality.

## Verification and certification

Each stability capability requires:

- analytic dispersion/eigenvalue cases;
- symmetry/nullspace/gauge handling;
- JVP/VJP duality;
- convergence under spatial/velocity/moment refinement;
- an independent code/literature reference where possible;
- explicit distinction between discrete and continuous spectral claims.

Nonlinear stability/saturation additionally requires a trajectory or
finite-amplitude theorem; a linear eigenvalue interval does not prove it.

Chris&AI
