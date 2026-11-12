# Stability, response and bifurcation

Status (2026-10-21): canonical analysis-consumer specification. Model-specific
MHD energy, linear-response and free-boundary Hessian consumers are implemented.
The free-boundary probe now retains each full cosine/sine harmonic pair and
checks its eigen-residual; an analytic cylinder external-kink limit exercises
the same Hessian and vacuum elimination. See
[ideal-free-boundary.md](ideal-free-boundary.md#harmonic-pair-stability-audit).
The first [S0 plasma-volume spectral slice](s0-spectrum.md) adds physical
inertia, sparse real/complex shift targeting and real-shift inertia diagnostics
for caller-supplied displacement spaces.
A positive shifted-iteration result does not establish an inertia count or
continuous stability. A general consumer covering the hierarchy below
remains planned; the plan for global 3D ideal-MHD stability is
[the P1e section](#global-3d-ideal-mhd-stability-p1e); no ITER n=3 continuous stability certificate is available.

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

MILCA should consume matrix-free actions \(Av\) and \(A^\dagger v\) generated
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
external parameters. MILCA keeps those as separate equilibria and analyzes the
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

## Global 3D ideal-MHD stability (P1e)

Status (2026-11-12): stages P1e.1 and P1e.2 implemented, P1e.3 and P1e.5 in
progress (see
[implementation status](#p1e-implementation-status)); later stages planned.
This section is the plan
for roadmap item P1e: a global ideal-MHD stability capability for 3D
equilibria, fixed and free boundary, that covers the published capabilities of
CAS3D, CASTOR3D, TERPSICHORE and DCON and goes beyond them in certified errors,
derivatives and cost at matched accuracy. Every source below is a publication
or an open-source repository; the provenance rule is "Source provenance
(strict)" in [AGENTS.md](../AGENTS.md).

### Reference codes as published

Each entry gives its evidence: **P** stated in the title or abstract of the
cited paper, **S** taken from a secondary description (another paper or a
conference abstract) because the full text has not been read yet, **?** not
established. Stage P1e.0 below replaces every S and ? by a full-text reading
or removes the entry. A closed code enters only through these publications and
their published inputs and results.

| Capability | CAS3D | CASTOR3D | TERPSICHORE | DCON / DCON3D |
|---|---|---|---|---|
| Equilibria | 3D nested (P) [1,2] | 2D and 3D tokamak (P) [5]; stellarator applications (S) | 3D nested (P) [7] | axisymmetric (P) [9]; 3D nested (P) [10] |
| Model | compressible ideal, kinetic energy with density (P) [2] | resistive MHD with viscosity and toroidal flow (S) [5]; two-fluid terms added later (S) | ideal (P) [7,19] | ideal, no kinetic energy (P) [9] |
| Output | unstable global modes, stable spectra, Alfven eigenmodes (P) [2] | growth rates over five orders, from no-wall to resistive-wall time scales (S) [5] | low-n global modes (P) [7] | stable/unstable decision without eigenvalues; plasma response matrix (P) [9,10] |
| Boundary | free boundary (P) [3] | ideal and resistive 3D thin walls, coils (S) [5,6] | ? | free boundary via plasma plus vacuum response matrix (P) [9]; 3D: ? |
| Walls and feedback | conducting structures and feedback design with STARWALL (P) [4] | STARWALL thin-wall model (S) [4,5] | ? | ? |
| Local criteria | ? | ? | Mercier (P) [7] | singular-surface crossing (P) [9] |
| Coordinates | Boozer (S) | ? | ? | straight field line flux coordinates (P) [9,10] |
| Derivatives, error bounds | none found so far | none found so far | none found so far | none found so far |
| Availability | closed | closed | Apache-2.0 [R1] | MIT within GPEC [R2,R3]; DCON3D: ? |

The P1e target is the union of the ideal rows, plus what no publication found
so far reports: an error estimate or bound, or an explicit unavailable status,
on every reported number; derivatives of every observable with respect to
equilibrium and design parameters; and measured cost at a declared accuracy.
The resistive, viscous, flow and two-fluid rows of CASTOR3D are not ideal MHD;
they enter only after the ideal gates, as stage P1e.11.

### Notation

B, J, p and rho are the equilibrium magnetic field, current density
(J = curl B / mu_0), pressure and mass density; mu_0 is the vacuum
permeability and gamma the adiabatic index. xi is the plasma displacement,
nu the outward unit normal of the plasma boundary, and n the toroidal mode
number. A perturbation proportional to exp(-i omega t) has eigenvalue
lambda = omega^2, and lambda < 0 means growth rate sqrt(-lambda). The unit of
frequency is the on-axis Alfven frequency omega_A = v_A0/R_0, with
v_A0 = B_0/sqrt(mu_0 rho_0) from the field and density on the magnetic axis
and R_0 the major radius of the axis. After discretization, A is the
Hermitian stiffness matrix (twice delta W) and M the mass matrix (twice the
kinetic energy), and A x = lambda M x.

### Model hierarchy

1. **M0 fixed boundary.** Static equilibrium without flow on nested 3D flux
   surfaces, potential energy [11,12]

   \[
   \delta W_p=\frac12\int\Big[\frac{|Q|^2}{\mu_0}
   -\xi^*\cdot(J\times Q)
   +\gamma p|\nabla\!\cdot\xi|^2
   +(\xi\cdot\nabla p)\,\nabla\!\cdot\xi^*\Big]dV,\qquad
   Q=\nabla\times(\xi\times B),
   \]

   which is real on an equilibrium because the force operator is
   self-adjoint, and kinetic energy \(\tfrac12\int\rho|\xi|^2dV\), with the
   boundary condition \(\nu\cdot\xi=0\). Three closures are required:
   gamma = 0 (no pressure perturbation from compression), gamma > 0
   (compressible, the CAS3D model [2]) and the incompressible limit. The
   perpendicular kinetic norm \(\rho|\xi_\perp|^2\) is kept as a second
   option, because published growth rates use either norm; the norm changes
   growth rates, not the unstable count or marginal points. Toroidal
   harmonics couple within one mode family (n fixed modulo the number of
   field periods).
2. **M1 free boundary.** M0 plus the vacuum energy
   \(\delta W_v=\int|\delta B_v|^2/(2\mu_0)\,dV\) outside the plasma and the
   surface term on the plasma boundary [11,12], with either no wall (exterior
   to infinity), an ideal conducting wall of arbitrary 3D shape, or
   fixed-current coils.
3. **M2 resistive thin wall.** M1 with a thin 3D wall of finite
   conductivity: the plasma stays ideal and the wall currents add a
   first-order time derivative, so growth rates span the ideal (no-wall) to
   resistive-wall time scales [4,6]. The wall time tau_w is the slowest decay
   time of free currents in the wall without plasma.
4. **M3 local criteria.** Mercier and ideal ballooning on each flux surface
   [7,13,14], as diagnostics and as the high-mode-number limit that the
   global result must approach.

M0 to M3 are P1e. Resistive plasma, viscosity, flow and two-fluid terms
(stage P1e.11) linearize the residual of
[nonlinear-resistive-mhd.md](nonlinear-resistive-mhd.md) instead of adding a
second stability model.

### Observables and acceptance tolerances

Frozen with this specification. A change needs the escalation rule in
[AGENTS.md](../AGENTS.md) and must not be made in the same step as an
implementation change. The internal column bounds the estimated discretization
error: the difference between the two finest levels of a refinement sequence
that shows the expected convergence rate, or an a posteriori estimate once one
is available. The reference column bounds the difference to a named
reference. Every tolerance is relative, with the absolute floor given where
the observable can cross zero.

| Id | Observable | Internal | Against references |
|---|---|---|---|
| O1 | number of eigenvalues below -delta, delta = 1e-4 omega_A^2, per mode family (inertia of A + delta M) | exact; floating-point count until stage P1e.9, certified after | exact where the reference states it |
| O2 | growth rate of the most unstable mode, units of omega_A | 1 % | analytic: 1e-4; published codes: 3 % for growth rates above 1e-2 omega_A, smaller ones recorded |
| O3 | marginal point: parameter value (beta, current, boundary coefficient) where the lowest eigenvalue crosses zero | 1 %, floor 1e-3 of the scan range | 2 %, same floor |
| O4 | discrete stable eigenvalues and continuum edges | 1e-3 | analytic: 1e-4; published codes: 1 % |
| O5 | spectral pollution: an eigenvalue that is not in the exact spectrum (where it is known) or that does not converge at the element's rate over three refinement levels | none | none |
| O6 | lowest eigenvalues of the plasma, vacuum and total boundary response matrices W_P, W_V, W_T, harmonics normalized as in [9] | 1e-3 | 2 % against DCON |
| O7 | Mercier criterion and ballooning eigenvalue per surface | 1e-3, floor 1e-3 of the largest magnitude on the profile | 1 %, same floor |
| O8 | resistive-wall-mode growth rate, units of 1/tau_w | 1 % | 3 % |
| O9 | derivatives of O2, O3 and O6 with respect to equilibrium and boundary parameters | tangent and adjoint agree to 1e-10; both agree with central finite differences (complex step on real-arithmetic paths) to 1e-6 | not applicable |
| O10 | wall time and peak memory to reach O2 at 1 % | measured | same hardware as the running comparators |

O1 counts below -delta because the continua accumulate at zero, and at
gamma = 0 the parallel displacements form a neutral subspace. When O7 shows a
locally unstable surface, the continuous problem has infinitely many unstable
modes and O1 grows under refinement; it is then reported with that status and
not as converged. A disagreement with a published or running code beyond
these bands is recorded as unresolved with its suspected cause; tolerances are
not widened to absorb it.

### Formulation decisions

- **Energy source.** One FortSym definition of M0 and M1 generates the
  primal, tangent and adjoint leaves. It is derived from [11,12] in MILCA and
  replaces the externally pinned generator behind
  `src/generated/ideal_mhd_*`. Independent checks: equivalence with the
  intuitive form of delta W [12], proved in FortSym and repeated in SymPy;
  reduction to Newcomb's cylinder Lagrangian [16]; agreement of the potential
  energy with the pinned leaves on the existing tests before the pin is
  removed. The kinetic leaf changes from the perpendicular norm of the pinned
  package ([ideal-xi.md](ideal-xi.md)) to the full norm of M0, with the
  perpendicular norm kept as the option above.
- **Discretization.** Radial finite elements aligned with flux surfaces,
  Fourier harmonics in the two angles, one mode family at a time. The radial
  pair generalizes the hybrid elements of [1,15]: the normal component of xi
  continuous of degree p, the tangential components discontinuous of degree
  p-1, designed so that the divergence and field-line bending terms are
  approximated compatibly and the continua stay unpolluted, as analysed in
  [15] for one- and two-dimensional equilibria; in 3D, O5 checks it. Degree,
  radial mesh and harmonic set adapt, with refinement at rational surfaces.
  Unaligned meshes are out of P1e because the published codes and the
  pollution analysis assume alignment.
- **Equilibria.** The native nested map equilibrium
  ([nested-equilibrium-input.md](nested-equilibrium-input.md)) and imported
  VMEC, GVEC and DESC states through an adapter, both transformed to
  straight field line coordinates. The equilibrium force residual and its
  error enter the error ledger of every stability result.
- **Vacuum.** The exterior energy uses the Steklov-Poincare map of the
  existing boundary-element operators
  ([ideal-free-boundary.md](ideal-free-boundary.md)): exact to infinity,
  with an ideal wall or coils as additional boundaries. The thin resistive
  wall is a surface current potential on the wall, following the
  formulation published in [4].
- **Spectrum.** A real-shift LDL factorization of A - sigma M gives the
  number of eigenvalues below sigma by Sylvester's law of inertia, as in
  [s0-spectrum.md](s0-spectrum.md); shift-invert Lanczos [8] gives the
  eigenpairs. W_P, W_V and W_T are Schur complements on the boundary
  harmonics.
- **Derivatives.** For a simple eigenvalue,
  \(\lambda'=x^H(A'-\lambda M')x/x^HMx\), where the prime is the
  derivative with respect to a parameter mu and A', M' come from the tangent
  leaves, including the shape derivative of the vacuum boundary-element map
  when the boundary moves. Derivatives through the equilibrium use the
  adjoint of the equilibrium residual, not the solver iterations. Marginal
  points are Newton roots of lambda(mu) = 0. Degenerate and crossing
  eigenvalues report a failed status instead of a derivative.
- **Certification.** The discrete claim is an interval bound on the inertia
  count and on isolated eigenvalues. Continuous claims are made only where a
  discretization bound exists; elsewhere the status is certified-discrete
  with an empirical convergence estimate, stated as such.

### P1e implementation status

**P1e.1 (done).** `gen/milca_symbolic_ideal.f90` holds the frozen-in energy
of [ideal-map-3d.md](ideal-map-3d.md); the map, xi, GLL and current-sheet
leaves regenerate from it byte for byte, and the external generator, its
pin and its coordinate-reduced and batched leaves are removed. The
coordinate-reduced leaves had no consumer besides their own test and were
removed without a comparison with the new leaves; the pinch test below
replaces that check. `gen/gen_ideal_dw.f90` emits, with tangent and adjoint
leaves, the Hermitian part of the Eulerian delta W (`dw_bernstein`), the
adiabatic Lagrangian Hessian (`dw_lagrange`) and the free-boundary vacuum and
surface densities
(`dw_vacuum`, `dw_surface`); `milca_ideal_dw` wraps them for real and complex
fields. FortSym proves 85 identities exactly: symmetry of both forms, the 81
entries of the Hessian of the adiabatic frozen-in energy at the identity, the
link between the two forms in a uniform medium, and the shear-Alfven and
magnetosonic characteristic polynomial of the plane-wave matrix. These
proofs do not involve the current and pressure-gradient terms. The test
`ideal_dw_pinch_equivalence_and_derivatives` therefore integrates both forms
for compactly supported fields on a force-balanced screw pinch with uniform
current: they agree to 2e-16 relative, differ by 1.5 % when the pressure
gradient is removed, and differ by 8 % when the sign of the current term is
flipped. The same test checks the gamma = 0 limit against the xi form,
Hermitian symmetry and the tangent and adjoint leaves. Still open from the
energy-source plan: the equivalence with the intuitive form of delta W, the
repetition of the identities in SymPy and the
reduction to Newcomb's Lagrangian, which rung R2 tests numerically.

**P1e.2 (done).** `milca_ideal_cylinder` assembles the stiffness (Eulerian
form) and mass of a periodic screw pinch for one pair exp(i(m theta +
kz z)). The unknowns are the contravariant densities X = r xi_r in the
B-spline space S_p (FortFEM `evaluate_bspline_basis`) and T = xi_theta,
Z = r xi_z in its derivative space S_(p-1); interior knots of multiplicity p
make this the hybrid pair, continuous X and discontinuous T and Z. X vanishes
on the axis and wall, Z on the axis, T on the axis unless |m| = 1. The
physical components come from a generated leaf. Results
(`ideal_cylinder_*`):

| Rung | Case | Result |
|---|---|---|
| R1 | uniform cylinder, m = 2, p = 3, 24 cells | fast and slow roots within 7e-8 (O4); no eigenvalue in the exact gaps (O5); Alfven point exact |
| R2 | fixed-wall m = 1 kink of [s0-spectrum.md](s0-spectrum.md), gamma = 0 | omega^2 within 1.3e-10 of the independent oracle at p = 3, 32 cells, error falling about as h^(2p); one eigenvalue below -1e-4 by sparse LDL^H inertia (O1); shift-invert Lanczos agrees |
| R2 | analytic limit Btheta = r, Bz -> 0, kz = -40, -80, -160 | extrapolated omega^2 = -3 within 1.1e-6 (O2) |

Compression (gamma = 5/3) raises the kink eigenvalue to -0.149, as it must.
The perpendicular kinetic norm is not offered yet: it leaves the parallel
displacement without inertia, which needs the field-aligned components of
P1e.4. Newcomb and Suydam (O7) belong to P1e.7.

**P1e.3 (in progress): general flux-coordinate space.** `milca_ideal_flux`
takes any nested chart x(rho, theta, zeta) with its first and second
derivatives, the toroidal flux per radian psi(rho), iota(rho), the jets of a
stream function lambda (B = grad psi x grad(theta + lambda) + iota grad zeta x
grad psi, so that VMEC, Boozer and PEST coordinates all fit) and p(rho). The
unknowns are the contravariant densities D^u = sqrt(g) xi^u of one mode
family on the hybrid radial pair; `gen/gen_ideal_flux.f90` emits B, grad B,
grad p and the linear map from density jets to (xi, grad xi), and the
coefficient matrix C of the Eulerian form, so that the stiffness integrand is
J^H T^T C T J. Harmonics couple through the angular Fourier coefficients of
T^T C T on a uniform angle grid. Checks (`ideal_flux_*`): on a
non-axisymmetric chart with lambda, B agrees with an independent dual-basis
formula to 1e-13, grad B and grad xi with central differences to 1e-7, and
div B = 0 and sqrt(g) div xi = d_u D^u hold to rounding (the general symbolic
proof of these two identities did not finish in 9 minutes and is replaced by
this check). On a straight-cylinder chart the m = 1 kink agrees with the
independent oracle to 1.3e-10, also with m = 0, 1, 2 coupled, and the uniform
cylinder with m = 1, 2, 3 gives the union of the single-harmonic spectra to
2e-11.

Toroidal input: `milca_ideal_flux_tokamak` turns a poloidal-flux grid
(G-EQDSK) into a chart, with rho = sqrt(psin/psin_wall), the PEST or Hamada
straight-field-line angle of `milca_vacuum_surfaces`, or the geometric angle
with lambda = theta_PEST - theta. R, Z, the toroidal shift and lambda are
Fourier series whose coefficients, divided by rho^k with the parity of m
(k <= 2), are fitted by quintic B-splines in s = rho^2. On an exact Solov'ev
equilibrium (`ideal_flux_tokamak_*`, elongation 1.5, q from 0.89 to 1.15,
fixed wall at psin = 0.8) all three charts reproduce the analytic field to
4e-7 and force balance to 3e-3 (straight charts 5e-4); the chart jets agree
with central differences to 2e-9. The fixed-boundary m = 1 internal kink
converges, with m = -10..28 and 32 cells, to omega^2 = -1.30328e-3 in PEST
and -1.30331e-3 in Hamada coordinates (2.6e-5 apart; the geometric chart
gives -1.294e-3 with the same harmonics, its basis converging more slowly);
the gate uses m = -8..18 and 16 cells, where the two straight charts agree to
9e-4. Force balance converges as h^3.5 in the spacing of the flux grid,
which therefore sets the equilibrium error of EQDSK input. Still open in
P1e.3: published PEST cases and DCON on the same input (rung R3), Boozer-file
input, and the sparse complex LDL^H of FortSparse, which takes 95 s for an
11k-row banded matrix with 4.7M nonzeros and limits the harmonic range of
the tests. Its real symmetric embedding [[Re, -Im], [Im, Re]] factors in
20 s with the real front-based LDL^T and twice the inertia; dense-front
complex factorization in FortSparse is part of P1e.10.

**P1e.5 (in progress): free boundary.** In the cylinder the edge density
X(a) is free and the vacuum adds 2 delta W_V = -(a/mu) F_a^2 G |xi_r(a)|^2 per
unit 4 pi^2 R0, with G = phi/phi' at the edge from modified Bessel functions
(conducting wall at b, or decay to infinity). On the uniform-current case of
[ideal-free-boundary.md](ideal-free-boundary.md#independent-cylindrical-external-kink)
(`ideal_cylinder_free`, q = 0.5): the rigid m = 1 displacement, exactly
representable in the hybrid space, has the analytic plasma-plus-vacuum
quotient q^2 - 1 + f (q - 1)^2 = -0.5060096 to 8e-13; the lowest eigenvalue
lies below it (-0.52682 at ka = 0.1, the Ritz bound) and tends to
2q(q - 1) = -0.5 as ka -> 0 (8e-6 at ka = 1e-3, O2). A conducting wall
stabilizes monotonically as it approaches, and the mode switches from
unstable to stable between b/a = 1.43 and 1.40 around the analytic marginal
wall sqrt(2) (O3). The toroidal vacuum (boundary elements) and W_P, W_V, W_T
follow.

### Verification ladder

Each rung is a frozen gate; a later rung does not start before the earlier
one passes or its failure is recorded.

| Rung | Case | Reference | Observables |
|---|---|---|---|
| R0 | algebraic identities of delta W: self-adjointness, equivalence of the two forms, the three gamma closures | FortSym proof, repeated in SymPy | - |
| R1 | homogeneous slab and cylinder | exact Alfven and slow continua and fast modes [12] | O4, O5 |
| R2 | screw pinch | Newcomb criterion [16], Suydam [12], external kink of the uniform-current cylinder, all already in [linear-mhd-oracles.md](linear-mhd-oracles.md) and [ideal-free-boundary.md](ideal-free-boundary.md) | O1, O2, O7 |
| R3 | axisymmetric torus | Bussac internal kink (existing oracle), published PEST cases [17], DCON in GPEC run on the same input [R2]; toroidal Alfven gap and its eigenmodes | O1, O2, O4, O5, O6, O7 |
| R4 | low-shear stellarator, low-n and Mercier modes | TERPSICHORE published cases [7] and runs on the same VMEC input [R1] | O1, O2, O7 |
| R5 | W7-X-type configurations | CAS3D published results [1,2,3] where the published inputs suffice | O1, O2, O3, O4 |
| R6 | 3D tokamak, vertical and resistive wall modes | CASTOR3D and STARWALL published results [4,5,6] | O2, O8 |
| R7 | compact quasi-axisymmetric stellarator | published TERPSICHORE results [18]; DCON3D [10] if P1e.0 finds it available | O1, O2, O6 |

A case whose published input data are incomplete is listed as not comparable,
not reconstructed from figures by fitting.

### Stages

| Stage | Content | Exit gate | Needs |
|---|---|---|---|
| P1e.0 | confirm every S and ? entry above from full text; collect published inputs for R3 to R7; establish DCON3D availability | matrix holds only P and full-text entries | - |
| P1e.1 | MILCA-owned delta W (M0, M1) in FortSym; retire the external pin | R0; potential energy agrees with the pinned leaves on the existing tests, then the pin is removed | P1e.0 |
| P1e.2 | radial hybrid element times Fourier space, fixed boundary, one mode family, with the existing shift-invert eigen solver | R1; R2 fixed boundary (O1 as floating-point count, O2) | P1e.1 |
| P1e.3 | axisymmetric toroidal equilibria in straight field line coordinates; stable spectrum with gap targeting | R3 fixed boundary (O1, O2, O4, O5) | P1e.2 |
| P1e.4 | 3D equilibria, native and imported; mode-family coupling | R4, R5 fixed boundary | P1e.3 |
| P1e.5 | free boundary: vacuum to infinity, ideal wall, coils; W_P, W_V, W_T | R2 external kink; R3 O6 against DCON; R5 free boundary and R7 | P1e.3; P1e.4 for R5 and R7 |
| P1e.6 | thin resistive 3D wall | R6 | P1e.4, P1e.5 |
| P1e.7 | Mercier and ballooning; global-to-local limit at high n | O7 on R2 and R3; on R4 | P1e.2; P1e.4 for R4 |
| P1e.8 | derivatives and marginal points | O9 on R2 and R3; O3 as beta scans on R3; R5 | P1e.3; P1e.4 and P1e.5 for R5 and free-boundary derivatives |
| P1e.9 | interval certificates for inertia and isolated eigenvalues | O1 certified on R2 and R3; on R4 and R5 | P1e.3; P1e.4 for R4 and R5 |
| P1e.10 | performance on the clusters at 1 % on R3 to R5; CUDA optional | O10 against TERPSICHORE and DCON at matched accuracy | P1e.4, P1e.5 |
| P1e.11 | resistive plasma, viscosity, flow, two-fluid terms (CASTOR3D rows) | published CASTOR3D cases [5] | P1e.6 and a linearization of the resistive MHD residual |

Library work follows step 3 of the P1e procedure in
[ROADMAP.md](../ROADMAP.md): interval LDL inertia and shift-invert Lanczos in
FortSparse and FortNum, the radial hybrid element in FortFEM, the straight
field line transform in libneo if it owns one, else in the MILCA adapter.

### Sources

Publications (S and ? entries above are still to be read in full):

1. C. Schwab, Phys. Fluids B 5, 3195 (1993).
2. C. Nührenberg, Phys. Plasmas 6, 137 (1999).
3. C. Nührenberg, Nucl. Fusion 56, 076010 (2016),
   doi:10.1088/0029-5515/56/7/076010.
4. P. Merkel and E. Strumberger, arXiv:1508.04911 (2015).
5. E. Strumberger and S. Günter, Nucl. Fusion 57, 016032 (2017),
   doi:10.1088/0029-5515/57/1/016032.
6. S. Günter, E. Strumberger, P. Merkel and M. Sempf, Phys. Plasmas 15,
   056110 (2008), doi:10.1063/1.2884579.
7. W. A. Cooper, U. Schwenn, G. Y. Fu, R. Gruber and D. V. Anderson,
   Phys. Fluids B 4, 1401 (1992).
8. N. Nakajima, M. Okamoto and J. Chen, Comput. Phys. Commun. 113, 1 (1998),
   doi:10.1016/S0010-4655(98)00078-2.
9. A. H. Glasser, Phys. Plasmas 23, 072505 (2016), doi:10.1063/1.4958328.
10. A. H. Glasser, Phys. Plasmas 32, 082504 (2025), doi:10.1063/5.0274154.
11. I. B. Bernstein, E. A. Frieman, M. D. Kruskal and R. M. Kulsrud,
    Proc. R. Soc. A 244, 17 (1958).
12. J. P. Goedbloed and S. Poedts, *Principles of Magnetohydrodynamics*,
    Cambridge University Press (2004); J. P. Freidberg, *Ideal MHD*,
    Cambridge University Press (2014).
13. C. Mercier, Nucl. Fusion Suppl. 2, 801 (1962).
14. J. W. Connor, R. J. Hastie and J. B. Taylor, Phys. Rev. Lett. 40, 396
    (1978).
15. R. Gruber and J. Rappaz, *Finite Element Methods in Linear Ideal
    Magnetohydrodynamics*, Springer (1985).
16. W. A. Newcomb, Ann. Phys. 10, 232 (1960).
17. R. C. Grimm, R. L. Dewar and J. Manickam, J. Comput. Phys. 49, 94
    (1983), doi:10.1016/0021-9991(83)90116-X.
18. W. A. Cooper et al., Comput. Phys. Commun. 141, 55 (2001),
    doi:10.1016/S0010-4655(01)00396-4.
19. D. V. Anderson, W. A. Cooper, R. Gruber, S. Merazzi and U. Schwenn,
    in *Scientific Computing on Supercomputers II*, Springer (1990),
    doi:10.1007/978-1-4613-0659-7_8.

Open-source comparators, to be studied and run, no code taken from them. The
commit is the one inspected on 2026-11-12; each comparison run replaces it
with the commit it ran.

- R1. TERPSICHORE, Apache-2.0, https://gitlab.epfl.ch/spc/public/terpsichore,
  master at 04dcf9af3ceea7ef36198508f2beb5fafd9cef20.
- R2. GPEC (Fortran, includes DCON), MIT,
  https://github.com/PrincetonUniversity/GPEC, develop at
  02f90d4a432a3e7f8586dbd1d22f0c2e84db2c20.
- R3. GPEC (Julia), MIT, https://github.com/OpenFUSIONToolkit/GPEC, develop
  at 052978b42404e53b2ba80d433ede3680f87fe995.

STARWALL is used only through [4]; its code is not public.

Chris&AI
