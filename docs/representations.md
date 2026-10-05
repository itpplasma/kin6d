# Distribution representations and particle methods

Status: canonical representation policy. The compact continuum laboratory
includes free streaming, Krook balance, single-mode Poisson, a coupled
source-balanced Poisson--FP family and bounded OU/FP equilibrium. Particle and
other carrier implementations remain staged capabilities. The current
continuous bounds cover declared analytical cases and arbitrary represented
sources for [damped transport](general-transport.md). General nonlinear
kinetic certification remains planned.

## Principle

KIN6D separates the **physical kinetic model** from the **finite-dimensional
representation** used to approximate its distribution. A full-kinetic,
gyrokinetic or drift-kinetic species must not be tied permanently to one carrier
such as an Eulerian mesh or PIC particles.

The principal independent axes are

\[
\boxed{
\begin{array}{ll}
\text{distribution form} & \text{full-f, delta-f, control variate},\\
\text{carrier} & \text{Eulerian, semi-Lagrangian, particles, blobs},\\
\text{field coupling} & \text{grid/FE, PIC, PIF, particle-in-basis},\\
\text{sampling} & \text{deterministic, QMC, randomized QMC, MC},\\
\text{collision realization} & \text{deterministic operator, stochastic process}.
\end{array}}
\]

A configuration is accepted only if these choices implement the same named
continuous model and have compatible measures, fields, boundaries and
conservation laws.

## Accuracy and result contract

Choose a named observable and accuracy scale before comparing carriers. The
practical workhorse target is approximately 1%; near-zero observables need an
absolute scale. The goal is the fastest measured total time at matched model,
accuracy and claim strength, including setup, solve, certification and requested
derivatives. Record memory and measurement variability as well as runtime.

Every reported state, field, observable and sensitivity carries an applicable
bound/status or explicitly states that its bound is unavailable. Account for
representation, quadrature, field, timestep, boundary and rounding errors as
applicable. A validated coarse rounding budget may suffice for the requested
accuracy; deterministic containment is required for a certified result, while
tighter audits are optional. Statistical and physical-model uncertainties remain
separate.

Plot full applicable intervals. An L2 state radius does not support pointwise
bands, and a sum of Dirac particles has no ordinary pointwise distribution
error. Use weak-observable bounds or name the smooth reconstruction and bound
its own error. A finite algebraic inclusion does not enclose the continuous
solution until truncation, reconstruction and stability premises also close.

## Quantum state representations

Quantum plasma uses a different state ontology where appropriate; it is not
forced into a classical distribution carrier.

Possible microscopic/continuum quantum states include:

- many-body wavefunction \(\Psi(x_1,\ldots,x_N)\) or density matrix \(\rho_N\)
  for tiny-\(N\) exact references;
- reduced density matrices/cumulants \(\rho_1,\rho_2,\ldots\);
- Hartree/Hartree--Fock/Kohn--Sham orbitals on ordinary 3D meshes/bases;
- QMC/PIMC walkers or imaginary-time paths for stochastic equilibrium
  references;
- Wigner quasidistributions \(f_W(x,p)\) on Eulerian/semi-Lagrangian
  phase-space representations;
- quantum-fluid/Fermi-EOS state variables on physical-space meshes.

These representations are not interchangeable by storage conversion alone.
Each corresponds to a named physical/statistical model with its own measures,
symmetry/antisymmetry, positivity or quasiprobability semantics and observable
map.

The same separation principle still applies: physical model and numerical
carrier are distinct. A Wigner model may be spectral, DG or semi-Lagrangian;
an orbital model may use FEM, plane waves or finite differences.

See [quantum-plasma.md](quantum-plasma.md).

## Full-f, delta-f and control variates

The full distribution may be represented directly,

\[
f=f_h,
\]

or split,

\[
f=F_0+\delta f.
\]

The reference \(F_0\) may be an analytic Maxwellian, imported equilibrium,
previous stationary state, fluid/gyrofluid reconstruction, DK/GK solution or
another explicitly named approximation.

Delta-f is not intrinsically a particle method. Both \(F_0\) and \(\delta f\)
may use continuum or particle carriers. The split must record which moments,
fields, collisions and sources act on each piece and how conservation is
restored or shared.

A reduced model may also act as a control variate for a more expensive model.
This changes estimator variance or solver conditioning, not the target
equations. The residual of the target model remains authoritative.

## Eulerian and semi-Lagrangian continuum

Eulerian representations store coefficients in a basis/grid and apply
transport, fields and collisions through weak/strong operator actions.

Semi-Lagrangian/characteristic methods follow the phase flow and reconstruct
the distribution. They are attractive for low-noise advection over large
timesteps, but are accepted only with explicit control of:

- mass/charge conservation;
- phase-space measure;
- interpolation/projection defects;
- positivity where claimed;
- long-time Hamiltonian/free-energy structure.

The project does not identify "semi-Lagrangian" with one interpolation scheme.
Conservative DG, characteristic Galerkin, variational maps and exact subflows
are legitimate realizations.

The stationary certified residual remains a mathematical differential/weak
problem. A time-marching interpolation algorithm does not define the steady
physics.

## Particles

A particle carrier represents

\[
f_h(z)=\sum_{p=1}^{N_p} w_p\delta(z-Z_p)
\]

or a regularized approximation. Particle coordinates may be chosen
deterministically, by QMC, or randomly.

A particle method is not automatically PIC and not automatically Monte Carlo.
Pure characteristic particles may never deposit to a field. They are useful
for orbit maps, loss observables, quadrature and validation and naturally
generalize SIMPLE/GORILLA-like workloads.

Particle storage and movement should be structure-of-arrays/AoSoA, tileable,
GPU-resident when accelerated, and migration-aware when distributed. A generic
particle runtime is not split out of KIN6D until repeated independent consumers
show which abstraction is actually stable.

### Physical particles versus superparticles

For the microscopic Maxwell--Lorentz reference, one represented particle is a
physical particle with its actual species charge and mass. Replacing many real
particles by one macro-particle changes finite-\(N\) noise, close encounters and
correlations.

FMM/tree, Ewald/spectral-Ewald and P3M are therefore preferred microscopic
accelerations: they compress the **interaction evaluation**, not the physical
particle state. A blob or PIF marker may still be used in the microscopic
branch only when its shape/weight is part of the explicitly modified model.

## PIC

Particle-in-cell couples particle distributions to local field bases. The
preferred long-range design is variational/geometric PIC with compatible FEEC
field spaces, so charge continuity and gauge structure arise from the same
discrete complex used by continuum Maxwell clients.

Required components include:

- finite particle shapes;
- charge/current deposition;
- charge-conserving current update;
- field gather;
- compatible Maxwell/Poisson solve;
- particle pusher;
- boundary/surface interaction;
- sorting/binning and distributed migration.

PIC is valuable both as a production representation in high-dimensional
collisionless regimes and as an independent oracle for continuum methods.
Sampling noise and grid/shape error must be reported separately.

The first teaching comparison uses only periodic 1X1V electrostatic
deposit/Poisson/gather. Declare its particle shape, pairing and tested
conservation properties. Maxwell current continuity, geometric long-time
structure, walls and distributed migration remain later admission gates rather
than claims of the initial baseline.

## PIF and general particle-in-basis

For Fourier field modes,

\[
\rho_k=\sum_p q_p w_p\widehat S(k)e^{-ik\cdot X_p},
\]

and the spectral field is solved and evaluated directly. This is
particle-in-Fourier (PIF).

The general abstraction is

\[
\rho_a=\sum_p q_pw_p\langle S_p,\psi_a\rangle
\]

for field basis functions \(\psi_a\). Fourier basis gives PIF; local FE/grid
functions give PIC/particle-FE.

PIF is especially attractive when only a moderate number of physical modes is
required. Direct \(O(N_pN_k)\) coupling and NUFFT-accelerated variants should
both be benchmarked. Never materialize the logical particle-mode matrix;
generate tiled map/reduce kernels that accumulate modes or particle forces
while tiles are resident.

## Blobs and moving local bases

A finite blob uses

\[
f_h(z)=\sum_p w_p S_p(z-Z_p;C_p),
\]

where \(C_p\) may encode width, covariance or orientation. Fixed PIC shapes are
one special case; adaptive Gaussian/compact packets are another.

Blobs are interesting because they interpolate between delta particles and a
smooth continuum representation. Their evolution, remapping and overlap must be
derived rather than treated as arbitrary smoothing.

Potential capabilities include:

- adaptive anisotropic covariance;
- deterministic moving quadrature;
- reconstruction/remapping when distortion grows;
- particle-to-continuum and continuum-to-particle transfer;
- hybrid smooth bulk plus sparse energetic tail.

Any compression/remapping step carries an explicit defect or statistical
status.

## Rarefied-particle and direct-correlation models

DSMC-like neutral/gas particles and direct classical N-body/MD particles reuse
parts of the same particle runtime but represent **different continuous/
microscopic models** from ordinary Vlasov PIC.

- DSMC approximates a rarefied Boltzmann collision process.
- MCC/binary-collision particles approximate selected stochastic collision/
  reaction operators.
- Relativistic Maxwell--Lorentz finite-(N) resolves explicit physical
  particles and their electromagnetic interaction. It serves both weak-coupling
  closure validation and strongly coupled classical reference problems.

Do not claim equivalence to FK+Landau merely because all three use arrays of
particles. Their convergence targets, invariants and statistical/correlation
errors differ.

## Deterministic sampling, QMC and Monte Carlo

Deterministic quadrature remains the preferred reference where practical.
Low-discrepancy/QMC sampling is important for high-dimensional integrals,
particle initialization and loss observables. Randomized QMC adds statistical
error estimation while retaining low discrepancy.

Ordinary Monte Carlo is first-class when appropriate, notably for:

- kinetic neutral histories and DSMC-like rarefied gas;
- stochastic reaction/wall events;
- binary particle collisions;
- rare-event and importance-sampling estimators;
- problems where full continuum resolution is prohibitively dimensional.

The rule is not "no MC". It is **no uncontrolled numerical noise**.

A deterministic certificate has a deterministic inclusion statement. A
stochastic calculation may instead provide a confidence statement

\[
\Pr(Q\in I)\ge 1-\delta.
\]

The two are never merged silently. Hybrid results state both deterministic and
statistical contributions.

## Stochastic collisions and reactions

Particle clients may use stochastic binary-collision or reaction algorithms.
These are alternative realizations of a specified collision/reaction model, not
automatically equivalent to the deterministic Landau/Beliaev--Budker operator.

Where a stochastic scheme approximates a Fokker--Planck limit, validate:

- conservation in expectation and/or per realization as claimed;
- equilibrium distribution;
- relaxation rates;
- variance and timestep effects;
- convergence to the named deterministic operator.

The deterministic operator remains especially important for low-noise
transport coefficients, Onsager checks and rigorous numerical certification.

## Hybrid representations

Different populations may use different carriers. Examples include:

- continuum thermal bulk + particle energetic alphas;
- gyrofluid/DK control variate + delta-f particle correction;
- deterministic neutral continuum in dense regions + MC histories in dilute
  long-mean-free-path regions;
- continuum GK bulk + FK particle tail.

A spatial switch between representations requires a conservative transfer or
overlap formulation, not a mask applied to unrelated states. Transfer error,
sampling error and model-reduction error are distinct.

## Local and temporal representation adaptation

Representation may vary by species, region, phase-space subset and time when a
derived transfer exists.

Examples include

- continuum bulk with particle energetic tail;
- deterministic neutral continuum with MC/DSMC rarefied regions;
- Eulerian kinetic state with a transient characteristic/particle patch;
- physical-particle microscopic near region coupled to a field/FMM far
  representation;
- classical state with a temporary quantum close-encounter microproblem.

A representation switch does not by itself change the named physical model;
however, macro-particle weighting, smoothing/blobs or reduced quantum/classical
state maps may change physical fluctuation/correlation content and must be
classified correctly.

Every transfer declares conserved moments, reconstruction assumptions,
statistical semantics, derivative behavior and a transfer defect. Demotion
cannot silently discard unresolved information. See
[adaptive-multiscale-uq.md](adaptive-multiscale-uq.md).

## Differentiability

Deterministic particle and blob equations are differentiated like any other
finite-dimensional variational system. PIF derivatives are especially simple
because particle-mode factors are analytic and can be generated by FortSym.

Generate primal/JVP/VJP from the first implementation, including particle
weights/positions, deposit/gather and implicit field response. Declare shape
regularity, cell boundaries and event surfaces, and verify derivatives on their
stated domain. The physical/discrete pairing defines the adjoint. Derivative
values carry their own error/status rather than inheriting the state bound.

For random algorithms, derivative semantics are declared explicitly:

- pathwise/reparameterized derivatives with fixed random variables;
- common-random-number finite differences;
- score/likelihood-ratio estimators;
- nondifferentiable event surfaces.

A noisy stochastic gradient is not interchangeable with the deterministic
implicit adjoint of the underlying continuum model. Optimization should prefer
the deterministic/control-variate path when it is available and validate with
stochastic simulations separately.

## Performance

The particle path follows the same performance principles as continuum
operators:

- no per-particle allocation;
- batch/tile interactions;
- sorting/binning only when reuse exceeds its cost;
- persistent device state;
- bounded migration buffers;
- compact 32-bit local metadata where safe;
- direct-vs-NUFFT/PIC crossover measured rather than assumed.

KeOps-style pairwise or particle-mode expressions are logical reductions. A
large interaction matrix is not stored.

Establish reusable upstream benchmark support before new hot paths. KIN6D
supplies fixed physical cases; generic basis mathematics stays in FortFEM,
generated mathematics in FortSym, and rigorous arithmetic/reductions in
FortNum. Extract a generic particle runtime only after repeated consumers
demonstrate a stable abstraction. Compare matched total result time, including
certification and requested derivative work, rather than only pusher speed.

## Verification ladder

Begin with small cases before the full nonlinear Vlasov--Poisson ladder:

1. Exact free streaming compares continuum and deterministic characteristic
   particles against analytical density modes and sensitivities. This is a
   transport comparison; particles alone are not PIC/PIF.
2. Neutral single-mode Poisson compares continuum projection, direct PIF and
   simple local-basis PIC against analytical field/energy. Account separately
   for quadrature, particle shape and field-space error. This coupling test does
   not establish a self-consistent trajectory.
3. Short self-consistent Vlasov--Poisson compares all three with the same
   initial distribution, physical parameters, boundaries and observation times.
   Declare the actual trajectory-bound status: certified individual field
   solves, conservation diagnostics and refinement agreement are insufficient
   for a continuous trajectory certificate.

Use deterministic sampling first, then add randomized variants with separate
statistical assumptions and uncertainty. Semi-Lagrangian methods, blobs,
delta-f/control variates and geometric/current-conserving PIC follow concrete
cases. Defer NUFFT, adaptive particle/remapping machinery, stochastic collisions
and GPU/MPI until the simple comparisons establish a need.

Compare physical observables, sensitivities, invariant/exchange ledgers and
appropriate bounds before ranking total time and memory at matched accuracy.
Later dispersion/growth cases retain the same shared physical operators and
earlier independent oracles. Cross-representation agreement is validation
evidence; its error status is stated explicitly.

Chris&AI

See [the numerical-method catalogue](numerical-methods.md) for shock/front,
AMR, DSMC, MD, radiation and low-rank methods that complement this distribution-focused document.
