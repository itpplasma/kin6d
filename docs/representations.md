# Distribution representations and particle methods

Status: canonical representation policy. Only the current compact 1X1V
continuum slice is implemented; the remaining representations are staged
capabilities.

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

## Deterministic sampling, QMC and Monte Carlo

Deterministic quadrature remains the preferred reference where practical.
Low-discrepancy/QMC sampling is important for high-dimensional integrals,
particle initialization and loss observables. Randomized QMC adds statistical
error estimation while retaining low discrepancy.

Ordinary Monte Carlo is first-class when appropriate, notably for:

- kinetic neutral histories;
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

## Differentiability

Deterministic particle and blob equations are differentiated like any other
finite-dimensional variational system. PIF derivatives are especially simple
because particle-mode factors are analytic and can be generated by FortSym.

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

## Verification ladder

The first cross-representation suite should use identical 1X1V
Vlasov--Poisson initial data and compare:

1. Eulerian continuum;
2. semi-Lagrangian continuum;
3. deterministic PIF;
4. PIC;
5. blobs;
6. delta-f variants;
7. MC/QMC variants where meaningful.

Compare dispersion/growth rate, fields, invariants, observables, memory, runtime
and statistical noise. Agreement between methods with very different dominant
errors is a major validation asset.

Chris&AI
