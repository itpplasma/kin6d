# Numerical-method catalogue

Status: canonical numerical option map. This is a library of admissible
approaches, not a requirement to implement all methods before science work.

## Principle

KIN6D chooses numerics to preserve the structure and asymptotics of the selected
physical model at the lowest time/memory cost for the requested observable.
No single discretization is optimal for every plasma regime.

The primary high-order kinetic path remains mapped tensor/FEEC/spectral where
it is advantageous, but the common framework must also admit shock/front,
rarefied, stochastic and strongly correlated regimes.

## Spatial methods

Supported design families include:

- mapped tensor-product B-spline/NURBS/IGA;
- continuous finite elements;
- FEEC-compatible H1/H(curl)/H(div)/L2 complexes;
- DG/HDG;
- Fourier/spectral/pseudospectral methods;
- finite-volume conservative methods;
- entropy-stable/high-resolution DG or finite-volume shock capturing;
- Riemann-solver, limiter and WENO-like reconstructions when hyperbolic shocks
  require them;
- cut/unstructured meshes and mortars;
- FCI/field-line-map operators;
- BEM/DtN/PML/open-boundary methods.

Use the simplest structure-preserving method that matches the PDE. Do not force
shock-dominated compressible flow into the kinetic tensor path or Maxwell FEEC
into a generic finite-volume representation merely for uniformity.

## Mesh and representation adaptivity

Admit:

- h-, p- and hp-adaptivity;
- block/patch AMR;
- velocity-domain and harmonic adaptation;
- local representation/model switching with conservative transfer;
- moving/ALE meshes when a concrete moving-boundary problem requires them.

AMR is especially relevant to streamers, sheaths, shocks, current sheets and
localized reaction/radiation fronts.

Goal-oriented error indicators should eventually rank adaptation by observable
benefit per cost.

## Kinetic continuum methods

Candidate families include:

- discrete velocity;
- spherical-harmonic/polynomial velocity bases;
- Hermite--Laguerre/gyromoment expansions;
- conservative DG in phase space;
- semi-Lagrangian/characteristic Galerkin;
- spectral collision operators;
- sparse grids/low-rank/tensor decompositions when they demonstrate real
  5D/6D savings.

Low-rank/tensor-train formats are optional representations with explicit
compression defects, not a universal state format.

## Particle and stochastic methods

The common particle infrastructure should support:

- deterministic characteristics/orbit following;
- PIC/charge-conserving geometric PIC;
- PIF/particle-in-basis;
- delta-f/control-variate particles;
- blobs/remapping;
- Monte-Carlo collision (MCC)/binary-collision operators;
- DSMC-like rarefied neutral gas;
- photon Monte Carlo;
- full relativistic Maxwell--Lorentz finite-(N) reference calculations.

Particle sorting, binning, migration and event handling are reusable runtime
services only after real consumers establish the abstraction.

## Long-range microscopic interaction methods

For physical-particle microscopic references, retain:

- direct \(O(N^2)\) force evaluation as the tiny gold standard;
- FMM/tree methods for controlled long-range pair interactions;
- Ewald/spectral-Ewald and P3M for periodic Coulomb systems;
- near/far decompositions that keep close encounters explicit while
  accelerating smooth distant fields;
- local full-Maxwell field evolution for large-\(N\) relativistic
  electrodynamics;
- direct retarded Lienard--Wiechert mutual-field evaluation for small-\(N\)
  reference cases;
- space-time/hierarchical retarded kernels only if they provide measured
  benefit at matched error.

The interaction matrix is logical and must not be stored densely by default.
FMM truncation tolerance is numerical error; macro-particle reduction is a
different physical model.

## Quantum numerical methods

Quantum plasma is not particle-only. The appropriate state determines the
numerics.

### Exact many-body state

A direct wavefunction \(\Psi(x_1,\ldots,x_N)\) or density matrix is a tiny-\(N\)
reference. In 3D its storage grows exponentially with particle number.
Use symmetry, antisymmetry and basis truncation, but do not present it as a
macroscopic solver.

### Reduced density matrices and correlation tensors

Quantum BBGKY/cumulant models may use low-rank separated forms, tensor trains,
hierarchical Tucker, sparse grids or center-of-mass/relative coordinates.
Every truncation/compression carries an explicit defect/status.

### Orbital and density methods

Hartree, Hartree--Fock, Kohn--Sham DFT and time-dependent variants evolve
ordinary 3D orbitals/densities. Candidate discretizations include:

- FEM/IGA;
- plane waves/Fourier;
- finite differences;
- spectral/wavelet bases.

These methods therefore reuse ordinary mesh, eigensolver, Poisson/Maxwell and
linear-algebra infrastructure rather than creating a particle-only quantum
runtime.

### Quantum Monte Carlo/path integrals

Variational/diffusion/path-integral Monte Carlo sample many-body
configuration/path space. Report statistical uncertainty, finite-size error,
imaginary-time discretization and sign-problem limitations separately.

### Wigner continuum

Wigner--Poisson/Maxwell is a phase-space continuum problem. Candidate methods
include Eulerian spectral/DG, semi-Lagrangian, tensor-product, sparse-grid and
dynamic low-rank representations. Its \(\hbar\to0\) Vlasov limit is a permanent
cross-model oracle.

### Quantum fluids and EOS closures

Fermi pressure, exchange/correlation and derived quantum-gradient terms live in
ordinary physical-space PDEs. FEM/FEEC/FV/DG can therefore propagate quantum
closure effects across macroscopic devices or stellar interiors without
resolving Fermi wavelengths.

### Planning scale envelope

These are design expectations, not current performance claims:

- exact 3D many-body wavefunction: typically 2--few particles;
- exact diagonalization: tiny basis-dependent systems;
- QMC/PIMC: tens--thousands of microscopic particles, strongly problem
  dependent;
- HF/DFT/TDDFT: roughly \(10^2\)--\(10^4\)-scale particles/orbitals as a
  plausible reference range depending on basis/hardware;
- Wigner: local kinetic boxes, with full 3X3P requiring the same tensor/low-rank
  discipline as classical FK;
- quantum-fluid/EOS: macroscopic whole-device or whole-star domains.

See [quantum-plasma.md](quantum-plasma.md).

## Time integration

The method catalogue includes:

- explicit Runge--Kutta/SSP methods;
- implicit BDF/DIRK/fully implicit methods;
- IMEX schemes;
- exponential/integrating-factor methods;
- operator splitting;
- multirate and subcycling schemes;
- symplectic/variational/geometric integrators;
- space--time methods when useful.

The timestepper is selected from the model's stiffness, invariants and error
target; it does not define the physical equations.

## Asymptotic-preserving and multiscale numerics

Scale-separated plasma models should support asymptotic-preserving or uniformly
accurate discretizations where they avoid resolving irrelevant fast scales
while retaining the correct limit.

Examples of relevant limits include:

- quasineutral/Debye limits;
- fluid limits of kinetic equations;
- highly collisional diffusion limits;
- drift/gyro reductions;
- radiation diffusion limits;
- stiff reaction equilibria.

AP behavior is a mathematical property to verify against both the parent and
limit models, not a label on an implicit solver.

## Nonlinear, linear and spectral solvers

The solver toolbox may include:

- matrix-free Newton--Krylov;
- nonlinear multigrid/FAS;
- geometric/algebraic multigrid;
- block/domain-decomposition/Schwarz preconditioning;
- static condensation/Schur complements;
- low-storage Krylov and pipelined variants;
- Arnoldi/Lanczos/Jacobi--Davidson;
- shift-invert/resolvent actions;
- continuation and bifurcation algorithms.

Physics-based reduced models are preferred preconditioners when they lower
total time and memory without changing the target residual.

## Optimization, inverse and uncertainty numerics

KIN6D supplies differentiable forward models and numerical/model uncertainty.
The wider workflow may use:

- tangent and adjoint optimization;
- deterministic interval/ball validation;
- Monte Carlo/QMC/randomized QMC;
- stochastic collocation/polynomial-chaos-like methods where useful;
- ensemble methods;
- reduced-order/surrogate models with explicit domain and error status.

TIAGO owns diagnostic likelihood/posterior inference in the current
architecture. Other inference/data-assimilation consumers must use the same
state/observable/JVP/VJP contract rather than redefine plasma physics.

## Reduced-order and surrogate acceleration

ROMs, response surfaces and learned surrogates are allowed as accelerators or
closures if they expose:

- training/calibration domain;
- conservation/structure constraints where relevant;
- derivative behavior;
- out-of-domain detection;
- a validation/error status.

A surrogate never silently replaces the named target model in a certified
claim.

## Radiation numerics

Radiation uses the same generic numerical families plus, where appropriate,

- diffusion/multigroup block solvers;
- moment closures;
- characteristic/discrete-angle transport;
- Monte-Carlo photons.

See [radiation.md](radiation.md).

## Performance selection

Every implemented method keeps a small independent reference and a realistic
performance case. Selection is empirical after mathematical admissibility:
time, peak memory, derivative cost and certificate cost all matter.

The permanent optimization order remains

\[
\text{model}
\rightarrow
\text{representation}
\rightarrow
\text{algorithm}
\rightarrow
\text{data movement}
\rightarrow
\text{parallel decomposition}
\rightarrow
\text{microarchitecture}.
\]

Chris&AI
