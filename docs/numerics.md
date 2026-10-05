# Numerical foundation

Status: canonical numerical architecture. FortFEM already contains substantial
tensor-product spline/FEEC machinery, but KIN6D does not yet have the final
matrix-free/tiled production backend described here.

## Ownership

FortFEM owns generic spatial mathematics. KIN6D owns plasma-specific equations,
constitutive laws and composition.

The preferred high-order kinetic path in this document is not the only
admissible numerical family. Shock/front, rarefied, stochastic, radiation and
strong-correlation regimes may require different methods. The complete option
map is [the numerical-method catalogue](numerical-methods.md); implement only
the method families demanded by prioritized vertical slices.

The intended spatial stack includes:

- mapped tensor-product B-splines/NURBS;
- FEEC-compatible \(H^1\), \(H({\rm curl})\), \(H({\rm div})\), \(L^2\) spaces;
- DG/HDG and broken spaces;
- Fourier directions where genuinely periodic;
- multipatch/mortar coupling;
- FCI parallel operators;
- cut/unstructured fallback where geometry requires it;
- FEM/BEM/DtN/PML boundary machinery.

Do not reimplement generic spatial kernels in KIN6D.

## Hyperbolic, shock and front-capable spatial path

Strongly compressible MHD/fluid, collisionless-shock fluid closures,
ionization fronts and some HEDP/arc regimes need conservative discontinuity
handling. Admit finite-volume or entropy-stable/high-resolution DG paths with
appropriate Riemann/limiter/WENO-like reconstruction when the selected PDE
requires them.

These methods should share geometry, fields, reaction sources, diagnostics and
certification contracts with the FE/FEEC path. Do not create a second plasma
state architecture merely because the spatial discretization differs.

Block/patch AMR and hp refinement are first-class when thin sheaths, streamers,
current sheets, shocks or reaction/radiation fronts dominate. Adaptivity must
carry conservative transfer and a measurable/bounded defect.

## Rarefied and strongly correlated paths

Neutral rarefaction may use deterministic Boltzmann or DSMC-like particle
transport. Particle collision methods may use MCC/binary-collision schemes with
explicit statistical status.

For strongly coupled classical reference problems, small direct N-body/MD
calculations are admissible as an independent representation/oracle. Their
pairwise execution should use tiled/reduced interactions where possible, but
the physical \(N\)-body model is distinct from mean-field FK plus a collision
operator.


## Prefer structured tensor geometry

"General 3D geometry" does not mean tetrahedralize everything.

The preferred path is

\[
\boxed{\text{mapped tensor patches}+\text{interfaces/escape hatches}.}
\]

This preserves regular data access, tensor contractions and high-order
accuracy. Use unstructured/cut elements only where the topology or wall
geometry makes them worthwhile.

For conventional fusion devices, a poloidal tensor patch times a Fourier or
periodic coordinate can be extremely efficient. Figure-8/G-frame and spherical
space-plasma charts use the same map abstraction.

## Compatible fields

For electromagnetics,

\[
E_h\in H({\rm curl}),\qquad
B_h\in H({\rm div})
\]

is preferred so the de Rham identities are represented structurally.

Scalar densities/temperatures and kinetic distribution coefficients may use
different spaces. Do not force all state variables into the same finite element
family for convenience.

## Tensor-product spline production path

The existing assembled FortFEM spline path is a reference. The production hot
path should add:

1. span-local basis evaluation returning only nonzero basis functions;
2. no allocation inside basis/quadrature application;
3. precomputed 1D basis/derivative tables;
4. sum factorization;
5. partial assembly of compact geometry/Hodge factors;
6. matrix-free operator application;
7. native many-RHS batching.

In 3D, never construct a full basis-by-quadrature tensor when three 1D
contractions suffice.

The same operator should retain an assembled small-problem reference for
validation/coarse/direct use.

## Many-RHS kinetic spatial actions

KIN6D frequently applies the same spatial operator to many velocity/harmonic
channels:

\[
Y=A_xX,
\qquad
X\in\mathbb R^{N_x\times N_{\rm channel}}.
\]

This is a first-class API, not a loop around a vector operator.

Geometry and basis data should be reused across large RHS batches, converting a
nominally bandwidth-bound FEM action into a denser tensor workload.

CPU and GPU schedules may choose different tile shapes while representing the
same action.

## Phase-space representation

A strong continuum FK reference is

\[
f_s(x,p,\Omega_p)
=
\sum_{\ell,m}
f_{\ell m}^s(x,p)Y_{\ell m}(\Omega_p).
\]

Each species has independent radial momentum extent/order and angular cutoff.
A local magnetic frame or restricted \(m\) is a later optimization requiring
its own geometry/derivative/certificate treatment.

Fourier modes are used where a true periodic coordinate exists. Derivatives in
that direction are diagonal; nonlinear products use measured direct-convolution
versus FFT/pseudospectral crossovers.

This representation is not mandatory for every regime. Beam/ring/loss-cone or
particle methods may use other velocity carriers.

## Sparse angular and collision structure

Angular coupling tensors are stored as sparse interaction graphs when sparse.
Do not store dense \(L^4\)-like Gaunt tensors.

At low harmonic count, direct sparse contractions may win. At larger count,
angular transforms/pseudospectral products may win. Benchmark the crossover.

For Rosenbluth/Landau radial solves, factors depending only on species/grid and
\(\ell\) should be constructed once and reused for all spatial/mode RHS. Use
batched banded/specialized solves when the discretization permits them.

## Particle/PIF/PIC numerics

The particle representation follows [representations.md](representations.md).

The generic hot patterns are:

- particle push;
- particle-field projection/gather;
- local PIC deposition;
- particle-mode PIF reduction;
- sorting/binning/migration;
- wall intersection;
- stochastic event processing.

PIF and other particle-in-basis couplings should use tiled reductions rather
than an explicit particle-by-mode matrix. Direct and NUFFT schedules are
benchmarked.

PIC deposition/current updates must respect the compatible field discretization
and charge continuity of the selected model.

Particles and continuum must be able to coexist without duplicating common
field, geometry, diagnostic or source machinery.

## KeOps/FlashAttention-style indexed execution

A logical operator

\[
y_i=\sum_jF(x_i,y_j,p)
\]

is represented as an indexed reduction. The schedule tiles \(i,j\), stages
reused data in cache/shared memory/registers, evaluates \(F\), and immediately
reduces. The full \(i\times j\) object is never written to global memory.

The same execution pattern applies to:

- PIF;
- boundary/pairwise kernels;
- some collision/reaction integrals;
- quadrature residuals;
- sparse graph interactions.

Do not replace a superior \(O(N)\) or separable physical algorithm by a
quadratic map/reduce merely because the scheduler can run it efficiently.

## Data layout

One canonical copy of the large kinetic state is preferred.

Use tiled AoSoA-like layouts so both spatial and velocity kernels can process
contiguous small blocks. When another orientation is essential, transpose a
bounded tile into scratch instead of maintaining a second full-state copy.

Local metadata may use 32-bit indices where safe; global identity/counts use
wider types when needed.

Resolution adaptation should initially use a small set of patch/block classes
rather than fully ragged per-cell dimensions. This preserves SIMD/GPU batching.

## Time integration and asymptotic preservation

The time consumer may choose explicit SSP/RK, implicit BDF/DIRK,
IMEX, exponential/integrating-factor, splitting, multirate/subcycling or
geometric/variational methods according to stiffness and structure.

When a parameter limit removes expensive fast scales, prefer an
asymptotic-preserving or uniformly accurate formulation over brute-force
resolution if it can be verified against both parent and limit models.
Relevant limits include quasineutral/Debye, collisional fluid, drift/gyro,
stiff-reaction and radiation-diffusion limits.

The time integrator does not define the physical model; a reduced target
equation and a reduced preconditioner for a fuller target remain distinct.

## Matrix-free nonlinear solves

The steady solver interface is

\[
R(U),\qquad J(U)v,\qquad J(U)^\dagger v,\qquad P^{-1}v.
\]

Do not assemble a global kinetic Jacobian.

The preconditioner may use cheaper physics without changing the target
equations: DK/GK/MHD blocks, mode-diagonal approximations, local collision
solves or coarse field solves are legitimate accelerators.

Krylov memory is a first-class budget. A state too large to hold 50 copies
requires small-restart or low-storage algorithms and stronger preconditioning.

## Static condensation/HDG

Where local interior unknowns can be eliminated algebraically, use static
condensation/Schur complements to reduce global communication and state.

Local elimination is exact with respect to the discrete model and therefore
certification-friendly, provided roundoff/solve error is bounded.

## Distributed-memory decomposition

The default kinetic decomposition is configuration space, with each owner
holding the complete local velocity/harmonic state. This keeps Coulomb
collisions local.

Change this only when memory or scaling data demonstrate the need.

Start halo exchange, compute local collision/reaction work while communication
is in flight, then complete spatial transport. Count global reductions
explicitly; they can dominate Krylov scaling.

## GPU residency

Once accelerated, the production state remains on device across residual/JVP
applications. CPU orchestration must not cause host/device ping-pong.

Use large batched kernels, not one launch per element/mode/species. Kernel
fusion is tuned against register pressure and occupancy; maximum fusion and
maximum occupancy are not objectives by themselves.

CUDA graphs, persistent kernels and asynchronous staged loads are later
schedule optimizations triggered by profiles.

## Mixed precision

Candidate solves may use FP32/TF32 or approximate local factors in
preconditioners if the strict FP64/rigorous residual and final certificate
remain independent.

The default safety pattern is

\[
\text{FP64 target state/residual}
+
\text{lower-precision preconditioner}
+
\text{rigorous final check}.
\]

Mixed precision does not weaken the mathematical target; it only changes how
the candidate is obtained.

## Adaptivity

Numerical and physical adaptivity are separate.

Numerical controls include

\[
h,\quad p,\quad N_p,\quad L,\quad M,\quad N_n,\quad \text{quadrature}.
\]

Physics controls include FK/GK/DK/gyrofluid/MHD, collision fidelity, neutral
model and field reduction.

The certifier/adjoint should eventually identify the dominant contribution to
the requested observable rather than blindly refining every axis.

## Reference implementations

For a hot mathematical operator, retain when appropriate:

\[
\boxed{
\text{simple/assembled reference}
\leftrightarrow
\text{optimized CPU}
\leftrightarrow
\text{optimized CUDA}.
}
\]

They must agree under independent tests. The slow path is valuable for
certification, debugging and cross-checking and should not be deleted simply
because a faster backend exists.

Chris&AI
