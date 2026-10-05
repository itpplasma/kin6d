# kin6d design

Status: canonical architecture, 2026-10-05.

## Mission

KIN6D is not one 6D solver and not a collection of compatibility modes. It is a
common plasma-physics engine in which the usual hierarchy of plasma models and
numerical representations can coexist, be compared one approximation at a
time, differentiated, optimized and, where the mathematics permits, certified.

The highest-fidelity charged-particle reference is relativistic 3X3P
Vlasov--Maxwell with multispecies nonlinear Coulomb collisions. Practical
calculations deliberately use cheaper models and representations when their
error relative to the scientific observable is acceptable.

The primary scientific target is **direct kinetic steady state**:
source/collision/transport/field-balanced non-Maxwellian states should be
computed from a stationary residual rather than by resolving irrelevant fast
waves for long physical times. Transient and turbulence consumers use the same
operators and are developed early because they expose structure and performance
defects quickly.

## Canonical equation interface

A selected model defines a residual and, when evolutionary,

\[
M(U)\dot U + R(U;\theta)=0.
\]

Stationary, forced-response, eigenvalue and trajectory problems are different
consumers of shared model operators. No subsystem is allowed to exist only as a
time-step routine if its mathematical residual can be stated directly.

Where established, closed reversible dynamics is written with a Poisson
operator and energy, irreversible dynamics with a metric/metriplectic
operator, and external interaction through ports:

\[
\dot U =
J(U)\frac{\delta H}{\delta U}
+
M_d(U)\frac{\delta S}{\delta U}
+
P_{\rm open}(U),
\qquad
J^\dagger=-J,\quad M_d^\dagger=M_d\succeq0.
\]

Claims of Hamiltonian, Poisson, symplectic, metriplectic, entropy or free-energy
structure require the corresponding discrete identities and degeneracies;
the notation is not decorative. Open-system sources, radiation, reactions,
walls and reservoirs carry an explicit exchange ledger instead of being forced
into an ordinary least-action principle.

## Orthogonal configuration axes

KIN6D configuration is intentionally factored.

### Physical model per species

\[
{\rm fixed/adiabatic/fluid/gyrofluid/gyromoment/PK/BA/DK/GK/FK}.
\]

A separate bulk model may aggregate explicitly owned thermal populations into
ideal/resistive/Hall/XMHD or multi-fluid state. Species and bulk ownership must
prevent double counting.

### Distribution form

\[
{\rm full}\text{-}f,\qquad
f=F_0+\delta f,\qquad
{\rm control\ variate/multilevel}.
\]

The reference \(F_0\) may be analytic, a previous steady solution or a cheaper
fluid/DK/GK state. The split is independent of whether the correction is stored
on a grid or carried by particles.

### Numerical carrier

\[
{\rm Eulerian},\quad
{\rm semi\!-\!Lagrangian/characteristic},\quad
{\rm particles},\quad
{\rm finite\ blobs},\quad
{\rm hybrid/remapped}.
\]

Particles may be deterministic quadrature points or stochastic samples.

### Particle-field coupling

\[
{\rm PIC},\quad
{\rm PIF},\quad
{\rm particle\!-\!FE},\quad
{\rm general\ particle\!-\!basis}.
\]

PIF is the Fourier specialization of a more general particle-in-basis
projection. PIC is a local-basis/grid specialization. Their physical equation
must match the continuum model selected.

### Sampling and stochasticity

\[
{\rm deterministic},\quad
{\rm QMC/randomized\ QMC},\quad
{\rm Monte\ Carlo}.
\]

Stochastic binary collisions, neutral histories, wall/reaction events and rare
event estimators are permitted. Statistical error is a separate result
component.

### Fields

\[
{\rm prescribed}
\rightarrow
{\rm Poisson/electrostatic}
\rightarrow
{\rm quasineutral\ EM}
\rightarrow
{\rm Darwin}
\rightarrow
{\rm Maxwell}.
\]

The selected species reductions must be compatible with the selected field
action. Background equilibrium, imposed drive and induced response remain
separate.

### Coupling roles

A population independently declares whether it evolves, contributes to fields
and feels induced response. Test populations, fixed reservoirs and active
species therefore do not require separate physics implementations.

## Model hierarchy

The intended reusable hierarchy is

\[
{\rm MHD}
\rightarrow
{\rm Hall/XMHD,\ multi-fluid}
\rightarrow
{\rm gyrofluid/Landau\ fluid}
\rightarrow
{\rm Hermite\!-\!Laguerre\ gyromoments}
\rightarrow
{\rm DK/GK}
\rightarrow
{\rm FK}.
\]

This is not a universal monotone accuracy ordering. Each reduction names its
ordering, closure, field action, phase-space measure and valid observables.

Gyrofluid is not a side code. Low-order gyrofluid and Landau-fluid models are
truncations/closures of a systematic gyromoment hierarchy wherever that
derivation is available. Increasing the moment set gives a controlled route
toward the selected gyrokinetic model, with moment-tail/closure errors kept
explicit.

MHD is a bulk aggregate, not a per-species moment label. Kinetic-MHD and
spatial hybrids require a derived pressure/current/exchange coupling and
matching conditions.

## Kinetic parent and collisions

The relativistic FK reference uses physical momentum \(p\),

\[
\partial_t f_s
+v_s\cdot\nabla_xf_s
+q_s(E+v_s\times B)\cdot\nabla_pf_s
=
\sum_{s'}C_{ss'}[f_s,f_{s'}]+S_s.
\]

A central bilinear collision primitive \(Q_{ab}(g,h)\) exposes test-particle,
field-particle, linearized and nonlinear pieces consistently. The production
ladder includes none, Krook/BGK, Lorentz, LBD/Dougherty, model FP, Landau and
Beliaev--Budker, with screened/matched/correlation models as later research.

For linearized collisional transport, declare the natural weighted pairing,
collision invariants, self-adjoint or Onsager--Casimir transformation,
dissipation form and constrained coercivity. NEO-2-style independent reciprocal
coefficients are permanent verification cases.

## Geometry

Physical coordinates and stored physical vectors are Cartesian. Computational
elements use maps

\[
x=X(\xi^1,\xi^2,\xi^3).
\]

Cylindrical, toroidal, spherical, G-frame/curve-following and arbitrary mapped
patches are specializations. The architecture may exploit field-period or
axisymmetry but may not require it.

FCI is a field-line map/operator service for anisotropic parallel dynamics, not
the fundamental geometry. Full kinetic transport remains valid when topology
changes or field lines hit a wall.

The Figure-8 stellarator is a deliberate geometry acceptance test because it
guards against accidental dependence on a conventional toroidal angle.

## Spatial and phase-space discretization

FortFEM owns generic spatial spaces/maps/operators. Compatible FEEC spaces are
preferred for Maxwell/MHD fields. Tensor-product B-splines/NURBS, Fourier,
DG/HDG and multipatch mortars form the high-performance regular path;
cut/unstructured methods handle regions that need them.

For continuum FK a strong reference representation is

\[
f_s(x,p,\Omega_p)=\sum_{\ell m}f^s_{\ell m}(x,p)Y_{\ell m}(\Omega_p),
\]

with a species-specific radial momentum basis/grid and Fourier representation
in genuinely periodic configuration directions. This is a backend, not a
universal requirement.

The production execution rule is to exploit tensor/separable structure before
forming operators. Matrix-free sum factorization, partial assembly and
many-RHS batching are preferred. No global kinetic Jacobian or pairwise kernel
matrix is part of the normal state.

## Particle and blob discretizations

A particle representation is

\[
f_h(z)=\sum_p w_p\delta(z-Z_p),
\]

and a blob representation replaces the delta by a finite shape
\(S_p(z-Z_p)\), possibly with adaptive width/covariance.

The particle-field projection is derived from the same selected field
variational principle where possible. Charge-conserving geometric PIC,
variational PIF and particle-in-FE therefore share field spaces and structure
with continuum clients rather than becoming a separate code.

Deterministic/QMC particle initialization is supported alongside random
sampling. Remapping may reconstruct a smooth distribution and redistribute
particles. Hybrid continuum/particle populations are allowed, especially for
energetic tails, when transfer and conservation defects are controlled.

## Differentiability

Each continuous project-owned operator has one retained symbolic mathematical
definition whenever practical. FortSym generates primal, JVP and VJP; HVP is
added for a concrete consumer. FortAD and independent directional/analytic
checks validate derivatives but do not define the production path.

Global solves are differentiated implicitly. For

\[
R(U,\theta)=0,
\quad
R_U\,\delta U=-R_\theta\,\delta\theta,
\]

and for scalar \(Q\),

\[
R_U^\dagger\lambda=Q_U^\dagger,\qquad
\frac{dQ}{d\theta}=Q_\theta-\lambda^\dagger R_\theta.
\]

No full sensitivity matrix is required. Physical pairings, not blindly
Euclidean transposes, define adjoints.

Stochastic algorithms expose derivative semantics separately. Fixed-randomness
pathwise derivatives, likelihood-ratio estimators or common-random-number
methods are optional tools and carry statistical uncertainty.

## Certification and uncertainty

The candidate solver is untrusted. A Fortran certifier using verified
interval/ball arithmetic checks the mathematical premises needed for the
claimed result. Lean is an optional theorem/checker layer after useful
numerical certificates exist.

Result components remain separate:

- numerical discretization/algebraic error;
- model/reduction error;
- physical-input/data uncertainty;
- stochastic sampling/confidence error.

For SPD/coercive operators, residual/energy bounds and equilibrated fluxes are
preferred because they provide tight inexpensive certificates. General
non-normal problems require constrained inverse/inf-sup/resolvent bounds.
Adjoints convert residual defects into goal-oriented observable estimates and
eventually certified model-error corrections.

## Performance architecture

The mathematical definition is separated from indexed operator structure and
hardware schedule. Indexed reductions such as

\[
y_i=\sum_jF(x_i,y_j)
\]

are logical objects, not matrices. KeOps/FlashAttention-style tiled
map/reduce kernels evaluate and reduce tiles while data are resident in
registers/cache/shared memory.

Production design targets:

- no per-call allocation in hot kernels;
- one canonical large-state layout plus bounded scratch;
- tile/block resolution classes rather than maximally ragged state;
- factor/geometry reuse;
- batched radial/sparse/angular solves;
- direct-vs-FFT and sparse-vs-transform crossovers selected by measurement;
- configuration-space distributed ownership by default;
- overlapped communication and local collision work;
- bounded-Krylov memory;
- reduced-model preconditioners and control variates;
- GPU-resident state when accelerators are used.

Fast candidate arithmetic may use mixed precision or relaxed execution when the
strict certifier remains independent.

## Physical scope

Long-range capabilities include equilibria, MHD/XMHD, gyrofluid/gyromoments,
DK/GK/FK, PIC/PIF/blobs, edge/SOL/divertor, deterministic and MC neutrals,
atomic/molecular reactions, sheaths/walls, radiation, RF/full-wave heating,
NBI, fusion reactions/products, vacuum/coils/circuits, Solar-System plasmas,
synthetic diagnostics, inference, design and control.

The architecture must remain one engine, but features are admitted only when a
real use case, derivation, independent oracle, differentiability statement,
certificate target and performance target exist.

## Immediate implementation boundary

The implemented 1X1V slice remains the authoritative first rung. Its compact
Fourier representation must not dictate the future 6D layout. The next work is
self-consistent 1X1V/representation comparisons and the 0X3V collision
laboratory, while an early HM/gyrofluid branch exercises nonlinear structure.

Do not split KIN6D into continuum/PIC/PIF/gyrofluid repositories. If a generic
particle engine or indexed-kernel scheduler acquires multiple independent
consumers, extract the stable generic layer upstream only after the abstraction
is demonstrated.

Chris&AI
