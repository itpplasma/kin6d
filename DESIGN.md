# kin6d design

Status: canonical architecture, 2026-10-05.

## Mission

KIN6D is not one 6D solver and not a collection of compatibility modes. It is a
**general classical plasma multiphysics engine** in which plasma, neutral,
radiation, surface and transport models plus their numerical representations
can coexist, be compared one approximation at a time, differentiated,
optimized and, where the mathematics permits, certified.

The intended application envelope covers fusion, Solar-System/natural,
industrial/low-temperature, beam/accelerator, atmospheric, dusty/non-neutral
and classical high-energy-density plasmas. The highest-fidelity
charged-particle reference is relativistic 3X3P Vlasov--Maxwell; practical
calculations deliberately use cheaper models, closures and representations
when their error relative to the scientific observable is acceptable.

The full design space is recorded in the capability and numerical-method
catalogues. They prevent architectural dead ends but do not determine
implementation order.

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

## Equilibria and stationary states

An equilibrium is a stationary solution of the selected equations,

\[
R(U;\theta)=0,
\]

not a special file type and not necessarily a unique mapping
\(\theta\mapsto U\). If two distinct stationary solutions exist, KIN6D stores
and treats them as two equilibria. Their initial states, continuation histories,
external inputs and solution identifiers are provenance; KIN6D does not need a
general probabilistic branch framework.

Axisymmetric MHD must acquire a native Grad--Shafranov path for fixed and free
boundary problems. EQDSK/gEQDSK input is an interchange representation: import,
normalize conventions, evaluate physical and equation residuals, optionally
re-solve with the native Grad--Shafranov model, and certify the named result.
A free-boundary re-solve also requires the external coil/current or equivalent
boundary data; an EQDSK file alone does not reconstruct all of that physics.

In three dimensions, global nested flux surfaces are an optional cheap
representation, not the equilibrium ontology. VMEC/DESC/GVEC-like states are
valuable candidates, initial guesses, reduced representations and independent
comparators. The general physical-space MHD formulation must remain meaningful
with islands, stochastic regions, separatrices and wall contact. This does not
remove the equations' own restrictions: for smooth static ideal MHD,
\(B\cdot\nabla p=0\), so an ergodic magnetic region cannot carry an
arbitrary smooth pressure gradient. HINT/PIES-like full-pressure relaxation and
SPEC/SPECTRE-like relaxed-region models are benchmark/model precedents, each
with its own physical assumptions.

Topology/history is specified by the physical model rather than by a universal
switch. Ideal evolution preserves its frozen-in constraints; MRxMHD-like
relaxation retains a different finite invariant set; resistive/XMHD evolution
may reconnect. A stationary result records which model, constraints, initial
state and continuation/evolution path selected it.

The same architecture extends beyond MHD. DK/GK/FK equilibria are stationary
solutions of their corresponding equations, and source/collision/wall-balanced
open systems are stationary states even when the word equilibrium is
inappropriate. Imported states may seed these solves but never bypass the
native residual/certificate contract.

See [equilibria and inference](docs/equilibria-inference.md).

## Orthogonal configuration axes

KIN6D configuration is intentionally factored.

### Physical model per species

\[
{\rm fixed/adiabatic/fluid/gyrofluid/gyromoment/PK/BA/DK/GK/FK}.
\]

A separate bulk model may aggregate explicitly owned populations into reacting
neutral/multispecies compressible flow, ideal/resistive/Hall/XMHD or
multi-fluid state. Species and bulk ownership must prevent double counting.
Reduced Boltzmann/EEDF models are kinetic reductions with their own angular/
energy representation and reaction-rate contract, not special industrial
solver modes.

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


### Interaction and correlation fidelity

The microscopic interaction model is another independent axis:

\[
\text{collisionless mean field}
\rightarrow
\text{model collisions}
\rightarrow
\text{Landau/Beliaev--Budker}
\rightarrow
\text{screened/mean-force nonideal closures}
\rightarrow
\text{direct classical correlated-particle reference}.
\]

Direct N-body/MD is reserved for small strongly coupled classical reference
problems; fully quantum many-body physics remains outside the classical core.

### Ionization, chemistry and neutral matter

The same species/reaction graph supports charge states, atoms, molecules,
negative ions, metastables and internal states. Heavy-particle fidelity ranges
from prescribed gas through reaction-diffusion and multispecies compressible
flow to rarefied kinetic/DSMC-like transport.

Electron fidelity for weakly ionized plasmas can progress from local-field/
electron-energy closures through multi-term Boltzmann/EEDF to FK electrons.
See [reacting and partially ionized plasmas](docs/reacting-partially-ionized.md).

### Radiation

Radiation is a state/closure axis:

\[
\text{optically thin}
\rightarrow
\text{opacity/escape}
\rightarrow
\text{grey diffusion}
\rightarrow
\text{multigroup/moment transport}
\rightarrow
\text{angular transport}
\rightarrow
\text{external high-fidelity non-LTE spectral provider}.
\]

Reduced radiation transport is native KIN6D physics; implementing the entire
line-by-line non-LTE transfer discipline is not a prerequisite for the common
plasma engine. See [radiation](docs/radiation.md).

### Coupling roles

A population independently declares whether it evolves, contributes to fields
and feels induced response. Test populations, fixed reservoirs and active
species therefore do not require separate physics implementations.

### Analysis consumer

The same physical model may be consumed as:

\[
\text{steady/equilibrium},\quad
\text{linear initial value},\quad
\text{eigenvalue/spectrum},\quad
\text{forced response/resolvent},\quad
\text{nonmodal growth},\quad
\text{continuation/bifurcation},\quad
\text{nonlinear stability},\quad
\text{transient/statistical evolution}.
\]

This is an orthogonal axis. Linear stability uses the same \(R_Uv\) and
\(R_U^\dagger v\) actions needed by Newton, adjoints and optimization.
Kink/ballooning/tearing/RWM/Alfvén and gyrokinetic microinstabilities therefore
do not require independent solver architectures.

### Macro transport and multiscale closure

A macro transport state may be 0D, 1D, 1.5D, 2D or 3D. Its conservative core is

\[
\partial_tu_a+\nabla\cdot\Gamma_a=S_a,
\]

while the closure provider may range from prescribed \(D,\chi,V\) through
neoclassical/Onsager and linear/quasilinear models to nonlinear flux-tube
fluxes, nonlocal orbit kernels, electromagnetic cross-spectra and direct
global gyrofluid/GK/FK fluxes.

Coupling strength ranges from offline tables and in-memory loose coupling to
iterated transport, separated-timescale micro/macro coupling, sparse global
microphysics sampling and fully coupled global evolution. Spatial locality and
timescale separation are independent assumptions.

Traditional equilibrium--transport splitting is an important cheap
configuration, not a legacy workaround. KIN6D should reproduce QL-Balance/
ASTRA-like integrated modelling early and then relax its approximations toward
2D/3D nonlocal and fully coupled evolution.

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

## Inference boundary

KIN6D owns forward physics, stationary/evolution solves, observables,
JVP/VJP/implicit-adjoint actions, numerical certificates and model-reduction
error estimates. TIAGO owns reconstruction/inference from experimental data:
likelihoods, diagnostic/calibration nuisance parameters, priors and posterior
uncertainty.

For a KIN6D state defined by \(R(U,\theta)=0\), TIAGO should be able to
consume matrix-free actions for the composed diagnostic map without knowing
Newton/Krylov details or the internal state layout. The VMEC++ reconstruction
already implemented in TIAGO is a retained provider/reference path; KIN6D is a
planned provider, not a reverse dependency.

If several KIN6D equilibria are compatible with the same measurements, TIAGO
may fit each separately and, when useful, report a finite mixture of local
Laplace approximations. The equilibria remain distinct physical solutions; do
not average topologically or physically distinct states into a fictitious
equilibrium merely to obtain one mean and covariance.

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

Long-range capabilities include equilibria, linear/nonlinear stability and
bifurcation, 0D--3D integrated/global transport, MHD/XMHD,
gyrofluid/gyromoments, DK/GK/FK, PIC/PIF/blobs, compressible/shock-capable
plasma, edge/SOL/divertor, deterministic/MC/rarefied neutrals, partial
ionization and large sparse chemistry/internal-state networks, drift-diffusion
and ionization fronts, conducting/dielectric/reactive surfaces, reduced native
radiation transport, RF/full-wave/wave-kinetic physics, energetic and fusion
products, vacuum/coils/circuits, Solar-System/natural, industrial,
beam/accelerator, atmospheric, dusty/non-neutral and classical HEDP
configurations, synthetic diagnostics, TIAGO inference integration, design and
control.

The architecture must remain one engine. Named application domains are thin
configurations over common primitives. Features are admitted only when a real
use case, derivation, independent oracle, differentiability statement,
certificate target and performance target exist.

See [the complete capability catalogue](docs/capability-catalogue.md) and
[numerical-method catalogue](docs/numerical-methods.md).

## Configuration, result and UX boundary

Fortran owns the typed configuration semantics and validation. The canonical
human-editable text syntax is Fortran namelist. A run preserves both the
literal input and a fully resolved namelist with defaults/automatic choices
materialized. Machine-readable schema output may drive documentation, AI
configuration proposals, Python convenience layers and a browser UI, but those
layers must not redefine model compatibility or numerical defaults.

Scientific result files use a NetCDF-4-compatible HDF5 subset through FortIO.
Every result embeds the original and resolved input plus source/dependency,
compiler/flags, hardware and external-input provenance. A KIN6D result from
which its effective input cannot be recovered is invalid. Large external inputs
may be referenced by normalized metadata and content hashes instead of copied
verbatim.

The command-line/Fortran path remains sufficient on a cluster. Optional Python
may provide notebook/server orchestration and TypeScript may provide a browser
front end, including scheduler submission and post-processing. Neither is a
scientific runtime dependency. Build/test support targets GNU, Intel LLVM,
NVIDIA and LLVM Flang families; compiler-specific optimization and strict
certifier flags are target-local, never hard-coded globally.

See [configuration, UX, I/O and reproducibility](docs/interface-reproducibility.md).

## Immediate implementation boundary and prioritization

The implemented 1X1V slice remains the authoritative first rung. Its compact
Fourier representation must not dictate the future 6D layout. The current
high-value path remains self-consistent 1X1V/representation comparisons, the
0X3V collision laboratory, early stability/gyrofluid work and cheap integrated
transport.

The broader catalogue is **not an implementation queue**. At each planning
cycle choose the vertical slice with the strongest combination of scientific
value, cross-domain reuse, dependency unlock, validation strength,
performance/certification leverage and acceptable cost/risk.

Do not split KIN6D into continuum/PIC/PIF/gyrofluid/industrial/space
repositories. If a generic particle engine, reaction runtime, radiation
transport primitive or indexed-kernel scheduler acquires multiple independent
consumers, extract the stable generic layer upstream only after the abstraction
is demonstrated.

Chris&AI
