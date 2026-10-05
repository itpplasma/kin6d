# kin6d design

Status: canonical architecture, 2026-10-05.

## Non-negotiable scientific principles

These are **project invariants, not scope items**. YAGNI applies to features,
models, algorithms and abstractions; it does not permit a vertical slice to
drop these properties merely to ship sooner. A slice may implement the
smallest useful instance of a principle, but must preserve the contract needed
to strengthen it later.

1. **One mathematical truth / DRY.** Each physical operator, convention,
   configuration meaning and conservation law has one authoritative production
   definition. Primal, derivative, enclosure and optimized specializations come
   from that definition wherever practical.
2. **Differentiability by construction.** Every smooth continuous production
   operator used by a differentiable consumer exposes consistent primal/JVP/VJP
   semantics from the same mathematical definition. Nonlinear/equilibrium/
   linear/eigen solves are differentiated at equation level, not through
   incidental iterations. Nonsmooth events and unavailable derivatives are
   explicit, never silently approximated.
3. **Error control and certification by construction.** A numerical answer is
   incomplete without an error/status contract. The fast solver is untrusted;
   strict residuals, independent checks and interval/ball or other rigorous
   bounds are developed wherever the mathematics permits. Numerical error,
   model-reduction error, input/reconstruction uncertainty and stochastic
   sampling error remain distinct.
4. **Structure preservation.** Conservation, symmetries, Hamiltonian/Poisson or
   metriplectic structure, entropy production, positivity and
   Onsager--Casimir relations are preserved or their defect is explicitly
   measured. Numerical convenience does not silently change the physics.
5. **Independent verification.** Important results retain analytic,
   manufactured, structural or cross-code oracles independent of the production
   path. Shared-code agreement alone is not evidence.
6. **Reproducibility and provenance.** Effective inputs, model assumptions,
   source/dependency/compiler revisions, transformations and result status are
   recoverable from the scientific artifact.
7. **Performance at controlled accuracy.** Performance is judged at equal
   physical observable accuracy/certificate width, with time and memory both
   first-class. Optimization may change implementation, never the target
   mathematics.
8. **Clean separation of policy and details.** Physics/model semantics are
   independent of replaceable solver, file-format, hardware, scheduler and UI
   details. New scope composes with the same contracts instead of spawning a
   second numerical universe.

A proposed feature that conflicts with these principles is rejected or remains
an explicitly non-production experiment, even if it would broaden the nominal
scope.

## Mission

KIN6D is not one 6D solver and not a collection of compatibility modes. It is a
**general plasma multiphysics engine** in which plasma, neutral,
radiation, surface and transport models plus their numerical representations
can coexist, be compared one approximation at a time, differentiated,
optimized and, where the mathematics permits, certified.

The intended application envelope covers fusion, Solar-System/natural,
industrial/low-temperature, beam/accelerator, atmospheric, dusty/non-neutral
and classical high-energy-density plasmas. The deepest classical
charged-particle reference is the multispecies relativistic Maxwell--Lorentz
many-body system. Relativistic 3X3P Vlasov--Maxwell/FK is the highest continuum
mean-field reference. Practical calculations deliberately use cheaper models,
closures and representations when their error relative to the scientific
observable is acceptable.

The full design space is recorded in the capability and numerical-method
catalogues. They prevent architectural dead ends but do not determine
implementation order. Quantum plasma fundamentals are included where they
clarify or supply plasma closures; the fundamental ownership boundary is QED.

The primary scientific target is **direct kinetic steady state**:
source/collision/transport/field-balanced non-Maxwellian states should be
computed from a stationary residual rather than by resolving irrelevant fast
waves for long physical times. Transient and turbulence consumers use the same
operators and are developed early because they expose structure and performance
defects quickly.

## Clean dependency rule

KIN6D freezes scientific policies and contracts earlier than implementation
details. Clean Architecture here is a dependency rule, not a prescribed folder
tree.

Physics/model definitions are the inner policy. Application consumers such as
steady solves, stability and transport compose that policy. External solvers,
file formats, accelerator runtimes, schedulers, CLI/web frontends and other
framework details live outside those layers and depend inward through narrow
contracts/adapters.

Do not introduce an abstraction merely to make a detail theoretically
replaceable. Use the simplest direct implementation for the current vertical
slice and extract a seam when a real second implementation/consumer, testing
boundary or costly external dependency demonstrates it. This keeps details
swappable without freezing hypothetical interfaces prematurely.

The repository-wide KISS/YAGNI/DRY/SRP/DIP/ISP and related rules are specified
in [software-engineering.md](docs/software-engineering.md).

## Microscopic classical parent

The continuum kinetic hierarchy is itself a reduction. Beneath it, KIN6D
recognizes the exact finite-\(N\) classical microscopic state: individual
relativistic charged particles coupled to the full Maxwell field.

For species \(s\), particle \(i\),

\[
\dot{\mathbf x}_{si}
=
\mathbf p_{si}/(\gamma_{si}m_s),
\qquad
\dot{\mathbf p}_{si}
=
q_s(\mathbf E+\mathbf v_{si}\times\mathbf B),
\]

with microscopic charge/current sources in Maxwell's equations.

The electrodynamic fidelity ladder is independent of development priority:

\[
\text{nonrelativistic Coulomb}
\rightarrow
\text{relativistic Coulomb}
\rightarrow
\text{Darwin}
\rightarrow
\text{full Maxwell--Lorentz}
\rightarrow
\text{direct retarded Lienard--Wiechert reference}
\rightarrow
\text{declared radiation-reaction extension}.
\]

Full Maxwell--Lorentz is the target endpoint. Coulomb and Darwin are reduced
benchmarks. Direct Lienard--Wiechert evaluation is a derived retarded,
field-eliminated reference under matching radiation conditions, not a more
fundamental model than the Maxwell initial-value problem.

Point-charge self-fields are singular, so mutual interaction, independent field
data and self-force/radiation-reaction modeling are kept distinct. Ordinary
fusion collision/correlation validation normally does not require a
self-radiation model.

This parent is a declared classical reference with an explicit self-field
prescription, not an already implemented singular point-charge solver. Early
teaching cases use nonrelativistic instantaneous Coulomb in a stated regime;
retardation and self-force are not prerequisites for collective examples.

The corresponding statistical reduction is

\[
F_N
\rightarrow
(f_1,g_2,g_3,\ldots)
\rightarrow
\left\{
\begin{array}{l}
\text{Vlasov/FK mean-field limit},\\
\text{Balescu--Lenard/Landau/BB effective weak correlations},\\
\text{explicit cumulant/mean-force closures}.
\end{array}
\right.
\]

See [microscopic relativistic many-body reference](docs/microscopic-many-body.md).
BBGKY supplies the unclosed marginal hierarchy; the displayed reductions are
different limits/closures rather than a sequence of collision operators.
Collective polarization and waves already occur in mean field. Pairwise
microscopic forces can generate higher connected correlations without a
fundamental three-body force term. Track mean response, finite-particle noise
and collisionality separately in the declared Debye-volume scan.

## Quantum fundamentals boundary

The classical hierarchy is not asserted beyond its domain. KIN6D also admits a
quantum-plasma branch that exposes the transition from classical ideal plasma
to degenerate/correlated quantum plasma.

The main regime coordinates are the degeneracy parameter

\[
\Theta=k_BT/E_F,
\]

the Coulomb-coupling parameter \(\Gamma\), and the relativistic Fermi parameter

\[
x_F=p_F/(m_ec).
\]

The reduction hierarchy may use finite-temperature Fermi gases,
Thomas--Fermi/Lindhard/RPA response, quantum hydrodynamic closures, Wigner
kinetics, Hartree/Hartree--Fock/DFT-like orbital references, quantum
BBGKY/cumulants and quantum collision/dielectric closures.

The numerical representation is independent of the word "quantum": exact
many-body wavefunctions/density matrices are tiny-\(N\) references; orbital
methods live on ordinary 3D meshes; Wigner models reuse phase-space
Eulerian/semi-Lagrangian/DG/spectral machinery; quantum-fluid/EOS closures use
ordinary FEM/FV/DG on macroscopic domains.

The native fundamentals branch stops at **QED**. General quantized photon
dynamics, pair creation/annihilation, strong-field QED cascades and radiative
QED corrections remain external. QED-derived rates or closures may enter
through provider interfaces.

See [quantum plasma fundamentals](docs/quantum-plasma.md).

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

The microscopic reference is not restricted to strongly coupled plasmas: small
weakly coupled boxes also validate Debye shielding, dielectric response and
collision closures. Direct (N^2), FMM/Ewald/P3M and field-mediated Maxwell
paths accelerate the same physical particles; macro-particles/blobs change the
finite-(N) microscopic model unless explicitly declared. Fully quantum
many-body physics remains outside the classical core.

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

Certification is available from the first supported case and covers every
reported scientific output. Each result states its pointwise, cellwise,
observable, norm, weak or statistical claim and any unavailable/failed bound.
Norm bounds are not pointwise plot bands; particle measures need declared
observable or reconstructed-field claims.

The practical default targets a 1% observable radius relative to a declared
physical scale and absolute floor. Initially, bound evaluation receives 1%
of that radius after stability or adjoint amplification. Workhorse arithmetic
uses fast centers and proved aggregate roundoff/remainder bounds; strict
outward arithmetic supplies audits and fallback. Both modes enclose the true
quantity under stated hypotheses. Tight last-bit endpoints are unnecessary.

General certificates check a continuous residual, reconstruction, boundaries
and truncation against a proved stability/inverse bound. Known solutions are
independent benchmark oracles. The first general case is implemented for
independently supplied Fourier/P1 periodic damped transport, followed by planned
forced FP with continuous flux/stability control. The manufactured Poisson/FP
laboratory is reference evidence,
not a general nonlinear certificate.

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
certifier remains independent and its arithmetic theorem accounts for its own
execution. The performance goal is the fastest measured total result time at
matched model and requested accuracy, including every reported band and
requested derivative. Upstream verified-arithmetic and generated-kernel
benchmarks precede optimization and broader physics; external comparisons are
required before claiming that this goal is achieved.

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
