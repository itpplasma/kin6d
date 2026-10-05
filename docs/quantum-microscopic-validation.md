# Microscopic quantum experiments and validation

Status: implementation plan, 2026-10-05. Audience: lecture students and
implementers who need to see microscopic quantum physics and measure the
approximations used to reach larger systems. Planned gates are not support
claims. The existing pair quench and ideal Fermi EOS are retained references;
the localized collision and degenerate interacting plasma below remain open.

## Physical parent and comparisons

Start with the antisymmetric finite-N Schroedinger equation or density
operator for spin-1/2 electrons with unscreened Coulomb interactions. Use
isolated encounters for scattering and neutralized periodic 3D cells for
electron gases. Declare mass, charge, hbar, initial preparation, spin sector,
geometry, background and boundary conditions in every case. Nuclear motion
may initially be classical or prescribed, with that approximation explicit.
Nonrelativistic electrostatics is the first controlled domain; record v/c,
x_F and the omitted electromagnetic/relativistic physics. QED is outside the
native scope. No screened potential is inserted to demonstrate emergent
screening.

A microscopic quantum reference evolves the interacting state. Hartree,
HF, wave-packet, Wigner and fluid closures remain separately identified
approximations. A classical finite-N ensemble is a different physical model;
its difference from quantum mechanics is the phenomenon being illustrated,
not the quantum solver's numerical error. BBGKY evolution and stellar
structure remain deferred. Observable reduced density matrices and connected
correlations do not require implementing a BBGKY solver.

Compare quantum and classical preparations at the same physical parameters
and imposed disturbance. For unexchanged Gaussian preparations, a positive
initial Wigner law supplies a useful classical phase-space ensemble. For
antisymmetrized or thermal degenerate states, record the preparation
difference: do not force classical Maxwell occupations to equal Fermi
occupations or invent a positive law for a negative Wigner function.

## QG0: upstream infrastructure and baseline

Before the new encounter solver, freeze small correctness and realistic
performance cases for complex Hermitian actions, evolution, residual
evaluation, probability-bin observables, requested derivatives and all bands.
Retain the current 125/343-mode pair reference and its numerical scope.
The reviewed 90.89-times solve-only quantum cost is an optimization baseline,
not an acceptable workhorse performance claim.

Measure assembled and structured/matrix-free schedules, batching, immutable
Hamiltonian reuse and memory traffic. Choose the next implementation from
these measurements; do not add all backends speculatively. Compare cold and
amortized results at equal observable accuracy, including certificate setup.

Generic gaps belong to their owners: FortNum for verified quadrature,
evolution/residual and spectral evaluation primitives; FortSym for shared
complex primal/JVP/VJP and enclosure definitions; FortFEM for generic spatial
bases/actions and singular integration; FortSparse for needed sparse actions;
FortPlot for interval rendering. Extend an owner only when this case needs
the feature. Preserve project presets and repair the documented Fo resident
CMake input-capture gap rather than dropping the feedback gate.

Gate: independent behavioral reproducer, owner-side focused tests, verified
revision pinned and consumer rechecked for each required upstream increment.
Retain solve-only, error-estimation and memory baselines before adding physics.

## QG1: localized two-electron encounter

Prepare two separated incoming Gaussian electron packets with declared
centroids, momentum covariances, widths and impact parameter. Run opposite-spin
and same-spin preparations; state the actual spin density/state and its
spatial symmetry. Singlet/triplet references separate exchange effects from
interaction correlations. Start with repulsive electron-electron scattering;
binding to an ion is a separate future encounter.

Exploit exact center-of-mass/relative-coordinate separation for equal-mass
Coulomb dynamics when the initial preparation factorizes in those coordinates.
Matched equal-width incoming Gaussians provide such a family. Relative
motion is still a 3D quantum problem. Retain every angular sector needed for
nonzero impact parameter and exchange; do not confine particles to a line.
The existing fixed-total-momentum quench sector is not a localized packet.
General nonfactorizing preparations need a later 6D or controlled low-rank
representation, with its separate budget.

Use the bare Coulomb operator with a declared singularity treatment and
isolated-domain approximation. Bound or report domain, initial-tail and
boundary errors. Coulomb has long-range scattering phases: finite-time
outgoing bins and true asymptotic cross sections are different observables.
An absorber requires an explicit lost-flux/boundary ledger.

Classical comparisons are (a) point trajectories from the packet's admitted
initial phase-space ensemble and (b) an optional declared classical extended
charge model. A classical packet is initially an ensemble of independent
two-particle realizations, not a swarm interacting as extra electrons.
Centroid-only point trajectories are an additional labeled illustration.
Charge-smearing changes the interaction law; its difference is model error.

Independent oracles: free Gaussian propagation; classical Rutherford
scattering; quantum Coulomb scattering with declared spin, angular and
long-range conventions; independently assembled small Hamiltonians. The
distinguishable plane-wave Coulomb angular intensity can agree with
Rutherford, so do not design a test that demands a difference. Packet
diffraction, exchange and conditional correlations supply other comparisons.

Gate: reproducible incoming/encounter/outgoing movies and bounded integrated
probabilities, independent derivative oracles, conservation checks, initial
and discretization budgets, resolution audit and full cost ledger. Missing
continuum/boundary bounds remain explicit; finite-basis success alone does not
close the continuum gate.

## QG2: coherent packet deformation and splitting

The microscopic reference retains the complex wavefunction through
deformation, diffraction and separated outgoing regions. Repulsive electrons
do not acquire a bound capture channel. A single Gaussian is an approximation
whose failure is measured against that reference.

First compare one Gaussian, then a coherent sum of Gaussians or a richer
adaptive spatial basis on the same two-electron case. A sum retains complex
coefficients and all interference cross terms. Splitting a packet is a basis
enrichment, not the creation of additional electrons. Fit/evolve the new
representation with a measured state-transfer defect and preserve normalization.
Basis removal and merging require an error budget too.

Multiple Gaussians per one-electron orbital do not supply arbitrary
interelectron entanglement. Add a joint relative-coordinate expansion or a
controlled sum of antisymmetrized configurations when the correlated reference
requires it. Monitor overlap conditioning and rank growth with explicit memory
and refusal limits. Orthogonalization alone does not remove physical branches.

Gate: enriched representation reproduces integrated outgoing probabilities
and interference where present; report shape, exchange and correlation errors
separately. An incoherent mixture is permitted only with an explicit
decoherence/coarse-graining model and information-loss status.

## QG3: finite-N degeneracy and exchange

Build noninteracting finite-N fermion references first: occupied orbitals,
antisymmetric states and finite-temperature density operators. Validate
canonical partition sums independently. Grand-canonical Fermi--Dirac formulas
and thermodynamic-limit EOS are separate references, not exact fixed-N
occupations for every small cell. Declare spin populations and shell effects.

Show how the same-spin exchange hole and Fermi pressure occur without Coulomb
interaction. Validate one-body and pair observables from the full Slater or
thermal state against independent ideal-gas identities. Compare a classical
gas at the same density and temperature, retaining its different preparation.

Gate: correct exchange symmetry, particle count and thermal normalization;
independently checked pressure/energy and pair probabilities; temperature,
orbital and finite-size scans. These are microscopic statistical references,
not a claim that noninteracting electrons constitute a self-consistent plasma.

## QG4: interacting microscopic degenerate plasma

Use a neutralized 3D electron cell with bare Coulomb interactions. Start at
N=2; attempt N=4 and then richer cells only within measured Hilbert-space
budgets. Record orbital count, sector dimension, basis error, temperature and
coupling; N alone is not a capacity specification. Compute ground/thermal
states without an imposed Hartree factorization, then apply a weak test charge
or density disturbance and evolve the interacting state.

Measure pair probabilities, connected correlations, energy, pressure, static
response and density-mode evolution. Screening emerges from the interacting
system. A tiny box demonstrates microscopic response; thermodynamic screening,
Friedel tails and collective dispersion require independent finite-size and
wavelength limits. A discrete finite-N spectrum has recurrences and must not
be advertised as irreversible collision damping.

Independent checks combine noninteracting limits, separate tiny-basis
diagonalization/evolution and suitable equilibrium QMC/PIMC references.
Declare statistical, fixed-node and sign-problem limitations. Ground-energy
upper bounds do not imply pressure or lower-energy bounds; thermal spectral
and orbital tails need their own argument. Equilibrium comparisons do not
validate real-time dynamics.

Gate: at least one genuinely interacting degenerate microscopic comparison
with the classical finite-N cell, full error/status ledger, conditional-density
and mode-response movies. Finite-size/model differences stay separate from
numerical error. BBGKY RHS work remains deferred.

## QG5: validated scalable electron packets and orbitals

Compare Hartree, antisymmetric HF and restricted Gaussian variants against
QG1--QG4, changing one approximation at a time. Hartree/HF define reduced
state manifolds; single-determinant exchange is not general correlation.
Scalable product packets with effective Pauli forces have an additional
declared exchange approximation. Width confinement or fitted correlation
terms are physical model changes requiring independent validation.

Accept a scalable model only for named observables and a documented regime
scan. Prefer full-parent residual or controlled asymptotic estimates; use
sparse references when rigorous model bounds remain unavailable. Multiple
packets per orbital and multiple configurations enrich different missing
physics and have different costs.

After correctness, benchmark approximately 100, 1000 and 10000 electron
packets with accelerated electrostatics. These are targets, not promised
capacities; deep degeneracy may require orbital/Fermi methods instead.
Explore 100000--1000000 only after measured accuracy, force-tolerance,
memory and total-cost gates. CPU batching/vectorization precedes GPU/MPI
unless evidence selects another schedule. No million-electron full-state claim.

## QG6: local quantum encounters in a classical plasma

Use QG1--QG2 to validate local diffraction/exchange scattering and, where
appropriate, reusable transition kernels. The classical bulk retains real
physical particles. Indicators involving relative wavelength, impact scale,
occupation and encounter overlap choose quantum refinement; they do not
themselves bound error. A universal distance-only switch is insufficient.

Specify initial-state lifting, external-field treatment, replacement of the
classical near interaction, and outgoing reconstruction. Budget matching,
truncation, sampling and lost coherence; avoid double counting near forces
and keep momentum/energy balances. Overlapping encounters or persistent
coherence trigger richer treatment or explicit refusal. Occupation-dependent
Pauli blocking requires a population description even when encounters are local.

Gate: hybrid versus full small-system scattering, conservation and transfer
oracles, plus a cost comparison with classical and all-quantum runs at matched
observable accuracy. Bulk degeneracy never reduces to pair-only promotion.

## QG7: collective quantum reductions

Only after the corresponding microscopic references, compare finite-T Fermi
screening, Thomas--Fermi, Lindhard/RPA, HF response, Wigner and derived fluid
closures. Ideal-gas Lindhard response, interacting RPA and a full interacting
cell discharge different claims. Retain Debye/classical limits, Friedel
structure, finite-k oscillation/dispersion and Pauli-blocked collision
comparisons where their regimes apply. Small perturbations and long-wavelength
closures have explicit applicability and remainder status.

Check Wigner's classical limit without imposing positivity; a reduced Wigner
equation needs a declared correlation closure. Fermi pressure alone does not
validate dynamical damping or dispersion. DFT/TDDFT and QMC are optional pinned
providers when they supply a required independent reference, not prerequisites
for QG1. Transfer/lifting differences are reported before comparing models.

## QG8: optional white-dwarf pressure

Retain the existing numerical nonrelativistic ideal Fermi EOS. Extend the
momentum integrals to relativistic energies only as a small independent
increment; benchmark classical, low-T nonrelativistic and ultrarelativistic
limits with tail and root/quadrature bounds. Show pressure versus density and
temperature and the associated occupations.

Interaction/ion corrections require separate microscopic or literature
validation. A pressure calculation is not a correlated stellar simulation.
QG8 can follow QG3 alongside other work if inexpensive; it must not delay the
localized encounter. Hydrostatic stellar structure and neutron-star nuclear
matter are outside this implementation sequence.

## Error estimates and promotion

Every requested observable, including correlations and requested gradients,
has a scale, an absolute floor and its own status. Default target: about 1%
of that scale; bound-evaluation roundoff receives a small fraction of the
radius. Do not demand a 1% full-state error to deliver a 1% observable, or
infer the latter from a small diagnostic residual alone.

Maintain separate budgets for numerical evolution; spatial/orbital and
domain truncation; Hamiltonian/force evaluation; representation transfer;
physical model reduction; input uncertainty; and sampling. Triangle bounds
may combine applicable deterministic contributions. Statistical intervals
remain separately labeled, with their sampling law and dependence assumptions.
Classical ensemble realizations, not individual correlated particles or movie
frames, determine the sampling units.

For a normalized approximate wavefunction in the domain of the full declared
self-adjoint H, use the candidate residual r=i*hbar*d_t(psi_tilde)-H*psi_tilde.
Unitary stability gives

\[
e(t)\le e(0)+\hbar^{-1}\int_0^t\|r(s)\|\,ds,
\qquad
|\Delta\langle A\rangle|\le2\|A\|e(t)
\]

for bounded A and normalized states. Enclose residual integration and
reconstruction, including singular Coulomb and omitted basis components.
Finite-basis residuals only establish finite-basis error. Unbounded energy or
pointwise density needs additional regularity/moment estimates; integrated
probability bins provide bounded observables. Seek tighter goal-specific
estimates when the norm bound is too pessimistic.

Apply this route to packet/HF model defects in the full parent when affordable;
do not assume its cost stays linear or that it stays sharp for long times.
An auxiliary DFT residual certifies its own equation, not the interacting
many-electron parent. Thermal-state errors use spectral/partition bounds or
appropriate density-operator estimates, not a pure-state formula by analogy.

Label each estimate as a rigorous bound, a reference comparison, a convergence
indicator, or a statistical interval. Small-system calibration cannot prove a
large-system error bound; two resolutions and agreement between closures do
not prove a remainder. Publish missing terms and failed certificates rather
than displaying a fabricated complete band. Production promotion requires
the declared numerical certificate and model applicability gate to close.

## Figures and animations

Each artifact has a one-sentence teaching message, reproducible data, physical
units and normalization. Prefer the following small, synchronized views:

| Case | Main animation | Quantitative companion |
|---|---|---|
| QG1 encounter | Quantum one-body probability projections beside classical ensemble projections; common time and axes | Relative separation/angle-bin probabilities and their bands |
| QG2 enrichment | Full reference, one Gaussian and coherent enriched density | Outgoing-bin errors and declared model-error estimates versus basis size |
| QG3 ideal degeneracy | Occupation distribution/Fermi sphere and exchange-hole profile | Pressure and pair-bin probabilities versus classical reference |
| QG4 interacting cell | Density disturbance and conditional pair density with stated conditioning | Density-mode history, connected pair probabilities and finite-size scan |
| QG5 scaling | Small validated reference versus scalable model | Error and time/memory versus electron count |
| QG6 hybrid encounter | Classical approach, admitted quantum region and outgoing distribution | Transfer/conservation defects, scattering probabilities and extra cost |
| QG7 screening/oscillation | Induced density or potential and coherent mode response | Debye/Fermi/Lindhard curves or dispersion with applicability stated |
| QG8 pressure | Occupations as density/temperature vary, if useful | Numerical pressure with classical and degenerate limiting curves |

Use spatial projections integrating the omitted coordinate; a slice needs its
own explicit normalization. Quantum density is probability, not individual
electron trajectories. Classical plots state ensemble size and each displayed
realization's role. No unexplained matched-pair markers or inference of
correlation from one-body density alone. Conditional densities specify their
conditioning bin and normalization; a small denominator needs its own bound.

Use equal scales and a shared clock for comparisons. Keep the main collision
movie visually simple; technical diagnostics belong in companion figures.
Show full numerical intervals and separately styled model/statistical bands.
Probability heatmaps have companion bin/projection bounds or uncertainty maps;
colors alone do not constitute a bound. Smooth connected envelopes require a
valid interpolation allowance, not subpixel stairs. Inspect actual first,
encounter and outgoing frames, grayscale and labels before sharing.

Retain GIF, MP4, snapshots, CSV values/endpoints and manifests outside Git.
Lecture captions name the experiment and observable rather than internal
solver machinery. Keep numerical/refinement plots available for validation.

## Execution and completion

Immediate sequence: QG0, QG1, QG2, QG3, QG4. QG5 uses both encounter and
degenerate references; QG6 can branch from validated QG1/QG2; QG7 uses the
matching microscopic and ideal-gas references. QG8 is the inexpensive side
branch described above. This gallery does not reorder the direct-steady
scientific spine. Implement one complete runnable increment at a time.

Use Fortran and retained FortSym generators for owned numerical leaves,
primal/JVP/VJP and enclosures. Parameter sensitivities cover impact parameter,
packet width, coupling, density or temperature where smooth. Spin choices,
basis changes and promotion events are discrete; declare event derivatives or
their unavailable status and budget recomputation/checkpoints explicitly.

At each gate update the benchmark atlas with implemented cases, independent
oracles, source/dependency revisions and actual numerical status. Freeze
review inputs by base commit and patch digest. Retain native Release and
checked Debug CTest evidence, meaningful end-to-end postprocessing, and
independent review for substantive mathematical certificates.

Report solve-only time, error-estimation time, total time, added seconds,
total/solve multiplier, peak memory and repeated-trial variability. Separate
requested derivatives, preparation/reuse, sampling and rendering. Compare
models and external codes at matched observable accuracy. No overhead target
or published particle count becomes a measured KIN6D performance claim.

## Primary references

- [Svensson et al.](https://arxiv.org/abs/2211.08560): restricted Gaussian
  trajectories, approximate Pauli interactions and near-linear scaling.
- [Angermeier and White](https://arxiv.org/abs/2104.03441): distinct packet-shape,
  exchange and correlation approximation errors.
- [Chabrier and Potekhin](https://arxiv.org/abs/physics/9807042): electron-ion
  thermodynamics from classical to degenerate regimes.
- [White](https://arxiv.org/abs/2410.23599): dynamical response limitations of
  an orbital-free density approximation.

Chris&AI
