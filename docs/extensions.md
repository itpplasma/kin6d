# Actuators, fields, diagnostics and further extensions

Status: reserved scientific directions, not implemented capabilities. This
map completes the current requirements without fixing their implementation
order. Each feature needs a concrete scientific question, derived compatible
coupling, independent analytic/manufactured and literature/cross-code evidence,
and a measurable certification/performance goal before production admission.

## Common mathematical engine

The state contains only the active blocks needed by the selected case:
species distributions/moments, an optional bulk aggregate, field unknowns,
and eventually neutral, surface, wave, radiation or circuit states. Algebraic
closures and constraints are distinguished from evolution blocks. Consumers
solve steady equations, forced response, eigenproblems or trajectories using
the same derived operators; the discrete evolution may be a mass-matrix DAE.

The mathematical decomposition reserves transport, collisions, reactions,
waves/heating, fields, sources and boundary/surface actions. This is not a
requirement to implement a plugin framework or all interfaces immediately.
Shared libraries can serve thin steady, transient, response, eigenproblem and
certifier drivers as real consumers arise. One engine need not mean one
executable with hundreds of switches. Unsupported compositions fail during
configuration validation, with no silent lower-fidelity fallback.

## Field-model axis

| Choice | Required definition |
|---|---|
| Prescribed fields | External field/source compatibility and interpolation bounds |
| Electrostatic Poisson | Charge equation, potential/gauge, boundary and loop-voltage treatment |
| Quasineutral electromagnetic | Specific ordering, constraint/Ohm law and compatible species response |
| Darwin | Explicit longitudinal/transverse split, gauge and boundary conditions |
| Full Maxwell | Displacement current, induction, Gauss constraints and compatible boundaries |

These are distinct approximations, not a universally nested ladder. In the
usual Darwin formulation the transverse displacement current is omitted
while the longitudinal contribution maintains charge consistency. Simply
deleting all displacement current is not a complete Darwin model. A named
quasineutral model is likewise more than substituting `rho=0` into Maxwell.
Reduced-species polarization/magnetization, charge continuity and energy must
match the selected fields. Reduced fields need their own model-error theorem
before their numerical certificate can describe a full-Maxwell solution.

AP formulations and physics preconditioning retain the requested equations
while addressing stiffness; selecting Darwin or quasineutrality changes the
model. Implicit stability does not prove accurate unresolved fast waves.
Certificate tails/temporal defects determine what a coarse representation can
claim. Multirate, exponential and AP algorithms are added when physical-scale
benchmarks demonstrate the need, not as unused day-one backends.

## Waves and heating

Reserve ECRH/ECCD, ICRH, LHCD and generic antenna sources through a wave/heating
axis: prescribed quasilinear diffusion, ray plus quasilinear response, and
full wave coupled to kinetic response. Space wave scattering uses the same
family. Each specifies frequency conventions, spectrum, resonance broadening,
coordinates/measure, polarization, boundary conditions and energy/momentum
exchange. Diffusion positivity does not establish correct deposited power
or current drive. Rays require a justified short-wavelength regime;
mode conversion or resonances can require a different model.

Resolved Maxwell fields and an imposed quasilinear operator need an explicit
scale/energy partition to avoid counting the same wave interaction twice.
A frozen dielectric computed from a background differs from self-consistent
nonlinear wave/particle feedback. Wave amplitudes, dielectric singularities,
resonance integration and unresolved spectrum enter certification premises.
KAMEL is a planned heritage/reference source only for documented overlapping
wave/linear-response equations; its full capabilities are not presumed here.

NBI reuses neutral trajectories, attenuation/ionization/CX and a specified
fast-ion birth distribution, then the selected FK/GK/DK transport and coupling
role. Gas puffing and pellets are external mass/energy/phase-space sources;
pellet ablation/deposition is a separate closure. Track injected particles,
charge, momentum, energy and any beam/pellet reservoir. A source is not
validated solely because its integrated heating power matches an input.

## Fusion reactions and radiation

Extend the reaction graph to D–T, both D–D branches, and D–helium-3 as needed.
For distribution-aware rates use the appropriate relative-flux/cross-section
measure. Identical reactants need the pair-counting factor `1/2` in a number
reaction rate. Include fuel depletion, branching ratios, reaction mass defects
and product distributions with matching source/sink bookkeeping.

D–T alpha birth has approximately 3.5 MeV in the conventional low-energy
center-of-mass reference. It is not always an isotropic monoenergetic source
in the laboratory frame: reactant energies and center-of-mass motion broaden
and shift products. Neutrons can initially escape to a recorded energy,
momentum, nuclei and diagnostic channel. Full neutron transport and nuclear
engineering are outside initial scope. Data provenance and uncertainties are
separate from integration/discretization error.

Radiation choices extend optically thin atomic losses to bremsstrahlung,
cyclotron/synchrotron emission and relativistic radiation reaction. Account
for emitted energy and force consistently; radiative force can alter phase-
space volume, as discussed in [space scope](space-plasmas.md). Photon transport
is deferred until opacity or a requested diagnostic requires it. Prescribed
radiation losses and resolved photon dynamics must not duplicate energy sinks.

## Conductors, circuits and equilibrium

Add vacuum domains, coils, finite-conductivity structures and external circuit
states when free-boundary evolution, resistive-wall response, vertical
stability or feedback requires them. A circuit model specifies resistance,
inductance/mutual coupling, imposed voltages, wall currents and plasma fluxes
in one consistent power balance. It must not duplicate induced fields already
represented in the field subsystem. Time-varying coil fields include their
compatible electric drive, not only a changing magnetic snapshot.

Extend equilibria beyond static isotropic ideal force balance to flowing,
anisotropic, energetic-particle and fully kinetic states as needed. Their
force/moment equations and boundary data differ. An invariant distribution
`F(H,P_phi,...)` is stationary only for actual invariants of the chosen dynamics;
magnetic moment is generally approximate for full orbits. Collision/source
balance and field self-consistency remain additional obligations. Imported
MHD equilibrium is a candidate/initialization, not a unique kinetic lift.

Positive masses and signed/zero charges are species parameters, without a
hard-coded electron/ion taxonomy. Pair-plasma benchmarks can therefore be
admitted without a special solver branch. Collision/reaction data must still
support the pair. Special-relativistic fluid/MHD reductions are deferred
models with their own equations; arbitrary masses do not implement them.
GR, QED, pair creation, detailed solid-state evolution and comprehensive
radiation hydrodynamics remain outside initial scope.

## Diagnostics, inverse problems and control

A synthetic diagnostic maps a reconstructed state to measurement space,
including instrument geometry, response, resolution and calibration. Targets
include interferometry, Thomson scattering, ECE, CXRS/FIDA, spectroscopy,
bolometry, magnetic/Mirnov probes, neutron cameras, fast-ion loss detectors,
reflectometry, Langmuir probes and IR wall loads. Some need waves, opacity,
atomic data or a sheath model; they are not all simple local moments.

Propagate numerical enclosures through the diagnostic, including bounded
quadrature and instrument mapping. Experimental noise/calibration and
physical-model uncertainty have separate labels. Compare with observations
in measurement space without interpreting numerical error as an experimental
confidence interval.

Inverse reconstruction and assimilation can later use forward diagnostic
operators, parameter constraints, observation covariance/likelihood and
regularization. A local optimizer or adjoint supplies a candidate inference;
it does not prove identifiability, uniqueness or a global optimum. Verified
parameter regions and numerical/model-error-aware inference are separate
research gates. Optimization and control reuse these operators only after
the forward model works. Coil, actuator, stellarator and divertor optimization
and feedback design are downstream applications; no hardware action is
implicitly authorized by a simulation configuration.

Use fortad derivatives where supported, with independently checked actions
and the fortsym-generated numerical-leaf policy. Non-smooth limiter, event,
model-switching and adaptive-mesh decisions need explicit derivative treatment;
AD of a floating-point execution is not automatically a derivative theorem
for the continuous solution or its certificate.

## Adaptive fidelity and difficult observables

Separate numerical refinement (`h`, spatial order, radial momentum basis/
extent, Fourier and spherical harmonics) from physical reduction per species
or spatial region. Data structures grow the smallest needed seam; no immutable
universal layout forbids later changes. Refine from rigorous error
contributions, retaining old uncompressed/fixed-resolution oracles.

Model adaptation requires a validated ordering/remainder and conservative
interfaces, transfers and population ownership. Small heuristic ordering
indicators alone are not model-error bounds. Non-additive interactions between
closures must be controlled before selecting the cheapest model automatically.

Rare-event/loss observables include wall loads, escape times, separatrix
crossings and loss fractions from deterministic distributions or specified
random laws. A norm-ball does not directly bound a discontinuous hitting
indicator. Grazing trajectories, event-location uncertainty, boundary trace
regularity and probability/flux measures need dedicated enclosures. Statistical
rare-event estimates remain statistical. Prescribed fluctuation samples define
conditional deterministic cases; claims over a stochastic field ensemble need
a specified law and their own probabilistic/robust theorem.

Topology-changing fields leave physical geometry valid. FCI maps are rebuilt
or invalidated when their field changes; never assume a smooth fixed map
through reconnection or magnetic nulls. Their approximation contributes to
the numerical certificate for any model that uses them.

## Units, provenance and the final error claim

Physical SI definitions are the conceptual contract. A run declares its
normalization, dimensional scales, charge/sign conventions, measures and
transformation rules. Check dimensional consistency and independent physical
invariants; test conversions with known physical cases, not only round-trip
identity. Version the model/normalization definition alongside the state.

Every result records code/generator/dependency/proof revisions; model and
coupling configuration; fields and geometry; all initial/boundary/source data;
collision, reaction and wall datasets; normalization; resolution; hardware/
compiler arithmetic assumptions; reconstruction; certificate and diagnostic
identifiers. Content hashes make external assets identifiable without
committing them. Failed or unavailable proof premises remain visible.

A combined model-plus-numerical observable enclosure may add error intervals
only when a proved comparison chain connects each to the same target quantity
under compatible assumptions. Paired simulations yield a certified difference
of two named models, not a remainder to physical reality. Correlated input
uncertainties and closure interactions cannot be independently added by habit.
No available theorem means no asserted bound for that contribution.

The performance metric is total time and peak memory to obtain a requested
certified observable tolerance, with solve/certification/postprocessing times
reported separately. Report failure/time limits and certificate strength;
uncertified fast runs are not a faster certified result.

Chris&AI
