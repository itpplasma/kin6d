# Long-range physics extensions

Status: scope map. These capabilities are admitted only when a concrete
scientific case reaches them; this file is not an implementation checklist to
complete before useful results.

## Field reductions

The field model is independent of the species model:

\[
{\rm prescribed}
\rightarrow
{\rm Poisson/electrostatic}
\rightarrow
{\rm quasineutral\ electromagnetic}
\rightarrow
{\rm Darwin}
\rightarrow
{\rm full\ Maxwell}.
\]

Each reduction states exactly which waves, displacement current, polarization
or quasineutral constraints are removed. Reduced fields may also be used as
preconditioners without changing a fuller target residual.

## RF and wave-particle physics

Progress from a prescribed quasilinear diffusion tensor

\[
\partial_tf
=
\partial_{J_i}\!\left(D_{ij}^{QL}\partial_{J_j}f\right)
\]

to ray/quasilinear coupling and eventually full-wave kinetic response.

Target applications include:

- ECRH/ECCD;
- ICRH;
- LHCD;
- generic wave-particle diffusion;
- space-plasma resonant scattering.

KAMEL/KiLCA/KIM provide valuable internal reference cases. Full-wave
electromagnetics should reuse FortFEM \(H({\rm curl})\) spaces and open-boundary
operators.

## NBI, pellets and external sources

Sources are phase-space operators \(S_s(x,p,t)\), not special solver modes.

NBI eventually includes neutral trajectories, attenuation, ionization/CX,
birth distributions and active/passive fast-ion evolution.

Pellets/gas puffing create localized plasma/neutral sources with their own
mass/energy/momentum ledgers.

## Fusion reactions

At minimum support

\[
D+T\rightarrow\alpha+n
\]

plus D-D and D-\(^3\)He.

For arbitrary distributions,

\[
R_{12}
=
\int f_1f_2\,\sigma(g)g\,d^3v_1d^3v_2.
\]

Fusion-product birth distributions become ordinary kinetic sources. Neutrons
are initially an escaping energy/diagnostic channel; complete neutronics is
outside scope.

## Relativistic electrons and radiation

The FK parent already uses physical momentum. Add when demanded:

- relativistic Beliaev--Budker collisions;
- runaway acceleration;
- synchrotron/cyclotron radiation;
- bremsstrahlung;
- radiation reaction;
- distribution-aware emission.

DREAM/CODE/CQL3D-like regimes are comparison targets, not named modes.

## Atomic, molecular and internal-state network

Represent reactions through a sparse stoichiometric graph

\[
\sum_a\nu^-_{ar}X_a\rightarrow\sum_a\nu^+_{ar}X_a,
\qquad
S_a=\sum_r\nu_{ar}R_r.
\]

Processes include ionization, recombination, excitation, charge exchange,
dissociation, MAR/MAD, photoionization and chemistry.

Fidelity spans:

1. tabulated local-field/Maxwellian effective rates;
2. electron-energy or reduced Boltzmann/EEDF-derived rates;
3. collisional-radiative/internal-state population equations;
4. distribution-aware rate integrals;
5. differential kinetic reaction kernels.

Large weakly ionized/reacting cases may contain hundreds of species/states and
thousands of reactions. Use generated sparse residual/JVP/VJP evaluation rather
than a domain-specific chemistry solver. See
[reacting and partially ionized plasmas](reacting-partially-ionized.md).

ADAS/OPEN-ADAS, AMJUEL, HYDHEL and similar sources remain external data with
version/hash provenance.

## Neutrals

Neutral transport supports:

\[
{\rm fixed}
\rightarrow
{\rm fluid}
\rightarrow
{\rm deterministic\ kinetic}
\leftrightarrow
{\rm Monte\ Carlo\ kinetic}.
\]

EIRENE is an early backend/oracle. A native deterministic solver is valuable
for differentiability and rigorous intervals; MC histories remain valuable in
long-mean-free-path and reaction-rich regimes.

Hybrid deterministic/MC neutral regions are admissible only with explicit
transfer and uncertainty accounting.

## Sheath and material wall

Keep

\[
\text{bulk plasma}\rightarrow\text{sheath}\rightarrow\text{material}
\]

as separate layers.

Sheath fidelity spans algebraic/Bohm transmission closures through kinetic
closures to locally resolved FK Debye-scale regions.

A wall patch is a boundary operator with state,

\[
\mathcal B_{\rm wall}:
\{f_s^{in}\},W\mapsto\{f_s^{out}\},\dot W.
\]

Processes include absorption, reflection, implantation, recycling,
neutralization, molecular formation/release, sputtering, secondary electrons,
desorption and pumping.

Detailed solid-state evolution remains external unless a specific problem
requires a bounded reduced material state.

## Radiation

Radiation is a native fidelity axis rather than a binary local-loss/full-transfer
choice:

\[
\text{optically thin}
\rightarrow
\text{opacity/escape closure}
\rightarrow
\text{grey diffusion}
\rightarrow
\text{multigroup/moment transport}
\rightarrow
\text{angular transport}
\rightarrow
\text{external high-fidelity non-LTE spectral provider}.
\]

This is sufficient to make reduced radiation hydrodynamics an ordinary KIN6D
configuration for dense/thermal plasma where justified. Detailed line-by-line
angle/frequency non-LTE transfer remains an external specialty provider until a
prioritized case justifies internalizing a subset.

Radiated energy/momentum and photo-reaction coupling appear explicitly in the
total ledger. See [radiation](radiation.md).

## Compressible, beam and high-energy-density classical extensions

The common fluid/Maxwell/reaction/radiation framework should admit
shock-capable compressible plasma, nonideal EOS/opacity closures and
laser/beam energy deposition when a prioritized case needs them. This covers
classical portions of arcs, propulsion, atmospheric-entry, Z-pinch,
magnetized-target/ICF-like and laser-produced plasmas without a separate HEDP
code architecture.

Relativistic Maxwell + FK/PIC already covers many beam/wakefield problems.
Moving-window, envelope or specialized source models are numerical/model
optimizations to add only when demonstrated useful.

Strongly coupled classical benchmarks may progress from screened/mean-force
closures to small direct N-body/MD reference calculations. Fully quantum dense
matter remains outside the classical core.

## Vacuum, conductors, coils and circuits

For MHD response/free-boundary/control problems, couple

\[
{\rm plasma}+{\rm vacuum}+{\rm conducting\ structures}+{\rm coils}+{\rm circuits}.
\]

This supports RMPs, RWMs, vertical stability, feedback and free-boundary
evolution.

External circuit states may obey ODEs such as

\[
L\dot I+RI+M\dot I_{\rm plasma}=V.
\]

Field coupling must avoid double counting induced electric/magnetic fields.

## Synthetic diagnostics

Every experimental diagnostic is a forward operator

\[
D_{\rm syn}[U].
\]

Targets include:

- interferometry;
- Thomson scattering;
- ECE;
- CXRS/FIDA;
- spectroscopy/bolometry;
- magnetic/Mirnov probes;
- neutron diagnostics;
- fast-ion loss detectors;
- reflectometry;
- Langmuir probes;
- IR heat flux.

Instrument response, calibration and experimental noise are distinct from
numerical error.

## Inverse problems and control

Once forward diagnostics and implicit adjoints exist, support:

- profile/source inference;
- equilibrium reconstruction;
- wall/recycling parameter inference;
- coil/stellarator/divertor optimization;
- actuator optimization;
- feedback-control design.

Optimization convergence does not prove identifiability or a global optimum.
Certified gradients and verified parameter regions are later, stronger claims.

## Solar-System plasma extensions

Reuse the same physics hierarchy with:

- gravity;
- rotating-frame forces;
- radiation pressure;
- open/inflow/outflow boundaries;
- CGL/six/ten-moment/Landau-fluid closures;
- pickup-ion distributions;
- charge exchange and ENA diagnostics;
- photoionization;
- mass loading;
- surface charging/photoelectrons.

Dust is architecturally possible but low priority.

## Adaptive physical fidelity

The long-range goal is to choose model and numerical fidelity by its effect on
the requested observable.

Possible adaptation axes include:

- MHD/gyrofluid/GK/FK by species or region;
- full-f versus delta-f/control variate;
- continuum versus particle/blob;
- collision operator fidelity;
- neutral fluid/deterministic/MC treatment;
- atomic-rate fidelity;
- field reduction.

A switch is guaranteed only when state transfer, interface and model remainder
bounds exist. Before that, adaptation is heuristic with explicit status.

## Rare-event methods

Wall losses, alpha confinement, escape times and tail observables deserve
first-class support.

Available methods may include:

- deterministic phase-space quadrature;
- discrepancy-controlled QMC;
- randomized QMC;
- Monte Carlo;
- importance sampling;
- splitting/rare-event methods.

The target certified alpha-loss observable is

\[
L(T)=\int 1_{\tau(z)\le T}f_0(z)\,dz,
\]

with deterministic or probabilistic bounds and explicit control of the
grazing/ambiguous set.

## Outside the initial scope

Do not silently absorb:

- general relativity/GRMHD;
- strong-field QED and pair creation/cascades;
- fully quantum dense-matter electronic structure;
- full reactor neutron transport/nuclear engineering;
- comprehensive bulk solid mechanics/material microstructure;
- comprehensive liquid/solid chemistry;
- the highest-fidelity line-by-line angle/frequency non-LTE radiation-transfer
  discipline.

Reduced radiation hydrodynamics through grey/multigroup/moment closures,
classical nonideal EOS/opacity models and conservative adjacent-media interfaces
remain inside the KIN6D envelope. External specialties may couple through narrow
reusable provider interfaces.

Chris&AI
