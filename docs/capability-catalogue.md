# Plasma-physics capability catalogue

Status: canonical long-range scope map, 2026-10-05. This document describes
the **available design space KIN6D intends to cover**, not implementation order
and not a promise that every Cartesian product is meaningful.

## Product goal

KIN6D aims to be a one-stop shop for **classical plasma physics** across

- magnetic and inertial/magnetized fusion regimes;
- Solar, heliospheric, magnetospheric and planetary plasmas;
- industrial/low-temperature/reacting plasmas;
- beam, accelerator and laser plasmas;
- atmospheric and terrestrial discharges;
- dusty, non-neutral, ultracold and other laboratory plasmas;
- classical high-energy-density and strongly coupled regimes where an
  appropriate classical/nonideal model remains valid.

Breadth comes from a reusable set of species, field, collision, reaction,
radiation, surface, transport, stability and numerical primitives. Application
names select configurations and validation cases; they do not define new solver
architectures.

The catalogue is deliberately broader than the active roadmap. The roadmap
selects a small number of vertical slices by scientific value and leverage.

## Orthogonal physics axes

### Species/state fidelity

Charged and neutral populations may use, as appropriate,

\[
\text{fixed/algebraic}
\rightarrow
\text{drift-diffusion or moment fluid}
\rightarrow
\text{multi-fluid/anisotropic fluid}
\rightarrow
\text{gyrofluid/Landau-fluid}
\rightarrow
\text{gyromoments}
\rightarrow
\text{reduced Boltzmann/EEDF}
\rightarrow
\text{orbit/bounce averaged, PK, DK, GK}
\rightarrow
\text{FK}.
\]

These are not a single monotone accuracy chain. Their phase spaces, orderings
and closures differ.

A separate bulk state can represent neutral/reacting gas dynamics or
ideal/resistive/Hall/XMHD plasma. Kinetic minorities and reduced species may
couple to the bulk with explicit population ownership.

### Distribution representation

For any compatible kinetic state:

- full-f, delta-f or control-variate formulation;
- Eulerian/spectral/DG/discrete-velocity continuum;
- semi-Lagrangian/characteristic;
- particles, PIC, PIF, particle-in-FE;
- finite blobs/remapped particles;
- deterministic quadrature, QMC or Monte Carlo initialization;
- hybrid continuum/particle representations.

Classical direct N-body/molecular-dynamics representations are reserved for
strong-correlation reference problems where mean-field kinetic closure is not
adequate.

### Electromagnetic-field fidelity

\[
\text{prescribed}
\rightarrow
\text{electrostatic/Poisson}
\rightarrow
\text{quasineutral EM}
\rightarrow
\text{Darwin}
\rightarrow
\text{full Maxwell}.
\]

Vacuum, coils, conducting structures, dielectrics, circuits and external drives
use the same field infrastructure.

### Interaction/correlation fidelity

\[
\text{collisionless mean field}
\rightarrow
\text{model collisions}
\rightarrow
\text{Landau/Beliaev--Budker}
\rightarrow
\text{screened/Balescu--Lenard-like}
\rightarrow
\text{mean-force/nonideal closures}
\rightarrow
\text{direct classical correlated-particle reference}.
\]

A higher interaction model is admitted only with a regime statement and
verification target. Fully quantum many-body physics is outside the classical
core.

### Ionization, chemistry and internal states

The same reaction graph supports:

- charge-state evolution;
- ionization/recombination/attachment/detachment;
- excitation/de-excitation and metastables;
- molecular dissociation/association;
- charge exchange;
- vibrational/rotational/electronic internal states;
- collisional-radiative reductions;
- distribution-dependent and differential reaction kernels.

Sparse networks may contain hundreds of species/states and thousands of
reactions. Chemistry is data-driven; no domain-specific chemistry solver is
hard-coded.

### Matter/ionization regime

KIN6D should span

\[
\text{fully ionized}
\leftrightarrow
\text{partially ionized/reacting}
\leftrightarrow
\text{weakly ionized gas},
\]

with consistent charged-neutral momentum/energy exchange. Heavy-particle gas
fidelity ranges from prescribed background through diffusion/reaction and
multispecies compressible flow to rarefied kinetic/DSMC-like descriptions.

See [reacting and partially ionized plasmas](reacting-partially-ionized.md).

### Radiation fidelity

Radiation is an independent closure/state axis:

\[
\text{optically thin emissivity/loss}
\rightarrow
\text{line/continuum opacity and escape closures}
\rightarrow
\text{grey diffusion}
\rightarrow
\text{multigroup diffusion/moments}
\rightarrow
\text{angular transport/moment closure}
\rightarrow
\text{external high-fidelity non-LTE spectral provider}.
\]

This keeps radiation native where reduced transport is sufficient while
avoiding a requirement to implement the entire non-LTE line-transfer discipline
inside KIN6D.

See [radiation](radiation.md).

### Surfaces and adjacent media

One surface/interface contract covers:

- absorbing/reflecting/recycling/pumping walls;
- kinetic sheaths;
- conducting electrodes and circuits;
- dielectrics and evolving surface charge;
- secondary, thermionic, field and photo-emission;
- sputtering, implantation, desorption;
- reactive surface coverages, deposition and etching;
- finite conducting bodies/spacecraft/moons;
- coupling to an external liquid/solid/material backend.

Detailed bulk solid/liquid chemistry and mechanics stay outside the plasma core;
the interface exchanges conserved fluxes and state variables.

### External forces and frames

Reusable force terms include gravity, rotation, radiation pressure, drag and
other declared body forces. The conservative phase-space measure and energy
ledger must be derived for the selected model.

## Bulk and reduced plasma hierarchies

### Fluid/MHD family

\[
\text{drift-diffusion/reactive fluid}
\rightarrow
\text{compressible multi-fluid}
\rightarrow
\text{ideal/resistive MHD}
\rightarrow
\text{Hall/XMHD}
\rightarrow
\text{kinetic-MHD}.
\]

Compressible/shock-capable variants serve arcs, propulsion, atmospheric/HEDP
plasmas and space shocks as well as fusion transients.

### Magnetized kinetic family

\[
\text{gyrofluid/Landau-fluid}
\rightarrow
\text{gyromoments}
\rightarrow
\text{DK/GK}
\rightarrow
\text{FK}.
\]

### Low-temperature electron family

\[
\text{local-field/tabulated closure}
\rightarrow
\text{electron-energy model}
\rightarrow
\text{two-/multi-term Boltzmann/EEDF}
\rightarrow
\text{full kinetic electrons}.
\]

This supplies transport and reaction rates without requiring full electron FK
where a reduced EEDF is sufficient.

### Energetic-particle family

\[
\text{test orbit}
\rightarrow
\text{focused/orbit-averaged transport}
\rightarrow
\text{QL wave diffusion}
\rightarrow
\text{active DK/GK/FK}.
\]

The same hierarchy applies to fusion fast ions, radiation-belt/SEP populations,
pickup ions and beam/accelerator populations when assumptions match.

## Analysis and dynamics

Every compatible model may be consumed as:

- equilibrium or open-system stationary state;
- linear initial-value problem;
- eigenvalue/spectral stability problem;
- forced response/resolvent;
- transient/nonmodal amplification;
- continuation and bifurcation;
- nonlinear finite-amplitude stability/saturation;
- transient evolution;
- statistically stationary/turbulent state.

Reconnection, shocks, edge transitions, disruptions and other named phenomena
are combinations of these consumers with the appropriate physical model, not
new top-level architectures.

## Transport and multiscale coupling

Macro transport may be 0D, 1D, 1.5D, 2D or 3D. Closure fidelity spans

\[
\text{prescribed}
\rightarrow
\text{neoclassical/local}
\rightarrow
\text{linear/QL}
\rightarrow
\text{nonlinear local}
\rightarrow
\text{nonlocal orbit/global kernel}
\rightarrow
\text{global nonlinear microphysics}
\rightarrow
\text{fully coupled evolution}.
\]

Microphysics may supply local coefficients, fluxes, Onsager matrices, response
kernels, mode data, wave/turbulence cross-spectra or direct global flux fields.

Spatial locality and timescale separation are independent assumptions.

## Waves and statistical closures

KIN6D should support

- prescribed quasilinear diffusion;
- ray/eikonal transport;
- wave-kinetic/statistical-spectrum evolution;
- full-wave Maxwell/kinetic response;
- self-consistent nonlinear field evolution.

Wave spectra can be external reservoirs or evolved states. Cross-spectra and
phase information are first-class when transport depends on them.

## Natural and industrial coverage through shared primitives

The catalogue should be able to compose, without new solver families:

- fusion core, edge/divertor, energetic particles, stability, RF and reactor
  scenario transport;
- solar wind/corona, magnetospheres, ionospheres, shocks, reconnection and
  energetic particles;
- semiconductor/plasma-processing reactors;
- Hall/ion/MPD thrusters and plasma sources;
- arcs, torches and circuit-breaker plasmas;
- atmospheric streamers/lightning-like discharges;
- plasma medicine/catalysis/combustion assistance;
- laser/beam/wakefield and non-neutral plasmas;
- dusty/complex and ultracold classical plasmas;
- classical HEDP/Z-pinch/magnetized-target/ICF-like regimes to the extent that
  the chosen EOS, opacity and radiation closures remain valid.

These examples test coverage; they do not receive dedicated core modes.

## Deliberate outer boundary

KIN6D may couple to but does not initially own:

- fully quantum dense-matter electronic structure/Wigner/TDDFT many-body
  physics beyond controlled reduced closures;
- strong-field QED pair cascades;
- general-relativistic gravity/plasma;
- complete reactor neutronics and nuclear engineering;
- detailed bulk solid mechanics/material microstructure;
- comprehensive liquid/solid chemistry;
- highest-fidelity line-by-line, angle- and frequency-resolved non-LTE
  radiative transfer.

Narrow reusable interfaces may be added when a concrete plasma application
needs them.

## Catalogue versus priorities

This document is a **capability catalogue, not a backlog**. New work is chosen
by the prioritization rules in [ROADMAP.md](../ROADMAP.md). A capability should
be implemented only when a vertical slice demonstrates scientific value and
clarifies the abstraction. The catalogue prevents architectural dead ends; it
does not justify speculative framework construction.

Chris&AI
