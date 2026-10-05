# Radiation hierarchy

Status: canonical radiation scope. Radiation is a native physical axis through
useful reduced transport models; the highest-fidelity non-LTE spectral problem
may remain an external provider.

## Why radiation is an axis

The relevant discriminator is optical depth and spectral coupling, not the
application label. A fusion core, solar corona, arc, photosphere and
laser-produced plasma can require very different radiation closures.

Radiation exchanges energy, momentum and sometimes reaction rates with the
plasma and material surfaces. The selected closure therefore belongs in the
same residual/source ledger as collisions, reactions and fields.

## Fidelity hierarchy

### R0 — optically thin local emission/loss

Use emissivity and radiated-power models when photons escape without important
reabsorption. This covers much magnetic-fusion core/coronal work.

Outputs may include line/continuum emissivity for synthetic diagnostics in
addition to the energy sink.

### R1 — opacity and escape closures

Add frequency/group-dependent opacity, absorption and escape-factor style
closures when self-absorption is important but a full transport solve is not.

These are physical models with stated geometry assumptions, not numerical
corrections.

### R2 — grey radiation diffusion

For optically thick regimes,

\[
\partial_t E_r + \nabla\cdot F_r = S_r,
\qquad
F_r\approx -D_r\nabla E_r
\]

with a flux limiter/closure only when mathematically specified.

Couple matter/radiation energy conservatively.

### R3 — multigroup diffusion and moments

Resolve broad spectral groups and, where necessary, radiation moments such as

\[
E_{r,g},\quad F_{r,g},\quad P_{r,g}.
\]

Closures may include diffusion, \(P_N\)-type or M1-like moments with verified
limits and causality/positivity properties appropriate to the model.

### R4 — angular transport

Where \(\tau\sim1\) and anisotropy matters, admit discrete-angle/characteristic
or Monte-Carlo photon transport as another transport representation.

This may be native when a focused case justifies it, using the same mapped
geometry, boundary and stochastic-certification infrastructure.

### R5 — high-fidelity non-LTE spectral provider

Detailed line-by-line, angle/frequency-dependent transfer with scattering,
frequency redistribution and large non-LTE atomic populations is a major
discipline of its own.

KIN6D should expose a provider interface for such a solver and may internalize
subsets later, but it does not require this capability to claim general
classical plasma coverage.

## Atomic populations and opacity

Radiation and collisional-radiative populations are coupled but distinct.
Population kinetics determines emissivity/opacity; the radiation field in turn
changes photoexcitation, photoionization and energy exchange.

The coupling may be loose, iterated or monolithic through the same multiscale
contract used elsewhere.

All atomic/opacity data carry units, provenance, validity ranges and uncertainty
status.

## Applications without dedicated modes

The hierarchy can cover, at the appropriate rung:

- optically thin fusion/core and coronal losses;
- selected edge/divertor line trapping;
- arcs/thermal plasmas;
- atmospheric-entry radiating gas;
- dense laser/Z-pinch/HEDP calculations;
- stellar interior diffusion;
- photospheric/chromospheric problems up to the fidelity supported natively.

Application names do not change the radiation state contract.

## Numerics

Candidate numerical families include:

- finite-element/finite-volume diffusion;
- implicit stiff matter-radiation energy coupling;
- multigroup block solves;
- moment closures;
- characteristic/discrete-angle transport;
- deterministic or Monte-Carlo photon packets.

The method is chosen from the general numerical catalogue rather than coupled
to one radiation fidelity level by convention.

## Certification and uncertainty

Separate:

- numerical transport error;
- group/angular truncation;
- opacity/emissivity interpolation;
- closure/model error;
- atomic-data uncertainty;
- Monte-Carlo sampling error.

A grey or escape-factor model may be certified numerically without claiming it
is an exact approximation to a line-resolved non-LTE plasma.

Chris&AI
