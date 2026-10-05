# Reacting and partially ionized plasma models

Status: canonical cross-domain scope. These capabilities serve fusion edge,
industrial/low-temperature, atmospheric and planetary plasmas through one
species/reaction/transport framework.

## Principle

Weakly and partially ionized plasmas are not a separate KIN6D application.
They are configurations in which charged species, neutral gas, internal states,
reactions, fields and surfaces have comparable importance.

The reusable state may contain

\[
\{f_e,f_i,f_n\},
\quad
\{n_a,\mathbf u_a,T_a\},
\quad
\text{internal-state populations},
\quad
E,B,
\quad
\text{surface state}.
\]

Each population chooses the cheapest adequate closure.

## Electron hierarchy

Low-temperature electron kinetics needs an explicit reduced path:

\[
\text{local-field/tabulated mobility/rates}
\rightarrow
\text{electron-energy equation}
\rightarrow
\text{two-/multi-term Boltzmann EEDF}
\rightarrow
\text{full kinetic electrons}.
\]

Reaction and transport coefficients are functionals of the declared EEDF, not
automatically Maxwellian functions of one temperature.

For a reaction \(r\),

\[
K_r[f_e]
=
\int f_e(\mathbf v)\,\sigma_r(E)\,v\,d^3v
\]

or the corresponding reduced-energy/angular form. The same cross-section data
must be usable by the reduced Boltzmann and FK paths.

## Heavy-particle/neutral hierarchy

Heavy particles span

\[
\text{prescribed background}
\rightarrow
\text{reaction-diffusion}
\rightarrow
\text{multi-temperature/multispecies fluid}
\rightarrow
\text{compressible Navier--Stokes}
\rightarrow
\text{rarefied kinetic/DSMC-like}.
\]

The charged-neutral exchange operator accounts for mass, momentum, energy,
charge exchange and reactions. A fluid neutral model and a kinetic neutral model
are two representations of named physical closures, not separate application
codes.

## Internal states and chemistry

Sparse reaction networks may include atoms, ions, molecules, negative ions,
metastables and vibrational/rotational/electronic states.

Possible reductions include:

- explicit state-to-state kinetics;
- grouped levels;
- collisional-radiative quasi-steady closures;
- effective vibrational/rotational temperatures;
- LTE where justified.

Reaction records carry thresholds, internal energies, product distributions,
units, validity ranges and source provenance.

FortSym-generated sparse residual/JVP/VJP evaluation is the preferred owned
implementation for large networks; generic sparse linear algebra belongs
upstream.

## Drift-diffusion and ionization fronts

A cheap charged-fluid branch has the generic form

\[
\partial_t n_s+\nabla\cdot
\left(n_s\mathbf u_s
-\mathbf D_s\nabla n_s\right)
=S_s,
\]

with mobility/drift supplied by the selected field/transport closure and
Poisson/Maxwell coupling as appropriate.

This supports ionization fronts and streamer-like problems when combined with
reaction and photon/source terms. Front resolution requires AMR or an
equivalent adaptive/shock/front-capable spatial method; a uniform fine mesh is
not the architecture.

Photoionization may begin as a nonlocal reduced source and later use the common
radiation transport hierarchy.

## Partial-ionization electrodynamics

Reduced one-/multi-fluid models may include, when derived for the selected
regime,

- ion-neutral friction;
- ambipolar diffusion;
- Hall/Pedersen conductivity;
- anisotropic resistivity;
- neutral-wind/electrodynamic coupling.

These are closures of the common species equations. They do not replace the
kinetic charged/neutral path and must expose their reduction assumptions.

## Multi-temperature thermodynamics

The state may distinguish electron, ion, neutral and internal-mode
temperatures, for example

\[
T_e,\quad T_i,\quad T_n,\quad T_{\rm vib},\quad T_{\rm rot}.
\]

Energy exchange between them is explicit. One-temperature LTE is a selectable
limit, not a universal assumption.

EOS, transport and reaction data may come from external tables with hashes,
units and uncertainty metadata.

## Electrodes, dielectrics and surfaces

Extend the common surface operator with:

- conducting/floating electrodes;
- dielectric permittivity and conductivity;
- evolving surface charge;
- secondary-electron emission;
- thermionic, field and photo-emission;
- adsorption/desorption and surface reaction coverages;
- sputtering, deposition and etching yields.

A dielectric surface may evolve

\[
\partial_t \sigma_s
=
J^{\rm plasma}_n-J^{\rm material}_n
\]

with the field boundary condition using the same surface charge.

Incoming energy-angle-resolved fluxes are passed to surface chemistry through a
generic interface. Detailed solid chemistry/morphology can remain an external
provider.

## Adjacent-media coupling

Plasma-liquid or detailed material applications use a conservative interface
port rather than importing an entire liquid/solid solver into KIN6D.

The interface exchanges declared particle/species fluxes, current, heat,
potential and reaction products. External media return boundary state and
fluxes with their own numerical/model uncertainty.

## Generic validation ladder

Start from simple reusable cases:

1. homogeneous reaction network with exact equilibrium/relaxation;
2. electron swarm or reduced-Boltzmann coefficient benchmark;
3. 1D drift-diffusion/Poisson front;
4. dielectric surface charging;
5. neutral-flow/reaction coupling;
6. kinetic-electron versus reduced-EEDF comparison;
7. cross-domain physical case.

Application-specific databases and geometries enter only after the primitive
works.

Chris&AI
