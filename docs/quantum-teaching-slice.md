# Quantum teaching slice

This document records the completed periodic classical/quantum pair quench
and ideal Fermi-gas pressure reference. The quench is a correlation experiment,
not incoming wave-packet scattering or a degenerate plasma demonstration.
The next teaching deliverables and their acceptance gates are defined in
[microscopic quantum validation](quantum-microscopic-validation.md).
Quantum BBGKY evolution and stellar structure remain deferred.

## Matched interacting experiment

Use two equal-mass, equal-charge electrons in a periodic three-dimensional
box with a uniform neutralizing background. The quantum preparation has one
spin-up and one spin-down electron in the zero-momentum orbital. Its initial
one-body Wigner distribution is positive: position is uniform and momentum is
zero. The classical ensemble uses the same position and momentum preparation.
Both start with an uncorrelated unlike-spin pair distribution, and the same
bare Coulomb coupling is switched on at time zero.

The quantum solver evolves the full interacting, antisymmetric finite-basis
many-body state. It does not impose a screened interaction, Hartree
factorization, or collision operator. The classical solver evolves actual
three-dimensional particle pairs with periodic Ewald Coulomb forces. An
ensemble supplies a distribution of trajectories; it does not turn each
physical particle into a PIC superparticle.

Compare the relative-x probability in sixteen equal bins, integrating over
the complete transverse torus. This is a three-dimensional slab average, not
a one-dimensional interaction model. Normalize each bin to the initial
uniform law, so both population distributions begin at g=1. The initial
classical histogram retains its finite-sample fluctuations; it is not
renormalized to hide them. Quantum densities are probabilities;
classical dots are trajectories. Do not draw trajectories for quantum
particles or interpret the quantum density as a classical cloud of markers.

Two opposite spins isolate interacting quantum diffraction/correlation from
same-spin exchange. The general fermionic Hamiltonian retains exchange; a
same-spin example can follow after this case is validated. A two-electron
quench is a microscopic reference, not a demonstration of thermodynamic Debye
shielding or a finite-temperature equilibrium plasma.

## General interfaces and explicit budgets

The Hamiltonian accepts supplied one- and two-body coefficients and arbitrary
fixed-N initial amplitudes. The physical constructor uses three-dimensional
plane-wave momenta, spin, box size, mass, hbar, and Coulomb coupling. The
uniform background removes the zero Fourier mode. Orbital cutoffs and the
occupation-bit backend's capacity are numerical representation limits.
They are not definitions of the physical parent model.

The teaching pair also has an exact total-momentum-zero, one-up/one-down
sector constructor. Its physical ordered states are |k up, -k down>;
this symmetry reduction avoids the occupation-bit limit and a dense
four-index interaction tensor. The default uses 125 momenta, with a
343-momentum convergence audit. Their late-time slab difference is about
0.0052 on the initial g=1 scale; this is an indicator, not a continuum bound.

The evolution interface accepts arbitrary complex Hermitian matrices. Primal
and directional evolution share generated actions, with equation-level
adjoints. Fixed storage and explicit Hilbert-space/memory/time-step budgets
prevent silent growth. Small assembled matrices are reference schedules;
large sparse or matrix-free production evolution requires separate measured
work. Targeted pair annihilation avoids materializing the full four-index
pair density when only a pair observable is needed.

Wick products and connected correlations are observables and optional
factorizations, not an enforced dynamics closure. Fixed-N pure-state tests
and fractional grand-canonical occupations have different contraction
identities. Do not apply a fixed-N identity to a finite-temperature ideal-gas
Wick state.

## Error band scope

A unitary residual argument bounds continuous time evolution for the declared
finite-basis Hermitian Hamiltonian. Strict reconstruction and observable
lifting account for numerical roundoff. Enclosed physical coefficient
assembly can contribute a separate Hamiltonian perturbation allowance.
Spatial curves need a full plot-cell interpolation allowance; node intervals
alone do not certify the connecting curve.

The orbital-cutoff difference is a convergence indicator until a continuum
projection remainder has been proved. It cannot be hidden in a time-step
certificate. Classical Ewald cutoff/time-step checks are convergence evidence
unless a deterministic remainder is supplied. Classical ensemble confidence
intervals describe statistical uncertainty and are never deterministic
numerical error bounds. Report solve-only and solve-plus-check cost separately
at the same input and accuracy target; exclude rendering from both timings.

The reviewed local Release run records 1.155 s for the matched quantum
candidate initialization/evolution and 104.967 s including certification and
all plotted bands: 90.89 times solve alone. The adaptive first result costs
157.35 s before rendering. Its maximum numerical display radius is 0.002757
on scale 1 at 125 modes and 16 steps per frame. Classical evolution costs
33.508 s plus 2.318 s for sampling/diagnostics, a factor of 1.069; that extra
work is not a deterministic trajectory certificate. These are single CPU
trials with retained ledgers, not a cross-code speed ranking. Dense interval
residual actions remain a measured optimization target; the fast production
workhorse goal is not yet achieved for this quantum reference.

## Ideal electron pressure

Solve the finite-temperature ideal spin-two Fermi occupation numerically at
fixed density. The chemical potential follows from the density equation;
pressure and compressibility follow from momentum integrals. Zero-temperature
and classical formulas are independent limiting-case references, not the
production answer or certificate input. Strict quadrature, tail bounds, and
root bracketing determine the enclosure. Implicit differentiation follows the
fixed-density equation, rather than solver iteration history.

An illustrative low-density stellar-envelope electron pressure is allowed
through the same SI normalization. Declare density, electron fraction,
temperature, pF/(m_e c), and the ideal nonrelativistic model. This is an
electron-pressure calculation, not a star model. Relativistic stellar cores,
ion pressure, interaction corrections, and hydrostatic structure are deferred.

At mass density 10^6 kg/m^3, electron fraction 0.5 and temperature 10^6 K,
the numerical ideal-electron pressure is 3.19695e16 Pa, with outward-rounded
enclosure [3.1931e16, 3.2008e16] Pa. Here Theta=0.05260 and pF/(m_e c)=0.08007.
The benchmark CSV retains the full endpoints. Central CODATA constants and the
declared binary64 inputs are used, without measurement-uncertainty propagation.

## Hot encounters and degeneracy are different questions

Bulk degeneracy depends on occupation: compare k_B T with E_F, or density
with the thermal wavelength. Quantum diffraction in an individual encounter
instead compares the wave-packet scale with spatial variation of the
Coulomb field. High temperature does not guarantee classical close encounters.
At 10 keV, the reduced thermal wavelength hbar/sqrt(m_e k_B T) is about
2.76 pm, while the Landau length e^2/(4 pi epsilon_0 k_B T) is about 0.144 pm.
These definitions are explicit; other wavelength/impact-parameter conventions
have different factors. They do not imply that a dilute tokamak plasma is
bulk degenerate or that every encounter has a large quantum correction.

A primary discussion of this short-range regime is
[Alharbi et al., PTEP 2026, 023J01](https://doi.org/10.1093/ptep/ptaf188).
Its mean-field model is not a certificate for the finite-N experiment here.
For plane-wave electronic Hamiltonians see
[Babbush et al., Physical Review X 8, 011044](https://doi.org/10.1103/PhysRevX.8.011044).
Fermi-integral conventions are given in
[NIST DLMF 25.12(iii)](https://dlmf.nist.gov/25.12#iii).
