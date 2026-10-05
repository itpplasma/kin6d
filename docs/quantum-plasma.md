# Quantum plasma fundamentals and dense-matter bridge

Status: canonical long-range scope. Quantum-plasma capabilities are reference,
teaching and closure-validation branches unless a prioritized scientific case
promotes them. They do not make every KIN6D calculation quantum.

## Scope boundary

KIN6D extends beyond the purely classical hierarchy far enough to expose and
test the quantum assumptions that matter in plasmas:

- Fermi--Dirac degeneracy and relativistic Fermi statistics;
- quantum diffraction/exchange corrections to close encounters;
- ideal and interacting electron gases;
- Thomas--Fermi and Lindhard/RPA screening;
- quantum hydrodynamic closures when derived;
- Wigner kinetic equations and their classical limit;
- Hartree/Hartree--Fock and reduced-density-matrix references;
- quantum BBGKY/cumulant hierarchies;
- quantum collision/dielectric closures including Pauli blocking;
- equilibrium quantum many-body reference methods such as exact
  diagonalization and quantum/path-integral Monte Carlo where practical;
- DFT/TDDFT or equivalent orbital-density methods as microscopic
  quantum-matter providers/references.

The deliberate fundamental boundary is **QED**. KIN6D does not attempt a
general quantum-field-theoretic electron/positron/photon solver. Pair creation,
strong-field QED cascades, radiative QED corrections and fully dynamical
quantized radiation fields are external specialties.

Relativistic Fermi equations of state, relativistic single-particle kinematics
and controlled semirelativistic/Dirac-Coulomb references may still be used
below that boundary when pair creation and quantum radiation are negligible.

## Regime coordinates

The quantum/classical and weak/strong-coupling distinctions are independent.

For mean interparticle spacing

\[
a=\left(\frac{3}{4\pi n}\right)^{1/3},
\]

define the Fermi energy

\[
E_F=
\frac{\hbar^2}{2m_e}(3\pi^2n_e)^{2/3}
\]

in the nonrelativistic limit and the degeneracy parameter

\[
\Theta=\frac{k_BT}{E_F}.
\]

Classical Coulomb coupling is measured by

\[
\Gamma=
\frac{e^2}{4\pi\epsilon_0ak_BT},
\]

while in a degenerate system a useful interaction ratio is

\[
\Gamma_q=
\frac{e^2/(4\pi\epsilon_0a)}{E_F}
=
\Gamma\Theta.
\]

Relativistic degeneracy is indicated by

\[
x_F=\frac{p_F}{m_ec}.
\]

This gives a clean regime map:

- \(\Theta\gg1,\Gamma\ll1\): classical ideal/weakly coupled plasma;
- \(\Theta\gg1,\Gamma\gtrsim1\): classical strongly coupled plasma;
- \(\Theta\lesssim1,\Gamma_q\ll1\): ideal/weakly interacting degenerate Fermi
  plasma;
- \(\Theta\lesssim1,\Gamma_q\gtrsim1\): correlated quantum plasma/warm dense
  matter.

"Ideal plasma" here means weak interaction relative to the relevant kinetic
energy. It is unrelated to the separate ideal-MHD approximation.

## Quantum statistical hierarchy

The microscopic quantum parent can be represented by a many-body wavefunction
or density operator,

\[
\Psi(x_1,\ldots,x_N,t),
\qquad
\rho_N.
\]

Reduced density matrices satisfy a quantum BBGKY hierarchy. The practical
reduction ladder is

\[
\boxed{
\rho_N
\rightarrow
(\rho_1,\rho_2,\ldots)
\rightarrow
\text{quantum cumulants}
\rightarrow
\text{Hartree/HF/RPA/kinetic closures}
\rightarrow
\text{Wigner or quantum-fluid models}.
}
\]

This parallels the classical

\[
F_N\rightarrow(f_1,g_2,g_3,\ldots)
\]

hierarchy and should reuse generic low-rank/tensor/cumulant infrastructure
where mathematically appropriate.

## Wigner kinetic bridge

The Wigner quasidistribution \(f_W(x,p,t)\) is the most direct quantum analogue
of the classical phase-space distribution. Its evolution replaces the
classical force/Poisson bracket by the corresponding Moyal/nonlocal quantum
operator.

The required limiting test is

\[
\hbar\rightarrow0
\quad\Longrightarrow\quad
\text{Wigner}\rightarrow\text{Vlasov}.
\]

This makes existing Eulerian, semi-Lagrangian, DG, spectral, tensor and
physical-space FEM infrastructure directly relevant to quantum plasma
fundamentals.

Wigner functions may be negative. Positivity semantics from classical
probability densities must not be copied blindly.

## Screening hierarchy

A particularly useful common observable is static screening. In a general
electron gas,

\[
k_s^2
=
\frac{e^2}{\epsilon_0}
\frac{\partial n}{\partial\mu}.
\]

The classical limit recovers Debye screening, while a degenerate ideal Fermi
gas gives Thomas--Fermi screening. The dynamic/nonlocal Lindhard/RPA response
adds wave-number and frequency dependence and produces Friedel oscillatory
structure beyond a simple exponential model.

The verification sequence is therefore

\[
\boxed{
\text{Debye}
\rightarrow
\text{finite-}T\text{ Fermi screening}
\rightarrow
\text{Thomas--Fermi}
\rightarrow
\text{Lindhard/RPA}
\rightarrow
\text{correlated quantum reference}.
}
\]

## Quantum collision hierarchy

Quantum effects can matter in collisions even when the bulk distribution is
nondegenerate. The hierarchy should distinguish:

- quantum-diffraction short-distance cutoff/cross sections;
- exchange symmetry;
- Pauli blocking;
- Uehling--Uhlenbeck-like kinetic corrections;
- quantum Lenard--Balescu/dielectric closures;
- mean-force/correlated quantum closures.

A hot tokamak may therefore be classically nondegenerate while still requiring
quantum diffraction in the short-distance collision physics.

## Numerical representations

Quantum plasma does **not** imply particle-only computation.

### Exact wavefunction/density matrix

Direct \(3N\)-dimensional \(\Psi\) or \(6N\)-dimensional density-matrix
representations are gold-standard tiny-\(N\) references. They scale
exponentially and are not macroscopic methods.

### Orbital methods

Hartree, Hartree--Fock, Kohn--Sham DFT and time-dependent variants evolve
ordinary 3D orbitals/densities. Finite elements, plane waves, spectral bases,
finite differences and wavelets are legitimate representations.

Initially these may be independent providers/oracles; internalize only the
subset demanded by a plasma case.

### Quantum Monte Carlo and path integrals

Variational/diffusion/path-integral Monte Carlo sample many-body configuration
or imaginary-time path space. They are valuable equilibrium EOS,
pair-correlation, structure-factor and free-energy references. Their
statistical/sign-problem status is explicit.

### Wigner continuum

Use Eulerian, spectral, DG, semi-Lagrangian, tensor-product, sparse-grid or
low-rank phase-space methods. This is the quantum branch most directly aligned
with KIN6D's Vlasov machinery.

### Quantum hydrodynamics and EOS closures

At macroscopic scales use Fermi pressure, exchange/correlation and derived
gradient/Bohm-like terms only under their declared closure assumptions. These
ordinary physical-space PDEs can use FEM/FEEC/FV/DG and can cover whole stars
or devices without resolving microscopic quantum wavelengths.

## Adaptive classical--quantum coupling

Quantum fidelity may be spatially and temporally local.

In a hot nondegenerate plasma, most dynamics may remain classical while rare
close encounters require quantum diffraction/scattering. An event-local
two-body quantum microproblem can be used when its validity indicator is met,
without promoting the entire plasma.

Degeneracy and exchange are nonlocal in occupation space. If
\(\Theta=T/T_F\lesssim O(1)\), pair-only promotion is generally inadequate;
promote an appropriate region/species/model to Fermi/Wigner/HF/quantum kinetic.

Demotion back to a classical state is allowed only with a declared
coarse-graining/decoherence/information-loss model and propagated uncertainty.
Persistent phase/coherence information cannot be silently discarded.

Quantum model error should preferentially be estimated through asymptotic
\(\hbar\)-expansions, classical-limit defects, richer-model residuals and sparse
quantum reference calculations. Probabilistic model discrepancy is a fallback,
not the default.

See [adaptive-multiscale-uq.md](adaptive-multiscale-uq.md).

## Expected computational scales

The numbers below are planning scales, not performance claims.

- direct real-space many-body wavefunction: roughly 2--few particles in 3D at
  useful grid resolution;
- exact diagonalization: tiny basis-dependent particle counts;
- quantum/path-integral Monte Carlo: typically tens to thousands of particles
  for microscopic equilibrium reference boxes, strongly problem dependent;
- HF/DFT/TDDFT: microscopic cells with \(10^2\)--\(10^4\)-scale
  particles/orbitals are plausible targets depending on method, basis and
  hardware;
- Wigner: local 1X1P/2X2P demonstrations are easy; 3X3P has the same
  high-dimensional pressure as classical FK and needs structure, low rank,
  reduced geometry or HPC;
- quantum fluid/Fermi-EOS models: macroscopic domains up to whole-star scales,
  because microscopic quantum physics enters through closures.

Do not turn these planning ranges into support claims before benchmarks exist.

## Stellar-interior connection

Stellar interiors use the same hierarchy by changing representation with
scale.

A low-cost whole-star model combines

\[
\text{self-gravity}
+
\text{compressible plasma}
+
\text{EOS}
+
\text{radiative/conductive/convective transport}
+
\text{nuclear reaction network}.
\]

Quantum plasma supplies the EOS/conductivity/screening corrections; it does
not require a star-sized wavefunction.

A staged EOS path is

\[
\text{classical ideal gas}
\rightarrow
\text{finite-}T\text{ Fermi gas}
\rightarrow
\text{Coulomb/exchange corrections}
\rightarrow
\text{partial ionization/nonideal EOS}.
\]

This covers main-sequence/solar-interior and white-dwarf-like plasma physics
within the non-QED/non-neutron-star regime. Neutron-star interiors require
nuclear-density matter, strong-interaction physics, superfluidity/
superconductivity and GR and remain outside KIN6D's intended ownership.

General stellar burning may extend the reaction graph beyond D--T/DD/D3He to
pp chains, CNO, triple-alpha and larger networks. Neutrino transport and
advanced late-stage weak-interaction networks are separate high-energy
specialties and need not be internalized merely to support ordinary stellar
plasma fundamentals.

## Teaching/verification gallery

High-value examples are:

1. Maxwell--Boltzmann versus Fermi--Dirac occupation;
2. Fermi sphere and \(E_F\);
3. classical thermal pressure versus nonrelativistic and relativistic
   degeneracy pressure;
4. the \((\Theta,\Gamma,x_F)\) regime diagram with tokamak, solar-core,
   warm-dense and white-dwarf-like examples;
5. Debye \(\rightarrow\) Thomas--Fermi \(\rightarrow\) Lindhard screening and
   Friedel oscillations;
6. Wigner evolution converging to Vlasov as \(\hbar\to0\);
7. quantum-degenerate plasma oscillation/dispersion;
8. Pauli blocking and quantum-diffraction effects on collision rates;
9. a two-electron exchange/correlation reference;
10. a reduced stellar-structure example showing thermal-to-degeneracy-pressure
    crossover.

These should use the same result/error/provenance contract as classical
examples.

## Relationship to the QED boundary

The hierarchy intentionally stops before general quantum field theory. If the
dominant physics requires creation/annihilation of particles, quantized
electromagnetic radiation, strong-field nonlinear QED or radiative QED loop
corrections, the problem is beyond KIN6D's native fundamentals branch.

The provider boundary may still pass QED-derived rates or closure data into a
classical/quantum plasma calculation when that is scientifically useful.

Chris&AI
