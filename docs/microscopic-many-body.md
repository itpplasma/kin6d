# Microscopic relativistic many-body reference

Status: canonical microscopic foundation. This is a reference hierarchy beneath
the continuum kinetic models; it is not the default production representation
for macroscopic plasma calculations.

## Classical parent model

The deepest classical model in KIN6D is the **multispecies relativistic
Maxwell--Lorentz many-body system**: individual particle worldlines coupled to
the electromagnetic field.

For particle \(i\) of species \(s\),

\[
\dot{\mathbf x}_{si}
=
\frac{\mathbf p_{si}}{\gamma_{si}m_s},
\qquad
\gamma_{si}
=
\sqrt{1+\frac{|\mathbf p_{si}|^2}{m_s^2c^2}},
\]

\[
\dot{\mathbf p}_{si}
=
q_s
\left[
\mathbf E(\mathbf x_{si},t)
+
\mathbf v_{si}\times\mathbf B(\mathbf x_{si},t)
\right],
\]

with microscopic sources

\[
\rho(\mathbf x,t)
=
\sum_{s,i}q_s\delta(\mathbf x-\mathbf x_{si}),
\qquad
\mathbf j(\mathbf x,t)
=
\sum_{s,i}q_s\mathbf v_{si}
\delta(\mathbf x-\mathbf x_{si}),
\]

and the full Maxwell equations.

This is the classical parent reference. Relativistic Vlasov--Maxwell/FK is the
highest **continuum mean-field** reference, not the deepest microscopic model.

## Covariant action and the self-field caveat

The formal covariant particle-field action is

\[
S
=
-\sum_{s,i}m_sc\int ds_{si}
-\sum_{s,i}q_s\int A_\mu\,dx^\mu_{si}
-\frac{1}{4\mu_0}\int F_{\mu\nu}F^{\mu\nu}\,d^4x.
\]

Variation yields Maxwell plus relativistic Lorentz dynamics. For point charges,
however, the self-field is singular. KIN6D must therefore distinguish:

- exact mutual interparticle electromagnetic interaction;
- independent/free electromagnetic field data;
- the selected self-force/radiation-reaction prescription.

Ordinary tokamak collision/correlation validation normally needs mutual
interaction but not radiation reaction. When radiation reaction matters,
Landau--Lifshitz or another explicitly named regularized/renormalized model is a
separate physical choice. Do not call a numerically softened point charge the
same microscopic model.

## Electrodynamic fidelity ladder

The microscopic hierarchy is

\[
\boxed{
\text{relativistic Coulomb/electrostatic N-body}
\rightarrow
\text{relativistic Darwin N-body}
\rightarrow
\text{full Maxwell--Lorentz N-body}
\rightarrow
\text{direct retarded Lienard--Wiechert reference}
\rightarrow
\text{radiation-reaction extension}.
}
\]

Coulomb and Darwin are reduced models and valuable benchmarks. The target
endpoint remains full relativistic Maxwell--Lorentz.

Direct Lienard--Wiechert evaluation is a derived field-eliminated formulation
under a stated retarded/no-incoming-radiation condition. It is not more
fundamental than the Maxwell initial-value problem, because the latter may also
contain independent radiation.

## Statistical-reduction hierarchy

For the exact microscopic empirical distribution,

\[
f_K(z,t)
=
\sum_i\delta(z-z_i(t)),
\]

the Klimontovich equation is exact. Ensemble averaging introduces fluctuation
correlations. Equivalently, the \(N\)-particle Liouville density generates the
BBGKY hierarchy.

Write

\[
f_2(1,2)
=
f_1(1)f_1(2)+g_2(1,2),
\]

and similarly higher cumulants \(g_3,\ldots\). The statistical reduction ladder
is therefore

\[
\boxed{
F_N
\rightarrow
(f_1,g_2,g_3,\ldots)
\rightarrow
\text{BBGKY/cumulant closures}
\rightarrow
\text{Balescu--Lenard/mean-force/Landau-like kinetic closures}
\rightarrow
\text{Vlasov/FK mean field}.
}
\]

This is not a perfectly linear accuracy ladder; different closures take
different asymptotic limits. The microscopic solver exists to test those
limits directly.

## Quantum-many-body connection

The classical correlation hierarchy may borrow numerical ideas from quantum
many-body physics, and the dedicated quantum branch uses the corresponding
reduced-density-matrix hierarchy directly:

- reduced density/marginal hierarchies;
- cumulant/cluster expansions;
- pair-correlation and Jastrow-like ansatzes;
- dynamic low-rank factorizations;
- tensor trains/hierarchical Tucker;
- Green-function/response formulations;
- center-of-mass/relative-coordinate localization.

For example,

\[
g_2(z_1,z_2)
\approx
\sum_{\alpha=1}^r
a_\alpha(z_1)b_\alpha(z_2)
\]

may make an explicit pair-correlation model tractable. Such compression carries
an explicit truncation defect; it does not magically preserve arbitrary
\(N\)-body information.

## Exact particles versus macro-particles

A computational marker with charge \(wq\) representing \(w\gg1\) physical
particles changes microscopic fluctuation amplitudes, close collisions and pair
correlations. Therefore

\[
\boxed{
\text{PIC/PIF macro-particles are not exact microscopic coarse graining.}
}
\]

Likewise a finite blob changes short-range interactions unless its form factor
is part of the declared microscopic model.

To preserve the physical finite-\(N\) model, retain the physical particles and
accelerate **force evaluation**, not particle count.

## Fast interaction algorithms

For electrostatic/periodic microscopic calculations retain a hierarchy:

1. direct \(O(N^2)\) pair sum as the gold-standard small-\(N\) oracle;
2. FMM/treecode with controlled force error;
3. Ewald/spectral-Ewald/P3M for periodic long-range Coulomb interaction;
4. near/far decompositions that keep exact/controlled near encounters while
   accelerating smooth far fields.

For full electrodynamics, large-\(N\) production should normally evolve Maxwell
fields locally rather than evaluate all retarded pairs. Direct retarded
Lienard--Wiechert interaction remains a small-\(N\) oracle. Space-time FMM,
hierarchical retarded kernels or exact near-pair plus field-mediated far
interaction are research paths when they beat the field solve at matched error.

An FMM approximation with a verified tolerance changes numerical error, not the
physical particle model. A macro-particle changes the physical finite-\(N\)
model.

## Multispecies from the beginning

The microscopic state is natively multispecies through \((m_s,q_s,N_s)\).
Electron-electron, electron-ion, ion-ion, impurity, pair-plasma and energetic
minority interactions require no separate collision formula.

External prescribed fields, walls or radiation may be added through the same
coupling contracts as the continuum models.

## Microscopic observables

Do not validate only temperature relaxation. Important microscopic observables
include:

- two-body deflection/cross sections;
- pair correlation \(g_{ab}(r)\);
- static structure factor \(S_{ab}(k)\);
- dynamic structure factor \(S_{ab}(k,\omega)\);
- dielectric response \(\epsilon(k,\omega)\);
- force/velocity autocorrelation functions;
- Green--Kubo diffusion, conductivity, viscosity and heat conductivity;
- slowing-down and pitch-angle diffusion;
- multispecies momentum/energy equilibration;
- magnetized collision and finite-Larmor corrections.

These observables connect microscopic dynamics to kinetic collision and
transport closures.

## Microscopic teaching/verification ladder

Retain especially transparent examples:

### M0 — two particles, no field

Relativistic two-body Coulomb scattering with controlled impact parameter.
Compare conserved quantities and the nonrelativistic Rutherford limit.

### M1 — two particles in prescribed electromagnetic field

Resolve gyro-motion plus mutual interaction. This exposes the transition from
an isolated binary encounter to a magnetized collision without invoking a
collision operator.

### M2 — retarded two-particle electrodynamics

For a small prescribed case, compare Maxwell-field evolution and direct
retarded Lienard--Wiechert mutual fields under matching radiation/boundary
conditions.

### M3 — microscopic plasma oscillation

Initialize many physical electrons/ions in a small periodic box and recover the
collective plasma-frequency response with no Vlasov equation imposed.

### M4 — real Debye shielding

Insert or hold a weak test charge in a thermal multispecies microscopic plasma.
Measure the ensemble/time-averaged response and compare, in the weak-coupling
linear regime, with

\[
\phi(r)
=
\frac{q}{4\pi\epsilon_0r}
\exp(-r/\lambda_D).
\]

The exponential screening must **emerge from the particles**. Do not impose a
Debye dielectric, collision operator or screened pair potential in this test.
Finite-box, periodic-image, finite-\(N_D\), averaging and test-charge
nonlinearity errors are reported separately.

### M5 — correlation and dielectric spectrum

Measure \(g(r)\), \(S(k)\), \(S(k,\omega)\) and dielectric response, then compare
against Vlasov/RPA/Balescu-type predictions in their valid limits.

### M6 — microscopic collision coefficients

Measure slowing down, pitch-angle diffusion, temperature/momentum relaxation
and Green--Kubo transport; compare against Landau, Beliaev--Budker,
Balescu--Lenard and mean-force closures.

### M7 — strong-coupling departure

Increase coupling until weak-correlation kinetic closures fail and compare
explicit pair-correlation/mean-force models against direct particles.

## Computational scope

The microscopic branch is intentionally local. It need not represent a whole
tokamak. Its purpose is to resolve a statistically meaningful volume over
microscopic lengths/times and connect those observables upward to continuum
closures.

For weakly coupled fusion-core parameters, one Debye-scale volume already
contains very many physical particles. Therefore the development ladder should
use:

- tiny direct two-/few-particle tests;
- controlled finite-\(N_D\) boxes;
- convergence in plasma parameter/coupling;
- FMM/Ewald/P3M or near/far acceleration;
- larger physical-particle boxes only when the observable requires them.

Never claim that a low-\(N_D\) surrogate is the same microscopic tokamak plasma.
Use it to test asymptotic scaling and closure convergence.

## Relationship to the quantum branch

Within **classical** plasma physics, Maxwell--Lorentz many-body dynamics is the
bottom of the hierarchy. KIN6D now continues into quantum plasma fundamentals
for degeneracy, exchange/diffraction, Wigner/quantum kinetic models,
reduced-density-matrix/cumulant methods and microscopic quantum references.

In hot fusion plasmas, quantum diffraction can set the short-distance cutoff
for rare close electron encounters even while collective Debye-scale dynamics
is classical. The classical microscopic branch therefore validates the
classical closure and its asymptotic structure; the quantum branch quantifies
where that closure fails.

The common native boundary is QED: pair creation/annihilation, quantized
radiation dynamics, strong-field QED cascades and radiative QED corrections are
external. See [quantum-plasma.md](quantum-plasma.md).

Chris&AI
