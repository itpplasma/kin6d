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

The electromagnetic fidelity hierarchy is independent of implementation
priority. Its first useful teaching limit is nonrelativistic instantaneous
Coulomb dynamics, followed by relativistic electrostatics, Darwin and the full
Maxwell model when the case needs them. Coulomb admission requires declared
small \(v/c\), small \(\omega L/c\), and negligible radiative loss on the
observation interval; these statements need a model-remainder status.

The microscopic hierarchy is

\[
\boxed{
\text{nonrelativistic Coulomb N-body}
\rightarrow
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

and similarly higher cumulants \(g_3,\ldots\), with a declared marginal
normalization (fixed-number normalized marginals also have finite-\(N\)
combinatorial corrections). These are a branching family of reductions:

\[
\boxed{
F_N
\rightarrow
(f_1,g_2,g_3,\ldots)
\rightarrow
\left\{
\begin{array}{l}
\text{mean-field limit: Vlasov/continuum FK},\\
\text{weak-correlation closures: Balescu--Lenard/Landau/BB},\\
\text{retained correlations: cumulant/mean-force models}.
\end{array}
\right.
}
\]

This is not a perfectly linear accuracy ladder; different closures take
different asymptotic limits. The microscopic solver exists to test those
limits directly.

Collective polarization, plasma oscillation and collisionless phase mixing
already occur in self-consistent Vlasov theory. They do not require a binary
collision operator. A resolved two-body encounter is a different observation
from collective mean-field motion. Landau describes accumulated small-angle
encounters in its kinetic limit; Balescu--Lenard also retains dielectric
dressing of weak correlations. Neither resolves arbitrary irreducible
three-/higher-particle correlations. BBGKY is the hierarchy, not itself a
particular effective collision operator. Do not add a collision operator to
resolved interacting particles without deriving a split and preventing double
counting.

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

The IDs below identify cases, not prerequisite order. First implement M0,
the Coulomb form of M1, and the Coulomb collective M3/M4 cases. M2 and the
relativistic/radiative extensions follow concrete validation needs. Every
microscopic case remains planned until its own numerical and physical gates
pass; the current free-streaming markers have no mutual interaction.

### M0 — two particles, no field

Begin with nonrelativistic repulsive two-body Coulomb scattering at a declared
impact parameter, using the Rutherford deflection and momentum/energy ledger
as independent oracles. Animate both trajectories in their actual 3D space.
Extend to relativistic kinematics later. Close-approach events require a
resolved/adaptive time and force error budget, not hidden numerical softening.

### M1 — two particles in prescribed electromagnetic field

Resolve gyro-motion plus mutual interaction. This exposes the transition from
an isolated binary encounter to a magnetized collision without invoking a
collision operator.

Start with a few interacting 3D3V particles, then add the prescribed-field
limit and compare it with the independent single-particle orbit reference.

### M2 — retarded two-particle electrodynamics

For a small prescribed case, compare Maxwell-field evolution and direct
retarded Lienard--Wiechert mutual fields under matching radiation/boundary
conditions.

### M3 — microscopic plasma oscillation

Start with physical electrons and an explicit uniform positive neutralizing
background (an electron one-component plasma, OCP) in a periodic 3D3V box.
Use a small nonzero spatial-mode displacement, recover the collective response
and compare with the matching Vlasov/continuum and analytical cold/thermal
limits. A uniform translation in a periodic OCP is not a charge-separation
oscillation. Mobile ions are a later distinct population model. No effective
collision operator is imposed.

For \(L=4\lambda_D\), the smallest nonzero wave number has
\(k\lambda_D=2\pi/4\), so a warm box is not the long-wavelength
\(\omega\simeq\omega_p\) limit. Compare its full finite-\(k\)
dielectric/initial-value response, or increase \(L/\lambda_D\), or declare
a separate cold/nonthermal model. Do not use the small shielding box's size
as evidence that the cold Langmuir oracle applies to a warm microscopic run.

### M4 — real Debye shielding

Begin with that OCP and a weak, fixed, repulsive impurity charge, compensating
its charge in the declared background to maintain periodic neutrality. This
avoids classical attractive electron--ion bound collapse. A positive finite
charged sphere is a separate finite-source benchmark; its radius and source
law belong to the physical model. Mobile positive point ions and attractive
close encounters require their own classical validity and equilibrium rules.

Insert or hold a weak test charge in a thermal microscopic plasma.
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

For a point impurity, the linear response is a far-field claim, outside the
nonlinear radius \(|qQ|/(4\pi\epsilon_0 k_BT)\). At finite periodic box size,
use the corresponding periodic screened-Poisson Green function as the
reference; the isolated Yukawa curve above is an overlay with a finite-box
model-difference status. A \(4\lambda_D\) box is not an infinite-domain
shielding experiment. Specify which species respond when defining
\(\lambda_D^{-2}=\sum_{\text{responsive }s}n_sq_s^2/(\epsilon_0 k_BT_s)\).
Electron-only and equilibrated electron/ion responses have different lengths.

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

### Debye-volume scan and separate limits

Use the explicit conventions

\[
\Lambda_s=n_s\lambda_D^3,\qquad
N_{D,s}=\frac{4\pi}{3}\Lambda_s,\qquad
N_s=\Lambda_s(L/\lambda_D)^3.
\]

The responsive species and Debye-length convention are fixed for each scan.
Start with \(\Lambda=100,300,1000,3000,10000\), at fixed normalized domain,
temperature ratios, mass ratios and observation times. Add smaller values
only as explicitly non-asymptotic correlation examples. Vary physical
density/temperature/charge consistently to obtain this coupling scan; assigning
larger weights to fewer markers is not the same finite-\(N\) experiment.

At \(L=4\lambda_D\), the counts per species are respectively
6400, 19200, 64000, 192000 and 640000. These are planning counts, not timing
evidence or a claim that direct summation is cheap. At 64000 particles,
an isolated unordered direct pair action already visits about
\(2.05\times10^9\) pairs per step, before periodic-image and checking work.
Calibrate force, near-encounter, step, observable and certificate costs before
choosing a default. Start with smaller boxes and lower counts for mechanics,
then promote the several-Debye-length case with a verified long-range action.

Maintain separate scan axes and acceptance gates:

| Gate | Measured quantity | Interpretation |
|---|---|---|
| Collective mean response | Mean potential, screening length, plasma frequency, dielectric response | Test convergence toward the matching continuum limit |
| Finite-particle fluctuations | Density variance, force fluctuations, normalized correlation/structure spectra | Preserve the declared finite-\(N\) physical dependence |
| Effective collisions | Drift/diffusion, short-time increments and correlation integrals | Compare normalized rates under the closure's hypotheses |
| Geometry and boundaries | \(L/\lambda_D\), periodic images, radial observation window | Separate finite-domain effects from the coupling scan |
| Estimator and averaging | Independent ensembles, autocorrelation, averaging time | Separate uncertainty/bias from deterministic numerical error |

The rough independent-count fluctuation scale \(\Lambda^{-1/2}\) is a
teaching estimate, not a rigorous confidence interval for an interacting
plasma. Temporal samples are correlated; radial bins do not automatically
provide independent observations. State the estimator and effective-sample
assumptions, with statistical uncertainty separate from certified force and
integration errors. The actual averaged potential, not a frame's particle
picture, is compared with the equilibrium screening prediction. How the
initial ensemble is prepared and whether an equilibration thermostat is used
must be recorded; a thermostat changes the dynamical model.

Mean response may converge well before fluctuation amplitudes or collision
rates match the target plasma. The schematic weak-coupling rate is
\(\nu/\omega_p\sim C\ln\mathcal L/\Lambda\), where the Coulomb logarithm
\(\ln\mathcal L=\ln(b_{\max}/b_{\min})\), its cutoffs and \(C\) depend on
the species/model/convention. Replacing it by \(\ln\Lambda\) is an
asymptotic scaling illustration, not a calibrated physical collision rate.
Likewise an \(O(\Lambda^{-1})\) finite-coupling correction must be derived
or checked for the chosen mean observable; fitting it is not a remainder
certificate. Normalized structure factors can retain a nontrivial screened
limit: do not claim that all correlations or all definitions of \(S\) vanish.

For context, the declared \(n_e=10^{20}\,\mathrm{m}^{-3}\),
\(k_BT_e=10\,\mathrm{keV}\) electron-only case has approximately
\(\lambda_{De}=74.3\,\mu\mathrm m\),
\(\Lambda_e=4.1\times10^7\) and \(N_{D,e}=1.7\times10^8\).
These are reference scales, not parameters of the small teaching box.
Low-\(\Lambda\) results can validate an asymptotic collective limit without
reproducing this case's collision rate or finite-particle noise.

### Spatial dimensionality and boundary gates

The small microscopic box is genuinely 3D3V. A planar projection, radial
histogram or spherical average is an observable of that simulation, not a
constraint on the particles' motion. Display a projection and a thin-slice
selection explicitly when those make the dancing particles legible.

A one-dimensional line or two-dimensional sheet of point charges interacting
through \(1/r\) is a different constrained model. A spherical shell or slab
reduction requires derived interactions, phase-space measures and a separate
model-error statement; do not call it the same 3D physical-particle plasma.

Periodic Coulomb needs the declared periodic Green function, neutralizing
background, potential gauge and electrostatic boundary convention. Minimum
image or abruptly cut off \(1/r\) is not that force law. The small periodic
oracle must use a converged/checked Ewald or equivalent sum before comparing
accelerated P3M/FMM actions. Numerical force tolerance is distinct from physical
softening, excluded self force, particle shape and finite-domain model errors.

### Promotion and animation evidence

Before implementing a microscopic case, add its exact interaction definitions,
generated primal/JVP/VJP, initial/boundary data and independent oracle. Pair
force antisymmetry, the energy-gradient identity and Rutherford scattering are
initial gates. A periodic action also needs checked reciprocal/real-space
tails, background and self terms. Pair forces may be tiled; never store the
logical all-pairs matrix. Close encounters, minimum separation and domain or
chart events define the derivative and certification domain explicitly.

A conserved numerical energy or an enclosed instantaneous force is not a
continuous trajectory certificate. Report the actual trajectory/observable
bound status; derive stability amplification and integration/force remainders
before drawing a deterministic trajectory band. Keep statistical estimator
intervals, finite-domain effects and continuum-closure differences separate.

Animations show the same named example at the physical-particle, mean-field
and effective-closure levels when available. Show individual trajectories and
the averaged density/field, with full applicable bands and explicit unavailable
statuses. Record time and memory for solve alone, derivative work, numerical
checking and statistical sampling separately, then report absolute overhead
and total multipliers. Larger counts and accelerated forces are promoted only
on measured accuracy/cost evidence, not hypothetical workstation capacity.

See the primary [Debye-response derivation](https://farside.ph.utexas.edu/teaching/plasma/lectures/node7.html),
which makes the responsive-species convention explicit, and
[LAMMPS's long-range boundary contract](https://docs.lammps.org/kspace_style.html)
for periodic image sums and force-accuracy settings. Pin any external
comparison implementation when its executable case is admitted.

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
