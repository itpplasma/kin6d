# Microscopic quantum experiments and validation

Status: current plan, 2026-10-06. Audience: lecture students and implementers
who need to see microscopic quantum physics and measure the approximations
used to reach larger systems. The campaign log in
[hierarchy-campaign.md](hierarchy-campaign.md) is the authoritative task
state; this document records the scientific contract, what each implemented
case claims, the measurements that changed the plan and the remaining order.
The classical counterpart is
[microscopic-many-body.md](microscopic-many-body.md#microscopic-teaching-and-verification-cases).

## Physical parent and comparisons

The parent is the antisymmetric finite-N Schroedinger equation or density
operator for spin-1/2 electrons with unscreened Coulomb interactions. Isolated
encounters serve scattering; neutralized periodic 3D cells serve electron
gases. Every case declares mass, charge, hbar, preparation, spin sector,
geometry, background and boundary conditions, all taken from the shared
species/unit record and the single symbolic Coulomb source
(`gen/kin6d_symbolic_physics.f90`). Nuclear motion is classical or prescribed
when present, with that approximation stated. Nonrelativistic electrostatics
is the controlled domain; v/c, x_F and the omitted electromagnetic physics are
recorded, and QED is outside the native scope. No screened potential is
inserted to demonstrate screening.

A microscopic quantum reference evolves the interacting state. Hartree, HF,
wave-packet, Wigner and fluid closures are separately identified
approximations. A classical finite-N ensemble is a different physical model:
its difference from quantum mechanics is the phenomenon being shown, not the
quantum solver's numerical error. Classical and quantum runs use the same
physical parameters and disturbance. For unexchanged Gaussian preparations
the positive initial Wigner law supplies the classical ensemble; for
antisymmetrized or degenerate thermal states the preparation difference is
recorded rather than hidden.

Bulk degeneracy and encounter diffraction answer different questions.
Degeneracy compares k_B T with E_F, or the density with the thermal
wavelength. Diffraction in one encounter compares the packet scale with the
spatial variation of the Coulomb field, measured by the Sommerfeld parameter
eta = g/(hbar v). High temperature makes encounters more quantum, not less:
at 10 keV the reduced thermal wavelength hbar/sqrt(m_e k_B T) is about
2.76 pm and the Landau length e^2/(4 pi eps0 k_B T) about 0.144 pm, while a
tokamak plasma is far from bulk degenerate. Other wavelength conventions
change these factors. A discussion of this short-range regime is
[Alharbi et al., PTEP 2026, 023J01](https://doi.org/10.1093/ptep/ptaf188);
its mean-field model is not a certificate for the finite-N cases here.

## Status vocabulary

Every reported output carries one of these labels, following AGENTS.md:

- *certified-discrete*: a rigorous enclosure for the declared finite basis or
  finite model, including arithmetic; continuum projection error is
  unavailable unless stated;
- *certified*: an enclosure of the named continuous model (for example the
  ideal Fermi EOS integrals and root);
- *reference comparison*: agreement with an independent exact result,
  literature value or richer model, with the measured difference;
- *indicator*: resolution, cutoff or estimator differences that are not a
  remainder bound;
- *statistical*: a sampling interval with a declared law (CLT, batch means,
  jackknife); never a deterministic certificate;
- *candidate*: a fast result without any of the above.

## Implemented cases

Gate IDs QG0--QG8 are the original plan names; Q/C/H IDs are the campaign case
names. Each row names its own document, which carries the full ledger.

| Gate | Case | Document | Claim and status |
|---|---|---|---|
| QG0 | Q0 periodic opposite-spin pair quench, 125-mode exact sector | [quantum-reference-gallery.md](quantum-reference-gallery.md), [classical-pair.md](classical-pair.md) | Certified-discrete continuous-time Hilbert residual and slab-bin bounds (bin radius 2.25e-3 at 1024 steps); 125 vs 343 modes differ by 5.2e-3, an indicator; classical pairs statistical |
| QG0 | Finite-basis evolution certificate (FortNum `fortnum_ode_residual_certificate`) | [quantum-evolution.md](quantum-evolution.md) | Certified-discrete for a supplied Hermitian matrix |
| QG1 | Q1 two-body encounter in the relative coordinate, partial waves on a GLL FEM-DVR, product/singlet/triplet and classical ensemble; Q1b Rydberg packet | [coulomb-encounter.md](coulomb-encounter.md) | Candidate with finite-space bounds; continuum error by reference comparison (hydrogen levels 9e-13, Coulomb phases 1e-9, Mott ratio 0.018, Rutherford fluence 1%, Kepler autocorrelation 1e-10) and indicators; bin JVP/VJP vs FD 1.6e-8 |
| QG2/QG5 | Q2 Gaussian packet models M2--M4 (WPMD, two-determinant, relative TDVP, vMCG K = 2/4/8) vs Q1 | [packet-models.md](packet-models.md) | Candidate models with a full-Hamiltonian Duhamel estimator measured against the Q1 discrete reference (held at every sample, vacuous above 2 at eta = 10) |
| QG3 | Q3 finite-N ideal fermions: canonical sums, exchange hole, density matrix, box EOS | [fermion-ensemble.md](fermion-ensemble.md) | Interval enclosures with a rigorous omitted-orbital remainder; thermodynamic-limit and classical limits as reference comparison |
| QG4 | Q4 interacting N = 2--4 jellium cell, sector Lanczos, 1/2/3-RDMs, weak-kick evolution | [jellium-cell.md](jellium-cell.md) | Finite-basis reference; kick band 2e-6 at t = 8 certified-discrete; cutoff differences indicators; Theta = 0.5 rows for N = 3, 4 unconverged indicators |
| QG5 | Q5 Hartree, RHF/UHF, Mermin, static TDHF/TDH in the Q4 basis | [mean-field-cell.md](mean-field-cell.md) | Reference comparison against Q4; full-parent residual with a Temple bound for N = 2 (gap input is a Lanczos candidate) |
| QG6 | Hybrid classical test-electron plasma with promotion of close quantum encounters to an l-local ring kernel | [hybrid-encounter.md](hybrid-encounter.md) | Candidate kernels with quadrature/truncation indicators; Monte Carlo rates statistical; explicit interface-defect entry |
| QG7 | Q7 Debye -> finite-T Fermi/Thomas--Fermi -> Lindhard static screening; 1X1V Wigner--Poisson -> Vlasov | [screening-wigner.md](screening-wigner.md) | Candidate with rigorous truncation tails and quadrature indicators; dispersion roots reference comparison |
| QG8 | Q8 relativistic ideal Fermi EOS; nonrelativistic finite-T Fermi gas | [fermi-gas.md](fermi-gas.md) | Certified enclosures (relativistic to x_F = 10 at 1e7 K; nonrelativistic to theta = 1e-4) |
| H1b | Quantum-classical matched finite-N correlation map vs (Gamma, Theta) | [hierarchy-campaign.md](hierarchy-campaign.md) | In progress on a worker branch; no result claimed |

All cases share the species record, Coulomb leaves, result I/O with embedded
configuration and provenance ([interface-reproducibility.md](interface-reproducibility.md))
and the figure style ([figure-style.md](figure-style.md)). Movies and figures
are listed in [chapter-gallery.md](chapter-gallery.md).

## Findings that change the plan

**Visible encounter quantum structure comes from exchange.** In the Q1
eta = 3, sigma = r0/2 movies the singlet and triplet relative densities show
standing-wave fringes in the encounter region (maxima 3.3--3.5/k apart, the
local wavelength being longer than lambda/2 = pi/k near the turning point)
and the triplet nodal plane z = 0. The distinguishable (spin-unresolved)
density shows no such fringes: on the axis it agrees with the classical
ensemble except for a tunnelling tail inside r0, and the singlet/triplet
average matches the classical ensemble within its standard error. Outgoing
singlet/triplet bins follow the plane-wave Mott oscillation
cos(eta ln tan^2(theta/2)) up to eta = 10, and at eta = 10 quantum and
classical distinguishable bins agree within 4 standard deviations. A quantum
correction to a classical encounter picture is therefore first an exchange
correction, then diffraction at eta <~ 1.

**Gaussian packet models fail for eta >= 1, and classical ensembles beat
them for eta >= 3.** At b = 0 and r0/sigma = 0.5 the state error of the best
model, M4 with K = 8 Gaussians, is 0.24/0.59/1.29/1.40 at
eta = 0.3/1/3/10 (distance of normalized states, at most 2); its angular-bin
total variation is 0.12/0.22/0.78/1.13, against 0.27/0.19/0.076/0.020 for the
classical Wigner ensemble M1. Exchange inside WPMD dynamics has no
measurable effect because the parity images never overlap, so these models
cannot produce the Mott interference. Scalable packet methods are admitted
only in the regime where they were validated (eta < 1 here); for eta >= 3 a
classical trajectory with a quantum encounter kernel is cheaper and more
accurate.

**The quantum Coulomb logarithm correction is universal.** Quantum minus
classical ln Lambda equals ln eta - Re psi(1 + i eta) for both a sharp cut
and Yukawa screening; the Yukawa-minus-cut offset ln 2 - 1/2 - gamma is the
same classically and quantum mechanically. The l-local ring kernel
chi_l = 2 atan(eta/(l+1)) is exact for momentum transfer and has a 10.6%
second-moment defect at eta << 1; unpolarized exchange adds -ln 2/2 at
eta << 1 (kernel only). In relaxation rates for T = 1e-3..100 E_h the
classical model errs by up to 15%, the hybrid interface defect stays below
3.5e-3, and the hybrid is 1.5 times cheaper at matched 1% statistics for
eta_u < 0.3. Plain Maxwellian velocity draws have infinite variance and
naive velocity rotation at the detection sphere double counts (+11%); both
are now excluded by construction.

**Open-shell Hartree--Fock breaks translation symmetry.** For N = 3, 4 in the
jellium cell the UHF minima break translation symmetry and are non-aufbau;
the plane-wave determinant is a saddle, so the Temple bound does not apply.
For closed-shell N = 2 the Temple bound overestimates the true HF energy
error by 1.25--2.2. Static TDHF tracks the exact N = 2 response (chi/chi0 at
r_s = 8: exact 0.69, TDHF 0.66, RPA 0.49), while RPA overscreens. Mermin
thermal energies differ from the exact canonical cell by a factor 2.7 at
N = 2, r_s = 0.5, dominated by the grand-canonical orbital ensemble. A
mean-field reference must therefore report its symmetry and ensemble, and
finite-cell comparisons need a canonical or thermodynamic-limit reference.

**Three-body cumulants are not small in the quantum cell.** The ratio
|lambda3|/|Gamma3| of the connected three-body RDM grows from 0.060
(r_s = 0.5) to 0.236 (r_s = 8) in the N = 3 ground multiplet; for N = 4 it is
0.29--0.41, dominated by multiplet mixing. Exact static response relative to
the same-cell Lindhard function falls from 0.97 (r_s = 0.5) to 0.69
(r_s = 8). A quantum BBGKY closure at the pair level is therefore an
approximation to be measured, not assumed.

**Smaller findings.** The alternating power-sum recursion for canonical
fermion sums cancels by up to 5e36 when cold; the subtraction-free recursion
is the production path. Wigner--Poisson converges to Vlasov as O(hbar^2)
(observed rate 1.99), with Wigner negativity 0.26 at t = 20 in the two-stream
run. The finite-size energy of ideal fermions behaves as e_N - e_TL ~ 1/N.

## Remaining plan

The direction is the quantum BBGKY/cumulant hierarchy and adaptive
correlation-order and classical/quantum switching
([adaptive-multiscale-uq.md](adaptive-multiscale-uq.md)). In order:

1. **H1b matched correlation map.** Finish the quantum-classical finite-N map
   (pair and three-body cumulants vs Gamma and Theta from Q4 and the
   classical H1a sampler); it is the reference for steps 2--3.
2. **Closures beyond the first YBG equation.** Self-consistent BGY--Kirkwood
   and the second YBG equation for g3, measured against the H1a Monte Carlo
   g2/g3; c3 resolution near contact.
3. **Quantum BBGKY with measured cumulants.** The first RDM equation with the
   Q4 two-RDM, then pair-level cumulant truncations, compared against the Q4
   kick dynamics and the H1b map. Report the omitted lambda3 effect as a
   model-defect indicator.
4. **Thermodynamic-limit references.** Finite-size and ensemble corrections
   of Q4/Q5 against PIMC/literature data for the uniform electron gas and the
   OCP; canonical versus grand-canonical comparisons stated explicitly.
5. **Lindhard with local-field corrections.** Static LFC and an interacting
   RPA/HF comparison against the Q4 response, with a rigorous Lindhard
   quadrature enclosure where affordable.
6. **Collision kernels and larger Lambda scans.** Exchange in the hybrid
   Monte Carlo; microscopic collision coefficients (classical side, C3/C4 at
   Lambda = 100..10000) compared with the QG6 kernels.
7. **Scalable packets only where validated.** Packet or orbital methods for
   100--10000 electrons only in regimes where Q2-type comparisons pass for the
   named observable; otherwise classical trajectories with quantum kernels.
8. **Adaptive switching.** Correlation-order and quantum promotion
   indicators (eta, Theta, c3 size, Kirkwood defect, lambda3 ratio) with
   lift/projection and interface-defect ledgers, first in the QG6 setting.

Generic pieces found on the way go upstream: Gauss--Lobatto and Legendre
rules are in FortNum; the Z function, digamma/Riccati--Bessel, graded Gauss
rules and series statistics are FortNum candidates; FortSym complex code
generation is needed for packet models.

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
probability bins provide bounded observables. Q2 shows the limit of this
route: the Duhamel integral grows with the number of phase cycles while the
true distance saturates, so its efficiency index rises from 2--4 at
eta <= 1 to 14--18 at eta = 10, where it exceeds the trivial bound.
Goal-specific estimators are needed there.

An auxiliary DFT or HF residual certifies its own equation, not the
interacting many-electron parent; Q5 evaluates the HF determinant in the
parent basis for that reason. Thermal-state errors use spectral/partition
bounds or density-operator estimates, not a pure-state formula by analogy.

Small-system calibration cannot prove a large-system error bound; two
resolutions and agreement between closures do not prove a remainder. Publish
missing terms and failed certificates rather than displaying a fabricated
complete band. Production promotion requires the declared numerical
certificate and model applicability gate to close.

## Figures and animations

Each artifact has a one-sentence teaching message, reproducible data,
physical units and normalization; the current list is in
[chapter-gallery.md](chapter-gallery.md) and the style in
[figure-style.md](figure-style.md). Spatial projections integrate the omitted
coordinate; a slice states its own normalization (Q1 normalizes each spin
sector by its norm on the y = 0 slice). Quantum density is probability, not
trajectories. Classical plots state ensemble size and each displayed
realization's role. Conditional densities specify their conditioning bin and
normalization.

Comparisons use equal scales and a shared clock; slow-motion windows are
labelled in the frame. Numerical intervals are shown in full, and model and
statistical bands are styled separately. Probability heatmaps need companion
bin or projection bounds; colors alone are not a bound. First, encounter and
outgoing frames, grayscale and labels are inspected before sharing. GIF,
MP4, snapshots, CSV and NetCDF results stay outside Git.

## Execution rules

This gallery does not reorder the direct-steady scientific spine. One
complete runnable increment at a time, in Fortran, with retained FortSym
generators for owned leaves, primal/JVP/VJP and enclosures. Parameter
sensitivities cover impact parameter, packet width, coupling, density or
temperature where smooth. Spin choices, basis changes and promotion events
are discrete; their derivatives are declared or marked unavailable.

At each gate record independent oracles, source and dependency revisions and
the actual numerical status. Review inputs are frozen by base commit and
patch digest. Release and checked Debug CTest evidence is retained. Cost
reports separate solve-only time, error-estimation time, derivatives,
preparation/reuse, sampling and rendering, and compare models at matched
observable accuracy. No published particle count becomes a KIN6D
performance claim.

## Primary references

- [Svensson et al.](https://arxiv.org/abs/2211.08560): restricted Gaussian
  trajectories, approximate Pauli interactions and near-linear scaling.
- [Angermeier and White](https://arxiv.org/abs/2104.03441): distinct
  packet-shape, exchange and correlation approximation errors.
- [Chabrier and Potekhin](https://arxiv.org/abs/physics/9807042):
  electron-ion thermodynamics from classical to degenerate regimes.
- [White](https://arxiv.org/abs/2410.23599): dynamical response limitations
  of an orbital-free density approximation.
- [Babbush et al., Phys. Rev. X 8, 011044](https://doi.org/10.1103/PhysRevX.8.011044):
  plane-wave electronic Hamiltonians.

Chris&AI
