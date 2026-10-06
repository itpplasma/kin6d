# Q1/Q1b: localized two-body quantum Coulomb encounter

Status: implemented candidate with finite-space error bounds, independent
oracles, equation-level derivatives and a classical comparison ensemble
(campaign cases Q1 and Q1b, docs/hierarchy-campaign.md). Continuum
(discretization) error is reported by reference comparison and convergence
indicators, not by a rigorous bound. See also
[quantum-microscopic-validation.md](quantum-microscopic-validation.md) (QG1).

## Model

Two nonrelativistic point particles with masses m1, m2 and bare Coulomb
coupling g = q1 q2/(4 pi eps0) (declared reduced units, hbar declared; g > 0
repulsive, g < 0 attractive). Centre-of-mass and relative motion separate
exactly. The centre of mass is a free Gaussian (analytic spreading). The
relative coordinate is a full 3D Schroedinger problem for the reduced mass mu
with V = g/r. No softening, screening, absorber or radiation; v/c,
relativistic, spin-orbit and radiative corrections are omitted and
unavailable.

Preparation (Q1): the relative packet
psi0 = (2 pi sigma^2)^(-3/4) exp(-|r - r_c|^2/(4 sigma^2)) exp(i k z),
r_c = (b, 0, -z0), i.e. |psi0|^2 has rms sigma per axis, mean wave vector k
along z and impact parameter b along x. Rotational invariance makes this
choice general for isotropic Gaussians. For two equal-mass electrons with equal
rms width sigma_e this is the exact factorization of two independent Gaussian
electrons: sigma = sqrt(2) sigma_e, sigma_cm = sigma_e/sqrt(2) = sigma/2.
Configuration fields: `encounter_config_t` (src/coulomb_encounter.f90).

Spin/exchange: Y_lm(-r) = (-1)^l Y_lm(r). The spatial singlet keeps even l,
the triplet odd l, the product (|up,down>, distinguishable) preparation all l.
Every sector observable is normalized by its own sector norm (reported).
Spin-unresolved detection of the antisymmetrized |up,down> state gives
(|psi(r)|^2 + |psi(-r)|^2)/2: no Mott cross term. Singlet/triplet carry the
cross term 2 Re f(theta) f*(pi-theta) ~ cos(eta ln tan^2(theta/2)).

Scales: lambdabar = 1/k, v = hbar k/mu, Sommerfeld eta = g/(hbar v),
head-on turning radius r0 = 2 eta lambdabar (kin6d_species and the generated
leaves; physics from gen/kin6d_symbolic_physics.f90).

## Representation

psi(r) = sum_{l, m>=0} u_lm(r)/r P_lm(cos theta) c_m cos(m phi): real
harmonics suffice because psi0 is even in y. H_l is independent of m.

Radial FEM-DVR (src/radial_dvr.f90): Gauss-Lobatto-Legendre Lagrange
elements, bridge functions, node r=0 removed (u(0)=0 built in) and node r=R
removed (hard wall). In the first element every retained basis function
vanishes at the origin, so chi_i chi_j/r and chi_i chi_j/r^2 are polynomials
of degree <= 2p-1 and the GLL rule integrates the Coulomb and centrifugal
matrix elements exactly: the singularity is handled by the basis, not
softened. The kinetic matrix is exact on every element. Remaining
approximation: diagonal (lumped) overlap/potential quadrature on smooth
elements, converging spectrally with order (measured below). Repulsive cases
use uniform elements sized for the packet's maximal wavenumber (default 6 GLL
nodes per local wavelength); attractive cases grade geometrically from the
Bohr radius.

Evolution is exact in the represented space: per l one banded symmetric
eigendecomposition (LAPACK dsbevx for small windows, dsbevd otherwise;
LAPACK is FortNum's declared external dependency) reused for all times,
a(t) = V exp(-i E t/hbar) V^T a(0). Only eigenpairs in an adaptive energy
window are kept; the discarded initial weight is measured.

Projection onto Y_lm: Gauss-Legendre in cos theta (l_est + 16 nodes) and
periodic trapezoid in phi (sized by the Bessel content b r/(2 sigma^2)).
l_est = k_hi (|b| + z0 + 8 sigma) + 8, then l_max/m_max are trimmed to the
smallest values whose discarded weight fits tolerance/3.

At matched nodes per wavelength the uniform
Cartesian alternative needs (2 R k_hi ppw/2 pi)^3 points; for the eta scan
below that is 7e7-1.8e9 complex unknowns with a singular point potential to
regularize, versus 1.5e4-1.8e5 partial-wave unknowns with exact singular
integrals (output/coulomb-eta-summary.csv).

| eta | sigma/lambdabar | r0/sigma | R/lambdabar | n_radial | l_max | kept states | setup s | e(t_final) bound |
|---|---|---|---|---|---|---|---|---|
| 0.3 | 10 | 0.06 | 154 | 208 | 73 | 2834 | 0.4 | 2.0e-6 |
| 1 | 10 | 0.2 | 155 | 208 | 73 | 3385 | 0.4 | 1.9e-6 |
| 3 | 12 | 0.5 | 185 | 241 | 86 | 5253 | 0.9 | 2.1e-6 |
| 10 | 40 | 0.5 | 583 | 637 | 278 | 63742 | 52 | 2.1e-6 |

(b = 0, z0 = 6 sigma, t_final = 2 z0/v, gfortran 15 Release, Apple M-series.)

## Observables

- `encounter_angular_bins`: probability that r points into a theta bin (phi
  integrated) with |r| >= r_cut, per sector. Exact Gauss-Legendre bin
  matrices. A finite-time observable: Coulomb trajectories approach their
  asymptotes as O(r0/r), and long-range Coulomb phases mean finite-time bins
  are not asymptotic cross sections.
- `encounter_separation_density`, norm, independent banded-action energy,
  wall-layer probability, point synthesis.
- `encounter_relative_projection`: (x,z) density integrated over y
  (visualization path; bilinear table interpolation, relative amplitude error
  about (k dtab)^2/8, not a bound). `encounter_lab_projection`: one-body
  density of either particle, x1 = X + (m2/M) r, x2 = X - (m1/M) r, by
  convolution with the analytic free CM Gaussian.
- `encounter_classical_ensemble`: independent point pairs sampled from the
  positive initial Wigner function (relative and CM), relative motion by the
  shared nbody/interaction-laboratory step with encounter-resolved step
  control, CM free. Repulsive only. Statistical intervals are binomial.
- Q1b (src/coulomb_kepler.f90): one partial wave, radial Gaussian at rest at
  the outer turning point of nbar; <r>(t), breathing amplitude per period,
  radial density, autocorrelation A(t) and its oracle with every identified
  bound eigenvalue replaced by E_n = -mu g^2/(2 hbar^2 n^2).

## Derivatives

`encounter_bin_jvp` / `encounter_bin_vjp`: sector-normalized bin probability
P = <a|Pi Q|a>/<a|Q|a> w.r.t. (b, sigma, g). Equation level in the full
represented space (complete eigenbasis per l, recomputed, O(n^2) scratch):
d a = U d a0 + dg V[(V^T dH/dg V o F) V^T a0] with the Duhamel
(Daleckii-Krein) kernel F_jk = -(i t/hbar) e^{-i(E_j+E_k)t/2hbar}
sinc((E_j-E_k)t/2hbar). d a0 is the projected FortSym JVP of the packet. The
VJP uses one backward propagation and the transposed contractions. The
discrete mesh, l_max and m_max are frozen design choices (not
differentiated); finite differences must freeze them too (`cfg%hmax`,
`cfg%lmax`, `cfg%mmax`).

## Error ledger (`encounter_ledger_t`)

| Term | Label |
|---|---|
| e(t) <= e0 + t rho/hbar, e0 = ‖a0 - V V^T a0‖, rho = sqrt(sum_lm (sum_j abs(c_j) ‖H v_j - E_j v_j‖)^2) | finite-space bound (inputs in binary64, combination by outward interval leaf) |
| sector bin radius from e(t) | bound for the represented model, same caveat |
| spectral window discarded weight | bound (exact identity, binary64) |
| l/m truncation (discarded, edge weights) | indicator (quadrature-measured) |
| projection defect 1 - ‖a0‖^2, pointwise reconstruction rms | indicator |
| initial packet weight outside R (chi-3 tail) | analytic bound |
| wall-layer probability at output times | indicator for hard-wall influence |
| FEM-DVR vs continuum | reference comparison (hydrogen, Coulomb phases, F_l, free Gaussian) |
| finite time vs asymptotic cross section | model/observable difference, reported, not an error |

The first three formulas require enclosed eigenvectors, residual norms and
orthogonality defects before they establish a deterministic certificate.
Their binary64 inputs are currently unenclosed: report these evaluated
quantities as indicators/finite-space formulas, not `certified-discrete` or
validated Tier A. The benchmark now times embedded eigen-residual and
reconstruction diagnostics separately through the ledger's
`estimator_seconds` and `indicator_seconds`. See
[the measured tier choices](error-control-strategy.md).

## Verification (test/test_coulomb_encounter.f90)

Independent oracles: GLL exactness and Legendre orthonormality; hydrogen
levels n <= 20, l <= 3; Coulomb phases sigma_l = arg Gamma(l+1+i eta) from box
quantization through the A&S 14.5 asymptotic expansion against a binary128
Stirling evaluation (repulsive and attractive); regular Coulomb F_l shape
against a binary128 power series; analytic free 3D Gaussian with b != 0;
exchange parity and the singlet+triplet = spin-unresolved product identity;
late-time Mott ratio and Rutherford fluence (eta = 1, sigma = 10 lambdabar);
eta = 10 quantum vs classical Wigner ensemble; JVP/VJP vs central differences
and the dot-product identity; Q1b autocorrelation vs exact levels, Kepler
period, collapse, half and full revival.

## Evidence

`build/coulomb_encounter_demo [--config FILE] [--no-movie | --movie-only]`
(namelist `&coulomb_encounter`: `etas`, `nclass` (b-scan pairs), `nclass_scan`
(eta-scan pairs, default 1e5), `seed`, `movie_seed`, `movie`) writes to output/
(not committed): coulomb-eta-summary.csv, coulomb-eta-bins.csv (1 degree bins),
coulomb-eta-classical-bins.csv (5 degree bins with binomial standard error),
coulomb-mott-eta-scan.png, coulomb-b-scan.{csv,png}, coulomb-kepler.{csv,png},
coulomb-kepler-summary.csv, coulomb-cost-ledger.csv, the movies
coulomb-encounter-{headon,offaxis}.mp4 with frame PNGs and per-frame CSVs
(coulomb-movie-{headon,offaxis}.csv: time, slow-motion factor, sector norms,
energy, wall probability), coulomb-turning-cut.{csv,png}, plus
`coulomb_encounter_demo.nc` (kin6d_result_io: eta/impact-parameter bin
ratios, per-eta state bounds, Kepler series and errors, bin VJP, cost
timings, literal/resolved configuration and provenance).

Movies (visualization path, no pointwise bound): eta = 3, sigma = 3 lambdabar
= r0/2, z0 = 20 lambdabar, kb = 0 (head-on) and kb = 2. Each frame shows the
slice y = 0 of the relative density |psi(x,0,z)|^2 k^-3 (each sector divided
by its own norm) for singlet, triplet and the spin-unresolved product, beside
the classical Wigner ensemble (8e5 pairs, spin-unresolved, slab |y| < 1/(2k)),
on one fixed four-decade log colour scale. The frame clock runs up to 4x
slower in a Gaussian window around and after closest approach (labelled in
the title). The panels show the classically forbidden hole r < r0; exchange
standing-wave fringes along z (spacing near lambda/2 = pi/k) in singlet and
triplet but not in the spin-unresolved product or the classical ensemble;
the triplet nodal plane z = 0 (psi_T(-r) = -psi_T(r), axially symmetric for
b = 0); weak self-interference arcs of the distinguishable pair after
reflection (t v/z0 ~ 1.2-1.5); the classical X-shaped outgoing caustic and
the Mott node pattern of the outgoing singlet/triplet waves.
coulomb-turning-cut.png quantifies this at closest approach: on the axis the
distinguishable quantum density agrees with the classical ensemble except for
the tunnelling tail inside r0 (no fringes: incoming and reflected parts of
one packet barely overlap for sigma = r0/2); at x = 2 r0 the singlet and
triplet fringes are complementary, their maxima are 3.3-3.5/k apart (the
local wavelength is longer than the asymptotic one near the turning region),
and their average is the fringe-free product, which matches the classical
ensemble within its standard error.

Teaching messages: (1) distinguishable particles follow Rutherford even
quantum mechanically, identical particles interfere (singlet doubles,
triplet vanishes at 90 degrees, oscillating with eta ln tan^2(theta/2));
(2) for eta >> 1 and sigma < r0 quantum bins approach the classical ensemble
except for the exchange term; (3) a Rydberg packet breathes at the Kepler
period, collapses and revives near T_rev/2 and T_rev = (2 nbar/3) T_K, with
the full revival degraded by cubic dispersion (~ (8 pi/3) dn^3/nbar);
(4) when lambdabar is comparable to r0 (eta = 3), exchange produces
lambda/2 fringes in the encounter region that the classical ensemble cannot
show, while the spin-unresolved density stays classical-like.

## Open issues

- GLL rule and normalized associated Legendre recurrence are private here;
  they are generic and belong in FortNum (next action: upstream with tests,
  then pin and delete the private copies).
- No rigorous continuum (FEM-DVR discretization) bound; the finite-space
  bound uses binary64-evaluated residual norms.
- The g-derivative recomputes complete eigenbases per l (O(l_max n^3));
  windowed or Krylov tangents are a measured future optimization.
- Classical comparison is repulsive only (interaction laboratory scope).
- Movie projections are a visualization path without bounds.
