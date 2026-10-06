# Q2: reduced wave-packet models versus the Q1 encounter reference

Status: implemented candidate models with a full-Hamiltonian a-posteriori
estimator against the Q1 partial-wave reference
([coulomb-encounter.md](coulomb-encounter.md)); campaign case Q2
([hierarchy-campaign.md](hierarchy-campaign.md)), gate QG2/QG5 of
[quantum-microscopic-validation.md](quantum-microscopic-validation.md).
These are measurements of model defect in one two-body encounter, not
validated scalable models.

## Models (one assumption changed at a time)

Same physical parameters, preparation, observables and bins as Q1
(`packet_spec_t`, `src/packet_models.f90`).

| ID | State manifold | Approximation added |
|---|---|---|
| M1 | positive Wigner ensemble of classical point pairs (Q1 `encounter_classical_ensemble`) | classical dynamics |
| M2 | product of two isotropic thawed electron Gaussians (WPMD, Heller/Klakow-Toepffer-Reinhard canonical form) | single Gaussian per electron; pair Coulomb seen through erf(q/sqrt(2(s1+s2)))/q; CM width tied to the relative width |
| M2S/M2T | the same pair symmetrized / antisymmetrized (two-determinant TDVP, exact overlap and exchange integrals) | as M2, with exchange in the dynamics |
| M3 | one thawed Gaussian in the relative coordinate under g/r, CM exactly free | single relative Gaussian (no CM tie) |
| M3S/M3T | M3 symmetrized / antisymmetrized | as M3, with exchange |
| M4 (K) | coherent sum of K relative Gaussians, complex amplitudes, all cross terms (variational multi-configuration Gaussians) | finite basis; K = 1 equals M3 |

Canonical TDVP (M2, M3): psi = (2 pi s)^(-3/4) exp(-(1/(4s) - i beta)|r-q|^2
+ i p.(r-q)/hbar + i gamma); with ps = 3 hbar beta the equations are
Hamilton's equations of <H>(q, p, s, ps) and hbar gammadot = p.qdot - s psdot
- <H> (Dirac-Frenkel phase). Energy, forces, JVP/VJP: FortSym leaves from
`gen/gen_packet_models.f90`; the Gaussian-averaged Coulomb definition is
`gaussian_coulomb_potential` in `gen/kin6d_symbolic_physics.f90` (series form
for q^2/(2s) < 0.02).

Metric engine (M2S/M2T/M3S/M3T/M4, cross-check for M2/M3): real parameters
(q, p, alpha, d) per Gaussian, McLachlan form Re(M) xdot = Im(h)/hbar
(equivalent to Dirac-Frenkel on this complex manifold), parity and CM ties
as linear parameter maps, Jacobi equilibration of populated directions,
eigendecomposition with Tikhonov eps (1e-10, 1e-6 for M4). Zero-population
Gaussians (K > 1 starts from one Gaussian, others with d = 0 at fixed
phase-space offsets) follow their own one-Gaussian TDVP in the metric null
space until the variational equations give them weight: basis enrichment,
initial state identical to Q1. Pair integrals (`src/gaussian_pair.f90`):
overlap, moments to fourth order, kinetic, harmonic and Coulomb elements via
the complex Boys function E_n(T) = int_0^1 u^2n e^{T(1-u^2)} du (Gauss-Legendre
for |T| <= 36, erfc endpoint asymptotics beyond). They are hand-written
complex closed forms because FortSym emits no complex arithmetic.

## Estimator and what is rigorous

The model state is projected onto the Q1 basis, a_m(t) = P psi_model(t)
(`encounter_project`, same quadrature as Q1's packet projection), together
with P dpsi/dt. With the Q1 banded Hamiltonian H_h,
r_h = i hbar da_m/dt - H_h a_m. P is linear and time independent, so

  ||a_m(t) - a_ref(t)|| <= ||a_m(0) - a_ref(0)|| + (1/hbar) int_0^t ||r_h|| ds

is a bound in the represented space (binary64; time integral by trapezoid
over the samples, a quadrature of a smooth integrand, not enclosed). The
initial term equals Q1's spectral-window defect (about 1.4e-6). For M2
the CM factor's free-Hamiltonian residual sqrt(6) s |c2| (closed form) is
added in quadrature (Dirac-Frenkel phases split so the cross term vanishes).
Bounded observables follow: |dP_bin| <= 2e + e^2 (sector radius formula).
Not bounded: the Q1 discretization error versus the continuum and the
projection defect of the model state (both reported as indicators). The
estimator includes regularization and integration error automatically.

## Results (b = 0, Q1 eta scan parameters, t = 2 z0/v)

Per cell: state error e (phase-sensitive, the bounded quantity) / phase-free
error / estimator B / angular-bin total variation TV (36 bins, |r| > z0/4).
Sector-normalized for S/T models. Efficiency index I = B/e. Errors are
distances of normalized states (trivially <= 2).

| model | eta = 0.3 | eta = 1 | eta = 3 | eta = 10 |
|---|---|---|---|---|
| M2 | 0.3 / 0.299 / 0.66 / 0.33 | 0.815 / 0.778 / 2.3 / 0.86 | 1.55 / 1.17 / 7.6 / 1.2 | 1.46 / 1.34 / 25.5 / 1.3 |
| M2S | 0.304 / 0.304 / 0.71 / 0.34 | 0.811 / 0.773 / 2.47 / 0.84 | 1.54 / 1.17 / 8.1 / 1.2 | 1.48 / 1.36 / 26.7 / 1.3 |
| M2T | 0.295 / 0.295 / 0.61 / 0.33 | 0.82 / 0.784 / 2.11 / 0.88 | 1.55 / 1.16 / 7 / 1.2 | 1.47 / 1.36 / 25.5 / 1.3 |
| M3 | 0.277 / 0.276 / 0.65 / 0.3 | 0.76 / 0.728 / 2.27 / 0.71 | 1.51 / 1.14 / 7.4 / 0.86 | 1.49 / 1.33 / 25 / 0.96 |
| M3S | 0.282 / 0.281 / 0.7 / 0.3 | 0.756 / 0.723 / 2.45 / 0.7 | 1.51 / 1.13 / 7.9 / 0.9 | 1.49 / 1.33 / 25.6 / 0.96 |
| M3T | 0.272 / 0.271 / 0.6 / 0.29 | 0.764 / 0.734 / 2.08 / 0.72 | 1.51 / 1.14 / 6.8 / 0.88 | 1.49 / 1.33 / 24.4 / 0.97 |
| M4 K=2 | 0.28 / 0.279 / 0.67 / 0.26 | 1.07 / 1.07 / 2.46 / 1.2 | 1.33 / 1.25 / 7 / 0.92 | 1.4 / 1.37 / 24.2 / 1 |
| M4 K=4 | 0.22 / 0.219 / 0.57 / 0.16 | 0.864 / 0.852 / 2.23 / 0.81 | 1.33 / 1.26 / 6.5 / 1.1 | 1.38 / 1.38 / 22 / 1.1 |
| M4 K=8 | 0.245 / 0.244 / 0.92 / 0.12 | 0.585 / 0.582 / 1.81 / 0.22 | 1.29 / 1.16 / 5.5 / 0.78 | 1.4 / 1.37 / 19.1 / 1.1 |
| M1 classical | TV 0.27 | TV 0.19 | TV 0.076 | TV 0.020 (sampling sum-sd 0.024) |

Efficiency index I: 2.1-3.7 at eta <= 1, 4.2-5.3 at eta = 3, 14-18 at
eta = 10 (the estimator then exceeds the trivial bound 2 and is uninformative;
the residual integral grows with the number of phase cycles k z0 while the
actual distance saturates). The bound held at every sample in every run.
Exchange in the dynamics changes nothing measurable here: the parity images
never overlap appreciably (M2S = M2T bins to plotting accuracy), so WPMD
exchange cannot produce the Mott interference visible in the reference.

Impact-parameter scan (eta = 1, `packet-b-scan.csv`): M2/M3 state error
0.82/0.76 at b = 0 falls to 0.44/0.44 at b = 20 lambdabar (2 sigma), where
the packet is mostly deflected as a whole. M4 K = 8 gains at b = 0 (0.59,
TV 0.22) but not at b = 5 and 10 lambdabar (0.75 and 0.71, TV 0.68 and
0.64); at b = 20 lambdabar K = 8 hit the RHS refusal limit (status recorded).

Cost (eta = 1, Apple M-series, Release): Q1 reference setup 0.38 s, 3385
kept states, 6 MB; M3 6 parameters 1e-4 s; M2S/M2T 0.016 s; M4 K = 8 64
parameters 3.2 s; each model's full-H diagnostic (61 projections + DVR
residuals) 0.7-2.8 s. At eta = 10: Q1 setup 53 s, 63742 kept
states, 330 MB; M3 1e-4 s, M4 K = 8 17.5 s; diagnostics 17-38 s per model (`packet-cost-ledger.csv`).

## Derivatives

`packet_bin_jvp`: JVP of a sector-normalized bin probability w.r.t. (b, sigma)
for M3 through the FortSym JVP of the canonical RHS (tangent ODE) and the
exact polarization identity for the quadratic bin functional; basis frozen.
Matches central differences to 5.6e-8 (test).

## Verification (test/test_packet_models.f90)

Complex Boys E_0..2 vs binary128 quadrature (63 points, |T| to 900, all
phases) 3.6e-15; overlap/moment, kinetic (gradient form) and Coulomb pair
elements vs direct 3D quadrature 1e-13; free TDVP width vs analytic
spreading 1.5e-10 and free Gaussian vs the Q1 g = 0 reference 6e-8;
harmonic oscillator (TDVP exact for quadratic H): centroid and width of M3,
the metric engine, M3S and M3T vs analytic motion 6e-11/3e-10; canonical vs
metric engine (M3; M2 with CM tie) 7e-10/4e-7; energy and norm conservation
(M4 K = 4, M2S, M2T, M2) 1e-9/1e-8 at eps = 1e-10; preparation equals Q1;
Duhamel bound >= actual error at every sample for six models; enrichment
K = 1, 2, 4 reduces the error (eta = 1, sigma = 3 lambdabar); M3 bin JVP vs
central differences.

## Physics

- Single Gaussians (M2, M3) carry the centroid on a mean-field trajectory
  with erf-smeared Coulomb: an effective smearing over the packet width, a
  model error and not numerics. At b = 0 the centroid stays on the axis and
  the packet cannot split, so all large-angle Rutherford/Mott weight is
  missing; the error does not decrease with eta in this scan because sigma
  grows with eta (r0/sigma <= 0.5).
- Width dynamics: the free width grows without bound (spreading) and the
  Coulomb width force only adds to it; M2's tie lets the Coulomb force widen
  the CM factor, so its two-body error exceeds the relative error.
- Exchange in the dynamics (M2S/M2T, M3S/M3T) changes the sector state only
  where the two parity images overlap; it cannot produce the Mott
  interference structure without splitting.
- K-Gaussian enrichment reduces the angular-bin error most (K = 8) but is
  not monotone in K with fixed guide offsets and eps = 1e-6 regularization.

## Evidence

`build/packet_models_demo [--config FILE]` (namelist `&packet_models`:
samples, nclass, movie_frames, etas, b_values, movie, b_scan) writes to
output/: packet-eta-errors.csv, packet-eta-bins.csv, packet-cost-ledger.csv,
packet-b-scan.csv, packet-movie-separation.csv, packet-error-vs-eta.png,
packet-angular-eta1-eta10.png, packet-density-{mid,final}.png,
packet-density.mp4, and packet_models_demo.nc (result_io provenance).

## Open issues

- Enrichment is not monotone in K with fixed guide offsets: K = 2 and 4 are
  worse than K = 1 at eta = 1 (sigma = 10 lambdabar) because the eps = 1e-6
  Tikhonov regularization (needed to keep small-amplitude vMCG directions
  non-stiff) is not variational; eps = 1e-8 makes runs exceed minutes.
  Next action: adaptive spawning where the projected residual is largest,
  with a measured transfer defect, instead of fixed guide offsets.
- No model here is acceptable for eta >= 3 at r0/sigma = 0.5, b = 0: the
  scattered wave is a broad angular shell that isotropic Gaussians on the
  axis cannot represent. M1 is the better reduced model there (TV 0.02 at
  eta = 10). Off-axis rings of Gaussians need m > 0 reference content.
- Estimator: the norm-residual Duhamel bound is 2-4 times the actual error
  at eta <= 1 but uninformative (> 2) at eta >= 3; a goal-oriented
  (bin-projected) estimator is the next step. The time integral is a
  trapezoid of samples, not an enclosure.
- Complex Gaussian integrals are hand-written (FortSym has no complex
  codegen); M2/M3 leaves are generated. Complex support in FortSym would
  remove the exception.
- The continuum statement still lacks a bound on the Q1 discretization error.
