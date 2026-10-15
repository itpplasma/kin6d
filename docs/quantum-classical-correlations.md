# Quantum versus classical finite-N correlations (H1b)

Status: finite-N reference map, 2026-10-06. Hierarchy-campaign item H1, the
quantum half of the correlation map. Planned items are not support claims.

## Question and setting

How large are diffraction, exchange and degeneracy effects on pair and
three-body correlations, measured against the classical parent with the
same N, periodic cube, uniform neutralizing background, density and
temperature? Electrons only, atomic units, `r_s`, `Theta = kT/E_F` and
`Gamma = g/(a kT)` from `kin6d_species` (`a = r_s`). On the grid
`Gamma = r_s/(1.8416 Theta)` and the thermal wavelength
`lambda_T = (2 pi hbar^2/(m kT))^(1/2)` gives `lambda_T/a = 1.847 Theta^(-1/2)`,
a function of `Theta` alone.

## Models (one cell, two mechanics)

- Classical parent (`src/classical/classical_cell.f90`): N point electrons with the
  shared Ewald energy (`kin6d_nbody`/`kin6d_ewald`). The energy is pairwise
  additive with the zero-mean Ewald pair potential `u` (Fourier `Q = 0`
  removed, as in the quantum Hamiltonian). Spin does not enter the classical
  Boltzmann weight; spin labels are spectators with probability 1/2 that a
  pair has equal spins.
  - N = 2: the centre of mass integrates out. `p(d) = exp(-beta u(d))/Z` on
    the torus is integrated by the periodic midpoint rule (Z, Fourier
    moments) and by spherical Gauss rules (shell bins). This is deterministic.
  - N = 3, 4: canonical Metropolis Monte Carlo with single-particle moves
    by the one periodic-cell sampler shared with H1a
    (`src/classical/ocp_monte_carlo.f90`, `mc_setup_cell`: incremental Ewald with a
    real image cube chosen by the `kin6d_ewald` cost model). The Kelbg and
    Kelbg+Pauli laws are a declared short-range addition
    (`src/pairs/pair_law.f90`) on top of the one Ewald Coulomb energy.
    Intervals are batch-means statistical errors (40 batches,
    `kin6d_series_statistics`). They are labelled statistical, never
    certified. The N = 2 quadrature uses the independent two-body
    `kin6d_nbody` pair potential and is the sampler's oracle.
- Exact quantum N = 2 (`src/quantum/quantum_pair_cell.f90`): the centre of mass and
  relative coordinate separate exactly on the torus. The relative wave
  function is periodic or antiperiodic per axis according to the parity of
  the centre-of-mass momentum component (twist classes t = 0..3). Even and
  odd relative parity correspond to the singlet and the triplet. All
  eigenstates of the relative plane-wave basis `|kv|^2 <= kmax2` (units
  `pi/L`) are kept. `kmax2` makes the cutoff kinetic energy at least
  `tau kT` at the largest Theta (tau = 36). The cutoff indicator repeats the
  calculation with tau = 25.
- Q4 cell for N = 3, 4 (`kin6d_jellium_cell`): every `(S_z, K)` sector is
  diagonalized densely (`jellium_full_spectrum`), so the canonical ensemble
  has no thermal truncation. The plane-wave cutoff (`|n|^2 <= c`) stays the
  limiting approximation.
- Ideal fermions: Q3 `kin6d_fermion_ensemble` with fixed-N_sigma states,
  combined over all `S_z` with weights `Z_{N_up} Z_{N-N_up}`. The same Q4
  cell at coupling 0 serves as the ideal reference for three-body moments.

## Common estimators

Both mechanics are measured by the same diagonal (position-space)
quantities. They are phase-space free, so the classical and quantum values
mean the same thing.

- Spin-summed shell-averaged pair function on equal shells inside
  `|r| < L/2`: `g_b = V/(N(N-1)) <sum_{i/=j} 1_b(r_ij)>/V_b`. Quantum spin
  channels `g_same`, `g_opp` use the Q4 normalization `2 (N/2)^2`. With
  spin-blind labels the classical channels are `g (N-1)/N`.
- Diagonal Fourier moments on `S = {Q = 2 pi n/L : 1 <= |n|^2 <= 3}`
  (26 vectors): `f2(Q) = <sum_{i/=j} e^{iQ.(x_i-x_j)}>/(N(N-1))` and
  `f3(Q1,Q2) = <sum_{i,j,k distinct} e^{i(Q1.x_i+Q2.x_j+Q3.x_k)}>/(N(N-1)(N-2))`
  with `Q3 = -Q1-Q2` in S (all ordered triples). For a homogeneous state at
  nonzero Q these are the Fourier coefficients of the Ursell (connected)
  pair and triplet cumulants of the diagonal densities. The norms are
  `c2 = rms_S f2` and `c3 = rms f3`. Classical Monte Carlo norms carry a
  jackknife bias correction. Quantum `f2`, `f3` come from the ensemble
  2-RDM and the per-state 3-RDM entries. The quantum RDM cumulant `lambda3`
  (Q4) removes the ideal exchange part and is a different quantity.
- Difference norm: `||Delta|| = (sum_b V_b Delta_b^2 / sum_b V_b)^(1/2)`.

## Decomposition

- total `Delta = g_q - g_cl`;
- (a) ideal part `g_{q,id} - 1`: exchange and degeneracy without
  interaction (the classical ideal finite-N gas has `g = 1` exactly);
- (b) interaction-induced part `(g_q - g_{q,id}) - (g_cl - 1)`: diffraction
  of the Coulomb hole plus interaction-exchange cross terms;
- (c) channels `g_same,q - g_cl/2` (exchange and diffraction) and
  `g_opp,q - g_cl/2` (diffraction only, N = 2 normalization).

## Semiclassical predictors (reduced models, model-defect statements)

Generated leaves are in `gen/gen_quantum_classical.f90`.

- Wigner--Kirkwood, leading order in hbar^2. Expanding the N-body Slater sum
  gives the configuration weight
  `w = 1 - (hbar beta)^2/(12 m) sum_i lap_i U + hbar^2 beta^3/(24 m) sum_i |grad_i U|^2`.
  This is checked against the exact harmonic-oscillator density: the error is
  O(beta^4) with a measured ratio of 16.0. For N = 2 it reduces to the relative
  coordinate with reduced mass `m/2`. With the Ewald potential,
  `lap_i U = 4 pi g (N-1)/V` away from contact, which is constant and cancels.
  The delta-function contact term multiplies `exp(-beta u) = 0`. The correction
  is therefore carried by the force term alone.
- Leading exchange: identical spin-1/2 pairs have equal spins with
  probability 1/2. The ideal same-spin factor is
  `1 - X(r)`, `X = sum_images exp(-|r + nL|^2/lambda^2)`,
  `lambda = hbar/(m kT)^(1/2)`. The spin-averaged weight is therefore
  `1 - sum_pairs X/2`.
- Predictor `g_pred` = classical ensemble reweighted by `w * (1 - sum X/2)`
  (quadrature for N = 2, Monte Carlo for N = 3, 4). The efficiency index is
  `||g_pred - g_cl|| / ||g_q - g_cl||` and the residual is
  `||g_pred - g_q|| / ||g_q - g_cl||`.
- Kelbg potential
  `u_K = g [1 - exp(-x^2) + sqrt(pi) x erfc(x)]/r`, `x = r/lambda`. It
  replaces the minimum-image Coulomb term, and every other image adds the
  Gaussian-decaying excess. The Kelbg+Pauli model adds
  `-kT ln(1 - X/2)` per pair. The model error is
  `||g_model - g_q|| / ||g_cl - g_q||` (1 means no better than classical).

## Verification (`quantum_classical_h1b_oracles`)

| Oracle | Result |
|---|---|
| Ewald energy pairwise additive (N = 3, different alpha/cutoffs) | 1.7e-15 |
| N = 2 Monte Carlo vs quadrature of `exp(-beta u)` (g, f2, dg/dbeta, reweighted g) | within 4.5 s.e. |
| Classical g depends on Gamma only (r_s = 4, Theta = 2 vs r_s = 1, Theta = 0.5) | < 1e-6 |
| Relative reduction in the Q4 pair set = Q4 S_z = 0 two-electron sector spectrum (K = 000, 110) | 1.4e-14 |
| Ideal relative reduction = Q3 canonical ideal fermions over all S_z | 3e-15 |
| WK weight vs exact harmonic density, error ratio for beta omega = 0.2/0.1 | 15.99 (O(beta^4)) |
| Ewald force and Laplacian (`4 pi g/V` off contact) vs central differences | 4e-7 rel, 1e-3 |
| Kelbg: Coulomb at 8 lambda, contact `sqrt(pi) g/lambda`, series switch, JVP vs FD | 1e-14 to 1e-8 |
| N = 2 spectral dg/dTheta vs central differences | 2.5e-8 rel |
| Q4 diagonal f2/f3 of a Slater determinant vs brute-force Wick sum | exact |
| Sparse 2-RDM shell bins (`fock_rdm2_entries`) vs dense `jellium_pair_bins` | 4e-15 |
| Dense full-sector spectrum (364 states) vs block-Lanczos ground state | 1e-10 rel |
| Shared sampler energy minus offset = independent two-body pair sum, N = 3 (Coulomb, Kelbg+Pauli; minimum image and image cube) | 1.3e-13 |
| OCP form of the shared sampler (a = 1, minimum image), N = 2, vs pair quadrature, Gamma = 2.17 | 2.0 s.e. |

## Results (`correlation_map_quantum_demo`, defaults)

Grid r_s = 0.5, 1, 2, 4, 8 and Theta = 0.5, 1, 2, 4, 10 (Gamma 0.03 to 8.7),
12 shells. The demo writes `output/qcmap-{n2,n2-g,many,derivative,cost}.csv`,
`output/correlation_map_quantum_demo.nc` and the figures
`output/qcmap-{pair-g,regime-pair,regime-three-body,wigner-kirkwood,kelbg}.png`.
A full run takes 40 min on one core; N = 4 accounts for 2/3 of that.

N = 2 (deterministic). The relative-basis cutoff indicator is at most
3.4e-4, below 2.5 % of the smallest difference on the grid.

- Total `||g_q - g_cl||` ranges from 0.47 (r_s = 8, Theta = 0.5) to 0.006
  (r_s = 4, Theta = 10). It falls monotonically with Theta at fixed r_s.
  Outside `r > lambda_T` the difference is at most 0.013 (zero bins exist
  for Theta <= 2, where `lambda_T > L/2`). The classical limit is recovered
  except inside the thermal wavelength. The quantum contact value
  `g_q(0)` stays at 0.02 to 0.74, while the classical value is
  `exp(-beta u)`, effectively 0.
- Ideal part: 0.02 (Theta = 0.5, singlet-dominated N = 2 ground state) and
  0.11 at Theta = 1, then 0.037 at Theta = 10. It depends on Theta only.
- The interaction part exceeds the total wherever the exchange hole and the
  Coulomb hole overlap. The interaction-exchange cross term is as large as
  either part. The opposite-spin (pure diffraction) channel dominates for
  Theta <= 1 (0.15 to 0.52). At Theta >= 4 the same-spin channel carries
  most of the difference at small r_s.
- Theta derivative at r_s = 2, Theta = 2: `d||Delta||/dTheta = -1.9053e-2`
  from the exact thermal covariances (quantum spectral, classical
  quadrature). Central differences (h = 2e-3) agree to 2.5e-6 relative.
- Wigner--Kirkwood + exchange: the residual is 0.14 at r_s = 0.5,
  Theta = 10 (Gamma = 0.03, `lambda_T/a = 0.58`) and 0.16 to 0.27 for
  r_s = 8 at Theta = 2 to 4. It is 0.4 to 1.2 for `lambda_T/a >= 1.3` at
  small r_s and fails (residual 0.7 to 8) for Theta = 0.5. The first-order
  reweighting overshoots near contact, where `hbar^2 beta^3 |grad u|^2` is
  not small. Useful regime: `lambda_T/a <~ 0.6`, or strong repulsion with
  `lambda_T` below the Coulomb-hole size (r_s = 8, Theta >= 2).
- Kelbg model error (ratio to the classical error). Kelbg alone gives
  0.016 to 0.09 at Theta = 0.5 for r_s <= 2 and 0.07 to 0.5 for r_s >= 4 up
  to Theta = 2. It is 1.2 to 2.3 for r_s <= 2 at Theta >= 2 because it
  omits exchange, which at high T is comparable to diffraction.
  Kelbg+Pauli reaches 0.03 to 0.08 at Theta = 10 for r_s <= 1. It is no
  better than classical (0.8 to 1.2) for r_s >= 4 at Theta >= 2, and gives
  0.39 to 0.61 at Theta = 0.5. There the spin-averaged Pauli term is wrong
  because the singlet dominates.

N = 3, 4 (Q4 cutoff 3 / 2; Monte Carlo 40000 sweeps). Q4 thermal rows are
computed only where the ideal canonical occupation outside the cutoff is
at most 5 %: Theta <= 1, plus the ground multiplet. Theta = 2 misses 7 %
(N = 3) and 22 % (N = 4) and is skipped. The cutoff indicator (c vs c-1)
is calibrated on N = 2 against the exact cell. There the Q4 c = 3 error is
0.005 to 0.024 and the c = 2 to 3 change 0.002 to 0.008, so the indicator
underestimates the error by up to a factor of 3.3.

- Pair: total 0.05 to 0.33 (N = 3), 0.04 to 0.30 (N = 4). Ideal part
  0.12 to 0.13, interaction part 0.10 to 0.44. N = 3 rows pass the reference
  test (indicator, ideal truncation, statistics below 0.2 of the
  difference) except r_s <= 1 at Theta = 1. No N = 4 row passes: the c = 2
  indicator is 0.015 to 0.13.
- Three-body connected moment `c3` at Theta = 0.5, N = 3: quantum
  0.0020/0.0036/0.0066/0.012/0.020 against classical
  0.0034/0.0082/0.018/0.036/0.064 (r_s = 0.5..8; statistical s.e.
  <= 4e-4). Quantum delocalization suppresses connected triplet
  correlations by a factor 1.7 to 3.2. The ideal exchange part is 4.5e-4
  (N = 3) and 1.3e-3 (N = 4). Ground-multiplet values: 0.0015 to 0.020
  (N = 3).
- Finite size: classical `c3(N=4)/c3(N=3)` = 0.78 to 0.70 at Theta = 0.5,
  and the quantum ratio is 0.75 to 1.26 with the unconverged N = 4 basis.

What tiny N can and cannot say. The cell holds every particle within about
one Wigner-Seitz radius of all the others, and S probes only the three
longest wavelengths (`|Q| a` of about 3 to 6). Classical and quantum values
share these finite-size constraints, so the ratio and difference at fixed N
are meaningful. They measure how strongly diffraction, exchange and
degeneracy reshape the correlations under identical constraints. The
absolute `c3`, its N dependence, and any Kirkwood or closure statement are
not thermodynamic-limit quantities. Two values of N cannot be extrapolated:
the 20 to 30 % N = 3/4 change is the finite-size indicator, not a
correction.

## Status and open items

Every quantum number above is a finite-basis candidate. N = 2 is close to
the basis limit (cutoff indicator <= 3.4e-4); N = 3, 4 rows are indicators
unless flagged `reference = 1`. Classical Monte Carlo intervals are
statistical, and the N = 2 classical rows are deterministic quadrature. No
row carries a certificate.

- No rigorous plane-wave remainder for the Coulomb cusp. The N = 2
  calibration shows that the Cauchy indicator underestimates the error.
- N = 4 thermal moments cost about 3 ms per state (sparse 3-RDM); Theta > 1
  needs cutoff >= 3 (316k states) and an ensemble (per-sector BLAS) RDM.
- One classical periodic sampler (`ocp_monte_carlo`) serves H1a and H1b.
  The pair-law addition is the same function in the sampler and in the
  N = 2 quadrature, so the quadrature checks sampling and the Ewald part
  independently but not the Kelbg/Pauli definition itself (that is the
  Kelbg leaf oracle).
- WK reweighting of the contact bin is heavy-tailed in Monte Carlo. For
  N = 3, 4 the predictor rows are qualitative.
- fortplot: a line clipped by `ylim` in one subplot leaks segments into the
  next subplot. The demo clips the data at the axis top as a workaround.
