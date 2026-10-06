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

- Classical parent (`src/classical_cell.f90`): N point electrons with the
  shared Ewald energy (`kin6d_nbody`/`kin6d_ewald`). The energy is pairwise
  additive with the zero-mean Ewald pair potential `u` (Fourier `Q = 0`
  removed, as in the quantum Hamiltonian). Spin does not enter the classical
  Boltzmann weight; spin labels are spectators with probability 1/2 that a
  pair has equal spins.
  - N = 2: the centre of mass integrates out. `p(d) = exp(-beta u(d))/Z` on
    the torus is integrated by the periodic midpoint rule (Z, Fourier
    moments) and by spherical Gauss rules (shell bins). This is deterministic.
  - N = 3, 4: canonical Metropolis Monte Carlo with single-particle moves.
    Intervals are batch-means statistical errors (40 batches). They are
    labelled statistical, never certified.
- Exact quantum N = 2 (`src/quantum_pair_cell.f90`): the centre of mass and
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
