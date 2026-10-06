# Classical correlation hierarchy of the 3D OCP (H1a)

Status: candidate evidence. Every interval below is statistical (one
standard error from autocorrelation, batch means or a jackknife over
batches). No deterministic enclosure is claimed. Campaign item H1
(docs/hierarchy-campaign.md), classical part.

## Model and units

The model is the classical one-component plasma: electrons, a uniform
neutralizing background and the bare Coulomb interaction in a periodic cube
with tin-foil Ewald sums. Units are Wigner-Seitz: `a = 1`,
`e^2/(4 pi eps0) = 1`, `n = 3/(4 pi)`, `L = (4 pi N/3)^(1/3)` and
`beta = Gamma`, so `beta v(r) = Gamma/r`. Gamma and the Debye wavenumber
`kappa a = sqrt(3 Gamma)` come from the shared `kin6d_species` scales and the
generated `ch_debye_wavenumber` leaf. Couplings are
Gamma = 0.1, 0.3, 1, 3, 10, 30 (fluid; freezing near 175).

## Sampling

- `src/ocp_monte_carlo.f90` runs canonical Metropolis Monte Carlo with
  single-particle moves, uniform in `[-d, d]^3`. The energy is the Ewald sum
  of `kin6d_coulomb_leaves`. Its parameters come from `kin6d_ewald` with the
  new `real_cutoff = 0` (minimum-image real space; every omitted image is at
  least `L/2` away, and the existing shell bound covers it). The reciprocal
  sum uses the cube `0 < |n|_inf <= nk`. A move costs O(N) for the real-space
  difference plus an O(modes) update of `rho_k`; the energy change is exact
  for the truncated Hamiltonian. At tolerance 1e-5 the reported truncation
  bound is about 6e-9 per pair (units 1/L).
- Preparation: a uniform random start, then 2000 equilibration sweeps. The
  step is tuned towards 40% acceptance in the first half only. Production is
  10000 sweeps in 20 batches. Samples (g2, g3, S(k)) are taken every 5
  sweeps, total forces every 10. Diagnostics per chain are the
  equilibration trace, acceptance, integrated autocorrelation time `tau`
  (Sokal window), effective sample size and the drift of the incremental
  energy against a full recompute.
- Finite size: N = 1000 and N = 500 at every Gamma, plus a linear 1/N
  extrapolation of `u`.
- Independent method: molecular dynamics at Gamma = 1, N = 500, from
  `kin6d_ocp` (PME, encounter-resolved splitting) at a quarter of its
  default step. An isokinetic thermostat acts during equilibration only;
  production is NVE.

## Definitions

Normalized distributions `g2(r)` and `g3(r12, r13, r23)`, with
`h_ij = g2(r_ij) - 1`. The histograms are normalized by `N(N-1)/V^2` and
`N(N-1)(N-2)/V^3`. Because `S(0) = 0` for the OCP, the leading canonical
1/N correction is this normalization alone. Tables therefore report
`g2 (N-1)/N` and `g3 (N-1)(N-2)/N^2`. The raw `g2` approaches `1 + 1/N` at
large r. Without the correction, the g2-route energy at N = 200,
Gamma = 10 is off by 10%.

- Kirkwood product `K = g12 g13 g23`. Cumulant-truncation reference
  `P = 1 + h12 + h13 + h23`.
- Connected triplet cumulant (task definition):
  `g3 = 1 + sum h + sum h h + h12 h13 h23 + c3`, that is,
  `c3 = g3 - K = K delta_K`.
- Kirkwood superposition defect `delta_K = g3/K - 1`.
- Ursell (BBGKY density) cumulant `U3 = g3 - P`. It is the part dropped by a
  hierarchy truncated after `g_2`.
- Reducible pair-product part `K - P`.
- Bin averages: g3 is a bin average over the exact ideal measure
  `d3r2 d3r3 = 8 pi^2 r12 r13 r23 dr12 dr13 dr23` on the triangle domain.
  K and P are averaged over the same measure, so `c3` and `U3` contain no
  first-order binning bias. Bins have width 0.2a, and all sides are < 3a.
- Sizes:
  - `||h2||^2 = n int_{r<3a} h^2 d3r`;
  - `||X||^2 = n^2 int int_D X^2 d3r2 d3r3`, where D is the triplet domain
    (all sides < 3a), for `X = c3, U3, K - P`;
  - rms defect `[int K V delta_K^2 / int K V]^(1/2)` over bins with
    `K > 1e-3`.

  Squared sizes are jackknife bias-corrected, which removes the
  `O(1/samples)` noise bias. A size is called resolved when it exceeds twice
  its standard error; otherwise the figures show an upper limit.
- Bin-discretization indicator: `||K_bin - K_centre||`, the bin average
  minus the centre product, reported separately from the statistical error.

## YBG closures (influence of three-body correlation)

The first YBG (BBGKY equilibrium) equation for the OCP, with the background
force included, reads

`-dw/dR = beta F(R) + n int beta F(r13) . R-hat [g3/g2 - 1] d3r3`,
`g2 = exp(-w)`.

Three closures are compared:

- (a) Exact. The right side equals the conditional mean total force
  `<F_1 . R-hat | r12 = R>` in the periodic box. This is sampled directly
  from the Ewald forces, with no binning of g3. Bins with fewer than 20 pair
  samples are excluded.
- (b) Kirkwood (`c3 = 0`). The integral becomes
  `n int F(r13) g13 h23 d3r3`.
- (c) First order (`g3/g2 -> 1 + h13 + h23`). The integral becomes
  `n int F(r13) h23 d3r3`, the mean field of particle 2's screening cloud.
  This is the Debye/Poisson-Boltzmann level. It agrees with `U3 = 0` to
  first order in h. The strict `U3 = 0` form divides by `g12` and diverges
  where `g12 -> 0`.

The closures take the MC `g2` as input (one-shot, not self-consistent).
Before use, `h` is shifted on `r > L/4` so that the table satisfies perfect
screening, `n int h = -1`. A non-neutral hole produces a spurious far-field
force; the size of the shift is reported as `neutral_mc + 1`. The integrals
use a midpoint rule in s and 32-point Gauss-Legendre in `mu` (leaf
`ch_ybg_force_moment`). The integration starts at `w(L/2) = 0`. Outputs per
closure: `||g - g_MC||` (n-weighted L2 over r < L/2) and the max norm; the
excess energy `beta u = int 2 pi n r^2 beta v h dr` (leaf
`ch_excess_energy_density`) plus an HNC tail beyond `L/2`; and the
Stillinger-Lovett moment `n int r^2 h d3r kappa^2/6 -> -1`.

Reduced models:

- HNC (`src/ocp_hnc.f90`): Ornstein-Zernike with the Ng (1974) Ewald split
  of the Coulomb tail, a DST-I through the FortNum FFT and damped Picard
  iteration. Kernels are generated leaves.
- Debye-Hueckel `exp(-Gamma e^{-kappa r}/r)` and its linearization.

## Oracles (CTest `correlation_hierarchy_h1a_oracles`)

- The MC Ewald energy agrees with the tiled direct Ewald sum of
  `kin6d_nbody` (different alpha, real cutoff 2) to 2e-11, below the 2.6e-9
  bound. After sweeps, the incremental energy matches a full recompute to
  1e-10. MC forces match central differences to 6e-9 (relative max).
- The N = 2 Metropolis chain gives `<U>` equal to direct quadrature of the
  periodic pair Boltzmann factor within statistics.
- The triplet bin measure matches the analytic two-ball overlap integral to
  2.4e-5. For an ideal gas, `g3 = 1` within 4 standard errors, and K = P = 1.
- The first-order YBG closure with `h = A exp(-r^2)` reproduces the
  analytic Gaussian-cloud force to 2e-6.
- HNC reproduces the Caillol-Gilles values at Gamma = 0.1 and 1
  (-0.0256885, -0.570455) to 1e-6 and 1.3e-6, and the Ortner value at
  Gamma = 0.01 to 1e-8. The DST matches a direct sine sum.
- The AR(1) autocorrelation time is 4.500 against an exact 4.5. The
  JVP leaf of the literature-fit derivative matches finite differences, and
  the Debye leaves match their closed forms.

## Literature

- Potekhin & Chabrier, Phys. Rev. E 62, 8554 (2000): u_ex fit of Caillol
  (1999) and earlier MC data (generated leaf `ch_ocp_energy_fit`; its JVP is
  the derivative oracle).
- Caillol, J. Chem. Phys. 111, 6538 (1999); Caillol & Gilles (HNC);
  Demyanov & Levashov, arXiv:2509.02390 (2025): tabulated MC u_ex
  (Gamma = 0.1: -0.0256975(35); 1: -0.571414(24); 10: -7.998170(16)) and the
  HNC values used above.
- Ng, J. Chem. Phys. 61, 2680 (1974): HNC long-range renormalization.
