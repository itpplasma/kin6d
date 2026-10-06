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
  bound is about 6e-9 per pair (units 1/L). The same sampler serves the
  H1b finite-N cell (`mc_setup_cell`: general units, a real image cube for
  small N, and the `kin6d_pair_law` Kelbg/Pauli addition;
  docs/quantum-classical-correlations.md). The OCP path is unchanged
  bit for bit.
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
- Discretization indicators, reported separately from the statistical
  error:
  - triplet measure quadrature: sub = 32 midpoint points per side in
    r12, r13, with r23 clipped exactly to the triangle range. The per-bin
    relative error is 9e-4. With sub = 8, the measure error alone produced a
    spurious floor `||c3||^2 ~ 1e-3`.
  - resolution: `||c3||` from the same counts merged into 0.6a bins
    (`c3_sq_coarse`).
  - core sensitivity: sizes with all sides >= 0.4a (`*_outer`).

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

## Results (N = 1000, `benchmarks/correlation-map/h1a-classical.nml`)

Run time: 12 chains plus MD on 13 threads, 12.4 min wall. Data are in
`output/correlation-map-*.csv` and `output/correlation_map_classical_demo.nc`.
Raw accumulators are in `output/correlation-map-raw.bin`; `reuse_samples`
re-runs the analysis in about 7 s.

Sampling diagnostics: acceptance is 0.90/0.77/0.50/0.41/0.39/0.42 from
Gamma = 0.1 to 30; at weak coupling the step is capped at L/2. `tau` of
`beta U/N` is 0.6 to 4.9 sweeps, and the ESS is 8300 down to 1000. The
incremental-energy drift is <= 4e-12 per particle in `k_B T`.

| Gamma | u (MC) | u (fit) | du/dGamma (fluct.) | (fit) |
|---|---|---|---|---|
| 0.1 | -0.02620(4) | -0.02567 | -0.3687(18) | -0.3687 |
| 0.3 | -0.12010(9) | -0.11960 | -0.5414(20) | -0.5414 |
| 1 | -0.57202(15) | -0.57183 | -0.7090(18) | -0.7089 |
| 3 | -2.11360(27) | -2.11317 | -0.8068(19) | -0.8059 |
| 10 | -7.99946(55) | -7.99808 | -0.8594(15) | -0.8583 |
| 30 | -25.44407(95) | -25.44333 | -0.8785(7) | -0.8798 |

- Finite size, linear in 1/N from N = 500 and 1000:
  - Gamma = 0.1: -0.025689(92), against MC -0.0256975(35)
    (Demyanov-Levashov) and -0.025672 (fit); the N = 1000 value is still 2%
    low (`lambda_D = 1.8a`, `L = 16a`).
  - Gamma = 1: -0.57193(37), against -0.571414(24).
  - Gamma = 10: -7.9999(14), against -7.998170(16).
- The fluctuation derivative agrees with the fit derivative (JVP leaf) to
  within 1.8 standard errors at every Gamma.
- MD at Gamma = 1, N = 500 (NVE, 1000/omega_p): `u = -0.57265(49)` at
  `T/T_ref = 0.99846(40)`. Corrected to the effective Gamma, this matches
  MC N = 500 (-0.5732) within 1 sigma, including the temperature uncertainty. The g2 difference to MC
  (n-weighted L2, 0.017) equals the MC noise level (0.015). The energy drift
  is 1e-2 `k_B T` per particle over the run.
- `S(k_min = 0.39/a)`, MC / DH / HNC: Gamma 0.1: 0.331(5) / 0.336 / 0.337;
  Gamma 1: 0.0486(5) / 0.0482 / 0.0490; Gamma 10: 0.00518(8) / 0.00504 /
  0.00518.
- First peak of g2: Gamma 10 at r = 1.674(10)a with g = 1.139 (HNC 1.687,
  1.105); Gamma 30 at 1.668(3)a with g = 1.440 (HNC 1.657, 1.335).

Sizes (sqrt of the jackknife-corrected squared norms; triplet domain with
sides < 3a, 0.2a bins):

| Gamma | h2 | K - P | U3 | c3 | rms delta_K | c3/(K-P) |
|---|---|---|---|---|---|---|
| 0.1 | 0.146 | 0.054 | 0.052 | < 0.017 (2 sigma) | noise | - |
| 0.3 | 0.289 | 0.209 | 0.197 | 0.019(2) | 0.0016 | 0.09 |
| 1 | 0.513 | 0.616 | 0.584 | 0.061(1) | 0.0052 | 0.10 |
| 3 | 0.756 | 1.228 | 1.191 | 0.137(1) | 0.0131 | 0.11 |
| 10 | 1.079 | 2.221 | 2.236 | 0.361(1) | 0.0294 | 0.16 |
| 30 | 1.554 | 4.031 | 4.203 | 1.060(2) | 0.0640 | 0.26 |

- Effective exponents (`output/correlation-map-scaling.csv`). The quoted
  errors are statistical only and ignore curvature; three to five points
  each:
  - Gamma <= 1: `||h2|| ~ Gamma^0.50`, `||K-P|| ~ Gamma^0.96`,
    `||U3|| ~ Gamma^0.91`;
  - c3 over Gamma 0.3..3: `~ Gamma^0.74(1)`; over 0.3..30:
    `~ Gamma^0.89`.

  These are local slopes on a finite domain with a nonlinear core, not
  asymptotic Debye-limit exponents.
- The Ursell cumulant `U3` (what a `g2`-truncated hierarchy drops) is as
  large as the pair products at every Gamma. The superposition-connected
  part `c3` grows faster than them, but `||c3||/||K-P||` is only 0.26 at
  Gamma = 30. A power-law extrapolation of the ratio reaches 1 near
  Gamma ~ 600, beyond freezing. Inside the fluid range c3 never becomes
  comparable to the pair products.
- Indicators:
  - N = 500 gives the same `||c3||` within errors;
  - the core-excluded values are equal from Gamma = 1;
  - 0.6a bins lower `||c3||^2` by 33% at Gamma = 10 and 30 (smoothing),
    which suggests the 0.2a values are low by about 4% (h^2 scaling).

Influence on the pair level (one-shot YBG with MC g2 input):

| Gamma | L2 g error: exact / Kirkwood / 1st order / HNC / DH | du/u: exact / Kirkwood / 1st / HNC / DH |
|---|---|---|
| 0.1 | 0.007 / 0.009 / 0.009 / 0.006 / 0.007 | 0.04 / 0.07 / 0.08 / 0.02 / -0.05 |
| 0.3 | 0.007 / 0.008 / 0.014 / 0.007 / 0.019 | 0.006 / -0.02 / 0.03 / 0.005 / -0.18 |
| 1 | 0.007 / 0.021 / 0.081 / 0.007 / 0.106 | 0.006 / -0.04 / 0.17 / 0.003 / -0.51 |
| 3 | 0.006 / 0.079 / 0.38 / 0.015 / 0.35 | 0.010 / -0.05 / 0.52 / 0.005 / -1.1 |
| 10 | 0.006 / 0.27 / 2.8 / 0.080 / 0.87 | 0.000 / -0.32 / 1.8 / 0.008 / -2.4 |
| 30 | 0.011 / 3.3 / 296 / 0.26 / 1.5 | -0.03 / 5.9 / 190 / 0.008 / -4.6 |

The MC noise floor of the L2 error is 0.006 to 0.01.

- The exact-force YBG route reproduces the histogram g2 at the noise level
  for all Gamma. This checks the conditional-force estimator and the
  periodic YBG identity. Its energy errors at Gamma = 0.1 and 30 come from
  the sparsely sampled core.
- The first-order (mean-field) closure fails from Gamma ~ 1. It
  over-counts the ion-sphere attraction, giving a spurious contact peak
  (g = 2.9 at Gamma = 10).
- Kirkwood is good to Gamma ~ 1 (2% in g, 4% in u). Error grows like the
  superposition-connected force `n int F(r13) c3/g12`; divided by small g12
  and exponentiated, it reaches 30% in u at Gamma = 10 and diverges at
  Gamma = 30. Three-body correlation therefore controls the pair structure
  through the YBG route from Gamma of a few on. The self-consistent closed
  HNC stays at 1% in u up to Gamma = 30.
- The Stillinger-Lovett moment is noise-dominated (r^4 weight) beyond
  Gamma = 1 and is not interpreted.

## Derivatives

`du/dGamma` uses the canonical fluctuation identity
`d<beta U>/dGamma = <U> - Gamma Var(U)`, an equation-level derivative of
the sampled measure (leaf `ch_energy_fluctuation_derivative`). The oracle is
the generated JVP of the literature fit. All generated leaves carry
JVP/VJP. No derivative of the YBG/HNC outputs is propagated yet.

## Open

- Self-consistent BGY-Kirkwood, and closures beyond the first YBG
  equation (the second YBG equation for g3), are not implemented.
- The quoted c3 is a 0.2a bin average (resolution indicator above). A
  finer, statistics-limited g3 near contact needs more sampling.
- Time-series statistics (`kin6d_series_statistics`) are generic and could
  move to FortNum. `ocp_shielding_demo` and the H1b batch means use them.
- FortSym emits `sqrt(1*3)` (an integer argument) for `sqrt(3*one)` with a
  rational one. The generator works around it with `sqrt(rat(3))`; this is a
  candidate FortSym folding fix.
- MD at the default OCP step fails closed (cluster check) at Gamma = 1
  within about 1e4 steps; the demo uses dt/4. The C3/C4 defaults target
  `Lambda >> 1`.

## Oracles (CTest `correlation_hierarchy_h1a_oracles`)

- The MC Ewald energy agrees with the tiled direct Ewald sum of
  `kin6d_nbody` (different alpha, real cutoff 2) to 2e-11, below the 2.6e-9
  bound. After sweeps, the incremental energy matches a full recompute to
  1e-10. MC forces match central differences to 6e-9 (relative max).
- The N = 2 Metropolis chain gives `<U>` equal to direct quadrature of the
  periodic pair Boltzmann factor within statistics.
- The triplet bin measure matches the analytic two-ball overlap integral to
  2.4e-5 in total, and to 9e-4 per bin (sub 32 against 96). For an ideal gas, `g3 = 1` within 4 standard errors, and K = P = 1.
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
