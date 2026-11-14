# Teaching discharge: operating limits and disruptions

Status: reduced teaching model. It runs natively and in the browser
(WebAssembly through LFortran and emscripten, `web/tokamak-limits/build.sh`)
and gives the same results in both. It is not a validated disruption
predictor and makes no claim beyond the checks listed below.

## Purpose

Students change the machine (major and minor radius, elongation,
triangularity, field, wall) and the plasma (current, central q, density,
temperature, impurity fraction), then ramp density, heating or current and
watch which limit the plasma reaches first and how it disrupts in time.

## Model

`src/continuum/tokamak_discharge.f90`, module `milca_tokamak_discharge`.

- Geometry: a straight cylinder of area-equivalent radius
  `a_e = a sqrt(kappa)`; q carries shaping and toroidicity through the ITER
  Physics Basis q95 factor, so q at the edge equals q95.
- State: temperature `T(x,t)` (electrons and ions equal) and enclosed current
  `I(x,t)` on 65 nodes, `x = r/a_e`.
- Heat: `3n dT/dt = div(2 n chi grad T) + eta j^2 + p_aux - n^2 L(T)`, backward
  Euler. `chi = chi0 (1 + 3x^2)` with `chi0` calibrated at the start and then
  scaled with the IPB98(y,2) confinement time, so confinement degrades with
  power as in the scaling. Radiation: bremsstrahlung and a carbon-like cooling
  curve (surrogate, peak 1.5e-31 W m^3 at 10 eV).
- Current: resistive diffusion with Spitzer resistivity. The edge holds `I_p`
  (feedback) until the thermal quench and then obeys
  `L_ext dI_p/dt = -2 pi R0 eta j`.
- Tearing: `r_s Delta'` of the 2/1 and 3/2 surfaces from the cylindrical outer
  equation (`src/mhd/cylinder_stability.f90`); Rutherford growth
  `dw/dt = 1.22 eta/mu0 Delta'` with a saturation surrogate
  `w_sat = r_s rDelta'/(rDelta' + 12)`; islands flatten T across their width.
- Ideal kink: cylindrical low-beta energy with a wall (Newcomb shooting),
  m = 1..4, n = 1; growth estimate `gamma tau_A = sqrt(-dW/2K)`. Only m = 2
  with `dW < -2e-3` (a tenth of the depth of its window) disrupts. The
  cylinder also opens windows m - 1 < q95 < m for m >= 3 once the current
  profile broadens; toroidal coupling closes them in real tokamaks, whose
  low-n external-kink limit is q95 = 2, so they are reported, not acted on.
- Ballooning: local s-alpha first boundary (`src/mhd/ballooning_s_alpha.f90`)
  tabulated for s = 0.1..3. Shaping enters through the surrogate factor
  `(1 + kappa^2)/2` on alpha. A disruption needs the core (`x <= 0.7`) with
  shear `s >= 0.4` beyond the first boundary, confirmed by the direct Newcomb
  test at the local `(s, alpha)` (second stability). The pedestal region and the
  low-shear core are reported only: the local high-n model is pessimistic
  there (finite Larmor radius, diamagnetic and sawtooth effects).
- Triggers of the thermal quench: island overlap (Chirikov) or an island
  reaching the edge; an ideal mode growing to 0.1 a; the radiating front
  reaching x = 0.8 (radiative collapse, the density limit of this model).
- Disruption: stochastic transport (chi x 3000), wall impurity influx
  (surrogate, carbon fraction raised to 3 %), current flattening a third of
  the way to uniform with poloidal-flux conservation (current spike), current
  quench, and a rigid vertical displacement once position control is lost,
  `gamma_v = (kappa - 1)/((2.3 - kappa) tau_w)`. Contact with the wall scrapes
  the plasma; the current outside the remaining radius is lost.
- Sawteeth: q0 < 0.97 flattens T and q inside 1.4 x(q=1) (surrogate).

Greenwald density and the Troyon beta_N = 2.8 appear only as empirical
reference lines; they never trigger a disruption.

## Checks

- `test/test_cylinder_stability.f90`: uniform-current kink energy against the
  closed form, including the window m - 1 < q_a < m and wall stabilisation;
  `r_s Delta'` of a smoothed uniform current channel against the step-current
  closed form (agreement 0.6 %); convergence in eps and in the grid.
- `test/test_tokamak_discharge.f90`: ITER-like operating point (q95 = 3.13
  from the shaping formula, Greenwald density 1.19e20 m^-3); density, power and
  current ramps and a safe ramp. Native gfortran and the WebAssembly build
  (`web/tokamak-limits/check.mjs`) give the same quench times and causes.

## Browser interface

`milca_td_sizes`, `milca_td_limits(par, out, prof, status)` and
`milca_td_run(par, ramp, maxrec, nrec, series, profiles, status)` are
`bind(C)`. `par` (20): R0, a, kappa, delta, B0, I_p [MA], q0, nbar [1e20],
density peaking, T0 [keV], temperature peaking, Zeff, impurity fraction, b/a,
wall time [s], H98, ion mass, T_edge [keV]. `ramp` (8): duration, final nbar,
final P_aux [MW] (negative keeps it), final I_p, end time, settling time.
`series` (28 per record): t, I_p, T0, W, nbar, n/n_G, beta_N, q0, q95, li, Z,
remaining minor radius, w_21, w_32, x_21, x_32, P_rad, P_ohm, P_aux, V_loop,
phase, cause, core ballooning alpha/alpha_c (the trigger quantity), ideal amplitude, its m, r_s Delta' (2/1, 3/2), dt.
`profiles` per record: T, n, q, j, alpha/alpha_c on the 65 nodes.

## Sources

Wesson, Tokamaks (4th ed.), sections 6.5 (kink), 6.8 (tearing), 7 (disruptions);
Furth, Rutherford and Selberg, Phys. Fluids 16, 1054 (1973); Rutherford, Phys.
Fluids 16, 1903 (1973); Connor, Hastie and Taylor, PRL 40, 396 (1978);
Greenwald, PPCF 44, R27 (2002); Troyon et al., PPCF 26, 209 (1984);
ITER Physics Basis, Nucl. Fusion 39 (1999), chapters 1 (q95, IPB98(y,2)) and 3
(disruptions). No non-public code was used.
