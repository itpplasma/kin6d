# Periodic microscopic OCP: plasma oscillation (C3) and Debye shielding (C4)

Status: candidate evidence. Statistical intervals are labelled statistical;
no deterministic trajectory or observable certificate is claimed.

## Model

Classical electron one-component plasma in a periodic cube, 3D3V, with a
uniform neutralizing background and the bare Coulomb interaction. No
softening, screened pair potential, collision operator or wall is imposed.
Plasma units of the responsive species (electrons only): `k_B T = m = 1`,
`lambda_De = 1`, `omega_pe = 1`, so the density is `Lambda = n lambda_D^3`
and the Coulomb constant is `1/(4 pi Lambda)`; `kin6d_ocp` checks these
against the shared `kin6d_species` scales. C4 adds a fixed repulsive point
impurity of charge `Q = -kappa 4 pi Lambda e` (same sign as the electrons),
compensated in the background. Its nonlinear radius is
`r_nl = |qQ|/(4 pi eps0 k_B T) = kappa lambda_D`. "Weak" means linear response
outside a small radius: `kappa = 0.05` (linearity gate `kappa = 0.1`). At
`Lambda = 100`, `|Q| = 63 e`. The impurity is an external source, not a mobile
ion.

## Numerics

- `src/pme.f90`: particle-mesh Ewald (PME). Real space uses a linked cell
  list with cutoff `< L/2`. The mesh part uses cardinal B-splines of order
  4/6/8 and an energy-optimal influence function over aliases `|M| <= 2`.
  Transforms are FortNum 1D FFTs applied along x, y and z. Forces are the
  exact gradient of the mesh energy, so the dynamics is Hamiltonian.
  `pme_select` takes a requested rms force tolerance relative to
  `F0 = |c| e^2/a^2`. Each component is estimated: the Kolafa-Perram real-space
  estimate and the exact mean-square reciprocal pair-force functional of the
  used influence function. Each must be `<= tol/(2 sqrt 2)`; the factor 2
  covers the mesh self force, which is not estimated. The cost model picks
  among the admissible (order, cutoff, mesh) choices. The direct Ewald sum
  of `kin6d_nbody` is the oracle.
- `src/near_encounter.f90`: encounter-resolved near/far splitting. Fixed-step
  KDK with bare Coulomb heats an OCP at `Lambda = 100` by about 20% in
  `10/omega_p` (`dt = 0.05` and `0.01`; measured). This is the reason for
  the splitting. The pair interaction is split exactly with a C2 switch
  (`near_switch_pair`, generated). One step is a FortNum Strang composition:
  far kick, near-pair flow, far kick. In the near-pair flow, clusters of
  candidate pairs take their own KDK substeps, sized by the two-body
  periapsis Kepler, crossing and switch-crossing times. All other particles
  drift exactly. Cluster completeness is verified after the flow and fails
  closed. The tangent uses the same flow with the substep structure frozen.
- `src/vlasov_response.f90` (independent oracles):
  - the plasma dispersion function from the entire Dawson series in
    quadruple precision;
  - the Landau root of the Maxwellian dielectric, found by Newton iteration;
  - the full finite-k linear Vlasov initial-value response of a density
    perturbation, solved as a Volterra equation (trapezoid rule plus
    Richardson extrapolation).
- `src/pair_correlation.f90`: estimators for `g(r)` (cell list) and for the
  shell-averaged `S(k)`.
- `kin6d_ocp`: preparation, mode amplitude, displacement and its tangent,
  and the radial histogram. It also holds the periodic screened-Poisson Green
  function. That function is the image sum of `exp(-r)/r` minus the zero mode
  `4 pi/(V kappa^2)`, with a shell truncation bound. Its leaves are generated
  in `gen/gen_ocp.f90`.

## Preparation

Electron positions are uniform; with an impurity, a draw is accepted with the
bare-impurity Boltzmann factor `exp(-kappa/r)`. Velocities are Maxwellian at
T = 1 with zero total momentum. Equilibration runs for `5/omega_p` under an
isokinetic velocity-rescaling thermostat. This is a different, declared
dynamics and is switched off before any measurement. Every member records its
seed, rejected draws and cost.

## Cases and gates

- C3 (`app/ocp_oscillation_demo.f90`, configs `benchmarks/ocp/*.nml`). The
  displacement `x_a += (eps/k) sin(k x_a)` is applied on all three axes,
  with `k = 2 pi/L` and `eps = 0.01`. From the same equilibrated state, a
  paired unperturbed run (common random numbers) gives
  `R(t) = (2/N) Re(rho_k^pert - rho_k^unpert)/(-eps)`, averaged over the axes.
  The step tangent `dR/d eps` at `eps = 0` uses the same runs. Comparisons:
  - the linear Vlasov response;
  - the Landau root mode;
  - the Bohm-Gross frequency;
  - the cold slab, `kin6d_slab_oscillation`.

  Two boxes:
  - "clean": `k lambda_D = 0.3`, `Lambda = 10`, N = 91 855;
  - "Landau": `k lambda_D = 0.5`, `Lambda = 20`, N = 39 478.

  The structure factor `S(k)` and `g(r)` come from unperturbed snapshots.
  Binary snapshots (`N, t, Lambda, L, x(3N), v(3N)`) are written to
  `output/` for C5/H1.
- C4 (`app/ocp_shielding_demo.f90`): `L = 4 lambda_D`. The estimator is the
  time- and ensemble-averaged shell density. The reference is `-kappa G(r)`
  with G the periodic Green function, averaged over each bin. The isolated
  Yukawa curve is an overlay. The fit window `[0.3, 1.8] lambda_D` excludes
  `r_nl`. Separate gates:
  - mean response: fitted screening length per member and for the mean
    profile;
  - fluctuations: `S(2 pi/L)` against Debye-Hueckel;
  - geometry: periodic minus isolated reference;
  - estimator: integrated autocorrelation time and effective samples.

  Numerical settings (PME estimate, step, energy defect) are reported
  separately.

## Results

See the log entry in `docs/hierarchy-campaign.md` and `output/*.csv`.

## Derivatives and limits

The step tangent is checked against central differences at the level of a
whole run (`test_ocp`, `test_pme`). The microscopic dynamics is chaotic. The
position-tangent norm grows at a measured Lyapunov-type rate of order
`omega_p`, so `dR/d eps` from a single trajectory is meaningful only for a
few plasma periods. Beyond that it is dominated by exponential divergence;
the paired finite-amplitude estimator remains usable. Ensemble-mean
sensitivities need a separate, non-trajectory method (for example a
linear-response or shadowing estimator) and are not claimed. Encounter
substep counts and cluster structure are held fixed in the derivative, so it
is piecewise smooth.

## Open

- No VJP of the split step.
- Quiet-start positions are not implemented.
- The mesh self force is measured, not estimated.
- The Debye fit uses a fixed amplitude `kappa`; the amplitude is not fitted.
